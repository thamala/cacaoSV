/*
 Copyright (C) 2021 Tuomas Hamala

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 For any other inquiries send an email to tuomas.hamala@gmail.com
 
 ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

 Program for compiling RNA- and DNA-seq counts for heterozygous SNPs found in heterozygous SVs

 Compiling: gcc SV2ase.c -o SV2ase -lm

 Usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -snp1 [file] A file with SNPs from haplotype 1 (produced by MUMmmer's show-snps tool)
 -snp2 [file] A file with SNPs from haplotype 2 (produced by MUMmmer's show-snps tool)
 -dna [file] A VCF file for DNA-seq based SNPs
 -rna [file] An mpile file for RNA-seq based SNPs
 -ind [int] Number defining which individual from the VCF files is analyzed
 -t [str] Either INS, DEL, DUP, TRA, or INV.
 -l [int]-[int] Minimum and maximum length of the SVs.

 Example:
 ./SV2ase -sv SV.vcf -snp1 haplo1.txt -snp2 haplo2.txt -dna SNPs.vcf -rna SNPs.mpile -ind 0 -t INV -l 1e4-1e6 > out.txt
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n\n"

typedef struct{
    int start, stop, ma;
    char chr[101];
}SV_s;

typedef struct{
    int pos, ma;
    char ref, alt, chr[101];
}Site_s;

typedef struct{
    int pos, maj_i, min_i;
    char maj_a, min_a, chr[101];
}Dna_s;

void openFiles(int argc, char *argv[]);
SV_s *readSV(FILE *sv_file, char **id, int ind, int *n, int *l, int type, double len[]);
Site_s *readSNP(FILE *snp_file, SV_s *svl, int sv_n, int sv_l, int hap, int *n);
Dna_s *readDNA(FILE *dna_file, Site_s *snp1, Site_s *snp2, int ind, int snp1_n, int snp2_n, int *n);
void readRNA(FILE *rna_file, Dna_s *dna, char id[], int ind, int dna_n);
int *splitInfo(char str[]);
int isNumeric(const char *s);
void lineTerminator(char *line);

int main(int argc, char *argv[]){
    
    int second=0, minute=0, hour=0;
    time_t timer=0;
    
    if(argc == 1)
        exit(EXIT_FAILURE);
    
    timer = time(NULL);
    
    openFiles(argc, argv);
    
    second = time(NULL) - timer;
    
    minute = second / 60;
    
    hour = second / 3600;
    
    fprintf(stderr,"\nDone!");
    if(hour > 0)
        fprintf(stderr,"\nElapsed time: %i h, %i min & %i sec\n\n", hour, minute-hour*60, second-minute*60);
    else if(minute > 0)
        fprintf(stderr,"\nElapset time: %i min & %i sec\n\n", minute, second-minute*60);
    else if(second > 5)
        fprintf(stderr,"\nElapsed time: %i sec\n\n", second);
    else
        fprintf(stderr,"\n\n");
    
    return 0;
}

void openFiles(int argc, char *argv[]){
    
    int i, j, ind=-1, type=-1, sv_n=0, sv_l=0, snp1_n=0, snp2_n=0, gene_n=0, dna_n=0;
    double len[2]={0};
    char *temp1=NULL, *temp2=NULL, *id=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    Site_s *snp1=NULL, *snp2=NULL;
    SV_s *svl=NULL;
    Dna_s *dna;
    FILE *sv_file=NULL, *snp1_file=NULL, *snp2_file=NULL, *rna_file=NULL, *dna_file=NULL;

    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((sv_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sv %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-snp1") == 0){
            if((snp1_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-snp1 %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-snp2") == 0){
            if((snp2_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-snp2 %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-rna") == 0){
            if((rna_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-rna %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-dna") == 0){
            if((dna_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-dna %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-ind") == 0){
            if(isNumeric(argv[++i]))
                ind = atoi(argv[i]);
            fprintf(stderr, "\t-ind %i\n", ind);
        }
        
        else if(strcmp(argv[i], "-t") == 0){
            argv[++i];
            for(j=0;j<5;j++){
                if(strcmp(argv[i], vars[j]) == 0)
                    type = j;
            }
            if(type == -1){
                fprintf(stderr,"\nERROR: Unknown parameter for -t\n");
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-t %s\n", vars[type]);
        }
        
        else if(strcmp(argv[i], "-l") == 0){
            temp1 = argv[++i];
            temp2 = strtok(temp1, "-");
            if(isNumeric(temp2))
                len[0] = atof(temp2);
            temp2 = strtok(NULL, "-");
            if(isNumeric(temp2))
                len[1] = atof(temp2);
            if(len[1] != 0){
                fprintf(stderr, "\t-l ");
                if(len[0] < 1e5)
                    fprintf(stderr,"%.0f-", len[0]);
                else
                    fprintf(stderr,"%.0e-", len[0]);
                if(len[1] < 1e5)
                    fprintf(stderr,"%.0f\n", len[1]);
                else
                    fprintf(stderr,"%.0e\n", len[1]);
            }
            else{
                fprintf(stderr,"\nERROR: Unknown parameter for -l\n");
                exit(EXIT_FAILURE);
            }
        }
        
        else{
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    
    if(sv_file==NULL || snp1_file==NULL || snp2_file==NULL || dna_file==NULL || rna_file==NULL || ind == -1){
        fprintf(stderr,"\nERROR: -sv [file] -snp1 [file] -snp2 [file] -dna [file] -rna [file] -ind [int] are needed\n");
        exit(EXIT_FAILURE);
    }
    
    while(1){
        svl = readSV(sv_file, &id, ind, &sv_n, &sv_l, type, len);
        if(sv_n == 0)
            break;
        snp1 = readSNP(snp1_file, svl, sv_n, sv_l, 1, &snp1_n);
        snp2 = readSNP(snp2_file, svl, sv_n, sv_l, 2, &snp2_n);
        if(snp1_n == 0 && snp2_n == 0){
            fprintf(stderr,"No heterozygous SNPs found\n");
            break;
        }
        dna = readDNA(dna_file, snp1, snp2, ind, snp1_n, snp2_n, &dna_n);
        if(dna_n == 0)
            break;
        readRNA(rna_file, dna, id, ind, dna_n);
        break;
    }
    
    if(id != NULL)
        free(id);
}

SV_s *readSV(FILE *sv_file, char **id, int ind, int *n, int *l, int type, double len[]){
    
    int i, j, k, char_i=0, maxchar=0, line_i=0, geno=9, *info;
    double p=0;
    char c, *temp=NULL, *line=NULL, *end=NULL, *id_cpy=NULL;
    SV_s *svl;
    
    while((c=fgetc(sv_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(sv_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((svl=malloc(line_i*sizeof(SV_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((id_cpy=malloc(maxchar*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while(fgets(line, maxchar+1, sv_file) != NULL){
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        p = 0;
        geno = 9;
        if(strcmp(temp,"#CHROM") == 0){
            while(temp != NULL){
                if(i > 9){
                    if(j == ind){
                        strcpy(id_cpy, temp);
                        *id = id_cpy;
                        break;
                    }
                    j++;
                }
                temp = strtok_r(NULL,"\t",&end);
                i++;
            }
            continue;
        }
        else if(temp[0] == '#')
            continue;
        while(temp != NULL){
            if(i == 1)
                strncpy(svl[*n].chr, temp, 100);
            else if(i == 2)
                svl[*n].start = atoi(temp);
            else if(i == 8){
                info = splitInfo(temp);
                if(type != -1 && info[0] != type)
                    break;
                if(len[1] > 0 && (info[2] < len[0] || info[2] > len[1]))
                    break;
                svl[*n].stop = info[1];
            }
            else if(i > 9){
                if(temp[0] == '0' && temp[2] == '0')
                    geno = 0;
                else if(temp[0] == '0' && temp[2] == '1')
                    geno = 1;
                else if(temp[0] == '1' && temp[2] == '0')
                    geno = 1;
                else if(temp[0] == '1' && temp[2] == '1')
                    geno = 2;
                else
                    geno = 9;
                if(j == ind){
                    if(geno != 1)
                        break;
                    else{
                        if(temp[0] == '0' && temp[2] == '1')
                            svl[*n].ma = 1; //first haplo carrying ref allele
                        else
                            svl[*n].ma = 2; //second haplo carrying ref allele
                    }
                }
                if(geno != 9)
                    p += geno;
                j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
            if(temp == NULL){
                if(p/((double)j*2)> 0.5){ //alt allele is maj
                    if(svl[*n].ma == 1)
                        svl[*n].ma = 2; //second haplo carrying maj allele
                    else
                        svl[*n].ma = 1; //first haplo carrying maj allele
                }
                *l = *l + info[2];
                *n = *n + 1;
            }
        }
    }
    
    fprintf(stderr,"\nFound %i heterozygous SVs\n", *n);
    
    free(line);
    
    return svl;
}

Site_s *readSNP(FILE *snp_file, SV_s *svl, int sv_n, int sv_l, int hap, int *n){
    
    int i=0, j=0;
    char *line=NULL, *temp=NULL;
    Site_s *snps;
    size_t len=0;
    ssize_t read;
    
    if((snps=malloc(sv_l*sizeof(Site_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, snp_file)) != -1){
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        temp = strtok(line,"\t");
        i = 0;
        while(temp != NULL){
            if(i == 0)
                snps[*n].pos = atoi(temp);
            else if(i == 1)
                snps[*n].ref = temp[0];
            else if(i == 2)
                snps[*n].alt = temp[0];
            else if(i == 8)
                strncpy(snps[*n].chr,temp,100);
            temp = strtok(NULL,"\t");
            i++;
        }
        while(j < sv_n){
            if(strcmp(snps[*n].chr, svl[j].chr) == 0){
                if(snps[*n].pos <= svl[j].stop && snps[*n].pos >= svl[j].start){
                    if(hap == 1){
                        if(svl[j].ma == 1) // haplo 1 carries major
                            snps[*n].ma = 1; //snp is from the major allele
                        else
                            snps[*n].ma = 0; //snp is from the minor allele
                    }
                    else{
                        if(svl[j].ma == 2) //haplo 2 carries major
                            snps[*n].ma = 1; //snp is from the major allele
                        else
                            snps[*n].ma = 0; //snp is from the minor allele
                    }
                    *n = *n + 1;
                    break;
                }
                else if(snps[*n].pos < svl[j].start)
                    break;
            }
            else if(strcmp(snps[*n].chr, svl[j].chr) < 0)
                break;
            j++;
        }
    }
    
    free(line);
    fclose(snp_file);
    
    return snps;
}

Dna_s *readDNA(FILE *dna_file, Site_s *snp1, Site_s *snp2, int ind, int snp1_n, int snp2_n, int *n){
    
    int i, j, snp1_i=0, snp2_i=0, ok=0, maj=-1, ref_i=0, alt_i=0;
    char ref, alt, *line=NULL, *temp=NULL, *temp2=NULL, *end=NULL, *end2=NULL, *end3=NULL;
    Dna_s *dna;
    size_t len=0;
    ssize_t read;
    
    if(snp1_n > snp2_n){
        if((dna=malloc(snp1_n*sizeof(Dna_s))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    else{
        if((dna=malloc(snp2_n*sizeof(Dna_s))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while((read = getline(&line, &len, dna_file)) != -1){
        if(line[0] == '\n' || line[0] == '#')
            continue;
        lineTerminator(line);
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        maj = -1;
        ok = 0;
        while(temp != NULL){
            if(i == 1)
                strncpy(dna[*n].chr, temp, 100);
            else if(i == 2)
                dna[*n].pos = atoi(temp);
            else if(i == 4)
                ref = temp[0];
            else if(i == 5){
                alt = temp[0];
                while(snp1_i < snp1_n){
                    if(strcmp(dna[*n].chr, snp1[snp1_i].chr) == 0){
                        if(dna[*n].pos == snp1[snp1_i].pos){
                            if(ref != snp1[snp1_i].ref || alt != snp1[snp1_i].alt){
                                fprintf(stderr,"Warning: Alleles differ at position %s %i\n", dna[*n].chr, dna[*n].pos);
                                break;
                            }
                            if(snp1[snp1_i].ma == 1){ //haplo1 carries maj allele
                                dna[*n].maj_a = snp1[snp1_i].alt;
                                dna[*n].min_a = snp1[snp1_i].ref;
                                maj = 1;
                            }
                            else{ //haplo1 carries min allele
                                dna[*n].maj_a = snp1[snp1_i].ref;
                                dna[*n].min_a = snp1[snp1_i].alt;
                                maj = 0;
                            }
                            ok = 1;
                            break;
                        }
                        else if(dna[*n].pos < snp1[snp1_i].pos)
                            break;
                    }
                    else if(strcmp(dna[*n].chr, snp1[snp1_i].chr) < 0)
                        break;
                    snp1_i++;
                }
                while(snp2_i < snp2_n){
                    if(strcmp(dna[*n].chr, snp2[snp2_i].chr) == 0){
                        if(dna[*n].pos == snp2[snp2_i].pos){
                            if(ok == 1){
                                ok = 0;
                                break;
                            }
                            if(ref != snp2[snp2_i].ref || alt != snp2[snp2_i].alt){
                                fprintf(stderr,"Warning: Alleles differ at position %s %i\n", dna[*n].chr, dna[*n].pos);
                                ok = 0;
                                break;
                            }
                            if(snp2[snp2_i].ma == 1){ //haplo2 carries maj allele
                                dna[*n].maj_a = snp2[snp2_i].alt;
                                dna[*n].min_a = snp2[snp2_i].ref;
                                maj = 1;
                            }
                            else{ //haplo2 carries min allele
                                dna[*n].maj_a = snp2[snp2_i].ref;
                                dna[*n].min_a = snp2[snp2_i].alt;
                                maj = 0;
                            }
                            ok = 1;
                            break;
                        }
                        else if(dna[*n].pos < snp2[snp2_i].pos)
                            break;
                    }
                    else if(strcmp(dna[*n].chr, snp2[snp2_i].chr) < 0)
                        break;
                    snp2_i++;
                }
                if(ok == 0)
                    break;
            }
            else if(i > 9){
                if(j == ind){
                    if(temp[0] != '0' || temp[2] != '1')
                        break;
                    temp2 = strtok_r(temp,":",&end2);
                    temp2 = strtok_r(NULL,":",&end2);
                    temp2 = strtok_r(NULL,":",&end2);
                    temp2 = strtok_r(NULL,":",&end2);
                    ref_i = atoi(strtok_r(temp2,",",&end3));
                    alt_i = atoi(strtok_r(NULL,",",&end3));
                    if(ref_i == 0 || alt_i == 0)
                        break;
                    if(maj == 0){
                        dna[*n].maj_i = ref_i;
                        dna[*n].min_i = alt_i;
                    }
                    else{
                        dna[*n].maj_i = alt_i;
                        dna[*n].min_i = ref_i;
                    }
                    *n = *n + 1;
                    break;
                }
                else
                    j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
        }
    }
    if(*n > 0)
        fprintf(stderr,"Found %i heterozygous SNPs\n", *n);
    else
        fprintf(stderr,"No heterozygous SNPs found\n");
    
    free(line);
    free(snp1);
    free(snp2);
    fclose(dna_file);
    
    return(dna);
}

void readRNA(FILE *rna_file, Dna_s *dna, char id[], int ind, int dna_n){
    
    int i=0, j=0, k=0, pos=0, alt_i=0, ref_i=0, maj=0, min=0, ok=0;
    char ref, chr[101], *line=NULL, *temp=NULL, *temp2=NULL, *end=NULL, *end2=NULL, *end3=NULL;
    size_t len=0;
    ssize_t read;
    
    while((read = getline(&line, &len, rna_file)) != -1){
        if(line[0] == '\n' || line[0] == '#')
            continue;
        lineTerminator(line);
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        ok = 0;
        while(temp != NULL){
            if(i == 1)
                strncpy(chr, temp, 100);
            else if(i == 2){
                pos = atoi(temp);
                while(k < dna_n){
                    if(strcmp(chr, dna[k].chr) == 0){
                        if(pos == dna[k].pos){
                            ok = 1;
                            break;
                        }
                        else if(pos < dna[k].pos)
                            break;
                    }
                    else if(strcmp(chr, dna[k].chr) < 0)
                        break;
                    k++;
                }
                if(ok == 0)
                    break;
            }
            else if(i == 4)
                ref = temp[0];
            else if(i > 9){
                if(j == ind){
                    temp2 = strtok_r(temp,":",&end2);
                    temp2 = strtok_r(NULL,":",&end2);
                    temp2 = strtok_r(NULL,":",&end2);
                    ref_i = atoi(strtok_r(temp2,",",&end3));
                    alt_i = atoi(strtok_r(NULL,",",&end3));
                    if(ref == dna[k].maj_a){
                        maj = ref_i;
                        min = alt_i;
                    }
                    else if(ref == dna[k].min_a){
                        maj = alt_i;
                        min = ref_i;
                    }
                    else
                        break;
                    if(maj > 0 || min > 0)
                        printf("%s %s %i %i %i %i %i\n", id, chr, pos, dna[k].maj_i, dna[k].min_i, maj, min);
                    break;
                }
                else
                    j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
        }
    }
    
    free(id);
    free(line);
    free(dna);
    fclose(rna_file);
}


int *splitInfo(char str[]){
    
    int i;
    static int info[3];
    char *temp1, *temp2, *end;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    
    temp1 = strtok_r(str,";",&end);
    temp2 = strtok(temp1,"=");
    temp2 = strtok(NULL,"=");
    for(i=0;i<5;i++){
        if(strcmp(temp2, vars[i]) == 0)
            info[0] = i;
    }
    temp1 = strtok_r(NULL,";",&end);
    temp2 = strtok(temp1,"=");
    temp2 = strtok(NULL,"=");
    info[1] = atoi(temp2);
    temp1 = strtok_r(NULL,";",&end);
    temp2 = strtok(temp1,"=");
    temp2 = strtok(NULL,"=");
    info[2] = atoi(temp2);
    
    return info;
}

int isNumeric(const char *s){
    
    char *p;
    
    if(s == NULL || *s == '\0' || isspace(*s))
        return 0;
    strtod(s, &p);
    return *p == '\0';
}

void lineTerminator(char *line){
    
    int i;
    
    for(i=0;line[i]!=0;i++){
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}
