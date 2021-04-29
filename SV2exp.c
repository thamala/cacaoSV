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

 Program for compiling RNA- and DNA-seq counts for genes within 5 Kb of SVs.

 Compiling: gcc sv2exp.c -o sv2exp -lm

 Usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -rna [file] A tab-delimited file with rna-seq read counts for each gene.
 -dna [file] A tab-delimited file with dna-seq read count for each gene.
 -genes [file] A tab-delemited file listing locations for each gene: chr start end strand (+ or -)
 -maf [double] Minor allele frequency cutoff for the SV.
 -t [str] Either INS, DEL, DUP, TRA, or INV
 -l [int]-[int] Minimum and maximum length of the SVs.

 Example:
 ./SV2der -sv SV.vcf -rna RNA.txt -dna DNA.txt -genes genes.txt -maf 0.15 -t DEL -l 500-1e5 > out.txt
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
    int type, start, stop, len, ma, gene_i, dist[1000], *geno;
    char chr[20], id[1000][50], **ind;
}SV_s;

typedef struct{
    int start, stop, len, ph;
    char chr[20], id[50];
}Gene_s;

typedef struct{
    double *exp;
    char id[50];
}Exp_s;

void openFiles(int argc, char *argv[]);
Gene_s *readGenes(FILE *gene_file, int *n);
SV_s *readVCF(FILE *vcf_file, Gene_s *genes, int *n, int *m, int gene_n, int type, double maf, double len[]);
Exp_s *readExp(FILE *exp_file, int sv_m, int *n);
void defStats(SV_s *svl, Exp_s *rna, Exp_s *dna, int sv_n, int sv_m, int rna_n, int dna_n);
int *splitInfo(char str[]);
int isNumeric(const char *s);
void lineTerminator(char *line);

int main(int argc, char *argv[]){
    
    int second=0, minute=0, hour=0;
    time_t timer=0;
    
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
    
    int i, j, type=-1, sv_n=0, sv_m=0, rna_n=0, dna_n=0, gene_n=0;
    double ol=0, maf=0, len[2]={0};
    char *temp1=NULL, *temp2=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    SV_s *svl;
    Gene_s *genes;
    Exp_s *rna, *dna;
    FILE *vcf_file, *rna_file, *dna_file, *gene_file;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((vcf_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sv %s\n", argv[i]);
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
        
        else if(strcmp(argv[i], "-genes") == 0){
            if((gene_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-genes %s\n", argv[i]);
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
        
        else if(strcmp(argv[i], "-maf") == 0){
            if(isNumeric(argv[++i]))
                maf = atof(argv[i]);
            fprintf(stderr, "\t-maf %.2f\n", maf);
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
    
    if(vcf_file == NULL | rna_file == NULL || dna_file == NULL || gene_file == NULL){
        fprintf(stderr,"\nERROR: The following are required: -sv [file] -rna [file] -dna [file] -genes [file]\n");
        exit(EXIT_FAILURE);
    }
    
    genes = readGenes(gene_file, &gene_n);
    svl = readVCF(vcf_file, genes, &sv_n, &sv_m, gene_n, type, maf, len);
    rna = readExp(rna_file, sv_n, &rna_n);
    dna = readExp(dna_file, sv_n, &dna_n);
    defStats(svl, rna, dna, sv_n, sv_m, rna_n, dna_n);
}

Gene_s *readGenes(FILE *gene_file, int *n){
    
    int char_i=0, line_i=0, maxchar=0;
    char c, *line=NULL, *temp=NULL;
    Gene_s *list;
    
    while((c=fgetc(gene_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(gene_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((list=malloc(line_i*sizeof(Gene_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while(fgets(line, maxchar+1, gene_file) != NULL){
        lineTerminator(line);
        strncpy(list[*n].chr,strtok(line,"\t"),19);
        list[*n].start = atoi(strtok(NULL,"\t"));
        list[*n].stop = atoi(strtok(NULL,"\t"));
        list[*n].len = list[*n].stop - list[*n].start + 1;
        strncpy(list[*n].id,strtok(NULL,"\t"),49);
        temp = strtok(NULL,"\t");
        if(temp[0] == '+')
            list[*n].ph = 0;
        else if(temp[0] == '-')
            list[*n].ph = 1;
        else if(isNumeric(temp))
            list[*n].ph = atoi(temp);
        *n = *n + 1;
    }
    
    free(line);
    fclose(gene_file);
    
    return list;
}

SV_s *readVCF(FILE *vcf_file, Gene_s *genes, int *n, int *m, int gene_n, int type, double maf, double len[]){
    
    int i, j, k, char_i=0, maxchar=0, tab_i=0, maxtab=0, line_i=0, id_s=0, ok=0, dist=0, *info;
    double p=0;
    char c, *temp=NULL, *line=NULL, *end=NULL, *gene=NULL, **id=NULL;
    SV_s *svl;
    
    while((c=fgetc(vcf_file)) != EOF){
        char_i++;
        if(c == '\t')
            tab_i++;
        if(c == '\n'){
            line_i++;
            if(tab_i > maxtab)
                maxtab = tab_i;
            tab_i = 0;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(vcf_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((svl=malloc(line_i*sizeof(SV_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<line_i;i++){
        if((svl[i].geno=malloc(maxtab*sizeof(int))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    if((svl[0].ind=malloc(maxtab*sizeof(char*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<maxtab;i++){
        if((svl[0].ind[i]=malloc(101*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while(fgets(line, maxchar+1, vcf_file) != NULL){
        lineTerminator(line);
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        p = 0;
        if(strcmp(temp,"#CHROM") == 0){
            while(temp != NULL){
                if(i > 9){
                    strncpy(svl[0].ind[j],temp,100);
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
                strncpy(svl[*m].chr, temp, 100);
            else if(i == 2)
                svl[*m].start = atoi(temp);
            else if(i == 8){
                info = splitInfo(temp);
                if(type != -1 && info[0] != type)
                    break;
                if(len[1] > 0 && (info[2] < len[0] || info[2] > len[1]))
                    break;
                svl[*m].type = info[0];
                svl[*m].stop = info[1];
                if(svl[*m].start > svl[*m].stop){
                    svl[*m].stop = svl[*m].start;
                    svl[*m].start = info[1];
                }
                svl[*m].len = info[2];
                svl[*m].gene_i = 0;
                for(k=0;k<gene_n;k++){
                    if(strcmp(svl[*m].chr,genes[k].chr) == 0){
                        strncpy(svl[*m].id[svl[*m].gene_i], genes[k].id, 49);
                        if(genes[k].start > svl[*m].stop){
                            if(genes[k].ph == 0)
                                svl[*m].dist[svl[*m].gene_i] = genes[k].start - svl[*m].stop;
                            else
                                svl[*m].dist[svl[*m].gene_i] = svl[*m].stop - genes[k].start;
                            if(abs(svl[*m].dist[svl[*m].gene_i]) <= 5000){
                                svl[*m].gene_i++;
                                if(svl[*m].gene_i == 1000){
                                    fprintf(stderr,"Warning: max SV amount of 1000 reached\n");
                                    break;
                                }
                            }
                        }
                        else if(genes[k].stop < svl[*m].start){
                            if(genes[k].ph == 1)
                                svl[*m].dist[svl[*m].gene_i] = svl[*m].start - genes[k].stop;
                            else
                                svl[*m].dist[svl[*m].gene_i] = genes[k].stop - svl[*m].start;
                            if(abs(svl[*m].dist[svl[*m].gene_i]) <= 5000){
                                svl[*m].gene_i++;
                                if(svl[*m].gene_i == 1000){
                                    fprintf(stderr,"Warning: max SV amount of 1000 reached\n");
                                    break;
                                }
                            }
                        }
                        else if(svl[*m].start <= genes[k].stop && svl[*m].stop >= genes[k].start){
                            svl[*m].dist[svl[*m].gene_i] = 0;
                            svl[*m].gene_i++;
                            if(svl[*m].gene_i == 1000){
                                fprintf(stderr,"Warning: max SV amount of 1000 reached\n");
                                break;
                            }
                        }
                    }
                    else if(strcmp(svl[*m].chr, genes[k].chr) < 0)
                        break;
                }
            }
            else if(i > 9){
                if(temp[0] == '0' && temp[2] == '0')
                    svl[*m].geno[j] = 0;
                else if(temp[0] == '0' && temp[2] == '1')
                    svl[*m].geno[j] = 1;
                else if(temp[0] == '1' && temp[2] == '0')
                    svl[*m].geno[j] = 1;
                else if(temp[0] == '1' && temp[2] == '1')
                    svl[*m].geno[j] = 2;
                else
                    svl[*m].geno[j] = 9;
                if(svl[*m].geno[j] != 9)
                    p += svl[*m].geno[j];
                j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
            if(temp == NULL){
                p /= j * 2;
                if(p < maf || 1 - p < maf)
                    break;
                if(p > 0.5)
                    svl[*m].ma = 1;
                else
                    svl[*m].ma = 0;
                if(*n == 0)
                    *n = j;
                *m = *m + 1;
            }
        }
    }
    
    free(line);
    
    return svl;
}

Exp_s *readExp(FILE *exp_file, int sv_n, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0;
    char c, *line=NULL, *temp=NULL;
    Exp_s *exp;
    
    while((c=fgetc(exp_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(exp_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((exp=malloc(line_i*sizeof(Exp_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<line_i;i++){
        if((exp[i].exp=malloc(sv_n*sizeof(double))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while(fgets(line, maxchar+1, exp_file) != NULL){
        if(line[0] == '#')
            continue;
        lineTerminator(line);
        strncpy(exp[*n].id, strtok(line,"\t"), 49);
        temp = strtok(NULL,"\t");
        i = 0;
        while(temp != NULL){
            exp[*n].exp[i] = atof(temp);
            temp = strtok(NULL,"\t");
            i++;
        }
        *n = *n + 1;
    }
    
    free(line);
    fclose(exp_file);
    
    return exp;
}

void defStats(SV_s *svl, Exp_s *rna, Exp_s *dna, int sv_n, int sv_m, int rna_n, int dna_n){
    
    int i, j, k, l, m, n=0, ok=0;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    
    printf("n\ttype\tsize\tdist\tid\tgene\trna\tdna\tgeno\n");
    for(i=0;i<sv_m;i++){
        for(j=0;j<svl[i].gene_i;j++){
            ok = 0;
            for(k=0;k<rna_n;k++){
                if(strcmp(svl[i].id[j], rna[k].id) == 0){
                    for(l=0;l<dna_n;l++){
                        if(strcmp(svl[i].id[j], dna[l].id) == 0){
                            ok = 1;
                            break;
                        }
                    }
                    break;
                }
            }
            if(ok == 0)
                continue;
            for(m=0;m<sv_n;m++){
                printf("%i\t%s\t%i\t%i\t%s\t%s\t%f\t%f\t",n,vars[svl[i].type],svl[i].len,svl[i].dist[j],svl[0].ind[m],rna[k].id,rna[k].exp[m],dna[l].exp[m]);
               if(svl[i].ma == 0){
                   if(svl[i].geno[m] == 1)
                       printf("1\n");
                   else if(svl[i].geno[m] == 2)
                       printf("2\n");
                   else
                       printf("0\n");
               }
               else{
                   if(svl[i].geno[m] == 1)
                       printf("1\n");
                   else if(svl[i].geno[m] == 0)
                       printf("2\n");
                   else
                       printf("0\n");
                }
            }
            n++;
        }
    }
    
    free(svl);
    free(rna);
    free(dna);
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
