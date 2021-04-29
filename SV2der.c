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

 Program for counting the numer of derived alleles in the major and minor SV arrangements

 compiling: gcc SV2der.c -o SV2der -lm

 usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -snp [file] SNP-vcf file for the ingroup, including variant and invariant sites.
 -anc [file] SNP-vcf file for the outgroup, including variant and invariant sites.
 -sites [file] A tab-delimited file with chr and pos for sites to use.
 -excl [file] A tab-delimited file with chr, start, and end of regions to exclude.
 -maf [double] Minor allele frequency cutoff for the SV.
 -t [str] Either INS, DEL, DUP, TRA, or INV
 -l [int]-[int] Minimum and maximum length of the SVs.

 example:
 ./SV2der -sv SV.vcf -snp SNP.vcf -anc ANC.vcf -sites 0fold.sites -t INV -l 1e4-1e6 > out.txt
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
    int type, start, stop, len, *geno;
    char chr[101], **id;
}SV_s;

typedef struct{
    int pos;
    char chr[101], nuc;
}Site_s;

typedef struct{
    int start, stop;
    char chr[101];
}Gene_s;

typedef struct{
    int pos, *geno;
    char chr[101];
}Snp_s;

void openFiles(int argc, char *argv[]);
Site_s *readSites(FILE *site_file, int *n);
Gene_s *readGenes(FILE *gene_file, int *n);
SV_s *readSV(FILE *sv_file, double maf, int *l, int *n, int *m, int type, double len[]);
Site_s *readAnc(FILE *anc_file, SV_s *svl, Site_s *sites, Gene_s *genes, int sv_l, int sv_n, int site_n, int gene_n, int *n);
Snp_s *readSNP(FILE *snp_file, Site_s *anc, int anc_n, int sv_m, int *n);
void defStat(SV_s *svl, Snp_s *snps, int ind, int sv_n, int sv_m, int snp_n);
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
    
    int i, j, ind=-1, type=-1, sv_l=0, sv_n=0, sv_m=0, site_n=0, gene_n=0, anc_n=0, snp_n=0;
    double maf=0, len[2]={0};
    char *temp1=NULL, *temp2=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    SV_s *svl=NULL;
    Site_s *sites=NULL, *anc=NULL;
    Gene_s *genes=NULL;
    Snp_s *snps=NULL;
    FILE *sv_file=NULL, *anc_file=NULL, *snp_file=NULL, *site_file=NULL, *gene_file=NULL;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((sv_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sv %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-anc") == 0){
            if((anc_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-anc %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-snp") == 0){
            if((snp_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-snp %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-sites") == 0){
            if((site_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sites %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-excl") == 0){
            if((gene_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-excl %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-i") == 0){
            ind = atoi(argv[++i]);
            fprintf(stderr,"\t-i %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-maf") == 0){
            maf = atof(argv[++i]);
            fprintf(stderr,"\t-maf %s\n", argv[i]);
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
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    
    fprintf(stderr,"\n");
    
    if(sv_file==NULL || anc_file==NULL || snp_file==NULL){
        fprintf(stderr,"\nERROR: -sv [file] -anc [file] -snp [file] are required\n");
        exit(EXIT_FAILURE);
    }
    
    if(site_file != NULL)
        sites = readSites(site_file, &site_n);
    if(gene_file != NULL)
        genes = readGenes(gene_file, &gene_n);
    svl = readSV(sv_file, maf, &sv_l, &sv_n, &sv_m, type, len);
    anc = readAnc(anc_file, svl, sites, genes, sv_l, sv_n, site_n, gene_n, &anc_n);
    snps = readSNP(snp_file, anc, anc_n, sv_m, &snp_n);
    defStat(svl, snps, ind, sv_n, sv_m, snp_n);
}

Site_s *readSites(FILE *site_file, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0;
    char c, *line=NULL;
    Site_s *list=NULL;
    
    while((c=fgetc(site_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(site_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((list=malloc(line_i*sizeof(Site_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while(fgets(line, maxchar+1, site_file) != NULL){
        lineTerminator(line);
        strncpy(list[*n].chr, strtok(line,"\t"), 100);
        list[*n].pos = atoi(strtok(NULL,"\t"));
        *n = *n + 1;
    }
    
    free(line);
    fclose(site_file);
    
    return list;
}

Gene_s *readGenes(FILE *gene_file, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0;
    char c, *line=NULL;
    Gene_s *list=NULL;
    
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
        strncpy(list[*n].chr, strtok(line,"\t"), 100);
        list[*n].start = atoi(strtok(NULL,"\t"));
        list[*n].stop = atoi(strtok(NULL,"\t"));
        *n = *n + 1;
    }
    
    free(line);
    fclose(gene_file);
    
    return list;
}

SV_s *readSV(FILE *sv_file, double maf, int *l, int *n, int *m, int type, double len[]){
    
    int i, j, k, char_i=0, maxchar=0, tab_i=0, maxtab=0, line_i=0, id_s=0, ok=0, *info;
    double p=0;
    char c, *temp=NULL, *line=NULL, *end=NULL;
    SV_s *svl=NULL;
    
    while((c=fgetc(sv_file)) != EOF){
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
    
    rewind(sv_file);
    
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
    if((svl[0].id=malloc(maxtab*sizeof(char*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<maxtab;i++){
        if((svl[0].id[i]=malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while(fgets(line, maxchar+1, sv_file) != NULL){
        lineTerminator(line);
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        p = 0;
        if(strcmp(temp,"#CHROM") == 0){
            while(temp != NULL){
                if(i > 9){
                    strcpy(svl[0].id[j], temp);
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
                svl[*n].type = info[0];
                svl[*n].stop = info[1];
                svl[*n].len = info[2];
            }
            else if(i > 9){
                if(temp[0] == '0' && temp[2] == '0')
                    svl[*n].geno[j] = 0;
                else if(temp[0] == '0' && temp[2] == '1')
                    svl[*n].geno[j] = 1;
                else if(temp[0] == '1' && temp[2] == '0')
                    svl[*n].geno[j] = 1;
                else if(temp[0] == '1' && temp[2] == '1')
                    svl[*n].geno[j] = 2;
                else
                    svl[*n].geno[j] = 9;
                if(svl[*n].geno[j] != 9)
                    p += svl[*n].geno[j];
                j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
            if(temp == NULL && (p/((double)j*2) >= maf && p/((double)j*2) <= (1-maf))){
                if(*m == 0)
                    *m = j;
                if(p/((double)j*2) > 0.5){
                    for(i=0;i<*m;i++){
                        if(svl[*n].geno[i] == 0)
                            svl[*n].geno[i] = 2;
                        else if(svl[*n].geno[i] == 2)
                            svl[*n].geno[i] = 0;
                    }
                }
                *l = *l + svl[*n].len;
                *n = *n + 1;
            }
        }
    }
    
    free(line);
    fclose(sv_file);
    
    return svl;
}

Site_s *readAnc(FILE *anc_file, SV_s *svl, Site_s *sites, Gene_s *genes, int sv_l, int sv_n, int site_n, int gene_n, int *n){
    
    int i, j, k=0, l=0, char_i=0, line_i=0, maxchar=0, pos=0, ok=0, ref_i=0, alt_i=0;
    char c, chr[20], ref, alt, *line=NULL, *temp=NULL;
    Site_s *anc=NULL;
    size_t len=0;
    ssize_t read;
    
    if(site_n > 0 && site_n < sv_l)
        sv_l = site_n;
    if((anc=malloc(sv_l*sizeof(Site_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, anc_file)) != -1){
        lineTerminator(line);
        temp = strtok(line,"\t");
        ok = 0;
        i = 0;
        j = 1;
        ref_i = 0;
        alt_i = 0;
        if(temp[0] != '#'){
            strncpy(chr, temp, 19);
            while(temp != NULL){
                if(j == 2){
                    pos =  atoi(temp);
                    while(k < sv_n){
                        if(strcmp(chr,svl[k].chr) == 0){
                            if(pos >= svl[k].start && pos <= svl[k].stop){
                                ok = 1;
                                break;
                            }
                            else if(pos < svl[k].start)
                                break;
                        }
                        else if(strcmp(chr,svl[k].chr) < 0)
                            break;
                        k++;
                    }
                    if(ok == 0)
                        break;
                    if(site_n > 0){
                        ok = 0;
                        while(l < site_n){
                            if(strcmp(chr,sites[l].chr) == 0){
                                if(pos == sites[l].pos){
                                    ok = 1;
                                    break;
                                }
                                else if(pos < sites[l].pos)
                                    break;
                            }
                            else if(strcmp(chr,sites[l].chr) < 0)
                                break;
                            l++;
                        }
                        if(ok == 0)
                            break;
                    }
                    else if(gene_n > 0){
                        ok = 0;
                        while(l < gene_n){
                            if(strcmp(chr,genes[l].chr) == 0){
                                if(pos >= genes[l].start && pos <= genes[l].stop){
                                    ok = 1;
                                    break;
                                }
                                else if(pos < genes[l].start)
                                    break;
                            }
                            else if(strcmp(chr,genes[l].chr) < 0)
                                break;
                            l++;
                        }
                        if(ok == 1)
                            break;
                    }
                }
                else if(j == 4)
                    ref = temp[0];
                else if(j == 5)
                    alt = temp[0];
                else if(j > 9){
                    if(temp[0] == '0' && temp[2] == '0')
                        ref_i += 2;
                    else if(temp[0] == '1' && temp[2] == '0'){
                        ref_i++;
                        alt_i++;
                    }
                    else if(temp[0] == '0' && temp[2] == '1'){
                        ref_i++;
                        alt_i++;
                    }
                    else if(temp[0] == '1' && temp[2] == '1')
                        alt_i += 2;
                    i++;
                }
                temp = strtok(NULL,"\t");
                j++;
                if(temp == NULL){
                    if((alt_i == 0 && ref_i > 0) || (ref_i == 0 && alt_i > 0)){
                        strcpy(anc[*n].chr, chr);
                        anc[*n].pos = pos;
                        if(ref_i > 0)
                            anc[*n].nuc = ref;
                        else
                            anc[*n].nuc = alt;
                        *n = *n + 1;
                    }
                }
            }
        }
    }
    
    if(site_n > 0)
        free(sites);
    if(gene_n > 0)
        free(genes);
    free(line);
    fclose(anc_file);
    
    return anc;
}

Snp_s *readSNP(FILE *snp_file, Site_s *anc, int anc_n, int sv_m, int *n){
    
    int i, j, k=0, ok=0;
    char ref, alt, *line=NULL, *temp=NULL;
    Snp_s *snps;
    size_t len=0;
    ssize_t read;
    
    if((snps=malloc(anc_n*sizeof(Snp_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<anc_n;i++){
        if((snps[i].geno=malloc(sv_m*sizeof(int))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while((read = getline(&line, &len, snp_file)) != -1){
        lineTerminator(line);
        temp = strtok(line,"\t");
        ok = 0;
        i = 0;
        j = 1;
        if(temp[0] != '#'){
            strncpy(snps[*n].chr, temp, 19);
            while(temp != NULL){
                if(j == 2){
                    snps[*n].pos =  atoi(temp);
                    while(k < anc_n){
                        if(strcmp(snps[*n].chr,anc[k].chr) == 0){
                            if(snps[*n].pos == anc[k].pos){
                                ok = 1;
                                break;
                            }
                            else if(snps[*n].pos < anc[k].pos)
                                break;
                        }
                        else if(strcmp(snps[*n].chr,anc[k].chr) < 0)
                            break;
                        k++;
                    }
                    if(ok == 0)
                        break;
                }
                else if(j == 4)
                    ref = temp[0];
                else if(j == 5){
                    alt = temp[0];
                    if(anc[k].nuc != ref && anc[k].nuc != alt)
                        break;
                }
                else if(j > 9){
                    if(temp[0] == '0' && temp[2] == '0'){
                        if(anc[k].nuc == ref)
                            snps[*n].geno[i] = 0;
                        else
                            snps[*n].geno[i] = 2;
                    }
                    else if(temp[0] == '1' && temp[2] == '0')
                        snps[*n].geno[i] = 1;
                    else if(temp[0] == '0' && temp[2] == '1')
                        snps[*n].geno[i] = 1;
                    else if(temp[0] == '1' && temp[2] == '1'){
                        if(anc[k].nuc == alt)
                            snps[*n].geno[i] = 0;
                        else
                            snps[*n].geno[i] = 2;
                    }
                    else
                        snps[*n].geno[i] = 9;
                    i++;
                }
                temp = strtok(NULL,"\t");
                j++;
                if(temp == NULL)
                    *n = *n + 1;
            }
        }
    }
    
    free(anc);
    free(line);
    fclose(snp_file);
    
    return snps;
}

void defStat(SV_s *svl, Snp_s *snps, int ind, int sv_n, int sv_m, int snp_n){
    
    int i, j, k;
    double geno[3][3], prob[3];
    
    memset(geno,0,sizeof(geno[0][0])*3*3);
    
    printf("sv\tid\tder\n");
    
    for(i=0;i<sv_n;i++){
        for(j=0;j<snp_n;j++){
            if(strcmp(svl[i].chr, snps[j].chr) == 0){
                if(svl[i].start <= snps[j].pos && svl[i].stop >= snps[j].pos){
                    for(k=0;k<sv_m;k++){
                        if((svl[i].geno[k] == 9 || snps[j].geno[k] == 9) || (ind != -1 && k != ind))
                            continue;
                        geno[svl[i].geno[k]][snps[j].geno[k]]++;
                        if(snps[j].geno[k] == 0){
                            printf("%i\t%s\t0\n", svl[i].geno[k], svl[0].id[k]);
                            printf("%i\t%s\t0\n", svl[i].geno[k], svl[0].id[k]);
                        }
                        else if(snps[j].geno[k] == 1){
                            printf("%i\t%s\t0\n", svl[i].geno[k], svl[0].id[k]);
                            printf("%i\t%s\t1\n", svl[i].geno[k], svl[0].id[k]);
                        }
                        else if(snps[j].geno[k] == 2){
                            printf("%i\t%s\t1\n", svl[i].geno[k], svl[0].id[k]);
                            printf("%i\t%s\t1\n", svl[i].geno[k], svl[0].id[k]);
                        }
                        
                    }
                }
                else if(svl[i].start < snps[j].pos)
                    break;
            }
            else if(strcmp(svl[i].chr, snps[j].chr) < 0)
                break;
        }
    }
    
    prob[0] = (2*geno[0][2]+geno[0][1])/((geno[0][0]+geno[0][1]+geno[0][2])*2);
    prob[1] = (2*geno[1][2]+geno[1][1])/((geno[1][0]+geno[1][1]+geno[1][2])*2);
    prob[2] = (2*geno[2][2]+geno[2][1])/((geno[2][0]+geno[2][1]+geno[2][2])*2);
    fprintf(stderr,"\tanc-anc\tanc-der\tder-der\tder%\n");
    fprintf(stderr,"major-major\t%.0f\t%.0f\t%.0f\t%f\n", geno[0][0], geno[0][1], geno[0][2], prob[0]);
    fprintf(stderr,"major-minor\t%.0f\t%.0f\t%.0f\t%f\n", geno[1][0], geno[1][1], geno[1][2], prob[1]);
    fprintf(stderr,"minor-minor\t%.0f\t%.0f\t%.0f\t%f\n", geno[2][0], geno[2][1], geno[2][2], prob[2]);
    
    free(svl);
    free(snps);
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
