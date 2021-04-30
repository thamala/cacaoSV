// ./sift2SV -g SV.geno -s Amaz-15_15_S19_SIFTannotations.xls -i 0 -t INV -l 500-inf > test.txt

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

 Program for compiling SIFT scores for SNPs in the major and minor SV arrangements

 compiling: gcc SV2sift.c -o SV2sift -lm

 usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -sift [file] A file produced by SIFT4G_Annotator.jar
 -ind [ind] Number defining which individual from the VCF files is analyzed
 -t [str] Either INS, DEL, DUP, TRA, or INV
 -l [int]-[int] Minimum and maximum length of the SVs.

 example:
 ./SV2sift -sv SV.vcf -sift SIFTannotations.xls -ind 0 -t INV -l 1e4-1e6 > out.txt
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
    int type, start, stop, len, ind;
    char chr[101];
}SV_s;

typedef struct{
    int pos, geno;
    double sift, id;
    char chr[101];
}Sift_s;

void openFiles(int argc, char *argv[]);
SV_s *readSV(FILE *sv_file, int *n, int ind, int type, double len[]);
Sift_s *readSift(FILE *sift_file, int *n);
void defStat(SV_s *svl, Sift_s *sift, int ind, int sv_n, int sift_n);
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
    
    int i, j, type=-1, ind=-1, sv_n=0, sift_n=0;
    double af=1, len[2]={0};
    char *temp1=NULL, *temp2=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    SV_s *svl;
    Sift_s *sift;
    FILE *sv_file=NULL, *sift_file=NULL;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((sv_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sv %s\n", argv[i]);
        }

        else if(strcmp(argv[i], "-sift") == 0){
            if((sift_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-sift %s\n", argv[i]);
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
    
    if(sv_file==NULL || sift_file==NULL || ind == -1){
        fprintf(stderr,"\nERROR: -sv [file] -sift [file] -ind [int] are needed\n");
        exit(EXIT_FAILURE);
    }
    
    svl = readSV(sv_file, &sv_n, ind, type, len);
    sift = readSift(sift_file, &sift_n);
    defStat(svl, sift, ind, sv_n, sift_n);
}

SV_s *readSV(FILE *sv_file, int *n, int ind, int type, double len[]){
    
    int i, j, k, char_i=0, maxchar=0, line_i=0, geno=0, *info;
    double geno_i=0;
    char c, *line, *temp, *end;
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
    
    while(fgets(line, maxchar+1, sv_file) != NULL){
        lineTerminator(line);
        if(line[0] == '#')
            continue;
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        k = 0;
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
                    geno = 0;
                else if(temp[0] == '0' && temp[2] == '1')
                    geno = 1;
                else if(temp[0] == '1' && temp[2] == '0')
                    geno = 1;
                else if(temp[0] == '1' && temp[2] == '1')
                    geno = 2;
                else
                    geno = 9;
                if(geno != 9)
                    k += geno;
                if(j == ind)
                    svl[*n].ind = geno;
                j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
            if(temp == NULL){
                j *= 2;
                if((double)k/(double)j == 1)
                    continue;
                else if((double)k/(double)j > 0.5){
                    if(svl[*n].ind == 0)
                        svl[*n].ind = 2;
                    else if(svl[*n].ind == 2)
                        svl[*n].ind = 0;
                }
                *n = *n + 1;
            }
        }        
    }
    
    free(line);
    fclose(sv_file);
    
    return svl;
}

Sift_s *readSift(FILE *sift_file, int *n){
    
    int i=0, char_i=0, maxchar=0, line_i=0;
    char c, *line, *temp;
    Sift_s *list;
    
    while((c=fgetc(sift_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(sift_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((list=malloc(line_i*sizeof(Sift_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    fgets(line, maxchar+1, sift_file);
    
    while(fgets(line, maxchar+1, sift_file) != NULL){
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        strncpy(list[*n].chr, strtok(line,"\t"), 100);
        list[*n].pos = atoi(strtok(NULL,"\t"));
        temp = strtok(NULL, "\t");
        i = 1;
        while(temp != NULL){
            if(i == 11 && strcmp(temp,"NA")!=0)
                list[*n].sift = atof(temp);
            else if(i == 12 && strcmp(temp,"NA")!=0){
                if(atof(temp) >= 2.75 && atof(temp) <= 3.5)
                    *n = *n + 1;
            }
            temp = strtok(NULL,"\t");
            i++;
        }
    }
    
    
    free(line);
    fclose(sift_file);
    
    return list;
}

void defStat(SV_s *svl, Sift_s *sift, int ind, int sv_n, int sift_n){
    
    int i, j, k=0, het_i=0, hom_i=0, tot_i=0, ok=0;
    double tot=0, het=0, hom=0;
    Sift_s *outs;
    
    if((outs=malloc(sift_n*sizeof(Sift_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<sv_n;i++){
        ok = 0;
        for(j=0;j<sift_n;j++){
            if(strcmp(svl[i].chr, sift[j].chr) == 0){
                if(svl[i].start <= sift[j].pos && svl[i].stop >= sift[j].pos){
                    if(svl[i].ind == 0){
                        strcpy(outs[k].chr, sift[j].chr);
                        outs[k].pos = sift[j].pos;
                        outs[k].sift = sift[j].sift;
                        outs[k].geno = 0;
                        k++;
                        if(sift[j].sift < 0.05)
                            tot++;
                        tot_i++;
                    }
                    else if(svl[i].ind == 1){
                        strcpy(outs[k].chr, sift[j].chr);
                        outs[k].pos = sift[j].pos;
                        outs[k].sift = sift[j].sift;
                        outs[k].geno = 1;
                        k++;
                        if(sift[j].sift < 0.05)
                            het++;
                        het_i++;
                    }
                    else if(svl[i].ind == 2){
                        strcpy(outs[k].chr, sift[j].chr);
                        outs[k].pos = sift[j].pos;
                        outs[k].sift = sift[j].sift;
                        outs[k].geno = 2;
                        k++;
                        if(sift[j].sift < 0.05)
                            hom++;
                        hom_i++;
                    }
                }
                else if(svl[i].start < sift[j].pos)
                    break;
            }
            else if(strcmp(svl[i].chr, sift[j].chr) < 0)
                break;
        }
    }
    
    if(het_i > 0 || hom_i > 0){
        printf("chr pos sift sv id\n");
        for(i=0;i<k;i++)
            printf("%s %i %f %i %i\n", outs[i].chr, outs[i].pos, outs[i].sift, outs[i].geno, ind);
    }
    
    tot /= (double)tot_i;
    het /= (double)het_i;
    hom /= (double)hom_i;
    
    fprintf(stderr,"\nmajor-major = %f n = %i\nmajor-minor = %f n = %i\nminor-minor = %f n = %i\n", tot, tot_i, het, het_i, hom, hom_i);
    
    free(svl);
    free(sift);
    free(outs);
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
