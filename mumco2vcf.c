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

 Program for combining SV calls from MUM&Co into a VCF file

 Compiling: gcc mumco2vcf.c -o mumco2vcf -lm

 Usage:
 -i [file] File listing samples to combine. Haplotypes are in one row, separated by a tab.
 -t [str] either INS, DEL, DUP, TRA, or INV. (optional)
 -l [int]-[int] minimum and maximum length of the SVs. (optional)
 -o [int] overlap in base pairs to combine variants. (optional)

 Examples:
 ./mumco2sv -i list.txt -t DEL -l 500-1e5 -o 250 > out.vcf
 ./mumco2sv -i list.txt -t DEL -l 500-1e5 -o 250 | sort -k1,1 -k2,2n > sorted.vcf
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#define merror "\nERROR: System out of memory\n\n"
#define version "2019.05.19"

typedef struct{
    char *haplo1, *haplo2, *names;
}List_s;

typedef struct{
    int type, start, stop, len, stat, het, ok;
    char chr[101];
}SV_s;

void openFiles(int argc, char *argv[]);
List_s *readList(FILE *list_file, int *n, int *dip);
SV_s *readSV(FILE *sv_file, int *n, int *l_ok, int type, double len[], double gb);
SV_s **defHet(SV_s **svp1, SV_s **svp2, int n, int m, double gb);
void defGeno(SV_s **svp, int n, double gb);
int isNumeric (const char *s);
void lineTerminator(char *line);
void printInfo(void);

int main(int argc, char *argv[]){
    
    int second=0, minute=0, hour=0;
    time_t timer=0;
    
    if(argc == 1){
        printInfo();
        exit(EXIT_FAILURE);
    }
    
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
    
    int i, j, type=-1, l_ok=0, list_n=0, hp=0, tmp_n=0, sv_n=0;
    double gb=50, len[2]={0};
    char *temp1=NULL, *temp2=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    List_s *list=NULL;
    SV_s **svp1=NULL, **svp2=NULL, **svp=NULL;
    FILE *list_file=NULL, *sv_file=NULL;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0){
            if((list_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-i %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0){
            argv[++i];
            for(j=0;j<5;j++){
                if(strcmp(argv[i], vars[j]) == 0)
                    type = j;
            }
            if(type == -1){
                fprintf(stderr,"\nERROR: Unknown parameter for -t\n");
                printInfo();
                exit(EXIT_FAILURE);
            }
            fprintf(stderr, "\t-t %s\n", vars[type]);
        }
        
        else if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "--length") == 0){
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
                printInfo();
                exit(EXIT_FAILURE);
            }
        }
        
        else if(strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--overlap") == 0){
            if(isNumeric(argv[++i]))
                gb = atof(argv[i]);
            fprintf(stderr, "\t-o %s\n", argv[i]);
            
        }
        
        else{
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n", argv[i]);
            printInfo();
            exit(EXIT_FAILURE);
        }
    }
    
    list = readList(list_file, &list_n, &hp);
    
    if((svp1=malloc(list_n*sizeof(SV_s*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if(hp == 1){
        if((svp2=malloc(list_n*sizeof(SV_s*))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    for(i=0;i<list_n;i++){
        if((sv_file=fopen(list[i].haplo1,"r"))==NULL){
            fprintf(stderr,"\nERROR:Cannot open file %s\n", list[i].haplo1);
            exit(EXIT_FAILURE);
        }
        svp1[i] = readSV(sv_file, &tmp_n, &l_ok, type, len, gb);
        if(tmp_n > sv_n)
            sv_n = tmp_n;
        tmp_n = 0;
        if(hp == 1){
            if((sv_file=fopen(list[i].haplo2,"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", list[i].haplo2);
                exit(EXIT_FAILURE);
            }
            svp2[i] = readSV(sv_file, &tmp_n, &l_ok, type, len, gb);
            if(tmp_n > sv_n)
                sv_n = tmp_n;
            tmp_n = 0;
        }
    }
    
    if(l_ok == 1)
        fprintf(stderr,"\nWarning: Chromosomes have long names. Make sure they are distinguishable from the first 100 characters.\n");
    
    printf("##fileformat=VCFv4.0\n");
    printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for(i=0;i<list_n;i++)
        printf("\t%s",list[i].haplo1);
    printf("\n");
    
    if(hp == 1)
        svp = defHet(svp1, svp2, list_n, sv_n, gb);
    else
        svp = svp1;
    
    defGeno(svp, list_n, gb);
}

List_s *readList(FILE *list_file, int *n, int *dip){
    
    int i, char_i=0, maxchar=0, line_i=0;
    char c, *line, *temp=NULL;
    List_s *list;
    
    while((c=fgetc(list_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(list_file);
    
    if((line = malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((list = malloc(line_i*sizeof(List_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<line_i;i++){
        if((list[i].haplo1 = malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        if((list[i].haplo2 = malloc(maxchar*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while(fgets(line, maxchar+1, list_file) != NULL){
        if(line[0] != '#' && line[0] != '\n' && line[0] != '\r'){
            lineTerminator(line);
            temp = strtok(line, "\t");
            strcpy(list[*n].haplo1, temp);
            temp = strtok(NULL, "\t");
            if(temp != NULL){
                strcpy(list[*n].haplo2, temp);
                *dip = 1;
                temp = strtok(NULL,"\t");
                if(temp != NULL){
                    fprintf(stderr,"\nWarning: List-file contains > 2 columns. Only the first 2 are used to define haplotypes.\n\n");
                }
            }
            *n = *n + 1;
        }
    }
    
    free(line);
    fclose(list_file);
    
    return list;
}

SV_s *readSV(FILE *sv_file, int *n, int *l_ok, int type, double len[], double gb){
    
    int i=0, j=0, k=0, char_i=0, maxchar=0, line_i=0, ok=0, tmp_n=0;
    char c, qchr[101], *line, *temp;
    const char *vars[5]={"insertion","deletion","duplication","transloc","inversion"};
    SV_s *svl, *svl_cp;
    
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
    if((svl=malloc((line_i+1)*sizeof(SV_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if(gb > 0){
        if((svl_cp=malloc((line_i+1)*sizeof(SV_s))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while(fgets(line, maxchar+1, sv_file) != NULL){
        if(line[0] == '\n' || line[0] == '\r')
            continue;
        lineTerminator(line);
        temp = strtok(line,"\t");
        if(strcmp(temp, "ref_chr") == 0)
            continue;
        if(strlen(temp) >= 100)
            *l_ok = 1;
        strncpy(svl[i].chr, temp, 100);
        temp = strtok(NULL, "\t");
        strncpy(qchr, temp, 100);
        svl[i].start = atoi(strtok(NULL, "\t"));
        svl[i].stop = atoi(strtok(NULL, "\t"));
        if(svl[i].start > svl[i].stop){
            svl[i].start = svl[i].start + svl[i].stop;
            svl[i].stop = svl[i].start - svl[i].stop;
            svl[i].start = svl[i].start - svl[i].stop;
        }
        if(i > 0 && (abs(svl[i].start-svl[i-1].start) <= 10 && abs(svl[i].stop-svl[i-1].stop) <= 10))
            continue;
        svl[i].len = atoi(strtok(NULL,"\t"));
        temp = strtok(NULL, "\t");
        for(j=0;j<5;j++){
            if(strcmp(temp, vars[j]) == 0)
                svl[i].type = j;
        }
        if(svl[i].type == 0 && svl[i].stop-svl[i].start > 100000)
            continue;
        for(j=0;j<3;j++)
            temp = strtok(NULL,"\t");
        if(temp != NULL){
            if(strcmp(temp,"complicated") == 0 || strcmp(temp,"double") == 0)
                continue;
        }
        if(svl[i].type == -1 || svl[i].len <= 0)
            continue;
        if(type != -1 && svl[i].type != type)
            continue;
        if(gb == 0){
            if(len[1] > 0 && (svl[i].len < len[0] || svl[i].len > len[1]))
                continue;
        }
        svl[i].het = 0;
        svl[i].ok = 0;
        i++;
    }
    svl[i].chr[0] = '\0';
    
    free(line);
    fclose(sv_file);
    
    if(gb > 0){
        for(i=0;svl[i].chr[0]!='\0';i++){
            if(svl[i].ok == 1)
                continue;
            for(j=i+1;svl[j].chr[0]!='\0';j++){
                if(svl[i].type != svl[j].type)
                    break;
                if(strcmp(svl[i].chr, svl[j].chr) == 0){
                    if(abs(svl[i].start-svl[j].start) <= gb && abs(svl[i].stop-svl[j].stop) <= gb){
                        if(svl[i].start > svl[j].start)
                            svl[i].start = svl[j].start;
                        if(svl[i].stop < svl[j].stop)
                            svl[i].stop = svl[j].stop;
                        if(svl[i].type != 0)
                            svl[i].len = svl[i].stop - svl[i].start;
                        else if(svl[j].len > svl[i].len)
                            svl[i].len = svl[j].len;
                        svl[j].ok = 1;
                    }
                    else if(svl[i].stop < svl[j].start)
                        break;
                }
                else if(strcmp(svl[i].chr, svl[j].chr) < 0)
                    break;
            }
            if(len[1] > 0 && (svl[i].len < len[0] || svl[i].len > len[1]))
                continue;
            strcpy(svl_cp[k].chr, svl[i].chr);
            svl_cp[k].type = svl[i].type;
            svl_cp[k].start = svl[i].start;
            svl_cp[k].stop = svl[i].stop;
            svl_cp[k].len = svl[i].len;
            svl_cp[k].het = 0;
            svl_cp[k].ok = 0;
            k++;
        }
        svl_cp[k].chr[0] = '\0';
        *n = k;
        
        free(svl);
        
        return svl_cp;
    }
    else
        *n = i;
    
    return svl;
}

SV_s **defHet(SV_s **svp1, SV_s **svp2, int n, int m, double gb){
    
    int i, j, k, l, ok=0;
    SV_s **svp;
    
    if((svp=malloc(n*sizeof(SV_s*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<n;i++){
        if((svp[i]=malloc(m*2*sizeof(SV_s))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    for(i=0;i<n;i++){
        l = 0;
        for(j=0;svp1[i][j].chr[0]!='\0';j++){
            ok = 0;
            for(k=0;svp2[i][k].chr[0]!='\0';k++){
                if(strcmp(svp1[i][j].chr, svp2[i][k].chr) == 0){
                    if(abs(svp1[i][j].start-svp2[i][k].start) <= gb && abs(svp1[i][j].stop-svp2[i][k].stop) <= gb){
                        if(svp1[i][j].type == svp2[i][k].type){
                            svp1[i][j].ok = 0;
                            if(svp1[i][j].start > svp2[i][k].start)
                                svp1[i][j].start = svp2[i][k].start;
                            if(svp1[i][j].stop < svp2[i][k].stop)
                                svp1[i][j].stop = svp2[i][k].stop;
                            if(svp1[i][j].type != 0)
                                svp1[i][j].len = svp1[i][j].stop-svp1[i][j].start;
                            else if(svp2[i][k].len > svp1[i][j].len)
                                svp1[i][j].len = svp2[i][k].len;
                            svp2[i][k].ok = 1;
                            ok = 1;
                        }
                        else{
                            svp1[i][j].ok = 1;
                            svp2[i][k].ok = 1;
                        }
                    }
                }
            }
            if(svp1[i][j].ok ==1)
                continue;
            strcpy(svp[i][l].chr, svp1[i][j].chr);
            svp[i][l].type = svp1[i][j].type;
            svp[i][l].start = svp1[i][j].start;
            svp[i][l].stop = svp1[i][j].stop;
            svp[i][l].len = svp1[i][j].len;
            if(ok == 1)
                svp[i][l].het = 0;
            else
                svp[i][l].het = 1;
            svp[i][l].ok = 0;
            l++;
        }
        for(j=0;svp2[i][j].chr[0]!='\0';j++){
            if(svp2[i][j].ok == 1)
                continue;
            strcpy(svp[i][l].chr, svp2[i][j].chr);
            svp[i][l].type = svp2[i][j].type;
            svp[i][l].start = svp2[i][j].start;
            svp[i][l].stop = svp2[i][j].stop;
            svp[i][l].len = svp2[i][j].len;
            svp[i][l].het = 2;
            svp[i][l].ok = 0;
            l++;
        }
    }
    
    free(svp1);
    free(svp2);
    
    return svp;
}

void defGeno(SV_s **svp, int n, double gb){
    
    int i, j, k, l, r=0, geno_i=0, *genos;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"}, nucs[]="ACGT";
    
    srand(time(NULL));
    
    if((genos=malloc(n*sizeof(int))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<n;i++)
        genos[i] = 0;
    
    for(i=0;i<n;i++){
        for(j=0;svp[i][j].chr[0]!='\0';j++){
            if(svp[i][j].ok == 1)
                continue;
            if(svp[i][j].het == 1)
                genos[i] = 1;
            else if(svp[i][j].het == 2)
                genos[i] = 2;
            else
                genos[i] = 3;
            for(k=i+1;k<n;k++){
                for(l=0;svp[k][l].chr[0]!='\0';l++){
                    if(svp[k][l].ok == 1)
                        continue;
                    if(strcmp(svp[i][j].chr, svp[k][l].chr) == 0){
                        if(abs(svp[i][j].start-svp[k][l].start) <= gb && abs(svp[i][j].stop-svp[k][l].stop) <= gb){
                            if(svp[i][j].type != svp[k][l].type)
                                continue;
                            if(svp[i][j].start > svp[k][l].start)
                                svp[i][j].start = svp[k][l].start;
                            if(svp[i][j].stop < svp[k][l].stop)
                                svp[i][j].stop = svp[k][l].stop;
                            if(svp[i][j].type != 0)
                                svp[i][j].len = svp[i][j].stop-svp[i][j].start;
                            else if(svp[k][l].len > svp[i][j].len)
                                svp[i][j].len = svp[k][l].len;
                            if(svp[k][l].het == 1)
                                genos[k] = 1;
                            else if(svp[k][l].het == 2)
                                genos[k] = 2;
                            else
                                genos[k] = 3;
                            svp[k][l].ok = 1;
                        }
                    }
                }
            }
            r = rand() % 4;
            printf("%s\t%i\t.\t%c\t<%s>\t.\tPASS\tSVTYPE=%s;END=%i;LEN=%i\tGT",svp[i][j].chr,svp[i][j].start,nucs[r],vars[svp[i][j].type],vars[svp[i][j].type],svp[i][j].stop,svp[i][j].len);
            for(k=0;k<n;k++){
                if(genos[k] == 1)
                    printf("\t1|0");
                else if(genos[k] == 2)
                    printf("\t0|1");
                else if(genos[k] == 3)
                    printf("\t1|1");
                else
                    printf("\t0|0");
                genos[k] = 0;
            }
            printf("\n");
        }
    }
    
    free(svp);
    free(genos);
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

void printInfo(void){
    
}
