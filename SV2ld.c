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

 Program for estimating LD within major and minor SV arrangements.

 Compiling: gcc SV2ld.c -o SV2ld -lm

 Usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -snp [file] SNP-vcf file.
 -mat [file] Kinship matrix.
 -t [str] Either INS, DEL, DUP, TRA, or INV.
 -l [int]-[int] Minimum and maximum length of the SVs.

 Example:
 ./SV2ld -sv SV.vcf -snp SNP.vcf -mat kinshp.covmat -t INV -l 1e4-1e6 > out.txt
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
    int type, start, stop, len, snp_n, sv_n, *geno, *pos, **snps;
    double p;
    char chr[101];
}SV_s;

typedef struct{
    int n;
    double r;
}Cor_s;

void openFiles(int argc, char *argv[]);
SV_s *readSV(FILE *sv_file, int *n, int *m, int *maxlen, int type, double len[]);
double **readMat(FILE *mat_file, int *n);
void readSNP(FILE *snp_file, SV_s *svl, int sv_n, int sv_m, int maxlen);
void compSites(SV_s *svl, double **mat, int sv_n, int sv_m);
Cor_s estR(double snp1[], double snp2[], double **V, int n);
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
    
    int i, j, type=-1, sv_n=0, sv_m=0, mat_n=0, maxlen=0;
    double len[2]={0}, **mat;
    char *temp1=NULL, *temp2=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    SV_s *svl;
    FILE *sv_file=NULL, *snp_file=NULL, *mat_file=NULL;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((sv_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sv %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-snp") == 0){
            if((snp_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-snp %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-mat") == 0){
            if((mat_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-mat %s\n", argv[i]);
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
    
    if(sv_file == NULL || snp_file == NULL || mat_file == NULL){
        fprintf(stderr,"\nERROR: The following are required: -sv [file] -snp [file] -mat [file]\n\n");
        exit(EXIT_FAILURE);
    }
    
    svl = readSV(sv_file, &sv_n, &sv_m, &maxlen, type, len);
    mat = readMat(mat_file, &mat_n);
    if(sv_m != mat_n){
        fprintf(stderr,"\nERROR: SV-file and matrix-file have different number of individuals\n\n");
        exit(EXIT_FAILURE);
    }
    
    readSNP(snp_file, svl, sv_n, sv_m, maxlen);
    compSites(svl, mat, sv_n, sv_m);
}

SV_s *readSV(FILE *sv_file, int *n, int *m, int *maxlen, int type, double len[]){
    
    int i, j, k, char_i=0, maxchar=0, tab_i=0, maxtab=0, line_i=0, id_s=0, ok=0, ref=0, alt=0, *info;
    double p=0;
    char c, *temp, *line, *end;
    SV_s *svl;
    
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
    
    while(fgets(line, maxchar+1, sv_file) != NULL){
        lineTerminator(line);
        temp = strtok_r(line,"\t",&end);
        i = 1;
        j = 0;
        p = 0;
        ref=0;
        alt=0;
        if(temp[0] == '#')
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
                if(svl[*n].geno[j] != 9){
                    p += svl[*n].geno[j];
                    if(svl[*n].geno[j] == 0)
                        ref++;
                    if(svl[*n].geno[j] == 2)
                        alt++;
                }
                j++;
            }
            temp = strtok_r(NULL,"\t",&end);
            i++;
            if(temp == NULL){
                if(*m == 0)
                    *m = j;
                p /= (double)j*2;
                if(p >= 0.15 && p <= 0.85){
                    if(p > 0.5){
                        for(i=0;i<j;i++){
                            if(svl[*n].geno[i] == 0)
                                svl[*n].geno[i] = 2;
                            else if(svl[*n].geno[i] == 2)
                                svl[*n].geno[i] = 0;
                        }
                        svl[*n].sv_n = ref;
                    }
                    else
                        svl[*n].sv_n = alt;
                    if(svl[*n].len > *maxlen)
                        *maxlen = svl[*n].len;
                    svl[*n].snp_n = 0;
                    *n = *n + 1;
                }
            }
        }
    }
    
    free(line);
    fclose(sv_file);
    
    return svl;
}

double **readMat(FILE *mat_file, int *n){
    
    int i, j, char_i=0, line_i=0, maxchar=0;
    double **mat;
    char c, *line=NULL, *temp=NULL;
    
    while((c=fgetc(mat_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(mat_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((mat=malloc(line_i*sizeof(double*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<line_i;i++){
        if((mat[i]=malloc(line_i*sizeof(double))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    while(fgets(line, maxchar+1, mat_file) != NULL){
        lineTerminator(line);
        temp = strtok(line,"\t");
        j = 0;
        while(temp != NULL){
            mat[*n][j] = atof(temp);
            temp = strtok(NULL,"\t");
            j++;
        }
        *n = *n + 1;
    }
    
    fclose(mat_file);
    
    return mat;
}

void readSNP(FILE *snp_file, SV_s *svl, int sv_n, int sv_m, int maxlen){
    
    int i, j, k, pos=0, sv_i=0, ok=0, *geno;
    double ac=0, at=0;
    char chr[20], *line=NULL, *temp=NULL;
    size_t len=0;
    ssize_t read;
    
    for(i=0;i<sv_n;i++){
        if((svl[i].pos=malloc(maxlen*sizeof(int))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        if((svl[i].snps=malloc(maxlen*sizeof(int*))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        for(j=0;j<maxlen;j++){
            if((svl[i].snps[j]=malloc(sv_m*sizeof(int))) == NULL){
                fprintf(stderr,merror);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    if((geno=malloc(sv_m*sizeof(int))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, snp_file)) != -1){
        if(line[0] == '\n' || line[0] == '#')
            continue;
        lineTerminator(line);
        temp = strtok(line,"\t");
        j = 1;
        k = -1;
        ac = 0;
        at = 0;
        strncpy(chr, temp, 19);
        while(temp != NULL){
            if(j == 2){
                pos =  atoi(temp);
                ok = 0;
                while(sv_i < sv_n){
                    if(strcmp(chr, svl[sv_i].chr) == 0){
                        if(pos <= svl[sv_i].stop && pos >= svl[sv_i].start){
                            ok = 1;
                            break;
                        }
                        else if(pos < svl[sv_i].start)
                            break;
                    }
                    else if(strcmp(chr, svl[sv_i].chr) < 0)
                        break;
                    sv_i++;
                }
                if(ok == 0)
                    break;
            }
            else if(j > 9){
                if(temp[0] == '0' && temp[2] == '0')
                    geno[++k] = 0;
                else if(temp[0] == '1' && temp[2] == '0')
                    geno[++k] = 1;
                else if(temp[0] == '0' && temp[2] == '1')
                   geno[++k] = 1;
                else if(temp[0] == '1' && temp[2] == '1')
                    geno[++k] = 2;
                else
                    geno[++k] = 9;
                if(geno[k] != 9){
                    ac += geno[k];
                    at += 2;
                }
            }
            temp = strtok(NULL,"\t");
            j++;
            if(temp == NULL){
                if(ac / at > 0.05 && ac / at < 0.95){
                    for(i=0;i<sv_n;i++){
                        if(strcmp(chr, svl[i].chr) == 0){
                            if(pos >= svl[i].start && pos <= svl[i].stop){
                                svl[i].pos[svl[i].snp_n] = pos;
                                for(j=0;j<sv_m;j++)
                                    svl[i].snps[svl[i].snp_n][j] = geno[j];
                                svl[i].snp_n++;
                            }
                            else if(pos < svl[i].start)
                                break;
                        }
                        else if(strcmp(chr, svl[i].chr) < 0)
                            break;
                    }
                }
            }
        }
    }

    free(line);
    fclose(snp_file);
}

void compSites(SV_s *svl, double **mat, int sv_n, int sv_m){
    
    int i, j, k, l, dist=0, c_n=0, s_n=0, c1_i=0, s1_i=0, c2_i=0, s2_i=0, *mat_c, *mat_s;
    double r2_c=0, r2_s=0, p=0, *snp1_c, *snp1_s, *snp2_c, *snp2_s, **kin_c, **kin_s;
    char *line=NULL, *temp=NULL;
    Cor_s r_c, r_s;
    size_t len=0;
    ssize_t read;
    
    if((snp1_c=malloc(sv_m*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((snp1_s=malloc(sv_m*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((snp2_c=malloc(sv_m*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((snp2_s=malloc(sv_m*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((mat_c=malloc(sv_m*sizeof(int))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((mat_s=malloc(sv_m*sizeof(int))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((kin_c=malloc(sv_m*sizeof(double*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((kin_s=malloc(sv_m*sizeof(double*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<sv_m;i++){
        if((kin_c[i]=malloc(sv_m*sizeof(double))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        if((kin_s[i]=malloc(sv_m*sizeof(double))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    
    printf("distance\tc.r\tc.n\ts.r\ts.n\n");
    
    for(i=0;i<sv_n;i++){
        if(svl[i].snp_n < 2)
            continue;
        for(j=0;j<svl[i].snp_n;j++){
            c1_i = -1;
            s1_i = -1;
            for(k=0;k<sv_m;k++){
                if(svl[i].geno[k] == 0 && c1_i < svl[i].sv_n)
                    snp1_c[++c1_i] = svl[i].snps[j][k];
                else if(svl[i].geno[k] == 2 && s1_i < svl[i].sv_n)
                    snp1_s[++s1_i] = svl[i].snps[j][k];
            }
            for(k=j+1;k<svl[i].snp_n;k++){
                dist = svl[i].pos[k] - svl[i].pos[j];
                c2_i = -1;
                s2_i = -1;
                for(l=0;l<sv_m;l++){
                    if(svl[i].geno[l] == 0 && c2_i < svl[i].sv_n){
                        snp2_c[++c2_i] = svl[i].snps[k][l];
                        mat_c[c2_i] = l;
                    }
                    else if(svl[i].geno[l] == 2 && s2_i < svl[i].sv_n){
                        snp2_s[++s2_i] = svl[i].snps[k][l];
                        mat_s[s2_i] = l;
                    }
                }
            }
            for(k=0;k<=c2_i;k++){
                for(l=0;l<=c2_i;l++)
                    kin_c[k][l] = mat[mat_c[k]][mat_c[l]];
            }
            for(k=0;k<=s2_i;k++){
                for(l=0;l<=s2_i;l++)
                    kin_s[k][l] = mat[mat_s[k]][mat_s[l]];
            }
            r_c = estR(snp1_c, snp2_c, kin_c, c2_i+1);
            r_s = estR(snp1_s, snp2_s, kin_s, s2_i+1);
            if(r_s.n > 2 && r_c.n > 2 && isnan(r_c.r) == 0 && isnan(r_s.r) == 0)
                printf("%i\t%f\t%i\t%f\t%i\n",dist,r_c.r,r_c.n,r_s.r,r_s.n);
        }
    }
}

Cor_s estR(double snp1[], double snp2[], double **V, int n){
    
    //From Mangin et al. 2012
    
    int i, j, k=0, m=0;
    double s1=0, s2=0, X1=0, X2=0, nom=0, denom1=0, denom2=0;
    Cor_s cor;
    
    m = n;
    
    for(i=0;i<m;i++){
        if(snp1[i] != 9 && snp2[i] != 9){
            s1 = 0;
            s2 = 0;
            for(j=0;j<m;j++){
                if(snp1[j] == 9 || snp2[j] == 9)
                    continue;
                s1 += V[i][j]*snp1[j];
                s2 += V[i][j]*snp2[j];
            }
            snp1[k] = s1;
            snp2[k] = s2;
            X1 += snp1[k];
            X2 += snp2[k];
            k++;
        }
        else
            n--;
    }
    X1 /= (double)n;
    X2 /= (double)n;
    for(i=0;i<n;i++){
        nom += (snp1[i]-X1)*(snp2[i]-X2);
        denom1 += (snp1[i]-X1)*(snp1[i]-X1);
        denom2 += (snp2[i]-X2)*(snp2[i]-X2);
    }
    cor.n = n;
    cor.r = fabs(nom/(sqrt(denom1)*sqrt(denom2)));
    return cor;
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
