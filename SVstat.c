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

 Program for computing basic statistics using SVs

 Compiling: gcc SVstat.c -o SVstat -lm

 Usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -fst [file] For estimating Fst. File listing individuals from a single population (one per line).
 -dxy[file] For estimating Dxy. File listing individuals from a single population (one per line).
 -sfs Produces a 1D allele frequency spectrum.
 -pca Estimates a kinship matrix for PCA (model by Patterson et al. 2006).
 -n [int] Number of permutation cycles for Fst or Dxy.
 -z Reports Z-values instead of P-values for permutation testing.
 -maf [double] Minimum minor allele frequency.
 -t [str] Either INS, DEL, DUP, TRA, or INV.
 -l [int]-[int] Minimum and maximum length of the SVs.
 -o [int] Name for the output file.

 Examples:
 ./SVstat -sv SV.vcf -fst pop1.txt -fst pop2.txt -n 100 -Z -sfs -pca -maf 0.05 -t DEL -l 500-1e5 -o test
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
    int type, start, stop, len, maj, *geno;
    char chr[101], **id;
}SV_s;

void openFiles(int argc, char *argv[]);
SV_s *readVCF(FILE *vcf_file, int *n, int *m, int type, double len[]);
char **readFst(FILE *fst_file);
void doSFS(SV_s *svl, char out[], int n, int m);
void doPCA(SV_s *svl, char out[], int n, int m, double maf);
void doFst(SV_s *svl, char **fst1, char **fst2, char out[], int n, int m, int zval, int dxy, double maf, double perm);
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
    
    int i, j, type=-1, sfs=0, pca=0, theta=0, zval=0, n=0, m=0, f_ok=0, dxy=0;
    double maf=0, perm=0, len[2]={0};
    char out[101]="SVstat", **fst1=NULL, **fst2=NULL, *temp1=NULL, *temp2=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    SV_s *svl;
    FILE *vcf_file=NULL, *fst_file=NULL;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((vcf_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sv %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-fst") == 0 || strcmp(argv[i], "-dxy") == 0){
            if(strcmp(argv[i], "-dxy") == 0)
                dxy = 1;
            if((fst_file=fopen(argv[++i],"r"))==NULL){
                fprintf(stderr,"\nERROR:Cannot open file %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            if(fst1 == NULL){
                fst1 = readFst(fst_file);
                fprintf(stderr, "\t-fst %s\n", argv[i]);
            }
            else if(fst2 == NULL){
                fst2 = readFst(fst_file);
                fprintf(stderr, "\t-fst %s\n", argv[i]);
            }
            else
                f_ok = 1;
        }
        
        else if(strcmp(argv[i], "-o") == 0){
            strncpy(out, argv[++i], 100);
            fprintf(stderr, "\t-o %s\n", out);
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
        
        else if(strcmp(argv[i], "-n") == 0){
            if(isNumeric(argv[++i]))
                perm = atof(argv[i]);
            fprintf(stderr, "\t-n %.0f\n", perm);
        }
        
        else if(strcmp(argv[i], "-z") == 0){
            zval = 1;
            fprintf(stderr, "\t-z\n");
        }
        
        else if(strcmp(argv[i], "-maf") == 0){
            if(isNumeric(argv[++i]))
                maf = atof(argv[i]);
            fprintf(stderr, "\t-maf %.2f\n", maf);
        }
        
        else if(strcmp(argv[i], "-sfs") == 0){
            sfs = 1;
            fprintf(stderr, "\t-sfs\n");
        }
        
        else if(strcmp(argv[i], "-pca") == 0){
            pca = 1;
            fprintf(stderr, "\t-pca\n");
        }
        
        else{
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    
    fprintf(stderr,"\n");
    
    if(vcf_file == NULL){
        fprintf(stderr,"\nERROR: -vcf [file] is needed\n");
        exit(EXIT_FAILURE);
    }
    
    if((fst1 != NULL && fst2 == NULL) || (fst1 == NULL && fst2 != NULL)){
        fprintf(stderr,"\nERROR: To estimate Fst, two lists are needed\n");
        exit(EXIT_FAILURE);
    }
    
    if(f_ok == 1)
        fprintf(stderr,"\nWarning: Only the first two Fst-lists are used\n");
    
    svl = readVCF(vcf_file, &n, &m, type, len);
    
    if(sfs == 1)
        doSFS(svl, out, n, m);
    
    if(pca == 1)
        doPCA(svl, out, n, m, maf);
    
    if(fst1 != NULL && fst2 != NULL){
        doFst(svl, fst1, fst2, out, n, m, zval, dxy, maf, perm);
        free(fst1);
        free(fst2);
    }
}

SV_s *readVCF(FILE *vcf_file, int *n, int *m, int type, double len[]){
    
    int i, j=0, char_i=0, maxchar=0, tab_i=0, maxtab=0, line_i=0, id_s=0, *info;
    double p=0;
    char c, *temp, *line, *end;
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
    
    while(fgets(line, maxchar+1, vcf_file) != NULL){
        if(line[0] == '\n')
            continue;
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
                svl[*m].len = info[2];
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
                if(*n == 0)
                    *n = j;
                if(p/((double)j*2) > 0.5)
                    svl[*m].maj = 1;
                else
                    svl[*m].maj = 0;
                *m = *m + 1;
            }
        }
    }
    
    free(line);
    
    return svl;
}

char **readFst(FILE *fst_file){
    
    int i=0, char_i=0, maxchar=0, line_i=0;
    char c, *line, **list;
    
    while((c=fgetc(fst_file)) != EOF){
        char_i++;
        if(c == '\n'){
            line_i++;
            if(char_i > maxchar)
                maxchar = char_i;
            char_i = 0;
        }
    }
    
    rewind(fst_file);
    
    if((line=malloc((maxchar+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    if((list=malloc((line_i+1)*sizeof(char*))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    for(i=0;i<line_i+1;i++){
        if((list[i]=malloc((maxchar+1)*sizeof(char))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
    }
    i = 0;
    while(fgets(line, maxchar+1, fst_file) != NULL){
        if(line[0] == '\n')
            continue;
        lineTerminator(line);
        strcpy(list[i], line);
        i++;
    }
    list[i][0] = '\0';
    
    free(line);
    fclose(fst_file);
    
    return list;
}

void doSFS(SV_s *svl, char out[], int n, int m){
    
    int i, j, count=0;
    double tot=0, *sfs;
    char name[60];
    FILE *out_file;
    
    n *= 2;
    
    sprintf(name,"%s.sfs", out);
    
    if((out_file=fopen(name,"w"))==NULL){
        fprintf(stderr,"\nERROR:Cannot create file %s\n", name);
        exit(EXIT_FAILURE);
    }
    
    if((sfs=malloc(n*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<n;i++)
        sfs[i] = 0;
    
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            if(svl[i].geno[j] == 1)
                count++;
            else if(svl[i].geno[j] == 2)
                count += 2;
        }
        sfs[count-1]++;
        count = 0;
    }
    
    for(i=0;i<n;i++){
        fprintf(out_file, "%.0f ", sfs[i]);
        tot += sfs[i];
    }
    fprintf(out_file,"\n");
    for(i=0;i<n;i++)
        fprintf(out_file, "%f ", sfs[i]/tot);
    fprintf(out_file,"\n");
    
    free(sfs);
    fclose(out_file);
}

void doPCA(SV_s *svl, char out[], int n, int m, double maf){
    
    int i, j, k;
    double g1=0, g2=0, sum=0, sum_i=0, *p;
    char name[60];
    FILE *out_file;
    
    sprintf(name,"%s.pca", out);
    
    if((out_file=fopen(name,"w"))==NULL){
        fprintf(stderr,"\nERROR:Cannot create file %s\n", name);
        exit(EXIT_FAILURE);
    }
    
    if((p=malloc(m*sizeof(double))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<m;i++){
        for(j=0;j<n;j++)
            p[i] += svl[i].geno[j];
        p[i] /= n*2;
    }
    
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            for(k=0;k<m;k++){
                if(p[k] == 1)
                    continue;
                else if(p[k] <= 0.5 && p[k] < maf)
                    continue;
                else if(p[k] > 0.5 && (1-p[k]) < maf)
                    continue;
                g1 = svl[k].geno[i];
                g2 = svl[k].geno[j];
                sum += (g1-2*p[k])*(g2-2*p[k])/(2*p[k]*(1-p[k]));
                sum_i++;
            }
            fprintf(out_file,"%f\t",sum/sum_i);
            sum=0;
            sum_i=0;
        }
        fprintf(out_file,"\n");
    }
    
    free(p);
    fclose(out_file);
}

void doFst(SV_s *svl, char **fst1, char **fst2, char out[], int n, int m, int zval, int dxy, double maf, double perm){
    
    int i, j, k, r_i=0, p1_i=0, p2_i=0, *plist1=NULL, *plist2=NULL, *rlist=NULL;
    double p=0, p1=0, p2=0, n1=0, n2=0, pi1=0, pi2=0, fst=0, hw=0, hb=0, r_p1=0, r_p2=0, r_fst=0, f_p=0, mean=0, sd=0, *flist=NULL;
    char name[60];
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    FILE *out_file;
    
    if(dxy == 1)
        sprintf(name,"%s.dxy", out);
    else
        sprintf(name,"%s.fst", out);
    
    if((out_file=fopen(name,"w"))==NULL){
        fprintf(stderr,"\nERROR:Cannot create file %s\n", name);
        exit(EXIT_FAILURE);
    }
    
    if((plist1=malloc(n*sizeof(int))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    if((plist2=malloc(n*sizeof(int))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    for(i=0;i<n;i++){
        for(j=0;fst1[j][0]!='\0';j++){
            if(strcmp(svl[0].id[i], fst1[j]) == 0){
                plist1[p1_i] = i;
                p1_i++;
            }
        }
        for(j=0;fst2[j][0]!='\0';j++){
            if(strcmp(svl[0].id[i], fst2[j]) == 0){
                plist2[p2_i] = i;
                p2_i++;
            }
        }
    }
    
    if(p1_i == 0 || p2_i == 0){
        fprintf(stderr,"ERROR: individuals in the Fst lists are not found in the VCF-file\n\n");
        exit(EXIT_FAILURE);
    }
    
    if(dxy == 1){
        if(perm > 0 && zval == 1)
            fprintf(out_file,"Chr\tType\tStart\tEnd\tLength\tD1\tD2\tDxy\tZ\tP1\tP2\n");
        else if(perm > 0 && zval == 0)
            fprintf(out_file,"Chr\tType\tStart\tEnd\tLength\tD1\tD2\tDxy\tP\tP1\tP2\n");
        else
            fprintf(out_file,"Chr\tType\tStart\tEnd\tLength\tD1\tD2\tDxy\tP1\tP2\n");
    }
    else{
        if(perm > 0 && zval == 1)
            fprintf(out_file,"Chr\tType\tStart\tEnd\tLength\tFst\tZ\tP1\tP2\n");
        else if(perm > 0 && zval == 0)
            fprintf(out_file,"Chr\tType\tStart\tEnd\tLength\tFst\tP\tP1\tP2\n");
        else
            fprintf(out_file,"Chr\tType\tStart\tEnd\tLength\tFst\tP1\tP2\n");
    }
     
    if(perm > 0){
        srand(time(NULL));
        if((rlist=malloc(n*sizeof(int))) == NULL){
            fprintf(stderr,merror);
            exit(EXIT_FAILURE);
        }
        if(zval == 1){
            if((flist=malloc(perm*sizeof(double))) == NULL){
                fprintf(stderr,merror);
                exit(EXIT_FAILURE);
            }
        }
    }
    
    for(i=0;i<m;i++){
        p1 = 0;
        p2 = 0;
        n1 = 0;
        n2 = 0;
        r_i = 0;
        f_p = 0;
        mean = 0;
        sd = 0;
        for(j=0;j<p1_i;j++){
            p1 += svl[i].geno[plist1[j]];
            n1 += 2;
            if(perm > 0){
                rlist[r_i] = svl[i].geno[plist1[j]];
                r_i++;
            }
        }
        for(j=0;j<p2_i;j++){
            p2 += svl[i].geno[plist2[j]];
            n2 += 2;
            if(perm > 0){
                rlist[r_i] = svl[i].geno[plist2[j]];
                r_i++;
            }
        }
        p = (p1+p2)/(n1+n2);
        if(p <= 0.5 && p < maf)
            continue;
        else if(p > 0.5 && (1-p) < maf)
            continue;
        if(dxy == 1){
            if(p1 > n1/2)
                pi1 = (2*(n1-p1)*(n1-1))/(n1*n1-1);
            else
                pi1 = (2*p1*(n1-1))/(n1*n1-1);
            if(p2 > n2/2)
                pi2 = (2*(n2-p2)*(n2-1))/(n2*n2-1);
            else
                pi2 = (2*p2*(n2-1))/(n2*n2-1);
            p1 /= n1;
            p2 /= n2;
            if(svl[i].maj == 1){
                p1 = 1-p1;
                p2 = 1-p2;
            }
            fst = p1*(1-p2)+p2*(1-p1);
        }
        else{
            p1 /= n1;
            p2 /= n2;
            if(svl[i].maj == 1){
                p1 = 1-p1;
                p2 = 1-p2;
            }
            fst = (pow((p1-p2),2)-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1))/(p1*(1-p2)+p2*(1-p1));
            hw += (pow((p1-p2),2)-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1));
            hb += (p1*(1-p2)+p2*(1-p1));
        }
        if(isnan(fst) || isinf(fst))
            continue;
        if(perm > 0){
            for(j=0;j<perm;j++){
                r_p1 = 0;
                r_p2 = 0;
                for(k=0;k<n1/2;k++)
                    r_p1 += rlist[rand()%r_i];
                for(k=0;k<n2/2;k++)
                    r_p2 += rlist[rand()%r_i];
                r_p1 /= n1;
                r_p2 /= n2;
                if(dxy == 1)
                    r_fst = r_p1*(1-r_p2)+r_p2*(1-r_p1);
                else
                    r_fst = (pow((r_p1-r_p2),2)-r_p1*(1-r_p1)/(n1-1)-r_p2*(1-r_p2)/(n2-1))/(r_p1*(1-r_p2)+r_p2*(1-r_p1));
                if(isnan(r_fst) || isinf(r_fst)){
                    j--;
                    continue;
                }
                if(zval == 1){
                    flist[j] = r_fst;
                    mean += r_fst;
                }
                else{
                    if(r_fst >= fst)
                        f_p++;
                }
            }
            if(zval == 1){
                mean /= (double)perm;
                for(j=0;j<perm;j++)
                    sd += pow((flist[j]-mean),2);
                sd = sqrt(sd/((double)perm-1));
                f_p = (fst-mean)/sd;
            }
            else
                f_p = (1+f_p)/(1+(double)perm);
            fprintf(out_file,"%s\t%s\t%i\t%i\t%i\t", svl[i].chr, vars[svl[i].type], svl[i].start, svl[i].stop, svl[i].len);
            if(dxy == 1)
                fprintf(out_file,"%f\t%f\t%f\t%f\t%f\t%f\n",pi1, pi2, fst, f_p, p1, p2);
            else
                fprintf(out_file,"%f\t%f\t%f\t%f\n",fst, f_p, p1, p2);
        }
        else{
            fprintf(out_file,"%s\t%s\t%i\t%i\t%i\t", svl[i].chr, vars[svl[i].type], svl[i].start, svl[i].stop, svl[i].len);
            if(dxy == 1)
                fprintf(out_file,"%f\t%f\t%f\t%f\t%f\n",pi1, pi2, fst, p1, p2);
            else
                fprintf(out_file,"%f\t%f\t%f\n",fst, p1, p2);
        }
    }
    
    if(dxy == 0){
        fprintf(stderr, "\nAverage weighted Fst = %f\n", hw/hb);
        if(isatty(fileno(stdout)) == 0)
            printf("%f\n", hw/hb);
    }
    
    
    if(perm > 0){
        free(rlist);
        if(zval == 1)
            free(flist);
    }
    fclose(out_file);
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
