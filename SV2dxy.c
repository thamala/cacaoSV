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

 Program for estimating Dxy and Fst between the major and minor SV arrangements.

 Compiling: gcc SV2dxy.c -o SV2dxy -lm

 Usage:
 -sv [file] A VCF file produced by mumco2vcf.
 -snp [file] SNP-vcf file.
 -t [str] Either INS, DEL, DUP, TRA, or INV.
 -l [int]-[int] Minimum and maximum length of the SVs.

 Example:
 ./SV2dxy -sv SV.vcf -snp SNP.vcf -t INV -l 1e4-1e6 > out.txt
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
    char chr[101];
}SV_s;

typedef struct{
    int pos;
    char chr[101];
}Site_s;

typedef struct{
    int start, stop;
    char chr[101];
}Gene_s;

void openFiles(int argc, char *argv[]);
Site_s *readSites(FILE *site_file, int *n);
Gene_s *readGenes(FILE *gene_file, int *n);
SV_s *readSV(FILE *sv_file, int *n, int *m, int type, double len[]);
void readSNP(FILE *snp_file, SV_s *svl, Site_s *sites, Gene_s *genes, int sv_n, int sv_m, int site_n, int gene_n);
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
    
    int i, j, type=-1, sv_n=0, sv_m=0, site_n=0, gene_n=0;
    double len[2]={0};
    char *temp1=NULL, *temp2=NULL, *sv_name, *snp_name, *site_name=NULL, *gene_name=NULL;
    const char *vars[5]={"INS","DEL","DUP","TRA","INV"};
    SV_s *svl;
    Site_s *sites;
    Gene_s *genes;
    FILE *sv_file=NULL, *snp_file=NULL, *site_file=NULL, *gene_file=NULL;
    
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
        
        else if(strcmp(argv[i], "-sites") == 0){
            if((site_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sites %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-genes") == 0){
            if((gene_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-genes %s\n", argv[i]);
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
    
    if(sv_file == NULL || snp_file == NULL){
        fprintf(stderr,"\nERROR: The following are required: -sv [file] -snp [file]\n\n");
        exit(EXIT_FAILURE);
    }
    
    if(site_name != NULL)
        sites = readSites(site_name, &site_n);
    if(gene_name != NULL)
        genes = readGenes(gene_name, &gene_n);
    svl = readSV(sv_name, &sv_n, &sv_m, type, len);
    readSNP(snp_name, svl, sites, genes, sv_n, sv_m, site_n, gene_n);
}

Site_s *readSites(FILE *site_file, int *n){
    
    int i, char_i=0, line_i=0, maxchar=0;
    char c, *line;
    Site_s *list;
    
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
    char c, *line;
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
        strncpy(list[*n].chr, strtok(line,"\t"), 100);
        list[*n].start = atoi(strtok(NULL,"\t"));
        list[*n].stop = atoi(strtok(NULL,"\t"));
        *n = *n + 1;
    }
    
    free(line);
    fclose(gene_file);
    
    return list;
}

SV_s *readSV(FILE *sv_file, int *n, int *m, int type, double len[]){
    
    int i, j, k, char_i=0, maxchar=0, tab_i=0, maxtab=0, line_i=0, id_s=0, ok=0, *info;
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
        if(temp[0] == '#')
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
                if(p/((double)j*2) >= 0.1 && p/((double)j*2) <= 0.9)
                    *m = *m + 1;
            }
        }
    }
    
    free(line);
    fclose(sv_file);
    
    return svl;
}

void readSNP(FILE *snp_file, SV_s *svl, Site_s *sites, Gene_s *genes, int sv_n, int sv_m, int site_n, int gene_n){
    
    int i, j, k=0, l=0, pos=0, ok=0, sv=-1, geno[3][2];
    double n1=0, n2=0, p1=0, p2=0, hw=0, hb=0, fst=0, hw_tot=0, hb_tot=0, hw_temp=0, hb_temp=0;
    double dxy=0, dxy_i=0, dxy_tot=0, dxy_tot_i=0;
    char chr[20], *line=NULL, *temp=NULL;
    size_t len=0;
    ssize_t read;
    
    printf("Fst\tDxy\n");
    
    while((read = getline(&line, &len, snp_file)) != -1){
        lineTerminator(line);
        temp = strtok(line,"\t");
        ok = 0;
        i = 0;
        j = 1;
        n1 = 0;
        n2 = 0;
        memset(geno,0,sizeof(geno[0][0])*3*2);
        if(temp[0] != '#'){
            strncpy(chr, temp, 19);
            while(temp != NULL){
                if(j == 2){
                    pos =  atoi(temp);
                    while(k < sv_m){
                        if(strcmp(chr,svl[k].chr) == 0){
                            if(pos >= svl[k].start && pos <= svl[k].stop){
                                if(sv == -1)
                                    sv = k;
                                else if(sv != k){
                                    hw_tot += hw;
                                    hb_tot += hb;
                                    dxy_tot += dxy;
                                    dxy_tot_i += dxy_i;
                                    fst = hw/hb;
                                    dxy /= dxy_i;
                                    if(isnan(fst) == 0 && isnan(dxy) == 0)
                                        printf("%f\t%f\n", fst, dxy);
                                    hw = 0;
                                    hb = 0;
                                    dxy = 0;
                                    dxy_i = 0;
                                    sv = k;
                                }
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
                else if(j > 9){
                    if(svl[k].geno[i] != 9){
                        if(svl[k].geno[i] == 0)
                            n1 += 2;
                        else if(svl[k].geno[i] == 2)
                            n2 += 2;
                        if(temp[0] == '0' && temp[2] == '0')
                            geno[svl[k].geno[i]][0] += 2;
                        else if(temp[0] == '1' && temp[2] == '0'){
                            geno[svl[k].geno[i]][0]++;
                            geno[svl[k].geno[i]][1]++;
                        }
                        else if(temp[0] == '0' && temp[2] == '1'){
                            geno[svl[k].geno[i]][0]++;
                            geno[svl[k].geno[i]][1]++;
                        }
                        else if(temp[0] == '1' && temp[2] == '1')
                            geno[svl[k].geno[i]][1] += 2;
                        i++;
                    }
                    else
                        i++;
                }
                temp = strtok(NULL,"\t");
                j++;
                if(temp == NULL){
                    if(n1 < 6 || n2 < 6)
                        break;
                    if(geno[0][0]+geno[2][0] > geno[0][1]+geno[2][1]){
                        p1 = geno[0][1] / n1;
                        p2 = geno[2][1] / n2;
                    }
                    else{
                        p1 = geno[0][0] / n1;
                        p2 = geno[2][0] / n2;
                    }
                    hw_temp = (p1-p2)*(p1-p2)-p1*(1-p1)/(n1-1)-p2*(1-p2)/(n2-1);
                    hb_temp = p1*(1-p2)+p2*(1-p1);
                    if(isnan(hb_temp) == 0){
                        dxy += hb_temp;
                        dxy_i++;
                    }
                        
                    if(isnan(hw_temp) == 0 && isnan(hb_temp) == 0){
                        hw += hw_temp;
                        hb += hb_temp;
                    }
                }
            }
        }
    }
    
    hw_tot += hw;
    hb_tot += hb;
    dxy_tot += dxy;
    dxy_tot_i += dxy_i;
    fst = hw/hb;
    dxy /= dxy_i;
    if(isnan(fst) == 0 && isnan(dxy) == 0)
        printf("%f\t%f\n", fst, dxy);
    
    fprintf(stderr,"Hw = %f\nHb = %f\nWeighted Fst = %f\n\n", hw_tot, hb_tot, hw_tot/hb_tot);
    fprintf(stderr,"Dxy = %f\nSites = %f\nMean Dxy = %f\n", dxy_tot, dxy_tot_i, dxy_tot/dxy_tot_i);
    
    if(site_n > 0)
        free(sites);
    free(line);
    free(svl);
    fclose(snp_file);
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
