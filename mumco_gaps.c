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

 Program for filtering SVs calls by MUM&co depending on their gap content

 Compiling: gcc mumco_gaps.c -o mumco_gaps -lm

 Usage:
 -sv [file] SV calls produced by MUM&Co
 -ref [file] fasta file for the reference assembly
 -que [file] fasta file for the query assembly
 -chr [int] length of the longest chromosome (or some bigger value)
 -gap [double] only keep SVs with lower propotion of gaps

 Examples:
 ./mumco_gaps -sv SV.tsv -ref ref.fa -que que.fa -chr 41236440 -gap 0.05 > out.tsv
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
    int ref_start, ref_end, que_start, que_end, size;
    double gaps;
    char ref_chr[50], que_chr[50], type[50];
}SV_s;

void openFiles(int argc, char *argv[]);
SV_s *readSV(FILE *sv_file, int *n);
void readFa(FILE *fa_file, SV_s *svl, int sv_n, int chr_l, int ref);
void countGaps(char nuc[], char chr[], int nuc_n, SV_s *svl, int sv_n, int ref);
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
    
    int i, sv_n=0, chr_l=0;
    double gap=0.1;
    SV_s *svl;
    FILE *sv_file=NULL, *ref_file=NULL, *que_file=NULL;
    
    fprintf(stderr,"\nParameters:\n");
    
    for(i=1;i<argc;i++){
        
        if(strcmp(argv[i], "-sv") == 0){
            if((sv_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-sv %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-ref") == 0){
            if((ref_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-ref %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-que") == 0){
            if((que_file = fopen(argv[++i], "r")) == NULL){
                fprintf(stderr,"\nERROR: Cannot open file %s\n\n", argv[i]);
                exit(EXIT_FAILURE);
            }
            fprintf(stderr,"\t-que %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-chr") == 0){
            chr_l = atoi(argv[++i]);
            fprintf(stderr, "\t-chr %s\n", argv[i]);
        }
        
        else if(strcmp(argv[i], "-gap") == 0){
            gap = atof(argv[++i]);
            fprintf(stderr, "\t-gap %s\n", argv[i]);
        }
        
        else{
            fprintf(stderr,"\nERROR: Unknown argument '%s'\n\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }
    
    fprintf(stderr,"\n");
    
    if(sv_file == NULL || ref_file == NULL || que_file == NULL || chr_l == 0){
        fprintf(stderr,"\nERROR: The following are required: -sv [file] -ref [file] -que [file] -chr [int]\n\n");
        exit(EXIT_FAILURE);
    }
    
    svl = readSV(sv_file, &sv_n);
    readFa(ref_file, svl, sv_n, chr_l, 1);
    readFa(que_file, svl, sv_n, chr_l, 0);
    
    printf("ref_chr\tquery_chr\tref_start\tref_stop\tsize\tSV_type\tquery_start\tquery_stop\n");
    for(i=0;i<sv_n;i++){
        if(strcmp(svl[i].type, "deletion") == 0)
            svl[i].gaps /= svl[i].ref_end - svl[i].ref_start + 1;
        else
            svl[i].gaps /= svl[i].que_end - svl[i].que_start + 1;
        if(svl[i].gaps < gap)
            printf("%s\t%s\t%i\t%i\t%i\t%s\t%i\t%i\n", svl[i].ref_chr, svl[i].que_chr, svl[i].ref_start, svl[i].ref_end, svl[i].size, svl[i].type, svl[i].que_start, svl[i].que_end);
    }
}

SV_s *readSV(FILE *sv_file, int *n){
    
    int i, char_i=0, maxchar=0, line_i=0, sv_l=0;
    char c, *line=NULL, *temp=NULL, *chr=NULL;
    SV_s *list;
    
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
    if((list=malloc((line_i+1)*sizeof(SV_s))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while(fgets(line, maxchar+1, sv_file) != NULL){
        if(line[0] == '\n' || line[0] == '\r')
            continue;
        lineTerminator(line);
        strncpy(list[*n].ref_chr, strtok(line,"\t"), 49);
        if(strcmp(list[*n].ref_chr, "ref_chr") == 0)
            continue;
        strncpy(list[*n].que_chr, strtok(NULL,"\t"), 49);
        list[*n].ref_start = atoi(strtok(NULL, "\t"));
        list[*n].ref_end = atoi(strtok(NULL, "\t"));
        list[*n].size = atoi(strtok(NULL, "\t"));
        strncpy(list[*n].type, strtok(NULL,"\t"), 49);
        if(strcmp(list[*n].type, "insertion") == 0 && list[*n].ref_end-list[*n].ref_start > 100000)
            continue;
        list[*n].que_start = atoi(strtok(NULL, "\t"));
        list[*n].que_end = atoi(strtok(NULL, "\t"));
        temp = strtok(NULL,"\t");
        if(temp != NULL){
            if(strcmp(temp,"complicated") == 0 || strcmp(temp,"double") == 0)
                continue;
        }
        list[*n].gaps = 0;
        *n = *n + 1;
    }
    
    free(line);
    fclose(sv_file);
    
    return list;
}

void readFa(FILE *fa_file, SV_s *svl, int sv_n, int chr_l, int ref){
    
    int i=0, j=0, k=0, ok=0;
    char chr[20]={0}, *line=NULL, *temp=NULL, *nuc=NULL;
    size_t len=0;
    ssize_t read;
    
    if((nuc=malloc((chr_l+1)*sizeof(char))) == NULL){
        fprintf(stderr,merror);
        exit(EXIT_FAILURE);
    }
    
    while((read = getline(&line, &len, fa_file)) != -1){
        if(line[0] == '\n' || line[0] == '\r')
            continue;
        else if(line[0] == '>'){
            lineTerminator(line);
            temp = strtok(line,">");
            if(chr[0] == '\0')
                strncpy(chr, temp, 19);
            else{
                countGaps(nuc, chr, j, svl, sv_n, ref);
                strncpy(chr, temp, 19);
                j = 0;
            }
            continue;
        }
        else{
            for(i=0;line[i]!=0;i++){
                if(line[i] == '\n' || line[i] == '\r')
                    continue;
                nuc[j] = line[i];
                j++;
                if(j == chr_l){
                    fprintf(stderr,"ERROR: fasta file contains scaffolds longer than %i char\n\n", chr_l);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    
    free(nuc);
    free(line);
    fclose(fa_file);
}

void countGaps(char nuc[], char chr[], int nuc_n, SV_s *svl, int sv_n, int ref){
    
    int i, j;
    
    if(ref == 1){
        for(i=0;i<sv_n;i++){
            if(strcmp(svl[i].type, "deletion") != 0)
                continue;
            if(strcmp(chr, svl[i].ref_chr) != 0)
                continue;
            for(j=svl[i].ref_start-1;j<svl[i].ref_end;j++){
                if(nuc[j] == 'N')
                    svl[i].gaps++;
            }
        }
    }
    else{
        for(i=0;i<sv_n;i++){
            if(strcmp(svl[i].type, "deletion") == 0)
                continue;
            if(strcmp(chr, svl[i].que_chr) != 0)
                continue;
            for(j=svl[i].que_start-1;j<svl[i].que_end;j++){
                if(nuc[j] == 'N')
                    svl[i].gaps++;
            }
        }
    }
}

void lineTerminator(char *line){
    
    int i;
    
    for(i=0;line[i]!=0;i++){
        if(line[i] == '\n' || line[i] == '\r')
            line[i] = '\0';
    }
}
