
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


struct Coord{
	int i;
	int j;
};

struct TrueMismatch{
	char* perfStr;
	char* queyStr;
	int start_ind;
	int end_ind;
};

void mfree_TrueMismatch(struct TrueMismatch tm){
	free(tm.perfStr);
	tm.perfStr = NULL;
	free(tm.queyStr);
	tm.queyStr=NULL;
}

void m_init(char* marr, int msize){
	int i;
	for (i=0; i< msize; i++){
		marr[i] = '\0';
	}
}

struct TrueMismatch* getMismatch(char* mmismatch, int mislen, int mismnum){
	int i,j,k,prep,curp;
	struct TrueMismatch * mmismatchlist;
	char* numstr;

	if (mismnum>0){
		//printf("%s\n", mmismatch);
		numstr = malloc(sizeof(char)*(50));
		m_init(numstr, 50);

		mmismatchlist = malloc(sizeof(struct TrueMismatch) * (mismnum));
		prep=0; i=0; k = 0;
		for (curp=0; curp<mislen; curp++){
			//printf("%d(%c) ", curp, mmismatch[curp]);
			if (mmismatch[curp]==':' || curp==mislen-1){
				if (curp==mislen-1){
					j = curp - prep + 1;
				}else{
					j = curp - prep;
				}
				if (i==0){
					mmismatchlist[k].perfStr=malloc(sizeof(char)*(j+1));
					m_init(mmismatchlist[k].perfStr, j+1);
					strncpy((mmismatchlist[k].perfStr), &mmismatch[prep], j);
					mmismatchlist[k].perfStr[j] = '\0';
				}else if (i==1){
					mmismatchlist[k].queyStr=malloc(sizeof(char)*(j+1));
					m_init(mmismatchlist[k].queyStr, j+1);
					strncpy((mmismatchlist[k].queyStr), &mmismatch[prep], j);
					mmismatchlist[k].queyStr[j] = '\0';	
				}else if (i==2){
					strncpy(numstr, &mmismatch[prep], j);
					numstr[j] = '\0';
					mmismatchlist[k].start_ind = atoi(numstr);
				}else if (i==3){
					strncpy(numstr, &mmismatch[prep], j);
					numstr[j] = '\0';
					mmismatchlist[k].end_ind = atoi(numstr);
				}else{
					printf("Error too much field (%d) in %s (%d)\n", i, mmismatch, curp);
				}
				
				prep = curp + 1;
				i = i + 1;
			}else if (mmismatch[curp]==';'){
				j = curp - prep;
				if (i==3){
					strncpy(numstr, &mmismatch[prep], j);
					numstr[j] = '\0';
					mmismatchlist[k].end_ind = atoi(numstr);
				}else{
					printf("Error wrong field (%d) in %s (%d)\n", i, mmismatch, curp);
				}
				
				prep = curp + 1;
				i = 0;
				k = k + 1;
			}
		}
		
		//printf("free in getMismatch numstr");
		free(numstr); numstr=NULL;
	}else{
		mmismatchlist = NULL;
	}

	return mmismatchlist;
}


struct AlignRes {
	char* correctedRead;
	char* perfRep;
	int size;
	int score;
};

void myfree(int lenRead, int ** alignmentMatrix, struct Coord ** ij){
	int i;

        for (i=0; i<lenRead+1; i++){
                free(alignmentMatrix[i]);
                alignmentMatrix[i] = NULL;
                free(ij[i]);
                ij[i] = NULL;
        }

        //printf("free in myfree alignmentMatrix");
        free(alignmentMatrix);
        //printf("free in myfree ij");
        free(ij);
}

//match=10; mismatch=-10; gap_in_perf=-1; gap_in_read=-10;
struct AlignRes _unsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum){
	int i,ii,j,jj,k,kk,diagonal,up,left; //, hasacgt;
	int ** alignmentMatrix;
	struct Coord ** ij;
	struct Coord maxscoreij;
	int maxscore;
	struct AlignRes ar;
	int isspecial;
	struct TrueMismatch * mmismatchlist;

	ar.size = 0;
	maxscore = INT_MIN;
	maxscoreij.i = -1;
	maxscoreij.j = -1;
	
	ar.correctedRead = malloc(sizeof(char) * (lenPerfRep+lenRead+1));
	m_init(ar.correctedRead, lenPerfRep+lenRead+1);
	ar.perfRep = malloc(sizeof(char) * (lenPerfRep+lenRead+1));
	m_init(ar.perfRep, lenPerfRep+lenRead+1);
	
	alignmentMatrix = malloc(sizeof(int*) * (lenRead+1));	
	ij = malloc(sizeof(struct Coord*) * (lenRead+1));

	mmismatchlist = getMismatch(mmismatch, mislen, mismnum);
	for (i=0; i<mismnum; i++){
		if (isprint) printf("two string %s %s, %d, %d\n", mmismatchlist[i].perfStr, mmismatchlist[i].queyStr, mmismatchlist[i].start_ind, mmismatchlist[i].end_ind);
	}

	for (i=0; i<lenRead+1; i++){
		alignmentMatrix[i] = malloc(sizeof(int) * (lenPerfRep+1));
		ij[i] = malloc(sizeof(struct Coord) * (lenPerfRep+1));
 		
		for (j=0; j<lenPerfRep+1; j++){
			if (i==0 || j==0){ 
				if (abs(i-j)>bandWidth && bandWidth>10){
					alignmentMatrix[i][j] = INT_MIN;
				}else{
					alignmentMatrix[i][j] = (i+j)*gap_before_after;
				}
			}
			else{ alignmentMatrix[i][j] = INT_MIN;}
			if (i==0){ ij[i][j].j = j-1; ij[i][j].i = i;}
			if (j==0){ ij[i][j].i = i-1; ij[i][j].j = j;}
			if (i==0 && j==0){ ij[i][j].i=-1; ij[i][j].j=-1;}
		}
	}

	diagonal=0; up=0; left=0;
	for (i=1; i<lenRead+1; i++){
		for (j=1; j<lenPerfRep+1; j++){
			if (abs(i-j)>bandWidth && bandWidth>10){
				continue;
			}

			if (alignmentMatrix[i-1][j-1]==INT_MIN){
				diagonal = INT_MIN;
			}else{
				diagonal = alignmentMatrix[i-1][j-1];
				if (perfectRepeat[j-1]==read[i-1]){ diagonal += match;}
				//else{ diagonal += mismatch;}
				else {
					isspecial = 0;
					for (k=0; k<mismnum; k++){
						kk=1;
						for (ii=mmismatchlist[k].start_ind, jj=0; ii<=mmismatchlist[k].end_ind; ii++, jj++){
							if (j-1+ii>-1 && j-1+ii<lenPerfRep && mmismatchlist[k].perfStr[jj]==perfectRepeat[j-1+ii] && i-1+ii>-1 && i-1+ii<lenRead && mmismatchlist[k].queyStr[jj]==read[i-1+ii] ){}
							else {kk=-1; break;}
						}
						if (kk==1){
							isspecial = 1; 
							for (ii=mmismatchlist[k].start_ind, jj=0; ii<=mmismatchlist[k].end_ind; ii++, jj++){
								//if (isprint) printf("%d %d %c:%c\n", ii, jj, perfectRepeat[j-1+ii], read[i-1+ii]);
							}
							break;
						}
					}
					if(isspecial==1){
						diagonal += match;
						//diagonal += mismatch;
					}else{ diagonal += mismatch;}
				}
			}

			if (alignmentMatrix[i-1][j]==INT_MIN){
				up = INT_MIN;
			}else{
				up = alignmentMatrix[i-1][j];
				up += ((j==lenPerfRep)?gap_before_after:gap_in_perf);
			}

			if (alignmentMatrix[i][j-1]==INT_MIN){
				left = INT_MIN;
			}else{
				left = alignmentMatrix[i][j-1];
				left += ((i==lenRead)?gap_before_after:gap_in_read);
			}

			/*
			diagonal = alignmentMatrix[i-1][j-1];
			up = alignmentMatrix[i-1][j];
			left = alignmentMatrix[i][j-1];
			
			if (perfectRepeat[j-1]==read[i-1]){ diagonal += match;}
			else{ diagonal += mismatch;}
			
			up += ((j==lenPerfRep)?gap_before_after:gap_in_perf);
			left += ((i==lenRead)?gap_before_after:gap_in_read);
			*/

			if (up >= diagonal && up >= left){
				alignmentMatrix[i][j] = up;
				ij[i][j].i = i-1; ij[i][j].j = j;
				
				//if (up>=maxscore){
					maxscore = up;
					maxscoreij.i = i; maxscoreij.j = j;
				//}
			}else 
			  if (diagonal >= up && diagonal >= left){
				alignmentMatrix[i][j] = diagonal;
				ij[i][j].i = i-1; ij[i][j].j = j-1;
				
				//if (diagonal>=maxscore){
					maxscore = diagonal;
					maxscoreij.i = i; maxscoreij.j = j;
				//} 
			}else 
			  /*if (up >= diagonal && up >= left){
				alignmentMatrix[i][j] = up;
				ij[i][j].i = i-1; ij[i][j].j = j;
			}else*/ 
			  if (left >= diagonal && left >= up){
				alignmentMatrix[i][j] = left;
				ij[i][j].i = i; ij[i][j].j = j-1;
				
				//if (left>=maxscore){
					maxscore = left;
					maxscoreij.i = i; maxscoreij.j = j;
				//}
			}else{
				printf("Error no largest vaues among (diagonal=%d, up=%d, left=%d) at (row=%d, column=%d)", diagonal, up, left, i,j);
			}
			//printf("diagonal=%d(%d), up=%d(%d) left=%d(%d) at (%d, %d) and maxscore=%d at ij=(%d, %d) \n", diagonal, alignmentMatrix[i-1][j-1], up, alignmentMatrix[i-1][j], left, alignmentMatrix[i][j-1], i, j, maxscore, maxscoreij.i, maxscoreij.j);
		}
	}

	//k = lenPerfRep+lenRead-1;
	k = lenPerfRep+lenRead;
	//printf("in:%d, %d, %d, %d\n", lenPerfRep, lenRead, lenPerfRep+lenRead+1, k);
	ar.correctedRead[k] = '\0'; ar.perfRep[k]='\0'; k--;
 	i = lenRead; j = lenPerfRep;
	ar.score = alignmentMatrix[i][j];
	//printf("in:%d, %d, %d, %d, %d, %d, %d\n", lenPerfRep, lenRead, lenPerfRep+lenRead+1, k, i, j, ar.score);
	
	i = maxscoreij.i; j = maxscoreij.j;
	ar.score = maxscore;

	//printf("%s\n%s\n", perfectRepeat, read);

	while (i>0 || j>0){
		if (alignmentMatrix[i][j]==INT_MIN){
			printf("Errors: current value at (%d, %d) is INT_MIN. Matrix is (%d, %d)\n", i, j, lenRead, lenPerfRep);
			printf("Errors: gap_before_after=%d, bandWidth=%d\n", gap_before_after, bandWidth);
			printf("Errors: maxscore at (%d, %d) is %d; INT_MIN+2=%d\n", maxscoreij.i, maxscoreij.j, maxscore, INT_MIN+2);
			myfree(lenRead, alignmentMatrix, ij);
			free(ar.correctedRead);
			free(ar.perfRep);

			exit(11);
		}
		
		if (ij[i][j].i==i-1 && ij[i][j].j == j-1){
			//printf("\t%d,%d,%c, %c, %d\n",i,j,read[i-1], perfectRepeat[j-1],k);
			ar.correctedRead[k] = read[i-1]; ar.perfRep[k]=perfectRepeat[j-1]; k--;
			i = i-1; j = j-1;
		}
		else if (ij[i][j].i==i-1 && ij[i][j].j == j){
			//printf("\t%d,%d,%c, %c, %d\n",i,j,read[i-1], '-',k);
			ar.correctedRead[k] = read[i-1]; ar.perfRep[k]='-'; k--;
			i = i-1; //j = j;
		}
		else if (ij[i][j].i==i && ij[i][j].j == j-1){
			//printf("\t%d,%d,%c, %c, %d\n",i,j,'-',perfectRepeat[j-1],k);
			ar.correctedRead[k] = '-'; ar.perfRep[k]=perfectRepeat[j-1]; k--;
			j = j-1; //i = i;
		}
	}
	//printf("%s\n%s %d\n", ar.perfRep, ar.correctedRead, k);
	i = 0; k++; 
	for (j=k; j<lenPerfRep+lenRead; j++){
		//printf("%d, %d, %d, %d, %c, %c\n", i, j, k, lenPerfRep+lenRead, ar.correctedRead[j], ar.perfRep[j]);
		ar.correctedRead[i] = ar.correctedRead[j];
		ar.perfRep[i] = ar.perfRep[j];

		if ((ar.correctedRead[i]=='\0' && ar.perfRep[i]!='\0') || (ar.correctedRead[i]!='\0' && ar.perfRep[i]=='\0')){
			printf("Errors: the size after alignment is not equal: \n%s, \n%s\n", ar.correctedRead, ar.perfRep);
		}

		i++;
	}
	//printf("%s\n%s\n", ar.perfRep, ar.correctedRead);

	ar.size = i;

	//printf("free in _unsymmetricPairWiseAlignment alignmentMatrix, ij");
	myfree(lenRead, alignmentMatrix, ij);
	if (mismnum>0){
		for (j=0; j<mismnum; j++){
			mfree_TrueMismatch(mmismatchlist[j]);
		}
		free(mmismatchlist);
	}
	/*
	for (i=0; i<lenRead+1; i++){
		free(alignmentMatrix[i]);
		alignmentMatrix[i] = NULL;
		free(ij[i]);
		ij[i] = NULL;
	}

	free(alignmentMatrix);
	free(ij);*/

	return ar;
}

char * unsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum){
	int j, scorelen;
	char * returnstr;
	char * scorestr;
	
	struct AlignRes ar = _unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint, mmismatch, mislen, mismnum);

	returnstr = NULL;
	scorelen = 20;

	returnstr = malloc(sizeof(char) * ((lenPerfRep+lenRead+2)*2+scorelen+1));
	m_init(returnstr, (lenPerfRep+lenRead+2)*2+scorelen+1);
	scorestr = malloc(sizeof(char) * scorelen);
	m_init(scorestr, scorelen);
	sprintf(scorestr, "%d",  ar.score);

	for(j=0; j<ar.size; j++){
                if (ar.correctedRead[j]=='\0'){
                        returnstr[j] = ';';
                }else{
                        returnstr[j] = ar.correctedRead[j];
                }
                if (ar.perfRep[j]=='\0'){
                        //returnstr[j+ar.size] = '\0';
								returnstr[j+ar.size] = ';';
                }else{
                        returnstr[j+ar.size] = ar.perfRep[j];
                }
        }

	for(j=0; j<scorelen; j++){
		returnstr[j+2*ar.size] = scorestr[j];
	}

   if (isprint) {
		printf("In unsymmetricPairWiseAlignment:\nCorrect=%s\nPerfect=%s\n", ar.correctedRead, ar.perfRep);
      printf("In unsymmetricPairWiseAlignment: all(len=%zu)=%s\n", strlen(returnstr), returnstr);
   }

	//printf("free in unsymmetricPairWiseAlignment ar.correctedRead");
   free(ar.correctedRead);
	//printf("free in unsymmetricPairWiseAlignment ar.perfRep");
   free(ar.perfRep);
	//printf("free in unsymmetricPairWiseAlignment scorestr");
	free(scorestr);

	return returnstr;
}

char* correctedByunsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum){
	char * newRead;
	int i, newi;

	struct AlignRes ar = _unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint, mmismatch, mislen, mismnum);

	newRead = malloc(sizeof(char) * (lenPerfRep+lenRead+1));
	m_init(newRead, lenPerfRep+lenRead+1);
	newi = 0;
        for (i=0; i<ar.size; i++){
                if(ar.correctedRead[i]!='-' && ar.perfRep[i]!='-' && ar.correctedRead[i]!='\0' && ar.perfRep[i]!='\0'){
                        newRead[newi] = ar.correctedRead[i]; newi++;
                }

                if ((ar.correctedRead[i]=='\0' && ar.perfRep[i]!='\0') || (ar.correctedRead[i]!='\0' && ar.perfRep[i]=='\0')){
                        printf("Errors: the size after alignment is not equal: \n%s, \n%s\n", ar.correctedRead, ar.perfRep);
                }
        }

        newRead[newi] = '\0'; newi++;

	if (isprint) {
                printf("In correctedByunsymmetricPairWiseAlignment:\nCorrect=%s\nPerfect=%s\nNew=%s\n", ar.correctedRead, ar.perfRep, newRead);
        }

	//printf("In correctedByunsymmetricPairWiseAlignment:\nCorrect=%s\nPerfect=%s\n    New=%s\n", ar.correctedRead, ar.perfRep, newRead);

	//printf("free in correctedByunsymmetricPairWiseAlignment ar.correctedRead %d", ar.size);
	//printf("test:%s\n", ar.perfRep);
	//printf("test:%s\n     %s\n     %s", perfectRepeat, read, ar.correctedRead);
	free(ar.correctedRead);
	//printf("free in correctedByunsymmetricPairWiseAlignment ar.perfRep");
   free(ar.perfRep);

	return newRead;

}

int main(int argc, const char * argv[]){
	int match, mismatch, gap_in_perf, gap_in_read, gap_before_after;
	char* perfectRepeat;
	char* read;
	int lenPerfRep, lenRead;
	char* newread;
	int isprint;
	int bandWidth;
	char* mmismatch;
	int mislen;
	int mismnum;

	mismnum = 2; mismnum = 1;
	mmismatch = "TG:TT:-1:0;CT:TT:0:1"; mmismatch = "TG:TT:-1:0";
	//           12345678901234567890
	mislen = 20; mislen = 10;

	match=10; mismatch=-9; gap_in_perf=-2; gap_in_read=-13; gap_before_after = -1;
	match=3; mismatch=-3; gap_in_perf=-2; gap_in_read=-10; gap_before_after = -1;
	bandWidth = 50;
	
	perfectRepeat = "CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
	read          = "CTGAGTGTTGAGCTGGCTGACTAGATTGCTTGGTATCGTGCTTTTGCGCTG";
	
	lenPerfRep = strlen(perfectRepeat); lenRead = strlen(read);
	isprint = 1;
	newread = NULL;
	newread = correctedByunsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint, mmismatch, mislen, mismnum);
	//newread = unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, isprint);

	printf("Main1: all(len=%zu)=%s\n", strlen(newread), newread);
	free(newread);

	newread = NULL;
	newread = unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint, mmismatch, mislen, mismnum);
	printf("Main2: all(len=%zu)=%s\n", strlen(newread), newread);

	free(newread);
	
	return 0;
}

