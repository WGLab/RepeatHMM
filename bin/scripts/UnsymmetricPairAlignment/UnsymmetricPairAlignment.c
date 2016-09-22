
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


struct Coord{
	int i;
	int j;
};

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

        free(alignmentMatrix);
        free(ij);
}

//match=10; mismatch=-10; gap_in_perf=-1; gap_in_read=-10;
struct AlignRes _unsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint){
	int i,j,k,diagonal,up,left; //, hasacgt;
	int ** alignmentMatrix;
	struct Coord ** ij;
	struct Coord maxscoreij;
	int maxscore;
	struct AlignRes ar;

	ar.size = 0;
	maxscore = INT_MIN;
	maxscoreij.i = -1;
	maxscoreij.j = -1;
	
	ar.correctedRead = malloc(sizeof(char) * (lenPerfRep+lenRead+1));
	ar.perfRep = malloc(sizeof(char) * (lenPerfRep+lenRead+1));
	
	alignmentMatrix = malloc(sizeof(int*) * (lenRead+1));	
	ij = malloc(sizeof(struct Coord*) * (lenRead+1));
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
				else{ diagonal += mismatch;}
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

	k = lenPerfRep+lenRead-1;
	ar.correctedRead[k] = '\0'; ar.perfRep[k]='\0'; k--;
 	i = lenRead; j = lenPerfRep;
	ar.score = alignmentMatrix[i][j];
	
	i = maxscoreij.i; j = maxscoreij.j;
	ar.score = maxscore;

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
			ar.correctedRead[k] = read[i-1]; ar.perfRep[k]=perfectRepeat[j-1]; k--;
			i = i-1; j = j-1;
		}
		else if (ij[i][j].i==i-1 && ij[i][j].j == j){
			ar.correctedRead[k] = read[i-1]; ar.perfRep[k]='-'; k--;
			i = i-1; //j = j;
		}
		else if (ij[i][j].i==i && ij[i][j].j == j-1){
			ar.correctedRead[k] = '-'; ar.perfRep[k]=perfectRepeat[j-1]; k--;
			j = j-1; //i = i;
		}
	}

	i = 0; k++; 
	for (j=k; j<lenPerfRep+lenRead; j++){
		ar.correctedRead[i] = ar.correctedRead[j];
		ar.perfRep[i] = ar.perfRep[j];

		if ((ar.correctedRead[i]=='\0' && ar.perfRep[i]!='\0') || (ar.correctedRead[i]!='\0' && ar.perfRep[i]=='\0')){
			printf("Errors: the size after alignment is not equal: \n%s, \n%s\n", ar.correctedRead, ar.perfRep);
		}

		i++;
	}

	ar.size = i;

	myfree(lenRead, alignmentMatrix, ij);
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

char * unsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint){
	int j;
	char * returnstr;
	char * scorestr;
	returnstr = NULL;
	
	struct AlignRes ar = _unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint);

	returnstr = malloc(sizeof(char) * ((lenPerfRep+lenRead+2)*2+11));
	scorestr = malloc(sizeof(char) * 11);
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

	for(j=0; j<11; j++){
		returnstr[j+2*ar.size] = scorestr[j];
	}

        if (isprint) {
		printf("In unsymmetricPairWiseAlignment: Correct=%s\nPerfect=%s\n", ar.correctedRead, ar.perfRep);
                printf("In unsymmetricPairWiseAlignment: all(len=%zu)=%s\n", strlen(returnstr), returnstr);
        }

        free(ar.correctedRead);
        free(ar.perfRep);
	free(scorestr);

	return returnstr;
}

char* correctedByunsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint){
	char * newRead;
	int i, newi;

	struct AlignRes ar = _unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint);

	newRead = malloc(sizeof(char) * (lenPerfRep+lenRead+1));
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
                printf("In correctedByunsymmetricPairWiseAlignment: Correct=%s\nPerfect=%s\nNew=%s\n", ar.correctedRead, ar.perfRep, newRead);
        }

	free(ar.correctedRead);
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

	match=10; mismatch=-9; gap_in_perf=-2; gap_in_read=-13; gap_before_after = -1;
	bandWidth = 50;
	
	perfectRepeat = "CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG";
	read          = "CTGAGTGTTGAGCTGGCTGACTAGATTGCTTGGTATCGTGCTTTTGCGCTG";
	
	lenPerfRep = strlen(perfectRepeat); lenRead = strlen(read);
	isprint = 1;
	newread = NULL;
	newread = correctedByunsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint);
	//newread = unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, isprint);

	printf("Main1: all(len=%zu)=%s\n", strlen(newread), newread);
	free(newread);

	newread = NULL;
	newread = unsymmetricPairWiseAlignment(perfectRepeat, lenPerfRep, read, lenRead, match, mismatch, gap_in_perf, gap_in_read, gap_before_after, bandWidth, isprint);
	printf("Main2: all(len=%zu)=%s\n", strlen(newread), newread);

	free(newread);
	
	return 0;
}

