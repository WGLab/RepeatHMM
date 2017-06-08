
%newobject unsymmetricPairWiseAlignment;
%newobject correctedByunsymmetricPairWiseAlignment;
%module UnsymmetricPairAlignment
%{
/* put header files here or function declarations like below*/

extern char* unsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum);
extern char * correctedByunsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum);
%}

extern char* unsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum);
extern char * correctedByunsymmetricPairWiseAlignment(char* perfectRepeat, int lenPerfRep, char* read, int lenRead, int match, int mismatch, int gap_in_perf, int gap_in_read, int gap_before_after, int bandWidth, int isprint, char* mmismatch, int mislen, int mismnum);

