
gcc -o upa UnsymmetricPairAlignment.c

valgrind --tool=memcheck upa

