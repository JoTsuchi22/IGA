#include <stdio.h>
// #include <string.h>
// #include <math.h>
#include <stdlib.h>

int main(){
    double **A;
    A = malloc(sizeof(double*) * 16);
    // A = (double **)malloc(sizeof(double*) * 16);
    // A[0] = 1.0;
    A[3][3] = 6.0;
    // printf("%le", A[0]);
    printf("%le", A[3][3]);
    free(A);
}