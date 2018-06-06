#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double dot_vectors(double *a, double *b, int n) {

    double tmp = 0.0;
    double tmpa, tmpb;
    int i;
    
    for (i=0; i<n; i++) {
        tmpa = *(a+i);
        tmpb = *(b+i);
        tmp += tmpa*tmpb;
    }
    
    return(tmp);
}
