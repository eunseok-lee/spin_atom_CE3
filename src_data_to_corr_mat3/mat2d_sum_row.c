#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double mat2d_sum_row(double *matA, int row_id, int n_col) {
    int i, j;
    double sum = 0.0;
    
    for (i=0;i<n_col;i++) {
        sum += *(matA+n_col*row_id+i);
    }
    
    return(sum);
}
