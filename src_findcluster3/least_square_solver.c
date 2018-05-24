#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cem.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void least_square_solver(double *matA, int n_row, int n_col, double *y, double *sol_x) {
    
    int i, j, k, tmp;
    double sum;
    gsl_matrix * A = gsl_matrix_alloc (n_row,n_col);
    gsl_matrix * V = gsl_matrix_alloc (n_col,n_col);
    gsl_vector * S = gsl_vector_alloc (n_col);
    gsl_vector * work = gsl_vector_alloc (n_col);
//    gsl_matrix * U = gsl_matrix_alloc (n_row-1,n_col_ids);
    gsl_vector * b = gsl_vector_alloc (n_row);
    gsl_vector * x = gsl_vector_alloc (n_col);
    
/*    //printf A and b
    FILE *fpA, *fpU, *fpb;
    fpA = fopen("A.dat","w");
    fpU = fopen("U.dat","w");
    fpb = fopen("b.dat","w");*/
    
    // need to copy matA to A, y to b!
    for (i=0;i<n_row;i++) {
        for (j=0; j<n_col; j++) {
            gsl_matrix_set (A, i, j, *(matA+i*n_col+j));
        }
        gsl_vector_set (b, i, y[i]);
    }
/*    for (i=0;i<n_row;i++) {
        for (j=0; j<n_col; j++) {
            fprintf(fpA,"%f ", gsl_matrix_get (A, i, j));
        }
        fprintf(fpA,"\n");
        fprintf(fpb,"%f\n", gsl_vector_get(b,i));
    }*/
    gsl_linalg_SV_decomp(A, V, S, work);    // on output, A is replaced by U!
/*    for (i=0;i<n_row;i++) {
        for (j=0; j<n_col; j++) {
            fprintf(fpU,"%f ", *(matA+i*n_col+j));
        }
        fprintf(fpU,"\n");
    }*/
    gsl_linalg_SV_solve (A, V, S, b, x);    // So, A should be used instead of U, here.
        
    for (i=0;i<n_col;i++) {
        sol_x[i] = gsl_vector_get(x,i);
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(b);
    gsl_vector_free(x);

}
