#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// convex hull is obtained as cvh(n_bin[0],...,n_bin[end])
// inputs are y = Ef(x). x:n_data by 4, y:n_data by 1
void obtain_convex_h_Nd(int n_data, double *ef2cvh, int *x, double *y, int Ndim, int *nCat_min, int *nCat_max, double *E0) {
    
    int i, j, k, l, nid, id1, id2, ctr;
    double tmp;
    int convex_h_ids[n_data];   //preallocated, not all filled
    double convex_h_x[n_data][Ndim];   //preallocated, not all filled
    double convex_h_y[n_data];  //preallocated, not all filled

    printf("Convex hull is created.\n");
    printf("Input data\n");
    printf("---------------------------\n");
    printf("  nL   nN   nM   nC   nV    Ef\n");
    for (i=0;i<n_data;i++) {
        for (j=0;j<Ndim;j++)
            printf("%4d ",x[Ndim*i+j]);
        printf("  % .4f\n",y[i]);
    }
    printf("---------------------------\n");

    // calculate the gradient and find Ndim base convex hull points
    double dEdx[n_data][Ndim];
    double dx[n_data][Ndim];
    for (i=0;i<n_data;i++) {
        printf("%d-data:%d,%d,%d,%d,%d, dE/dx=[",i,x[i*Ndim],x[i*Ndim+1],x[i*Ndim+2],x[i*Ndim+3],x[i*Ndim+4]);
        for (j=0;j<Ndim;j++) {
            dx[i][j] = (double) (nCat_max[j]-x[i*Ndim+j])*1.0/(nCat_max[j]-nCat_min[j]);
            if (dx[i][j] == 0.0)  // case of convex hull
                dEdx[i][j] = 0.0;
            else
                dEdx[i][j] = (E0[j]-y[i]) / dx[i][j];
            printf("% .2e ",dEdx[i][j]);
        }
        printf("]\n");
    }
    // check if the point is the convex hull
    int ctr_set[Ndim];
    nid = 0;
    for (i=0;i<n_data;i++) {
        if (dx[i][j] == 0.0) {
            convex_h_ids[nid] = i;
            for (j=0;j<Ndim;j++)
                convex_h_x[nid][j] = (double) x[i*Ndim+j]*1.0;
            convex_h_y[nid] = y[i];
            nid++;
            continue;
        }
        for (k=0;k<Ndim;k++) {
            ctr_set[k] = 1;
            for (j=0;j<n_data;j++) {
                if (j==i)
                    continue;
                if (dEdx[j][k] < dEdx[i][k]) {
                    ctr_set[k] = 0;
                    break;
                }
            }
        }
        ctr = 0;
        for (k=0;k<Ndim;k++)
            ctr += ctr_set[k];
        if (ctr > 0) {
            convex_h_ids[nid] = i;
            for (j=0;j<Ndim;j++)
                convex_h_x[nid][j] = (double) x[i*Ndim+j]*1.0;
            convex_h_y[nid] = y[i];
            nid++;
        }
    }
    printf("%d cvh points were identified.\n",nid);
    for (i=0;i<nid;i++) {
        printf("%d ",convex_h_ids[i]);
        for (j=0;j<Ndim;j++)
            printf("%f ",convex_h_x[i][j]);
        printf("%f\n",convex_h_y[i]);
    }

    // A*c = b
    gsl_matrix * A = gsl_matrix_alloc (Ndim,Ndim);
    gsl_matrix * V = gsl_matrix_alloc (Ndim,Ndim);
    gsl_vector * S = gsl_vector_alloc (Ndim);
    gsl_vector * work = gsl_vector_alloc (Ndim);
    gsl_vector * b = gsl_vector_alloc (Ndim);
    gsl_vector * c = gsl_vector_alloc (Ndim);

    double dij2[nid];
    int ids[nid];
    double ef_cvh;
    for (i=0;i<n_data;i++) {
        for (j=0;j<nid;j++) {
            dij2[j] = 0.0;
            for (k=0;k<Ndim;k++) {
                tmp = (double) x[i*Ndim+k]*1.0;
                dij2[j] += pow(tmp - convex_h_x[j][k],2);
            }
            ids[j] = j;
        }
        sort_array_ids(dij2,ids,nid);
        if (dij2[0]==0.0)
            ef2cvh[i] = 0.0;
        else {
            for (j=0;j<Ndim;j++) {
                for (k=0;k<Ndim;k++)
                    gsl_matrix_set (A,k,j,convex_h_x[ids[j]][k]);
                gsl_vector_set(b,j,x[i*Ndim+j]*1.0);
            }
            gsl_linalg_SV_decomp(A, V, S, work);
            gsl_linalg_SV_solve (A, V, S, b, c);
            ef_cvh = 0.0;
            for (j=0;j<Ndim;j++)
                ef_cvh += gsl_vector_get(c,j)*convex_h_y[ids[j]];
            ef2cvh[i] = y[i]-ef_cvh;
        }
    }
    
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(b);
    gsl_vector_free(c);

}













