#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define qh_QHimport
#include <libqhull_r/qhull_ra.h>

void print_summary(qhT *qh);
double findDelaunay(qhT *qh, int dim, coordT *TargetCoord);
double findNearestCentrum(qhT *qh, int dim, coordT *TargetCoord);
double findProjectedVertices(qhT *qh, int dim, coordT *TargetCoord);

void obtain_convex_h_qhull(int n_data, double *Ef2cvh, int *x, double *y, int NCatDim, int *nCat_min, int *nCat_max, double *E0, char *qhullflags) {
    
    int i, j, k, l, nid, Ndim;
    double tmp;
    
    Ndim = NCatDim - 1;     // because nV is dependent on other nCats
    printf("Convex hull is created using qhull.\n");
    printf("Input data\n");
    printf("---------------------------\n");
    printf("  nL   nN   nM   nC   nV    Ef\n");
    for (i=0;i<n_data;i++) {
        for (j=0;j<NCatDim;j++)
            printf("%4d ",x[NCatDim*i+j]);
        printf("  % .4f\n",y[i]);
    }
    printf("---------------------------\n");

    // Convex Hull by calling qhull
    int DIM;
    int SIZEcube, SIZEdiamond, TOTpoints;
    int dim;            /* dimension of points */
    int numpoints;      /* number of points */
    DIM = Ndim + 1;     /* because of ef as one additional dimension */
    dim = DIM;
    TOTpoints = n_data;
    numpoints= n_data;

//    coordT points[(DIM+1)*TOTpoints]; /* array of coordinates for each point */
    coordT points[DIM*TOTpoints]; /* array of coordinates for each point */
    coordT *rows[TOTpoints];
    boolT ismalloc= False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
    char flags[250];          /* option flags for qhull, see qh-quick.htm */
    FILE *outfile= stdout;    /* output from qh_produce_output()
                               use NULL to skip qh_produce_output() */
    FILE *errfile= stderr;    /* error messages from qhull code */
    int exitcode;             /* 0 if no error from qhull */
    facetT *facet;            /* set by FORALLfacets */
    int curlong, totlong;     /* memory remaining after qh_memfreeshort */

    qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
    qhT *qh= &qh_qh;
    
    QHULL_LIB_CHECK
    
    qh_zero(qh, errfile);

    printf( "\ncompute %d-d convexhull, options:%s\n", dim, qhullflags);
//    sprintf(flags, "qhull o Qx ");     // Error generated
    sprintf(flags, "qhull o %s ",qhullflags);
//    sprintf(flags, "qhull o QJ ");
//    sprintf(flags, "qhull o C-0 ");     // Error generated
//    sprintf(flags, "qhull o C-n ");     // Error generated
//    sprintf(flags, "qhull o Cn ");     // Error generated
//    sprintf(flags, "qhull o Qt ");     // Error generated
//    sprintf(flags, "qhull o Qx C-n ");     // Error generated
    for (i=0;i<n_data;i++) {    //size of points is n_data by (Ndim+1)
        for (j=0;j<Ndim;j++)    //size of x is n_data by Ndim
            points[dim*i+j] = x[NCatDim*i+j];
        points[dim*i+Ndim] = y[i];
    }
    for (i=0;i<numpoints;i++)   //size of rows is n_data x (Ndim+1)
        rows[i]= points+dim*i;
//    for (i=numpoints; i--; )
//        rows[i]= points+dim*i;
    qh_printmatrix(qh, outfile, "input", rows, numpoints, dim);
    exitcode= qh_new_qhull(qh, dim, numpoints, points, ismalloc,
                           flags, outfile, errfile);
    if (!exitcode) {                  /* if no error */
        print_summary(qh);
        printf( "\nThe relative Ef w.r.t. the convex hull:\n", dim);
        exitcode= setjmp(qh->errexit);
        if (!exitcode) {
            /* Trap Qhull errors in findDelaunay().  Without the setjmp(), Qhull
             will exit() after reporting an error */
            qh->NOerrexit= False;
            for (i=0;i<numpoints;i++) {
//                Ef2cvh[i] = findDelaunay(qh, DIM, rows[i]);
//                Ef2cvh[i] = findNearestVertices(qh, DIM, rows[i]);
//                Ef2cvh[i] = findNearestCentrum(qh, DIM, rows[i]);
                Ef2cvh[i] = findProjectedVertices(qh, DIM, rows[i]);
                printf("Ef2cvh(%d) = % .4e\n", i, Ef2cvh[i]);
            }
        }
        qh->NOerrexit= True;
    }
    qh_freeqhull(qh, !qh_ALL);                 /* free long memory */
    qh_memfreeshort(qh, &curlong, &totlong);  /* free short memory and memory allocator */
    if (curlong || totlong)
        fprintf(errfile, "qhull internal warning (user_eg, #2): did not free %d bytes of long memory (%d pieces)\n", totlong, curlong);
}

void print_summary(qhT *qh) {
    facetT *facet;
    vertexT *vertex;
    int k;
    
    printf("\n%d facets with normals:\n", qh->num_facets);
    FORALLfacets {
        for (k=0; k < qh->hull_dim; k++) {
            printf("%6.2g ", facet->normal[k]);
//            printf("", facet->vertices);
        }
        printf("\n");
    }
    printf("\n%d vertices:\n",qh->num_vertices);
    FORALLvertices {
        qh_printvertex(qh, qh->fout, vertex);
    }
}

double findDelaunay(qhT *qh, int dim, coordT *TargetCoord) {
    int j, k;
    coordT tmppoint[ 100];
    boolT isoutside;
    realT bestdist;
    facetT *facet;
    vertexT *vertex, **vertexp;
    
    for (k= 0; k < dim; k++)
        tmppoint[k]= TargetCoord[k];
    qh_setdelaunay(qh, dim+1, 1, tmppoint);
    facet= qh_findbestfacet(qh, tmppoint, qh_ALL, &bestdist, &isoutside);
    if (facet->tricoplanar) {
        fprintf(stderr, "findDelaunay: not implemented for triangulated, non-simplicial Delaunay regions (tricoplanar facet, f%d).\n",
                facet->id);
        qh_errexit(qh, qh_ERRqhull, facet, NULL);
    }
    
    for (k= 0; k < dim-1; k++)
        printf("% .3f ",tmppoint[k]);
    printf("% .3f, ",tmppoint[dim-1]);
    printf("facet id: %d\n",facet->id);
/*    FOREACHvertex_(facet->vertices) {
        for (k=0; k < dim; k++)
            printf("%5.3f ", vertex->point[k]);
        printf("vertex id: %d, point id: %d\n", vertex->id, qh_pointid(qh, vertex->point));
    }*/
    FOREACHvertex_(facet->vertices) {
        printf("vertex id: %d, point id: %d\n", vertex->id, qh_pointid(qh, vertex->point));
    }
    
    gsl_matrix * A = gsl_matrix_alloc (dim,dim);
    gsl_matrix * V = gsl_matrix_alloc (dim,dim);
    gsl_vector * S = gsl_vector_alloc (dim);
    gsl_vector * work = gsl_vector_alloc (dim);
    gsl_vector * b = gsl_vector_alloc (dim);
    gsl_vector * x = gsl_vector_alloc (dim);
    
    j= 0;
    FOREACHvertex_(facet->vertices) {
        for (k= 0; k < dim; k++)
            gsl_matrix_set (A,k,j, vertex->point[k]);
        gsl_vector_set(b,j,tmppoint[j]);
        j++;
    }
    double ef[dim];
    // dismiss the last row corresponding to Ef
    // & apply the constriction on the coefficients sum_x_i = 1
    for (k= 0; k < dim; k++) {
        ef[k] = gsl_matrix_get (A,dim-1,k);
        gsl_matrix_set (A,dim-1,k,1.0);
    }
    gsl_vector_set(b,dim-1,2.0);

    gsl_linalg_SV_decomp(A, V, S, work);
    gsl_linalg_SV_solve (A, V, S, b, x);
    
    double ef_cvh = 0.0;
    for (k= 0; k < dim; k++)
        ef_cvh += gsl_vector_get(x,k)*ef[k];

    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(b);
    gsl_vector_free(x);

    return(TargetCoord[dim-1]-ef_cvh);
}

double findNearestCentrum(qhT *qh, int dim, coordT *TargetCoord) {
    int i, j, k;
    int nf, nv;
    double ef_cvh;
    double tmp;
    double tmppoint[dim];
    boolT isoutside;
    realT bestdist;
    facetT *facet;
    vertexT *vertex, **vertexp;
    pointT *centrum;
    
    nf = qh->num_facets;
    nv = qh->num_vertices;
    
    for (k= 0; k < dim; k++)
        tmppoint[k]= TargetCoord[k];
    
    double facet_normal[nf][dim];
    double facet_center[nf][dim];
    int facet_id[nf];
    int num_lowerhullfacets;

    j = 0;
    FORALLfacets {
        if (facet->normal[dim-1] < 0.0) {//is the slope in ef direction negative?
            for (k= 0; k < dim; k++) {
                facet_normal[j][k] = facet->normal[k];
                if (qh->CENTERtype == qh_AScentrum)
                    centrum= facet->center;
                else
                    centrum= qh_getcentrum(qh, facet);
                facet_center[j][k] = centrum[k];
                facet_id[j] = facet->id;
            }
            j++;
        }
    }
    num_lowerhullfacets = j;

    double dij2[num_lowerhullfacets];
    int ids[num_lowerhullfacets];
    for (j= 0; j< num_lowerhullfacets; j++) {
        tmp = 0.0;
        for (k= 0; k < dim-1; k++) {    //ef is excluded for configurational coordinate
            tmp += pow(facet_center[j][k] - tmppoint[k],2);
        }
        dij2[j] = tmp;
        ids[j] = j;
    }
    sort_array_ids(dij2,ids,num_lowerhullfacets);

    int projectedfacet = facet_id[ids[0]];
    double vertices_point[dim][dim];
    FORALLfacets {
        if (facet->id == projectedfacet) {
            j = 0;
            FOREACHvertex_(facet->vertices) {
                for (k=0; k < dim; k++)
                    vertices_point[j][k] = vertex->point[k];
                j++;
            }
        }
    }

    gsl_matrix * A = gsl_matrix_alloc (dim,dim);
    gsl_matrix * V = gsl_matrix_alloc (dim,dim);
    gsl_vector * S = gsl_vector_alloc (dim);
    gsl_vector * work = gsl_vector_alloc (dim);
    gsl_vector * b = gsl_vector_alloc (dim);
    gsl_vector * x = gsl_vector_alloc (dim);
    
    for (j= 0; j< dim; j++) {
        for (k= 0; k < dim; k++)
            gsl_matrix_set (A,k,j,vertices_point[j][k]);
        gsl_vector_set(b,j,tmppoint[j]);
    }
    for (k= 0; k < dim; k++) {      //dismiss the last row corresponding to ef
        gsl_matrix_set (A,dim-1,k,1.0);
    }
    gsl_vector_set(b,dim-1,1.0);    //dismiss the last row corresponding to ef
    
    printf("A:\n");
    for (j= 0; j< dim; j++) {
        for (k= 0; k < dim; k++)
            printf("%g ", gsl_matrix_get (A,j,k));
        printf("\n");
    }
    printf("b:\n");
    for (j= 0; j< dim; j++)
        printf("%g\n", gsl_vector_get(b,j));
    
    gsl_linalg_SV_decomp(A, V, S, work);
    gsl_linalg_SV_solve (A, V, S, b, x);
    printf("solution:\n");
    for (j= 0; j< dim; j++)
        printf("%g\n", gsl_vector_get(x,j));
    // need to check if the components of x is in [0,2]
    for (j= 0; j< dim; j++) {
        tmp = gsl_vector_get(x,j);
        if (tmp < -0.1 || tmp > 1.0) {
            printf("\n%d-component of x is %f, out of range\n", j, tmp);
            exit(-1);
        }
    }
    tmp = 0.0;
    for (j= 0; j< dim; j++)
        tmp += gsl_vector_get(x,j)*vertices_point[j][dim-1];
    ef_cvh = tmppoint[dim-1]-tmp;
    
    gsl_matrix_free(A);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(b);
    gsl_vector_free(x);
    
    return(ef_cvh);
}

double findProjectedVertices(qhT *qh, int dim, coordT *TargetCoord) {
    int i, j, k;
    int nf, nv;
    int ctr;
    double ef_cvh;
    double tmp;
    double tmppoint[dim];
    double Ef_store[dim];
    boolT isoutside;
    realT bestdist;
    facetT *facet;
    vertexT *vertex, **vertexp;
    pointT *centrum;

    gsl_matrix * A = gsl_matrix_alloc (dim,dim);
    gsl_matrix * V = gsl_matrix_alloc (dim,dim);
    gsl_vector * S = gsl_vector_alloc (dim);
    gsl_vector * work = gsl_vector_alloc (dim);
    gsl_vector * b = gsl_vector_alloc (dim);
    gsl_vector * x = gsl_vector_alloc (dim);
    gsl_matrix * A0 = gsl_matrix_alloc (dim,dim);

    nf = qh->num_facets;
    nv = qh->num_vertices;
    
    for (k= 0; k < dim; k++)
        tmppoint[k]= TargetCoord[k];
    
    FORALLfacets {
        if (facet->normal[dim-1] < 0.0) {
            j = 0;
            FOREACHvertex_(facet->vertices) {
                for (k=0; k < dim; k++) {
                    gsl_matrix_set (A,k,j,vertex->point[k]);    //(j,k) to (k,j)
                    gsl_matrix_set (A0,k,j,vertex->point[k]);   //store original A
                }
                gsl_vector_set(b,j,tmppoint[j]);
                j++;
            }
            if (j != dim) {
                printf("[E] j is not equal to dim(=%d)\n",dim);
                exit(-1);
            }
            for (k= 0; k < dim; k++) {      //dismiss the last row corresponding to ef
                Ef_store[k] = gsl_matrix_get (A,dim-1,k);
                gsl_matrix_set (A,dim-1,k,1.0);
            }
            gsl_vector_set(b,dim-1,1.0);    //dismiss the last row corresponding to ef
            
            gsl_linalg_SV_decomp(A, V, S, work);
            gsl_linalg_SV_solve (A, V, S, b, x);
            ctr = 1;
            for (k= 0; k < dim; k++) {
                tmp = gsl_vector_get(x,k);
                if (tmp < -1.0e-8 || tmp > 1.0)
                    ctr = ctr*0;
            }
            if (ctr == 1) {     //matching case!
                if (0) {
                    printf("Coordate:\n");
                    for (j= 0; j< dim; j++)
                        printf("%g ",tmppoint[j]);
                    printf("\n");
                    printf("Matching Facet id:%d\n",facet->id);
                    printf("A & b:\n");
                    for (j= 0; j< dim-1; j++) {
                        for (k= 0; k < dim; k++)
                            printf("%g ", gsl_matrix_get (A0,j,k));
                        printf("%.4f\n", gsl_vector_get(b,j));
                    }
                    printf("x:\n");
                    for (j= 0; j< dim; j++)
                        printf("%.4f ",gsl_vector_get(x,j));
                    printf("\n");
                }
                tmp = 0.0;
                for (j= 0; j< dim; j++)
                    tmp += gsl_vector_get(x,j)*Ef_store[j];
                ef_cvh = tmppoint[dim-1]-tmp;
                return(ef_cvh);
            }
        }
    }
    printf("[E] No matching facet\n");
    exit(-1);
}










