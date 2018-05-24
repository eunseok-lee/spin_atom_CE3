#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cem.h>
#include <mpi.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void randperm(int*, int);
void load_double_mat(double*, int, int, int, char*);
void load_int_mat(int*, int, int, int, char*);
void load_double_array(double*, int, char*);
void load_int_array(int*, int, char*);
void load_double_value(int,char);
void load_int_value(int,char);
int mat_size_row(int*);
int find_value_in_array(int*,int,int,int);
void mat2d_prod(double*, int, int, double*, int, int, double*, int, int);
double obtainerr2_par(int,int,double*, int, int, int*, int, double*);
void least_square_solver(double*, int, int, double*, double*);

#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

int main(int argc, char **argv)
{
    int n1, n2, n3;
    int np, nfu;
    double kv[3][3], pk[3][3], pp[8][3];
    double *rp, *rpL, *rpN, *rpO;
    int *map_to_cluster1, *map_to_cluster2, *map_to_cluster3;
    int *nlist;
    int neighbor_num;
    int ncluster, ncluster2, ncluster3, ncorr_col;
    int c2start, c3start;
    int data;
    double *E, *Ef;
    int ndata;
    int i, j, k, nj, nk, ctn;
    int iter, max_iter;
    double kT=0.0517;
    int howmanycluster;
    int dispfreq;
    
    int construct_corr_mat = 0;

    // mpi parameters and initialization
    int numprocs, rank, mtype;
    int row_dist_size;
    int row_ini, row_end, row_offset;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // Get # processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      // Get my rank (id)

    if (numprocs>1)
        printf("Parallel computing: the input files are read by all nodes. Multiple messages are not due to error!\n",numprocs);

    // Construct conrrelation matrix. It can also be loaded.
    
    if (construct_corr_mat) {
        printf("Not activated yet.\n");
/*        // Task: make a function to load values from clusterlist run, later
        np = 48;            // for 3x2x2 cell
        nfu = np/2;         // for 3x2x2 cell
        neighbor_num = 47;  // for 3x2x2 cell
        ncluster2 = 10;
        ncluster3 = 62;
        
        // load lattice and cluster data
        load_double_mat(kv,3,3,1,"kv.dat");
        load_double_mat(rp,np,3,1,"rp.dat");
        load_double_mat(rpL,npL,3,1,"rpL.dat");
        load_double_mat(rpN,npN,3,1,"rpN.dat");
        load_double_mat(rpO,npO,3,1,"rpO.dat");
        load_double_array(cluster2,"cluster2.dat");
        load_double_array(cluster3,"cluster3.dat");
        load_int_array(map_to_cluster1,"map_to_cluster1.dat");
        load_int_array(map_to_cluster2,"map_to_cluster2.dat");
        load_int_array(map_to_cluster3,"map_to_cluster3.dat");
        load_int_array(nlist,"nlist.dat");
        
        ncluster  = ncluster2+ncluster3+1+1;
        ncorr_col = 1 + 2*np + 3*ncluster2 + 4*ncluster3;
        c2start = 3;    //2*length(cluster1)+1. Don't need -1: cluster2 idx starts from 0
        c3start = 3*ncluster2 + c2start;

        // load LNO occupation variable data
        load_int_array(data,"LNC_occu_data.dat");
        load_double_array(E,"E.dat");
        load_double_array(Ef,"Ef.dat");
        ndata = mat_size_row(data);
        int nL[ndata];
        int nC[ndata];
        int nN[ndata];
        
        for (i=0;i<ndata;i++) {
            nL[i] = find_value_in_array(data,i,0);
            nC[i] = find_value_in_array(data,i,1);
            nN[i] = find_value_in_array(data,i,-1);
        }
        
        // initialize the correlation matrix
        double corr_mat[ndata][ncorr_col];
        obtain_corr_mat(corr_mat, ndata, ncorr_col, data, c2start, c3start);

        // post data processing: excluding degenerate cases
        int nondegenerate[ndata];
        for (i=0;i<ndata;i++)
            nondegenerate[i] = 1;
        for (i=0;i<ndata;i++) {
            for (j=0;j<i;j++) {
                ctn = 1;
                for (l=0;l<ncorr_col;l++)
                    if (corr_mat[j][l]~=corr_mat[i][l]) {
                        ctn = 0;
                        break;
                    }
                if (ctn == 1)
                    if (Ef[j] < Ef[i])
                        nondegenerate[i] = 0;
                    else
                        nondegenerate[j] = 0;
            }
            for (j=i+1;j<ndata;j++) {
                ctn = 1;
                for (l=0;l<ncorr_col;l++)
                    if (corr_mat[j][l]~=corr_mat[i][l]) {
                        ctn = 0;
                        break;
                    }
                if (ctn == 1)
                    if (Ef[j] < Ef[i])
                        nondegenerate[i] = 0;
                    else
                        nondegenerate[j] = 0;
            }
        }
        
        int n_nondegenerate = 0;
        for (i=0;i<ndata;i++)
            if (nondegenerate[i] == 1)
                n_nondegenerate++;
        double corr_mat_u[n_nondegenerate][ncorr_col];
        double Ef_u[n_nondegenerate];
        double nC_r[n_nondegenerate];
        l = 0;
        for (i=0;i<ndata;i++)
            if (nondegenerate[i] == 1) {
                for (j=0;j<ncorr_col;j++)
                    corr_mat_u[l][j] = corr_mat[i][j];
                Ef_u[l] = Ef[i];
                nC_u[l] = nC[i];
                l++;
            }
        
        // select data sets near convex hull
        int n_corr_mat_ug_row = 0;
        int near_convh[n_nondegenerate];
        for (i=0;i<n_nondegenerate;i++)
            if (Ef_u[i] < ef_conv(nC[i]) + near_convh_cutoff) {
                near_convh[i] = 1;
                n_corr_mat_ug_row++;
            }
            else
                near_convh[i] = 0;
        
        double corr_mat_ug[n_corr_mat_ug_row][ncorr_col];
        double Ef_ug[n_corr_mat_ug_row];
        double nC_ug[n_corr_mat_ug_row];
        l = 0;
        for (i=0;i<n_corr_mat_ug_row)
            if (near_convh[i] == 1) {
                for (j=0;j<ncorr_col;j++)
                    corr_mat_ug[l][j] = corr_mat_u[i][j];
                Ef_ug[l] = Ef_u[i];
                nC_ug[l] = nC_u[i];
                l++;
            }
        
        // select the non-singular columns only
        int n_non_singular_col = 0;
        int col_singularity[ncorr_col];
        for (i=0;i<ncorr_col;i++) {
            ctn = 0;
            for (j=1;j<n_corr_mat_ug_row;j++)
                if (corr_mat_ug[j][i] ~= corr_mat_ug[0][i]) {
                    ctn = 1;
                    break;
                }
            if (ctn == 1) {
                col_singularity[i] = 1;  // 0:singular, 1:non-singular
                n_non_singular_col++;
            }
            else
                col_singularity[i] = 0;
        }
        double corr_mat_ugs[n_corr_mat_ug_row][n_non_singular_col];
        int non_singular_col_id[n_non_singular_col];
        l = 0;
        for (i=0;i<n_non_singular_col)
            if (col_singularity[i] == 1) {
                for (j=0;j<n_corr_mat_ug_row;j++)
                    corr_mat_ugs[j][l] = corr_mat_ug[j][l];
                non_singular_col_id[l] = i;
                l++;
            }
 */
    }
//    else {
    int n_corr_mat_ug_row;
    int n_non_singular_col;
    FILE *fp;
    char paramfilename[100]="param_findcluster.dat";
    char corr_mat_ugs_filename[200];
    fp = fopen(paramfilename, "r");
    fscanf(fp, "%s %d %d %d %d %lf %d",corr_mat_ugs_filename, &n_corr_mat_ug_row, &n_non_singular_col, &howmanycluster, &max_iter, &kT, &dispfreq);
    fclose(fp);
    printf("n_corr_mat_ug_row=%d, n_non_singular_col=%d, howmanycluster=%d, max_iter=%d, dispfreq=%d\n",n_corr_mat_ug_row, n_non_singular_col, howmanycluster, max_iter, dispfreq);
    double *corr_mat_ugs;
    corr_mat_ugs = (double*) malloc(n_corr_mat_ug_row*n_non_singular_col*sizeof(double));
    load_double_mat(corr_mat_ugs,n_corr_mat_ug_row,n_non_singular_col,1,corr_mat_ugs_filename);
    
    int usefulcorr_col[n_non_singular_col];
    load_int_array(usefulcorr_col,n_non_singular_col,"usefulcorr_col.dat");
    
    double Ef_ug[n_corr_mat_ug_row];
    load_double_array(Ef_ug,n_corr_mat_ug_row,"Ef_ug.dat");
    
    int nC_ug[n_corr_mat_ug_row];
    load_int_array(nC_ug,n_corr_mat_ug_row,"nC_ug.dat");
    
//    }
    printf("step1\n");
    
    row_offset = n_corr_mat_ug_row/numprocs;
    double err2_sum_from_master;
    double err2_sum_from_worker;
    int cluster_set1[howmanycluster];
    
    /**** Start of MC simulation ****/
    if (rank == MASTER) {
        int cluster_set1_old[howmanycluster], cluster_set1_min[howmanycluster], cluster_set2[n_non_singular_col-howmanycluster];
        double cvs, cvs_old, cvs_min, cvs_tol;
        int select1, select2, target, candidate;
        
        cvs_tol = 0.01;
        for (i=0;i<howmanycluster;i++)
            cluster_set1[i] = i;
        for (i=0;i<n_non_singular_col-howmanycluster;i++)
            cluster_set2[i] = i+howmanycluster;
        
        mtype = FROM_MASTER;
        for (i=1;i<numprocs;i++) {
            MPI_Send(&cluster_set1, howmanycluster, MPI_INT, i, mtype, MPI_COMM_WORLD);
        }
        row_ini = rank * row_offset;
        row_end = row_ini +row_offset;
        err2_sum_from_master = obtainerr2_par(row_ini,row_end,corr_mat_ugs,n_corr_mat_ug_row,n_non_singular_col,cluster_set1,howmanycluster,Ef_ug);
        printf("initial err2 = %f\n",err2_sum_from_master);
        for (i=1;i<numprocs;i++) {
            mtype = FROM_WORKER;
            MPI_Recv(&err2_sum_from_worker, 1, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
            err2_sum_from_master += err2_sum_from_worker;
        }
        printf("initial err2 = %f\n",err2_sum_from_master);
        cvs = sqrt(err2_sum_from_master/n_corr_mat_ug_row);
        cvs_min = cvs;
        for (i=0;i<howmanycluster;i++)
            cluster_set1_min[i] = cluster_set1[i];
        printf("initial cvs = %f\n",cvs);
        for (iter=0;iter<max_iter;iter++) {
            select1 = (int) floor(drand48()*howmanycluster);
            select2 = (int) floor(drand48()*(n_non_singular_col-howmanycluster));
            target = cluster_set1[select1];
            candidate = cluster_set2[select2];
            cvs_old = cvs;
            cluster_set1[select1] = candidate;
            mtype = FROM_MASTER;
            for (i=1;i<numprocs;i++) {
                MPI_Send(&cluster_set1, howmanycluster, MPI_INT, i, mtype, MPI_COMM_WORLD);
            }
            row_ini = rank * row_offset;
            row_end = row_ini +row_offset;
            err2_sum_from_master = obtainerr2_par(row_ini,row_end,corr_mat_ugs,n_corr_mat_ug_row,n_non_singular_col,cluster_set1,howmanycluster,Ef_ug);
            for (i=1;i<numprocs;i++) {
                mtype = FROM_WORKER;
                MPI_Recv(&err2_sum_from_worker, 1, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
                err2_sum_from_master += err2_sum_from_worker;
            }
            cvs = sqrt(err2_sum_from_master/n_corr_mat_ug_row);
            if (cvs < cvs_old)
                cluster_set2[select2] = target;
            else if (exp(-(cvs-cvs_old)/kT) > drand48())
                cluster_set2[select2] = target;
            else {
                cluster_set1[select1] = target;
                cvs = cvs_old;
            }
            if (cvs < cvs_min) {
                for (i=0;i<howmanycluster;i++) 
                    cluster_set1_min[i] = cluster_set1[i];
                cvs_min = cvs;
            }
            if (cvs < cvs_tol)
                break;
            if (iter%dispfreq==0)
                printf("iter=%8d, cvs=%f, cvs_min=%f, cvs_min/fu=%f\n",iter,cvs,cvs_min,cvs_min/72);
        }
        
        // print the result
        printf("End of MC iteration. The selected clusters are...\n");
        printf("ID in non-singular | ID in entire\n");
        for (i=0;i<howmanycluster;i++)
            printf("    %4d                %4d\n",cluster_set1_min[i],usefulcorr_col[cluster_set1_min[i]]);
        
        double err_pred, rms;
        double Ef_pred[n_corr_mat_ug_row];
        double *corr_mat_ugs_reduced;
        corr_mat_ugs_reduced = (double*) malloc(n_corr_mat_ug_row*howmanycluster*sizeof(double));
        double eci[howmanycluster];
        for (i=0;i<n_corr_mat_ug_row;i++) {
            for (j=0;j<howmanycluster;j++)
                *(corr_mat_ugs_reduced+howmanycluster*i+j) = *(corr_mat_ugs+n_non_singular_col*i+cluster_set1_min[j]);
        }
        least_square_solver(corr_mat_ugs_reduced,n_corr_mat_ug_row,howmanycluster,Ef_ug,eci);
        mat2d_prod(Ef_pred,n_corr_mat_ug_row,1,corr_mat_ugs_reduced,n_corr_mat_ug_row,howmanycluster,eci,howmanycluster,1);
        err_pred = 0.0;
        for (i=0;i<n_corr_mat_ug_row;i++)
            err_pred += (Ef_ug[i] - Ef_pred[i])*(Ef_ug[i] - Ef_pred[i]);
        rms = sqrt(err_pred/n_corr_mat_ug_row);

        // Write result into files
        FILE *fp2;
        char cluster_set_filename[200];
        sprintf(cluster_set_filename, "cluster_set_magN%d.dat",howmanycluster);
        fp2 = fopen(cluster_set_filename,"w");
        for (i=0;i<howmanycluster;i++)
            fprintf(fp2,"%4d\t%4d\t%.4f\n",cluster_set1_min[i],usefulcorr_col[cluster_set1_min[i]],eci[i]);
    }
    else {
        for (iter=0;iter<(max_iter+1);iter++) {
            mtype = FROM_MASTER;
            MPI_Recv(&cluster_set1, howmanycluster, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            row_ini = rank * row_offset;
            row_end = row_ini + row_offset;
            err2_sum_from_worker = obtainerr2_par(row_ini,row_end,corr_mat_ugs,n_corr_mat_ug_row,n_non_singular_col,cluster_set1,howmanycluster,Ef_ug);
            mtype = FROM_WORKER;
            MPI_Send(&err2_sum_from_worker, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
        }
    }
    MPI_Finalize();
    
}

void load_double_mat(double *A, int Arow, int Acol, int Apgs, char *datfilename) {
    
    FILE *fp;
    int i, j, k;
    double tmp;
    
    printf("filename: %s\n", datfilename);
    fp = fopen(datfilename, "r");
    for (k=0;k<Apgs;k++)
        for (i=0;i<Arow;i++)
            for (j=0;j<Acol;j++) {
                fscanf(fp, "%lf", &tmp);
                *(A+Arow*Acol*k+Acol*i+j) = tmp;
            }
    
    fclose(fp);
    
}

void load_int_mat(int *A, int Arow, int Acol, int Apgs, char *datfilename) {
    
    FILE *fp;
    int i, j, k, tmp;
    
    printf("filename: %s\n", datfilename);
    fp = fopen(datfilename, "r");
    for (k=0;k<Apgs;k++)
        for (i=0;i<Arow;i++)
            for (j=0;j<Acol;j++) {
                fscanf(fp, "%d", &tmp);
                *(A+Arow*Acol*k+Acol*i+j) = tmp;
            }
    
    fclose(fp);
    
}

void load_int_array(int *array, int arraysize, char *datfilename) {
    
    FILE *fp;
    int i;
    
    printf("filename: %s\n", datfilename);
    fp = fopen(datfilename, "r");
    //    printf("step2\n");
    for (i=0;i<arraysize;i++) {
        fscanf(fp, "%d", &array[i]);
        //        printf("%d\n",array[i]);
    }
    
    fclose(fp);
    //    printf("step3\n");
}

void load_double_array(double *array, int arraysize, char *datfilename) {
    
    FILE *fp;
    int i;
    
    printf("filename: %s\n", datfilename);
    fp = fopen(datfilename, "r");
    //    printf("step4\n");
    for (i=0;i<arraysize;i++) {
        fscanf(fp, "%lf", &array[i]);
        //        printf("%f\n",array[i]);
    }
    
    fclose(fp);
    //    printf("step5\n");
}

void randperm(int *a, int np) {
    
    int i, id1, id2, id3, tmpr, num_perm;
    
    num_perm = 100*np;
    
    for (i=0;i<np;i++)
        *(a+i) = i;
    
    for (i=0;i<num_perm;i++) {
        id1 = rand() % np;
        id3 = rand() % np;
        id2 = (id1 + id3) % np;
        if (id1 != id2) {
            tmpr = *(a+id2);
            *(a+id2) = *(a+id1);
            *(a+id1) = tmpr;
        }
    }
}

int find_value_in_array(int *A, int target_rowA, int colA, int target_value) {
    int i;
    
    for (i=0;i<colA;i++)
        if (*(A+colA*target_rowA+i)==target_value)
            return(i);
    if (i==colA)
        return(-1);
}


























