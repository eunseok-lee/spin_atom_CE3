#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <mpi.h>
//#include <cem.h>

void randperm(int*, int);
void load_double_mat(double*, int, int, int, char*);
void load_int_mat(int*, int, int, int, char*);
void load_double_array(double*, int, char*);
void load_int_array(int*, int, char*);
void load_double_value(int,char);
void load_int_value(int,char);
int mat_size_row(int*);
int find_value_in_array(int*,int,int,int);
int count_value_in_array(int*,int,int,int);
void mat2d_prod(double*, int, int, double*, int, int, double*, int, int);
double mat2d_sum_row(double*, int, int);
void sort_array(double*, int);
void sort_int_array(int*, int);
void sort_array_ids(double*, int*, int);
void obtain_corr_mat_mag_par2(int, int, double*, int, int, int*, int*, int, int, int*, int*, int*, int, int);
void obtain_convex_h_xy(double*, int, int*, double*, int);
void eval_anomaly(double*, double*, int, int, double*, int, int);

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
    int ncluster, ncluster1, ncluster2, ncluster3, ncorr_col;
    int c2start, c3start;
    int *data, *magmom;
    double *E, *Ef;
    int ndata;
    double *corr_mat;
    double near_convh_cutoff;
    int i, j, k, l, nj, nk, ctn;
    char buff_line[200], dummy[200], param_name[100];
    int param_name_len;

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
    else
        printf("Single CPU computing\n");

    // anomaly detection parameters
    int use_anomaly_detection;
    double anomaly_score_crit;
    int knn;

    FILE *fp;
    char dirname[100];
    char paramfilename[100];
    char loadfilename[100];
    strcpy(dirname,"dir_inputs");
    sprintf(paramfilename,"%s/param.dat",dirname);
    fp = fopen(paramfilename, "r");
    while(fgets(buff_line,sizeof(buff_line),fp) != NULL) {
        if (buff_line[0] == '#') {
            //            printf("Comment line: %s",buff_line);
            continue;
        }
        strcpy(param_name,"np");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&np);
            if (rank == MASTER)
                printf("np = %d\n",np);
        }
        strcpy(param_name,"neighbor_num");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&neighbor_num);
            if (rank == MASTER)
                printf("neighbor_num = %d\n",neighbor_num);
        }
        strcpy(param_name,"ncluster1");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&ncluster1);
            if (rank == MASTER)
                printf("ncluster1 = %d\n",ncluster1);
        }
        strcpy(param_name,"ncluster2");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&ncluster2);
            if (rank == MASTER)
                printf("ncluster2 = %d\n",ncluster2);
        }
        strcpy(param_name,"ncluster3");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&ncluster3);
            if (rank == MASTER)
                printf("ncluster3 = %d\n",ncluster3);
        }
        strcpy(param_name,"ndata");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&ndata);
            if (rank == MASTER)
                printf("ndata = %d\n",ndata);
        }
        strcpy(param_name,"near_convh_cutoff");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %lf",dummy,&near_convh_cutoff);
            if (rank == MASTER)
                printf("near_convh_cutoff = %f\n",near_convh_cutoff);
        }
        strcpy(param_name,"use_anomaly_detection");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&use_anomaly_detection);
            if (rank == MASTER)
                printf("use_anomaly_detection = %d\n",use_anomaly_detection);
        }
        strcpy(param_name,"knn");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&knn);
            if (rank == MASTER)
                printf("knn = %d\n",knn);
        }
        strcpy(param_name,"anomaly_score_crit");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %lf",dummy,&anomaly_score_crit);
            if (rank == MASTER)
                printf("anomaly_score_crit = %f\n",anomaly_score_crit);
        }
    }
    fclose(fp);

    nlist           = (int*) malloc(np*neighbor_num*sizeof(int));
    map_to_cluster1 = (int*) malloc(np*sizeof(int));
    map_to_cluster2 = (int*) malloc(np*neighbor_num*sizeof(int));
    map_to_cluster3 = (int*) malloc(np*neighbor_num*neighbor_num*sizeof(int));
    sprintf(loadfilename,"%s/nlist.dat",dirname);
    load_int_mat(nlist, np, neighbor_num, 1, loadfilename);
    sprintf(loadfilename,"%s/map_to_cluster1.dat",dirname);
    load_int_array(map_to_cluster1, np, loadfilename);
    sprintf(loadfilename,"%s/map_to_cluster2.dat",dirname);
    load_int_mat(map_to_cluster2, np, neighbor_num, 1, loadfilename);
    sprintf(loadfilename,"%s/map_to_cluster3.dat",dirname);
    load_int_mat(map_to_cluster3, np, neighbor_num, neighbor_num, loadfilename);
    
    nfu       = np/2;
    ncluster  = 1+ncluster1+ncluster2+ncluster3;    //1 is for 0-cluster
    ncorr_col = 1 + 8*ncluster1 + 30*ncluster2 + 80*ncluster3;
    c2start = 8*ncluster1 + 1;      // Don't need -1: cluster2 idx starts from 0
    c3start = 30*ncluster2 + c2start;
    
    // load LNO occupation variable data
    
    data   = (int*) malloc(ndata*np*sizeof(int));
    magmom = (int*) malloc(ndata*np*sizeof(int));
    E      = (double*) malloc(ndata*sizeof(double));
    Ef     = (double*) malloc(ndata*sizeof(double));
    sprintf(loadfilename,"%s/data_orig.dat",dirname);
    load_int_mat(data,ndata,np,1,loadfilename);
    sprintf(loadfilename,"%s/magmom_orig.dat",dirname);
    load_int_mat(magmom,ndata,np,1,loadfilename);
    sprintf(loadfilename,"%s/E_orig.dat",dirname);
    load_double_array(E,ndata,loadfilename);
    sprintf(loadfilename,"%s/Ef_orig.dat",dirname);
    load_double_array(Ef,ndata,loadfilename);
    
    int nL[ndata];
    int nN[ndata];
    int nM[ndata];
    int nC[ndata];
    int nV[ndata];
    for (i=0;i<ndata;i++) {
        nL[i] = count_value_in_array(data,i,np,2);
        nN[i] = count_value_in_array(data,i,np,1);
        nM[i] = count_value_in_array(data,i,np,0);
        nC[i] = count_value_in_array(data,i,np,-1);
        nV[i] = count_value_in_array(data,i,np,-2);
    }
    // Update Ef, based on the end states: LNO and LMO
    int nM_min, nM_max;
    double E_nM_min, E_nM_max;
    nM_min = 0;
    nM_max = nfu;
    for (i=0;i<ndata;i++)
        if (nM[i] == nM_min) {
            E_nM_min = E[i];
            break;
        }
    for (i=0;i<ndata;i++)
        if (nM[i] == nM_max) {
            E_nM_max = E[i];
            break;
        }
    for (i=0;i<ndata;i++) {
        if (nM[i]==nM_min && E[i] < E_nM_min)
            E_nM_min = E[i];
        if (nM[i]==nM_max && E[i] < E_nM_max)
            E_nM_max = E[i];
    }
    for (i=0;i<ndata;i++)
        *(Ef+i) = *(E+i) - nN[i]*1.0/nfu*E_nM_min - nM[i]*1.0/nfu*E_nM_max;
    
    corr_mat = (double*) malloc(ndata*ncorr_col*sizeof(double));
    // obtain the correlation matrix through parallel computing
    row_offset = (int) ceil(np*1.0/numprocs);
    double tmp_corr_mat_trial[numprocs][ncorr_col];
    double tmp_corr_mat_trial_task[ncorr_col];
    double corr_mat_trial[ncorr_col];
    int data_trial[np], magmom_trial[np];
    
    if (rank == MASTER) {
        for (i=0;i<ndata;i++) {
            for (j=0;j<np;j++) {
                data_trial[j] = *(data + i*np+j);
                magmom_trial[j] = *(magmom + i*np+j);
            }
            row_ini = rank * row_offset;
            row_end = row_ini +row_offset;
            printf("ri=%d, re=%d, ncorr_col=%d, np=%d, neighbor_num=%d, c2start=%d, c3start=%d\n",row_ini,row_end,ncorr_col,np,neighbor_num,c2start, c3start);
            obtain_corr_mat_mag_par2(row_ini, row_end, corr_mat_trial, 1, ncorr_col, data_trial, magmom_trial, c2start, c3start, nlist, map_to_cluster2, map_to_cluster3, np, neighbor_num);
            for (j=1;j<numprocs;j++) {
                mtype = FROM_WORKER;
                MPI_Recv(&tmp_corr_mat_trial[j][0], ncorr_col, MPI_DOUBLE, j, mtype, MPI_COMM_WORLD, &status);
            }
            for (j=1;j<ncorr_col;j++) {     //corr_mat for empty cluster from workers must be excluded to avoid accumulation over workers
                for (k=1;k<numprocs;k++)
                    corr_mat_trial[j] += tmp_corr_mat_trial[k][j];
                *(corr_mat + i*ncorr_col + j) = corr_mat_trial[j];
            }
            *(corr_mat + i*ncorr_col + 0) = 1.0;
            printf("%d-th data set is converted to corr_mat.\n",i);
        }
        
        strcpy(dirname,"dir_result");
        struct stat st = {0};
        if (stat(dirname, &st) == -1) {
            mkdir(dirname,0777);
            printf("Created a new directory for storing result.\n");
        }

        char tmpwritefilename[200];
        sprintf(tmpwritefilename,"%s/tmp_corr_mat.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<ndata;i++) {
            for (j=0;j<ncorr_col;j++)
                fprintf(fp,"%f ",*(corr_mat + i*ncorr_col + j));
            fprintf(fp,"\n");
        }
        fclose(fp);

        // post data processing: excluding degenerate cases
        int nondegenerate[ndata];
        for (i=0;i<ndata;i++)
            nondegenerate[i] = 1;
        for (i=0;i<ndata;i++) {
            for (j=0;j<i;j++) {
                ctn = 1;
                for (l=0;l<ncorr_col;l++)
                    if (*(corr_mat+j*ncorr_col+l) != *(corr_mat+i*ncorr_col+l)) {
                        ctn = 0;
                        break;
                    }
                if (ctn == 1)
                    if (*(Ef+j) < *(Ef+i))
                        nondegenerate[i] = 0;
                    else
                        nondegenerate[j] = 0;
            }
        }
        sprintf(tmpwritefilename,"%s/tmp_nondegenerate.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<ndata;i++)
            fprintf(fp,"%d\n",nondegenerate[i]);
        fclose(fp);

        int n_nondegenerate = 0;
        for (i=0;i<ndata;i++)
            if (nondegenerate[i] == 1)
                n_nondegenerate++;
        double *corr_mat_u;
        corr_mat_u = (double*) malloc(n_nondegenerate*ncorr_col*sizeof(double));
        double Ef_u[n_nondegenerate];
        int nM_u[n_nondegenerate];
        l = 0;
        for (i=0;i<ndata;i++)
            if (nondegenerate[i] == 1) {
                for (j=0;j<ncorr_col;j++)
                    *(corr_mat_u+l*ncorr_col+j) = *(corr_mat+i*ncorr_col+j);
                Ef_u[l] = *(Ef+i);
                nM_u[l] = nM[i];
                l++;
            }
        sprintf(tmpwritefilename,"%s/tmp_Ef_u.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_nondegenerate;i++)
            fprintf(fp,"%f\n",Ef_u[i]);
        fclose(fp);
        sprintf(tmpwritefilename,"%s/tmp_nM_u.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_nondegenerate;i++)
            fprintf(fp,"%d\n",nM_u[i]);
        fclose(fp);
        printf("Denegrated data were filtered out; n_nondegenerate: %d\n\n",n_nondegenerate);
        
        // Apply anomaly detection
        double anomaly_score[n_nondegenerate];
        int anomaly_status[n_nondegenerate];
        int n_normal = 0;
        double max_anomaly_score = 0.0;
        for (i=0;i<n_nondegenerate;i++)
            anomaly_status[i] = 0;
        eval_anomaly(anomaly_score, corr_mat_u, n_nondegenerate, ncorr_col, Ef_u, 1, knn);
        for (i=0;i<n_nondegenerate;i++) {
            printf("anomaly_score[%d] = %f\n",i,anomaly_score[i]);
            if (anomaly_score[i] > max_anomaly_score)
                max_anomaly_score = anomaly_score[i];
        }
        printf("max_anomaly_score = %f\n", max_anomaly_score);
        if (use_anomaly_detection == 0)
            anomaly_score_crit = max_anomaly_score*1.1;
        for (i=0;i<n_nondegenerate;i++) {
            if (anomaly_score[i] < anomaly_score_crit) {
                anomaly_status[i] = 1;
                n_normal++;
            }
        }
        double *corr_mat_ua;
        corr_mat_ua = (double*) malloc(n_normal*ncorr_col*sizeof(double));
        double Ef_ua[n_normal];
        int nM_ua[n_normal];
        l = 0;
        for (i=0;i<n_nondegenerate;i++)
            if (anomaly_status[i] == 1) {
                for (j=0;j<ncorr_col;j++)
                    *(corr_mat_ua+l*ncorr_col+j) = *(corr_mat_u+i*ncorr_col+j);
                Ef_ua[l] = Ef_u[i];
                nM_ua[l] = nM_u[i];
                l++;
            }
        sprintf(tmpwritefilename,"%s/tmp_Ef_ua.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_normal;i++)
            fprintf(fp,"%f\n",Ef_ua[i]);
        fclose(fp);
        sprintf(tmpwritefilename,"%s/tmp_nC_ua.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_normal;i++)
            fprintf(fp,"%d\n",nM_ua[i]);
        fclose(fp);
        sprintf(tmpwritefilename,"%s/tmp_corr_mat_ua.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_normal;i++) {
            for (j=0;j<ncorr_col;j++)
                fprintf(fp,"%f ",*(corr_mat_ua+i*ncorr_col));
            fprintf(fp,"\n");
        }
        fclose(fp);
        printf("Anomaly data sets were filtered out; n_normal: %d\n\n",n_normal);
        
        // Generate convex hull
        int n_bin = np/2+1;
        double ef_conv[n_bin];
        obtain_convex_h_xy(ef_conv,n_bin,nM_ua,Ef_ua,n_normal);
        sprintf(tmpwritefilename,"%s/ef_conv_hull.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_bin;i++)
            fprintf(fp,"%d %f\n",i,ef_conv[i]);
        fclose(fp);
        
        // select data sets near convex hull
        int n_corr_mat_ug_row = 0;
        int near_convh[n_normal];
        for (i=0;i<n_normal;i++)
            if (Ef_ua[i] - ef_conv[nM_ua[i]] < near_convh_cutoff*nfu) {
                near_convh[i] = 1;
                n_corr_mat_ug_row++;
            }
            else
                near_convh[i] = 0;
        sprintf(tmpwritefilename,"%s/tmp_near_convh.dat",dirname);
        fp = fopen(tmpwritefilename, "w");
        for (i=0;i<n_normal;i++)
            fprintf(fp,"%d\n",near_convh[i]);
        fclose(fp);
        
        double *corr_mat_ug;
        corr_mat_ug = (double*) malloc(n_corr_mat_ug_row*ncorr_col*sizeof(double));
        double Ef_ug[n_corr_mat_ug_row];
        int nM_ug[n_corr_mat_ug_row];
        l = 0;
        for (i=0;i<n_normal;i++)
            if (near_convh[i] == 1) {
                for (j=0;j<ncorr_col;j++)
                    *(corr_mat_ug+l*ncorr_col+j) = *(corr_mat_ua+i*ncorr_col+j);
                Ef_ug[l] = Ef_ua[i];
                nM_ug[l] = nM_ua[i];
                l++;
            }
        printf("States beyond tol+convex hull were filterd out; n_corr_mat_ug_row: %d\n\n",n_corr_mat_ug_row);

        // select the non-singular columns only
        int n_non_singular_col;
        int col_singularity[ncorr_col];
        col_singularity[0] = 1;
        n_non_singular_col = 1;
        for (i=1;i<ncorr_col;i++) {
            ctn = 0;
            for (j=1;j<n_corr_mat_ug_row;j++)
                if ( *(corr_mat_ug+j*ncorr_col+i) != *(corr_mat_ug+0*ncorr_col+i) ) {
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

        double *corr_mat_ugs;
        corr_mat_ugs = (double*) malloc(n_corr_mat_ug_row*n_non_singular_col*sizeof(double));
        int usefulcorr_col[n_non_singular_col];
        l = 0;
        for (i=0;i<ncorr_col;i++)
            if (col_singularity[i] == 1) {
                for (j=0;j<n_corr_mat_ug_row;j++)
                    *(corr_mat_ugs+j*n_non_singular_col+l) = *(corr_mat_ug+j*n_non_singular_col+i);
                usefulcorr_col[l] = i;
                l++;
            }
        printf("non_singular_col were extracted; n_non_singular_col: %d\n\n",n_non_singular_col);

        // Write result into files
        FILE *fp1, *fp2, *fp3, *fp4;
        char corr_mat_ugs_filename[200];
        char Ef_ug_filename[200];
        char usefulcorr_col_filename[200];
        sprintf(corr_mat_ugs_filename, "%s/corr_mat_ugs.dat",dirname);
        sprintf(Ef_ug_filename, "%s/Ef_ug.dat",dirname);
        sprintf(usefulcorr_col_filename, "%s/usefulcorr_col.dat",dirname);
        sprintf(paramfilename,"%s/result_param.dat",dirname);
        fp1 = fopen(corr_mat_ugs_filename,"w");
        fp2 = fopen(Ef_ug_filename,"w");
        fp3 = fopen(usefulcorr_col_filename,"w");
        fp4 = fopen(paramfilename,"w");
        for (i=0;i<n_corr_mat_ug_row;i++) {
            for (j=0;j<n_non_singular_col;j++)
                fprintf(fp1,"%f ",*(corr_mat_ugs+i*n_non_singular_col+j));
            fprintf(fp1,"\n");
            fprintf(fp2,"%f\n",Ef_ug[i]);
        }
        for (i=0;i<n_non_singular_col;i++)
            fprintf(fp3,"%d\n",usefulcorr_col[i]);
        fprintf(fp4,"n_corr_mat_ug_row\t%d\n",n_corr_mat_ug_row);
        fprintf(fp4,"n_non_singular_col\t%d\n",n_non_singular_col);
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);

        printf("The Correlation Matrix, the Ef, and the usefulcorr_col were regularized...\n\n");

    }
    else {
        for (i=0;i<ndata;i++) {
            row_ini = rank * row_offset;
            row_end = row_ini + row_offset;
            if (rank == numprocs-1)
                row_end = np;  // last processor till the last row, -1 is taken inside obtainerr2_par
            obtain_corr_mat_mag_par2(row_ini, row_end, tmp_corr_mat_trial_task, 1, ncorr_col, data_trial, magmom_trial, c2start, c3start, nlist, map_to_cluster2, map_to_cluster3, np, neighbor_num);
            mtype = FROM_WORKER;
            MPI_Send(&tmp_corr_mat_trial_task, ncorr_col, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
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

int find_value_in_array(int *A, int target_rowA, int ncolA, int target_value) {
    int i;
    
    for (i=0;i<ncolA;i++)
        if (*(A+ncolA*target_rowA+i)==target_value)
            return(i);
    if (i==ncolA)
        return(-1);
}

int count_value_in_array(int *A, int target_rowA, int ncolA, int target_value) {
    int i;
    int count=0;
    
    for (i=0;i<ncolA;i++)
        if (*(A+ncolA*target_rowA+i)==target_value)
            count++;

    return(count);
}


























