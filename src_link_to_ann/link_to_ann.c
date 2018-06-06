#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <mpi.h>

void feed_forward(double*, double*, double*, double*, double*, int, int, int);
void load_double_mat(double*, int, int, int, char*);
double sigmoidf(double);
double dot_vectors(double*, double*, int);

#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

int main(int argc, char **argv) {

    int i, j, k;
    double err2;
    int nfu, disp_freq;

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

    // load inputs
    FILE *fp, *fp2, *fp3, *tmp_fp;
    char buff_line[200], dummy[200], param_name[100];
    int param_name_len;
    char dirname[100]="dir_inputs";
    char paramfilename[100];
    sprintf(paramfilename,"%s/param.dat",dirname);
    fp = fopen(paramfilename, "r");
    char datainput_filename[200];
    char targets_filename[200];
    int datainput_rowsize, datainput_colsize;
    int targets_rowsize, targets_colsize;
    int iter, max_iter, num_hidden_neuron, num_output_neuron, data_idx;

    while(fgets(buff_line,sizeof(buff_line),fp) != NULL) {
        if (buff_line[0] == '#') {
            //            printf("Comment line: %s",buff_line);
            continue;
        }
        strcpy(param_name,"datainput_filename");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %s",dummy,datainput_filename);
            if (rank == MASTER)
                printf("datainput_filename = %s\n",datainput_filename);
        }
        strcpy(param_name,"datainput_rowsize");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&datainput_rowsize);
            if (rank == MASTER)
                printf("datainput_rowsize = %d\n",datainput_rowsize);
        }
        strcpy(param_name,"datainput_colsize");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&datainput_colsize);
            if (rank == MASTER)
                printf("datainput_colsize = %d\n",datainput_colsize);
        }
        strcpy(param_name,"targets_filename");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %s",dummy,targets_filename);
            if (rank == MASTER)
                printf("targets_filename = %s\n",targets_filename);
        }
        strcpy(param_name,"targets_rowsize");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&targets_rowsize);
            if (rank == MASTER)
                printf("targets_rowsize = %d\n",targets_rowsize);
        }
        strcpy(param_name,"targets_colsize");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&targets_colsize);
            if (rank == MASTER)
                printf("targets_colsize = %d\n",targets_colsize);
        }
        strcpy(param_name,"max_iter");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&max_iter);
            if (rank == MASTER)
                printf("max_iter = %d\n",max_iter);
        }
        strcpy(param_name,"num_hidden_neuron");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&num_hidden_neuron);
            if (rank == MASTER)
                printf("num_hidden_neuron = %d\n",num_hidden_neuron);
        }
        strcpy(param_name,"disp_freq");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&disp_freq);
            if (rank == MASTER)
                printf("disp_freq = %d\n",disp_freq);
        }
        strcpy(param_name,"nfu");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&nfu);
            if (rank == MASTER)
                printf("nfu = %d\n",nfu);
        }
    }
    fclose(fp);
    num_output_neuron = targets_colsize;

    char loadfilename[100];
    double *datainput;
    datainput = (double*) malloc(datainput_rowsize*datainput_colsize*sizeof(double));
    sprintf(loadfilename,"%s/%s",dirname,datainput_filename);
    load_double_mat(datainput,datainput_rowsize,datainput_colsize,1,loadfilename);

    double *targets;
    targets = (double*) malloc(targets_rowsize*targets_colsize*sizeof(double));
    sprintf(loadfilename,"%s/%s",dirname,targets_filename);
    load_double_mat(targets,targets_rowsize,targets_colsize,1,loadfilename);

    // initialize the hidden layers and output layer
    double *hidden_layer, *output_layer;
    hidden_layer = (double*) malloc(num_hidden_neuron*(datainput_colsize+1)*sizeof(double));
    for (i=0; i<num_hidden_neuron*(datainput_colsize+1); i++)
        hidden_layer[i] = (double) rand()*1.0/RAND_MAX;
    output_layer = (double*) malloc(num_output_neuron*(num_hidden_neuron+1)*sizeof(double));
    for (i=0; i<num_output_neuron*(num_hidden_neuron+1); i++)
        output_layer[i] = (double) rand()*1.0/RAND_MAX;

/*    tmp_fp = fopen("inter_result.dat","w");
    fprintf(tmp_fp, "Hidden Layer\n");
    for (i=0; i<num_hidden_neuron*(datainput_colsize+1); i++)
        fprintf(tmp_fp, "%f\n", hidden_layer[i]);
    fprintf(tmp_fp, "Output Layer\n");
    for (i=0; i<num_output_neuron*(num_hidden_neuron+1); i++)
        fprintf(tmp_fp, "%f\n", output_layer[i]);
*/
    double input_vector[datainput_colsize];
    double input_vector_w_bias[datainput_colsize+1];
    double target_vector[targets_colsize];
    double output_deltas[targets_colsize];
    double hidden_outputs[num_hidden_neuron];
    double hidden_deltas[num_hidden_neuron];
    double hidden_outputs_w_bias[num_hidden_neuron+1];
    double outputs[num_output_neuron];
    double hidden_output[num_hidden_neuron];
    double output_layer_vector[num_output_neuron];

    // backpropagation
    for (iter=0; iter<max_iter; iter++) {
        err2 = 0.0;
        for (data_idx=0; data_idx<datainput_rowsize; data_idx++) {
            // assign input_vector and target_vector
            for (i=0; i<datainput_colsize; i++) {
                input_vector[i] = datainput[datainput_colsize*data_idx+i];
                input_vector_w_bias[i] = input_vector[i];
            }
            input_vector_w_bias[datainput_colsize] = 1.0;
            for (i=0; i<targets_colsize; i++) {
                target_vector[i] = targets[targets_colsize*data_idx+i]/nfu;     // in case of very large target values, *1/nfu
            }

            feed_forward(hidden_layer, output_layer, input_vector_w_bias, hidden_outputs, outputs, num_hidden_neuron, datainput_colsize+1, num_output_neuron);

            for (i=0; i<targets_colsize; i++) {     //num_output_neuron=targets_colsize
                output_deltas[i] = outputs[i]*(1-outputs[i])*(outputs[i]-target_vector[i]);
                err2 += pow(outputs[i]-target_vector[i],2.0);
            }

            // adjust weights for the output_layer
            for (i=0; i<num_hidden_neuron; i++)
                hidden_outputs_w_bias[i] = hidden_outputs[i];
            hidden_outputs_w_bias[num_hidden_neuron] = 1.0;
            for (i=0; i<num_output_neuron; i++)
                for (j=0; j<=num_hidden_neuron; j++)
                    *(output_layer + (num_hidden_neuron+1)*i + j) -= output_deltas[i]*hidden_outputs_w_bias[j];

            // back-propagate error to the hidden layer
            for (i=0; i<num_hidden_neuron; i++) {
                for (j=0; j<num_output_neuron; j++)
                    output_layer_vector[j] = *(output_layer + (num_hidden_neuron+1)*j + i);
                hidden_deltas[i] = hidden_output[i]*(1-hidden_output[i])*dot_vectors(output_deltas,output_layer_vector,num_output_neuron);
            }

            // adjust weights for the hidden layer
            for (i=0; i<num_hidden_neuron; i++)
                for (j=0; j<(datainput_colsize+1); j++)
                    *(hidden_layer + (datainput_colsize+1)*i + j) -= hidden_deltas[i]*input_vector_w_bias[j];
        }
        if (iter%disp_freq == 0)
            printf("iter=%d, error_mse=%e\n",iter,sqrt(err2/datainput_rowsize));
    }

    // write result into files
    char hidden_layer_filename[200];
    char output_layer_filename[200];
    strcpy(dirname, "dir_result");
    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        mkdir(dirname,0777);
        printf("Created a new directory for storing result.\n");
    }
    sprintf(hidden_layer_filename, "%s/hidden_layer.dat", dirname);
    sprintf(output_layer_filename, "%s/output_layer.dat", dirname);
    fp2 = fopen(hidden_layer_filename, "w");
    fp3 = fopen(output_layer_filename, "w");
    for (i=0; i<num_hidden_neuron; i++) {
        for (j=0; j<(datainput_colsize+1); j++)
            fprintf(fp2, "%e ", *(hidden_layer + (datainput_colsize+1)*i + j));
        fprintf(fp2, "\n");
    }
    for (i=0; i<num_output_neuron; i++) {
        for (j=0; j<(num_hidden_neuron+1); j++)
            fprintf(fp3, "%e ", *(output_layer + (num_hidden_neuron+1)*i + j));
        fprintf(fp3, "\n");
    }
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
