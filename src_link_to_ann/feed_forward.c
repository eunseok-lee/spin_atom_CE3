#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double sigmoidf(double);
double dot_vectors(double*, double*, int);

void feed_forward(double *hidden_layer, double *output_layer, double *input_vector_w_bias, double *hidden_outputs, double *outputs, int num_hidden_neuron, int hidden_layer_colsize, int num_output_neuron) {

    /*
     hidden_layer: num_hidden_neuron x (datainput_colsize+1)
     output_layer: num_output_neuron x (num_hidden_neuron+1)
     input_vector_w_bias: (datainput_colsize+1)
     hidden_outputs: num_hidden_neuron
     outputs: num_output_neuron
     */

    int i,j ;
    double tmp_hidden_layer[hidden_layer_colsize];

    // process for the hidden_layer
    for (i=0; i<num_hidden_neuron; i++) {
        for (j=0; j<hidden_layer_colsize; j++)
            tmp_hidden_layer[j] = *(hidden_layer + hidden_layer_colsize*i + j);
        *(hidden_outputs + i) = sigmoidf(dot_vectors(tmp_hidden_layer, input_vector_w_bias, hidden_layer_colsize));
//        printf("%d-hidden_neuron, dot_vectors = %f, outputs = %f\n",i,dot_vectors(tmp_hidden_layer, input_vector_w_bias, hidden_layer_colsize),*(hidden_outputs + i));
    }

    // process for the output_layer
    double input_for_output_layer[num_hidden_neuron+1];
    double tmp_output_layer[num_hidden_neuron+1];
    for (i=0; i<num_hidden_neuron; i++)
        input_for_output_layer[i] = *(hidden_outputs + i);
    input_for_output_layer[num_hidden_neuron] = 1.0;
    for (i=0; i<num_output_neuron; i++) {
        for (j=0; j<(num_hidden_neuron+1); j++)
            tmp_output_layer[j] = *(output_layer + (num_hidden_neuron+1)*i + j);
        *(outputs + i) = sigmoidf(dot_vectors(tmp_output_layer, input_for_output_layer, num_hidden_neuron+1));
//        printf("%d-output_neuron, dot_vectors = %f, outputs = %f\n",i,dot_vectors(tmp_output_layer, input_for_output_layer, num_hidden_neuron+1), *(outputs + i));
    }
}
