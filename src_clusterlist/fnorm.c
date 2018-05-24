#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double fnorm(double *a)
{
    double a1, a2, a3, a_norm;
    a1 = *(a+0);
    a2 = *(a+1);
    a3 = *(a+2);
    a_norm = sqrt(a1*a1 + a2*a2 + a3*a3);

    return(a_norm);
}
