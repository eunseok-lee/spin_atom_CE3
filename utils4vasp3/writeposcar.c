#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define sp_Li 2
#define sp_Va -2
#define sp_Mn 1
#define sp_Ni -1

int main(int argc, char *argv[]) {
    
    FILE *fp, *fp1, *fp2, *fp3, *fp4;
    int i, np, howmanyLi, howmanyVa, howmanyMn, howmanyNi;
    double *rp, *rpO, *kv;
    double tmp;
    char datfilename[100];
    char poscarfilename[100];
    char buff_line[200];
    
    fp = fopen("rp.dat","r");
    np = 0;
    while(fgets(buff_line,sizeof(buff_line),fp) != NULL) {
        if (buff_line[0] == '\n') {
            continue;
        }
        else
            np++;
    }
    fclose(fp);
    printf("np = %d\n", np);
    
    rp = (double*) calloc (np*3,sizeof(double));
    rpO= (double*) calloc (np*3,sizeof(double));
    kv = (double*) calloc (3*3,sizeof(double));
    int array[np];
    
    printf("rp data is read\n");
    fp1 = fopen("rp.dat","r");
    for (i=0;i<np;i++) {
        fscanf(fp1, "%lf", &tmp);
        *(rp+3*i) = tmp;
        fscanf(fp1, "%lf", &tmp);
        *(rp+3*i+1) = tmp;
        fscanf(fp1, "%lf", &tmp);
        *(rp+3*i+2) = tmp;
        printf("%d\t%f %f %f\n",i,*(rp+3*i),*(rp+3*i+1),*(rp+3*i+2));
    }
    fclose(fp1);
    
    printf("rpO data is read\n");
    fp2 = fopen("rpO.dat","r");
    for (i=0;i<np;i++) {
        fscanf(fp2, "%lf", &tmp);
        *(rpO+3*i) = tmp;
        fscanf(fp2, "%lf", &tmp);
        *(rpO+3*i+1) = tmp;
        fscanf(fp2, "%lf", &tmp);
        *(rpO+3*i+2) = tmp;
        printf("%d\t%f %f %f\n",i,*(rpO+3*i),*(rpO+3*i+1),*(rpO+3*i+2));
    }
    fclose(fp2);
    
    sprintf(datfilename,"%s",argv[1]);
    printf("Reading the datafile, %s\n",datfilename);
    fp3 = fopen(datfilename,"r");
    for (i=0;i<np;i++) {
        fscanf(fp3, "%d", &array[i]);
        printf("%d\n",array[i]);
    }
    fclose(fp3);
    
    // Analyze the loaded array
    for (i=0;i<np;i++) {
        switch (array[i]) {
            case sp_Li:
                howmanyLi++;
                break;
            case sp_Mn:
                howmanyMn++;
                break;
            case sp_Ni:
                howmanyNi++;
                break;
            case sp_Va:
                howmanyVa++;
                break;
            default:
                printf("[E] %d-element: %d doesn't match to Li/Va/Mn/Ni\n",i,array[i]);
        }
    }
    
    printf("kv data is read\n");
    fp4 = fopen("kv.dat","r");
    for (i=0;i<3;i++) {
        fscanf(fp4, "%lf", &tmp);
        *(kv+3*i) = tmp;
        fscanf(fp4, "%lf", &tmp);
        *(kv+3*i+1) = tmp;
        fscanf(fp4, "%lf", &tmp);
        *(kv+3*i+2) = tmp;
        printf("%d\t%f %f %f\n",i,*(kv+3*i),*(kv+3*i+1),*(kv+3*i+2));
    }
    fclose(fp4);

    sprintf(poscarfilename,"POSCAR");
    fp = fopen(poscarfilename,"w");
    fprintf(fp,"LixMnyNizO2\n");
    fprintf(fp,"1.0\n");
    fprintf(fp,"%.4f %.4f %.4f\n",*(kv+3*0),*(kv+3*0+1),*(kv+3*0+2));
    fprintf(fp,"%.4f %.4f %.4f\n",*(kv+3*1),*(kv+3*1+1),*(kv+3*1+2));
    fprintf(fp,"%.4f %.4f %.4f\n",*(kv+3*2),*(kv+3*2+1),*(kv+3*2+2));
    if (howmanyLi > 0)
        fprintf(fp, "Li ");
    if (howmanyMn > 0)
        fprintf(fp, "Mn ");
    if (howmanyNi > 0)
        fprintf(fp, "Ni ");
    fprintf(fp,"O\n");
    if (howmanyLi > 0)
        fprintf(fp, "%d ", howmanyLi);
    if (howmanyMn > 0)
        fprintf(fp, "%d ", howmanyMn);
    if (howmanyNi > 0)
        fprintf(fp, "%d ", howmanyNi);
    fprintf(fp,"%d\n", np);

    fprintf(fp,"Direct\n");

    if (howmanyLi > 0)
        for (i=0;i<np;i++)
            if (array[i]==sp_Li)
                fprintf(fp, "%.4f %.4f %.4f Li\n",*(rp+3*i),*(rp+3*i+1),*(rp+3*i+2));
    
    if (howmanyMn > 0)
        for (i=0;i<np;i++)
            if (array[i]==sp_Mn)
                fprintf(fp, "%.4f %.4f %.4f Mn\n",*(rp+3*i),*(rp+3*i+1),*(rp+3*i+2));

    if (howmanyNi > 0)
        for (i=0;i<np;i++)
            if (array[i]==sp_Ni)
                fprintf(fp, "%.4f %.4f %.4f Ni\n",*(rp+3*i),*(rp+3*i+1),*(rp+3*i+2));
    
    for (i=0;i<np;i++)
        fprintf(fp, "%.4f %.4f %.4f O\n",*(rpO+3*i),*(rpO+3*i+1),*(rpO+3*i+2));
    
    fclose(fp);
}
