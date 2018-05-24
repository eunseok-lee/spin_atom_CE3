#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define sLi 2
#define sNi 1
#define sMn 0
#define sCo -1
#define sVa -2

int main(int argc, char **argv) {
    
    FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9, *fp10;
    int i, j, k, np, npL, npN, npM, npC;
    int setn, ctr;
    int ndata, data_ini_id;
    double *rp;
    double tmp;
    int tmp1;
    double tmp2, tmp3, tmp4, tmp5;
    char dirname[100];
    char datfilename[100];
    char buff_line[200];

    double max_dist_cutoff = 0.08;
    double ind_dist_cutoff = 0.001;
    double dist_cutoff;
    double dij[3];
    double dist;
    
    double magLi = atof(argv[1]);
    double magNi = atof(argv[2]);
    double magMn = atof(argv[3]);
    double magCo = atof(argv[4]);

    for (i=1;i<100000;i++) {
        sprintf(dirname,"data%d",i);
        struct stat st = {0};
        if (stat(dirname, &st) == 0) {
            ndata++;
            if (ndata == 1) {
                data_ini_id = i;
                printf("data_ini_id = %d\n",data_ini_id);
            }
        }
        else if (stat(dirname, &st) == -1 && ndata > 0)
            break;
    }
    
    printf("ndata = %d\n",ndata);
    fp1 = fopen("rp.dat","r");
    np = 0;
    while(fgets(buff_line,sizeof(buff_line),fp1) != NULL) {
        if (buff_line[0] == '\n') {
            continue;
        }
        else
            np++;
    }
    fclose(fp1);
    printf("np = %d\n", np);
    printf("magLi = %f\n",magLi);
    printf("magNi = %f\n",magNi);
    printf("magMn = %f\n",magMn);
    printf("magCo = %f\n",magCo);

    rp = (double*) calloc (np*3,sizeof(double));
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
    printf("----------------------------------------\n");
    
    int data[ndata][np];
    int magmom[ndata][np];
    int nL[ndata];
    int nM[ndata];
    int nN[ndata];
    int nC[ndata];
    int nV[ndata];
    double E[ndata];
    double Ef[ndata];
    double *rpL, *rpM, *rpN, *rpC;
    double *mL, *mM, *mN, *mC;
    int *Lid, *Mid, *Nid, *Cid;
    
    // fill up data and magmom
    for (i=0;i<ndata;i++)
        for (j=0;j<np;j++) {
            data[i][j] = sVa; //default, also value for Va
            magmom[i][j] = 0;   //default, also value for Va
        }

    for (setn=0;setn<ndata;setn++) {
        printf("data%d is being processed\n",setn+data_ini_id);
        sprintf(datfilename,"data%d/posLi.dat",setn+data_ini_id);
        fp2 = fopen(datfilename,"r");
        sprintf(datfilename,"data%d/posNi.dat",setn+data_ini_id);
        fp3 = fopen(datfilename,"r");
        sprintf(datfilename,"data%d/posMn.dat",setn+data_ini_id);
        fp4 = fopen(datfilename,"r");
        sprintf(datfilename,"data%d/posCo.dat",setn+data_ini_id);
        fp5 = fopen(datfilename,"r");

        npL = 0;
        while(fgets(buff_line,sizeof(buff_line),fp2) != NULL) {
            if (buff_line[0] == '\n') {
                continue;
            }
            else
                npL++;
        }
        printf("npL = %d\n", npL);
        npN = 0;
        while(fgets(buff_line,sizeof(buff_line),fp3) != NULL) {
            if (buff_line[0] == '\n') {
                continue;
            }
            else
                npN++;
        }
        printf("npN = %d\n", npN);
        npM = 0;
        while(fgets(buff_line,sizeof(buff_line),fp4) != NULL) {
            if (buff_line[0] == '\n') {
                continue;
            }
            else
                npM++;
        }
        printf("npM = %d\n", npM);
        npC = 0;
        while(fgets(buff_line,sizeof(buff_line),fp5) != NULL) {
            if (buff_line[0] == '\n') {
                continue;
            }
            else
                npC++;
        }
        printf("npC = %d\n", npC);
        
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
        fclose(fp5);

        rpL = (double*) calloc(npL*3,sizeof(double));
        rpN = (double*) calloc(npN*3,sizeof(double));
        rpM = (double*) calloc(npM*3,sizeof(double));
        rpC = (double*) calloc(npC*3,sizeof(double));

        Lid = (int*) malloc(npL*sizeof(int));
        Nid = (int*) malloc(npN*sizeof(int));
        Mid = (int*) malloc(npM*sizeof(int));
        Cid = (int*) malloc(npC*sizeof(int));
        for (i=0;i<npL;i++)
            *(Lid+i) = -1;
        for (i=0;i<npM;i++)
            *(Mid+i) = -1;
        for (i=0;i<npN;i++)
            *(Nid+i) = -1;
        for (i=0;i<npC;i++)
            *(Cid+i) = -1;

        printf("posL data is read\n");
        sprintf(datfilename,"data%d/posLi.dat",setn+data_ini_id);
        fp2 = fopen(datfilename,"r");
        for (i=0;i<npL;i++) {
            fscanf(fp2, "%lf", &tmp);
            *(rpL+3*i) = tmp;
            fscanf(fp2, "%lf", &tmp);
            *(rpL+3*i+1) = tmp;
            fscanf(fp2, "%lf", &tmp);
            *(rpL+3*i+2) = tmp;
//            printf("%d\t%f %f %f\n",i,*(rpL+3*i),*(rpL+3*i+1),*(rpL+3*i+2));
        }
        fclose(fp2);

        printf("posM data is read\n");
        sprintf(datfilename,"data%d/posMn.dat",setn+data_ini_id);
        fp3 = fopen(datfilename,"r");
        for (i=0;i<npM;i++) {
            fscanf(fp3, "%lf", &tmp);
            *(rpM+3*i) = tmp;
            fscanf(fp3, "%lf", &tmp);
            *(rpM+3*i+1) = tmp;
            fscanf(fp3, "%lf", &tmp);
            *(rpM+3*i+2) = tmp;
//            printf("%d\t%f %f %f\n",i,*(rpM+3*i),*(rpM+3*i+1),*(rpM+3*i+2));
        }
        fclose(fp3);

        printf("posN data is read\n");
        sprintf(datfilename,"data%d/posNi.dat",setn+data_ini_id);
        fp4 = fopen(datfilename,"r");
        for (i=0;i<npN;i++) {
            fscanf(fp4, "%lf", &tmp);
            *(rpN+3*i) = tmp;
            fscanf(fp4, "%lf", &tmp);
            *(rpN+3*i+1) = tmp;
            fscanf(fp4, "%lf", &tmp);
            *(rpN+3*i+2) = tmp;
//            printf("%d\t%f %f %f\n",i,*(rpN+3*i),*(rpN+3*i+1),*(rpN+3*i+2));
        }
        fclose(fp4);

        printf("posC data is read\n");
        sprintf(datfilename,"data%d/posCo.dat",setn+data_ini_id);
        fp5 = fopen(datfilename,"r");
        for (i=0;i<npC;i++) {
            fscanf(fp5, "%lf", &tmp);
            *(rpC+3*i) = tmp;
            fscanf(fp5, "%lf", &tmp);
            *(rpC+3*i+1) = tmp;
            fscanf(fp5, "%lf", &tmp);
            *(rpC+3*i+2) = tmp;
//            printf("%d\t%f %f %f\n",i,*(rpN+3*i),*(rpN+3*i+1),*(rpN+3*i+2));
        }
        fclose(fp5);

        mL = (double*) calloc(npL,sizeof(double));
        mM = (double*) calloc(npM,sizeof(double));
        mN = (double*) calloc(npN,sizeof(double));
        mC = (double*) calloc(npC,sizeof(double));
        
        printf("magL data is read\n");
        sprintf(datfilename,"data%d/magLi.dat",setn+data_ini_id);
        fp6 = fopen(datfilename,"r");
        for (i=0;i<npL;i++) {
            fscanf(fp6, "%d %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
            *(mL+i) = tmp5;
        }
        fclose(fp6);
        
        printf("magM data is read\n");
        sprintf(datfilename,"data%d/magMn.dat",setn+data_ini_id);
        fp7 = fopen(datfilename,"r");
        for (i=0;i<npM;i++) {
            fscanf(fp7, "%d %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
            *(mM+i) = tmp5;
        }
        fclose(fp7);
        
        printf("magN data is read\n");
        sprintf(datfilename,"data%d/magNi.dat",setn+data_ini_id);
        fp8 = fopen(datfilename,"r");
        for (i=0;i<npN;i++) {
            fscanf(fp8, "%d %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
            *(mN+i) = tmp5;
        }
        fclose(fp8);
        
        printf("magC data is read\n");
        sprintf(datfilename,"data%d/magCo.dat",setn+data_ini_id);
        fp9 = fopen(datfilename,"r");
        for (i=0;i<npC;i++) {
            fscanf(fp9, "%d %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
            *(mC+i) = tmp5;
        }
        fclose(fp9);
        
        for (i=0;i<npL;i++) {
            dist_cutoff = ind_dist_cutoff;
            ctr = 0;
            while (dist_cutoff < max_dist_cutoff && Lid[i] < 0) {
                for (j=0;j<np;j++) {
                    dij[0] = pow(fmod(*(rp+3*j+0) - *(rpL+3*i+0) + 0.5, 1.0) - 0.5, 2.0);
                    dij[1] = pow(fmod(*(rp+3*j+1) - *(rpL+3*i+1) + 0.5, 1.0) - 0.5, 2.0);
                    dij[2] = pow(fmod(*(rp+3*j+2) - *(rpL+3*i+2) + 0.5, 1.0) - 0.5, 2.0);
                    dist = sqrt(dij[0] + dij[1] + dij[2]);
                    if (dist < dist_cutoff) {
                        Lid[i] = j;
                        data[setn][j] = sLi;
                        if (*(mL+i) > magLi)
                            magmom[setn][j] = 1;
                        else if (*(mL+i) < -1*magLi)
                            magmom[setn][j] = -1;
                        else
                            magmom[setn][j] = 0;
                        ctr = 1;
                        break;
                    }
                    dist_cutoff = dist_cutoff+ind_dist_cutoff;
                }
            }
            if (ctr == 0)
                printf("[E] Lid[%d] not identified\n",i);
            else
                printf("Lid[%d] = %d\n",i,Lid[i]);
        }
        for (i=0;i<npM;i++) {
            dist_cutoff = ind_dist_cutoff;
            ctr = 0;
            while (dist_cutoff < max_dist_cutoff && Mid[i] < 0) {
                for (j=0;j<np;j++) {
                    dij[0] = pow(fmod(*(rp+3*j+0) - *(rpM+3*i+0) + 0.5, 1.0) - 0.5, 2.0);
                    dij[1] = pow(fmod(*(rp+3*j+1) - *(rpM+3*i+1) + 0.5, 1.0) - 0.5, 2.0);
                    dij[2] = pow(fmod(*(rp+3*j+2) - *(rpM+3*i+2) + 0.5, 1.0) - 0.5, 2.0);
                    dist = sqrt(dij[0] + dij[1] + dij[2]);
                    if (dist < dist_cutoff) {
                        Mid[i] = j;
                        data[setn][j] = sMn;
                        if (*(mM+i) > magMn)
                            magmom[setn][j] = 1;
                        else if (*(mM+i) < -1*magMn)
                            magmom[setn][j] = -1;
                        else
                            magmom[setn][j] = 0;
                        ctr = 1;
                        break;
                    }
                    dist_cutoff = dist_cutoff+ind_dist_cutoff;
                }
            }
            if (ctr == 0)
                printf("[E] Mid[%d] not identified\n",i);
            else
                printf("Mid[%d] = %d\n",i,Mid[i]);
        }
        for (i=0;i<npN;i++) {
            dist_cutoff = ind_dist_cutoff;
            ctr = 0;
            while (dist_cutoff < max_dist_cutoff && Nid[i] < 0) {
                for (j=0;j<np;j++) {
                    dij[0] = pow(fmod(*(rp+3*j+0) - *(rpN+3*i+0) + 0.5, 1.0) - 0.5, 2.0);
                    dij[1] = pow(fmod(*(rp+3*j+1) - *(rpN+3*i+1) + 0.5, 1.0) - 0.5, 2.0);
                    dij[2] = pow(fmod(*(rp+3*j+2) - *(rpN+3*i+2) + 0.5, 1.0) - 0.5, 2.0);
                    dist = sqrt(dij[0] + dij[1] + dij[2]);
                    if (dist < dist_cutoff) {
                        Nid[i] = j;
                        data[setn][j] = sNi;
                        if (*(mN+i) > magNi)
                            magmom[setn][j] = 1;
                        else if (*(mN+i) < -1*magNi)
                            magmom[setn][j] = -1;
                        else
                            magmom[setn][j] = 0;
                        ctr = 1;
                        break;
                    }
                    dist_cutoff = dist_cutoff+ind_dist_cutoff;
                }
            }
            if (ctr == 0)
                printf("[E] Nid[%d] not identified\n",i);
            else
                printf("Nid[%d] = %d\n",i,Nid[i]);
        }
        for (i=0;i<npC;i++) {
            dist_cutoff = ind_dist_cutoff;
            ctr = 0;
            while (dist_cutoff < max_dist_cutoff && Cid[i] < 0) {
                for (j=0;j<np;j++) {
                    dij[0] = pow(fmod(*(rp+3*j+0) - *(rpC+3*i+0) + 0.5, 1.0) - 0.5, 2.0);
                    dij[1] = pow(fmod(*(rp+3*j+1) - *(rpC+3*i+1) + 0.5, 1.0) - 0.5, 2.0);
                    dij[2] = pow(fmod(*(rp+3*j+2) - *(rpC+3*i+2) + 0.5, 1.0) - 0.5, 2.0);
                    dist = sqrt(dij[0] + dij[1] + dij[2]);
                    if (dist < dist_cutoff) {
                        Cid[i] = j;
                        data[setn][j] = sCo;
                        if (*(mC+i) > magCo)
                            magmom[setn][j] = 1;
                        else if (*(mC+i) < -1*magCo)
                            magmom[setn][j] = -1;
                        else
                            magmom[setn][j] = 0;
                        ctr = 1;
                        break;
                    }
                    dist_cutoff = dist_cutoff+ind_dist_cutoff;
                }
            }
            if (ctr == 0)
                printf("[E] Cid[%d] not identified\n",i);
            else
                printf("Cid[%d] = %d\n",i,Cid[i]);
        }

        nL[setn] = npL;
        nN[setn] = npN;
        nM[setn] = npM;
        nC[setn] = npC;
        nV[setn] = np - (npL+npN+npM+npC);
        
        sprintf(datfilename,"data%d/energy.dat",setn+data_ini_id);
        fp10 = fopen(datfilename,"r");
        fscanf(fp10, "%d %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
        E[setn] = tmp2;
        printf("E[%d] = %f\n",setn,E[setn]);
        fclose(fp10);
    }
    
    // identify the energy at end states
    double E0L, E0N, E0M, E0C, E0V;
    E0L = 0.0;
    E0N = 0.0;
    E0M = 0.0;
    E0C = 0.0;
    E0V = 0.0;
    for (i=0;i<ndata;i++) {
        if (nL[i] == np && E[i] < E0L)
            E0L = E[i];
        if (nN[i] == np && E[i] < E0N)
            E0N = E[i];
        if (nM[i] == np && E[i] < E0M)
            E0M = E[i];
        if (nC[i] == np && E[i] < E0C)
            E0C = E[i];
        if (nV[i] == np && E[i] < E0V)
            E0V = E[i];
    }
    printf("----------------------------------------\n");
    printf("E0L=%.4e, E0N=%.4e, E0M=%.4e, E0C=%.4e, E0V=%.4e\n",E0L,E0N,E0M,E0C,E0V);
    for (i=0;i<ndata;i++)
        Ef[i] = E[i] - (nL[i]*E0L + nN[i]*E0N + nM[i]*E0M + nC[i]*E0C + nV[i]*E0V)/np;

    printf("----------------------------------------\n");
    printf("Data Table\n");
    for (i=0;i<ndata;i++) {
        for (j=0;j<np;j++)
            printf("%2d ",data[i][j]);
        printf("\n");
    }
    printf("----------------------------------------\n");
    printf("Magmom Table\n");
    for (i=0;i<ndata;i++) {
        for (j=0;j<np;j++)
            printf("%2d ",magmom[i][j]);
        printf("\n");
    }
    printf("----------------------------------------\n");
    printf("E    \tEf\n");
    for (i=0;i<ndata;i++)
        printf("%.4e\t%.4e\n",E[i],Ef[i]);
    
    FILE *fp11, *fp12, *fp13, *fp14, *fp15;
    fp11 = fopen("data_orig.dat","w");
    fp12 = fopen("magmom_orig.dat","w");
    fp13= fopen("E_orig.dat","w");
    fp14= fopen("Ef_orig.dat","w");
    fp15= fopen("nCat_orig.dat","w");
    for (i=0;i<ndata;i++) {
        for (j=0;j<np;j++) {
            fprintf(fp11, "%d ", data[i][j]);
            fprintf(fp12, "%d ", magmom[i][j]);
        }
        fprintf(fp11, "\n");
        fprintf(fp12, "\n");
        fprintf(fp13, "%.6f\n", E[i]);
        fprintf(fp14, "%.6f\n", Ef[i]);
        fprintf(fp15, "%d %d %d %d\n",nL[i],nN[i],nM[i],nC[i]);
    }
    fclose(fp11);
    fclose(fp12);
    fclose(fp13);
    fclose(fp14);
    fclose(fp15);
}


















