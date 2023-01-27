#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <matrix.h>

#define DT 0.001		
#define NUMSIM 100000
#define MAXT 20

float noise(){		
 static float invcum[100];
 static int virgin = 1;
 if (virgin){
  
  float norm[600], cum[600];
  int i;
  float target;
  int current = 0;

  virgin = 0;
  for (i=0; i<600; i++)
   norm[i] = 1.0 / sqrt (2 * 3.14159) * exp ( -pow (0.01*i - 3, 2) / 2 );
  cum[0] = norm[0] * 0.01;
  for (i=1; i<600; i++)
   cum[i] = cum[i-1] + norm[i] * 0.01;

  for (i=0; i<100; i++){
   target = 0.01 * i;
   while (cum[current] < target && current < 600)
    current ++;
   invcum[i] = 0.01*current - 3;
  }
 }
 return invcum[rand() % 100] + 0.023;
}

double AverageRandom(){
    return rand()/(double)RAND_MAX;
}


void LBA(int N,double *B,double *C0,double *Ame,double *Astd,double T0,double rt[],double choice[]){
    int i;
    int iter;
    double *A=(double *)malloc(sizeof(double)*N);
    double *C=(double *)malloc(sizeof(double)*N);
    double *rt_indv=(double *)malloc(sizeof(double)*N);
    
    for(iter=0;iter<NUMSIM;iter++){
        
        
        /*% sample drift rate*/
        for(i=0;i<N;i++){
            A[i]=noise()*Astd[i]+Ame[i];
            C[i]=AverageRandom()*C0[i]; /*sample starting point*/
            rt_indv[i]=(B[i]-C[i])/A[i];
        }
        rt[iter]=rt_indv[iter%4];
        choice[iter]=(iter%4)+1;
       /* rt[iter]=MAXT;
        choice[iter]=0;
        
        for(i=0;i<N;i++){
            if(rt_indv[i]<rt[iter] && rt_indv[i]>0){
                
                rt[iter]=rt_indv[i];
                choice[iter]=(double)i+1;
            }
        }*/
        rt[iter]+=T0;
    }
}

/* interface between Matlab and C*/
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *rt;
    double *choice;
    int N;
    double *B;
    double *C0;
    double *Ame;
    double *Astd;
    double T0;
    
    /* get parameter values from matlab*/
    N =*mxGetPr (prhs[0]);	
    B = mxGetPr (prhs[1]);	
    C0 = mxGetPr (prhs[2]);	
    Ame = mxGetPr (prhs[3]);	
    Astd = mxGetPr (prhs[4]);	
    T0 = *mxGetPr (prhs[5]);	

	
    plhs[0] = mxCreateDoubleMatrix (1,NUMSIM, mxREAL);
    plhs[1] = mxCreateDoubleMatrix (1,NUMSIM, mxREAL);
    
    rt = mxGetPr (plhs[0]);	
    choice = mxGetPr (plhs[1]);
    
    srand(14);
    /* run the model*/
    LBA(N,B,C0,Ame,Astd,T0,rt,choice);
}

