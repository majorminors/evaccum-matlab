#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <matrix.h>
#include <time.h>
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


void LBA(int N,double *B,double *C0,double *Ame,double *Astd,double T0,double rt[],double choice[],double act_noavail[], double *cond_ratio){
    int i;
    int iter;
    double base,tag;
    double *A=(double *)malloc(sizeof(double)*N);
    double *C=(double *)malloc(sizeof(double)*N);
    double *rt_indv=(double *)malloc(sizeof(double)*N);
    
    for(iter=0;iter<NUMSIM;iter++){
        
        base=cond_ratio[0];
        tag=AverageRandom();
        for(i=0;i<N;i++){
            if (tag>=base){
                base+=cond_ratio[i+1];
                continue;
            }
            else{
                act_noavail[iter]=i;
                break;
            }
        }

        /*% sample drift rate
        for(i=0;i<N;i++){
            A[i]=-1;
            C[i]=-1;
            while (A[i]<0){
                A[i]=noise()*Astd[i]+Ame[i];
            }
            while (B[i]<C[i]){
                C[i]=AverageRandom()*C0[i]; 
            }
            rt_indv[i]=(B[i]-C[i])/A[i];
        }*/
        
        /*% sample drift rate*/
        
        for(i=0;i<N;i++){
            /* get drift rate*/
            A[i]=noise()*Astd[i]+Ame[i];
            /* get starting point*/
            C[i]=B[i]+1;
            while (B[i]<C[i]){
                C[i]=AverageRandom()*C0[i]; /*sample starting point*/
            }
            
            /* the decision time to reach boundary for each finger*/
            rt_indv[i]=(B[i]-C[i])/A[i];
        }
        
        rt[iter]=MAXT;
        choice[iter]=0;
        
        for(i=0;i<N;i++){
            /* loop through all alternatives to find the choice that has the minimum decision time */
            if(rt_indv[i]<rt[iter] && rt_indv[i]>0){                
                /* model predicted decision time and response*/
                rt[iter]=rt_indv[i];
                choice[iter]=(double)i+1;
            }
        }
        /* model predicted response time (decision time + non-decision time*/
        rt[iter]+=T0;
    }
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    double *rt;
    double *choice;
    int N;
    double *B;
    double *C0;
    double *Ame;
    double *Astd;
    double T0;
    double *cond_ratio; 
    double *act_noavail;
    
    N =*mxGetPr (prhs[0]);	
    B = mxGetPr (prhs[1]);	
    C0 = mxGetPr (prhs[2]);	
    Ame = mxGetPr (prhs[3]);	
    Astd = mxGetPr (prhs[4]);	
    T0 = *mxGetPr (prhs[5]);	
    cond_ratio = mxGetPr (prhs[6]);
    
    plhs[0] = mxCreateDoubleMatrix (1,NUMSIM, mxREAL);
    plhs[1] = mxCreateDoubleMatrix (1,NUMSIM, mxREAL);
    plhs[2] = mxCreateDoubleMatrix (1,NUMSIM, mxREAL);
    rt = mxGetPr (plhs[0]);	
    choice = mxGetPr (plhs[1]);
    act_noavail = mxGetPr (plhs[2]);
    
    /*srand(time(NULL));*/
    srand(14);
    LBA(N,B,C0,Ame,Astd,T0,rt,choice,act_noavail,cond_ratio);
}

