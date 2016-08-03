/***********************************************/
/* WORM ALGORITHM FOR BOSONS AT FINITE DENSITY */
/* BASED ON CHARACTER EXPANSION OF Z ***********/
/* Character_worm_v2.c *************************/
/* AUTHOR: MICHAEL G. ENDRES *******************/
/* LAST MODIFIED: 09/14/06 *********************/
/***********************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "bessi.h"
#include "ran2.h"
#include "ran3.h"
#include "irbit2.h"
#include "mod.h"
#include "parameters.h"

#define RNG ran3
#define IRNG irbit2 

/*************/
/* FUNCTIONS */
/*************/
int floatToInt(float x, int m) // GENERATE AN INTEGER ON [0,m) GIVEN A FLOAT ON (0,1)
{
  return (int)floor(x*m);
};

void vecEq(int v1[3], int v2[3]) // ASSIGN v1=v2
{
  int i;
  for (i=0;i<3;i++) {
    v1[i]=v2[i];
  };
};

int vecComp(int v1[3], int v2[3]) // COMPUTE v1==v2
{
  int i;
  int ans=1;
  for (i=0;i<3;i++) {
    ans=ans&&(v1[i]==v2[i]);
  };
  return ans;
};

void printfVec(int v[3]) // PRINT VECTOR
{
  int i;
  for (i=0;i<3;i++) {
    printf("%d",v[i]);
    printf(" ");
  };
  printf("\n");
};

/***********************/
/* PROBABILITY WEIGHTS */
/***********************/
float wS(int j, int k, int n, float s[NT][NS][NS], int l[NT][NS][NS][3], float ds)
{
  float ans=0.0;
  float sTrial;

  sTrial=s[j][k][n]+ds;
  if (sTrial<=0.0) {
    return ans;
  } else {
    ans += sTrial/s[j][k][n];
    ans *= exp(-B*(2.0+1.0/(B*B)+MSq/2.0)*(sTrial*sTrial-s[j][k][n]*s[j][k][n]));
    ans *= exp(-B*(LAMBDA/16.0)*(sTrial*sTrial*sTrial*sTrial-s[j][k][n]*s[j][k][n]*s[j][k][n]*s[j][k][n]));
    ans *= bessi(l[j][k][n][0],sTrial*s[mod(j+1,NT)][k][n]/B);
    ans /= bessi(l[j][k][n][0],s[j][k][n]*s[mod(j+1,NT)][k][n]/B);
    ans *= bessi(l[j][k][n][1],sTrial*s[j][mod(k+1,NS)][n]*B);
    ans /= bessi(l[j][k][n][1],s[j][k][n]*s[j][mod(k+1,NS)][n]*B);
    ans *= bessi(l[j][k][n][2],sTrial*s[j][k][mod(n+1,NS)]*B);
    ans /= bessi(l[j][k][n][2],s[j][k][n]*s[j][k][mod(n+1,NS)]*B);
    ans *= bessi(l[mod(j-1,NT)][k][n][0],sTrial*s[mod(j-1,NT)][k][n]/B);
    ans /= bessi(l[mod(j-1,NT)][k][n][0],s[j][k][n]*s[mod(j-1,NT)][k][n]/B);
    ans *= bessi(l[j][mod(k-1,NS)][n][1],sTrial*s[j][mod(k-1,NS)][n]*B);
    ans /= bessi(l[j][mod(k-1,NS)][n][1],s[j][k][n]*s[j][mod(k-1,NS)][n]*B);
    ans *= bessi(l[j][k][mod(n-1,NS)][2],sTrial*s[j][k][mod(n-1,NS)]*B);
    ans /= bessi(l[j][k][mod(n-1,NS)][2],s[j][k][n]*s[j][k][mod(n-1,NS)]*B);
  }
  return ans;
};

float wBulkS(int j, int k, int n, float s[NT][NS][NS], int l[NT][NS][NS][3], float ds)
{
  float ans=0.0;
  float sTrial;

  sTrial=s[j][k][n]+ds;
  if (sTrial<=0.0) {
    return ans;
  } else {
    ans += sTrial/s[j][k][n];
    ans *= exp(-B*(2.0+1.0/(B*B)+MSq/2.0)*(sTrial*sTrial-s[j][k][n]*s[j][k][n]));
    ans *= exp(-B*(LAMBDA/16.0)*(sTrial*sTrial*sTrial*sTrial-s[j][k][n]*s[j][k][n]*s[j][k][n]*s[j][k][n]));
    ans *= bessi(l[j][k][n][0],sTrial*s[j+1][k][n]/B);
    ans /= bessi(l[j][k][n][0],s[j][k][n]*s[j+1][k][n]/B);
    ans *= bessi(l[j][k][n][1],sTrial*s[j][k+1][n]*B);
    ans /= bessi(l[j][k][n][1],s[j][k][n]*s[j][k+1][n]*B);
    ans *= bessi(l[j][k][n][2],sTrial*s[j][k][n+1]*B);
    ans /= bessi(l[j][k][n][2],s[j][k][n]*s[j][k][n+1]*B);
    ans *= bessi(l[j-1][k][n][0],sTrial*s[j-1][k][n]/B);
    ans /= bessi(l[j-1][k][n][0],s[j][k][n]*s[j-1][k][n]/B);
    ans *= bessi(l[j][k-1][n][1],sTrial*s[j][k-1][n]*B);
    ans /= bessi(l[j][k-1][n][1],s[j][k][n]*s[j][k-1][n]*B);
    ans *= bessi(l[j][k][n-1][2],sTrial*s[j][k][n-1]*B);
    ans /= bessi(l[j][k][n-1][2],s[j][k][n]*s[j][k][n-1]*B);
  }
  return ans;
};

/***************/
/* OBSERVABLES */
/***************/
float s2(float s[NT][NS][NS]) {
  float ans=0.0;
  int j;
  int k;
  int n;

  for (j=0;j<NT;j++) {
    for (k=0;k<NS;k++) {
      for (n=0;n<NS;n++) {
        ans += s[j][k][n]*s[j][k][n];
      };
    };
  };
  ans /= NT*NS*NS;
  return ans;
};

int Q(int t, int l[NT][NS][NS][3]) {
  int ans=0;
  int k;
  int n;

  for (k=0;k<NS;k++) {
    for (n=0;n<NS;n++) {
      ans += l[t][k][n][0];
    };
  };
  return ans;
};

/*
 *  Simulation
 */

int main(void)
{
  int L[3]={NT,NS,NS};
  int i;  /*  MC step counter  */
  int j;  /*  Site counter for t direction */
  int k;  /*  Site counter for x direction  */
  int n;  /*  Site counter for y direction  */
  long seed;  /*  Random number generator seed  */
  long iseed;  /*  Random bit generator seed  */
  float s[NT][NS][NS];  /*  Site variables  */
  int l[NT][NS][NS][3];  /*  Link variables  */

  int head[3];
  int tail[3];
  int trialHead[3];
  int stepBv;
  int stepDir;
  float fug[3]={exp(MU),1.0,1.0};
  float invFug[3]={exp(-MU),1.0,1.0};

  float ds;  /*  Trial change in site variable  */
  int rAcc;  /*  Acceptance rate counter  */
  FILE *binaryp;

  /********************/
  /* INITIALIZE SEEDS */
  /********************/
  seed = -time((time_t *)NULL);  /* Random generator */
  iseed = -time((time_t *)NULL);  /* Bit generator */

 /******************************/
  /* INITIALIZE STATE OF SYSTEM */
  /******************************/
  binaryp=fopen("system.conf","rb");
  if (binaryp==NULL) {
    /******************/
    /* USE COLD START */
    /******************/
    for (j=0;j<NT;j++) {
      for (k=0;k<NS;k++) {
        for (n=0;n<NS;n++) {
          s[j][k][n]=SINIT;
          l[j][k][n][0]=0;
          l[j][k][n][1]=0;
          l[j][k][n][2]=0;
        };
      };
    };
  }else{
    /******************************/
    /* USE PREVIOUS CONFIGURATION */
    /******************************/
    fread(s, sizeof (float), NT*NS*NS, binaryp);
    fread(l, sizeof (int), NT*NS*NS*3, binaryp);
    fclose(binaryp);
  };

/***  BEGIN MC SIMULATIONS  ***/
  for (i=0;i<NIT;i++) {
    rAcc=0;  /* Accceptance rate counter  */
/***  UPDATE SITE VARIABLES  ***/
    for (j=0;j<NT;j++) {  /*  Sum over t  */
      for (k=0;k<NS;k++) {  /*  Sum over x  */
        for (n=0;n<NS;n++) {  /*  Sum over y  */
          if (j==0||k==0||n==0||j==(NT-1)||k==(NS-1)||n==(NS-1)) {
            /***  IF ON THE ENDGE OF THE BOX, DO THIS  ***/
            ds=D*(2.0*RNG(&seed)-1.0);  /*  Generate a trial step  */
            if (wS(j,k,n,s,l,ds)>RNG(&seed)) {  /* Update site variable  */
              s[j][k][n] += ds;
              rAcc++;
            };
          } else {
            /***  OTHERWISE, DO THIS  ***/
            ds=D*(2.0*RNG(&seed)-1.0);  /*  Generate a trial step  */
            if (wBulkS(j,k,n,s,l,ds)>RNG(&seed)) {  /* Update site variable  */
              s[j][k][n] += ds;
              rAcc++;
            };
          };
        };  /* End of n for loop  */
      };  /*  End of k for loop  */
    };  /*  End of j for loop  */

/****************************/
/*** A WORMY ALTERNATIVE  ***/
/****************************/
    for (j=0;j<NIT_WORM;j++) {
      for (k=0;k<3;k++) { 
        head[k]=floatToInt(RNG(&seed),L[k]);
      };
      vecEq(tail,head);
      do {
        stepBv=floatToInt(RNG(&seed),3);
        stepDir=2*IRNG(&iseed)-1;
        vecEq(trialHead,head);
        trialHead[stepBv]=mod(trialHead[stepBv]+stepDir,L[stepBv]);
        if (stepDir>0) {
          if (fug[stepBv]*bessi(l[head[0]][head[1]][head[2]][stepBv]+1,s[head[0]][head[1]][head[2]]*s[trialHead[0]][trialHead[1]][trialHead[2]])\
          /bessi(l[head[0]][head[1]][head[2]][stepBv],s[head[0]][head[1]][head[2]]*s[trialHead[0]][trialHead[1]][trialHead[2]])>RNG(&seed)) {
            l[head[0]][head[1]][head[2]][stepBv]+=1;
            vecEq(head,trialHead);
          };
        } else {
          if (invFug[stepBv]*bessi(l[trialHead[0]][trialHead[1]][trialHead[2]][stepBv]-1,s[head[0]][head[1]][head[2]]*s[trialHead[0]][trialHead[1]][trialHead[2]])\
          /bessi(l[trialHead[0]][trialHead[1]][trialHead[2]][stepBv],s[head[0]][head[1]][head[2]]*s[trialHead[0]][trialHead[1]][trialHead[2]])>RNG(&seed)) {
            l[trialHead[0]][trialHead[1]][trialHead[2]][stepBv]-=1;
            vecEq(head,trialHead);
          };
        };
      } while (vecComp(head,tail)!=1);
    };
/***  PRINT AVERAGES  ***/
    if ((i>NEQ)&&(i%TCORR==0)) {
      printf("%d\n",Q(0,l));  /*  Print total charge corresponding to field configuration  */
//      printf("%f\n",s2(s));  /*  Print average |\phi|^2  */
//      printf("%f\n",rAcc*1.0/(NT*NS*NS));  /*  Print acceptance rate  */
    };
  };  /* End of i for loop */

  /******************************/
  /* WRITE SYSTEM STATE To DISK */
  /******************************/
  binaryp=fopen("system.conf","wb");
  fwrite(s, sizeof (float), NT*NS*NS, binaryp);
  fwrite(l, sizeof (int), NT*NS*NS*3, binaryp);
  fclose(binaryp);
  /************************/
  /* END OF MC SIMULATION */
  /************************/
  return 0;
}
