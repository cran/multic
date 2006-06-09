#ifndef MAX_FUN_FUN_H
#define MAX_FUN_FUN_H

#include <stdio.h>
#include "verS.h"

#define MAXFUN_NP 45
#define MAXFUN_NPV MAXFUN_NP
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))

typedef struct{  
  FILE *iout, *idet;
  double thin[MAXFUN_NP],thl[MAXFUN_NP],thu[MAXFUN_NP],
    stpin[MAXFUN_NP],epsd,yota,epst,epsc1,epsc2,epsc3;
  int istin[MAXFUN_NP],nt,maxit,method,ixvc,ihit,lprt;
  char label[MAXFUN_NP][29];
} Maxfun_exp;


class Maxfun_fun{
 public:
  virtual void eval_fun(double *theta,double *f,int *nfe,int *lex)=0;
  virtual void dep_fun(double *tr, int *lex)=0;
  virtual Maxfun_exp getExp()=0;
};

#endif
