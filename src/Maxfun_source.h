#ifndef MAX_FUN_SOURCE_H
#define MAX_FUN_SOURCE_H

#include "Maxfun_fun.h"

class Maxfun_source {
 private:
  Maxfun_exp Exp;
  Maxfun_fun* mfun;
  double EPSILON;

 public:
  Maxfun_source(Maxfun_fun* mfun1);
  void maxfun(double *theta,double *f,int *nfe,int *lfl);

  void bsrch(double theta[],double *f,int *nfe,int *ifl);
  void nsrch(double theta[],double *f,int *nfe,double dth[],
	     int *ifl);
  void psrch(double theta[],double *f,int *nfe,double dth[],
	     int *ifl);
  void comb(int *ntotal,int *nchoos,int *nsw,int list[]);
  void nrstep(double theta[],double *f,int *nfe,int *ifl);
  void binit(double theta[],double *f,int *nfe,int *ihess);
  void bupdt(double gpr[],int *lex);
  void direct(int *lex);
  void lsrch(double theta[],double *f,int *nfe,int *ifirst,
	     int *ifl);
  void prepd(double theta[],double *f,int *nfe,int *lex);
  void fixbnd(double theta[],double *f,int *nfe,int *lex);
  void bndchk(double *th,double *thb,int *lex);
  void bcnvch(double theta[],int *lex);
  void vcnvch(double theta[],int *lex);
  void itest(double theta[],int *i,double *eloc,int *lex);
  void endit(double theta[],double *f);
  void vcmx(double theta[],double *f,int *nfe,int *ih);
  void augv(double theta[],int *ih);
  void deriv1(double theta[],double *f,int *nfe);
  void deriv2(double theta[],double *f,int *nfe,int *ih,int *lex);
  void fitder(double *d1,double *d2,double *d3,double *dd,
	      double *prmu,int *lex);
  double dfn(double *ath,double *sf);
  double dfninv(double *ath,double *delth);
  double efn(double *ath,double *eloc);
  void mxnvrt(double *a,int *mra,int *m,double *b,int *mrb,
	      int *lex);
  void mmult(double *a,int *mra,int *ma,int *namb,double *b,
	     int *mrb,int *nb,double *c,int *mrc);
  // After adding #include <S.h> to Maxfun_fun.h, the compiler said FILE was
  // ambiguous.  So the :: were addded here (and below).
  void mout(double *a,int *mra,int *m,int *n, ::FILE *ipr);
  void lbd(double theta[],double *f);
  void lbv(double theta[],double *f,int *lin);
  void litd(double theta[],double *f,int *nfe);
  void litv(double theta[],double *f,int *nfe);
  void lf(double theta[],double *f,int *nfe);
  void copyr(::FILE *pout);
};

#endif
