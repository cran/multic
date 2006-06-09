#ifndef INIT_VALUE_PAR_H
#define INIT_VALUE_PAR_H

#include "multic.h"
#include "TraitMarkerCov_par.h"
#include <fstream>

class InitValue_par {
 private:
  // from TraitMarkerCov_par.h
  int itraits,iloci,icovs,iinitvcnum,env_vcnum, irepeatmst;
  // output
  int init_mu_flag, init_s_flag,  init_m1_flag, init_m2_flag,
    init_t_flag, init_c_flag,  init_p_flag,  init_q_flag;
  double init_x_mu[MAXNUMTRAIT], init_x_S[INITVCNUM], 
    init_x_M1[INITVCNUM], init_x_M2[INITVCNUM],
    init_x_T[INITVCNUM],    init_x_beta[MAXNUMTRAIT][MAXNUMCOV],
    init_x_c[INITVCNUM],    init_x_p[INITVCNUM], 
    init_x_q[INITVCNUM];
  int N1;
  int convertPrecision(char *);
 public:
  InitValue_par(TraitMarkerCov_par *tmc, double *initalValues,
		char **precisions, int maxIterations, double *initialBetas);
  ~InitValue_par();
  void GetInitValues(double *initalValues, char **precision);
  void GetXValues(int, double[], char[]);
  void GetInitCovValues(char[]);
  int getmuflag();
  int getsflag();
  int getm1flag();
  int getm2flag();
  int gettflag();
  int getcflag();
  int getpflag();
  int getqflag();
  int getiterationnum();
  double* getinitmu();
  double* getinitS();
  double* getinitM1();
  double* getinitM2();
  double* getinitT();
  double* getinitc();
  double* getinitp();
  double* getinitq();
  void getinitbeta(double beta[MAXNUMTRAIT][MAXNUMCOV]);
};

#endif
