#ifndef LEAST_H
#define LEAST_H

/* least1.c is designed for multi-traits variance component analysis and 
using the overall mean of each trait as the mean of the trait values of all 
fathers, mothers and sibs.  Also, least1.c includes least squares estimation 
for boundary */

#include "Estimator.h"
#include "ReadFamilyData.h"
#include "multic.h"
#include <fstream>

class Least : public Estimator{
  ReadFamilyData *readfam;
  ShareRelation *shareArray;
  int relationSize;
  int relationIndex;

  ofstream outfp;
  ifstream pt_loci;
  int familyCount;
  int* familySizes;
  int m; //total number of individual
  int itraits;
  int* t_start;
  int* t_len;
  //	double  y[MAXNUMTRAIT][MAXNUMFAM][FAMMEMNUM];
  //	double	a[MAXNUMFAM][FAMMEMNUM][FAMMEMNUM]; 
  double  ***y, ***a, ***r_mat;
  double	ybar[MAXNUMTRAIT], var_g[MAXNUMTRAIT][MAXNUMTRAIT], 
    var_G[MAXNUMTRAIT][MAXNUMTRAIT], var_e[MAXNUMTRAIT];
  double  A11, A12, A22;
  int **missing_flag;
  int header;
  Family fam;

 public:
  Least(TraitMarkerCov_par *tmc, ShareRelation *shareArray, int relationSize,
	FortData *fortArray, int dataSize);
  virtual ~Least();
  void CloseFiles();
  void getresponse();
  void getcoefficient();
  void getestimator();
  void least_main();
  void getmajorgene(double   init_x_M1[INITVCNUM]);
  void getpolygene(double    init_x_S[INITVCNUM]);
  void getenvironment(double init_x_T[INITVCNUM]);
  double getmajorgene1();
  double getmajorgene2();
  double getpolygene();
  double getenvironment();
  double getC();
  double getP();
  double getQ();
};

#endif
