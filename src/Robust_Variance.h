#ifndef ROBUST_VARIANCE_H
#define ROBUST_VARIANCE_H

#include <fstream>
#include "Model.h"
#include "multic.h"
#include "Composite.h"

class Robust_Variance{
  int env_vcnum;
  int iloci,itraits,icovs, iinitvcnum, imubeta_dim, ismtcpq_dim, itraitcovnum, imubetacovnum;
  int  nfam;    // number of family
  int* N;    // member of each family

  Trait*     trait_array;
  Marker*    marker_array;
  Cov*       cov_array;
  // Composite composite;
  // double** D;
  double  ybar;
  double  **y;
  int *missing_flag;
  int header;

  Family fam;
  ifstream      fp_para;
  ifstream      fp_data;
  ifstream      fp_loci;
  ifstream      fp_share;    

 public:
  Robust_Variance();
  ~Robust_Variance();
  void main_fun();
  //void fix_effect();
  void buildModel(Composite* composite, ofstream* outfp,
		  TraitMarkerCov_par *tmc, FortData *fortArray,
		  int dataSize);
  //void add(double**,double**,int,int,double**,int,int,double*,int);
  void getY(ShrinkData*,double*,int);
  void getybar(ReadFamilyData*);
  //assume files are opened already
  //void OpenFiles();
  //void CloseFiles();
};

#endif
