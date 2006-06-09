#ifndef TRAIT_MARKER_COV_PAR_H
#define TRAIT_MARKER_COV_PAR_H

#include "verS.h"
#include "ShrinkData.h"
// As per the trend in C++ programming, this next include directive
// should read '#include <fstream>'  Unforutnately that directive
// causes Splus to not find the header file.  Using '#include <fstream.h>'
// eliminates this problem
#include <fstream>
using namespace std;

class TraitMarkerCov_par{
 private:
  int itraits,iloci, icovs,iinitvcnum,env_vcnum;
  int irepeatmst, datatype;
  char **traitNames;
  double missingValue;
  int total_trait_values, total_cov_values;
  int t_start[MAXNUMTRAIT],t_len[MAXNUMTRAIT];
  int m_start[MAXNUMMARKER],m_len[MAXNUMMARKER];
  int c_start[TOTALCOV],c_len[TOTALCOV];
  Trait trait_array[MAXNUMTRAIT];
  Marker marker_array[MAXNUMMARKER];
  Cov cov_array[TOTALCOV];
  int marker_type[MAXNUMMARKER];
  ifstream *fp_para;
 public:
  TraitMarkerCov_par(Sint, Sint, Sint, char **, double/*, Sint **/);
  ~TraitMarkerCov_par();
  void GetParameters(/*Sint **/);
  void GetFromFormat(Sint *);
  void GetMarkerType(int,int*,int);
  void TraitLociCov();
  int gettraitnum();
  int getmarkernum();
  int getcovnum();
  int getinitvcnum();
  int getiterationum();
  Cov* getcov_array();
  int* get_t_start();
  int* get_t_len();
  int* get_m_start();
  int* get_m_len();
  int* get_c_start();
  int* get_c_len();
  int getirepeatmst();
  Trait* gettrait_array();
  int getinitenvvcnum();
  Marker* getmarker_array();
  int gettotal_cov_values();
  int gettotal_trait_values();
  void TraitMarkerCovParameter(int, int*);
};

#endif
