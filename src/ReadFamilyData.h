#ifndef READ_FAMILY_DATA_H
#define READ_FAMILY_DATA_H

#include "InitValue_par.h"
#include "multic.h"
#include <fstream>

#define MAXNUMFAM 10000

class ReadFamilyData {
  int familyCount;
  int* familySizes; 
  FortData *fortArray;
  int dataSize;
  int dataIndex;

  Family* fam;
  //  ifstream* fp_para;
  int env_vcnum;
  int iloci,itraits,icovs, iinitvcnum, imubeta_dim, ismtcpq_dim, itraitcovnum, imubetacovnum;
  int total_missing_trait[MAXNUMTRAIT];
  int total_missing_marker[MAXNUMMARKER];
  int total_trait_values, total_cov_values;
  int* t_start;
  int* t_len; /* For multiple traits.       */
  int* m_start;
  int* m_len;/* For multiple markers.      */
  int* c_start;
  int* c_len;   /* For multiple covariates.   */
  Trait*     trait_array;
  Marker*    marker_array;
  Cov*       cov_array;

 public:
  ReadFamilyData(TraitMarkerCov_par *tmc, int col, Family* fam1,
		 FortData *fa, int ds);
  ~ReadFamilyData();
  void getdata(int);
  int gettrait();
  int getmarker();
  int getcov();
  int gett_start();

  int getFamilyCount();
  int* getFamilySizes();
  void resetDataIndex();
  FortData *getFortArray();
};

#endif
