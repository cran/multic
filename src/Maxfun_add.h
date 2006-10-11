#ifndef MAX_FUN_ADD_H
#define MAX_FUN_ADD_H

#include "verS.h"
//#include "Maxfun_source.h"
#include "TraitMarkerCov_par.h"
#include "InitValue_par.h"
#include "multic.h"

class Maxfun_add {
  double theta[45];
  int lfl;  
  double f; 
  int nfe; 
 public:
  void run(TraitMarkerCov_par*, InitValue_par*,
	   FortData *, int, ShareRelation *, int);
};

#endif
