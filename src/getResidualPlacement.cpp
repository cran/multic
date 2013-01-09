#include <iostream>

#include "verS.h"
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif
using namespace std;

extern "C" {

/*********************
Title: getResidualPlacement
Description: allocates and populates a vector by repeating hasAllValues per
             family lociCount times
Input: s_object *familySizes - integer vector with length equal to the number
                               of families in the study, where each index is
                               the size of a family 
       s_object *hasAllValues - logical vector with length equal to the
                                number of people in the study, where each 
                                index is T if and that person has no missing
                                values for any trait or covariate
       s_object *lociCount - integer value representing the number of markers
                             analyzed during the multic calculations plus one
                             for the polygene analysis
Output: s_object *residualPlacement - logical vector with length lociCount *
                                      total people count where each family's
                                      hasAllValues' data is repeated
                                      lociCount times
Side Effects: NONE
Author: Eric Lunde, 12/17/2004
*********************/
s_object *getResidualPlacement(s_object *familySizes, s_object *hasAllValues,
			       s_object *lociCount) {
  // I'm not sure what S_EVALUATOR does, but I'm guessing we don't need it now
  // that we're dropping support for S-PLUS (PV 1/7/2009)
  //S_EVALUATOR
 
  Sint *familySizesPtr = INTEGER_POINTER(familySizes);
  Sint *hasAllValuesPtr = LOGICAL_POINTER(hasAllValues);
  Sint lociCountValue = INTEGER_VALUE(lociCount);

  Sint familyCount = LENGTH(familySizes);
  Sint hasAllValuesOffset = 0;

  Sint peopleCount = 0;
  for(int i = 0; i < familyCount; i++) {
    peopleCount += familySizesPtr[i];
  }

  s_object *residualPlacement = NEW_LOGICAL(peopleCount * lociCountValue);
  Sint *residualPlacementPtr = LOGICAL_POINTER(residualPlacement);
  int residualPlacementIndex = 0;

  // For each family, repeat their hasAllValues data lociCount times.
  for(int i = 0; i < familyCount; i++) {
    for(int j = 0; j < lociCountValue; j++) {
      for(int k = 0; k < familySizesPtr[i]; k++) {
	residualPlacementPtr[residualPlacementIndex++] =
	  hasAllValuesPtr[hasAllValuesOffset + k];
      }
    }
    hasAllValuesOffset += familySizesPtr[i];
  }

  return residualPlacement;
}

}
