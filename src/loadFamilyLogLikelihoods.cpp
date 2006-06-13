/******************************************************************************
Title: loadFamilyLogLikelihoods
Description: loadFamilyLogLikelihoods is designed to be called from Splus in
             the multic program.  It reads the file "fam.lik" and saves the
             log likelihood information per family per locus in
             Splus-managed memory.
Input: s_object *familyCount - integer value specifying the number of number
                               of families fam.lik contains information for.
       s_object *hypothesisCount - hypothesisCount is an integer value
                                   specifying the number of hypotheses
                                   (null and all alternative) fam.lik will
                                    contain.
Output: s_object *familyLogLiks - numeric array of dimensions familyCount by
                                  2 by hypothesisCount.  
Side Effects: This function reads from "fam.lik".  This is not a user-defined
              name.  It is hardcoded.  It must be in the same directory as
              the Splus session is occuring.  If it does not exist, an empty
              numeric vector will be returned and an error message will be
              printed.
Author: Eric Lunde, 10-17-03
Updates: (Date, Modifyied By, Modification Description)
         04-08-2004, Eric Lunde, char **lociNames, long *familyIds, and double
         *logLiks were all removed as input parameters.  The desired output of
        this function changed and since I had more knowledge about this type
        of Splus/C++ object modification, I made the switch to build the
        entire object in C++ rather than some here and the rest in Splus. 
        Also, alternativeCount was replaced with hypothesisCount and the
        expected value increased by one to account for the null hypothesis.
******************************************************************************/
#include <iostream>
#include <fstream>
#include <cstring>
#include "multicString.h"

#include <S.h>
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif
using namespace std;

extern "C" {

s_object *loadFamilyLogLikelihoods(s_object *familyCount,
				    s_object *hypothesisCount) {
  long familyCountValue = INTEGER_VALUE(familyCount);
  long hypothesisCountValue = INTEGER_VALUE(hypothesisCount);
  char tempString[1024];
  int logLikIndex = 0;
  double zero = 0.0;
  char file[] = "fam.lik";

  ifstream logLikelihoodFile(file);
  if(logLikelihoodFile.fail()) {
    cerr << "The file " << file << " could not be opened for reading" << endl
	 << "loadFamilyLogLikelihoods.cpp key 51" << endl;
    return NEW_NUMERIC(0);
  }

  // Read the family ids and store them in s_object *familyIds
  s_object *familyIds = NEW_CHARACTER(familyCountValue);
  for(int i = 0; i < familyCountValue; i++) {
    logLikelihoodFile >> tempString;
    multic_SET_STRING_ELT(familyIds, i, tempString);
  }

  // Read the rest of the data, each line beginning with the ibd name or
  // null (to be stored in lociNames)  followed by the family log
  // likelihoods for each family (to be stored in familyLogliks).
  s_object *lociNames = NEW_CHARACTER(hypothesisCountValue);
  
  s_object *familyLogLiks = NEW_NUMERIC(familyCountValue
					* 2
					* hypothesisCountValue);
  double *familyLogLiksPtr = NUMERIC_POINTER(familyLogLiks);

  for(int i = 0; i < hypothesisCountValue; i++) {
    // Read the ibd name
    logLikelihoodFile >> tempString;
    multic_SET_STRING_ELT(lociNames, i, tempString);
    
    // Read the family log likelihood data
    for(int j = 0; j < familyCountValue; j++) {
      logLikIndex = i * familyCountValue * 2 + j;
      // Assign the log likelihood
      logLikelihoodFile >> familyLogLiksPtr[logLikIndex];
      if(i == 0) {
	// The null hypothesis will not have a lod score so generate a
	// NaN which will be translated to NA by Splus
	familyLogLiksPtr[logLikIndex + familyCountValue] = zero / zero;
      } else {
	// Calculate the lod score
	familyLogLiksPtr[logLikIndex + familyCountValue]
	  = -2
	  * (familyLogLiksPtr[j] - familyLogLiksPtr[logLikIndex])
	  / 4.6;
      }
    }
  }

  logLikelihoodFile.close();

  // Generate the 3rd dimension names
  s_object *likelihoodNames = NEW_CHARACTER(2);
  multic_SET_STRING_ELT(likelihoodNames, 0, "log.lik");
  multic_SET_STRING_ELT(likelihoodNames, 1, "lod.score");
  
  // Set the dimensions and dimnames for the returned object.
  int dimensionCount = 3;
  s_object *dim = NEW_INTEGER(dimensionCount);
#ifdef USING_R
  int *
#else
  long *
#endif
    dimPtr = INTEGER_POINTER(dim);
  dimPtr[0] = familyCountValue;
  dimPtr[1] = 2;
  dimPtr[2] = hypothesisCountValue;
  SET_DIM(familyLogLiks, dim);

  s_object *dimnames = NEW_LIST(dimensionCount);
  SET_ELEMENT(dimnames, 0, familyIds);
  SET_ELEMENT(dimnames, 1, likelihoodNames);
  SET_ELEMENT(dimnames, 2, lociNames);
  SET_DIMNAMES(familyLogLiks, dimnames);

  return (familyLogLiks);
}

}
