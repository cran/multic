/******************************************************************************
Title: loadVarSandwich
Description: loadVarSandwich is responsible for reading the file 'varCovar.log'.
             varCovar.log is an output file from multic and represents the
             variance covariance matrix calculated in
             Robust_Variance::buildModel.  This function reads that information
             and "returns" the information back to Splus.
Input: double *varCovarMatrix - varCovarMatrix is an array allocated in Splus
                                which will be filled with the variance 
                                covariance matrix from the null hypothesis and
                                the variance covariance matrices from the
                                alternative hypotheses from the last run of
                                multic.
       Sint *randomEffectsCount - randomEffectsCount is a singular value that
                                  indicates how many polygenic and major gene
                                  values there are.
       Sint *altHypCount - altHypCount is a singular value that indicates the
                           number of alternative hypotheses were calculated
                           during the last run of multic.  Keep in mind that
                           this value does not account for the null hypothesis.
Ouput: NONE
Side Effects: If 'varCovar.log' can not be opend for reading, an appropriate
              error message will be printed to standard error and the program
              will terminate.
Author: Eric Lunde, 10-24-03
Updates: (Date, Modified By, Modification Description)
03-31-2005, Eric Lunde, The number of environmental variance (covariance)
            values no longer is dependent upon whether multic is run in
            multivariate or longitudinal mode (it now always operates as if
            in multivariate mode).  Therefore I've removed the
            environmentCount variable and replaced it with the
            randomEffectsCount which now always is the correct value.
******************************************************************************/
#include <fstream>
#include <iostream>
#include "verS.h"
using namespace std;

extern "C" {

void loadVarSandwich(double *varCovarMatrix,
		     Sint *randomEffectsCount,
		     Sint *altHypCount,
		     Sint *nVC,
		     Sint *isPolyFixed,
		     Sint *isMg1Fixed) {
  // Open "varCovar.log" and test for failure
  ifstream varCovar("varCovar.log");
  if(varCovar.fail()) {
    PROBLEM "The file varCovar.log could not be opened for reading.\nloadVarSandwich.cpp key 53\n"
      RECOVER(NULL_ENTRY);
  }

  // Although the data is gouped in two dimensional pieces, our data structure
  // is linear, nextIndex will keep track of our position in that structure
  int nextIndex = 0;

  if(*isMg1Fixed) {
    // Extract from varCovar.log the variance covariance matrix associated
    // with all of the hypotheses.
    int length = (int) (*nVC * *randomEffectsCount)
      * (*nVC * *randomEffectsCount) * (*altHypCount + 1);
    for(int i=0; i<length; i++) {
      varCovar >> varCovarMatrix[nextIndex++];
    }
  } else {
    // Fill the upper third of the variance covariance matrix under the null
    // hypothesis
    if(!*isPolyFixed) {
      for(int i=0; i<*randomEffectsCount; i++) {
	for(int j=0; j<*randomEffectsCount; j++) {
	  varCovar >> varCovarMatrix[nextIndex++];
	}
	for(int j=0; j<*randomEffectsCount; j++) {
	  varCovarMatrix[nextIndex++] = -9.0;
	}
	// *vVC - 2 because if we are here, then both isMg1Fixed and
	// isPolyFixed are true (equal to one each)
	for(int j=0; j<*nVC - 2; j++) {
	  varCovar >> varCovarMatrix[nextIndex++];
	}
      }
    }

    // Fill the middle third of null matrix
    for(int i=0; i<*randomEffectsCount; i++) {
      for(int j=0; j<*nVC; j++) {
	varCovarMatrix[nextIndex++] = -9.0;
      }
    }

    // Fill bottom third of null matrix
    for(int i=0; i<*nVC - !*isPolyFixed - 1; i++) {
      if(!*isPolyFixed) {
	for(int j=0; j<*randomEffectsCount; j++) {
	  varCovar >> varCovarMatrix[nextIndex++];
	}
      }
      for(int j=0; j<*randomEffectsCount; j++) {
	varCovarMatrix[nextIndex++] = -9.0;
      }
      // We know we have the -1 because isMg1Fixed is True (1)
      for(int j=0; j<*nVC - !*isPolyFixed - 1; j++) {
	varCovar >> varCovarMatrix[nextIndex++];
      }
    }

    // Extract from varCovar.log the variance covariance matrix associated with
    // all of the alternative hypotheses.
    int altLength = (int) (*nVC * *randomEffectsCount)
      * (*nVC * *randomEffectsCount) * (*altHypCount);
    for(int i=0; i<altLength; i++) {
      varCovar >> varCovarMatrix[nextIndex++];
    }
  }

  // Close varCovar.log
  varCovar.close();
}

}
