/******************************************************************************
Title: loadInvExpSecDerRandom
Description: loadInvExpSecDerRandom is responsible for reading the file
             'invExpSecDerRandom.log'.  invSecDerRandom.log is an output file
             from multic and represents the inverse of the expected second
             derivative matrices calculated in Calculate.cpp.  This function
             reads that information and "returns" the information back to
             Splus.
Input: double *invExpSecDerRandomMatrix - invExpSecDerRandomMatrix is an array
                                          allocated in Splus which will be
                                          filled with the expected second
                                          deriviative matrix from the null
                                          hypothesis from the last run of
                                          multic followed by the expected
                                          second derivative matrices from all
                                          of the alternative hypotheses from
                                          the last run of multic.
       Sint *randomEffectsCount - randomEffectsCount is a singular value that
                                  indicates how many polygenic and major gene
                                  values there are.
       Sint *altHypCount - altHypCount is a singular value that indicates the
                           number of alternative hypotheses were calculated
                           during the last run of multic.  Keep in mind that
                           this value does not account for the null hypothesis.
Ouput: NONE
Side Effects: If 'invExpSecDerRandom.log' can not be opend for reading, an
              appropriate error message will be printed to standard error and
              the program will terminate.
Author: Eric Lunde, 10-24-03
Update: (Date, Modified By, Modification Description)
03/08/2004, Eric Lunde, Since the format of the invExpSecDerRandom.log file
            has changed, so has the reading of it.  Now it has to read the 
            ibd file name that corresponds to that matrix and all the matrices
            are the same size.  If there are internal fixes or constrained 
            effects, they are represented in the file with -9's.
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

extern "C" {

void loadInvExpSecDerRandom(double *invExpSecDerRandomMatrix,
			    Sint *randomEffectsCount,
			    Sint *altHypCount) {
  // Open "invExpSecDerRandom.log" and test for failure
  std::ifstream invExpSecDerRandom("invExpSecDerRandom.log");
  if(invExpSecDerRandom.fail()) {
    PROBLEM "The file invExpSecDerRandom.log could not be opened for reading.\ninvExpSecDerRandom.cpp key 55\n"
      RECOVER(NULL_ENTRY);
  }

  // Although the data is gouped in two dimensional pieces, our data structure
  // is linear, nextIndex will keep track of our position in that structure
  int nextIndex = 0;
  char ibdFileName[256];

  // Extract from invExpSecDerRandom.log the expected second derivative
  // matrices associated with all of the alternative hypothees.
  int totalMatrices = (int) (*altHypCount + 1);
  int matrixSize = (int) (3 * *randomEffectsCount)
    * (3 * *randomEffectsCount);

  for(int i = 0; i < totalMatrices; i++) {
    invExpSecDerRandom >> ibdFileName;
    for(int j = 0; j < matrixSize; j++) {
      invExpSecDerRandom >> invExpSecDerRandomMatrix[nextIndex++];
    }
  }

  // Close invExpSecDerRandom
  invExpSecDerRandom.close();
}

}
