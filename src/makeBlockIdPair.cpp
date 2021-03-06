#include <iostream>
#include <cstring>
#include "multicString.h"

#include "verS.h"
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif
#include "Rostream.h"
#include "Rstreambuf.h"

using namespace Rcpp;

extern "C" {

/****************************************************************************
 * uniqueIds is a character vector of the famid-id pairs from a pedigree
 *           where each index i is a concatenation of person i's famid, '-',
 *           and id.
 * familySizes is an integer vector where each index i is the number of 
 *             people in family i
 ***************************************************************************/
s_object *makeBlockIdPair(s_object *uniqueIds, s_object *familySizes) {
  S_EVALUATOR

  // cumFamilySize will offset us to the correct location in the uniqueIds
  // array
  Sint cumFamilySize = 0;
  Sint currentFamilySize = 0;


  Rcout << "Calculating idPairLength" << std::endl;
  // Compute how many id-pairs there will be.
  Sint idPairLength = 0;
  for(Sint i = 0; i < LENGTH(familySizes); i++) {
    currentFamilySize = INTEGER_POINTER(familySizes)[i];
    idPairLength += (currentFamilySize - 1) * currentFamilySize / 2;
    //cout << "currentFamilySize = " << currentFamilySize << std::endl;
    //cout << "idPairLength = " << idPairLength << std::endl;
  }
  Sint idPairIndex = 0;

  //cout << "After loop to compute # of id pairs" << std::endl;

  s_object *idPairs = NEW_CHARACTER(idPairLength * 2);

  //cout << "Before last big loop" << std::endl;


  // PROBLEM WITH sw2mloci HAPPENS IN THIS LOOP
  // i iterates over the families
  for(Sint i = 0; i < LENGTH(familySizes); i++) {
    currentFamilySize = INTEGER_POINTER(familySizes)[i];

    //cout << "\ni = " << i << std::endl;
    //cout << "currentFamilySize = " << currentFamilySize << std::endl;

    // j iterates over the first familySizes[i] - 1 persons of each family
    // (This will be the first of an id pair and the last person in the
    // family will not be the first of the pair)
    for(Sint j = 0; j < currentFamilySize - 1; j++) {
      // k iterates over the last familySizes[i] - 1 persons of each family
      // (This will be the second of an id pair and the first person in the
      // family will not be the second of the pair)
      for(Sint k = j + 1; k < currentFamilySize; k++) {
	// The first assignment sets the first of the id-pair,
	// the second assignment sets the second.
	
	multic_SET_STRING_ELT(idPairs, idPairIndex,
			      multic_STRING_ELT(uniqueIds, 
						cumFamilySize + j));
	multic_SET_STRING_ELT(idPairs, idPairIndex + idPairLength,
			      multic_STRING_ELT(uniqueIds,
						cumFamilySize + k));

	idPairIndex++;
      }
    }
    cumFamilySize += currentFamilySize;
  }

  //cout << "After last big loop" << std::endl;

  s_object *dim = NEW_INTEGER(2);
  INTEGER_POINTER(dim)[0] = idPairLength;
  INTEGER_POINTER(dim)[1] = 2;

  SET_DIM(idPairs, dim);

  Rcout << "Returning from makeBlockIdPair" << std::endl;

  return idPairs;
}

}
