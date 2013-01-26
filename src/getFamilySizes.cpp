#include <iostream>
#include "multicString.h"

#include "verS.h"
#ifdef USING_R
#include <Rinternals.h>
#include <Rdefines.h>
#endif
#include "Rostream.h"
#include "Rstreambuf.h"

using namespace Rcpp;

extern "C" {

/****************************************************************************
Title: getFamilySizes
Description: Determines the family sizes and returns them in the same order
             given.  This function fills a hole that S-PLUS's functions could
             not fill easily.  table would rearrage the order of the families
             and would join families with the same famid into one logical
             family.
Input: s_object *famids - character vector specifying an individual's family
                          identifier
       s_object *ids - character vector specifying an individual's personal
                       identifier (must be unique within the family)
Output: s_object *familySizes - integer vector specifying the family sizes
                                that is consistant with the order that famids
                                was given.
Side Effects: NONE
Author: Eric Lunde, 2005-08-30
****************************************************************************/
s_object *getFamilySizes(s_object *famids, s_object *ids) {
  /*
    If we have found a new family, they will either have a different famid
    than the saved famid or they will have the same famid as the saved famid 
    and the same id as the saved id.  In only the case with a different famid
    will we want to save over the saved famid and id.

    if old fam != new fam || new id == saved id
      if old fam != new fam
        save fam
	save id
        point to new family
      increase family size
  */

  // Do a little length checking on famids and ids
  Sint famidsLength = LENGTH(famids);
  Sint idsLength = LENGTH(ids);

  if(famidsLength != idsLength) {
    Rcerr << "The length of famids (" << famidsLength << ") does not match "
	  << "the length of ids (" << idsLength << ")." << std::endl
	  << "getFamilySizes.cpp key 41" << std::endl;
  }

  // Create the storage space for the famid sizes.  The length of familySizes
  // and familyLables is originally the max size of famids.  In the largest
  // scenario, each person would be their own family, and thus we have
  // famidsLength amount of families.
  s_object *familySizes = NEW_INTEGER(famidsLength);

  /*#ifdef USING_R
  int
#else
  Sint
#endif
  */
  Sint *familySizesPtr = INTEGER_POINTER(familySizes);
  Sint familySizesIndex = 0;
  familySizesPtr[0] = 1;

  // Save the first family's famid snd first person's id
  const char *savedFamid = multic_STRING_ELT(famids, 0);
  const char *savedId = multic_STRING_ELT(ids, 0);

  // Create the storage space for the famid labels for the famid sizes
  s_object *familyLabels = NEW_CHARACTER(famidsLength);

  multic_SET_STRING_ELT(familyLabels, familySizesIndex, savedFamid);

  for(Sint i = 1; i < famidsLength; i++) {
    if(strcmp(savedFamid, multic_STRING_ELT(famids, i)) != 0
       || strcmp(savedId, multic_STRING_ELT(ids, i)) == 0) {
      // Since we've found a new family, increase the familySizesIndex
      familySizesIndex++;
      // famidSizesPtr should be set to 1 (because we now know that there is 
      // at least one person in this family), but the ++ operater below does
      // that for us.
      familySizesPtr[familySizesIndex] = 0;

      // Save the famid and id for later comparison
      savedFamid = multic_STRING_ELT(famids, i);
      savedId = multic_STRING_ELT(ids, i);

      multic_SET_STRING_ELT(familyLabels, familySizesIndex, savedFamid);
    }

    familySizesPtr[familySizesIndex]++;
  }

  // Remove any trailing 0's
  SET_LENGTH(familySizes, familySizesIndex + 1);
  SET_LENGTH(familyLabels, familySizesIndex + 1);

  SET_NAMES(familySizes, familyLabels);

  return familySizes;
}

}
