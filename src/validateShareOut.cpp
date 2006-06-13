#include <fstream>
#include <iostream>
#include <S.h>
using namespace std;

extern "C" {

/****************************************************************************
 * validateShareOut is a routine that makes sure that the ids that are about
 * to be used to create fort.12 are in fact the same ids in the named
 * shareOutName (most likely share.out) and represent the same ordering that
 * has already been defined in shareOutName (most likely share.out).
 ****************************************************************************/
void validateShareOut(long *familySizes, long *familiesCount, char **id,
		      long *idsCount, char **shareOutName,
		      long *passedValidation)
{
  // Perform some bounds checking before using familySizes, familiesCount,
  // id, and idsCount
  if( *familiesCount < 1 ) {
    PROBLEM "familiesCount was passed a value of '%ld'.\nThis value must be positive.\nvalidateShareOut.cpp key 14\n",
      *familiesCount RECOVER(NULL_ENTRY);
  }

  for(long l = 0; l < *familiesCount; l++) {
    if(familySizes[l] < 1) {
      PROBLEM "familySizes[%ld] was passed a value of '%ld'.\nThis value must be positive.\nvalidateShareOut.cpp key 23\n",
	l, familySizes[l] RECOVER(NULL_ENTRY);
    }
  }

  if( *idsCount < 1 ) {
    PROBLEM "idsCount was passed a value of '%ld'.\nThis value must be positive.\nvalidateShareOut.cpp key 30\n",
      *idsCount RECOVER(NULL_ENTRY);
  }

  /* Commented out on 11/9/2004 by Eric Lunde because id changed from type
     long * to char **
  for(long l = 0; l < *idsCount; l++) {
    if(id[l] < 1) {
    PROBLEM "id[%d] was passed a value of '%d'.\nThis value must be positive.\nvalidateShareOut.cpp key 38\n",
    l, id[l] RECOVER(NULL_ENTRY);
    }
  }
  */

  ifstream shareOut(shareOutName[0]);
  if( !shareOut.is_open() ) {
    PROBLEM "The file '%s' could not be opened.\nvalidateShareOut.cpp key 46\n",
      shareOutName[0] RECOVER(NULL_ENTRY);
  }
    
  // Since we are dealing with ids in an array format, firstIdIndex and
  // secondIdIndex represent the indices of those two values.
  long firstIdIndex, secondIdIndex;
  
  // i, j, and k are used only as for loop iteration variables.
  long i, j, k;

  // firstId and secondId are the actual id values specified in the array
  // 'id'.  shareOutFirstId and shareOutSecondId are the actual id values
  // taken from the file shareOutName (most likely share.out).
  //  long firstId, secondId, shareOutFirstId, shareOutSecondId;
  char *firstId, *secondId, shareOutFirstId[32], shareOutSecondId[32];

  // shareOutLine Number is a variable to keep track of what line number of
  // shareOutName (most likely, share.out) we are currently extracting data
  // from.  This is used to aid the user when an inconsistancy has been found.
  long shareOutLineNumber = 0;

  // firstIdIndex begins at -2 because when switching families in the linear
  // data, we need to shift the firstIdIndex two places (once for it to index
  // the last person in the family, but since they are the last person in the
  // family, we need to increment it again to be the correct index for the
  // first person in the next family.
  firstIdIndex = -2;

  // For all families, process the current family.
  for (i = 0; i < *familiesCount; i++) {
    // Increment the index to point to the last person in the previous family
    firstIdIndex++;

    // Processing means make sure for a family that the order of the ids is
    // the same order represented by the sequence of id combinations in
    // shareOutName.
    for (j = familySizes[i]-1; j > 0; j--) {
      // Increment the index to point to the next person in the family (in 
      // general, when switching families, we actually increment the index
      // to point to the first person of the next family.
      firstIdIndex++;

      // Since shareOutName only lists combinations between two different
      // people, we need to have the secondIdIndex be one position greater
      // than the firstIdIndex.
      secondIdIndex = firstIdIndex + 1;

      // Get the id of the first person.
      firstId = id[firstIdIndex];
      
      // Read a line of shareOut and compare ids.  Reporting and returning
      // if a mismatch is found.
      for(k = 0; k < j; k++) {
	// Get the id of the second person.
	secondId = id[secondIdIndex++];
	
	// Read the id data from shareOutName.
	shareOutLineNumber++;
	shareOut >> shareOutFirstId >> shareOutSecondId;
	shareOut.ignore(512, '\n');

	// If an error is found, do not exit the program, just print an error
	// message, set the passedValidation value to 0 (false), and return.
	// I think it is better to not quit Splus, just the function call.
	//	if(firstId != shareOutFirstId || secondId != shareOutSecondId) {
	if(strcmp(firstId, shareOutFirstId) != 0
	   || strcmp(secondId, shareOutSecondId) != 0) {
	  cerr << "The id combination, " << firstId << " and " << secondId 
	       << ", does not match the" << endl << "file's id combination, "
	       << shareOutFirstId << " and " << shareOutSecondId
	       << ", on line number: " << shareOutLineNumber << " of " << endl
	       << shareOutName[0] << endl << "Data validation failed."
	       << endl;
	  *passedValidation = 0;
	  return;
	}
      }
    }
  }

  shareOut.close();

  // Set the passedValidation variable to 1 (true) and return.
  *passedValidation = 1;
}

}
