#include <fstream>
#include <iostream>
#include <cstring>
#include "multicString.h"

#include "verS.h"
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif

extern "C" {

bool hasPreviousEntry(const char *id, 
		      s_object *shareOutOrder,
		      int shareOutOrderIndex,
		      int depth);

s_object *getShareOutOrder(s_object *fileName, s_object *peopleCount, 
			   s_object *usingMloci) {
  S_EVALUATOR

  const char *fileNameValue = CHARACTER_VALUE(fileName);
  Sint peopleCountValue = INTEGER_VALUE(peopleCount);
  bool usingMlociValue = INTEGER_VALUE(usingMloci);
  char firstIdBuffer[1024];
  char secondIdBuffer[1024];
  char famidBuffer[1024];
  char famid[1024];
  char FAMILY_ID_DELIMITER[2] = "-";
  Sint shareOutOrderIndex = 0L;
  int currentFamilyCount = 0;

  if(peopleCountValue < 1) {
    PROBLEM "peopleCountValue (%ld) cannot be less than 1.\ngetShareOutOrder.cpp key 18\n",
      peopleCountValue RECOVER(NULL_ENTRY);
  }

  s_object *shareOutOrder = NEW_CHARACTER(peopleCountValue);
  
  std::ifstream shareOut(fileNameValue);
  if(shareOut.fail()) {
    PROBLEM "The file '%s' could not be opened for reading.\ngetShareOutOrder.cpp key 29\n",
      fileNameValue RECOVER(NULL_ENTRY);
  }

  // Read first and second person and save, ignore the rest
  // Save famid
  // Repeat until done with file
  // Save first
  // If first is not equal to the previous person
  //   If famid is not equal, append second person to vector, save famid
  //   Append first person to vector
  //   Save second person
  // Ignore rest, eat '\n'

  // Read first line if using mloci.out
  if(usingMlociValue) {
    shareOut.ignore(1024, '\n');
  }

  // Get the first unique id from the family
  shareOut >> firstIdBuffer;
  multic_SET_STRING_ELT(shareOutOrder, shareOutOrderIndex,
			firstIdBuffer);
  shareOutOrderIndex++;
  currentFamilyCount++;

  // Save the family id
  strcpy(famid, strtok(firstIdBuffer, FAMILY_ID_DELIMITER));

  // Get the second unique id from the family
  shareOut >> secondIdBuffer;

  // Read the rest of the line
  shareOut.ignore(1024, '\n');

  while(shareOutOrderIndex+1 != peopleCountValue) {
    shareOut >> firstIdBuffer;

    if(strcmp(firstIdBuffer,
	      multic_STRING_ELT(shareOutOrder, shareOutOrderIndex - 1)
	      ) != 0) {
      // Compare famids
      strcpy(famidBuffer, firstIdBuffer);
      strtok(famidBuffer, FAMILY_ID_DELIMITER);
      
      if(strcmp(famid, famidBuffer) != 0
	 || (hasPreviousEntry(firstIdBuffer,
			      shareOutOrder,
			      shareOutOrderIndex - 1,
			      currentFamilyCount))
	 ) {
	// Save the second id, because we are about to change families
	multic_SET_STRING_ELT(shareOutOrder, shareOutOrderIndex,
			      secondIdBuffer);
	shareOutOrderIndex++;
	currentFamilyCount = 0;

	// Save new famid
	strcpy(famid, famidBuffer);
      }
      
      // Save the first id, it is a new person
      multic_SET_STRING_ELT(shareOutOrder, shareOutOrderIndex,
			    firstIdBuffer);
      shareOutOrderIndex++;
      currentFamilyCount++;

      shareOut >> secondIdBuffer;
    }
    
    // Read the rest of the line
    shareOut.ignore(1024, '\n');
  }

  // Add the last remaining secondIdBuffer to the array.
  multic_SET_STRING_ELT(shareOutOrder, shareOutOrderIndex,
			secondIdBuffer);
  shareOutOrderIndex++;
  currentFamilyCount++;

  shareOut.close();

  return shareOutOrder;
}

bool hasPreviousEntry(const char *id, 
		      s_object *shareOutOrder,
		      int shareOutOrderIndex,
		      int depth) {
  /* Starting at sharOutOrderIndex going to shareOutOrderIndex - depth,
     determine if id is in shareOutOrderPtr */
  bool found = false;

  for(int i = shareOutOrderIndex; i > shareOutOrderIndex - depth; i--) {
    if(strcmp(id, multic_STRING_ELT(shareOutOrder, i)) == 0) {
      found = true;
      break;
    }
  }

  return found;
}

}
