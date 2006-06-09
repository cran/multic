#include <iostream>
#include <cstring>
#include "verS.h"
using namespace std;

extern "C" {

void isSpouse(char **x, char **y, Sint *idsLength, char **dadid, char **momid,
	      Sint *dadidLength, Sint *familySizes, Sint *familyLength,
	      Sint *isSpouse) {

  Sint *familyStartingIndex = new Sint[*familyLength + 1];
  familyStartingIndex[0] = 0;
  for(Sint i = 1; i < *familyLength + 1; i++) {
    familyStartingIndex[i] = familyStartingIndex[i - 1] + familySizes[i - 1];
  }

  Sint familyIndex = 0;
  Sint numberOfValidIds =
	(familySizes[familyIndex] - 1) * familySizes[familyIndex] / 2;
  Sint numberOfIdsUsed = 0;

  for(Sint i = 0; i < *idsLength; i++) {
    isSpouse[i] = 0;
    
    for(Sint j = familyStartingIndex[familyIndex];
	j < familyStartingIndex[familyIndex + 1]; j++) {
      if( (strcmp(dadid[j], x[i]) == 0 && strcmp(momid[j], y[i]) == 0)
	  || (strcmp(dadid[j], y[i]) == 0 && strcmp(momid[j], x[i]) == 0) ) {
	isSpouse[i] = 1;
	break;
      }
    }

    numberOfIdsUsed++;

    if(numberOfIdsUsed == numberOfValidIds) {
      familyIndex++;
      numberOfValidIds =
	(familySizes[familyIndex] - 1) * familySizes[familyIndex] / 2;
      numberOfIdsUsed = 0;
    }
  }
  
  delete [] familyStartingIndex;
}

}
