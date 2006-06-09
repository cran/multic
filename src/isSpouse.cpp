#include <iostream>
#include <cstring>
using namespace std;

extern "C" {

void isSpouse(char **x, char **y, long *idsLength, char **dadid, char **momid,
	      long *dadidLength, long *familySizes, long *familyLength,
	      long *isSpouse) {

  long *familyStartingIndex = new long[*familyLength + 1];
  familyStartingIndex[0] = 0;
  for(long i = 1; i < *familyLength + 1; i++) {
    familyStartingIndex[i] = familyStartingIndex[i - 1] + familySizes[i - 1];
  }

  long familyIndex = 0;
  long numberOfValidIds =
	(familySizes[familyIndex] - 1) * familySizes[familyIndex] / 2;
  long numberOfIdsUsed = 0;

  for(long i = 0; i < *idsLength; i++) {
    isSpouse[i] = 0;
    
    for(long j = familyStartingIndex[familyIndex];
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
