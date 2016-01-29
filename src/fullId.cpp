#include <iostream>
#include <cstring>
#include "verS.h"

extern "C" {

void fullId(char **id, Sint *idLength, Sint *familySizes,
	    Sint *familySizesLength, Sint *loci, char **fullId,
	    Sint *fullIdLength, Sint *idMaxLength) {

  Sint fullIdIndex = 0;
  Sint baseFamilyIndex = 0;
  Sint familySize = 0;

  if(*idLength * *loci != *fullIdLength) {
    PROBLEM "*idLength * *loci = %d, *fullIdLength = %d\nThe values need to be equal.\nfullId.cpp key 18\n",
      (int)(*idLength * *loci), (int)*fullIdLength RECOVER(NULL_ENTRY);
  }

  for(Sint i = 0; i < *familySizesLength; i++) {
    familySize = familySizes[i];
    for(Sint j = 0; j < *loci; j++) {
      for(Sint k = 0; k < familySize; k++) {
        //      fullId[fullIdIndex] = id[baseFamilyIndex + k];
        //      cout << "i = " << i << ", j = " << j << ", k = " << k << endl;
        strncpy(fullId[fullIdIndex], id[baseFamilyIndex + k], *idMaxLength);
        fullIdIndex++;
      }
    }
    baseFamilyIndex += familySize;
  }

  if(*fullIdLength != fullIdIndex) {
    PROBLEM "*fullIdLength = %d, fullIdIndex = %d\nThe values need to be equal.\nfullId.cpp key 29\n",
      (int)*fullIdLength, (int)fullIdIndex RECOVER(NULL_ENTRY);
  }
  if(*idLength != baseFamilyIndex) {
    PROBLEM "*idLength = %d, baseFamilyIndex = %d\nThe values need to be equal.\nfullId.cpp key 35\n",
      (int)*idLength, (int)baseFamilyIndex RECOVER(NULL_ENTRY);
  }
}

}
