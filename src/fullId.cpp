#include <iostream>
#include <cstring>
#include <S.h>
using namespace std;

extern "C" {

void fullId(char **id, long *idLength, long *familySizes,
	    long *familySizesLength, long *loci, char **fullId,
	    long *fullIdLength, long *idMaxLength) {

  long fullIdIndex = 0;
  long baseFamilyIndex = 0;
  long familySize = 0;

  if(*idLength * *loci != *fullIdLength) {
    PROBLEM "*idLength * *loci = %ld, *fullIdLength = %ld\nThe values need to be equal.\nfullId.cpp key 18\n",
      (*idLength * *loci), *fullIdLength RECOVER(NULL_ENTRY);
  }

  for(long i = 0; i < *familySizesLength; i++) {
    familySize = familySizes[i];
    for(long j = 0; j < *loci; j++) {
      for(long k = 0; k < familySize; k++) {
	//	fullId[fullIdIndex] = id[baseFamilyIndex + k];
	//	cout << "i = " << i << ", j = " << j << ", k = " << k << endl;
	strncpy(fullId[fullIdIndex], id[baseFamilyIndex + k], *idMaxLength);
	fullIdIndex++;
      }
    }
    baseFamilyIndex += familySize;
  }

  if(*fullIdLength != fullIdIndex) {
    PROBLEM "*fullIdLength = %ld, fullIdIndex = %ld\nThe values need to be equal.\nfullId.cpp key 29\n",
      *fullIdLength, fullIdIndex RECOVER(NULL_ENTRY);
  }
  if(*idLength != baseFamilyIndex) {
    PROBLEM "*idLength = %ld, baseFamilyIndex = %ld\nThe values need to be equal.\nfullId.cpp key 35\n",
      *idLength, baseFamilyIndex RECOVER(NULL_ENTRY);
  }
}

}
