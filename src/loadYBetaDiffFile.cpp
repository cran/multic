#include <fstream>
#include <iostream>
#include <cstring>

#include "verS.h"
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif
using namespace std;

extern "C" {

/*********************
Title: loadYBetaDiffFile
Description: reads the file name specified in fileName, generates a list of
             numeric vectors, and returns that list
Input: s_object *fileName - character value specifying the name of the file
                            to be read
       s_object *familyCount - integer value specifying the number of
                               families in the givin file name
Output: s_object *list - a list of numeric vectors
Side Effects: NONE
Author: Eric Lunde, 12/12/2004
*********************/
s_object *loadYBetaDiffFile(s_object *fileName, s_object *familyCount) {
  S_EVALUATOR

  const char *fileNameValue = CHARACTER_VALUE(fileName);
  Sint familyCountValue = INTEGER_VALUE(familyCount);

  s_object *list = NEW_LIST(familyCountValue);
 
  ifstream fin(fileNameValue);
  if(fin.fail()) {
    PROBLEM "The file '%s' could not be opened for reading.\nloadYBetaDiff.cpp key 20\n",
      fileNameValue RECOVER(NULL_ENTRY);
  }

  char holder[1024];
  Sint familySize;
  Sint traitCount;
  s_object **yBetaDiffs;
  double *yBetaDiffsPtr;

  yBetaDiffs = (s_object **) Salloc(familyCountValue, s_object *);

  for(int fam = 0; fam < familyCountValue; fam++) {
    fin >> holder >> holder >> holder >> familySize >> holder >> traitCount;
    
    yBetaDiffs[fam] = NEW_NUMERIC(familySize * traitCount);
    yBetaDiffsPtr = NUMERIC_POINTER(yBetaDiffs[fam]);
    for(int i = 0; i < familySize * traitCount; i++) {
      fin >> yBetaDiffsPtr[i];
    }
    
    SET_ELEMENT(list, fam, yBetaDiffs[fam]);
  }

  fin.close();

  return list;
}
  
}

