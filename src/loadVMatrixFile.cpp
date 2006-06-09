#include <fstream>
#include <iostream>
#include <cstring>

#include <S.h>
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif
using namespace std;

extern "C" {

void swap(double *matrix, int i, int j, int familySize);
void transpose(double *matrix, int familySize);

/*********************
Title: loadVMatrixFile
Description: reads the file name specified in fileName, generates a list of
             numeric matrices, and returns that list
Input: s_object *fileName - character value specifying the name of the file
                            to be read
       s_object *familyCount - integer value specifying the number of
                               families in the givin file name
Output: s_object *list - a list of numeric matrices
Side Effects: NONE
Author: Eric Lunde, 12/12/2004
*********************/
s_object *loadVMatrixFile(s_object *fileName, s_object *familyCount) {
  S_EVALUATOR

  char *fileNameValue = CHARACTER_VALUE(fileName);
  long familyCountValue = INTEGER_VALUE(familyCount);

  s_object *list = NEW_LIST(familyCountValue);
 
  ifstream fin(fileNameValue);
  if(fin.fail()) {
    PROBLEM "The file '%s' could not be opened for reading.\nloadVMatrix.cpp key 23\n",
      fileNameValue RECOVER(NULL_ENTRY);
  }

  char holder[1024];
  long familySize;
  long traitCount;
  long matrixRowCount;
  long matrixElementCount;
  s_object **vMatrix;
  double *vMatrixPtr;
  s_object **vMatrixDim;
#ifdef USING_R
  int *
#else
  long *
#endif
    vMatrixDimPtr;

  vMatrix = (s_object **) Salloc(familyCountValue, s_object *);
  vMatrixDim = (s_object **) Salloc(familyCountValue, s_object *);
  for(int i = 0; i < familyCountValue; i++) {    
    vMatrixDim[i] = NEW_INTEGER(2);
  }

  for(int fam = 0; fam < familyCountValue; fam++) {
    //    cout << "fam = " << fam << endl;
    fin >> holder >> holder >> holder >> familySize >> holder >> traitCount;
    matrixRowCount = familySize * traitCount;
    matrixElementCount = matrixRowCount * matrixRowCount;
    
    vMatrix[fam] = NEW_NUMERIC(matrixElementCount);
    vMatrixPtr = NUMERIC_POINTER(vMatrix[fam]);
    for(int i = 0; i < matrixElementCount; i++) {
      fin >> vMatrixPtr[i];
    }
    
    transpose(vMatrixPtr, matrixRowCount);
    
    vMatrixDimPtr = INTEGER_POINTER(vMatrixDim[fam]);
    vMatrixDimPtr[0] = matrixRowCount;
    vMatrixDimPtr[1] = matrixRowCount;
    SET_DIM(vMatrix[fam], vMatrixDim[fam]);
    
    SET_ELEMENT(list, fam, vMatrix[fam]);
  }

  fin.close();

  return list;
}
  
/*********************
Title: transpose
Description: (shhh!  it doesn't really transpose a matrix) given a (logical)
             square matrix (actually an array), swap the values along the
             major diagonal
Input: double *matrix - an array representing a sqpare matrix
       int rowCount - the number of rows (and columns for that matter)
Output: NONe
Side Effects: NONE
Author: Eric Lunde, 12/12/2004
*********************/
void transpose(double *matrix, int rowCount) {
  for(int i = 0; i < rowCount; i++) {
    for(int j = i+1; j < rowCount; j++) {
      swap(matrix, i, j, rowCount);
    }    
  }
}

void swap(double *matrix, int i, int j, int rowCount) {
  if(i == j) {
    return;
  }
  double temp = matrix[j + i * rowCount];
  matrix[j + i * rowCount] = matrix[i + j * rowCount];
  matrix[i + j * rowCount] = temp;
}

}

