#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "multicString.h"

#include <S.h>
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif
using namespace std;

int *getVarianceIndices(long traitCountValue);
int *getCovarianceIndices(long traitCountValue);
s_object *getCorrelations(s_object *values, int *varianceIndices,
			  int *covarianceIndices, long traitCountValue);
void setDim(s_object *object, int rows, int columns);
void setDimNames(s_object *object, s_object *values, int traitCountValue);
char *strdupS(char *string);
const int MISSING_VALUE = -9;

extern "C" {

/*********************************
Title: calculateCorrelations
Description: calculate the correlations by dividing the covariance value by 
             the square root of the product of its two variances
Input: s_object *values - numeric matrix of variance and covariance values
                          where each column in a new hypothesis and each row
                          is a variance or covariance for that hypothesis
       s_object *traitCount - the number of traits so we know the dimensions
                              of values (this value could be extracted from
                              values directly, but I learned this midway
                              through programming and chose not to fix it)
Output: s_object *cors - the correlation calculations with dimensions and
                         dimnames set
Side Effects: NONE
Author: Eric Lunde, 01/27/05
*********************************/
s_object *calculateCorrelations(s_object *values, s_object *traitCount) {
  S_EVALUATOR
    
  // Get the actual value out of the s_object
  long traitCountValue = INTEGER_VALUE(traitCount);

  // Get the length of values
  long valuesLength = LENGTH(values);

  // Based on traitCountValue, determine how much of the values each
  // hypothesis takes.
  int hypothesisLength = traitCountValue * (traitCountValue + 1)/ 2;
  int hypothesisCount = valuesLength / hypothesisLength;
  int covarianceLength = hypothesisLength - traitCountValue;

  // Then, based on hypothesisLength, determine which indices point to
  // variance values and which point to covariance values and set
  // varianceIndices and covarianceIndices to their apporpriate values.
  int *varianceIndices = getVarianceIndices(traitCountValue);
  int *covarianceIndices = getCovarianceIndices(traitCountValue);

  // Send off the work to getCorrelations 
  s_object *cors = getCorrelations(values, varianceIndices,
				   covarianceIndices, traitCountValue);
  setDim(cors, hypothesisCount, covarianceLength);
  setDimNames(cors, values, traitCountValue);

  // Clean up the finished memory
  delete [] varianceIndices;
  delete [] covarianceIndices;

  return cors;
}

}

int *getVarianceIndices(long traitCountValue) {
  // Create the array to hold the index values
  int *varianceIndices = new int[traitCountValue];

  // Have a variable to hold the same value as traitCountValue, becuase we
  // want to keep the value traitCountValue holds and alter it as well.
  long indexIncrease = traitCountValue;
  
  // Set the first (and many times only) index to 0;
  varianceIndices[0] = 0;

  // For the rest, take the previous value and add indexIncrease to it,
  // while deceasing indexIncrease every iteration.
  for(int i = 1; i < traitCountValue; i++) {
    varianceIndices[i] = varianceIndices[i - 1] + indexIncrease--;
  }

  // Do a little error checking during development.
  if(indexIncrease != 1) {
    PROBLEM "indexIncrease = %ld\nError occurred calculating variance indices.\ncalculateCorrelations.cpp key 84\n",
      indexIncrease RECOVER(NULL_ENTRY);
  }

  return varianceIndices;
}

int *getCovarianceIndices(long traitCountValue) {
  // Get the variane indices from our other function, instead of passing it
  // as a parameter.
  int *varianceIndices = getVarianceIndices(traitCountValue);

  int lastIntegerValue = varianceIndices[traitCountValue - 1];

  // Create the array to hold the index values.
  int *covarianceIndices = new int[lastIntegerValue - traitCountValue + 1];

  int covarianceIndicesIndex = 0;
  int varianceIndicesIndex = 0;

  // Taking i from 0 to lastIntegerValue, if i is not in varianceIndices, 
  // save it in covarianceIndices.
  for(int i = 0; i < lastIntegerValue + 1; i++) {
    if(i != varianceIndices[varianceIndicesIndex]) {
      covarianceIndices[covarianceIndicesIndex++] = i;
    }else {
      varianceIndicesIndex++;
    }
  }

  delete [] varianceIndices;

  return covarianceIndices;
}

s_object *getCorrelations(s_object *values, int *varianceIndices,
			  int *covarianceIndices, long traitCountValue) {

  double *valuesPtr = 
#ifdef USING_R
    NUMERIC_POINTER(values)
#else
    NUMERIC_POINTER(GET_DATA(values))
#endif
    ;
  long valuesLength = LENGTH(values);

  // These index values are used to access the varianceIndices and
  // covarianceIndices arrays.  The produce other indices to be used to
  // access the values object.
  int covarIndex = 0, correctHypothesis = 0;

  // The calculations to get the covariance, variance1, variance2,
  // and correlation are complicated.  I'll save their values in these temp
  // variables to separate calculations.
  double covariance, variance1, variance2, correlation;

  int varianceIndicesLength = traitCountValue;
  int hypothesisLength = traitCountValue * (traitCountValue + 1)/ 2;
  int covarianceIndicesLength = hypothesisLength - traitCountValue;
  long hypothesesCount = valuesLength / hypothesisLength;

  // Create a correlations s_object to hold the correlations.  It will be a
  // matrix of the number of total hypotheses rows (calculated by 
  // hypothesesCount) by the number of covariances (or correlations, same
  //  value) rows (calculated by covarianceIndicesLength).
  s_object *correlations = NEW_NUMERIC(covarianceIndicesLength
				       * hypothesesCount);
  double *correlationsPtr = NUMERIC_POINTER(correlations);

  // For each hypothesis, calculated all of the correlations within that 
  // hypthesis.  This entails for each unique variance combination, calculate
  // their covariance correlation and then place that value in its correct 
  // place in correlations.
  for(int i = 0; i < hypothesesCount; i++) {
    covarIndex = 0;
    correctHypothesis = i * hypothesisLength;
    for(int j = 0; j < varianceIndicesLength - 1; j++) {
      for(int k = j + 1; k < varianceIndicesLength; k++) {
	// correctHypothesis determines the hypothesis and the
	// (co)varianceIndices determine the offset inside the hypothesis.
	covariance = valuesPtr[correctHypothesis
			       + covarianceIndices[covarIndex]];
	variance1 = valuesPtr[correctHypothesis + varianceIndices[j]];
        variance2 = valuesPtr[correctHypothesis + varianceIndices[k]];
	
	if(variance1 == 0 || variance2 == 0) {
	  correlation = MISSING_VALUE;
	} else {
	  correlation = covariance / sqrt(variance1 * variance2);
	}

	// hypothesisCount equals the number of rows in a single column, so
	// covarIndex * hypothesisCount calculates the correct (covarIndex)
	// column and i offsets to the correct correlation.
	correlationsPtr[covarIndex * hypothesesCount + i] = correlation;
	covarIndex++;
      }
    }
    if(covarIndex != covarianceIndicesLength) {
      PROBLEM "covarIndex != covarianceIndicesLength (%d != %d)\ncalculateCorrelations.cpp key 135\n",
	covarIndex, covarianceIndicesLength RECOVER(NULL_ENTRY);
    }    
  }
  
  return correlations;
}

void setDim(s_object *object, int rows, int columns) {
  s_object *objectDim = NEW_INTEGER(2);
#ifdef USING_R
  int *
#else
  long *
#endif
    objectDimPtr = INTEGER_POINTER(objectDim);

  objectDimPtr[0] = rows;
  objectDimPtr[1] = columns;

  SET_DIM(object, objectDim);
}

void setDimNames(s_object *object, s_object *values, int traitCountValue) {
  // Get the correlation (also called covariate) indices so we know which
  // names to retrieve.  
  int *correlationIndices = getCovarianceIndices(traitCountValue);
  int correlationIndicesLength = (traitCountValue - 1) * traitCountValue / 2;
  int correlationIndicesIndex = 0;

  // Create a new character vector to hold the covariance (correlation) 
  // names.
  s_object *correlationNames = NEW_CHARACTER(correlationIndicesLength);

  // Copy the deisred correlation names to the vector.  These rowNames will
  // become the column names of the correlation object
  s_object *rowNames = 
#ifdef USING_R
    GET_ROWNAMES(GET_DIMNAMES(values));
#else
    GET_ROWNAMES(values);
#endif

  for(int i = 0; i < LENGTH(rowNames); i++) {
    if(i == correlationIndices[correlationIndicesIndex]) {
      multic_SET_STRING_ELT(correlationNames, correlationIndicesIndex,
			    multic_STRING_ELT(rowNames, i));
      correlationIndicesIndex++;
    }
  }
    
  s_object *dimNames = NEW_LIST(2);
  // Make the row names the same as the column names of the input
  s_object *colNames = 
#ifdef USING_R
    GET_COLNAMES(GET_DIMNAMES(values))
#else
    GET_COLNAMES(values)
#endif
    ;
  SET_ELEMENT(dimNames, 0, colNames);
  // Make the column names the names of the covariances of the input.  The 
  // correlationNames become the column names of the correlation object.
  SET_ELEMENT(dimNames, 1, correlationNames);
  delete [] correlationIndices;

  SET_DIMNAMES(object, dimNames);
}

char *strdupS(char *string) {
  char *newString = Salloc(strlen(string) + 1, char);
  strcpy(newString, string);
  return newString;
}
