/******************************************************************************
Title: loadEffects
Description: loadEffects is designed to be called from Splus during the
             execution of the multic program.  It opens and reads from
             "summary.log" and saves its information in Splus-allocated memory.
Input: char **lociNames - Each (row) index of lociNames will be filled with
                          the name of the ibd file that multic used as an
                          alternative hypothesis or some string specifying that
                          the information is from the null hypothesis.
       char **variableNames - Each (row) index of variableNames will contain
                              the name of the variable that each
                              estimate/standard error represent.
       double *fixedEffects - Each index of fixedEffects will contain an
                              estimate or standard error for one variable at
                              one locus.  The ordering of fixedEffects is as
                              follows: the overall mean estimate, each trait's
                              mean, overall mean standard error, and each
                              trait's standard errors.  This pattern repeats
                              for each different alternative hypothesis.
       double *polygene - Each index of polygene contains an estimate or
                          standard error for one variable at one locus.  The
                          ordering of polygene is as follows: all polygene
                          variable estimates for one locus followed by all
                          polygene variable standard errors for one locus.
                          This pattern repeats for all alternative hypotheses.
       double *majorGene1 - Similar description as polygene except this will
                            hold major gene informaion.
       double *environment -  Similar desription as polygene except this will
                              hold environment informaion.
       double *logLikelihoods - Each index of logLikelihoods will contain the
                                log likelihood for the entire locus (or the
                                sum of the family log likelihoods for that
                                locus).
       char **logLikStatus - Each (row) index of logLikStatus will contain a
                             string representing whether the loglikelihood
                             represnts a convervgent or non-convergent value.
       long *traitCount - An integer value specifying the number of traits that
                          were used in the multic calculations.
       long *covariateCount - An integer value specifying the number of
                              covariates that were used in the multic
                              calculations.
       long *lociCount - An integer value specifying the number of loci that
                         were used in the multic calculations.  This value does
                         include the null hypothesis calculation.
Output: NONE
Side Effects: This function reads from "summary.log". It is not user-defined
              name.  It is hardcoded.  It must be in the same directory as the
              Splus session is occuring.  If it does not exist, and error
              message will be printed and the program will quit.
Author: Eric Lunde 10-17-03
Update: (Date, Modified By, Modification Description)
03-31-2005, Eric Lunde, The number of environmental variance (covariance) 
            values no longer is dependent upon whether multic is run in
            multivariate or longitudinal mode (it now always operates as if
            in multivariate mode).  Therefore I've removed the
            environmentCount variable and replaced it with the
            randomEffectsCount which now always is the correct value.
******************************************************************************/
#include <fstream>
#include <iostream>
#include <S.h>
using namespace std;

extern "C" {

void loadEffects(char **lociNames,
		 char **variableNames,
		 double *fixedEffects,
		 double *polygene,
		 double *majorGene1,
		 double *environment,
		 double *siblingSibling,
		 double *parentParent,
		 double *parentOffspring,
		 double *logLikelihoods,
		 char **logLikStatus,
		 long *fixedEffectsCountLong,
		 long *randomEffectsCountLong,
		 long *lociCount) {

  // Open for reading "summary.log" and test for failure
  ifstream summaryLog("summary.log");
  if(summaryLog.fail()) {
    PROBLEM "The file summary.log could not be opened for reading.\nloadEffects.cpp key 33\n"
      RECOVER(NULL_ENTRY);
  }

  // newLociIndicator will hold the special character "#" seperating the
  // logical "multic.out"s (seperate alternative hypothesis information) in
  // the file summary.log.
  char newLociIndicator[8];
  // The file "summary.log" contains more information than multic currently
  // needs.  Instead of not writing those values to "summary.log" from multic,
  // we will just not read them back into multic.  placeHoldingDouble is the
  // space we use to continue parsing the file without having to save the
  // informaion.  placeHoldingString serves the same purpose.
  double placeHoldingDouble;
  char placeHoldingString[128];

  // fixedEffectsCount is the total number of fixed effects (overall mean and
  // covariate means) per alternative hypothesis.
  int fixedEffectsCount = //(int)(*traitCount * (1 + *covariateCount));
    (int) *fixedEffectsCountLong;
  int fixedEffectsIndex = 0;
  // randomEffectsCount is the total number of random effects per effect per
  // locus.
  int randomEffectsCount = (int) *randomEffectsCountLong;

  // Since ibd file names, polygene, major gene, and environment information
  // is dispursed throughout "summary.log", we need an independent means of
  // storing them in their own locations.  These variables grant us this
  // ability.
  int variableNamesIndex = 0;
  int polyGeneEffectsIndex = 0;
  int majorGeneEffectsIndex = 0;
  int environmentEffectsIndex = 0;
  int siblingSiblingEffectsIndex = 0;
  int parentParentEffectsIndex = 0;
  int parentOffspringEffectsIndex = 0;

  // Note: The logic behind code like this: 
  // summaryLog >> fixedEffects[fixedEffectsIndex]
  //            >> fixedEffects[fixedEffectsIndex + fixedEffectsCount];
  // is to store all of the necessary information in the following order.
  // First save all of the estimates for fixed effects, or polygene, or etc for
  // one locus in contiguous memory, followed by the standard errors (also in
  // contiguous memory) for those variables directly after the group of
  // estimates.  After saving that information, update fixedEffectsIndex with
  // the statemtent: fixedEffectsIndex += fixedEffectsCount; to reposition the
  // fixedEffectsIndex value.
  for(int i=0; i<*lociCount; i++) {
    summaryLog >> newLociIndicator >> lociNames[i];

    // Read fixed effects data
    for(int j=0; j<fixedEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> fixedEffects[fixedEffectsIndex]
		 >> fixedEffects[fixedEffectsIndex + fixedEffectsCount];
      fixedEffectsIndex++;
    }
    fixedEffectsIndex += fixedEffectsCount;

    // Read Polygene data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> polygene[polyGeneEffectsIndex]
		 >> polygene[polyGeneEffectsIndex + randomEffectsCount];
      polyGeneEffectsIndex++;
    }
    polyGeneEffectsIndex += randomEffectsCount;
    
    // Read MajorGene1 data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> majorGene1[majorGeneEffectsIndex]
		 >> majorGene1[majorGeneEffectsIndex + randomEffectsCount];
      majorGeneEffectsIndex++;
    }
    majorGeneEffectsIndex += randomEffectsCount;
    
    // Read (but ignore) MajorGene2 data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> placeHoldingString;
      }
      summaryLog >> placeHoldingDouble >> placeHoldingDouble;
    }

    // Read Environment data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> environment[environmentEffectsIndex]
		 >> environment[environmentEffectsIndex
				+ randomEffectsCount];
      environmentEffectsIndex++;
    }
    environmentEffectsIndex += randomEffectsCount;

    // Read Shared Shiship data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> siblingSibling[siblingSiblingEffectsIndex]
		 >> siblingSibling[siblingSiblingEffectsIndex
				   +  randomEffectsCount];
      siblingSiblingEffectsIndex++;
    }
    siblingSiblingEffectsIndex += randomEffectsCount;
    
    // Read Shared Spouse data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> parentParent[parentParentEffectsIndex]
		 >> parentParent[parentParentEffectsIndex
				 + randomEffectsCount];
      parentParentEffectsIndex++;
    }
    parentParentEffectsIndex += randomEffectsCount;
    
    // Read Shared Parent-Offspring data
    for(int j=0; j<randomEffectsCount; j++) {
      if(i == 0) {
	summaryLog >> variableNames[variableNamesIndex++];
      }
      summaryLog >> parentOffspring[parentOffspringEffectsIndex]
		 >> parentOffspring[parentOffspringEffectsIndex
				 + randomEffectsCount];
      parentOffspringEffectsIndex++;
    }
    parentOffspringEffectsIndex += randomEffectsCount;

    // Read locus log likelihood and convergence status
    summaryLog >> logLikelihoods[i];
    summaryLog >> logLikStatus[i];
  }

  summaryLog.close();
}

}
