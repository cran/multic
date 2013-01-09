/******************************************************************************
File: multic.cpp
$Id: multic.cpp,v 1.22 2013/01/22 22:57:20 m015733 Exp $
Description: This file defines the void function multic.  In the vast majority
             of cases, multic will be called from an Splus routine.  multic is
             a program that performs linkage analysis.  It uses variance
             component analysis to calculate likelihoods of gene markers
             causing quantitative traits.
Author: Eric Lunde, 6-13-03
Updates: (Date, Modified By, Modification Description)
7-10-03, Eric Lunde, Just after executing the uppermost Robust_Variance section
         (circa line 160) the program exited, which killed the Splus session.
         I changed that exit to a return statement so the Splus session remains
         after the program is completed.
7-11-03, Eric Lunde, I added code to read the first line of loci.out which
         contains the name of the (m)ibd file that genereated it.  This
         information is passed to the Calculate object so the (m)ibd file name
         can be added to the output file summary.Log.
7-25-03, Eric Lunde, I altered the second half of this file dramatically.  When
         building the Robust_Variance object we use other objects like PP, PS,
         and MajorGene1.  All of these latter objects read from files, but the
         files were opened and closed in multic.  I thought that this was poor
         design, so I moved the file opening and closing to within the objects
         to promote encapsulation.
7-30-03, Eric Lunde, Many of the subclasses of Model read from the file
         share.out.  Since we have already created the internal structure of
         share.out, I now have them referencing that array to save time.
7-31-03, Eric Lunde, Just like share.out, fort.12 was being read many many
         times.  I have written code to represent all of fort.12 in memory.
         What remains is to alter the sections that read from fort.12 and have
         them now read from the array.
8-29-03, Eric Lunde, I added a small amount of code to alter the screen output
         during execution.  I didn't like all the information output during
         runtime.  I thought it would be better to let the user know what the
         program was doing and where it was at for execution.  I will probably
         add things like time of execution and number of iterations before
         (non)convergence.
9-02-03, Eric Lunde, I added code to test if fort.12 contains information
         detailing which participants are probands.  I changed the code in
         readFort12 and printFortData.
Update logs are now kept in CVS.
******************************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "multic.h"
#include "Least.h"
#include "Maxfun_add.h"
#include "Calculate.h"
#include "Robust_Variance.h"
#include "Model.h"
#include "Composite.h"
#include <ctime>
#include <cstring>
#include "verS.h"
#include "Lib.h"
#include "Rostream.h"
#include "Rstreambuf.h"

using namespace Rcpp;

int ascert_ch;
int ascert_flag;

int runopt_flag;
int run_flag;
int need_ibd_flag;

double epsilon;

char fixtoboundary_ch;
int fixtoboundary_flag;

int algorithm_flag;

void readShareOut(const char *shareFileName, ShareRelation **sa, int *rc);
void readFort12(FortData **fa, int traitCount,
		int covariateCount, int repeatCount, char **famid, char **id,
		char **dadid, char **momid, Sint *sex, double **traitMatrix,
		double **covariateMatrix, Sint *proband,
		int observationsCount);
void freeFortData(FortData *fortArray, int observations);

extern "C" {

void multic(char **famid,
	    char **id,
	    char **dadid,
	    char **momid,
	    Sint *sex,
	    double *trait,
	    double *covariate,
	    Sint *proband,
	    Sint *observationsCount,
	    Sint *asc_ch,
	    Sint *runopt_f,
	    double *epsilonValue,
	    Sint *fixtoboundary_c,
	    Sint *algorithm_f,
	    Sint *traitCount,
	    Sint *covariateCount,
	    Sint *repeatCount,
	    char **traitNames,
	    double *missingValue,
	    double *initialValues,
	    char **constraints,
	    Sint *maxIterations,
	    double *initialBetas,
	    char **shareFileName,
	    Sint *printProgress,
	    Sint *calculateResiduals,
	    Sint *familySizes,
	    char **uniqueFamilies,
	    Sint *familyCount
	    )
{
  // Initialize the global variables in multicParameters.h  
  ascert_ch = (int)*asc_ch;
  ascert_flag = (ascert_ch ? YES : NO);
  
  runopt_flag = *runopt_f;
  if ((runopt_flag == 1) || (runopt_flag == 3) || (runopt_flag == 4)) {
    need_ibd_flag = YES;
  }
  if (runopt_flag == 2) {
    need_ibd_flag = NO ;
  }

  epsilon = *epsilonValue;

  fixtoboundary_ch = *fixtoboundary_c;
  fixtoboundary_flag = (fixtoboundary_ch ? YES : NO);

  algorithm_flag = *algorithm_f;

  /*
    cout << std::endl << "multic.cpp variables:" << std::endl
    << "ascert_ch = " << ascert_ch << std::endl
    << "ascert_flag = " << ascert_flag << std::endl
    << "runopt_flag = " << runopt_flag << std::endl
    << "need_ibd_flag = " << need_ibd_flag << std::endl
    << "epsilon = " << epsilon << std::endl
    << "fixtoboundary_ch = char(" << (int)fixtoboundary_ch << ")" << std::endl
    << "fixtoboundary_flag = " << fixtoboundary_flag << std::endl
    << "algorithm_flag = " << algorithm_flag << std::endl << std::endl;
  */
  
  // TraitMarkerCov_par is a class that stores the values for trait, loci,
  // covariate, and repeat values.  It also stores the names, missing values,
  // and lociations in fort.12 of those traits, loci, and covariates.
  TraitMarkerCov_par tmc(*traitCount, *covariateCount,
			 *repeatCount, traitNames, *missingValue);
  
  int itraits = tmc.gettraitnum();
  int irepeatmst = tmc.getirepeatmst();

  // InitValue_par is a class that stores the starting values for the
  // numerical algorithm.
  InitValue_par init(&tmc, initialValues, constraints, (int)*maxIterations,
		     initialBetas);

  ShareRelation *shareArray = NULL;
  int relationCount = 0;
  readShareOut(shareFileName[0], &shareArray, &relationCount);
  /*
    for(int i=0 ; i<relationCount; i++) {
    printShareRelation(&shareArray[i]);
    }
  */

  // 1) alter readFort12 function
  // 2) alter function naming and documentation
  // 3) alter multic.h
  // 4) alter multic.q and functions used to create fort.12

  double **traitMatrix =
    Lib::dVecToMat(trait, 0, *observationsCount - 1,
		   0, (*traitCount * *repeatCount) - 1);
  double **covariateMatrix =
    Lib::dVecToMat(covariate, 0, *observationsCount - 1,
		   0, (*covariateCount * *repeatCount) - 1);

  FortData *fortArray = NULL;
  readFort12(&fortArray, (int)*traitCount, (int)*covariateCount,
	     (int)*repeatCount, famid, id, dadid, momid, sex,
	     traitMatrix, covariateMatrix, proband, (int)*observationsCount);

  Lib::free_dmatrix(traitMatrix, 0, *observationsCount - 1,
		    0, (*traitCount * *repeatCount) - 1);
  Lib::free_dmatrix(covariateMatrix, 0, *observationsCount - 1,
		    0, (*covariateCount * *repeatCount) - 1);
  /*
    for(int i=0 ; i<*observationsCount; i++) {
    printFortData(&fortArray[i], (int)*traitCount, (int)*covariateCount,
    (int)*repeatCount);
    }
  */
  
  
  if(algorithm_flag == 4){
    if(irepeatmst!=1){
      PROBLEM "Sorry, Maxfun doesn't work for Longitutal data so far.\nmultic.cpp key 202\n"
	RECOVER(NULL_ENTRY);
    }
    Rcout << "enter maxfun" << std::endl;
    Maxfun_add maxfun;
    maxfun.run(&tmc, &init, fortArray, *observationsCount, shareArray,
	       relationCount);
    Rcout << "exit maxfun" << std::endl;
  }
  
  // This Calculate object will be initialized later, but it must
  // be defined here because it will be used outside the scope of
  // where it was given value.
  Calculate *c = NULL;

  // We need the name of the (m)ibd file that generated this loci.out
  // The next 10 lines were added by Eric Lunde on 7-11-03
  char ibdFileName[BUF];  
  if(runopt_flag == 1) {
    ifstream loci_out(loci_file);
    if(loci_out.fail()) {
      PROBLEM "The file %s could not be opened for reading.\nmultic.cpp key 100\n",
	loci_file RECOVER(NULL_ENTRY);
    }
    loci_out.getline(ibdFileName, BUF);
    loci_out.close();
  }else if(runopt_flag == 2) {
    strncpy(ibdFileName, "# null", BUF);
  }

  //assume least square works for Longititual data
  if(algorithm_flag != 1 && algorithm_flag != 3) {
    // The code '&ibdFileName[2]' needs some explaining.  In loci.out, the
    // first line is a pound sign followed by a space and then the name of the
    // (m)ibd file that generated it.  We want that name inside the Calculate
    // object so it can print it to file.  So the code segment (approximately
    // line 100) reads in the first line of loci.out.  We want to discard the
    // first two characters.  That is why we have '&ibdFileName[2]'.
    c = new Calculate(NULL, &ibdFileName[2], epsilon, shareArray,
		      relationCount, fortArray, *observationsCount,
		      initialValues, initialBetas, *printProgress,
		      *calculateResiduals, familySizes, uniqueFamilies,
		      (int) *familyCount);
  }else {
    if(irepeatmst!=1){
      PROBLEM "Sorry, Least Square option doesn't work for Longitutal data so far.\nmultic.cpp key 247"
	RECOVER(NULL_ENTRY);
    }
    Rcout << "Running Least Square ..." << std::endl;
    Least least(&tmc, shareArray, relationCount, fortArray,
		*observationsCount);
    least.least_main();
    Rcout << "Exit Least Square." << std::endl;
    
    if(itraits!=1 && algorithm_flag==1){
      // There should be an error message printed here, but as of 4-25-03,
      // Eric Lunde (284-5630) does not know why this is an error.  Once he
      // does, he will add an appropriate message.
      PROBLEM "Unknown problem occured.\nmultic.cpp key 260"
	RECOVER(NULL_ENTRY);
    }

    // algorithm_flag == 1 means least square estimation method
    if(itraits == 1 && algorithm_flag == 1) {
      //      cout << "1) Entering Robust Variance Component." << std::endl;
      Composite composite;
      
      MajorGene1 m1(&least);
      
      PolyGene p1(&least, shareArray);

      Environment p2(&least);

      composite.add(&m1);
      composite.add(&p1);
      composite.add(&p2);

      ofstream outfp("least.log", ios::app);
      if(outfp == NULL) {
	PROBLEM "The file least.log could not be opened for appending.\nmultic.cpp key 281"
	  RECOVER(NULL_ENTRY);
      }

      Robust_Variance rv;
      rv.buildModel(&composite, &outfp, &tmc, fortArray, *observationsCount);
      outfp.close();

      freeFortData(fortArray, (int) *observationsCount);
      free(shareArray);
      
      //      cout << "Exiting Robust Variance Component." << std::endl;
      return;
    }
    
    // The code '&ibdFileName[2]' needs some explaining.  In loci.out, the
    // first line is a pound sign followed by a space and then the name of the
    // (m)ibd file that genereated it.  We want that name inside the Calculate
    // object so it can print it to file.  So the code segment (approximately
    // line 100) reads in the first line of loci.out.  We want to discard the
    // first two characters.  That is why we have '&ibdFileName[2]'.
    c = new Calculate(&least, &ibdFileName[2], epsilon, shareArray,
		      relationCount, fortArray, *observationsCount,
		      initialValues, initialBetas, *printProgress,
		      *calculateResiduals, familySizes, uniqueFamilies,
		      (int) *familyCount);
  }

  //  cout << "Entering multic." << std::endl;

  // Calculate current time to display to user
  struct tm *timePtr;
  time_t localTime;

  localTime = time(NULL);
  timePtr = localtime(&localTime);

  if(runopt_flag == 1) {
    // This first newline is in response to the format of outputting the
    // iteration number.  Since after printing the iteration number, I didn't
    // print a newline, we need to print one here before we start the next 
    // ibd file.  Eric Lunde 01-08-04
    if(*printProgress) {
      Rcout << "Calculating likelihoods for locus: '" << &ibdFileName[2]
	   << "' on " << asctime(timePtr);
    }
  }else if(runopt_flag == 2) {
    if(*printProgress) {
      Rcout << "Calculating likelihoods under null hypothesis on "
	   << asctime(timePtr);
    }
  }

  c->main_fun(&tmc, &init);

  if(irepeatmst==1 && itraits==1){
    Composite composite;

    MajorGene1 *m1 = NULL;
    MajorGene2 *m2 = NULL;
    if(runopt_flag != 2) {
      m1 = new MajorGene1(c);
      if(init.getm1flag() != FIXED) {
	composite.add(m1);
      }
      
      m2 = new MajorGene2(c);
      if(init.getm2flag() != FIXED) {
	composite.add(m2);
      }
    }      

    PolyGene p1(c, shareArray);
    if(init.getsflag() != FIXED) {
      composite.add(&p1);
    }

    Environment p2(c);
    if(init.gettflag()!=FIXED){
      composite.add(&p2);
    }

    SS ss(c, shareArray);
    if(init.getcflag()!=FIXED){
      composite.add(&ss);
    }

    PP pp(c, shareArray);
    if(init.getpflag()!=FIXED){
      composite.add(&pp);
    }

    PS ps(c, shareArray);
    if(init.getqflag()!=FIXED){
      composite.add(&ps);
    }

    ofstream outfp("multic.log", ios::app);
    if(outfp == NULL) {
      PROBLEM "The file multic.log could not be opened for appending.\nmultic.cpp key 246\n"
	RECOVER(NULL_ENTRY);
    }

    Robust_Variance rv;
    rv.buildModel(&composite, &outfp, &tmc, fortArray, *observationsCount);
    outfp.close();

    if(runopt_flag != 2) {
      delete m1;
      delete m2;
    }
  }

  delete c;

  freeFortData(fortArray, (int) *observationsCount);
  free(shareArray);
}

}

/*****************************
File name: multic.cpp
Function name: readShareOut
Description: Create an array of ShareRelation structs.  Fill it with all the
             data from share.out.
Input: const char *shareFileName - the name of the file used to provide the
                                   share data.  Most likely this will be 
                                   'share.out', but the directory may
                                   something other than the current directory.
       ShareRelation **sa - return the created array in this pointer.  It
                            stands for share array.
       int *rc - return the size of the array in this pointer.  It stands for
                 relation count.
Output: NONE.
Side Effects: NONE.
Author: Eric Lunde, 7-29-03
******************************/
void readShareOut(const char *shareFileName, ShareRelation **sa, int *rc) {
  // Create the internal storage space for representing share.out in main
  // memory.  First read the number of lines (number of slots for our array).
  // Second read in the actual data.  relationCount begins at -1 because
  // fp_share.good() is false after reading the very last line of the file
  // which contains no data (""). - Eric Lunde, 7/29/03
  ifstream fp_share(shareFileName);
  if(fp_share.fail()) {
    PROBLEM "The file %s could not be opened for reading.\nmultic.cpp key 251\n",
      shareFileName RECOVER(NULL_ENTRY);
  }
  int relationCount = -1;
  char line[BUF];
  while(fp_share.good()) {
    fp_share.getline(line, BUF);
    relationCount++;
  }

  // Reset our file pointer so we can read the real data. - Eric Lunde, 7/29/03
  fp_share.clear();
  fp_share.seekg(0, ios::beg);

  // Create the storage space
  ShareRelation *shareArray;
  shareArray = (ShareRelation *) malloc( relationCount
					 * sizeof(ShareRelation) );
  if(!shareArray) {
    PROBLEM "Error during memory allocation for the shareArray.\nmultic.cpp key 293\n"
      RECOVER(NULL_ENTRY);
  }

  for(int i=0; i<relationCount; i++) {
    fp_share >> shareArray[i].seqId1 >> shareArray[i].seqId2
	     >> shareArray[i].geneticSimilarity >> shareArray[i].areSiblings
	     >> shareArray[i].areSpouses >> shareArray[i].areParentChild;
  }

  fp_share.close();

  *rc = relationCount;
  *sa = shareArray;
}

void printShareRelation(ShareRelation *sr) {
  Rcout << sr->seqId1 << '\t' << sr->seqId2 << '\t' // << setw(5)
       << sr->geneticSimilarity << '\t' << sr->areSiblings << '\t'
       << sr->areSpouses << '\t' << sr->areParentChild << std::endl;
}

/*****************************
File name: multic.cpp
Function name: readFort12
Description: Create an array of FortData structs.  Fill it with all the data
             from fort.12.
Input: FortData **fa - return the created array in this pointer.  It stands for
                       fort array.
       int traitCount - the number of traits to be analyzed
       int covariateCount - the number of covariates to be analyzed
       int repeatCount - the number of times the longitudinal measurements 
                         (this number is 1 when we do a multivariate analysis)
Output: NONE.
Side Effects: NONE.
Author: Eric Lunde, 7-30-03
******************************/
void readFort12(FortData **fa, int traitCount,
		int covariateCount, int repeatCount, char **famid, char **id,
		char **dadid, char **momid, Sint *sex, double **traitMatrix,
		double **covariateMatrix, Sint *proband,
		int observationsCount)
{
  // Create the internal storage space for representing fort.12 in main
  // memory.  First, create the space (observationsCount).
  // Second read in the actual data.  
  // - Eric Lunde, 12/16/04
  
  // Create the storage space
  FortData *fortArray;
  fortArray = (FortData *) malloc( observationsCount * sizeof(FortData) );
  if(!fortArray) {
    PROBLEM "Error during memory allocation for the fortArray.\nmultic.cpp key 350\n"
      RECOVER(NULL_ENTRY);
  }  

  int j;
  for(int i = 0; i < observationsCount; i++) {
    fortArray[i].familyId = (char *) malloc((strlen(famid[i]) + 1)
					    * sizeof(char));
    strcpy(fortArray[i].familyId, famid[i]);

    fortArray[i].seqId = (char *) malloc((strlen(id[i]) + 1)
					 * sizeof(char));
    strcpy(fortArray[i].seqId, id[i]);

    fortArray[i].fatherId = (char *) malloc((strlen(dadid[i]) + 1)
					    * sizeof(char));
    strcpy(fortArray[i].fatherId, dadid[i]);

    fortArray[i].motherId = (char *) malloc((strlen(momid[i]) + 1)
					    * sizeof(char));
    strcpy(fortArray[i].motherId, momid[i]);

    fortArray[i].sex = (int) sex[i];
   
    fortArray[i].traits = (double *) malloc(traitCount * repeatCount
					    * sizeof(double));
    for(j=0; j<traitCount * repeatCount; j++) {      
      fortArray[i].traits[j] = traitMatrix[i][j];      
    }

    fortArray[i].covariates = (double *) malloc(covariateCount * repeatCount
						* sizeof(double));
    for(j=0; j<covariateCount * repeatCount; j++) {
      fortArray[i].covariates[j] = covariateMatrix[i][j];
    }

    if(ascert_flag == YES) {      
      fortArray[i].proband = (int) proband[i];
    }
  }
    

  *fa = fortArray;
}

void printFortData(FortData *fd, int traitCount, int covariateCount,
		   int repeatCount) {
  Rcout << fd->familyId << ' ' << fd->seqId << ' ' << fd->fatherId << ' '
       << fd->motherId << ' ' << fd->sex;
  for(int j=0; j<traitCount * repeatCount; j++) {
    Rcout << ' ' << fd->traits[j];
  }
  for(int j=0; j<covariateCount * repeatCount; j++) {
    Rcout << ' ' << fd->covariates[j];
  }
  if(ascert_flag == YES) {
    Rcout << ' ' << fd->proband;
  }
  Rcout << std::endl;
}

void freeFortData(FortData *fortArray, int observations) {
  for(int i = 0; i < observations; i++) {
    free(fortArray[i].familyId);
    free(fortArray[i].seqId);
    free(fortArray[i].fatherId);
    free(fortArray[i].motherId);
    free(fortArray[i].traits);
    free(fortArray[i].covariates);
  }
  free(fortArray);
}
