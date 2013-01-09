/******************************************************************************
File: InitValue_par.cpp
Description: This class originally was responsible for parsing part of the
             'multic.par' file.  It obtained the initial values for the
             variance component analysis.  Those values were stored in the
             arrays with names like init_x_M1.  As of writing this
             documentation, it still receives and stores the initial values,
             but it receives them via function parameters generated in
             multic.s (Splus proram).  In addition to storing the initial
             values, it also stores the 'precisions'.  These precisions tell
             the program how to treat each of the initial parameters.  It can
             treat them as values to begin estimation, as fixed, or as
             constrained.
Author: Eric Lunde, 6-12-03
Updates: (Date, Modified By, Modification Description)
******************************************************************************/
#include "InitValue_par.h"
#include "multic.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <iostream>
#include <fstream>

/*****************************
Class name: InitValue_par
Method name: InitValue_par (constructor)
Description: Construct an instance of an InitValue_par object and initialize
             some of the member variables.
Input: TraitMarkerCov_par *in2 - This object holds information such as number
                                 of traits, loci, and covariates.  We need
                                 this object to get at these values.
       double *initialValues - An array of doubles that contain all the
                               initial values for the numerical algorithm.  It
                               also contains the maximum number of iterations
                               for the algorithm.  Finally, if there are
			       covariates,  it contains the covariates
			       coeffecient.  This will be passed to
                               GetInitValues without being modified here.
       char **precisions - An array of strings containing how the initial
                           values are going to be used.  This will be passed
                           to GetInitValues without being modified here.
Output: NONE.
Side Effects: NONE.
Author: Eric Lunde, 6-12-03
******************************/
InitValue_par::InitValue_par(TraitMarkerCov_par *in2, double *initialValues,
			     char **precisions, int maxIterations,
			     double *initialBetas){
  itraits = in2->gettraitnum();
  iloci  = in2->getmarkernum();
  icovs  = in2->getcovnum();
  irepeatmst = in2->getirepeatmst();
  iinitvcnum = in2->getinitvcnum();
  env_vcnum  = in2->getinitenvvcnum();
  int count = 0;
  for ( count = 0; count < iinitvcnum; count++ ) {
    init_x_S[count] =0.0;
    init_x_M1[count]=0.0;
    init_x_M2[count]=0.0;
    init_x_c[count] =0.0;
    init_x_p[count] =0.0;
    init_x_q[count] =0.0;
  }
  for(count=0;count<env_vcnum;count++){
    init_x_T[count] =0.0;
  }
  GetInitValues(initialValues, precisions);

  N1 = maxIterations;

  if(icovs > 0) {
    int betasIndex = 0;
    for(int i=0; i<itraits; i++) {
      for(int j=0; j<icovs; j++) {
	init_x_beta[i][j] = initialBetas[betasIndex++];
      }
    }
  }

}

InitValue_par::~InitValue_par() {
}

/* Get the initial values of init_x_mu, init_x_S, init_x_M1, init_x_M2, init_x_T,
   init_x_c, init_x_p, init_x_q   from the strings.  */
void InitValue_par::GetXValues(int numval, double *x_val, char *x_str){
  char str[STRING], temp_str[STRING];
         
  int j  = 0;
  int kv = 0;   
  int len = strlen(x_str);

  for (int i = 0; i < len; i++) {
    // for the first itraits-1  value(s) of x_mu. 
    if (kv < numval-1) {
      if (x_str[i] != ',') {
	str[j]   = x_str[i];
	str[j+1] = '\0';
	j++;
      }else {
	strcpy(temp_str, str);
	x_val[kv] = atof(temp_str);
	j = 0;
	kv++;
      }
    }
    /* for the last value of x_mu.  */
    else {
      str[j] = x_str[i];
      str[j+1] = '\0';
      j++;
    }
  }
  strncpy(temp_str, str, STRING-1);
  x_val[kv] = atof(temp_str);
}


/* Get the initial values of covariates 'init_x_beta' from the string.  */
void InitValue_par::GetInitCovValues(char x_str[]){
  int    i,j,kv;
  int    len;
  double *cov_val;
  char   str[STRING], temp_str[STRING];
         
        
  /* Initialize the size of the pointers.    5-12-98; E.Y.  */
  cov_val = (double *) malloc((itraits*icovs) * sizeof(double));

  j  = 0;
  kv = 0;
  len = strlen(x_str);

  for (i = 0; i < len; i++) {
    /* for the first itraits*icovs  value(s) of beta.  */
    if ( kv < (itraits*icovs -1) ) {
      if (x_str[i] != ',') {
	str[j]   = x_str[i];
	str[j+1] = '\0';
	j++;
      }   
      else {   
	strcpy(temp_str, str);
	cov_val[kv] = atof(temp_str);
	j = 0;
	kv++;
      }
    }
    /* for the last value of x_cov.  */
    else {
      str[j] = x_str[i];
      str[j+1] = '\0';
      j++;
    }
  }
  strcpy(temp_str, str);
  cov_val[kv] = atof(temp_str);  
               
  /* put the values of cov_val[] to array beta[][].     */
  kv = 0;    
  for (i = 0; i < itraits; i++) {
    for (j = 0; j < icovs; j++) {
      init_x_beta[i][j] = cov_val[kv];
      kv++;
    }
  }
  free(cov_val);
}

/*****************************
Class name: InitValue_par
Method name: GetInitValues
Description: Given an array of doubles and strings, set the appropriate
             internal init_x_* arrays to the initial values calculated from
             Splus (or the multic.par file) and assign the precisions (or
             variable usage) to the init_*_flag variables.
Input: double *initialValues - An array of doubles that contain all the
                               initial values for the numerical algorithm.  It
                               also contains the maximum number of iterations
                               for the algorithm.  Finally, if there are
			       covariates,  it contains the covariates
			       coeffecient.
       char **precisions - An array of strings containing how the initial
                           values are going to be used.
Output: NONE.
Side Effects: NONE.
Other: Logically, initialValues is a matrix of doubles.  But since it would be
       ragged, it is stored as a single vector of doubles.
Author: Eric Lunde, 6-12-03
******************************/
void InitValue_par::GetInitValues(double *initialValues, char **precisions)
{
  int valuesIndex = 0;

  init_mu_flag = convertPrecision(precisions[0]);
  for(int i=0; i<itraits; i++) {
    init_x_mu[i] = initialValues[valuesIndex++];
    //    cout << "init_x_mu[" << i << "]  = " << init_x_mu[i] << std::endl;
  }

  init_s_flag = convertPrecision(precisions[1]);
  for(int i=0; i<iinitvcnum; i++) {
    init_x_S[i] = initialValues[valuesIndex++];
    //    cout << "init_x_S[" << i << "]  = " << init_x_S[i] << std::endl;
  }

  init_m1_flag = convertPrecision(precisions[2]);
  for(int i=0; i<iinitvcnum; i++) {
    init_x_M1[i] = initialValues[valuesIndex++];
    //    cout << "init_x_M1[" << i << "]  = " << init_x_M1[i] << std::endl;
  }

  init_m2_flag = convertPrecision(precisions[3]);
  for(int i=0; i<iinitvcnum; i++) {
    init_x_M2[i] = initialValues[valuesIndex++];
    //    cout << "init_x_M2[" << i << "]  = " << init_x_M2[i] << std::endl;
  }

  init_t_flag = convertPrecision(precisions[4]);
  for(int i=0; i<env_vcnum; i++) {
    init_x_T[i] = initialValues[valuesIndex++];
    //    cout << "init_x_T[" << i << "]  = " << init_x_T[i] << std::endl;
  }

  init_c_flag = convertPrecision(precisions[5]);
  for(int i=0; i<iinitvcnum; i++) {
    init_x_c[i] = initialValues[valuesIndex++];
    //    cout << "init_x_c[" << i << "]  = " << init_x_c[i] << std::endl;
  }

  init_p_flag = convertPrecision(precisions[6]);
  for(int i=0; i<iinitvcnum; i++) {
    init_x_p[i] = initialValues[valuesIndex++];
    //    cout << "init_x_p[" << i << "]  = " << init_x_p[i] << std::endl;
  }

  init_q_flag = convertPrecision(precisions[7]);
  for(int i=0; i<iinitvcnum; i++) {
    init_x_q[i] = initialValues[valuesIndex++];
    //    cout << "init_x_q[" << i << "]  = " << init_x_q[i] << std::endl;
  }
}

/*****************************
Class name: InitValue_par
Method name: convertPrecision
Description: Given a string representation of how an initial value is to be
             used, return the appropriate interger representation.
Input: char *precision - A string specifying how an initial value will be
                         fixed, estimated, or constrained.
Output: int - This int is the integer representation of a variable's desired
              use.
Side Effects: If precision is not of the correct form, an error message is
              printed to standard error and the program exits with an error
              level of 1.
Author: Eric Lunde, 6-12-03
******************************/
int InitValue_par::convertPrecision(char *precision) {
  if(strcmp(precision, "F")==0 || strcmp(precision, "f")==0 ) {
    return FIXED;
  }else if(strcmp(precision, "E")==0 || strcmp(precision, "e")==0 ) {
    return ESTIMATE;
  }else if(strcmp(precision, "C")==0 || strcmp(precision, "c")==0 ) {
    return CONSTRAINT;
  }else {
    PROBLEM "%s%s%s%s%s%s%s%s",
      "Error in InitValue_par::GetInitValues.\n",
      "Precision values must be F, f, E, e, C, or c.\n",
      "Your precision was ", precision, "\n",
      "If you cannot fix the problem, contact programmer ",
      "Eric Lunde @ 507-284-5630",
      "InitValue_par.cpp line 277\n"
      RECOVER(NULL_ENTRY);
    // Adding a return -1; makes the compiler happy.
    return -1;
  }
}

int InitValue_par::getmuflag(){
  return init_mu_flag;
}

int InitValue_par::getsflag(){
  return init_s_flag;
}
            
int InitValue_par::getm1flag(){
  return init_m1_flag;
}
            
int InitValue_par::getm2flag(){
  return init_m2_flag;
}
            
int InitValue_par::gettflag(){
  return init_t_flag;
}
            
int InitValue_par::getcflag(){
  return init_c_flag;
}
            
int InitValue_par::getpflag(){
  return init_p_flag;
}

int InitValue_par::getqflag(){
  return init_q_flag;
}
            
int InitValue_par::getiterationnum(){
  return N1;
}

double* InitValue_par::getinitmu(){
  return init_x_mu;
}

double* InitValue_par::getinitS(){
  return init_x_S;
}
         
double* InitValue_par::getinitM1(){
  return init_x_M1;
}
         
double* InitValue_par::getinitM2(){
  return init_x_M2;
}
         
double* InitValue_par::getinitT(){
  return init_x_T;
}
         
double* InitValue_par::getinitc(){
  return init_x_c;
}
         
double* InitValue_par::getinitp(){
  return init_x_p;
}
         
double* InitValue_par::getinitq(){
  return init_x_q;
}
         

void InitValue_par::getinitbeta(double beta[MAXNUMTRAIT][MAXNUMCOV]){
  for(int i=0;i<itraits;i++)
    for(int j=0;j<icovs;j++){
      beta[i][j] = init_x_beta[i][j];
    }
}
