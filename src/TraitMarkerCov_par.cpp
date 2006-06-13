#include "TraitMarkerCov_par.h"
#include "Lib.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <S.h>
using namespace std;

TraitMarkerCov_par::TraitMarkerCov_par(long traitCount,
				       long covariateCount,
				       long repeatCount,
				       char **trait_names,
				       double missing_value
				       //,long *format
				       )
{
  itraits = traitCount;
  iloci = 1;
  icovs = covariateCount;
  irepeatmst = repeatCount;
  traitNames = trait_names;
  missingValue = missing_value;
  GetParameters();
}

TraitMarkerCov_par::~TraitMarkerCov_par() {
}


/*
Get the position(s) of trait(s), the position of marker and marker's type
from parameter file --- 'multic.par' or 'longic.par'.
This function also calls functions: 
    TraitLociCov,
    TraitMarkerCovParameter,
    GetFromFormat, 
*/
void TraitMarkerCov_par::GetParameters(/*long *format*/) {
  int  flag;             /* flag of trait or marker.            */

  /* Get initial G_VAL 'itraits', 'iinitvcnum', 'iloci', 'icovs', 
     'imubetacovnum',  'env_vcnum'  values here and NO CHANGE.   */ 
  TraitLociCov();        /* Type of Analysis record.             */

  int currentTraitMarkerCovariateIndex = 0;

  /* Trait dependent parameters.          */
  if (itraits >0) { 
    flag = TRAIT;
    TraitMarkerCovParameter(flag, &currentTraitMarkerCovariateIndex);
  }

  /* Marker dependent parameters.         */
  if (iloci > 0) {
    flag = MARKER;
    TraitMarkerCovParameter(flag, &currentTraitMarkerCovariateIndex);
  }
  /* Covariates.                          */
  if (icovs > 0) {
    flag = COV;
    TraitMarkerCovParameter(flag, &currentTraitMarkerCovariateIndex);
  }

  /* Get initial G_VAL 't_start', 't_len', 'm_start', 'm_len',
     'c_start', 'c_len' and 'marker_type' values here and NO CHANGE.  */

  //    GetFromFormat(format);            // Format For Family Data File Records.
}

/*
Get the name and missing value for each trait, marker or covariate from
parameter file 'multic.par' or 'longic.par'. 
Take care of missing values part.
*/
void
TraitMarkerCov_par::TraitMarkerCovParameter(int flag,
					    int *currentTraitMarkerCovariateIndex)
{
  char *str;
  char temp[STRING];
   
  int r = 0;
  int  count_t_m_c = -1;
  if (flag == TRAIT) {
    count_t_m_c = total_trait_values;
  }else if (flag == MARKER) {
    count_t_m_c = iloci;
  }else if (flag == COV) {
    count_t_m_c = total_cov_values;
  }
  
  for(int i = 0; i < count_t_m_c; i++, r++, (*currentTraitMarkerCovariateIndex)++) {
    // get the name from Trait, Marker or Covariate Dependent Parameters
    str = traitNames[*currentTraitMarkerCovariateIndex];

    int k = 0;
    int length = strlen(str);
    for ( int j= 0; j < length ;j++ ) {
      if (!isspace(str[j]) ) {
	if (isupper(str[j]))
	  str[j] = tolower(str[j]);
	temp[k] = str[j]; 
	temp[k+1] = '\0';
	k++;
      }
    }
    if (flag == TRAIT) {
      strcpy(trait_array[r].t_name,temp);
    }else if (flag == MARKER) {
      strcpy(marker_array[r].m_name,temp);
    }else if (flag == COV) {
      strcpy(cov_array[r].c_name,temp);
    }

    // get the missing value from Trait, Marker or Covariate Dependent Parameters
    if (flag == TRAIT) {
      trait_array[r].t_misval = missingValue;
    }else if (flag == MARKER) {
      sprintf(temp, "%f", missingValue);
      strcpy(marker_array[r].m_misval, temp);
    }else if (flag == COV) {
      cov_array[r].c_misval = missingValue;
    }
  }
}

/*
Get the format for family data file records.
There are two lines reserved for it.
This module calls   
	GetMarkerType  
function.
Take care of position & length part. (..,..,..,.....)
*/
void TraitMarkerCov_par::GetFromFormat(long *format){
  // store strings between separator(',') from Format For Family Data File Records
  // in file 'multic.par'.
  //  char s[NUMOFSUM][STRING];

  int i = 0;

  int m, n, p;
  m = n = p = 0;
  int num_items = total_trait_values + iloci + total_cov_values; 

  // count traits, marker loci and covariates
  int count_items   = 0;

  // In this for loop we have i+=2.  This is because we take the information
  // in pairs of column and length.
  for (i = 0; count_items < num_items; i+=2) {
    if (count_items < total_trait_values) {
      t_start[m] = format[i];
      t_len[m]   = format[i+1];
      m++;
      count_items++;
    }else if ( count_items < (total_trait_values+iloci) ) {
      m_start[n] = format[i];
      m_len[n]   = format[i+1];
      GetMarkerType(m_len[n], marker_type, n);
      n++;
      count_items++;
    }else {
      c_start[p] = format[i];
      c_len[p]   = format[i+1];
      p++;
      count_items++;
    }
    // }// closes the biggest if block
  }// closes the for loop
}

/* get marker type from the format.       */
void TraitMarkerCov_par::GetMarkerType(int m_len1, int *marker_type1, int n){
  switch (m_len1){
  case 2:
  case 3:
    marker_type1[n] = PHENOTYPE;
    break;
  default:
    marker_type1[n] = GENOTYPE;
    break;
  }
}

/*
Get the Type of Analysis recode from parameter file 'multic.par' or
'longic.par'.
    itraits     number of traits.
    iloci       number of marker loci.
    icovs       number of covariates. 
    irepeatmst      number of time-frame.  (for longitudinal data only)
Take care of third line of the parameter file.
*/
void TraitMarkerCov_par::TraitLociCov(){

  //  *fp_para >> itraits >> iloci >> icovs >> irepeatmst;

  // Read the rest of the line (even if it is only the newline char).
  //  fp_para->getline(str, BUF);

  if(irepeatmst==1) {
    datatype  = MULTIVARIATE;
  }else if(irepeatmst>1) {
    datatype  = LONGITUDINAL;
  }else {
    // This else clause was added by Eric Lunde 4-23-03
    PROBLEM "%s%s%d%s%s%s%s",
      "Multivariate/Longitudinal paramater was not correctly\n",
      "specified.  Your entry was ", irepeatmst, ".\n",
      "Please make sure this parameter is a positive value.\n",
      "TraitMarkerCov_par.cpp key 316\n",
      "Programmer Contact: Eric Lunde 284-5630\n"
      RECOVER(NULL_ENTRY);
  }

  /* assign the env_vcnum value here according the datatype.
     added 9/29/97.                                         */
  /*  changed by chen, 10/1/99 to remove one par file. */
  if (datatype == MULTIVARIATE) {
    iinitvcnum   = (itraits*(itraits+1))/2;
    env_vcnum    = iinitvcnum;
    total_trait_values = itraits; 
    total_cov_values   = icovs; 
  }else if (datatype == LONGITUDINAL) {
    iinitvcnum   = ((itraits*irepeatmst)*((itraits*irepeatmst)+1))/2;
    env_vcnum    = iinitvcnum;
    //env_vcnum = itraits*irepeatmst;
    total_trait_values = itraits*irepeatmst; 
    total_cov_values   = icovs*irepeatmst; 
  }
}


int TraitMarkerCov_par::gettraitnum(){
  return itraits;
}

int TraitMarkerCov_par::getmarkernum(){
  return iloci;
}

int TraitMarkerCov_par::getcovnum(){
  return icovs;
}

int TraitMarkerCov_par::getinitvcnum(){
  return iinitvcnum;
}
    
int TraitMarkerCov_par::getinitenvvcnum(){
  return env_vcnum;
}
    
int TraitMarkerCov_par::getirepeatmst(){
  return irepeatmst;
}

int TraitMarkerCov_par::gettotal_trait_values(){
  return total_trait_values;
}

int TraitMarkerCov_par::gettotal_cov_values(){
  return total_cov_values;
}

Trait* TraitMarkerCov_par::gettrait_array(){
  return trait_array;
}

Marker* TraitMarkerCov_par::getmarker_array(){
  return marker_array;
}

Cov* TraitMarkerCov_par::getcov_array(){
  return cov_array;
}

int* TraitMarkerCov_par::get_t_start(){
  return t_start;
}

int* TraitMarkerCov_par::get_t_len(){
  return t_len;
}
        
int* TraitMarkerCov_par::get_m_start(){
  return m_start;
}
        
int* TraitMarkerCov_par::get_m_len(){
  return m_len;
}
        
int* TraitMarkerCov_par::get_c_start(){
  return c_start;
}
        
int* TraitMarkerCov_par::get_c_len(){
  return c_len;
}
