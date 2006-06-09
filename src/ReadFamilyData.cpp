/******************************************************************************
File: ReadFamilyData.cpp
Description: ReadFamilyData is a class that provides methods that access the
             data contained in fort.12.  
Author: Eric Lunde, 8-6-03 (This file was not writen by Eric, but this is the
        date he took over and began modifying it.)
Updates: (Date, Modified By, Modification Description)
8-6-03, Eric Lunde, Originally, ReadFamilyData inherited from Getheader.  The
        main purpose of Getheader was to test for the presence of a header on
        the file fort.12.  A header was the number of families followed by the
        size of each family.  The header came before the real data in fort.12.
        Since multic.s (Splus) creates fort.12 for us now and it does not
        genereate a header, Getheader was no longeer necessary.  I decided to
        remove it for simplifing the class heirarchy.  Now ReadFamilyData
        contains a pointer to a FortData struct (defined in multic.h).  Each
        slot of this array contains all the information from one line of the
        file fort.12.  This way, when we want to read from the file, we can
        just access the array.  I have deleted all of the code referring to the
        file fort.12.  Via the getdata method, we read the file.  Some code
        sections used this getdata method and others just read from the file
        directly.  Because of this, I've had to implement the array to have
        both ablilities.  I keep an internal int dataIndex which keeps track
        of where in the array we are.  There is a related resetDataIndex method
        to reset that variable back to 0 (akin to reseting the file get
        pointer).  For the sections that read from the file itself, I've
        written a method to return the array, so the calling code can read the
        array just like it had a handle to the file.
******************************************************************************/
#include "ReadFamilyData.h"
#include "multic.h"
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

ReadFamilyData::ReadFamilyData(TraitMarkerCov_par *tmc, int col,
			       Family *fam1, FortData *fa, int ds)
{
  fam = fam1;

  itraits = tmc->gettraitnum();
  iloci = tmc->getmarkernum();
  icovs = tmc->getcovnum();
  iinitvcnum = tmc->getinitvcnum();
  env_vcnum = tmc->getinitenvvcnum();
  total_trait_values = tmc->gettotal_trait_values();
  total_cov_values = tmc->gettotal_cov_values();

  trait_array = tmc->gettrait_array();
  marker_array = tmc->getmarker_array();
  cov_array = tmc->getcov_array();
  t_start = tmc->get_t_start();
  t_len = tmc->get_t_len();
  m_start = tmc->get_m_start();
  m_len = tmc->get_m_len();   
  c_start = tmc->get_c_start();
  c_len = tmc->get_c_len();

  fortArray = fa;
  dataSize = ds;
  dataIndex = 0;
  
  // Test to see if there is at least 1 family.
  if(dataSize <= 0) {
    cerr << "fort.12 contains no family data" << endl
	 << "ReadFamilyData.cpp line 37" << endl;
  }else {
    familyCount = 1;
    for(int i=1; i<dataSize; i++) {
      if(strcmp(fortArray[i-1].familyId, fortArray[i].familyId) != 0) {
	familyCount++;
      }
    }
    familySizes = (int *) malloc( familyCount * sizeof(int) );
    
    // We konw that there has to be at least one person in each family.
    for(int i=0; i<familyCount; i++) {
      familySizes[i] = 1;
    }
    
    int currentFamily = 0;
    for(int i=1; i<dataSize; i++) {
      if(strcmp(fortArray[i-1].familyId, fortArray[i].familyId) == 0) {
	familySizes[currentFamily]++;
      }else {
	currentFamily++;
      }
    }
  }
}

ReadFamilyData::~ReadFamilyData() {
  free(familySizes);
}

void ReadFamilyData::getdata(int n_fam_mem){
  int  i,j;

  for(i=0;i<n_fam_mem;i++){
    fam->person[i].missing_traitflag=0;
  }

  for ( j = 0; j < n_fam_mem; j++, dataIndex++) {
    for ( i = 0; i < total_trait_values; i++) {
      fam->person[j].trait[i] = fortArray[dataIndex].traits[i];
      if (fam->person[j].trait[i] == trait_array[i].t_misval) {
	fam->person[j].missing_traitflag = 1;
      }
    }

    for ( i = 0; i < total_cov_values; i++) {
      fam->person[j].cov[i] = fortArray[dataIndex].covariates[i];
    }
  }
}

int ReadFamilyData::gettrait(){
  return itraits;
}

int ReadFamilyData::getmarker(){
  return iloci;
}
                                
int ReadFamilyData::getcov(){
  return icovs;
}
                                
int ReadFamilyData::gett_start(){
  return t_start[0];
}

int ReadFamilyData::getFamilyCount(){
  return familyCount;
}

int *ReadFamilyData::getFamilySizes(){
  return familySizes;
}

void ReadFamilyData::resetDataIndex() {
  dataIndex = 0;
}

FortData *ReadFamilyData::getFortArray() {
  return fortArray;
}
