/******************************************************************************
File: Model.cpp
Description: Model.cpp defines many of the utility classes used in the multic
             program.  The MajorGene series is used to represent data
             corresponding to the gene that we are analyzing.  The class PP
             represents a parent-parent (or spouse) relationship.  SS
             represents a sibling-sibling (often called sib-sib) relationship,
             and PS represents a parent-sibling (or parent-offspring)
             relationship.
Author: Eric Lunde, 7-30-03
Updates: (Date, Modified By, Modification Description)
7-30-03, Eric Lunde, Many of the subclasses of Model (SS, PP, PO) have
         previously read from the file share.out.  Since much of the work
         recently has been storing share.out in main memory, I now have them
         referencing that array to save time.
******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <S.h>
#include "Model.h"
#include "Lib.h"
#include "Calculate.h"
using namespace std;

MajorGene1::MajorGene1(Estimator *est){
  char restOfLine[BUF];
  estimator = est;
  fp_loci = new ifstream(loci_file);
  if(fp_loci->fail()) {
    PROBLEM "The file %s could not be opened for reading.\nModel.cpp key 15\n",
      loci_file RECOVER(NULL_ENTRY);
  }
  // This next line was added on 7-11-03 by Eric Lunde to read in the new first
  // line of the file, which is now '# ' followed by the name of the (m)ibd
  // file that generated it.
  fp_loci->getline(restOfLine, BUF);
  sigma = estimator->getmajorgene1();
  char str[] = "MajorGene1";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}

MajorGene1::~MajorGene1() {
  fp_loci->close();
  delete fp_loci;
  free(name);
}

double MajorGene1::getbeta_hat(){
  return sigma;
}

void MajorGene1::getv(double** v, int n){
  double temp1, temp2, temp3;
  char restOfLine[BUF];
  for(int i=0;i<n;i++){ 
    v[i][i] = 1.0;
    for(int j=i+1;j<n;j++){
      *fp_loci >> temp1 >> temp2 >> temp3;
      fp_loci->getline(restOfLine, BUF);
      if(temp3 == non_value) {
	temp3 = 0.5;
      }
      v[i][j] = temp3;
      v[j][i] = temp3;
    }
  }
}

char *MajorGene1::getname(){
  return name;
}

MajorGene2::MajorGene2(Estimator *est){
  char restOfLine[BUF];
  estimator = est;
  fp_loci = new ifstream(loci_file);
  if(fp_loci->fail()) {
    PROBLEM "The file %s could not be opened for reading.\nModel.cpp key 61\n",
      loci_file RECOVER(NULL_ENTRY);
  }
  // This next line was added on 7-11-03 by Eric Lunde to read in the new first
  // line of the file, which is now '# ' followed by the name of the (m)ibd
  // file that generated it.
  fp_loci->getline(restOfLine, BUF);
  sigma = estimator->getmajorgene2();
  char str[] = "MajorGene2";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}

MajorGene2::~MajorGene2() {
  fp_loci->close();
  delete fp_loci;
  free(name);
}

double MajorGene2::getbeta_hat(){
  return sigma;
}
         
void MajorGene2::getv(double** v, int n){
  double temp1, temp2, temp3, temp4, temp5, temp6;
  char restOfLine[BUF];
  for(int i=0;i<n;i++){
    v[i][i] = 1.0;
    for(int j=i+1;j<n;j++){
      *fp_loci >> temp1 >> temp2 >> temp3 >> temp4 >> temp5 >> temp6;
      fp_loci->getline(restOfLine, BUF);
      if(temp3 == non_value) {
	temp3 = 0.5;
      }
      if(temp6 == non_value) {
	temp6 = 0.5;
      }
      v[i][j] = temp6;
      v[j][i] = temp6;
    }
  }
}

char *MajorGene2::getname(){
  return name;
}

Environment::Environment(Estimator *est){
  estimator = est;
  sigma = estimator->getenvironment();
  char str[] = "Environment";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}
        
double Environment::getbeta_hat(){
  return sigma;
}

void Environment::getv(double** v, int n){
  for(int i=0;i<n;i++){
    v[i][i] = 1.0;
  }
}

char *Environment::getname(){
  return name;
}

Environment::~Environment() {
  free(name);
}

PolyGene::PolyGene(Estimator *est, ShareRelation *sa){
  estimator = est;
  shareArray = sa;
  relationIndex = 0;
  sigma = estimator->getpolygene();
  char str[] = "PolyGene";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}

PolyGene::~PolyGene() {
  free(name);
}

double PolyGene::getbeta_hat(){
  return sigma;
}

void PolyGene::getv(double** v, int n){
  for(int i=0;i<n;i++){
    v[i][i] = 1.0;
    for(int j=i+1;j<n;j++) {
      v[i][j] = v[j][i] = shareArray[relationIndex++].geneticSimilarity;
    }
  }
}

char *PolyGene::getname(){
  return name;
}

SS::SS(Estimator *est, ShareRelation *sa){
  estimator = est;
  shareArray = sa;
  relationIndex = 0;
  sigma = estimator->getC();
  char str[] = "Sib-Sib";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}

SS::~SS() {
  free(name);
}
                        
double SS::getbeta_hat(){
  return sigma;
}
         
void SS::getv(double** v, int n){
  for(int i=0;i<n;i++){
    v[i][i] = 1.0;
    for(int j=i+1;j<n;j++){
      v[i][j] = v[j][i] = shareArray[relationIndex++].areSiblings;
    }
  }
}

char *SS::getname(){
  return name;
}

PP::PP(Estimator *est, ShareRelation *sa){
  estimator = est;
  shareArray = sa;
  relationIndex = 0;
  sigma = estimator->getP();
  char str[] = "Parent-Parent";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}

PP::~PP() {
  free(name);
}
                        
double PP::getbeta_hat(){
  return sigma;
}
         
void PP::getv(double** v, int n){
  for(int i=0;i<n;i++){
    v[i][i] = 1.0;
    for(int j=i+1;j<n;j++){
      v[i][j] = v[j][i] = shareArray[relationIndex++].areSpouses;
    }
  }
}

char *PP::getname(){
  return name;
}

PS::PS(Estimator *est, ShareRelation *sa){
  estimator = est;
  shareArray = sa;
  relationIndex = 0;
  sigma = estimator->getQ();
  char str[] = "Parent-Sib";
  name = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(name, str);
}

PS::~PS() {
  free(name);
}
                        
double PS::getbeta_hat(){
  return sigma;
}
         
void PS::getv(double** v, int n){
  for(int i=0;i<n;i++){
    v[i][i] = 1.0;
    for(int j=i+1;j<n;j++){
      v[i][j] = v[j][i] = shareArray[relationIndex++].areParentChild;
    }
  }
}

char *PS::getname(){
  return name;
}
