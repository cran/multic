#include "Composite.h"
#include "Model.h"
#include "Lib.h"
#include <cstdlib>
#include <iostream>
#include <cstring>
#include "Rostream.h"
#include "Rstreambuf.h"
using namespace Rcpp;


Composite::Composite(){
  //  vector = new Vector();
}

Composite::~Composite() {
}

void Composite::add(Model *model){
  vector.add(model);
}     

int Composite::getElementCount(){
  return vector.getcount();
}

double Composite::getbeta_hat(){
  Rcout<<"Error found at composite: getbeta_hat()"<<std::endl;
  return 0;
}

void Composite::getv(double **v, int n) {
  Rcout << "Call wrong method" << std::endl;
}

int Composite::getbeta_hat(double* beta, int n){
  if(n!=vector.getcount()) return -1;
  for(int i=0;i<n;i++){
    Model* temp = (Model*) vector.elementAt(i);
    beta[i] = temp->getbeta_hat();
  }
  return 1;
}

Model* Composite::getElementAt(int i){
  return (Model*) vector.elementAt(i);
}

// get v matrix which is stack of the upper triangle each gene component
void Composite::getvsum(double** vorig, double** vstack, int n,
			ShrinkData *sd)
{
  int c1,c2;

  //for(int j=0;j<n;j++){
  //  printf("trait=%f\n",fam.person[j].trait[0]);
  //  printf("missing=%d\n",fam.person[j].missing_traitflag);
  //}
  // can't create an object this way
  //ShrinkData sd(fam);

  int count = sd->getfamrealcount(n);

  double** tempv = Lib::dmatrix(0,n-1,0,n-1);

  for(int i=0;i<vector.getcount();i++){
    Model* temp = (Model*) vector.elementAt(i);
    //initialize tempv matrix
    for(c1=0;c1<n;c1++)
      for(c2=0;c2<n;c2++)
	tempv[c1][c2] = 0.0;
    //get individual matrix

    temp->getv(tempv,n);

    double localestimator = temp->getbeta_hat();
    //cout<<"beta="<<localestimator<<std::endl;
    sd->Changedata(tempv,n);
    int kk = 0;
    //stack into one column of V matrix
    for(c1=0;c1<count;c1++){
      for(c2=c1;c2<count;c2++){
	vstack[kk++][i] = tempv[c1][c2];
	//cout<<tempv[c1][c2]<<"  ";
      }
      //cout<<"next"<<std::endl;
    }
    for(c1=0;c1<count;c1++){
      for(c2=0;c2<count;c2++){
	vorig[c1][c2] += localestimator*tempv[c1][c2];
      }
    }
  }
  Lib::free_dmatrix(tempv, 0, n-1, 0, n-1);
}

char* Composite::getname(){
  char str[] = "Composite";
  char *returnStr = (char *) malloc( (strlen(str) + 1) * sizeof(char));
  strcpy(returnStr, str);
  return returnStr;
}
