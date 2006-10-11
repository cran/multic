//#include <sys/times.h>
#include <limits.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "Robust_Variance.h"
#include "Calculate.h"
#include "Lib.h"
#include "Model.h"
#include "multic.h"
#include <cstdlib>
#define singular 1.0e-10

using namespace std;

Robust_Variance::Robust_Variance(){
  ybar = 0;
}

Robust_Variance::~Robust_Variance() {
}

void Robust_Variance::main_fun(){
  //OpenFiles();
  //buildModel();
  //fix_effect();
  //CloseFiles();
}


void Robust_Variance::buildModel(Composite *composite, ofstream *outfp,
				 TraitMarkerCov_par *tmc,
				 FortData *fortArray,
				 int dataSize){
  int c1,c2;
/*
  fp_para.open(multic_para_file);
  if (fp_para == NULL) {
  PROBLEM "The file %s could not be opened for reading.\nRobust_Variance.cpp key 28\n",
  multic_para_file RECOVER(NULL_ENTRY);
  }


  fp_data.open(data_file);
  if (fp_data == NULL) {
  PROBLEM "The file %s" could not be opened for reading.\nRobust_Variance.cpp key 36\n",
  data_file RECOVER(NULL_ENTRY);
  }
*/
  ReadFamilyData readfam(tmc, 1, &fam, fortArray, dataSize);
  nfam = readfam.getFamilyCount();
  N = readfam.getFamilySizes();

  itraits = readfam.gettrait();
  icovs = readfam.getcov();   
  iloci = readfam.getmarker();

  missing_flag = Lib::ivector(0, nfam-1);
  for(c1=0;c1<nfam;c1++)
    missing_flag[c1] = 0;

  y = Lib::dmatrix(0,nfam-1,0,FAMMEMNUM-1);
  for(c1=0;c1<nfam;c1++)
    for(c2=0;c2<FAMMEMNUM;c2++)
      y[c1][c2] = 0.0;

  getybar(&readfam);
  //cout<<"ybar="<<ybar<<endl;

  int comp_count = composite->getElementCount();

  double** D1D;
  double** midd;
  D1D  = Lib::dmatrix(0,comp_count-1,0,comp_count-1);
  midd = Lib::dmatrix(0,comp_count-1,0,comp_count-1);
  for(c1=0;c1<comp_count;c1++){
    for(c2=0;c2<comp_count;c2++){
      D1D[c1][c2]  = 0.0;
      midd[c1][c2] = 0.0;
    }
  }

  readfam.resetDataIndex();

  for(int i=0;i<nfam;i++){
    int n = N[i];
    //    cout<<"N["<<i<<"]="<<N[i]<<endl;

    readfam.getdata(n);

    ShrinkData sd(&fam);
    int count = sd.getfamrealcount(n);

    int dim = count*(count+1)/2;

    //get X matrix (correspond to D in my notation)
    double** X = Lib::dmatrix(0,dim-1,0,comp_count-1);

    double** orig_v = Lib::dmatrix(0,count-1,0,count-1);
    for(c1=0;c1<count;c1++)
      for(c2=0;c2<count;c2++)
	orig_v[c1][c2] = 0;

    for(c1=0;c1<dim;c1++)
      for(c2=0;c2<comp_count;c2++) 
	X[c1][c2] = 0;

    composite->getvsum(orig_v, X, n, &sd);
		
    double** W = Lib::dmatrix(0,dim-1,0,dim-1);
    int irow=0;
    int icol=0;
    for(c1=0; c1<count; c1++){
      for(c2=c1; c2<count; c2++){
	icol=0;
	for(int l=0; l<count; l++){
	  for(int m=l; m<count; m++){
	    W[irow][icol]=orig_v[c1][l]*orig_v[c2][m]+
	      orig_v[c1][m]*orig_v[c2][l];
	    icol++;
	  }
	}
	irow++;
      }
    }

    double** W_inv = Lib::dmatrix(0,dim-1,0,dim-1);
    double return_val = Lib::InverseOfMatrix(W_inv,dim,W);

    if(return_val==BAD){
      cout << "Singular matrix in Robust Variance calculation" << endl;	
      // << "Robust_Varience.cpp line 185" << endl;
      // This is not necessarily an error.  It appears to not be a major
      // problem to cause multic to exit(1).  Eric Lunde, 2005-09-07
      return;
    }

    //X1 = X'
    double** X1 = Lib::dmatrix(0,comp_count-1,0,dim-1);
    Lib::Transpose(X1,X,dim,comp_count);

    //X1W_inv = X'*W_inv
    double** X1W_inv = Lib::dmatrix(0,comp_count-1,0,dim-1);
    Lib::Multi_Matrices(X1W_inv,X1,comp_count,dim,W_inv,dim,dim);

    double** X1X = Lib::dmatrix(0,comp_count-1,0,comp_count-1);
    Lib::Multi_Matrices(X1X,X1W_inv,comp_count,dim,X,dim,comp_count);

    for(c1=0;c1<comp_count;c1++){
      for(c2=0;c2<comp_count;c2++){
	D1D[c1][c2] += X1X[c1][c2];
      }
    }

    double* Y = Lib::dvector(0,dim-1);
    //shrink data
    //getY(sd,Y,n);

    int countt=0;
    for(c1=0;c1<count;c1++){
      for(c2=c1;c2<count;c2++){
	Y[countt++] = y[i][c1]*y[i][c2];
      }   
    }

    /*
      cout<<"orignal Y vector"<<endl;   
      for(c1=0;c1<dim;c1++)
      cout<<Y[c1]<<endl;
    */

    double* Y_hat = Lib::dvector(0,dim-1);
    double* beta_hat = Lib::dvector(0,comp_count-1);

    composite->getbeta_hat(beta_hat,comp_count);

    //cout<<"Y hat=DB=E(Y) vector"<<endl;

    for(c1=0;c1<dim;c1++){
      Y_hat[c1] = 0.0;
      for(int c2=0;c2<comp_count;c2++)
	Y_hat[c1] += X[c1][c2]*beta_hat[c2];
      //cout<<Y_hat[c1]<<endl;
    }

    //now, Y = Y - XB (E)
    //cout<<"Y=Y-E(Y)"<<endl;
    for(c1=0;c1<dim;c1++){
      Y[c1] = Y[c1] - Y_hat[c1];
      //cout<<Y[c1]<<"  ";
    }
    //cout<<endl;

    //newY = X'[Y-E(Y)]
    //cout<<"D'[Y-E(Y)]"<<endl;
    double* newY = Lib::dvector(0,comp_count-1);
    for(c1=0;c1<comp_count;c1++){
      newY[c1] = 0.0;
      for(c2=0;c2<dim;c2++)
	newY[c1] += X1W_inv[c1][c2]*Y[c2];
      //cout<<newY[c1]<<endl;
    }
                
    //get newY*newY'
    //           -1	                  -1
    //cout<<"D'(W) [Y-E(Y)][Y-E(Y)]'(W) D"<<endl;
    for(c1=0;c1<comp_count;c1++){
      for(c2=0;c2<comp_count;c2++){
	midd[c1][c2] += newY[c1]*newY[c2];
	//cout<<newY[c1]*newY[c2]<<"  ";
      }
      //cout<<endl;
    }

    //add(result,X,dim,comp_count,V,dim,dim,Y,dim);
    Lib::free_dmatrix(X, 0, dim-1, 0, comp_count-1);
    Lib::free_dmatrix(orig_v, 0, count-1, 0, count-1);
    Lib::free_dmatrix(W, 0, dim-1, 0, dim-1);
    Lib::free_dmatrix(W_inv, 0, dim-1, 0, dim-1);
    Lib::free_dmatrix(X1, 0, comp_count-1, 0, dim-1);
    Lib::free_dmatrix(X1W_inv, 0, comp_count-1, 0, dim-1);
    Lib::free_dmatrix(X1X, 0, comp_count-1, 0, comp_count-1);

    Lib::free_dvector(Y, 0, dim-1);
    Lib::free_dvector(Y_hat, 0, dim-1);
    Lib::free_dvector(beta_hat, 0, comp_count-1);
    Lib::free_dvector(newY, 0, comp_count-1);
  }

  int nomissingcount = 0;
  int *nomissing = Lib::ivector(0,comp_count-1);
  for(c1=0;c1<comp_count;c1++){
    nomissing[c1] = -1;
    for(int c2=0;c2<comp_count;c2++){
      if(D1D[c1][c2]!=0) nomissing[c1] = 0;
    }
    if(nomissing[c1]==0) nomissingcount++;
  }

  double** D1D_new = Lib::dmatrix(0,nomissingcount-1,0,nomissingcount-1);
  double** midd_new = Lib::dmatrix(0,nomissingcount-1,0,nomissingcount-1);
  int row=0,column;
  for(c1=0;c1<comp_count;c1++){
    column = 0;
    if(nomissing[c1]==0){
      for(c2=0;c2<comp_count;c2++)
	if(nomissing[c2]==0){
	  D1D_new[row][column]=D1D[c1][c2];
	  midd_new[row][column]=midd[c1][c2];
	  column++;
	}
      row++;
    }
  }

  double** inv = Lib::dmatrix(0,nomissingcount-1,0,nomissingcount-1);
  double return_val = Lib::InverseOfMatrix(inv,nomissingcount,D1D_new);

  if(return_val==BAD){
    cout << "Singular matrix in Robust Variance calculation" << endl;
    // << "Robust_Varience.cpp line 263" << endl;
    // This is not necessarily an error.  It appears to not be a major
    // problem to cause multic to exit(1).  Eric Lunde, 2005-09-07
    return;
  }
  double** temp2 = Lib::dmatrix(0,nomissingcount-1,0,nomissingcount-1);
  double** result = Lib::dmatrix(0,nomissingcount-1,0,nomissingcount-1);
 
  Lib::Multi_Matrices(temp2,nomissingcount,inv,midd_new);
  Lib::Multi_Matrices(result,nomissingcount,temp2,inv);

  *outfp << endl << "Variance Covariance Matrix:" << endl;
  for(c1=0;c1<comp_count;c1++){
    if(nomissing[c1]==0){
      Model* temp = composite->getElementAt(c1);
      *outfp << setw(15) << temp->getname() << "  ";
    }
  }
  *outfp << endl;

  ofstream varCovar("varCovar.log", ios::app);
  if(varCovar.fail()) {
    PROBLEM "The file varCovar.log could not be opened for appending.\nRobust_Variance.cpp key 261\n"
      RECOVER(NULL_ENTRY);
  }

  //cout<<"result"<<endl;
  for(c1=0;c1<nomissingcount;c1++){
    for(int c2=0;c2<nomissingcount;c2++){
      //cout<<result[c1][c2]<<"  ";
      *outfp << setw(15) << result[c1][c2] << "  ";
      /// add file pint statements here - eric lunde 10-22-03
      varCovar << result[c1][c2] << endl;
    }
    *outfp << endl;
    //cout<<endl;
  }

  varCovar.close();
  // Width 0 and precision 6 are the default values
  outfp->width(0);
  outfp->precision(6);
  *outfp << endl;

  //  fp_para.close();
  //  fp_data.close();

  Lib::free_ivector(missing_flag, 0, nfam-1);
  Lib::free_dmatrix(y, 0, nfam-1, 0, FAMMEMNUM-1);
  Lib::free_dmatrix(D1D, 0, comp_count-1, 0, comp_count-1);
  Lib::free_dmatrix(midd, 0, comp_count-1, 0, comp_count-1);
  Lib::free_dmatrix(D1D_new, 0, nomissingcount-1, 0, nomissingcount-1);
  Lib::free_dmatrix(midd_new, 0, nomissingcount-1, 0, nomissingcount-1);
  Lib::free_dmatrix(inv, 0, nomissingcount-1, 0, nomissingcount-1);
  Lib::free_dmatrix(temp2, 0, nomissingcount-1, 0, nomissingcount-1);
  Lib::free_dmatrix(result, 0, nomissingcount-1, 0, nomissingcount-1);
  Lib::free_ivector(nomissing, 0, comp_count-1);
}


// get Y of the second moment model
void Robust_Variance::getY(ShrinkData *sd, double* Y, int n){
  int c1,c2;
  int count = sd->getfamrealcount(n);
  double* old_y;
  old_y = Lib::dvector(0, count-1);
  int i = 0;
  int j=0;
  for(j=0;j<n;j++){
    if (fam.person[j].missing_traitflag != 1) {
      old_y[i] = fam.person[j].trait[0]-ybar;
      i++;
    }
  }

  cout<<"old_y"<<endl;
  for(j=0;j<count;j++)
    cout<<old_y[j]<<endl;

  double** s;
  s = Lib::dmatrix(0,count-1,0,count-1);

  /*
    cout<<"S matrix"<<endl;
    for(c1=0;c1<count;c1++){
    for(c2=0;c2<count;c2++){
    s[c1][c2] = old_y[c1]*old_y[c2];
    cout<<s[c1][c2]<<"  ";
    }
    cout<<endl;
    }
  */

  //y = (double *) malloc((count*(count+1)/2) * sizeof(double));
  i = 0;
  for(c1=0;c1<count;c1++)
    for(c2=c1;c2<count;c2++)
      Y[i++] = s[c1][c2];

  Lib::free_dmatrix(s, 0,count-1,0,count-1);
  Lib::free_dvector(old_y, 0, count-1);
}

// get overall mean
void Robust_Variance::getybar(ReadFamilyData *readfam) {
  int i=0;
  int count = 0;
  int ttcount = 0;
  int index = 0;
  ifstream data(data_file);
  //get ride of header
  if(header == 1){
    int temp;
    for(int i=0;i<=nfam;i++){
      data >> temp;
    } 
    while(data.get() != (int)('\n'));
    // delete the newline character.
  }

  for (i=0; i<nfam; i++){
    readfam->getdata(N[i]);
    index = 0;
    for(int j=0;j<N[i];j++){
      ttcount++;
      if (fam.person[j].missing_traitflag != 1) {
	ybar += fam.person[j].trait[0];
	y[i][index] = fam.person[j].trait[0];
	//cout<<y[i][index]<<endl;
	missing_flag[i]++;
	count++;
	index++;
      }
    }
  }
  //cout<<"ttcount="<<ttcount<<endl;
  //cout<<"count="<<count<<endl;
  ybar = ybar/count;
  //cout<<"y"<<endl;
  for (i=0; i<nfam; i++){
    for(int j=0;j<N[i];j++){
      y[i][j] -= ybar;
      //cout<<y[i][j]<<"  ";
    }
    //cout<<endl;
  }
  data.close();
}
