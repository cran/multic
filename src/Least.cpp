/******************************************************************************
File: Least.cpp
Description: Least.cpp is designed for multi-traits variance component
             analysis and using the overall mean of each trait as the mean of
             the trait values of all fathers, mothers and sibs.  Also,
             Least.cpp includes least squares estimation for boundary.
Author: Eric Lunde, 7-9-03 (This is the date Eric began documnetation of this
                           file)
Updates: (Date, Modified By, Modification Description)
7-9-03, Eric Lunde, One of the resposiblities of the method 'void
        Least::getresponse()' is to retreive the data from fort.12.  The
        internal values were not matching the file itself.  The reason for this
        was the file was being closed at the end of the constructor, which was
        before the method was called.  Therefore when it tried to get the data,
        the file was closed.  I moved the CloseFiles() method to the destructor
        to prevent the premature file close.
7-11-03, Eric Lunde, The format of loci.out has changed to have the name of
         the (m)ibd file at the top of the file.  I added multiple lines of
         code to simply read past that title line.
8-6-03, Eric Lunde, I removed the code in getresponse that read directly from
        fort.12.  It has been replaced with the necessary ReadFamilyData
        methods to get the same data.

        I also added ShareRelation *shareArray, int relationSize, and int
        relationIndex to the member variable list.  The method getcoeffecient
        was reading from share.out.  Since we have that information already
        stored in main memory, we should be reading from that instead of the
        file.
******************************************************************************/

#include "Least.h"
#include "Lib.h"
#include "Calculate.h"
#include "TraitMarkerCov_par.h"
#include "Lib.h"
#include "ReadFamilyData.h"
#include "multic.h"
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <S.h>
using namespace std;

Least::Least(TraitMarkerCov_par *tmc, ShareRelation *sa, int rs,
	     FortData *fortArray, int dataSize) :
  outfp("least.log", ios::app),
  pt_loci(loci_file)
{

  shareArray = sa;
  relationSize = rs;
  relationIndex = 0;

  m = 0;
  A11=0.0;
  A12=0.0;
  A22=0.0;
  int i,j;
  for(i=0; i<MAXNUMTRAIT; i++){
    ybar[i] = 0;
    var_e[i]=0;
    for(j=0;j<MAXNUMTRAIT;j++){
      var_g[i][j]=0;
      var_G[i][j]=0;
    }
  }

  //  outfp("least.log", ios::app);
  if(!outfp) {
    PROBLEM "The file least.log could not be opened for appending.\nLeast.cpp key 28\n"
      RECOVER(NULL_ENTRY);
  }

  if(pt_loci == NULL) {
    PROBLEM "The file %s could not be opened for reading.\nLeast.cpp key 40\n",
      loci_file RECOVER(NULL_ENTRY);
  }

  itraits = tmc->gettraitnum();
  t_start = tmc->get_t_start();
  t_len = tmc->get_t_len();

  readfam = new ReadFamilyData(tmc, 1, &fam, fortArray, dataSize);
  familyCount = readfam->getFamilyCount();
  familySizes = readfam->getFamilySizes();

  missing_flag = Lib::imatrix(0, familyCount-1, 0, FAMMEMNUM-1);
  for(i=0; i<familyCount; i++)
    for(j=0; j<FAMMEMNUM; j++)
      missing_flag[i][j] = 1;

  readfam->resetDataIndex();

  for(i = 0; i < familyCount; i++) {
    readfam->getdata(familySizes[i]);
    for(j = 0; j < familySizes[i]; j++)
      if (fam.person[j].missing_traitflag != 1) 
	missing_flag[i][j] = 0;
  }
}

Least::~Least() {
  for(int i=0; i<familyCount; i++){
    Lib::free_dmatrix(a[i], 0, FAMMEMNUM-1, 0, FAMMEMNUM-1);
    Lib::free_dmatrix(r_mat[i], 0, FAMMEMNUM-1, 0, FAMMEMNUM-1);
  }
  for(int i=0; i<itraits; i++) {
    Lib::free_dmatrix(y[i], 0, familyCount-1, 0, FAMMEMNUM-1);
  }
  Lib::free_imatrix(missing_flag, 0, familyCount-1, 0, FAMMEMNUM-1);
  delete readfam;
  CloseFiles();
}

void Least::getresponse(){
  /* Read fort.12 file to obtain yij and ybar--the overall mean of yij. */
  /* contribute N, n[i], y and ybar to the class*/
  int i, j, k, yu1;

  outfp << endl << endl << " **************************************************";
  outfp << endl << "Total number of familys: " << familyCount << ' ';

  y = (double ***)malloc(itraits * sizeof(double **));
  a = (double ***)malloc(familyCount * sizeof(double **));
  r_mat = (double ***)malloc(familyCount * sizeof(double **));

  for(i=0; i<familyCount; i++){
    a[i] = Lib::dmatrix(0,FAMMEMNUM-1,0,FAMMEMNUM-1);
    r_mat[i] = Lib::dmatrix(0,FAMMEMNUM-1,0,FAMMEMNUM-1);
  }
  for(i=0; i<itraits; i++)
    y[i] = Lib::dmatrix(0, familyCount-1, 0, FAMMEMNUM-1);

  FortData *fortArray = readfam->getFortArray();
  int dataIndex = 0;

  for(i = 0; i < familyCount; i++) {
    for(j = 0; j < familySizes[i]; j++, dataIndex++){
      for ( yu1 = 0; yu1 < itraits; yu1++) {
	y[yu1][i][j] = fortArray[dataIndex].traits[yu1];
      }

      /* trait values in fort.12 are long float. */
      for(k=0; k < itraits; k++){
	if (missing_flag[i][j] != 1) {
	  ybar[k] = ybar[k] + y[k][i][j];
	  // printf("%.5f \n",y[k][i][j]);
	}
      }
      if (missing_flag[i][j] != 1) m++;
    }
  }
  outfp << endl << "Total number of individuals: " << m << ' ' << endl;
  for(k=0; k<itraits; k++) {
    ybar[k] = ybar[k]/m;
    outfp << endl << "Mean of trait " << k+1 << ": " << ybar[k] << endl;
  }
}

void Least::CloseFiles() {
  outfp.close();
  pt_loci.close();
}

void Least::getcoefficient(){
  /* Read loci.out file to obtain number of aij=0, aij!=0 and */
  /* calculate coeficient A11, A12, A22.  */
  int i, j, k;
  double temp1, temp2;
  char str[BUF];
  // This next line is to read past the new (7-11-03) title line at the
  // beginning of loci.out - Eric Lunde
  pt_loci.getline(str, BUF);

  for(i=0; i<familyCount; i++) {
    for(j=0; j<familySizes[i]-1; j++) {
      for(k=j+1; k<familySizes[i]; k++) {
	pt_loci >> temp1 >> temp2 >> a[i][j][k];
	pt_loci.getline(str, BUF);

	r_mat[i][j][k] = shareArray[relationIndex++].geneticSimilarity;
/*	fp_share >> r_mat[i][j][k] >> temp1 >> temp2
		 >> temp3 >> temp4 >> temp5;
	fp_share.getline(str, BUF);
*/

	a[i][k][j] = a[i][j][k];
	r_mat[i][k][j] = r_mat[i][j][k];
	//if(a[i][j][k]==0.0 && missing_flag[i][j] != 1) { 
	//  a0 = a0 + 1;
	//}else if(a[i][j][k]!=0.0 && missing_flag[i][j] != 1) {
	if(missing_flag[i][j] != 1) {
	  A11 = A11 + 2*a[i][j][k]*a[i][j][k];
	  A12 = A12 + 2*r_mat[i][j][k]*a[i][j][k];
	  A22 = A22 + 2*r_mat[i][j][k]*r_mat[i][j][k];
	  //					a1 = a1 + 1;
	}
      }
    }
  }


  //		A22 = 0.5*(a1 + a0 - N);
  // printf("\nThe numbers of aij=0, and !=0 are %d, and %d. \n", a0, a1);
  // printf("A11=%f A12=%f A22=%f \n",A11,A12,A22);
  // fclose(pt_loci);
}


void Least::getestimator(){
  int i, j, k, l, h;
  double  b1[MAXNUMTRAIT][MAXNUMTRAIT], b2[MAXNUMTRAIT][MAXNUMTRAIT],
    b3[MAXNUMTRAIT][MAXNUMTRAIT];

  /* Calculate b1, b2 and b3 */
  for(l=0; l<itraits; l++) {
    for(h=l; h<itraits; h++) {
      b1[l][h] = 0;
      b2[l][h] = 0;
      b3[l][h] = 0;
      /***************************
      The following algorithm was detailed by Eric Lunde on 6-19-03.
      for (each family) {
        for (each family member, j) {
          for (each family member, k) {
            if (family member j in famliy i has a missing value) {
              if (family member j is family member k) {
		b3[l][h] += (y[l][i][j]-ybar[l]) * (y[h][i][k]-ybar[h]);
              }else if (family member i and k are the first two, but different,
                         members of the household [I don't know if this
                         necessarily means that they are the parents or not]) {
		b1[l][h] = b1[l][h];
		b2[l][h] = b2[l][h];
              }else {
		b1[l][h] += a[i][j][k]
		  * (y[l][i][j]-ybar[l]) * (y[h][i][k]-ybar[h]);
		b2[l][h] += r_mat[i][j][k]
		  * (y[l][i][j]-ybar[l]) * (y[h][i][k]-ybar[h]);
              }
            }
          }
        }
      }
      ***************************/
      for(i=0; i<familyCount; i++) {
	for(j=0; j<familySizes[i]; j++) {
	  for(k=0; k<familySizes[i]; k++) {
	    if(missing_flag[i][j]!=1){
	      if(j==k) {
		b3[l][h] += (y[l][i][j]-ybar[l]) * (y[h][i][k]-ybar[h]); 
	      }else if( (j==0 && k==1)||(j==1 && k==0)) {
		b1[l][h] = b1[l][h];
		b2[l][h] = b2[l][h];
	      }else {
		b1[l][h] += a[i][j][k]
		  * (y[l][i][j]-ybar[l]) * (y[h][i][k]-ybar[h]);
		b2[l][h] += r_mat[i][j][k]
		  * (y[l][i][j]-ybar[l]) * (y[h][i][k]-ybar[h]);
	      }
	    }
	  }   
	}      
      }

      /* Calculate var_g, var_G, var_e */  
      if(l == h) {
	/// assigning var_g
	var_g[l][h] = (b1[l][h]*A22 - b2[l][h]*A12) / (A11*A22 - A12*A12);
	var_G[l][h] = (b2[l][h]*A11 - b1[l][h]*A12)/(A11*A22 - A12*A12);
	var_e[l] = b3[l][h]/m - var_g[l][h] - var_G[l][h];
	if(var_g[l][h] < 0) {
	  /// assigning var_g
	  var_g[l][h] = 0;
	  var_G[l][h] = b2[l][h]/A22;
	  if(var_G[l][h] < 0) {
	    var_G[l][h] = 0;
	  }
	  var_e[l] = b3[l][h]/m - var_G[l][h];
	}
	if(var_G[l][h] < 0) {
	  var_G[l][h] = 0;
	  /// assigning var_g
	  var_g[l][h] = b1[l][h]/A11;
	  if(var_g[l][h] < 0) {
	    var_g[l][h] = 0;
	  }
	  var_e[l] = b3[l][h]/m - var_g[l][h];
	}
	if(var_e[l] < 0) {   
	  var_e[l] = 0;   
	  if(var_g[l][h] < 0) {
	    /// assigning var_g
	    var_g[l][h] = 0;
	    var_G[l][h] = (b3[l][h]+b2[l][h])/(m+A22);
	  }
	  if(var_G[l][h] < 0) {
	    var_G[l][h] = 0;
	    /// assigning var_g
	    var_g[l][h] = (b3[l][h]+b1[l][h])/(m+A11);
	  }
	}
      }
      if(l < h) {
	/// assigning var_g
	var_g[l][h] = ((b3[l][h]+b1[l][h])*(m+A22)-(b3[l][h]+b2[l][h])*(m+A12))/((m+A11)*(m+A22)-(m+A12)*(m+A12));
	var_G[l][h] = ((m+A11)*(b3[l][h]+b2[l][h])-(m+A12)*(b3[l][h]+b1[l][h]))/((m+A11)*(m+A22)-(m+A12)*(m+A12));
      }
    }
  }
  Lib::free_imatrix(missing_flag, 0, familyCount-1, 0, FAMMEMNUM-1);
  free(y);
  free(a);
  free(r_mat);
}


void Least::least_main(){
  int h,l;
  double time_secs;
  time_t t1, t2;
  clock();
  time(&t1);
  getresponse();
  getcoefficient();
  getestimator();
  time_secs=clock()*1.0/CLOCKS_PER_SEC;

  time(&t2);
        
  for(l=0; l < itraits; l++) {
    for(h=l; h < itraits; h++) {
      outfp << endl << "The polygene s" << l+1 << h+1 << " = "
	    << setprecision(3) << var_G[l][h] << endl;
      outfp << endl << "The major gene m" << l+1 << h+1 << " = "
	    << setprecision(3) << var_g[l][h] << endl;
      if(l == h) {
	outfp << endl << "The environment t" << l+1 << h+1 << " = "
	      << setprecision(3) << var_e[l] << endl;
      }
    }
  }
  
  outfp << endl << "The CPU time: " << time_secs << " seconds." << endl;
  outfp << endl << "The real time: " << (t2-t1) << " seconds." << endl;
  outfp.close();
}


void Least::getmajorgene(double init_x_M1[INITVCNUM]){
  int k=0;
  for(int i=0;i<itraits;i++)
    for(int j=i;j<itraits;j++)
      init_x_M1[k++] = var_g[i][j];
}


void Least::getpolygene(double init_x_S[INITVCNUM]){
  int k = 0;
  for(int i=0;i<itraits;i++)
    for(int j=i;j<itraits;j++)
      init_x_S[k++] = var_G[i][j];
}


void Least::getenvironment(double init_x_T[INITVCNUM]){
  int k=0;
  for(int i=0;i<itraits;i++)
    init_x_T[k++] = var_e[i];
}


double Least::getmajorgene1(){
  return var_g[0][0];
}

//possible extension of current least square version
double Least::getmajorgene2(){
  return 0;
}

double Least::getpolygene(){
  return var_G[0][0];
}

double Least::getenvironment(){
  return var_e[0];
}

//possible extension of current least square version
double Least::getC(){
  return 0;
}

//possible extension of current least square version
double Least::getP(){
  return 0;
}

//possible extension of current least square version
double Least::getQ(){
  return 0;
}
