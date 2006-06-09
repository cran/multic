#include "ShrinkData.h"

ShrinkData::ShrinkData(Family *famd){
  famdata = famd;
}

int ShrinkData::getfamrealcount(int n_fam_member) {
  int count = 0;
  for (int k = 0; k < n_fam_member; k++) {
    if (famdata->person[k].missing_traitflag != 1) {   
      count++;
    }
  }
  return count;
}

void ShrinkData::Changedata(double **g_mat, int n_fam_member) {
  int k;
  int prev_count;
  int n_mb;
  prev_count = 0;
  n_mb       = n_fam_member;

  for (k = 0; k < n_fam_member; k++) {
    if (famdata->person[k].missing_traitflag == 1) {
      Shrinkdata(g_mat, k-prev_count, n_mb);
      prev_count++;
      n_mb--;
    }
  }
}


void ShrinkData::Changedata(double *g_vec, int n_fam_member) {
  int k;
  int prev_count;
  int n_mb;
  prev_count = 0;
  n_mb       = n_fam_member;

  for (k = 0; k < n_fam_member; k++) {
    if (famdata->person[k].missing_traitflag == 1) {
      Shrinkdata(g_vec, k-prev_count, n_mb);
      prev_count++;
      n_mb--;
    }
  }
}


void ShrinkData::Shrinkdata(double **g_mat, int n_order, int n_member){
  int i,j,ii,jj;
  //    int numofcov;

  for(i = 0,ii = 0; i < n_member-1; i++,ii++) {
    for(j = 0,jj = 0; j < n_member-1; j++,jj++) {
      if ( (i == n_order) && (j == n_order) ) {
	ii = i+1;
	jj = j+1;
      }
      else if (i == n_order)
	ii = i+1;
      else if (j == n_order)
	jj = j+1;
      g_mat[i][j] = g_mat[ii][jj];
    }
  }
}


void ShrinkData::Shrinkdata(double *Y1, int n_order, int n_member){
  int i,ii; 
  //    int numofcov;

  for(i = 0,ii = 0; i < n_member-1; i++,ii++) {
    if (i == n_order) ii = i+1;
    Y1[i] = Y1[ii];
  }
}
