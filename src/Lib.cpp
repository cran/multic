#include <iomanip>
#include <cstdlib>
#include <iostream>
#include "Lib.h"
#include <stdlib.h>
#include <math.h>
#include "Rostream.h"
#include "Rstreambuf.h"

using namespace Rcpp;
#define NR_END 1

void Lib::PrintOneVector(int n, double *vec) {
  int i;

  for( i = 0; i < n; i++) {
    //cout << std::setw(25) << setprecision(15) << vec[i] << ' ' << std::endl;
    Rcout << std::setw(15) << vec[i] << ' ' << std::endl;
  }
}

void Lib::PrintOneVector(int n, int *vec) {
  int i;

  for( i = 0; i < n; i++) {
    Rcout << std::setw(15) << vec[i] << ' ' << std::endl;
  }
}

/* Print each element of a matrix. 
*/
void Lib::PrintOneMatrix(int nr, int nc, double **a) {
  int i, j;

  for( i = 0; i < nr; i++) {
    for( j = 0; j < nc; j++) {
      Rcout << std::setw(15) << a[i][j] << ' ';
    }
    Rcout << std::endl;
  }
}


/* get the determinant of an LU decomposed matrix is just the product of the
   diaginal elements.
*/
double Lib::DetermOfMatrix(int n,  double **a){
  int i,j, *indx;
  double  rt_flag;
  double d;

  d = 1.0;
  indx = ivector(0,n-1);

  for(i=0;i<n;i++) 
    indx[i] = 0;

  rt_flag = ludcmp(a, n, indx, &d);

  if (rt_flag == BAD) {
    return BAD;
  }

  for (j = 0; j < n; j++) {
    d *= a[j][j];
  }

  Lib::free_ivector(indx, 0, n-1);

  return d;
}


/* copy one matrix (a_in) to another one, and return the copy.
*/
void Lib::CopyMatrix(double **a_out, int n, double **a_in){
  int i, j;

  for( i = 0; i < n; i++) {
    for( j = 0; j < n; j++) {
      a_out[i][j] = a_in[i][j];
    }
  }
}

/* find the inverse of a matrix column by column, a is the input and 
   return the inverse matrix of a.   */
double Lib::InverseOfMatrix(double **y, int n, double **a){
  int i,j, *indx;
  double rt_flag;
  double d, *col;

  indx = ivector(0,n-1);
  col  = dvector(0,n-1);

  for(i=0;i<n;i++)
    indx[i] = 0;

  rt_flag = ludcmp(a,n,indx, &d);
  if (rt_flag == BAD)
    return BAD;

  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++)
      col[i] = 0.0;
    col[j] = 1.0;
    lubksb(a,n,indx,col);
    for (i = 0; i < n; i++)
      y[i][j] = col[i];
  }

  Lib::free_ivector(indx, 0, n-1);
  Lib::free_dvector(col, 0, n-1);

  return GOOD;
}

/* Multiply two square matrices, mat1 and mat2 are input, return mat_out
   matrix.   
*/
void Lib::Multi_Matrices(double **mat_out, int n, double **mat1, double **mat2){
  int i,j,k;

  for (i = 0; i < n; i++) {
    for ( j = 0; j < n; j++) {
      mat_out[i][j] = 0;
      for (k = 0; k < n; k++) {
	mat_out[i][j] = mat_out[i][j] + mat1[i][k]*mat2[k][j];
      }
    }
  }
}

/* transpose a matrix */

void Lib::Transpose(double **mat_out,double **mat1,int c1,int c2){
  for(int i=0;i<c1;i++)
    for(int j=0;j<c2;j++)
      mat_out[j][i] = mat1[i][j];
}

/* Multiply any two matrices, mat1 and mat2 are input, return mat_out
   matrix.
*/
int Lib::Multi_Matrices(double **mat_out,double **mat1,int c1,int c2,double **mat2,int c3,int c4){
  int i,j,k;

  //check dimension agreement
  if(c2!=c3) return -1;

  for (i = 0; i < c1; i++) {
    for ( j = 0; j < c4; j++) {
      mat_out[i][j] = 0;
      for (k = 0; k < c2; k++) {
	mat_out[i][j] += mat1[i][k]*mat2[k][j];  
      } 
    }
  }
  return 0;
}


/* Multiply two vectors and return a double value.
*/
double Lib::Multi_Vectors(int n, double *vec1, double *vec2){
  int    i;
  double *vec;
  double rv;

  vec = dvector(0,n-1);

  rv = 0.0;
  for (i = 0; i < n; i++) {
    rv += vec1[i] * vec2[i];	
  }

  Lib::free_dvector(vec, 0, n-1);

  return rv; 
}


/* Multiply matrix and vector, mat and vec are input, 
   vec_out is the output vector.
*/
void Lib::Multi_Matrix_Vector(double *vec_out, int nr, int nc, double **mat,
			      double *vec){
  int    i,j;

  /* a matrix multiplies a vector.    */
  for ( i = 0; i < nr; i++) {
    vec_out[i] = 0;
    for ( j = 0; j < nc; j++) {
      vec_out[i] += mat[i][j] * vec[j];
    }
    /*
      printf("nr=%d,nc=%d\nvec_out[%d]=%lf\n",nr,nc,i,vec_out[i]);
    */
  }
}


/* Multiply vectors and matrix, vec1, mat and vec2 are input, return a double
   value.
*/
double Lib::Multi_Vectors_Matrix(int n, double *vec1, double **mat,
				 double *vec2){
  int    i,j;
  double *temp_vec;
  double rv;

  /* Initialize the size of the pointer.   5-12-98; E.Y.  */
  temp_vec = dvector(0, n-1);

  /* a vector multiplies a matrix.    */ 
  for ( j = 0; j < n; j++) {
    temp_vec[j] = 0;
    for ( i = 0; i < n; i++) {
      temp_vec[j] += vec1[i] * mat[i][j];  
    }
  }

  /* a vector multiplies a vector.    */
  rv = 0.0;
  for (i = 0; i < n; i++) {
    rv += temp_vec[i] * vec2[i];	
  }

  Lib::free_dvector(temp_vec, 0, n-1);

  return rv; 
}

/* Get the trace value of a matrix, return a double value.
*/
double Lib::Trace_Matrix(int n, double **mat){
  int i;
  double rv;

  rv = 0.0;
  for (i = 0; i < n; i++) {
    rv += mat[i][i];	
  }
  return rv;
}


/* Given an nxn matrix a[1..n][1..n], this routine replaces it by the LU 
   decomposition of a rowwise permutation of itself. a and n are input. a is
   output, arranged as in equation; indx[1..n] is an output vector which
   records the row permutation effected by the partial pivoting; d is output
   as +(-)1 depending on whether the number os row interchanges was even or
   odd, respectively. THis routine is used in combination with lubksb to
   solve linear equations or invert a matrix.
*/
double Lib::ludcmp(double **a, int n, int *indx, double *d){
  int     i, imax = 0, j, k;
  double  big, dum, sum, temp;
  double  *vv;        
  /* vv stores the implicit scaling of each row. */
  vv = dvector(0,n-1);
  *d = 1.0;                 /* No row interchanges yet.    */
  for (i = 0; i < n; i++) { /* Loop over rows to get the implicit scaling
			       information.                */
    big = 0.0;
    for (j = 0; j < n; j++) {
      if ((temp = fabs(a[i][j])) > big) 
	big = temp;
    }
    if (big == 0.0) {      /* No nonzero largest element. */
      return BAD;        /* Singular matrix here.   changed 7/31/96.  */

      /*
	nrerror("Singular matrix in routine LUDCMP");
      */
    }

    vv[i] = 1.0/big;      /* Save the scaling.           */
  }

  for (j = 0; j < n; j++) { /* This is the loop over columns of Crout's
			       method.   */
    for (i = 0; i < j; i++) {
      sum = a[i][j];
      for (k = 0; k < i; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;	
    }
    big = 0.0;            /* Initialize for the search for largest pivot
			     element.   */
    for (i = j; i < n; i++) {
      sum = a[i][j];
      for (k = 0; k < j; k++)
	sum -= a[i][k]*a[k][j];
      a[i][j] = sum;	
      if ( (dum = vv[i]*fabs(sum)) >= big) {
	big  = dum;
	imax = i;
      }
    }
    if (j != imax) {      /* Do we need to interchange rows?
			     Yes, do so ...        */
      for (k = 0; k < n; k++) {
	dum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k]    = dum;
      }
      *d = -(*d);       /* ... and change the parity of d.    */
      vv[imax] = vv[j]; /* Also interchange the scale factor. */
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)  a[j][j] = TINY;
    if (j != (n-1)) {       /* Now, finally, divide by the pivot element.*/
      dum = 1.0/(a[j][j]);
      /* If the pivot element is zero the matrix is singular (at least
	 to the precision of the algorithm).
	 For some applications on singular matrices, it is desirable to
	 substitute TINY for zero. */
      for (i = j+1; i < n; i++)
	a[i][j] *= dum;
    }
  }       /* Go back for the next column in the reduction.        */

  Lib::free_dvector(vv, 0, n-1);

  return GOOD;
}

/* Solves the set of n linear equations A*X = B. Here a[1..n][1..n] is input,
   not as the matrix A but rather as its LU decomposition, determined by the
   routine ludcmp. indx[1..n] is input as the permutation vector returned by
   ludcmp. b[1..n] is input as the right-hand side vector B, and returns with
   the solution vector X. a,n, and indx are not modified by this routine and
   can be left in place for successive calls with different right-hand sides
   b. This routine takes into account the possibility that b will begin with
   many zero elements, so it is efficient for use in matrix inversion.
*/
void Lib::lubksb(double **a, int n, int *indx, double b[]){
  int i, ii=-1, ip, j;
  double sum;

  for (i = 0; i< n; i++) {
    ip   = indx[i];
    sum  = b[ip];
    b[ip]= b[i];
    if (ii+1) {
      for (j =ii; j <=i-1; j++)
	sum -=a[i][j]*b[j];
    }
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--) {
    sum = b[i];
    for (j = i+1; j < n; j++) 
      sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}

/* Standard error handler.
*/
void Lib::nrerror(char *error_text){
  PROBLEM "Run-time error...\n%s\n...now exiting to system...\nLib.cpp key 365\n",
    error_text RECOVER(NULL_ENTRY);
}

/* Get the items for first derivatives, return a double value.
*/
double Lib::GetFderItem(int t_flag, int n, double **mat_v, double **mat_x,
			double *y_vec) {
  //    int i,j;
  double rv = 0;
  double **temp_mat1, **temp_mat2, **mat_out;

  temp_mat1 = dmatrix(0, n-1, 0, n-1);
  temp_mat2 = dmatrix(0, n-1, 0, n-1);
  mat_out   = dmatrix(0, n-1, 0, n-1);

  switch(t_flag) {
  case 1:
    Multi_Matrices(mat_out, n, mat_v,     mat_x);
    rv = (-1.0/2.0) * Trace_Matrix(n, mat_out);
    break;
  case 2:
    Multi_Matrices(temp_mat1, n, mat_v,     mat_x);
    Multi_Matrices(mat_out, n, temp_mat1, mat_v);
    rv = (1.0/2.0) * Multi_Vectors_Matrix(n, y_vec, mat_out, y_vec);
    break;
  default:
    PROBLEM "The variable 'int t_flag' needs to have the value 1 or 2.\nLib.cpp key 368\n"
      RECOVER(NULL_ENTRY);
  }
  free_dmatrix(temp_mat1, 0, n-1, 0, n-1);
  free_dmatrix(temp_mat2, 0, n-1, 0, n-1);
  free_dmatrix(mat_out, 0, n-1, 0, n-1); 
  return rv;
}

/* Get the items for second derivatives according to the 
   scoring-algorithm, return a double value.
*/
double Lib::GetExSder(int n, double **mat_v, double **mat_x, double **mat_y) {
  //    int i,j;
  double rv;
  double **temp_mat1, **temp_mat2, **mat_out;

  temp_mat1 = dmatrix(0, n-1, 0, n-1);
  temp_mat2 = dmatrix(0, n-1, 0, n-1);
  mat_out   = dmatrix(0, n-1, 0, n-1);

  Multi_Matrices(temp_mat1, n, mat_v,     mat_x);
  Multi_Matrices(temp_mat2, n, temp_mat1, mat_v);
  Multi_Matrices(mat_out, n, temp_mat2, mat_y);
  rv = (1.0/2.0) * Trace_Matrix(n, mat_out);

  free_dmatrix(temp_mat1, 0, n-1, 0, n-1);
  free_dmatrix(temp_mat2, 0, n-1, 0, n-1);
  free_dmatrix(mat_out, 0, n-1, 0, n-1);
  return rv;
}

/*****************************************************************************/
double **Lib::dVecToMat(double *Yvec,  
			Sint nrl,
			Sint nrh,
			Sint ncl,
			Sint nch){
  /**************************************************************************
   * allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]     * 
   * and copy the values of the vector Yvec into the allocated double matrix.* 
   **************************************************************************/
  Sint i,j,k;
  double **Y;

  Y=dmatrix(nrl,nrh,ncl,nch);
  k=0;
  for (j=ncl;j<=nch;j++){
    for (i=nrl;i<=nrh;i++){
      Y[i][j]=Yvec[k];
      k++;
    }
  }
  return Y;
}
/*****************************************************************************/
void Lib::dMatToVec(double **Ymat,
		     double *Yvec,  
		     Sint nrl,
		     Sint nrh,
		     Sint ncl,
		     Sint nch){
  /**************************************************************
* copy the values of a double matrix with subscript range     *
* Ymat[nrl..nrh][ncl..nch] into a double vector Yvec.         *
**************************************************************/
  Sint i,j,k;

  k=0;
  for (j=ncl;j<=nch;j++){
    for (i=nrl;i<=nrh;i++){
      Yvec[k]=Ymat[i][j];
      k++;
    }
  }
}
/*****************************************************************************/
int **Lib::iVecToMat(int *Yvec,  
		     Sint nrl,
		     Sint nrh,
		     Sint ncl,
		     Sint nch){
  /*************************************************************************
   * allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]       *
   * and copy the values of the vector Yvec into the allocated int matrix.  *
   *************************************************************************/
  Sint i,j,k;
  int **Y;

  Y=imatrix(nrl,nrh,ncl,nch);
  k=0;
  for (j=ncl;j<=nch;j++){
    for (i=nrl;i<=nrh;i++){
      Y[i][j]=Yvec[k];
      k++;
    }
  }
  return Y;
}
/*****************************************************************************/
void Lib::iMatToVec(int **Ymat,
		     int *Yvec,  
		     Sint nrl,
		     Sint nrh,
		     Sint ncl,
		     Sint nch){
  /**************************************************************
* copy the values of a int matrix with subscript range        *
* Ymat[nrl..nrh][ncl..nch] into a int vector Yvec.            *
**************************************************************/
  Sint i,j,k;

  k=0;
  for (j=ncl;j<=nch;j++){
    for (i=nrl;i<=nrh;i++){
      Yvec[k]=Ymat[i][j];
      k++;
    }
  }
}

/* Allocates a double matrix with range [nrl..nrh][ncl..nch].
*/
double** Lib::dmatrix(int nrl, int nrh, int ncl, int nch){
    int i;
    double **m;

    /* Allocate pointers to rows.   */
    m = (double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) {
      char errorString[] = "allocation failure for rows in dmatrix()\n";
      nrerror(&errorString[0]);
    }
    m -= nrl;

    /* Allocate rows and set pointers to them.   */
    for (i = nrl; i <= nrh; i++) {
	m[i] = (double *) malloc((unsigned)(nch-ncl+1)*sizeof(double));
	if (!m[i]) {
	  char errorString[] = "allocation failure for columns in dmatrix()\n";
	  nrerror(&errorString[0]);
        }
        m[i] -= ncl;
    }
    /* Return pointer to array of pointers to rows.   */
    return m;
}

/* Allocates an int matrix with range [nrl..nrh][ncl..nch].
*/
int** Lib::imatrix(int nrl, int nrh, int ncl, int nch){
    int i;
    int **m;

    /* Allocate pointers to rows.   */
    m = (int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if (!m) {
      char errorString[] = "allocation failure for rows in imatrix()\n";
      nrerror(&errorString[0]);
    }
    m -= nrl;

    /* Allocate rows and set pointers to them.   */
    for (i = nrl; i <= nrh; i++) {
	m[i] = (int *) malloc((unsigned)(nch-ncl+1)*sizeof(int));
	if (!m[i]) {
	  char errorString[] = "allocation failure for columns in imatrix()\n";
	  nrerror(&errorString[0]);
        }
        m[i] -= ncl;
    }
    /* Return pointer to array of pointers to rows.   */
    return m;
}

/* Allocates a int vector with range [nl..nh].
*/
int* Lib::ivector(int nl, int nh){
    int  *v;

    v = (int *)malloc((unsigned)(nh-nl+1)*sizeof(int));
    if(!v) {
      char errorString[] = "allocation failure in ivector()";
      nrerror(&errorString[0]);
    }
    return v-nl;
}


/* Allocates a double vector with range [nl..nh].
*/
double* Lib::dvector(int nl, int nh){
    double *v;

    v = (double *)malloc((unsigned)(nh-nl+1)*sizeof(double));
    if (!v) {
      char errorString[] = "allocation failure in dvector()";
      nrerror(&errorString[0]);
    }
    return v-nl;
}
/* Frees a int vector allocated by ivector().
*/

void Lib::free_ivector(int *v, int nl, int nh){
    free((char*)(v+nl));
}


/* Frees a double vector allocated by dvector().
*/
void Lib::free_dvector(double *v, int nl, int nh){
    free((char*)(v+nl));
}

/* Frees a matrix allocated with dmatrix.
*/
void Lib::free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch){
    int i;

    for (i = nrh; i >=nrl; i--)
	free((char *)(m[i]+ncl));
    free((char *)(m+nrl));
}

/* Frees a matrix allocated with imatrix.
*/
void Lib::free_imatrix(int **m, int nrl, int nrh, int ncl, int nch){
    int i;

    for (i = nrh; i >=nrl; i--)
	free((char *)(m[i]+ncl));
    free((char *)(m+nrl));
}
