#ifndef LIB_H
#define LIB_H

#define GOOD 1.0
#define BAD  0.0
#define TINY 1.0e-20

class Lib {
 public:
  // Statistics part
  //static double mean(int data[]);
  //static double mean(double data[]);

  // Linear algebra part
  static void   PrintOneVector(int n, double *vec);
  static void   PrintOneVector(int n, int *vec);
  static void   PrintOneMatrix(int nr, int nc, double **a);
  static double DetermOfMatrix(int n,  double **a);
  static void   CopyMatrix(double **a_out, int n, double **a_in);
  static double InverseOfMatrix(double **y, int n, double **a);
  static void   Transpose(double **mat_out,double **mat1,int c1,int c2);
  static void   Multi_Matrices(double **mat_out, int n, double **mat1, double **mat2);
  static int    Multi_Matrices(double**, double**, int, int, double**, int, int);
  static double Multi_Vectors(int n, double *vec1,double *vec2);
  static void   Multi_Matrix_Vector(double *vec_out, int nr, int nc, double **mat,
				    double *vec);
  static double Multi_Vectors_Matrix(int n, double *vec1, double **mat,
				     double *vec2);
  static double Trace_Matrix(int n, double **mat);
  //add by chen
  static double Gausslud(double **a, int n, int *iperm);

  static double ludcmp(double **a, int n, int *indx, double *d);
  static void   lubksb(double **a, int n, int *indx, double b[]);

  static double GetFderItem(int t_flag, int n, double **mat_v, double **mat_x, 
			    double *y_vec);
  static double GetExSder(int n, double **mat_v, double **mat_x, double **mat_y); 
  static void   nrerror(char error_text[]);

/******************************************************************************
When calling externally defined C functions within S-Plus there is a need to
pass data from S-Plus to the C function.

You can only pass vector data structures as arguments from S-Plus to a 
externally defined C routine.

A developer might need to perform the following tasks on the data type(s) - 
double, int, float, long.

1) Dynamically allocate a vector.
2) Dynamically allocate a matrix.
3) Copy the contents of a vector to a matrix. 
4) Copy contents of a matrix to a vector.

We have created a series of functions useful for allocating and manipulating
vector and matrix data structures in C.

The vectorS family of functions return an instance of a vector.
-- dvectorS(),ivectorS(),fvectorS(),lvectorS()

The matrixS family of functions return an instance of a matrix.
-- dmatrixS(),imatrixS(),fmatrixS(),lmatrixS()

The VecToMat family of functions returns an instance of a matrix and copies 
the values of a vector into the allocated  matrix.
-- dVecToMat(),iVecToMat(),fVecToMat(),lVecToMat()

The MatToVec family of functions returns void and copies the values of a 
matrix into a vector.
-- dMatToVec(),iMatToVec(),fMatToVec(),lMatToVec()


NOTE: MatToVec family functions return void and VecToMat family functions 
return a matrix.  

EXAMPLES:

Create int vector i
--------------------
int *i;
i=ivectorS(0,9);

Create long matrix Y
---------------------
float **Y;
Y=fmatrixS(0,4,0,5);

Copy contents of double vector d to double matrix D
----------------------------------------------------
double *d,**D;
D=dVecToMat(d,0,4,0,1); 

Copy contents of float matrix F to float vector Fvec
-----------------------------------------------------
float *Fvec,**F;
fMatToVec(F,Fvec,0,4,0,1);


******************************************************************************/

  static double *dvectorS(long nl, long nh);
  /* allocate a double vector with subscript range v[nl..nh] */

  static float *fvectorS(long nl, long nh);
  /* allocate a float vector with subscript range v[nl..nh] */

  static int *ivectorS(long nl, long nh);
  /* allocate an int vector with subscript range v[nl..nh] */

  static long *lvectorS(long nl, long nh);
  /* allocate an long vector with subscript range v[nl..nh] */

  static double **dmatrixS(long nrl, long nrh, long ncl, long nch);
  /* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */

  static float **fmatrixS(long nrl, long nrh, long ncl, long nch);
  /* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */

  static int **imatrixS(long nrl, long nrh, long ncl, long nch);
  /* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */

  static long **lmatrixS(long nrl, long nrh, long ncl, long nch);
  /* allocate a long matrix with subscript range m[nrl..nrh][ncl..nch] */

  static double **dVecToMat(double *Yvec, long nrl,long nrh,long ncl,long nch);
  /***************************************************************************
   * allocate a double matrix with subscript range m[nrl..nrh][ncl..nch]      *
   * and copy the values of the vector Yvec into the allocated double matrix. *
   ***************************************************************************/

  static long **lVecToMat(long *Yvec, long nrl,long nrh,long ncl,long nch);
  /*************************************************************************
   * allocate a long matrix with subscript range m[nrl..nrh][ncl..nch]      *
   * and copy the values of the vector Yvec into the allocated long matrix. *
   *************************************************************************/

  static float **fVecToMat(float *Yvec, long nrl,long nrh,long ncl,long nch);
  /**************************************************************************
   * allocate a float matrix with subscript range m[nrl..nrh][ncl..nch]      * 
   * and copy the values of the vector Yvec into the allocated float matrix. *
   **************************************************************************/

  static int **iVecToMat(int *Yvec, long nrl,long nrh,long ncl,long nch);
  /*************************************************************************
   * allocate a int matrix with subscript range m[nrl..nrh][ncl..nch]       *
   * and copy the values of the vector Yvec into the allocated int matrix.  *
   *************************************************************************/

  static void dMatToVec(double **Ymat,double *Yvec,long nrl,long nrh,long ncl,
		  long nch);
  /**************************************************************
   * copy the values of a double matrix with subscript range     *
   * Ymat[nrl..nrh][ncl..nch] into a double vector Yvec.         *
   **************************************************************/

  static void fMatToVec(float **Ymat,float *Yvec,long nrl,long nrh,long ncl,
		  long nch);
  /**************************************************************
   * copy the values of a float matrix with subscript range      *
   * Ymat[nrl..nrh][ncl..nch] into a float vector Yvec.          *
   **************************************************************/

  static void iMatToVec(int **Ymat,int *Yvec,long nrl,long nrh,long ncl,
		  long nch);
  /**************************************************************
   * copy the values of a int matrix with subscript range        *
   * Ymat[nrl..nrh][ncl..nch] into a int vector Yvec.            *
   **************************************************************/

  static void lMatToVec(long **Ymat,long *Yvec,long nrl,long nrh,long ncl,
		  long nch);
  /**************************************************************
   * copy the values of a long matrix with subscript range       *
   * Ymat[nrl..nrh][ncl..nch] into a long vector Yvec.           *
   **************************************************************/

static int*   ivector(int nl, int nh);
static int**  imatrix(int nrl, int nrh, int ncl, int nch);
static double* dvector(int nl, int nh);
static double** dmatrix(int nrl, int nrh, int ncl, int nch);
static void   free_ivector(int *v, int nl, int nh);
static void   free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
static void   free_dvector(double *v, int nl, int nh);
static void   free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

};

#endif

