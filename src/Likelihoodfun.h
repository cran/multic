#ifndef LIKELIHOODFUN_H
#define LIKELIHOODFUN_H

#include <fstream>
#include "Estimator.h"
#include "multic.h"
using namespace std;

class Likelihoodfun : public Estimator{
 protected:
  // Added by Eric Lunde on 7-31-03 to represent fort.12 internally in memory
  FortData *fortArray;
  int dataSize;
  int dataIndex;

  // Added by Eric Lunde on 7/29/03 to represent share.out internally in memory
  ShareRelation *shareArray;
  int relationSize;
  int relationIndex;

  int datatype;   // normal or longitudinal data type.
  int env_vcnum;
  int iloci,itraits,icovs, iinitvcnum, imubeta_dim, ismtcpq_dim,
    itraitcovnum, imubetacovnum;
  int irepeatmst; // number of repeat measurement for longitudinal datatype.
  int total_trait_values, total_cov_values; // added for multic and longic data
  int* t_start;
  int* t_len; // For multiple traits.
  int* m_start;
  int* m_len;  // For multiple markers.
  int* c_start;
  int* c_len;   // For multiple covariates.
  int  nfam;    // number of family
  int* N;    // member of each family
  int N1;    // number of iteration N1.
  int init_mu_flag, init_s_flag,  init_m1_flag, init_m2_flag, init_t_flag;
  int init_c_flag,  init_p_flag,  init_q_flag;
  /* The flags are for fixed or estimated  
     initial values, static varables.    */
  int mu_flag, s_flag,  m1_flag, m2_flag, t_flag; /* working varables.*/
  int c_flag,  p_flag,  q_flag;                   /* working varables.*/
  int s_flag_array[INITVCNUM],  m1_flag_array[INITVCNUM],
    m2_flag_array[INITVCNUM], t_flag_array[INITVCNUM];
  int c_flag_array[INITVCNUM],  p_flag_array[INITVCNUM],
    q_flag_array[INITVCNUM];
  /* working varables.
     The flags are for fixed or estimated
     parameters for all variances and
     covariances components.             */
  int s_fix_flag_array[INITVCNUM],  m1_fix_flag_array[INITVCNUM],
    m2_fix_flag_array[INITVCNUM], t_fix_flag_array[INITVCNUM];
  int c_fix_flag_array[INITVCNUM],  p_fix_flag_array[INITVCNUM],
    q_fix_flag_array[INITVCNUM];
  /* working varables.
     The flags are for fixed internally
     or externally.                      */
  int s_vcnum,  m1_vcnum, m2_vcnum, t_vcnum; // number of variance components.
  int c_vcnum,  p_vcnum,  q_vcnum;           // number of variance components.
  double* init_x_mu;
  double* init_x_S;
  double* init_x_M1;
  double* init_x_M2;
  double* init_x_T;
  double* init_x_c;
  double* init_x_p;                        
  double* init_x_q;
  double  init_x_beta[MAXNUMTRAIT][MAXNUMCOV];
  /* initial values of MU,SIGMA,TAU,MG,C,P,Q
     and covariates for multiple traits. */
  double x_mu[MAXNUMTRAIT], x_S[INITVCNUM], x_M1[INITVCNUM], x_M2[INITVCNUM],
    x_T[INITVCNUM],    x_beta[MAXNUMTRAIT][MAXNUMCOV];
  double x_c[INITVCNUM],    x_p[INITVCNUM], x_q[INITVCNUM]; 
  /* working values of MU,SIGMA,TAU,MG,C,P,Q
     and covariates for multiple traits. */ 


  int header;
  int printfirstfam_flag;  // print first family's data into lookup log file.
  //  int ascert_flag;         // Ascertainment flag; 1 = ascert, 0 = unascert.

  Trait*     trait_array;
  Marker*    marker_array;
  Cov*       cov_array;
        
  Family fam;
  ifstream  fp_loci;                   /* pi values from ibd part   */
  ofstream  fp_lookuplog;              /* the intermediate values*/

 public:
  Likelihoodfun();
  virtual ~Likelihoodfun();
  // common block for constructing likelihood function
  //              void GetFromPara();
  //              void PrintMarkerErrorMg(fpos_t *curfam_ptr, fpos_t *nextfam_ptr);
  //              void PrintMatrix(int n_fam_member, double **g_mat, double **z1_mat,
  //                      double **z2_mat,  double **i_mat);
  void PrintResults(int debug_flag);
  void ChangeMatrices(double **g_mat, double **z1_mat,  double **z2_mat,
		      double **i_mat, double **c_mat,   double **p_mat,
		      double **q_mat, int n_fam_member, int  *n_member);
  void ShrinkMatrices(double **g_mat, double **z1_mat, double **z2_mat,
		      double **i_mat, double **c_mat, double **p_mat,
		      double **q_mat, int n_order,    int n_member);
  void GetVGZICPQMatrix(double **V_mat,   double **G_mat[], double **Z1_mat[],
			double **Z2_mat[], double **I_mat[], double **C_mat[],
			double **P_mat[],  double **Q_mat[], double **g_mat,
			double **z1_mat,  double **z2_mat,   double **i_mat,
			double **c_mat,   double **p_mat,   double **q_mat,
			int n_fam_member);
  void GetSmtcpqMatrix(double **Smtcpq_mat[], double **G_mat[],
		       double **Z1_mat[], double **Z2_mat[], double **I_mat[],
		       double **C_mat[], double **P_mat[], double **Q_mat[],
		       int n_fam_member);
  void GetXValues(int numval, double x_val[], char x_str[]);
  void InitMissingTraitflag(int n_fam_member);
  void InitMat(int n_fam_member,double **G_mat[],double **Z1_mat[],
	       double **Z2_mat[], double **I_mat[],double **g_mat,
	       double **z1_mat,   double **z2_mat, double **i_mat); 
  void InitVec(int n_fam_member, double *i_vec,double *ydb_vec,
	       double *MuBeta_vec[]);
  void InitCPQMat(int n_fam_member,double **C_mat[],double **P_mat[],
		  double **Q_mat[], double **c_mat, double **p_mat,
		  double **q_mat);
  void GetStartMubeta_SMT();
  void GetStartCPQ();
  void GetInitValues();
  void GetST_flag_array();
  void GetM_flag_array();
  void GetCPQ_flag_array();
  void GetSMT_vcnum();
  void GetCPQ_vcnum();
  void GetInitCovValues(char x_str[]);
  void Build_beta_Vector(int nt, double *beta_vec);
  void GetData(int n_fam_member);
  void SaveProband(int family_member_number);
  void Get_z_mat(int n_fam_member, double **z1_mat, double **z2_mat);
  void Get_gcpq_mat(int n_fam_member,
		    double **g_mat, double **c_mat, double **p_mat,
		    double **q_mat);
  void AssignInitValues(); 

  // longic
  void LGetMuBetaVector(double *MuBeta_vec[], double *i_vec, double ***ld_vec,
			int n_member);
  void Build_lD_Matrices(int n_fam_member, double ***D_mat, double ***d_vec);
  void LChangeVec(double ***ld_vec,int n_fam_member);

  virtual double getmajorgene1()=0;
  virtual double getmajorgene2()=0;
  virtual double getpolygene()=0;
  virtual double getenvironment()=0;
  virtual double getC()=0;
  virtual double getP()=0;
  virtual double getQ()=0;
};	

#endif
