#ifndef CALCULATE_H
#define CALCULATE_H

#include <fstream>
#include "Likelihoodfun.h"
#include "Least.h"
#include "TraitMarkerCov_par.h"
#include "InitValue_par.h"
#include "multic.h"

static char multic_para_file[] = "multic.par";
static char data_file[] = "fort.12";     /* family data file name*/
static char loci_file[] = "loci.out";    /* marker loci file name*/
static char share_file[] = "share.out";   /* G,C,PP,PO file name  */
static char out_file[]  = "multic.out";
static char mu_file[]   = "multic.mu";
static char s_file[]    = "multic.poly";
static char m1_file[]   = "multic.mg1";
static char m2_file[]   = "multic.mg2";
static char t_file[]    = "multic.env";
static char sib_file[]  = "multic.sib";
static char pp_file[]   = "multic.pp";
static char po_file[]   = "multic.po";
static char lik_file[]  = "multic.lik";
static char log_file[]  = "multic.log";
static char lookuplog_file[]  = "multiclookup.log";
static char summary_file[] = "summary.log";
static char iteration_file[] = "iterations.log";

// SUMMARY_LOG_FIELD_WIDTH was added to centraliize the value used to specify
// field with of the file summary.log.  Eric Lunde 8-22-03
#define SUMMARY_LOG_FIELD_WIDTH 25

// FAMILY_ID_MAX_LENGTH was added to centraliize the value used to specify
// the maximum length of the strings used to hold the family ids.  Family and 
// person ids were converted to chars from ints by Eric Lunde on 11-9-2004
#define FAMILY_ID_MAX_LENGTH 32

class Calculate : public Likelihoodfun {
 protected:
  // These values are added to search for NaN values, so they are not output
  // to the v matrix temporary files. - EML 09-29-2004  
  // These values were replaced by float.h's DBL_MAX and DBL_MIN
  // - EML 07-04-2005
  //int MAX_DOUBLE_BINARY, MIN_DOUBLE_BINARY;
  //double MAX_DOUBLE, MIN_DOUBLE;

  // savedVMatrix and savedYBetaDiff were added so we can report to the user
  // the values in the V Matrix and the difference between a person's trait 
  // values and the estimated betas at the time of convergence.  Since at the 
  // location of the calculation of such values, we do not know whether or not
  // convergence will occur, we save the values in the addressed by these
  // variables.  Then after convergence, we can write the values to the files
  // necessary.  Eric Lunde 01-08-04
  double ***savedVMatrix;
  // savedYBetaDiff stands for the difference between the Y's (individual
  // person's trait value) and the beta's (estimated values)
  double ***savedYBetaDiff;
  // Since the amount of trait valued members per family changes with each
  // family, we also need a storage of the number of trait valued members for
  // each family.  savedTraitValuedMembers saves that information.
  int *savedTraitValuedMembers;

  // estimates was added because we want to save the estimates from the null
  // hypothesis back in multic.s to be used as initial values for the
  // alternative hypotheses.  We will pass the address of that array all the
  // way into Calculate.cpp and store the results here.
  double *estimates;
  // The correct order for the values in estimates is the same order that the
  // values are printed to multic.log.  So we just keep an index into the
  // array, and append all new values.
  int estimatesIndex;

  // coefficients was added because we want to save the coefficients from the
  // null hypothesis back in multic.s to be used as initial values for the
  // alternative hypotheses.  We will pass the address of that array all the
  // way into Calculate.cpp and store the results here.
  double *coefficients;
  // The correct order for the values in coefficients is the same order that
  // the values are printed to multic.log.  So we just keep an index into the
  // array, and append all new values.
  int coefficientsIndex;
  
  // familyLoglikelihoods was created to be able to capture the loglikelihood
  // for each of the families
  double *familyLoglikelihoods;

  // uniqueFamilyIds is an char * array to hold the unique family identifers
  // for quick lookup in writeSavedVMatrix
  char **uniqueFamilyIds;

  // ibdFileName was added by Eric Lunde on 7-11-03 so we can print the name of
  // the (m)ibd file that genereated the output
  char ibdFileName[BUF];
  Least *least;
  /* -------- Global variables for the program  ------- */
  int total_missing_trait[MAXNUMTRAIT];
  int total_missing_marker[MAXNUMMARKER];

  /* The parameters saved for PrintSummary 
     under Alternative hypotheses H(A).  */
  int  s_A_flag_array[INITVCNUM],  m1_A_flag_array[INITVCNUM],
    m2_A_flag_array[INITVCNUM], t_A_flag_array[INITVCNUM];
  int  c_A_flag_array[INITVCNUM],  p_A_flag_array[INITVCNUM],
    q_A_flag_array[INITVCNUM];

  /* The parameters saved for PrintSummary 
     under Null hypotheses H(0).         */
  int  s_0_flag_array[INITVCNUM],  m1_0_flag_array[INITVCNUM],
    m2_0_flag_array[INITVCNUM], t_0_flag_array[INITVCNUM];
  int  c_0_flag_array[INITVCNUM],  p_0_flag_array[INITVCNUM],
    q_0_flag_array[INITVCNUM];

  /* The parameters saved for PrintSummary 
     under Alternative hypotheses H(A).  */
  int  s_A_fix_flag_array[INITVCNUM],  m1_A_fix_flag_array[INITVCNUM],
    m2_A_fix_flag_array[INITVCNUM],  t_A_fix_flag_array[INITVCNUM];
  int  c_A_fix_flag_array[INITVCNUM],  p_A_fix_flag_array[INITVCNUM],
    q_A_fix_flag_array[INITVCNUM];

  /* The parameters saved for PrintSummary 
     under Null hypotheses H(0).         */
  int  s_0_fix_flag_array[INITVCNUM],  m1_0_fix_flag_array[INITVCNUM],
    m2_0_fix_flag_array[INITVCNUM],  t_0_fix_flag_array[INITVCNUM];
  int  c_0_fix_flag_array[INITVCNUM],  p_0_fix_flag_array[INITVCNUM],
    q_0_fix_flag_array[INITVCNUM];

  //  int runopt_flag;                      /* run the hypotheses with no major gene or not. */ 
  //  int run_flag;                         /* the flag for first or second run.   */ 
  int final_flag;                       /* final value flag.                   */ 
  int converg_flag;                     /* convergent flag for the first run.  */ 
  //  int need_ibd_flag;                    /* the flag  for using IBD's pi values or not.  */ 
  int need_mg1_flag, need_mg2_flag;     /* the flags for using IBD's pi values of mg1 and mg2  */ 
  //  int fixtoboundary_flag;
                                        /* the flag  for fixing parameters to
  					   boundary or not.                    */ 
  //  int algorithm_flag;
                                       /* specify algorithm used: (added by chen 10/1/99)
					 1. Least square
					 2. Multic;
					 3. Multic with Least square as initial values
					 4. Robust variance;
					 5. Maxfun;
				      */
  int breakstephalf_flag;               /* break the step halfing within the 
					   program.                            */ 
  int breakstephalf_run1;               /* break the step halfing within the 
					   program for the first run.          */ 
  int breakstephalf_run2;               /* break the step halfing within the 
					   program for the second run.         */ 

  double **inv_Exsder_mubeta_mat, **inv_Exsder_smtcpq_mat;   /* No ascertainment.      */
  double **inv_Exsder_mat;                                /* Ascertainment.         */   

  double x_mu_A[MAXNUMTRAIT], x_S_A[INITVCNUM], x_M1_A[INITVCNUM], x_M2_A[INITVCNUM], 
    x_T_A[INITVCNUM],    x_beta_A[MAXNUMTRAIT][MAXNUMCOV];   
  double x_c_A[INITVCNUM],    x_p_A[INITVCNUM], x_q_A[INITVCNUM]; 
  /* The values of MU,SIGMA,TAU,MG,C,P,Q  
     covariates under Alternative
     hypotheses.                         */
  double x_mu_0[MAXNUMTRAIT], x_S_0[INITVCNUM], x_M1_0[INITVCNUM], x_M2_0[INITVCNUM], 
    x_T_0[INITVCNUM],    x_beta_0[MAXNUMTRAIT][MAXNUMCOV];   
  double x_c_0[INITVCNUM],    x_p_0[INITVCNUM], x_q_0[INITVCNUM]; 
  /* The values of MU,SIGMA,TAU,MG,C,P,Q 
     covariates under Null hypotheses.   */
  double inc_mu[MAXNUMTRAIT], inc_beta[INITCOVNUM], inc_s[INITVCNUM], 
    inc_m1[INITVCNUM],   inc_m2[INITVCNUM],    inc_t[INITVCNUM];          
  double inc_c[INITVCNUM],    inc_p[INITVCNUM],     inc_q[INITVCNUM]; 
  /* increased values of MU,SIGMA,TAU,MG,C,P,Q
     before CalculateValues function.    */
  double new_inc_mu[MAXNUMTRAIT], new_inc_beta[INITCOVNUM], new_inc_s[INITVCNUM], 
    new_inc_m1[INITVCNUM],   new_inc_m2[INITVCNUM],    new_inc_t[INITVCNUM];          
  double new_inc_c[INITVCNUM],    new_inc_p[INITVCNUM],     new_inc_q[INITVCNUM]; 
  /* increased values of MU,SIGMA,TAU,MG,C,P,Q
     after CalculateValues function.     */
  double LogLikbreakval;
  double loglikelihood_A, loglikelihood_0;

  double *mubeta_A_sd, *smtcpq_A_sd;              /* standard error vector under
						     H(A).                  */
  double *mubeta_0_sd, *smtcpq_0_sd;              /* standard error vector under
						     H(0).                  */

  ofstream fp_out,                   /* output file              */   
    fp_mu,                    /* output file for mu       */
    fp_s,                     /* output file for x_S(s)   */
    fp_m1,                     /* output file for x_M1(m1) */
    fp_m2,                     /* output file for x_M2(m2) */
    fp_t,                     /* output file for x_T(t)   */
    fp_sib,                   /* output file for x_c(sib) */
    fp_pp,                    /* output file for x_p(pp)  */
    fp_po,                    /* output file for x_q(po)  */
    fp_lik,                   /* log likelihood file      */
    fp_log,                    /* summary for the program  */
    summaryLog; // summary.log file

  long printProgress;
  long calculateResiduals;

  // markerName was added to when multic prints the iteration count to file
  // it can also print the marker that produced that many iterations.
  char *markerName;

 public:
  Calculate(Least *, char *, double, ShareRelation *, int, FortData *, int,
	    double *est, double *coeff, long print_progress,
	    long calculateResiduals, long *familySizes,
	    char **uniqueFamilies, int famliyCount);
  ~Calculate();
  void main_fun(TraitMarkerCov_par *tmc, InitValue_par *init);
  void OpenFiles();
  void CloseFiles();
  void PrintTitle(ostream *fp);
  void PrintParaInfo();
  void PrintInputFiles();
  void SetDataType();
  void GetIbd_flag();
  void PrintSummary(struct tm *now, double time_secs);
  void PrintSummaryP1();
  void PrintSummaryP2();
  void PrintSummaryP3(int H_Aor0_flag);
  void PrintSummaryP3_Result_1(int H_Aor0_flag);
  void PrintSummaryP3_Result_2(int H_Aor0_flag);

  double ConvergCriterion1(double prev_val);
  int    Fixtoboundary();
  void   CheckBoundary();
  void   CalculateDiffOfPara(double *prev_mu,  double **prev_beta,
			     double *prev_s,   double *prev_m1, double *prev_m2,
			     double *prev_t,   double *prev_c,  double *prev_p,
			     double *prev_q, int *converg1_flag,int stephalf_flag);
  double CalculateLik(double *prev_mu, double **prev_beta, double *prev_s,
		      double *prev_m1, double *prev_m2, double *prev_t,
		      double *prev_c,  double *prev_p,  double *prev_q,
		      double prev_ln_func, int    *converg1_2_flag);
  void CheckFixTraitcov();
  void SCondOfInc();
  void MCondOfInc();
  void SMCPQ_Inc(int x_flag, int x_flag_array[], double temp_val[], double x_str[],
		 double inc_x[]);
  //		void GetIbd();
  //		void SaveMeanAndFlag(int run_flag);
  void SaveMeanAndFlag();
  //		void Get_SE(int run_flag, int ascert_flag);
  void Get_SE();
  void Run();
  void SRun();
  void FixCovAndAssignInitValues();
  void Build_D_Matrix(int n_fam_member, double **D_mat, double *d_vec[]);
  void GetMuBetaVector(double *MuBeta_vec[], double *i_vec, double *d_vec[],
		       int n_member);
  void ChangeVec(double *d_vec[],int n_fam_member);

  // common block for constructing likelihood function
  ////		void GetFromPara();
  void PrintMarkerErrorMg(int curfam_ptr, int nextfam_ptr);
  ////		void PrintMatrix(int n_fam_member, double **g_mat, double **z1_mat, 
  ////		 	double **z2_mat,  double **i_mat);
  void PrintResults(int debug_flag);
  //		void ChangeMatrices(double **g_mat, double **z1_mat,  double **z2_mat,   
  //			double **i_mat, double **c_mat,   double **p_mat,   
  //		        double **q_mat, int n_fam_member, int  *n_member);
  //		void ShrinkMatrices(double **g_mat, double **z1_mat, double **z2_mat, double **i_mat, 
  //                        double **c_mat, double **p_mat,  double **q_mat, 
  //                        int n_order,    int n_member);
  //		void GetVGZICPQMatrix(double **V_mat,   double **G_mat[], double **Z1_mat[], double **Z2_mat[], 
  //		      double **I_mat[], double **C_mat[], double **P_mat[],  double **Q_mat[],
  //		      double **g_mat,   double **z1_mat,  double **z2_mat,   double **i_mat,
  //		      double **c_mat,   double **p_mat,   double **q_mat, int n_fam_member);
  //		void GetSmtcpqMatrix(double **Smtcpq_mat[], double **G_mat[], double **Z1_mat[], 
  //			double **Z2_mat[], double **I_mat[], double **C_mat[], double **P_mat[],
  //			double **Q_mat[], int n_fam_member);
  //		void GetXValues(int numval, double x_val[], char x_str[]);
  //		void InitMissingTraitflag(int n_fam_member);
  //		void InitMat(int n_fam_member,double **G_mat[],double **Z1_mat[], double **Z2_mat[], 
  //	     		double **I_mat[],double **g_mat,  double **z1_mat,   double **z2_mat, 
  //	     		double **i_mat);
  //		void InitVec(int n_fam_member, double *i_vec,double *ydb_vec, 
  //	     		double *MuBeta_vec[]);
  //		void InitCPQMat(int n_fam_member,double **C_mat[],double **P_mat[], double **Q_mat[], 
  //			double **c_mat, double **p_mat, double **q_mat);
  //		void GetStartMubeta_SMT();
  //		void GetStartCPQ();
  //		void GetInitValues();
  //		void GetST_flag_array();
  //		void GetM_flag_array();
  //		void GetCPQ_flag_array();
  //		void GetSMT_vcnum();
  //		void GetCPQ_vcnum();
  //		void GetInitCovValues(char x_str[]);
  //		void Build_beta_Vector(int nt, double *beta_vec);
  //		void GetData(int n_fam_member);
  //		void SaveProband(int num_member, char str[]);
  //		void Get_z_mat(int n_fam_member, double **z1_mat, double **z2_mat); 
  //		void Get_gcpq_mat(int n_fam_member, double **g_mat, double **c_mat, 
  //			double **p_mat,   double **q_mat);
  //		void AssignInitValues();

  void GetDim();
  void Getismtcpq_dim();

  // ascertainment
  void AS_SGetExSderMatrix(double **mat, double exsder_mubeta[], double exsder_smtcpq[],
			   double exsder_cfmubeta[], double exsder_cfmusmtcpq, double exsder_cfsmtcpq);
  void AS_SGetFderVector(double *vec, double fder_mubeta[], double fder_s[], 
			 double fder_m1[], double fder_m2[],double fder_t[], double fder_cfmubeta[], 
			 double fder_cfsmtcpq, double fder_c[],double fder_p[],double fder_q[]); 
  double AS_CalculateValues();
  void AS_GetIncValues(double *vec);
  //		void AS_GetIncValues(double *vec,    double inc_mu[], double inc_beta[], 
  //		     	double inc_s[], double inc_m1[], double inc_m2[], 
  //		     	double inc_t[], double inc_c[],  double inc_p[], double inc_q[]);

  // no ascertainment
  void SGetExSderMatrix(int mu_smtcpqflag, double **mat, double exsder_array[]);
  void SGetFderVector(double *vec, double fder_s[], double fder_m1[], double fder_m2[],
		      double fder_t[], double fder_c[], double fder_p[], double fder_q[]);
  double NOAS_CalculateValues();
  //		void GetIncValues(double *mubeta_vec, double *smtcpq_vec, double inc_mu[], 
  //			double inc_beta[], double inc_s[], double inc_m1[], double inc_m2[], 
  //			double inc_t[], double inc_c[], double inc_p[], double inc_q[]);
  void GetIncValues(double *mubeta_vec, double *smtcpq_vec);
  void PrintErrMsg(char *str_mat, char *str_routin, int mat_dim, double **x_mat);

  // longic
  //		void LGetMuBetaVector(double *MuBeta_vec[], double *i_vec, double ***ld_vec, int n_member);
  //		void Build_lD_Matrices(int n_fam_member, double ***D_mat, double ***d_vec);
  //		void LChangeVec(double ***ld_vec,int n_fam_member);

  //new interface added for robust variance
  double getmajorgene1();
  double getmajorgene2();
  double getpolygene();
  double getenvironment();
  //shared sibship, offspring-offspring
  double getC();
  //shared spouse, parent-parent
  double getP();
  //shared parent-offspring
  double getQ();
  void writeInvExpSecDer();
  void printIterationNumber(int iterationNumber);

  // V matrix functions
  double ***allocateSavedVMatrix(int familyCount, int totalTraitValues,
				 int *familySizes);
  void saveVMatrix(double ***localSavedVMatrix, double **V_mat,
		   int familyIndex, int totalTraitValues, int familySize);
  void writeSavedVMatrixToFile(double ***localSavedVMatrix, int familyCount,
			       int totalTraitCount, int *familySizes,
			       char *localIbdFileName);
  void freeSavedVMatrix(double ***localSavedVMatrix, int familyCount,
			int totalTraitValues, int *familySizes);

  // (y-beta) vector functions
  double ***allocateSavedYBetaDiff(int familyCount, int totalTraitValues,
				   int *familySizes);
  void saveYBetaDiff(double ***localSavedYBetaDiff,
		     double **y_vec, int familyIndex,
		     int totalTraitValues, int familySize);
  void writeSavedYBetaDiffToFile(double ***localSavedYBetaDiff,
				 int familyCount,
				 int totalTriatCount,
				 int *familySizes,
				 char *localIbdFileName);
  void freeSavedYBetaDiff(double ***localSavedYBetaDiff, int familyCount,
			  int totalTraitValues);
};
	
#endif
