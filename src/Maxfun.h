#ifndef MAX_FUN_H
#define MAX_FUN_H

#include "Maxfun_fun.h"
#include "Calculate.h"
#include "multic.h"
#include <fstream>
using namespace std;

class Maxfun : public Maxfun_fun, public Likelihoodfun {
  int need_ibd_flag;
  int method_flag, proband_flag;
  ifstream maxfunfile;
  double x_deta[INITVCNUM], x_gama[INITVCNUM];
  Maxfun_exp Exp;
  void init(TraitMarkerCov_par *tmc, InitValue_par *init);
 public:
  Maxfun(TraitMarkerCov_par *tmc, InitValue_par *init, FortData *fa, int ds,
	 ShareRelation *sa, int rs);
  virtual ~Maxfun();
  virtual Maxfun_exp getExp();
  virtual void eval_fun(double *TR,double *FQE,int *NFE,int *LEX);
  virtual void dep_fun(double *tr,int *lex);
  //scorefun
  void   OpenFiles();
  void   CloseFiles();
  void   getamp(int *a_flag, int *m_flag, int *p_flag);
  double Like_func();
  //The new model component for Maxfun only 
  void Get_z_mat(int n_fam_member,double **z1_mat,double **z2_mat,double **z3_mat);
  void ChangeMatrices(double **g_mat, double **z1_mat,double **z2_mat,double **z3_mat,
		      double **i_mat,double **c_mat, double **p_mat, double **q_mat,
		      int n_fam_member, int *n_member);
  void ShrinkMatrices(double **g_mat, double **z1_mat,double **z2_mat,double **z3_mat,
		      double **i_mat,double **c_mat,double **p_mat,double **q_mat,
		      int n_order,int n_member);
  double SRun();
  void   PrintErrMsg(const char *str_mat,const char *str_routine,int mat_dim,
		     double **x_mat);
  //ascert 
  double Correct_factor();
  double erf(double x);
  double dens(double x);
  double cump(double x);

  double CF_1or2prob();
  void GetProb(int n_fam_member);
  void GetProb_z_mat(int n_fam_member, double **z1_mat, double **z2_mat);
  void GetProb_gcpq_mat(int n_fam_member,  double **g_mat, 
			double **c_mat, double **p_mat, double **q_mat);

  void Build_D_Matrix(int n_fam_member, double **D_mat,
		      double *d_vec[]);
  void GetMuBetaVector(double *MuBeta_vec[],double *i_vec,
		       double *d_vec[],int n_member);
  //new interface added for robust variance
  virtual double getmajorgene1();
  virtual double getmajorgene2();
  virtual double getpolygene();
  virtual double getenvironment();
  virtual double getC();
  virtual double getP();
  virtual double getQ();
};

#endif
