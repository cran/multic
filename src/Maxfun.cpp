#include "Maxfun.h"
#include "Lib.h"
#include "multic.h"
#include <iostream>
#include <cctype>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "Rostream.h"
#include "Rstreambuf.h"

using namespace Rcpp;

extern int ascert_ch;
extern int ascert_flag;

Maxfun::Maxfun(TraitMarkerCov_par *tmc, InitValue_par *initVal,
	       FortData *fa, int ds, ShareRelation *sa, int rs){
  OpenFiles();
  //need_ibd_flag = YES;

  fortArray = fa;
  dataSize = ds;
  dataIndex = 0;

  shareArray = sa;
  relationSize = rs;
  relationIndex = 0;

  init(tmc, initVal);
}

Maxfun::~Maxfun(){
  free(N);
  CloseFiles();
}

Maxfun_exp Maxfun::getExp() {
  return Exp;
}


void Maxfun::init(TraitMarkerCov_par *tmc, InitValue_par *init) {
  int i = 0;
  //read from multic.par
  //  TraitMarkerCov_par tmc(&fp_para);
  //  InitValue_par init(&fp_para,tmc);
        
  itraits = tmc->gettraitnum();
  iloci = tmc->getmarkernum(); 
  icovs = tmc->getcovnum();    
  iinitvcnum = tmc->getinitvcnum();
  env_vcnum = tmc->getinitenvvcnum();
  irepeatmst = tmc->getirepeatmst(); 
  total_trait_values = tmc->gettotal_trait_values();
  total_cov_values = tmc->gettotal_cov_values();
        
  trait_array = tmc->gettrait_array();
  marker_array = tmc->getmarker_array();
  cov_array = tmc->getcov_array();

  t_start = tmc->get_t_start();
  t_len = tmc->get_t_len();
  c_start = tmc->get_c_start();
  c_len = tmc->get_c_len();

  init_mu_flag = init->getmuflag(); 
  init_s_flag = init->getsflag();
  init_m1_flag = init->getm1flag();
  init_m2_flag = init->getm2flag();
  init_t_flag = init->gettflag();
  init_c_flag = init->getcflag();
  init_p_flag = init->getpflag(); 
  init_q_flag = init->getqflag();
  N1 = init->getiterationnum();
  init_x_mu = init->getinitmu();
  init_x_S = init->getinitS();
  init_x_M1 = init->getinitM1();
  init_x_M2 = init->getinitM2();
  init_x_T = init->getinitT();
  init_x_c = init->getinitc();
  init_x_p = init->getinitp();
  init_x_q = init->getinitq();

  init->getinitbeta(init_x_beta);

  if(irepeatmst == 1)
    datatype  = MULTIVARIATE;
  else datatype = LONGITUDINAL;

  nfam = 1;
  for(i=1; i<dataSize; i++) {
    if(strcmp(fortArray[i-1].familyId, fortArray[i].familyId) != 0) {
      nfam++;
    }
  }
  N = (int *) malloc( nfam * sizeof(int) );
  
  // We konw that there has to be at least one person in each family.
  for(i=0; i<nfam; i++) {
    N[i] = 1;
  }
  
  int currentFamily = 0;
  for(i=1; i<dataSize; i++) {
    if(fortArray[i-1].familyId == fortArray[i].familyId) {
      N[currentFamily]++;
    }else {
      currentFamily++;
    }
  }
  
  GetStartMubeta_SMT();
  GetStartCPQ();

  /*
   *--READ IN BASIC CONTROL INPUT FOR MAXFUN
   */
    
  Exp.lprt = 1;

  //read first line
  /* Get flags for ascert_flag:  ascertainment   (1=YES,    2=NO);
     mt_flag:  method          (1=Elston, 2=Hopper);
     p_flag:  proband's number(1 or 2).
  */
  char ascert_ch, str[BUF], str1[BUF];
  maxfunfile.getline(str,BUF);
  maxfunfile >> ascert_ch >> method_flag >> proband_flag >> str1;
  maxfunfile.getline(str,BUF);
  if ((ascert_ch == 'Y') || (ascert_ch == 'y'))
    ascert_flag = 1;
  else if ((ascert_ch == 'N') || (ascert_ch == 'n'))
    ascert_flag = 2;
        
  Exp.maxit = N1;
  Exp.epsc1 = Exp.epsc2 = Exp.epsc3 = Exp.epsd = Exp.yota = Exp.epst=0;
  //	Exp.nt = 8;
  //	Exp.method = 6;
  //	Exp.ixvc = Exp.ihit = 1;
  maxfunfile >> Exp.nt >> Exp.method >> Exp.ixvc >> Exp.ihit >> str1;
  maxfunfile.getline(str,BUF);

  //get lower bound of parameters for maxfun
  int ddd = 0;
  for(i=0;i<Exp.nt;i++){
    maxfunfile >> ddd;
    Exp.thl[i] = ddd;
    //printf("Exp.thl[%d]==%f\n",i,Exp.thl[i]);
  }

  //get upper bound of parameters for maxfun 
  for(i=0;i<Exp.nt;i++){
    maxfunfile >> ddd;
    Exp.thu[i] = ddd;
    //printf("Exp.thu[%d]==%f\n",i,Exp.thu[i]); 
  }

  //get status of parameters for maxfun
  for(i=0; i<MAXFUN_NP; i++) Exp.istin[i] = 4;
  if(init_mu_flag == ESTIMATE)
    Exp.istin[0] = 1;
  else if (init_mu_flag == FIXED)
    Exp.istin[0] = 4;
  if(init_s_flag == ESTIMATE)
    Exp.istin[1] = 1;
  else if (init_s_flag == FIXED)
    Exp.istin[1] = 4;
  if(init_m1_flag == ESTIMATE)
    Exp.istin[2] = 1;
  else if (init_m1_flag == FIXED)
    Exp.istin[2] = 4;
  if(init_m2_flag == ESTIMATE)
    Exp.istin[3] = 1;
  else if (init_m2_flag == FIXED)
    Exp.istin[3] = 4;
  if(init_t_flag == ESTIMATE)
    Exp.istin[4] = 1;
  else if (init_t_flag == FIXED)
    Exp.istin[4] = 4;
  if(init_c_flag == ESTIMATE)
    Exp.istin[5] = 1;
  else if (init_c_flag == FIXED)
    Exp.istin[5] = 4;
  if(init_p_flag == ESTIMATE)
    Exp.istin[6] = 1;
  else if (init_p_flag == FIXED)
    Exp.istin[6] = 4;
  if(init_q_flag == ESTIMATE)
    Exp.istin[7] = 1;
  else if (init_q_flag == FIXED)
    Exp.istin[7] = 4;

  for(i=0;i<icovs;i++)
    Exp.istin[8+i] = 1;

  for(i=0; i<MAXFUN_NP; i++) Exp.stpin[i] = 0.01;

  //get initial value for maxfun
  for(i=0; i<MAXFUN_NP; i++) Exp.thin[i] = 0;
  Exp.thin[0] = init_x_mu[0];
  Exp.thin[1] = init_x_S[0];
  Exp.thin[2] = init_x_M1[0];
  Exp.thin[3] = init_x_M2[0];
  Exp.thin[4] = init_x_T[0];
  Exp.thin[5] = init_x_c[0];
  Exp.thin[6] = init_x_p[0];
  Exp.thin[7] = init_x_q[0];

  for(i=0;i<icovs;i++){
    Exp.thin[i+8] = init_x_beta[0][i];
    //printf("%f\n",Exp.thin[i+8]);
  }

  //   for(i=0; i<Exp.nt; i++)
  //        printf("i=%d %f\n",i, Exp.thin[i]);

  if ((Exp.iout = fopen("fort.12", "w")) == NULL) {
    PROBLEM "The file fort.12 could not be opened for writing.\nMaxFun.cpp key 182\n"
      RECOVER(NULL_ENTRY);
  }

  if ((Exp.idet = fopen("details.txt", "w")) == NULL) {
    PROBLEM "The file details.txt could not be opened for writing.\nMaxFun.cpp key 189\n"
      RECOVER(NULL_ENTRY);
  }
}

void Maxfun::dep_fun(double tr[],int *lex) {
  *lex = 0;
  return;
}

void Maxfun::eval_fun(double TR[], double *FQE, int *NFE, int *LEX){
  int    i;
  //int    a_flag, mt_flag, p_flag;


  for (i = 0; i < 8+icovs; i++)
    Rcout << "TR[" << i << "]=" << TR[i] << std::endl;

  x_mu[0]      = TR[0];
  x_S[0]       = TR[1];
  x_M1[0]      = TR[2];
  x_M2[0]      = TR[3];
  x_T[0]       = TR[4];
  x_c[0]       = TR[5];
  x_p[0]       = TR[6];
  x_q[0]       = TR[7];
  for(i=0;i<icovs;i++){
    x_beta[0][i] = TR[8+i];
  }

  /*
    printf("ascert_flag = %d\n", ascert_flag);
    printf("method_flag = %d\n", method_flag);
    printf("proband_flag = %d\n", proband_flag);
  */
  *FQE  = SRun();

  //printf("FQE=%lf\n",*FQE);

  *NFE += 1;
  *LEX  = 0;
    
  return;
}

/* Get flags for a_flag:  ascertainment   (1=YES,    2=NO);
                 mt_flag:  method          (1=Elston, 2=Hopper);
		 p_flag:  proband's number(1 or 2).
*/
void Maxfun::getamp(int *a_flag, int *mt_flag, int *p_flag2){
  int  aflag = 0, mflag, pflag;

  if ((ascert_ch == 'Y') || (ascert_ch == 'y'))
    aflag = 1;
  else if ((ascert_ch == 'N') || (ascert_ch == 'n'))
    aflag = 2;

  mflag = 2;
  pflag = 2;

  *a_flag = aflag;
  *mt_flag = mflag;
  *p_flag2 = pflag;
}

void Maxfun::OpenFiles()
{
  fp_loci.open(loci_file);
  if (fp_loci.fail()) {
    PROBLEM "The file %s could not be opened for reading.\nMaxFun.cpp key 280\n",
      loci_file RECOVER(NULL_ENTRY);
  }
  maxfunfile.open("maxfun.par");
  if (maxfunfile.fail()) {
    PROBLEM "The file maxfun.par could not be opened for reading.\nMaxFun.cpp key 292\n"
      RECOVER(NULL_ENTRY);
  }
}

void Maxfun::CloseFiles()
{
  fp_loci.close();
  maxfunfile.close();
}

/* Get z1_mat, z2_mat, z3_mat matrices from file fp_loci.
*/
void Maxfun::Get_z_mat(int n_fam_member, double **z1_mat, double **z2_mat,double **z3_mat) {
  int    i,j;
  double mg1_f1, mg1_f2, mg1_pi;
  double mg2_f1, mg2_f2, mg2_pi;
  char   str[BUF];

  for ( i = 0; i < n_fam_member; i++) {
    for ( j = i; j < n_fam_member; j++) {
      if (j > i) {
	fp_loci.getline(str, BUF);
	if (iloci == 1) {
	  sscanf(str,"%lf %lf %lf ", &mg1_f1,&mg1_f2,&mg1_pi);
	  if (mg1_pi == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_pi = 0.5;
	  }
	  if (mg1_f2 == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_f2 = 0.5;
	  }
	  z1_mat[i][j] = mg1_pi;
	  z1_mat[j][i] = mg1_pi;
	  z3_mat[i][j] = mg1_f2;
	  z3_mat[j][i] = mg1_f2;
	}
	else if (iloci == 2) {
	  sscanf(str,"%lf %lf %lf %lf %lf %lf ",
		 &mg1_f1,&mg1_f2,&mg1_pi,
		 &mg2_f1,&mg2_f2,&mg2_pi);
	  if (mg1_pi == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_pi = 0.5;
	  }
	  if (mg1_f2 == non_value) {   
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_f2 = 0.5;
	  }
	  z1_mat[i][j] = mg1_pi;
	  z1_mat[j][i] = mg1_pi;
	  if (mg2_pi == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg2_pi = 0.5;
	  }
	  z2_mat[i][j] = mg2_pi;   
	  z2_mat[j][i] = mg2_pi;
	  z3_mat[i][j] = mg1_f2;
	  z3_mat[j][i] = mg1_f2;
	}
	else {
	  Rcerr << "*** Get_z_mat debug ***" << std::endl <<
	    "Error: can't handle more than two loci," << std::endl <<
	    "please change the parameter file." << std::endl <<
	    "(Maxfun.cpp line 343)" << std::endl;
	}
      }
    }           
  }
}


void Maxfun::ChangeMatrices(double **g_mat, double **z1_mat,double **z2_mat, double **z3_mat,
			    double **i_mat, double **c_mat, double **p_mat, double **q_mat,
			    int n_fam_member, int  *n_member) {
  int k;
  int prev_count;
  int n_mb;

  prev_count = 0;
  n_mb       = n_fam_member;

  //printf("n_fam=%d\n",n_mb);

  for (k = 0; k < n_fam_member; k++) {
    if (fam.person[k].missing_traitflag == 1) {
      //printf("k=%d missing\n",k);
      ShrinkMatrices(g_mat, z1_mat, z2_mat, z3_mat, i_mat, c_mat, p_mat, q_mat,
		     k-prev_count, n_mb);
      prev_count++;
      n_mb--;
    }
  }
  *n_member = n_mb;
}


void Maxfun::ShrinkMatrices(double **g_mat, double **z1_mat,double **z2_mat,double **z3_mat,
			    double **i_mat, double **c_mat, double **p_mat, double **q_mat,
			    int n_order,    int n_member) {
  int i,j,ii,jj;
  //    int numofcov;
  /*
printf("n_order=%d, n_member=%d\n", n_order, n_member);
*/

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
      z1_mat[i][j]= z1_mat[ii][jj];
      z2_mat[i][j]= z2_mat[ii][jj];
      z3_mat[i][j]= z3_mat[ii][jj];
      i_mat[i][j] = i_mat[ii][jj];
      c_mat[i][j] = c_mat[ii][jj];
      p_mat[i][j] = p_mat[ii][jj];
      q_mat[i][j] = q_mat[ii][jj];
    }
  }
}

/* This function is modified from AS_Calculate and NOAS_Calculate functions and 
   subtract redundant parts from the two functions.
*/
double Maxfun::Like_func() {
  int    i,j,k,kk,Ni;
  int    count_smt;
  double sing_flag;              /* singularity flag of the matrix.  */
  double **temp_mat1;
  double **V_mat, **V_mat_inv, **temp_mat;
  double ***G_mat, ***Z1_mat, ***Z2_mat, ***I_mat;
  double ***C_mat, ***P_mat, ***Q_mat;
  double **g_mat, **z1_mat, **z2_mat, **i_mat, **z3_mat;
  double **c_mat, **p_mat,  **q_mat;
  double **D_mat, **d_vec;
  double *ydb_vec,*beta_vec;
  double *i_vec,  **y_vec, *Y_vec, **I_vec;
  //    double **Exsder_mat, *fder_vec, *inc_vec;
  double ***Smtcpq_mat;
  double **MuBeta_vec;
  double V1, V2, V3, ln_func;
  char str[BUF];

  sing_flag      = GOOD;          /* set default to non_singularity.  */
	
  ln_func   =0.0;

  count_smt = iinitvcnum*3+env_vcnum+iinitvcnum*3;

  G_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z1_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z2_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  I_mat  = (double ***)malloc(env_vcnum  * sizeof(double **));
  C_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  P_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Q_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  // haven't considered about longitute data
  d_vec  = (double **) malloc(icovs      * sizeof(double *));
  y_vec  = (double **) malloc(total_trait_values    * sizeof(double *));
  Smtcpq_mat  = (double ***)malloc(count_smt          * sizeof(double **));
  MuBeta_vec  = (double **) malloc(((icovs+1)*total_trait_values)* sizeof(double *));

  fp_loci.clear();
  fp_loci.seekg(0, ios::beg);
  fp_loci.getline(str, BUF);
  dataIndex = 0;
  relationIndex = 0;

  /* get the traits and marker loci which need to be calculated 
     for each family. */
  for (i = 0; i < nfam; i++) {
    V_mat      = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    V_mat_inv  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    temp_mat   = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    temp_mat1  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    g_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    z1_mat     = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    z2_mat     = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    z3_mat     = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    i_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    c_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    p_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    q_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    D_mat      = Lib::dmatrix(0,N[i]-1,0,icovs);
    i_vec      = Lib::dvector(0,N[i]-1); 
    ydb_vec    = Lib::dvector(0,N[i]-1); 
    beta_vec   = Lib::dvector(0,icovs); 
    Y_vec      = Lib::dvector(0,total_trait_values*N[i]-1);

    for (j = 0; j < env_vcnum; j++) {
      I_mat[j]  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }

    for (j = 0; j < iinitvcnum; j++) {
      G_mat[j]  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Z1_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Z2_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      C_mat[j]  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      P_mat[j]  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Q_mat[j]  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }

    for (j = 0; j < total_trait_values; j++) {
      y_vec[j] = Lib::dvector(0,N[i]-1);
    }

    I_vec = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);

    //need change? 
    //for (j = 0; j <total_trait_values; j++) {
    //    I_vec[j] = Lib::dvector(0,total_trait_values*N[i]-1); 
    //}

    for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++)
      Smtcpq_mat[j]=Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);

    for (j = 0; j < (icovs+1)*total_trait_values; j++)
      MuBeta_vec[j]=Lib::dvector(0,total_trait_values*N[i]-1);

    // need change? for Sint data
    for (j = 0; j < icovs; j++) {
      d_vec[j]  = Lib::dvector(0,N[i]-1);
    }

    InitMat   (N[i], G_mat, Z1_mat, Z2_mat, I_mat, g_mat, z1_mat, z2_mat, i_mat);
    InitCPQMat(N[i], C_mat, P_mat,  Q_mat,  c_mat, p_mat, q_mat);
    InitVec   (N[i], i_vec, ydb_vec, MuBeta_vec);
    InitMissingTraitflag(N[i]);

    /* Get initial G_VAL 'fam' value here and it may be changed within 
       this function.                                             */                  
    GetData(N[i]);

    Build_D_Matrix(N[i], D_mat, d_vec);

    for (j = 0; j <total_trait_values; j++) {
      Build_beta_Vector(j, beta_vec);
      Lib::Multi_Matrix_Vector(ydb_vec, N[i], icovs+1, D_mat, beta_vec); 
      for(k = 0,kk = 0; k < N[i]; k++) {
	if (fam.person[k].missing_traitflag != 1) {
	  y_vec[j][kk] = fam.person[k].trait[j] - ydb_vec[k];
	  kk++;
	}
      }
    }

    fam.n_sib    = N[i] - 2;
    fam.sib_pair = fam.n_sib*(fam.n_sib-1)/2;

    /* changed 4/22/97;     E.Y.    */ 
    Get_z_mat(N[i],    z1_mat, z2_mat, z3_mat);
    Get_gcpq_mat(N[i], g_mat,  c_mat, p_mat, q_mat);

    /* After Change_gzi_Matrix_d_vec function, we got 'Ni' for new family member number.
       So we change N[i] to Ni from now until before free_dmatrix.  */ 
    ChangeMatrices(g_mat,z1_mat,z2_mat,z3_mat,i_mat,c_mat,p_mat,q_mat,N[i],&Ni);

    for (j = 0; j <total_trait_values; j++) {
      for(k = 0; k < Ni; k++) {
	Y_vec[j*Ni+k] = y_vec[j][k];
	I_vec[j][j*Ni+k] = i_vec[k];
      }
    }

    GetVGZICPQMatrix(V_mat,G_mat,Z1_mat,Z2_mat,I_mat,C_mat,P_mat,Q_mat,
		     g_mat,z1_mat,z2_mat,i_mat,c_mat,p_mat,q_mat,Ni); 

 
    Lib::CopyMatrix(temp_mat,total_trait_values*Ni, V_mat);
    GetSmtcpqMatrix(Smtcpq_mat, G_mat, Z1_mat, Z2_mat, I_mat, C_mat, P_mat, Q_mat,Ni);
    GetMuBetaVector(MuBeta_vec,  i_vec, d_vec, Ni);

    Lib::CopyMatrix(temp_mat1,total_trait_values*Ni, V_mat);
    sing_flag = Lib::DetermOfMatrix(total_trait_values*Ni, temp_mat1);
    if (sing_flag == BAD) {
      PrintErrMsg("V_mat","AS_CalculateValues",total_trait_values*Ni,V_mat);
    }

    Lib::InverseOfMatrix(V_mat_inv, total_trait_values*Ni, V_mat);
    /*
      PrintOneMatrix(itraits*Ni, temp_mat);
    */

    V1 = Lib::DetermOfMatrix(total_trait_values*Ni, temp_mat);
    V2 = log(V1);
    V3 = Lib::Multi_Vectors_Matrix(total_trait_values*Ni, Y_vec, V_mat_inv, Y_vec);
    /*
      #ifdef DEBUG 
    */
    //    printf("V1 = %lf\n", V1);
    //    printf("V2 = %lf\n", V2);
    //    printf("V3 = %lf\n", V3);
    Rcout << "ln_func = " << ln_func << std::endl;
    /*
      #endif
    */
    ln_func +=(-1.0/2.0)* (V2 + V3);

    Lib::free_dmatrix(V_mat,      0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(V_mat_inv , 0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(temp_mat,   0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(temp_mat1,  0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(g_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z1_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z2_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z3_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(i_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(c_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(p_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(q_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(D_mat, 0, N[i]-1, 0, icovs );
    Lib::free_dvector(i_vec,      0, N[i]-1); 
    Lib::free_dvector(ydb_vec,    0, N[i]-1); 
    Lib::free_dvector(beta_vec,   0, icovs ); 
    Lib::free_dvector(Y_vec,      0, total_trait_values*N[i]-1); 
    for (j = 0; j < total_trait_values; j++) {
      Lib::free_dvector(y_vec[j], 0,N[i]-1);
    }
    for (j = 0; j < env_vcnum; j++) {
      Lib::free_dmatrix(I_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }
    /*
      for (j = 0; j <total_trait_values; j++) {
      Lib::free_dvector(I_vec[j], 0,total_trait_values*N[i]-1); 
      }
    */
    Lib::free_dmatrix(I_vec,0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    for (j = 0; j < iinitvcnum; j++) {
      Lib::free_dmatrix(G_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Lib::free_dmatrix(Z1_mat[j],0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);	
      Lib::free_dmatrix(Z2_mat[j],0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Lib::free_dmatrix(C_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Lib::free_dmatrix(P_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Lib::free_dmatrix(Q_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }
    for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++) {
      Lib::free_dmatrix(Smtcpq_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }
    for (j = 0; j < (icovs+1)*total_trait_values; j++)
      Lib::free_dvector(MuBeta_vec[j], 0,total_trait_values*N[i]-1);
    for (j = 0; j < icovs; j++) {
      Lib::free_dvector(d_vec[j], 0,N[i]-1);
    }
  }

  free(G_mat);
  free(Z1_mat);
  free(Z2_mat);
  free(I_mat);
  free(C_mat);
  free(P_mat);
  free(Q_mat);
  free(d_vec);
  free(y_vec);
  free(Smtcpq_mat);
  free(MuBeta_vec);

#ifdef DEBUG 
  printf("Loglik  = %lf\n", ln_func);
#endif
  return ln_func;
}

/*
Once we got everything from the parameter file, it is the time to do 
the analysis and to get the new MU,SIGMA,MG and TAU.
This module calls following functions:
        Like_func;
        Correct_factor.
*/
double Maxfun::SRun()
{
  double ln_func, c_ln_func, logcf;
  //    double sum_mubeta, sum_vc;
  //    double temp_a[iinitvcnum], temp_c[iinitvcnum], temp_b[itraits];
   
#ifdef DEBUG
  Rprintf("This is the start values after Check&Fix Covariates\n");
  for (i = 0; i < itraits; i++) {
    Rprintf("%10f ", x_mu[i]);
    for (k = 0; k < icovs; k++)
      Rprintf("%10f ", x_beta[i][k]);
  }
  Rprintf("\n");
  //PrintResults(DEBUG);
#endif
       
  ismtcpq_dim = 0;
  GetSMT_vcnum();
  GetCPQ_vcnum();

  ismtcpq_dim = s_vcnum  + m1_vcnum + m2_vcnum + t_vcnum
    + c_vcnum  + p_vcnum  + q_vcnum;

  itraitcovnum = ismtcpq_dim*(ismtcpq_dim+1)/2;

  //printf("ascert_flag=%d\n",ascert_flag);

  /* 1 = YES .   */
  if (ascert_flag == 1) { 
    ln_func = Like_func();
    logcf = Correct_factor();
    c_ln_func = ln_func - logcf;
    return c_ln_func;
  }
  /* 2 = NO.   */
  else if (ascert_flag == 2) {
    ln_func = Like_func();
    return ln_func;
  } 
  return 0;
}

/* Print the error message and the matrix if it is singular matrix.
*/
void Maxfun::PrintErrMsg(const char *str_mat, const char *str_routine,
			 int mat_dim, double **x_mat)
{
  PROBLEM "This is the singular %s matrix at routine %s.\nMaxfun.cpp key 724\n",
    str_mat, str_routine RECOVER(NULL_ENTRY);
}


double Maxfun::Correct_factor()
{
  //    int     i, j;
  double  logcf; 
  double  S_val, Z_val, D_val, G_val;

  logcf = 0.0;

  S_val = sqrt(x_S[0]+x_M1[0]+x_M2[0]+x_T[0]+x_c[0]+x_p[0]+x_q[0]); 

  /* Using Elston method.              */
  if (method_flag == 1) {
    if (proband_flag == 1) {
      Z_val = (trait_array[0].Th_mu - x_mu[0])/S_val;
      D_val = dens(Z_val);
      G_val = (1. - cump(Z_val));

      //printf("S_val=%lf;   Z_val=%lf;   D_val=%lf\n", S_val, Z_val, D_val);
      //printf("G_val=%lf\n", G_val);
      //printf("\n %10.2f %10.2f %10.2f \n", Z_val, D_val, G_val); 

      logcf = nfam*log(G_val);
    }
    else {
      PROBLEM "Not available now.\nMaxFun.cpp key 753\n" RECOVER(NULL_ENTRY);
    }
  }
  /* Using Hopper method.              */
  if (method_flag == 2) {
    logcf = CF_1or2prob();
  }
  return logcf;
}

/* This function is modified from Like_func function, can calculate log likelihood for 
   1 or 2 probands. 
*/
double Maxfun::CF_1or2prob() {
  int    i,j,k,Ni;
  int    count_smt;
  double sing_flag;              /* singularity flag of the matrix.  */
  double **temp_mat1;
  double **V_mat, **V_mat_inv, **temp_mat;
  double ***G_mat, ***Z1_mat, ***Z2_mat, ***I_mat;
  double ***C_mat, ***P_mat, ***Q_mat;
  double **g_mat, **z1_mat, **z2_mat, **i_mat;
  double **c_mat, **p_mat,  **q_mat;
  double **D_mat, **d_vec;
  double *ydb_vec,*beta_vec;
  double *i_vec,  **y_vec, *Y_vec, **I_vec;
  double ***Smtcpq_mat; 
  double **MuBeta_vec;
  double V1, V2, V3, ln_func; 
  char str[BUF];

  sing_flag      = GOOD;          /* set default to non_singularity.  */

  ln_func   =0.0;
    
  Ni = proband_flag;              /* we use 2 probands.               */

  fp_loci.clear();
  fp_loci.seekg(0, ios::beg);
  fp_loci.getline(str, BUF);
  dataIndex = 0;
  relationIndex = 0;
    
  count_smt = iinitvcnum*3+env_vcnum+iinitvcnum*3;
  G_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z1_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z2_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  I_mat  = (double ***)malloc(env_vcnum  * sizeof(double **));
  C_mat  = (double ***)malloc(iinitvcnum * sizeof(double **)); 
  P_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Q_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  // haven't considered about longitute data
  d_vec  = (double **) malloc(icovs      * sizeof(double *));
  y_vec  = (double **) malloc(total_trait_values    * sizeof(double *));
  Smtcpq_mat  = (double ***)malloc(count_smt          * sizeof(double **));
  MuBeta_vec  = (double **) malloc(((icovs+1)*itraits)* sizeof(double *));

  V_mat      = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  V_mat_inv  = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  temp_mat   = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  temp_mat1  = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  g_mat      = Lib::dmatrix(0,Ni-1,0,Ni-1);
  z1_mat     = Lib::dmatrix(0,Ni-1,0,Ni-1);
  z2_mat     = Lib::dmatrix(0,Ni-1,0,Ni-1);
  i_mat      = Lib::dmatrix(0,Ni-1,0,Ni-1);
  c_mat      = Lib::dmatrix(0,Ni-1,0,Ni-1);
  p_mat      = Lib::dmatrix(0,Ni-1,0,Ni-1);
  q_mat      = Lib::dmatrix(0,Ni-1,0,Ni-1);
  D_mat      = Lib::dmatrix(0,Ni-1,0,icovs);
  i_vec      = Lib::dvector(0,Ni-1);
  ydb_vec    = Lib::dvector(0,Ni-1);
  beta_vec   = Lib::dvector(0,icovs);
  Y_vec      = Lib::dvector(0,total_trait_values*Ni-1);
  for (j = 0; j < env_vcnum; j++) {
    I_mat[j] = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  }
  for (j = 0; j < total_trait_values; j++) {   
    y_vec[j] = Lib::dvector(0,Ni-1);
  }

  I_vec = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  //need change?
  //for (j = 0; j < total_trait_values; j++) {
  //    I_vec[j] = Lib::dvector(0,total_trait_values*Ni-1);
  //}
  for (j = 0; j < iinitvcnum; j++) {
    G_mat[j] =  Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Z1_mat[j] = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Z2_mat[j] = Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    C_mat[j] =  Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    P_mat[j] =  Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Q_mat[j] =  Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  }

  for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++)
    Smtcpq_mat[j]=Lib::dmatrix(0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
        
  for (j = 0; j < (icovs+1)*total_trait_values; j++)
    MuBeta_vec[j]=Lib::dvector(0,total_trait_values*Ni-1);

  for (j = 0; j < icovs; j++) {
    d_vec[j]  = Lib::dvector(0,Ni-1);
  }

  InitMat   (proband_flag, G_mat, Z1_mat, Z2_mat, I_mat, g_mat, z1_mat, z2_mat, i_mat);
  InitCPQMat(proband_flag, C_mat, P_mat,  Q_mat,  c_mat, p_mat, q_mat);
  InitVec   (proband_flag, i_vec, ydb_vec, MuBeta_vec);
  //assume proband no missing, otherwise we are in trouble.
  //InitMissingTraitflag(proband_flag);

  /* get the traits and marker loci which need to be calculated 
     for each family. */
  for (i = 0; i < nfam; i++) {

    //cout<<"i="<<i<<std::endl;

    GetProb(N[i]);

    Build_D_Matrix(Ni, D_mat, d_vec);

    for (j = 0; j < total_trait_values; j++) {
      Build_beta_Vector(j, beta_vec);
      Lib::Multi_Matrix_Vector(ydb_vec,Ni,icovs+1,D_mat,beta_vec); 
      for(k = 0; k < Ni; k++) {
	y_vec[j][k] = fam.proband[k].trait_p - ydb_vec[k];
      }
    }

    // need to change?
    if (proband_flag == 2){
      GetProb_z_mat(N[i], z1_mat, z2_mat);
      GetProb_gcpq_mat(N[i],g_mat,c_mat,p_mat,q_mat);
    }

    /*
      printf("This is the z_mat\n");
      PrintOneMatrix(Ni, z_mat);
    */
    /* PrintMatrix(Ni, g_mat, z_mat, i_mat); */ 

    /*
      printf("Ni=%d\n", Ni);
    */
    /* PrintMatrix(Ni, g_mat, z_mat, i_mat);   */ 
    for (j = 0; j <total_trait_values; j++) {
      for(k = 0; k < Ni; k++) {
	Y_vec[j*Ni+k] = y_vec[j][k];
	I_vec[j][j*Ni+k] = i_vec[k];
      }
    }

    GetVGZICPQMatrix(V_mat,G_mat,Z1_mat,Z2_mat,I_mat,C_mat,P_mat,Q_mat,
		     g_mat,z1_mat,z2_mat,i_mat,c_mat,p_mat,q_mat,Ni);

    for(j=0;j<total_trait_values*Ni;j++){
      for(k=0;k<total_trait_values*Ni;k++)
	Rprintf("%.3f ",V_mat[j][k]);
      Rprintf("\n");
    }

    Lib::CopyMatrix(temp_mat,total_trait_values*Ni, V_mat);

    Lib::CopyMatrix(temp_mat1,total_trait_values*Ni, V_mat);
    sing_flag = Lib::DetermOfMatrix(total_trait_values*Ni, temp_mat1);
    if (sing_flag == BAD) {
      PrintErrMsg("V_mat","Like_func",total_trait_values*Ni,V_mat);
    }

    Lib::InverseOfMatrix(V_mat_inv, total_trait_values*Ni, V_mat);
    /*
      Lib::PrintOneMatrix(itraits*Ni, temp_mat);
    */

    V1 = Lib::DetermOfMatrix(total_trait_values*Ni, temp_mat);
    V2 = log(V1);
    V3 = Lib::Multi_Vectors_Matrix(total_trait_values*Ni, Y_vec, V_mat_inv, Y_vec);
    /*
      #ifdef DEBUG 
      printf("V1 = %lf\n", V1);
      printf("V2 = %lf\n", V2);
      printf("V3 = %lf\n", V3);
      printf("ln_func = %lf\n", ln_func);
      #endif
    */
    if (proband_flag == 1)
      ln_func += -(1.0/2.0)*log(2*3.1415926535897932) + (-1.0/2.0)* (V2 + V3);
    else if (proband_flag == 2)
      ln_func += -log(2*3.1415926535897932) + (-1.0/2.0)* (V2 + V3);
  }

  Lib::free_dmatrix(V_mat,      0, total_trait_values*Ni-1, 0, total_trait_values*Ni-1);
  Lib::free_dmatrix(V_mat_inv , 0, total_trait_values*Ni-1, 0, total_trait_values*Ni-1);
  Lib::free_dmatrix(temp_mat,   0, total_trait_values*Ni-1, 0, total_trait_values*Ni-1);
  Lib::free_dmatrix(temp_mat1,  0, total_trait_values*Ni-1, 0, total_trait_values*Ni-1);
  Lib::free_dmatrix(g_mat,      0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(z1_mat,     0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(z2_mat,     0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(i_mat,      0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(c_mat,      0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(p_mat,      0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(q_mat,      0, Ni-1, 0, Ni-1);
  Lib::free_dmatrix(D_mat,      0, Ni-1, 0, icovs );
  Lib::free_dvector(i_vec,      0, Ni-1); 
  Lib::free_dvector(ydb_vec,    0, Ni-1); 
  Lib::free_dvector(beta_vec,   0, icovs ); 
  Lib::free_dvector(Y_vec,      0, total_trait_values*Ni-1); 
  for (j = 0; j < total_trait_values; j++) {
    Lib::free_dvector(y_vec[j], 0,Ni-1);
  }
  for (j = 0; j < env_vcnum; j++) {
    Lib::free_dmatrix(I_mat[j], 0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  }
  Lib::free_dmatrix(I_vec,0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  for (j = 0; j < iinitvcnum; j++) {
    Lib::free_dmatrix(G_mat[j],  0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Lib::free_dmatrix(Z1_mat[j], 0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Lib::free_dmatrix(Z2_mat[j], 0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Lib::free_dmatrix(C_mat[j],  0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Lib::free_dmatrix(P_mat[j],  0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
    Lib::free_dmatrix(Q_mat[j],  0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  }
  for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++) {
    Lib::free_dmatrix(Smtcpq_mat[j], 0,total_trait_values*Ni-1,0,total_trait_values*Ni-1);
  }
  for (j = 0; j < (icovs+1)*total_trait_values; j++)
    Lib::free_dvector(MuBeta_vec[j], 0,total_trait_values*Ni-1);

  for (j = 0; j < icovs; j++) {
    Lib::free_dvector(d_vec[j], 0,Ni-1);
  }

  free(G_mat);
  free(Z1_mat);
  free(Z2_mat);
  free(I_mat);
  free(C_mat);
  free(P_mat);
  free(Q_mat);
  free(d_vec);
  free(y_vec);
  free(Smtcpq_mat);
  free(MuBeta_vec);

#ifdef DEBUG 
  Rprintf("\nLoglik = %lf\n", ln_func);
#endif
  return ln_func;
}


/*
Get the proband(s) for each family from family data 
file 'fort.12'.
*/
void Maxfun::GetProb(int n_fam_member)
{
  int  i,j,k,r,s;
  int  n_len, prob_sign, prob_num;
  //    int  digit_flag, space_flag;    /* for searching phenotype alleles.    */
  //    int  right_flag;                /* checking the phenotype in datafile. */
  int  t_end[MAXNUMTRAIT];
  char str[BUF];
  char trait_str[MAXNUMTRAIT][BUF];
  char prob_str[STRING];

  for ( i = 0; i < total_trait_values; i++) {
    t_end[i] = t_start[i] + t_len[i] -1;
  }
  prob_num = 0;
  for ( j = 0; j < n_fam_member; j++) {
    /*  get the trait, marker and covariate from each record.  */
    //    fp_data.getline(str,BUF);
    n_len = strlen(str);

    s = 0;
    r = 0;

    //Proband located in last 2 columns (default)

    for ( k = (n_len-2); k < n_len; k++) {
      if (!isspace(str[k]) ) {
	prob_str[r]   = str[k]; 
	prob_str[r+1] = '\0';
	r++;
      }		
    }
    prob_sign = atoi(prob_str);
    if ((prob_sign == 1) || (prob_sign == 2)) {
      for (i = 0; i < total_trait_values; i++) {
	r = 0;
	for (k = t_start[i]-1; k < t_end[i]; k++) {
	  if (!isspace(str[k]) ) {
	    trait_str[i][r] = str[k]; 
	    trait_str[i][r+1] = '\0';
	    r++;
	  }
	}		
	fam.proband[prob_num].trait_p     = atof(trait_str[i]);
	fam.proband[prob_num].prob_mem = j; 
	prob_num++;
	/*
	  printf("proband string is:\n%s\n", str);
	*/
      }
    }
  }
}


/* Get z1_mat, z2_mat matrices from file fp_loci.
   changed 9/17/1997.
*/
void Maxfun::GetProb_z_mat(int n_fam_member, double **z1_mat, double **z2_mat) {
  int    i,j;
  double mg1_f1, mg1_f2, mg1_pi;
  double mg2_f1, mg2_f2, mg2_pi;
  char   str[BUF];
  int    zi, zj;

  zi = 1;
  zj = 0;  

  if(iloci==0) return;

  for ( i = 0; i < n_fam_member; i++) {
    for ( j = i; j < n_fam_member; j++) {
      if (j > i) {
	fp_loci.getline(str, BUF);
	if (iloci == 1) {
	  sscanf(str,"%lf %lf %lf ", &mg1_f1,&mg1_f2,&mg1_pi);
	  if (mg1_pi == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_pi = 0.5;
	  }
	  if ((i == fam.proband[0].prob_mem) &&
	      (j == fam.proband[1].prob_mem) ){
	    z1_mat[zi][zj] = mg1_pi;
	    z1_mat[zj][zi] = mg1_pi;
	  }
	}
	else if (iloci == 2) {
	  sscanf(str,"%lf %lf %lf %lf %lf %lf ",
		 &mg1_f1,&mg1_f2,&mg1_pi,
		 &mg2_f1,&mg2_f2,&mg2_pi);
	  if (mg1_pi == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_pi = 0.5;
	  }

	  if (mg2_pi == non_value) {
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg2_pi = 0.5;
	  }
	  if ((i == fam.proband[0].prob_mem) &&
	      (j == fam.proband[1].prob_mem) ){
	    z1_mat[zi][zj] = mg1_pi;
	    z1_mat[zj][zi] = mg1_pi;
	    z2_mat[zi][zj] = mg2_pi;
	    z2_mat[zj][zi] = mg2_pi;
	  }
	}
	else {
	  Rcerr << "Error: can't handle more than two loci," << std::endl;
	  Rcerr << "please change the parameter file." << std::endl << std::endl;
	}
      }
    }
  }
}

/* Get g_mat, c_mat, p_mat, q_mat matrices from file fp_share.
   changed 9/11/97.
*/
void Maxfun::GetProb_gcpq_mat(int n_fam_member,  double **g_mat, double **c_mat,
			      double **p_mat,    double **q_mat) {
  int    i,j;
  int    zi, zj;

  zi = 1;
  zj = 0;

  for ( i = 0; i < n_fam_member; i++) {
    for ( j = i; j < n_fam_member; j++) {
      if (j > i) {
	// Variables shareOutTemp1 and shareOutTemp2 were added here because
	// the format of share.out has changed.  These variables are unused in
	// the program.
	if ((i == fam.proband[0].prob_mem) &&
	    (j == fam.proband[1].prob_mem) ){
	  g_mat[i][j] = g_mat[j][i] = shareArray[relationIndex].geneticSimilarity;
      
	  c_mat[i][j] = c_mat[j][i] = shareArray[relationIndex].areSiblings;
	  
	  p_mat[i][j] = p_mat[j][i] = shareArray[relationIndex].areSpouses;
	  
	  q_mat[i][j] = q_mat[j][i] = shareArray[relationIndex].areParentChild;
	  
	  relationIndex++;
	}
      }
    }
  }
}


/*
 Error function (related simply to probability integral).
 It calculates the area under the standard normal curve from -x to x
 */
double Maxfun::erf(double  x)
{
  double  t, z;

  if (x < 0.) return -erf(-x);

  t = 1./(1.+x*.3275911);
  z = ((((t*1.061405429 - 1.453152027)*t+1.421413741)*t-.284496736)*t
       +.254829592)*t;

  return (1.-exp(-x*x)*z);
}

/*
 It determines the f(x) where f(.) is the standard normal density
 */
double Maxfun::dens(double  x)
{
  double  t;
     
  t = exp(-(x*x)/2.)/sqrt(2.*3.141592654);

  return(t);
} 

/*
 Cummulative probability integral.
 It determines the area under normal curve from -infinty to x
 */
double Maxfun::cump(double  x)
{
  static double rt = .7071067812;

  return (.5*(1.+erf(x*rt)));
}


/* Get (itraits*n_fam_member) MuBeta_vec vector by using  i_vec, d_vec vectors.
*/
void Maxfun::GetMuBetaVector(double *MuBeta_vec[],double *i_vec,double *d_vec[],int n_member) {
  int i,j,k,kv;

  kv = 0;
     
  for (i = 0; i < itraits; i++) {
    for (j = 0; j < (icovs+1); j++) {
      if (j == 0) {
	for(k = 0; k < n_member; k++) {
	  MuBeta_vec[kv][i*n_member+k] = i_vec[k]; 
	}
      }
      else {
	for(k = 0; k <n_member ; k++) {
	  MuBeta_vec[kv][i*n_member+k] = d_vec[j-1][k];
	}
      }
      kv++;
    }
  }
}


/*
Build D matrix and the covariate columns d_vec of D matrix.
*/
void Maxfun::Build_D_Matrix(int n_fam_member, double **D_mat, double *d_vec[]) {
  int i,j;
            
  /* Build D matrix.         */
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < (icovs+1); j++) {
      if (j == 0)
	D_mat[i][j] = 1;
      else
	D_mat[i][j] = fam.person[i].cov[j-1];
    }
  }

  /* Build the covariate columns d_vec of D matrix.    */
  for (j = 0; j < icovs; j++) {
    for (i = 0; i < n_fam_member; i++) {
      d_vec[j][i] = fam.person[i].cov[j];
    }
  }
}


double Maxfun::getmajorgene1(){
  return x_M1[0];
}
        
double Maxfun::getmajorgene2(){
  return x_M2[0];
}
        
double Maxfun::getpolygene(){
  return x_S[0];
}
        
double Maxfun::getenvironment(){
  return x_T[0];
}
        
double Maxfun::getC(){
  return x_c[0];
}

double Maxfun::getP(){
  return x_p[0];
}
        
double Maxfun::getQ(){
  return x_q[0];
}
