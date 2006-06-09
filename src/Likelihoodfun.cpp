/******************************************************************************
File: Likelihoodfun.cpp
Description:  Methods of a Likelihoodfun object performs many matrix operations
              such as shrinking and initalization from values located in files.
Author: Eric Lunde, 7-28-03 (This file was not writen by Eric, but this is the
        date he took over and began modifying it.)
Updates: (Date, Modified By, Modification Description)
7-28-03, Eric Lunde, About 2-3 months ago, the file multic.par was removed from
         the multic process.  I removed the ifstream fp_para object that was
         meant to reference that file.

         The method GetInitValues was used to read values from multic.par.
         Since that file is no longer used, I removed this method.  It is still
         located in version 6 and before.

         One of the problems I encountered in this process was that all of a
         sudden, multic would not compete execution.  This was because I
         removed code in Calculate.cpp (a subclass of Likelihood.cpp) that
         opened, closed, and reset the get pointer to a file named share.out
         because no where else in that code was share.out being read.  The
         actual reading was taking place in Likelihood.cpp.  I moved the
         opening and closing of fp_share, fp_data, and fp_lookuplog into this
         file to promote encapsulation.  I did not move the opening and closing
         of fp_loci because, it is opened inside of a conditional statement.  I
         thought that it belonged in Calculate.cpp.
7-29-03, Eric Lunde, I modified Get_gcpq_mat to read from the ShareRelation
         array instead of the file share.out.  I then commeneted the previous
         (read from share.out) code to test the timing of the program.  I then
8-6-03, Eric Lunde, In the method GetData, I have removed all of the code that
        read from fort.12.  Instead, it now reads form the FortData array
        fortArray.
8-22-03, Eric Lunde, When reading from loci.out, we only want to discard the
         first narrative line.  The method Get_z_mat was discarding the next
         line (first line if it is the first time called) every time it is
         called.  This incorrect algorithm has been corrected.
******************************************************************************/
#include "Likelihoodfun.h"
#include "Calculate.h"
#include "Lib.h"
#include "multic.h"
#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <S.h>

using namespace std;

extern int ascert_flag;

Likelihoodfun::Likelihoodfun() {
  fp_lookuplog.open(lookuplog_file, ios::app);
  if(fp_lookuplog == NULL) {
    PROBLEM "The file %s could not be opened for reading.\nLikelihoodfun.cpp key 44\n",
      lookuplog_file RECOVER(NULL_ENTRY);
  }
}

Likelihoodfun::~Likelihoodfun() {
  fp_lookuplog.fill('#');
  fp_lookuplog.width(70);
  fp_lookuplog << endl << "#" << endl;
  fp_lookuplog.close();
}


/* Change g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat if there is 
   any missing trait value in the family members.
   It calls the function:
   ShrinkMatrices.
*/
void Likelihoodfun::ChangeMatrices(double **g_mat, double **z1_mat,
				   double **z2_mat,double **i_mat,
				   double **c_mat, double **p_mat,
				   double **q_mat, int n_fam_member,
				   int  *n_member) {
  int k;
  int prev_count;
  int n_mb;

  prev_count = 0;
  n_mb       = n_fam_member;

  //printf("n_fam=%d\n",n_mb);
  for (k = 0; k < n_fam_member; k++) {
    if (fam.person[k].missing_traitflag == 1) {
      //printf("k=%d missing\n",k);
      ShrinkMatrices(g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat, 
		     k-prev_count, n_mb);
      prev_count++;
      n_mb--;
    }
  }
  *n_member = n_mb;
}

/* Shrink g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat matrices 
   at one time.
*/
void Likelihoodfun::ShrinkMatrices(double **g_mat, double **z1_mat,double **z2_mat, double **i_mat, 
				   double **c_mat, double **p_mat, double **q_mat, 
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
      }else if (i == n_order) {
	ii = i+1;
      }else if (j == n_order) {
	jj = j+1;
      }
      g_mat[i][j] = g_mat[ii][jj];
      z1_mat[i][j]= z1_mat[ii][jj];
      z2_mat[i][j]= z2_mat[ii][jj];
      i_mat[i][j] = i_mat[ii][jj];
      c_mat[i][j] = c_mat[ii][jj];
      p_mat[i][j] = p_mat[ii][jj];
      q_mat[i][j] = q_mat[ii][jj];
    }
  }
}

/* Get (n_fam_member*total_trait_values) dimension V_mat by using formula: 
       V_mat = (x_S  * g_mat)  + (x_M1 * z1_mat) + (x_M2 * z2_mat) + (x_T  * i_mat)
             + (x_c * c_mat)   + (x_p * p_mat)   + (x_q * q_mat)
   and G_mat, Z1_mat, Z2_mat, I_mat, C_mat, P_mat, Q_mat by using 
   g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat.
*/
void Likelihoodfun::GetVGZICPQMatrix(double **V_mat, double **G_mat[],
				     double **Z1_mat[], double **Z2_mat[],
				     double **I_mat[], double **C_mat[],
				     double **P_mat[],  double **Q_mat[],
				     double **g_mat, double **z1_mat,
				     double **z2_mat, double **i_mat,
				     double **c_mat, double **p_mat,
				     double **q_mat, int n_fam_member) {
  int i,j,k,kk,kv;
  double ***vm;

  /* Initialize the size of the pointers.    5-12-98; E.Y.  */
  vm = (double ***)malloc(iinitvcnum * sizeof(double **));

  for (k = 0; k < iinitvcnum; k++) {
    vm[k] = Lib::dmatrix(0, n_fam_member-1, 0, n_fam_member-1);
  }

  for (k = 0; k < iinitvcnum; k++) {
    for (i = 0; i < n_fam_member; i++) {
      for (j = 0; j < n_fam_member; j++) {
	vm[k][i][j] = 0.0;
      }
    }
  }

  /*
  cout << endl << "n_fam_member = " << n_fam_member << endl;
  printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",x_S[0],x_M1[0],x_M2[0],x_T[0],x_c[0],x_p[0],x_q[0]);
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ", g_mat[i][j]);
    printf("\n");
  }
  printf("\n");
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ",z1_mat[i][j]);
    printf("\n");
  }   
  printf("\n");
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ",z2_mat[i][j]);
    printf("\n");
  }   
  printf("\n");
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ",i_mat[i][j]);
    printf("\n");
  }   
  printf("\n");
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ",c_mat[i][j]);
    printf("\n");
  }   
  printf("\n");
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ",p_mat[i][j]);
    printf("\n");
  }   
  printf("\n");
  for (i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++)  
      printf("%.3f ",q_mat[i][j]);
    printf("\n");
  }   
  printf("\n");
  */
  
  /*  vm = x_S *g_mat + x_M1*z1_mat + x_M2*z2_mat + x_c*c_mat + x_p*p_mat + x_q*q_mat. */ 
  for (k = 0; k < iinitvcnum; k++) {
    /*
    cout<<x_S[k]<<endl;
    cout<<x_M1[k]<<endl;
    cout<<x_M2[k]<<endl;
    cout<<x_c[k]<<endl;
    cout<<x_p[k]<<endl;
    cout<<x_q[k]<<endl;
    */
    
    for (i = 0; i < n_fam_member; i++) {
      for (j = 0; j < n_fam_member; j++) {
	vm[k][i][j] += x_S[k]  * g_mat[i][j]
	  + x_M1[k] * z1_mat[i][j]
	  + x_M2[k] * z2_mat[i][j]
	  + x_c[k]  * c_mat[i][j]
	  + x_p[k] * p_mat[i][j]
	  + x_q[k]  * q_mat[i][j];
      }
    }
  }

  /* vm = vm(previous) + x_T*i_mat  .    */
/*
  if (datatype == MULTIVARIATE) {
    k_prev = 0;
    for (k = 0, kk = total_trait_valules; k < total_trait_values; k++, kk--) {
      for (i = 0; i < n_fam_member; i++) {
	for (j = 0; j < n_fam_member; j++) {
	  vm[k_prev][i][j] +=  x_T[k] * i_mat[i][j]; 
	}
      }
      k_prev = k_prev + kk; 
    }
  }else if (datatype == LONGITUDINAL) {
*/
    for (k = 0; k < env_vcnum; k++) {
      for (i = 0; i < n_fam_member; i++) {
	for (j = 0; j < n_fam_member; j++) {
	  vm[k][i][j] += x_T[k] * i_mat[i][j]; 
	}
      }
    }
    //  }

  /* build upper corner of V_mat by using (total_trait_values X total_trait_values) parts of vm and   
       build upper corner of G_mat,Z1_mat,Z2_mat,I_mat, C_mat, P_mat, Q_mat by using 
       g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat.    */
  kv = 0; 
  for (k = 0; k < total_trait_values; k++) {
    for (kk=k; kk < total_trait_values; kk++) {
      for (i = 0; i < n_fam_member; i++) {
	for (j = 0; j < n_fam_member; j++) {

	  V_mat[k*n_fam_member+i][kk*n_fam_member+j]      = vm[kv][i][j];
	  V_mat[kk*n_fam_member+i][k*n_fam_member+j]      = vm[kv][i][j];

	  G_mat[kv][k*n_fam_member+i][kk*n_fam_member+j]  = g_mat[i][j];
	  G_mat[kv][kk*n_fam_member+i][k*n_fam_member+j]  = g_mat[i][j];

	  Z1_mat[kv][k*n_fam_member+i][kk*n_fam_member+j] = z1_mat[i][j];
	  Z1_mat[kv][kk*n_fam_member+i][k*n_fam_member+j] = z1_mat[i][j];

	  Z2_mat[kv][k*n_fam_member+i][kk*n_fam_member+j] = z2_mat[i][j];
	  Z2_mat[kv][kk*n_fam_member+i][k*n_fam_member+j] = z2_mat[i][j];
/*
	  if (datatype == MULTIVARIATE) {
	    if ( kk == k)
	      I_mat[k][k*n_fam_member+i][k*n_fam_member+j]= i_mat[i][j];
	  }else if (datatype == LONGITUDINAL) {
*/
	    I_mat[kv][k*n_fam_member+i][kk*n_fam_member+j]  = i_mat[i][j];
	    I_mat[kv][kk*n_fam_member+i][k*n_fam_member+j]  = i_mat[i][j];
	    //	  }

	  C_mat[kv][k*n_fam_member+i][kk*n_fam_member+j]  = c_mat[i][j];
	  C_mat[kv][kk*n_fam_member+i][k*n_fam_member+j]  = c_mat[i][j];

	  P_mat[kv][k*n_fam_member+i][kk*n_fam_member+j]  = p_mat[i][j];
	  P_mat[kv][kk*n_fam_member+i][k*n_fam_member+j]  = p_mat[i][j];

	  Q_mat[kv][k*n_fam_member+i][kk*n_fam_member+j]  = q_mat[i][j];
	  Q_mat[kv][kk*n_fam_member+i][k*n_fam_member+j]  = q_mat[i][j];
	}
      }
      kv++;
    }
  }

  /*
        for(j=0;j<total_trait_values*n_fam_member;j++){
                for(k=0;k<total_trait_values*n_fam_member;k++)
                        printf("%.3f ",V_mat[j][k]);
                printf("\n");
                }
	*/
  for (k = 0; k < iinitvcnum; k++)
    Lib::free_dmatrix(vm[k], 0, n_fam_member-1, 0, n_fam_member-1);
  free(vm);
}

/* Get (total_trait_values*n_fam_member) dimension Smtcpq_mat by using  G_mat, Z1_mat, Z2_mat 
   and I_mat. 
*/
void Likelihoodfun::GetSmtcpqMatrix(double **Smtcpq_mat[], double **G_mat[],
				    double **Z1_mat[], double **Z2_mat[], 
				    double **I_mat[], double **C_mat[],
				    double **P_mat[],  double **Q_mat[],
				    int n_fam_member) {
  int i,j,ii,jj,k,kv;
  //    int loop_num;   

  kv = 0;

  for (k = 0; k < iinitvcnum; k++) {
    if (s_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = G_mat[k][i][j];
	}
      }
      kv++;
    }
  }


  for (k = 0; k < iinitvcnum; k++) {
    if (m1_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = Z1_mat[k][i][j];
	}
      }
      kv++;
    }
  }

  for (k = 0; k < iinitvcnum; k++) {
    if (m2_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = Z2_mat[k][i][j];
	}
      }
      kv++;
    }
  }


  /* changed 9/17/96.     */
  for (k = 0; k < env_vcnum; k++) {
    if (t_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = I_mat[k][i][j];
	}
      }
      kv++;
    }
  }

  /* added 7/8/97.       */
  for (k = 0; k < iinitvcnum; k++) {
    if (c_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = C_mat[k][i][j];
	}
      }
      kv++;
    }
  }

  for (k = 0; k < iinitvcnum; k++) {
    if (p_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = P_mat[k][i][j];
	}
      }
      kv++;
    }
  }
  for (k = 0; k < iinitvcnum; k++) {
    if (q_flag_array[k] == ESTIMATE) {
      for (i = 0; i < total_trait_values*n_fam_member; i++) {
	for (j = 0; j < total_trait_values*n_fam_member; j++) {
	  Smtcpq_mat[kv][i][j] = Q_mat[k][i][j];
	}
      }
      kv++;
    }
  }

  for (kv = 0; kv < ismtcpq_dim; kv++) {
    for (i = 0; i < total_trait_values*n_fam_member; i++) {
      for (j = 0; j < total_trait_values*n_fam_member; j++) {
	if (j > i) {
	  ii = i;
	  jj = j;
	  Smtcpq_mat[kv][jj][ii] = Smtcpq_mat[kv][i][j];
	}
      }
    }
  }
}

                  
 
/* Get the initial values of init_x_mu, init_x_S, init_x_M1, init_x_M2, init_x_T,
   init_x_c, init_x_p, init_x_q   from the strings.  */
void Likelihoodfun::GetXValues(int numval, double x_val[], char x_str[]) {
  int i,j,kv;
  int len;
  char str[STRING], temp_str[STRING];

  j  = 0;
  kv = 0;
  len = strlen(x_str);

  for (i = 0; i < len; i++) {
    /* for the first itraits-1  value(s) of x_mu.  */
    if (kv < numval-1) { 
      if (x_str[i] != ',') { 
	str[j]   = x_str[i];
	str[j+1] = '\0';
	j++;
      }
      else { 
	strcpy(temp_str, str);
	x_val[kv] = atof(temp_str);
	j = 0;
	kv++;
      }
    }
    /* for the last value of x_mu.  */
    else {
      str[j] = x_str[i];
      str[j+1] = '\0';
      j++;
    }
  }
  strcpy(temp_str, str);
  x_val[kv] = atof(temp_str);
}


/* Set Initial missing_traitflag = 0 to avoid the previous family's effect.
*/
void Likelihoodfun::InitMissingTraitflag(int n_fam_member) {
  int k;

  for(k = 0; k < n_fam_member; k++)
    fam.person[k].missing_traitflag = 0;
}

/* Initiate the values in G_mat, Z1_mat, Z2_mat, I_mat, g_mat, 
   z1_mat, z2_mat and i_mat matrices. 
*/
void Likelihoodfun::InitMat(int n_fam_member, double **G_mat[],
			    double **Z1_mat[], double **Z2_mat[], 
			    double **I_mat[], double **g_mat,
			    double **z1_mat,  double **z2_mat, 
			    double **i_mat) {
  int i,j,k;

  /* initiate G_mat, Z1_mat, Z2_mat.      */         
  for (k = 0; k < iinitvcnum; k++) {
    for (i = 0; i < total_trait_values*n_fam_member; i++) {
      for (j = 0; j < total_trait_values*n_fam_member; j++) {
	G_mat[k][i][j]  = 0.0;
	Z1_mat[k][i][j] = 0.0;
	Z2_mat[k][i][j] = 0.0;
      } 
    }
  }

  /* initiate I_mat.             */         
  for (k = 0; k < env_vcnum; k++) {
    for (i = 0; i < total_trait_values*n_fam_member; i++) {
      for (j = 0; j < total_trait_values*n_fam_member; j++) {
	I_mat[k][i][j] = 0.0;
      } 
    }
  }

  for(i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++) {
      if (j == i) {
	g_mat[i][j]  = 1.0;
	z1_mat[i][j] = 1.0;
	z2_mat[i][j] = 1.0;
	i_mat[i][j]  = 1.0;
      }else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
	g_mat[i][j]  = 0.0;
	z1_mat[i][j] = 0.0;
	z2_mat[i][j] = 0.0;
	i_mat[i][j]  = 0.0;
      }else {
	g_mat[i][j]  = 0.5;
	z1_mat[i][j] = 0.5;
	z2_mat[i][j] = 0.5;
	i_mat[i][j]  = 0.0;
      }
    }
  }
}


/* Initiate the values in  i_vec, ydb_vec, Mubeta_vec vectors.
*/
void Likelihoodfun::InitVec(int n_fam_member, double *i_vec,double *ydb_vec, 
			    double *MuBeta_vec[]) {
  int i,k;

  /* initiate MuBeta_vec.      */         
  for (k = 0; k < (icovs+1)*itraits; k++) {
    for (i = 0; i < total_trait_values*n_fam_member; i++) {
      MuBeta_vec[k][i] = 0.0;
    }
  }

  for(i = 0; i < n_fam_member; i++) {
    i_vec[i]   = 1.0;
    ydb_vec[i] = 0.0;
  }
}


/* Initiate the values in C_mat, P_mat, Q_mat, c_mat, p_mat and q_mat matrices. 
*/
void Likelihoodfun::InitCPQMat(int n_fam_member, double **C_mat[],
			       double **P_mat[], double **Q_mat[], 
			       double **c_mat, double **p_mat,
			       double **q_mat) {
  int i,j,k;

  /* initiate C_mat, P_mat, Q_mat.      */         
  for (k = 0; k < iinitvcnum; k++) {
    for (i = 0; i < total_trait_values*n_fam_member; i++) {
      for (j = 0; j < total_trait_values*n_fam_member; j++) {
	C_mat[k][i][j] = 0.0;
	P_mat[k][i][j] = 0.0;
	Q_mat[k][i][j] = 0.0;
      } 
    }
  }

  for(i = 0; i < n_fam_member; i++) {
    for (j = 0; j < n_fam_member; j++) {
      if (j == i) {
	c_mat[i][j] = 1.0;
	p_mat[i][j] = 1.0;
	q_mat[i][j] = 1.0;
      }else if ((i == 0 && j == 1) || (i == 1 && j == 0)) {
	c_mat[i][j] = 0.0;
	p_mat[i][j] = 1.0;
	q_mat[i][j] = 0.0;
      }else {
	c_mat[i][j] = 0.0;
	p_mat[i][j] = 0.0;
	q_mat[i][j] = 0.0;
      }
    }
  }
}


/* Get start values of G_VAL  
   'mu_flag',   's_flag', 'm1_flag', 'm2_flag', 't_flag', 
   'x_mu',      'x_S',    'x_M1',    'x_M2',    'x_T' and
   'x_beta' in this function. */
void Likelihoodfun::GetStartMubeta_SMT() {
  int i,j;

  mu_flag = init_mu_flag;
  s_flag  = init_s_flag ;
  m1_flag = init_m1_flag;
  m2_flag = init_m2_flag;
  t_flag  = init_t_flag;

  for (i = 0; i < iinitvcnum; i++) {
    x_S[i]       =  init_x_S[i] ;
    x_M1[i]      =  init_x_M1[i];
    x_M2[i]      =  init_x_M2[i];
  }
  for (i = 0; i < env_vcnum; i++) {
    x_T[i]       =  init_x_T[i];
  }
  for (i = 0; i < itraits; i++) {
    x_mu[i]      =  init_x_mu[i];
    for (j = 0; j < icovs; j++) 
      x_beta[i][j] =  init_x_beta[i][j]; 
  }

  GetST_flag_array();
  GetM_flag_array();
}

/* Get start values of G_VAL  
   'c_flag',   'p_flag','q_flag', 
   'x_c',      'x_p',   'x_q',
   in this function.                */
void Likelihoodfun::GetStartCPQ() {
  int i;

  c_flag =  init_c_flag;
  p_flag =  init_p_flag;
  q_flag =  init_q_flag;

  for (i = 0; i < iinitvcnum; i++) {
    x_c[i]      =  init_x_c[i];
    x_p[i]      =  init_x_p[i];
    x_q[i]      =  init_x_q[i];
  }

  GetCPQ_flag_array();
}


/* Get G_VAL 's_flag_array', 't_flag_array' values 
   here and they may be changed within the iteration in SRun.   
*/
void Likelihoodfun::GetST_flag_array() {
  int i,j,k;
    
  k = total_trait_values;

  /* Assign values to s_flag_array.   */
  if (s_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      s_flag_array[i] = ESTIMATE;
      s_fix_flag_array[i] = non_value;
    } 
  }
  else if (s_flag == CONSTRAINT) {
    for (i = 0, j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	s_flag_array[i] = ESTIMATE;
	s_fix_flag_array[i] = non_value;
	j = i;
	if (i != 0)
	  k--;
      }
      else { 
	s_flag_array[i] = CONSTRAINT;
	s_fix_flag_array[i] = non_value;
      }
    }
  }
  else if (s_flag == FIXED) {
    for (i = 0; i < iinitvcnum; i++) {
      s_flag_array[i]     = FIXED;
      s_fix_flag_array[i] = EXTERNALLY;
    }
  }

  k = total_trait_values;
  /* Assign values to t_flag_array.   */
  if (t_flag == ESTIMATE) {
    for (i = 0; i < env_vcnum; i++) {
      t_flag_array[i] = ESTIMATE;
      t_fix_flag_array[i] = non_value;
    }
  }else if (t_flag == CONSTRAINT) {
    //    if (datatype == LONGITUDINAL) {
      for (i = 0, j = 0; i < env_vcnum; i++) {
	if ((i == 0) || (i == j+k)) {
	  t_flag_array[i] = ESTIMATE;
	  t_fix_flag_array[i] = non_value;
	  j = i;
	  if (i != 0)
	    k--;
	}else { 
	  t_flag_array[i] = CONSTRAINT;
	  t_fix_flag_array[i] = non_value;
	}
      }
      //    }
  }else if (t_flag == FIXED) {
    for (i = 0; i < env_vcnum; i++) {
      t_flag_array[i]     = FIXED;
      t_fix_flag_array[i] = EXTERNALLY;
    }
  }
}

/* Get G_VAL 'm1_flag_array' , 'm2_flag_array' values 
   here and they may be changed within the iteration in SRun.   
*/
void Likelihoodfun::GetM_flag_array() {
  int i,j,k;
    
  k = total_trait_values;

  /* Assign values to m1_flag_array.   */
  if (m1_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      m1_flag_array[i] = ESTIMATE;
      m1_fix_flag_array[i] = non_value;
    }
  }else if (m1_flag == CONSTRAINT) {
    for (i = 0, j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	m1_flag_array[i] = ESTIMATE;
	m1_fix_flag_array[i] = non_value;
	j = i;
	if (i != 0)
	  k--;
      }else { 
	m1_flag_array[i] = CONSTRAINT;
	m1_fix_flag_array[i] = non_value;
      }
    }
  }else if (m1_flag == FIXED) {
    for (i = 0; i < iinitvcnum; i++) {
      m1_flag_array[i]     = FIXED;
      m1_fix_flag_array[i] = EXTERNALLY;
    }
  }

  /* Assign values to m2_flag_array.   */
  if (m2_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      m2_flag_array[i] = ESTIMATE;
      m2_fix_flag_array[i] = non_value;
    }
  }
  else if (m2_flag == CONSTRAINT) {
    for (i = 0, j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	m2_flag_array[i] = ESTIMATE;
	m2_fix_flag_array[i] = non_value;
	j = i;
	if (i != 0)
	  k--;
      }
      else { 
	m2_flag_array[i] = CONSTRAINT;
	m2_fix_flag_array[i] = non_value;
      }
    }
  }
  else if (m2_flag == FIXED) {
    for (i = 0; i < iinitvcnum; i++) {
      m2_flag_array[i]     = FIXED;
      m2_fix_flag_array[i] = EXTERNALLY;
    }
  }
}

/* Get G_VAL 'c_flag_array', 'p_flag_array', 'q_flag_array' values 
   here and they may be changed within the iteration in SRun.   
*/
void Likelihoodfun::GetCPQ_flag_array() {
  int i,j,k;
    
  k = total_trait_values;

  /* Assign values to c_flag_array.   */
  if (c_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      c_flag_array[i] = ESTIMATE;
      c_fix_flag_array[i] = non_value;
    } 
  }
  else if (c_flag == CONSTRAINT) {
    for (i = 0, j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	c_flag_array[i] = ESTIMATE;
	c_fix_flag_array[i] = non_value;
	j = i;
	if (i != 0)
	  k--;
      }
      else { 
	c_flag_array[i] = CONSTRAINT;
	c_fix_flag_array[i] = non_value;
      }
    }
  }
  else if (c_flag == FIXED) {
    for (i = 0; i < iinitvcnum; i++) {
      c_flag_array[i]     = FIXED;
      c_fix_flag_array[i] = EXTERNALLY;
    }
  }


  /* Assign values to p_flag_array.   */
  if (p_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      p_flag_array[i] = ESTIMATE;
      p_fix_flag_array[i] = non_value;
    } 
  }
  else if (p_flag == CONSTRAINT) {
    for (i = 0, j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	p_flag_array[i] = ESTIMATE;
	p_fix_flag_array[i] = non_value;
	j = i;
	if (i != 0)
	  k--;
      }
      else { 
	p_flag_array[i] = CONSTRAINT;
	p_fix_flag_array[i] = non_value;
      }
    }
  }
  else if (p_flag == FIXED) {
    for (i = 0; i < iinitvcnum; i++) {
      p_flag_array[i]     = FIXED;
      p_fix_flag_array[i] = EXTERNALLY;
    }
  }


  /* Assign values to q_flag_array.   */
  if (q_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      q_flag_array[i] = ESTIMATE;
      q_fix_flag_array[i] = non_value;
    } 
  }
  else if (q_flag == CONSTRAINT) {
    for (i = 0, j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	q_flag_array[i] = ESTIMATE;
	q_fix_flag_array[i] = non_value;
	j = i;
	if (i != 0)
	  k--;
      }
      else { 
	q_flag_array[i] = CONSTRAINT;
	q_fix_flag_array[i] = non_value;
      }
    }
  }
  else if (q_flag == FIXED) {
    for (i = 0; i < iinitvcnum; i++) {
      q_flag_array[i]     = FIXED;
      q_fix_flag_array[i] = EXTERNALLY;
    }
  }
}


/* Get estimated variance components G_VAL 's_vcnum', 'm1_vcnum', 'm2_vcnum', 
   't_vcnum' values here. They may be changed whenever this function is called.
*/
void Likelihoodfun::GetSMT_vcnum() {
  int i;
    
  s_vcnum = m1_vcnum = m2_vcnum = t_vcnum = 0;

  /* Assign values to s_vcnum.   */
  if (s_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      if (s_flag_array[i] == ESTIMATE)
	s_vcnum++;
    }
  }else if (s_flag == CONSTRAINT) {
    for (i = 0; i < iinitvcnum; i++) {
      if (s_flag_array[i] == ESTIMATE)
	s_vcnum++;
    }
  }

  /* Assign values to m1_vcnum.   */
  if (m1_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      if (m1_flag_array[i] == ESTIMATE)
	m1_vcnum++;
    }
  }else if (m1_flag == CONSTRAINT) {
    for (i = 0; i < iinitvcnum; i++) {
      if (m1_flag_array[i] == ESTIMATE)
	m1_vcnum++;
    }
  }

  /* Assign values to m2_vcnum.   */
  if (m2_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      if (m2_flag_array[i] == ESTIMATE)
	m2_vcnum++;
    }
  }else if (m2_flag == CONSTRAINT) {
    for (i = 0; i < iinitvcnum; i++) {
      if (m2_flag_array[i] == ESTIMATE)
	m2_vcnum++;
    }
  }

  /* Assign values to t_vcnum.   */
  if (t_flag == ESTIMATE) {
    for (i = 0; i < env_vcnum; i++) {
      if (t_flag_array[i] == ESTIMATE)
	t_vcnum++;
    }
  }else if (t_flag == CONSTRAINT) {
    //    if (datatype == LONGITUDINAL) {
      for (i = 0; i < env_vcnum; i++) {
	if (t_flag_array[i] == ESTIMATE)
	  t_vcnum++;
      }
      //    }
  }
}

/* Get estimated variance components G_VAL 'c_vcnum', 'p_vcnum', 'q_vcnum' 
   values here. They may be changed whenever this function is called.
*/
void Likelihoodfun::GetCPQ_vcnum() {
  int i;
    
  c_vcnum = p_vcnum = q_vcnum = 0;

  /* Assign values to c_vcnum.   */
  if (c_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      if (c_flag_array[i] == ESTIMATE)
	c_vcnum++;
    }
  }
  else if (c_flag == CONSTRAINT) {
    for (i = 0; i < iinitvcnum; i++) {
      if (c_flag_array[i] == ESTIMATE)
	c_vcnum++;
    }
  }

  /* Assign values to p_vcnum.   */
  if (p_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      if (p_flag_array[i] == ESTIMATE)
	p_vcnum++;
    }
  }
  else if (p_flag == CONSTRAINT) {
    for (i = 0; i < iinitvcnum; i++) {
      if (p_flag_array[i] == ESTIMATE)
	p_vcnum++;
    }
  }

  /* Assign values to q_vcnum.   */
  if (q_flag == ESTIMATE) {
    for (i = 0; i < iinitvcnum; i++) {
      if (q_flag_array[i] == ESTIMATE)
	q_vcnum++;
    }
  }
  else if (q_flag == CONSTRAINT) {
    for (i = 0; i < iinitvcnum; i++) {
      if (q_flag_array[i] == ESTIMATE)
	q_vcnum++;
    }
  }
}

/* Get the initial values of covariates 'init_x_beta' from the string.  */
void Likelihoodfun::GetInitCovValues(char x_str[]) {
  int    i,j,kv;
  int    len;
  double *cov_val;
  char   str[STRING], temp_str[STRING];


  /* Initialize the size of the pointers.    5-12-98; E.Y.  */
  cov_val = (double *) malloc((itraits*icovs) * sizeof(double));

  j  = 0;
  kv = 0;
  len = strlen(x_str);

  for (i = 0; i < len; i++) {
    /* for the first itraits*icovs  value(s) of beta.  */
    if ( kv < (itraits*icovs -1) ) { 
      if (x_str[i] != ',') { 
	str[j]   = x_str[i];
	str[j+1] = '\0';
	j++;
      }
      else { 
	strcpy(temp_str, str);
	cov_val[kv] = atof(temp_str);
	j = 0;
	kv++;
      }
    }
    /* for the last value of x_cov.  */
    else {
      str[j] = x_str[i];
      str[j+1] = '\0';
      j++;
    }
  }
  strcpy(temp_str, str);
  cov_val[kv] = atof(temp_str);

  /* put the values of cov_val[] to array beta[][].     */
  kv = 0;
  for (i = 0; i < itraits; i++) {
    for (j = 0; j < icovs; j++) {
      init_x_beta[i][j] = cov_val[kv]; 
      kv++;
    }
  }
  free(cov_val);
}


/* 
Build beta_vec vector by using beta[][] values.
*/
void Likelihoodfun::Build_beta_Vector(int nt, double *beta_vec) {
  int i,j;

  /* Build beta_vec  vector.         */
  j = 0;
  for (i = 0; i < icovs+1; i++) {
    if (i == 0) 
      beta_vec[i] = x_mu[nt];
    else {
      beta_vec[i] = x_beta[nt][j];
      j++;
    }
    /*
      printf("beta_vec[%d]=%lf\n", i, beta_vec[i]);
    */
  }
}

void Likelihoodfun::SaveProband(int family_member_number) {
  int    i;
  int    prob_sign;           /* for probands info.             */ 
  double prob_trait;

  prob_sign = 0;

  FortData *fd = &fortArray[dataIndex];

  // if the person is a proband
  //  if ((prob_sign == 1) || (prob_sign == 2)) { 
  if((fd->proband == 1) || (fd->proband == 2)) {
    for (i = 0; i < itraits; i++) {
      prob_sign = fd->proband;
      prob_trait = fd->traits[i]; 
      // . . . if trait is missing, => ERROR . . .
      if (prob_trait == trait_array[i].t_misval) {
	PROBLEM "The proband's trait value is missing, please fix your data file first.\nLikelihoodfun.cpp key 1065\n"
	  RECOVER(NULL_ENTRY);
      }
      // . . . save the trait value and number of the family member
      if (prob_sign == 1) {
	fam.proband[0].trait_p     = prob_trait;
	fam.proband[0].prob_mem    = family_member_number;
      }
      else if (prob_sign == 2) {
	fam.proband[1].trait_p     = prob_trait;
	fam.proband[1].prob_mem    = family_member_number;
      }
    }
  }
}

/*
Get the trait(s) and marker from each individual from family data 
file 'fort.12'.
*/
void Likelihoodfun::GetData(int n_fam_member) {
  int  i,j;

  for ( j = 0; j < n_fam_member; j++, dataIndex++) {
    /* print the first family's data into lookup log file and
       change G_VAL 'printfirstfam_flag = NO' here.            */
    if (printfirstfam_flag == YES) { 
      if (j == 0) {
	fp_lookuplog << endl << "FIRST FAMILY's DATA:" << endl;
      }
      fp_lookuplog << fortArray[dataIndex].familyId << ' '
		   << fortArray[dataIndex].seqId << ' '
		   << fortArray[dataIndex].fatherId << ' '
		   << fortArray[dataIndex].motherId << ' '
		   << fortArray[dataIndex].sex << ' ';
      for(i=0; i<itraits; i++) {
	fp_lookuplog << fortArray[dataIndex].traits[i] << ' ';
      }
      for(i=0; i<icovs; i++) {
	fp_lookuplog << fortArray[dataIndex].covariates[i] << ' ';
      }
      fp_lookuplog << endl;
    }

    if (ascert_flag == YES) {
      SaveProband(j);
    }

    for ( i = 0; i < total_trait_values; i++) {
      fam.person[j].trait[i] = fortArray[dataIndex].traits[i];
      if(fam.person[j].trait[i] == trait_array[i].t_misval) {
	fam.person[j].missing_traitflag = 1;
      }
      //      cout<<"trait="<<i<<" val="<<fam.person[j].trait[i]<<endl;      
    }
    for ( i = 0; i < total_cov_values; i++) {
      fam.person[j].cov[i] = fortArray[dataIndex].covariates[i];
      //      cout<<"cov"<<i<<"="<<fam.person[j].cov[i]<<endl;
    }
  } // end of 'for ( j = 0; j < n_fam_member; j++) {'

  // add by cjf Aug9, 2001 to adjust the missing value by the mean of nomissing from this family
  double fammean = 0.0;
  int count = 0;
  for ( j = 0; j < n_fam_member; j++) {
    for ( i = 0; i < total_cov_values; i++) {
      if(fam.person[j].cov[i] != cov_array[i].c_misval){
	//printf("ind %d cov=%f\n",j,fam.person[j].cov[i]);
	fammean += fam.person[j].cov[i];
	count++; 
      }
    }
  }
  fammean = fammean/count;
  for ( j = 0; j < n_fam_member; j++) {
    for ( i = 0; i < total_cov_values; i++) {
      if(fam.person[j].cov[i] == cov_array[i].c_misval){
	//printf("ind %d cov=%f\n",j,fammean);
	fam.person[j].cov[i] = fammean;
      }
    }
  }

  printfirstfam_flag = NO;
}

/* Get z1_mat, z2_mat matrices from file fp_loci.
   changed 9/17/1997.
*/
void Likelihoodfun::Get_z_mat(int n_fam_member, double **z1_mat, double **z2_mat) {
  int    i,j;
  // stringHolder replaces mg1_f1, mg1_f2, mg2_f1, and mg2_f2
  char stringHolder[1024];
  double /* mg1_f1, mg1_f2, */ mg1_pi;
  double /* mg2_f1, mg2_f2, */ mg2_pi;
  char   str[BUF];

  for ( i = 0; i < n_fam_member; i++) {
    for ( j = i; j < n_fam_member; j++) {
      if (j > i) {
	fp_loci.getline(str, BUF);
	if (iloci == 0) {
	  PROBLEM "multic.par must specify the number of loci to be 1 or 2.\nLikelihoodfun.cpp key 1300\n"
	    RECOVER(NULL_ENTRY);
	}else if (iloci == 1) {
	  sscanf(str,"%s %s %lf ", stringHolder, stringHolder, &mg1_pi);
	  if (mg1_pi == non_value) { 
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_pi = 0.5;                
	  }
	  z1_mat[i][j] = mg1_pi;
	  z1_mat[j][i] = mg1_pi;
	}
	else if (iloci == 2) {
	  sscanf(str,"%s %s %lf %s %s %lf ", 
		 stringHolder, stringHolder, &mg1_pi, 
		 stringHolder, stringHolder, &mg2_pi);
	  if (mg1_pi == non_value) { 
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg1_pi = 0.5;                
	  }
	  z1_mat[i][j] = mg1_pi;
	  z1_mat[j][i] = mg1_pi;

	  if (mg2_pi == non_value) { 
	    /* changed 7/31/96;  if marker is missing for the sib.  */
	    mg2_pi = 0.5;                
	  }
	  z2_mat[i][j] = mg2_pi;
	  z2_mat[j][i] = mg2_pi;
	}
	else {
	  cerr << "*** Get_z_mat debug ***" << endl <<
	    "Error: can't handle more than two loci," << endl <<
	    "please change the parameter file." << endl <<
	    "(Likelihoodfun.cpp line 1327)" << endl;
	}
      }
    }
  }
}

/* Get g_mat, c_mat, p_mat, q_mat matrices from file fp_share.
   changed 9/11/97.

Updated:
7-30-03, Eric Lunde, This method previously read from share.out many many
         times.  To speed up the program, I stored all of share.out internally
         in an array.  Now, this method reads from that array.
*/
void Likelihoodfun::Get_gcpq_mat(int n_fam_member,
				 double **g_mat, double **c_mat,
				 double **p_mat, double **q_mat) {
  int    i,j;

  for ( i = 0; i < n_fam_member; i++) {
    for ( j = i+1; j < n_fam_member; j++) {
      g_mat[i][j] = g_mat[j][i] = shareArray[relationIndex].geneticSimilarity;
      
      c_mat[i][j] = c_mat[j][i] = shareArray[relationIndex].areSiblings;
      
      p_mat[i][j] = p_mat[j][i] = shareArray[relationIndex].areSpouses;
      
      q_mat[i][j] = q_mat[j][i] = shareArray[relationIndex].areParentChild;

      relationIndex++;
    }
  }
}

/* Assign the initial values to all the parameters except for the major gene
   components which already be fixed internally.
   10-30-96; Added by Emily.
*/
void Likelihoodfun::AssignInitValues() {
  int i,j;

  mu_flag = init_mu_flag;
  s_flag  = init_s_flag;
  m1_flag = init_m1_flag;
  m2_flag = init_m2_flag;
  t_flag  = init_t_flag;
  c_flag  = init_c_flag;
  p_flag  = init_p_flag;
  q_flag  = init_q_flag;

  for (i = 0; i < iinitvcnum; i++) {
    x_S[i]       = init_x_S[i];
    if (m1_fix_flag_array[i] != INTERNALLY) {
      x_M1[i]       = init_x_M1[i];
    }
    if (m2_fix_flag_array[i] != INTERNALLY) {
      x_M2[i]       = init_x_M2[i];
    }
  }
  for (i = 0; i < env_vcnum; i++) {
    x_T[i]       = init_x_T[i];
  }
  for (i = 0; i < itraits; i++) {
    x_mu[i]      = init_x_mu[i];
    for (j = 0; j < icovs; j++) {
      x_beta[i][j] = init_x_beta[i][j]; 
    }
  }

  GetST_flag_array();
  GetCPQ_flag_array();
}


/* 
Build lD matrices and the covariate columns ld_vec of lD matrices.
*/
void Likelihoodfun::Build_lD_Matrices(int n_fam_member, double ***lD_mat,
				      double ***ld_vec) {
  int i,j,k,r;
  r = 0;
  int num_rptmst;

    
  num_rptmst = 0;
  /* Build lD matrices.         */
  for (i = 0; i < irepeatmst; i++) {
    for (j = 0; j < n_fam_member; j++) {
      r = num_rptmst;
      for (k = 0; k < (icovs+1); k++) {
	if (k == 0)
	  lD_mat[i][j][k] = 1;
	else { 
	  lD_mat[i][j][k] = fam.person[j].cov[r];
	  r++;
	}
      }
    }
    num_rptmst = r;
  }
  num_rptmst = 0;
  /* Build the covariate columns ld_vec of lD matrices.    */
  for (i = 0; i < irepeatmst; i++) {
    r = num_rptmst;
    for (k = 0; k < icovs; k++) {
      for (j = 0; j < n_fam_member; j++) { 
	ld_vec[i][k][j] = fam.person[j].cov[r]; 
      }
      r++;
    }
    num_rptmst = r;
  }
}

/* Change ld_vec if there is any missing trait value in the family members.
*/
void Likelihoodfun::LChangeVec(double ***ld_vec,int n_fam_member) {
  int i,j,k,kk,numofcov;
  int prev_count;
  int n_mb, n_order;

  prev_count = 0;
  n_mb       = n_fam_member;

  for (k = 0; k < n_fam_member; k++) {
    if (fam.person[k].missing_traitflag == 1) {
      n_order = k - prev_count;
      /*
	printf("n_order=%d, n_mb=%d\n", n_order, n_mb);
      */
   
      kk = 0;
      for(i = 0; i < n_mb-1; i++) {
	if (i == n_order) 
	  kk = i+1;
	for ( numofcov=0; numofcov < icovs; numofcov++) {
	  for (j = 0; j < irepeatmst; j++)
	    ld_vec[j][numofcov][i] = ld_vec[j][numofcov][kk]; 
	}
	kk++;
      }
      prev_count++;
      n_mb--;
    }
  }
}

/* Get (itraits*n_fam_member) MuBeta_vec vector by using  i_vec, ld_vec vectors. 
*/
void Likelihoodfun::LGetMuBetaVector(double *MuBeta_vec[], double *i_vec, double ***ld_vec, int n_member) {
  int i,j,k,kv;

  kv = 0;

  for (j = 0; j < (icovs+1); j++) {
    for (i = 0; i < irepeatmst; i++) {
      if (j == 0) {
	for(k = 0; k < n_member; k++) {
	  MuBeta_vec[kv][i*n_member+k] = i_vec[k];
	}
      }
      else {
	for(k = 0; k <n_member ; k++) {
	  MuBeta_vec[kv][i*n_member+k] = ld_vec[i][j-1][k];
	}
      }
    }
    kv++;
  }
}
