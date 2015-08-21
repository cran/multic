/*
 *-Produced by PROMULA.FORTRAN to C V7.16 on 09/24/99 at 11:31:25
 *-With heavy modification by John Venier
 */

#include "Maxfun_source.h"
#include <cmath>
#include <cstring>
using namespace std;

/*
Maxfun_exp Exp = {
    stderr, NULL,
    { 0.0 }, { 0.0 }, { 0.0 }, { 0.0 },
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    { 0 }, 0, 0, 0, 0, 0, 1,
    { "" }
};
*/

/*
 * File Static Scalars and Arrays (private global variables)
 */
static const char Xstatus[10][37]={"IND-FN, MAY VARY",
                                   "IND, MAY VARY",
                                   "DEPENDENT",
                                   "FIXED EXTERNALLY",
                                   "IND-FN, FIXED BY MAXFUN AT BOUND",
                                   "IND, FIXED BY MAXFUN AT BOUND",
                                   "IND-FN, FIXED BY MAXFUN NEAR BOUN\
D",
                                   "IND, FIXED BY MAXFUN NEAR BOUND",
                                   "IND-FN, FIXED BY MAXFUN NOT NEAR \
BND",
                                   "IND, FIXED BY MAXFUN NOT NEAR BOU\
ND"};
static double Xstp[MAXFUN_NP],Xg[MAXFUN_NP],Xthpr[MAXFUN_NP],Xfpr,
  Xfch,Xdifmax,Xv[MAXFUN_NP][MAXFUN_NP],Xpdir[MAXFUN_NP],Xptg,
  Xtstep,Xh[MAXFUN_NP][MAXFUN_NP],Xult[MAXFUN_NP][MAXFUN_NP],
  Xdiag[MAXFUN_NP],Xcth[MAXFUN_NP],Xerm,Xav[MAXFUN_NP][MAXFUN_NP],
  Xstde[MAXFUN_NP],Xgtg;
static int Xist[MAXFUN_NP],Xne,Xnd,Xni,Xnv,Ximpbnd,Xit,Xnsurf2,Xigfl,
  Xivfl,Xigage,Xivage,Xidif,Xnb;

/*
 * Exported routine
 */

/*
 * *******************************************************************
 */

Maxfun_source::Maxfun_source(Maxfun_fun* mfun1){
  mfun = mfun1;
  Exp = mfun->getExp();
}

void Maxfun_source::maxfun(double *theta, double *f,
			   int *nfe, int *lfl){

#define PRECIS pow(2.e0,(double) - 45)
#define TAU 0.5e0

  /*
 * Local Scalars
 */
  static double si,thi;
  static int i,ifirst,ifl,iflb,ih,ihess,isti,iupdt,l,lex;
  static int K1 = 0;
  static int K2 = 1;

  /*
   * Local Arrays
   */
  static double dth[MAXFUN_NP],gpr[MAXFUN_NPV];

  /*
 * Format strings
 */
  const char* F9000 = "\n------------------------------------------\
----------------------------\n---------------------------------------\
-------------------------------\n";
  const char* F9100 = "MAXFUN              FEB 1997\n";
  const char* F9200 = "\nUNUSED Exp.iout =%#X\n";
  const char* F9300 = "\n";
  const char* F9400 = "UNUSED Exp.idet =%#X\n";
  const char* F9500 = "INVALID Exp.method =%4d\n";
  const char* F9600 = "INVALID Exp.nt =%4d\n";
  const char* F9700 = "INVALID Exp.istin(%2d) =%4d\n";
  const char* F9800 = "ALL PARAMETERS FIXED\n";
  const char* F9900 = "ALL PARAMETERS FIXED OR DEPENDENT\n";
  const char* F10000 = "\nMAXIMIZATION NOT ATTEMPTED DUE TO IMPROPE\
R CONTROL INPUT\n";
  const char* F10100 = "ITERATION DETAILS WILL BE WRITTEN IN FILE C\
ORRESPONDING TO\n   UNIT %#X\n";
  const char* F10200 = "TOTAL # PARAMETERS IN THE MODEL NT =\
%3d\n   # PARAMETERS TO BE ESTIMATED NE =%3d\n   # DEPENDENT PARAMETE\
RS       ND =%3d\n   MAXIMUM # ITERATIONS MAXIT =%4d\n\nCONTROL VALUE\
S:    IXVC  =%2d;         IHIT =%2d\n   EPSD  =%#10.3E; YOTA  =\
%#10.3E; EPST  =%#10.3E\n   EPSC1 =%#10.3E; EPSC2 =%#10.3E; EPSC3 =\
%#10.3E\n------------------------------------------------------------\
----------\n";
  const char* F10300 = "\nMAXIMIZATION NOT ATTEMPTED:\n   INITIAL E\
STIMATES NOT IN DOMAIN OF THE FUNCTION\n";
  const char* F10400 = "\nMETHOD 1:  DIRECT SEARCH METHOD\n";
  const char* F10500 = "\nSTOPPED AT%4d ITERATIONS (REACHED ITERATI\
ON LIMIT)\n";
  const char* F10600 = "\n2**(NI)-TRIAL SEARCH PRODUCED SIGNIFICAN\
T IMPROVEMENT TWICE,\n   NOT TRIED A THIRD TIME;\n   DIRECT SEARCH IT\
ERATION METHOD MAY NOT BE APPROPRIATE\n";
  const char* F10700 = "\nMETHOD 2:  DIRECT SEARCH METHOD, WITHOUT \
2**(NI)-TRIAL SEARCH\n";
  const char* F10800 = "\nMETHOD 3:  NEWTON-RAPHSON METHOD, WITHOU\
T RECOMPUTATION OF\n           VARIANCE-COVARIANCE MATRIX\n";
  const char* F10900 = "\nMETHOD 4:  NEWTON-RAPHSON METHOD\n";
  const char* F11000 = "\nMETHOD 5:  VARIABLE METRIC METHOD, USING \
NEGATIVE IDENTITY FOR INITIAL\n           HESSIAN\n";
  const char* F11100 = "\nSTOPPED AFTER ITERATION%4d BECAUSE SEARC\
H DIRECTION IS NOT UPWARDS\n";
  const char* F11300 = "-------------------------------------------\
---------------------------\nSWITCH TO CENTRAL DIFFERENCE IN GRADIEN\
T COMPUTATION\n------------------------------------------------------\
----------------\n";
  const char* F11200 = "\nSTOPPED AFTER ITERATION%4d BECAUSE ACCUMU\
LATION OF ROUNDING ERRORS\n   PREVENTS FURTHER PROGRESS\n";
  const char* F11400 = "\nMETHOD 6:  VARIABLE METRIC METHOD, COMPUT\
ING INITIAL HESSIAN\n";
  const char* F11500 = "\nVARIANCE-COVARIANCE MATRIX WILL BE COMPUT\
ED TO CORRESPOND TO FINAL\n   ESTIMATES\n";
    
  /*
 * Begin executable statements
 */

  /*
   *--INITIALIZE EXIT FLAG
   */
  *lfl = 0;

  /*
 *--PRINT INITIAL MESSAGES AND CHECK INPUT
 */
  if(Exp.iout == NULL) {
    (void)Rprintf(F9000);
    (void)Rprintf(F9100);
    (void)Rprintf(F9200,Exp.iout);
  } else {
    (void)fprintf(Exp.iout,"%s",F9000);
    (void)fprintf(Exp.iout,"%s",F9100);
    copyr(Exp.iout);
    (void)fprintf(Exp.iout,"%s",F9300);
  }

  if(Exp.idet == NULL) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9400,Exp.idet);
    else (void)Rprintf(F9400,Exp.idet);
  }
  else if(Exp.idet != Exp.iout) {
    (void)fprintf(Exp.idet,"%s",F9000);
    (void)fprintf(Exp.idet,"%s",F9100);
    copyr(Exp.idet);
    (void)fprintf(Exp.idet,"%s",F9300);
  }

  if(Exp.method < 1 || Exp.method > 6) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9500,Exp.method);
    *lfl = 13;
  }

  if(Exp.nt > 0 && Exp.nt <= MAXFUN_NP) goto S10;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9600,Exp.nt);
  *lfl = 13;
  goto S70;

 S10:
  Xne = Xnd = 0;

  for(i=1; i<=Exp.nt; i++) {

    //printf("i=%d %f\n",i, Exp.thin[i-1]);

    isti = Exp.istin[i - 1];

    switch(isti) {
    case 1: goto S30;
    case 2: goto S30;
    case 3: goto S20;
    case 4: goto S40;
    default: break;
    }

    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9700,i,isti);
    *lfl = 13;
    goto S40;

  S20:
    Xnd += 1;

  S30:
    Xne += 1;

  S40:;

  }

  if(Xne > 0) goto S50;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9800);
  *lfl = 13;
  goto S70;

 S50:
  if(Xnd < Xne) goto S60;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9900);
  *lfl = 13;
  goto S70;
   
 S60:
  if(*lfl == 0) goto S80;
   
 S70:
  Xigfl = 2;
  if(Exp.iout == NULL) return;
  (void)fprintf(Exp.iout,"%s",F10000);
  (void)fprintf(Exp.iout,"%s",F9000);
  return;

 S80:
  /*
 *--SUBSTITUTE DEFAULTS IF NECESSARY
 */
  if(Exp.ixvc < 0 || Exp.ixvc > 2) Exp.ixvc = 0;
  if(Exp.ihit < 0 || Exp.ihit > 1) Exp.ihit = 0;
  if(Exp.epsc1 <= 0.e0 || Exp.epsc1 >= TAU) Exp.epsc1 = 1.e-3;
  if(Exp.epsc2 <= 0.e0) Exp.epsc2 = 1.e-15;
  Exp.epsc3 = MAX(0.e0,Exp.epsc3);
  if(Exp.epsd <= 1.e-9 || Exp.epsd >= TAU) {
    if(Exp.epsc1 > 1.e-9) Exp.epsd = Exp.epsc1;
    else Exp.epsd = 1.e-3;
  }
  if(Exp.yota <= 0.e0 || Exp.yota >= 1.e0)
    Exp.yota = 10.e0 * sqrt(PRECIS);
  if(Exp.epst <= 0.e0) Exp.epst = 10.e0 * sqrt(Exp.yota);
  for(i=1; i<=Exp.nt; i++) {
    si = Exp.stpin[i - 1];
    if(si <= 0.e0 || si >= 1.e0) Exp.stpin[i - 1] = .1e0;
  }

  /*
 *--PRINT OUT CONTROL VALUES
 */
  if(Exp.iout != NULL) {
    if(Exp.idet != NULL) (void)fprintf(Exp.iout,F10100,Exp.idet);
    (void)fprintf(Exp.iout,F10200,Exp.nt,Xne,Xnd,Exp.maxit,
		  Exp.ixvc,Exp.ihit,Exp.epsd,Exp.yota,Exp.epst,
		  Exp.epsc1,Exp.epsc2,Exp.epsc3);
  }

  //here
  //printf("lfl=%d\n",*lfl);

  /*
 *--GENERAL INITIALIZATION
 */
  Xigage = Xivage = -1;
  Xigfl = 0;
  Xivfl = 5;
  ih = Exp.ihit;
  Xnsurf2 = Ximpbnd = Xit = *nfe = 0;
  Xni = Xne - Xnd;
  for(i=1; i<=Exp.nt; i++) {
    isti = Exp.istin[i - 1];
    Xist[i - 1] = isti;
    Xstp[i - 1] = Exp.stpin[i - 1];
    thi = Exp.thin[i - 1];
    if(isti != 3) {
      if(thi < Exp.thl[i - 1] || thi > Exp.thu[i - 1]) {
	//printf("i=%d\n",i);
	//printf("thi=%d thl=%d thu=%d\n",thi,Exp.thl[i-1],Exp.thu[i-1]);
	*lfl = 12;
      }
    }
    theta[i - 1] = thi;
  }

  //printf("lfl=%d\n",*lfl);

  if(*lfl > 0) goto S90;

  //printf("lfl=%d\n",*lfl);

  mfun->eval_fun(theta,f,nfe,&lex);

  if(lex <= 0) goto S100;

  *lfl = 12;

 S90:
  /*
   *--INITIAL ESTIMATES INVALID
   */
  Xigfl = 2;
  if(Exp.iout == NULL) return;
  (void)fprintf(Exp.iout,"%s",F10300);
  lbv(theta,f,&K1);
  (void)fprintf(Exp.iout,"%s",F9000);
  return;

 S100:
  /*
 *--CHECK FOR POSSIBLE ZERO-ITERATION RUN
 */
  if(Exp.maxit > 0) goto S120;

  prepd(theta,f,nfe,&lex);

  if(lex <= 0) goto S110;

  Xigfl = 2;
  goto S730;

 S110:
  if(Exp.ixvc <= 0) goto S700;

  goto S680;

 S120:
  /*
 *--SELECT PATH FOR CHOSEN METHOD
 */
  switch(Exp.method) {
  case 1: goto S130;
  case 2: goto S200;
  case 3: goto S250;
  case 4: goto S360;
  case 5: goto S460;
  case 6: goto S670;
  default: break;
  }

 S130:
  /*
   *--METHOD 1:  DIRECT SEARCH METHOD, INCLUDING BASIC SEARCH AND
   *--2**(NI)-TRIAL SEARCH
   *
   *--PRINT METHOD IDENTIFICATION AND INITIAL VALUES
   */
  if(Exp.idet != NULL && Exp.idet != Exp.iout)
    (void)fprintf(Exp.idet,"%s",F10400);
  if(Exp.iout != NULL) {
    (void)fprintf(Exp.iout,"%s",F10400);
    lbd(theta,f);
  }

 S140:
  /*
   *--START BASIC ITERATION PROCESS
   */
  bsrch(theta,f,nfe,&iflb);

  /*
 *--IF REACHED MAXIMUM # ITERATIONS, EXIT
 */
  if(iflb < 3) goto S150;

  *lfl = 4;
  goto S690;

 S150:
  /*
   *--SKIP OTHER TYPES OF SEARCH IF ONLY ONE PARAMETER TO VARY
   */
  if(Xni <= 1) goto S180;

  /*
 *--SEARCH 2(NI) "NEIGHBORING" PARAMETER SETS
 */
  nsrch(theta,f,nfe,dth,&ifl);
  if(ifl > 0) goto S160;

  /*
   *--SIGNIFICANT IMPROVEMENT; IF < MAXIMUM # ITERATIONS, GO BACK TO
   *--BASIC ITERATION PROCESS
   */
  if(Xit < Exp.maxit) goto S140;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F10500,Xit);
  *lfl = 4;
  goto S690;

 S160:
  /*
 *--CHECK WHETHER HAVE DONE 2**(NI)-TRIAL SEARCH MORE THAN ONCE
 *--ALREADY
 */
  if(Xnsurf2 <= 1) goto S170;

  Xnsurf2 = 3;
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F10600);
  goto S180;

 S170:
  /*
 *--SEARCH 2**(NI) NEIGHBORING PARAMETER SETS;
 *--IF SIGNIFICANT IMPROVEMENT, GO BACK TO BASIC SEARCH PROCESS
 */
  psrch(theta,f,nfe,dth,&ifl);
  if(ifl <= 0) goto S140;

 S180:
  /*
 *--FINISH UP
 */
  *lfl = 1;
  if(iflb > 1) *lfl = 3;
  prepd(theta,f,nfe,&lex);
  if(lex <= 0) goto S190;

  Xigfl = 2;
  goto S690;

 S190:
  if(Exp.ixvc <= 0) goto S690;

  goto S680;

 S200:

  /*
   *--METHOD 2:  DIRECT SEARCH METHOD, SKIPPING 2**(NI)-TRIAL SEARCH
   *
   *--PRINT METHOD IDENTIFICATION AND INITIAL VALUES
   */

  if(Exp.idet != NULL && Exp.idet != Exp.iout)
    (void)fprintf(Exp.idet,"%s",F10700);

  if(Exp.iout != NULL) {
    (void)fprintf(Exp.iout,"%s",F10700);
    lbd(theta,f);
  }
   
 S210:
  /*
   *--START BASIC ITERATION PROCESS
   */
  bsrch(theta,f,nfe,&iflb);

  /*
 *--IF REACHED MAXIMUM # ITERATIONS, EXIT
 */
  if(iflb < 3) goto S220;

  *lfl = 4;
  goto S690;

 S220:
  /*
   *--SKIP OTHER TYPE OF SEARCH IF ONLY ONE PARAMETER TO VARY
   */
  if(Xni <= 1) goto S230;

  /*
 *--SEARCH 2(NI) "NEIGHBORING" PARAMETER SETS
 */
  nsrch(theta,f,nfe,dth,&ifl);
  if(ifl > 0) goto S230;

  /*
   *--SIGNIFICANT IMPROVEMENT; IF < MAXIMUM # ITERATIONS, GO BACK TO
   *--BASIC ITERATION PROCESS
   */
  if(Xit < Exp.maxit) goto S210;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F10500,Xit);
  *lfl = 4;
  goto S690;

 S230:
  /*
 *--FINISH UP
 */
  *lfl = 1;
  if(iflb > 1) *lfl = 3;
  prepd(theta,f,nfe,&lex);
  if(lex <= 0) goto S240;

  Xigfl = 2;
  goto S690;

 S240:
  if(Exp.ixvc <= 0) goto S690;

  goto S680;

 S250:
  /*
   *--METHOD 3:  NEWTON-RAPHSON METHOD WITHOUT RECOMPUTATION OF
   *--VARIANCE-COVARIANCE MATRIX
   *
   *--PRINT METHOD IDENTIFICATION AND INITIAL VALUES
   */
  if(Exp.idet != NULL && Exp.idet != Exp.iout)
    (void)fprintf(Exp.idet,"%s",F10800);
  if(Exp.iout != NULL) {
    (void)fprintf(Exp.iout,"%s",F10800);
    lbv(theta,f,&K2);
  }

  /*
 *--PREPARE FOR NEWTON-RAPHSON ITERATION
 */
  prepd(theta,f,nfe,&lex);
  if(lex > 0) {
    *lfl = lex + 9;
    Xigfl = 2;
    return;
  }

  /*
 *--USE CENTRAL DIFFERENCE FOR GRADIENT CALCULATIONS
 */
  Xidif = 2;

 S260:
  /*
   *--COMPUTE NEW VARIANCE-COVARIANCE MATRIX FIRST TIME OR IF PARAMETERS
   *--BECOME FIXED
   */
  vcmx(theta,f,nfe,&Exp.ihit);
  if(Xivfl < 3) goto S270;

  *lfl = 9;
  goto S690;

 S270:
  /*
 *--COMPUTE GRADIENT VECTOR
 */
  deriv1(theta,f,nfe);
  if(Xigfl <= 0) goto S280;

  *lfl = 8;
  goto S690;

 S280:
  /*
 *--DO ONE NEWTON-RAPHSON ITERATION
 */
  nrstep(theta,f,nfe,&ifl);

  /*
 *--GO BACK TO COMPUTE NEW G UNLESS MET SOME STOPPING CRITERION
 */
  switch(ifl + 1) {
  case 1: goto S290;
  case 2: goto S310;
  case 3: goto S320;
  case 4: goto S310;
  case 5: goto S350;
  default: break;
  }

 S290:
  /*
   *--CHECK FOR CONVERGENCE TO BOUNDS BEFORE NEXT ITERATION
   */
  fixbnd(theta,f,nfe,&lex);

  switch(lex + 1) {
  case 1: goto S270;
  case 2: goto S260;
  case 3: goto S300;
  default: break;
  }

 S300:
  *lfl = 10;
  goto S690;

 S310:
  /*
 *--CONVERGED; IF ITERATION DONE, CHECK FOR CONVERGENCE TO BOUNDS
 *--BEFORE QUITTING
 */
  fixbnd(theta,f,nfe,&lex);

 S320:
  *lfl = ifl;

  switch(Exp.ixvc + 1) {
  case 1: goto S690;
  case 2: goto S340;
  case 3: goto S330;
  default: break;
  }

 S330:
  if(Exp.ihit <= 0) goto S680;

 S340:
  if(Exp.ixvc + Xivage <= 2) goto S690;

  goto S680;

 S350:
  /*
 *--REACHED MAXIMUM # ITERATIONS
 */
  *lfl = 4;
  goto S690;

 S360:
  /*
 *--METHOD 4:  NEWTON-RAPHSON METHOD WITH RECOMPUTATION OF VARIANCE-
 *--COVARIANCE MATRIX AT EACH ITERATION
 *
 *--PRINT METHOD IDENTIFICATION AND INITIAL VALUES
 */
  if(Exp.idet != NULL && Exp.idet != Exp.iout)
    (void)fprintf(Exp.idet,"%s",F10900);
  if(Exp.iout != NULL) {
    (void)fprintf(Exp.iout,"%s",F10900);
    lbv(theta,f,&K2);
  }

  /*
 *--PREPARE FOR NEWTON-RAPHSON ITERATION
 */
  prepd(theta,f,nfe,&lex);
  if(lex > 0) {
    *lfl = lex + 9;
    Xigfl = 2;
    return;
  }

  /*
 *--USE CENTRAL DIFFERENCE FOR GRADIENT CALCULATIONS
 */
  Xidif = 2;

 S370:
  /*
   *--COMPUTE NEW VARIANCE-COVARIANCE MATRIX AND VECTOR OF 1ST PARTIAL
   *--DERIVATIVES
   */
  vcmx(theta,f,nfe,&Exp.ihit);
  if(Xivfl < 3) goto S380;

  *lfl = 9;
  goto S690;

 S380:
  deriv1(theta,f,nfe);
  if(Xigfl <= 0) goto S390;

  *lfl = 8;
  goto S690;

 S390:
  /*
   *--DO ONE NEWTON-RAPHSON ITERATION
   */
  nrstep(theta,f,nfe,&ifl);

  /*
 *--GO BACK TO COMPUTE NEW V AND G UNLESS MET SOME STOPPING CRITERION
 */
  switch(ifl + 1) {
  case 1: goto S400;
  case 2: goto S410;
  case 3: goto S420;
  case 4: goto S410;
  case 5: goto S450;
  default: break;
  }

 S400:
  /*
   *--CHECK FOR CONVERGENCE TO BOUNDS BEFORE NEXT ITERATION
   */
  fixbnd(theta,f,nfe,&lex);
  if(lex <= 1) goto S370;

  *lfl = 10;
  goto S690;

 S410:
  /*
 *--CONVERGED; IF ITERATION DONE, CHECK FOR CONVERGENCE TO BOUNDS
 *--BEFORE QUITTING
 */
  fixbnd(theta,f,nfe,&lex);

 S420:
  *lfl = ifl;

  switch(Exp.ixvc + 1) {
  case 1: goto S690;
  case 2: goto S440;
  case 3: goto S430;
  default: break;
  }

 S430:
  if(Exp.ihit <= 0) goto S680;

 S440:
  if(Exp.ixvc + Xivage <= 2) goto S690;

  goto S680;

 S450:
  /*
 *--REACHED MAXIMUM # ITERATIONS
 */
  *lfl = 4;
  goto S690;

 S460:
  /*
 *--METHOD 5:  VARIABLE METRIC METHOD WITH INITIAL B = IDENTITY
 *
 *--PRINT METHOD IDENTIFICATION
 */
  if(Exp.idet != NULL && Exp.idet != Exp.iout)
    (void)fprintf(Exp.idet,"%s",F11000);
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F11000);

  /*
 *--SET SWITCH FOR INITIAL B:  DON'T COMPUTE HESSIAN
 */
  ihess = 0;

 S470:
  /*
   *--METHOD 6 JOINS PATH AT THIS POINT
   *
   *--PRINT INITIAL VALUES
   */
  if(Exp.iout != NULL) lbv(theta,f,&K2);

  /*
 *--PREPARE FOR DERIVATIVE COMPUTATION
 */
  prepd(theta,f,nfe,&lex);
  if(lex > 0) {
    *lfl = lex + 9;
    Xigfl = 2;
    return;
  }

 S480:
  /*
   *--INITIALIZE B MATRIX (ULT, DIAG) AND FLAGS
   */
  binit(theta,f,nfe,&ihess);
  Xidif = ifirst = 1;

 S490:
  iupdt = 0;

 S500:
  /*
 *--START ITERATION
 *
 *--COMPUTE NEW GRADIENT G, FIRST SAVING OLD ONE FOR BUPDT
 *
 */
  for(l=0; l<Xnv; l++) gpr[l] = Xg[l];
  deriv1(theta,f,nfe);
  if(Xigfl <= 0) goto S510;

  *lfl = 8;
  goto S690;

 S510:
  /*
   *--IF APPROPRIATE, UPDATE B (ULT, DIAG)
   */
  if(iupdt <= 0) goto S520;

  bupdt(gpr,&lex);
  if(lex <= 0) goto S520;

  *lfl = 7;
  goto S660;

 S520:
  /*
   *--COMPUTE NEW DIRECTION OF SEARCH
   */
  direct(&lex);

  switch(lex + 1) {
  case 1: goto S550;
  case 2: goto S530;
  case 3: goto S540;
  default: break;
  }

 S530:
  *lfl = 2;
  goto S660;

 S540:
  if(Xidif == 1) goto S650;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F11100,Xit);
  *lfl = 6;
  goto S660;

 S550:
  /*
   *--DO LINE SEARCH IN CHOSEN DIRECTION TO COMPLETE ONE ITERATION
   */
  lsrch(theta,f,nfe,&ifirst,&ifl);

  switch(ifl + 2) {
  case 1: goto S560;
  case 2: goto S570;
  case 3: goto S600;
  case 4: goto S610;
  case 5: goto S630;
  case 6: goto S640;
  default: break;
  }

 S560:
  /*
   *--ITERATION COMPLETED, BUT T < INITIAL T AND < TMIN
   */
  iupdt = 0;
  if(Xidif == 2) goto S580;

  if(Exp.idet != NULL) (void)fprintf(Exp.idet,"%s",F11300);
  Xidif = 2;
  goto S580;

 S570:
  /*
   *--ITERATION SUCCESSFULLY COMPLETED  WITH T >= INITIAL T AND/OR TMIN
   */
  iupdt = 1;

 S580:
  /*
   *--CHECK FOR CONVERGENCE TO BOUNDS BEFORE NEXT ITERATION
   */
  fixbnd(theta,f,nfe,&lex);

  switch(lex + 1) {
  case 1: goto S500;
  case 2: goto S480;
  case 3: goto S590;
  default: break;
  }

 S590:
  *lfl = 10;
  goto S690;

 S600:
  /*
 *--CONVERGED BY STANDARD TEST
 */
  *lfl = 1;
  goto S620;

 S610:
  /*
 *--NEGLIGIBLE FUNCTION CHANGE
 */
  *lfl = 3;

 S620:
  /*
   *--CHECK FOR CONVERGENCE TO BOUNDS BEFORE QUITTING
   */
  fixbnd(theta,f,nfe,&lex);
  goto S660;

 S630:
  /*
 *--REACHED MAXIMUM # ITERATIONS
 */
  *lfl = 4;
  goto S690;

 S640:
  /*
 *--COULD NOT COMPLETE ITERATION
 */
  if(Xidif == 1) goto S650;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F11200,Xit);
  *lfl = 5;
  goto S660;

 S650:
  /*
 *--SWITCH TO CENTRAL DIFFERENCE COMPUTATION OF GRADIENT
 *--AND RETRY THIS ITERATION
 */
  if(Exp.idet != NULL) (void)fprintf(Exp.idet,"%s",F11300);
  Xidif = 2;
  goto S490;

 S660:
  /*
   *--FINISH UP
   */
  if(Exp.ixvc > 0) goto S680;

  goto S690;

 S670:
  /*
 *--METHOD 6:  VARIABLE METRIC METHOD WITH INITIAL B = -H
 *--           COMPUTED AT INITIAL ESTIMATES
 *
 *--PRINT METHOD IDENTIFICATION
 */
  if(Exp.idet != NULL && Exp.idet != Exp.iout)
    (void)fprintf(Exp.idet,"%s",F11400);
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F11400);

  /*
 *--SET SWITCH TO COMPUTE INITIAL HESSIAN
 */
  ihess = 1;

  /*
 *--JOIN METHOD 5 PATH
 */
  goto S470;

 S680:
  /*
   *--FINAL OUTPUT AND TERMINATION
   *
   *--OBTAIN VARIANCE-COVARIANCE MATRIX FOR FINAL ESTIMATES
   */
  if(Exp.idet != NULL) (void)fprintf(Exp.idet,"%s",F11500);
  if(Exp.ixvc >= 2) ih = 1;
  if(ih > Exp.ihit) {
    for(i=1; i<=Exp.nt; i++) Xstp[i - 1] = Exp.epsd;
  }
  vcmx(theta,f,nfe,&ih);

 S690:
  /*
 *--AUGMENT AND PRINT FINAL VARIANCE-COVARIANCE MATRIX (IF ONE IS
 *--AVAILABLE) AND COMPUTE STANDARD DEVIATIONS
 */
  if(Xivfl < 3) augv(theta,&ih);

 S700:
  if(Xigage > 0) goto S720;

  if(Xigage == 0) goto S730;

  if(Xigfl > 0) goto S730;

 S720:
  Xidif = 2;
  deriv1(theta,f,nfe);

 S730:
  /*
 *--PRINT FINAL VALUES
 */
  if(Exp.iout != NULL) lf(theta,f,nfe);
  return;

  /*
   * Undefine Local Preprocessor Macros
   */
#undef PRECIS
#undef TAU
}

/*
 * File static routines
 */

/*
 * *******************************************************************
 */

void Maxfun_source::bsrch(double theta[],double *f,int *nfe,int *ifl) {
  /*
   *====================================================================
   *
   *--BASIC ITERATION PROCESS OF DIRECT SEARCH
   *
   *--RETURNS FLAG IFL:
   *--   1:  CONVERGED BY CRITERION 1
   *--   2:  CONVERGED BY NEGLIGIBLE FUNCTION CHANGE (CRITERION 3)
   *--   3:  REACHED MAXIMUM # ITERATIONS
   *
   *====================================================================
   */
    
  /*
   * Local Scalars
   */
  static double T1;
  static double athi,dbthi,den,dthi,ei8,fc,fm,fp,fy,si,thi,thib,
    thic,thim,thip,thmax,thmin,tpn,umtpn;
  static int i,ifix,ii,isti,lex,lret,n,ndecrm,ndecrp;

  /*
   * Local Arrays
   */
  static double thm[MAXFUN_NP],thp[MAXFUN_NP],thy[MAXFUN_NP];
  static int mex[MAXFUN_NP];

  /*
   * Equivalence
   */
  static double *ths = thm;

  /*
   * Format strings
   */
  const char* F9000 = "\nPARAMETER%3d IS TEMPORARILY FIXED AT BOUND\
%#17.8G\n";
  const char* F9100 = "\nPARAMETER%3d IS TEMPORARILY FIXED NEAR BOU\
ND%#17.8G\n";
  const char* F9200 = "\nWHEN PARAMETER%3d =%#17.8G\n   AND STEPSIZ\
E FACTOR =%#17.8G\n";
  const char* F9300 = "LOWER TRIAL VALUES ONLY WERE USED IN THIS IT\
ERATION BECAUSE UPPER\n   TRIAL VALUE MADE FUNCTION UNDEFINED\n";
  const char* F9400 = "UPPER TRIAL VALUES ONLY WERE USED IN THIS IT\
ERATION BECAUSE LOWER\n   TRIAL VALUE MADE FUNCTION UNDEFINED\n";
  const char* F9500 = "\nPARAMETER%3d IS TEMPORARILY FIXED\n";
  const char* F9600 = "AFTER THIS ITERATION,%3d IMPROVED VALUES FOU\
ND IN SAME DIRECTION.\n   NEW FUNCTION VALUE%#18.10E\n---------------\
-------------------------------------------------------\n";

  /*
   * Executable Statements
   */

  for(i=0; i<Exp.nt; i++) mex[i] = 0;

 S10:
  /*
   *--EACH ITERATION BEGINS HERE
   */
  for(i=0; i<Exp.nt; i++) Xthpr[i] = theta[i];
  Xfpr = *f;
    
  /*
   *--INITIALIZE FOR THIS ITERATION
   */
  Ximpbnd = 0;
  Xdifmax = 0.e0;
    
  /*
   *--LOOP THROUGH INDEPENDENT PARAMETERS THAT ARE STILL ALLOWED TO VARY
   */
  for(i=1; i<=Exp.nt; i++) {
    isti = Xist[i - 1];
    if(isti <= 2) {
            
      /*
       *--INITIALIZE WORKING VALUES FOR THIS PARAMETER
       */
      ifix = 0;
      fy = *f;
      thi = theta[i - 1];
      athi = ABS(thi);
      si = Xstp[i - 1];
      thmax = Exp.thu[i - 1];
      thmin = Exp.thl[i - 1];
            
      /*
       *--FIRST SELECT 3 APPROPRIATE TRIAL VALUES
       *
       *--CHECK FOR CLOSENESS TO A BOUND
       */
      bndchk(&thi,&thmax,&lex);
      if(lex <= 0) goto S20;

      if(si > efn(&athi,&Exp.epsd)) goto S50;

      thib = thmax;
      goto S30;

    S20:
      bndchk(&thi,&thmin,&lex);
      if(lex <= 0) goto S110;

      if(si > efn(&athi,&Exp.epsd)) goto S70;

      thib = thmin;
            
    S30:
      /*
       *--FIX THIS PARAMETER FOR THE REST OF BASIC ITERATION PROCESS
       */
      for(ii=0; ii<Exp.nt; ii++) thy[ii] = theta[ii];
      thy[i - 1] = thib;
      mfun->eval_fun(thy,&fy,nfe,&lex);
      if(lex > 0 || fy < *f) goto S40;

      Xist[i - 1] = isti + 4;
      if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,i,thib);
      if(Exp.idet != NULL && Exp.idet != Exp.iout)
	(void)fprintf(Exp.idet,F9000,i,thib);
      goto S320;
            
    S40:
      Xist[i - 1] = isti + 6;
      if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9100,i,thib);
      if(Exp.idet != NULL && Exp.idet != Exp.iout)
	(void)fprintf(Exp.idet,F9100,i,thib);
      goto S340;
            
    S50:
      /*
       *--CLOSE TO UPPER BOUND
       */
      lret = 3;
      thip = thmax;

      /*
       *--COME BACK TO THIS POINT (FOR LRET = 3) IF HAD TO DECREASE SI
       */
    S60:
      T1 = ABS(thmax);
      dbthi = dfn(&T1,&si);
      thim = thip - dbthi;
      goto S90;

    S70:
      /*
       *--CLOSE TO LOWER BOUND
       */
      lret = 2;
      thim = thmin;

      /*
       *--COME BACK TO THIS POINT (FOR LRET = 2) IF HAD TO DECREASE SI
       */
    S80:
      T1 = ABS(thmin);
      dbthi = dfn(&T1,&si);
      thip = thim + dbthi;

    S90:
      /*
       *--TWO CLOSE-TO-BOUND CASES JOIN AT THIS POINT
       */
      thi = 0.5e0 * (thip + thim);
      athi = ABS(thi);
      dthi = 0.5e0 * dbthi;
      si = dfninv(&athi,&dthi);
      ei8 = efn(&athi,&Exp.epsd) / 8.e0;

      for(ii=0; ii<Exp.nt; ii++) thy[ii] = theta[ii];
      thy[i - 1] = thi;
      mfun->eval_fun(thy,&fy,nfe,&lex);
      if(lex > 0) goto S100;

      Xdifmax = MAX(Xdifmax,ABS(fy - Xfpr));
      goto S160;

    S100:
      Ximpbnd = 1;
      si = 0.5e0 * si;
      if(si >= ei8) {
	switch(lret - 1) {
	case 1: goto S80;

	case 2: goto S60;

	default: break;
	}
      }
      ifix = 1;
      goto S330;

    S110:
      /*
       *--NOT CLOSE TO EITHER BOUND
       */
      lret = 1;
      ei8 = efn(&athi,&Exp.epsd) / 8.e0;
      if(mex[i - 1] > 0) goto S130;

      if(mex[i - 1] == 0) goto S140;

      si = 0.5e0 * si;
      if(si >= ei8) goto S140;

      /*
       *--SI HAS NOW BECOME SMALL IN NORMAL ITERATION PROCESS (NOT DUE TO
       *--PROBLEMS), SO GO AHEAD AND USE IT BEFORE FIXING PARAMETER
       */
      ifix = 1;
      goto S140;

    S130:
      si = MIN(si + si,0.5e0);

    S140:
      ndecrp = ndecrm = 0;

      /*
       *--COME BACK TO THIS POINT (FOR LRET = 1) IF HAD TO DECREASE SI
       */
    S150:
      dthi = dfn(&athi,&si);
      thip = thi + dthi;
      thim = thi - dthi;

    S160:
      /*
       *--CLOSE-TO-BOUND CASES JOIN NORMAL CASE AT THIS POINT
       *
       *--CHECK THAT ALL TRIAL VALUES ARE IN BOUNDS
       */
      if(thim < thmin || thip > thmax) goto S200;

      /*
       *--GET FUNCTION VALUES FOR UPPER AND LOWER TRIAL VALUES
       */
      for(ii=0; ii<Exp.nt; ii++) thp[ii] = theta[ii];
      thp[i - 1] = thip;
      mfun->eval_fun(thp,&fp,nfe,&lex);
      if(lex > 0) goto S170;

      Xdifmax = MAX(Xdifmax,ABS(fp - Xfpr));
      if(lret <= 1 && fp > fy) goto S220;

      goto S180;

    S170:
      Ximpbnd = 1;
      if(lret > 1) goto S200;

      ndecrp += 1;
      if(ndecrp <= 3) goto S200;

      thmax = thi;
      if(Exp.idet == NULL) goto S50;

      (void)fprintf(Exp.idet,F9200,i,thi,si);
      (void)fprintf(Exp.idet,"%s",F9300);
      goto S50;

    S180:
      for(ii=0; ii<Exp.nt; ii++) thm[ii] = theta[ii];
      thm[i - 1] = thim;
      mfun->eval_fun(thm,&fm,nfe,&lex);
      if(lex > 0) goto S190;

      Xdifmax = MAX(Xdifmax,ABS(fm - Xfpr));
      goto S210;

    S190:
      Ximpbnd = 1;
      if(lret > 1) goto S200;

      ndecrm += 1;
      if(ndecrm <= 3) goto S200;

      thmin = thi;
      if(Exp.idet == NULL) goto S70;

      (void)fprintf(Exp.idet,F9200,i,thi,si);
      (void)fprintf(Exp.idet,"%s",F9400);
      goto S70;

    S200:
      /*
       *--REDUCE SI AND TRY AGAIN
       */
      si = 0.5e0 * si;
      if(si >= ei8) {
	switch(lret) {
	case 1: goto S150;

	case 2: goto S80;

	case 3: goto S60;

	default: break;
	}
      }
      ifix = 1;
      goto S310;

    S210:
      /*
       *--HAVE THREE TRIAL VALUES AND CORRESPONDING FUNCTION VALUES;
       *--FIND BEST ESTIMATES
       */
      if(fy >= fp && fy >= fm) goto S250;

      /*
       *--EITHER THIP OR THIM IS BEST
       */
      if(fm > fp) goto S230;

      /*
       *--THIP IS BEST
       */
      if(fp >= *f) goto S220;

      mex[i - 1] = -1;
      goto S330;

    S220:
      mex[i - 1] = 1;
      *f = fp;
      goto S300;

    S230:
      /*
       *--THIM IS BEST
       */
      if(fm >= *f) goto S240;

      mex[i - 1] = -1;
      goto S330;

    S240:
      mex[i - 1] = 1;
      *f = fm;
      for(ii=0; ii<Exp.nt; ii++) theta[ii] = thm[ii];
      goto S330;

    S250:
      /*
       *--THI IS BEST (OR AS GOOD)
       */
      mex[i - 1] = -1;
      den = 2.e0 * (fm + fp - fy - fy);
      /* if(den == 0.e0) goto S310; */
      if(ABS(den) < EPSILON) goto S310;

      /*
       *--TRY PARABOLIC APPROXIMATION TO GET A BETTER ESTIMATE
       *
       *--(USE ARRAY THP ALREADY PARTLY PREPARED)
       */
      thic = thi - (fp - fm) * dthi / den;
      thp[i - 1] = thic;
      mfun->eval_fun(thp,&fc,nfe,&lex);
      if(lex <= 0) goto S260;

      Ximpbnd = 1;
      goto S280;

    S260:
      Xdifmax = MAX(Xdifmax,ABS(fc - Xfpr));
      if(fc - fy < 0.e0) goto S280;

      if(fc - fy > 0.e0) goto S290;

      /*
       *--THIC EQUALLY GOOD; FIX PARAMETER
       */
      ifix = 1;
      goto S310;

    S280:
      /*
       *--THI IS BETTER; RETRY WITH SMALLER STEPSIZE
       */
      dthi = MIN(0.5e0 * si,ABS(thi - thic));
      si = dfninv(&athi,&dthi);
      if(si < 8.e0 * ei8) goto S310;
            
      goto S150;

    S290:
      /*
       *--THIC IS BETTER
       */
      if(*f > fc) goto S330;

      *f = fc;

    S300:
      for(ii=0; ii<Exp.nt; ii++) theta[ii] = thp[ii];
      goto S330;

    S310:
      /*
       *--CHOOSE BETTER OF THI, THETA(I) (DIFFERENT ONLY IF LRET > 1)
       */
      if(lret == 1 || *f > fy) goto S330;

    S320:
      /*
       *--SET PARAMETER TO THI
       */
      *f = fy;
      for(ii=0; ii<Exp.nt; ii++) theta[ii] = thy[ii];

    S330:
      /*
       *  SAVE CURRENT STEPSIZE FACTOR
       */
      Xstp[i - 1] = si;

      /*
       *--FIX PARAMETER IF INDICATED
       */
      if(ifix > 0) {
	Xist[i - 1] = isti + 8;
	if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9500,i);
	if(Exp.idet != NULL && Exp.idet != Exp.iout)
	  (void)fprintf(Exp.idet,F9500,i);

	/*
	 *--END OF LOOP
	 */
      }
    }

  S340:;
        
  }

  /*
   *--FINISH UP THIS ITERATION
   */
  endit(theta,f);
  Xdifmax = MAX(Xdifmax,Xfch);

  /*
   *--WRITE OUT DETAILS OF THIS ITERATION IF DESIRED
   */
  if(Exp.idet != NULL) litd(theta,f,nfe);

  /*
   *--CHECK FOR CONVERGENCE
   */
  bcnvch(theta,&lex);
  switch(lex + 1) {
  case 1: goto S350;

  case 2: goto S380;

  case 3: goto S380;

  case 4: goto S420;

  default: break;
  }

 S350:
  /*
   *--DO "SUPPLEMENTARY ITERATION PROCESS" BEFORE DOING ANOTHER
   *--BASIC ITERATION
   *
   *--INITIALIZE
   */
  for(i=1; i<=Exp.nt; i++) ths[i - 1] = theta[i - 1];
  n = 1;
  tpn = 2.e0;

 S360:
  /*
   *--GO FARTHER IN SAME DIRECTION
   */
  umtpn = 1.e0 - tpn;
  for(i=1; i<=Exp.nt; i++) {
    if(Xist[i - 1] > 2) thy[i - 1] = ths[i - 1];
    else {
      thy[i - 1] = tpn * ths[i - 1] + umtpn * Xthpr[i - 1];
      if(thy[i - 1] < Exp.thl[i - 1] ||
	 thy[i - 1] > Exp.thu[i - 1])
	goto S370;

    }
  }
  mfun->eval_fun(thy,&fy,nfe,&lex);
  if(lex > 0 || fy <= *f) goto S370;

  /*
   *--THY CONTAINS IMPROVED ESTIMATES
   */
  for(i=1; i<=Exp.nt; i++) theta[i - 1] = thy[i - 1];
  *f = fy;

  n += 1;
  tpn += tpn;
  goto S360;

 S370:
  /*
   *--NO MORE IMPROVEMENT
   */
  if(n <= 1) goto S10;

  n -= 1;
  if(Exp.idet == NULL) goto S10;

  (void)fprintf(Exp.idet,F9600,n,*f);
  goto S10;

 S380:
  /*
   *--HAVE CONVERGED;
   *--BEFORE RETURNING, UNFIX ANY TEMPORARILY FIXED PARAMETERS
   */
  for(i=1; i<=Exp.nt; i++) {
    switch(Xist[i - 1]) {
    case  1: goto S410;

    case  2: goto S410;

    case  3: goto S410;

    case  4: goto S410;

    case  5: goto S390;

    case  6: goto S400;

    case  7: goto S390;

    case  8: goto S400;

    case  9: goto S390;

    case 10: goto S400;

    default: break;
    }

  S390:
    Xist[i - 1] = 1;
    goto S410;

  S400:
    Xist[i - 1] = 2;

  S410:;

  }

 S420:
  /*
   *--HAVE CONVERGED OR REACHED MAXIMUM # ITERATIONS; EXIT
   */
  *ifl = lex;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::nsrch(double theta[],double *f,int *nfe,double dth[],
			  int *ifl) {
  /*
 *====================================================================
 *
 *--TRY 2*(NI) "NEIGHBORING" PARAMETER SETS TO FIND A BETTER ONE
 *
 *--RETURNS FLAG IFL:
 *--   0:  SIGNIFICANT CHANGE; CONTINUE ITERATING
 *--   1:  NEGLIGIBLE OR NO CHANGE
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static double athi,fs,fy,si,thi,thiy;
  static int i,ii,lex;

  /*
   * Local Arrays
   */
  static double ths[MAXFUN_NP],thy[MAXFUN_NP];

  /*
 * Format Strings
 */
  const char* F9000 = "\nNO IMPROVEMENT IN 2(NI)-TRIAL SEARCH\n";
  const char* F9100 = "\nIMPROVED VALUES(S) FOUND IN 2(NI)-TRIAL SE\
ARCH\n";
  const char* F9200 = "\n2(NI)-TRIAL ITERATION:\n";
  const char* F9300 = "\nNEGLIGIBLE CHANGE IN FUNCTION VALUE:\
%#20.10E\n   IN 2(NI)-TRIAL SEARCH\n";

  /*
 * Executable Statements
 */

  /*
   *--SAVE INCOMING ESTIMATES AND ESTABLISH NEW STEPSIZES
   */
  for(i=0; i<Exp.nt; i++) {
    ths[i] = thy[i] = theta[i];
    athi = ABS(thy[i]);
    si = efn(&athi,&Exp.epsd);
    dth[i] = dfn(&athi,&si);
    Xstp[i] = si;
  }
  fs = *f;

  /*
   *--INITIALIZE BOUNDARY FLAG
   */
  Ximpbnd = 0;

  /*
 *--LOOP THROUGH PARAMETERS TO CHECK NEIGHBORING VALUES
 */
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) {
      thi = thy[i];
      thiy = thi + dth[i];
      if(thiy <= Exp.thu[i]) {
	thy[i] = thiy;
	mfun->eval_fun(thy,&fy,nfe,&lex);
	if(lex > 0) Ximpbnd = 1;
	else {
	  if(fy > *f) {
	    for(ii=0; ii<Exp.nt; ii++)
	      theta[ii] = thy[ii];
	    *f = fy;
	  }
	}
      }
      thiy = thi - dth[i];
      if(thiy >= Exp.thl[i]) {
	thy[i] = thiy;
	mfun->eval_fun(thy,&fy,nfe,&lex);
	if(lex > 0) Ximpbnd = 1;
	else {
	  if(fy > *f) {
	    for(ii=0; ii<Exp.nt; ii++)
	      theta[ii] = thy[ii];
	    *f = fy;
	  }
	}
      }
      thy[i] = theta[i];
    }
  }

  /*
 *--IF ANY IMPROVEMENT HAS OCCURRED, CONSIDER THIS ANOTHER ITERATION
 */
  if(*f <= fs) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9000);
    *ifl = 1;
    return;

  }
  for(i=0; i<Exp.nt; i++) Xthpr[i] = ths[i];
  Xfpr = fs;
  endit(theta,f);
  Xdifmax = Xfch;

  /*
 *--PRINT ITERATION DETAILS IF DESIRED
 */
  if(Exp.iout != NULL && Exp.iout != Exp.idet)
    (void)fprintf(Exp.iout,"%s",F9100);
  if(Exp.idet != NULL) {
    (void)fprintf(Exp.idet,"%s",F9200);
    litd(theta,f,nfe);
  }

  /*
 *--CHECK FOR CONVERGENCE
 */
  if(Xfch <= Exp.epsc1 * Exp.epsc1 || Xfch <= Exp.epsc3) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9300,Xfch);
    *ifl = 1;
    return;

  }
    
  /*
 *--SIGNIFICANT CHANGE
 */
  for(i=0; i<Exp.nt; i++) Xstp[i] = 2.e0 * Xstp[i];
  *ifl = 0;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::psrch(double theta[],double *f,int *nfe,double dth[],
			  int *ifl) {
  /*
 *====================================================================
 *
 *--SEARCH 2**(NI) NEIGHBORING PARAMETER SETS FOR BETTER ESTIMATES
 *
 *--RETURNS FLAG IFL:
 *--   0:  SIGNIFICANT CHANGE; CONTINUE ITERATING
 *--   1:  NEGLIGIBLE OR NO CHANGE
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static double fs,fy,thiy;
  static int i,k,lex,np1,nplus,nsw;

  /*
   * Local Arrays
   */
  static double ths[MAXFUN_NP],thy[MAXFUN_NP];
  static int list[(MAXFUN_NPV + 1)];

  /*
   * Format Strings
   */
  const char* F9000 = "\nNO IMPROVEMENT IN 2**(NI)-TRIAL SEARCH\n";
  const char* F9100 = "\nIMPROVED VALUE(S) FOUND IN 2**(NI)-TRIAL S\
EARCH\n";
  const char* F9200 = "\n2**(NI)-TRIAL ITERATION:\n";
  const char* F9300 = "\nNEGLIGIBLE CHANGE IN FUNCTION VALUE:\
%#20.10E IN 2**(NI)-TRIAL SEARCH\n";

  /*
 * Executable Statements
 */

  /*
   *--SAVE INCOMING VALUES
   */
  for(i=0; i<Exp.nt; i++) ths[i] = theta[i];
  fs = *f;

  /*
   *--INITIALIZE BOUNDARY FLAG
   */
  Ximpbnd = 0;

  /*
 *--CONSIDER TRIAL PARAMETER SETS WITH NPLUS INDEPENDENT PARAMETERS
 *--CHANGED IN + DIRECTION (AND OTHERS IN - DIRECTION),
 *--FOR NPLUS = 0 TO NI
 *
 *--NPLUS = 0 CASE:  NO PARAMETERS TO CHANGE IN + DIRECTION
 *
 */
  nplus = list[0] = 0;

 S10:
  /*
   *--SET UP TRIAL PARAMETER SET WITH NPLUS PARAMETERS CHANGING IN
   *-- + DIRECTION (WHICH ONES ARE INDICATED IN LIST ARRAY)
   */
  np1 = 1;
  k = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) goto S20;

    thy[i] = ths[i];
    goto S30;

  S20:
    k += 1;
    if(k == list[np1 - 1]) {
      thy[i] = ths[i] + dth[i];
      np1 += 1;
    }
    else thy[i] = ths[i] - dth[i];
    if(thy[i] > Exp.thu[i] || thy[i] < Exp.thl[i]) goto S40;

  S30:;
  }

  /*
 *--CHECK OUT THIS TRIAL SET
 */
  mfun->eval_fun(thy,&fy,nfe,&lex);
  if(lex > 0) Ximpbnd = 1;
  else {
    if(fy > *f) {
      /*
       *--HAVE AN IMPROVEMENT
       */
      for(i=0; i<Exp.nt; i++) theta[i] = thy[i];
      *f = fy;
    }
  }
 S40:
  /*
   *--FINISHED WITH THIS TRIAL SET; GET NEXT ONE
   */
  if(nplus > 0) goto S60;

 S50:
  /*
   *--DONE WITH THIS NPLUS; INCREASE # PARAMETERS TO CHANGE IN
   *-- + DIRECTION
   */
  nplus += 1;
  if(nplus > Xni) goto S70;

  nsw = 0;

 S60:
  /*
   *--GET LIST OF NPLUS PARAMETERS TO CHANGE IN + DIRECTION
   */
  comb(&Xni,&nplus,&nsw,list);
  if(nsw > 1) goto S50;

  goto S10;

 S70:
  /*
   *--HAVE TRIED EVERY DIRECTION OF CHANGE; ANY IMPROVEMENT?
 */
  if(*f <= fs) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9000);
    *ifl = 1;
    return;

  }

  /*
 *--INCREMENT # TIMES 2**(NI)-TRIAL SEARCH MADE IMPROVEMENT
 */
  Xnsurf2 += 1;

  /*
 *--TRY GOING FARTHER IN BEST DIRECTION
 */
  for(i=0; i<Exp.nt; i++) {
    thiy = 2.e0 * theta[i] - ths[i];
    if(Xist[i] > 2 || thiy < Exp.thl[i] || thiy > Exp.thu[i])
      thy[i] = theta[i];
    else thy[i] = thiy;
  }
  mfun->eval_fun(thy,&fy,nfe,&lex);
  if(lex > 0) Ximpbnd = 1;
  else {
    if(fy > *f) {
      for(i=0; i<Exp.nt; i++) theta[i] = thy[i];
      *f = fy;
    }
  }

  /*
 *--COMPLETE THIS "ITERATION"
 */
  for(i=0; i<Exp.nt; i++) Xthpr[i] = ths[i];
  Xfpr = fs;
  endit(theta,f);
  Xdifmax = Xfch;

  /*
   *--PRINT ITERATION DETAILS IF DESIRED
   */
  if(Exp.iout != NULL && Exp.iout != Exp.idet)
    (void)fprintf(Exp.iout,"%s",F9100);
  if(Exp.idet != NULL) {
    (void)fprintf(Exp.idet,"%s",F9200);
    litd(theta,f,nfe);
  }

  /*
 *--CHECK FOR CONVERGENCE
 */
  if(Xfch <= Exp.epsc1 * Exp.epsc1 || Xfch <= Exp.epsc3) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9300,Xfch);
    *ifl = 1;
    return;

  }

  /*
 *--SIGNIFICANT CHANGE
 */
  for(i=0; i<Exp.nt; i++) Xstp[i] = Exp.epsd;
  *ifl = 0;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::comb(int *ntotal,int *nchoos,int *nsw,int list[]) {
  /*
   *====================================================================
   *
   *--CHOOSE NCHOOS ELEMENTS OUT OF NTOTAL; INDICATE WHICH ONES IN LIST
   *--(GENERATES COMBINATIONS TO BE USED IN THE 2**(NI)-TRIAL SEARCH:
   *--INDICATE WHICH PARAMETERS TO CHANGE IN + DIRECTION)
   *
   *====================================================================
   */

  /*
   * Local Scalars and Arrays
   */
  static int k,kv,listkv,notch,maxlst[MAXFUN_NPV];

  /*
   * Executable Statements
   */
  if(*nsw > 0) goto S10;

  /*
   *--START WITH A NEW NCHOOS
   */
  notch = *ntotal - *nchoos;
  for(k=1; k<=*nchoos; k++) maxlst[k - 1] = notch + k;
  *nsw = 1;
  list[*nchoos] = 0;
  kv = list[0] = 1;
  goto S30;

 S10:
  /*
   *--IN PROGRESS WITH THIS NCHOOS
   */
  listkv = list[kv - 1] + 1;
  list[kv - 1] = listkv;
  if(listkv <= maxlst[kv - 1]) return;

 S20:
  kv -= 1;
  if(kv <= 0) {
    *nsw = 2;
    return;

  }
  list[kv - 1] += 1;
  if(list[kv - 1] > maxlst[kv - 1]) goto S20;

 S30:
  /*
   *--HAVE ESTABLISHED KV'TH CHOSEN ELEMENT; FILL IN REST OF LIST WITH
   *--LOWEST POSSIBLE VALUES
   */
  if(kv < *nchoos) {
    for(k=kv + 1; k<=*nchoos; k++) list[k - 1] = list[k - 2] + 1;
  }
  kv = *nchoos;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::nrstep(double theta[],double *f,int *nfe,int *ifl) {
  /*
   *====================================================================
   *
   *--PERFORM ONE NEWTON-RAPHSON ITERATION
   *
   *--RETURNS FLAG IFL:
   *--   0:  CONTINUE ITERATING
   *--   1:  CONVERGED BY CRITERION 1 (IMPLIED TEST)
   *--   2:  CONVERGED BY CRITERION 2:  PTG <= EPSC2
   *--   3:  CONVERGED BY NEGLIGIBLE FUNCTION CHANGE (CRITERION 3)
   *--   4:  REACHED MAXIMUM # ITERATIONS
   *
   *====================================================================
   */

  /*
   * Parameter
   */
  static int K3 = 1;

  /*
   * Local Scalars
   */
  static int T1,T2,T4;
  static double fs,pl,tt;
  static int i,l,lex;

  /*
   * Local Arrays
   */
  static double ths[MAXFUN_NP];

  /*
   * Format Strings
   */
  const char* F9000 = "\nCONVERGED BY CRITERION 2 AFTER ITERATION\
%4d:\n   NORMALIZED GRADIENT%#18.10E WITHIN SPECIFIED TOLERANCE\n";
  const char* F9100 = "\nAFTER ITERATION%4d,\n   UNDEFINED FUNCTIO\
N WITH STEP SIZE%#18.10E\n";
  const char* F9200 = "\nAFTER ITERATION%4d,\n   DECREASING FUNCTIO\
N WITH STEP SIZE%#18.10E\n";

  /*
   * Executable Statements
   */

  /*
   *--SAVE INCOMING VALUES
   */
  for(i=0; i<Exp.nt; i++) ths[i] = theta[i];
  fs = *f;

  /*
   *--COMPUTE DIRECTION OF CHANGE
   */
  T1 = MAXFUN_NPV;
  T2 = MAXFUN_NPV;
  T4 = MAXFUN_NPV;
  /*   mmult(Xv,&T1,&Xnv,&Xnv,Xg,&T2,&K3,Xpdir,&T4); */
  mmult(&Xv[0][0],&T1,&Xnv,&Xnv,&Xg[0],&T2,&K3,&Xpdir[0],&T4);

  /*
   *--REDUCE AMOUNT OF CHANGE IF NECESSARY TO PRESERVE BOUNDS ON
   *--INDEPENDENT PARAMETERS
   */
  Xptg = 0.e0;
  tt = 1.e0;
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] > 2) goto S40;

    l += 1;
    pl = Xpdir[l - 1];
    /* if(pl == 0.e0) goto S40; */
    if(ABS(pl) < EPSILON) goto S40;

    else if(pl > 0.e0) {
      tt = MIN(tt,(Exp.thu[i] - ths[i]) / pl);

    S10:
      if(ths[i] + tt * pl <= Exp.thu[i]) goto S30;

      tt = (1.e0 - Exp.yota) * tt;
      goto S10;

    }
    else {
      tt = MIN(tt,(Exp.thl[i] - ths[i]) / pl);

    S20:
      if(ths[i] + tt * pl >= Exp.thl[i]) goto S30;

      tt = (1.e0 - Exp.yota) * tt;
      goto S20;
            
    }

  S30:
    Xptg += (pl * Xg[l - 1]);

  S40:;

  }

  /*
   *--SEE IF PTG SMALL ENOUGH TO STOP
   */
  if(ABS(Xptg) <= Exp.epsc2) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,Xit,Xptg);
    *ifl = 2;
    return;

  }

  /*
   *--INITIALIZE BOUNDARY FLAG
   */
  Ximpbnd = 0;

 S50:
  /*
   *--COMPUTE NEW ESTIMATES
   */
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) {
      l += 1;
      theta[i] = ths[i] + tt * Xpdir[l - 1];
    }
  }

  /*
   *--GET FUNCTIONAL VALUE IF POSSIBLE
   */
  mfun->eval_fun(theta,f,nfe,&lex);
  if(lex <= 0) goto S60;

  if(Exp.idet != NULL) (void)fprintf(Exp.idet,F9100,Xit,tt);
  Ximpbnd = 1;
  goto S70;

 S60:
  /*
   *--CHECK FOR DECREASING VALUE
   */
  if(*f >= fs) goto S80;

  if(Exp.idet != NULL) (void)fprintf(Exp.idet,F9200,Xit,tt);

 S70:
  /*
   *--RAN INTO TROUBLE; CUT DOWN STEP SIZE AND TRY AGAIN
   */
  tt = 0.5e0 * tt;
  goto S50;
    
 S80:
  /*
   *--OK TO COMPLETE THIS ITERATION
   */
  Xtstep = tt;
  for(i=0; i<Exp.nt; i++) Xthpr[i] = ths[i];
  Xfpr = fs;
  endit(theta,f);

  /*
   *--INCREMENT NUMBER OF ITERATIONS WITH THIS VARIANCE-COVARIANCE
   *--MATRIX AND THIS GRADIENT VECTOR
   */
  Xivage += 1;
  Xigage += 1;

  /*
   *--WRITE OUT DETAILS OF THIS ITERATION IF DESIRED
   */
  if(Exp.idet != NULL) litv(theta,f,nfe);

  /*
   *--CHECK FOR CONVERGENCE
   */
  vcnvch(theta,&lex);

  /*
   *--SET FLAG AND RETURN
   */
  if(lex <= 1) *ifl = lex;
  else *ifl = lex + 1;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::binit(double theta[],double *f,int *nfe,int *ihess) {
  /*
   *====================================================================
   *
   *--OBTAIN INITIAL B MATRIX (APPROXIMATING OPPOSITE OF HESSIAN AT
   *--OPTIMAL ESTIMATES) IN FACTORIZED FORM B = ULT*DIAG*ULT'
   *--(NOTE:  OPPOSITES ARE NECESSARY SINCE UPDATING PROCEDURE REQUIRES
   *--POSITIVE DEFINITE B)
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static double bl,ultl2l;
  static int l,l1,l2,lex;

  /*
   * Local Arrays
   */
  static double b[MAXFUN_NPV][MAXFUN_NPV];

  /*
   * Format Strings
   */
  const char* F9000 = "\nCOMPUTING INITIAL B FROM HESSIAN AT INITIA\
L ESTIMATES\n";
  const char* F9100 = "FACTORIZATION OF B MATRIX FAILED\n";
  const char* F9200 = "\nUSING IDENTITY FOR INITIAL B MATRIX\n";
  const char* F9300 = "--------------------------------------------\
--------------------------\n";

  /*
   * Executable Statements
   */

  /*
   *--SELECT TYPE OF INITIAL MATRIX TO BE USED
   */
  if(*ihess <= 0) goto S20;

  /*
   *--COMPUTE MATRIX OF SECOND PARTIAL DERIVATIVES FOR INITIAL ESTIMATES
   *--OF PARAMETERS
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9000);
  deriv2(theta,f,nfe,&Exp.ihit,&lex);
  if(lex > 1) goto S20;

  /*
   *--INITIALIZE B TO -H
   */
  for(l1=0; l1<Xnv; l1++) {
    for(l2=0; l2<Xnv; l2++) b[l2][l1] = -Xh[l2][l1];
  }

  /*
   *--DO CHOLESKY FACTORIZATION
   */
  for(l=1; l<=Xnv; l++) {
    bl = b[l - 1][l - 1];
    if(bl > 0.e0) goto S10;

    if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9100);
    goto S20;

  S10:
    Xdiag[l - 1] = bl;
    Xult[l - 1][l - 1] = 1.e0;
    if(l >= Xnv) goto S30;

    for(l1=l + 0; l1<Xnv; l1++)
      Xult[l - 1][l1] = b[l - 1][l1] / bl;
    for(l2=l + 1; l2<=Xnv; l2++) {
      ultl2l = Xult[l - 1][l2 - 1];
      for(l1=l2 - 1; l1<Xnv; l1++)
	b[l2 - 1][l1] -= (ultl2l * Xult[l - 1][l1] * bl);
    }
  }
  goto S30;

 S20:
  /*
   *--USE IDENTITY MATRIX
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9200);
  for(l1=0; l1<Xnv; l1++) {
    for(l2=1; l2<=l1; l2++) Xult[l2 - 1][l1] = 0.e0;
    Xult[l1][l1] = Xdiag[l1] = 1.e0;
  }

 S30:
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,"%s",F9300);
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::bupdt(double gpr[],int *lex) {
  /*
   *====================================================================
   *
   *--UPDATE MATRIX B = ULT*DIAG*ULT'
   *
   *--RETURNS FLAG LEX:
   *--   0:  UPDATE SUCCESSFUL
   *--   1:  UPDATE FAILED:  LOST ALL SIGNIFICANT DIGITS IN OPTIMAL
   *--       CONDITIONING
   *
   *--USING PREVIOUS GRADIENT (BEFORE ITERATION IT) IN GPR AND CURRENT
   *--GRADIENT (AFTER ITERATION IT) IN G
   *
   *--------------------------------------------------------------------
   *- NOTE THAT G, GPR, Y ARE DEFINED OPPOSITE TO QUANTITIES DESCRIBED -
   *- IN GEMINI DOCUMENTATION, AND FORMULAS ARE CHANGED ACCORDINGLY    -
   *--------------------------------------------------------------------
   *- NOTE ALSO THAT VARIABLES TAU1, TAU2, AND PHL REFER TO QUANTITIES -
   *- THET1, THET2, AND THL IN GEMINI; THESE DO NOT REFER TO THE       -
   *- PARAMETER ESTIMATE VECTOR THETA                                  -
   *--------------------------------------------------------------------
   *
   *====================================================================
   */
    
  /*
   * Local Scalars
   */
  static double aa,al,alph1,alph2,alpha,bb,bettl,bl,cc,del1,del2,
    eta,fbcd,gamtl,lambl,lambl2,mul,nu,phl,rdiv,sa,sb,sbc,sc,sd,
    sdb,sdbc,sden,tau1,tau2,ytp;
  static int iswup,l,ll;

  /*
   * Local Arrays
   */
  static double dp[MAXFUN_NPV],s[MAXFUN_NPV],w[MAXFUN_NPV],
    wpl[MAXFUN_NPV],wplp1[MAXFUN_NPV],wtil[MAXFUN_NPV],
    wtplp1[MAXFUN_NPV],y[MAXFUN_NPV],z[MAXFUN_NPV],
    zpl[MAXFUN_NPV],zplp1[MAXFUN_NPV],ztil[MAXFUN_NPV],
    ztplp1[MAXFUN_NPV];

  /*
   * Format String
   */
  const char* F9000 = "\nSTOPPED AFTER ITERATION%5d BECAUSE ALL SIG\
NIFICANT DIGITS LOST\n   THROUGH CANCELLATION IN OPTIMAL CONDITIONING\
\n";
    
  /*
   * Executable Statements
   */

  /*
   *--CHECK THAT DIRECTION OF CHANGE IS ONE OF DECREASING GRADIENT;
   *--IF NOT, DON'T UPDATE
   */
  ytp = 0.e0;
  for(l=0; l<Xnv; l++) {
    y[l] = Xg[l] - gpr[l];
    ytp += (y[l] * Xpdir[l]);
  }
  if(ytp >= 0.e0) goto S60;

  if(Xnv > 1) goto S10;

  /*
   *--FOR NV = 1:
   */
  Xdiag[0] = -(y[0] * Xdiag[0] / (Xtstep * gpr[0]));
  goto S60;

 S10:
  /*
   *--FOR NV > 1 (FORMULA NUMBERS REFER TO GEMINI DOCUMENTATION):
   *
   *--------------------------------------------------------------------
   *
   *--GOLDFARB METHOD FOR VARIABLE METRIC B-MATRIX UPDATING;
   *--UPDATE OF FACTORIZED B-MATRIX BY RANK TWO MODIFICATION
   *--IN REAL PRODUCT FORM WITH FORMULA (3.25) OR (3.28) OF GEMINI
   *--DOCUMENTATION, DEPENDING ON WHETHER NONE OR AT LEAST ONE
   *--LAMBA(L)**2 IS GREATER THAN 10.
   *--OPTIMAL CONDITIONING OF DAVIDON ALSO INCORPORATED.
   *
   *--------------------------------------------------------------------
   *--SOLVE (3.20) AND (3.21) FOR WTIL AND ZTIL
   */
  wtil[0] = gpr[0];
  ztil[0] = -y[0];
  for(l=1; l<Xnv; l++) {
    wtil[l] = gpr[l];
    ztil[l] = -y[l];
    for(ll=0; ll<l; ll++) {
      wtil[l] -= (Xult[ll][l] * wtil[ll]);
      ztil[l] -= (Xult[ll][l] * ztil[ll]);
    }
  }

  /*
   *--------------------------------------------------------------------
   *--GET SCALARS B, C AND D, DENOTED SB, SC, SD HERE (3.27)
   */
  sb = sc = sd = 0.e0;
  for(l=0; l<Xnv; l++) {
    sb += (Xtstep * ztil[l] * wtil[l] / Xdiag[l]);
    sc += (Xtstep * Xtstep * wtil[l] * wtil[l] / Xdiag[l]);
    sd += (ztil[l] * ztil[l] / Xdiag[l]);
  }

  /*
   *--------------------------------------------------------------------
   *--OPTIMAL CONDITIONING OF UPDATE, ACCORDING TO THEOREM 3 OF DAVIDON
   *--(1975) -- MAY BE DISABLED BY "goto S60;"
   *---------------------
   *
   *--PREVENTING UNDERFLOW:
   *--IF ANY OF SB, SC, SD IS SMALLER THAN 1.D-10, USE SR1-UPDATE FOR B
   */
  fbcd = MIN(MIN(sb,sc),sd);
  if(fbcd <= 1.e-10) goto S20;

  fbcd = 2.e0 * sc * (sd / sb) / (sc + sd);
  if(fbcd < 1.e0) goto S20;

  /*
   *---------------------
   *--RANK TWO UPDATE
   *
   *--GET ALPHA
   *
   */
  aa = sb / sc - 2.e0 * (sd / sb) + sd / sc;
  bb = sb / sc - 1.e0;
  cc = 1.e0 - sb / sd;
  del2 = bb * bb - aa * cc;

  /*
   *--IF DISCRIMINANT NEGATIVE OR EQUALS ZERO, TAKE ALPHA EQUAL TO ZERO
   */
  if(del2 <= 1.e-8) goto S30;

  /*
   *----------------
   */
  del1 = sqrt(del2);
  alph1 = (-bb + del1) / aa;
  alph2 = (-bb - del1) / aa;

  /*
   *--FOR NOW, ALWAYS CHOOSE SOLUTION OF SMALLEST MODULUS
   */
  if(ABS(alph1) <= ABS(alph2)) alpha = alph1;
  else alpha = alph2;

  /*
   *--IF ALPHA VERY SMALL, ALPHA TAKEN EQUAL TO ZERO
   */
  if(ABS(alpha) <= 1.e-5) goto S30;

  /*
   *--GET SA
   */
  sa = (alpha + 1.e0) * (alpha + 1.e0) + sc / sb - alpha * alpha *
    (sc / sb) * (sd / sb) - 1.e0 + alpha * alpha * (sd / sb);
  if(sa <= 0.e0) sa = 0.e0;
  else sa = 1.e0 / (sqrt(sa) * sb);

  /*
   *--GET TAU1 AND TAU FOR NON-TRIVIAL ALPHA
   */
  rdiv = 1.e0 / (alpha * alpha * sd + 2.e0 * alpha * sb + sc);
  tau1 = -((sa * (alpha * (alpha * sd + sb)) + 1.e0) * rdiv);
  tau2 = sa + (alpha * sa * (sc + alpha * sb) - alpha) * rdiv;
  goto S40;

 S20:
  /*
   *---------------------
   *--SR1 (SYMMETRIC RANK-ONE) UPDATE
   */
  alpha = -1.e0;
  sbc = sb - sc;
  /*   if(sbc == 0.e0) goto S80; */
  if(ABS(sbc) < EPSILON) goto S80;

  sdb = sd - sb;
  /*   if(sdb == 0.e0) goto S80; */
  if(ABS(sdb) < EPSILON) goto S80;

  sden = sqrt(ABS(sbc)) * sqrt(ABS(sdb));
  /*   if(sden == 0.e0) goto S80; */
  if(ABS(sden) < EPSILON) goto S80;

  sa = 1.e0 / sden;
  sdbc = sd - 2.e0 * sb + sc;
  /*   if(sdbc == 0.e0) goto S80; */
  if(ABS(sdbc) < EPSILON) goto S80;

  tau1 = -((sdb * sa + 1.e0) / sdbc);
  tau2 = sa + (sa * sbc + 1.e0) / sdbc;
  goto S40;

 S30:
  /*
   *---------------------
   *--ALPHA ZERO:  DFP-UPDATE FOR B
   *
   */
  alpha = 0.e0;
  sa = 1.e0 / (sqrt(sb) * sqrt(sc));
  tau1 = -(1.e0 / sc);
  tau2 = sa;

 S40:
  /*
   *--------------------------------------------------------------------
   *--START WITH SWITCH SET TO DO (3.25) TYPE UPDATE
   */
  iswup = 1;

  /*
   *--SAVE ULT LOWER TRIANGLE IN UPPER TRIANGLE OF ARRAY ULT, IN CASE A
   *--SWITCH TO (3.28) SHOULD BE REQUIRED
   */
  for(l=1; l<Xnv; l++) {
    for(ll=0; ll<l; ll++) Xult[l][ll] = Xult[ll][l];
  }
    
  /*
   *--GET W AND Z AS GIVEN IN (3.22), (3.23)
   */
  for(l=0; l<Xnv; l++) {
    w[l] = Xtstep * wtil[l] + alpha * ztil[l];
    z[l] = Xtstep * tau1 * wtil[l] + tau2 * ztil[l];
  }

  /*
   *--------------------------------------------------------------------
   *--READY TO APPLY GOLDFARB RECURRENCE 3 TO SOLVE CONCURRENTLY (3.25),
   *--(3.24) BEING ALSO SOLVED SIMULTANEOUSLY
   *--(OR (3.28) AND CORRESPONDING ANALOG OF (3.24))
   *
   *--GET S FIRST (P.802 (TOP) OF ?)
   */
  s[Xnv - 1] = 0.e0;
  for(l=Xnv - 1; l>=1; l--)
    s[l - 1] = s[l] + w[l] * w[l] / Xdiag[l];

 S50:
  /*
   *--------------------------------------------------------------------
   *--INITIALIZE NU, ETA
   */
  nu = 1.e0;
  eta = 0.e0;

  /*
   *--------------------------------------------------------------------
   *--INITIALIZE WTPLP1, ZTPLP1 OR WPLP1, ZPLP1
   */
  if(iswup <= 1) {
    for(l=0; l<Xnv; l++) {
      wtplp1[l] = gpr[l];
      ztplp1[l] = -y[l];
    }
  }
  else {
    for(l=0; l<Xnv; l++) {
      wplp1[l] = Xtstep * Xg[l] - alpha * y[l];
      zplp1[l] = Xtstep * tau1 * Xg[l] - tau2 * y[l];
    }

    /*
     *--------------------------------------------------------------------
     */
  }
  for(l=0; l<Xnv - 1; l++) {

    /*
     *--RECURRENCE FORMULA FROM (3.24) OR ANALOG
     */
    if(iswup <= 1) {
      for(ll=l + 1; ll<Xnv; ll++) {
	wtplp1[ll] -= (wtil[l] * Xult[l][ll]);
	ztplp1[ll] -= (ztil[l] * Xult[l][ll]);
      }
    }
    else {
      for(ll=l + 1; ll<Xnv; ll++) {
	wpl[ll] = wplp1[ll];
	zpl[ll] = zplp1[ll];
	wplp1[ll] = wpl[ll] - w[l] * Xult[l][ll];
	zplp1[ll] = zpl[ll] - z[l] * Xult[l][ll];
      }

      /*
       *--RECURRENCE 3 (?) TO GET AL, BL, ETC...
       */
    }
    al = nu * z[l] - eta * w[l];
    phl = 1.e0 + al * w[l] / Xdiag[l];
    lambl2 = phl * phl + al * al * s[l] / Xdiag[l];

    /*
     *--SWITCH TO (3.28) UPDATE IF ANY LAMBL2 GREATER THAN 4
     */
    if(iswup <= 1 && lambl2 > 4.e0) goto S70;

    /*
     *--COMPUTE L'TH ELEMENT OF D+ (NEW DIAG)
     */
    dp[l] = Xdiag[l] * lambl2;

    /*
     *--TAKES SIGN OF LAMBDA OPPOSITE TO THAT OF THETA
     */
    lambl = sqrt(lambl2);
    if(phl > 0.e0) lambl = -lambl;
        
    mul = phl - lambl;
    bl = phl * w[l] + al * s[l];
        
    /*
     *--NOTE:  GAMTL AND BETTL STAND FOR GAMMA TILDE AND BETA TILDE
     *--RESPECTIVELY IN (3.26)
     */
    gamtl = bl * nu / (lambl2 * Xdiag[l]);
    bettl = (al - bl * eta) / (lambl2 * Xdiag[l]);

    /*
     *--UPDATE L'TH COLUMN OF ULT
     */
    if(iswup <= 1) {
      for(ll=l + 1; ll<Xnv; ll++)
	Xult[l][ll] += (Xtstep * (bettl + tau1 * gamtl) *
			wtplp1[ll] + (alpha * bettl + tau2 *
				      gamtl) * ztplp1[ll]);
    }
    else {
      for(ll=l + 1; ll<Xnv; ll++)
	Xult[l][ll] = Xult[l][ll] /
	  lambl2 + bettl * wpl[ll] + gamtl * zpl[ll];

      /*
       *--UPDATE NU, ETA
       */
    }
    nu = -(nu / lambl);
    eta = -((eta + al * al / (mul * Xdiag[l])) / lambl);
  }
    
  /*
   *--GET LAST LAMBDA FOR D+
   */
  al = nu * z[Xnv - 1] - eta * w[Xnv - 1];
  lambl = 1.e0 + al * w[Xnv - 1] / Xdiag[Xnv - 1];

  dp[Xnv - 1] = Xdiag[Xnv - 1] * lambl * lambl;

  /*
   *--UPDATE DIAG FROM D+
   */
  for(l=0; l<Xnv; l++) Xdiag[l] = dp[l];

 S60:
  /*
   *--NORMAL EXIT
   */
  *lex = 0;
  return;

 S70:
  /*
   *--SWITCH TO UPDATING PROCEDURE 2
   */
  iswup = 2;
  for(l=1; l<Xnv; l++) {
    for(ll=0; ll<l; ll++) Xult[ll][l] = Xult[l][ll];
  }
  goto S50;

 S80:
  /*
   *--ERROR EXIT IN CASE IF DIVISION BY ZERO WOULD OCCUR IN OPTIMAL
   *--CONDITIONING
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,Xit);
  *lex = 1;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::direct(int *lex) {
  /*
   *====================================================================
   *
   *--COMPUTE DIRECTION OF SEARCH (SOLVING (-ULT*DIAG*ULT')*PDIR = -G
   *--FOR PDIR) AND NORMALIZED GRADIENT
   *
   *--RETURNS FLAG LEX:
   *--   0:  PDIR SUCCESSFULLY DETERMINED; CONTINUE
   *--   1:  CONVERGED BY CRITERION 2:  PTG <= EPSC2
   *--   2:  PDIR NOT IN ASCENT DIRECTION
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static int l,ll;

  /*
   * Local Array
   */
  static double q[MAXFUN_NPV];

  /*
   * Format String
   */
  const char* F9000 = "\nCONVERGED BY CRITERION 2 AFTER ITERATION\
%4d:\n   NORMALIZED GRADIENT%#18.10E WITHIN SPECIFIED TOLERANCE\n";

  /*
   * Executable Statements
   */

  /*
   *--SEPARATE PATHS FOR ONE OR MORE PARAMETERS
   */
  if(Xnv <= 1) {
    /*
     *--FOR NV = 1:
     */
    Xpdir[0] = Xg[0] / Xdiag[0];
  }
  else {
    /*
     *--FOR NV > 1:
     *
     *--FIRST SOLVE   ULT*Q = G   FOR Q
     *--(USE G RATHER THAN -G TO COMPENSATE FOR THE FACT THAT
     *--ULT*DIAG*ULT' CORRESPONDS TO OPPOSITE OF ACTUAL HESSIAN MATRIX)
     */
    q[0] = Xg[0];
    for(l=1; l<Xnv; l++) {
      q[l] = Xg[l];
      for(ll=0; ll<l; ll++) q[l] -= (Xult[ll][l] * q[ll]);
    }

    /*
     *--THEN SOLVE   ULT'*PDIR = DIAG**(-1)*Q   FOR PDIR
     */
    Xpdir[Xnv - 1] = q[Xnv - 1] / Xdiag[Xnv - 1];
    for(l=Xnv - 1 - 1; l>=0; l--) {
      Xpdir[l] = q[l] / Xdiag[l];
      for(ll=l + 1; ll<Xnv; ll++) Xpdir[l] -= (Xult[l][ll] *
					       Xpdir[ll]);
    }
  }

  /*
   *--INITIALIZE EXIT FLAG
   */
  *lex = 0;

  /*
   *--COMPUTE NORMALIZED GRADIENT PDIR'G
   */
  Xptg = 0.e0;
  for(l=0; l<Xnv; l++) Xptg += (Xpdir[l] * Xg[l]);

  /*
   *--MAKE SURE PTG > 0 (DIRECTION WILL INCREASE FUNCTION)
   */
  if(Xptg <= 0.e0) {
    *lex = 2;
    return;

  }

  /*
   *--SEE IF PTG SMALL ENOUGH TO STOP
   */
  if(Xptg > Exp.epsc2) return;

  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,Xit,Xptg);
  *lex = 1;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::lsrch(double theta[],double *f,int *nfe,int *ifirst,
			  int *ifl) {
  /*
 *====================================================================
 *
 *--LINE SEARCH IN CHOSEN DIRECTION (PDIR)
 *
 *--RETURN FLAG IFL:
 *--  -1:  ITERATION COMPLETE, BUT TSTEP < TIN AND < TMIN
 *--   0:  ITERATION COMPLETE, WITH TSTEP >= TIN AND/OR TMIN
 *--   1:  ITERATION COMPLETE; CONVERGED BY IMPLIED TEST
 *--   2:  ITERATION COMPLETE; NEGLIGIBLE FUNCTION CHANGE
 *--   3:  ITERATION COMPLETE; REACHED MAXIMUM # ITERATIONS
 *--   4:  COULD NOT COMPLETE ITERATION
 *
 *====================================================================
 */

  /*
   * Define Local Preprocessor Macro
   */
#define CURV .75e0

  /*
 * Local Scalars
 */
  static double T1;
  static double diff1,diff2,fs,ft,ft2,scal,tbnd,thpm,tin,tmin,tt,
    tt2;
  static int i,l,lex;

  /*
   * Local Arrays
   */
  static double ths[MAXFUN_NP],tht[MAXFUN_NP];

  /*
 * Format Strings
 */
  const char* F9000 = "   INITIAL STEP SIZE%#13.5E\n   MINIMUM STE\
P SIZE%#13.5E  MAXIMUM STEP SIZE%#13.5E\n";
  const char* F9100 = "\nACTIVE BOUNDARY CONSTRAINT ON INDEPENDENT \
PARAMETERS\n";

  /*
 * Executable Statements
 */

  /*
   *--SET TIN (INITIAL TT), TMIN (SIMPLE INITIAL VALUES IF THIS IS FIRST
   *--ITERATION WITH THIS NV)
   */
  if(*ifirst > 0) {
    /*
     *--FIRST TIME
     */
    tmin = 0.e0;
    tin = .1e0;
  }
  else {
    /*
     *--NOT FIRST TIME
     */
    thpm = 1.e30;
    l = 0;
    for(i=0; i<Exp.nt; i++) {
      if(Xist[i] <= 2) {
	l += 1;
	/* if(Xpdir[l - 1] != 0.e0) { */
	if(ABS(Xpdir[l - 1]) > EPSILON) {
	  T1 = ABS(theta[i]);
	  thpm = MIN(thpm,dfn(&T1,&Exp.yota) / 
		     ABS(Xpdir[l - 1]));
	}
      }
    }
    tmin = .5e0 * thpm / Exp.epst;
        
    tin = 2.e0 * Xfch / Xptg;
    if(tin <= 0.e0) tin = 1.e0;
    tin = MIN(tin,1.e0);
    if(Xidif > 1) tin = MAX(tin,1.e0);
  }

  /*
 *--ESTABLISH MAXIMUM TSTEP (TBND) TO SATISFY BOUNDARY RESTRICTIONS ON
 *--INDEPENDENT PARAMETERS
 */
  tbnd = 1.e5;
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) {
      l += 1;
      if(Xpdir[l - 1] > 0.e0) {
	tbnd = MIN(tbnd,
		   (Exp.thu[i] - theta[i]) / Xpdir[l - 1]);

      S10:
	if(theta[i] + tbnd * Xpdir[l - 1] <= Exp.thu[i])
	  goto S30;

	tbnd = (1.e0 - Exp.yota) * tbnd;
	goto S10;

      }
      else if(Xpdir[l - 1] < 0.e0) {
	tbnd = MIN(tbnd,
		   (Exp.thl[i] - theta[i]) / Xpdir[l - 1]);

      S20:
	if(theta[i] + tbnd * Xpdir[l - 1] >= Exp.thl[i])
	  goto S30;

	tbnd = (1.e0 - Exp.yota) * tbnd;
	goto S20;

      }
    }

  S30:;

  }

  /*
 *--MODIFY TIN IF NECESSARY
 */
  if(tin * (2.e0 + Exp.yota) >= tbnd)
    tin = tbnd * (.5e0 - Exp.yota);
  tt = tin;

  /*
 *--WRITE POSSIBLE TSTEP VALUES IN DETAIL FILE, IF ANY
 */
  if(Exp.idet != NULL) (void)fprintf(Exp.idet,F9000,tin,tmin,tbnd);

  /*
 *--BEGIN LINE SEARCH IN CHOSEN DIRECTION
 *
 *--SAVE INCOMING VALUES
 */
  for(i=0; i<Exp.nt; i++) ths[i] = theta[i];
  fs = *f;

  /*
   *--INITIALIZE IMPLIED BOUNDARY FLAG
   */
  Ximpbnd = 0;

 S40:
  /*
   *--COMPUTE THETA-1
   */
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) {
      l += 1;
      tht[i] = ths[i] + tt * Xpdir[l - 1];
    }
    else tht[i] = ths[i];
  }

  /*
 *--GET CORRESPONDING FUNCTIONAL VALUE F-1
 */
  mfun->eval_fun(tht,&ft,nfe,&lex);
  if(lex > 0) goto S130;

  if(ft <= *f) goto S110;

 S50:
  /*
   *--IMPROVING; TRY GOING FARTHER
   */
  tt2 = tt + tt;
  if(tt2 < tbnd) goto S60;

  /*
   *--CAN'T GO ANY FARTHER, SO STOP AT THETA-1
   */
  if(Exp.idet != NULL) (void)fprintf(Exp.idet,"%s",F9100);
  Xtstep = tt;
  *f = ft;
  for(i=0; i<Exp.nt; i++) theta[i] = tht[i];
  goto S140;

 S60:
  /*
   *--COMPUTE THETA-2 (FIXED PARAMETERS ALREADY ESTABLISHED; DEPENDENT
   *--ONES WILL BE COMPUTED IN DEPAR VIA FUN)
   */
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    theta[i] = tht[i];
    if(Xist[i] <= 2) {
      l += 1;
      tht[i] = theta[i] + tt * Xpdir[l - 1];
    }
  }

  /*
 *--GET CORRESPONDING FUNCTIONAL VALUE F-2
 */
  mfun->eval_fun(tht,&ft2,nfe,&lex);
  if(lex <= 0) goto S70;

  Ximpbnd = 1;
  goto S80;

 S70:
  if(ft2 > ft) goto S90;

  Ximpbnd = 0;

 S80:
  /*
   *--NO LONGER IMPROVING; STOP AT THETA-1
   */
  Xtstep = tt;
  *f = ft;
  goto S140;

 S90:
  /*
   *--STILL IMPROVING; HOW FAST?
 */
  tt = tt2;
  diff1 = ft2 - ft;
  diff2 = ft - *f;
  if(diff1 > diff2) goto S100;

  /*
   *--SLOWER; TOO SLOW?
 */
  if(diff1 / diff2 >= CURV) goto S100;

  /*
 *--TOO SLOW; STOP AT THETA-2
 */
  Xtstep = tt2;
  *f = ft2;
  for(i=0; i<Exp.nt; i++) theta[i] = tht[i];
  goto S140;

 S100:
  /*
 *--IMPROVING FASTER OR AT LEAST FAST ENOUGH;
 *--PROCEED WITH THETA-2 AS NEW THETA-1
 */
  ft = ft2;
  goto S50;

 S110:
  /*
 *--NOT IMPROVING; TRY NOT GOING SO FAR
 */
  Ximpbnd = 0;
  if(tt <= tmin) {
    /*
     *--CAN'T CUT ANY FARTHER; QUIT WITHOUT COMPLETING ITERATION
     */
    *ifl = 4;
    return;

  }

  /*
 *--HALVE TT AND COMPUTE THETA-2 (FIXED PARAMETERS ALREADY
 *--ESTABLISHED; DEPENDENT ONES WILL BE COMPUTED IN DEPAR VIA FUN)
 */
  tt = .5e0 * tt;
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) {
      l += 1;
      tht[i] = ths[i] + tt * Xpdir[l - 1];
    }
  }

  /*
 *--GET CORRESPONDING FUNCTIONAL VALUE F-2
 */
  mfun->eval_fun(tht,&ft2,nfe,&lex);
  if(lex > 0) goto S130;

  if(ft2 <= *f) goto S120;

  /*
 *--IMPROVED THIS TIME; STOP AT THETA-2
 */
  Xtstep = tt;
  *f = ft2;
  for(i=0; i<Exp.nt; i++) theta[i] = tht[i];
  goto S140;

 S120:
  /*
 *--STILL NOT IMPROVING; REDUCE TT FURTHER BY SCALING FACTOR
 */
  diff1 = ft + *f - ft2 - ft2;
  scal = .1e0;
  if(diff1 < 0.e0) scal = MAX(scal,1.e0+.5e0 * (*f - ft) / diff1);
  tt = scal * tt;

  /*
   *--COMPUTE NEW THETA-2 (FIXED PARAMETERS ALREADY ESTABLISHED;
   *--DEPENDENT ONES WILL BE COMPUTED IN DEPAR VIA FUN)
   */
  l = 0;
  for(i=0; i<Exp.nt; i++) {
    if(Xist[i] <= 2) {
      l += 1;
      tht[i] = ths[i] + tt * Xpdir[l - 1];
    }
  }

  /*
 *--GET CORRESPONDING FUNCITONAL VALUE F-2
 */
  mfun->eval_fun(tht,&ft,nfe,&lex);
  if(lex > 0) goto S130;

  /*
   *--IF STILL NOT IMPROVED, GO BACK TO CUT TT MORE, WITH THETA-2 AS NEW
   *--THETA-1
   */
  if(ft <= *f) goto S110;

  /*
 *--IMPROVED THIS TIME; STOP AT THETA-2
 */
  Xtstep = tt;
  *f = ft;
  for(i=0; i<Exp.nt; i++) theta[i] = tht[i];
  goto S140;

 S130:
  /*
 *--IF RUN INTO UNDEFINED FUNCTION BEFORE ACHIEVING ANY IMPROVEMENT,
 *--HALVE TT AND START OVER
 */
  tt = .5e0 * tt;
  Ximpbnd = 1;
  goto S40;

 S140:
  /*
   *--CHECK WHETHER TSTEP < TMIN
   */
  if(Xtstep < tmin && Xtstep < tin) *ifl = -1;
  else {
    /*
     *--NORMAL TERMINATION
     */
    *ifl = 0;
  }

  /*
 *--FINISH UP THIS ITERATION
 */
  for(i=0; i<Exp.nt; i++) Xthpr[i] = ths[i];
  Xfpr = fs;
  endit(theta,f);
  *ifirst = 0;

  /*
   *--INCREMENT NUMBER OF ITERATIONS WITH THIS GRADIENT VECTOR
   */
  Xigage += 1;

  /*
 *--WRITE OUT DETAILS OF THIS ITERATION IF DESIRED
 */
  if(Exp.idet != NULL) litv(theta,f,nfe);

  /*
 *--CHECK FOR CONVERGENCE
 */
  vcnvch(theta,&lex);

  /*
 *--SET FLAG (IF APPROPRIATE) AND RETURN
 */
  if(lex > 0) *ifl = lex;
  return;

  /*
   * Undefine Local Preprocessor Macro
   */
#undef CURV
}

/*
 * *******************************************************************
 */

void Maxfun_source::prepd(double theta[],
			  double *f,int *nfe,int *lex) {
  /*
 *====================================================================
 *
 *--PREPARE FOR DERIVATIVE COMPUTATION:
 *--CHECK FOR CONVERGENCE TO BOUNDS AND WHETHER # PARAMETERS OK,
 *--INITIALIZE INCREMENT STEPSIZE FACTORS FOR DERIVATIVE COMPUTATION
 *
 *--RETURNS FLAG LEX:
 *--   0:  EVERYTHING OK
 *--   1:  ALL INDEPENDENT PARAMETERS CONVERGED TO BOUNDS
 *--   2:  TOO MANY PARAMETERS TO COMPUTE DERIVATIVES
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static double sinit;
  static int i,lexb;

  /*
   * Format String
   */
  const char* F9000 = "\nPROGRAM DIMENSIONS TOO SMALL TO COMPUTE DE\
RIVATIVES FOR%3d PARAMETERS;\n   MAXIMUM NI IS%3d\n";

  /*
 * Executable Statements
 */

  /*
   *--INITIALIZE COUNTER OF PARAMETERS CONVERGED TO BOUNDS
   */
  Xnb = 0;

  /*
 *--CHECK FOR CONVERGENCE TO BOUNDS
 */
  fixbnd(theta,f,nfe,&lexb);
  if(lexb >= 2) {
    *lex = 1;
    return;

  }

  /*
 *--CHECK FOR TOO MANY PARAMETERS
 */
  if(Xni > MAXFUN_NPV) {
    if(Exp.iout != NULL)
      (void)fprintf(Exp.iout,F9000,Xni,MAXFUN_NPV);
    *lex = 2;
    return;

  }

  /*
 *--INITIALIZE INCREMENT STEPSIZE FACTORS FOR 2ND DERIVATIVE
 *--COMPUTATION
 */
  if(Exp.ihit > 0) sinit = Exp.epsd;
  else sinit = Exp.yota;
  for(i=0; i<Exp.nt; i++) Xstp[i] = sinit;

  *lex = 0;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::fixbnd(double theta[],
			   double *f,int *nfe,int *lex) {
  /*
 *====================================================================
 *
 *--FIX PARAMETERS WHICH ARE CLOSE TO A BOUND
 *
 *--RETURNS FLAG LEX:
 *--   0:  NOTHING NEW FIXED
 *--   1:  FIXED AT LEAST ONE INDEPENDENT PARAMETER
 *--   2:  ALL INDEPENDENT PARAMETERS NOW FIXED
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static double fs,thib,thiz;
  static int i,isti,lexb,nbs;

  /*
   * Format strings
   */
  const char* F9000 = "\nPARAMETER%3d FIXED AT BOUND%#17.8G\n";
  const char* F9100 = "\nPARAMETER%3d FIXED NEAR BOUND%#17.8G\n";
  const char* F9200 = "\nAFTER ITERATION%4d, ALL INDEPENDENT PARAME\
TERS CONVERGED TO BOUNDS\n";

  /*
 * Executable Statements
 */

  /*
   *--SAVE INCOMING VALUE AND INITIALIZE EXIT FLAG
   */
  nbs = Xnb;

  /*
 *--CHECK EACH VARYING INDEPENDENT PARAMETER
 */
  for(i=1; i<=Exp.nt; i++) {
    isti = Xist[i - 1];
    if(isti <= 2) {
      fs = *f;
      thiz = theta[i - 1];
      thib = Exp.thu[i - 1];
      bndchk(&thiz,&thib,&lexb);
      if(lexb > 0) goto S10;

      thib = Exp.thl[i - 1];
      bndchk(&thiz,&thib,&lexb);
      if(lexb <= 0) goto S20;
            
    S10:
      /*
       *--CLOSE TO BOUND; FIX
       */
      Xnb += 1;
      theta[i - 1] = thib;
      mfun->eval_fun(theta,f,nfe,&lexb);
      if(lexb <= 0 && *f >= fs) {

	Xist[i - 1] = isti + 4;
	if(Exp.iout != NULL)
	  (void)fprintf(Exp.iout,F9000,i,thib);
	if(Exp.idet != NULL && Exp.idet != Exp.iout)
	  (void)fprintf(Exp.idet,F9000,i,thib);
      }
      else {

	Xist[i - 1] = isti + 6;
	if(Exp.iout != NULL)
	  (void)fprintf(Exp.iout,F9100,i,thib);
	if(Exp.idet != NULL && Exp.idet != Exp.iout)
	  (void)fprintf(Exp.idet,F9100,i,thib);
	theta[i - 1] = thiz;
	*f = fs;

      }
    }

  S20:;

  }

  /*
 *--FIX UP DEPENDENT PARAMETERS TO GO WITH FINAL THETA
 *--(NEED NOT CHECK RETURN FLAG AS ALREADY TESTED THIS COMBINATION)
 */
  if(Xnd > 0) mfun->dep_fun(theta,&lexb);

  /*
 *--SET # PARAMETERS THAT MAY VARY
 */
  Xnv = Xni - Xnb;

  /*
 *--SEE IF ANY PARAMETERS HAVE BEEN FIXED
 */
  if(Xnb <= nbs) {
    *lex = 0;
    return;

  }

  /*
 *--SEE IF ANY LEFT UNFIXED
 */
  if(Xnv > 0) {
    *lex = 1;
    return;
        
  }
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9200,Xit);
  *lex = 2;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::bndchk(double *th,double *thb,int *lex) {
  /*
   *====================================================================
   *
   *--TEST FOR CLOSENESS TO A BOUND
   *
   *--RETURNS FLAG LEX:
   *--   0:  NOT CLOSE TO BOUND
   *--   1:  CLOSE TO BOUND
   *
   *====================================================================
   */

  /*
   * Local Scalar
   */
  static double rth;

  /*
   * Executable Statements
   */

  /*
   *--PERFORM APPROPRIATE FORM OF TEST
   */
  if(ABS(*th) >= Exp.epsd) goto S10;

  if(ABS(*th - *thb) <= Exp.epsd * Exp.epsd) goto S30;

  goto S20;

 S10:
  rth = *thb/ *th;
  if(rth >= 1.e0 - Exp.epsd && rth <= 1.e0 + Exp.epsd) goto S30;

 S20:
  /*
   *--NOT CLOSE TO BOUND
   */
  *lex = 0;
  return;

 S30:
  /*
   *--CLOSE TO BOUND
   */
  *lex = 1;
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::bcnvch(double theta[],int *lex) {
  /*
   *====================================================================
   *
   *--CHECK FOR CONVERGENCE -- IMPLIED AND REGULAR TESTS AND OTHER STUFF
   *--APPROPRIATE FOR DIRECT SEARCH (BASIC)
   *
   *--RETURNS FLAG LEX:
   *--   0:  FAILED CONVERGENCE TEST; CONTINUE ITERATING
   *--   1:  PASSED CONVERGENCE TEST (CRITERION 1)
   *--   2:  NEGLIGIBLE FUNCTION CHANGE (CONVERGENCE CRITERION 3)
   *--   3:  FAILED CONVERGENCE TEST, BUT REACHED MAXIMUM # ITERATIONS
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static double T1;
  static int D3;
  static double eloc;
  static int i;

  /*
   * Format Strings
   */
  const char* F9000 = "\nBASIC ITERATION PROCESS OF DIRECT SEARCH C\
ONVERGED BY CRITERION 1 AT\n   ITERATION%4d\n";
  const char* F9100 = "(USING EPSILON%#13.6G)\n";
  const char* F9200 = "\nBASIC ITERATION PROCESS STOPPED AFTER ITER\
ATION%4d DUE TO NEGLIGIBLE\n   DIFFERENCE IN COMPUTED FUNCTION VALUES\
%#20.10E\n   (CONVERGENCE CRITERION 3)\n";
  const char* F9300 = "\nSTOPPED WITHOUT CONVERGENCE AT%4d ITERATIO\
NS\n   (REACHED ITERATION LIMIT)\n";

  /*
   * Executable Statements
   */

  /*
   *--INITIALIZE INTERNAL EPSILON
   */
  eloc = Exp.epsc1;

 S10:
  /*
   *--LOOP THROUGH PARAMETERS, USING IMPLIED TEST FOR DEPENDENT ONES
   */
  for(i=1,D3=(Exp.nt - i + 1); D3>0; D3--,i+=1) {
    switch(Xist[i - 1]) {
    case  1: goto S20;

    case  2: goto S20;

    case  3: goto S30;

    case  4: goto S40;

    case  5: goto S40;

    case  6: goto S40;

    case  7: goto S40;

    case  8: goto S40;

    case  9: goto S40;

    case 10: goto S40;

    default: break;
    }

  S20:
    if(Xstp[i - 1] > efn((T1= ABS(Xthpr[i - 1]),&T1),&eloc))
      goto S50;

    goto S40;

  S30:
    itest(theta,&i,&eloc,lex);
    if(*lex <= 0) goto S50;

  S40:;

  }

  /*
   *--PASSED BASIC CONVERGENCE TEST
   */
  if(Exp.iout != NULL) {
    (void)fprintf(Exp.iout,F9000,Xit);
    if(eloc > Exp.epsc1) (void)fprintf(Exp.iout,F9100,eloc);
  }
  *lex = 1;
  return;

 S50:
  /*
   *--FAILED STANDARD CONVERGENCE TEST
   *
   *--TRY INCREASING EPSILON IF NOT ALREADY TRIED AND DIFMAX SMALL
   *--ENOUGH
   */
  if(Exp.epsc1 <= 0.e0 || eloc > Exp.epsc1 || Xdifmax > eloc * eloc)
    goto S60;

  eloc *= 10.e0;
  goto S10;

 S60:
  /*
   *--COMPARE DIFMAX WITH USER-SPECIFIED EPSC3
   */
  if(Xdifmax < Exp.epsc3) {
    if(Exp.iout != NULL)
      (void)fprintf(Exp.iout,F9200,Xit,Xdifmax);
    *lex = 2;
    return;

  }

  /*
   *--GIVE UP ON CONVERGENCE; CHECK FOR MAXIMUM # ITERATIONS
   */
  if(Xit >= Exp.maxit) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9300,Xit);
    *lex = 3;
    return;

  }

  /*
   *--NOT THAT EITHER; CONTINUE ITERATING
   */
  *lex = 0;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::vcnvch(double theta[],int *lex) {
  /*
   *====================================================================
   *
   *--CHECK FOR CONVERGENCE -- IMPLIED TEST APPROPRIATE FOR
   *--NEWTON-RAPHSON OR VARIABLE METRIC METHODS
   *--RETURNS FLAG LEX:
   *--   0:  FAILED CONVERGENCE TEST; CONTINUE ITERATING
   *--   1:  PASSED (IMPLIED) CONVERGENCE TEST (CRITERION 1)
   *--   2:  NEGLIGIBLE FUNCTION CHANGE (CONVERGENCE CRITERION 3)
   *--   3:  FAILED CONVERGENCE TEST, BUT REACHED MAXIMUM # ITERATIONS
   *====================================================================
   *
   */

  /*
   * Local Scalars
   */
  static int D2;
  static int i;

  const char* F9000 = "\nCONVERGED BY CRITERION 1 AT ITERATION%4d\n\
";
  const char* F9100 = "\nSTOPPED AFTER ITERATION%4d DUE TO NEGLIGIB\
LE FUNCTION\n   CHANGE%#18.10E (CONVERGENCE CRITERION 3)\n";
  const char* F9200 = "\nSTOPPED WITHOUT CONVERGENCE AT%4d ITERATIO\
NS\n   (REACHED ITERATION LIMIT)\n";

  /*
   * Executable Statements
   */

  /*
   *--DO IMPLIED TEST FOR EACH VARYING PARAMETER
   */
  for(i=1,D2=(Exp.nt - i + 1); D2>0; D2--,i+=1) {
    if(Xist[i - 1] > 3) goto S10;

    itest(theta,&i,&Exp.epsc1,lex);
    if(*lex <= 0) goto S20;

  S10:;

  }

  /*
   *--CONVERGED
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,Xit);
  return;

 S20:
  /*
   *--DID NOT CONVERGE BY IMPLIED TEST
   *--CHECK FOR NEGLIGIBLE CHANGE IN FUNCTION
   */
  if(Xfch < Exp.epsc3) {
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9100,Xit,Xfch);
    *lex = 2;
    return;

  }

  /*
   *--DID NOT CONVERGE; CHECK FOR TOO MANY ITERATIONS
   */
  if(Xit < Exp.maxit) return;

  /*
   *--STOP AT MAXIMUM # ITERATIONS
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9200,Xit);
  *lex = 3;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::itest(double theta[],int *i,double *eloc,int *lex) {
  /*
   *====================================================================
   *
   *--CHECK FOR CONVERGENCE ON ONE PARAMETER -- IMPLIED TEST
   *
   *--RETURNS FLAG LEX:
   *--   0:  FAILED IMPLIED CONVERGENCE TEST
   *--   1:  PASSED IMPLIED CONVERGENCE TEST
   *
   *====================================================================
   */

  /*
   * Local Scalar
   */
  static double rth;

  /*
   * Executable Statements
   */

  if(ABS(Xthpr[*i - 1]) > *eloc) goto S10;

  if(ABS(Xcth[*i - 1]) > *eloc * *eloc) goto S30;

  goto S20;

 S10:
  rth = theta[*i - 1] / Xthpr[*i - 1];
  if(rth < 1.e0 - *eloc || rth > 1.e0 + *eloc) goto S30;

 S20:
  /*
   *--PASSED TEST
   */
  *lex = 1;
  return;

 S30:
  /*
   *--FAILED TEST
   */
  *lex = 0;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::endit(double theta[],double *f) {
  /*
   *====================================================================
   *
   *--FINISH UP AN ITERATION
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static double er;
  static int i;

  /*
   * Executable Statements
   */

  Xerm = 0.e0;
  for(i=0; i<Exp.nt; i++) {
    er = theta[i] - Xthpr[i];
    Xcth[i] = er;
    er = ABS(er);
    Xerm = MAX(er,Xerm);
  }

  Xfch = *f - Xfpr;

  Xit += 1;

  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::vcmx(double theta[],
			 double *f,int *nfe,int *ih) {
  /*
 *====================================================================
 *
 *--COMPUTE VARIANCE-COVARIANCE MATRIX
 *
 *--RETURNS FLAG IVFL:
 *--   0:  NO PROBLEM WITH H
 *--   1:  ROUND-OFF ERROR IN H
 *--   3:  H COULD NOT BE INVERTED
 *--   4:  H COULD NOT BE COMPUTED
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static int T1,T2,T3;
  static int l1,l2,lex;

  /*
   * Local Array
   */
  static double hn[MAXFUN_NPV][MAXFUN_NPV];

  /*
 * Format Strings
 */
  const char* F9000 = "\nAFTER ITERATION%4d MATRIX OF SECOND PARTIA\
L DERIVATIVES CANNOT BE\n   INVERTED\n";
  const char* F9100 = "\nVARIANCE-COVARIANCE MATRIX\nBEFORE ADDING \
ROWS AND COLUMNS FOR PARAMETERS THAT ARE DEPENDENT OR\n   CONVERGED T\
O A BOUND\n\n";
  const char* F9200 = "\nINVERSE CHECK:  SHOULD BE IDENTITY MATRIX\
\n\n";
  const char* F9300 = "\n\n";

  /*
   * Executable Statements
   */

  /*
   *--INITIALIZE MATRIX COMPUTATION STATUS FLAG
   */
  Xivfl = 0;

  /*
 *--GET SECOND PARTIAL DERIVATIVES
 */
  deriv2(theta,f,nfe,ih,&lex);
  if(lex - 1 < 0) goto S20;

  if(lex - 1 != 0) {
    Xivfl = 4;
    return;

  }
  Xivfl = 1;

 S20:
  /*
 *--PUT NEGATIVE OF MATRIX OF 2ND PARTIALS IN ANOTHER ARRAY
 */
  for(l1=0; l1<Xnv; l1++) {
    for(l2=0; l2<Xnv; l2++) hn[l2][l1] = -Xh[l2][l1];
  }

  /*
 *--INVERT TO GET VARIANCE-COVARIANCE MATRIX
 */
  T1 = MAXFUN_NPV;
  T2 = MAXFUN_NPV;
  mxnvrt(&hn[0][0],&T1,&Xnv,&Xv[0][0],&T2,&lex);
  if(lex > 0) {
    Xivfl = 3;
    if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,Xit);
    Xivage = -1;
    return;

  }

  /*
 *--WRITE OUT MATRIX AND INVERSE CHECK IN DETAIL FILE, IF ANY
 */
  if(Exp.idet != NULL) {

    (void)fprintf(Exp.idet,"%s",F9100);
    T1 = MAXFUN_NPV;
    mout(&Xv[0][0],&T1,&Xnv,&Xnv,Exp.idet);

    T1 = MAXFUN_NPV;
    T2 = MAXFUN_NPV;
    T3 = MAXFUN_NPV;
    /*      mmult(Xh,&T1,&Xnv,&Xnv,Xv,&T2,&Xnv,hn,&T3); */
    mmult(&Xh[0][0],&T1,&Xnv,&Xnv,&Xv[0][0],&T2,&Xnv,&hn[0][0],
	  &T3);
    for(l1=0; l1<Xnv; l1++) {
      for(l2=0; l2<Xnv; l2++) hn[l2][l1] = -hn[l2][l1];
    }
    (void)fprintf(Exp.idet,"%s",F9200);
    T1 = MAXFUN_NPV;
    mout(&hn[0][0],&T1,&Xnv,&Xnv,Exp.idet);
    (void)fprintf(Exp.idet,"%s",F9300);
  }

  /*
 *--NOW HAVE NEW V; INDICATE THAT IT MATCHES CURRENT THETA
 */
  Xivage = 0;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::augv(double theta[],int *ih) {
  /*
   *====================================================================
   *
   *--AUGMENT THE VARIANCE-COVARIANCE MATRIX TO INCLUDE INDEPENDENT
   *--PARAMETERS THAT HAVE CONVERGED TO BOUNDS AND DEPENDENT
   *--PARAMETERS, PRINT IT, AND COMPUTE STANDARD DEVIATIONS
   *
   *--NOTE THAT VALUES COMPUTED HERE CORRESPOND TO THE FINAL ESTIMATES
   *--THETA, WHILE THE REST OF THE MATRIX MAY NOT
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static double T1;
  static int T2;
  static double dthi,s,thi;
  static int i,ii,isti,istii,j,jj,jm,jmm,k,kk,l,lex,ll,m,mm,n,newit;

  /*
   * Local Arrays
   */
  static double df[MAXFUN_NPV][MAXFUN_NPV],prmuv[MAXFUN_NPV],
    thm[MAXFUN_NP],thp[MAXFUN_NP],u[3][MAXFUN_NPV],
    vn[MAXFUN_NPV][MAXFUN_NPV];
  static int jdp[MAXFUN_NPV],newv[MAXFUN_NPV];

  /*
   * Format Strings
   */
  const char* F9000 = "\nROUNDING ERRORS MAY BE AFFECTING THE DERIV\
ATIVE OF%3dTH FUNCTIONALLY\n   DEPENDENT PARAMETER WITH RESPECT TO PA\
RAMETER%3d\nTHE 3 VALUES ARE%#16.8G%#16.8G%#16.8G\n";
  const char* F9100 = "\nVARIANCE-COVARIANCE MATRIX (NOT INCLUDING \
USER-FIXED PARAMETERS):\n\n";
  const char* F9200 = "\nTHIS MATRIX CORRESPONDS TO PARAMETER ESTIM\
ATES OF%4d ITERATIONS AGO\n";
  const char* F9300 = "\nTHIS MATRIX MAY BE AFFECTED BY ROUND-OFF E\
RRORS IN COMPUTATION OF 2ND\n   PARTIAL DERIVATIVES\n";

  /*
   * Executable Statements
   */

  /*
   *--FIRST CONSIDER INDEPENDENT PARAMETERS THAT HAVE CONVERGED TO
   *--BOUNDS
   */
  if(Xnb <= 0) {
    /*
     *--NO INDEPENDENT PARAMETERS CONVERGED TO BOUNDS (NI = NV); JUST
     *--COPY V TO VN AND GO ON
     */
    for(kk=1; kk<=Xni; kk++) {
      for(k=1; k<=Xni; k++)
	vn[kk - 1][k - 1] = Xv[kk - 1][k - 1];
    }
  }
  else {
    /*
     *--I INDEXES ALL NT PARAMETERS
     *--K INDEXES NI INDEPENDENT PARAMETERS
     *--L INDEXES NV ITERABLE (INDEPENDENT, NOT NEAR BOUND) PARAMETERS
     */
    k = l = 0;
    for(i=1; i<=Exp.nt; i++) {
      isti = Xist[i - 1];
      if(isti < 3 || isti > 4) {

	k += 1;

	if(isti > 2) {
	  /*
	   *--PARAMETER I HAS CONVERGED TO A BOUND; INSERT ZEROS
	   */
	  for(kk=1; kk<=Xni; kk++)
	    vn[kk - 1][k - 1] = vn[k - 1][kk - 1] = 0.e0;
	}
	else {
	  /*
	   *--PARAMETER I IS INDEPENDENT AND VARYING
	   */
	  l += 1;

	  /*
	   *--II, KK, LL PROVIDE SIMILAR INDEXING SCHEME FOR 2ND PARAMETER
	   */
	  kk = ll = 0;
	  for(ii=0; ii<i; ii++) {
	    istii = Xist[ii];
	    if(istii < 3 || istii > 4) {

	      kk += 1;

	      if(istii <= 4) {
		/*
		 *--PARAMETER II IS INDEPENDENT AND NOT FIXED
		 */
		ll += 1;
		vn[kk - 1][k - 1] = Xv[ll - 1][l - 1];
		vn[k - 1][kk - 1] = Xv[l - 1][ll - 1];
	      }
	    }
	  }

	}
      }
    }
  }
  /*
   *--NOW CONSIDER DEPENDENT PARAMETERS
   */
  if(Xnd <= 0) {
    /*
     *--NO DEPENDENT PARAMETERS (NE = NI); JUST COPY VN TO AV AND GO ON TO
     *--COMPUTE STANDARD DEVIATIONS
     */
    for(jj=1; jj<=Xne; jj++) {
      for(j=1; j<=Xne; j++)
	Xav[jj - 1][j - 1] = vn[jj - 1][j - 1];
    }
  }
  else {
    /*
     *--THERE ARE DEPENDENT PARAMETERS TO BE TAKEN CARE OF
     *
     *--FIRST MAJOR LOOP IS THROUGH INDEPENDENT PARAMETERS TO OBTAIN
     *--DERIVATIVES OF ALL DEPENDENT PARAMETERS WITH RESPECT TO THEM
     *
     *--I INDEXES ALL NT PARAMETERS
     *--J INDEXES ALL NE NON-FIXED (BY USER) PARAMETERS
     *--K INDEXES NI INDEPENDENT PARAMETERS
     *--M INDEXES ND DEPENDENT PARAMETERS
     */
    j = k = m = 0;

    for(i=1; i<=Exp.nt; i++) {
      isti = Xist[i - 1];
      if(isti != 4) {

	j += 1;

	if(isti == 3) {
	  /*
	   *--STORE INDEX (IN 1 TO NE SYSTEM) TO MTH DEPENDENT PARAMETER
	   */
	  m += 1;
	  jdp[m - 1] = j;
	}
	else {
	  /*
	   *--THIS IS AN INDEPENDENT PARAMETER
	   */
	  k += 1;

	  if(isti != 1) {
	    /*
	     *--FOR INDEPENDENT PARAMETERS WHICH HAVE CONVERGED TO A BOUND OR ARE
	     *--NOT INVOLVED IN FUNCTIONAL RELATIONSHIPS, SET TO ZERO ALL
	     *--DERIVATIVES WITH RESPECT TO THEM
	     */
	    for(mm=1; mm<=Xnd; mm++)
	      df[k - 1][mm - 1] = 0.e0;
	  }
	  else {
	    /*
	     *--PREPARE TO COMPUTE DERIVATIVES WITH RESPECT TO PARAMETER I
	     */
	    thi = theta[i - 1];

	    for(mm=1; mm<=Xnd; mm++) {
	      prmuv[mm - 1] = 999999.e0;
	      newv[mm - 1] = 1;
	    }
	    T1 = ABS(thi);
	    dthi = dfn(&T1,&Xstp[i - 1]);

	    n = 1;
                        
	  S10:
	    /*
	     *--OBTAIN NTH APPROXIMATION TO DERIVATIVE FOR EACH DEPENDENT
	     *--PARAMETER WHICH STILL NEEDS ANOTHER ITERATION
	     */
	    for(ii=0; ii<Exp.nt; ii++)
	      thp[ii] = thm[ii] = theta[ii];
	    thp[i - 1] = thi + dthi;
	    thm[i - 1] = thi - dthi;
	    mfun->dep_fun(thp,&lex);
	    if(lex > 0) goto S20;

	    mfun->dep_fun(thm,&lex);
	    if(lex > 0) goto S20;

	    /*
	     *--II INDEXES ALL NT PARAMETERS
	     *--MM INDEXES ND DEPENDENT PARAMETERS
	     */
	    mm = 0;
	    for(ii=0; ii<Exp.nt; ii++) {
	      if(Xist[ii] == 3) {
		mm += 1;
		if(newv[mm - 1] > 0)
		  u[n - 1][mm - 1] = 
		    (thp[ii] - thm[ii]) /
		    (dthi + dthi);
	      }
	    }

	    /*
	     *--PREPARE TO DO ANOTHER APPROXIMATION IF NECESSARY
	     */
	    if(*ih <= 0) goto S50;

	    n += 1;
	    if(n > 3) goto S30;

	  S20:
	    /*
	     *--PREPARE FOR ITERATION WITH SMALLER STEPSIZE
	     */
	    dthi *= 0.5e0;
	    goto S10;

	  S30:
	    /*
	     *--FOR EACH DEPENDENT PARAMETER STILL BEING WORKED ON, FIT DERIVATIVE
	     *--USING 3 (LATEST) APPROXIMATIONS
	     *
	     *--MM INDEXES ND DEPENDENT PARAMETERS
	     */
	    newit = 0;
	    for(mm=1; mm<=Xnd; mm++) {
	      if(newv[mm - 1] > 0) {
		fitder(&u[0][mm - 1],&u[1][mm - 1],
		       &u[2][mm - 1],
		       &df[k - 1][mm - 1],
		       &prmuv[mm - 1],&lex);
		if(lex - 1 >= 0) {
		  if(lex - 1 != 0) {
		    /*
		     *--POSSIBLE ROUND-OFF ERROR DETECTED IN COMPUTATION OF DERIVATIVE
		     */
		    if(Exp.iout != NULL)
		      (void)fprintf(Exp.iout,
				    F9000,mm,i,
				    u[0][mm-1],
				    u[1][mm-1],
				    u[2][mm-1]);
		  }
		  /*
		   *--GO DIRECTLY HERE IF DERIVATIVE SUCCESSFULLY APPROXIMATED
		   */
		  newv[mm - 1] = 0;
		}
		else {
		  /*
		   *--ANOTHER ITERATION NECESSARY FOR THIS DERIVATIVE
		   */
		  newit = 1;
		}
	      }
	    }

	    if(newit == 0) goto S60;

	    n = 3;
	    goto S20;

	  S50:
	    /*
	     *--STORE 1ST APPROXIMATION (WHEN ONLY DOING ONE) FOR EACH DEPENDENT
	     *--PARAMETER
	     */
	    for(mm=1; mm<=Xnd; mm++)
	      df[k - 1][mm - 1] = u[0][mm - 1];
	  }
                    
	S60:
	  /*
	   *--COPY ROW, COLUMN FOR THIS INDEPENDENT PARAMETER TO OUTPUT
	   *--VARIANCE-COVARIANCE MATRIX
	   *
	   *--II INDEXES ALL NT PARAMETERS
	   *--JJ INDEXES ALL NE NON-FIXED PARAMETERS
	   *--KK INDEXES NI INDEPENDENT PARAMETERS
	   */
	  jj = kk = 0;
	  for(ii=0; ii<i; ii++) {
	    istii = Xist[ii];
	    if(istii != 4) {

	      jj += 1;
                            
	      if(istii != 3) {
                                
		kk += 1;
		Xav[jj - 1][j - 1] =
		  vn[kk - 1][k - 1];
		Xav[j - 1][jj - 1] =
		  vn[k - 1][kk - 1];
	      }
	    }
	  }
                    
	}
      }
    }

    /*
     *--SECOND MAJOR LOOP IS THROUGH DEPENDENT PARAMETERS TO COMPUTE AND
     *--STORE CORRESPONDING VARIANCES, COVARIANCES
     *
     *--M INDEXES ND DEPENDENT PARAMETERS
     *--JM GIVES CORRESPONDING INDEX (IN 1 TO NE SYSTEM) FOR OUTPUT MATRIX
     */
    for(m=1; m<=Xnd; m++) {
      jm = jdp[m - 1];

      /*
       *--FIRST DO COVARIANCES BETWEEN THIS DEPENDENT PARAMETER AND EACH
       *--INDEPENDENT PARAMETER
       *
       *--I INDEXES ALL NT PARAMETERS
       *--J INDEXES NE NON-FIXED PARAMETERS
       *--K INDEXES NI INDEPENDENT PARAMETERS
       */
      j = k = 0;
      for(i=1; i<=Exp.nt; i++) {
	isti = Xist[i - 1];
	if(isti != 4) {
                    
	  j += 1;
                    
	  if(isti != 3) {
                        
	    k += 1;
	    s = 0.e0;
	    for(kk=1; kk<=Xni; kk++)
	      s += (df[kk - 1][m - 1] *
		    vn[k - 1][kk - 1]);
	    Xav[j - 1][jm - 1] = Xav[jm - 1][j - 1] = s;
                        
	  }
	}
      }

      /*
       *--THEN DO COVARIANCES BETWEEN THIS DEPENDENT PARAMETER AND OTHER
       *--DEPENDENT PARAMETERS (INCLUDING ITS OWN VARIANCE)
       *
       *--MM INDEXES ND DEPENDENT PARAMETERS
       *--JMM GIVES CORRESPONDING INDEX (IN 1 TO NE SYSTEM) FOR OUTPUT
       *--MATRIX
       */
      for(mm=1; mm<=m; mm++) {
	jmm = jdp[mm - 1];
	s = 0.e0;
	for(kk=1; kk<=Xni; kk++) {
	  for(k=1; k<=Xni; k++)
	    s += (df[kk - 1][m - 1] * df[k - 1][mm - 1] *
		  vn[k - 1][kk - 1]);
	}
	Xav[jmm - 1][jm - 1] = Xav[jm - 1][jmm - 1] = s;
      }
    }
  }
    
  /*
   *--PRINT OUT VARIANCE-COVARIANCE MATRIX AND RELATED INFORMATION
   */
  if(Exp.iout != NULL) {
    (void)fprintf(Exp.iout,"%s",F9100);
    T2 = MAXFUN_NP;
    mout(&Xav[0][0],&T2,&Xne,&Xne,Exp.iout);
    (void)fprintf(Exp.iout,F9200,Xivage);
    if(Xivfl > 0) (void)fprintf(Exp.iout,"%s",F9300);
  }

  /*
   *--COMPUTE STANDARD DEVIATIONS WHERE POSSIBLE
   */
  j = 0;
  for(i=1; i<=Exp.nt; i++) {
    if(Xist[i - 1] == 4) Xstde[i - 1] = 0.e0;
    else {
      j += 1;
      s = Xav[j - 1][j - 1];
      if(s < 0.e0) {
	Xstde[i - 1] = s;
	Xivfl = 2;
      }
      else Xstde[i - 1] = sqrt(s);
    }
  }

  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::deriv1(double theta[],double *f,int *nfe) {
  /*
   *====================================================================
   *
   *--COMPUTE GRADIENT (VECTOR OF FIRST PARTIAL DERIVATIVES) AND ITS
   *--NORM
   *
   *--RETURNS FLAG IGFL:
   *--   0:  DERIVATIVES SUCCESSFULLY COMPUTED
   *--   1:  DERIVATIVES COULD NOT BE COMPUTED DUE TO UNDEFINED
   *--       FUNCTION
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static double T1;
  static double dthi,fm,fp,si,thiz,thiy;
  static int i,l,lex;

  /*
   * Local Array
   */
  static double thy[MAXFUN_NP];

  /*
   * Format String
   */
  const char* F9000 = "\nAFTER ITERATION%4d,\n    1 OR MORE 1ST PAR\
TIAL DERIVATIVES COULD NOT BE COMPUTED, AS\n       PARAMETER%4d CLOS\
E TO EXPLICIT OR IMPLIED BOUNDARY\n";
    
  /*
   * Executable Statements
   */

  /*
   *--INITIALIZE
   */
  Xigage = -1;
  for(i=1; i<=Exp.nt; i++) thy[i - 1] = theta[i - 1];
  Xgtg = 0.e0;

  /*
   *--LOOP THROUGH NV INDEPENDENT VARYING PARAMETERS
   */
  l = 0;
  for(i=1; i<=Exp.nt; i++) {
    if(Xist[i - 1] > 2) goto S70;

    l += 1;
    thiz = thy[i - 1];
    si = Exp.yota;
        
  S10:
    T1 = ABS(thiz);
    dthi = dfn(&T1,&si);
    thiy = thiz + dthi;
    if(thiy > Exp.thu[i - 1]) goto S30;

    thy[i - 1] = thiy;
    mfun->eval_fun(thy,&fp,nfe,&lex);
    if(lex <= 0) goto S40;

  S20:
    /*
     *--RUNNING INTO IMPLIED BOUNDARY
     */
    Ximpbnd = 1;

  S30:
    /*
     *--RUNNING INTO TROUBLE; DECREASE INCREMENT
     */
    si = 0.5e0 * si;
    if(si < Exp.epsd * Exp.epsd / 8.e0) goto S80;

    goto S10;

  S40:
    /*
     *--DIVIDE FORWARD, CENTRAL DIFFERENCE PATHS
     */
    if(Xidif > 1) goto S50;

    /*
     *--FORWARD DIFFERENCE
     */
    Xg[l - 1] = (fp - *f) / dthi;
    goto S60;

  S50:
    /*
     *--CENTRAL DIFFERENCE
     */
    thiy = thiz - dthi;
    if(thiy < Exp.thl[i - 1]) goto S30;

    thy[i - 1] = thiy;
    mfun->eval_fun(thy,&fm,nfe,&lex);
    if(lex > 0) goto S20;

    Xg[l - 1] = (fp - fm) / (dthi + dthi);
        
  S60:
    Xgtg += pow(Xg[l - 1],2.0);
    thy[i - 1] = thiz;

  S70:;
        
  }

  /*
   *--FINISH NORM OF GRADIENT
   */
  Xgtg = sqrt(Xgtg);

  /*
   *--INDICATE HAVE G CORRESPONDING TO CURRENT THETA
   */
  Xigfl = Xigage = 0;
    
  return;

 S80:
  /*
   *--ERROR EXIT; PARAMETER I STUCK
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9000,Xit,i);
  Xigfl = 1;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::deriv2(double theta[],
			   double *f,int *nfe,int *ih,int *lex) {
  /*
 *====================================================================
 *
 *--COMPUTE 2ND PARTIAL DERIVATIVES OF THE FUNCTION
 *
 *--RETURNS FLAG LEX:
 *--   0:  NO PROBLEM WITH H
 *--   1:  ROUND-OFF ERROR IN H
 *--   2:  COULD NOT COMPUTE H
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static int T1;
  static double athi,athii,dthi,dthii,e2d8,fm,fmm,fmp,fp,fpm,fpp,
    prmuh,si,sii,thi,thii,thim,thip;
  static int i,ii,isti,istii,l,lexf,ll,n,nk,nl,nm;

  /*
   * Local Arrays
   */
  static double hn[3][MAXFUN_NPV][MAXFUN_NPV],thy[MAXFUN_NP];
  static int lrh[MAXFUN_NPV][MAXFUN_NPV],nbd[4];

  /*
   * Format Strings
   */
  const char* F9000 = "\nCOMPUTING SECOND PARTIAL DERIVATIVES USIN\
G ITERATIVE APPROXIMATION\n";
  const char* F9100 = "\nCOMPUTING SECOND PARTIAL DERIVATIVES USIN\
G A SINGLE APPROXIMATION\n";
  const char* F9200 = "\nMATRIX OF SECOND PARTIAL DERIVATIVES\nNOT \
COUNTING PARAMETERS THAT ARE FIXED, DEPENDENT,\n   OR CONVERGED TO A \
BOUND\n\n";
  const char* F9300 = "\nROUNDING ERRORS MAY BE AFFECTING THE FOLLO\
WING ELEMENTS IN THE MATRIX\n   OF SECOND PARTIAL DERIVATIVES:\n\n  \
L LL       D1(L,LL)        D2(L,LL)        D3(L,LL)\n";
  const char* F9400 = "%3d%3d%#16.8G%#16.8G%#16.8G\n";
  const char* F9500 = "\nAFTER ITERATION%4d,\n   1 OR MORE 2ND PART\
IAL DERIVATIVES COULD NOT BE COMPUTED,\n      AS PARAMETER%4d CLOSE T\
O EXPLICIT OR IMPLIED BOUNDARY\n";

  /*
 * Executable Statements
 */

  /*
   *--INDICATE TYPE OF APPROXIMATION USED FOR SECOND DERIVATIVES
   */
  if(Exp.idet != NULL) {
    if(*ih > 0) (void)fprintf(Exp.idet,"%s",F9000);
    else (void)fprintf(Exp.idet,"%s",F9100);
  }
    
  /*
 *--INITIALIZE
 */
  *lex = 0;
  e2d8 = Exp.epsd * Exp.epsd / 8.e0;
  for(i=1; i<=Exp.nt; i++) thy[i - 1] = theta[i - 1];
  for(ll=1; ll<=Xnv; ll++) {
    for(l=1; l<=Xnv; l++) lrh[ll - 1][l - 1] = 0;
  }

  /*
 *--LOOP THROUGH INDEPENDENT, VARYING PARAMETERS
 */
  l = 0;
  for(i=1; i<=Exp.nt; i++) {
    isti = Xist[i - 1];
    if(isti > 2) goto S400;

    l += 1;
    thi = thy[i - 1];
    athi = ABS(thi);
    si = Xstp[i - 1];
    dthi = dfn(&athi,&si);

  S10:
    /*
     *--CHECK WHETHER NEIGHBORING VALUES ARE OK WITH RESPECT TO BOUNDS,
     *--DEPENDENT PARAMETERS
     */
    thip = thi + dthi;
    thim = thi - dthi;
    if(thip <= Exp.thu[i - 1] && thim >= Exp.thl[i - 1]) goto S20;

    si *= 0.5e0;
    if(si < e2d8) goto S430;

    dthi = 0.5e0 * dthi;
    goto S10;

  S20:
    if(isti == 2) goto S50;

  S30:
    thy[i - 1] = thip;
    mfun->dep_fun(thy,&lexf);
    if(lexf > 0) goto S40;

    thy[i - 1] = thim;
    mfun->dep_fun(thy,&lexf);
    if(lexf <= 0) goto S50;

  S40:
    Ximpbnd = 1;
    si *= 0.5e0;
    if(si < e2d8) goto S430;

    dthi = 0.5e0 * dthi;
    thim = thi - dthi;
    thip = thi + dthi;
    goto S30;

  S50:
    /*
     *--RESTORE THY(I) AND SAVE CURRENT VALUE OF STEPSIZE FACTOR
     */
    Xstp[i - 1] = si;
    thy[i - 1] = thi;

    /*
     *--COMPUTE 2ND PARTIAL DERIVATIVES FOR (I,II) AND (II,I) PAIRS
     *--WHERE II < I
     */
    if(i == 1) goto S310;

    ll = 0;
    for(ii=1; ii<=i - 1; ii++) {
      istii = Xist[ii - 1];
      if(istii > 2) goto S300;

      ll += 1;
            
      thii = thy[ii - 1];
      athii = ABS(thii);
      sii = Xstp[ii - 1];
      dthii = dfn(&athii,&sii);
      si = Xstp[i - 1];
      dthi = dfn(&athi,&si);

      /*
       *--IF BOTH INDEPENDENT PARAMETERS ARE INVOLVED IN FUNCTIONAL
       *--RELATIONSHIPS, CHECK ALL COMBINATIONS OF NEIGHBORING PARAMETER
       *--PAIRS WITH RESPECT TO DEPENDENT PARAMETERS
       */
      if(isti == 2 || istii == 2) goto S210;

    S60:
      nl = 1;

    S70:
      for(nk=1; nk<=4; nk++) nbd[nk - 1] = 0;
      thy[i - 1] = thi + dthi;
      thy[ii - 1] = thii + dthii;
      nm = 0;
      nk = 1;

    S80:
      mfun->dep_fun(thy,&lexf);
      if(lexf > 0) {
	Ximpbnd = nbd[nk - 1] = 1;
	nm += 1;
      }
      switch(nk) {
      case 1: goto S90;

      case 2: goto S100;

      case 3: goto S110;

      case 4: goto S120;

      default: break;
      }

    S90:
      thy[ii - 1] = thii - dthii;
      nk = 2;
      goto S80;

    S100:
      thy[i - 1] = thi - dthi;
      thy[ii - 1] = thii + dthii;
      nk = 3;
      goto S80;

    S110:
      thy[ii - 1] = thii - dthii;
      nk = 4;
      goto S80;

    S120:
      if(nm <= 0) goto S210;

      switch(nm) {
      case 1: goto S130;

      case 2: goto S170;

      case 3: goto S180;

      case 4: goto S180;

      default: break;
      }

    S130:
      switch(nl) {
      case 1: goto S140;

      case 2: goto S150;

      case 3: goto S160;

      case 4: goto S140;

      default: break;
      }

    S140:
      si *= 0.5e0;
      if(si < e2d8) goto S430;

      dthi = 0.5e0 * dthi;
      nl = 2;
      goto S70;

    S150:
      sii *= 0.5e0;
      if(sii < e2d8) goto S420;

      dthii = 0.5e0 * dthii;
      si += si;
      dthi += dthi;
      nl = 3;
      goto S70;

    S160:
      si *= 0.5e0;
      if(si < e2d8) goto S430;

      dthi = 0.5e0 * dthi;
      nl = 4;
      goto S70;

    S170:
      if(nbd[0] * nbd[1] > 0 || nbd[2] * nbd[3] > 0) goto S200;

      if(nbd[0] * nbd[2] > 0 || nbd[1] * nbd[3] > 0) goto S190;

    S180:
      si *= 0.5e0;
      if(si < e2d8) goto S430;

      dthi = 0.5e0 * dthi;

    S190:
      sii *= 0.5e0;
      if(sii < e2d8) goto S420;

      dthii = 0.5e0 * dthii;
      goto S60;

    S200:
      si *= 0.5e0;
      if(si < e2d8) goto S430;

      dthi = 0.5e0 * dthi;
      goto S60;

    S210:
      /*
       *--READY TO GO AHEAD WITH DERIVATIVE ESTIMATION
       */
      prmuh = 999999.e0;
      n = 1;

    S220:
      /*
       *--COMPUTE NTH APPROXIMATION TO 2ND PARTIAL DERIVATIVE
       */
      thy[i - 1] = thi + dthi;
      thy[ii - 1] = thii + dthii;
      mfun->eval_fun(thy,&fpp,nfe,&lexf);
      if(lexf > 0) goto S230;

      thy[i - 1] = thi - dthi;
      mfun->eval_fun(thy,&fmp,nfe,&lexf);
      if(lexf > 0) goto S230;

      thy[ii - 1] = thii - dthii;
      mfun->eval_fun(thy,&fmm,nfe,&lexf);
      if(lexf > 0) goto S230;

      thy[i - 1] = thi + dthi;
      mfun->eval_fun(thy,&fpm,nfe,&lexf);
      if(lexf <= 0) goto S240;

    S230:
      Ximpbnd = 1;
      if(n > 1) goto S410;

      si *= 0.5e0;
      if(si < e2d8) goto S410;

      sii *= 0.5e0;
      if(sii < e2d8) goto S410;

      dthi *= 0.5e0;
      dthii *= 0.5e0;
      goto S220;
            
    S240:
      hn[n - 1][ll - 1][l - 1] =
	(fpp - fmp - fpm + fmm) / (4.e0 * dthi * dthii);

      /*
       *--STORE FIRST APPROXIMATION IF ONLY DOING ONE
       */
      if(*ih > 0) goto S250;

      Xh[ll - 1][l - 1] = hn[0][ll - 1][l - 1];
      goto S290;

    S250:
      /*
       *--PREPARE TO DO ANOTHER APPROXIMATION IF APPROPRIATE
       */
      n += 1;
      if(n <= 3) goto S270;

      /*
       *--FIT DERIVATIVES FROM 3 APPROXIMATIONS
       */
      fitder(&hn[0][ll - 1][l - 1],&hn[1][ll - 1][l - 1],
	     &hn[2][ll - 1][l - 1],&Xh[ll - 1][l - 1],&prmuh,
	     &lexf);
      if(lexf - 1 > 0) goto S280;

      if(lexf - 1 == 0) goto S290;

      /*
       *--NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE
       */
      n = 3;

    S270:
      si *= 0.5e0;
      dthi *= 0.5e0;
      sii *= 0.5e0;
      dthii *= 0.5e0;
      goto S220;

    S280:
      /*
       *--ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING
       *--2ND DERIVATIVE
       */
      lrh[ll - 1][l - 1] = *lex = 1;

    S290:
      /*
       *--COPY 2ND DERIVATIVE TO ELEMENT (LL,L) OF SYMMETRIC MATRIX
       */
      Xh[l - 1][ll - 1] = Xh[ll - 1][l - 1];

      /*
       *--FINISH UP LOOP
       *
       *--RESTORE ORIGINAL VALUES OF ITH, IITH PARAMETERS
       */
      thy[i - 1] = thi;
      thy[ii - 1] = thii;
            
    S300:;

    }

  S310:
    /*
     *--TAKE CARE OF (I,I) PAIR
     */
    si = Xstp[i - 1];
    dthi = dfn(&athi,&si);

    prmuh = 999999.e0;
    n = 1;

  S320:
    /*
     *--START HERE TO WORK ON NTH APPROXIMATION
     */
    thy[i - 1] = thi + dthi;
    mfun->eval_fun(thy,&fp,nfe,&lexf);
    if(lexf > 0) goto S330;

    thy[i - 1] = thi - dthi;
    mfun->eval_fun(thy,&fm,nfe,&lexf);
    if(lexf <= 0) goto S340;

  S330:
    Ximpbnd = 1;
    if(n > 1) goto S430;

    si *= 0.5e0;
    if(si < e2d8) goto S430;

    dthi *= 0.5e0;
    goto S320;
        
  S340:
    hn[n - 1][l - 1][l - 1] =
      (fp - *f - *f + fm) / (dthi * dthi);

    /*
     *--STORE FIRST APPROXIMATION IF ONLY DOING ONE
     */
    if(*ih > 0) goto S350;

    Xh[l - 1][l - 1] = hn[0][l - 1][l - 1];
    goto S390;

  S350:
    /*
     *--PREPARE TO DO ANOTHER APPROXIMATION IF NECESSARY
     */
    n += 1;
    if(n <= 3) goto S370;

    /*
     *--FIT DERIVATIVES FROM 3 APPROXIMATIONS
     */
    fitder(&hn[0][l - 1][l - 1],&hn[1][l - 1][l - 1],
	   &hn[2][l - 1][l - 1],&Xh[l - 1][l - 1],&prmuh,&lexf);
    if(lexf - 1 > 0) goto S380;

    if(lexf - 1 == 0) goto S390;

    /*
     *--NEED TO DO ANOTHER ITERATION FOR 2ND DERIVATIVE
     */
    n = 3;
        
  S370:
    si *= 0.5e0;
    dthi = 0.5e0 * dthi;
    goto S320;

  S380:
    /*
     *--ITERATION FINISHED BUT THERE IS ROUND-OFF ERROR IN COMPUTING
     *--2ND DERIVATIVE
     */
    lrh[l - 1][l - 1] = *lex = 1;

  S390:
    /*
     *--RESTORE ORIGINAL VALUE OF I'TH PARAMETER
     */
    thy[i - 1] = thi;

  S400:;

  }

  /*
 *--PRINT H AND WARNINGS ABOUT ROUND-OFF ERRORS IN DETAIL FILE, IF ANY
 */
  if(Exp.idet == NULL) return;

  (void)fprintf(Exp.idet,"%s",F9200);
  T1 = MAXFUN_NPV;
  mout(&Xh[0][0],&T1,&Xnv,&Xnv,Exp.idet);

  if(*lex <= 0) return;

  /*
 *--PRINT WARNINGS ABOUT 2ND PARTIAL DERIVATIVES
 */
  (void)fprintf(Exp.idet,"%s",F9300);
  for(l=1; l<=Xnv; l++) {
    for(ll=1; ll<=l; ll++) {
      if(lrh[ll - 1][l - 1] > 0)
	(void)fprintf(Exp.idet,F9400,l,ll,hn[0][ll-1][l-1],
		      hn[1][ll-1][l-1],hn[2][ll-1][l-1]);
    }
  }

  return;

 S410:
  /*
 *--ERROR EXIT; PARAMETERS II AND I STUCK
 */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9500,Xit,ii);
  goto S430;

 S420:
  /*
 *--ERROR EXIT; PARAMETER II STUCK
 */
  i = ii;

 S430:
  /*
   *--ERROR EXIT; PARAMETER I STUCK
   */
  if(Exp.iout != NULL) (void)fprintf(Exp.iout,F9500,Xit,i);

  *lex = 2;
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::fitder(double *d1,double *d2,double *d3,double *dd,
			   double *prmu,int *lex) {
  /*
 *====================================================================
 *
 *--FITDER COMPUTES AN IMPROVED VALUE FOR A NUMERICAL DERIVATIVE
 *--USING 3 ESTIMATES WHICH HAVE BEEN COMPUTED USING 3 DIFFERENT
 *--STEP SIZES.
 *
 *--RETURNS FLAG LEX:
 *--   0:  NEED ANOTHER ITERATION
 *--   1:  SUCCESSFUL DERIVATIVE ESTIMATION
 *--   2:  ROUND-OFF ERROR
 *
 *====================================================================
 */

  /*
   * Define Local Preprocessor Macro
   */
#define ZER 0.1e-07

  /*
 * Local Scalars
 */
  static double r,rden,rmu,teps;

  /*
 * Executable Statements
 */

  /*
   *--SEE IF ALL ESTIMATES < 10**-8
   */
  if(!(ABS(*d1) >= ZER || ABS(*d2) >= ZER || ABS(*d3) >= ZER)) {
    *dd = 0.e0;
    *lex = 1;
    return;

  }

  /*
 *--COMPUTE RATIO R
 */
  rden = *d3 + *d3 - *d2;
  /*   if(rden == 0.e0) r = 0.e0; */
  if(ABS(rden) < EPSILON) r = 0.e0;
  else {
    r = (*d2 + *d2 - *d1) / rden;
  }

  /*
 *--IF R "CLOSE" TO 1, FINISHED
 */
  teps = 10.e0 * Exp.epsd;
  if(r >= 1.e0 - teps && r <= 1.e0 + teps) goto S10;

  /*
   *--UNLESS R GETTING FARTHER FROM 1, DO ANOTHER ITERATION
   */
  rmu = ABS(r - 1.e0);
  if(rmu <= *prmu) {
    *prmu = rmu;
    *d1 = *d2;
    *d2 = *d3;
    *lex = 0;
    return;
        
  }
  *dd = *d1;
  *lex = 2;
  return;

 S10:
  *dd = (*d1 - 6.0e0 * *d2 + 8.0e0 * *d3) / 3.0e0;
  *lex = 1;
  return;

  /*
   * Undefine Local Preprocessor Macro
   */
#undef ZER
}

/*
 * *******************************************************************
 */

double Maxfun_source::dfn(double *ath,double *sf) {
  /*
   *====================================================================
   *
   *--COMPUTES INCREMENT ("DELTA THETA") AS A FUNCTION OF ATH AND
   *--STEPSIZE FACTOR SF, WITH PARAMETER TAU
   *
   *====================================================================
   */

  /*
   * Define Local Preprocessor Macro
   */
#define TAU 0.5e0

  /*
   * Return Value
   */
  static double return_val;

  /*
   * Executable Statements
   */
  if(*ath > TAU) return_val = *ath;
  else return_val = TAU;
  return_val = *sf * return_val;
  return return_val;

  /*
   * Undefine Local Preprocessor Macro
   */
#undef TAU
}

/*
 * *******************************************************************
 */

double Maxfun_source::dfninv(double *ath,double *delth) {
  /*
   *====================================================================
   *
   *--COMPUTES INCREMENT STEPSIZE FACTOR AS A FUNCTION OF ATH AND
   *--INCREMENT DELTH, WITH PARAMETER TAU
   *
   *====================================================================
   */

  /*
   * Define Local Preprocessor Macro
   */
#define TAU 0.5e0

  /*
   * Return Value
   */
  static double return_val;

  /*
   * Executable Statements
   */

  if(*ath > TAU) return_val = *delth/ *ath;
  else return_val = *delth / TAU;
  return return_val;
    
  /*
   * Undefine Local Preprocessor Macro
   */
#undef TAU
}

/*
 * *******************************************************************
 */

double Maxfun_source::efn(double *ath,double *eloc) {
  /*
   *====================================================================
   *
   *--COMPUTES MINIMAL VALUE OF INCREMENT STEPSIZE FACTOR (DELTA) AS A
   *--FUNCTION OF ATH, WITH PARAMETERS ELOC AND TAU
   *
   *====================================================================
   */

  /*
   * Define Local Preprocessor Macro
   */
#define TAU 0.5e0

  /*
   * Return Value
   */
  static double return_val;

  /*
   * Executable Statements
   */

  return_val = *eloc;
  if(*ath > TAU) return return_val;

  return_val /= TAU;
  if(*ath > *eloc) return_val *= *ath;
  else return_val *= *eloc;
  return return_val;

  /*
   * Undefine Local Preprocessor Macro
   */
#undef TAU
}

/*
 * *******************************************************************
 */

void Maxfun_source::mxnvrt(double *a,int *mra,int *m,double *b,int *mrb,
			   int *lex) {
  /*
 *====================================================================
 *
 *--INVERT POSITIVE DEFINITE MATRIX A (M BY M) AND PLACE INVERSE
 *--MATRIX IN B (M BY M)
 *
 *--GAUSS-JORDAN ALGORITHM IS PERFORMED WITHOUT PERMUTATION OF ROWS:
 *--A ZERO PIVOT ELEMENT AT ANY STEP CAUSES ERROR EXIT
 *
 *--RETURNS FLAG LEX:
 *--   0:  NO ERROR
 *--   1:  MATRIX CANNOT BE INVERTED BY THIS SUBROUTINE
 *--   2:  INVALID ARGUMENT: ORDER OF MATRIX <= 0
 *
 *--NOTE THAT INDICES I, J, K ARE LOCAL TO THIS SUBROUTINE AND MAY NOT
 *--CORRESPOND TO INDICES OF THE SAME NAMES IN OTHER PROGRAM UNITS
 *
 *====================================================================
 */

  /*
   * Local Scalars
   */
  static double d,p;
  static int i,j,k;

  /*
   * Executable Statements
   */

  if(*m - 1 > 0) goto S20;

  if(*m - 1 != 0) {
    /*
     *--INVALID ORDER OF MATRIX
     */
    *lex = 2;
    return;

  }

  /*
 *--SCALAR CASE:  M = 1
 */
  d = *(a + 0 + 0 * *mra);
  /*   if(d == 0.e0) goto S40; */
  if(ABS(d) < EPSILON) goto S40;

  *(b + 0 + 0 * *mrb) = 1.e0 / d;
  goto S30;

 S20:
  /*
   *--USUAL CASE:  M > 1
   *
   *--INITIALIZE
   */
  for(j=0; j<*m; j++) {
    for(i=0; i<*m; i++) *(b + i + j * *mrb) = *(a + i + j * *mra);
  }

  /*
 *--LOOP THROUGH PIVOT ELEMENTS
 */
  for(k=1; k<=*m; k++) {

    /*
     *--GET PIVOT ELEMENT
     */
    p = *(b + k - 1 + (k - 1) * *mrb);

    /*
     *--CHECK FOR ZERO PIVOT ELEMENT
     */
    /*      if(p == 0.e0) goto S40; */
    if(ABS(p) < EPSILON) goto S40;

    /*
     *--INVERT PIVOT ELEMENT
     */
    p = 1.e0 / p;

    /*
     *--PROCESS K'TH COLUMN (EXCEPT K'TH ROW ELEMENT)
     */
    for(i=1; i<=*m; i++) {
      if(i != k) *(b + i - 1 + (k - 1) * *mrb) *= p;
    }

    /*
     *--LOOP THROUGH OTHER COLUMNS
     */
    for(j=1; j<=*m; j++) {
      if(j != k) {
	/*
	 *--PROCESS J'TH COLUMN (EXCEPT K'TH ROW ELEMENT)
	 */
	for(i=1; i<=*m; i++) {
	  if(i != k)
	    *(b + i - 1 + (j - 1) * *mrb) -=
	      (*(b + i - 1 + (k - 1) * *mrb) *
	       *(b + k - 1 + (j - 1) * *mrb));
	}
                
      }
    }

    /*
     *--NOW TAKE CARE OF K'TH ROW
     */
    *(b + k - 1 + (k - 1) * *mrb) = p;
        
    p = -p;
    for(j=1; j<=*m; j++) {
      if(j != k)
	*(b + k - 1 + (j - 1) * *mrb) =
	  p * *(b + k - 1 + (j - 1) * *mrb);
    }
  }

 S30:
  /*
   *--NORMAL RETURN
   */
  *lex = 0;
  return;

 S40:
  /*
 *--ERROR EXIT IF MATRIX CANNOT BE INVERTED
 */
  *lex = 1;
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::mmult(double *a,int *mra,int *ma,int *namb,double *b,
			  int *mrb,int *nb,double *c,int *mrc) {
  /*
 *====================================================================
 *
 *--MULTIPLY MATRIX A (MA BY NAMB) BY MATRIX B (NAMB BY NB) AND STORE
 *--THE RESULT IN MATRIX C (MA BY NB)
 *
 *--NOTE THAT INDICES I, J, K ARE LOCAL TO THIS GENERAL SUBROUTINE AND
 *--MAY NOT CORRESPOND TO INDICES OF THE SAME NAMES IN OTHER PROGRAM
 *--UNITS
 *
 *====================================================================
 */
    
  /*
   * Local Scalars
   */
  static double z;
  static int i,j,k;

  /*
   * Executable Statements
   */

  for(i=0; i<*ma; i++) {
    for(j=0; j<*nb; j++) {
      z = 0.e0;
      for(k=0; k<*namb; k++)
	z += (*(a + i + k * *mra) * *(b + k + j * *mrb));
      *(c + i + j * *mrc) = z;
    }
  }

  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::mout(double *a,int *mra,int *m,int *n,FILE *ipr) {
  /*
   *====================================================================
   *
   *--WRITE (IN IPR) MATRIX A (M BY N) WHOSE REAL FIRST DIMENSION SIZE
   *--IS MRA
   *
   *--NOTE THAT INDICES I, J ARE LOCAL TO THIS GENERAL SUBROUTINE AND
   *--MAY NOT CORRESPOND TO INDICES OF THE SAME NAMES IN OTHER PROGRAM
   *--UNITS
   *
   *====================================================================
   */

  /*
   * N_FULL_ROW is the maximum number of 14 wide floats to be printed on
   * a single output line.  For 132 column output, it is 9; for 80
   * column output, it is 5.
   */

  /*
   * Define Local Preprocessor Macro
   */
#define N_FULL_ROW 5

  /*
   * Local Scalars
   */
  static int i,j,jf,ji;
    
  /*
   * Format Strings
   */
  const char* F9010 = "%#14.6E";
  const char* F9020 = "\n";
  const char* F9100 = "\n\n";

  /*
   * Executable Statements
   */

  ji = 1;
  if(*n <= N_FULL_ROW) goto S20;

  jf = N_FULL_ROW;
  goto S30;

 S20:
  jf = *n;

 S30:
  for(i=0; i<*m; i++) {
    for(j=ji - 1; j<jf; j++)
      (void)fprintf(ipr,F9010,*(a + i + j * *mra));
    (void)fprintf(ipr,"%s",F9020);
  }
  if(jf >= *n) return;

  (void)fprintf(ipr,"%s",F9100);
  ji += N_FULL_ROW;
  jf += N_FULL_ROW;
  if(jf > *n) goto S20;

  goto S30;
    
  /*
   * Undefine Local Preprocessor Macro
   */
#undef N_FULL_ROW
}

/*
 * *******************************************************************
 */

void Maxfun_source::lbd(double theta[],double *f) {
  /*
   *====================================================================
   *
   *--LIST (IN IOUT) INITIAL CONDITIONS FOR DIRECT SEARCH METHOD
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static int i,isti;

  /*
   * Format Strings
   */
  const char* F9000 = "\nINITIAL CONDITIONS\n\nPARAMETE\
R                                                      STATUS\
\n       LOWER BOUND      INITIAL EST      UPPER BOUND     IN STEP FA\
CT\n";
  const char* F9100 = "\n%2d  %28s  %36s\n  %#17.8G%#17.8G%#17.8G\
%#17.8G\n";
  const char* F9200 = "  INITIAL FUNCTION VALUE%#18.10E\n";
  const char* F9300 = "--------------------------------------------\
--------------------------\n";

  /*
   * Executable Statements
   */

  (void)fprintf(Exp.iout,"%s",F9000);
  for(i=1; i<=Exp.nt; i++) {
    isti = Exp.istin[i - 1];
    (void)fprintf(Exp.iout,F9100,i,Exp.label[i - 1],
		  Xstatus[isti - 1],Exp.thl[i - 1],theta[i - 1],
		  Exp.thu[i - 1],Exp.stpin[i - 1]);
  }
  if(Exp.lprt > 0) (void)fprintf(Exp.iout,F9200,*f);
  (void)fprintf(Exp.iout,"%s",F9300);
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::lbv(double theta[],double *f,int *lin) {
  /*
   *====================================================================
   *
   *--LIST (IN IOUT) INITIAL CONDITIONS FOR NEWTON-RAPHSON OR VARIABLE
   *--METRIC METHOD
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static int i,isti;

  /*
   * Format Strings
   */
  const char* F9000 = "\nINITIAL CONDITIONS\n\nPARAMETE\
R                                                      STATUS\
\n                        LOWER BOUND      INITIAL EST      UPPER BOU\
ND\n";
  const char* F9100 = "\n%2d  %28s  %36s\n                   \
%#17.8G%#17.8G%#17.8G\n";
  const char* F9200 = "  INITIAL FUNCTION VALUE%#18.10E\n";
  const char* F9300 = "--------------------------------------------\
--------------------------\n";

  /*
   * Executable Statements
   */

  (void)fprintf(Exp.iout,"%s",F9000);
  for(i=1; i<=Exp.nt; i++) {
    isti = Exp.istin[i - 1];
    (void)fprintf(Exp.iout,F9100,i,Exp.label[i - 1],
		  Xstatus[isti - 1],Exp.thl[i - 1],theta[i - 1],
		  Exp.thu[i - 1]);
  }
  if(*lin > 0 && Exp.lprt > 0) (void)fprintf(Exp.iout,F9200,*f);
  (void)fprintf(Exp.iout,"%s",F9300);
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::litd(double theta[],double *f,int *nfe) {
  /*
   *====================================================================
   *
   *--LIST (IN IDET) DETAILS OF ONE DIRECT-SEARCH ITERATION
   *
   *====================================================================
   */

  /*
   * Local Scalar
   */
  static int i;

  /*
   * Format Strings
   */
  const char* F9000 = "ITERATION%5d\nPAR.     STEP FAC\
T              OLD              NEW           CHANGE\n";
  const char* F9100 = "%2d%#17.8G%#17.8G%#17.8G%#17.4E\n";
  const char* F9200 = "   MAX ABSOLUTE CHANGE%#18.10E\nOLD FUNCTN V\
ALUE%#18.10E  NEW FUNCTN VALUE%#18.10E\nCHANGE          %#18.10E  MA\
X F DIF       %#18.10E\n   NUMBER OF FUNCTION EVALUATIONS%5d\n-------\
---------------------------------------------------------------\n";
    
  /*
   * Executable Statements
   */

  (void)fprintf(Exp.idet,F9000,Xit);
  for(i=1; i<=Exp.nt; i++)
    (void)fprintf(Exp.idet,F9100,i,Xstp[i - 1],Xthpr[i - 1],
		  theta[i - 1],Xcth[i - 1]);
  (void)fprintf(Exp.idet,F9200,Xerm,Xfpr,*f,Xfch,Xdifmax,*nfe);
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::litv(double theta[],double *f,int *nfe) {
  /*
   *====================================================================
   *
   *--LIST (IN IDET) DETAILS OF ONE NEWTON-RAPHSON OR VARIABLE METRIC
   *--ITERATION
   *
   *====================================================================
   */

  /*
   * Local Scalars
   */
  static int i,l;

  /*
   * Format Strings
   */
  const char* F9000 = "ITERATION%5d\n\n   PARAMETER        GRADIEN\
T       SEARCH DIR\n                                OLD              \
NEW           CHANGE\n";
  const char* F9100 = "\n   %2d       %#17.8E%#17.8E\
\n                   %#17.8G%#17.8G%#17.4E\n";
  const char* F9200 = "\n   %2d          --------------   ---------\
-----\n                   %#17.8G%#17.8G%#17.4E\n";
  const char* F9300 = "   NORM OF G%#18.10E  P'G%#18.10E\n   STEP S\
IZE%#18.10E   MAX ABSOLUTE CHANGE%#18.10E\n   OLD FUNCTION VALUE\
%#18.10E\n   NEW FUNCTION VALUE%#18.10E       CHANGE%#18.10E\n   NUMB\
ER OF FUNCTION EVALUATIONS%5d\n--------------------------------------\
--------------------------------\n";

  /*
   * Executable Statements
   */

  (void)fprintf(Exp.idet,F9000,Xit);
  l = 0;
  for(i=1; i<=Exp.nt; i++) {
    if(Xist[i - 1] <= 2) {
      l += 1;
      (void)fprintf(Exp.idet,F9100,i,Xg[l - 1],Xpdir[l - 1],
		    Xthpr[i - 1],theta[i - 1],Xcth[i - 1]);
    }
    else (void)fprintf(Exp.idet,F9200,i,Xthpr[i - 1],theta[i - 1],
		       Xcth[i - 1]);
  }
  (void)fprintf(Exp.idet,F9300,Xgtg,Xptg,Xtstep,Xerm,Xfpr,*f,Xfch,
		*nfe);
  return;
    
}

/*
 * *******************************************************************
 */

void Maxfun_source::lf(double theta[],double *f,int *nfe) {
  /*
   *====================================================================
   *
   *--LIST (IN IOUT) FINAL CONDITIONS
   *
   *====================================================================
   */

  /*
   * Parameter
   */
  const char* Undef = "  ----------    ";

  /*
   * Local Scalars ..
   */
  static double sei;
  static int i,isti,l;
  static char star[2],cg[17],cse[17];

  /*
   * Format Strings
   */
  const char* F9000 = "\nFINAL CONDITIONS\n\nPARAMETE\
R                                                      STATUS\
\n                          FINAL EST          STD DEV      FIRST DER\
IV\n";
  const char* F9100 = "%#16.8G";
  const char* F9200 = "%12.8f    ";
  const char* F9300 = "\n%2d  %28s  %36s\n                   \
%#17.8G %16s%1s%16s\n";
  const char* F9400 = "  FINAL FUNCTION VALUE%#18.10E\n";
  const char* F9500 = "  NUMBER OF FUNCTION EVALUATIONS%5d\n";
  const char* F9600 = " !WARNING:  ITERATION STOPPED AT OR NEAR BOU\
ND OF DEPENDENT PARAMETER!\n !OR IMPLIED BOUNDARY;  FUNCTION NOT NECE\
SSARILY AT A LOCAL MAXIMUM  !\n";
  const char* F9700 = " *NEGATIVE VALUE IN S. D. COLUMN IS COMPUTE\
D \"VARIANCE\"\n";
  const char* F9800 = "\n------------------------------------------\
----------------------------\n---------------------------------------\
-------------------------------\n";

  /*
   * Executable Statements
   */

  (void)fprintf(Exp.iout,"%s",F9000);

  l = 0;
  for(i=1; i<=Exp.nt; i++) {
    (void)strcpy(star," ");
    isti = Xist[i - 1];
    if(isti <= 2) l += 1;
        
    if(Xivfl <= 2) {
      sei = Xstde[i - 1];
      (void)sprintf(cse,F9100,sei);
      if(sei < 0.e0) (void)strcpy(star,"*");
    }
    else (void)strcpy(cse,Undef);

    if(Xigage == 0 && isti <= 2)
      (void)sprintf(cg,F9200,Xg[l - 1]);
    else (void)strcpy(cg,Undef);

    (void)fprintf(Exp.iout,F9300,i,Exp.label[i - 1],
		  Xstatus[isti - 1],theta[i - 1],cse,star,cg);
  }
    
  if(Exp.lprt > 0) (void)fprintf(Exp.iout,F9400,*f);
  (void)fprintf(Exp.iout,F9500,*nfe);

  if(Ximpbnd > 0) (void)fprintf(Exp.iout,"%s",F9600);

  if(Xivfl == 2) (void)fprintf(Exp.iout,"%s",F9700);
  (void)fprintf(Exp.iout,"%s",F9800);
  return;

}

/*
 * *******************************************************************
 */

void Maxfun_source::copyr(FILE *pout) {
  /*
   *====================================================================
   *
   * COPYR VERSION 2.0:
   * SUBROUTINE TO WRITE COPYRIGHT MESSAGE TO SPECIFIED OUTPUT FILE
   *
   * ALEXA J. M. SORANT AND ALEXANDER F. WILSON 1-APR-1994
   *
   *====================================================================
   */

  /*
   * Format String
   */
  const char* F9000 = "COPYRIGHT (C) 1994 BY R. C. ELSTON, INC.\n";

  /*
   * Executable Statements
   */

  (void)fprintf(pout,"%s",F9000);
  return;
    
}
