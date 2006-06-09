/******************************************************************************
File: Calculate.cpp
Description: The main function of Calculate.cpp is to output the results of
             multic.  It outputs estimates and standard errors such as mu and
             mg1 to files named multic.mu and multic.mg1 respectively.
Author: Eric Lunde, 6-12-03 (This file was not writen by Eric, but this is the
        date he took over and began modifying it.)
Updates: (Date, Modified By, Modification Description)
6-12-03, Eric Lunde, I found a chain of an unused variable.  char str_lik[]
         was a local variable of main_fun and was passed to Run, SRun, and
         PrintSummary.  The only time it was given value was in SRun, which it
         received the string representation of another variable, ln_func.  It
         never was used for its contents.
7-23-03, Eric Lunde, One of the problems I encountered in this process was that
         all of a sudden, multic would not compete execution.  This was because
         I removed code in this file that opened, closed, and reset the get
         pointer to a file named share.out because no where else in this file
         was share.out being read.  The actual reading was taking place in
         Likelihood.cpp (Calculate's parent class).  I moved the opening and
         closing of fp_share, fp_data, and fp_lookuplog into this file to
         promote encapsulation.  I did not move the opening and closing of
         fp_loci because, it is opened inside of a conditional statement.  I
         thought that it belonged in Calculate.cpp.
7-29-03, Eric Lunde, I decided to implement representing share.out internally
         as an array of ShareRelation structs.  I added a ShareRelation* and
         int to the constructor of Calculate.  When the program would reset the
         file get pointer, it know resets the relationIndex (position in the
         array) back to 0.
7-30-03, Eric Lunde, I deleted the code that reset the get pointers for the
         file share.out.
8-6-03, Eric Lunde, I've added fortArray and dataSize to the parameter list for
        the constructor.  Those variables are members of Likelihoodfun.  Also
        I've removed Calculate's use of the class Getheader.  We were creating
        an instance of Getheader to acquire the number of families and their
        respective sizes.  Now main_fun gets that information from simple
        algorithms on the fortArray.
8-27-03, Eric Lunde, I added code to make two new files, fam.lik0 and fam.likA.
         These files contain the family loglikelihoods for each family.  The
         structure of the file is that each line contains one hypothesis (null
         or alternative) and each line holds one value per family.  I create
         the array to hold the values in main_fun after calculating the
         variable nfam.  I store values in the array in (NO)AS_CalulateValues.
         The values are then printed from the method SRun after the final
         calulations have been made.
9-09-03, Eric Lunde, We need to save the results (estimates) from the null
         hypothesis to be used as initial values for the alternative
         hypotheses.  This will be done by saving those values in the variable
         double *estimates.  This will contain the address of the
         initial.values parameter that came all the way from the
         .C("multic", ...) call.  I added a variable 'int coeffecientIndex' in
         PrintSummaryP3_Result_1.  Since the coeffecients are stored at the end
         of the array, we need to calculate the appropriate offset and store it
         in coeffecientIndex.
9-29-03, Eric Lunde, At the end of AS_CalulateValues, we were free'ing memory
         that had already been free'd.  This was causing the program to crash.
         It's funny, I've spent to time dealing with AS_CalculateValues, but
         along the development, this code that I deleted, was added in.
10-13-03, Eric Lunde, When reading the file 'summary.log' back into Splus, we
          need to have like types in the columns of summary.log.  This is not
          the case currently as when we write "(FIXED_EXTERNALLY)" or
          "(FIXED_INTERNALLY)".  We need to write to summary.log some value
          that idicates these states.  This value should be decided by Beth
          Atkinson and Mariza de Andrade.  But for now I will use 69 to
          indicate "(FIXED_EXTERNALLY)" and 73 to indicate
          "(FIXED_INTERNALLY)".  The reasoning is that 69 and 73 are the ASCII
          values for 'E' and 'I' respectively.
01-08-04, Eric Lunde, I've added code to print the iteration number in a more
          presentable manner (with commas and new lines).  I've also added
          code to save the V matrix and the (y-beta) vector.  I've created the
          variables and allocation/free methods to manipulate the initial
          variables, methods to save the last iteration's V matrix and (y-beta)
          vector, and methods to write those structures to their appropriate
          files.  Those files are named with the name of the ibd file that
          created that alternative hypothesis's loci.out.
03-18-2004, Eric Lunde, I added code to print the ibd file name to mark each 
            inv.exp.sec.der.random matrix written to file.  Also, the values
            are written in a more "matrix" looking format.  Finally, whenever
            an effect becomes fixed or is constrained, -9's will be printed in 
            their absence to keep the sizes the same for each alternative
            hypothesis.
09-20-2004, Eric Lunde, In the 9-9-03 entry, I talk about saving the
            estimates and coefficients in the same array (hence the need to
            calculate the coefficeintsIndex.  Since then, the estimates and
            coefficients have been split up in to their own separate arrays.
            I neglected to reflect this in how I operated on the estimates
            array (I was writing the coefficients off the end the estimates
            array).  I have created a second coefficients array to serve the
            same purpose as estimates served.  Hopefully this will end multic
            hanging up for no apparent reason just before returning the
            completed multic object.
******************************************************************************/
#include <sys/times.h>
#include <ctime>
#include <limits.h>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include "Calculate.h"
#include "Likelihoodfun.h"
#include "Lib.h"
#include "multic.h"
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <new>
#include <cfloat>
using namespace std;

extern int ascert_flag, run_flag, runopt_flag, need_ibd_flag,
  fixtoboundary_flag, algorithm_flag;
extern double epsilon;

// char *ibdFN was added by Eric Lunde on 7-11-03 to initialize ibdFileName
// ShareRelation *shareArr was added by Eric Lunde on 7-30-03 to initialize
//   the shareArray variable in super-class Likelihoodfun.
// int relationCount was added by Eric Lunde on 7-30-03 to initialize the
//   relationSize variable in super-class Likelihoodfun.
// FortData *fortArr was added by Eric Lunde on 7-31-03 to initialize
//   the fortArray variable in super-class Likelihoodfun.
// int dataCount was added by Eric Lunde on 7-31-03 to initialize the
//   dataSize variable in super-class Likelihoodfun.
// double *est was added by Eric Lunde on 9-09-03
// Sint *familySizes was added by Eric Lunde on 2005-09-07
// char **uniqueFamilies was added by Eric Lunde on 2005-09-07
// int familyCount was added by Eric Lunde on 2005-09-07
Calculate::Calculate(Least *lea, char *ibdFN, double epsilon,
		     ShareRelation *shareArr, int relationCount,
                     FortData *fortArr, int dataCount, double *est,
		     double *coeff, Sint print_progress,
		     Sint calculate_residuals, Sint *familySizes,
		     char **uniqueFamilies, int familyCount) :
  fp_out(out_file, ios::app),
  fp_mu(mu_file, ios::app),
  fp_s(s_file, ios::app),
  fp_m1(m1_file, ios::app),
  fp_m2(m2_file, ios::app),
  fp_t(t_file, ios::app),
  fp_sib(sib_file, ios::app),
  fp_pp(pp_file, ios::app),
  fp_po(po_file, ios::app),
  fp_lik(lik_file, ios::app),
  fp_log(log_file, ios::app),
  summaryLog(summary_file, ios::app)
{
  least = lea;

  LogLikbreakval = epsilon;

  shareArray = shareArr;
  relationSize = relationCount;
  relationIndex = 0;

  fortArray = fortArr;
  dataSize = dataCount;
  dataIndex = 0;

  // ibdFileName was added on 7-11-03 by Eric Lunde so we can print the name of
  // the (m)ibd file that genereated the output
  strncpy(ibdFileName, ibdFN, BUF);

  estimates = est;
  estimatesIndex = 0;
  coefficients = coeff;
  coefficientsIndex = 0;

  // printProgress determines whether or not to print the iteration number,
  // this is to be turned OFF when calling multic to get a base values for
  // calculating the proportion of variance
  printProgress = print_progress;

  // when calculateResiduals is 0, do not print the V matrices nor the y beta
  // diff files.
  calculateResiduals = calculate_residuals;

  // familySizes, uniqueFamilies, and familyCount are being calulated
  // correctly in get.family.sizes, get.unique.families, and get.family.count
  // (S code) so, we just reference that data here.  Eric Lunde, 2005-09-07
  nfam = familyCount;

  N = (int *) malloc(nfam * sizeof(int) );
  for(int i = 0; i < nfam; i++) {
    N[i] = (int) familySizes[i];
  }  

  uniqueFamilyIds = uniqueFamilies;
}

Calculate::~Calculate() {
  free(N);
  free(familyLoglikelihoods);
  CloseFiles();
}

/*  This is the main program which calls
    OpenFiles, PrintTitle, PrintParaInfo,
    GetFlagValues, GetFromPara, GetStartMubeta_SMT,
    GetDim, GetIbd_flag, Run, PrintSummary,
    CloseFiles. functions.

    main_fun:
    1 - Prints copyright information and parameters received from 'multic.par'
    2 - Opens files such as 'multic.log' and 'multic.mu'
    3 - Copies data from the TraitMarkerCov_par and InitValue_par objects into
        it's own local memory
*/
void Calculate::main_fun(TraitMarkerCov_par *tmc, InitValue_par *init){
  int    i;
  int    mubeta_dim, smtcpq_dim;
  int    init_need_ibd_flag;
  clock_t mytime;
  struct tms tbuffer;
  time_t tim;
  struct tm *now;
  double time_secs;
  char line[BUF];

  clock();

  need_ibd_flag= YES;  /* set default G_VAL 'need_ibd_flag' here and
			  it may be changed within this function.   */
  need_mg1_flag = need_mg2_flag = YES;

  // These next three function calls were commented to alter the screen output
  // of multic during execution.  The statements did not describe where the
  // program was in execution and that seemed to be important for a program
  // that could take a Sint time to run. - Eric Lunde 10-20-03
  // PrintTitle(&cout);
  // PrintParaInfo();
  // PrintInputFiles();

  OpenFiles();

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

  SetDataType();

  // added by chen, maybe further change in case least square
  // include c,p,q component in the model.
  if(algorithm_flag==3 && least!=NULL){
    least->getmajorgene(init_x_M1);
    least->getpolygene(init_x_S);
    least->getenvironment(init_x_T);
  }

  // Test to see if there is at least 1 family, and initialize the nfam, N, and
  // familyLoglikelihoods variables.
  if(dataSize <= 0) {
    PROBLEM "There is no family data\nCalculate.cpp line 185"
      RECOVER(NULL_ENTRY);
  }

  familyLoglikelihoods = (double *) malloc( nfam * sizeof(double) );
  for(i=0; i<nfam; i++) {
    familyLoglikelihoods[i] = 0.0;
  }

  savedVMatrix = allocateSavedVMatrix(nfam, total_trait_values, N);
  savedYBetaDiff = allocateSavedYBetaDiff(nfam, total_trait_values, N);
  savedTraitValuedMembers = (int *) malloc(nfam * sizeof(int));

  GetStartMubeta_SMT();
  GetStartCPQ();

  // If the user chose without major gene run (null hypothesis) only.
  if (runopt_flag == 2) {
    m1_flag = FIXED;
    m2_flag = FIXED;
    for (i = 0; i < iinitvcnum; i++){
      x_M1[i]              = 0.0;
      x_M2[i]              = 0.0;
      m1_flag_array[i]     = FIXED;
      m2_flag_array[i]     = FIXED;
      m1_fix_flag_array[i] = EXTERNALLY;
      m2_fix_flag_array[i] = EXTERNALLY;
    }
  }

  /* Get initial G_VAL 'imubeta_dim' and 'ismtcpq_dim' here and they may be
     changed if "runopt == 3" or some parameters are fixed.         */
  //  cout << "GetDim - 1 key 277" << endl;
  GetDim();

  mubeta_dim = imubeta_dim;
  smtcpq_dim = ismtcpq_dim;

  mubeta_A_sd   = Lib::dvector(0,mubeta_dim-1);
  mubeta_0_sd   = Lib::dvector(0,mubeta_dim-1);
  smtcpq_A_sd   = Lib::dvector(0,smtcpq_dim-1);
  smtcpq_0_sd   = Lib::dvector(0,smtcpq_dim-1);

  GetIbd_flag();

  // added one line here;  3/10//97;   E.Y.
  init_need_ibd_flag = need_ibd_flag;

  if (init_need_ibd_flag == YES) {
    fp_loci.open(loci_file);
    if (fp_loci == NULL) {
      PROBLEM "The file %s could not be opened for reading\nCalculate.cpp key 192",
	loci_file RECOVER(NULL_ENTRY);
    }
    fp_loci.getline(line, BUF);
  } else {
    strcpy(line, "# null");
  }
  markerName = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(markerName, &line[2]);
  dataIndex = 0;

  Run();
  // Write the V matrix and (y-beta) vector to file and release the memory.
  // Eric Lunde, 01-09-04
  // This is commented to test whatn happens when we use
  // familyNonMissingSizes instead of N
  if(calculateResiduals) {
    writeSavedVMatrixToFile(savedVMatrix, nfam, total_trait_values,
			    savedTraitValuedMembers,
			    ibdFileName);
    writeSavedYBetaDiffToFile(savedYBetaDiff, nfam, total_trait_values,
			      savedTraitValuedMembers, ibdFileName);
  }


  freeSavedVMatrix(savedVMatrix, nfam, total_trait_values, N);
  freeSavedYBetaDiff(savedYBetaDiff, nfam, total_trait_values);
  free(savedTraitValuedMembers);

  if (init_need_ibd_flag == YES) {
    fp_loci.close();
  }
  free(markerName);

  mytime = times(&tbuffer);

  time(&tim);
  now = localtime(&tim);
  // time_secs=clock()/CLOCKS_PER_SEC;
  time_secs = (tbuffer.tms_utime
	       +tbuffer.tms_stime
	       +tbuffer.tms_cutime
	       +tbuffer.tms_cstime)*1.0/CLOCKS_PER_SEC/*CLK_TCK*/;
  PrintSummary(now,time_secs);

  /* Get CPU time used by the program. 8/14/96.
     system("ps -ef | grep $USER | grep multic >> cputime.log");
  */

  Lib::free_dvector(mubeta_A_sd,   0, mubeta_dim-1);
  Lib::free_dvector(mubeta_0_sd,   0, mubeta_dim-1);
  Lib::free_dvector(smtcpq_A_sd,   0, smtcpq_dim-1);
  Lib::free_dvector(smtcpq_0_sd,   0, smtcpq_dim-1);
}

double Calculate::getmajorgene1(){
  return x_M1[0];
}

double Calculate::getmajorgene2(){
  return x_M2[0];
}

double Calculate::getpolygene(){
  return x_S[0];
}

double Calculate::getenvironment(){
  return x_T[0];
}

double Calculate::getC(){
  return x_c[0];
}

double Calculate::getP(){
  return x_p[0];
}

double Calculate::getQ(){
  return x_q[0];
}


/* Open files fp_para, fp_data, fp_share, fp_out, fp_mu, fp_s,
   fp_m1, fp_m2, fp_t, fp_sib,  fp_pp,    fp_po, fp_lik, fp_log.
*/
void Calculate::OpenFiles(){
  if(fp_out.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 292\n",
      out_file RECOVER(NULL_ENTRY);
  }

  if(fp_mu.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 299\n",
      mu_file RECOVER(NULL_ENTRY);
  }

  if(fp_s.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 306\n",
      s_file RECOVER(NULL_ENTRY);
  }

  if(fp_m1.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 313\n",
      m1_file RECOVER(NULL_ENTRY);
  }

  if(fp_m2.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 321\n",
      m2_file RECOVER(NULL_ENTRY);
  }

  if(fp_t.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 328\n",
      t_file RECOVER(NULL_ENTRY);
  }

  if(fp_sib.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 335\n",
      sib_file RECOVER(NULL_ENTRY);
  }

  if(fp_pp.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 342\n",
      pp_file RECOVER(NULL_ENTRY);
  }

  if(fp_po.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 349\n",
      po_file RECOVER(NULL_ENTRY);
  }

  if(fp_lik.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 356\n",
      lik_file RECOVER(NULL_ENTRY);
  }

  if(fp_log.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 363\n",
      log_file RECOVER(NULL_ENTRY);
  }

  if(summaryLog.fail()) {
    PROBLEM "The file %s could not be opened for appending\nCalculate.cpp key 370\n",
      summary_file RECOVER(NULL_ENTRY);
  }
}

/* Close files fp_para, fp_data, fp_share, fp_out, fp_mu, fp_s,   fp_m,
	       fp_t,    fp_sib,  fp_pp,    fp_po, fp_lik, fp_log.
*/
void Calculate::CloseFiles() {
  fp_log.fill('#');
  fp_log.width(70);
  fp_log << endl << "#" << endl;

  fp_out.close();
  fp_mu.close();
  fp_s.close();
  fp_m1.close();
  fp_m2.close();
  fp_t.close();
  fp_sib.close();
  fp_pp.close();
  fp_po.close();
  fp_lik.close();
  fp_log.close();
  summaryLog.close();
}


/* Print the title, update time and copyright for the program.
*/
void Calculate::PrintTitle(ostream *fp) {
  *fp << endl
      << "\tAnalyze Multivariate Traits with Covariance Program\t(MULTIC)"
      << endl;
  *fp << "\tACT:\tAnalysis for Complex Traits Package" << endl;
  *fp << "\tRevision " << REVISION << " (" << DISTRIBDATE << ")" << endl;
  *fp << "\tCopyright(C) 1997" << endl;
  *fp << "\tDepartment of Epidemiology" << endl;
  *fp << "\tUT M.D. Anderson Cancer Center" << endl;
  *fp << "\tAll rights reserved." << endl << endl << endl;
}

/* Print the paremeters information for the program.
*/
void Calculate::PrintParaInfo() {
  cout << "\t(Convergence value :  " << BREAKVAL << ")" << endl << endl;
  cout << "\tHere are the maximum numbers for your data file," << endl;
  cout << "\tplease check them before you continue to run the program." <<endl;
  cout << "\t Maximum number of traits:               " << MAXNUMTRAIT
       << " (multivariate data)" << endl;
  cout << "\t                                   or    " << 1
       << " (longitudinal data)" << endl;
  cout << "\t Maximum number of unlinked marker loci: " << MAXNUMMARKER <<endl;
  cout << "\t Maximum number of covariates:           " << MAXNUMCOV << endl;
  cout << "\t Maximum number of repeat measurement (longitudinal data only): "
       << MAXREPEATMST << endl;
  cout << "\t Maximum number of alleles at a locus:   " << NUMALLELE << endl;
  cout << "\t Maximum number of each family member:   " << FAMMEMNUM << endl;
}

/****************************
   Print the three input files (and their purpose) to the screen
*****************************/
void Calculate::PrintInputFiles() {
  cout << "\nYou should have the following input files:\n";
  cout << "\t'fort.12'     --- family data file;\n";
  cout << "\t'loci.out'    --- marker information file;\n";
  cout << "\t'share.out'   --- relationship information file.\n";
  cout << "\n\t The program is running, please wait ... \n\n";
}

/***************************
   Set the datatype variable based on irepeatmst
****************************/
void Calculate::SetDataType() {
  if(irepeatmst==1)
    datatype  = MULTIVARIATE;
  if(irepeatmst>1)
    datatype  = LONGITUDINAL;
}


/* Get the flag for running IBD part or not.
*/
void Calculate::GetIbd_flag() {
  int i;
  int count_mg1, count_mg2;

  count_mg1 = count_mg2 = 0;
  for (i = 0; i < iinitvcnum; i++) {
    if (x_M1[i] == 0.0)
      count_mg1++;
    if (x_M2[i] == 0.0)
      count_mg2++;
  }
  if ((count_mg1 == iinitvcnum) && (m1_flag == FIXED)){
    need_mg1_flag = NO;
  }
  if ((count_mg2 == iinitvcnum) && (m2_flag == FIXED)){
    need_mg2_flag = NO;
  }

  if ((need_mg1_flag == NO) && (need_mg2_flag == NO)){
    need_ibd_flag= NO;
  }
}



/* Print the summary in the log file.
*/
void Calculate::PrintSummary(struct tm *now, double time_secs) {
  char *run_date;
  run_date = asctime(now);
  fp_log << "\t\t" << run_date << endl;
  PrintTitle(&fp_log);
  /*
    printf("\t(It is running at %d-%d-%d)\n",now->tm_mon +1,
    now->tm_mday,now->tm_year);
  */

  fp_log << "\tThe program used cpu time: "<< time_secs << "  seconds"
	 << endl << endl << endl;
  fp_log << "\t\t =====================" << endl;
  fp_log << "\t\t| SUMMARY OF ANALYSIS |" << endl;
  fp_log << "\t\t =====================" << endl;

  /* Print the data type here, longitudinal or multivariate data.
     added 12-16-97.
  */
  if (datatype == LONGITUDINAL)
    fp_log << "\t\t ( longitudinal data )" << endl;
  else if (datatype == MULTIVARIATE)
    fp_log << "\t\t ( multivariate data )" << endl;

  PrintSummaryP1();
  PrintSummaryP2();
  if (converg_flag == YES) {    
    PrintSummaryP3(FIRST);
    fp_log << endl << "\t(3). Log Likelihood after convergence:" << endl;
    if (runopt_flag == 1) {
      fp_log << "\t        L(A) =    " << loglikelihood_A << endl;
      summaryLog << setprecision(10)
		 << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_A
		 << setw(SUMMARY_LOG_FIELD_WIDTH) << "converg" << endl;
    }else {
      fp_log << "\t        L(0) =    " << loglikelihood_0 << endl;
      summaryLog << setprecision(10)
		 << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_0
		 << setw(SUMMARY_LOG_FIELD_WIDTH) << "converg" << endl;
    }
  }else {
    if (breakstephalf_run1 == NO) {
      /* Second time use G_VAL 'N1' in the program.       */
      fp_log << endl << endl
	     << "\t***THERE ARE NO CONVERGENT VALUES AFTER " << N1
	     << " ITERATIONS." << endl;
      fp_log << endl << "\t***THE LAST OUTPUTS ARE:";
      PrintSummaryP3(FIRST);
      fp_log << endl << "\t(3). Log Likelihood (non-convergent value):"
	     << endl;
      if (runopt_flag == 1) {
	fp_log << "\t        L(A) (non-convergent)  =    " << loglikelihood_A
	       << endl;
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_A
		   << setw(SUMMARY_LOG_FIELD_WIDTH) << "non-converg" << endl;
      }else {
	fp_log << "\t        L(0) (non-convergent)  =    " << loglikelihood_0
	       << endl;
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_0
		   << setw(SUMMARY_LOG_FIELD_WIDTH) << "non-converg" << endl;
      }
    }else if (breakstephalf_run1 == YES) {
      fp_log << endl << endl << "\t*$*THERE ARE NO CONVERGENT VALUES AFTER "
	     << BREAKSTEPHALF <<" STEP-HALFING." << endl;
      fp_log << endl << "\t***THE LAST OUTPUTS ARE:";
      PrintSummaryP3(FIRST);
      fp_log << endl << "\t(3). Log Likelihood (non-convergent value):"
	     << endl;
      if (runopt_flag == 1) {
	fp_log << "\t        L(A) (non-convergent)  =    " << loglikelihood_A
	       << endl;
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_A
		   << setw(SUMMARY_LOG_FIELD_WIDTH) << "non-converg" << endl;
      }else {
	fp_log << "\t        L(0) (non-convergent)  =    " << loglikelihood_0
	       << endl;
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_0
		   << setw(SUMMARY_LOG_FIELD_WIDTH) << "non-converg" << endl;
      }
    }
  }
  if ((runopt_flag == 3) || (runopt_flag == 4)) {
    if (final_flag == YES) {
      PrintSummaryP3(SECOND);
      fp_log << endl << "\t(3). Log Likelihood under the hypothesis of"
	     << endl;
      fp_log << "\t     with major gene component(s):" << endl;
      fp_log << "\t        L(A) =    " << loglikelihood_A << endl;
      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_A
		 << setw(SUMMARY_LOG_FIELD_WIDTH) << "converg" << endl;

      /* added LRT ; 5/13/97.  */
      fp_log << endl << "\t        LRT = -2*(L(0)-L(A)) =    "
	     << -2*(loglikelihood_0 - loglikelihood_A) << endl;
    }else {
      if (breakstephalf_run2 == NO) {
	fp_log << endl << endl
	       << "\t***THERE ARE NO CONVERGENT VALUES AFTER " << N1
	       << " ITERATIONS." << endl;
	fp_log << endl << "\t***THE LAST OUTPUTS ARE:";
	PrintSummaryP3(SECOND);
	fp_log << endl
	       << "\t(3). Log Likelihood (non-convergent value) under" << endl;
	fp_log << "\t     the hypothesis of with major gene component(s):"
	       << endl;
	fp_log << "\t        L(A) (non-convergent)  =    "
	       << loglikelihood_A << endl;
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_A
		   << setw(SUMMARY_LOG_FIELD_WIDTH) << "non-converg" <<  endl;
      }else if (breakstephalf_run2 == YES) {
	fp_log << endl << endl
	       << "\t*$*THERE ARE NO CONVERGENT VALUES AFTER " << BREAKSTEPHALF
	       << " STEP-HALFING." << endl;
	fp_log << endl << "\t*$*THE LAST OUTPUTS ARE:";
	PrintSummaryP3(SECOND);
	fp_log << endl
	       << "\t(3). Log Likelihood (non-convergent value) under" << endl;
	fp_log << "\t     the hypothesis of with major gene component(s):"
	       << endl;
	fp_log << "\t        L(A) (non-convergent)  =    " << loglikelihood_A
	       << endl;
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << loglikelihood_A
		   << setw(SUMMARY_LOG_FIELD_WIDTH) << "non-converg" <<  endl;
      }

      /* added (non-converg) LRT ; 2/20/98.  */
      fp_log << endl << "\t       (non-converg)LRT = -2*(L(0)-L(A)) =    "
	     << -2*(loglikelihood_0 - loglikelihood_A) << endl;
    }
  }
  //  summaryLog << endl;
}


/* Print the input file names in the log file.
*/
void Calculate::PrintSummaryP1() {
  int i, total_individual;
  // int i, N[MAXNUMFAM], nfam, total_individual;
  char para_str2[BUF], para_str3[BUF];

  dataIndex = 0;
  total_individual = 0;

  sprintf(para_str2, "%d %d %f %d %d", ascert_flag, runopt_flag,
	  epsilon, fixtoboundary_flag, algorithm_flag);

  sprintf(para_str3, "%5d%5d%5d%5d", itraits, iloci, icovs, irepeatmst);

  /* get the total number of families and the total number of
     individuals.  Added 3-3-98; E.Y.  */
  // fscanf(fp_data,"%d", &nfam);  by chen
  for (i = 0; i < nfam; i++){
    // fscanf(fp_data,"%d", &N[i]); by chen
    total_individual += N[i];
  }

  fp_log << endl << endl << "\t";
  int originalFill = fp_log.fill('-');
  int originalWidth = fp_log.width(11);
  fp_log << "-" << endl;
  fp_log.fill(originalFill);
  fp_log.width(originalWidth);

  fp_log << "\tINPUT FILES" << endl;
  fp_log << "\t";
  originalFill = fp_log.fill('-');
  originalWidth = fp_log.width(11);
  fp_log << "-" << endl;
  fp_log.fill(originalFill);
  fp_log.width(originalWidth);

  //    if (datatype == MULTIVARIATE)
  fp_log << "\t(1). Parameter file:   '" << multic_para_file << "'" << endl;
  //    else if (datatype == LONGITUDINAL)
  //        fp_log << "\t(1). Parameter file:   '" << longic_para_file << "'" << end;
  fp_log << "\t     --------------" << endl;
  //  fp_log << "\t     " << para_str1 << endl;
  fp_log << "\t     " << para_str2 << endl;
  fp_log << "\t     " << para_str3 << endl;
  fp_log << "\t     ..." << endl;
  fp_log << "\t     --------------" << endl << endl;
  fp_log << "\t(2). Family data file: '" << data_file << "'" << endl;
  fp_log << "\t     (Total number of families   : " << nfam << endl;
  fp_log << "\t      Total number of individuals: " << total_individual
	 << ")" << endl;
}

/* Print the parameters of the analysis in the log file.
*/
void Calculate::PrintSummaryP2() {
  int i;

  fp_log << endl << endl << "\t";
  int originalFill = fp_log.fill('-');
  int originalWidth = fp_log.width(27);
  fp_log << "-" << endl;
  fp_log.fill(originalFill);
  fp_log.width(originalWidth);

  fp_log << "\tPARAMETERS FOR THE ANALYSIS" << endl;
  fp_log << "\t";
  originalFill = fp_log.fill('-');
  originalWidth = fp_log.width(27);
  fp_log << "-" << endl;
  fp_log.fill(originalFill);
  fp_log.width(originalWidth);

  fp_log << "\t(1). " << total_trait_values << " Trait(s):     ";
  for (i = 0; i < total_trait_values; i++) {
    fp_log << trait_array[i].t_name;
    if ((total_trait_values > 1) && (i != (total_trait_values-1)))
      fp_log << ", ";
  }
  fp_log << "." << endl;

  fp_log << "\t(2). " << iloci << " Marker(s):    ";
  for (i = 0; i < iloci; i++) {
    fp_log << marker_array[i].m_name;
    if ((iloci > 1) && (i != (iloci-1)))
      fp_log << ", ";
  }
  fp_log << "." << endl;

  fp_log << "\t(3). " << icovs << " Covariate(s): ";
  for (i = 0; i < icovs; i++) {
    fp_log << cov_array[i].c_name;
    if ((icovs > 1) && (i != (icovs-1)))
      fp_log << ", ";
  }
  fp_log << "." << endl;

  fp_log << "\t(4). Ascertainment:  ";
  if (ascert_flag == YES)
    fp_log << "YES." << endl;
  else if (ascert_flag == NO)
    fp_log << "NO." << endl;
}

/* Print the results of the analysis in the log file.
*/
void  Calculate::PrintSummaryP3(int H_Aor0_flag) {
  fp_log << endl << endl << "\t";
  int originalFill = fp_log.fill('-');
  int originalWidth = fp_log.width(31);
  fp_log << "-" << endl;
  fp_log.fill(originalFill);
  fp_log.width(originalWidth);

  fp_log << "\tRESULTS FOR THE ANALYSIS  ";

  if ( (H_Aor0_flag == SECOND)
       || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
    fp_log << "(H_A)" << endl;
    summaryLog << "# " << ibdFileName << endl;
  }else if ((H_Aor0_flag == FIRST)
	    &&((runopt_flag == 2) || (runopt_flag == 3) || (runopt_flag == 4))){
    fp_log << "(H_0)" << endl;
    summaryLog << "# null " << endl;
  }

  fp_log << "\t";
  originalFill = fp_log.fill('-');
  originalWidth = fp_log.width(31);
  fp_log << "-" << endl;
  fp_log.fill(originalFill);
  fp_log.width(originalWidth);

  PrintSummaryP3_Result_1(H_Aor0_flag);
  PrintSummaryP3_Result_2(H_Aor0_flag);
}

/* Print the first part of results of the analysis in the log file.
*/
void Calculate::PrintSummaryP3_Result_1(int H_Aor0_flag) {
  int   i,j,kv;
  char  temp_str[STRING+5];

  kv = 0;
  fp_log << "\t(1). Covariate coefficients:" << endl;
  fp_log << "\t                             Estimate             S.E." << endl;
  for (i = 0; i < itraits; i++) {
    fp_log << "\tTrait " << i+1 << " (" << trait_array[i].t_name << "):" << endl;
    for (j = 0; j < icovs+1; j++) {
      if (j == 0) {
	if ((H_Aor0_flag == SECOND)
	    || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	  sprintf(temp_str,"mean(A)%d",i+1);
	  fp_log << setiosflags(ios::left) << "\t        " << setw(13)
		 << temp_str << "=" << setw(16) << setiosflags(ios::right)
		 << x_mu_A[i] << '\t';
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_mu_A[i];
	  fp_log.unsetf(ios::left | ios::right);
	  if (init_mu_flag == ESTIMATE) {
	    fp_log << setw(15) << mubeta_A_sd[kv] << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << mubeta_A_sd[kv]
		       << endl;
	  }else if (init_mu_flag == FIXED) {
	    fp_log << setw(15) << "(Fixed)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << "(Fixed)" << endl;
	  }
	}else if ((H_Aor0_flag == FIRST)
		  &&
		  ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	  sprintf(temp_str,"mean(0)%d",i+1);
	  fp_log << setiosflags(ios::left) << "\t        " << setw(13)
		 << temp_str << "=" << setw(16) << setiosflags(ios::right)
		 << x_mu_0[i] << '\t';
	  fp_log.unsetf(ios::left | ios::right);
	  summaryLog << setw(15) << trait_array[i].t_name;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_mu_0[i];
	  // esitmates is the array passed from multic.s, it is used to give
	  // the alternative hypothesis starting values, Eric Lunde, 9-09-03
	  estimates[estimatesIndex++] = x_mu_0[i];
	  if (init_mu_flag == ESTIMATE) {
	    fp_log << setw(15) << mubeta_0_sd[kv] << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << mubeta_0_sd[kv]
		       << endl;
	  }else if (init_mu_flag == FIXED) {
	    fp_log << setw(15) << "(Fixed)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << "(Fixed)" << endl;
	  }
	}
      }
      else {
	if ((H_Aor0_flag == SECOND)
	    || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	  sprintf(temp_str,"%s(A)%d%d", cov_array[j-1].c_name,i+1,j);
	  fp_log << setiosflags(ios::left) << "\t        " << setw(13)
		 << temp_str << "=" << setw(16) << setiosflags(ios::right)
		 << x_beta_A[i][j-1] << '\t';
	  fp_log.unsetf(ios::left | ios::right);
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_beta_A[i][j-1];
	  if (init_mu_flag == ESTIMATE) {
	    fp_log << setw(15) << mubeta_A_sd[kv] << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << mubeta_A_sd[kv]
		       << endl;
	  }else if (init_mu_flag == FIXED) {
	    fp_log << setw(15) << "(Fixed)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << "(Fixed)" << endl;
	  }
	}
	else if ((H_Aor0_flag == FIRST)
		 && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	  sprintf(temp_str,"%s(0)%d%d", cov_array[j-1].c_name,i+1,j);
	  fp_log << setiosflags(ios::left) << "\t        " << setw(13)
		 << temp_str << "=" << setw(16) << setiosflags(ios::right)
		 << x_beta_0[i][j-1] << '\t';
	  fp_log.unsetf(ios::left | ios::right);
	  summaryLog << setw(15) << cov_array[j-1].c_name;
	  // This if statement is used to only output two numeric values when
	  // those values are different
	  /*
	    This comment is to make the output for the fixed effects show two
	    numbers.  Earlier, I made the program only print the number 2 for
            covariance between the second effect and itself.  That is not the
	    same as writing 22 here.  Eric Lunde 01-16-04.
	  if(i+1 != j) {
	    summaryLog << j;
	  }
	  */
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_beta_0[i][j-1];
	  coefficients[coefficientsIndex] = x_beta_0[i][j-1];
	  coefficientsIndex++;
	  if (init_mu_flag == ESTIMATE) {
	    fp_log << setw(15) << mubeta_0_sd[kv] << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << mubeta_0_sd[kv]
		       << endl;
	  }else if (init_mu_flag == FIXED) {
	    fp_log << setw(15) << "(Fixed)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << "(Fixed)" << endl;
	  }
	}
      }
      kv++;
    }
    fp_log << endl;
  }
  fp_log << endl;
}

/* Print the second part of  results of the analysis in the log file.
*/
void Calculate::PrintSummaryP3_Result_2(int H_Aor0_flag) {
  int   i,j,kv,ks;

  kv = ks = 0;
  fp_log << "\t(2). Variance components:" << endl;
  fp_log << "\t                             Estimate             S.E." << endl;
  fp_log << "\t\tPolygenic:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if ((H_Aor0_flag == SECOND)
	  || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	fp_log << "\t        s(A)" << i+1 << j+1 << "       = " << setw(15)
	       << x_S_A[ks] << setw(0) << "\t";
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_S_A[ks];
	if (s_A_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		     << endl;
	  kv++;
	}else if (s_A_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (s_A_flag_array[ks] == FIXED) {
	  if (s_A_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (s_A_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      else if ((H_Aor0_flag == FIRST)
	       && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	fp_log << "\t        s(0)" << i+1 << j+1 << "       = " << setw(15)
	       << x_S_0[ks] << setw(0) << "\t";
	summaryLog << setw(15) << "s" << i+1;
	// This if statement is used to only output two numeric values when
	// those values are different
	if(i != j) {
	  summaryLog << j+1;
	}
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_S_0[ks];
	estimates[estimatesIndex++] = x_S_0[ks];
	if (s_0_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		     << endl;
	  kv++;
	}else if (s_0_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (s_0_flag_array[ks] == FIXED) {
	  if (s_0_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (s_0_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      ks++;
    }
  }
  fp_log << endl;

  ks = 0;
  fp_log << "\t\tFirst Major gene:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if ((H_Aor0_flag == SECOND)
	  || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	fp_log << "\t        m1(A)" << i+1 << j+1 << "      = " << setw(15)
	       << x_M1_A[ks] << "\t";
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_M1_A[ks];
	if (m1_A_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		     << endl;
	  kv++;
	}else if (m1_A_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (m1_A_flag_array[ks] == FIXED) {
	  if (m1_A_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (m1_A_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      else if ((H_Aor0_flag == FIRST)
	       && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	fp_log << "\t        m1(0)" << i+1 << j+1 << "      = " << setw(15)
	       << x_M1_0[ks] << setw(0) << "\t";
	summaryLog << setw(15) << "mg" << i+1;
	// This if statement is used to only output two numeric values when
	// those values are different
	if(i != j) {
	  summaryLog << j+1;
	}
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_M1_0[ks];
	estimates[estimatesIndex++] = x_M1_0[ks];
	if (m1_0_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		     << endl;
	  kv++;
	}else if (m1_0_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (m1_0_flag_array[ks] == FIXED) {
	  if (m1_0_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (m1_0_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      ks++;
    }
  }
  fp_log << endl;

  ks = 0;
  fp_log << "\t\tSecond Major gene:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if ((H_Aor0_flag == SECOND)
	  || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	fp_log << "\t        m2(A)" << i+1 << j+1 << "      = " << setw(15)
	       << x_M2_A[ks] << setw(0) << "\t";
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_M2_A[ks];
	if (m2_A_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		     << endl;
	  kv++;
	}else if (m2_A_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (m2_A_flag_array[ks] == FIXED) {
	  if (m2_A_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (m2_A_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      else if ((H_Aor0_flag == FIRST)
	       && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	fp_log << "\t        m2(0)" << i+1 << j+1 << "      = " << setw(15)
	       << x_M2_0[ks] << setw(0) << "\t";
	summaryLog << setw(15) << "m2" << i+1;
	// This if statement is used to only output two numeric values when
	// those values are different
	if(i != j) {
	  summaryLog << j+1;
	}
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_M2_0[ks];
	estimates[estimatesIndex++] = x_M2_0[ks];
	if (m2_0_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		     << endl;
	  kv++;
	}else if (m2_0_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (m2_0_flag_array[ks] == FIXED) {
	  if (m2_0_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (m2_0_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      ks++;
    }
  }
  fp_log << endl;

  ks = 0;
  fp_log << "\t\tEnvironment:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
/*
      if (datatype == MULTIVARIATE) {
	if (j == i) {
	  if ((H_Aor0_flag == SECOND)
	      || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	    fp_log << "\t        t(A)" << i+1 << j+1 << "       = " << setw(15)
		   << x_T_A[i] << setw(0) << "\t";
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_T_A[i];
	    if (t_A_flag_array[i] == ESTIMATE) {
	      fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
			 << endl;
	      kv++;
	    }else if (t_A_flag_array[i] == FIXED) {
	      if (t_A_fix_flag_array[ks] == INTERNALLY) {
		fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
		summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			   << (int)'I' << endl;
	      }else if (t_A_fix_flag_array[ks] == EXTERNALLY) {
		fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
		summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			   << (int)'E' << endl;
	      }
	    }
	  }else if ((H_Aor0_flag == FIRST)
		   && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	    fp_log << "\t        t(0)" << i+1 << j+1 << "       = " << setw(15)
		   << x_T_0[i] << setw(0) << "\t";
	    summaryLog << setw(15) << "e" << i+1;
	    // This if statement is used to only output two numeric values when
	    // those values are different
	    if(i != j) {
	      summaryLog << j+1;
	    }
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_T_0[i];
	    estimates[estimatesIndex++] = x_T_0[i];
	    if (t_0_flag_array[i] == ESTIMATE) {
	      fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
			 << endl;
	      kv++;
	    }else if (t_0_flag_array[i] == FIXED) {
	      if (t_0_fix_flag_array[ks] == INTERNALLY) {
		fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
		summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			   << (int)'I' << endl;
	      }else if (t_0_fix_flag_array[ks] == EXTERNALLY) {
		fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
		summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			   << (int)'E' << endl;
	      }
	    }
	  }
	  ks++;
	  }
      }else if (datatype == LONGITUDINAL) {
*/
	if ((H_Aor0_flag == SECOND)
	    || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	  fp_log << "\t        t(A)" << i+1 << j+1 << "       = " << setw(15)
		 << x_T_A[ks] << '\t';
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_T_A[ks];
	  if (t_A_flag_array[ks] == ESTIMATE) {
	    fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		       << endl;
	    kv++;
	  }else if (t_A_flag_array[ks] == FIXED) {
	    if (t_A_fix_flag_array[ks] == INTERNALLY) {
	      fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			 << (int)'I' << endl;
	    }else if (t_A_fix_flag_array[ks] == EXTERNALLY) {
	      fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			 << (int)'E' << endl;
	    }
	  }
	}else if ((H_Aor0_flag == FIRST)
		 && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	  fp_log << "\t        t(0)" << i+1 << j+1 << "       = " << setw(15)
		 << x_T_0[ks] << '\t';
	  summaryLog << setw(15) << "e" << i+1;
	  // This if statement is used to only output two numeric values when
	  // those values are different
	  if(i != j) {
	    summaryLog << j+1;
	  }
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_T_0[ks];
	  estimates[estimatesIndex++] = x_T_0[ks];
	  if (t_0_flag_array[ks] == ESTIMATE) {
	    fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		       << endl;
	    kv++;
	  }else if (t_0_flag_array[ks] == FIXED) {
	    if (t_0_fix_flag_array[ks] == INTERNALLY) {
	      fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			 << (int)'I' << endl;
	    }else if (t_0_fix_flag_array[ks] == EXTERNALLY) {
	      fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	      summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
			 << (int)'E' << endl;
	    }
	  }
	}
	ks++;
//      }
    }
  }
  fp_log << endl;

  fp_log << "\t(3). Shared common environmental variance components:" << endl;
  ks = 0;
  fp_log << "\t                             Estimate             S.E." << endl;
  fp_log << "\t\tShared Sibship:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if ((H_Aor0_flag == SECOND) || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	fp_log << "\t        sib(A)" << i+1 << j+1 << "     = " << setw(15)
	       << x_c_A[ks] << setw(0) << "\t";
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_c_A[ks];
	if (c_A_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		     << endl;
	  kv++;
	}else if (c_A_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (c_A_flag_array[ks] == FIXED) {
	  if (c_A_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (c_A_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      else if ((H_Aor0_flag == FIRST)
	       && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	fp_log << "\t        sib(0)" << i+1 << j+1 << "     = " << setw(15)
	       << x_c_0[ks] << setw(0) << "\t";
	summaryLog << setw(15) << "sib" << i+1;
	// This if statement is used to only output two numeric values when
	// those values are different
	if(i != j) {
	  summaryLog << j+1;
	}
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_c_0[ks];
	estimates[estimatesIndex++] = x_c_0[ks];
	if (c_0_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		     << endl;
	  kv++;
	}else if (c_0_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (c_0_flag_array[ks] == FIXED) {
	  if (c_0_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (c_0_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      ks++;
    }
  }
  fp_log << endl;

  ks = 0;
  fp_log << "\t\tShared Spouse:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if ((H_Aor0_flag == SECOND)
	  || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	fp_log << "\t        p(A)" << i+1 << j+1 << "       = " << setw(15)
	       << x_p_A[ks] << setw(0) << "\t";
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_p_A[ks];
	if (p_A_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		     << endl;
	  kv++;
	}else if (p_A_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (p_A_flag_array[ks] == FIXED) {
	  if (p_A_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
		       << (int)'I' << endl;
	  }else if (p_A_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH)
		       << (int)'E' << endl;
	  }
	}
      }else if ((H_Aor0_flag == FIRST)
	       && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	fp_log << "\t        p(0)" << i+1 << j+1 << "       = " << setw(15)
	       << x_p_0[ks] << setw(0) << "\t";
	summaryLog << setw(15) << "p" << i+1;
	// This if statement is used to only output two numeric values when
	// those values are different
	if(i != j) {
	  summaryLog << j+1;
	}
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_p_0[ks];
	estimates[estimatesIndex++] = x_p_0[ks];
	if (p_0_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		     << endl;
	  kv++;
	}else if (p_0_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (p_0_flag_array[ks] == FIXED) {
	  if (p_0_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (p_0_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      ks++;
    }
  }
  fp_log << endl;

  ks = 0;
  fp_log << "\t\tShared Parent-Offspring:" << endl;
  for (i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if ((H_Aor0_flag == SECOND)
	  || ((H_Aor0_flag == FIRST) && (runopt_flag == 1))) {
	fp_log << "\t        q(A)" << i+1 << j+1 << "       = " << setw(15)
	       << x_q_A[ks] << setw(0) << "\t";
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_q_A[ks];
	if (q_A_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_A_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_A_sd[kv]
		     << endl;
	  kv++;
	}else if (q_A_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (q_A_flag_array[ks] == FIXED) {
	  if (q_A_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (q_A_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }else if ((H_Aor0_flag == FIRST)
	       && ((runopt_flag==2) || (runopt_flag==3) || (runopt_flag==4))) {
	fp_log << "\t        q(0)" << i+1 << j+1 << "       = " << setw(15)
	       << x_q_0[ks] << setw(0) << "\t";
	summaryLog << setw(15) << "q" << i+1;
	// This if statement is used to only output two numeric values when
	// those values are different
	if(i != j) {
	  summaryLog << j+1;
	}
	summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << x_q_0[ks];
	estimates[estimatesIndex++] = x_q_0[ks];
	if (q_0_flag_array[ks] == ESTIMATE) {
	  fp_log << setw(15) << smtcpq_0_sd[kv] << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << smtcpq_0_sd[kv]
		     << endl;
	  kv++;
	}else if (q_0_flag_array[ks] == CONSTRAINT) {
	  fp_log << setw(15) << "(Constraint)" << endl;
	  summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'C'
		     << endl;
	}else if (q_0_flag_array[ks] == FIXED) {
	  if (q_0_fix_flag_array[ks] == INTERNALLY) {
	    fp_log << setw(15) << "(FIXED INTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'I'
		       << endl;
	  }else if (q_0_fix_flag_array[ks] == EXTERNALLY) {
	    fp_log << setw(15) << "(FIXED EXTERNALLY)" << endl;
	    summaryLog << setw(SUMMARY_LOG_FIELD_WIDTH) << (int)'E'
		       << endl;
	  }
	}
      }
      ks++;
    }
  }
  fp_log << endl;
}

double Calculate::ConvergCriterion1(double prev_val) {
  double converg1_val;

  converg1_val = EPSILON1 *(fabs(prev_val) + EPSILON2);
  return converg1_val;
}

/* Fix all the parameter(s) to the BOUNDARY if the it(they) close to the BOUNDARY.
   And save the intermediate values to lookuplog file.
*/
int Calculate::Fixtoboundary() {
  int    i,j,k, jv;
  int    reachboundary_flag;
  double temp_mean;
  double temp_se;

  reachboundary_flag = NO;

  /* Check all ESTIMATE parameters using BOUNDARY.   */

  k = total_trait_values;
  for ( i = 0,jv = 0; i < iinitvcnum; i++) {
    if ((i == 0) || (i == jv+k)) {
      if (s_flag_array[i] == ESTIMATE) {
	if (fabs(x_S[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+i;
	    temp_mean = x_S[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }else {
	    j = i;
	    temp_mean = x_S[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      s(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      s(A)" << i << "     S.E." << endl;
	  }else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      s(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      s(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " << temp_se
		       << endl;

	  x_S[i]              = BOUNDARY;
	  s_flag_array[i]     = FIXED;
	  s_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag  = YES;
	}
      }

      if (m1_flag_array[i] == ESTIMATE) {
	if (fabs(x_M1[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+s_vcnum+i;
	    temp_mean = x_M1[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }else {
	    j = s_vcnum+i;
	    temp_mean = x_M1[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      m1(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      m1(A)" << i << "     S.E." << endl;
	  }else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      m1(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      m1(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " << temp_se
		       << endl;

	  x_M1[i]              = BOUNDARY;
	  m1_flag_array[i]     = FIXED;
	  m1_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag   = YES;
	}
      }

      if (m2_flag_array[i] == ESTIMATE) {
	if (fabs(x_M2[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+s_vcnum+m1_vcnum+i;
	    temp_mean = x_M2[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }else {
	    j = s_vcnum+m1_vcnum+i;
	    temp_mean = x_M2[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      m2(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      m2(A)" << i << "     S.E." << endl;
	  }else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      m2(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      m2(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " << temp_se
		       << endl;

	  x_M2[i]              = BOUNDARY;
	  m2_flag_array[i]     = FIXED;
	  m2_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag   = YES;
	}
      }
      jv = i;
      if ( i != 0)
	k--;
    }
  }

/*
  if (datatype == MULTIVARIATE) {
    for (i = 0; i < env_vcnum; i++) {
      if (t_flag_array[i] == ESTIMATE) {	
	//cout << "x_T["<<i<<"]  = " << x_T[i] << endl;
	if (fabs(x_T[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+s_vcnum+m1_vcnum+m2_vcnum+i;
	    temp_mean = x_T[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }
	  else {
	    j = s_vcnum+m1_vcnum+m2_vcnum+i;
	    temp_mean = x_T[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      t(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      t(A)" << i << "     S.E." << endl;
	  }
	  else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      t(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      t(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " << temp_se
		       << endl;

	  x_T[i]              = BOUNDARY;
	  t_flag_array[i]     = FIXED;
	  t_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag  = YES;
	}
      }
    }
  }else if (datatype == LONGITUDINAL) {
*/
    k = total_trait_values;
    for ( i = 0,jv = 0; i < env_vcnum; i++) {
      if ((i == 0) || (i == jv+k)) {
	if (t_flag_array[i] == ESTIMATE) {
	  if (fabs(x_T[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	    if (ascert_flag == YES) {
	      j = imubeta_dim+s_vcnum+m1_vcnum+m2_vcnum+i;
	      temp_mean = x_T[i];
	      temp_se   = sqrt(inv_Exsder_mat[j][j]);
	    }
	    else {
	      j = s_vcnum+m1_vcnum+m2_vcnum+i;
	      temp_mean = x_T[i];
	      temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	    }
	    if (run_flag == FIRST) {
	      if ((runopt_flag == 3) || (runopt_flag == 4))
		fp_lookuplog << "      t(0)" << i << "     S.E." << endl;
	      else
		fp_lookuplog << "      t(A)" << i << "     S.E.xz" << endl;
	    }
	    else {
	      if ((runopt_flag == 3) || (runopt_flag == 4))
		fp_lookuplog << "      t(A)" << i << "     S.E." << endl;
	      else
		fp_lookuplog << "      t(0)" << i << "     S.E." << endl;
	    }

	    fp_lookuplog << "        " << temp_mean << "       " << temp_se
			 << endl;

	    x_T[i]              = BOUNDARY;
	    t_flag_array[i]     = FIXED;
	    t_fix_flag_array[i] = INTERNALLY;
	    reachboundary_flag  = YES;
	  }
	}
	jv = i;
	if ( i != 0)
	  k--;
      }
    }
    //  }

  k = total_trait_values;
  for ( i = 0,jv = 0; i < iinitvcnum; i++) {
    if ((i == 0) || (i == jv+k)) {
      if (c_flag_array[i] == ESTIMATE) {
	if (fabs(x_c[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+s_vcnum+m1_vcnum+m2_vcnum+t_vcnum+i;
	    temp_mean = x_c[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }else {
	    j = s_vcnum+m1_vcnum+m2_vcnum+t_vcnum+i;
	    temp_mean = x_c[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      c(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      c(A)" << i << "     S.E." << endl;
	  }else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      c(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      c(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " << temp_se
		       << endl;

	  x_c[i]              = BOUNDARY;
	  c_flag_array[i]     = FIXED;
	  c_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag  = YES;
	}
      }

      if (p_flag_array[i] == ESTIMATE) {
	if (fabs(x_p[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+s_vcnum+m1_vcnum+m2_vcnum+t_vcnum+c_vcnum+i;
	    temp_mean = x_p[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }else {
	    j = s_vcnum+m1_vcnum+m2_vcnum+t_vcnum+c_vcnum+i;
	    temp_mean = x_p[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      p(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      p(A)" << i << "     S.E." << endl;
	  }else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      p(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      p(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " <<  temp_se
		       << endl;

	  x_p[i]              = BOUNDARY;
	  p_flag_array[i]     = FIXED;
	  p_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag  = YES;
	}
      }

      if (q_flag_array[i] == ESTIMATE) {
	if (fabs(x_q[i] - BOUNDARY) <= CLOSETOBOUNDARY) {
	  if (ascert_flag == YES) {
	    j = imubeta_dim+s_vcnum+m1_vcnum+m2_vcnum+t_vcnum+c_vcnum+p_vcnum+i;
	    temp_mean = x_q[i];
	    temp_se   = sqrt(inv_Exsder_mat[j][j]);
	  }else {
	    j = s_vcnum+m1_vcnum+m2_vcnum+t_vcnum+c_vcnum+p_vcnum+i;
	    temp_mean = x_q[i];
	    temp_se   = sqrt(inv_Exsder_smtcpq_mat[j][j]);
	  }
	  if (run_flag == FIRST) {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      q(0)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      q(A)%d     S.E." << endl;
	  }else {
	    if ((runopt_flag == 3) || (runopt_flag == 4))
	      fp_lookuplog << "      q(A)" << i << "     S.E." << endl;
	    else
	      fp_lookuplog << "      q(0)" << i << "     S.E." << endl;
	  }

	  fp_lookuplog << "        " << temp_mean << "       " << temp_se
		       << endl;

	  x_q[i]              = BOUNDARY;
	  q_flag_array[i]     = FIXED;
	  q_fix_flag_array[i] = INTERNALLY;
	  reachboundary_flag  = YES;
	}
      }

      jv = i;
      if ( i != 0)
	k--;
    }
  }
  if (reachboundary_flag == YES) {
    FixCovAndAssignInitValues();
    mu_flag = init_mu_flag;
  }
  return reachboundary_flag;
}

/*****************************
Class name: Calculate
Method name: CheckBoundary
Description:  If user chose 'runopt_flag = 4' and 'fixtoboundary_flag = y'
              (which means the program will fix the parameter to the
              BOUNDARY(0) internally if the parameter is close to the
              BOUNDARY), then the parameter has to be re-initiate to the
              original values and flag according to the parameter file.
              changed 2-4-98; E.Y.
	      Check all ESTIMATE parameters and see if they are fixed to
	      the BOUNDARY internally or not.
Input: NONE.
Output: NONE.
Side Effects: NONE.
Author: Eric Lunde, 6-12-03
******************************/
void Calculate::CheckBoundary() {
  int    i;

  for ( i = 0; i < iinitvcnum; i++) {
    if ( s_fix_flag_array[i]  == INTERNALLY )
      x_S[i]    = init_x_S[i];

    if ( m1_fix_flag_array[i] == INTERNALLY )
      x_M1[i]   = init_x_M1[i];

    if ( m2_fix_flag_array[i] == INTERNALLY )
      x_M2[i]   = init_x_M2[i];
  }

  for (i = 0; i < env_vcnum; i++) {
    if ( t_fix_flag_array[i] == INTERNALLY )
      x_T[i]    = init_x_T[i];
  }

  for ( i = 0; i < iinitvcnum; i++) {
    if ( c_fix_flag_array[i] == INTERNALLY )
      x_c[i]    = init_x_c[i];

    if ( p_fix_flag_array[i] == INTERNALLY )
      x_p[i]    = init_x_p[i];

    if ( q_fix_flag_array[i] == INTERNALLY )
      x_q[i]    = init_x_q[i];
  }
}


/* Calculate the difference of the parameters and assign the value to
   converg1_flag which is the 1st criterion.
*/
void Calculate::CalculateDiffOfPara(double *prev_mu, double **prev_beta,
				    double *prev_s, double *prev_m1,
				    double *prev_m2, double *prev_t,
				    double *prev_c,  double *prev_p,
				    double *prev_q, int *converg1_flag,
				    int stephalf_flag)
{
  int    i,j;
  double sum_mubeta;
  double *diff_s, *diff_m1, *diff_m2, *diff_t;
  double *diff_mu,**diff_beta;
  double *diff_c, *diff_p,  *diff_q;

  /* Initialize the size of the pointers.         changed 5-12-98; E.Y.   */
  diff_s    = (double *) malloc(iinitvcnum * sizeof(double));
  diff_m1   = (double *) malloc(iinitvcnum * sizeof(double));
  diff_m2   = (double *) malloc(iinitvcnum * sizeof(double));
  diff_t    = (double *) malloc(env_vcnum * sizeof(double));
  diff_mu   = (double *) malloc(itraits * sizeof(double));
  diff_beta = (double **)malloc(itraits    * sizeof(double *));
  for (i = 0; i < itraits; i++)
    diff_beta[i] = (double *) malloc(icovs * sizeof(double));
  diff_c    = (double *) malloc(iinitvcnum * sizeof(double));
  diff_p    = (double *) malloc(iinitvcnum * sizeof(double));
  diff_q    = (double *) malloc(iinitvcnum * sizeof(double));

  *converg1_flag = YES;       /* set it to YES, if one of the parameters
				 doesn't satisfy the 1st criterion then
				 set this to NO.                      */

  /* Fix mu and beta values if the sum_mubeta  < BREAKVAL .           */
  sum_mubeta = 0.0;
  if (mu_flag == ESTIMATE) {
    for (i = 0; i < itraits; i++) {
      diff_mu[i]  = fabs(x_mu[i]-prev_mu[i]);
      sum_mubeta += diff_mu[i];
      if (diff_mu[i] > ConvergCriterion1(prev_mu[i]))
	*converg1_flag = NO;

      for (j = 0; j < icovs; j++) {
	diff_beta[i][j] = fabs(x_beta[i][j] - prev_beta[i][j]) ;
	sum_mubeta     += diff_beta[i][j] ;
	if (diff_beta[i][j] > ConvergCriterion1(prev_beta[i][j]))
	  *converg1_flag = NO;
      }
    }
    //    cout << "(sum_mubeta < BREAKVAL) = " << (sum_mubeta < BREAKVAL) << " Calc.cpp key 1887" << endl;
    //    cout << "(stephalf_flag == 0) = " << (stephalf_flag == 0) << " Calc.cpp key 1887" << endl;
    if ((sum_mubeta < BREAKVAL) && (stephalf_flag == 0)) {
      /* Change G_VAL 'mu_flag' here.              */
#ifdef DEBUG
      cout << "sum_mubeta < BREAKVAL" << endl;
#endif
      mu_flag = FIXED;
      for (i = 0; i < imubeta_dim; i++) {
	if (ascert_flag == YES) {
	  if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1)))
	    mubeta_A_sd[i] = sqrt(inv_Exsder_mat[i][i]);
	  else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
					   || (runopt_flag == 4)))
	    mubeta_0_sd[i] = sqrt(inv_Exsder_mat[i][i]);
	}else {
	  if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1)))
	    mubeta_A_sd[i] = sqrt(inv_Exsder_mubeta_mat[i][i]);
	  else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
					   || (runopt_flag == 4)))
	    mubeta_0_sd[i] = sqrt(inv_Exsder_mubeta_mat[i][i]);
	}
      }
    }
  }

  /* check x_S, x_M1, x_M2, x_T, x_c, x_p, x_q convergent criterion.  */
  for (i = 0; i < env_vcnum; i++) {
    if (t_flag_array[i] == ESTIMATE) {
      diff_t[i] = fabs(x_T[i] - prev_t[i]);
      if (diff_t[i] > ConvergCriterion1(prev_t[i])) {
	*converg1_flag = NO;
      }
    }
  }

  for ( i = 0; i < iinitvcnum; i++) {
    if (s_flag_array[i] == ESTIMATE) {
      diff_s[i] = fabs(x_S[i] - prev_s[i]);
      if (diff_s[i] > ConvergCriterion1(prev_s[i]))
	*converg1_flag = NO;
    }
    if (m1_flag_array[i] == ESTIMATE) {
      diff_m1[i] = fabs(x_M1[i] - prev_m1[i]);
      if (diff_m1[i] > ConvergCriterion1(prev_m1[i]))
	*converg1_flag = NO;
    }
    if (m2_flag_array[i] == ESTIMATE) {
      diff_m2[i] = fabs(x_M2[i] - prev_m2[i]);
      if (diff_m2[i] > ConvergCriterion1(prev_m2[i]))
	*converg1_flag = NO;
    }
    if (c_flag_array[i] == ESTIMATE) {
      diff_c[i] = fabs(x_c[i] - prev_c[i]);
      if (diff_c[i] > ConvergCriterion1(prev_c[i]))
	*converg1_flag = NO;
    }
    if (p_flag_array[i] == ESTIMATE) {
      diff_p[i] = fabs(x_p[i] - prev_p[i]);
      if (diff_p[i] > ConvergCriterion1(prev_p[i]))
	*converg1_flag = NO;
    }
    if (q_flag_array[i] == ESTIMATE) {
      diff_q[i] = fabs(x_q[i] - prev_q[i]);
      if (diff_q[i] > ConvergCriterion1(prev_q[i]))
	*converg1_flag = NO;
    }
  }
  free(diff_s);
  free(diff_m1);
  free(diff_m2);
  free(diff_t);
  free(diff_mu);
  for (i = 0; i < itraits; i++)
    free(diff_beta[i]);
  free(diff_beta);
  free(diff_c);
  free(diff_p);
  free(diff_q);
}

/*
   Get the new MU,SIGMA,MG,TAU and the loglikelihood.
*/
double Calculate::CalculateLik(double *prev_mu, double **prev_beta,
			       double *prev_s, double *prev_m1,
			       double *prev_m2, double *prev_t,
			       double *prev_c, double *prev_p,
			       double *prev_q, double prev_ln_func,
			       int *converg1_2_flag) {

  int    i,j;
  int    stephalf_flag;
  // int    nfam;
  // int    N[MAXNUMFAM];
  int    converg1_flag = YES, converg2_flag = YES;
  int    reachboundary_flag;
  double denom_half, stephalf_val;
  double ln_func = 0.0;
  double *temp_inc_mu,   *temp_inc_beta;
  double *temp_inc_s, *temp_inc_m1, *temp_inc_m2, *temp_inc_t;
  double *temp_inc_c, *temp_inc_p,  *temp_inc_q;

  /* Initialize the size of the pointers.         changed 5-12-98; E.Y.   */
  temp_inc_mu   = (double *) malloc(itraits         * sizeof(double));
  temp_inc_beta = (double *) malloc((itraits*icovs) * sizeof(double));
  temp_inc_s    = (double *) malloc(iinitvcnum      * sizeof(double));
  temp_inc_m1   = (double *) malloc(iinitvcnum      * sizeof(double));
  temp_inc_m2   = (double *) malloc(iinitvcnum      * sizeof(double));
  temp_inc_t    = (double *) malloc(env_vcnum       * sizeof(double));
  temp_inc_c    = (double *) malloc(iinitvcnum      * sizeof(double));
  temp_inc_p    = (double *) malloc(iinitvcnum      * sizeof(double));
  temp_inc_q    = (double *) malloc(iinitvcnum      * sizeof(double));

  /* Save all the previous increment parameters.                    */
  for (i = 0; i < itraits; i++) {
    temp_inc_mu[i] = new_inc_mu[i];
    inc_mu[i]      = new_inc_mu[i];
  }
  for (i = 0; i < env_vcnum; i++) {
    temp_inc_t[i]  = new_inc_t[i];
    inc_t[i]       = new_inc_t[i];
  }
  for (i = 0; i < iinitvcnum; i++) {
    temp_inc_s[i]  = new_inc_s[i];
    inc_s[i]       = new_inc_s[i];
    //cout << "inc_s["<<i<<"] = " << inc_s[i] << endl;
    temp_inc_m1[i] = new_inc_m1[i];
    inc_m1[i]      = new_inc_m1[i];
    temp_inc_m2[i] = new_inc_m2[i];
    inc_m2[i]      = new_inc_m2[i];
    temp_inc_c[i]  = new_inc_c[i];
    inc_c[i]       = new_inc_c[i];
    temp_inc_p[i]  = new_inc_p[i];
    inc_p[i]       = new_inc_p[i];
    temp_inc_q[i]  = new_inc_q[i];
    inc_q[i]       = new_inc_q[i];
  }
  for (i = 0; i < (itraits*icovs); i++) {
    temp_inc_beta[i]  = new_inc_beta[i];
    inc_beta[i]       = new_inc_beta[i];
  }
  stephalf_flag    = 0;
  *converg1_2_flag = NO;
  reachboundary_flag = NO;

  do {
    if (stephalf_flag > 0) {
      stephalf_val = 1.0 * stephalf_flag;
      denom_half   = pow(2.0,stephalf_val);

      if (mu_flag == ESTIMATE) {
	for (i = 0; i < itraits; i++) {
	  x_mu[i] = prev_mu[i];
	  for (j = 0; j < icovs; j++) {
	    x_beta[i][j] = prev_beta[i][j];
	  }
	}
      }
      for (i = 0; i < env_vcnum; i++) {
	if (t_flag_array[i] == ESTIMATE) {
	  x_T[i]  = prev_t[i];
	}
      }
      for (i = 0; i < iinitvcnum; i++) {
	if (s_flag_array[i] == ESTIMATE) {
	  x_S[i]  = prev_s[i];
	}
	if (m1_flag_array[i] == ESTIMATE) {
	  x_M1[i] = prev_m1[i];
	}
	if (m2_flag_array[i] == ESTIMATE) {
	  x_M2[i] = prev_m2[i];
	}
	if (c_flag_array[i] == ESTIMATE) {
	  x_c[i]  = prev_c[i];
	}
	if (p_flag_array[i] == ESTIMATE) {
	  x_p[i]  = prev_p[i];
	}
	if (q_flag_array[i] == ESTIMATE) {
	  x_q[i]  = prev_q[i];
	}
      }

      /* step halfing increase values for all the parameters. */
      for (i = 0; i < itraits; i++) {
	inc_mu[i] = temp_inc_mu[i]/denom_half;
      }
      for (i = 0; i < env_vcnum; i++) {
	inc_t[i]  = temp_inc_t[i]/denom_half;
      }
      for (i = 0; i < iinitvcnum; i++) {
	inc_s[i]  = temp_inc_s[i]/denom_half;
	inc_m1[i] = temp_inc_m1[i]/denom_half;
	inc_m2[i] = temp_inc_m2[i]/denom_half;
	inc_c[i]  = temp_inc_c[i]/denom_half;
	inc_p[i]  = temp_inc_p[i]/denom_half;
	inc_q[i]  = temp_inc_q[i]/denom_half;
      }
      for (i = 0; i < (itraits*icovs); i++) {
	inc_beta[i]  = temp_inc_beta[i]/denom_half;
      }
    }

    /* Change G_VAL 'x_S', 'x_M1', 'x_M2', 'x_T' and
       'x_c', 'x_p', 'x_p' values here.   */
    SCondOfInc();

    /* Change G_VAL 'x_mu', 'x_beta' values here.      */
    MCondOfInc();

    CheckFixTraitcov();

    CalculateDiffOfPara(prev_mu, prev_beta, prev_s, prev_m1, prev_m2, prev_t,
			prev_c, prev_p, prev_q, &converg1_flag,
			stephalf_flag);
    if (fixtoboundary_flag == YES) {
      reachboundary_flag = Fixtoboundary();
    }

    /* Get new G_VAL 'imubeta_dim' and 'ismtcpq_dim' here.             */
    GetDim();

    /* Get initial G_VAL 'imubetacovnum' and 'itraitcovnum' value here and
       they may be changed if mu_flag or smtcpq_flag changed.          */
    imubetacovnum= imubeta_dim*(imubeta_dim+1)/2;
    itraitcovnum = ismtcpq_dim*(ismtcpq_dim+1)/2;

    if (ascert_flag == YES) {
      ln_func = AS_CalculateValues();
    }else if (ascert_flag == NO) {
      ln_func = NOAS_CalculateValues();
    }

    if (reachboundary_flag == YES) {
      /* changed 10-10-96.
	 ln_func = prev_ln_func;
      */
      break;
    }
    stephalf_flag++;

    /* break step-halfing look if stephalf_flag >= BREAKSTEPHALF; 10-10-96. */
    if (stephalf_flag >= BREAKSTEPHALF) {
      breakstephalf_flag = YES;
      break;
    }
    if ((fabs((ln_func - prev_ln_func)/prev_ln_func)) <= LogLikbreakval){
      if (ln_func >= prev_ln_func) {
	converg2_flag = YES;
	break;
      }
    }
  }while(ln_func < prev_ln_func); // end 'do {'

  if ((converg1_flag == YES) && (converg2_flag == YES))
    *converg1_2_flag = YES;

  free(temp_inc_mu);
  free(temp_inc_beta);
  free(temp_inc_s);
  free(temp_inc_m1);
  free(temp_inc_m2);
  free(temp_inc_t);
  free(temp_inc_c);
  free(temp_inc_p);
  free(temp_inc_q);

  return ln_func;
}


// Check and fix covariate's value according to the following
//   condition: [-cov(i)*cov(j)] <=  cov(ij) <= [cov(i)*cov(j)]
void Calculate::CheckFixTraitcov(){
  int i,j,k;
  double Va1, Va2, Vc1, Vc2;
  double **temp_s, **temp_m1, **temp_m2, **temp_t;

  temp_s  = (double **) malloc(total_trait_values * sizeof(double *));
  temp_m1 = (double **) malloc(total_trait_values * sizeof(double *));
  temp_m2 = (double **) malloc(total_trait_values * sizeof(double *));
  temp_t =  (double **) malloc(total_trait_values * sizeof(double *));

  for (i = 0; i < total_trait_values; i++) {
    temp_s[i]  = (double *) malloc(total_trait_values * sizeof(double));
    temp_m1[i] = (double *) malloc(total_trait_values * sizeof(double));
    temp_m2[i] = (double *) malloc(total_trait_values * sizeof(double));
    temp_t[i]  = (double *) malloc(total_trait_values * sizeof(double));
  }

  /* put x_S and x_M1, x_M2 to two dimension arrays.  */
  k = 0;
  for ( i = 0; i < total_trait_values; i++) {
    for ( j = i; j < total_trait_values; j++) {
      temp_s[i][j]  = x_S[k];
      temp_m1[i][j] = x_M1[k];
      temp_m2[i][j] = x_M2[k];
      //      if (datatype == LONGITUDINAL)
	temp_t[i][j]  = x_T[k];
      k++;
    }
  }

  for ( i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if (j > i) {
	/* Check x_S (sigma) covariates .  */
	Va1 = fabs(temp_s[i][j]);
	Va2 = sqrt(temp_s[i][i]*temp_s[j][j]);
	if (s_flag == ESTIMATE) {
	  if ( Va1 > Va2 ) {
	    /*
	      printf("changing the x_S values of covariates\n");
	    */
	    if (temp_s[i][j] < 0)
	      temp_s[i][j] = -Va2;
	    if (temp_s[i][j] > 0)
	      temp_s[i][j] = Va2;
	  }
	}else if (s_flag == CONSTRAINT) {
	  if (temp_s[i][j] < 0)
	    temp_s[i][j] = -Va2;
	  if (temp_s[i][j] > 0)
	    temp_s[i][j] = Va2;
	}

	/* Check x_M1 (major gene) covariates .  */
	Vc1 = fabs(temp_m1[i][j]);
	Vc2 = sqrt(temp_m1[i][i]*temp_m1[j][j]);
	if (m1_flag == ESTIMATE) {
	  if ( Vc1 > Vc2 ) {
	    /*
	      printf("changing the x_M1 values of covariates\n");
	    */
	    if (temp_m1[i][j] < 0)
	      temp_m1[i][j] = -Vc2;
	    if (temp_m1[i][j] > 0)
	      temp_m1[i][j] = Vc2;
	  }
	}else if (m1_flag == CONSTRAINT) {
	  if (temp_m1[i][j] < 0)
	    temp_m1[i][j] = -Vc2;
	  if (temp_m1[i][j] > 0)
	    temp_m1[i][j] = Vc2;
	}

	/* Check x_M2 (second major gene) covariates .  */
	Vc1 = fabs(temp_m2[i][j]);
	Vc2 = sqrt(temp_m2[i][i]*temp_m2[j][j]);
	if (m2_flag == ESTIMATE) {
	  if ( Vc1 > Vc2 ) {
	    /*
	      printf("changing the x_M2 values of covariates\n");
	    */
	    if (temp_m2[i][j] < 0)
	      temp_m2[i][j] = -Vc2;
	    if (temp_m2[i][j] > 0)
	      temp_m2[i][j] = Vc2;
	  }
	}else if (m2_flag == CONSTRAINT) {
	  if (temp_m2[i][j] < 0)
	    temp_m2[i][j] = -Vc2;
	  if (temp_m2[i][j] > 0)
	    temp_m2[i][j] = Vc2;
	}
	//	if (datatype == LONGITUDINAL) {
	  /* Check x_T (environment) covariates .  */
	  Vc1 = fabs(temp_t[i][j]);
	  Vc2 = sqrt(temp_t[i][i]*temp_t[j][j]);
	  if (t_flag == ESTIMATE) {
	    if ( Vc1 > Vc2 ) {
	      /*
		printf("changing the x_T values of covariates\n");
	      */
	      if (temp_t[i][j] < 0)
		temp_t[i][j] = -Vc2;
	      if (temp_t[i][j] > 0)
		temp_t[i][j] = Vc2;
	    }
	  }else if (t_flag == CONSTRAINT) {
	    if (temp_t[i][j] < 0)
	      temp_t[i][j] = -Vc2;
	    if (temp_t[i][j] > 0)
	      temp_t[i][j] = Vc2;
	  }
	  //	}
      }
    }
  }

  /* restore the values of x_S and x_M1, x_M2. */
  k = 0;
  for ( i = 0; i < total_trait_values; i++) {
    for ( j = i; j < total_trait_values; j++) {
      x_S[k]  = temp_s[i][j];
      x_M1[k] = temp_m1[i][j];
      x_M2[k] = temp_m2[i][j];
      //      if (datatype == LONGITUDINAL)
	x_T[k] = temp_t[i][j];
      k++;
    }
  }

  for (i = 0; i < total_trait_values; i++) {
    free(temp_s[i]);
    free(temp_m1[i]);
    free(temp_m2[i]);
    free(temp_t[i]);
  }
  free(temp_s);
  free(temp_m1);
  free(temp_m2);
  free(temp_t);
}




/* This is the increasement for s or m or c or p or q.
 */
void Calculate::SMCPQ_Inc(int x_flag, int x_flag_array[], double temp_val[],
			  double x_str[], double inc_x[]) {
  int i,j,k;

  k = total_trait_values;
  for ( i = 0,j = 0; i < iinitvcnum; i++) {
    if ((i == 0) || (i == j+k)) {
      /* Calculate the new values of x_str.       */
      if (x_flag_array[i] == ESTIMATE) {
	if ( x_str[i] == 0.0) {
	  if ( inc_x[i] < 0)  {
	    ;
	  } else {
	    x_str[i] += inc_x[i];
	  }
	}else if ( (x_str[i]+inc_x[i]) < 0) {
	  while ( (temp_val[i]+inc_x[i]) < 0)  {
	    inc_x[i] = inc_x[i]/2.0;
	    x_str[i]   = temp_val[i] + inc_x[i] ;
	  }
	}else {
	  x_str[i] += inc_x[i];
	}
      }

      j = i;
      if ( i != 0) {
	k--;
      }
    }else {
      if (x_flag_array[i] == ESTIMATE) {
	x_str[i] += inc_x[i];
      }
    }
  }
}

/* This is the condition of the increasement for s,m,t,c,p,q.
   It calls function:
	SMCPQ_Inc.
 */
void Calculate::SCondOfInc() {
  int i,j,k;
  double *temp_s, *temp_m1,*temp_m2, *temp_t;
  double *temp_c, *temp_p, *temp_q;

  temp_s  = (double *) malloc(iinitvcnum * sizeof(double));
  temp_m1 = (double *) malloc(iinitvcnum * sizeof(double));
  temp_m2 = (double *) malloc(iinitvcnum * sizeof(double));
  temp_t  = (double *) malloc(env_vcnum  * sizeof(double));
  temp_c  = (double *) malloc(iinitvcnum * sizeof(double));
  temp_p  = (double *) malloc(iinitvcnum * sizeof(double));
  temp_q  = (double *) malloc(iinitvcnum * sizeof(double));

  /* set values to the temp initiate values.   */
  for ( i = 0; i < env_vcnum; i++)
    temp_t[i] = x_T[i];

  for ( i = 0; i < iinitvcnum; i++) {
    temp_s[i]  = x_S[i];
    temp_m1[i] = x_M1[i];
    temp_m2[i] = x_M2[i];
    temp_c[i]  = x_c[i];
    temp_p[i]  = x_p[i];
    temp_q[i]  = x_q[i];
  }

  /// I'm currently working on finding the difference with s, but
  /// afterwards, I may want to also do thea sam with M1, and the others
  SMCPQ_Inc(s_flag,  s_flag_array,  temp_s,  x_S,  inc_s);
  SMCPQ_Inc(m1_flag, m1_flag_array, temp_m1, x_M1, inc_m1);
  SMCPQ_Inc(m2_flag, m2_flag_array, temp_m2, x_M2, inc_m2);
  SMCPQ_Inc(c_flag,  c_flag_array,  temp_c,  x_c,  inc_c);
  SMCPQ_Inc(p_flag,  p_flag_array,  temp_p,  x_p,  inc_p);
  SMCPQ_Inc(q_flag,  q_flag_array,  temp_q,  x_q,  inc_q);

/*
  if (datatype == MULTIVARIATE) {
    for ( i = 0; i < total_trait_values; i++) {
      // about new values of x_T.
      if (t_flag_array[i] == ESTIMATE) {
	if ( x_T[i] == 0.0) {
	  if (inc_t[i] < 0) {
	    ;
	  }else {
	    x_T[i] += inc_t[i];
	  }
	}else if( (x_T[i]+inc_t[i]) < 0) {
	  while ( (temp_t[i]+inc_t[i]) < 0)  {
	    //cout << "inc_T["<<i<<"] = " << inc_t[i] << endl;
	    //cout << "temp_t["<<i<<"] = " << temp_t[i] << endl;
	    inc_t[i] = inc_t[i]/2.0;
	    x_T[i]   = temp_t[i] + inc_t[i] ;
	  }
	}else {
	  //cout << "inc_T["<<i<<"] = " << inc_t[i] << endl;
	  x_T[i] += inc_t[i];
	}
      }
    }
  }else if (datatype == LONGITUDINAL) {
*/
    k = total_trait_values;
    for ( i = 0,j = 0; i < iinitvcnum; i++) {
      if ((i == 0) || (i == j+k)) {
	// Calculate the new values of x_T.
	if (t_flag_array[i] == ESTIMATE) {
	  if ( x_T[i] == 0.0) {
	    if ( inc_t[i] < 0) {
	      ;
	    }else {
	      x_T[i] += inc_t[i];
	    }
	  }else if ( (x_T[i]+inc_t[i]) < 0) {
	    while ( (x_T[i]+inc_t[i]) < 0)  {
	      inc_t[i] = inc_t[i]/2.0;
	      x_T[i] = temp_t[i] + inc_t[i] ;
	    }
	  }else {
	    x_T[i] += inc_t[i];
	  }
	}

	j = i;
	if ( i != 0) {
	  k--;
	}
      }else {
	if (t_flag_array[i] == ESTIMATE) {
	  x_T[i] += inc_t[i];
	}
      }
    }
    //  }
  free(temp_s);
  free(temp_m1);
  free(temp_m2);
  free(temp_t);
  free(temp_c);
  free(temp_p);
  free(temp_q);
}

/* This is the condition of the increasement for mu and beta.
 */
void Calculate::MCondOfInc() {
  int i,j,kv;

  kv = 0;
  /* use the new values of x_mu.          */
  if (mu_flag == ESTIMATE) {
    for (i = 0; i < itraits; i++) {
      x_mu[i] += inc_mu[i];

      for (j = 0; j < icovs; j++) {
	x_beta[i][j] += inc_beta[kv];
	kv++;
      }
    }
  }
}

/* Save estimated values for MU, SIGMA, MG, TAU
   and C,P,Q.
*/
//void Calculate::SaveMeanAndFlag(int run_flag) {

void Calculate::SaveMeanAndFlag() {
  int i,j;

  for (i = 0; i < itraits; i++) {
    /* MuBeta.         */
    for (j = 0; j < icovs+1; j++) {
      if ( j == 0) {
	if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1)))
	  x_mu_A[i] = x_mu[i];
	else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
					 || (runopt_flag == 4)))
	  x_mu_0[i] = x_mu[i];
      }
      else {
	if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1)))
	  x_beta_A[i][j-1]  = x_beta[i][j-1];
	else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
					 || (runopt_flag == 4)))
	  x_beta_0[i][j-1]  = x_beta[i][j-1];
      }
    }
  }

  for (i = 0; i < env_vcnum; i++) {
    /* TAU.            */
    if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1))) {
      x_T_A[i]              = x_T[i];
      t_A_flag_array[i]     = t_flag_array[i];
      t_A_fix_flag_array[i] = t_fix_flag_array[i];
    }
    else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
				     || (runopt_flag == 4))) {
      x_T_0[i]              = x_T[i];
      t_0_flag_array[i]     = t_flag_array[i];
      t_0_fix_flag_array[i] = t_fix_flag_array[i];
    }
  }

  for (i = 0; i < iinitvcnum; i++) {
    /* SIGMA, MG and C,P,Q.    */
    if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1))) {
      x_S_A[i] = x_S[i];
      x_M1_A[i]= x_M1[i];
      x_M2_A[i]= x_M2[i];
      x_c_A[i] = x_c[i];
      x_p_A[i] = x_p[i];
      x_q_A[i] = x_q[i];
      s_A_flag_array[i]  = s_flag_array[i];
      m1_A_flag_array[i] = m1_flag_array[i];
      m2_A_flag_array[i] = m2_flag_array[i];
      c_A_flag_array[i]  = c_flag_array[i];
      p_A_flag_array[i]  = p_flag_array[i];
      q_A_flag_array[i]  = q_flag_array[i];
      s_A_fix_flag_array[i]  = s_fix_flag_array[i];
      m1_A_fix_flag_array[i] = m1_fix_flag_array[i];
      m2_A_fix_flag_array[i] = m2_fix_flag_array[i];
      c_A_fix_flag_array[i]  = c_fix_flag_array[i];
      p_A_fix_flag_array[i]  = p_fix_flag_array[i];
      q_A_fix_flag_array[i]  = q_fix_flag_array[i];
    } else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
				     || (runopt_flag == 4))) {
      x_S_0[i]  = x_S[i];
      x_M1_0[i] = x_M1[i];
      x_M2_0[i] = x_M2[i];
      x_c_0[i]  = x_c[i];
      x_p_0[i]  = x_p[i];
      x_q_0[i]  = x_q[i];
      s_0_flag_array[i]  = s_flag_array[i];
      m1_0_flag_array[i] = m1_flag_array[i];
      m2_0_flag_array[i] = m2_flag_array[i];
      c_0_flag_array[i]  = c_flag_array[i];
      p_0_flag_array[i]  = p_flag_array[i];
      q_0_flag_array[i]  = q_flag_array[i];
      s_0_fix_flag_array[i]  = s_fix_flag_array[i];
      m1_0_fix_flag_array[i] = m1_fix_flag_array[i];
      m2_0_fix_flag_array[i] = m2_fix_flag_array[i];
      c_0_fix_flag_array[i]  = c_fix_flag_array[i];
      p_0_fix_flag_array[i]  = p_fix_flag_array[i];
      q_0_fix_flag_array[i]  = q_fix_flag_array[i];
    }
  }
}

/* Get standard errors for MU, SIGMA, MG, TAU
   and C, P, Q.
*/
//void Calculate::Get_SE(int run_flag, int ascert_flag) {
void Calculate::Get_SE() {
  int i,kv;
  
  //  cout << "GetDim - 3 key 2567" << endl;
  GetDim();      /* Get new 'imubeta_dim' and 'ismtcpq_dim'.   */
  
  if (ascert_flag == YES) {
    kv = 0;
    if (mu_flag == ESTIMATE) {
      for (i = 0; i < imubeta_dim; i++) {
	if ((run_flag == SECOND)
	    || ((run_flag == FIRST) && (runopt_flag == 1))) {
	  // Added code to alert the user if inv_Exsder_mat is negative
	  // Eric Lunde, 2005-08-26
	  if(inv_Exsder_mat[kv][kv] < 0) {
	    PROBLEM "inv_Exsder_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2657\n",
	      kv, kv, inv_Exsder_mat[kv][kv] RECOVER(NULL_ENTRY);
	  }
	  mubeta_A_sd[i] = sqrt(inv_Exsder_mat[kv][kv]);
	} else if ((run_flag == FIRST)
		   && ((runopt_flag == 2) || (runopt_flag == 3)
		       || (runopt_flag == 4))) {
	  // Added code to alert the user if inv_Exsder_mat is negative
	  // Eric Lunde, 2005-08-26
	  if(inv_Exsder_mat[kv][kv] < 0) {
	    PROBLEM "inv_Exsder_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2688\n",
	      kv, kv, inv_Exsder_mat[kv][kv] RECOVER(NULL_ENTRY);
	  }
	  mubeta_0_sd[i] = sqrt(inv_Exsder_mat[kv][kv]);
	  kv++;
	}
      }
    }
    for (i = 0; i < ismtcpq_dim; i++) {
      if ((run_flag == SECOND) ||
	  ((run_flag == FIRST) && (runopt_flag == 1))) {	
	// Added code to alert the user if inv_Exsder_mat is negative
	// Eric Lunde, 2005-08-26
	if(inv_Exsder_mat[kv][kv] < 0) {
	  PROBLEM "inv_Exsder_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2640\n",
	    kv, kv, inv_Exsder_mat[kv][kv] RECOVER(NULL_ENTRY);
	}
	smtcpq_A_sd[i] = sqrt(inv_Exsder_mat[kv][kv]);
      } else if ((run_flag == FIRST)
		 && ((runopt_flag == 2) || (runopt_flag == 3)
		     || (runopt_flag == 4)))
	// Added code to alert the user if inv_Exsder_mat is negative
	// Eric Lunde, 2005-08-26
	if(inv_Exsder_mat[kv][kv] < 0) {
	  PROBLEM "inv_Exsder_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2650\n",
	    kv, kv, inv_Exsder_mat[kv][kv] RECOVER(NULL_ENTRY);
	}
      smtcpq_0_sd[i] = sqrt(inv_Exsder_mat[kv][kv]);
      kv++;
    }
  } else {
    if (mu_flag == ESTIMATE) {
      for (i = 0; i < imubeta_dim; i++) {
	if ((run_flag == SECOND)
	    || ((run_flag == FIRST) && (runopt_flag == 1))) {
	  // Added code to alert the user if inv_Exsder_mubeta_mat is
	  // negative
	  // Eric Lunde, 2005-08-26
	  if(inv_Exsder_mubeta_mat[i][i] < 0) {
	    PROBLEM "inv_Exsder_mubeta_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2708\n",
	      i, i, inv_Exsder_mubeta_mat[i][i] RECOVER(NULL_ENTRY);
	  }
	  mubeta_A_sd[i] = sqrt(inv_Exsder_mubeta_mat[i][i]);
	}else if ((run_flag == FIRST)
		  && ((runopt_flag == 2) || (runopt_flag == 3)
		      || (runopt_flag == 4))) {
	  // Added code to alert the user if inv_Exsder_mubeta_mat is
	  // negative
	  // Eric Lunde, 2005-08-26
	  if(inv_Exsder_mubeta_mat[i][i] < 0) {
	    PROBLEM "inv_Exsder_mubeta_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2718\n",
	      i, i, inv_Exsder_mubeta_mat[i][i] RECOVER(NULL_ENTRY);
	  }
	  mubeta_0_sd[i] = sqrt(inv_Exsder_mubeta_mat[i][i]);
	}
      }
    }
    
    for (i = 0; i < ismtcpq_dim; i++) {
      if ((run_flag == SECOND)
	  || ((run_flag == FIRST) && (runopt_flag == 1))) {
	// Added code to alert the user if inv_Exsder_smtcpq_mat is
	// negative
	// Eric Lunde, 2005-08-26
	if(inv_Exsder_smtcpq_mat[i][i] < 0) {
	  PROBLEM "inv_Exsder_smtcpq_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2732\n",
	    i, i, inv_Exsder_smtcpq_mat[i][i] RECOVER(NULL_ENTRY);
	}
	smtcpq_A_sd[i]   = sqrt(inv_Exsder_smtcpq_mat[i][i]);
      }else if ((run_flag == FIRST)
		&& ((runopt_flag == 2) || (runopt_flag == 3)
		    || (runopt_flag == 4))) {
	// Added code to alert the user if inv_Exsder_smtcpq_mat is
	// negative
	// Eric Lunde, 2005-08-26
	if(inv_Exsder_smtcpq_mat[i][i] < 0) {
	  PROBLEM "inv_Exsder_smtcpq_mat[%d][%d] (%f) < 0.\nCalculate.cpp key 2732\n",
	    i, i, inv_Exsder_smtcpq_mat[i][i] RECOVER(NULL_ENTRY);
	}
	smtcpq_0_sd[i]   = sqrt(inv_Exsder_smtcpq_mat[i][i]);
      }
    }
  }
}

/*
This module calls the function:
	SRun.
*/



/*****************************
Class name: Calculate
Method name: Run
Description: I'm not sure about this method at all.  So I have no description.
             But here is the "documenation" when I started.

             "Once we got everything from the parameter file, it is the time to
             do the analysis under the hypotheses of with and without major
             gene components."

Input: NONE.
Output: NONE.
Side Effects:
Author: Eric Lunde, 6-12-03
******************************/
void Calculate::Run() {
  int    i;
  int    mubeta_dim, smtcpq_dim, i_bigdim;

  /* Get initial G_VAL 'run_flag' value here and it may be
     changed within this function.                              */
  run_flag     = FIRST;

  /* Get initial G_VAL 'converg_flag' values here and it may be
     changed within function "SRun".                            */
  converg_flag = NO;

  /* Set initial G_VAL 'printfirstfam_flag' value here and it
     may be changed within function "GetData".                  */
  printfirstfam_flag = YES;

  /* Set the initial G_VAL 'breakstephalf_run1' and 'breakstephalf_run2' here
     and they may be changed when run SRun.                     */
  breakstephalf_run1 = NO;
  breakstephalf_run2 = NO;

  do {
    if(run_flag == FIRST) {
      if ((runopt_flag == 3) || (runopt_flag == 4)) {
	GetStartMubeta_SMT();
	GetStartCPQ();

	for (i = 0; i < iinitvcnum; i++) {
	  x_M1[i] = 0.0;
	  x_M2[i] = 0.0;
	}
	m1_flag  = FIXED;
	m2_flag  = FIXED;

	GetST_flag_array();
	GetM_flag_array();
	GetCPQ_flag_array();

	need_ibd_flag = NO;
      }
    }else if (run_flag == SECOND) {
      if (runopt_flag == 3) {
	GetStartMubeta_SMT();
	GetStartCPQ();

	GetST_flag_array();
	GetM_flag_array();
	GetCPQ_flag_array();

	need_ibd_flag = YES;
      }
      if (runopt_flag == 4) {
	/* start with the values which come from the null hypotheses except
	   using the initial values for major gene in the parameter file. */
	for (i = 0; i < iinitvcnum; i++) {
	  x_M1[i] = init_x_M1[i];
	  x_M2[i] = init_x_M2[i];
	}

	/* check the parameter fix flag and see whether the parameter has been
	   fixed to BOUNDARY internally or not if user chose
	   'fixtoboundary_flag =y' option.
	   added 2-4-98; */
	if (fixtoboundary_flag == YES)
	  CheckBoundary();

	/* get 'mu_flag', 's_flag','m1_flag','m2_flag','t_flag' values from initial
	   values here.  */
	mu_flag = init_mu_flag;
	s_flag  = init_s_flag ;
	m1_flag = init_m1_flag;
	m2_flag = init_m2_flag;
	t_flag  = init_t_flag ;
	c_flag  = init_c_flag ;
	p_flag  = init_p_flag ;
	q_flag  = init_q_flag ;

	GetST_flag_array();
	GetM_flag_array();
	GetCPQ_flag_array();

	need_ibd_flag = YES;
      }
    }

    GetDim();
    mubeta_dim = imubeta_dim;
    smtcpq_dim = ismtcpq_dim;
    i_bigdim   = mubeta_dim + smtcpq_dim;

    if (ascert_flag == YES) {
      inv_Exsder_mat = Lib::dmatrix(0,i_bigdim-1, 0, i_bigdim-1);
    }else {
      inv_Exsder_mubeta_mat = Lib::dmatrix(0,mubeta_dim-1,0,mubeta_dim-1);
      inv_Exsder_smtcpq_mat = Lib::dmatrix(0,smtcpq_dim-1,0,smtcpq_dim-1);
    }

    SRun();

    /* Save G_VAL 'breakstephalf_run1' and 'breakstephalf_run2' for
       PrintSummary function.                                 */
    if (run_flag == FIRST)
      breakstephalf_run1 = breakstephalf_flag;
    else if (run_flag == SECOND)
      breakstephalf_run2 = breakstephalf_flag;

    //	SaveMeanAndFlag(run_flag);
    SaveMeanAndFlag();

    writeInvExpSecDer();

    /* catch standard error arrays here.                      */
    Get_SE();

    if(ascert_flag == YES) {
      Lib::free_dmatrix(inv_Exsder_mat, 0,i_bigdim-1, 0, i_bigdim-1);
    }else {
      Lib::free_dmatrix(inv_Exsder_mubeta_mat, 0,mubeta_dim-1,0,mubeta_dim-1);
      Lib::free_dmatrix(inv_Exsder_smtcpq_mat, 0,smtcpq_dim-1,0,smtcpq_dim-1);
    }

    if (run_flag == FIRST) {
      /* Change G_VAL 'run_flag' value here.                */
      run_flag = SECOND;
    }else if (run_flag == SECOND) {
      /* Change G_VAL 'run_flag' value here.                */
      run_flag = END;
    }

  } while ((run_flag == SECOND) && ((runopt_flag == 3) ||(runopt_flag == 4)));
}

/*****************************
Class name: Calculate
Method name: SRun
Description: Once we got everything from the parameter file, it is the time to
             do the analysis and to get the new MU,SIGMA,MG and TAU.
             (Description written before Eric began)
Input: NONE.
Output: NONE.
Side Effects:
Author: Eric Lunde, 6-12-03
******************************/
void Calculate::SRun() {
  int    i,j,k;
  // int    nfam,  N[MAXNUMFAM];
  int    break_flag = NO;   /* break flag for the convergent point.  */
  int    converg1_2_flag;
  int    satisfiedcount;
  double ln_func;
  double prev_ln_func;
  double *prev_s,  *prev_m1, *prev_m2, *prev_t;
  double *prev_mu, **prev_beta;
  double *prev_c,  *prev_p, *prev_q;

  /* Initialize the size of pointers.    5-12-98; E.Y.      */
  prev_s    = (double *) malloc(iinitvcnum * sizeof(double));
  prev_m1   = (double *) malloc(iinitvcnum * sizeof(double));
  prev_m2   = (double *) malloc(iinitvcnum * sizeof(double));
  prev_t    = (double *) malloc(env_vcnum  * sizeof(double));
  prev_mu   = (double *) malloc(itraits    * sizeof(double));
  prev_beta = (double **)malloc(itraits    * sizeof(double *));
  for (i = 0; i < itraits; i++)
    prev_beta[i] = (double *) malloc(icovs * sizeof(double));
  prev_c    = (double *) malloc(iinitvcnum * sizeof(double));
  prev_p    = (double *) malloc(iinitvcnum * sizeof(double));
  prev_q    = (double *) malloc(iinitvcnum * sizeof(double));

  /* Get initial G_VAL 'final_flag' value here and it may be
     changed within this function .                                  */
  final_flag = NO;

  /* Get initial G_VAL 'inc_mu', 'inc_t', 'inc_s', 'inc_m1', 'inc_m2','inc_beta',
     'inc_c', 'inc_p', 'inc_q' values here and they may be changed
     within the iterations.                                          */
  for (i = 0; i < itraits; i++) {
    inc_mu[i] = 0.0;
  }
  for (i = 0; i < env_vcnum; i++) {
    inc_t[i]  = 0.0;
  }
  for (i = 0; i < iinitvcnum; i++) {
    inc_s[i]  = 0.0;
    inc_m1[i] = 0.0;
    inc_m2[i] = 0.0;
    inc_c[i]  = 0.0;
    inc_p[i]  = 0.0;
    inc_q[i]  = 0.0;
  }
  for (i = 0; i < (itraits*icovs); i++) {
    inc_beta[i]  = 0.0;
  }

  satisfiedcount     = 0;
  breakstephalf_flag = NO;
  ln_func = 0.0; // add by Chen to initialize ln_func value

  /* Do N1 iteration by using the same family data file.      */
  /* First time use G_VAL 'N1' in the program.                */

  if(printProgress) {
    cout << "Iteration Number:\n";
  }
  for (j = 0; j < N1; j++) {
    // Print iteration count for user to know that the program is working and 
    // intervene if the amount of iterations is suspisious.
    // Eric Lunde 01-08-04
    printIterationNumber(j + 1);
    
    /* Save all the previous parameters. */
    for (i = 0; i < itraits; i++) {
      prev_mu[i] = x_mu[i];
      for (k = 0; k < icovs; k++) {
	prev_beta[i][k] = x_beta[i][k];
      }
    }
    for (i = 0; i < env_vcnum; i++) {
      prev_t[i]  = x_T[i];
    }
    for (i = 0; i < iinitvcnum; i++) {
      prev_s[i]  = x_S[i];
      prev_m1[i] = x_M1[i];
      prev_m2[i] = x_M2[i];
      prev_c[i]  = x_c[i];
      prev_p[i]  = x_p[i];
      prev_q[i]  = x_q[i];
    }
    prev_ln_func = ln_func;

    /* Calculate Log likelihood here.                               */
    if (j == 0) {
      if (fixtoboundary_flag == YES)
	Fixtoboundary();

      /* Get new G_VAL 'imubeta_dim' and 'ismtcpq_dim' here.         */
      GetDim();

      /* Get initial G_VAL 'imubetacovnum' and 'itraitcovnum' value here and
	 they may be changed if mu_flag or smtcpq_flag changed.      */
      imubetacovnum= imubeta_dim*(imubeta_dim+1)/2;
      itraitcovnum = ismtcpq_dim*(ismtcpq_dim+1)/2;

      if (ascert_flag == YES) {
	ln_func = AS_CalculateValues();
      }else if (ascert_flag == NO) {
	ln_func = NOAS_CalculateValues();
      }
    }else {
      ln_func = CalculateLik(prev_mu, prev_beta, prev_s, prev_m1,
			     prev_m2, prev_t,    prev_c, prev_p,
			     prev_q,  prev_ln_func,  &converg1_2_flag);
      if (breakstephalf_flag == YES) {
	break;
      }
    }

    if ( (j > 0) && (converg1_2_flag == YES)) {
      satisfiedcount++;

      if (satisfiedcount == SATISFIEDTIMES) {
	break_flag = YES;
      }
    }

    /* break the loop.         */
    if (break_flag == YES) {
      /* Change G_VAL 'final_flag' value here.          */
      final_flag = YES;

      if (run_flag == FIRST) {
	/* Change G_VAL 'converg_flag' value here.     */
	converg_flag = YES;

	if ((runopt_flag == 1) || (runopt_flag == 2)) {
	  for (i = 0; i < itraits; i++) {
	    fp_out << setw(10) << x_mu[i];
	    fp_mu << setw(10) << x_mu[i];
	    for (k = 0; k < icovs; k++) {
	      fp_out << setw(10) << x_beta[i][k];
	      fp_mu << setw(10) << x_beta[i][k];
	    }
	  }
	  fp_out << endl;
	  fp_mu << endl;

	  /* print the parameters to the output files.  */
	  PrintResults(NONDEBUG);
	}
      }else if (run_flag == SECOND) {
	if ((runopt_flag == 3) || (runopt_flag == 4)) {
	  for (i = 0; i < itraits; i++) {
	    fp_out << setw(10) <<  x_mu[i];
	    fp_mu << setw(10) << x_mu[i];
	    for (k = 0; k < icovs; k++) {
	      fp_out << setw(10) << x_beta[i][k];
	      fp_mu << setw(10) << x_beta[i][k];
	    }
	  }
	  fp_out << endl;
	  fp_mu << endl;

	  /* print the parameters to the output files.  */
	  PrintResults(NONDEBUG);
        }
      }
      dataIndex = 0;
      break;
    } // Closes 'if (break_flag == YES) {'
    dataIndex = 0;
  } // Closes 'for (j = 0; j < N1; j++) {'  
  cout << endl;

  // Write the number of iterations to a file for Splus to read back in
  ofstream iter(iteration_file, ios::app);
  iter << markerName << endl;
  iter << (j + 1) << endl;
  iter.close();

  // Write to the appropriate file each family's loglikelihood
  if(run_flag == FIRST && runopt_flag != 1) {
    ofstream fam_lik("fam.lik");
    if(fam_lik.fail()) {
      PROBLEM "The file fam.lik could not be opened for writing\nCalculate.cpp key 2786\n"
	RECOVER(NULL_ENTRY);
    }
    // Write family ids to fam.lik0
    for(int i=0; i<nfam; i++) {
      fam_lik << " " /* setw(15) */ << uniqueFamilyIds[i];
    }
    fam_lik << endl;
    fam_lik << "null";
    for(int i=0; i<nfam; i++) {
      fam_lik << " " /* setw(15) */ << familyLoglikelihoods[i];
    }
    fam_lik << endl;
    fam_lik.close();
  }else if(run_flag == SECOND
	   || (run_flag == FIRST && runopt_flag == 1) ) {
    ofstream fam_lik("fam.lik", ios::app);
    if(fam_lik.fail()) {
      PROBLEM "The file fam.lik could not be opened for appending\nCalculate.cpp key 2799\n"
	RECOVER(NULL_ENTRY);
    }
    fam_lik << ibdFileName;
    for(int i=0; i<nfam; i++) {
      fam_lik << " " /* setw(15) */ << familyLoglikelihoods[i];
    }
    fam_lik << endl;
    fam_lik.close();
  }

  if ((run_flag == SECOND) || ((run_flag == FIRST) && (runopt_flag == 1))) {
    fp_lik << ln_func << endl;
    loglikelihood_A = ln_func;
  }
  else if ((run_flag == FIRST) && ((runopt_flag == 2) || (runopt_flag == 3)
				   || (runopt_flag == 4))) {
    loglikelihood_0 = ln_func;
  }

  dataIndex = 0;

  free(prev_s);
  free(prev_m1); 
  free(prev_m2);
  free(prev_t);
  free(prev_mu);
  for (i = 0; i < itraits; i++)
    free(prev_beta[i]);
  free(prev_beta);
  free(prev_c);
  free(prev_p);
  free(prev_q);
}



/* If you fixed some paremters internally, you have to fix covariance internally as well.
   If you fixed one of the major gene then you have to use the initial values for all other
   parameters.
   10-24-96; Added  by Emily.
*/
void Calculate::FixCovAndAssignInitValues() {
  int     i,j,k;
  int     fix_mg_flag = NO;
  double  **temp_s, **temp_m1, **temp_m2,**temp_t;
  int  **temp_s_flag_array, **temp_m1_flag_array,**temp_m2_flag_array,**temp_t_flag_array;
  int  **temp_s_fix_flag_array, **temp_m1_fix_flag_array,
    **temp_m2_fix_flag_array,**temp_t_fix_flag_array;
  double  **temp_c, **temp_p, **temp_q;
  int  **temp_c_flag_array, **temp_p_flag_array, **temp_q_flag_array;
  int  **temp_c_fix_flag_array, **temp_p_fix_flag_array, **temp_q_fix_flag_array;

  /* Initialize the size of the pointers .    5-12-98;  E.Y.  */
  temp_s  = (double **) malloc(total_trait_values * sizeof(double *));
  temp_m1 = (double **) malloc(total_trait_values * sizeof(double *));
  temp_m2 = (double **) malloc(total_trait_values * sizeof(double *));
  temp_t  = (double **) malloc(total_trait_values * sizeof(double *));
  temp_c  = (double **) malloc(total_trait_values * sizeof(double *));
  temp_p  = (double **) malloc(total_trait_values * sizeof(double *));
  temp_q  = (double **) malloc(total_trait_values * sizeof(double *));

  temp_s_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_m1_flag_array  = (int **) malloc(total_trait_values * sizeof(int *));
  temp_m2_flag_array  = (int **) malloc(total_trait_values * sizeof(int *));
  temp_t_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_c_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_p_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_q_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));

  temp_s_fix_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_m1_fix_flag_array  = (int **) malloc(total_trait_values * sizeof(int *));
  temp_m2_fix_flag_array  = (int **) malloc(total_trait_values * sizeof(int *));
  temp_t_fix_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_c_fix_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_p_fix_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));
  temp_q_fix_flag_array   = (int **) malloc(total_trait_values * sizeof(int *));

  for (i = 0; i < total_trait_values; i++) {
    temp_s[i]  = (double *) malloc(total_trait_values * sizeof(double));
    temp_m1[i] = (double *) malloc(total_trait_values * sizeof(double));
    temp_m2[i] = (double *) malloc(total_trait_values * sizeof(double));
    temp_t[i]  = (double *) malloc(total_trait_values * sizeof(double));
    temp_c[i]  = (double *) malloc(total_trait_values * sizeof(double));
    temp_p[i]  = (double *) malloc(total_trait_values * sizeof(double));
    temp_q[i]  = (double *) malloc(total_trait_values * sizeof(double));

    temp_s_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_m1_flag_array[i] = (int *) malloc(total_trait_values * sizeof(int));
    temp_m2_flag_array[i] = (int *) malloc(total_trait_values * sizeof(int));
    temp_t_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_c_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_p_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_q_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));

    temp_s_fix_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_m1_fix_flag_array[i] = (int *) malloc(total_trait_values * sizeof(int));
    temp_m2_fix_flag_array[i] = (int *) malloc(total_trait_values * sizeof(int));
    temp_t_fix_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_c_fix_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_p_fix_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
    temp_q_fix_flag_array[i]  = (int *) malloc(total_trait_values * sizeof(int));
  }


  /* put x_S, x_M1, x_M2, x_c, x_p, x_q to two dimension arrays.  */
  k = 0;
  for ( i = 0; i < total_trait_values; i++) {
    for ( j = i; j < total_trait_values; j++) {
      temp_s[i][j]  = x_S[k];
      temp_m1[i][j] = x_M1[k];
      temp_m2[i][j] = x_M2[k];
      temp_c[i][j]  = x_c[k];
      temp_p[i][j]  = x_p[k];
      temp_q[i][j]  = x_q[k];
      temp_s_flag_array[i][j]  =  s_flag_array[k];
      temp_m1_flag_array[i][j] =  m1_flag_array[k];
      temp_m2_flag_array[i][j] =  m2_flag_array[k];
      temp_c_flag_array[i][j]  =  c_flag_array[k];
      temp_p_flag_array[i][j]  =  p_flag_array[k];
      temp_q_flag_array[i][j]  =  q_flag_array[k];
      temp_s_fix_flag_array[i][j]  =  s_fix_flag_array[k];
      temp_m1_fix_flag_array[i][j] =  m1_fix_flag_array[k];
      temp_m2_fix_flag_array[i][j] =  m2_fix_flag_array[k];
      temp_c_fix_flag_array[i][j]  =  c_fix_flag_array[k];
      temp_p_fix_flag_array[i][j]  =  p_fix_flag_array[k];
      temp_q_fix_flag_array[i][j]  =  q_fix_flag_array[k];

      //      if (datatype == LONGITUDINAL) {
	temp_t[i][j]  = x_T[k];
	temp_t_flag_array[i][j]  =  t_flag_array[k];
	temp_t_fix_flag_array[i][j]  =  t_fix_flag_array[k];
	//      }

      if ( (m1_fix_flag_array[k] == INTERNALLY) ||
	   (m2_fix_flag_array[k] == INTERNALLY) )
	fix_mg_flag = YES;
      k++;
    }
  }

  for ( i = 0; i < total_trait_values; i++) {
    for (j = i; j < total_trait_values; j++) {
      if (j > i) {
	/* Check x_S (sigma) covariates .  */
	if ((temp_s_fix_flag_array[i][i] == INTERNALLY) ||
	    (temp_s_fix_flag_array[j][j] == INTERNALLY)) {
	  if ((temp_s[i][i] == BOUNDARY) ||
	      (temp_s[j][j] == BOUNDARY)) {
	    temp_s[i][j] = BOUNDARY;
	    temp_s_flag_array[i][j] = FIXED;
	    temp_s_fix_flag_array[i][j] = INTERNALLY;
	  }
	}

	/* Check x_M1 (sigma) covariates .  */
	if ((temp_m1_fix_flag_array[i][i] == INTERNALLY) ||
	    (temp_m1_fix_flag_array[j][j] == INTERNALLY)) {
	  if ((temp_m1[i][i] == BOUNDARY) ||
	      (temp_m1[j][j] == BOUNDARY)) {
	    temp_m1[i][j] = BOUNDARY;
	    temp_m1_flag_array[i][j] = FIXED;
	    temp_m1_fix_flag_array[i][j] = INTERNALLY;
	  }
	}

	/* Check x_M2 (sigma) covariates .  */
	if ((temp_m2_fix_flag_array[i][i] == INTERNALLY) ||
	    (temp_m2_fix_flag_array[j][j] == INTERNALLY)) {
	  if ((temp_m2[i][i] == BOUNDARY) ||
	      (temp_m2[j][j] == BOUNDARY)) {
	    temp_m2[i][j] = BOUNDARY;
	    temp_m2_flag_array[i][j] = FIXED;
	    temp_m2_fix_flag_array[i][j] = INTERNALLY;
	  }
	}

	//	if (datatype == LONGITUDINAL) {
	  /* Check x_T (TAU) covariates .  */
	  if ((temp_t_fix_flag_array[i][i] == INTERNALLY) ||
	      (temp_t_fix_flag_array[j][j] == INTERNALLY)) {
	    if ((temp_t[i][i] == BOUNDARY) ||
		(temp_t[j][j] == BOUNDARY)) {
	      temp_t[i][j] = BOUNDARY;
	      temp_t_flag_array[i][j] = FIXED;
	      temp_t_fix_flag_array[i][j] = INTERNALLY;
	    }
	  }
	  //	}

	/* Check x_c (sibship) covariates .  */
	if ((temp_c_fix_flag_array[i][i] == INTERNALLY) ||
	    (temp_c_fix_flag_array[j][j] == INTERNALLY)) {
	  if ((temp_c[i][i] == BOUNDARY) ||
	      (temp_c[j][j] == BOUNDARY)) {
	    temp_c[i][j] = BOUNDARY;
	    temp_c_flag_array[i][j] = FIXED;
	    temp_c_fix_flag_array[i][j] = INTERNALLY;
	  }
	}

	/* Check x_p (parent-parent) covariates .  */
	if ((temp_p_fix_flag_array[i][i] == INTERNALLY) ||
	    (temp_p_fix_flag_array[j][j] == INTERNALLY)) {
	  if ((temp_p[i][i] == BOUNDARY) ||
	      (temp_p[j][j] == BOUNDARY)) {
	    temp_p[i][j] = BOUNDARY;
	    temp_p_flag_array[i][j] = FIXED;
	    temp_p_fix_flag_array[i][j] = INTERNALLY;
	  }
	}

	/* Check x_q (parent-offspring) covariates .  */
	if ((temp_q_fix_flag_array[i][i] == INTERNALLY) ||
	    (temp_q_fix_flag_array[j][j] == INTERNALLY)) {
	  if ((temp_q[i][i] == BOUNDARY) ||
	      (temp_q[j][j] == BOUNDARY)) {
	    temp_q[i][j] = BOUNDARY;
	    temp_q_flag_array[i][j] = FIXED;
	    temp_q_fix_flag_array[i][j] = INTERNALLY;
	  }
	}
      }
    }
  }
  /* restore x_S, x_M1, x_M2, x_c, x_p, x_q from two dimension arrays.  */
  k = 0;
  for ( i = 0; i < total_trait_values; i++) {
    for ( j = i; j < total_trait_values; j++) {
      x_S[k]  =  temp_s[i][j];
      x_M1[k] =  temp_m1[i][j];
      x_M2[k] =  temp_m2[i][j];
      x_c[k]  =  temp_c[i][j];
      x_p[k]  =  temp_p[i][j];
      x_q[k]  =  temp_q[i][j];

      s_flag_array[k]  =  temp_s_flag_array[i][j];
      m1_flag_array[k] =  temp_m1_flag_array[i][j];
      m2_flag_array[k] =  temp_m2_flag_array[i][j];
      c_flag_array[k]  =  temp_c_flag_array[i][j];
      p_flag_array[k]  =  temp_p_flag_array[i][j];
      q_flag_array[k]  =  temp_q_flag_array[i][j];

      s_fix_flag_array[k]  =  temp_s_fix_flag_array[i][j];
      m1_fix_flag_array[k] =  temp_m1_fix_flag_array[i][j];
      m2_fix_flag_array[k] =  temp_m2_fix_flag_array[i][j];
      c_fix_flag_array[k]  =  temp_c_fix_flag_array[i][j];
      p_fix_flag_array[k]  =  temp_p_fix_flag_array[i][j];
      q_fix_flag_array[k]  =  temp_q_fix_flag_array[i][j];

      //      if (datatype == LONGITUDINAL) {
	x_T[k] =  temp_t[i][j];
	t_flag_array[k] = temp_t_flag_array[i][j];
	t_fix_flag_array[k] = temp_t_fix_flag_array[i][j];
	//      }

      k++;
    }
  }
  /* Assign the initial values for all the parameters except for major gene.*/
  if (fix_mg_flag == YES) {
    AssignInitValues();
  }

  for (i = 0; i < total_trait_values; i++) {
    free(temp_s[i]);
    free(temp_m1[i]);
    free(temp_m2[i]);
    free(temp_t[i]);
    free(temp_c[i]);
    free(temp_p[i]);
    free(temp_q[i]);

    free(temp_s_flag_array[i]);
    free(temp_m1_flag_array[i]);
    free(temp_m2_flag_array[i]);
    free(temp_t_flag_array[i]);
    free(temp_c_flag_array[i]);
    free(temp_p_flag_array[i]);
    free(temp_q_flag_array[i]);

    free(temp_s_fix_flag_array[i]);
    free(temp_m1_fix_flag_array[i]);
    free(temp_m2_fix_flag_array[i]);
    free(temp_t_fix_flag_array[i]);
    free(temp_c_fix_flag_array[i]);
    free(temp_p_fix_flag_array[i]);
    free(temp_q_fix_flag_array[i]);
  }

  free(temp_s);
  free(temp_m1);
  free(temp_m2);
  free(temp_t);
  free(temp_c);
  free(temp_p);
  free(temp_q);

  free(temp_s_flag_array);
  free(temp_m1_flag_array);
  free(temp_m2_flag_array);
  free(temp_t_flag_array);
  free(temp_c_flag_array);
  free(temp_p_flag_array);
  free(temp_q_flag_array);

  free(temp_s_fix_flag_array);
  free(temp_m1_fix_flag_array);
  free(temp_m2_fix_flag_array);
  free(temp_t_fix_flag_array);
  free(temp_c_fix_flag_array);
  free(temp_p_fix_flag_array);
  free(temp_q_fix_flag_array);
}


/* If L(fam.likehd) == 0 means this family data marker
   information is wrong, then print the error message.
   Added 10-21-96; Emily.   */
void Calculate::PrintMarkerErrorMg(int curfam_ptr, int nextfam_ptr) {
  char curfam_str[BUF];

  dataIndex = curfam_ptr;

  fp_log << endl << "***ERROR: this family data marker information is wrong:"
	 << endl;
  fp_log << "> ";
  fp_lookuplog << fortArray[dataIndex].familyId << ' '
	       << fortArray[dataIndex].seqId << ' '
	       << fortArray[dataIndex].fatherId << ' '
	       << fortArray[dataIndex].motherId << ' '
	       << fortArray[dataIndex].sex << ' ';
  for(int i=0; i<itraits; i++) {
    fp_lookuplog << fortArray[dataIndex].traits[i] << ' ';
  }
  for(int i=0; i<icovs; i++) {
    fp_lookuplog << fortArray[dataIndex].covariates[i] << ' ';
  }
  cout << endl << endl;

  cout << endl << "***ERROR: this family data marker information is wrong:"
       << endl;
  cout << "> " << curfam_str << endl << endl;
  cout << "Please correct this family data." << endl << endl;
  dataIndex = nextfam_ptr;
}

/*
Build D matrix and the covariate columns d_vec of D matrix.
*/
void Calculate::Build_D_Matrix(int n_fam_member, double **D_mat,
			       double *d_vec[]) {
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

/* Get (itraits*n_fam_member) MuBeta_vec vector by using  i_vec, d_vec vectors.
*/
void Calculate::GetMuBetaVector(double **MuBeta_vec, double *i_vec, double **d_vec, int n_member) {
  int i,j,k,kv;

  kv = 0;

  for (i = 0; i < itraits; i++) {
    for (j = 0; j < (icovs+1); j++) {
      if (j == 0) {
	for(k = 0; k < n_member; k++) {
	  MuBeta_vec[kv][i*n_member+k] = i_vec[k];
	}
      }else {
	for(k = 0; k <n_member ; k++) {
	  MuBeta_vec[kv][i*n_member+k] = d_vec[j-1][k];
	}
      }
      kv++;
    }
  }
}

/* Change d_vec if there is any missing trait value in the family members.
*/
void Calculate::ChangeVec(double *d_vec[],int n_fam_member) {
  int i,k,kk,numofcov;
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
	  d_vec[numofcov][i] = d_vec[numofcov][kk];
	}
	kk++;
      }
      prev_count++;
      n_mb--;
    }
  }
}

/* Print the results of x_S, x_M1, x_M2, x_T to output file.
*/
void Calculate::PrintResults(int debug_flag) {
  int i;

  for (i = 0; i < iinitvcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw(10) << x_S[i] << ' ';
    else {
      fp_out << setw(10) << x_S[i] << ' ';
      fp_s << setw(10) << x_S[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << "debug_flag" << endl;
  else {
    fp_out << endl;
    fp_s << endl;
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw(10) <<  x_M1[i] << ' ';
    else {
      fp_out << setw(10) << x_M1[i] << ' ';
      fp_m1 << setw(10) << x_M1[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << endl;
  else {
    fp_out << endl;
    fp_m1 << endl;
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw(10) << x_M2[i] << ' ';
    else {
      fp_out << setw(10) << x_M2[i] << ' ';
      fp_m2 << setw(10) << x_M2[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << endl;
  else {
    fp_out << endl;
    fp_m2 << endl;
  }

  for (i = 0; i < env_vcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw(10) << x_T[i] << ' ';
    else {
      fp_out << setw(10) << x_T[i] << ' ';
      fp_t << setw(10) << x_T[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << endl;
  else {
    fp_out << endl;
    fp_t << endl;
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw(10) <<  x_c[i] << ' ';
    else {
      fp_out << setw(10) << x_c[i] << ' ';
      fp_sib << setw(10) << x_c[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << endl;
  else {
    fp_out << endl;
    fp_sib << endl;
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw (10) << x_p[i] << ' ';
    else {
      fp_out << setw(10) << x_p[i] << ' ';
      fp_pp << setw(10) << x_p[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << endl;
  else {
    fp_out << endl;
    fp_pp << endl;
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (debug_flag == 1)       /* 1= DEBUG  */
      cout << setw(10) << x_q[i] << ' ';
    else {
      fp_out << setw(10) << x_q[i] << ' ';
      fp_po << setw(10) << x_q[i] << ' ';
    }
  }
  if (debug_flag == 1)           /* 1= DEBUG  */
    cout << endl;
  else {
    fp_out << endl;
    fp_po << endl;
  }
}

/*****************************
Class name: Calculate
Method name: GetDim
Description: Set imubeta_dim to the size of the mubeta matrix and many other
             variables (via function call).
Input: NONE.
Output: NONE.
Side Effects: imubeta_dim and the other variable are global.  Assigning them
              values is a side effect.
Author: Eric Lunde, 6-12-03
******************************/
void Calculate::GetDim() {
  // get mubeta matrix dimension.
  if (mu_flag == ESTIMATE)
    imubeta_dim = (icovs+1)*itraits;
  else
    imubeta_dim = 0;

  Getismtcpq_dim();
}

/* Get the dimension of smtcpq matrix or vector which is for estimating variance
   components.
       ismtcpq_dim = # of estimating variance components.
*/
void Calculate::Getismtcpq_dim() {
  GetSMT_vcnum();
  GetCPQ_vcnum();

  ismtcpq_dim = s_vcnum  + m1_vcnum + m2_vcnum + t_vcnum
    + c_vcnum  + p_vcnum  + q_vcnum;
}

/* Calculate the increasements and new values for the four variables.
   This module calls the following functions:
   dmatrix, dvector, free_dmatrix, free_vector,
   InitMat, InitCPQMat, InitVec, InitMissingTraitflag,
   GetData, Nuclear_fam_marker,
   Fam_matrix_marker,
   ChangeMatrices,
   ChangeVec,
   GetVGZICPQMatrix,
   CopyMatrix,
   GetSmtcpqMatrix,
   InverseOfMatrix,
   DetermOfMatrix,
   Multi_Vectors_Matrix,
   GetFderItem,
   GetExSder,
   AS_SGetExSderMatrix,
   AS_SGetFderVector,
   Multi_Matrix_Vector,
   AS_GetIncValues.
*/
double Calculate::AS_CalculateValues() {
  int    i,j,k,kk,kv,kc,jc, Ni, l;
  int    i_bigdim;
  // int    nfam, N[MAXNUMFAM];
  int    count_smt;
  double sing_flag;              /* singularity flag of the matrix.  */
  double **temp_mat1, **temp_mat2;
  double **V_mat, **V_mat_inv, **temp_mat;
  double ***G_mat, ***Z1_mat, ***Z2_mat, ***I_mat;
  double ***C_mat, ***P_mat, ***Q_mat;
  double **g_mat, **z1_mat, **z2_mat, **i_mat;
  double **c_mat, **p_mat,  **q_mat;
  double **D_mat = NULL, **d_vec = NULL;
  double ***lD_mat = NULL, ***ld_vec = NULL;
  double *ydb_vec,*beta_vec;
  double *i_vec,  **y_vec, *Y_vec, **y_vecPlusMissing;
  double **Exsder_mat, *fder_vec, *inc_vec;
  double ***Smtcpq_mat;
  double **MuBeta_vec;
  double V1, V2, V3, ln_func;
  double *fder_mu,            *exsder_mu;
  double *fder_beta, *fder_mubeta;
  double *fder_s, *fder_m1, *fder_m2, *fder_t;
  double *fder_c, *fder_p,  *fder_q;
  double *exsder_smtcpq,    *exsder_mubeta;
  double *fder_cfmubeta, fder_cfsmtcpq;
  double *exsder_cfmubeta, exsder_cfmusmtcpq, exsder_cfsmtcpq;
  double S_val;
  double c_ln_func, logcf;
  double **drp_vec, *YDrp_vec;
  char str[BUF];

  /* Initialize the size of the pointers.         changed 5-12-98; E.Y.   */
  count_smt = iinitvcnum*3+env_vcnum+iinitvcnum*3;

  G_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z1_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z2_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  I_mat  = (double ***)malloc(env_vcnum  * sizeof(double **));
  C_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  P_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Q_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));

  if (datatype == MULTIVARIATE)
    d_vec  = (double **) malloc(icovs      * sizeof(double *));
  else if (datatype == LONGITUDINAL) {
    lD_mat = (double ***) malloc(total_trait_values  * sizeof(double **));
    ld_vec = (double ***) malloc(total_trait_values  * sizeof(double **));
    for (i = 0; i < total_trait_values; i++)
      ld_vec[i]  = (double **) malloc(icovs      * sizeof(double *));
  }
  y_vec       = (double **) malloc(total_trait_values * sizeof(double *));
  y_vecPlusMissing = (double **) malloc(total_trait_values * sizeof(double *));
  Smtcpq_mat  = (double ***)malloc(count_smt          * sizeof(double **));
  MuBeta_vec  = (double **) malloc(((icovs+1)*itraits)* sizeof(double *));
  fder_mu     = (double *) malloc(itraits             * sizeof(double));
  exsder_mu   = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_beta   = (double *) malloc((itraits*icovs)     * sizeof(double));
  fder_mubeta = (double *) malloc(((icovs+1)*itraits) * sizeof(double));
  fder_s      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_m1     = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_m2     = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_t      = (double *) malloc(env_vcnum           * sizeof(double));
  fder_c      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_p      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_q      = (double *) malloc(iinitvcnum          * sizeof(double));
  exsder_smtcpq = (double *) malloc(itraitcovnum      * sizeof(double));
  exsder_mubeta = (double *) malloc(imubetacovnum     * sizeof(double));
  fder_cfmubeta   = (double *) malloc(((icovs+1)*itraits) * sizeof(double));
  exsder_cfmubeta = (double *) malloc(imubetacovnum       * sizeof(double));
  drp_vec         = (double **)malloc(total_trait_values        * sizeof(double *));
  YDrp_vec        = (double *) malloc(total_trait_values        * sizeof(double));

  sing_flag      = GOOD;          /* set default to non_singularity.  */
  i_bigdim       = imubeta_dim + ismtcpq_dim;

  Exsder_mat     = Lib::dmatrix(0, i_bigdim-1, 0, i_bigdim-1);
  temp_mat2      = Lib::dmatrix(0, i_bigdim-1, 0, i_bigdim-1);
  fder_vec       = Lib::dvector(0, i_bigdim-1);
  inc_vec        = Lib::dvector(0, i_bigdim-1);

  fder_cfsmtcpq     = 0.0;
  exsder_cfmusmtcpq = 0.0;
  exsder_cfsmtcpq   = 0.0;

  for (i = 0; i < itraits; i++)
    fder_mu[i] =0.0;
  for (i = 0; i < total_trait_values; i++)
    drp_vec[i] = Lib::dvector(0,icovs);
  for (i = 0; i < imubeta_dim; i++) {
    fder_mubeta[i]   =0.0;
    fder_cfmubeta[i] =0.0;
  }
  for (i = 0; i < env_vcnum; i++)
    fder_t[i]  =0.0;
  for (i = 0; i < iinitvcnum; i++) {
    fder_s[i]    =0.0;
    fder_m1[i]   =0.0;
    fder_m2[i]   =0.0;
    fder_c[i]    =0.0;
    fder_p[i]    =0.0;
    fder_q[i]    =0.0;
    exsder_mu[i] =0.0;
  }
  for (i = 0; i < itraitcovnum; i++)
    exsder_smtcpq[i]  =0.0;
  for (i = 0; i < imubetacovnum; i++) {
    exsder_mubeta[i]    =0.0;
    exsder_cfmubeta[i]  =0.0;
  }

  ln_func   =0.0;
  logcf     =0.0;
  S_val     = x_S[0] + x_M1[0] + x_M2[0] + x_T[0] + x_c[0] + x_p[0] + x_q[0];

  /*
    cout << "x_mu[0]=" << x_mu[0] << endl;
    cout << "x_S[0]=" << x_S[0] << endl;
    cout << "x_M1[0]=" << x_M1[0] << endl;
    cout << "x_M2[0]=" << x_M2[0] << endl;
    cout << "x_T[0]=" << x_T[0] << endl;
  */

  dataIndex = 0;
  if (need_ibd_flag == YES) {
    fp_loci.clear();
    fp_loci.seekg(0, ios::beg);
    // This next line is to read past the new (7-11-03) title line at the
    // beginning of loci.out - Eric Lunde
    fp_loci.getline(str, BUF);
  }

  // We want to reset the index to the beginning of the ShareRelation array.
  // We only want to do this once, so that is why this command is outside the
  // loop.  It needs to be reset for the Get_gcpq_mat method.
  // - Eric Lunde, 7-28-03
  relationIndex=0;

  /* get the traits and marker loci which need to be calculated
     for each family. */
  for (i = 0; i < nfam; i++) {
    V_mat      = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
    V_mat_inv  = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
    temp_mat   = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
    temp_mat1  = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
    g_mat      = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    z1_mat     = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    z2_mat     = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    i_mat      = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    c_mat      = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    p_mat      = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    q_mat      = Lib::dmatrix(0, N[i]-1, 0, N[i]-1);
    if (datatype == MULTIVARIATE)
      D_mat      = Lib::dmatrix(0, N[i]-1, 0, icovs);
    else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++)
	lD_mat[j]  = Lib::dmatrix(0, N[i]-1, 0, icovs);
    }
    i_vec      = Lib::dvector(0, N[i]-1);
    ydb_vec    = Lib::dvector(0, N[i]-1);
    beta_vec   = Lib::dvector(0, icovs);
    Y_vec      = Lib::dvector(0, itraits*N[i]-1);
    for (j = 0; j < total_trait_values; j++) {
      y_vec[j] = Lib::dvector(0, N[i]-1);
      y_vecPlusMissing[j] = Lib::dvector(0, N[i]-1);
    }
    for (j = 0; j < env_vcnum; j++) {
      I_mat[j] = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
    }
    for (j = 0; j < iinitvcnum; j++) {
      G_mat[j]  = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Z1_mat[j] = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Z2_mat[j] = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
      C_mat[j]  = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
      P_mat[j]  = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Q_mat[j]  = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);
    }
    for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++)
      Smtcpq_mat[j] = Lib::dmatrix(0, itraits*N[i]-1, 0, itraits*N[i]-1);

    for (j = 0; j < (icovs+1)*itraits; j++)
      MuBeta_vec[j]=Lib::dvector(0, itraits*N[i]-1);

    if (datatype == MULTIVARIATE) {
      for (k = 0; k < icovs; k++) {
	d_vec[k]  = Lib::dvector(0, N[i]-1);
      }
    }else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++) {
	for (k = 0; k < icovs; k++) {
	  ld_vec[j][k] = Lib::dvector(0, N[i]-1);
	}
      }
    }

    InitMat(N[i], G_mat, Z1_mat, Z2_mat, I_mat, g_mat, z1_mat, z2_mat, i_mat);
    InitCPQMat(N[i], C_mat, P_mat,  Q_mat,  c_mat, p_mat, q_mat);
    InitVec(N[i], i_vec, ydb_vec, MuBeta_vec);
    InitMissingTraitflag(N[i]);

    /* Get initial G_VAL 'fam' value here and it may be changed within
       this function.                                             */
    GetData(N[i]);

    if (datatype == MULTIVARIATE) {
      Build_D_Matrix(N[i], D_mat, d_vec);
    }else if (datatype == LONGITUDINAL) {
      Build_lD_Matrices(N[i], lD_mat, ld_vec);
    }

    /* Get value of drp_vec which is the proband vector;added 2-11-97.*/
    for (j = 0; j < total_trait_values; j++) {
      drp_vec[j][0] = 1.0;
      for (k = 0; k < icovs; k++) {
	if (datatype == MULTIVARIATE){
	  drp_vec[j][k+1] = d_vec[k][fam.proband[0].prob_mem];
	}else if (datatype == LONGITUDINAL) {
	  drp_vec[j][k+1] = ld_vec[j][k][fam.proband[0].prob_mem];
	}
      }
    }

    for (j = 0; j < total_trait_values; j++) {
      if (datatype == MULTIVARIATE) {
	Build_beta_Vector(j, beta_vec);
        Lib::Multi_Matrix_Vector(ydb_vec, N[i], icovs+1, D_mat, beta_vec);
      }else if (datatype == LONGITUDINAL) {
	Build_beta_Vector(0, beta_vec);
	Lib::Multi_Matrix_Vector(ydb_vec, N[i], icovs+1, lD_mat[j], beta_vec);
      }
      for(k = 0,kk = 0; k < N[i]; k++) {
	if (fam.person[k].missing_traitflag != 1) {
	  y_vec[j][kk] = fam.person[k].trait[j] - ydb_vec[k];
	  kk++;
	  y_vecPlusMissing[j][k] = fam.person[k].trait[j] - ydb_vec[k];
	}else {
	  y_vecPlusMissing[j][k] = -999;
	}
      }
      YDrp_vec[j] = fam.proband[0].trait_p -
	Lib::Multi_Vectors(icovs+1, drp_vec[j], beta_vec);
    }

    fam.n_sib    = N[i] - 2;
    fam.sib_pair = fam.n_sib*(fam.n_sib-1)/2;

    /* changed 4/22/97;     E.Y.    */
    if (need_ibd_flag == YES)
      Get_z_mat(N[i],    z1_mat, z2_mat);

    Get_gcpq_mat(N[i], g_mat,  c_mat, p_mat, q_mat);

    /* After ChangeMatrices function, we got 'Ni' for new family member number.
       So we change N[i] to Ni from now until before free_dmatrix.  */
    // N[i] represents the number of people in a family, whereas Ni (after
    // returning from ChangeMatrices) represents the number of people in a
    // family with trait values.  I don't know yet if it means the number of
    // people with values for all traits or just the number of people in a
    // family without missing any trait values.  Eric Lunde 01-07-04
    ChangeMatrices(g_mat,z1_mat,z2_mat,i_mat,c_mat,p_mat,q_mat,N[i],&Ni);

    for (j = 0; j < total_trait_values; j++) {
      for(k = 0; k < Ni; k++) {
	Y_vec[j*Ni+k] = y_vec[j][k];
      }
    }

    GetVGZICPQMatrix(V_mat, G_mat, Z1_mat, Z2_mat, I_mat, C_mat, P_mat, Q_mat,
		     g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat, Ni);

    saveVMatrix(savedVMatrix, V_mat, i, total_trait_values, Ni);

    double **traitValuedMembersYBetaDiff
      = Lib::dmatrix(0, total_trait_values - 1, 0, Ni - 1);
    int traitValuedMembersYBetaDiffIndex = 0;
    for(k = 0; k < N[i]; k++) {
      if(fam.person[k].missing_traitflag != 1) {
	for(l = 0; l < total_trait_values; l++) {
	  traitValuedMembersYBetaDiff[l][traitValuedMembersYBetaDiffIndex]
	    = y_vecPlusMissing[l][k];
	}
	traitValuedMembersYBetaDiffIndex++;
      }
    }
    saveYBetaDiff(savedYBetaDiff, traitValuedMembersYBetaDiff, 
		  /*y_vecPlusMissing,*/ i, total_trait_values,
		  Ni);
    Lib::free_dmatrix(traitValuedMembersYBetaDiff, 0, total_trait_values - 1,
		      0, Ni - 1);

    // Here is were we finally assign values to the familyNonMissingSizes
    // array. - EML 09-30-2004
    savedTraitValuedMembers[i] = Ni;

    Lib::CopyMatrix(temp_mat,itraits*Ni, V_mat);
    GetSmtcpqMatrix(Smtcpq_mat, G_mat, Z1_mat, Z2_mat, I_mat, C_mat, P_mat,
		    Q_mat, Ni);
    if (datatype == MULTIVARIATE)
      GetMuBetaVector(MuBeta_vec,  i_vec, d_vec, Ni);
    else if (datatype == LONGITUDINAL)
      LGetMuBetaVector(MuBeta_vec,  i_vec, ld_vec, Ni);

    Lib::CopyMatrix(temp_mat1,itraits*Ni, V_mat);
    sing_flag = Lib::DetermOfMatrix(itraits*Ni, temp_mat1);
    if (sing_flag == BAD) {
      char mat[] = "V_mat";
      char routine[] = "AS_CalculateValues";
      PrintErrMsg(mat, routine, itraits*Ni, V_mat);
    }

    Lib::InverseOfMatrix(V_mat_inv, itraits*Ni, V_mat);

    V1 = Lib::DetermOfMatrix(itraits*Ni, temp_mat);
    V2 = log(V1);
    V3 = Lib::Multi_Vectors_Matrix(itraits*Ni, Y_vec, V_mat_inv, Y_vec);
    /*
      #ifdef DEBUG
      printf("V1 = %lf\n", V1);
      printf("V2 = %lf\n", V2);
      printf("V3 = %lf\n", V3);
      printf("ln_func = %lf\n", ln_func);
      #endif
    */
    ln_func +=(-1.0/2.0)* (V2 + V3);

    if(ln_func > DBL_MAX || ln_func < -DBL_MAX) {
      /// maybe this is not by itself a bad thing. It only causes the program
      // to crash later on with it tried to do sqrt of a negative number
      // (when the [major, output to the screen] iteration completes).
      // It may be inappropriate to exit the program now based only on this
      // info. - Eric Lunde, 2005-09-07
      if(ln_func > DBL_MAX) {
	cerr << "ln_func > " << DBL_MAX << endl;
      } else if(ln_func < -DBL_MAX) {
	cerr << "ln_func < " << -DBL_MAX << endl;
      } else {
	cerr << "This should never be printed to the screen" << endl;
      }
      cerr << "total_trait_values = " << total_trait_values << endl;
      cerr << "Ni = " << Ni << endl;
      cerr << "N[" << i << "] = " << N[i] << endl;
      Lib::PrintOneMatrix(total_trait_values*Ni, total_trait_values*Ni,
			  temp_mat);
      
      cerr << "V1 = " << V1 << endl;
      cerr << "V2 = " << V2 << endl;
      cerr << "V3 = " << V3 << endl;
      cerr << "ln_func = " << ln_func << endl;
      PROBLEM "ln_func is not a number\nCalculate.cpp key 4034\n"
	RECOVER(NULL_ENTRY);
    }

    familyLoglikelihoods[i] = (-1.0/2.0)* (V2 + V3);

    logcf   +=(-1.0/2.0)*log(S_val) +
      (-1.0/2.0)*(YDrp_vec[0]*YDrp_vec[0])/S_val;

    /* This code is commented because, all previous references to fam_lik
       or fam.lik were removed.  These calculations may be desired later, so
       I don't want to delete the formulas.
       fam_lik << "fam# = " << i
       << "  loglikelihood = "
       << setprecision(3) << (-1.0/2.0)*(V2+V3)
       << "  logcf = "
       << (-1.0/2.0)*log(S_val)+(-1.0/2.0)*(YDrp_vec[0]*YDrp_vec[0])/S_val
       << " diff="
       << (-1.0/2.0)*(V2+V3)-(-1.0/2.0)*log(S_val)-(-1.0/2.0)*(YDrp_vec[0]*YDrp_vec[0])/S_val << endl;
    */

    kc = 0;
    for (j = 0; j < itraits; j++) {
      if (mu_flag == ESTIMATE) {
	for (jc = 0; jc < (icovs+1); jc++) {
	  fder_mubeta[kc]  +=
	    Lib::Multi_Vectors_Matrix(itraits*Ni,MuBeta_vec[kc],V_mat_inv,Y_vec);
	  fder_cfmubeta[kc]+= (YDrp_vec[0])*drp_vec[j][jc]/S_val;
          kc++;
        }
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (s_flag_array[j] == ESTIMATE) {
	fder_s[kv] +=
	  Lib::GetFderItem(1,itraits*Ni,V_mat_inv,G_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,G_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (m1_flag_array[j] == ESTIMATE) {
	fder_m1[kv] +=
	  Lib::GetFderItem(1,itraits*Ni,V_mat_inv,Z1_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,Z1_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (m2_flag_array[j] == ESTIMATE) {
	fder_m2[kv] +=
	  Lib::GetFderItem(1,itraits*Ni,V_mat_inv,Z2_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,Z2_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < env_vcnum; j++) {
      if (t_flag_array[j] == ESTIMATE) {
	fder_t[kv]+= Lib::GetFderItem(1,itraits*Ni,V_mat_inv,I_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,I_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (c_flag_array[j] == ESTIMATE) {
	fder_c[kv] +=
	  Lib::GetFderItem(1,itraits*Ni,V_mat_inv,C_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,C_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (p_flag_array[j] == ESTIMATE) {
	fder_p[kv] +=
	  Lib::GetFderItem(1,itraits*Ni,V_mat_inv,P_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,P_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (q_flag_array[j] == ESTIMATE) {
	fder_q[kv] +=
	  Lib::GetFderItem(1,itraits*Ni,V_mat_inv,Q_mat[j],Y_vec)
	  + Lib::GetFderItem(2,itraits*Ni,V_mat_inv,Q_mat[j],Y_vec);
	kv++;
      }
    }
    fder_cfsmtcpq += - 1.0/(2.0*S_val)
      + YDrp_vec[0] * YDrp_vec[0] / (2.0*S_val*S_val);

    kc = 0;
    if (mu_flag == ESTIMATE) {
      for (j = 0; j < ((icovs+1)*itraits); j++) {
	for (k = j; k < ((icovs+1)*itraits); k++) {
	  exsder_mubeta[kc]   +=
	    Lib::Multi_Vectors_Matrix(itraits * Ni, MuBeta_vec[j], V_mat_inv,
				      MuBeta_vec[k]);
	  exsder_cfmubeta[kc] += drp_vec[0][j]*drp_vec[0][k]/S_val;
	  kc++;
	}
      }
    }

    kv = 0;
    for (j = 0; j < ismtcpq_dim; j++) {
      for (k = j; k < ismtcpq_dim; k++) {
	exsder_smtcpq[kv] +=
	  Lib::GetExSder(itraits*Ni,V_mat_inv,Smtcpq_mat[j],Smtcpq_mat[k]);
	kv++;
      }
    }
    exsder_cfsmtcpq += 1.0/(2.0*S_val*S_val);

    Lib::free_dmatrix(V_mat,      0, itraits*N[i]-1, 0, itraits*N[i]-1);
    Lib::free_dmatrix(V_mat_inv,  0, itraits*N[i]-1, 0, itraits*N[i]-1);
    Lib::free_dmatrix(temp_mat,   0, itraits*N[i]-1, 0, itraits*N[i]-1);
    Lib::free_dmatrix(temp_mat1,  0, itraits*N[i]-1, 0, itraits*N[i]-1);
    Lib::free_dmatrix(g_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z1_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z2_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(i_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(c_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(p_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(q_mat,      0, N[i]-1, 0, N[i]-1);
    if (datatype == MULTIVARIATE)
      Lib::free_dmatrix(D_mat,      0, N[i]-1, 0, icovs);
    else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++)
	Lib::free_dmatrix(lD_mat[j],  0, N[i]-1, 0, icovs);
    }
    Lib::free_dvector(i_vec,      0, N[i]-1);
    Lib::free_dvector(ydb_vec,    0, N[i]-1);
    Lib::free_dvector(beta_vec,   0, icovs );
    Lib::free_dvector(Y_vec,      0, itraits*N[i]-1);
    for (j = 0; j < total_trait_values; j++) {
      Lib::free_dvector(y_vec[j], 0, N[i]-1);
      Lib::free_dvector(y_vecPlusMissing[j], 0, N[i]-1);
    }
    for (j = 0; j < env_vcnum; j++) {
      Lib::free_dmatrix(I_mat[j], 0, itraits*N[i]-1, 0, itraits*N[i]-1);
    }
    for (j = 0; j < iinitvcnum; j++) {
      Lib::free_dmatrix(G_mat[j],  0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Lib::free_dmatrix(Z1_mat[j], 0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Lib::free_dmatrix(Z2_mat[j], 0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Lib::free_dmatrix(C_mat[j],  0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Lib::free_dmatrix(P_mat[j],  0, itraits*N[i]-1, 0, itraits*N[i]-1);
      Lib::free_dmatrix(Q_mat[j],  0, itraits*N[i]-1, 0, itraits*N[i]-1);
    }
    for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++) {
      Lib::free_dmatrix(Smtcpq_mat[j], 0, itraits*N[i]-1, 0, itraits*N[i]-1);
    }
    for (j = 0; j < (icovs+1)*itraits; j++)
      Lib::free_dvector(MuBeta_vec[j], 0, itraits*N[i]-1);

    if (datatype == MULTIVARIATE) {
      for (j = 0; j < icovs; j++) {
	Lib::free_dvector(d_vec[j], 0,N[i]-1);
      }
    }else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++) {
	for (k = 0; k < icovs; k++) {
	  Lib::free_dvector(ld_vec[j][k], 0, N[i]-1);
	}
      }
    }
  }

  /*
    cout << "after:  exsder_cfmubeta=" exsder_cfmubeta[0] << endl;
    cout << "after:  exsder_cfmusmtcpq=" << exsder_cfmusmtcpq << endl;
    cout << "after:  exsder_cfsmtcpq=" <<  exsder_cfsmtcpq << endl;
    cout << "after:  fder_cfmubeta=" << fder_cfmubeta[0] << endl;
    cout << "after:  fder_cfsmtcpq=" << fder_cfsmtcpq << endl;
  */

  AS_SGetExSderMatrix(Exsder_mat, exsder_mubeta, exsder_smtcpq,
		      exsder_cfmubeta, exsder_cfmusmtcpq, exsder_cfsmtcpq);

  AS_SGetFderVector(fder_vec, fder_mubeta, fder_s, fder_m1, fder_m2, fder_t,
		    fder_cfmubeta, fder_cfsmtcpq, fder_c, fder_p, fder_q);

  Lib::CopyMatrix(temp_mat2, i_bigdim, Exsder_mat);
  sing_flag = Lib::DetermOfMatrix(i_bigdim, temp_mat2);
  if (sing_flag == BAD) {
    char mat[] = "Exsder_mat";
    char routine[] = "AS_CalculateValues";
    PrintErrMsg(mat, routine,i_bigdim, Exsder_mat);
  }

  /* Get initial G_VAL 'inv_Exsder_mat' value here and NO CHANGE.              */
  Lib::InverseOfMatrix(inv_Exsder_mat, i_bigdim, Exsder_mat);

  if (final_flag == YES) {
    cout << endl << "This is the inverse of Expected Second Derivative Matrix"
	 << endl;
    cout << "for mu, covariates and variance components:" << endl << endl;
    Lib::PrintOneMatrix(i_bigdim, i_bigdim, inv_Exsder_mat);
  }
  Lib::Multi_Matrix_Vector(inc_vec,    i_bigdim, i_bigdim,
			   inv_Exsder_mat, fder_vec);

  /* Change G_VAL 'new_inc_mu', 'new_inc_beta', 'new_inc_s', 'new_inc_m1', 'new_inc_m2',
     'new_inc_t', 'new_inc_c', 'new_inc_p',   'new_inc_q'  values here.    */
  //    AS_GetIncValues(inc_vec,   new_inc_mu, new_inc_beta, new_inc_s, new_inc_m1, new_inc_m2,
  //		    new_inc_t, new_inc_c,  new_inc_p,    new_inc_q);

  AS_GetIncValues(inc_vec);

  Lib::free_dmatrix(Exsder_mat, 0, i_bigdim-1, 0, i_bigdim-1);
  Lib::free_dmatrix(temp_mat2,  0, i_bigdim-1, 0, i_bigdim-1);
  Lib::free_dvector(fder_vec,   0, i_bigdim-1);
  Lib::free_dvector(inc_vec,    0, i_bigdim-1);
  // itraits was commented out on 3/11/2004 by Eric Lunde,
  // when it is allocated, drp_vec is with total_trait_values slots
  for (i = 0; i < total_trait_values/*itraits*/; i++)
    Lib::free_dvector(drp_vec[i], 0, icovs);

  free(G_mat);
  free(Z1_mat);
  free(Z2_mat);
  free(I_mat);
  free(C_mat);
  free(P_mat);
  free(Q_mat);

  free(y_vec);
  free(y_vecPlusMissing);
  free(Smtcpq_mat);
  free(MuBeta_vec);
  free(fder_mu);
  free(exsder_mu);
  free(fder_beta);
  free(fder_mubeta);
  free(fder_s);
  free(fder_m1);
  free(fder_m2);
  free(fder_t);
  free(fder_c);
  free(fder_p);
  free(fder_q);
  free(exsder_smtcpq);
  free(exsder_mubeta);
  free(fder_cfmubeta);
  free(exsder_cfmubeta);
  free(YDrp_vec);
  free(drp_vec);

  if (datatype == MULTIVARIATE)
    free(d_vec);
  else if (datatype == LONGITUDINAL) {
    for (i = 0; i < total_trait_values; i++)
      free(ld_vec[i]);
    free(ld_vec);
    free(lD_mat);
  }

  c_ln_func = ln_func - logcf;

  return c_ln_func;
}

/* Get the matrix of expected second derivatives with ascertainment correction.
*/
void Calculate::AS_SGetExSderMatrix(double **mat, double exsder_mubeta[],
				    double exsder_smtcpq[],
				    double exsder_cfmubeta[],
				    double exsder_cfmusmtcpq,
				    double exsder_cfsmtcpq) {
  int i,j,ii,jj,kv;
  int mubeta_dim = 0, i_bigdim;

  if (mu_flag == ESTIMATE)
    mubeta_dim = imubeta_dim;
  else if (mu_flag == FIXED)
    mubeta_dim = 0;

  i_bigdim = mubeta_dim + ismtcpq_dim;

  /* initiate mat matrix.    */
  for (i = 0; i < i_bigdim; i++) {
    for (j = 0; j < i_bigdim; j++) {
      mat[i][j] = 0.0;
    }
  }

  kv = 0;
  for (i = 0; i < mubeta_dim; i++) {
    for ( j = i; j < mubeta_dim; j++) {
      mat[i][j] = exsder_mubeta[kv] - exsder_cfmubeta[kv];
      kv++;
    }
  }
  for (i = 0; i < mubeta_dim; i++) {
    for ( j = mubeta_dim; j < i_bigdim; j++) {
      mat[i][j] -=  exsder_cfmusmtcpq;
    }
  }
  kv = 0;
  for (i = mubeta_dim; i < i_bigdim; i++) {
    for (j = i; j < i_bigdim; j++) {
      mat[i][j] = exsder_smtcpq[kv] - exsder_cfsmtcpq;
      kv++;
    }
  }

  /* Get the lower corner of mat because symmetry.  */
  for (i = 0; i < i_bigdim; i++) {
    for (j = 0; j < i_bigdim; j++) {
      if (j > i) {
	ii = i;
	jj = j;
	mat[jj][ii] = mat[i][j];
      }
    }
  }
}

/* Get the vector for the first derivatives with ascertainment correction.
*/
void Calculate::AS_SGetFderVector(double *vec, double fder_mubeta[], double fder_s[],
				  double fder_m1[], double fder_m2[],double fder_t[],
				  double fder_cfmubeta[], double fder_cfsmtcpq,
				  double fder_c[], double fder_p[], double fder_q[]) {
  int i,kv;

  kv = 0;
  if (mu_flag == ESTIMATE) {
    for (i = 0; i < imubeta_dim; i++) {
      vec[kv] = fder_mubeta[i] - fder_cfmubeta[i];
#ifdef DEBUG
      printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_mubeta[i],fder_cfmubeta[i]);
#endif
      kv++;
    }
  }
  for (i = 0; i < s_vcnum; i++) {
    vec[kv] = fder_s[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_s[i], fder_cfsmtcpq);
#endif
    kv++;
  }
  for (i = 0; i < m1_vcnum; i++) {
    vec[kv] = fder_m1[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_m1[i] , fder_cfsmtcpq);
#endif
    kv++;
  }
  for (i = 0; i < m2_vcnum; i++) {
    vec[kv] = fder_m2[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_m2[i] , fder_cfsmtcpq);
#endif
    kv++;
  }
  for (i = 0; i < t_vcnum; i++) {
    vec[kv] = fder_t[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_t[i], fder_cfsmtcpq);
#endif
    kv++;
  }
  for (i = 0; i < c_vcnum; i++) {
    vec[kv] = fder_c[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_c[i] , fder_cfsmtcpq);
#endif
    kv++;
  }
  for (i = 0; i < p_vcnum; i++) {
    vec[kv] = fder_p[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_p[i] , fder_cfsmtcpq);
#endif
    kv++;
  }
  for (i = 0; i < q_vcnum; i++) {
    vec[kv] = fder_q[i] - fder_cfsmtcpq;
#ifdef DEBUG
    printf("vec[%d](%lf)=%lf-%lf\n",kv, vec[kv], fder_q[i] , fder_cfsmtcpq);
#endif
    kv++;
  }
}

/* Get the increased values for MU, SIGMA, MG, TAU, C, P, Q with ascertainment correction.
*/
void Calculate::AS_GetIncValues(double *vec) {
  /*
    void Calculate::AS_GetIncValues(double *vec, double new_inc_mu[], double new_inc_beta[],
    double new_inc_s[],  double new_inc_m1[], double new_inc_m2[],
    double new_inc_t[],  double new_inc_c[],  double new_inc_p[],
    double new_inc_q[]) {
  */
  int i,j,kb,kv;

  kb = kv = 0;
  if (mu_flag == ESTIMATE) {
    for (i = 0; i < itraits; i++) {
      for ( j = 0; j < (icovs+1); j++) {
	if ( j == 0)
	  new_inc_mu[i] = vec[kv];
	else  {
	  new_inc_beta[kb] = vec[kv];
	  kb++;
	}
#ifdef DEBUG
	cout << "new_inc_mu[" << i << "]=" << new_inc_mu[i] << endl;
#endif
	kv++;
      }
    }
  }
  for (i = 0; i < iinitvcnum; i++) {
    if (s_flag_array[i] == ESTIMATE) {
      new_inc_s[i] = vec[kv];
      kv++;
    }
  }
  for (i = 0; i < iinitvcnum; i++) {
    if (m1_flag_array[i] == ESTIMATE) {
      new_inc_m1[i] = vec[kv];
      kv++;
    }
  }
  for (i = 0; i < iinitvcnum; i++) {
    if (m2_flag_array[i] == ESTIMATE) {
      new_inc_m2[i] = vec[kv];
      kv++;
    }
  }
  for (i = 0; i < env_vcnum; i++) {
    if (t_flag_array[i] == ESTIMATE) {
      new_inc_t[i] = vec[kv];
      kv++;
    }
  }
  for (i = 0; i < iinitvcnum; i++) {
    if (c_flag_array[i] == ESTIMATE) {
      new_inc_c[i] = vec[kv];
      kv++;
    }
  }
  for (i = 0; i < iinitvcnum; i++) {
    if (p_flag_array[i] == ESTIMATE) {
      new_inc_p[i] = vec[kv];
      kv++;
    }
  }
  for (i = 0; i < iinitvcnum; i++) {
    if (q_flag_array[i] == ESTIMATE) {
      new_inc_q[i] = vec[kv];
      kv++;
    }
  }
#ifdef DEBUG
  cout << "new_inc_s[0]=" << new_inc_s[0] << endl;
  cout << "new_inc_m1[0]=" << new_inc_m1[0] << endl;
  cout << "new_inc_m2[0]=" << new_inc_m2[0] << endl;
  cout << "new_inc_t[0]=" << new_inc_t[0] << endl;
#endif
}




/******************************************************************************
noasfun.c  ----- the functions for NO_ASERTAINMENT data.
******************************************************************************/

/* Calculate the increasements and new values for the four variables.
   This module calls the following functions:
   dmatrix,
   dvector,
   free_dmatrix,
   free_vector,
   InitMat,
   InitCPQMat,
   InitVec,
   InitMissingTraitflag,
   GetData,
   Nuclear_fam_marker,
   Fam_matrix_marker,
   ChangeMatrices,
   ChangeVec,
   GetVGZICPQMatrix,
   CopyMatrix,
   GetSmtcpqMatrix,
   InverseOfMatrix,
   DetermOfMatrix,
   Multi_Vectors_Matrix,
   GetFderItem,
   GetExSder,
   SGetExSderMatrix,
   SGetFderVector,
   Multi_Matrix_Vector,
   GetIncValues.
*/
double Calculate::NOAS_CalculateValues() {
  int    i,j,k,kk,kv,kc,jc, Ni, l;
  int    mubeta_dim, smtcpq_dim;
  // int    nfam, N[MAXNUMFAM];
  int    count_smt;
  double sing_flag;              /* singularity flag of the matrix.  */
  double **temp_mat1, **temp_mat2 = NULL, **temp_mat3;
  double **V_mat, **V_mat_inv, **temp_mat;
  double ***G_mat, ***Z1_mat, ***Z2_mat, ***I_mat;
  double ***C_mat, ***P_mat, ***Q_mat;
  double **g_mat, **z1_mat, **z2_mat, **i_mat;
  double **c_mat, **p_mat,  **q_mat;
  double **D_mat = NULL, **d_vec = NULL;
  double ***lD_mat = NULL, ***ld_vec = NULL;
  double *ydb_vec,*beta_vec;
  double *i_vec,  **y_vec, *Y_vec, **y_vecPlusMissing;
  double **Exsder_mubeta_mat = NULL;
  double *inc_mubeta_vec = NULL;
  double **Exsder_smtcpq_mat;
  double *fder_smtcpq_vec, *inc_smtcpq_vec;
  double ***Smtcpq_mat;
  double **MuBeta_vec;
  double V1, V2, V3, ln_func;
  double *fder_mubeta, *fder_s, *fder_m1, *fder_m2, *fder_t;
  double *fder_c,      *fder_p, *fder_q;
  double *fder_smtcpq, *exsder_smtcpq,    *exsder_mubeta;
  char str[BUF];

  /* Initialize the size of the pointers.         changed 5-12-98; E.Y.   */
  count_smt = iinitvcnum*2+env_vcnum+iinitvcnum*3;

  //printf("total=%d datatype=%d\n",total_trait_values, datatype);

  G_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z1_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  Z2_mat = (double ***)malloc(iinitvcnum * sizeof(double **));
  I_mat  = (double ***)malloc(env_vcnum  * sizeof(double **));
  C_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  P_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));
  Q_mat  = (double ***)malloc(iinitvcnum * sizeof(double **));

  if (datatype == MULTIVARIATE)
    d_vec  = (double **) malloc(icovs      * sizeof(double *));
  else if (datatype == LONGITUDINAL) {
    lD_mat = (double ***) malloc(total_trait_values  * sizeof(double **));
    ld_vec = (double ***) malloc(total_trait_values  * sizeof(double **));
    for (i = 0; i < total_trait_values; i++)
      ld_vec[i]  = (double **) malloc(icovs      * sizeof(double *));
  }
  y_vec  = (double **) malloc(total_trait_values    * sizeof(double *));
  y_vecPlusMissing = (double **) malloc(total_trait_values    * sizeof(double *));
  Smtcpq_mat  = (double ***)malloc(count_smt          * sizeof(double **));
  MuBeta_vec  = (double **) malloc(((icovs+1)*itraits)* sizeof(double *));
  fder_mubeta = (double *) malloc(((icovs+1)*itraits) * sizeof(double));
  fder_s      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_m1     = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_m2     = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_t      = (double *) malloc(env_vcnum           * sizeof(double));
  fder_c      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_p      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_q      = (double *) malloc(iinitvcnum          * sizeof(double));
  fder_smtcpq = (double *) malloc(count_smt           * sizeof(double));
  exsder_smtcpq = (double *) malloc(itraitcovnum      * sizeof(double));
  exsder_mubeta = (double *) malloc(imubetacovnum     * sizeof(double));

  sing_flag             = GOOD;          /* set default to non_singularity.  */
  mubeta_dim            = imubeta_dim;

  if (mu_flag == ESTIMATE) {
    Exsder_mubeta_mat     = Lib::dmatrix(0,mubeta_dim-1,0,mubeta_dim-1);
    temp_mat2             = Lib::dmatrix(0,mubeta_dim-1,0,mubeta_dim-1);
    inc_mubeta_vec        = Lib::dvector(0,mubeta_dim-1);
  }

  smtcpq_dim            = ismtcpq_dim;

  //printf("smtcpq_dim=%d %d\n",smtcpq_dim,ismtcpq_dim);
  Exsder_smtcpq_mat     = Lib::dmatrix(0,smtcpq_dim-1,0,smtcpq_dim-1);
  temp_mat3             = Lib::dmatrix(0,smtcpq_dim-1,0,smtcpq_dim-1);
  fder_smtcpq_vec       = Lib::dvector(0,smtcpq_dim-1);
  inc_smtcpq_vec        = Lib::dvector(0,smtcpq_dim-1);

  for (i = 0; i < imubeta_dim; i++) {
    fder_mubeta[i] =0.0;
  }
  for (i = 0; i < env_vcnum; i++) {
    fder_t[i]  =0.0;
  }
  for (i = 0; i < iinitvcnum; i++) {
    fder_s[i]   =0.0;
    fder_m1[i]  =0.0;
    fder_m2[i]  =0.0;
    fder_c[i]   =0.0;
    fder_p[i]   =0.0;
    fder_q[i]   =0.0;
  }
  for (i = 0; i < itraitcovnum; i++)
    exsder_smtcpq[i]  =0.0;
  for (i = 0; i < imubetacovnum; i++)
    exsder_mubeta[i]  =0.0;

  ln_func = 0.0;

  dataIndex = 0;
  if (need_ibd_flag == YES) {
    fp_loci.clear();
    fp_loci.seekg(0, ios::beg);
    // This next line is to read past the new (7-11-03) title line at the
    // beginning of loci.out - Eric Lunde
    fp_loci.getline(str, BUF);
  }

  // We want to reset the index to the beginning of the ShareRealation array.
  // We only want to do this once, so that is why this command is outside the
  // loop.  It needs to be reset for the Get_gcpq_mat method.
  // - Eric Lunde, 7-28-03
  relationIndex=0;

  /* get the traits and marker loci which need to be calculated
     for each family. */

  for (i = 0; i < nfam; i++) {
    //cout << "i = " << i << endl;
    // Although the V_mat is allocated for enough space for
    // (total_trait_values*N[i])^2 elements.  They are not all used.  The
    // V_mat is only populated by (total_trait_values*Ni)^2 elements.
    V_mat      = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    V_mat_inv  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    temp_mat   = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    temp_mat1  = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    g_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    z1_mat     = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    z2_mat     = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    i_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    c_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    p_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);
    q_mat      = Lib::dmatrix(0,N[i]-1,0,N[i]-1);

    if (datatype == MULTIVARIATE)
      D_mat      = Lib::dmatrix(0,N[i]-1,0,icovs);
    else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++) {
	lD_mat[j]  = Lib::dmatrix(0,N[i]-1,0,icovs);
      }
    }
    i_vec      = Lib::dvector(0,N[i]-1);
    ydb_vec    = Lib::dvector(0,N[i]-1);
    beta_vec   = Lib::dvector(0,icovs);
    Y_vec      = Lib::dvector(0,total_trait_values*N[i]-1);
    for (j = 0; j < total_trait_values; j++) {
      y_vec[j] = Lib::dvector(0,N[i]-1);
    }
    for (j = 0; j < total_trait_values; j++) {
      y_vecPlusMissing[j] = Lib::dvector(0,N[i]-1);
    }
    for (j = 0; j < env_vcnum; j++) {
      I_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }
    for (j = 0; j < iinitvcnum; j++) {
      G_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Z1_mat[j]= Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Z2_mat[j]= Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      C_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      P_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
      Q_mat[j] = Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }
    for (j = 0; j < (iinitvcnum*2+env_vcnum+iinitvcnum*3); j++)
      Smtcpq_mat[j]=Lib::dmatrix(0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);

    for (j = 0; j < (icovs+1)*itraits; j++)
      MuBeta_vec[j]=Lib::dvector(0,total_trait_values*N[i]-1);

    if (datatype == MULTIVARIATE) {
      for (j = 0; j < icovs; j++) {
	d_vec[j]  = Lib::dvector(0,N[i]-1);
      }
    }else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++) {
	for (k = 0; k < icovs; k++) {
      	  ld_vec[j][k] = Lib::dvector(0,N[i]-1);
	}
      }
    }

    InitMat   (N[i], G_mat,  Z1_mat, Z2_mat, I_mat, g_mat, z1_mat, z2_mat, i_mat);
    InitCPQMat(N[i], C_mat,  P_mat,  Q_mat,  c_mat, p_mat, q_mat);
    InitVec   (N[i], i_vec, ydb_vec, MuBeta_vec);
    InitMissingTraitflag(N[i]);

    /* Get initial G_VAL 'fam' value here and it may be changed within
       this function.                                             */
    GetData(N[i]);

    if (datatype == MULTIVARIATE)
      Build_D_Matrix(N[i], D_mat, d_vec);
    else if (datatype == LONGITUDINAL)
      Build_lD_Matrices(N[i], lD_mat, ld_vec);

    for (j = 0; j < total_trait_values; j++) {
      if (datatype == MULTIVARIATE) {
	Build_beta_Vector(j, beta_vec);
	// This multiplication may be what mariza is looking for.
	Lib::Multi_Matrix_Vector(ydb_vec, N[i], icovs+1, D_mat, beta_vec);
      }else if (datatype == LONGITUDINAL) {
	/* assume have only one trait for longitudinal data.   6-98; E.Y. */
	Build_beta_Vector(0, beta_vec);
	Lib::Multi_Matrix_Vector(ydb_vec, N[i], icovs+1, lD_mat[j], beta_vec);
      }
      for(k = 0,kk = 0; k < N[i]; k++) {
	if (fam.person[k].missing_traitflag != 1) {
	  y_vec[j][kk] = fam.person[k].trait[j] - ydb_vec[k];
	  kk++;
	  y_vecPlusMissing[j][k] = fam.person[k].trait[j] - ydb_vec[k];
	}else {
	  y_vecPlusMissing[j][k] = -999;
	}
      }
    }

    fam.n_sib    = N[i] - 2;
    fam.sib_pair = fam.n_sib*(fam.n_sib-1)/2;

    /* changed 4/22/97;    E.Y.    */
    if (need_ibd_flag == YES) {
      Get_z_mat(N[i],    z1_mat, z2_mat);
    }

    Get_gcpq_mat(N[i], g_mat,  c_mat, p_mat, q_mat);

    /* After ChangeMatrices function, we got 'Ni' for new family member number.
       So we change N[i] to Ni from now until before free_dmatrix.  */
    // N[i] represents the number of people in a family, whereas Ni (after
    // returning from ChangeMatrices) represents the number of people in a
    // family with trait values.  I don't know yet if it means the number of
    // people with values for all traits or just the number of people in a
    // family without missing any trait values.  Eric Lunde 01-07-04
    ChangeMatrices(g_mat, z1_mat, z2_mat, i_mat, c_mat, p_mat, q_mat, N[i],
		   &Ni);

    if (datatype == MULTIVARIATE)
      ChangeVec(d_vec, N[i]);
    else if (datatype == LONGITUDINAL)
      LChangeVec(ld_vec, N[i]);

    for (j = 0; j < total_trait_values; j++) {
      for(k = 0; k < Ni; k++) {
	Y_vec[j*Ni+k] = y_vec[j][k];
      }
    }

    GetVGZICPQMatrix(V_mat,G_mat,Z1_mat,Z2_mat,I_mat,C_mat,P_mat,Q_mat,
		     g_mat,z1_mat,z2_mat,i_mat,c_mat,p_mat,q_mat,Ni);

    /*
    saveVMatrix(savedVMatrix, V_mat, i, total_trait_values, N[i]);
    saveYBetaDiff(savedYBetaDiff, y_vecPlusMissing, i, total_trait_values,
    N[i]);
    */
    saveVMatrix(savedVMatrix, V_mat, i, total_trait_values, Ni);

    double **traitValuedMembersYBetaDiff
      = Lib::dmatrix(0, total_trait_values - 1, 0, Ni - 1);
    int traitValuedMembersYBetaDiffIndex = 0;
    for(k = 0; k < N[i]; k++) {
      if(fam.person[k].missing_traitflag != 1) {
	for(l = 0; l < total_trait_values; l++) {
	  traitValuedMembersYBetaDiff[l][traitValuedMembersYBetaDiffIndex]
	    = y_vecPlusMissing[l][k];
	}
	traitValuedMembersYBetaDiffIndex++;
      }
    }
    saveYBetaDiff(savedYBetaDiff, traitValuedMembersYBetaDiff,
		  /*y_vecPlusMissing,*/ i, total_trait_values,
		  Ni);
    Lib::free_dmatrix(traitValuedMembersYBetaDiff, 0, total_trait_values - 1,
		      0, Ni - 1);

    // Here is were we finally assign values to the familyNonMissingSizes
    // array. - EML 09-30-2004
    savedTraitValuedMembers[i] = Ni;

    Lib::CopyMatrix(temp_mat,total_trait_values*Ni, V_mat);
    GetSmtcpqMatrix(Smtcpq_mat, G_mat, Z1_mat, Z2_mat, I_mat, C_mat, P_mat, Q_mat, Ni);
    if(datatype == MULTIVARIATE)
      GetMuBetaVector(MuBeta_vec,  i_vec, d_vec, Ni);
    else if (datatype == LONGITUDINAL)
      LGetMuBetaVector(MuBeta_vec,  i_vec, ld_vec, Ni);

    Lib::CopyMatrix(temp_mat1, total_trait_values*Ni, V_mat);
    sing_flag = Lib::DetermOfMatrix(total_trait_values*Ni, temp_mat1);
    if (sing_flag == BAD) {
      char mat[] = "V_mat";
      char routine [] = "NOAS_CalculateValues";
      PrintErrMsg(mat,routine,total_trait_values*Ni, V_mat);
    }
    Lib::InverseOfMatrix(V_mat_inv, total_trait_values*Ni, V_mat);
    
    V1 = Lib::DetermOfMatrix(total_trait_values*Ni, temp_mat);
    V2 = log(V1);
    V3 = Lib::Multi_Vectors_Matrix(total_trait_values*Ni, Y_vec, V_mat_inv, Y_vec);

    ln_func +=(-1.0/2.0)* (V2 + V3) ;

    if(ln_func > DBL_MAX || ln_func < -DBL_MAX) {
      /// maybe this is not by itself a bad thing. It only causes the program
      // to crash later on with it tried to do sqrt of a negative number
      // (when the [major, output to the screen] iteration completes).
      // It may be inappropriate to exit the program now based only on this
      // info. - Eric Lunde, 2005-09-07
      if(ln_func > DBL_MAX) {
	cerr << "ln_func > " << DBL_MAX << endl;
      } else if(ln_func < -DBL_MAX) {
	cerr << "ln_func < " << -DBL_MAX << endl;
      } else {
	cerr << "This should never be printed to the screen" << endl;
      }
      cerr << "total_trait_values = " << total_trait_values << endl;
      cerr << "Ni = " << Ni << endl;
      cerr << "N[" << i << "] = " << N[i] << endl;
      Lib::PrintOneMatrix(total_trait_values*Ni, total_trait_values*Ni,
			  temp_mat);
      
      cerr << "V1 = " << V1 << endl;
      cerr << "V2 = " << V2 << endl;
      cerr << "V3 = " << V3 << endl;
      cerr << "ln_func = " << ln_func << endl;
      PROBLEM "ln_func is not a number\nCalculate.cpp key 4859\n"
	RECOVER(NULL_ENTRY);
    }
    //cout << "ln_func = " << ln_func << endl;
    familyLoglikelihoods[i] = (-1.0/2.0)* (V2 + V3);

    kc = 0;
    for (j = 0; j < itraits; j++) {
      if (mu_flag == ESTIMATE) {
	for (jc = 0; jc < (icovs+1); jc++) {
	  fder_mubeta[kc] +=
	    Lib::Multi_Vectors_Matrix(total_trait_values*Ni,MuBeta_vec[kc],V_mat_inv,Y_vec);
	  kc++;
	}
      }
    }

    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (s_flag_array[j] == ESTIMATE) {
	fder_s[kv] +=
	  Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,G_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,G_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (m1_flag_array[j] == ESTIMATE) {
	fder_m1[kv] +=
	  Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,Z1_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,Z1_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (m2_flag_array[j] == ESTIMATE) {
	fder_m2[kv] +=
	  Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,Z2_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,Z2_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < env_vcnum; j++) {
      if (t_flag_array[j] == ESTIMATE) {
	fder_t[kv]+= Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,I_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,I_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (c_flag_array[j] == ESTIMATE) {
	fder_c[kv] +=
	  Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,C_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,C_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (p_flag_array[j] == ESTIMATE) {
	fder_p[kv] +=
	  Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,P_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,P_mat[j],Y_vec);
	kv++;
      }
    }
    kv = 0;
    for (j = 0; j < iinitvcnum; j++) {
      if (q_flag_array[j] == ESTIMATE) {
	fder_q[kv] +=
	  Lib::GetFderItem(1,total_trait_values*Ni,V_mat_inv,Q_mat[j],Y_vec)
	  + Lib::GetFderItem(2,total_trait_values*Ni,V_mat_inv,Q_mat[j],Y_vec);
	kv++;
      }
    }

    kc = 0;

    if (mu_flag == ESTIMATE) {
      for (j = 0; j < ((icovs+1)*itraits); j++) {
	for (k = j; k < ((icovs+1)*itraits); k++) {
	  exsder_mubeta[kc] +=
	    Lib::Multi_Vectors_Matrix(total_trait_values*Ni, MuBeta_vec[j],
				      V_mat_inv, MuBeta_vec[k]);
	  kc++;
	}
      }
    }

    kv = 0;
    for (j = 0; j < ismtcpq_dim; j++) {
      for (k = j; k < ismtcpq_dim; k++) {
	exsder_smtcpq[kv] +=
	  Lib::GetExSder(total_trait_values*Ni,V_mat_inv,Smtcpq_mat[j],Smtcpq_mat[k]);
	kv++;
      }
    }

    Lib::free_dmatrix(V_mat,      0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(V_mat_inv , 0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(temp_mat,   0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(temp_mat1,  0, total_trait_values*N[i]-1, 0, total_trait_values*N[i]-1);
    Lib::free_dmatrix(g_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z1_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(z2_mat,     0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(i_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(c_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(p_mat,      0, N[i]-1, 0, N[i]-1);
    Lib::free_dmatrix(q_mat,      0, N[i]-1, 0, N[i]-1);
    if (datatype == MULTIVARIATE)
      Lib::free_dmatrix(D_mat,      0, N[i]-1, 0, icovs );
    else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++)
	Lib::free_dmatrix(lD_mat[j],  0,N[i]-1,0,icovs);
    }
    Lib::free_dvector(i_vec,      0, N[i]-1);
    Lib::free_dvector(ydb_vec,    0, N[i]-1);
    Lib::free_dvector(beta_vec,   0, icovs );
    Lib::free_dvector(Y_vec,      0, total_trait_values*N[i]-1);
    for (j = 0; j < total_trait_values; j++) {
      Lib::free_dvector(y_vec[j], 0,N[i]-1);
    }
    for (j = 0; j < total_trait_values; j++) {
      Lib::free_dvector(y_vecPlusMissing[j], 0,N[i]-1);
    }
    for (j = 0; j < env_vcnum; j++) {
      Lib::free_dmatrix(I_mat[j], 0,total_trait_values*N[i]-1,0,total_trait_values*N[i]-1);
    }
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
    for (j = 0; j < (icovs+1)*itraits; j++)
      Lib::free_dvector(MuBeta_vec[j], 0,total_trait_values*N[i]-1);

    if (datatype == MULTIVARIATE) {
      for (j = 0; j < icovs; j++) {
	Lib::free_dvector(d_vec[j], 0,N[i]-1);
      }
    }
    else if (datatype == LONGITUDINAL) {
      for (j = 0; j < total_trait_values; j++) {
	for (k = 0; k < icovs; k++) {
	  Lib::free_dvector(ld_vec[j][k], 0,N[i]-1);
	}
      }
    }
  }

  if (mu_flag == ESTIMATE) {
    SGetExSderMatrix(MUBETA, Exsder_mubeta_mat, exsder_mubeta);
    Lib::CopyMatrix(temp_mat2, mubeta_dim, Exsder_mubeta_mat);
    sing_flag = Lib::DetermOfMatrix(mubeta_dim, temp_mat2);

    if (sing_flag == BAD) {
      char mat[] = "Exsder_mubeta_mat";
      char routine[] = "NOAS_CalculateValues";
      PrintErrMsg(mat, routine, mubeta_dim, Exsder_mubeta_mat);
    }

    /* Get initial G_VAL 'inv_Exsder_mubeta_mat' value here and NO CHANGE. */
    Lib::InverseOfMatrix(inv_Exsder_mubeta_mat, mubeta_dim, Exsder_mubeta_mat);

    if (final_flag == YES) {
      cout << endl << "This is the inverse of Expected Second Derivative Matrix";
      cout << "for mu and covariates:" << endl << endl;
      Lib::PrintOneMatrix(mubeta_dim, mubeta_dim, inv_Exsder_mubeta_mat);
    }
  }

  SGetExSderMatrix(SMTCPQ, Exsder_smtcpq_mat, exsder_smtcpq);

  SGetFderVector(fder_smtcpq_vec,fder_s, fder_m1, fder_m2, fder_t,
		 fder_c, fder_p, fder_q);

  Lib::CopyMatrix(temp_mat3, smtcpq_dim, Exsder_smtcpq_mat);
  sing_flag = Lib::DetermOfMatrix(smtcpq_dim, temp_mat3);
  if (sing_flag == BAD) {
    char mat[] = "Exsder_smtcpq_mat";
    char routine[] = "NOAS_CalculateValues";
    PrintErrMsg(mat, routine, smtcpq_dim, Exsder_smtcpq_mat);
  }

  /* Get initial G_VAL 'inv_Exsder_smtcpq_mat' value here and NO CHANGE.    */
  Lib::InverseOfMatrix(inv_Exsder_smtcpq_mat, smtcpq_dim, Exsder_smtcpq_mat);

  if (final_flag == YES) {
    cout << endl <<"This is the inverse of Expected Second Derivative Matrix"
	 << "for variance components:" << endl << endl;
    Lib::PrintOneMatrix(smtcpq_dim, smtcpq_dim, inv_Exsder_smtcpq_mat);
  }

  if (mu_flag == ESTIMATE) {
    Lib::Multi_Matrix_Vector(inc_mubeta_vec,    mubeta_dim, mubeta_dim,
			     inv_Exsder_mubeta_mat, fder_mubeta);
  }

  Lib::Multi_Matrix_Vector(inc_smtcpq_vec,       smtcpq_dim,    smtcpq_dim,
			   inv_Exsder_smtcpq_mat,    fder_smtcpq_vec);

  GetIncValues(inc_mubeta_vec, inc_smtcpq_vec);

  if (mu_flag == ESTIMATE) {
    Lib::free_dmatrix(Exsder_mubeta_mat,   0,mubeta_dim-1,    0,mubeta_dim-1);
    Lib::free_dmatrix(temp_mat2,           0,mubeta_dim-1,    0,mubeta_dim-1);
    Lib::free_dvector(inc_mubeta_vec,      0,mubeta_dim-1);
  }
  Lib::free_dmatrix(Exsder_smtcpq_mat,       0,smtcpq_dim-1,    0,smtcpq_dim-1);
  Lib::free_dmatrix(temp_mat3,               0,smtcpq_dim-1,    0,smtcpq_dim-1);
  Lib::free_dvector(fder_smtcpq_vec,         0,smtcpq_dim-1);
  Lib::free_dvector(inc_smtcpq_vec,          0,smtcpq_dim-1);

  free(G_mat);
  free(Z1_mat);
  free(Z2_mat);
  free(I_mat);
  free(C_mat);
  free(P_mat);
  free(Q_mat);

  free(y_vec);
  free(y_vecPlusMissing);
  free(Smtcpq_mat);
  free(MuBeta_vec);
  free(fder_mubeta);
  free(fder_s);
  free(fder_m1);
  free(fder_m2);
  free(fder_t);
  free(fder_c);
  free(fder_p);
  free(fder_q);
  free(fder_smtcpq);
  free(exsder_smtcpq);
  free(exsder_mubeta);

  if (datatype == MULTIVARIATE)
    free(d_vec);
  else if (datatype == LONGITUDINAL) {
    for (i = 0; i < total_trait_values; i++)
      free(ld_vec[i]);
    free(ld_vec);
    free(lD_mat);
  }

#ifdef DEBUG
  cout << endl << "Loglik = " << setw(20) << setprecision(10) << ln_func
       << endl;
#endif

  return ln_func;
}

/* Get the matrix of expected second derivatives.
*/
void Calculate::SGetExSderMatrix(int mu_smtcpqflag, double **mat,
				 double exsder_array[]) {

  int i,j,kv;
  int n_dim;

  if (mu_smtcpqflag == MUBETA)
    n_dim = (icovs+1)*itraits;
  else
    n_dim = ismtcpq_dim;

  /* initiate mat matrix.    */
  for (i = 0; i < n_dim ; i++) {
    for (j = 0; j < n_dim; j++) {
      mat[i][j] = 0.0;
    }
  }

  kv = 0;
  for (i = 0; i < n_dim; i++) {
    for ( j = i; j < n_dim; j++) {
      mat[i][j] = exsder_array[kv];
      mat[j][i] = mat[i][j];  /* because of symmetry. */
      kv++;
    }
  }
}

/* Get the vector for the first derivatives.
*/
void Calculate::SGetFderVector(double *vec,      double fder_s[],
			       double fder_m1[], double fder_m2[],
			       double fder_t[],  double fder_c[],
			       double fder_p[],  double fder_q[]) {
  int i,kv;

  kv = 0;
  for (i = 0; i < s_vcnum; i++) {
    vec[kv] = fder_s[i];
    kv++;
  }
  for (i = 0; i < m1_vcnum; i++) {
    vec[kv] = fder_m1[i];
    kv++;
  }
  for (i = 0; i < m2_vcnum; i++) {
    vec[kv] = fder_m2[i];
    kv++;
  }
  for (i = 0; i < t_vcnum; i++) {
    vec[kv] = fder_t[i];
    kv++;
  }
  for (i = 0; i < c_vcnum; i++) {
    vec[kv] = fder_c[i];
    kv++;
  }
  for (i = 0; i < p_vcnum; i++) {
    vec[kv] = fder_p[i];
    kv++;
  }
  for (i = 0; i < q_vcnum; i++) {
    vec[kv] = fder_q[i];
    kv++;
  }
}



/* Get the increased values for MU,BETA, SIGMA, MG and TAU.
*/
void Calculate::GetIncValues(double *mubeta_vec,double *smtcpq_vec){
  /*
    void Calculate::GetIncValues(double *mubeta_vec,double *smtcpq_vec, double new_inc_mu[],
    double new_inc_beta[], double new_inc_s[],  double new_inc_m1[],
    double new_inc_m2[],   double new_inc_t[],
    double new_inc_c[],    double new_inc_p[],  double new_inc_q[]) {
  */
  int i,j,kb,kv;

  kb = kv = 0;
  if (mu_flag == ESTIMATE) {
    for (i = 0; i < itraits; i++) {
      for ( j = 0; j < (icovs+1); j++) {
	if ( j == 0)
	  new_inc_mu[i] = mubeta_vec[kv];
	else  {
	  new_inc_beta[kb] = mubeta_vec[kv];
	  kb++;
	}
	kv++;
      }
    }
  }

  kv = 0;
  for (i = 0; i < iinitvcnum; i++) {
    if (s_flag_array[i] == ESTIMATE) {
      new_inc_s[i] = smtcpq_vec[kv];
      kv++;
    }
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (m1_flag_array[i] == ESTIMATE) {
      new_inc_m1[i] = smtcpq_vec[kv];
      kv++;
    }
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (m2_flag_array[i] == ESTIMATE) {
      new_inc_m2[i] = smtcpq_vec[kv];
      kv++;
    }
  }

  for (i = 0; i < env_vcnum; i++) {
    if (t_flag_array[i] == ESTIMATE) {
      new_inc_t[i] = smtcpq_vec[kv];
      kv++;
    }
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (c_flag_array[i] == ESTIMATE) {
      new_inc_c[i] = smtcpq_vec[kv];
      kv++;
    }
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (p_flag_array[i] == ESTIMATE) {
      new_inc_p[i] = smtcpq_vec[kv];
      kv++;
    }
  }

  for (i = 0; i < iinitvcnum; i++) {
    if (q_flag_array[i] == ESTIMATE) {
      new_inc_q[i] = smtcpq_vec[kv];
      kv++;
    }
  }

}


/* Print the error message and the matrix if it is singular matrix.
*/
void Calculate::PrintErrMsg(char *str_mat,char *str_routin,int mat_dim,double **x_mat) {
  cout << "This is the singular " << str_mat << " matrix at routine "
       << str_routin << ":" << endl << endl;
  Lib::PrintOneMatrix(mat_dim, mat_dim, x_mat);
  PROBLEM "Please fix the error first and then Run the program.\n"
    RECOVER(NULL_ENTRY);
}

void Calculate::printIterationNumber(int iterationNumber) {
    if(iterationNumber != 1) {
      cout << ',';
    }
    if( ( (iterationNumber - 1) % 10 ) == 0 && iterationNumber != 1) {
      cout << endl;
    }
    cout << setw(4) << iterationNumber;
    cout.flush();
}

void Calculate::writeInvExpSecDer() {
  ofstream fixedInvExpSecDer("invExpSecDerFixed.log", ios::app);
  if(fixedInvExpSecDer.fail()) {
    PROBLEM "The file invExpSecDerFixed.log could not be opened for appending\nCalculate.cpp key 4924\n"
      RECOVER(NULL_ENTRY);
  }

  ofstream randomInvExpSecDer("invExpSecDerRandom.log", ios::app);
  if(randomInvExpSecDer.fail()) {
    PROBLEM "The file invExpSecDerRandom.log could not be opened for appending\nCalculate.cpp key 5291\n"
      RECOVER(NULL_ENTRY);
  }

  // est* are the respective indices into the estimated array.  Each index of estimated is a 1 if 
  // the variance component represented by that slot is estimated, a 0 otherwise.
  int estI = 0, estJ = 0;
  
  // data* are the respective indices into the inverse of the expected second derivative matrix
  // for the fixed effects.
  int dataI = 0, dataJ = 0;

  // estimated is an array used for holding (in a linear fashion) whether or not the variance 
  // component was estimated.  In the first iinitvcnum/3 of it, it holds whether or not the 
  // polygene random effects were estimated or not.  In the second iinitvcnum/3, it holds the same
  // information for the major gene.  And the last iinitvcnum/3, it hold the same info for the
  // the environmental effects.
  int *estimated;

  estimated = new int[3 * iinitvcnum];
  for(int i = 0; i < iinitvcnum; i++) {
    estimated[i] = (int)(s_flag_array[i] == ESTIMATE);
    estimated[iinitvcnum + i] = (int)(m1_flag_array[i] == ESTIMATE);
    estimated[2 * iinitvcnum + i] = (int)(t_flag_array[i] == ESTIMATE);
  }
  
  //fixedInvExpSecDer << ibdFileName << endl;
  randomInvExpSecDer << ibdFileName;

  if(ascert_flag == YES) {
    for(int i=0; i<imubeta_dim; i++) {
      for(int j=0; j<imubeta_dim; j++) {
	fixedInvExpSecDer << inv_Exsder_mat[j][i] << endl;//setw(15) << inv_Exsder_mat[j][i];
      }
      //fixedInvExpSecDer << endl;
    }

    dataI = -1;
    for (estI = 0; estI < 3 * iinitvcnum; estI++) {
      if (estimated[estI]) {
	dataI++;
      }
      dataJ = 0;
      randomInvExpSecDer << endl;
      for (estJ = 0; estJ < 3 * iinitvcnum; estJ++) {
	if (estimated[estI] && estimated[estJ]) {
	  // Since Splus will use the array command to shape this data into matrices,
	  // we need to print the transposed matrix because Splus will read it in 
	  // column-major order.
	  randomInvExpSecDer << setw(15)
			     << inv_Exsder_mat[imubeta_dim + dataJ++][imubeta_dim + dataI];
	} else {
	  randomInvExpSecDer << setw(15) << -9;
	}
      }
    }
  }else {
    for(int i=0; i<itraits*(1 + icovs); i++) {
      for(int j=0; j<itraits*(1 + icovs); j++) {
	fixedInvExpSecDer << inv_Exsder_mubeta_mat[j][i] << endl; //setw(15) << inv_Exsder_mubeta_mat[j][i];
      }
      //fixedInvExpSecDer << endl;
    }
    
    dataI = -1;
    for (estI = 0; estI < 3 * iinitvcnum; estI++) {
      if (estimated[estI]) {
	dataI++;
      }
      dataJ = 0;
      randomInvExpSecDer << endl;
      for (estJ = 0; estJ < 3 * iinitvcnum; estJ++) {
	if (estimated[estI] && estimated[estJ]) {
	  // Since Splus will use the array command to shape this data into matrices,
	  // we need to print the transposed matrix because Splus will read it in 
	  // column-major order.
	  randomInvExpSecDer << setw(15) << inv_Exsder_smtcpq_mat[dataJ++][dataI];
	} else {
	  randomInvExpSecDer << setw(15) << -9;
	}
      }
    }
  }

  //fixedInvExpSecDer << endl;
  randomInvExpSecDer << endl << endl;

  delete [] estimated;

  fixedInvExpSecDer.close();
  randomInvExpSecDer.close();
}

/*****************************************************************************
Title: allocateSavedVMatrix
Class: Calculate
Description: allocateSavedVMatrix is a method to allocate the memory used to
             store the V matrix for each family during an iteration of the 
             numerical algorithm in SRun.
Input: int familyCount - familyCount is the integer value representing the
                         amount of families
       int totalTraitCount - totalTraitCount is the interger value representing
                             the total amount of traits (for longitudinal
                             analysis this value is the number of traits
                             multiplied by the number of time points)
       int *familySizes - the array familySizes is use to store the number of
                          members of each family.  The family numbered i has
                          familySizes[i] members.
Output: A reference to a 3-dimensional array.  The first dimension represents
        the family index, the second represents the first (row, vertical)
        family member, and the third dimension represents the second (column,
        horizontal) family member.
Side Effects: If any of the input variables are NULL or have data out of
              expected range, the program prints an error message and
              terminates execution.
Author: Eric Lunde 01-12-04
*****************************************************************************/
double ***Calculate::allocateSavedVMatrix(int familyCount,
					  int totalTraitCount,
					  int *familySizes) {
  if(familyCount < 0) {
    PROBLEM "The variable 'int familyCount' is %d it should be non-negative.\nCalculate.cpp key 5120\n",
      familyCount RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d it should be non-negative.\nCalculate.cpp key 5126\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }
  if(familySizes == NULL) {
    PROBLEM "The variable 'int *familySizes' is NULL, it should have value.\nCalculate.cpp key 5132\n"
      RECOVER(NULL_ENTRY);
  }
  for(int i=0; i<familyCount; i++) {
    if(familySizes[i] < 1) {
      PROBLEM "The variable 'int familySize[%d]' is %d, it should be positive\n.Calculate.cpp key 5139\n",
      i, familySizes[i] RECOVER(NULL_ENTRY);
    }
  }

  // Allocate enough memory for each family to have a matrix of size
  // total trait values * size of the family by
  // total trait values * size of the family
  double ***localSavedVMatrix;
  
  localSavedVMatrix = (double ***) malloc(familyCount * sizeof(double **));
  
  for(int i=0; i<familyCount; i++) {
    localSavedVMatrix[i] = (double **) malloc(totalTraitCount
					 * familySizes[i]
					 * sizeof(double *));
    
    for(int j=0; j<totalTraitCount * familySizes[i]; j++) {
      localSavedVMatrix[i][j] = (double *) malloc(totalTraitCount
					     * familySizes[i]
					     * sizeof(double));
    }
  }

  return localSavedVMatrix;
}


/*****************************************************************************
Title: saveVMatrix
Class: Calculate
Description: saveVMatrix is a method to transfer the data stored in the
             variable V_mat (V matrix for a family at a given locus at some
             iteration) to the class variable savedVMatrix.
Input: double ***localSavedVMatrix - localSavedVMatrix is a local reference to
                                     the global variable savedVMatrix.  During
                                     this method, this data structure will be 
                                     assigned new values.
       double **V_mat - V_mat is a reference to the real data.  It is a matrix
                        of size total traits * family size by total traits
                        * family size.  It's data will be transferred to
                        localSavedVMatrix.
       int familyIndex - familyIndex represents the integer for a particular
                         family.  It let's us know which family's V matrix
                         we are about to save.
       int totalTraitCount - totalTraitCount represents the total amount of
                             traits that we have values for in the file fort.12
                             (in a longitudinal analysis, this value is equal
                             to the number of traits multiplied by the amount
                             of time points).
       int familySize - familySize represents the number of total family
                        members for this particular family.
Output: NONE.
Side Effects: If any of the array variables are NULL or the scalar variables
              are out of range, an error message will be printed and the
              program will terminate. 
Author: Eric Lunde, 01-12-04
*****************************************************************************/
void Calculate::saveVMatrix(double ***localSavedVMatrix, double **V_mat,
			    int familyIndex, int totalTraitCount,
			    int familySize) {
  if(localSavedVMatrix == NULL || localSavedVMatrix[0] == NULL
     || localSavedVMatrix[0][0] == NULL) {
    PROBLEM "The variable 'double ***localSavedVMatrix' is NULL, it should have value\n.Calculate.cpp key 5173\n"
      RECOVER(NULL_ENTRY);
  }
  if(V_mat == NULL || V_mat[0] == NULL) {
    PROBLEM "The variable 'double **V_mat' is NULL, it should have value.\nCalculate.cpp key 5179\n"
      RECOVER(NULL_ENTRY);
  }
  if(familyIndex < 0 || familyIndex >= nfam) {
    PROBLEM "The variable 'int familyIndex' is %d, it is out of the valid range: 0 - %d\nCalculate.cpp key 5185\n",
      familyIndex, nfam-1 RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5191\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }
  if(familySize < 0) {
    PROBLEM "The variable 'int familySize' is %d, it should be non-negative.\nCalculate.cpp key 5197\n",
      familySize RECOVER(NULL_ENTRY);
  }

  for(int i=0; i<totalTraitCount * familySize; i++) {
    for(int j=0; j<totalTraitCount * familySize; j++) {
      localSavedVMatrix[familyIndex][i][j] = V_mat[i][j];
    }
  }
}

/*****************************************************************************
Title: writeSavedVMatrixToFile
Class: Calculate
Description: writeSavedVMatrixToFile is a method to write the variable V_mat
             to its designated v.matrix.*.log file (where the * is the names of
             the ibd file that generated a particular loci.out in the file
             mloci.out).
Input: double ***localSavedVMatrix - localSavedVMatrix is a local reference to
                                     the global variable savedVMatrix.  This
                                     array will hold each family's V matrix for
                                     a given iteration in SRun.
       int familyCount - familyCount is the amount of families in the input
                         file fort.12
       int totalTraitCount - totalTraitCount is the amount of traits (in a
                             longitudinal analysis, this value is the number of
                             traits multiplied by the number of time points).
       int *familySizes - the array familySizes holds at each index, the size 
                          of the family (the family numbered i will have 
                          familySizes[i] family members).
       char *localIbdFileName - localIbdFileName is a local reference to the
                                class variable ibdFileName.  Since we want to 
                                observe the V matrix for each family at each
                                locus, we need to write the data to seperate
                                files.  This is accomplished by using this
                                variable to open a unique file for each locus
                                operated on through the course of the program's
                                execution.
Output: NONE.
Side Effects: If any of the array variables are NULL or the scalar variables
              are out of range, an error message will be printed and the
              program will terminate. 
Author: Eric Lunde, 01-12-04
*****************************************************************************/
void Calculate::writeSavedVMatrixToFile(double ***localSavedVMatrix,
					int familyCount, int totalTraitCount,
					int *familySizes,
					char *localIbdFileName) {
  if(localSavedVMatrix == NULL || localSavedVMatrix[0] == NULL
     || localSavedVMatrix[0][0] == NULL) {
    PROBLEM "The variable 'double ***localSavedVMatrix' is NULL, it should have value.\nCalculate.cpp key 5221\n"
      RECOVER(NULL_ENTRY);
  }
  if(familyCount < 0) {
    PROBLEM "The variable 'int familyCount' is %d, it should be non-negative.\nCalculate.cpp key 5227\n",
      familyCount RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5233\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }
  if(familySizes == NULL) {
    PROBLEM "The variable 'int *familySizes' is NULL, it should have value.\nCalculate.cpp key 5239\n"
      RECOVER(NULL_ENTRY);
  }
  for(int i=0; i<familyCount; i++) {
    if(familySizes[i] < 0) {
      PROBLEM "The variable 'int familySize[%d]' is %d, it should be non-negative.\nCalculate.cpp key 5246\n",
	i, familySizes[i] RECOVER(NULL_ENTRY);
    }
  }

  char compositeFileName[256];
  snprintf(compositeFileName, 256, "v.matrix.%s.log", localIbdFileName);

  ofstream v_matrix_log(compositeFileName);

  // Check to see if the file opened correctly.  Eric Lunde, 01-07-04
  if(v_matrix_log.fail()) {
    PROBLEM "The file %s could not be opened for writing.\nCalculate.cpp key 5260\n",
      compositeFileName RECOVER(NULL_ENTRY);
  }
  
  // Simply write the number of the family, the size of the family, and that 
  // family's V matrix.  Eric Lunde, 01-07-04
  for(int i=0; i<familyCount; i++) {
    v_matrix_log << "family: " << /*i*/uniqueFamilyIds[i] << " familySize: "
		 << familySizes[i]
		 << " traitCount: " << totalTraitCount << endl;
    for(int j=0; j<totalTraitCount * familySizes[i]; j++) {
      for(int k=0; k<totalTraitCount * familySizes[i]; k++) {
	// We round to zero the values that are essentially zero.
	if(abs(localSavedVMatrix[i][j][k]) < 1.0e-100) {
	  localSavedVMatrix[i][j][k] = 0;
	}
	// The only "value" that is not (less than the max or greater than
	// the min) is NaN (-NaN also).  This NaN most likely was generated
	// becasuse of a 0.0 / 0.0 (0 / 0 {integer} would cause a
	// core dump, 0.0 / 0.0 {double} causes an NaN)
	if(!(localSavedVMatrix[i][j][k] <= DBL_MAX
	     || localSavedVMatrix[i][j][k] >= -DBL_MAX)) {
	  cout << "localSavedVMatrix[" << i << "][" << j << "][" << k << "] = "
	       << localSavedVMatrix[i][j][k] << endl;
	}
	// The explicit printing of the ' ' is to ensure whitespace around
	// the values so a simple fscanf will delimit the values.
	v_matrix_log << ' ' << setw(14) << localSavedVMatrix[i][j][k];
      }
      v_matrix_log << endl;
    }
    v_matrix_log << endl;
  }
  
  // Close the file.
  v_matrix_log.close();
}

/*****************************************************************************
Title: freeSavedVMatrix
Class: Calculate
Description: freeSavedVMatrix is responsible for using the C++ memory
             management command 'free' to release the memory allocated for the
             storage of the V matrix for each family during each iteration in 
             SRun.
Input: double ***localSavedVMatrix - localSavedVMatrix is a local reference to
                                     the global variable savedVMatrix.  This
                                     array will hold each family's V matrix for
                                     a given iteration in SRun.
       int familyCount - familyCount is the amount of families in the input
                         file fort.12
       int totalTraitCount - totalTraitCount is the amount of traits (in a
                             longitudinal analysis, this value is the number of
                             traits multiplied by the number of time points).
       int *familySizes - the array familySizes holds at each index, the size 
                          of the family (the family numbered i will have 
                          familySizes[i] family members).
Output: NONE.
Side Effects: If any of the array variables are NULL or the scalar variables
              are out of range, an error message will be printed and the
              program will terminate. 
Author: Eric Lunde, 01-12-04
*****************************************************************************/
void Calculate::freeSavedVMatrix(double ***localSavedVMatrix,
				 int familyCount, int totalTraitCount,
				 int *familySizes) {
  if(localSavedVMatrix == NULL || localSavedVMatrix[0] == NULL
     || localSavedVMatrix[0][0] == NULL) {
    PROBLEM "The variable 'double ***localSavedVMatrix' is NULL, it should have value.\nCalculate.cpp key 5291\n"
      RECOVER(NULL_ENTRY);
  }
  if(familyCount < 0) {
    PROBLEM "The variable 'int familyCount' is %d, it should be non-negative.\nCalculate.cpp key 5297\n",
      familyCount RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5303\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }
  if(familySizes == NULL) {
    PROBLEM "The variable 'int *familySizes' is NULL, it should have value.\nCalculate.cpp key 5309\n"
      RECOVER(NULL_ENTRY);
  }
  for(int i=0; i<familyCount; i++) {
    if(familySizes[i] < 1) {
      PROBLEM "The variable 'int familySize[%d]' is %d, it should be positive.\nCalculate.cpp key 5316\n",
	i, familySizes[i] RECOVER(NULL_ENTRY);
    }
  }

  for(int i=0; i<familyCount; i++) {
    for(int j=0; j<totalTraitCount * familySizes[i]; j++) {
      free(localSavedVMatrix[i][j]);
    }
    free(localSavedVMatrix[i]);
  }  
  free(localSavedVMatrix);
}

/*****************************************************************************
Title: allocateSavedYBetaDiff
Class: Calculate
Description: allocateSavedYBetaDiff is a method to allocate the memory used to
             store the differences between a family's trait values and the beta
             estimates.
Input: int familyCount - familyCount is the integer value representing the
                         amount of families
       int totalTraitCount - totalTraitCount is the interger value representing
                             the total amount of traits (for longitudinal
                             analysis this value is the number of traits
                             multiplied by the number of time points)
       int *familySizes - the array familySizes is use to store the number of
                          members of each family.  The family numbered i has
                          familySizes[i] members.
Output: A reference to a 3-dimensional array.  The first dimension represents
        the family index, the second represents each trait studied, and the
        third dimension represents each trait-valued family member.
Side Effects: If any of the input variables are NULL or have data out of
              expected range, the program prints an error message and
              terminates execution.
Author: Eric Lunde 01-12-04
*****************************************************************************/
double ***Calculate::allocateSavedYBetaDiff(int familyCount,
					    int totalTraitCount,
					    int *familySizes) {
  if(familyCount < 0) {
    PROBLEM "The variable 'int familyCount' is %d, it should be non-negative.\nCalculate.cpp key 5336\n",
      familyCount RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5342\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }
  if(familySizes == NULL) {
    PROBLEM "The variable 'int *familySizes' is NULL, it should have value.\nCalculate.cpp key 5348\n"
      RECOVER(NULL_ENTRY);
  }
  for(int i=0; i<familyCount; i++) {
    if(familySizes[i] < 1) {
      PROBLEM "The variable 'int familySize[%d]' is %d, it should be positive.\nCalculate.cpp key 5355\n",
	i, familySizes[i] RECOVER(NULL_ENTRY);
    }
  }

  double ***localSavedYBetaDiff;

  localSavedYBetaDiff = (double ***) malloc(familyCount * sizeof(double **));
  for(int i=0; i<familyCount; i++) {
    localSavedYBetaDiff[i] = (double **) malloc(totalTraitCount
						* sizeof(double *));
    
    for(int j=0; j<totalTraitCount; j++) {
      localSavedYBetaDiff[i][j] = (double *) malloc(familySizes[i]
						    * sizeof(double));
    }
  }

  return localSavedYBetaDiff;
}

/*****************************************************************************
Title: saveYBetaDiff
Class: Calculate
Description: saveYBetaDiff is a method to transfer the data stored in the
             variable y_vec (family trait values - beta estimates) to the
             class variable savedYBetaDiff.
Input: double ***localSavedYBetaDiff - localSavedYBetaDiff is a local reference
                                       to the class variable savedYBetaDiff.
                                       At the time of development,
                                       localSavedYBetaDiff always refers to the
                                       same memory as savedYBetaDiff.  But in 
                                       the scheme of using less global
                                       variables, I have made it a parameter 
                                       rather then accessing it globally.  This
                                       memory will be given the values stored
                                       in y_vec.
       double **y_vec - y_vec is an array that holds the difference between
                        each trait valued family member's trait value and the
                        estimated beta.  The first direction is the number of
                        trait and the second direction is actual difference
       int familyIndex - familyIndex is the integer value representing the
                         number of the current family being saved.
       int totalTraitCount - totalTraitCount is the interger value representing
                             the total amount of traits (for longitudinal
                             analysis this value is the number of traits
                             multiplied by the number of time points)
       int familySize - familySize is use to store the number of members of a
                        family.
Output: NONE.
Side Effects: If any of the input variables are NULL or have data out of
              expected range, the program prints an error message and
              terminates execution.
Author: Eric Lunde 01-12-04
*****************************************************************************/
void Calculate::saveYBetaDiff(double ***localSavedYBetaDiff,
			      double **y_vec, int familyIndex,
			      int totalTraitCount, int familySize) {

  if(localSavedYBetaDiff == NULL || localSavedYBetaDiff[0] == NULL
     || localSavedYBetaDiff[0][0] == NULL) {
    PROBLEM "The variable 'double ***localSavedYBetaDiff' is NULL, it should have value.\nCalculate.cpp key 5385\n"
      RECOVER(NULL_ENTRY);
  }
  if(y_vec == NULL || y_vec[0] == NULL) {
    PROBLEM "The variable 'double **y_vec' is NULL, it should have value.\nCalculate.cpp key 5396\n"
      RECOVER(NULL_ENTRY);
  }
  if(familyIndex < 0) {
    PROBLEM "The variable 'int familyIndex' %d, it should be non-negative.\nCalculate.cpp key 5402\n",
      familyIndex RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5408\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }

  if(familySize < 0) {
    PROBLEM "The variable 'int familySize' is %d, it should non-negative.\nCalculate.cpp key 5414\n",
      familySize RECOVER(NULL_ENTRY);
  }

  // localSavedTraitValuedMembers[familyIndex] = traitValuedMembers;

  for(int i=0; i<totalTraitCount; i++) {
    for(int j=0; j<familySize; j++) {
      localSavedYBetaDiff[familyIndex][i][j] = y_vec[i][j];
    }
  }
}

/*****************************************************************************
Title: writeYBetaDiffToFile
Class: Calculate
Description: writeYBetaDiffToFile is a method to write the variable y_vec
             (familys trait values - beta estimates) to its designated
             y.beta.diff.*.log file (where the * is the names of the ibd file
             that generated a particular loci.out in the file mloci.out).
Input: double ***localSavedYBetaDiff - localSavedYBetaDiff is a local reference
                                       to the class variable savedYBetaDiff.
                                       At the time of development,
                                       localSavedYBetaDiff always refers to the
                                       same memory as savedYBetaDiff.  But in 
                                       the scheme of using less global
                                       variables, I have made it a parameter 
                                       rather then accessing it globally.  The
                                       data stored in this array will be
                                       written to the appropriate output file.
       int familyCount - familyCount is the integer value representing the
                         amount of families
       int totalTraitCount - totalTraitCount is the interger value representing
                             the total amount of traits (for longitudinal
                             analysis this value is the number of traits
                             multiplied by the number of time points)
       int *familySizes - the array familySizes is use to store the number of
                          size of each family.  The family numbered i has
                          familySizes[i] members.
       char *localIbdFileName - localIbdFileName is a local reference to the
                                class variable ibdFileName.  Since we want to 
                                observe the difference between the actual trait
                                values and the estimated betas for each locus,
                                we need to write the data to seperate files.
                                This is accomplished by using this variable to
                                open a unique file for each locus operated on
                                through the course of the program's execution.
Output: NONE.
Side Effects: If any of the input variables are NULL or have data out of
              expected range, the program prints an error message and
              terminates execution.
Author: Eric Lunde 01-12-04
*****************************************************************************/
void Calculate::writeSavedYBetaDiffToFile(double ***localSavedYBetaDiff,
					  int familyCount,
					  int totalTraitCount,
					  int *familySizes,
					  char *localIbdFileName) {
  if(localSavedYBetaDiff == NULL || localSavedYBetaDiff[0] == NULL
     || localSavedYBetaDiff[0][0] == NULL) {
    PROBLEM "The variable 'double ***localSavedYBetaDiff' is NULL, it should have value.\nCalculate.cpp key 5443\n"
      RECOVER(NULL_ENTRY);
  }
  if(familyCount < 0) {
    PROBLEM "The variable 'int familyCount' %d, it should be non-negative.\nCalculate.cpp key 5449\n",
      familyCount RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5455\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }
  if(familySizes == NULL) {
    PROBLEM "The variable 'int *familySizes' is NULL, it should have value.\nCalculate.cpp key 5461\n"
      RECOVER(NULL_ENTRY);
  }
  for(int i=0; i<familyCount; i++) {
    if(familySizes[i] < 0) {
      PROBLEM "The variable 'int familySizes[%d] is %d, it should be non-negative.\nCalculate.cpp key 5468\n",
	i, familySizes[i] RECOVER(NULL_ENTRY);
    }
  }

  char compositeFileName[256];
  snprintf(compositeFileName, 256, "y.beta.diff.%s.log", localIbdFileName);

  ofstream yBetaDiffLog(compositeFileName);

  // Check to see if the file opened correctly.  Eric Lunde, 01-07-04
  if(yBetaDiffLog.fail()) {
    PROBLEM "The file %s could not be opened for writting.\nCalculate.cpp key 5482\n",
      compositeFileName RECOVER(NULL_ENTRY);
  }
  
  // Simply write the number of the family and that family's difference between
  // their trait values and the betas.  Eric Lunde, 01-07-04
  for(int i=0; i<familyCount; i++) {
    yBetaDiffLog << "family: " << /*i*/uniqueFamilyIds[i]
		 << " familySize: " << familySizes[i]
		 << " totalTraitCount: " << totalTraitCount
		 << endl;
    
    for(int j=0; j<totalTraitCount; j++) {
      for(int k=0; k<familySizes[i]; k++) {
	// The explicit printing of the ' ' is to ensure whitespace around
	// the values so a simple fscanf will delimit the values.
	yBetaDiffLog << ' ' << setw(14) << localSavedYBetaDiff[i][j][k];
      }
      if(familySizes[i] != 0) {
	yBetaDiffLog << endl;
      }
    }
    yBetaDiffLog << endl;
  }

  // Close the file
  yBetaDiffLog.close();
}

/*****************************************************************************
Title: freeSavedYBetaDiff
Class: Calculate
Description: freeSavedYBetaDiff is responsible for using the C++ memory
             management command 'free' to release the memory allocated for the
             storage of the numerical differences between the y values (actual
             individual trait values) and the estimated beta values.
Input: double ***localSavedYBetaDiff - localSavedYBetaDiff is a local reference
                                       to the class variable savedYBetaDiff.
                                       At the time of development,
                                       localSavedYBetaDiff always refers to the
                                       same memory as savedYBetaDiff.  But in 
                                       the scheme of using less global
                                       variables, I have made it a parameter 
                                       rather then accessing it globally.  This
                                       is the memory that will be free'd during
                                       the execution of this method.
       int familyCount - familyCount is the integer value representing the
                         amount of families
       int totalTraitCount - totalTraitCount is the interger value representing
                             the total amount of traits (for longitudinal
                             analysis this value is the number of traits
                             multiplied by the number of time points)
Output: NONE.
Side Effects: If any of the input variables are NULL or have data out of
              expected range, the program prints an error message and
              terminates execution.
Author: Eric Lunde, 01-12-04
*****************************************************************************/
void Calculate::freeSavedYBetaDiff(double ***localSavedYBetaDiff,
				   int familyCount,
				   int totalTraitCount) {
  if(localSavedYBetaDiff == NULL || localSavedYBetaDiff[0] == NULL
     || localSavedYBetaDiff[0][0] == NULL) {
    PROBLEM "The variable 'double ***localSavedYBetaDiff' is NULL, it should have value.\nCalculate.cpp key 5516\n"
      RECOVER(NULL_ENTRY);
  }
  if(familyCount < 0) {
    PROBLEM "The variable 'int familyCount' %d, it should be non-negative.\nCalculate.cpp key 5522\n",
      familyCount RECOVER(NULL_ENTRY);
  }
  if(totalTraitCount < 1) {
    PROBLEM "The variable 'int totalTraitCount' is %d, it should be positive.\nCalculate.cpp key 5528\n",
      totalTraitCount RECOVER(NULL_ENTRY);
  }

  for(int i=0; i<familyCount; i++) {    
    for(int j=0; j<totalTraitCount; j++) {
      free(localSavedYBetaDiff[i][j]);
    }

    free(localSavedYBetaDiff[i]);
  }

  free(localSavedYBetaDiff);
}
