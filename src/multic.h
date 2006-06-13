#ifndef MULTIC_H
#define MULTIC_H

/*********************************************************************
 multic.h  ----- define global variables
 ********************************************************************/

#define REVISION        "5.09" 
#define DISTRIBDATE     "08-22-03"
#define INITIALVALUE    10 
#define MULTIVARIATE    1
#define LONGITUDINAL    2 
#define ASCINFO         8 

#ifdef _DEBUG
#define DEBUG           1
#endif
#define NONDEBUG        0

/* ------ define for the condition of iteration.    */
#define BREAKVAL        0.00001 
#define EPSILON1        0.0001 
#define EPSILON2        0.001 
#define SATISFIEDTIMES  1 
#define BOUNDARY        0.0 
#define CLOSETOBOUNDARY 0.00001 
#define BREAKSTEPHALF   200 

/* ------ define for the program  ------- */
#define MAXNUMFAM      10000 
#define MAXCH          80 
#define FAMMEMNUM      200
#define BUF            256 
#define STRING         32 
#define MAXNUMTRAIT    5 
#define MAXNUMMARKER   2 
#define MAXNUMCOV      10 
#define MAXREPEATMST   MAXNUMTRAIT  
#define TOTALCOV       (MAXNUMCOV*MAXREPEATMST)
#define MAXNUMPROB     5 
#define INITCOVNUM     (MAXNUMCOV*MAXNUMTRAIT)
#define INITVCNUM      ((MAXNUMTRAIT*(MAXNUMTRAIT+1))/2)
#define NUMOFSUM       ((MAXNUMTRAIT+TOTALCOV+MAXNUMMARKER)*2)
#define NUMALLELE      30 
#define NUMPHENOTYPE   ((NUMALLELE+1)*NUMALLELE/2)
#define NUMLIKEHD      ((FAMMEMNUM-2)*(FAMMEMNUM-3)/2)
#define NUMTRAITPAIR   NUMLIKEHD
#define PHENOTYPE      0
#define GENOTYPE       1
#define SETON          0
#define SETOFF         1
#define MUBETA         1 
#define SMTCPQ         2 
#define YES            1 
#define NO             0   
#define GOOD           1.0 
#define BAD            0.0   
#define RIGHT          1 
#define BEGIN          0   
#define FIRST          1   
#define SECOND         2   
#define END            3   
#define non_value     -1 
#define MULTIHYPHEN    "--------------------------------------------------"
#define MULTIPOUND     "##################################################"
#define INTERNALLY     0
#define EXTERNALLY     1

enum {TRAIT, MARKER,   COV };
enum {FIXED, ESTIMATE, CONSTRAINT };

typedef struct {
    char   a_name[STRING]; 
    double a_fq;
} Allele;

typedef struct {
    char pheno_tp[STRING]; 
    char pheno_set[MAXCH];
    char pheno_marker[STRING];
} Phenotype;

typedef struct {
    char   t_name[STRING]; 
    double t_misval;
    double Th_mu;
} Trait;

typedef struct {
    char   c_name[STRING]; 
    double c_misval;
} Cov;

typedef struct {
    char      m_name[STRING]; 
    char      m_misval[STRING];
    int       n_allele;
    int       n_pheno;
    Allele    allele[NUMALLELE];
    Phenotype pheno[NUMPHENOTYPE];
} Marker;

typedef struct {
    char mk[STRING]; 
    char alle1[STRING];
    char alle2[STRING];
} Pmarker;

typedef struct {
    double  trait[MAXNUMTRAIT];
    double  cov[MAXNUMCOV];
    Pmarker marker[MAXNUMMARKER];
    Pmarker smarker;
    int     missing_traitflag; 
} Person;

typedef struct {
    double trait_1[NUMTRAITPAIR];
    double trait_2[NUMTRAITPAIR];
} SibTrait;

typedef struct {
     int    prob_mem;
     double trait_p;
} Proband;

typedef struct {
    int      n_sib;
    int      sib_pair;
    double   likehd;
    double   likehd_1[NUMLIKEHD];
    double   likehd_2[NUMLIKEHD];
    Proband  proband[MAXNUMPROB];
    SibTrait sibtrait[MAXNUMTRAIT];
    Person   person[FAMMEMNUM];
} Family;

/********
ShareRelation was created to hold one line of text from share.out.
seqIdX stands for Xth person's sequential identifier. - Eric Lunde 7/29/03
********/
#define MAX_ID_LENGTH 31
typedef struct {
  char seqId1[MAX_ID_LENGTH + 1];
  char seqId2[MAX_ID_LENGTH + 1];
  double geneticSimilarity;
  double areSiblings;
  double areSpouses;
  double areParentChild;
} ShareRelation;

void printShareRelation(ShareRelation *sr);

/********
FortData was created to hold one line of text from fort.12.
seqId stands for sequential identification. - Eric Lunde 7/29/03
********/
typedef struct {
  char *familyId;
  char *seqId;
  char *fatherId;
  char *motherId;
  int sex;
  double *traits;
  double *covariates;
  int proband;
} FortData;

void printFortData(FortData *fd, int, int, int);

#endif
