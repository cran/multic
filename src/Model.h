#ifndef MODEL_H
#define MODEL_H

#include "Least.h"
#include "multic.h"
#include <fstream>

class Object{
};

/*
struct linklist {
  Object *obj;
  struct linklist *next;
};

typedef struct linklist *LinkList;

class Vector{
  LinkList head;
  LinkList tail;
  int count;
 public:
  Vector();
  ~Vector();
  void add(Object *);
  int getcount();
  Object* elementAt(int);
  //contains() and containKey()
};
*/

class Model : public Object{
 protected:
  Estimator* estimator;
  char *name;
 public:
  virtual ~Model() {};
  virtual double getbeta_hat() = 0;
  virtual void getv(double** v, int n) = 0;
  virtual char* getname()=0;
};

// for one marker case
class MajorGene1 : public Model{
  ifstream *fp_loci;
  double sigma;
 public:
  MajorGene1(Estimator*);
  virtual ~MajorGene1();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};

// for two markers case
class MajorGene2 : public Model{
  ifstream* fp_loci;
  double sigma;
 public:
  MajorGene2(Estimator*);
  virtual ~MajorGene2();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};      


class PolyGene : public Model{
  ShareRelation *shareArray;
  int relationIndex;
  double sigma;
 public:
  PolyGene(Estimator*, ShareRelation *);
  virtual ~PolyGene();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};

class SS : public Model{
  ShareRelation *shareArray;
  int relationIndex;
  double sigma;
 public:
  SS(Estimator*, ShareRelation *);
  virtual ~SS();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};

class PS : public Model{
  ShareRelation *shareArray;
  int relationIndex;
  double sigma;
 public:
  PS(Estimator*, ShareRelation *);
  virtual ~PS();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};

class PP : public Model{
  ShareRelation *shareArray;
  int relationIndex;
  double sigma;
 public:
  PP(Estimator*, ShareRelation *);
  virtual ~PP();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};

class Environment : public Model{
  double sigma;
 public:
  Environment(Estimator*);
  virtual ~Environment();
  virtual double getbeta_hat();
  virtual void getv(double** v, int n);
  virtual char* getname();
};

/*
class Composite : public Model {
  Vector vector;
 public:
  Composite();
  virtual ~Composite();
  void add(Model *);
  int getElementCount();
  void getvsum(double** vorig, double** vstack,int n,ShrinkData *sd);
  int getbeta_hat(double* beta, int n);
  Model* getElementAt(int);
  virtual double getbeta_hat();
  virtual void getv(double** v, int n){cout << "call wrong method!\n";};
  virtual char* getname();
};
*/

#endif // end ifndef MODEL_H
