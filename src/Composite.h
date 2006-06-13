#ifndef COMPOSITE_H
#define COMPOSITE_H

#include "Model.h"
#include "Vector.h"

class Composite : public Model {
 private:
  Vector vector;
 public:
  Composite();
  virtual ~Composite();

  virtual void getv(double** v, int n);
  char* getname();
  double getbeta_hat();
  void add(Model *);
  void getvsum(double** vorig, double** vstack,int n,ShrinkData *sd);
  int getbeta_hat(double* beta, int n);
  int getElementCount();
  Model* getElementAt(int);

};

#endif
