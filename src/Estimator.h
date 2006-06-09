#ifndef ESTIMATOR_H
#define ESTIMATOR_H

class Estimator{
 public:
  virtual double getmajorgene1()=0;
  virtual double getmajorgene2()=0;
  virtual double getpolygene()=0;
  virtual double getenvironment()=0;
  virtual double getC()=0;
  virtual double getP()=0;
  virtual double getQ()=0;
};

#endif
