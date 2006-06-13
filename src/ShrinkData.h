#ifndef SHRINK_DATA_H
#define SHRINK_DATA_H

#include "multic.h"

class ShrinkData {
  Family *famdata;
 public:
  ShrinkData(Family *famd);
  int  getfamrealcount(int);
  void Changedata(double**,int);
  void Changedata(double*,int);
  void Shrinkdata(double**,int,int);
  void Shrinkdata(double*,int,int);
};

#endif
