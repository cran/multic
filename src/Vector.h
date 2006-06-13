#ifndef VECTOR_H
#define VECTOR_H

#include "Model.h"

typedef struct linklist {
  Object *obj;
  struct linklist *next;
} LinkList;

class Vector{
  LinkList *head;
  LinkList *tail;
  int count;
 public:
  Vector();
  ~Vector();
  void add(Object *);
  int getcount();
  Object* elementAt(int);
  //contains() and containKey()
};

#endif
