#include "Vector.h"
#include <cstdlib>

Vector::Vector(){
  head = NULL;
  tail = NULL;
  count = 0;
}

Vector::~Vector(){
  LinkList *temp;
  while(head!=NULL){
    temp = head;
    head = head->next;
    free(temp);
  }
  count = 0;
}

int Vector::getcount(){
  return count;
}

void Vector::add(Object *object){
  LinkList *temp = (LinkList *) malloc( sizeof(LinkList) );
  temp->obj = object;
  temp->next = NULL;
  if(head==NULL){
    head = temp;
    tail = head;
  } else{
    tail->next = temp;
    tail = tail->next;
  }
  count++;
}

Object* Vector::elementAt(int index){
  LinkList *temp = head;
  for(int i=index;i>0;i--)
    temp = temp->next;
  return temp->obj;
}
