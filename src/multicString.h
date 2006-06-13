#ifndef MULTIC_STRING_H
#define MULTIC_STRING_H

#include <S.h>
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif

void multic_SET_STRING_ELT(s_object *stringArray, int index,
			   const char *element);

char *multic_STRING_ELT(s_object *stringArray, int index);

// Function to safely do a strcpy
char *multic_strncpy(char *dest, const char *src, int n);

#endif
