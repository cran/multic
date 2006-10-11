#ifndef VER_S
#define VER_S

#include <S.h>

/*
** The next part is to make it easier to have common code with Splus and R, where
** "integers" are int.  In S they are Sint.  In R they are int.
** The definition of Sint also addresses issues with 64-bit architectures.
*/
#ifdef USING_R
#include <R.h>
typedef int Sint;
#else
typedef long Sint;
#endif
#endif
