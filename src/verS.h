#ifndef VER_S
#define VER_S

/* The following line is to avoid using the new Splus io classes
 * introduced in Splus 8.  It causes namespace conflicts and we
 * don't need the added GUI IO capabilities.
 * An alternative sugggested by Insightful is to make sure that
 * S.h is the last header file included by each source file.
 * I (PV) chose not to do that at this time because it would mean
 * editing all the C++ source files rather than just this one.
 */
#define S_NEWREDEF_H

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
