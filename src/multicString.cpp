#include <cstring>

#include <S.h>
#ifdef USING_R
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#endif

using namespace std;

// multic_SET_STRING_ELT abstracts the concept of putting a string into 
// an S character vector.  This is currently only used for added string to 
// unoccupied vector elements.  So if this is used to overwrite some string
// value in a vector, I do not worry about free'ing that string's memory
// space.  Eric Lunde 2005-09-02
void multic_SET_STRING_ELT(s_object *stringArray, int index,
			   const char *element) {
#ifdef USING_R
  SET_STRING_ELT(stringArray, index, mkChar(element));
#else
  CHARACTER_POINTER(stringArray)[index] = Salloc(strlen(element) + 1, char);
  strcpy(CHARACTER_POINTER(stringArray)[index], element);
#endif
}

// multic_SET_STRING abstracts the concept of getting a string out of
// an S character vector.  Eric Lunde 2005-09-02
char *multic_STRING_ELT(s_object *stringArray, int index) {
  char *element = 
#ifdef USING_R
    CHAR(STRING_ELT(stringArray, index))
#else
    CHARACTER_POINTER(stringArray)[index]
#endif
    ;
  return element;
}

/**
 * I found myself having to do both of these commands whenever I wanted a
 * safe string copy.  So I put them together here.
 */
char *multic_strncpy(char *dest, const char *src, int n) {
  char *result = strncpy(dest, src, n);
  result[n] = '\0';

  return result;
}
