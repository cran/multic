#include <iostream>

#include "verS.h"
#ifdef USING_R
#include <Rinternals.h>
#include <Rdefines.h>
#endif

extern "C" {
  /**
   * sortIdPairs - Accepts 2 R vectors as input.  The vectors are
   * subject id pairs.  This function makes sure that id1 <= id2.
   * This sorting is necessary to help merlin2mloci.  The merlin
   * files often have the id pairs where id1 > id2.
   */
  void sortIdPairs(s_object *id1, s_object *id2) {
    //cout << "inside sortIdPairs" << endl;

    // ADD WARNING IF id1 and id2 ARE DIFFERENT LENGTHS

    Sint *pid1 = INTEGER_POINTER(id1);
    Sint *pid2 = INTEGER_POINTER(id2);

    Sint temp;
    for (Sint i=0; i < LENGTH(id1); i++) {
      //cout << "i = " << i << endl;
      //cout << "id[" << i << "] = " << pid1[i] << endl;

      if (pid1[i] > pid2[i]) {
	temp = pid1[i];
	pid1[i] = pid2[i];
	pid2[i] = temp;
      }
    }

  }
} // end extern "C"
