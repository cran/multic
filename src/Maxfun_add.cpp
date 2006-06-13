#include "Maxfun.h"
#include "Maxfun_add.h"
#include "Maxfun_source.h"
#include "TraitMarkerCov_par.h"
#include "InitValue_par.h"
#include "multic.h"

void Maxfun_add::run(TraitMarkerCov_par *tmc, InitValue_par *init,
		     FortData *fa, int ds, ShareRelation *sa, int rs) {
  Maxfun* m = new Maxfun(tmc, init, fa, ds, sa, rs);
  Maxfun_source m1(m);
  m1.maxfun(theta,&f,&nfe,&lfl);
  delete m;
}
