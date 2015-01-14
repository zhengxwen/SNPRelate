//
// R_SNPRelate.c: link to R_GDS2.h in the gdsfmt package
//

#include <R_GDS.h>

// do not modify this file, R_GDS.c/R_GDS2.h is from the gdsfmt package

#if (defined(GDSFMT_R_VERSION) && (GDSFMT_R_VERSION>=0x010103))
#   include <R_GDS2.h>
#else
#   include <R_GDS.c>
#endif
