/* Added by J. Fox 2019-09-06 */

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */

extern void F77_NAME(bfac)(void *, void *, void *, void *, void *);
extern void F77_NAME(chisq)(void *);
extern void F77_NAME(chol2)(void *, void *, void *, void *, void *);
extern void F77_NAME(chols)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ctrsc)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(emn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gamm)(void *);
extern void F77_NAME(gauss)(void);
extern void F77_NAME(gtmc)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(gtoc)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(is1n)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(is2n)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(initn)(void *, void *);
extern void F77_NAME(invtrn)(void *, void *, void *, void *);
extern void F77_NAME(lasts)(void *, void *, void *, void *);
extern void F77_NAME(layers)(void *, void *, void *, void *);
extern void F77_NAME(lobsn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lprin)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mkpsi)(void *, void *);
extern void F77_NAME(mmn)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(moden)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ninvwn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(nmons)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ph2thn)(void *, void *, void *, void *);
extern void F77_NAME(ps1n)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ps2n)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rangen)(void *);
extern void F77_NAME(rngs)(void *);
extern void F77_NAME(sigex)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(sjn)(void *, void *, void *, void *);
extern void F77_NAME(stvaln)(void *, void *, void *, void *);
extern void F77_NAME(swp)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(swpobs)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tobsmn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tobsn)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"bfac",   (DL_FUNC) &F77_NAME(bfac),    5},
    {"chisq",  (DL_FUNC) &F77_NAME(chisq),   1},
    {"chol2",  (DL_FUNC) &F77_NAME(chol2),   5},
    {"chols",  (DL_FUNC) &F77_NAME(chols),   6},
    {"ctrsc",  (DL_FUNC) &F77_NAME(ctrsc),   6},
    {"emn",    (DL_FUNC) &F77_NAME(emn),    20},
    {"gamm",   (DL_FUNC) &F77_NAME(gamm),    1},
    {"gauss",  (DL_FUNC) &F77_NAME(gauss),   0},
    {"gtmc",   (DL_FUNC) &F77_NAME(gtmc),    7},
    {"gtoc",   (DL_FUNC) &F77_NAME(gtoc),    7},
    {"is1n",   (DL_FUNC) &F77_NAME(is1n),   16},
    {"is2n",   (DL_FUNC) &F77_NAME(is2n),   16},
    {"initn",  (DL_FUNC) &F77_NAME(initn),   2},
    {"invtrn", (DL_FUNC) &F77_NAME(invtrn),  4},
    {"lasts",  (DL_FUNC) &F77_NAME(lasts),   4},
    {"layers", (DL_FUNC) &F77_NAME(layers),  4},
    {"lobsn",  (DL_FUNC) &F77_NAME(lobsn),  14},
    {"lprin",  (DL_FUNC) &F77_NAME(lprin),  10},
    {"mkpsi",  (DL_FUNC) &F77_NAME(mkpsi),   2},
    {"mmn",    (DL_FUNC) &F77_NAME(mmn),     9},
    {"moden",  (DL_FUNC) &F77_NAME(moden),   9},
    {"ninvwn", (DL_FUNC) &F77_NAME(ninvwn), 10},
    {"nmons",  (DL_FUNC) &F77_NAME(nmons),   6},
    {"ph2thn", (DL_FUNC) &F77_NAME(ph2thn),  4},
    {"ps1n",   (DL_FUNC) &F77_NAME(ps1n),   12},
    {"ps2n",   (DL_FUNC) &F77_NAME(ps2n),   19},
    {"rangen", (DL_FUNC) &F77_NAME(rangen),  1},
    {"rngs",   (DL_FUNC) &F77_NAME(rngs),    1},
    {"sigex",  (DL_FUNC) &F77_NAME(sigex),   7},
    {"sjn",    (DL_FUNC) &F77_NAME(sjn),     4},
    {"stvaln", (DL_FUNC) &F77_NAME(stvaln),  4},
    {"swp",    (DL_FUNC) &F77_NAME(swp),     7},
    {"swpobs", (DL_FUNC) &F77_NAME(swpobs),  7},
    {"tobsmn", (DL_FUNC) &F77_NAME(tobsmn), 15},
    {"tobsn",  (DL_FUNC) &F77_NAME(tobsn),  11},
    {NULL, NULL, 0}
};

void R_init_norm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
