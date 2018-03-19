#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP palmA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                  SEXP, SEXP, SEXP);
extern SEXP palmB(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                  SEXP, SEXP);
extern SEXP palmC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                  SEXP, SEXP);
extern SEXP palmIP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP);
extern SEXP palmT(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simB(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simIP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP simThom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP smplxA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP smplxB(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP, SEXP);
extern SEXP smplxC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                   SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP smplxIP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP smplxThom(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                      SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"palmA",     (DL_FUNC) &palmA,     13},
    {"palmB",     (DL_FUNC) &palmB,     12},
    {"palmC",     (DL_FUNC) &palmC,     12},
    {"palmIP",    (DL_FUNC) &palmIP,    12},
    {"palmT",     (DL_FUNC) &palmT,     10},
    {"simA",      (DL_FUNC) &simA,       9},
    {"simB",      (DL_FUNC) &simB,       9},
    {"simC",      (DL_FUNC) &simC,      10},
    {"simIP",     (DL_FUNC) &simIP,      8},
    {"simThom",   (DL_FUNC) &simThom,    7},
    {"smplxA",    (DL_FUNC) &smplxA,    16},
    {"smplxB",    (DL_FUNC) &smplxB,    14},
    {"smplxC",    (DL_FUNC) &smplxC,    15},
    {"smplxIP",   (DL_FUNC) &smplxIP,   15},
    {"smplxThom", (DL_FUNC) &smplxThom, 12},
    {NULL, NULL, 0}
};

void R_init_NScluster(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
