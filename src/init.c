#include "regF77.h"
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Fortran calls */

static const R_FortranMethodDef FortEntries[] = {
    {"xqgausip",  (DL_FUNC) &F77_NAME(xqgausip),  14},
    {"palmt",     (DL_FUNC) &F77_NAME(palmt),     12},
    {"xqgausa",   (DL_FUNC) &F77_NAME(xqgausa),   15},
    {"palmb",     (DL_FUNC) &F77_NAME(palmb),     14},
    {"palmc",     (DL_FUNC) &F77_NAME(palmc),     14},
    {"smplxip",   (DL_FUNC) &F77_NAME(smplxip),   23},
    {"smplxthom", (DL_FUNC) &F77_NAME(smplxthom), 20},
    {"smplxa",    (DL_FUNC) &F77_NAME(smplxa),    24},
    {"smplxb",    (DL_FUNC) &F77_NAME(smplxb),    22},
    {"smplxc",    (DL_FUNC) &F77_NAME(smplxc),    23},
    {"simip",     (DL_FUNC) &F77_NAME(simip),     15},
    {"simthom",   (DL_FUNC) &F77_NAME(simthom),   14},
    {"sima",      (DL_FUNC) &F77_NAME(sima),      16},
    {"simb",      (DL_FUNC) &F77_NAME(simb),      22},
    {"simc",      (DL_FUNC) &F77_NAME(simc),      23},
    {NULL, NULL, 0}
};

void attribute_visible R_init_NScluster(DllInfo *dll) 
{
    R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
