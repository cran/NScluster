#include <R.h>
#include <Rdefines.h>
#include "NScluster.h"

extern void F77_NAME(palmtf)(double*, double*, int*, double*, double*, double*,
                             double*, double*, int*, int*, double*, double*);

SEXP palmT(SEXP x, SEXP y, SEXP np, SEXP delta, SEXP ty, SEXP amu, SEXP anu,
           SEXP v, SEXP m, SEXP jmax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3;

    SEXP ans = R_NilValue, palm = R_NilValue, palm1 = R_NilValue;
    double *xpalm, *xpalm1 = NULL;

    int jm, mm;
    int i;

    d1 = NUMERIC_POINTER(x);
    d2 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(np);
    d3 = NUMERIC_POINTER(delta);
    d4 = NUMERIC_POINTER(ty);
    d5 = NUMERIC_POINTER(amu);
    d6 = NUMERIC_POINTER(anu);
    d7 = NUMERIC_POINTER(v);
    i2 = INTEGER_POINTER(m);
    i3 = INTEGER_POINTER(jmax);

    jm = *i3;
    mm = (*i2) * jm;

    PROTECT(ans = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ans, 0, palm = allocVector(REALSXP, jm));
    SET_VECTOR_ELT(ans, 1, palm1 = allocVector(REALSXP, mm));

    d8 = NUMERIC_POINTER(palm);
    d9 = NUMERIC_POINTER(palm1);

    F77_CALL(palmtf) (d1,d2,i1,d3,d4,d5,d6,d7,i2,i3,d8,d9);

    xpalm = REAL(palm);
    xpalm1 = REAL(palm1);

    for(i=0; i<jm; i++) xpalm[i] = d8[i];
    for(i=0; i<mm; i++) xpalm1[i] = d9[i];

    UNPROTECT(1);

    return ans;
}
