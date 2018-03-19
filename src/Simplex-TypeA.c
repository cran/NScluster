#include <R.h>
#include <Rdefines.h>
#include "NScluster.h"

extern void F77_NAME(smplxaf)(double*, double*, int*, int*, double*, double*,
                              double*, double*, double*, double*, double*,
                              double*, int*, int*, int*, double*, double*,
                              double*, double*, double*, int*, int*, int*, int*);

SEXP smplxA(SEXP x, SEXP y, SEXP np, SEXP skip, SEXP ty, SEXP sclmu, SEXP sclnu,
            SEXP scla, SEXP scls1, SEXP scls2, SEXP x2, SEXP eps, SEXP itmax,
            SEXP itmax1, SEXP ipmax, SEXP ipflg)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14,*d15;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9;

    SEXP ans = R_NilValue, fn = R_NilValue, mple = R_NilValue, xx = R_NilValue,
  std = R_NilValue, f = R_NilValue, itr = R_NilValue, nip = R_NilValue,
  ipr = R_NilValue;
    double *xfn, *xmple, *xxx, *xstd, *xf = NULL;
    int  *xitr, *xnip, *xipr = NULL;

    int ipm, itm1;
    int i;

    d1 = NUMERIC_POINTER(x);
    d2 = NUMERIC_POINTER(y);
    i1 = INTEGER_POINTER(np);
    i2 = INTEGER_POINTER(skip);
    d3 = NUMERIC_POINTER(ty);
    d4 = NUMERIC_POINTER(sclmu);
    d5 = NUMERIC_POINTER(sclnu);
    d6 = NUMERIC_POINTER(scla);
    d7 = NUMERIC_POINTER(scls1);
    d8 = NUMERIC_POINTER(scls2);
    d9 = NUMERIC_POINTER(x2);
    d10 = NUMERIC_POINTER(eps);
    i3 = INTEGER_POINTER(itmax);
    i4 = INTEGER_POINTER(itmax1);
    i5 = INTEGER_POINTER(ipmax);
    i9 = INTEGER_POINTER(ipflg);

    ipm = *i5;
    itm1 = *i4;

    PROTECT(ans = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(ans, 0, fn = allocVector(REALSXP, ipm));
    SET_VECTOR_ELT(ans, 1, mple = allocVector(REALSXP, ipm*5));
    SET_VECTOR_ELT(ans, 2, xx = allocVector(REALSXP, itm1*5));
    SET_VECTOR_ELT(ans, 3, std = allocVector(REALSXP, itm1));
    SET_VECTOR_ELT(ans, 4, f = allocVector(REALSXP, itm1));
    SET_VECTOR_ELT(ans, 5, itr = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 6, nip = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 7, ipr = allocVector(INTSXP, ipm));

    d11 = NUMERIC_POINTER(fn);
    d12 = NUMERIC_POINTER(mple);
    d13 = NUMERIC_POINTER(xx);
    d14 = NUMERIC_POINTER(std);
    d15 = NUMERIC_POINTER(f);
    i6 = INTEGER_POINTER(itr);
    i7 = INTEGER_POINTER(nip);
    i8 = INTEGER_POINTER(ipr);

    F77_CALL(smplxaf) (d1,d2,i1,i2,d3,d4,d5,d6,d7,d8,d9,d10,i3,i4,i5,d11,d12,d13,d14,d15,i6,i7,i8,i9);

    xfn = REAL(fn);
    xmple = REAL(mple);
    xxx = REAL(xx);
    xstd = REAL(std);
    xf = REAL(f);
    xitr = INTEGER(itr);
    xnip = INTEGER(nip);
    xipr = INTEGER(ipr);

    for(i=0; i<ipm; i++) xfn[i] = d11[i];
    for(i=0; i<ipm*5; i++) xmple[i] = d12[i];
    for(i=0; i<itm1*5; i++) xxx[i] = d13[i];
    for(i=0; i<itm1; i++) xstd[i] = d14[i];
    for(i=0; i<itm1; i++) xf[i] = d15[i];
    *xitr = *i6;
    *xnip = *i7;
    for(i=0; i<ipm; i++) xipr[i] = i8[i];

    UNPROTECT(1);

    return ans;
}
