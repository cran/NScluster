#include <R.h>
#include <Rdefines.h>
#include "NScluster.h"

extern void F77_NAME(simthomf)(int*, double*, double*, double*, double*, int*,
                               int*, double*, double*, double*, double*, int*,
                               int*, int*);

SEXP simThom(SEXP ix, SEXP ty, SEXP amu, SEXP anu, SEXP sig, SEXP pmax,
             SEXP omax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8;
    int *i1,*i2,*i3,*i4,*i5,*i6;

    SEXP ans = R_NilValue, npts = R_NilValue, ncl = R_NilValue, x = R_NilValue,
  y = R_NilValue, xcl = R_NilValue, ycl = R_NilValue, ier = R_NilValue;
    double *xx, *yy, *xxcl, *yycl = NULL;
    int  *xnpts, *xncl, *nier = NULL;

    int npmax, nomax;
    int i;

    i1 = INTEGER_POINTER(ix);
    d1 = NUMERIC_POINTER(ty);
    d2 = NUMERIC_POINTER(amu);
    d3 = NUMERIC_POINTER(anu);
    d4 = NUMERIC_POINTER(sig);
    i4 = INTEGER_POINTER(pmax);
    i5 = INTEGER_POINTER(omax);

    npmax = *i4;
    nomax = npmax * (*i5);

    PROTECT(ans = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ans, 0, npts = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, ncl = allocVector(INTSXP, npmax));
    SET_VECTOR_ELT(ans, 2, x = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 3, y = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 4, xcl = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 5, ycl = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 6, ier = allocVector(INTSXP, 1));

    i2 = INTEGER_POINTER(npts);
    i3 = INTEGER_POINTER(ncl);
    d5 = NUMERIC_POINTER(x);
    d6 = NUMERIC_POINTER(y);
    d7 = NUMERIC_POINTER(xcl);
    d8 = NUMERIC_POINTER(ycl);
    i6 = INTEGER_POINTER(ier);

    F77_CALL(simthomf) (i1,d1,d2,d3,d4,i2,i3,d5,d6,d7,d8,i4,i5,i6);

    xnpts = INTEGER(npts);
    xncl = INTEGER(ncl);
    xx = REAL(x);
    yy = REAL(y);
    xxcl = REAL(xcl);
    yycl = REAL(ycl);
    nier = INTEGER(ier);

    *xnpts = *i2;
    for(i=0; i<npmax; i++) xncl[i] = i3[i];
    for(i=0; i<npmax; i++) xx[i] = d5[i];
    for(i=0; i<npmax; i++) yy[i] = d6[i];
    for(i=0; i<nomax; i++) xxcl[i] = d7[i];
    for(i=0; i<nomax; i++) yycl[i] = d8[i];
    *nier = *i6;

    UNPROTECT(1);

    return ans;
}
