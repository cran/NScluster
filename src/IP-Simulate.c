#include <R.h>
#include <Rdefines.h>
#include "NScluster.h"

extern void F77_NAME(simipf)(int*, int*, int*, double*, double*, double*, double*, double*, int*, int*, double*, double*, double*, double*, int*, int*, int*);

SEXP simIP(SEXP ix, SEXP iy, SEXP iz, SEXP ty, SEXP amu, SEXP anu, SEXP p, SEXP c, SEXP pmax, SEXP omax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8;

    SEXP ans = R_NilValue, npts = R_NilValue, ncl = R_NilValue, x = R_NilValue, y = R_NilValue, xcl = R_NilValue, ycl = R_NilValue, ier = R_NilValue;
    double *xx, *yy, *xxcl, *yycl = NULL;
    int  *xnpts, *xncl, *nier = NULL;

    int npmax, nomax;
    int i;

    i1 = INTEGER_POINTER(ix);
    i2 = INTEGER_POINTER(iy);
    i3 = INTEGER_POINTER(iz);
    d1 = NUMERIC_POINTER(ty);
    d2 = NUMERIC_POINTER(amu);
    d3 = NUMERIC_POINTER(anu);
    d4 = NUMERIC_POINTER(p);
    d5 = NUMERIC_POINTER(c);
    i6 = INTEGER_POINTER(pmax);
    i7 = INTEGER_POINTER(omax);

    npmax = *i6;
    nomax = npmax * (*i7);

    PROTECT(ans = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ans, 0, npts = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, ncl = allocVector(INTSXP, npmax));
    SET_VECTOR_ELT(ans, 2, x = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 3, y = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 4, xcl = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 5, ycl = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 6, ier = allocVector(INTSXP, 1));

    i4 = INTEGER_POINTER(npts);
    i5 = INTEGER_POINTER(ncl);
    d6 = NUMERIC_POINTER(x);
    d7 = NUMERIC_POINTER(y);
    d8 = NUMERIC_POINTER(xcl);
    d9 = NUMERIC_POINTER(ycl);
    i8 = INTEGER_POINTER(ier);

    F77_CALL(simipf) (i1,i2,i3,d1,d2,d3,d4,d5,i4,i5,d6,d7,d8,d9,i6,i7,i8);

    xnpts = INTEGER(npts);
    xncl = INTEGER(ncl);
    xx = REAL(x);
    yy = REAL(y);
    xxcl = REAL(xcl);
    yycl = REAL(ycl);
    nier = INTEGER(ier);

    *xnpts = *i4;
    for(i=0; i<npmax; i++) xncl[i] = i5[i];
    for(i=0; i<npmax; i++) xx[i] = d6[i];
    for(i=0; i<npmax; i++) yy[i] = d7[i];
    for(i=0; i<nomax; i++) xxcl[i] = d8[i];
    for(i=0; i<nomax; i++) yycl[i] = d9[i];
    *nier = *i8;

    UNPROTECT(1);

    return ans;
}
