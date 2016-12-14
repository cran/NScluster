#include <R.h>
#include <Rdefines.h>
#include "NScluster.h"

extern void F77_NAME(simbf)(int*, int*, int*, double*, double*, double*, double*, double*, double*, int*,  int*, double*, double*, double*, double*, int*,  int*, double*, double*, double*, double*, int*, int*, int*);

SEXP simB(SEXP ix, SEXP iy, SEXP iz, SEXP ty, SEXP amu, SEXP anu, SEXP a, SEXP sig1, SEXP sig2, SEXP pmax, SEXP omax)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12,*d13,*d14;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10;

    SEXP ans = R_NilValue, m1 = R_NilValue, ncl1 = R_NilValue, x1 = R_NilValue, y1 = R_NilValue, xx1 = R_NilValue, yy1 = R_NilValue, m2 = R_NilValue, ncl2 = R_NilValue, x2 = R_NilValue, y2 = R_NilValue, xx2 = R_NilValue, yy2 = R_NilValue, ier = R_NilValue;
    double *x10, *y10, *xx10, *yy10, *x20, *y20, *xx20, *yy20 = NULL;
    int  *xm1, *xncl1, *xm2, *xncl2, *nier = NULL;

    int npmax, nomax;
    int i;

    i1 = INTEGER_POINTER(ix);
    i2 = INTEGER_POINTER(iy);
    i3 = INTEGER_POINTER(iz);
    d1 = NUMERIC_POINTER(ty);
    d2 = NUMERIC_POINTER(amu);
    d3 = NUMERIC_POINTER(anu);
    d4 = NUMERIC_POINTER(a);
    d5 = NUMERIC_POINTER(sig1);
    d6 = NUMERIC_POINTER(sig2);
    i8 = INTEGER_POINTER(pmax);
    i9 = INTEGER_POINTER(omax);

    npmax = *i8;
    nomax = npmax * (*i9);

    PROTECT(ans = allocVector(VECSXP, 13));
    SET_VECTOR_ELT(ans, 0, m1 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 1, ncl1 = allocVector(INTSXP, npmax));
    SET_VECTOR_ELT(ans, 2, x1 = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 3, y1 = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 4, xx1 = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 5, yy1 = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 6, m2 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 7, ncl2 = allocVector(INTSXP, npmax));
    SET_VECTOR_ELT(ans, 8, x2 = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 9, y2 = allocVector(REALSXP, npmax));
    SET_VECTOR_ELT(ans, 10, xx2 = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 11, yy2 = allocVector(REALSXP, nomax));
    SET_VECTOR_ELT(ans, 12, ier = allocVector(INTSXP, 1));

    i4 = INTEGER_POINTER(m1);
    i5 = INTEGER_POINTER(ncl1);
    d7 = NUMERIC_POINTER(x1);
    d8 = NUMERIC_POINTER(y1);
    d9 = NUMERIC_POINTER(xx1);
    d10 = NUMERIC_POINTER(yy1);
    i6 = INTEGER_POINTER(m2);
    i7 = INTEGER_POINTER(ncl2);
    d11 = NUMERIC_POINTER(x2);
    d12 = NUMERIC_POINTER(y2);
    d13 = NUMERIC_POINTER(xx2);
    d14 = NUMERIC_POINTER(yy2);
    i10 = INTEGER_POINTER(ier);

    F77_CALL(simbf) (i1,i2,i3,d1,d2,d3,d4,d5,d6,i4,i5,d7,d8,d9,d10,i6,i7,d11,d12,d13,d14,i8,i9,i10);

    xm1 = INTEGER(m1);
    xncl1 = INTEGER(ncl1);
    x10 = REAL(x1);
    y10 = REAL(y1);
    xx10 = REAL(xx1);
    yy10 = REAL(yy1);
    xm2 = INTEGER(m2);
    xncl2 = INTEGER(ncl2);
    x20 = REAL(x2);
    y20 = REAL(y2);
    xx20 = REAL(xx2);
    yy20 = REAL(yy2);
    nier = INTEGER(ier);

    *xm1 = *i4;
    for(i=0; i<npmax; i++) xncl1[i] = i5[i];
    for(i=0; i<npmax; i++) x10[i] = d7[i];
    for(i=0; i<npmax; i++) y10[i] = d8[i];
    for(i=0; i<nomax; i++) xx10[i] = d9[i];
    for(i=0; i<nomax; i++) yy10[i] = d10[i];
    *xm2 = *i6;
    for(i=0; i<npmax; i++) xncl2[i] = i7[i];
    for(i=0; i<npmax; i++) x20[i] = d11[i];
    for(i=0; i<npmax; i++) y20[i] = d12[i];
    for(i=0; i<nomax; i++) xx20[i] = d13[i];
    for(i=0; i<nomax; i++) yy20[i] = d14[i];
    *nier = *i10;

    UNPROTECT(1);

    return ans;
}
