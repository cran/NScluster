#include <R.h>
#include <Rdefines.h>
#include "NScluster.h"

extern void F77_NAME(setnthreadf)(int*, int*);

SEXP setnthread(SEXP nthread)
{
    int *i1,*i2;

    SEXP ans = R_NilValue, mthread = R_NilValue;
    int  *xmthread = NULL;

    i1 = INTEGER_POINTER(nthread);

    PROTECT(ans = allocVector(VECSXP, 1));
    SET_VECTOR_ELT(ans, 0, mthread = allocVector(INTSXP, 1));

    i2 = INTEGER_POINTER(mthread);

    F77_CALL(setnthreadf) (i1,i2);

    xmthread = INTEGER(mthread);

    *xmthread = *i2;

    UNPROTECT(1);

    return ans;
}
