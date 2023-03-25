/*
*  NScluster : 
*  Simulation and estimation of the Neyman-Scott Type spatial cluster models
*  Copyright (C) 2010    The Institute of Statistical Mathematics
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*
*  ismrp at grp.ism.ac.jp
*/

#include <R.h>
#include <Rinternals.h>
#include <libintl.h>

#define _(String) (String)

/* Fortran : */

void F77_NAME(xqgausip)(double *x, double *y, int *np, double *delta, double *ty, double *x2,
     double *amu, double *anu, double *p, double *c, int *m, int *jmax,
     double *palm, double *palm1);
                     
void F77_NAME(palmt)(double *x, double *y, int *np, double *delta, double *ty, double *amu,
     double *anu, double *v, int *m, int *jmax, double *palm, double *palm1);
                 
void F77_NAME(xqgausa)(double *x, double *y, int *np, double *delta, double *ty, double *x2,
     double *amu, double *anu, double *a, double *s1, double *s2, int *m,
     int *jmax, double *palm, double *palm1);

void F77_NAME(palmb)(double *x, double *y, int *np, double *delta, double *ty, double *amu,
     double *anu, double *a, double *s1, double *s2, int *m, int *jmax,
     double *palm, double *palm1);

void F77_NAME(palmc)(double *x, double *y, int *np, double *delta, double *ty, double *alam,
     double *anu1, double *a, double *s1, double *s2, int *m, int *jmax,
     double *palm, double *palm1);
					 

void F77_NAME(smplxip)(double *x, double *y, int *np, int *skip, double *ty, double *sclmu,
     double *sclnu, double *sclp, double *sclc, double *x2, double *eps,
     int *itmax, int *itmax1, int *ipmax, double *fn, double *mple, double *xx,
     double *std, double *f, int *itr, int *nip, int *ipr, int *ipflg);
							  
void F77_NAME(smplxthom)(double *x, double *y, int *np, double *ty, double *sclmu, double *sclnu,
     double *scls, double *eps, int *itmax, int *itmax1, int *ipmax, double *fn,
     double *mple, double *xx, double *std, double *f, int *itr, int *nip,
     int *ipr, int *ipflg);

void F77_NAME(smplxa)(double *x, double *y, int *np, int *skip, double *ty, double *sclmu,
     double *sclnu, double *scla, double *scls1, double *scls2, double *x2,
     double *eps, int *itmax, int *itmax1, int *ipmax, double *fn, double *mple,
     double *xx, double *std, double *f, int *itr, int *nip, int *ipr,
     int *ipflg);

void F77_NAME(smplxb)(double *x, double *y, int *np, double *ty, double *mu1, double *mu2,
     double *nu, double *s1, double *s2, double *eps, int *itmax, int *itmax1,
     int *ipmax, double *fn, double *mple, double *xx, double *std, double *f,
     int *itr, int *nip, int *ipr, int *ipflg);
					  
void F77_NAME(smplxc)(double *x, double *y, int *np, double *ty, double *sclmu1, double *sclmu2,
     double *sclnu1, double *sclnu2, double *scls1, double *scls2, double *eps,
     int *itmax,  int *itmax1, int *ipmax, double *fn, double *mple, double *xx,
     double *std, double *f, int *itr, int *nip, int *ipr, int *ipflg);


void F77_NAME(simip)(int *ix, double *ty, double *amu, double *anu, double *p, double *c,
     int *npts, int *ncl, double *x, double *y, double *xcl, double *ycl,
     int *pmax, int *omax, int *ier);

void F77_NAME(simthom)(int *ix, double *ty, double *amu, double *anu, double *sig, int *npts,
     int *ncl, double *x, double *y, double *xcl, double *ycl, int *pmax,
     int *omax, int *ier);

void F77_NAME(sima)(int *ix, double *ty, double *amu, double *anu, double *a, double *sig1,
     double *sig2, int *npts, int *ncl, double *x, double *y, double *xcl,
     double *ycl, int *pmax, int *omax, int *ier);

void F77_NAME(simb)(int *ix, double *ty, double *amu, double *anu, double *a, double *sig1,
     double *sig2, int *m1, int *ncl1, double *x1, double *y1, double *xx1,
     double *yy1, int *m2, int *ncl2, double *x2, double *y2, double *xx2,
     double *yy2, int *pmax, int *omax, int *ier);
		  
void F77_NAME(simc)(int *ix, double *ty, double *amu1, double *amu2, double *anu1, double *anu2,
     double *sig1, double *sig2, int *m1,  int *ncl1, double *x1, double *y1,
     double *xx1, double *yy1, int *m2, int *ncl2, double *x2, double *y2,
     double *xx2, double *yy2, int *pmax, int *omax, int *ier);

