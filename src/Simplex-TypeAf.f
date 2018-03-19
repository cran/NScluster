      subroutine smplxAf(x, y, np, iskip1, ty, sclmu1, sclnu1, scla1,
     &    scls11, scls22, x22, eps, itmax, itmax1, ipmax, fn, mples,
     &    xinit, eps1, f, iter, nip, ipri, ipflag)
c
      include 'NScluster_f.h'
c
c simplx:  simplex minimization subroutine.
c minmax:  called by subroutine simplx. 
c first:             "
c centor:            "
c newsim:            "
c update:            "
c reduce:            "
c epslon:            "
c
c define parameter.
c maxh:    max dimension of vector x.
c
c procedure to use this program.
c 1. save this mail to a file. (save file-name.f)
c 2. delete the comments. (vi file name.f)
c 3. change the define parameter and function.
c
c program starts. ------------------------------------------------------
c
cx      implicit real * 8 (a-h,o-z)
cc      parameter   (maxh=6, maxh5=maxh+5)
      parameter   (n=5)
c
      integer :: np, iskip1, itmax, itmax1, ipmax, iter, nip, 
     1           ipri(ipmax), ipflag
      real(8) :: x(np), y(np), ty, sclmu1, sclnu1, scla1, scls11, 
     1           scls22, x22, eps, fn(ipmax), mples(ipmax,n), 
     2           xinit(n,itmax1), eps1(itmax1), f(itmax1)
      integer :: iskip
      real(8) :: scla, scls1, scls2, sclnu, sclmu, fmin, x2
      real(8) :: dist, rr(np**2), tx
      common/paramscl/scla, scls1, scls2, sclnu, sclmu
      common /fnmin/ fmin
cc      common /fname/filea
cc      character*50 filea
cc      common / sizes / tx,ty
      common/skip/iskip
      common/interval/x2
cc      real*8 scla,scls1,scls2,sclnu,sclmu
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
cc      external   funct
      external   afunctMP
c
      fmin = 1.d10
***************************
cc      open(13, file='TypeA.param')
cc      read(13,*) ix,iy,iz,iskip
      tx = 1.0d0
cc      read(13,*) ty
cc      read(13,*) sclmu, sclnu, scla, scls1, scls2
cc      read(13,1) filea
cc    1 format(a)
cc      read(13,*) x2
      iskip = iskip1
      sclmu = sclmu1
      sclnu = sclnu1
      scla = scla1
      scls1 = scls11
      scls2 = scls22
      x2 = x22
***************************
cc      call input
      call input(x, y, np, tx, ty, rr, nn)
cc      n = 5
cc      if(iskip.eq.1) open(8, FILE='TypeA.simplex.print')
cc      if(iskip.eq.1000) open(8, FILE='TypeA1000.simplex.print')
cc      write(8,9) '            -log L            mu            nu ',
cc     &'           a             s1            s2'
cx    9 format(2a)
cc      close(8)
c
cc      xinit(1) = 1.0d0
cc      xinit(2) = 1.0d0
cc      xinit(3) = 1.0d0
cc      xinit(4) = 1.0d0
cc      xinit(5) = 1.0d0
      xinit(1,1) = 1.0d0
      xinit(2,1) = 1.0d0
      xinit(3,1) = 1.0d0
      xinit(4,1) = 1.0d0
      xinit(5,1) = 1.0d0
c
      dist = 0.1d0
cc      eps = 1.0d-3
c
cc      call simplx(xinit, n, funct, dist, eps, f)
      nip = 1
      call simplx(xinit, n, rr, nn, afunctMP, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip = nip - 1
c
cc      stop
      return
      end
c--------------------------------------------------------------------- c
cc      subroutine funct(n,b,fn)
      subroutine afunctMP(n, b, fn, r, nn, nip, jpri, ffn, mples,
     & ipmax, ipflag)
c---------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c---------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
      integer :: n, nn, nip, ipmax, jpri(ipmax), ipflag
      real(8) :: b(n), fn, r(nn), ffn(ipmax), mples(ipmax,n)
c
      integer :: np, iskip
      real(8) :: ff, aic, rmin, rmax, scla, scls1, scls2, sclnu, sclmu,
     1           a, s1,s2, fmin
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
      common/ddd/ff, aic
      common/range/rmin, rmax
      common/paramscl/scla, scls1, scls2, sclnu, sclmu
      common/param/a, s1, s2
      common/events/np
      common /fnmin/ fmin
      common/skip/iskip
c
      real(8) :: pi, nu, mu, lambda, nu2pi, sum, dFr, Frmax, f, aKrmax
c
      data pi/3.14159265358979d0/
cc      dimension b(5),g(5),h(5)
      pi = 3.14159265358979d0
cc      amu=b(1)**2 *sclmu
cc      anu=b(2)**2 *sclnu
      mu = b(1)**2 *sclmu
      nu = b(2)**2 *sclnu
      a = b(3)**2 *scla
      s1 = b(4)**2 *scls1
      s2 = b(5)**2 *scls2 
c
      ier = 0
      sum=0.0d0
      lambda = nu*mu
      nu2pi = nu/2/pi
cc      do 30 i=1,nn,iskip
!$omp parallel do private(dFr) reduction(+:sum)
      do 30 i = 1, nn
c     if(mod(i,10000).eq.0) write(6,*) 'i=',i, '/',nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
c
cc      call power(rr(i),dFr,Frmax)
      call apowerMP(r(i), Frmax, dFr)
cc      ramdai=(amu*anu)+anu/2/pi/rr(i)*dFr
      f = lambda + nu2pi*dFr/r(i)
c
cc      if(ramdai.le.0.0) go to 190
cc      f1=f1+log(ramdai)
      if(f .le. 0.0) then
         ier = -1
      else
         sum = sum + log(f)
      end if
   30 continue
!$omp end parallel do
      if(ier .eq. -1) go to 190
c
cc      rmax=1.0d0/2
cc      call power(rmax,dFr,Frmax)
      call apowerMP(rmax, Frmax, dFr)
c
cc      aKrmax=pi*(rmax**2)+anu*Frmax/(amu*anu)
      aKrmax = pi*(rmax**2) + Frmax/mu
c
cc      fn=iskip*f1-(amu*anu)*(aKrmax)*np
      fn = iskip*sum - lambda*(aKrmax)*np
      fn = -fn
      if(fmin .gt. fn) then
      fmin = fn
      ipri = 1
      else
      ipri = 0
      end if
cc      if(iskip.eq.1) 
cc     &  open(8,FILE='TypeA.simplex.print',position='APPEND')
cc      if(iskip.eq.1000) 
cc     &  open(8, FILE='TypeA1000.simplex.print',position='APPEND')
cc      if(ipri.eq.0) write(8,2) 'testfn =',fn,amu,anu,a,s1,s2
cc      if(ipri.eq.1) write(8,2) 'update =',fn,amu,anu,a,s1,s2
cc      close(8)
cx    3 format(1h ,110x,d18.10)
cx    1 format(1h ,7d18.10)
cx    2 format(1h , a, d18.10,5d14.7)
      ffn(nip) = fn
cc      mples(nip,1) = amu
cc      mples(nip,2) = anu
      mples(nip,1) = mu
      mples(nip,2) = nu
      mples(nip,3) = a
      mples(nip,4) = s1
      mples(nip,5) = s2
      if((ipflag.eq.0) .or. (ipflag.eq.2)) return
      if(ipri.eq.0) jpri(nip) = -1
      if(ipri.eq.1) jpri(nip) = 1
      nip = nip + 1
      return
  190 continue
      fn = 1.d30
cc      write(8,2) 'fn190 =',fn,amu,anu,a,s1,s2
      return
      end
c
c--------------------------------------------------------------------- c
c--------------------------------------------------------------------- c
cc      subroutine power(ri,dFr,Fr)
      subroutine apowerMP(ri, Fr, dFr)
c
cx      implicit real*8(a-h,o-z)
C     driver for routine qgaus
cxx      common/distance/r0
cxx      common/case/kk
      real(8) :: ri, Fr, dFr
      real(8) :: x2, a, s1, s2
      common/interval/x2
***************************************************************** 
      common/param/a, s1, s2
***************************************************************** 
      real(8) :: r0
      integer :: kk
      common/distancep/r0
      common/casep/kk
!$omp threadprivate(/distancep/)
!$omp threadprivate(/casep/)
      real(8) :: x1, ss, tt, uu, Freps1, Freps2,
     1           xxm, xxr, xdx, xss, yy1, yy2, yy3, hMP1, hMP2
      common/param2/x1, ss, tt, uu, Freps1, Freps2
      common/param4/xxm, xxr, xdx, xss, yy1, yy2, yy3, hMP1, hMP2
!$omp threadprivate(/param2/)
!$omp threadprivate(/param4/)
c
cx      INTEGER NVAL
      integer :: NVAL
      real(8) :: delta, eps
c      REAL X1,X2
c      PARAMETER(X1=r0/2,X2=1.0,NVAL=10)
cx      INTEGER i
c      REAL dx,func,ss,x
cc      EXTERNAL func
      EXTERNAL afuncMP
c
      r0 = ri
c
cx   10 continue
      delta = 0.001d0
c----------------------------------------------------------
      x1 = r0/2
      nval = 100
        kk = 1
cc        call quad2d(X1,X2,ss)
        call qgausxMP(afuncMP, x1, x2)
        ss = xss
        kk = 2
cc        call quad2d(0.0d0,x1,tt)
        call qgausxMP(afuncMP, 0.0d0, x1)
        tt = xss
        kk = 3
cc        call quad2d(0.0d0,x1,uu)
        call qgausxMP(afuncMP, 0.0d0, x1)
        uu = xss
      Fr = 2*(ss+tt+uu)
c----------------------------------------------------------
      eps = 0.00001d0
*     r0=delta*j+eps
      r0 = ri + eps
      x1 = r0/2
      nval = 100
        kk = 1
cc        call quad2d(X1,X2,ss)
        call qgausxMP(afuncMP, X1, X2)
        ss = xss
        kk = 2
cc        call quad2d(0.0d0,x1,tt)
        call qgausxMP(afuncMP, 0.0d0, x1)
        tt = xss
        kk = 3
cc        call quad2d(0.0d0,x1,uu)
        call qgausxMP(afuncMP, 0.0d0, x1)
        uu = xss
      Freps1 = 2*(ss+tt+uu)
c----------------------------------------------------------
*     r0=delta*j-eps
      r0 = ri - eps
      x1 = r0/2
      nval = 100
        kk = 1
cc        call quad2d(X1,X2,ss)
        call qgausxMP(afuncMP, X1, X2)
        ss = xss
        kk = 2
cc        call quad2d(0.0d0,x1,tt)
        call qgausxMP(afuncMP, 0.0d0, x1)
        tt = xss
        kk = 3
cc        call quad2d(0.0d0,x1,uu)
        call qgausxMP(afuncMP, 0.0d0, x1)
        uu = xss
      Freps2 = 2*(ss+tt+uu)
      if (r0 .eq. 0) then
      Freps2 = 0
      end if   
      dFr = (Freps1-Freps2)/(2*eps)
      return
      END

cc      real*8 FUNCTION func(x,y)
cx      real*8 FUNCTION afuncMP(x,y)
      DOUBLE PRECISION FUNCTION afuncMP(x, y)
cx      implicit real*8(a-h,o-z)
cd      real(8) :: x, y
      real(8) :: x, y
      integer :: kk
      real(8) :: r0, a, s1, s2, qx, qy
cxx      common/distance/r0
cxx      common/case/kk
c      common/param/p,c
      common/distancep/r0
      common/casep/kk
!$omp threadprivate(/distancep/)
!$omp threadprivate(/casep/)
***************************************************************** 
      common/param/a, s1, s2
***************************************************************** 
      common/param3/qx, qy
!$omp threadprivate(/param3/)
c      REAL x,y
      real(8) :: pi, xyr0
      pi = 3.14159265358979d0
c     p=1.5d0
c     c=0.005d0
C      ak=(ap-1)*(ac**(ap-1)) 
c      write(6,*) ak
       qx = a/((s1)**2)*x*exp(-x**2/(2*(s1)**2))
     &     + (1-a)/((s2)**2)*x*exp(-x**2/(2*(s2)**2))
       qy = a/((s1)**2)*y*exp(-y**2/(2*(s1)**2))
     &     + (1-a)/((s2)**2)*y*exp(-y**2/(2*(s2)**2))
c
*      if (kk.le.2) func=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
*
*     &     *(ak/(x+ac)**ap)
*
*     &     *(ak/(y+ac)**ap)
*
*      if (kk.eq.3) func=1
*
*     &                  *(ak/(x+ac)**ap)
*
*     &                  *(ak/(y+ac)**ap)
*****************************************************************
***Mixed Gaussian***
***************************************************************** 
cc      if (kk.le.2) func=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
cx      if (kk.le.2) afuncMP=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
cx     &     *qx
cx     &     *qy
      if (kk.le.2) then
         xyr0 = (x**2 + y**2 - r0**2)/(2*x*y)
         if (abs(xyr0) .le. 1.0d0) then
            afuncMP = (1/pi)*acos(xyr0)*qx*qy
         else
            afuncMP = 0
         end if
      end if
cc      if (kk.eq.3) func=1
      if (kk.eq.3) afuncMP = 1
     &                  *qx
     &                  *qy
*****************************************************************
c      write(6,*) func,x,y
      return
      END
