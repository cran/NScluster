      subroutine smplxIP(x, y, np, iskip1, ty, sclmu1, sclnu1, sclp1,
     &   sclc1, x22, eps, itmax, itmax1, ipmax, fn, mples, xinit, eps1,
     &   f, iter, nip, ipri, ipflag)
c
      include 'NScluster.h'
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
      parameter   (n=4)
c
      integer iskip
      double precision sclp,  sclc,  sclnu,  sclmu, fmin, x2
ccx      common/paramscl/sclp, sclc, sclnu, sclmu
      common/paramip/sclp, sclc, sclnu, sclmu
      common /fnmin/ fmin
cc      common / sizes / tx,ty
      common /skip/iskip
cc      common /fname/filea
      common/interval/x2
cc      character*50 filea
cx      real*8 sclp,sclc,sclnu,sclmu
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
cc      external   funct
cx      real*8     xinit(n,itmax1), dist, eps, f(itmax1)
      integer np, iskip1, itmax, itmax1, ipmax, iter, nip,
     1        ipri(ipmax), ipflag
      double precision x(np), y(np), ty, sclmu1, sclnu1, sclp1, sclc1,
     1                 x22, eps, fn(ipmax), mples(ipmax,n),
     2                 xinit(n,itmax1), eps1(itmax1), f(itmax1)
      double precision dist, rr(np**2), tx
      external   ipfunctMP
c
cx      dimension  x(np), y(np), rr(np**2)
cx      dimension  eps1(itmax1)
cx      dimension  ipri(ipmax), fn(ipmax), mples(ipmax,n)
c
      fmin = 1.d10
***************************
cc      open(13, file='IP.param')
cc      read(13,*) ix,iy,iz,iskip
cc      read(13,*) ty
cc      read(13,*) sclmu, sclnu, sclp, sclc
cc      read(13,1) filea
cc      read(13,*) x2
cc    1 format(a)
      tx = 1.0d0
      iskip = iskip1
      sclmu = sclmu1
      sclnu = sclnu1
      sclp = sclp1
      sclc = sclc1
      x2 = x22
***************************
cc      call input
      call input(x, y, np, tx, ty, rr, nn)
c
cc      n = 4
cc      if(iskip.eq.1) open(8, FILE='IP.simplex.print')
cc      if(iskip.eq.1000) open(8, FILE='IP1000.simplex.print')
c     write(8,*) '                 -log L           p             c  ',
c    &'           mu            nu'
cc      write(8,*) '           -log L                 mu            nu ',
cc     &'           p             c'
cc      close(8)
c
cc      xinit(1) = 1.0d0
cc      xinit(2) = 1.0d0
cc      xinit(3) = 1.0d0
cc      xinit(4) = 1.0d0
      xinit(1,1) = 1.0d0
      xinit(2,1) = 1.0d0
      xinit(3,1) = 1.0d0
      xinit(4,1) = 1.0d0
c
      dist = 0.1d0
cc      eps = 1.0d-3
c
cc      call simplx(xinit, n, funct, dist, eps, f)
      nip = 1
      call simplx(xinit, n, rr, nn, ipfunctMP, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip = nip - 1
c
cc      stop
      return
      end
c--------------------------------------------------------------------- c
cc      subroutine funct(n,b,fn)
      subroutine ipfunctMP(n, b, fn, r, nn, nip, jpri, ffn, mples,
     & ipmax, ipflag)
c-----------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c-----------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
      integer n, nn, nip, ipmax, jpri(ipmax), ipflag
      double precision b(n), fn, r(nn), ffn(ipmax), mples(ipmax,n)
c
      integer np, iskip
      double precision ff, aic, rmin, rmax, sclp, sclc, sclnu, sclmu,
     1                 ap, ac, fmin
      common/ddd/ff, aic
      common/range/rmin, rmax
ccx      common/paramscl/sclp, sclc, sclnu, sclmu
cxx      common/param/ap, ac
      common/paramip/sclp, sclc, sclnu, sclmu
      common/param1/ap, ac
      common/events/np
      common /fnmin/ fmin
      common /skip/iskip
c
      double precision pi, nu, mu, lambda, nu2pi, sum, dFr, Frmax, f,
     1                 aKrmax
c
      data pi/3.14159265358979d0/
cc      dimension b(4),g(4),h(4)
      pi = 3.14159265358979d0
      ap = b(1)**2 *sclp
      ac = b(2)**2 *sclc
cc      anu=b(3)**2 *sclnu
cc      amu=b(4)**2 *sclmu
      nu = b(3)**2 *sclnu
      mu = b(4)**2 *sclmu
c
cc      f1=0.0
      sum = 0.0
      lambda = mu*nu
      nu2pi = nu/2/pi
      ier = 0
cc      do 30 i=1,nn,iskip
!$omp parallel do private(dFr, f) reduction(+:sum)
      do 30 i = 1, nn
c     if(mod(i,10000).eq.0) write(6,*) 'i=',i, '/',nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
c
cc      call power(rr(i),dFr,Frmax)
cc      ramdai=(amu*anu)+anu/2/pi/rr(i)*dFr
      call ippowerMP(r(i), Frmax, dFr)
      f = lambda + nu2pi*dFr/r(i)
c
cc      if(ramdai.le.0.0) go to 190
cc      f1=f1+log(ramdai)
      if(f. le. 0.0) then
         ier = -1
      else
         sum = sum + log(f)
      end if
   30 continue
!$omp end parallel do
      if(ier. eq. -1) go to 190
c
cc      rmax=1.0d0/2
cc      call power(rmax,dFr,Frmax)
      call ippowerMP(rmax, Frmax, dFr)
c
cc      aKrmax=pi*(rmax**2)+anu*Frmax/(amu*anu)
      aKrmax = pi*(rmax**2) + Frmax/mu
c
cc      fn=iskip*f1-(amu*anu)*aKrmax*np
      fn = iskip*sum - lambda*aKrmax*np
      fn = -fn
      if(fmin .gt. fn) then
      fmin = fn
      ipri = 1
      else
      ipri = 0
      end if
cc      if(iskip.eq.1) open(8,FILE='IP.simplex.print',position='APPEND')
cc      if(iskip.eq.1000) 
cc     & open(8,FILE='IP1000.simplex.print',position='APPEND')
c     if(ipri.eq.0) write(8,2) 'testfn =',fn,ap,ac,amu,anu
c     if(ipri.eq.1) write(8,2) 'update =',fn,ap,ac,amu,anu
cc      close(8)
cx    3 format(1h ,110x,d18.10)
cx    1 format(1h ,7d18.10)
cx    2 format(1h , a, d23.15,4d14.7)
      ffn(nip) = fn
      mples(nip,1) = mu
      mples(nip,2) = nu
      mples(nip,3) = ap
      mples(nip,4) = ac
      if((ipflag. eq. 0) .or. (ipflag.eq.2)) return
      if(ipri .eq. 0) jpri(nip) = -1
      if(ipri .eq. 1) jpri(nip) = 1
      nip = nip + 1
      return
  190 continue
      fn = 1.d30
cc      write(6,2) 'fn190 =',fn,ap,ac,anu,amu
      return
      end
c--------------------------------------------------------------------- c
c--------------------------------------------------------------------- c
cc      subroutine power(ri,dFr,Fr)
      subroutine ippowerMP(ri, Fr, dFr)
c
cx      implicit real*8(a-h,o-z)
C     driver for routine qgaus

cx      implicit real*8(a-h,o-z)
C     driver for routine qgaus
      double precision ri, Fr, dFr
c      common/distance/r0
c      common/case/kk
      double precision x2, ap, ac
      common/interval/x2
cxx      common/param/ap, ac
      common/param1/ap, ac
      integer kk
      double precision r0
      common/distancep/r0
      common/casep/kk
!$omp threadprivate(/distancep/)
!$omp threadprivate(/casep/)
      double precision x1, ss, tt, uu, Freps1, Freps2,
     1                 xxm, xxr, xdx, xss, yy1, yy2, yy3, hMP1, hMP2
      common/param2/x1, ss, tt, uu, Freps1, Freps2
      common/param4/xxm, xxr, xdx, xss, yy1, yy2, yy3, hMP1, hMP2
!$omp threadprivate(/param2/)
!$omp threadprivate(/param4/)
c
cx      INTEGER NVAL
      integer NVAL
      double precision eps, pi, delta
c      REAL X1,X2
c      PARAMETER(X1=r0/2,X2=1.0,NVAL=10)
cx      INTEGER i
c      REAL dx,func,ss,x
cc      EXTERNAL func
      EXTERNAL ipfuncMP
c
      r0 = ri
c
cx   10 continue
c      write(6,*) 'input r0'
c      read(5,*) r0
c     write(6,*) 'input x2'
c     read(5,*) x2    
c     x2=0.3d0
      delta=0.001d0
c     open(1, FILE='power.plot.txt')
c     open(2, FILE='power.plot.eps0.txt')
c----------------------------------------------------------
c      do 15 j=1, 500
c     r0=delta*j
      x1 = r0/2
      nval = 100
        kk = 1
cc        call quad2d(X1,X2,ss)
        call qgausxMP(ipfuncMP, X1, X2)
        ss = xss
        kk = 2
cc        call quad2d(0.0d0,x1,tt)
        call qgausxMP(ipfuncMP, 0.0d0, x1)
        tt = xss
        kk = 3
cc        call quad2d(0.0d0,x1,uu)
        call qgausxMP(ipfuncMP, 0.0d0, x1)
        uu = xss
c        write(6,*) 'r=', r0
c        write(6,*) 'integral=', x1, x2, ss, tt, uu
c        write(6,*) 'F(r)=', 2*(ss+tt+uu)
cx11    continue
      Fr = 2*(ss+tt+uu)
c----------------------------------------------------------
      eps = 0.001d0
*     r0=delta*j+eps
      r0 = ri + eps
      x1 = r0/2
      nval = 100
        kk = 1
cc        call quad2d(X1,X2,ss)
        call qgausxMP(ipfuncMP, X1, X2)
        ss = xss
        kk = 2
cc        call quad2d(0.0d0,x1,tt)
        call qgausxMP(ipfuncMP, 0.0d0, x1)
        tt = xss
        kk = 3
cc        call quad2d(0.0d0,x1,uu)
        call qgausxMP(ipfuncMP, 0.0d0, x1)
        uu = xss
c        write(6,*) 'r=', r0
c        write(6,*) 'integral=', x1, x2, ss, tt, uu
c        write(6,*) 'F(r)=', 2*(ss+tt+uu)
      Freps1 = 2*(ss+tt+uu)
c----------------------------------------------------------
*     r0=delta*j-eps
      r0 = ri - eps
      if (r0. ne. 0) then
      x1 = r0/2
      nval = 100
        kk = 1
cc        call quad2d(X1,X2,ss)
        call qgausxMP(ipfuncMP, X1, X2)
        ss = xss
        kk = 2
cc        call quad2d(0.0d0,x1,tt)
        call qgausxMP(ipfuncMP, 0.0d0, x1)
        tt = xss
        kk = 3
cc        call quad2d(0.0d0,x1,uu)
        call qgausxMP(ipfuncMP, 0.0d0, x1)
        uu = xss
c        write(6,*) 'r=', r0
c        write(6,*) 'integral=', x1, x2, ss, tt, uu
c        write(6,*) 'F(r)=', 2*(ss+tt+uu)
      Freps2 = 2*(ss+tt+uu)
      end if
      if (r0 .eq. 0) then
      Freps2 = 0
      end if   
      dFr = (Freps1-Freps2)/(2*eps)
c     if(ap.gt.1.6) then
c     write(6,*) 'kkkk', dFr,Fr,Freps1,Freps2,ap,ri
c     endif
      pi = 3.14159265358979d0
*     nu=30
*     mu=1500
*     nu:the mean number of points per cluster
*     r.v. the number of points per cluster which is Poisson distributed  
*     aKr:Ripley's K-function
*     mu*aKr=mu*pi*(r**2)+nu*Fr
c     aKr=pi*(r**2)+nu*Fr/mu
c     daKr=2*pi*r+nu*dFr/mu
c     Palm intensity:alambda
c     alambda=mu*dKr/(2*pi*r) 
c     write(1,*) delta*j, Fr
c     write(2,*) delta*j, dFr 
c      write(6,*) delta*j, dFr
c15    continue
c      go to 10
c     close(1)
c     close(2)
      return
      END


cc      real*8 FUNCTION func(x,y)
cx      real*8 FUNCTION ipfuncMP(x,y)
      DOUBLE PRECISION FUNCTION ipfuncMP(x,y)
cx      implicit real*8(a-h,o-z)
cxx      common/distance/r0
cxx      common/case/kk
cd      double precision x, y
      double precision x, y
      integer kk
      double precision r0, ap, ac, qx, qy
      common/distancep/r0
      common/casep/kk
!$omp threadprivate(/distancep/)
!$omp threadprivate(/casep/)
c      common/param/p,c
cxx      common/param/ap, ac
      common/param1/ap, ac
      common/param3/qx, qy
!$omp threadprivate(/param3/)
c
c      REAL x,y
      double precision pi, ak, xyr0
      pi = 3.14159265358979d0
c     p=1.5d0
c     c=0.005d0
      ak = (ap-1) * (ac**(ap-1))
cc      if (kk.le.2) func=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
cx      if (kk.le.2) ipfuncMP=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
cx     &     *(ak/(x+ac)**ap)
cx     &     *(ak/(y+ac)**ap)
      ipfuncMP = 0
      if (kk. le. 2) then
         xyr0 = (x**2 + y**2 - r0**2)/(2*x*y)
         if (abs(xyr0) .le. 1.0d0) then
            ipfuncMP = (1/pi)*acos(xyr0)*(ak/(x+ac)**ap)*(ak/(y+ac)**ap)
         else
            ipfuncMP = 0
         end if
      end if
cc      if (kk.eq.3) func=1
      if (kk .eq. 3) ipfuncMP = 1
     &                  *(ak/(x+ac)**ap)
     &                  *(ak/(y+ac)**ap)
c      write(6,*) func,x,y
      return
      END
