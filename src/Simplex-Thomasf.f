      subroutine smplxThomf(x, y, np, ty, sclmu1, sclnu1, scls1,
     &   eps, itmax, itmax1, ipmax, fn, mples, xinit, eps1, f, iter,
     &   nip, ipri, ipflag)
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
c program starts. ----------------------------------------------------
c
cc      parameter   (maxh=6, maxh5=maxh+5)
      parameter   (n=3)
c
cx      implicit real*8 (a-h, o-z)
      integer :: np, itmax, itmax1, ipmax, iter, nip, ipri(ipmax),
     1           ipflag
      real(8) :: x(np), y(np), ty, sclmu1, sclnu1, scls1, eps,
     1           fn(ipmax), mples(ipmax,n), xinit(n,itmax1),
     2           eps1(itmax1), f(itmax1)
      real(8) :: sclmu, sclnu, scls, fmin
ccx      common/paramscl/sclmu, sclnu, scls
      common/tparam/sclmu, sclnu, scls
cc      common / sizes / tx,ty
      common /fnmin/ fmin
cc      common /fname/filea
cc      character*50 filea
cx      real*8 sclmu,sclnu,scls
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
cc      external   funct
cx      real*8     xinit(n,itmax1), dist, eps, f(itmax1)
      external   functMP
c
cx      dimension  x(np), y(np), rr(np**2)
cx      dimension  eps1(itmax1)
cx      dimension  ipri(ipmax), fn(ipmax), mples(ipmax,n)
      integer :: iskip
      common /skip/iskip
      real(8) :: dist, rr(np**2), tx
c
      fmin = 1.d10
***************************
cc      open(13, file='Thomas.param')
cc      read(13,*) 
      tx = 1.0d0
cc      read(13,*) ty
cc      read(13,*) sclmu, sclnu, scls
      sclmu = sclmu1
      sclnu = sclnu1
      scls = scls1
cc      read(13,1) filea
cc    1 format(a)
***************************
c
cc      call input
      iskip = 1
      call input(x, y, np, tx, ty, rr, nn)
c
cc      n = 3
c
cc      open(8, FILE='Thomas.simplex.print')
cc      write(8,*) '           -log L           mu           nu     ',
cc     &'     sigma'
cc      close(8) 
c
      dist = 0.1d0
cc      eps = 1.0d-3
c
cc      xinit(1) = 1.0d0
cc      xinit(2) = 1.0d0
cc      xinit(3) = 1.0d0
      xinit(1,1) = 1.0d0
      xinit(2,1) = 1.0d0
      xinit(3,1) = 1.0d0
cc      call simplx(xinit, n, funct, dist, eps, f)
      nip = 1
      call simplx(xinit, n, rr, nn, functMP, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip = nip - 1
c
cc      stop
      return
      end
c--------------------------------------------------------------------- c
cc      subroutine funct(n,b,fn)
      subroutine functMP(n, b, fn, r, nn, nip, jpri, ffn, mples, ipmax,
     & ipflag)
c---------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c---------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
      integer :: n, nn, nip, ipmax, jpri(ipmax), ipflag
      real(8) :: b(n), fn, r(nn), ffn(ipmax), mples(ipmax,n)
c
      integer :: np
      real(8) :: ff, aic, rmin, rmax, sclmu, sclnu, scls, fmin
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
      common/ddd/ff, aic
      common/range/rmin, rmax
cxx      common/paramscl/sclmu, sclnu, scls
      common/tparam/sclmu, sclnu, scls
      common/events/np
      common /fnmin/ fmin
c
      real(8) :: pi, mu, nu, s, sum, lambda, s24, nupis, f, fs1, ainteg
      data pi/3.14159265358979d0/
cc      dimension b(3),g(3),h(3)
      pi = 3.14159265358979d0
cc      eps=0.1d-5
cc      amu=b(1)**2 * sclmu
cc      anu=b(2)**2 * sclnu
      mu = b(1)**2 * sclmu
      nu = b(2)**2 * sclnu
      s = b(3)**2 * scls
c
      ier = 0
      sum = 0.0d0
      lambda = mu*nu
      s24 = 4*s**2
      nupis = nu/(pi*s24)
!$omp parallel do private(f) reduction(+:sum)
      do 30 i = 1, nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
cc      ramdai=(amu*anu)+(exp(-((rr(i))**2)/(4*((s)**2))))
cc     &      *(anu/(4*pi*((s)**2)))
      f = lambda + nupis * (exp(-(r(i)**2)/s24))
cc      if(ramdai.le.0.0) go to 190
cc         f1=f1+log(ramdai)
      if(f .le. 0.0) then
         ier = -1
      else
         sum = sum + log(f)
      end if
   30 continue
!$omp end parallel do
      if(ier .eq. -1) go to 190
c
cc      fs1=pi*(rmax**2)*(amu*anu)
cc      ainteg=anu*(1-(exp(-(rmax**2)/(4*((s)**2)))))
cc      fn=f1-(fs1+ainteg)*np
      fs1 = pi*(rmax**2)*lambda
      ainteg = nu*(1-(exp(-(rmax**2)/(4*(s**2)))))
      fn = sum - (fs1+ainteg)*np
      fn = -fn
      ff = fn
      if(fmin .gt. fn) then
      fmin = fn
      ipri = 1
      else
      ipri = 0
      end if
cc      open(8, FILE='Thomas.simplex.print'
cc     &      ,position='APPEND')
cc      if(ipri.eq.0) write(8,2) 'testfn =',fn,amu,anu,s
cc      if(ipri.eq.1) write(8,2) 'update =',fn,amu,anu,s
cc      close(8)
      ffn(nip) = fn
cc      mples(nip,1) = amu
cc      mples(nip,2) = anu
      mples(nip,1) = mu
      mples(nip,2) = nu
      mples(nip,3) = s
      if((ipflag.eq.0) .or. (ipflag.eq.2)) return
      if(ipri.eq.0) jpri(nip) = -1
      if(ipri.eq.1) jpri(nip) = 1
      nip = nip+1
      return
cc      write(6,*) 'h(3)=', h(3)
cx    3 format(1h ,110x,d18.10)
cx    2 format(1h , a, d18.10,5d12.5)
cx    1 format(1h ,7d18.10)
c
  190 continue
      fn = 1.0d30
      return
      end
