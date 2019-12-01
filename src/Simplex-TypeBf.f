      subroutine smplxBf(x, y, np, ty1, mu1, mu2, nu1, scls11,
     &    scls22, eps, itmax, itmax1, ipmax, fn, mples, xinit, eps1,
     &    f, iter, nip, ipri, ipflag)
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
      integer :: np, itmax, itmax1, ipmax, iter, nip, ipri(ipmax),
     1            ipflag
      real(8) :: x(np), y(np), ty1, mu1, mu2, nu1, scls11, scls22,
     1           eps, fn(ipmax), mples(ipmax,n), xinit(n,itmax1),
     2           eps1(itmax1), f(itmax1)
      real(8) :: sclmu, sclnu, scla1, scls1, scls2, fmin
      integer :: iskip
      real(8) :: dist, rr(np**2), tx, ty
c
cxx      common/paramscl/sclmu,sclnu,scla1,scls1,scls2
      common/bparam/sclmu,sclnu,scla1,scls1,scls2
      common /fnmin/ fmin
cc      common / sizes / tx,ty
cc      common /fname/filea
      common/skip/iskip
cc      character*50 filea
cx      real*8 sclmu,sclnu,scla1,scls1,scls2
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
cc      external   funct
      external   bfunctMP
      fmin = 1.d10
***************************
cc      open(13, file='TypeB.param')
cc      read(13,*) 
      tx = 1.0d0
cc      read(13,*) ty
cc      read(13,*) amu1, amu2, anu, scls1, scls2
      ty = ty1
      scls1 = scls11
      scls2 = scls22
cc      sclmu=amu1+amu2
cc      sclnu=anu1
cc      scla1=amu1/(amu1+amu2)
      sclmu = mu1 + mu2
      sclnu = nu1
      scla1 = mu1/(mu1+mu2)
cc      read(13,1) filea
cc    1 format(a)
***************************
cc      call input
      iskip = 1
      call input(x, y, np, tx, ty, rr, nn)
c
cc      n = 5
cc      open(8, FILE='TypeB.simplex.print')
cc      write(8,2) '            -log L            mu          nu ',
cc     &'         a           s1          s2'
cx    2 format(2a)
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
      call simplx(xinit, n, rr, nn, bfunctMP, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip  =  nip - 1
c
cc      stop
      return
      end
c -------------------------------------------------------------------- c
c
cc      subroutine funct(n,b,fn)
      subroutine bfunctMP(n, b, fn, r, nn, nip, jpri, ffn, mples,
     & ipmax, ipflag)
c-----------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c-----------------------------------------------------------------------
cc      implicit real * 8 (a-h,o-z)
      integer :: n, nn, nip, ipmax, jpri(ipmax), ipflag
      real(8) :: b(n), fn, r(nn), ffn(ipmax), mples(ipmax,n)
c
      integer :: np
      real(8) :: ff, aic, rmin, rmax, sclmu, sclnu, scla1, scls1,
     1           scls2, fmin
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
cx      dimension rr(nn)
      common/ddd/ff, aic
      common/range/rmin, rmax
cxx      common/paramscl/sclmu, sclnu, scla1, scls1, scls2
      common/bparam/sclmu, sclnu, scla1, scls1, scls2
      common/events/np
      common /fnmin/ fmin
c
      real(8) :: pi, eps, mu, nu, a1, s1, s2, se, sum, lambda, f,
     1           s124, s224, fs1, ainteg1, ainteg2
      data pi/3.14159265358979d0/
cc      dimension b(5),g(5),h(5)
      pi = 3.14159265358979d0
      eps = 0.1d-5
cc      amu=b(1)**2 * sclmu
cc      anu=b(2)**2 * sclnu
      mu = b(1)**2 * sclmu
      nu = b(2)**2 * sclnu
      a1 = b(3)**2 * scla1
      s1 = b(4)**2 * scls1
      s2 = b(5)**2 * scls2
      se = s2 - eps
c
      ier = 0
      sum = 0.0d0
      lambda = mu*nu
      s124 = 4*(s1**2)
      s224 = 4*(s2**2)
!$omp parallel do private(f) reduction(+:sum)
      do 30 i = 1, nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
cc      ramdai = (amu*anu)
cc     &        +(
cc     &          a1*(exp(-((rr(i))**2)/(4*((s1)**2))))/((s1)**2)
cc     &         +(1-a1)*(exp(-((rr(i))**2)/(4*((s2)**2))))/((s2)**2)
cc     &          )
cc     &           *(anu/(4*pi))
      f = lambda
     &       + (a1*(exp(-(r(i)**2)/s124))/((s1)**2)
     &       + (1-a1)*(exp(-(r(i)**2)/s224))/((s2)**2))*(nu/(4*pi))
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
cc      fs1=pi*(rmax**2)*(amu*anu)
      fs1 = pi*(rmax**2)*lambda
c
cc      ainteg1 = a1*anu*(1-(exp(-(rmax**2)/(4*((s1)**2)))))
cc      ainteg2 = (1-a1)*anu*(1-(exp(-(rmax**2)/(4*((s2)**2)))) 
      ainteg1 = a1*nu*(1-(exp(-(rmax**2)/s124)))
      ainteg2 = (1-a1)*nu*(1-(exp(-(rmax**2)/s224)))
      fn = sum-(fs1+ainteg1+ainteg2)*np
      fn = -fn
      ff = fn
      if(fmin .gt. fn) then
      fmin = fn
      ipri = 1
      else
      ipri = 0
      end if
cc      open(8, FILE='TypeB.simplex.print',position='APPEND')
cc      if(ipri.eq.0) write(8,2) 'testfn =',fn,amu,anu,a1,s1,s2
cc      if(ipri.eq.1) write(8,2) 'update =',fn,amu,anu,a1,s1,s2
cc      close(8)
cc      write(6,*) 'h(3)=', h(3)
      ffn(nip) = fn
cc      mples(nip,1) = amu
cc      mples(nip,2) = anu
      mples(nip,1) = mu
      mples(nip,2) = nu
      mples(nip,3) = a1
      mples(nip,4) = s1
      mples(nip,5) = s2
      if((ipflag.eq.0) .or. (ipflag.eq.2)) return
      if(ipri.eq.0) jpri(nip) = -1
      if(ipri.eq.1) jpri(nip) = 1
      nip = nip + 1
      return
cx    3 format(1h ,110x,d18.10)
cx    2 format(1h , a, d18.10,5d12.5)
cx    1 format(1h ,7d18.10)
c
  190 continue
      fn = 1.0d30
      return
      end
