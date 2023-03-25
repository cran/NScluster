      subroutine smplxC(x, y, np, ty1, mu1, mu2, nu1, nu2, scls11,
     &   scls22, eps, itmax, itmax1, ipmax, fn, mples, xinit, eps1, f,
     &   iter, nip, ipri, ipflag)
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
cc      parameter   (maxh=6, maxh5=maxh+5)
      parameter   (n=5)
c
cx      implicit real * 8 (a-h,o-z)
      integer np, itmax, itmax1, ipmax, iter, nip, ipri(ipmax), ipflag
      double precision x(np), y(np), ty1, mu1, mu2, nu1, nu2, scls11,
     1                 scls22, eps, fn(ipmax), mples(ipmax,n),
     2                 xinit(n,itmax1), eps1(itmax1), f(itmax1)
c
      integer iskip
      double precision scllam, scla, sclnu1, scls1, scls2, fmin, tx, ty
cxx      common/paramscl/scllam, scla, sclnu1, sclnu2, scls1, scls2
      common/cparam/scllam, scla, sclnu1, scls1, scls2
      common /fnmin/ fmin
      common / sizes / tx, ty
cc      common /fname/filea
      common /skip/iskip
c
cc      character*50 filea
cx      real*8 scllam,scla,sclnu1,sclnu2,scls1,scls2
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
      double precision dist, rr(np**2)
cc      external   funct
      external   cfunctMP
c
      fmin = 1.d10
***************************
cc      open(13, file='TypeC.param')
cc      read(13,*) 
      tx = 1.0d0
cc      read(13,*) ty
cc      read(13,*) amu1,amu2,anu1,anu2,scls1,scls2
      ty = ty1
      scls1 = scls11
      scls2 = scls22
c
cc      scllam = amu1*anu1+amu2*anu2
cc      scla = amu1*anu1/scllam
cc      sclnu1=anu1
      scllam = mu1*nu1 + mu2*nu2
      scla = mu1*nu1/scllam
      sclnu1 = nu1
cc      read(13,1) filea
cc    1 format(a)
***************************
cc      call input
      iskip = 1
      call input(x, y, np, tx, ty, rr, nn)
c
cc      n = 5
cc      open(8, FILE='TypeC.simplex.print')
cc      write(8,2) '            -log L            lambda      nu1',
cc     &'         a           s1          s2'
cx    2 format(2a)

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
      call simplx(xinit, n, rr, nn, cfunctMP, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip = nip - 1
c
cc      close(8)
cc      stop
      return
      end
c -------------------------------------------------------------------- c
cc      subroutine funct(n,b,fn)
      subroutine cfunctMP(n, b, fn, r, nn, nip, jpri, ffn, mples, ipmax,
     & ipflag)
c-----------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c-----------------------------------------------------------------------
cx      implicit real * 8 (a-h,o-z)
      integer n, nn, nip, ipmax, jpri(ipmax), ipflag
      double precision b(n), fn, r(nn), ffn(ipmax), mples(ipmax,n)
c
      integer np
      double precision ff, aic, rmin, rmax, scllam, scla, sclnu1,
     1                 scls1, scls2, fmin
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
      common/ddd/ff, aic
      common/range/rmin, rmax
cxx      common/paramscl/scllam, scla, sclnu1, sclnu2, scls1, scls2
      common/cparam/scllam, scla, sclnu1, scls1, scls2
      common/events/np
      common /fnmin/ fmin
c
      double precision pi, eps, alam, a, a1, a2, nu1, nu2, s1, s2, se,
     1                 sum, s124, s224, exp1rmax, exp2rmax, exp1i,
     2                 exp2i, f, fs1, ainteg1, ainteg2
      data pi/3.14159265358979d0/
cc      dimension b(6),g(6),h(6)
      pi = 3.14159265358979d0
      eps = 0.1d-5
      alam = b(1)**2 * scllam
      a = b(2)**2 * scla
cc      anu1=b(3)**2 * sclnu1
      nu1 = b(3)**2 * sclnu1
      s1 = b(4)**2 * scls1
      s2 = b(5)**2 * scls2 
cc      anu2=anu1 * s2/s1
      nu2 = nu1*s2/s1
      se = s2 - eps
c
      ier = 0
      sum = 0.0d0
      s124 = 4*(s1**2)
      s224 = 4*(s2**2)
      exp1rmax = 2*(s1**2)*(1 - exp(-(rmax**2)/s124))
      exp2rmax = 2*(s2**2)*(1 - exp(-(rmax**2)/s224))
      a1 = (a)*(1/(s1**2))*(nu1/(4*pi))
      a2 = (1-a)*(1/(s2**2))*(nu2/(4*pi))
!$omp parallel do private(exp1i,exp2i,f) reduction(+:sum)
      do 30 i = 1, nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
cc      exp1i = exp(-((rr(i))**2)/(4*((s1)**2)))
cc      exp2i = exp(-((rr(i))**2)/(4*((s2)**2)))
cc      exp1rmax = 2*((s1)**2)*(1 - exp(-(rmax**2)/(4*((s1)**2))))
cc      exp2rmax = 2*((s2)**2)*(1 - exp(-(rmax**2)/(4*((s2)**2))))
cc      ramdai = (alam) + (a)*
cc     &         (1/((s1)**2))*(anu1/(4*pi))*exp1i+
cc     &         (1-a)*(1/((s2)**2))*(anu2/(4*pi))*
cc     &         exp2i       
      exp1i = exp(-(r(i)**2)/s124)
      exp2i = exp(-(r(i)**2)/s224)
      f = alam + a1*exp1i + a2*exp2i
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
       fs1 = pi*(rmax**2)*(alam)
cc      ainteg1 = 2*pi*(a)*
cc     &          (1/((s1)**2))*(anu1/(4*pi))*exp1rmax
cc      ainteg2 = 2*pi*(1-a)*
cc     &          (1/((s2)**2))*(anu2/(4*pi))*exp2rmax
cc      fn=f1-(fs1+ainteg1+ainteg2)*np
      ainteg1 = 2*pi*a1*exp1rmax
      ainteg2 = 2*pi*a2*exp2rmax
      fn = sum - (fs1+ainteg1+ainteg2)*np
      fn = -fn
      ff = fn
      if(fmin .gt. fn) then
      fmin = fn
      ipri = 1
      else
      ipri = 0
      end if
cc      open(8, FILE='TypeC.simplex.print',position='APPEND')
cc      if(ipri.eq.0) write(8,2) 'testfn =',fn,alam,anu1,a,s1,s2
cc      if(ipri.eq.1) write(8,2) 'update =',fn,alam,anu1,a,s1,s2
cc      close(8)
      ffn(nip) = fn
      mples(nip,1) = alam
cc      mples(nip,2) = anu1
      mples(nip,2) = nu1
      mples(nip,3) = a
      mples(nip,4) = s1
      mples(nip,5) = s2
      if((ipflag.eq.0) .or. (ipflag.eq.2)) return
      if(ipri.eq.0) jpri(nip) = -1
      if(ipri.eq.1) jpri(nip) = 1
      nip = nip + 1
      return
cc      write(6,*) 'h(6)=', h(6)
cx    3 format(1h ,110x,d18.10)
cx    2 format(1h , a, d18.10,6d12.5)
cx    1 format(1h ,7d18.10)
c
  190 continue
      fn = 1.0d30
      return
      end
