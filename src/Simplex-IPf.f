      subroutine smplxIPf(x, y, np, iskip1, ty, sclmu1, sclnu1, sclp1,
     &   sclc1, x22, eps, itmax, itmax1, ipmax, fn, mples, xinit, eps1,
     &   f, iter, nip, ipri, ipflag)
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
      implicit real * 8 (a-h,o-z)
cc      parameter   (maxh=6, maxh5=maxh+5)
      parameter   (n=4)
c
      common/paramscl/sclp,sclc,sclnu,sclmu
      common /fnmin/ fmin
cc      common / sizes / tx,ty
      common /skip/iskip
cc      common /fname/filea
      common/interval/x2
cc      character*50 filea
      real*8 sclp,sclc,sclnu,sclmu
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, f
cc      external   funct
      real*8     xinit(n,itmax1), dist, eps, f(itmax1)
      external   ipfunct
c
      dimension  x(np), y(np), rr(np**2)
      dimension  eps1(itmax1)
      dimension  ipri(ipmax), fn(ipmax), mples(ipmax,n)
c
      fmin=1.d10
***************************
cc      open(13, file='IP.param')
cc      read(13,*) ix,iy,iz,iskip
cc      read(13,*) ty
cc      read(13,*) sclmu, sclnu, sclp, sclc
cc      read(13,1) filea
cc      read(13,*) x2
cc    1 format(a)
      tx=1.0d0
      iskip=iskip1
      sclmu=sclmu1
      sclnu=sclnu1
      sclp=sclp1
      sclc=sclc1
      x2=x22
***************************
cc      call input
      call input(x,y,np,tx,ty,rr,nn)
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
      call simplx(xinit, n, rr, nn, ipfunct, dist, eps, f,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
      if( (ipflag.eq.1) .or. (ipflag.eq.3) ) nip = nip-1
c
cc      stop
      return
      end
c -------------------------------------------------------------------- 
c
cc      subroutine funct(n,b,fn)
      subroutine ipfunct(n,b,fn,rr,nn,nip,jpri,ffn,mples,
     &                       ipmax,ipflag)
c-----------------------------------------------------------------------
c     likelihood function of the inverse power poisson process
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
      dimension rr(nn)
      common/ddd/ff,aic
      common/range/rmin,rmax
      common/paramscl/sclp,sclc,sclnu,sclmu
cc      common/param/ap,ac
      common/param/ap,ac,ak
      common/events/np
      common /fnmin/ fmin
      common /skip/iskip
c
      data pi/3.14159265358979d0/
cc      dimension b(4),g(4),h(4)
      dimension b(n),g(n),h(n)
c
      real*8   ffn(ipmax), mples(ipmax,n)
      integer  jpri(ipmax)
c
      pi=3.14159265358979d0
      ap=b(1)**2 *sclp
      ac=b(2)**2 *sclc
      anu=b(3)**2 *sclnu
      amu=b(4)**2 *sclmu
c-----
      ak=(ap-1)*(ac**(ap-1)) 
      amuanu=amu*anu
      anu2pi=anu/2/pi
c-----
c
      f1=0.0
cc      do 30 i=1,nn,iskip
      do 30 i=1,nn
c     if(mod(i,10000).eq.0) write(6,*) 'i=',i, '/',nn
cc      if(rr(i).le.rmin.or.rr(i).ge.rmax) go to 30
c
cc      call power(rr(i),dFr,Frmax)
cc      ramdai=(amu*anu)+anu/2/pi/rr(i)*dFr
      call ippower(rr(i),dFr,Frmax)
      ramdai=amuanu+anu2pi/rr(i)*dFr
c
      if(ramdai.le.0.0) go to 190
      f1=f1+log(ramdai)
   30 continue
c
cc      rmax=1.0d0/2
cc      call power(rmax,dFr,Frmax)
      call ippower(rmax,dFr,Frmax)
c
cc      aKrmax=pi*(rmax**2)+anu*Frmax/(amu*anu)
      aKrmax=pi*(rmax**2)+Frmax/amu
c
      fn=iskip*f1-(amu*anu)*aKrmax*np
      fn=-fn
      if(fmin.gt.fn) then
      fmin=fn
      ipri=1
      else
      ipri=0
      end if
cc      if(iskip.eq.1) open(8,FILE='IP.simplex.print',position='APPEND')
cc      if(iskip.eq.1000) 
cc     & open(8,FILE='IP1000.simplex.print',position='APPEND')
c     if(ipri.eq.0) write(8,2) 'testfn =',fn,ap,ac,amu,anu
c     if(ipri.eq.1) write(8,2) 'update =',fn,ap,ac,amu,anu
cc      close(8)
    3 format(1h ,110x,d18.10)
    1 format(1h ,7d18.10)
    2 format(1h , a, d23.15,4d14.7)
      ffn(nip) = fn
      mples(nip,1) = amu
      mples(nip,2) = anu
      mples(nip,3) = ap
      mples(nip,4) = ac
      if((ipflag.eq.0) .or. (ipflag.eq.2)) return
      if(ipri.eq.0) jpri(nip) = -1
      if(ipri.eq.1) jpri(nip) = 1
      nip = nip+1
      return
  190 continue
      fn=1.d30
cc      write(6,2) 'fn190 =',fn,ap,ac,anu,amu
      return
      end
c--------------------------------------------------------------------- c
c--------------------------------------------------------------------- c
cc      subroutine power(ri,dFr,Fr)
      subroutine ippower(ri,dFr,Fr)
c
      implicit real*8(a-h,o-z)
C     driver for routine qgaus
      common/distance/r0
      common/case/kk
      common/param/ap,ac,ak
      common/interval/x2
      INTEGER NVAL
c      REAL X1,X2
c      PARAMETER(X1=r0/2,X2=1.0,NVAL=10)
      INTEGER i
c      REAL dx,func,ss,x
cc      EXTERNAL func
      EXTERNAL ipfunc
c
      r0=ri
c
   10 continue
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
      x1=r0/2
      nval=100
        kk=1
cc        call quad2d(X1,X2,ss)
        call quad2d(ipfunc,X1,X2,ss)
        kk=2
cc        call quad2d(0.0d0,x1,tt)
        call quad2d(ipfunc,0.0d0,x1,tt)
        kk=3
cc        call quad2d(0.0d0,x1,uu)
        call quad2d(ipfunc,0.0d0,x1,uu)
c        write(6,*) 'r=', r0
c        write(6,*) 'integral=', x1, x2, ss, tt, uu
c        write(6,*) 'F(r)=', 2*(ss+tt+uu)
11    continue
      Fr=2*(ss+tt+uu)
c----------------------------------------------------------
      eps=0.001d0
*     r0=delta*j+eps
      r0=ri+eps
      x1=r0/2
      nval=100
        kk=1
cc        call quad2d(X1,X2,ss)
        call quad2d(ipfunc,X1,X2,ss)
        kk=2
cc        call quad2d(0.0d0,x1,tt)
        call quad2d(ipfunc,0.0d0,x1,tt)
        kk=3
cc        call quad2d(0.0d0,x1,uu)
        call quad2d(ipfunc,0.0d0,x1,uu)
c        write(6,*) 'r=', r0
c        write(6,*) 'integral=', x1, x2, ss, tt, uu
c        write(6,*) 'F(r)=', 2*(ss+tt+uu)
      Freps1=2*(ss+tt+uu)
c----------------------------------------------------------
*     r0=delta*j-eps
      r0=ri-eps
      x1=r0/2
      nval=100
        kk=1
cc        call quad2d(X1,X2,ss)
        call quad2d(ipfunc,X1,X2,ss)
        kk=2
cc        call quad2d(0.0d0,x1,tt)
        call quad2d(ipfunc,0.0d0,x1,tt)
        kk=3
cc        call quad2d(0.0d0,x1,uu)
        call quad2d(ipfunc,0.0d0,x1,uu)
c        write(6,*) 'r=', r0
c        write(6,*) 'integral=', x1, x2, ss, tt, uu
c        write(6,*) 'F(r)=', 2*(ss+tt+uu)
      Freps2=2*(ss+tt+uu)
      if (r0.eq.0) then
      Freps2=0
      end if   
      dFr=(Freps1-Freps2)/(2*eps)
c     if(ap.gt.1.6) then
c     write(6,*) 'kkkk', dFr,Fr,Freps1,Freps2,ap,ri
c     endif
      pi=3.14159265358979d0
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
      real*8 FUNCTION ipfunc(x,y)
      implicit real*8(a-h,o-z)
      common/distance/r0
      common/case/kk
c      common/param/p,c
cc      common/param/ap,ac
      common/param/ap,ac,ak
c      REAL x,y
      pi=3.14159265358979d0
c     p=1.5d0
c     c=0.005d0
cc      ak=(ap-1)*(ac**(ap-1)) 
cc      if (kk.le.2) func=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
      if (kk.le.2) ipfunc=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
     &     *(ak/(x+ac)**ap)
     &     *(ak/(y+ac)**ap)
cc      if (kk.eq.3) func=1
      if (kk.eq.3) ipfunc=1
     &                  *(ak/(x+ac)**ap)
     &                  *(ak/(y+ac)**ap)
c      write(6,*) func,x,y
      return
      END

