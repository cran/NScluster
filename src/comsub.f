      include 'NScluster.h'
********************************************************************************
* common subroutine
c
c     Pois     ---   simThomas, simA, simB, simC, simIP
ccc     random   ---   simThomas, simA, simB, simC, simIP
c     init   ---   simThomas, simA, simB, simC, simIP
c     genrand   ---   simThomas, simA, simB, simC, simIP
c
c     input    ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     simplx   ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     minmax   ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     first    ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     center   ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     epslon(epsln)   ---    simplexThomas, simplexIP, simplexA, simplexB,
c                            simplexC
c     newsim   ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     update   ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
c     reduce   ---    simplexThomas, simplexIP, simplexA, simplexB, simplexC
ccc     quad2d   ---    simplexIP, simplexA, xqgausIP, xqgausA
ccc     qgausx   ---    simplexIP, simplexA, xqgausIP, xqgausA
ccc     qgausy   ---    simplexIP, simplexA, xqgausIP, xqgausA
c     quad2d   ---    xqgausIP, xqgausA
c     qgausx   ---    xqgausIP, xqgausA
c     qgausy   ---    xqgausIP, xqgausA

c
c-------------------------------------------------------------------------------
* common function
c
c     f   ---   simplexIP, simplexA
c     h   ---   simplexIP, simplexA
c     y1  ---   simplexIP, simplexA
c     y2  ---   simplexIP, simplexA
c     y3  ---   simplexIP, simplexA
c
********************************************************************************
*      *** Poisson random number ***
cc       subroutine Pois(ram,m)
cxx       subroutine Pois(ram,m,ix,iy,iz)
       subroutine Pois(ram,m)
cx       implicit real*8(a-h,o-z)
cxx      integer m, ix, iy, iz
      integer m
      double precision ram, alogu, u, random
cc       common ix,iy,iz
       alogu=ram
       m=0
cxx 10    u=random(ix,iy,iz)
 10    u=random()       
       alogu=alogu+log(u)
       if (alogu.gt.0) then
           m=m+1
           go to 10
       end if
c
       return
c
       end
********************************************************************************
      subroutine init(seed)
      integer seed, v(8)
      if ( seed .ge. 0 ) then
          call init_genrand64( seed )
      else
         call date_and_time( values=v )
         call init_genrand64( sum(v) )
         seed = sum(v)
      endif
      return
      end
********************************************************************************
cc      subroutine input
      subroutine input(x,y,n,tx,ty,rr,nn)
cx      implicit real * 8 (a-h,o-z)
      integer n, nn
      double precision x(n), y(n), tx, ty, rr(n**2)
cc      common/datpar/ nn
cc      common/xyod/rr(9234567),th(9234567)
cx      dimension rr(n**2)
      double precision ff, aic, rmin, rmax
      common/ddd/ff,aic
      common/range/rmin,rmax
cc      common / sizes / tx,ty
cc      dimension x(15000),y(15000)
cx      dimension x(n),y(n)
      integer npoi, iskip
      common/events/npoi
      common/skip/iskip
      double precision PI, t1, XX, YY, R2
      DATA PI /3.14159265358979D0/
c correlation lag band (Palm intensity)
      rmin=0.0d0
      rmax=1.d0/2
      t1=rmax
cc      CALL inlvy(X,Y,N,TX,TY,T1)
      npoi=n
c periodic boundary
      NN=0
      if(iskip.ne.1) go to 30
      DO 10 I=1,N
      DO 20 J=1,N
      if(J.eq.I) go to 20
      XX=X(J)-X(I)
      IF(XX.GT.TX/2) XX=-(TX-XX)
      IF(XX.LT.-TX/2) XX=(TX+XX)
      YY=Y(J)-Y(I)
      IF(YY.GT.TY/2) YY=-(TY-YY)
      IF(YY.LT.-TY/2) YY=(TY+YY)
      IF(ABS(XX).GT.T1.OR.ABS(YY).GT.T1) GO TO 20
      R2=SQRT(XX**2+YY**2)
      IF(R2.GT.T1) GO TO 20
c-----
      if(r2.le.rmin.or.r2.ge.rmax) go to 20
c-----
      NN=NN+1
cc      IF(NN.GT.9234567) STOP 229
      RR(NN)=R2
   20 CONTINUE
   10 CONTINUE
      return
c
   30 continue
      MM=0
      DO 11 I=1,N
      DO 21 J=1,N
      if(J.eq.I) go to 21
      XX=X(J)-X(I)
      IF(XX.GT.TX/2) XX=-(TX-XX)
      IF(XX.LT.-TX/2) XX=(TX+XX)
      YY=Y(J)-Y(I)
      IF(YY.GT.TY/2) YY=-(TY-YY)
      IF(YY.LT.-TY/2) YY=(TY+YY)
      IF(ABS(XX).GT.T1.OR.ABS(YY).GT.T1) GO TO 21
      R2=SQRT(XX**2+YY**2)
      IF(R2.GT.T1) GO TO 21
      NN=NN+1
cc      IF(NN.GT.9234567) STOP 229
cc      RR(NN)=R2
c-----
      if(r2.le.rmin.or.r2.ge.rmax) go to 21
c-----
      IF(MOD(NN,iskip).eq.1) THEN
         MM=MM+1
         RR(MM)=R2
      END IF
   21 CONTINUE
   11 CONTINUE
      NN=MM
c
cc      write(6,*) 'NN=', nn
c     prob=500**2*(0.5**2-0.1**2)*3.14159265
cc      write(6,1) 'nn,rmin,rmax ',nn,rmin,rmax
cc    1 format(1h ,a,1x,i10,3f15.5)
      return
      end
c
c
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
cc      subroutine simplx(xinit, n, funct, dist, eps, fr)
      subroutine simplx(xinit, n, rr, nn, funct, dist, eps, fr,
     & itmax, itmax1, iter, eps1, ipmax, nip, ipri, fn, mples, ipflag)
c -------------------------------------------------------------------- c
c -- parameter: ------------------------------------------------------ c
c ---- i/o xinit:      vector of initial and final x ----------------- c
c ---- i   n:          dimension of the vector x --------------------- c
c ---- i   funct:      minimize function ----------------------------- c
c ---- i   dist:       size of first simplex ------------------------- c
c ---- i   eps:        epsilon --------------------------------------- c
c ---- o   fr:         final value ----------------------------------- c
c -------------------------------------------------------------------- c
c -- this subroutine call next subroutine and function. -------------  c
c --------- first, minmax, newsim, center, reduce, update, epslon ---- c
c -------------------------------------------------------------------- c
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      external   funct
c
cc      integer    n
cc      real*8     xinit(maxh), dist, eps, fr
      integer n, nn, itmax, itmax1, iter, ipmax, nip, ipri(ipmax),
     1        ipflag
      double precision xinit(n,itmax1), rr(nn), dist, eps, fr(itmax1),
     1                 eps1(itmax1),fn(ipmax), mples(ipmax,n)
cx      real*8     xinit(n,itmax1), dist, eps, fr(itmax1)
c
cx      integer    iter, xh, xs, xl, x0, xr, xe, xc, j
cc      real*8     f(maxh5), x(maxh5,maxh)
cc      real*8     epslon
cc      real*8     eps1, alpha, beta, gamma
cx      real*8     f(n+5), x(n+5,n)
cx      real*8     epsln
cx      real*8     eps1(itmax1), alpha, beta, gamma
      integer xh, xs, xl, x0, xr, xe, xc, j
      double precision f(n+5), x(n+5,n), epsln, alpha, beta, gamma
c
cx      real*8     rr(nn)
cx      real*8     fn(ipmax), mples(ipmax,n)
cx      integer    ipri(ipmax)
c
      data       alpha/-1.0d0/, beta/0.5d0/, gamma/2.0d0/
c
      x0 = n + 2
      xr = n + 3
      xe = n + 4
      xc = n + 5

      x(1:(n+5),1:n) = 0
c
cc      call first(n, f, x, funct, xinit, dist)
      call first(n, f, x, rr, nn, funct, xinit, dist,
     &           ipmax, nip, ipri, fn, mples, ipflag)
c -------------------------------------------------------------------- c
cc      write(6,*)
cc      write(6,*) '   #    -log L     standard error  x(1)        x(2)',
cc     &'        x(3)'
      iter = 0
      it1 = 1
c
1     continue
c
c                                           ****   min max search   ****
      call minmax(n, f, xh, xs, xl)
c
cc      eps1 = epslon(n, f)
cc      write(6,31) iter, f(xl), eps1,(x(xl,i), i=1, n)
cc   31 format(i5, f15.3, 6d12.4)
cx   31 format(i5, f15.3, 6e12.4)
      eps1(it1) = epsln(n, f)
      fr(it1) = f(xl)
      do 10 j = 1, n
        xinit(j,it1) = x(xl, j)
10      continue
c
cc      if(iter .gt. 2 .and. eps1 .lt. eps) go to 2
      if(iter .gt. 2 .and. eps1(it1) .lt. eps) go to 2
c
c                                               ****   center set   ****
      call center(n, x, xh, x0)
c
c                                               ****   reflection   ****
cc      call newsim(n, f, x, funct, x0, xh, alpha, xr)
      call newsim(n, f, x, rr, nn, funct, x0, xh, alpha, xr,
     &               ipmax, nip, ipri, fn, mples, ipflag)
c
      if(f(xr) .le. f(xs)) then
c
        if(f(xr) .lt. f(xl)) then
c
c                                                ****   expansion   ****
cc          call newsim(n, f, x, funct, x0, xr, gamma, xe)
        call newsim(n, f, x, rr, nn, funct, x0, xr, gamma, xe,
     &                ipmax, nip, ipri, fn, mples, ipflag)
c

          if(f(xe) .lt. f(xr)) then

            call update(n, f, x, xe, xh)

          else

            call update(n, f, x, xr, xh)

            end if
c
        else
          call update(n, f, x, xr, xh)
          end if
c
      else
c
        if(f(xr) .lt. f(xh)) call update(n, f, x, xr, xh)
c
c                                              ****   contraction   ****
cc        call newsim(n, f, x, funct, x0, xh, beta, xc)
        call newsim(n, f, x, rr, nn, funct, x0, xh, beta, xc,
     &                 ipmax, nip, ipri, fn, mples, ipflag)
c
        if(f(xc) .lt. f(xh)) then
          call update(n, f, x, xc, xh)
c
        else
c
c                                                   ****   reduce   ****
cc          call reduce(n, f, x, funct, xl)
          call reduce(n, f, x, rr, nn, funct, xl,
     &                 ipmax, nip, ipri, fn, mples, ipflag)
c
          end if
c
        end if
c
      iter = iter + 1
      if(ipflag .ge. 2) it1 = it1+1
c
cc      if(iter .gt. 1000 ) go to 2
      if(iter .gt. itmax ) go to 2
c
c                                            ****   epsilon check   ****
      go to 1
c
2     continue
c
cc      do 10 j = 1, n
cc        xinit(j) = x(xl, j)
cc10      continue
cc      fr = f(xl)
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
      subroutine minmax(n, f, xh, xs, xl)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      integer    n, xh, xs, xl
cc      real*8     f(maxh5)
cx      real*8     f(n+5)
      double precision f(n+5)
c
      integer    i
cx      real*8     fmin1, fmax1, fmax2
      double precision fmin1, fmax1, fmax2
c
      if(f(1) .gt. f(2)) then
        fmax1 = f(1)
        xh = 1
        fmax2 = f(2)
        xs = 2
        fmin1 = f(2)
        xl = 2
      else
        fmax1 = f(2)
        xh = 2
        fmax2 = f(1)
        xs = 1
        fmin1 = f(1)
        xl = 1
        end if
c
      do 10 i = 3, n+1
c
        if(f(i) .gt. fmax1) then
          fmax2 = fmax1
          xs = xh
          fmax1 = f(i)
          xh = i
c
        else if(f(i) .gt. fmax2) then
          fmax2 = f(i)
          xs = i
          end if
c
        if(f(i) .lt. fmin1) then
          fmin1 = f(i)
          xl = i
          end if
c
10      continue
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
cc      subroutine first(n, f, x, funct, xinit, dist)
      subroutine first(n, f, x, rr, nn, funct, xinit, dist,
     &                    ipmax, nip, ipri, fn, mples, ipflag)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      external   funct
c
cx      integer    n
cc      real*8     f(maxh5), x(maxh5,maxh), xinit(maxh), xx(maxh), dist
cx      real*8     f(n+5), x(n+5,n), xinit(n), xx(n), dist
      integer n, nn,  ipmax, nip, ipri(ipmax), ipflag
      double precision f(n+5), x(n+5,n),  rr(nn), xinit(n), dist,
     1                 fn(ipmax), mples(ipmax,n), xx(n)
c
      integer    i, j
c
cx      real*8     rr(nn)
cx      real*8     fn(ipmax), mples(ipmax,n)
cx      integer    ipri(ipmax)
c
      do 10 i = 1, n + 1
c
        do 20 j = 1, n
          x(i,j) = xinit(j)
          xx(j) = x(i,j)
20        continue
c
        if(i .ne. 1 ) then
          x(i,i-1) = x(i,i-1) + dist
          xx(i-1) = x(i,i-1)
          end if
c
cc        call funct(n, xx, f(i))
        call funct(n, xx, f(i), rr, nn, nip, ipri, fn, mples,
     &              ipmax, ipflag)
c
10      continue
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
      subroutine center(n, x, xh, x0)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      integer    n, xh, x0
cc      real*8     x(maxh5,maxh)
cx      real*8     x(n+5,n)
      double precision x(n+5,n)
c
      integer    i, j
c
      do 10 j = 1, n
        x(x0,j) = 0.0d0
        do 20 i = 1, n+1
          if(i .ne. xh) x(x0,j) = x(x0,j) + x(i,j)
20        continue
        x(x0,j) = x(x0,j) / dble(n)
10      continue
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
cc      function   epslon(n, f)
cx      real*8 function   epsln(n, f)
      double precision function epsln(n, f)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      integer    n
cc      real*8     f(maxh5)
cx      real*8     f(n+5)
      double precision f(n+5)
c
      integer    i
cx      real*8     epslon, emean
      double precision epslon, emean
c
      emean = 0.0d0
      do 10 i = 1, n+1
        emean = emean + f(i)
10      continue
      emean = emean / dble(n+1)
c
      epslon = 0.0d0
      do 20 i = 1, n+1
        epslon = epslon + (f(i) - emean) ** 2
20      continue
cc      epslon = dsqrt(epslon) / dfloat(n+1)
      epsln = dsqrt(epslon) / dble(n+1)
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
cd  call newsim(n, f, x, rr, nn, funct, x0, xh, alpha, xr,ipmax, nip, ipri,
cd              fn, mples, ipflag)
cc      subroutine newsim(n, f, x, funct, x0, in, prm, out)
      subroutine newsim(n, f, x, rr, nn, funct, x0, in, prm, out,
     &                        ipmax, nip, ipri, fn, mples, ipflag)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      external   funct
c
cx      integer    n, x0, in, out
cc      real*8     f(maxh5), x(maxh5,maxh), xx(maxh), prm
cx      real*8     f(n+5), x(n+5,n), xx(n), prm
      integer  n, nn, x0, in, out, ipmax, nip, ipri(ipmax), ipflag
      double precision f(n+5), x(n+5,n), rr(nn), prm, fn(ipmax),
     1                 mples(ipmax,n), xx(n)
c
      integer    j
c
cx      real*8     rr(nn)
cx      real*8     fn(ipmax), mples(ipmax,n)
cx      integer    ipri(ipmax)
c
      do 10 j = 1, n
        x(out,j) = prm * x(in,j) + (1.0 - prm) * x(x0,j)
        xx(j) = x(out,j)
10      continue
cc      call funct(n, xx, f(out))
      call funct(n, xx, f(out), rr, nn, nip, ipri, fn, mples,
     &             ipmax, ipflag)
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
      subroutine update(n, f, x, in, out)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      integer    n, in, out
cc      real*8     f(maxh5), x(maxh5,maxh)
cx      real*8     f(n+5), x(n+5,n)
      double precision f(n+5), x(n+5,n)
c
      integer    j
c
      do 10 j = 1 , n
        x(out,j) = x(in,j)
10      continue
c
      f(out) = f(in)
c
      return
      end
c -------------------------------------------------------------------- c
c -------------------------------------------------------------------- c
cc      subroutine reduce(n, f, x, funct, xl)
      subroutine reduce(n, f, x, rr, nn, funct, xl,
     &                        ipmax, nip, ipri, fn, mples, ipflag)
c
cc      parameter   (maxh=6, maxh5=maxh+5)
c
      external   funct
c
cx      integer    n, xl
cc      real*8     f(maxh5), x(maxh5,maxh), xx(maxh)
cx      real*8     f(n+5), x(n+5,n), xx(n)
      integer n, nn, xl, ipmax, nip, ipri(ipmax), ipflag
      double precision f(n+5), x(n+5,n), rr(nn), fn(ipmax),
     &                 mples(ipmax,n)
c
      integer    i, j
c
cx      real*8     rr(nn)
cx      real*8     fn(ipmax), mples(ipmax,n)
cx      integer    ipri(ipmax)
      double precision xx(n)
c
      do 10 i = 1, n+1
        if(i .ne. xl) then
          do 20 j = 1, n
            x(i,j) = (x(i,j) + x(xl,j)) / 2.0d0
            xx(j) = x(i,j)
20          continue
c
cc          call funct(n, xx, f(i))
          call funct(n, xx, f(i), rr, nn, nip, ipri, fn, mples,
     &                ipmax, ipflag)
          end if
10      continue
c
      return
      end

cc      SUBROUTINE quad2d(x1,x2,ss)
      SUBROUTINE quad2d(func,x1,x2,ss)
cx      implicit real*8(a-h,o-z)
      double precision x1, x2, ss, r0
      common/distance/r0
c      REAL ss,x1,x2,h
cx      EXTERNAL h
      double precision, EXTERNAL :: h
CU    USES h,qgausx
cc      call qgausx(h,x1,x2,ss)
cx      EXTERNAL func
      double precision, EXTERNAL :: func
      call qgausx(func,h,x1,x2,ss)
      return
      END

cc      SUBROUTINE qgausx(func,a,b,ss)
      SUBROUTINE qgausx(func,h,a,b,ss)
cx      implicit real*8(a-h,o-z)
      double precision a, b, ss, r0
      common/distance/r0
c      REAL a,b,ss,func
cx      EXTERNAL func
      double precision, EXTERNAL :: func, h
      INTEGER j
c      REAL dx,xm,xr,w(5),x(5)
cx      dimension w(5), x(5)
      double precision w(5), x(5), xm, xr, dx
      SAVE w,x
      DATA w/.2955242247d0,.2692667193d0,.2190863625d0,.1494513491d0,
     *.0666713443d0/
      DATA x/.1488743389d0,.4333953941d0,.6794095682d0,.8650633666d0,
     *.9739065285d0/
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      ss=0
      do 11 j=1,5
        dx=xr*x(j)
cc        ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
        ss=ss+w(j)*(h(func,xm+dx)+h(func,xm-dx))
11    continue
      ss=xr*ss
      return
      END

      SUBROUTINE qgausy(func,a,b,ss)
cx      implicit real*8(a-h,o-z)
      double precision a, b, ss, r0
      common/distance/r0
c      REAL a,b,ss,func
cx      EXTERNAL func
      double precision, EXTERNAL :: func
      INTEGER j
c      REAL dx,xm,xr,w(5),x(5)
cx      dimension w(5), x(5)
      double precision f, w(5), x(5), xm, xr, dx
      SAVE w,x
      DATA w/.2955242247d0,.2692667193d0,.2190863625d0,.1494513491d0,
     *.0666713443d0/
      DATA x/.1488743389d0,.4333953941d0,.6794095682d0,.8650633666d0,
     *.9739065285d0/
      xm=0.5d0*(b+a)
      xr=0.5d0*(b-a)
      ss=0
      do 11 j=1,5
        dx=xr*x(j)
cc        ss=ss+w(j)*(f(xm+dx)+f(xm-dx))
        ss=ss+w(j)*(f(func,xm+dx)+f(func,xm-dx))
11    continue
      ss=xr*ss
      return
      END

cc      real*8 FUNCTION f(yy)
cx      real*8 FUNCTION f(func,yy)
      double precision FUNCTION f(func,yy)
cx      implicit real*8(a-h,o-z)
      double precision yy, r0, x, y, z, func
      common/distance/r0
c      REAL f,yy,func,x,y,z
      COMMON /xyz/ x,y,z
CU    USES func
      y=yy
      f=func(x,y)
      return
      END

cc      real*8 FUNCTION h(xx)
cx      real*8 FUNCTION h(func,xx)
      double precision FUNCTION h(func,xx)
cx      implicit real*8(a-h,o-z)
      double precision xx, r0, x, y, z, ss
      integer kk
      common/distance/r0
      common/case/kk 
c      REAL h,xx,func,y1,y2,x,y,z
cx      EXTERNAL func
      double precision, EXTERNAL :: func
      COMMON /xyz/ x,y,z
CU    USES g,qgausy,y1,y2
      double precision y1, y2, y3
c      REAL ss
      x=xx
      if (kk.eq.1) call qgausy(func,y1(x),y2(x),ss)
      if (kk.eq.2) call qgausy(func,y3(x),y2(x),ss)
      if (kk.eq.3) call qgausy(func,y1(x),y3(x),ss)
c      write(6,*) kk, ss,y1(x),y2(x),y3(x)
      h=ss
      return
      END

cx      real*8 FUNCTION y1(x)
      double precision FUNCTION y1(x)
cx      implicit real*8(a-h,o-z)
      double precision x, r0
      common/distance/r0
c      REAL x,y
      y1=x
      return
      END

cx      real*8 FUNCTION y2(x)
      double precision FUNCTION y2(x)
cx      implicit real*8(a-h,o-z)
      double precision x, r0
      common/distance/r0
c      REAL x,y
      y2=x+r0
      return
      END

cx      real*8 FUNCTION y3(x)
      double precision FUNCTION y3(x)
cx      implicit real*8(a-h,o-z)
      double precision x, r0
      common/distance/r0
c      REAL x,y
      y3=-x+r0
      return
      END
