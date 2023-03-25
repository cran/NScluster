cc      PROGRAM xqgaus
      subroutine xqgausa(x,y,np,delta,ty1,x2,amu,anu,aa,ss1,ss2,m,jmax,
     $                   palm,palm1)
c
      include 'NScluster.h'
c
cx      implicit real*8(a-h,o-z)
C     driver for routine qgaus
      integer np,m, jmax
      double precision x(np), y(np), delta, ty1, x2, amu(m), anu(m),
     1                 aa(m), ss1(m), ss2(m), palm(jmax), palm1(jmax,m)
      integer kk 
      double precision r0, a, s1, s2, tx, ty
      common/distance/r0
      common/case/kk
      common/av/a,s1,s2
cc      common/events/np
      common / sizes / tx,ty
cc      dimension  x(2000),y(2000), RR(4000000)
cc      dimension  nc(1000),palm(1000),palm1(1000,10)
cx      dimension  x(np),y(np), RR(np*np)
cx      dimension  nc(jmax),palm(jmax),palm1(jmax,m)
cx      dimension  amu(m),anu(m),aa(m),ss1(m),ss2(m)
cc      character*50 fname
cx      INTEGER NVAL
c     PARAMETER(X1=r0/2,X2=1.0,NVAL=10)
cx      INTEGER i
      integer nc(jmax), NVAL, i
      double precision RR(np*np), pi, t, r, x1, ss, tt, uu, Fr, eps,
     1                 Freps1, Freps2, dFr
c      EXTERNAL func
      EXTERNAL pafunc
cc      open(2,file='TypeAparam.palm')
cc      read(2,2) fname
cc    2 format(a)
cc      read(2,*) delta,ty
      ty=ty1
cc      read(2,*) x2
      tx=1.d0
cc      j=1
cc      open(1, FILE=fname)
cc 10   read(1, *, end=30) x(j), y(j) 
cc           j=j+1
cc           go to 10
cc 30   continue
cc      close(1)
cc      n=j-1
cc      np = n
c
cc      call  bdry(RR,NN,x,y)
      call  bdry(RR,NN,x,y,np)
c

      pi=3.14159265358979d0
      t = np
c
cc      do 15 k=1,1000
      do 15 k=1,jmax
          nc(k)=0
 15   continue  
c
      do 25 i=1, NN
cx          id=RR(i)/delta+1
          id=int(RR(i)/delta)+1
cc          nc(id)=nc(id)+1
          if(id.le.jmax)  nc(id)=nc(id)+1
 25   continue       
cc      do 35 i=1, 500
      do 35 i=1, jmax
         r = delta*i
         Palm(i)=((nc(i)/t)/((pi*((r+delta)**2))-(pi*((r)**2))))/t
 35   continue
cc      m=1
cc  120 continue
cc      read(2,*,end=110) amu,anu,a,s1,s2
      do 110 i=1,m
      a=aa(i)
      s1=ss1(i)
      s2=ss2(i)
c----------------------------------------------------------
cc      do 115 j=1, 500
      do 115 j=1, jmax
      r0=delta*j
      x1=r0/2
      nval=100
        kk=1
cc        call quad2d(X1,X2,ss)
        call quad2d(pafunc,X1,X2,ss)
        kk=2
cc        call quad2d(0.0d0,x1,tt)
        call quad2d(pafunc,0.0d0,x1,tt)
        kk=3
cc        call quad2d(0.0d0,x1,uu)
        call quad2d(pafunc,0.0d0,x1,uu)
      Fr=2*(ss+tt+uu)
c----------------------------------------------------------
      eps=0.001d0
      r0=delta*j+eps
      x1=r0/2
      nval=100
        kk=1
cc        call quad2d(X1,X2,ss)
        call quad2d(pafunc,X1,X2,ss)
        kk=2
cc        call quad2d(0.0d0,x1,tt)
        call quad2d(pafunc,0.0d0,x1,tt)
        kk=3
cc        call quad2d(0.0d0,x1,uu)
        call quad2d(pafunc,0.0d0,x1,uu)
      Freps1=2*(ss+tt+uu)
c----------------------------------------------------------
      r0=delta*j-eps
      if (r0.ne.0) then
      x1=r0/2
      nval=100
        kk=1
cc        call quad2d(X1,X2,ss)
        call quad2d(pafunc,X1,X2,ss)
        kk=2
cc        call quad2d(0.0d0,x1,tt)
        call quad2d(pafunc,0.0d0,x1,tt)
        kk=3
cc        call quad2d(0.0d0,x1,uu)
        call quad2d(pafunc,0.0d0,x1,uu)
      Freps2=2*(ss+tt+uu)
      end if
      if (r0.eq.0) then
      Freps2=0
      end if   
      dFr=(Freps1-Freps2)/(2*eps)
cc         palm1(j,m) = (amu * anu + anu*dFr/2/pi/(delta*j))/(amu*anu)
         palm1(j,i) = 
     &    (amu(i) * anu(i) + anu(i)*dFr/2/pi/(delta*j))/(amu(i)*anu(i))
115    continue
cc       m=m+1
cc       go to 120
110    continue
cc       m=m-1
cc       open(10, FILE='Palm-TypeA.txt')
cc       do 55 j=1,500
cc       r = delta*j
cc       write(10,1) r,palm(j),(palm1(j,i),i=1,m)
cc    1  format(f10.3,11e12.5)
cc 55    continue
cc      close(1)
cc      close(2)
cc      stop
      return
      END
c
cc      real*8 FUNCTION func(x,y)
cx      real*8 FUNCTION pafunc(x,y)
      DOUBLE PRECISION FUNCTION pafunc(x,y)
cx      implicit real*8(a-h,o-z)
      double precision x, y
      integer kk
      double precision r0, a, s1, s2
      common/distance/r0
      common/case/kk
      common/av/a,s1,s2
      double precision pi, qx, qy, xyr0
      pi=3.14159265358979d0
c
       qx = a/((s1)**2)*x*exp(-x**2/(2*(s1)**2))
     &     +(1-a)/((s2)**2)*x*exp(-x**2/(2*(s2)**2))
       qy = a/((s1)**2)*y*exp(-y**2/(2*(s1)**2))
     &     +(1-a)/((s2)**2)*y*exp(-y**2/(2*(s2)**2))
c
cc      if (kk.le.2) func=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
cx      if (kk.le.2) pafunc=(1/pi)*acos((x**2+y**2-(r0)**2)/(2*x*y))
cx     &                  *qx
cx     &                  *qy
      pafunc=0
      if (kk.le.2) then
         xyr0=(x**2+y**2-(r0)**2)/(2*x*y)
         if (abs(xyr0).le.1.0d0) then
            pafunc=(1/pi)*acos(xyr0)*qx*qy
         else
            pafunc=0
         end if
      end if
cc      if (kk.eq.3) func=1
      if (kk.eq.3) pafunc=1
     &                  *qx
     &                  *qy
      return
      END

