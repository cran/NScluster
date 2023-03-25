      include 'NScluster.h'
********************************************************************************
* common subroutine
c
c     qgausxMP   ---    simplexIP, simplexA
c     qgausyMP   ---    simplexIP, simplexA
c
********************************************************************************
cc      SUBROUTINE qgausx(func,a,b,ss)
      SUBROUTINE qgausxMP(func,a,b)
cx      implicit real*8(a-h,o-z)
      double precision, EXTERNAL :: func
      double precision a, b
      double precision r0
cxx      common/distance/r0
      common/distancep/r0

c      REAL a,b,ss,func
cx      EXTERNAL func
cx      INTEGER j
c      REAL dx,xm,xr,w(5),x(5)
cx      dimension w(5), x(5)
c
cxx      common/case/kk
cxx      COMMON /xyz/ xx,yy,zz
      integer j, kk
      double precision xx, yy, zz, xxm, xxr, xdx, xss, yy1, yy2, yy3,
     1                 hMP1, hMP2, yxm, yxr, ydx, yss, fMP1, fMP2
      common/casep/kk
      COMMON /xyzp/ xx,yy,zz
      common/param4/xxm,xxr,xdx,xss,yy1,yy2,yy3,hMP1,hMP2
      common/param5/yxm,yxr,ydx,yss,fMP1,fMP2
c
!$omp threadprivate(/distancep/)
!$omp threadprivate(/casep/)
!$omp threadprivate(/xyzp/)
!$omp threadprivate(/param4/)
!$omp threadprivate(/param5/)
c
      double precision w(5), x(5), ss
      SAVE w,x
      DATA w/.2955242247d0,.2692667193d0,.2190863625d0,.1494513491d0,
     *.0666713443d0/
      DATA x/.1488743389d0,.4333953941d0,.6794095682d0,.8650633666d0,
     *.9739065285d0/
c
cc      xm=0.5d0*(b+a)
cc      xr=0.5d0*(b-a)
      xxm=0.5d0*(b+a)
      xxr=0.5d0*(b-a)
      ss=0
      do 11 j=1,5
cc        dx=xr*x(j)
        xdx=xxr*x(j)
cc        ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
cc        ss=ss+w(j)*(hMP(func,xm+dx)+hMP(func,xm-dx))
        xx=xxm+xdx
        yy1=xx
        yy2=xx+r0
        yy3=-xx+r0
        if (kk.eq.1) call qgausyMP(func,yy1,yy2)
        if (kk.eq.2) call qgausyMP(func,yy3,yy2)
        if (kk.eq.3) call qgausyMP(func,yy1,yy3)
        hMP1=yss
c
        xx=xxm-xdx
        yy1=xx
        yy2=xx+r0
        yy3=-xx+r0
        if (kk.eq.1) call qgausyMP(func,yy1,yy2)
        if (kk.eq.2) call qgausyMP(func,yy3,yy2)
        if (kk.eq.3) call qgausyMP(func,yy1,yy3)
        hMP2=yss
        ss=ss+w(j)*(hMP1+hMP2)
11    continue
cc      ss=xr*ss
      xss=xxr*ss
      return
      END


cc      SUBROUTINE qgausy(func,a,b,ss)
      SUBROUTINE qgausyMP(func,a,b)
cx      implicit real*8(a-h,o-z)
      double precision, EXTERNAL :: func
      double precision a, b
      double precision r0
cxx      common/distance/r0
      common/distancep/r0
c      REAL a,b,ss,func
cx      EXTERNAL func
      INTEGER j
c      REAL dx,xm,xr,w(5),x(5)
cx      dimension w(5), x(5)
c
cxx      COMMON /xyz/ xx,yy,zz
      double precision xx, yy, zz, xxm, xxr, xdx, xss, yy1, yy2, yy3,
     1                 hMP1, hMP2, yxm, yxr, ydx, yss, fMP1, fMP2
      COMMON /xyzp/ xx,yy,zz
      common/param4/xxm,xxr,xdx,xss,yy1,yy2,yy3,hMP1,hMP2
      common/param5/yxm,yxr,ydx,yss,fMP1,fMP2
c
!$omp threadprivate(/distancep/)
!$omp threadprivate(/xyzp/)
!$omp threadprivate(/param4/)
!$omp threadprivate(/param5/)
c
      double precision w(5), x(5), ss
      SAVE w,x
      DATA w/.2955242247d0,.2692667193d0,.2190863625d0,.1494513491d0,
     *.0666713443d0/
      DATA x/.1488743389d0,.4333953941d0,.6794095682d0,.8650633666d0,
     *.9739065285d0/
c
      yxm=0.5d0*(b+a)
      yxr=0.5d0*(b-a)
      ss=0
      do 11 j=1,5
cc        dx=xr*x(j)
        ydx=yxr*x(j)
cc        ss=ss+w(j)*(f(xm+dx)+f(xm-dx))
cc        ss=ss+w(j)*(fMP(func,xm+dx)+fMP(func,xm-dx))
        yy=yxm+ydx
        fMP1=func(xx,yy)
        yy=yxm-ydx
        fMP2=func(xx,yy)
        ss=ss+w(j)*(fMP1+fMP2)
11    continue
cc      ss=xr*ss
      yss=yxr*ss
      return
      END
