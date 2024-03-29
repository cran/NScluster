ccx      subroutine simCf(ix,iy,iz,ty,amu1,amu2,anu1,anu2,s1,s2,
      subroutine simC(ix,ty,amu1,amu2,anu1,anu2,s1,s2,m1,ncl1,x1,y1,
     & xx1,yy1,m2,ncl2,x2,y2,xx2,yy2,mmax,nmax,ier)
c
      include 'NScluster.h'
c
cx      implicit real*8 (a-h, o-z)
cc      common ix,iy,iz
cc      dimension  x1(3000), y1(3000)
cc      dimension  x2(3000), y2(3000)
cc      dimension  xcl1(3000,5000), ycl1(3000,5000)
cc      dimension  xcl2(3000,5000), ycl2(3000,5000)
cc      dimension xx1(15000),yy1(15000)
cc      dimension xx2(15000),yy2(15000) 
cx      dimension  x1(mmax), y1(mmax)
cx      dimension  x2(mmax), y2(mmax)
cx      dimension  xcl1(mmax,nmax), ycl1(mmax,nmax)
cx      dimension  xcl2(mmax,nmax), ycl2(mmax,nmax)
cx      dimension xx1(mmax*nmax), yy1(mmax*nmax)
cx      dimension xx2(mmax*nmax), yy2(mmax*nmax)
cx      dimension ncl1(mmax), ncl2(mmax)
ccx      integer ix, iy, iz, m1, m2, mmax, nmax, ncl1(mmax), ncl2(mmax),
ccx     1           ier
      integer ix, m1, m2, mmax, nmax, ncl1(mmax), ncl2(mmax), ier
      double precision ty, amu1, amu2, anu1, anu2, s1, s2, x1(mmax),
     1                 y1(mmax), xx1(mmax*nmax), yy1(mmax*nmax),
     2                 x2(mmax), y2(mmax), xx2(mmax*nmax),
     3                 yy2(mmax*nmax)
      double precision xcl1(mmax,nmax), ycl1(mmax,nmax),
     1                 xcl2(mmax,nmax), ycl2(mmax,nmax), pi, r1, r2,
     2                 theta1, theta2, random, xcl1ij, ycl1ij, xcl2ij,
     2                 ycl2ij
c
      pi=3.14159265358979d0
c
cc      open(10, FILE='TypeC-offspring.xy')
cc      open(11,FILE='TypeC-parents.xy')
cc      open(12,file='TypeC.param')
cc      read(12,*) ix,iy,iz
cc      read(12,*) ty
cc      read(12,*) amu1,amu2,anu1,anu2,s1,s2
c
cc      write(6,*) 'a =', amu1*anu1/(amu1*anu1+amu2*anu2)
c
cc      call  Pois(amu1,m1)
ccx      call  Pois(amu1,m1,ix,iy,iz)
      call init(ix)
      call  Pois(amu1,m1)
c---
        ier=0
        if( m1 > mmax ) then
           ier=-1
           return
        end if
c---
c
      do 25 i=1, m1
ccx        x1(i)=random(ix,iy,iz)
ccx        y1(i)=random(ix,iy,iz)*ty
        x1(i)=random()
        y1(i)=random()*ty
cc        write(11,*) x1(i), y1(i)
 25   continue
cc        write(6,*) '#(parents)=', m1
c
      np=0
      np1=0
      do 40 i=1, m1
cc          call  Pois(anu1,ncl1)
ccx          call  Pois(anu1,ncl1(i),ix,iy,iz)
          call  Pois(anu1,ncl1(i))
cc          write(6,*) '#(offspring)=', i, ncl1
c---
          if( ncl1(i) > nmax ) then
             ier=-11
             return
          end if
c---
cc        do 45 j=1, ncl1
        do 45 j=1, ncl1(i)
ccx          r1=sqrt(-2*log(random(ix,iy,iz)))
ccx          theta1=2*pi*(random(ix,iy,iz))
          r1=sqrt(-2*log(random()))
          theta1=2*pi*(random())
          np1=np1+1
          np=np+1
          xcl1(i,j)=x1(i)+r1*cos(theta1)*s1
          ycl1(i,j)=y1(i)+r1*sin(theta1)*s1
          jx1=INT(xcl1(i,j))  
cx          jy1=ycl1(i,j)/ty
          jy1=int(ycl1(i,j)/ty)
c
          if (xcl1(i,j).le.0) then
             xcl1(i,j)=xcl1(i,j)+(1-jx1)
          end if
c
          if (ycl1(i,j).le.0) then
             ycl1(i,j)=ycl1(i,j)+(1-jy1)*ty
          end if
c
          if (xcl1(i,j).ge.1) then
             xcl1(i,j)=xcl1(i,j)-jx1
          end if 
c
          if (ycl1(i,j).ge.ty) then
             ycl1(i,j)=ycl1(i,j)-jy1*ty
          end if
c
           xcl1ij = xcl1(i,j)
           ycl1ij = ycl1(i,j)
           xx1(np)=xcl1ij
           yy1(np)=ycl1ij
cc          write(10,*) xx1(np), yy1(np)
 45       continue
 40     continue
c----------------------------------------------------------------------
cc      call  Pois(amu2,m2)
ccx      call  Pois(amu2,m2,ix,iy,iz)
      call  Pois(amu2,m2)
c---
        ier=0
        if( m2 > mmax ) then
           ier=-2
           return
        end if
c---
        do 35 i=1, m2
ccx          x2(i)=random(ix,iy,iz)
ccx          y2(i)=random(ix,iy,iz)*ty
          x2(i)=random()
          y2(i)=random()*ty
cc          write(11,*) x2(i), y2(i)
 35     continue
cc        write(6,*) '#(parents)=', m2
c
        np2=0
        do 50 i=1, m2
cc          call  Pois(anu2,ncl2)
ccx          call  Pois(anu2,ncl2(i),ix,iy,iz)
          call  Pois(anu2,ncl2(i))
c---
          if( ncl2(i) > nmax ) then
             ier=-22
             return
          end if
c---
cc          write(6,*) '#(offspring)=', i, ncl2
cc          do 55 j=1, ncl2
         do 55 j=1, ncl2(i)
ccx            r2=sqrt(-2*log(random(ix,iy,iz)))
ccx            theta2=2*pi*(random(ix,iy,iz))
            r2=sqrt(-2*log(random()))
            theta2=2*pi*(random())
            np2=np2+1
            np=np+1
            xcl2(i,j)=x2(i)+r2*cos(theta2)*s2
            ycl2(i,j)=y2(i)+r2*sin(theta2)*s2
            jx2=int(xcl2(i,j))
cx            jy2=ycl2(i,j)/ty
            jy2=int(ycl2(i,j)/ty)
               if (xcl2(i,j).le.0) then
               xcl2(i,j)=xcl2(i,j)+(1-jx2)
               end if
               if (ycl2(i,j).le.0) then
               ycl2(i,j)=ycl2(i,j)+(1-jy2)*ty
               end if
               if (xcl2(i,j).ge.1) then
               xcl2(i,j)=xcl2(i,j)-jx2
               end if 
               if (ycl2(i,j).ge.ty) then
               ycl2(i,j)=ycl2(i,j)-jy2*ty
               end if
           xcl2ij = xcl2(i,j)
           ycl2ij = ycl2(i,j)
cc           xx2(np)=xcl2ij
cc           yy2(np)=ycl2ij
           xx2(np2)=xcl2ij
           yy2(np2)=ycl2ij
cc          write(10,*) xx2(np), yy2(np)
 55       continue
 50     continue
cc        write(6,*) '#(total offspring)=', np1,np2,np
cc        close(10)
cc        close(11)
cc        close(12)
c
cc      stop
      return
c
      end
