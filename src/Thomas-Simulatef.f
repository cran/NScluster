      subroutine simThomf(ix,iy,iz,ty,amu,anu,sig,
     &                   npts,ncl,x,y,xcl,ycl,mmax,nmax,ier)
c
      include 'NScluster_f.h'
c
       implicit real*8 (a-h, o-z)
cc       common ix,iy,iz
cc       dimension  x(100), y(100)
cc       dimension  xcl(100,5000), ycl(100,5000)
       dimension  x(mmax), y(mmax)
       dimension  xcl(mmax,nmax), ycl(mmax,nmax)
       dimension  ncl(mmax)
c
       pi=3.14159265358979d0
c
cc       open(10, FILE='Thomas-offspring.xy')
cc       open(11, FILE='Thomas-parents.xy')
c
cc       open(12,file='Thomas.param')
cc       read(12,*) ix,iy,iz
cc       read(12,*) ty
cc       read(12,*) amu,anu,sig
c
        amu=amu*ty
cc        call Pois(amu,npts)
        call Pois(amu,npts,ix,iy,iz)
c---
        np=0
        ier=0
        if( npts > mmax ) then
           ier=-1
           return
        end if
c---
        do 15 i=1, npts
          x(i)=random(ix,iy,iz)
          y(i)=random(ix,iy,iz)*ty
cc          write(11,*) x(i), y(i)               
 15     continue  
cc        write(6,*) '#(parents)=', npts
c   
        do 25 i=1, npts
cc          call  Pois(anu,ncl)
          call  Pois(anu,ncl(i),ix,iy,iz)
cc          write(6,*) '#(offspring)=', i, ncl
c---
          if( ncl(i) > nmax ) then
             ier=-2
             return
          end if
c---
cc          do 35 j=1, ncl
          do 35 j=1, ncl(i)
 	    r=sqrt(-2*log(random(ix,iy,iz)))
	    theta=2*pi*(random(ix,iy,iz))	
            np=np+1
	    xclij=x(i)+r*cos(theta)*sig     
	    yclij=y(i)+r*sin(theta)*sig
            jx=INT(xclij)  
            jy=yclij/ty
               if (xclij.le.0) then
                 xclij=xclij+(1-jx)
               end if
               if (yclij.le.0) then
                 yclij=yclij+(1-jy)*ty
               end if
               if (xclij.ge.1) then
                 xclij=xclij-jx
               end if 
               if (yclij.ge.ty) then
                 yclij=yclij-jy*ty
               end if
c
              xcl(i,j)=xclij
              ycl(i,j)=yclij
cc           write(10,*) xcl(i,j), ycl(i,j)
 35       continue
 11	format(' ',2F10.6)
 25     continue
cc        write(6,*) '#(total offspring)=', np
c
cc        close(10)
cc        close(11)
cc        close(12)
c
cc        stop
        return
c
        end
