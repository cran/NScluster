      subroutine simAf(ix,iy,iz,ty,amu,anu,a,sig1,sig2,
     &                    npts,ncl,x,y,xcl,ycl,m,n,ier)
c
      include 'NScluster_f.h'
c
      implicit real*8 (a-h, o-z)
cc      common ix,iy,iz
cc      dimension  x(1000), y(1000)
cc      dimension  xcl(1000,5000), ycl(1000,5000)
      dimension  x(m), y(m)
      dimension  xcl(m,n), ycl(m,n)
      dimension  ncl(m)
c
      pi=3.14159265358979d0
c

cc      open(10, FILE='TypeA-offspring.xy')
cc      open(11, FILE='TypeA-parents.xy')
cc      open(12,file='TypeA.param')
cc      read(12,*) ix,iy,iz
cc      read(12,*) ty
cc      read(12,*) amu,anu,a,sig1,sig2
c
cc      call Pois(amu,npts)
      call Pois(amu,npts,ix,iy,iz)
c---
        np=0
        ier=0
        if( npts > m ) then
           ier=-1
           return
        end if
c---
c
      do 15 i=1, npts
        x(i)=random(ix,iy,iz)
        y(i)=random(ix,iy,iz)*ty               
cc        write(11,*) x(i), y(i)    
 15     continue
cc        write(6,*) '#(parents)=', npts
c
      do 25 i=1, npts
cc          call  Pois(anu,ncl)
          call  Pois(anu,ncl(i),ix,iy,iz)
cc          write(6,*) '#(offspring)=', i, ncl
c---
          if( ncl(i) > n ) then
             ier=-2
             return
          end if
c---
c
cc          do 35 j=1, ncl
          do 35 j=1, ncl(i)
 	    r=sqrt(-2*log(random(ix,iy,iz)))
	    theta=2*pi*(random(ix,iy,iz))	
            np=np+1
	    xclij=x(i)+r*cos(theta)*sig1     
	    yclij=y(i)+r*sin(theta)*sig1
	    xclij2=x(i)+r*cos(theta)*sig2     
	    yclij2=y(i)+r*sin(theta)*sig2
            jx=int(xclij)  
            jy=yclij/ty
            jx2=int(xclij2)  
            jy2=yclij2/ty
             if (xclij.le.0) then
               xclij=xclij+(1-jx)
             end if
             if (xclij2.le.0) then
               xclij2=xclij2+(1-jx2)
             end if
             if (yclij.le.0) then
               yclij=yclij+(1-jy)*ty
             end if
             if (yclij2.le.0) then
               yclij2=yclij2+(1-jy2)*ty
             end if
             if (xclij.ge.1) then
               xclij=xclij-jx
             end if 
             if (xclij2.ge.1) then
               xclij2=xclij2-jx2
             end if 
             if (yclij.ge.ty) then
               yclij=yclij-jy*ty
             end if
             if (yclij2.ge.ty) then
               yclij2=yclij2-jy2*ty
             end if
          choice=random(ix,iy,iz)               
          if(choice.le.a) then
            xcl(i,j)=xclij
            ycl(i,j)=yclij
          else
            xcl(i,j)=xclij2
            ycl(i,j)=yclij2
          endif
cc          write(10,*) xcl(i,j), ycl(i,j)
 35   continue
 11   format(' ',2F10.6)
 25   continue
cc        write(6,*) '#(total offspring)=', np
cc      close(10)
cc      close(11)
cc      close(12)
c
cc      stop
      return
c
      end
