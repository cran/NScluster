      subroutine palmTf(x,y,np,delta,ty1,amu,anu,v,m,jmax,palm,palm1)
c
      include 'NScluster_f.h'
c
       implicit real*8(a-h,o-z)
cc       dimension  x(2000),y(2000), RR(4000000)
cc       dimension  nc(1000),palm(1000),palm1(1000,10)
cc       common/events/np
       common / sizes / tx,ty
cc       character*50 fname
       dimension  x(np),y(np), RR(np*np)
       dimension  amu(m),anu(m),v(m)
       dimension  nc(jmax),palm(jmax),palm1(jmax,m)
cc       open(2,file='Thomasparam.palm')
cc       read(2,2) fname
cc  2    format(a)
       tx=1.d0
cc       read(2,*) delta,ty
       ty=ty1
cc       j=1
cc       open(1, FILE=fname)
cc 10    read(1, *, end=30) x(j), y(j) 
cc             j=j+1
cc             go to 10
cc 30    continue
cc       close(1)
cc       n=j-1
cc       np = n
c
cc       call  bdry(RR,NN,x,y)
       call  bdry(RR,NN,x,y,np)
c
       pi=3.14159265358979d0
       t = np
c
cc       do 15 k=1,1000
       do 15 k=1,jmax
          nc(k)=0
 15    continue  
c
       do 25 i=1, NN
          id=RR(i)/delta+1
cc          nc(id)=nc(id)+1
          if(id.le.jmax) nc(id)=nc(id)+1
 25    continue       
       pi=3.14159265358979d0       
       delta=0.001d0
c
cc       do 35 i=1, 500
       do 35 i=1, jmax
          r = delta*i
          Palm(i)=((nc(i)/t)/((pi*((r+delta)**2))-(pi*((r)**2))))/t
 35    continue
c
cc       m=1
cc 40    continue
cc       read(2,*,end=50) amu,anu,v
cc       do 45 j = 1, 500
       do 50 i=1,m
       do 45 j = 1, jmax
         r = delta*j
cc         ae = exp(-(r**2)/(4*(v**2)))
cc         Palm1(j,m) = (amu*anu + anu*ae/(4*pi*(v**2)))/(amu*anu)
         ae = exp(-(r**2)/(4*(v(i)**2)))
         Palm1(j,i) =
     &    (amu(i)*anu(i) + anu(i)*ae/(4*pi*(v(i)**2)))/(amu(i)*anu(i))
 45    continue
cc       m=m+1
cc       go to 40
 50    continue
cc       m=m-1
cc       open(10, FILE='Palm-Thomas.txt')
cc       do 55 j=1,500
cc       r = delta*j
cc       write(10,1) r,palm(j),(palm1(j,i),i=1,m)
cc    1  format(f10.3,11e12.5)
cc 55    continue
c
cc       close(10)
cc       close(2)
cc       stop
       return
       end
*************************
cc        subroutine bdry(RR,NN,x,y)
        subroutine bdry(RR,NN,x,y,np)
c
        implicit real*8(a-h,o-z)
cc        common/events/np
        common / sizes / tx,ty
cc        dimension  x(2000), y(2000), RR(4000000)
        dimension  x(np), y(np), RR(np*np)
        t1=0.5d0
c periodic boundary
        N = np
        NN=0
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
        NN=NN+1
cc        IF(NN.GT.4000000) STOP 229
        RR(NN)=R2
 20     CONTINUE
 10     CONTINUE
c
        return
        end        
