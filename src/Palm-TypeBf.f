       subroutine palmBf(x,y,np,delta,ty1,amu,anu,a,s1,s2,m,jmax,
     &                   palm,palm1)
c
      include 'NScluster_f.h'
c
       implicit real*8(a-h,o-z)
cc       dimension  x(2000),y(2000), RR(4000000)
cc       dimension  nc(1000),palm(1000),palm1(1000,10)
cc       common/events/np
       dimension  x(np),y(np), RR(np*np)
       dimension  nc(jmax),palm(jmax),palm1(jmax,m)
       dimension  amu(m), anu(m), a(m), s1(m), s2(m)
       common / sizes / tx,ty
cc       character*50 fname
cc       open(2,file='TypeBparam.palm')
cc       read(2,2) fname
cc  2    format(a)
       tx=1.d0
cc       read(2,*) delta,ty
      ty=ty1
cc       j=1
cc       open(1, FILE=fname)
cc 10    read(1, *, end=30) x(j), y(j) 
cc           j=j+1
cc           go to 10
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
          if(id.le.jmax)  nc(id)=nc(id)+1
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
cc       read(2,*,end=50) amu,anu,a,s1,s2
       do 50 i = 1, m
cc       alam=amu*anu
       alam=amu(i)*anu(i)
cc       do 45 j = 1, 500
       do 45 j = 1, jmax
         r = delta*j
cc         ae1 = a * exp(-(r**2)/(4*((s1)**2)))/((s1)**2)
cc         ae2 = (1-a) * exp(-(r**2)/(4*((s2)**2)))/((s2)**2) 
cc         Palm1(j,m) = (alam + anu*(ae1 + ae2)/(4*pi))/alam
         ae1 = a(i) * exp(-(r**2)/(4*((s1(i))**2)))/((s1(i))**2)
         ae2 = (1-a(i)) * exp(-(r**2)/(4*((s2(i))**2)))/((s2(i))**2) 
         Palm1(j,i) = (alam + anu(i)*(ae1 + ae2)/(4*pi))/alam
 45    continue
cc       m=m+1
cc       go to 40
 50    continue
cc       m=m-1
cc       open(10, FILE='Palm-TypeB.txt')
cc       do 55 j=1,500
cc       r = delta*j
cc       write(10,1) r,palm(j),(palm1(j,i),i=1,m)
cc       write(*,*) r,palm(j),(palm1(j,i),i=1,m)
cc    1  format(f10.3,11e12.5)
cc 55    continue
c
cc       close(10)
cc       close(2)
cc       stop
       return
       end



