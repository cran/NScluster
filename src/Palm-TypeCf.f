       subroutine palmCf(x,y,np,delta,ty1,alam,anu1,a,s1,s2,m,jmax,
     & palm,palm1)
c
      include 'NScluster_f.h'
c
       implicit real*8(a-h,o-z)
cc       dimension  x(2000),y(2000), RR(4000000)
cc       dimension  nc(1000),palm(1000),palm1(1000,10)
       dimension  x(np),y(2000), RR(np*np)
       dimension  nc(jmax),palm(jmax),palm1(jmax,m)
       dimension  alam(m),anu1(m),a(m),s1(m),s2(m)
cc       common/events/np
       common / sizes / tx,ty
cc       character*50 fname
cc       open(2,file='TypeCparam.palm')
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
cc       read(2,*,end=50) alam,anu1,a,s1,s2
       do 50 i = 1,m
cc       anu2=anu1*(s2/s1)
cc       do 45 j = 1, 500
       anu2=anu1(i)*(s2(i)/s1(i))
       do 45 j = 1, jmax
         r = delta*j
cc         e1j=(a * anu1/(s1**2))*exp(-(r**2)/(4*(s1**2)))
cc         e2j=((1-a) * anu2/(s2**2))*exp(-(r**2)/(4*(s2**2)))
cc         Palm1(j,m) = (alam + (e1j+e2j)/4/pi)/alam
         e1j=(a(i) * anu1(i)/(s1(i)**2))*exp(-(r**2)/(4*(s1(i)**2)))
         e2j=((1-a(i)) * anu2/(s2(i)**2))*exp(-(r**2)/(4*(s2(i)**2)))
         Palm1(j,i) = (alam(i) + (e1j+e2j)/4/pi)/alam(i)
 45    continue
cc       m=m+1
cc       go to 40
 50    continue
cc       m=m-1
cc       open(10, FILE='Palm-TypeC.txt')
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
