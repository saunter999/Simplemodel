!=================================!
!Calculating the band structure of graphene along the k path Gamma-K'-M-k-Gamma in the first BZ.(paramters:a=1,t=1,see related pdf)           2014.12.14
!=================================!
program band_struc
implicit none
real Gam(2),K_prime(2),M(2),K(2),dk(2),kp(2)
real E,E_ub,E_db,pi
parameter(pi=3.141592)
integer i,j,N
!==================!
!N=number of k points sampled along legs
!Gam(2),K_prime(2),M(2),K(2) special k points in 2d WS-BZ
!dk(2) increment of k points along legs
!kp(2) vectors to store the coordinates of k points along legs
!==================!
open(unit=7,file="upperband.txt")
open(unit=8,file="downband.txt")
N=100
Gam(1)=0
Gam(2)=0
K_prime(1)=4*pi/3*(1/(2*3**0.5))
K_prime(2)=4*pi/3*0.5
K(1)=-K_prime(1)
K(2)=K_prime(2)
! we don't need to define the M point since it is the middle point of K'-K.
!============== along the Gamma-K' leg===================!
do i=1,2
dk(i)=(K_prime(i)-Gam(i))/real(N)
kp(i)=Gam(i)
!write(*,*)dk(i),kp(i)
enddo
do i=1,N
E=(3+4*cos(0.5*3**0.5*kp(1))*cos(1.5*kp(2))+2*cos(3**0.5*kp(1)))**0.5
write(7,*)E
write(8,*)-E
	do j=1,2
	kp(j)=kp(j)+dk(j)
	enddo
!if(i==N)then
!write(*,*)kp(1)-K_prime(1),kp(2)-k_prime(2)
!endif
enddo
!================along the K'-K leg=====================!
do i=1,2
dk(i)=(K(i)-K_prime(i))/real(N)
kp(i)=K_prime(i)
enddo
do i=1,N
E=(3+4*cos(0.5*3**0.5*kp(1))*cos(1.5*kp(2))+2*cos(3**0.5*kp(1)))**0.5
write(7,*)E
write(8,*)-E
        do j=1,2
        kp(j)=kp(j)+dk(j)
        enddo
enddo
!==============along the K-Gamma leg====================!
do i=1,2
dk(i)=(Gam(i)-K(i))/real(N)
kp(i)=K(i)
enddo
do i=1,N
E=(3+4*cos(0.5*3**0.5*kp(1))*cos(1.5*kp(2))+2*cos(3**0.5*kp(1)))**0.5
write(7,*)E
write(8,*)-E
        do j=1,2
        kp(j)=kp(j)+dk(j)
        enddo
enddo

end program band_struc

