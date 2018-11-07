module grapheneband 
  implicit none
  real*8 :: t=1.
  integer*4::d=2
  real*8,allocatable ::Eupk(:,:),Ednk(:,:)
  real*8::pi=3.1415926
  integer*4 :: Nk=5000,Ne=500
  real*8::delta=1.0e-2
end module


program threed_dos
use grapheneband
implicit none
real*8 Emin,Emax
real*8 E,incr_E
real*8 gE,dos
complex*16 gloc,idel
integer*4 i,j
open(unit=7,file="dos.txt")
open(unit=8,file="gloc.txt")
call kmesh()

Emin=-3*t
Emax=3*t
incr_E=(Emax-Emin+2)/NE
idel=dcmplx(0.,delta)
E=Emin-1

do while (E<=(Emax+1) )
  dos=0.
  gE=0.
  gloc=(0.0,0.0)
     do j=1,Nk
	  do i=1,Nk
		     !##########1.using Gloc to compute  dos######!
		     gloc=gloc+(  E-Eupk(i,j)+idel  )**-1+( E-Ednk(i,j)+idel )**-1
		     !##########2.using Direct counting method to compute dos####!
		     if(Eupk(i,j)>=E.and.Eupk(i,j)<E+incr_E)then
			gE=gE+1
		     endif
		     if(Ednk(i,j)>=E.and.Ednk(i,j)<E+incr_E)then
			gE=gE+1
		     endif
	  enddo
     enddo
  dos=gE/(incr_E*Nk**d)
  gloc=gloc/Nk**d
  write(7,*)E,dos
  write(8,*)E,-1/pi*imag(gloc)
  E=E+incr_E
enddo
end program 


subroutine kmesh()
    use grapheneband
    implicit none
    integer*4::i,j,m,n
    real*8::k(d),G1(d),G2(d),dG1(d),dG2(d)
    allocate(Eupk(Nk,Nk))
    allocate(Ednk(Nk,Nk))

    G1(1)=4.*pi/3.*sqrt(3.)/2.  !first unit vector G1:G1(1) in kx direction;G1(2) in ky direction
    G1(2)=4.*pi/3.*0.5
    G2(1)=0                    !second unit vector G2
    G2(2)=4.*pi/3.
    do m=1,d
    dG1(m)=G1(m)/(Nk-1.)       !Determine k increment along G1 and G2(the conventional BZ) 
    dG2(m)=G2(m)/(Nk-1.)
    enddo

    do j=1,Nk 
      do i=1,Nk
	  do m=1,d
	    k(m)=(i-1.)*dG1(m)+(j-1.)*dG2(m)    !Generating kmesh along G1 and G2(the conventional BZ) 
	  enddo
	  Eupk(i,j)=t*sqrt( 3.+4.*cos( sqrt(3.)/2.*k(1) )*cos( 3./2.*k(2) )+2.*cos( sqrt(3.)*k(1) )  )
          Ednk(i,j)=-Eupk(i,j)
      enddo
    enddo
    return
end subroutine
