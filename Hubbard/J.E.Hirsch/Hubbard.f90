!03/07/2015  QH

module band
    real*8::t,E_max
    real*8,allocatable:: klist(:,:),qlist(:,:)
    real*8,allocatable:: Ek(:),Eom(:)
    complex*16,allocatable ::susc_qom(:),susc_om(:)
    real*8:: occup
    real*8:: pi=3.1415926,incr_k,incr_E,del
    integer*4:: Nk,Ne,d,beta,dim_k,dim_kom
contains
 
 subroutine input_params(bet,dimen,n_k,emax,n_e,dlt,occ,t1)
     implicit none
     integer*4  n_k,n_e,dimen,bet
     real*8 dlt,occ,t1,emax
     integer*4 i
     beta=bet
     d=dimen
     Nk=n_k
     E_max=emax
     Ne=n_e
     del=dlt
     occup=occ
     t=t1
     incr_k=2*pi/(Nk-1.)
     incr_E=E_max/(Ne-1.)
     dim_k=(Nk-1)**d
     dim_kom=dim_k*Ne
     allocate(susc_qom(dim_kom))
     allocate(Eom(Ne))
     allocate(susc_om(Ne))
     return
  end subroutine 
  
  subroutine kmesh()
    implicit none
    integer*4 i,j,ind
    real*8 kx,ky
    allocate(klist(dim_k,2))
    ind=1
      do j=1,Nk-1
        ky=-pi+(j-1.)*incr_k
        do i=1,Nk-1
          kx=-pi+(i-1.)*incr_k
          klist(ind,1)=kx
          klist(ind,2)=ky
          ind=ind+1
        enddo
      enddo
    return
  end subroutine
  
  subroutine qmesh()
     implicit none
     real*8 qx,qy
     integer i,j,ind
    allocate(qlist(dim_k,2))
     ind=1
     do j=1,Nk-1
       qy=-pi+(j-1.)*incr_k
       do i=1,Nk-1
         qx=-pi+(i-1.)*incr_k
         qlist(ind,1)=qx
         qlist(ind,2)=qy
         ind=ind+1
       enddo
      enddo
   end subroutine


  subroutine Emesh()
     implicit none
     integer i
     real*8 E
     do i=1,Ne
       Eom(i)=0.+incr_E*(i-1.)
     enddo
  end subroutine  
  
  subroutine ek_bz()
    implicit none
    integer*4 i,j 
    open(unit=10,FILE="Band_BZ.out")
    allocate(Ek(dim_k))
    j=1
    do i=1,dim_k
      Ek(i)=Ek_eval(klist(i,1),klist(i,2))
      if(mod(j,int((Nk-1)/10.))==0)then
        write(10,*)klist(i,1),klist(i,2),Ek(i)
      endif 
      if(mod(i,Nk-1)==0)then
         write(10,*)
      endif
    j=j+1
    enddo
    return
  end subroutine

  function Ek_eval(kx,ky)
    real*8 Ek_eval
    real*8 kx,ky
    Ek_eval=-2.*t*(cos(kx)+cos(ky)) 
  end function
  
  
  subroutine dos_mu(mu,n)
    implicit none
    complex*16 gloc
    complex*16 idel
    integer*4 i,j,k,n
    real*8 E,Emax,Emin,incr_E
    real*8 dos,mu
    character(len=1024)::filename
   ! open(unit=8,File="dos.out")
    write (filename,"(A7,I2)") "dos.txt",n 
    open(unit=8,File=filename)
    idel=dcmplx(0.,del)
    Emax=12.
    Emin=-12.
    incr_E=(Emax-Emin+1.)/NE
    E=Emin-0.5

    do while(E<=(Emax+0.5))
	dos=0.
	gloc=(0.,0.)
        do i=1,dim_k
		gloc=gloc+( E+mu-Ek(i)+idel )**-1
        enddo
        E=E+incr_E
	dos=-1/pi*imag(gloc)/dim_k
	write(8,*)E,dos
    enddo
    close(8)
    return
  end subroutine

  subroutine mu_fix(murng,mu) 
    implicit none
    real*8 murng(2)
    real*8,intent(out)::mu
    real*8 num(2),num_mu
    integer*4 i
   ! write(*,*)"Before convergence"
    do i=1,2
    call num_count(murng(i),num(i))  !In module,inter-subroutine call is checked!
   ! write(*,*)d,murng(i),num(i)!,num_p,num_d
    enddo
    if(num(1)<occup.and.num(2)>occup)then
    else
    write(*,*)"wrong inital guess for the lower bound and upper bound for the chemical potential"
    call exit(0)
    endif
    call bisect_mu(murng,num)
    mu=(murng(1)+murng(2))/2
    call num_count(mu,num_mu)
    write(*,*)mu,num_mu
    return
  end subroutine

  subroutine bisect_mu(mu,num)
    implicit none
    real*8 mu(2),num(2)
    real*8 mu_mid,num_mid
    integer*4 it,N_it
    real*8 eps
    eps=0.000002
    it=1
    N_it=100000
    do while(abs(mu(2)-mu(1))>eps.and.it<N_it)
    mu_mid=(mu(1)+mu(2))/2
    call num_count(mu_mid,num_mid)
        if(num_mid>occup)then
          mu(2)=mu_mid
          num(2)=num_mid
        else
          mu(1)=mu_mid
          num(1)=num_mid
        endif
    it=it+1
    enddo
    return
  end subroutine

  subroutine num_count(mu,num)
    implicit none
    real*8 mu,num
    integer*4 i
    num=0.
    do i=1,dim_k
         num=num+fermi(Ek(i),mu)
    enddo
    num=num*2./dim_k
    return
  end subroutine
   
  subroutine band_struc(mu) 
     implicit none
     real*8 mu
     real*8 kx,ky,Ekpath
     integer i,N
     open(unit=195,FILE="band_kpath.out")
     N=1000
     do i=1,N
       kx=pi*(i-1.)/N
       ky=0.
       Ekpath=Ek_eval(kx,ky)
       write(195,*)i,Ekpath
     enddo 
     do i=1,N
       kx=pi
       ky=pi*(i-1.)/N
       Ekpath=Ek_eval(kx,ky)
       write(195,*)i+N,Ekpath
     enddo 
     do i=1,N
       kx=pi-pi*(i-1.)/N
       ky=pi-pi*(i-1.)/N
       Ekpath=Ek_eval(kx,ky)
       write(195,*)i+2*N,Ekpath
     enddo 
  end subroutine
  
  subroutine fermi_surface(mu,m)
     implicit none
     real*8 mu
     integer m
     real*8 kx,ky,dk,Ekpath
     real*8 eps
     integer i,j,N
     character(len=1024):: filename
     write(filename,"(A6,I2)")"FS.out",m
     !open(unit=9,FILE="Fs.out")
     open(unit=9,FILE=filename)
     N=1000
     eps=0.02        !this parameter and N need to be tuned.
     dk=2.*pi/(N-1.)
     do i=1,N
       kx=-pi+dk*(i-1.)
       do j=1,N
         ky=-pi+dk*(j-1.)
         Ekpath=Ek_eval(kx,ky)
         if(abs(Ekpath-mu)<eps)then
         write(9,*)kx,ky
         else
         endif
         if(j==N)then
         write(9,*)
         endif
       enddo
     enddo
     close(9)
  end subroutine
  
  function fermi(E,mu)
     implicit none
     real*8 fermi
     real*8 E,mu
     fermi= 1./( exp( beta*(E-mu) )+1.)  
  end function

   function susc_kqom(kx,ky,qx,qy,omega,mu)
      implicit none
      complex*16 susc_kqom
      real*8 kx,ky,qx,qy,omega,mu
      real*8 kmqx,kmqy
      real*8 Ekval,Ekqval
      complex*16 idel
      idel=dcmplx(0,del)
      kmqx=kx-qx
      kmqy=ky-qy
      Ekval=Ek_eval(kx,ky)
      Ekqval=Ek_eval(kmqx,kmqy)
      susc_kqom=-( fermi(Ekqval,mu)-fermi(Ekval,mu) )/(omega-(Ekval-Ekqval)+idel)
   end function
   
   subroutine suscpt_screen(mu)
      implicit none
      real*8 mu
      real*8 kx,ky,qx,qy,omega
      integer i,j,k,ind
      open(unit=20,File="susc_sta_q.out")  
      open(unit=21,File="susc_loc_om.out")  
      ind=1
      ! omega loop
      do k=1,Ne
         omega=Eom(k)
         susc_om(k)=0.
         !! q loop.
         do j=1,dim_k
            qx=qlist(j,1)
            qy=qlist(j,2)
            susc_qom(ind)=0.
            !!! k sum looop
            do i=1,dim_k
               susc_qom(ind)=susc_qom(ind)+susc_kqom(kx,ky,qx,qy,omega,mu)
            enddo
            !!!
            susc_qom(ind)=susc_qom(ind)/dim_k
            if(omega==0.)then
               write(20,*)qx,qy,real(susc_qom(ind)),imag(susc_qom(ind))
               if (mod(j,Nk-1)==0)then
                   write(20,*)
               endif
            endif  
            susc_om(k)=susc_om(k)+susc_qom(ind)
            ind=ind+1
         enddo
         !!
         susc_om(k)=susc_om(k)/dim_k
         write(21,*)omega,real(susc_om(k)),imag(susc_om(k)) 
      enddo
      !
      return
   end subroutine


 
end module
