!03/01/2015  QH

module bubble
    real*8:: E_p,E_d,t_pp,t_dd,t_pd,E_max
    real*8,allocatable:: klist(:,:),qlist(:,:)
    real*8,allocatable:: Epk(:),Edk(:),tpdk(:),Eupk(:),Ednk(:),Eom(:)
    complex*16,allocatable :: bub_dp_qom(:),bub_dp_om(:)
    complex*16,allocatable :: bub_pp_qom(:),bub_pp_om(:)
    complex*16,allocatable :: bub_dd_qom(:),bub_dd_om(:)
    complex*16,allocatable :: u11_scr(:),u12_scr(:),u21_scr(:),u22_scr(:)
    complex*16 :: Prest(4,4),uscr(4,4),uscr_loc(4,4)
    complex*16 :: Pdd(2,2),wscr(2,2),udd_mat(2,2)
    real*8 :: U(4,4),iden(4,4)
    real*8:: U_pp,U_dd,U_pd,mu
    real*8:: pi=3.1415926,incr_k,incr_E,del
    integer*4:: Nk,Ne,d,beta
    integer*4 ::dim_k,dim_kom
contains
 
 subroutine input_params(bet,dimen,n_k,emax,n_e,dlt,Udd,Upp,Upd,Umatx,Pmatx,chmpot,Ep,Ed,tpp,tdd,tpd)
     implicit none
     integer*4 :: n_k,emax,n_e,dimen,bet
     real*8 :: dlt,Udd,Upp,Upd,chmpot,Ep,Ed,tpp,tdd,tpd
     real*8 :: Umatx(4,4),Pmatx(4,4)
     integer i,j
     beta=bet
     d=dimen
     Nk=n_k
     E_max=emax
     Ne=n_e
     del=dlt
     U_dd=Udd
     U_pp=Upp
     U_pd=Upd
     U=Umatx
     Prest=Pmatx
     mu=chmpot
     E_p=Ep
     E_d=Ed
     t_pp=tpp
     t_dd=tdd
     t_pd=tpd
     incr_k=2*pi/(Nk-1.)
     incr_E=E_max/(Ne-1.)
     dim_k=(Nk-1)**d
     dim_kom=((Nk-1)**d)*Ne
     allocate(klist(dim_k,2)) !The first entry stores the index of the kmesh,while the second entry stores (kx,ky)
     allocate(qlist(dim_k,2))
     allocate(Eom(Ne))
     allocate(bub_dp_qom(dim_kom)) 
     allocate(bub_dp_om(Ne))
     allocate(bub_pp_qom(dim_kom)) 
     allocate(bub_pp_om(Ne))
     allocate(bub_dd_qom(dim_kom)) 
     allocate(bub_dd_om(Ne))
     allocate(u11_scr(Ne),u12_scr(Ne),u21_scr(Ne),u22_scr(Ne))
     do i=1,4
        do j=1,4
        if(i==j)then
          iden(i,j)=1.
        else
          iden(i,j)=0.
        endif
        enddo
     enddo
    ! print *,U(:,:)
    ! print *,iden(:,:)
     return
  end subroutine
 
  subroutine kmesh()
     implicit none
     real*8 kx,ky
     integer i,j,ind
     ind=1
     do j=1,Nk-1                !Nk is equivalent to 1
       ky=-pi+(j-1.)*incr_k
       do i=1,Nk-1
         kx=-pi+(i-1.)*incr_k
         klist(ind,1)=kx
         klist(ind,2)=ky
         ind=ind+1
       enddo
     enddo
  end subroutine
  
  subroutine qmesh()
     implicit none
     real*8 qx,qy
     integer i,j,ind
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
  
  subroutine eps_k()
     implicit none
     integer i
     open(unit=8,FILE="kmesh.out")
     allocate(Epk(dim_k)) 
     allocate(Edk(dim_k)) 
     allocate(tpdk(dim_k)) 
     allocate(Eupk(dim_k)) 
     allocate(Ednk(dim_k))
     do i=1,dim_k
       Epk(i)=Epk_val(klist(i,1),klist(i,2))
       Edk(i)=Edk_val(klist(i,1),klist(i,2))
       tpdk(i)=tpdk_val(klist(i,1),klist(i,2))
       Eupk(i)=Eupk_val(klist(i,1),klist(i,2))!( Edk(i)+Epk(i)+sqrt( (Edk(i)-Epk(i))**2+4*tpdk(i)**2 ) )/2.    !upband dispersion  
       Ednk(i)=Ednk_val(klist(i,1),klist(i,2))!( Edk(i)+Epk(i)-sqrt( (Edk(i)-Epk(i))**2+4*tpdk(i)**2 ) )/2.  !lower band dispersion  
     enddo 

      do i=1,dim_k
       write(8,*)klist(i,1),klist(i,2),Epk(i)
      enddo
  end subroutine
  
  function Epk_val(kx,ky)
    implicit none
    real*8 Epk_val
    real*8 kx,ky
    Epk_val=E_p-2*t_pp*( cos(kx)+cos(ky) )
  end function

  function Edk_val(kx,ky)
    implicit none
    real*8 Edk_val
    real*8 kx,ky
    Edk_val=E_d-2*t_dd*( cos(kx)+cos(ky) )
  end function
  
  function tpdk_val(kx,ky)
    implicit none
    real*8 tpdk_val
    real*8 kx,ky
    tpdk_val=-t_pd*( sin(kx/2.)+sin(ky/2.) )
  end function
  
  function Eupk_val(kx,ky)
    implicit none
    real*8 Eupk_val
    real*8 kx,ky
    Eupk_val=( Edk_val(kx,ky)+Epk_val(kx,ky)+sqrt( (Edk_val(kx,ky)-Epk_val(kx,ky))**2+4*tpdk_val(kx,ky)**2 ) )/2.    !upband dispersion  
  end function

  function Ednk_val(kx,ky)
    implicit none
    real*8 Ednk_val
    real*8 kx,ky
    Ednk_val=( Edk_val(kx,ky)+Epk_val(kx,ky)-sqrt( (Edk_val(kx,ky)-Epk_val(kx,ky))**2+4*tpdk_val(kx,ky)**2 ) )/2.    !dnband dispersion  
  end function

  function fermi(E)
     implicit none
     real*8 fermi
     real*8 E
     fermi= 1./( exp( beta*(E-mu) )+1.) 
  end function
   
  function gf_dp_kom_rtd(kx,ky,omega,sig)
     implicit none
     complex*16 gf_dp_kom_rtd
     real*8 kx,ky
     real*8 omega,sig
     complex*16 om
     om=dcmplx(omega,sig*del)
     gf_dp_kom_rtd=tpdk_val(kx,ky)/(  (om-Edk_val(kx,ky))*(om-Epk_val(kx,ky))-tpdk_val(kx,ky)**2     )    
  end function

  function bubble_dp_kqom(kx,ky,qx,qy,omega)
     implicit none
     complex*16 bubble_dp_kqom
     real*8 kx,ky,qx,qy,omega 
     real*8 kmqx,kmqy
     real*8 tpdkval,Eupkval,Ednkval,tpdkq,Eupkq,Ednkq
     real*8 om1,om2,om3,om4,sig1,sig2
     kmqx=kx-qx
     kmqy=ky-qy
     tpdkval=tpdk_val(kx,ky)
     Eupkval=Eupk_val(kx,ky)
     Ednkval=Ednk_val(kx,ky)
     tpdkq=tpdk_val(kmqx,kmqy)
     Eupkq=Eupk_val(kmqx,kmqy)
     Ednkq=Ednk_val(kmqx,kmqy)
     om1=Eupkval-omega
     om2=Ednkval-omega
     sig1=-1.0
     sig2=1.0 
     om3=Eupkq+omega
     om4=Ednkq+omega
     bubble_dp_kqom=tpdkval/(Eupkval-Ednkval)* &
      ( fermi(Eupkval)*gf_dp_kom_rtd(kmqx,kmqy,om1,sig1)-fermi(Ednkval)*gf_dp_kom_rtd(kmqx,kmqy,om2,sig1) ) &
                    +tpdkq/(Eupkq-Ednkq)*&
      ( fermi(Eupkq)*gf_dp_kom_rtd(kx,ky,om3,sig2)-fermi(Ednkq)*gf_dp_kom_rtd(kx,ky,om4,sig2) ) 
  end function

 
  function gf_pp_kom(kx,ky,omega,sig)
     implicit none
     complex*16 gf_pp_kom
     real*8 kx,ky
     real*8 omega,sig
     complex*16 om
     om=dcmplx(omega,sig*del)
     gf_pp_kom=( om-Edk_val(kx,ky) )/(  (om-Edk_val(kx,ky))*(om-Epk_val(kx,ky))-tpdk_val(kx,ky)**2     )    
  end function

  function bubble_pp_kqom(kx,ky,qx,qy,omega)
     implicit none
     complex*16 bubble_pp_kqom
     real*8 kx,ky,qx,qy,omega 
     real*8 kmqx,kmqy
     real*8 Eupkval,Ednkval,Edkval,Eupmdk,Ednmdk
     real*8 Eupkq,Ednkq,Edkq,Eupmdkq,Ednmdkq
     real*8 om1,om2,om3,om4,sig1,sig2
     kmqx=kx-qx
     kmqy=ky-qy
     Eupkval=Eupk_val(kx,ky)
     Ednkval=Ednk_val(kx,ky)
     Edkval=Edk_val(kx,ky)
     Eupmdk=Eupkval-Edkval
     Ednmdk=Ednkval-Edkval
     Eupkq=Eupk_val(kmqx,kmqy)
     Ednkq=Ednk_val(kmqx,kmqy)
     Edkq=Edk_val(kmqx,kmqy)
     Eupmdkq=Eupkq-Edkq
     Ednmdkq=Ednkq-Edkq
     om1=Eupkval-omega
     om2=Ednkval-omega
     om3=Eupkq+omega
     om4=Ednkq+omega
     sig1=-1.0
     sig2=1.0 
     bubble_pp_kqom=1./(Eupkval-Ednkval)* &
   (fermi(Eupkval)*Eupmdk*gf_pp_kom(kmqx,kmqy,om1,sig1)-fermi(Ednkval)*Ednmdk*gf_pp_kom(kmqx,kmqy,om2,sig1)) &
                    +1./(Eupkq-Ednkq)*&
   (fermi(Eupkq)*Eupmdkq*gf_pp_kom(kx,ky,om3,sig2)-fermi(Ednkq)*Ednmdkq*gf_pp_kom(kx,ky,om4,sig2) ) 
  end function

  function gf_dd_kom(kx,ky,omega,sig)
     implicit none
     complex*16 gf_dd_kom
     real*8 kx,ky
     real*8 omega,sig
     complex*16 om
     om=dcmplx(omega,sig*del)
     gf_dd_kom=( om-Epk_val(kx,ky) )/(  (om-Edk_val(kx,ky))*(om-Epk_val(kx,ky))-tpdk_val(kx,ky)**2     )    
  end function

  function bubble_dd_kqom(kx,ky,qx,qy,omega)
     implicit none
     complex*16 bubble_dd_kqom
     real*8 kx,ky,qx,qy,omega 
     real*8 kmqx,kmqy
     real*8 Eupkval,Ednkval,Epkval,Eupmpk,Ednmpk
     real*8 Eupkq,Ednkq,Epkq,Eupmpkq,Ednmpkq
     real*8 om1,om2,om3,om4,sig1,sig2
     kmqx=kx-qx
     kmqy=ky-qy
     Eupkval=Eupk_val(kx,ky)
     Ednkval=Ednk_val(kx,ky)
     Epkval=Epk_val(kx,ky)
     Eupmpk=Eupkval-Epkval
     Ednmpk=Ednkval-Epkval
     Eupkq=Eupk_val(kmqx,kmqy)
     Ednkq=Ednk_val(kmqx,kmqy)
     Epkq=Epk_val(kmqx,kmqy)
     Eupmpkq=Eupkq-Epkq
     Ednmpkq=Ednkq-Epkq
     om1=Eupkval-omega
     om2=Ednkval-omega
     om3=Eupkq+omega
     om4=Ednkq+omega
     sig1=-1.0
     sig2=1.0 
     bubble_dd_kqom=1./(Eupkval-Ednkval)* &
   (fermi(Eupkval)*Eupmpk*gf_dd_kom(kmqx,kmqy,om1,sig1)-fermi(Ednkval)*Ednmpk*gf_dd_kom(kmqx,kmqy,om2,sig1)) &
                    +1./(Eupkq-Ednkq)*&
   (fermi(Eupkq)*Eupmpkq*gf_dd_kom(kx,ky,om3,sig2)-fermi(Ednkq)*Ednmpkq*gf_dd_kom(kx,ky,om4,sig2) ) 
  end function
  

  subroutine u_screened()
    implicit none
    real*8 kx,ky,qx,qy,omega
    complex*16 dp_bub,pp_bub
    complex*16 dys_inv(4,4)
    integer i,j,k,m,n,ind
    open(unit=9,FILE="bub_dp.out")
    open(unit=10,FILE="staticq_bub_dp.out")
    open(unit=11,FILE="bub_pp.out")
    open(unit=12,FILE="staticq_bub_pp.out")
    open(unit=13,FILE="uscr_dd.out")
    open(unit=14,FILE="uscr_stcq_dd.out")
    open(unit=15,FILE="bub_dd.out")
    open(unit=16,FILE="staticq_bub_dd.out")
    write(13,*)"#","omega","real(u_12)"   ,"imag(u_12)"  ,"real(u_21)"   ,"imag(u_21)"  ,&
    &"real(u_11)"   ,"imag(u_11)"   ,"real(u_22)"   ,"imag(u_22)"
    ind=1
    ! omega loop
    do k=1,Ne
      omega=Eom(k)
      bub_dp_om(k)=0.
      bub_pp_om(k)=0.
      bub_dd_om(k)=0.
      do m=1,4
         do n=1,4
            uscr_loc(m,n)=0.
         enddo
      enddo
      !! q loop
      do j=1,dim_k
         qx=qlist(j,1)
         qy=qlist(j,2)
	 bub_dp_qom(ind)=0.
	 bub_pp_qom(ind)=0.
	 bub_dd_qom(ind)=0.
         !!! k loop
         do i=1,dim_k
            kx=klist(i,1)
            ky=klist(i,2) 
            bub_dp_qom(ind)=bub_dp_qom(ind)+bubble_dp_kqom(kx,ky,qx,qy,omega) 
            bub_pp_qom(ind)=bub_pp_qom(ind)+bubble_pp_kqom(kx,ky,qx,qy,omega) 
            bub_dd_qom(ind)=bub_dd_qom(ind)+bubble_dd_kqom(kx,ky,qx,qy,omega) 
         enddo
         !!!
         bub_dp_qom(ind)=bub_dp_qom(ind)/dim_k
         bub_pp_qom(ind)=bub_pp_qom(ind)/dim_k
         bub_dd_qom(ind)=bub_dd_qom(ind)/dim_k
         dp_bub=bub_dp_qom(ind)
         pp_bub=bub_pp_qom(ind)
         call prest_init(dp_bub,pp_bub)
         dys_inv=iden-MATMUL(Prest,U)
         call Mat_inverse(dys_inv,4)
         uscr=MATMUL(U,dys_inv)
         if (omega==0.)then
           write(10,*)qx,qy,real(bub_dp_qom(ind)),imag(bub_dp_qom(ind))
           write(12,*)qx,qy,real(bub_pp_qom(ind)),imag(bub_pp_qom(ind))
           write(14,*)qx,qy,real(uscr(1,2)),imag(uscr(1,2))
           write(16,*)qx,qy,real(bub_dd_qom(ind)),imag(bub_dd_qom(ind))
           if (mod(j,Nk-1)==0)then
             write(10,*)
             write(12,*)
             write(14,*)
             write(16,*)
           endif 
         endif 
         bub_dp_om(k)=bub_dp_om(k)+bub_dp_qom(ind)
         bub_pp_om(k)=bub_pp_om(k)+bub_pp_qom(ind)
         bub_dd_om(k)=bub_dd_om(k)+bub_dd_qom(ind)
         uscr_loc=uscr_loc+uscr
         ind=ind+1
      enddo
      !!
      bub_dp_om(k)=bub_dp_om(k)/dim_k
      bub_pp_om(k)=bub_pp_om(k)/dim_k
      bub_dd_om(k)=bub_dd_om(k)/dim_k
      uscr_loc=uscr_loc/dim_k
      u11_scr(k)=uscr_loc(1,1)
      u12_scr(k)=uscr_loc(1,2)
      u21_scr(k)=uscr_loc(2,1)
      u22_scr(k)=uscr_loc(2,2)
      print *,ind
      write(9,*)omega,real(bub_dp_om(k)),imag(bub_dp_om(k))
      write(11,*)omega,real(bub_pp_om(k)),imag(bub_pp_om(k))
      write(13,*)omega,real(u12_scr(k)),imag(u12_scr(k)),real(u21_scr(k)),imag(u21_scr(k)),&
      &real(u11_scr(k)),imag(u11_scr(k)),real(u22_scr(k)),imag(u22_scr(k))
      write(15,*)omega,real(bub_dd_om(k)),imag(bub_dd_om(k))
    enddo
    !
  end subroutine

  subroutine prest_init(bub_dp,bub_pp)
     implicit none
     complex*16 bub_dp,bub_pp
     Prest(1,3)=bub_dp
     Prest(2,4)=bub_dp
     Prest(3,1)=bub_dp
     Prest(3,3)=bub_pp
     Prest(4,2)=bub_dp
     Prest(4,4)=bub_pp
     return
  end subroutine

  subroutine full_screened_w()
    implicit none
    integer i,j,k
    complex*16 :: dys_inv(2,2)
    complex*16 :: check(2,2)
    real*8 :: iden_mat(2,2)
    open(unit=29,FILE="Fullscr_dd.out")
    do i=1,2
       do j=1,2
        if(i==j)then
          iden_mat(i,j)=1.
        else
          iden_mat(i,j)=0.
        endif
       enddo
    enddo
    do k=1,Ne
       call Udd_mat_init(u11_scr(k),u12_scr(k),u21_scr(k),u22_scr(k)) 
       call Pdd_mat_init(bub_dd_om(k))
       dys_inv=iden_mat-MATMUL(Pdd,udd_mat)
       call Mat_inverse(dys_inv,2)
      ! check=MATMUL(Pdd,udd_mat)
        !wscr=MATMUL(udd_mat,MATMUL(check,check))
      !  wscr=udd_mat+MATMUL(udd_mat,check)+MATMUL(udd_mat,check*check)
       !MATMUL(udd_mat,check)
       wscr=MATMUL(udd_mat,dys_inv)
       !write(29,*)Eom(k),wscr(1,2) 
       write(29,*)Eom(k),real(wscr(1,2)),imag(wscr(1,2)) 
    enddo
  end subroutine 
  
  subroutine Udd_mat_init(u11,u12,u21,u22)
    implicit none
    complex*16 u11,u12,u21,u22
    udd_mat(1,1)=u11    !!!
    udd_mat(1,2)=u12
    udd_mat(2,1)=u21
    udd_mat(2,2)=u22     !!!
    return
  end subroutine
 
  subroutine Pdd_mat_init(p11)
    implicit none
    complex*16 p11
    Pdd(1,1)=p11
    Pdd(2,2)=p11
    Pdd(1,2)=0.
    Pdd(2,1)=0.
    return
  end subroutine
  subroutine Mat_inverse(Mat,n)
    implicit none
    integer n
    complex*16 :: Mat(n,n)
    complex*16 :: work(n)
    integer :: ipiv(n)
    integer info
    open(unit=999,FILE="info_inverse.out")
   ! allocate(Mat(n,n),work(n),ipiv(n))
    call ZGETRF(n,n,Mat,n,ipiv,info)
    if(info .eq. 0)then
      write(999,*)"succeded"
    else
      write(999,*)"failed"
    endif
    call ZGETRI(n,Mat,n,ipiv,work,n,info)
    if(info .eq. 0)then
      write(999,*)"succeded"
    else
      write(999,*)"failed"
    endif
    return
  end subroutine


end module
