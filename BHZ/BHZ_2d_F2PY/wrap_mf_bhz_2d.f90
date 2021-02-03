subroutine mf_bhz_2d(uloc,mh,res)
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  real(8),intent(in)                      :: uloc,mh
  integer,intent(out)                     :: res
  integer                                 :: Nparams
  integer,parameter                       :: Norb=2,Nspin=2,Nso=Nspin*Norb
  !
  character(len=20)                       :: Finput
  integer                                 :: Nkx
  real(8)                                 :: Jhratio
  integer 				  :: MaxIter
  real(8)                                 :: lambda,beta
  real(8)                                 :: wmix,it_error,sb_field
  !
  real(8)                                 :: Jh
  integer                                 :: Nky,Nktot
  real(8),dimension(:,:),allocatable      :: kgrid
  complex(8),dimension(:,:,:),allocatable :: Hk
  integer 				  :: Iter
  real(8)                                 :: Sz,Tz,error
  real(8)                                 :: E0,Ez,Ex
  logical                                 :: iexist,converged
  complex(8),dimension(Nso,Nso)           :: Gamma1,Gamma2,Gamma5,GammaS
  complex(8),dimension(Nso,Nso)           :: GammaN,GammaTz,GammaSz,GammaRz
  complex(8),dimension(Nso,Nso)           :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:),allocatable        :: params,params_prev
  real(8)                                 :: threshold


  call parse_cmd_variable(Finput,"FINPUT",default="inputBHZ.conf")
  call parse_input_variable(Nparams,"NPARAMS",Finput,default=6)
  call parse_input_variable(nkx,"NKX",Finput,default=25)
  call parse_input_variable(Jhratio,"JhRatio",Finput,default=0.25d0)
  call parse_input_variable(lambda,"LAMBDA",Finput,default=0.3d0)
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",Finput,default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD",Finput,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=1000)
  !

  select case(Nparams)
  case (1,2,5,6)
  case default
     stop "Wrong NPARAMS != [1,2,5,6]"
  end select
  !
  Jh = JhRatio*Uloc
  !
  Nky  = Nkx
  Nktot= Nkx*Nky
  !
  allocate( params(Nparams),params_prev(Nparams) )

  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaS=kron_pauli( pauli_sigma_z, pauli_tau_0)
  !

  gammaN=kron_pauli( pauli_sigma_0, pauli_tau_0 )
  gammaTz=kron_pauli( pauli_sigma_0, pauli_tau_z )
  gammaSz=kron_pauli( pauli_sigma_z, pauli_tau_0 )
  gammaRz=kron_pauli( pauli_sigma_z, pauli_tau_z )
  !
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron_pauli( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron_pauli( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron_pauli( pauli_sigma_z, pauli_tau_x )

  !Setup the k-space lattice basis:
  call TB_set_bk([pi2,0d0],[0d0,pi2])

  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)
  params_prev= 0d0
  params     = sb_field      ![Tz,Sz,S0,Ex,Ey,Ez]

  ! open(100,file="params_U"//str(Uloc)//"_M"//str(Mh))

  write(*,"(2F12.7)",advance="no")Uloc,Mh
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call solve_MF_bhz(iter,params)
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     converged = check_error(params,it_error,1,maxiter,error) 
  end do
  !
  threshold=1d-2
  res=1
  Tz = params(1)
  select case(Nparams)
  case(2)
     Sz = params(2)
     if(abs(Sz)>threshold)res=2
  case(5)
     Ex = sqrt(dot_product(params(5:6),params(5:6)))
     E0 = sqrt(dot_product(params(3:4),params(3:4)))
     if(abs(E0)>=threshold)res=4
     if(abs(Ex)>=threshold)res=3
  case(6)
     Sz = params(2)
     Ex = sqrt(dot_product(params(5:6),params(5:6)))
     E0 = sqrt(dot_product(params(3:4),params(3:4)))
     if(abs(Sz)>threshold)res=2
     if(abs(E0)>=threshold)res=4
     if(abs(Ex)>=threshold)res=3
  end select

  write(*,"(I4,4F21.12,I4,ES15.7)")iter,Tz,Sz,E0,Ex,res,error
  inquire(file="params.run",exist=iexist)
  if(iexist)then
     open(100,file="params.run",status="old", position="append", action="write")
  else
     open(100,file="params.run",status="new", action="write")
  endif
  write(100,*)Uloc,Mh,Tz,Sz,E0,Ex,res,jhratio
  close(100)

contains

  subroutine solve_MF_bhz(iter,a)
    real(8),dimension(:),intent(inout) :: a
    real(8),dimension(2)               :: kvec
    complex(8),dimension(Nso,Nso)      :: Hk,Hkmf
    real(8),dimension(Nso)             :: Ek,rhoDiag
    complex(8),dimension(Nso,Nso)      :: rhoHk,rhoH
    real(8)                            :: N,Tz,Sz,Rz
    real(8)                            :: E0,Ez,Ex,Ey
    integer                            :: ik,iter
    !
    rewind(100)
    !
    N  = 0d0
    Rz = 0d0
    Tz = 0d0
    Sz = 0d0
    E0 = 0d0
    Ez = 0d0
    Ex = 0d0
    Ey = 0d0
    !
    Hkmf=mf_hk_correction(a)
    !
    rhoH = 0d0
    do ik=1,Nktot
       kvec = kgrid(ik,:)             ![kx,ky]
       Hk   = hk_model(kvec,Nso)+Hkmf !H(k)
       !
       call eigh(Hk,Ek)       !diag Hk --> Ek
       !
       rhoDiag = fermi(Ek,beta)
       rhoHk   = matmul( Hk, matmul(diag(rhoDiag),conjg(transpose(Hk))) )
       rhoH    = rhoH + rhoHk/Nktot
    enddo
    !
    N  = N  + sum( GammaN*rhoH )
    Rz = Rz + sum( GammaRz*rhoH )
    Tz = Tz + sum( GammaTz*rhoH )
    Sz = Sz + sum( GammaSz*rhoH )
    E0 = E0 + sum( GammaE0*rhoH )
    Ez = Ez + sum( GammaEz*rhoH )
    Ex = Ex + sum( GammaEx*rhoH )
    Ey = Ey - dimag(sum( GammaEy*rhoH))
    !
    select case(Nparams)
    case(1)
       ! write(*,"(I4,1F21.12)")iter,Tz
       a = [Tz]
    case(2)
       ! write(*,"(I4,2F21.12)")iter,Tz,Sz
       a = [Tz,Sz]
    case(5)
       ! write(*,"(I4,6F21.12)")iter,Tz,E0,Ez,Ex,Ey
       a = [Tz,E0,Ez,Ex,Ey]
    case(6)
       ! write(*,"(I4,6F21.12)")iter,Tz,Sz,E0,Ez,Ex,Ey
       a = [Tz,Sz,E0,Ez,Ex,Ey]
    end select
    return
  end subroutine solve_MF_bhz


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    !
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
  end function hk_model

  function mf_Hk_correction(a) result(HkMF)
    real(8),dimension(:)          :: a
    complex(8),dimension(Nso,Nso) :: HkMF
    select case(Nparams)
    case(1)
       HkMF = -a(1)*(Uloc-5d0*Jh)/4d0*Gamma5
    case(2)
       HkMF = -a(1)*(Uloc-5d0*Jh)/4d0*Gamma5 &
            -a(2)*(Uloc+Jh)/4d0*GammaS
    case(5)
       HkMF = -a(1)*(Uloc-5d0*Jh)/4d0*Gamma5 &
            -a(2)*(Uloc-3d0*Jh)/4d0*GammaE0&
            -a(3)*(Uloc-3d0*Jh)/4d0*GammaEz&
            -a(4)*(Uloc-2d0*Jh)/4d0*GammaEx&
            -a(5)*(Uloc-2d0*Jh)/4d0*GammaEy
    case(6)
       HkMF = -a(1)*(Uloc-5d0*Jh)/4d0*Gamma5 &
            -a(2)*(Uloc+Jh)/4d0*GammaS &
            -a(3)*(Uloc-3d0*Jh)/4d0*GammaE0&
            -a(4)*(Uloc-3d0*Jh)/4d0*GammaEz&
            -a(5)*(Uloc-2d0*Jh)/4d0*GammaEx&
            -a(6)*(Uloc-2d0*Jh)/4d0*GammaEy
    end select
  end function mf_Hk_correction



  function check_error(Xnew,eps,N1,N2,error) result(convergence)
    real(8),intent(in)            :: Xnew(:)
    real(8)                       :: eps
    integer                       :: N1,N2
    integer                       :: Msize1
    logical                       :: convergence  
    real(8)                       :: err,error
    real(8),dimension(size(Xnew)) :: Verror
    real(8),save,allocatable      :: Xold(:,:)
    integer,save                  :: success=0,check=1
    Msize1=size(Xnew)
    if(.not.allocated(Xold))then
       allocate(Xold(1,Msize1))
       Xold=0.d0
    endif
    Verror=abs(Xnew-Xold(1,:))
    if(check==1)Verror=1d0
    err=sum(Verror)/dble(size(Verror))
    Xold(1,:)=Xnew
    if(err < eps)then
       success=success+1
    else
       success=0
    endif
    convergence=.false.
    if(success > N1)convergence=.true.
    ! if(convergence)then
    !    write(*,"(A,ES15.7,I8)")bold_green("    error="),err
    ! else
    !    if(err < eps)then
    !       write(*,"(A,ES15.7,I8)")bold_yellow("    error="),err
    !    else
    !       write(*,"(A,ES15.7,I8)")bold_red("    error="),err
    !    endif
    ! endif
    error=err
    ! if(check>=N2)convergence=.true.
    check=check+1

  end function check_error

end subroutine mf_bhz_2d


