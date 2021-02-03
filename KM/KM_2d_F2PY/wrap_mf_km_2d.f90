subroutine mf_km_2d(t2,uloc,res)
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  real(8),intent(in)                      :: t2,uloc
  integer,intent(out)                     :: res
  integer                                 :: Nparams
  integer,parameter                       :: Norb=1,Nspin=2,Nlat=2,Nlso=Nlat*Nspin*Norb
  !
  real(8),dimension(2)                    :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                    :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                    :: d1,d2,d3
  real(8),dimension(2)                    :: a1,a2,a3
  real(8),dimension(2)                    :: bklen
  !
  complex(8),dimension(:,:,:),allocatable :: Hk
  !
  character(len=20)                       :: Finput
  real(8)                                 :: t1
  real(8)                                 :: phi,Mh
  real(8)                                 :: wmix,lambda
  integer                                 :: Nkx
  integer 				  :: Iter,MaxIter,Nsuccess=2
  real(8)                                 :: beta
  real(8)                                 :: it_error,sb_field
  !
  real(8)                                 :: error
  integer                                 :: Nky,Nktot
  real(8),dimension(:,:),allocatable      :: kgrid
  real(8)                                 :: Sx,Sy,Sz,Rx,Ry,Rz,absS,absR
  logical                                 :: iexist,converged
  complex(8),dimension(Nlso,Nlso)         :: Gamma0,GammaX,GammaY,GammaZ,Gamma5
  complex(8),dimension(Nlso,Nlso)         :: GammaSz,GammaSx,GammaSy
  complex(8),dimension(Nlso,Nlso)         :: GammaRz,GammaRx,GammaRy
  real(8),dimension(:),allocatable        :: params,params_prev
  real(8)                                 :: threshold


  call parse_cmd_variable(Finput,"FINPUT",default="inputKM.conf")
  call parse_input_variable(Nparams,"NPARAMS",Finput,default=2,comment="2=AFMz,4=AFMxy,6=AFMxyz")
  call parse_input_variable(nkx,"NKX",Finput,default=25)
  call parse_input_variable(t1,"T1",finput,default=2d0)
  call parse_input_variable(mh,"MH",Finput,default=0d0)
  call parse_input_variable(lambda,"LAMBDA",Finput,default=0.3d0)
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",Finput,default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD",Finput,default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=1000)
  !

  select case(Nparams)
  case (2,4,6)
  case default
     stop "Wrong NPARAMS != [1,2,5,6]"
  end select
  !
  Nky  = Nkx
  Nktot= Nkx*Nky
  !
  allocate( params(Nparams),params_prev(Nparams) )

  !SETUP THE GAMMA MATRICES:
  !we must use the basis \Gamma_ab = \tau_a \circ \sigma_b
  ! tau_a   -> lattice
  ! sigma_b -> spin
  ! \psi = [A_up, A_dw; B_up, B_dw]^T
  !This convention is dictated by the use of DMFT_TOOLS functions
  gamma0=kron_pauli( pauli_tau_0, pauli_sigma_0) !G_00
  gammaZ=kron_pauli( pauli_tau_z, pauli_sigma_z) !G_33
  gammaX=kron_pauli( pauli_tau_x, pauli_sigma_0) !G_10
  gammaY=kron_pauli( pauli_tau_y, pauli_sigma_0) !G_20
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z) !G_03
  !
  gammaSx=kron_pauli( pauli_tau_0, pauli_sigma_x )
  gammaSy=kron_pauli( pauli_tau_0, pauli_sigma_y )
  gammaSz=kron_pauli( pauli_tau_0, pauli_sigma_z )
  gammaRx=kron_pauli( pauli_tau_z, pauli_sigma_x )
  gammaRy=kron_pauli( pauli_tau_z, pauli_sigma_y )
  gammaRz=kron_pauli( pauli_tau_z, pauli_sigma_z )

  !Lattice basis (a=1; a0=sqrt3*a) is:
  !e_1 = a0 [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
  !e_2 = a0 [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
  e1 = 3d0/2d0*[1d0, 1d0/sqrt(3d0)]
  e2 = 3d0/2d0*[1d0,-1d0/sqrt(3d0)]

  !LATTICE BASIS: nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]

  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d1-d3                    !3/2*a[1,1/sqrt3]
  a2 = d2-d3                    !3/2*a[1,-1/sqrt3]
  a3 = d1-d2

  !RECIPROCAL LATTICE VECTORS:
  bklen=2d0*pi/3d0
  bk1=bklen*[ 1d0, sqrt(3d0)] 
  bk2=bklen*[ 1d0,-sqrt(3d0)]
  call TB_set_bk(bkx=bk1,bky=bk2)

  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)



  params_prev= 0d0
  params     = sb_field      ![Svec,Rvec]

  write(*,"(2F12.7)",advance="no")Uloc,t2
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call solve_MF_km(iter,params)
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     converged = check_error(params,it_error,1,maxiter,error) 
  end do
  !
  threshold=1d-2
  res=1
  select case(Nparams)
  case(2)
     Sz = params(1)
     Rz = params(2)
     if(abs(Rz)>threshold)res=2
     write(*,"(I4,2F21.12,I4,ES15.7)")iter,params(1:2),res,error
  case(4)
     absS = sqrt(dot_product(params(1:2),params(1:2)))
     absR = sqrt(dot_product(params(3:4),params(3:4)))
     if(abs(absR)>=threshold)res=3
     if(abs(absS)>=threshold)res=4
     write(*,"(I4,4F21.12,I4,ES15.7)")iter,params(1:4),res,error
  case(6)
     Sz = params(1)
     Rz = params(2)
     absS = sqrt(dot_product(params(3:4),params(3:4)))
     absR = sqrt(dot_product(params(5:6),params(5:6)))
     if(abs(Rz)>threshold)res=2
     if(abs(absR)>=threshold)res=3
     if(abs(absS)>=threshold)res=4
     write(*,"(I4,6F21.12,I4,ES15.7)")iter,params(1:6),res,error
  end select


  inquire(file="params.run",exist=iexist)
  if(iexist)then
     open(100,file="params.run",status="old", position="append", action="write")
  else
     open(100,file="params.run",status="new", action="write")
  endif
  write(100,*)Uloc,t2,params,res
  close(100)

contains

  subroutine solve_MF_km(iter,a)
    real(8),dimension(:),intent(inout) :: a
    real(8),dimension(2)               :: kvec
    complex(8),dimension(Nlso,Nlso)    :: Hk,Hkmf
    real(8),dimension(Nlso)            :: Ek,rhoDiag
    complex(8),dimension(Nlso,Nlso)    :: rhoHk,rhoH
    real(8)                            :: Sx,Sy,Sz,Rx,Ry,Rz
    integer                            :: ik,iter
    !
    rewind(100)
    !
    Sx = 0d0; Rx = 0d0
    Sy = 0d0; Ry = 0d0
    Sz = 0d0; Rz = 0d0
    !
    Hkmf=mf_hk_correction(a)
    !
    rhoH = 0d0
    do ik=1,Nktot
       kvec = kgrid(ik,:)             ![kx,ky]
       Hk   = hk_model(kvec,Nlso)+Hkmf !H(k)
       !
       call eigh(Hk,Ek)       !diag Hk --> Ek
       !
       rhoDiag = fermi(Ek,beta)
       rhoHk   = matmul( Hk, matmul(diag(rhoDiag),conjg(transpose(Hk))) )
       rhoH    = rhoH + rhoHk/Nktot
    enddo
    !
    Sx = Sx + sum( GammaSx*rhoH )
    Sy = Sy - sum( GammaSy*rhoH )
    Sz = Sz + sum( GammaSz*rhoH )
    Rx = Rx + sum( GammaRx*rhoH )
    Ry = Ry - sum( GammaRy*rhoH )
    Rz = Rz + sum( GammaRz*rhoH )
    !
    !
    select case(Nparams)
    case(2)
       a = [Sz,Rz]
    case(4)
       a = [Sx,Sy,Rx,Ry]
    case(6)
       a = [Sx,Sy,Sz,Rx,Ry,Rz]
    end select
    return
  end subroutine solve_MF_km


  !Using \psi=[A_up,A_dw;B_up,B_dw]
  function hk_model(kpoint,N) result(hk)    
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    complex(8),dimension(N,N) :: hk
    real(8)                   :: h0,hx,hy,hz
    real(8)		      :: kdote1, kdote2
    !
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    hx = t1*( cos(kdote1) + cos(kdote2) + 1)
    hy = t1*( sin(kdote1) + sin(kdote2) )
    hz = 2*t2*( sin(kdote1) - sin(kdote2) - sin(kdote1-kdote2) )
    !
    Hk = hx*GammaX + hy*GammaY + hz*GammaZ + Mh*Gamma5
  end function hk_model


  function mf_Hk_correction(a) result(HkMF)
    real(8),dimension(:)            :: a
    complex(8),dimension(Nlso,Nlso) :: HkMF
    select case(Nparams)
    case(2)
       HkMF = a(1)*GammaSz + a(2)*GammaRz
    case(4)
       HkMF = a(1)*GammaSx + a(2)*GammaSy + a(3)*GammaRx + a(4)*GammaRy
    case(6)
       HkMF = a(1)*GammaSx + a(2)*GammaSy &
            + a(3)*GammaRx + a(4)*GammaRy &
            + a(5)*GammaSz + a(6)*GammaRz
    end select
    HkMF=-HkMf*Uloc/4d0
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

end subroutine mf_km_2d


