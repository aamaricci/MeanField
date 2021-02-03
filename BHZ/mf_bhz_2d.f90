program mf_bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  ! integer,parameter                             :: Nparams=6 ![Sz,Tz,S0,Tx,Ty,Tz]
  integer :: Nparams
  integer,parameter                             :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                       :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                       :: Nky,Nlat,Nx,Ny
  integer                                       :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                       :: ilat,jlat
  integer                                       :: ix,iy,iz
  real(8)                                       :: kx,ky,kz
  real(8),dimension(:,:),allocatable            :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable            :: Links
  complex(8),dimension(:,:,:),allocatable       :: Hk
  complex(8),dimension(:,:,:,:),allocatable     :: Hlat
  real(8),dimension(:),allocatable              :: Wtk
  integer 					:: Iter,MaxIter,Nsuccess=2
  real(8)                                       :: chern,z2,Uloc,Jh,JU,Sz,Tz,Rz,Ntot,Jx,Jp
  real(8)                                       :: mh,rh,lambda,delta
  real(8)                                       :: xmu,beta,eps
  real(8)                                       :: wmix,it_error,sb_field
  complex(8)                                    :: Hloc(Nso,Nso),Hmf_glob(Nso,Nso)
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gmats,Greal
  character(len=20)                             :: Finput
  logical                                       :: iexist,converged,withgf
  complex(8),dimension(Nso,Nso)                 :: Gamma1,Gamma2,Gamma5,GammaS
  complex(8),dimension(Nso,Nso)                 :: GammaN,GammaTz,GammaSz,GammaRz
  complex(8),dimension(Nso,Nso)                 :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:),allocatable              :: params,params_prev

  call parse_cmd_variable(Finput,"FINPUT",default="inputBHZ.conf")
  call parse_input_variable(Nparams,"NPARAMS",Finput,default=2)
  call parse_input_variable(nkx,"NKX",Finput,default=25)
  call parse_input_variable(nkpath,"NKPATH",Finput,default=500)
  call parse_input_variable(L,"L",Finput,default=2048)
  call parse_input_variable(Uloc,"ULOC",Finput,default=1d0)
  call parse_input_variable(Jh,"Jh",Finput,default=0.25d0)
  call parse_input_variable(Jx,"Jx",Finput,default=0.d0)
  call parse_input_variable(Jp,"Jp",Finput,default=0.d0)
  call parse_input_variable(mh,"MH",Finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",Finput,default=0.3d0)
  call parse_input_variable(xmu,"XMU",Finput,default=0.d0)
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",Finput,default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD",Finput,default=0.1d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=100)
  call parse_input_variable(withgf,"WITHGF",Finput,default=.false.)
  !
  call print_input(trim(Finput))
  call save_input_file(trim(Finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")

  select case(Nparams)
  case (1)
     write(*,*)"Solving BHZ with Orbital polarization: Tz"
  case (2)
     write(*,*)"Solving BHZ with Hartree terms: [Tz,Sz]"
  case(5)
     write(*,*)"Solving BHZ with non-magnetic Excitonic params: [Tz,E0,Ex,Ey,Ez] "
  case(6)
     write(*,*)"Solving BHZ with full Excitonic order and magnetic params: [Tz,Sz,E0,Ex,Ey,Ez] "
  case default
     stop "Wrong NPARAMS != [1,2,5,6]"
  end select

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

  !
  if(withgf)sb_field=0d0
  !
  Hmf_glob = zero
  params   = sb_field      ![Tz,Sz,S0,Ex,Ey,Ez]
  inquire(file="params.restart",exist=iexist)
  if(iexist)then
     call read_array("params.restart",params)
     !try to break other symmetry anyway:
     if(Nparams>1)params(2:)=params(2:)+sb_field
  endif
  call save_array("params.init",params) 	!FP: verifying the initial params


  if(withgf)then
     !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
     Hmf_glob =  mf_hk_correction(params)
     call TB_write_hloc(Hmf_glob)
     !
     allocate(Hk(Nso,Nso,Nktot))
     call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx],wdos=.false.)
     allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
     allocate(Greal(Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_matsubara(Hk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
     !
     call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=1)
     call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
     stop
  endif





  open(100,file="order_parameters_tz_sz_e0_ex_ey_ez.dat")
  converged=.false. ; iter=0
  do while(.not.converged.AND.iter<maxiter)
     iter=iter+1
     call start_loop(iter,maxiter,"MF-loop")
     !
     call solve_MF_bhz(iter,params)
     if(iter>1)params = wmix*params + (1d0-wmix)*params_prev
     params_prev = params
     !
     converged = check_convergence_local(params,it_error,nsuccess,maxiter) 
     !
     call end_loop
  end do
  call save_array("params.restart",params)
  close(100)
  !
  !Update Global Mean-Field Hamiltonian correction:
  Hmf_glob =  mf_hk_correction(params)
  call TB_write_hloc(Hmf_glob)

  !SOLVE ALONG A PATH IN THE BZ.
  Npts=5
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_X1
  kpath(2,:)=kpoint_Gamma
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_X1
  kpath(5,:)=kpoint_Gamma
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: 'X', 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")

  !GET Z2 INVARIANT:
  allocate(ktrims(2,4))
  ktrims=reshape( [ [0d0,0d0] , [0d0,pi] , [pi,0d0] , [pi,pi] ] , shape(ktrims))
  call get_z2_number(ktrims,[2,4],z2)


contains

  subroutine solve_MF_bhz(iter,a)
    real(8),dimension(:),intent(inout) :: a
    real(8),dimension(2)               :: kvec
    complex(8),dimension(Nso,Nso)      :: Hk,Hkmf
    real(8),dimension(Nso)             :: Ek,rhoDiag
    real(8),dimension(Nso,Nso)         :: rhoHk,rhoH
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
       !
       ! N  = N  + sum( GammaN*rhoHk )/Nktot
       ! Rz = Rz + sum( GammaRz*rhoHk )/Nktot
       ! !
       ! Tz = Tz + sum( GammaTz*rhoHk )/Nktot
       ! Sz = Sz + sum( GammaSz*rhoHk )/Nktot
       ! E0 = E0 + sum( GammaE0*rhoHk )/Nktot
       ! Ez = Ez + sum( GammaEz*rhoHk )/Nktot
       ! Ex = Ex + sum( GammaEx*rhoHk )/Nktot
       ! Ey = Ey - dimag(sum( GammaEy*rhoHk ))/Nktot
       !
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
       write(*,"(I4,1F21.12)")iter,Tz
       write(100,"(1F21.12)")Tz
       a = [Tz]
    case(2)
       write(*,"(I4,2F21.12)")iter,Tz,Sz
       write(100,"(2F21.12)")Tz,Sz
       a = [Tz,Sz]
    case(5)
       write(*,"(I4,6F21.12)")iter,Tz,E0,Ez,Ex,Ey
       write(100,"(6F21.12)")Tz,E0,Ez,Ex,Ey
       a = [Tz,E0,Ez,Ex,Ey]
    case(6)
       write(*,"(I4,6F21.12)")iter,Tz,Sz,E0,Ez,Ex,Ey
       write(100,"(6F21.12)")Tz,Sz,E0,Ez,Ex,Ey
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
    !
    Hk  = Hk + Hmf_glob
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


  subroutine get_z2_number(ktrims,band_indices,z2)
    real(8),dimension(:,:),intent(in)       :: ktrims
    integer,dimension(:),intent(in)         :: band_indices
    complex(8),dimension(:,:,:),allocatable :: Htrims
    real(8),dimension(:,:),allocatable      :: Etrims
    complex(8),dimension(:),allocatable     :: Delta
    real(8)                                 :: z2
    integer                                 :: i,j,Ntrim,itrim,Nocc,unit
    !
    Ntrim=size(Ktrims,2)
    Nocc = size(band_indices)
    allocate(Htrims(Nso,Nso,Ntrim),Etrims(Nocc,Ntrim))
    allocate(Delta(Ntrim))
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_model(Ktrims(:,itrim),Nso)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    z2=product(Delta(:))
    if(z2>0)then
       z2=0d0
    else
       z2=1d0
    end if
    open(free_unit(unit),file="z2_invariant.dat")
    write(unit,*) z2
    close(unit)
  end subroutine get_z2_number



end program mf_bhz_2d


