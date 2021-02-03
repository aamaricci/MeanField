program mf_bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  integer                                       :: Nparams
  integer,parameter                             :: Norb=1,Nspin=2,Nlat=2,Nlso=Nlat*Nspin*Norb
  integer                                       :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                       :: Nky,Nx,Ny
  integer                                       :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                       :: ilat,jlat
  integer                                       :: ix,iy,iz
  real(8)                                       :: kx,ky,kz
  real(8),dimension(:,:),allocatable            :: kgrid,kpath,ktrims,Rgrid
  real(8),dimension(2)                          :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                          :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: bklen

  complex(8),dimension(:,:,:),allocatable       :: Hk
  complex(8),dimension(:,:,:,:),allocatable     :: Hlat
  integer 					:: Iter,MaxIter,Nsuccess=2
  real(8)                                       :: Uloc,Jh,JU,Sz,Tz
  !
  !variables for the model:
  character(len=32)                             :: finput
  real(8)                                       :: t1,t2,phi,Mh,wmixing,lambda
  real(8)                                       :: xmu,beta,eps
  real(8)                                       :: wmix,it_error,sb_field
  complex(8)                                    :: Hmf_glob(Nlso,Nlso)
  complex(8),dimension(:,:,:,:,:,:),allocatable :: Gmats,Greal
  character(len=20)                             :: file
  logical                                       :: iexist,converged,withgf,withbnd
  complex(8),dimension(Nlso,Nlso)               :: Gamma0,GammaX,GammaY,GammaZ,Gamma5
  complex(8),dimension(Nlso,Nlso)               :: GammaSz,GammaSx,GammaSy
  complex(8),dimension(Nlso,Nlso)               :: GammaRz,GammaRx,GammaRy
  real(8),dimension(:),allocatable              :: params,params_prev

  call parse_cmd_variable(Finput,"FINPUT",default="inputKM.conf")
  call parse_input_variable(Nparams,"NPARAMS",Finput,default=2,comment="2=AFMz,4=AFMxy,6=AFMxyz")
  call parse_input_variable(nkx,"NKX",Finput,default=25)
  call parse_input_variable(nkpath,"NKPATH",Finput,default=500)
  call parse_input_variable(L,"L",Finput,default=2048)
  call parse_input_variable(Uloc,"ULOC",Finput,default=1d0)
  call parse_input_variable(t1,"T1",finput,default=2d0,comment='NN hopping, fixes noninteracting bandwidth')
  call parse_input_variable(t2,"T2",finput,default=0d0,comment='Haldane-like NNN hopping-strenght, corresponds to lambda_SO in KM notation')
  call parse_input_variable(mh,"MH",Finput,default=0d0)
  call parse_input_variable(xmu,"XMU",Finput,default=0.d0)
  call parse_input_variable(eps,"EPS",Finput,default=4.d-2)
  call parse_input_variable(beta,"BETA",Finput,default=1000.d0)
  call parse_input_variable(wmix,"WMIX",Finput,default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD",Finput,default=0.1d0)
  call parse_input_variable(it_error,"IT_ERROR",Finput,default=1d-5)
  call parse_input_variable(maxiter,"MAXITER",Finput,default=1000)
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
  case default
     stop "Wrong NPARAMS != [1,2,5,6]"
  case (2)
     write(*,*)"Solving BHZ with Z-magnetization: Sz,Rz"
  case (4)
     write(*,*)"Solving BHZ with XY-magnetization: [Sx,Sy,Rx,Ry]"
  case (6)
     write(*,*)"Solving BHZ with XYZ-magnetization: [Sx,Sy,Sz,Rx,Ry,Rz] "
  end select

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

  !
  if(withgf)sb_field=0d0
  !
  Hmf_glob = zero
  params   = sb_field      ![Sx,Sy,Sz]
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
     allocate(Hk(Nlso,Nlso,Nktot))
     call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nkx])
     allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_matsubara(Hk,Gmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hk,Greal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     !
     call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=4)
     call dmft_print_gf_realaxis(Greal,"Greal",iprint=4)
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
  allocate(Kpath(4,2))
  KPath(1,:)=[0,0]
  KPath(2,:)=[2*pi/3, 2*pi/3/sqrt(3d0)]
  Kpath(3,:)=[2*pi/3,-2*pi/3/sqrt(3d0)]
  KPath(4,:)=[0d0,0d0]
  call TB_Solve_model(hk_model,Nlso,KPath,Nkpath,&
       colors_name=[red,red,blue,blue],&
       points_name=[character(len=10) :: "{/Symbol G}","K","K`","{/Symbol G}"],&
       file="Eigenbands.nint",iproject=.false.)


contains

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
    !
    Hk  = Hk + Hmf_glob
  end function hk_model



  subroutine solve_MF_km(iter,a)
    real(8),dimension(:),intent(inout) :: a
    real(8),dimension(2)               :: kvec
    complex(8),dimension(Nlso,Nlso)    :: Hk,Hkmf
    real(8),dimension(Nlso)            :: Ek,rhoDiag
    complex(8),dimension(Nlso,Nlso)    :: rhoHk,rhoH
    real(8)                            :: Sx,Sy,Sz,Rx,Ry,Rz
    real(8)                            :: na(2),nb(2),Sza,Szb
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
       write(*,"(I4,2F21.12)")iter,Sz,Rz
       write(100,"(2F21.12)")Sz,Rz
       a = [Sz,Rz]
    case(4)
       write(*,"(I4,4F21.12)")iter,Sx,Sy,Rx,Ry
       write(100,"(4F21.12)")Sx,Sy,Rx,Ry
       a = [Sx,Sy,Rx,Ry]
    case(6)
       write(*,"(I4,12F21.12)")iter,Sx,Sy,Sz,Rx,Ry,Rz
       write(100,"(6F21.12)")Sx,Sy,Sz,Rx,Ry,Rz
       a = [Sx,Sy,Sz,Rx,Ry,Rz]
    end select
    return
  end subroutine solve_MF_km



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




end program mf_bhz_2d


