program mf_bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

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
  integer                                       :: Iter,MaxIter,Nsuccess=2
  real(8)                                       :: chern,z2,Uloc,Jh,JU,Sz,Tz,Rz,Ntot
  real(8)                                       :: mh,rh,lambda,delta
  real(8)                                       :: xmu,beta,eps,Eout(2)
  real(8)                                       :: dens(Nso),wmix,it_error,sb_field
  complex(8)                                    :: Hloc(Nso,Nso)
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gmats,Greal
  character(len=20)                             :: file
  logical                                       :: iexist,converged,withgf
  complex(8),dimension(Nso,Nso)                 :: Gamma1,Gamma2,Gamma5,GammaS
  real(8),dimension(Nso)                        :: params,params_prev,global_params


  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(L,"L","inputBHZ.conf",default=2048)
  call parse_input_variable(Uloc,"ULOC","inputBHZ.conf",default=1d0)
  call parse_input_variable(Jh,"Jh","inputBHZ.conf",default=0.25d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=1d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(wmix,"WMIX","inputBHZ.conf",default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD","inputBHZ.conf",default=0.1d0)
  call parse_input_variable(it_error,"IT_ERROR","inputBHZ.conf",default=1d-5)
  call parse_input_variable(maxiter,"MAXITER","inputBHZ.conf",default=100)
  call parse_input_variable(withgf,"WITHGF","inputBHZ.conf",default=.false.)
  !
  call save_input_file("inputBHZ.conf")
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")


  Nky  = Nkx
  Nktot= Nkx*Nky
  !

  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)
  gammaS=kron_pauli( pauli_tau_z, pauli_sigma_0)

  !Setup the k-space lattice basis:
  call TB_set_bk([pi2,0d0],[0d0,pi2])

  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)


  global_params = 0d0
  params        = [sb_field,sb_field,0d0,2d0]   ![Sz,Tz,Rz,N]
  !
  inquire(file="params.restart",exist=iexist)
  if(iexist)call read_array("params.restart",params)
  params(1)=params(1)+sb_field

  open(100,file="sz.dat")
  open(101,file="tz.dat")
  open(102,file="dens.dat")
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
  global_params = params
  close(100)
  close(101)
  close(102)

  if(withgf)then
     !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
     write(*,*) "Using Nk_total="//txtfy(Nktot)
     allocate(Hk(Nso,Nso,Nktot))
     allocate(Wtk(Nktot))
     call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx],wdos=.false.)
     Wtk = 1d0/Nktot
     allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
     allocate(Greal(Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
     !
     call dmft_print_gf_matsubara(Gmats,"Gmats",iprint=1)
     call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)

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
  endif






  !GET Z2 INVARIANT:
  allocate(ktrims(2,4))
  ktrims=reshape( [ [0d0,0d0] , [0d0,pi] , [pi,0d0] , [pi,pi] ] , shape(ktrims))
  call get_z2_number(ktrims,[2,4],z2)
  print*,z2


contains

  subroutine solve_MF_bhz(iter,a)
    real(8),dimension(Nso),intent(inout) :: a
    real(8),dimension(2)                 :: kvec
    complex(8),dimension(Nso,Nso)        :: Hk
    real(8),dimension(Nso)               :: Ek,rhoDiag
    real(8),dimension(Nso,Nso)           :: rhoHk
    real(8)                              :: dens(2,2)
    real(8)                              :: Tz,Sz,N
    integer                              :: ik,iter,iorb,ispin
    !
    Sz = a(1) 
    Tz = a(2)
    N  = a(4)
    !
    dens = 0d0
    do ik=1,Nktot
       kvec = kgrid(ik,:)        ![kx,ky]
       Hk   = hk_model(kvec,Nso) !H(k)
       Hk   = Hk + mf_hk_correction(a)
       !
       call eigh(Hk,Ek)       !diag Hk --> Ek
       !
       rhoDiag = fermi(Ek,beta)
       rhoHk   = matmul( Hk, matmul(diag(rhoDiag),conjg(transpose(Hk))) )
       !
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb+(ispin-1)*Norb
             dens(ispin,iorb) = dens(ispin,iorb) + rhoHk(io,io)
          enddo
       enddo
       !
    enddo
    dens = dens/Nktot
    !
    a(1) = 0.5d0*sum(dens(1,:)) - 0.5d0*sum(dens(2,:)) !N_up - N_dw
    a(2) = 0.5d0*sum(dens(:,1)) - 0.5d0*sum(dens(:,2)) !N_1  - N_2
    print*,iter,a(1),a(2)
    !
    rewind(100);rewind(101);rewind(102)
    write(100,*)a(1)
    write(101,*)a(2)
    write(102,*)(dens(1,iorb),iorb=1,Norb),(dens(2,iorb),iorb=1,Norb)
  end subroutine solve_MF_bhz


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8) :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
    !
    Hk = Hk + mf_hk_correction(global_params)
  end function hk_model

  function mf_Hk_correction(a) result(HkMF)
    real(8),dimension(Nso)        :: a
    complex(8),dimension(Nso,Nso) :: HkMF
    HkMF = -a(2)*(Uloc - 5d0*Jh)/2d0*Gamma5 - a(1)*(Uloc+Jh)/2d0*GammaS
    ! HkMF = HkMF + a(4)*(3*Uloc - 5*Jh)/4d0*diag(ones(Nso))
    ! HkMF = HkMF + (Uloc+Jh)/2d0*a(1)**2 + (Uloc-5d0*Jh)/2d0*a(2)**2
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


