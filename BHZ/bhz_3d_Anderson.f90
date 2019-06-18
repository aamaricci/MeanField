program bhz_3d_disorder
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                             :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                       :: Nk,Nktot,Nkx,Nky,Nkz,L
  integer                                       :: Nlat,Nx,Ny,Nz
  integer                                       :: Nkpath,Npts,z2(4)
  integer                                       :: i,j,k,ik
  integer                                       :: io,jo,ispin
  integer                                       :: iorb,ilat,jorb,jlat
  integer                                       :: ix,iy,iz
  real(8)                                       :: kx,ky,kz
  real(8),dimension(:),allocatable              :: Wtk,Evals,rhoDiag
  real(8),dimension(:,:),allocatable            :: kpath,ktrims,Kgrid,Rgrid,Ldens
  integer,dimension(:,:),allocatable            :: Links
  complex(8),dimension(:,:,:),allocatable       :: Hij
  complex(8),dimension(:,:,:,:),allocatable     :: Hlat
  integer                                       :: Iter,MaxIter,Nsuccess=2,Idum
  real(8)                                       :: chern,Uloc,Jh,JU,Sz,Tz,Rz,Ntot
  real(8)                                       :: ez,mh,rh,lambda,delta,lz,Wdis
  real(8)                                       :: xmu,beta,eps,Ekin,Eloc
  real(8)                                       :: wmix,it_error,sb_field
  real(8),allocatable,dimension(:)              :: erandom
  complex(8),dimension(:,:,:,:,:,:),allocatable :: GLmats,GLreal
  character(len=20)                             :: file
  logical                                       :: iexist,converged,withgf,bool
  complex(8),dimension(Nso,Nso)                 :: Gamma5,GammaX,GammaY,GammaZ,GammaS,Gamma0
  real(8),dimension(3)                          :: vecK,vecRi,vecRj
  complex(8),dimension(:,:),allocatable         :: rhoH
  real(8),dimension(:,:),allocatable            :: params,params_prev ![Nlat,Nso]==densities
  integer        :: disorder_type !0=scalar disorder, 1=magnetic_disorder

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(L,"L","inputBHZ.conf",default=2048)
  call parse_input_variable(Wdis,"WDIS","inputBHZ.conf",default=0d0)
  call parse_input_variable(idum,"IDUM","inputBHZ.conf",default=1234567)
  call parse_input_variable(disorder_type,"DISORDER_TYPE","inputBHZ.conf",default=0)
  call parse_input_variable(ez,"ez","inputBHZ.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(wmix,"WMIX","inputBHZ.conf",default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD","inputBHZ.conf",default=0.01d0)
  call parse_input_variable(it_error,"IT_ERROR","inputBHZ.conf",default=1d-5)
  call parse_input_variable(maxiter,"MAXITER","inputBHZ.conf",default=100)
  call parse_input_variable(withgf,"WITHGF","inputBHZ.conf",default=.false.)
  call save_input_file("inputBHZ.conf")
  call print_input()
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  !< Momentum space:
  Nky = Nkx
  Nkz = Nkx
  Nktot=Nkx*Nkx*Nkx
  !
  !< Real space
  Nx = Nkx
  Ny = Nkx
  Nz = Nkx
  Nlat = Nx*Ny*Nz

  !SETUP THE GAMMA MATRICES:
  gamma0=kron_pauli( pauli_tau_0, pauli_sigma_0)
  gammaX=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gammaY=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gammaZ=kron_pauli( pauli_tau_x, pauli_sigma_x)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)
  gammaS=kron_pauli( pauli_tau_z, pauli_sigma_0)


  !< start using TB_procedures:
  !< 1st set up the direct and momentum space lattice basis
  call TB_set_ei([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])


  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  Hlat=zero
  !
  allocate(Links(6,3))
  Links(1,:) = [1,0,0]
  Links(2,:) = [0,1,0]
  Links(3,:) = [-1,0,0]
  Links(4,:) = [0,-1,0]
  Links(5,:) = [0,0,1]
  Links(6,:) = [0,0,-1]
  print*,"Building up Hlat model"
  call TB_build_model(Hlat,ts_model,Nso,[Nkx,Nky,Nkz],Links,pbc=.true.)

  !< Build up disorder:
  allocate(erandom(Nlat))
  call mersenne_init(idum)  
  call mt_random(erandom)
  erandom=(2d0*erandom-1d0)*Wdis/2d0
  inquire(file='erandom.restart',exist=bool)
  if(bool)then
     if(file_length('erandom.restart')/=Nlat)stop "bhz_2d_anderson error: found erandom.restart with length different from Nlat"
     call read_array('erandom.restart',erandom)
  endif
  call save_array("erandom.used",erandom)

  allocate(Hij(Nlat*Nso,Nlat*Nso,1))
  Hij = zero
  !
  select case(disorder_type)
  case(0)
     do ilat=1,Nlat
        Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*Gamma0
     enddo
  case(1)
     do ilat=1,Nlat
        Hlat(:,:,ilat,ilat) = Hlat(:,:,ilat,ilat) + erandom(ilat)*GammaS
     enddo
  case default
     stop "disorder_type not [0-1]"
  end select
  !
  do ilat=1,Nlat
     do jlat=1,Nlat
        do io=1,Nso
           do jo=1,Nso
              i = io + (ilat-1)*Nso
              j = jo + (jlat-1)*Nso
              Hij(i,j,1) = Hlat(io,jo,ilat,jlat)
           enddo
        enddo
     enddo
  enddo

  open(100,file="sz.dat")
  open(101,file="tz.dat")
  open(102,file="dens.dat")
  !
  call solve_Anderson_bhz()
  !;
  close(100);close(101);close(102)


  !< BUILD THE LOCAL GF
  if(withgf)then
     allocate(GLmats(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(GLreal(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_matsubara(Hij,[1d0],GLmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hij,[1d0],GLreal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_print_gf_matsubara(GLmats,"Gloc",iprint=1)
     call dmft_print_gf_realaxis(GLreal,"Gloc",iprint=1)
  endif




contains


  subroutine solve_Anderson_bhz()
    real(8),dimension(Nlat)                   :: Tzii,Szii
    complex(8),dimension(Nlat*Nso,Nlat*Nso)   :: Hmat
    real(8),dimension(Nlat*Nso)               :: Evals,rhoDiag
    real(8),dimension(Nlat*Nso,Nlat*Nso)      :: rhoH
    real(8)                                   :: dens(Nlat,Nspin,Norb)
    integer                                   :: iter,iorb,ispin
    !
    dens = 0d0
    !
    Hmat = Hij(:,:,1)
    !
    call my_eigh(Hmat,Evals)       !diag Hij -> Ea_i
    !
    rhoDiag = fermi(Evals,beta)
    rhoH    = matmul(Hmat , matmul(diag(rhoDiag), conjg(transpose(Hmat))) )
    !
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             io = iorb+(ispin-1)*Norb+(ilat-1)*Nspin*Norb
             dens(ilat,ispin,iorb) = rhoH(io,io)
          enddo
       enddo
    enddo
    !
    !< get observables:
    do ilat=1,Nlat
       Szii(ilat) = 0.5d0*sum(dens(ilat,1,:)) - 0.5d0*sum(dens(ilat,2,:)) !N_up - N_dw
       Tzii(ilat) = 0.5d0*sum(dens(ilat,:,1)) - 0.5d0*sum(dens(ilat,:,2)) !N_1  - N_2
    enddo
    !
    rewind(100);rewind(101);rewind(102)
    do ilat=1,Nlat
       write(100,*)ilat,Szii(ilat)
       write(101,*)ilat,Tzii(ilat)
       write(102,*)ilat,(dens(ilat,1,iorb),iorb=1,Norb),(dens(ilat,2,iorb),iorb=1,Norb)
    enddo
    !
  end subroutine solve_Anderson_bhz



  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*Gamma5 
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaX
    case (2) !UP HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaY
    case (3) !LEFT HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaX
    case (4) !DOWN HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaY
    case (5) !TOP HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*GammaZ
    case (6) !BOTTOM HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*GammaZ
    case default 
       stop "ts_model ERROR: link index in {0..6}"
    end select
  end function ts_model







  subroutine my_eigh(A,W,jobz,uplo,vl,vu,il,iu,tol)
    complex(8),dimension(:,:),intent(inout)    :: A ! M v = E v/v(i,j) = ith component of jth vec.
    real(8),dimension(size(A,2)),intent(inout) :: W ! eigenvalues
    character(len=1),optional                  :: jobz,uplo
    character(len=1)                           :: jobz_,uplo_,range
    real(8),optional                           :: vl,vu,tol
    integer,optional                           :: il,iu
    real(8)                                    :: vL_,vU_,tol_
    integer                                    :: iL_,iU_
    integer                                    :: i,j,n,lda,info,ldz
    integer                                    :: lwork,liwork,lrwork,m
    integer                                    :: guess_iwork(1)
    real(8)                                    :: guess_rwork(1)
    complex(8)                                 :: guess_work(1)
    logical                                    :: boolV,boolI
    complex(8),dimension(:,:),allocatable      :: Z
    complex(8),dimension(:),allocatable        :: work
    complex(8),dimension(:),allocatable        :: rwork
    integer,dimension(:),allocatable           :: Isuppz
    integer,dimension(:),allocatable           :: Iwork
    !
    !
    jobz_='V' ;if(present(jobz))jobz_=jobz
    uplo_='U' ;if(present(uplo))uplo_=uplo
    vl_  = 1d0;if(present(vL))vL_=vL
    vu_  = 1d0;if(present(vU))vU_=vU
    iL_  = 1  ;if(present(iL))iL_=iL
    iU_  = 1  ;if(present(iU))iU_=iU
    tol_ = 0d0;if(present(tol))tol_=tol
    !
    range='A'
    boolV=present(vL).AND.present(vU)
    boolI=present(iL).AND.present(iU)
    if(boolV.and.boolI)stop "vX and iX arguments both present. Can not set RANGE"
    if(boolV)range='V'
    if(boolI)range='I'
    !
    lda   = max(1,size(A,1))
    n     = max(1,size(A,2))
    m     = n
    !
    allocate(Isuppz(2*m))
    !
    ldz   = n
    allocate(Z(ldz,n))
    !
    call zheevr('V','A','U',N,A,lda,vl_,vu_,iL_,iU_,tol_,M,W,Z,ldz,Isuppz,&
         guess_work,-1,&
         guess_rwork,-1,&
         guess_iwork,-1,info)
    lwork = guess_work(1)
    lrwork= guess_rwork(1)
    liwork= guess_iwork(1)
    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))
    call zheevr('V','A','U',N,A,lda,vl_,vu_,iL_,iU_,tol_,M,W,Z,ldz,Isuppz,&
         work,lwork,&
         rwork,lrwork,&
         iwork,liwork,info)
    !<copy the Evecs from Z to the input matrix A
    A = Z
  end subroutine my_eigh



end program bhz_3d_disorder







! !Build up the real-space Hamiltonian thru FT:"
! allocate(Hij(Nlat*Nso,Nlat*Nso,1)) ; Hij = zero
! call start_timer
! do ilat=1,Nlat
!    vecRi = Rgrid(ilat,:)
!    do jlat=1,Nlat
!       vecRj = Rgrid(jlat,:)
!       !
!       Htmp = zero
!       do ik=1,Nktot
!          vecK = Kgrid(ik,:)
!          arg=dot_product(vecK,vecRj-vecRi)
!          Htmp(:,:)= Htmp(:,:) + exp(xi*arg)*hk_model(vecK,Nso)/Nktot
!       enddo
!       !
!       do io=1,Nso
!          i = io + (ilat-1)*Nso
!          do jo=1,Nso
!             j = jo + (jlat-1)*Nso
!             !
!             Hij(i,j,1) = Htmp(io,jo)
!             !
!          enddo
!       enddo
!       !
!    enddo
!    call eta(ilat,Nlat)
! enddo
! where(abs(Hij)<1.d-6)Hij=zero
! call stop_timer


! print*,Nlat*Nso
! allocate(Evals(Nlat*Nso))
! allocate(rhoH(Nlat*Nso,Nlat*Nso))
! call eigh(Hij(:,:,1),Evals)
! rhoDiag = fermi(Evals,beta)
! rhoH    = matmul(Hij(:,:,1) , matmul(diag(rhoDiag), conjg(transpose(Hij(:,:,1)))) )

! ilat=1
! do io=1,Nso
!    dens(io) = dreal(rhoH(io+(ilat-1)*Nso,io+(ilat-1)*Nso))
! enddo
! write(*,"(A,10F14.9)")"Occupations =",(dens(io),io=1,Nso),sum(dens)
! open(10,file="Robservables.nint")
! write(10,"(10F20.12)")(dens(io),io=1,Nso),sum(dens)
! close(10)


