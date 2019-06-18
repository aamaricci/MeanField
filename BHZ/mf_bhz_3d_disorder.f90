program bhz_3d
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
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable       :: Hk,Hij
  complex(8),dimension(:,:,:,:),allocatable     :: Hlat
  integer                                       :: Iter,MaxIter,Nsuccess=2
  real(8)                                       :: chern,Uloc,Jh,JU,Sz,Tz,Rz,Ntot
  real(8)                                       :: ez,mh,rh,lambda,delta,lz
  real(8)                                       :: xmu,beta,eps,Ekin,Eloc
  real(8)                                       :: n(Nso),arg,dens(Nso),wmix,it_error,sb_field
  complex(8)                                    :: Hloc(Nso,Nso),Htmp(Nso,Nso)
  complex(8),dimension(:,:,:,:,:),allocatable   :: Gmats,Greal
  complex(8),dimension(:,:,:,:,:,:),allocatable :: GLmats,GLreal
  character(len=20)                             :: file
  logical                                       :: iexist,converged,withgf
  complex(8),dimension(Nso,Nso)                 :: Gamma5,GammaX,GammaY,GammaZ,GammaS
  real(8),dimension(3)                          :: vecK,vecRi,vecRj
  complex(8),dimension(:,:),allocatable         :: rhoH
  real(8),dimension(Nso)                      :: params,params_prev,global_params



  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(L,"L","inputBHZ.conf",default=2048)
  call parse_input_variable(Uloc,"ULOC","inputBHZ.conf",default=1d0)
  call parse_input_variable(Jh,"Jh","inputBHZ.conf",default=0.125d0)
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
  gammaX=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gammaY=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gammaZ=kron_pauli( pauli_tau_x, pauli_sigma_x)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)
  gammaS=kron_pauli( pauli_tau_z, pauli_sigma_0)


  !< start using TB_procedures:
  !< 1st set up the direct and momentum space lattice basis
  call TB_set_ei([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])


  !BUILD H(I,J) from Fourier Transforming the H(k)
  !< first step we define the grids in K and R space
  allocate(Kgrid(Nktot,3))
  call TB_build_kgrid([Nkx,Nky,Nkz],Kgrid)
  !
  allocate(Rgrid(Nlat,3))
  call TB_build_Rgrid([Nx,Ny,Nz],Rgrid)


  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  allocate(Links(6,3))
  Links(1,:) = [1,0,0]
  Links(2,:) = [0,1,0]
  Links(3,:) = [-1,0,0]
  Links(4,:) = [0,-1,0]
  Links(5,:) = [0,0,1]
  Links(6,:) = [0,0,-1]
  call TB_build_model(Hlat,ts_model,Nso,[Nkx,Nky,Nkz],Links,pbc=.true.)


  Hij = zero
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

  call eigh(Hij(:,:,1),Evals)
  rhoDiag = fermi(Evals,beta)
  rhoH    = matmul(Hij(:,:,1) , matmul(diag(rhoDiag), conjg(transpose(Hij(:,:,1)))) )

  ilat=1
  do io=1,Nso
     dens(io) = dreal(rhoH(io+(ilat-1)*Nso,io+(ilat-1)*Nso))
  enddo
  write(*,"(A,10F14.9)")"TOccupations =",(dens(io),io=1,Nso),sum(dens)
  open(10,file="Tobservables.nint")
  write(10,"(10F20.12)")(dens(io),io=1,Nso),sum(dens)
  close(10)





  !< BUILD THE LOCAL GF
  if(withgf)then
     allocate(GLmats(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(GLreal(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_matsubara(Hij,[1d0],GLmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hij,[1d0],GLreal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
     call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb+(ispin-1)*Norb
           n(io) = fft_get_density(GLmats(1,ispin,ispin,iorb,iorb,:),beta)
        enddo
     enddo
     open(10,file="observables.nint")
     write(10,"(10F20.12)")(n(io),io=1,Nso),sum(n)
     close(10)
     write(*,"(A,10F14.9)")"Occupations =",(n(io),io=1,Nso),sum(n)
  endif




contains



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



end program bhz_3d







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


