program bhz_3d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                             :: L=2048,Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                       :: Nk,Nktot,Nkx,Nky,Nkz
  integer                                       :: Nlat,Nx,Ny,Nz
  integer                                       :: Nkpath,Npts,z2(4)
  integer                                       :: i,j,k,ik
  integer                                       :: io,jo,ispin
  integer                                       :: iorb,ilat,jorb,jlat
  integer                                       :: ix,iy,iz
  real(8)                                       :: kx,ky,kz
  real(8),dimension(:),allocatable              :: Wtk,Evals,rhoDiag
  real(8),dimension(:,:),allocatable            :: kpath,ktrims,Kgrid,Rgrid
  complex(8),dimension(:,:,:),allocatable       :: Hk,Hij
  complex(8),dimension(:,:,:,:),allocatable     :: Hlat
  real(8)                                       :: ez,mh,rh,lambda,delta,lz
  real(8)                                       :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                          :: wm,wr
  real(8)                                       :: n(Nso),arg,dens(Nso)
  complex(8)                                    :: w,Hloc(Nso,Nso)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,L) :: Gmats,Greal
  character(len=20)                             :: file,nkstring
  logical                                       :: iexist,ibool,dcflag,iener
  complex(8),dimension(Nso,Nso)                 :: Gamma5,GammaX,GammaY,GammaZ
  real(8),dimension(3)                          :: vecK,vecRi,vecRj
  complex(8),dimension(:,:),allocatable         :: rhoH

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(ez,"ez","inputBHZ.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call save_input_file("inputBHZ.conf")
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


  !< start using TB_procedures:
  !< 1st set up the direct and momentum space lattice basis
  call TB_set_ei([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])



  write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nky,Nkz],wdos=.false.)
  Wtk = 1d0/Nktot

  call TB_write_hk(Hk,trim(file),Nso,&
       Nd=Norb,Np=1,Nineq=1,&
       Nkvec=[Nkx,Nky,Nkz])


  !GET LOCAL PART OF THE HAMILTONIAN

  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)


  !< SOLVE ALONG A PATH IN THE BZ.
  Npts=8
  allocate(kpath(Npts,3))
  kpath(1,:)=[0,1,0]!X
  kpath(2,:)=[0,0,0]!G
  kpath(3,:)=[1,1,0]!M
  kpath(4,:)=[1,1,1]!R
  kpath(5,:)=[0,0,1]!Z
  kpath(6,:)=[1,0,1]!A
  kpath(7,:)=[0,0,0]!G
  kpath(8,:)=[0,0,1]!Z
  kpath     =kpath*pi
  call TB_solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=10) :: "X","G","M","R","Z","A","G","Z"],&
       file="Eigenband.nint")



  !TO BE REMOVED:
  !< BUILD THE LOCAL GF
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"Gloc",1)
  call dmft_print_gf_realaxis(Greal,"Gloc",1)
  !
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb+(ispin-1)*Norb
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo
  open(10,file="observables.nint")
  write(10,"(10F20.12)")(n(io),io=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(n(io),io=1,Nso),sum(n)



  !ADD HERE THE FOURIER TRANSORM TO H(I,J)

  !< first step we define the grids in K and R space

  allocate(Kgrid(Nktot,3))
  call TB_build_kgrid([Nkx,Nky,Nkz],Kgrid)
  allocate(Rgrid(Nlat,3))
  call TB_build_Rgrid([Nx,Ny,Nz],Rgrid)




  !Build up the real-space Hamiltonian thru FT:"
  print*,"REAL SPACE:       H_bhz(i,j)=FT^-1(H_bhz(kx,ky))"

  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  Hlat=zero
  call start_timer
  do ik=1,Nktot
     vecK = Kgrid(ik,:)
     !
     do ilat=1,Nlat
        vecRi = Rgrid(ilat,:)
        do jlat=1,Nlat
           vecRj = Rgrid(jlat,:)
           !
           arg=dot_product(vecK,vecRj-vecRi)
           !
           Hlat(:,:,ilat,jlat)= Hlat(:,:,ilat,jlat) + exp(xi*arg)*Hk(:,:,ik)/Nktot!hk_model(vecK,Nso)/Nktot
        enddo
     enddo
     call eta(ik,Nktot)
  enddo
  where(abs(Hlat)<1.d-6)Hlat=zero
  call stop_timer



  allocate(Hij(Nlat*Nso,Nlat*Nso,1))
  do ilat=1,Nlat
     do io=1,Nso
        i = io + (ilat-1)*Nso
        do jlat=1,Nlat
           do jo=1,Nso
              j = jo + (jlat-1)*Nso
              !
              Hij(i,j,1) = Hlat(io,jo,ilat,jlat)
              !
           enddo
        enddo
     enddo
  enddo

  allocate(Evals(Nlat*Nso))
  call eigh(Hij(:,:,1),Evals)

  allocate(rhoH(Nlat*Nso,Nlat*Nso))

  rhoDiag = fermi(Evals,beta)
  rhoH    = matmul(Hij(:,:,1) , matmul(diag(rhoDiag), conjg(transpose(Hij(:,:,1)))) )

  do io=1,Nso
     dens(io) = dreal(rhoH(io,io))
  enddo
  open(10,file="Robservables.nint")
  write(10,"(10F20.12)")(dens(io),io=1,Nso),sum(dens)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(dens(io),io=1,Nso),sum(dens)



  open(10,file="DensError.nint")
  write(10,"(10F20.12)")(abs(dens(io)-n(io)),io=1,Nso)
  close(10)
  write(*,"(A,10F14.9)")"Delta_dens =",(abs(dens(io)-n(io)),io=1,Nso)


contains




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - cos(kx) - cos(ky) - ez*cos(kz))*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y
    Hk(3:4,3:4) = conjg( &
         (Mh-cos(-kx) - cos(-ky) - ez*cos(-kz))*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y)
    Hk(1:2,3:4) = lambda*sin(kz)*pauli_tau_x
    Hk(3:4,1:2) = lambda*sin(kz)*pauli_tau_x
  end function hk_model



end program bhz_3d


