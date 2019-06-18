program mf_hm_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                             :: Norb=1,Nspin=2,Nlat=2,Nlso=Nlat*Nspin*Norb
  integer                                       :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                       :: Nky,Nx,Ny
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
  real(8)                                       :: ts
  real(8)                                       :: xmu,beta,eps
  real(8)                                       :: wmix,it_error,sb_field
  complex(8),dimension(:,:,:,:,:,:),allocatable :: Gmats,Greal
  character(len=20)                             :: file
  logical                                       :: iexist,converged,withgf
  real(8),dimension(1)                          :: params,params_prev,global_params


  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(L,"L","inputHM.conf",default=2048)
  call parse_input_variable(ts,"TS","inputHM.conf",default=1d0)
  call parse_input_variable(Uloc,"ULOC","inputHM.conf",default=1d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmix,"WMIX","inputHM.conf",default=0.5d0)
  call parse_input_variable(sb_field,"SB_FIELD","inputHM.conf",default=0.d0)
  call parse_input_variable(it_error,"IT_ERROR","inputHM.conf",default=1d-5)
  call parse_input_variable(maxiter,"MAXITER","inputHM.conf",default=100)
  call parse_input_variable(withgf,"WITHGF","inputHM.conf",default=.false.)
  !
  call save_input_file("inputHM.conf")
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
  !>Reciprocal lattice basis vector  
  call TB_set_bk(pi*[1,-1],2*pi*[0,1])

  allocate(kgrid(Nktot,2))      !Nktot=# tot kpoints, 2= 2D
  call TB_build_kgrid([Nkx,Nky],kgrid)


  global_params = 0d0
  params        = [sb_field]   ![Sz,Tz,Rz,N]
  !
  inquire(file="params.restart",exist=iexist)
  if(iexist)call read_array("params.restart",params)
  params(1)=params(1)+sb_field

  open(100,file="sz.dat")
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
  call save_array("params.out",params)
  global_params = params
  close(100)
  close(101)
  close(102)

  if(withgf)then
     !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
     write(*,*) "Using Nk_total="//txtfy(Nktot)
     allocate(Hk(Nlso,Nlso,Nktot))
     allocate(Wtk(Nktot))
     call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nky],wdos=.false.)
     Wtk = 1d0/Nktot
     allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,L))
     allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
     call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
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
     call TB_Solve_model(Hk_model,Nlso,kpath,Nkpath,&
          colors_name=[red1,blue1,red2,blue2],&
          points_name=[character(len=20) :: 'X', 'G', 'M', 'X', 'G'],&
          file="Eigenband.nint")
  endif



contains

  subroutine solve_MF_bhz(iter,a)
    real(8),dimension(1),intent(inout) :: a
    real(8),dimension(2)                  :: kvec
    complex(8),dimension(Nlso,Nlso)       :: Hk
    real(8),dimension(Nlso)               :: Ek,rhoDiag
    real(8),dimension(Nlso,Nlso)          :: rhoHk
    real(8)                               :: dens(Nlat,Nspin)
    integer                               :: ik,iter,iorb,ispin
    !
    dens = 0d0
    do ik=1,Nktot
       kvec = kgrid(ik,:)                               ![kx,ky]
       Hk   = hk_model(kvec,Nlso) + mf_hk_correction(a) !H(k)
       !
       call eigh(Hk,Ek)       !diag Hk --> Ek
       !
       rhoDiag = fermi(Ek,beta)
       rhoHk   = matmul( Hk, matmul(diag(rhoDiag),conjg(transpose(Hk))) )
       !
       do ilat=1,Nlat
          do ispin=1,Nspin
             io = ispin+(ilat-1)*Nspin
             dens(ilat,ispin) = dens(ilat,ispin) + rhoHk(io,io)
          enddo
       enddo
       !
    enddo
    dens = dens/Nktot
    !
    a(1) = 0.5d0*dens(1,1) - 0.5d0*dens(1,2) !N_1up - N_1dw = -(N_2up-N_2dw)
    print*,iter,a(1)
    !
    rewind(100);rewind(102)
    write(100,*)a(1)
    write(102,*)(dens(1,ispin),ispin=1,Nspin),(dens(2,ispin),ispin=1,Nspin)
  end subroutine solve_MF_bhz


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    complex(8),dimension(2,2) :: tts
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    hk=zero
    !
    tts  = -ts*pauli_0
    !
    hk(1:2,3:4) = tts*(1d0 + exp(xi*2*kx) + exp(xi*(kx+ky)) + exp(xi*(kx-ky)))
    hk(3:4,1:2) = tts*(1d0 + exp(-xi*2*kx)+ exp(-xi*(kx+ky))+ exp(-xi*(kx-ky)))
  end function hk_model



  function mf_Hk_correction(a) result(HkMF)
    real(8),dimension(1)          :: a
    complex(8),dimension(Nlso,Nlso) :: HkMF
    HkMF = -Uloc*a(1)/3d0*diag([1d0,-1d0,-1d0,1d0]) 
  end function mf_Hk_correction



end program mf_hm_2d


