program bdg_disorder_mpi
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE MPI
  implicit none
  integer                                         :: Nside,Nlat
  integer                                         :: Nloop,Nblock
  integer                                         :: Lmats,Lreal
  integer                                         :: Nsuccess
  real(8)                                         :: Uloc
  real(8)                                         :: ts
  real(8)                                         :: beta
  real(8)                                         :: Wdis
  real(8)                                         :: wmixing
  real(8)                                         :: nsc
  real(8)                                         :: deltasc
  real(8)                                         :: wini,wfin
  real(8)                                         :: bdg_error

  real(8),allocatable,dimension(:)                :: erandom
  real(8),dimension(:),allocatable                :: nii,pii ![Nlat]
  real(8),dimension(:),allocatable                :: nii_prev,pii_prev ![Nlat]
  logical                                         :: converged,bool,igetgf,boolN,boolP
  integer                                         :: i,is,iloop,nrandom,idum
  real(8),dimension(:),allocatable                :: wm,wr
  complex(8),allocatable,dimension(:,:,:,:)       :: Hk ![2][Nlat][Nlat][1]
  real(8),dimension(:,:),allocatable              :: H0              ![Nlat][Nlat]
  character(len=50)                               :: finput
  integer,dimension(:),allocatable                :: icol,irow
  integer,dimension(:,:),allocatable              :: ij2site
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal ![2][Nlat][1][1][1][1][L]

  logical :: master

  call init_BLACS()
  master = get_master_BLACS()

  ! READ INPUT FILES !
  call parse_cmd_variable(finput,"FINPUT",default="inputBDG.conf")
  call parse_input_variable(Nside,"NSIDE",finput,default=10)
  call parse_input_variable(Nblock,"Nblock",finput,default=4)
  call parse_input_variable(ts,"TS",finput,default=0.5d0)
  call parse_input_variable(uloc,"ULOC",finput,default=-1d0,comment="Values of the local interaction")
  call parse_input_variable(beta,"BETA",finput,default=1000d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(Wdis,"WDIS",finput,default=0.d0)
  call parse_input_variable(idum,"IDUM",finput,default=1234567)
  call parse_input_variable(igetgf,"IGETGF",finput,default=.false.)
  call parse_input_variable(nsc,"NSC",finput,default=0.5d0,comment="Value of the initial density per spin (1/2=half-filling).")
  call parse_input_variable(deltasc,"DELTASC",finput,default=0.02d0,comment="Value of the SC symmetry breaking term.")
  call parse_input_variable(nloop,"NLOOP",finput,default=500,comment="Max number of iterations.")
  call parse_input_variable(Lmats,"LMATS",finput,default=2000,comment="Number of Matsubara frequencies.")
  call parse_input_variable(Lreal,"LREAL",finput,default=2000,comment="Number of real-axis frequencies.")
  call parse_input_variable(wini,"WINI",finput,default=-5.d0,comment="Smallest real-axis frequency")
  call parse_input_variable(wfin,"WFIN",finput,default=5.d0,comment="Largest real-axis frequency")
  call parse_input_variable(bdg_error,"BDG_ERROR",finput,default=0.00001d0,comment="Error threshold for the convergence")
  call parse_input_variable(nsuccess,"NSUCCESS",finput,default=1,comment="Number of successive iterations below threshold for convergence")  
  if(master)call save_input_file(finput)
  if(master)call print_input()

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(0d0,"XMU")
  call add_ctrl_var(wini,"WINI")
  call add_ctrl_var(wfin,"WFIN")
  call add_ctrl_var(0.1d0,"EPS")

  ! SET THE COMPRESSION THRESHOLD TO 1Mb (1024Kb)
  call set_store_size(1024)


  ! GET THE ACTUAL NUMBER OF LATTICE SITES !
  Nlat = Nside*Nside

  ! ALLOCATE ALL THE REQUIRED VARIABLES !
  allocate(nii(Nlat),nii_prev(Nlat))
  allocate(pii(Nlat),pii_prev(Nlat))

  ! ALLOCATE MATSUBARA AND REAL FREQUENCIES !
  allocate(wm(Lmats),wr(Lreal))
  wini  =wini-Wdis
  wfin  =wfin+Wdis
  wr = linspace(wini,wfin,Lreal)
  wm(:)  = pi/beta*real(2*arange(1,Lmats)-1,8)


  ! GET RANDOM ENERGIES
  allocate(erandom(Nlat))
  call mersenne_init(idum)
  call mt_random(erandom)
  erandom=(2.d0*erandom-1.d0)*Wdis/2.d0
  inquire(file='erandom.restart',exist=bool)
  if(bool)then
     if(file_length('erandom.restart')/=Nlat)stop "ed_ahm_disorder error: found erandom.restart with length different from Nlat"
     call read_array('erandom.restart',erandom)
  endif

  ! GET THE TIGHT BINDING HAMILTONIAN FOR THE SQUARE LATTICE 
  allocate(H0(Nlat,Nlat))
  H0 = Htb_square_lattice(Nrow=Nside,Ncol=Nside,ts=ts)
  H0 = H0 + diag(Erandom)


  nii=nsc;
  pii=deltasc;
  inquire(file="nVSisite.restart",exist=boolN)
  if(boolN)call read_array("nVSisite.restart",nii)
  inquire(file="phiVSisite.restart",exist=boolP)
  if(boolP)call read_array("phiVSisite.restart",pii)

  if(igetgf)then
     if(.not.boolN.OR..not.boolP)stop "Can not read the input Nii, Pii in .restart files"
     allocate(Gmats(2,Nlat,1,1,1,1,Lmats),Smats(2,Nlat,1,1,1,1,Lmats))
     allocate(Greal(2,Nlat,1,1,1,1,Lreal),Sreal(2,Nlat,1,1,1,1,Lreal))
     Smats=zero ; Sreal=zero
     Gmats=zero ; Greal=zero
     forall(i=1:Nlat)
        Smats(1,i,1,1,1,1,:)= Uloc*(nii(i)-0.5d0)
        Sreal(1,i,1,1,1,1,:)= Uloc*(nii(i)-0.5d0)
        Smats(2,i,1,1,1,1,:)= Uloc*pii(i)
        Sreal(2,i,1,1,1,1,:)= Uloc*pii(i)
     end forall
     !
     allocate(Hk(2,Nlat,Nlat,1))
     Hk(1,:,:,1)    =     H0
     Hk(2,:,:,1)    =    -H0
     !
     call dmft_gloc_matsubara(Hk,[1d0],Gmats,Smats)
     call dmft_gloc_realaxis(Hk,[1d0],Greal,Sreal)
     if(master)then
        call dmft_print_gf_matsubara(Gmats(1,:,:,:,:,:,:),"Gmats",iprint=4)
        call dmft_print_gf_matsubara(Gmats(2,:,:,:,:,:,:),"Fmats",iprint=4)
        call dmft_print_gf_realaxis(Greal(1,:,:,:,:,:,:),"Greal",iprint=4)
        call dmft_print_gf_realaxis(Greal(2,:,:,:,:,:,:),"Freal",iprint=4)
     endif
  endif


  iloop=0;
  converged=.false.;
  !+-------------------------------------+!
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"BdG-loop")
     !
     call BdG_Solve(nii,pii,H0)
     !     
     nii_prev = nii
     pii_prev = pii
     if(iloop>1)then
        nii = wmixing*nii + (1d0-wmixing)*nii_prev
        pii = wmixing*pii + (1d0-wmixing)*pii_prev
     endif
     !
     if(master)then
        converged = check_convergence_local(pii,bdg_error,nsuccess,nloop,file="error.err")
        call print_sc_out(converged)
     endif
     call bcast_MPI(MPI_COMM_WORLD,converged)

     if(master)call end_loop()
  enddo
  !+-------------------------------------+!

  if(master)then
     call save_array("nVSisite.restart",nii)
     call save_array("phiVSisite.restart",pii)
     call save_array('erandom.restart',erandom)
  endif


  call finalize_BLACS()

contains





  subroutine BdG_Solve(nii,pii,H)
    real(8),dimension(Nlat),intent(inout)   :: nii
    real(8),dimension(Nlat),intent(inout)   :: pii
    real(8),dimension(Nlat,Nlat),intent(in) :: H
    real(8),dimension(Nlat)                 :: Sigma_HFB
    real(8),dimension(Nlat)                 :: Self_HFB
    real(8),dimension(2*Nlat,2*Nlat)        :: Hnambu
    real(8),dimension(2*Nlat)               :: Enambu
    real(8),dimension(2*Nlat)               :: RhoDiag  !diagonal density matrix 
    real(8),dimension(2*Nlat,2*Nlat)        :: RhoNambu !Nambu density matrix
    integer                                 :: ilat
    Sigma_HFB(:) =  Uloc*(nii-0.5d0)
    Self_HFB(:)  =  Uloc*pii(:)
    Hnambu(1:Nlat,1:Nlat)               =  H + diag(Sigma_HFB)
    Hnambu(1:Nlat,Nlat+1:2*Nlat)        =    + diag(Self_HFB)
    Hnambu(Nlat+1:2*Nlat,1:Nlat)        =    + diag(Self_HFB)
    Hnambu(Nlat+1:2*Nlat,Nlat+1:2*Nlat) = -H - diag(Sigma_HFB)
    !
    call p_eigh(HNambu,ENambu,Nblock)
    RhoDiag  = fermi(ENambu,beta)
    RhoNambu = (HNambu.px.(diag(RhoDiag))).px.(transpose(HNambu))
    forall(ilat=1:Nlat)
       nii(ilat) = RhoNambu(ilat,ilat)
       pii(ilat) = RhoNambu(ilat,ilat+Nlat)
    end forall
    return
  end subroutine BdG_Solve










  ! !******************************************************************
  ! !******************************************************************




  subroutine print_sc_out(converged)
    integer                        :: i,j,is,row,col,unit
    real(8)                        :: nimp,phi,ccdw
    real(8),dimension(Nlat)        :: cdwii
    real(8),dimension(Nside,Nside) :: nij,pij
    real(8),dimension(Nside)       :: grid_x,grid_y
    real(8)                        :: mean,sdev,var,skew,kurt
    real(8),dimension(2,Nlat)      :: data_covariance
    real(8),dimension(2,2)         :: covariance_nd
    real(8),dimension(2)           :: data_mean,data_sdev,Eout
    logical                        :: converged
    character(len=50)              :: suffix
    nii=nii*2
    suffix=".bdg"
    !Get CDW "order parameter"
    do is=1,Nlat
       cdwii(is) = (-1.d0)**(irow(is)+icol(is))*(nii(is)-1.d0)
    enddo
    nimp = sum(nii(:))/dble(Nlat)
    phi  = sum(pii(:))/dble(Nlat)
    ccdw = sum(cdwii)/dble(Nlat)
    print*,"<nimp>  =",nimp
    print*,"<phi>   =",phi
    print*,"<ccdw>  =",ccdw
    open(unit=free_unit(unit),file="n_phi_cdwVSiloop.bdg",access='append')
    write(unit,*)iloop,nimp,phi,ccdw
    close(unit)
    call save_array("nVSisite"//trim(suffix),nii(:))
    call save_array("phiVSisite"//trim(suffix),pii(:))
    call save_array("cdwVSisite"//trim(suffix),cdwii)
    !
    !WHEN CONVERGED IS ACHIEVED PLOT ADDITIONAL INFORMATION:
    if(converged)then
       do row=1,Nside
          grid_x(row)=row
          grid_y(row)=row
          do col=1,Nside
             i            = ij2site(row,col)
             nij(row,col) = nii(i)
             pij(row,col) = pii(i)
          enddo
       enddo
       call splot3d("3d_nVSij"//trim(suffix),grid_x,grid_y,nij)
       call splot3d("3d_phiVSij"//trim(suffix),grid_x,grid_y,pij)
       call get_moments(nii(:),mean,sdev,var,skew,kurt)
       open(unit=free_unit(unit),file="statistics.n"//trim(suffix))
       write(unit,*)mean,sdev,var,skew,kurt
       close(unit)
       data_mean(1)=mean
       data_sdev(1)=sdev
       !
       call get_moments(pii(:),mean,sdev,var,skew,kurt)
       open(unit=free_unit(unit),file="statistics.phi"//trim(suffix))
       write(unit,*)mean,sdev,var,skew,kurt
       close(unit)
       data_mean(2)=mean
       data_sdev(2)=sdev
       !
       call get_moments(cdwii,mean,sdev,var,skew,kurt)
       open(unit=free_unit(unit),file="statistics.cdwn"//trim(suffix))
       write(unit,*)mean,sdev,var,skew,kurt
       close(unit)
       !
       data_covariance(1,:)=nii(:)
       data_covariance(2,:)=pii(:)
       covariance_nd = get_covariance(data_covariance,data_mean)
       open(10,file="covariance_n.phi"//trim(suffix))
       do i=1,2
          write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
       enddo
       close(10)
       !
       covariance_nd=0d0
       do i=1,2
          do j=1,2
             if(data_sdev(i)/=0d0.AND.data_sdev(j)/=0d0)then
                covariance_nd(i,j) = covariance_nd(i,j)/(data_sdev(i)*data_sdev(j))
             endif
          enddo
       enddo
       open(10,file="correlation_n.phi"//trim(suffix))
       do i=1,2
          write(10,"(2f24.12)")(covariance_nd(i,j),j=1,2)
       enddo
       close(10)
    end if
    nii=nii/2d0
  end subroutine print_sc_out





  ! !******************************************************************
  ! !******************************************************************



  function Htb_square_lattice(Nrow,Ncol,pbc_row,pbc_col,ts) result(H0)
    integer                                :: Nrow
    integer                                :: Ncol
    logical,optional                       :: pbc_row,pbc_col
    logical                                :: pbc_row_,pbc_col_
    real(8),optional                       :: ts
    real(8)                                :: ts_
    real(8),dimension(Nrow*Ncol,Nrow*Ncol) :: H0
    integer                                :: i,jj,row,col,link(4),j
    integer                                :: unit
    !
    pbc_row_=.true. ; if(present(pbc_row)) pbc_row_=pbc_row
    pbc_col_=.true. ; if(present(pbc_col)) pbc_col_=pbc_col
    ts_=0.5d0;if(present(ts))ts_=ts
    !
    H0 = 0.d0
    unit=free_unit()
    !+- 2D LATTICE (NROW x NCOL) -+!
    if(Nlat /= Nrow*Ncol) stop "get_lattice_hamiltonian error: Nlat != Nrow*Ncol"
    !THESE ARE STILL GLOBAL VARIABLES...
    allocate(icol(Nlat),irow(Nlat))
    allocate(ij2site(Nrow,Ncol))
    do row=0,Nrow-1
       do col=0,Ncol-1
          i=col+ 1 + row*Ncol
          !
          irow(i)=row+1
          icol(i)=col+1
          ij2site(row+1,col+1)=i
          !
          !right hop
          link(1)= i + 1     
          if((col+1)==Ncol) then
             link(1)=0  
             if(pbc_col_)link(1)=1+row*Ncol  
          end if
          !left  hop
          link(3)= i - 1    
          if((col-1)<0)     then
             link(3)=0  
             if(pbc_col_)link(3)=Ncol+row*Ncol
          end if
          !up    hop
          link(2)= i + Ncol 
          if((row+1)==Nrow) then
             link(2)=0  
             if(pbc_row_)link(2)=col+1
          end if
          !down  hop
          link(4)= i - Ncol 
          if((row-1)<0)     then
             link(4)=0  
             if(pbc_row_)link(4)=col+1+(Nrow-1)*Ncol
          end if
          !
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=-ts_ !! ts must be negative.
          enddo
          !
       enddo
    enddo
    open(unit,file='Htb_square_lattice.bdg')
    do i=1,Nlat
       write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,Nlat)
    enddo
    close(unit)
  end function Htb_square_lattice



end program bdg_disorder_mpi
