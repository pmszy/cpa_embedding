c**************************************************************************** 
c***************************************************************************
c***************************************************************************
c***************************************************************************
c***************************************************************************
c***************************************************************************
c***************************************************************************
	module Global
c****************************************************************************
c       This module should appear in all the main blocks of the code.  
c       It replaces the common block and declares and passes all the 
c       global variables and parameters (those globally accessible).
c       Note that many of the vaiables are allocatable.  They are
c       allocated/deallocated in the calls below for 1 and 2 particle
c       variables.
c****************************************************************************
        use mod_tprof
	        implicit none
	        save
            include 'mpif.h'     
c**************************************************************************** 
c       set some parameters.
c**************************************************************************** 
c	Character Parameters
c**************************************************************************** 
	character*4 :: cluster
c**************************************************************************** 
c       Logical parameters.
c****************************************************************************
        logical, parameter :: isw = .false. ! if true, use the inhomogeneous grid else use homogeneous grid
        logical, parameter :: histo_enable = .false. !Set true if you want to calculate histogram else, always false. It is very expensive
        logical, parameter :: binary = .false.
        logical, parameter :: box=.true.
        logical, parameter :: gauss=.false.
        logical, parameter :: boltz=.false.
        logical, parameter ::  ph_symm = .false.
          logical, parameter ::  my = .true. ! switch for SD, case of inverted Gscript & MB
          logical, parameter ::   MB1 = .true.
          logical, parameter ::   Gfav=.true. !use Green function in DCA
          logical, parameter ::   exact=.false. !evaluate all the disorder configurations for binary disorder
          logical, parameter ::   Gscript = .true.!using Gscript to instead of Delta
          logical, parameter ::   Gscript_pole = .false. !pole procedure only worked for nb=1
          logical, parameter ::   init_sigma=.true.
          logical, parameter ::   embedding1=.false.
          logical, parameter ::   embedding2=.false.
          logical, parameter ::   liz_embedding=.true.
          logical, parameter ::   free_standing=.false.
          logical, parameter ::   mb_tmt= .true. ! MB-has to be on
           logical, parameter ::   mb_tmt_ver1= .true. !Hanna's 
           logical, parameter ::   mb_tmt_ver2= .false. !rho_tot_typ basis invariant 
            logical, parameter ::   mb_tmt_ver6= .false. !Yi's new trial
            logical, parameter ::   mb_tmt_ver7= .false. !local ansatz
            logical, parameter ::   mb_tmt_ver8= .false. !matrix version
            logical, parameter ::   mb_tmt_ver9= .false. !matrix version+rho_loc_typ
          !	Interpolation flags. Works only if interpolation is true and case either 1,2, or 3.
        logical, parameter :: interpolation = .false. !Interpolation works only if case is either 1, 2, or 3. Set below in interpolation_type
!        logical, parameter :: spline_interpolation = .false. ! Set true calculate error bar
!        logical, parameter :: star_interpolation = .false. ! Set true calculate error bar
!        logical, parameter :: trilinear_interpolation = .true. ! Set true calculate error bar
c**************************************************************************** 
c       Integer parameters.
c**************************************************************************** 
	integer, parameter :: kind=8
        integer, parameter :: nwnm=1200,ntm=50,neta=4,ndim=3,nl=10,e_div=200
        integer, parameter :: ichim=0   ! if 1, measure chim (expensive)
        integer, parameter :: ichit=0   ! if 1, measure chit
        integer, parameter :: ibetts=0  ! If 1, use Betts Lattice
        !integer, parameter :: tables_type = 3 ! 1=ibetts,2=conventional,3=cubic
	integer, parameter :: lwork=100**2  ! zgetri and zgetrf parameter 
	integer, parameter :: lud=42,iprint=19 
	integer, parameter :: imeas=3000 ! For histogram. Should be larged enough
        integer :: interpolation_type = 2 ! 1=trilinear_interpolation,2=star_interpolation,3=spline_interpolation


c**************************************************************************** 
c       Real parameters.
c**************************************************************************** 
        real(kind), parameter :: tmax=30.d0, dpi=dacos(-1.d0)
        real(kind), parameter ::  epsilon=1.0d-05
        real(kind), parameter ::  depsi=1.0d-06
        real(kind), parameter ::  zeror=0.0d+00
        real(kind), parameter ::  halfr=0.5d+00
        real(kind), parameter ::   oner=1.0d+00
        real(kind), parameter ::   twor=2.0d+00
        real(kind), parameter :: threer=3.0d+00
        real(kind), parameter ::  fourr=4.0d+00
        real(kind), parameter ::     pi=3.141592653589793, hbar=1.0545717260d-34
	    real(kind), parameter :: estart=-1.51d0, eend=1.51d0
	    real(kind), parameter :: afact = 0.2d0 ! Adjust according when using the pole procedure
	    real(kind), parameter :: dim = 3 ! Dimension of the cluster
c**************************************************************************** 
c       Complex parameters.
c****************************************************************************
        complex(kind), parameter :: ii=(0.0d+00,1.0d+00)
        complex(kind), parameter :: eta=(0.0d+00,0.0d-02) ! Broadening parameter
        complex(kind), parameter :: zeroc=(0.d0,0.d0)
        complex(kind), parameter :: onec=(1.d0,0.d0)
c**************************************************************************** 
c       Integer Variables Global
c****************************************************************************
        integer :: iseed,iter,niter,nover,nwn,Nc,nt,run,nrun,nf,
     &		meas,nmeas,nslr,isi,iso,isd,isf,isv,isg,tables_type,
     &		ncw,a(3,3),ngroup,icR0,icK0,icKpi,ntw,ick,icr,
     &		N_sp,ntot,nwsc,ic00,icpipi,igroup,x_pole,
     &		a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z,nproc,threads
        integer :: nb=1, liz, lizsize
c**************************************************************************** 
c       Integer Variables for Broyden
c****************************************************************************
	integer :: vlen,broylen
	integer, allocatable  :: ipiv(:)
c****************************************************************************
c        integer, allocatable  
c****************************************************************************

        integer, allocatable  :: neighbor(:,:),
     &		ickdeg(:),Rc(:,:),ict(:,:,:),
     &		ickmap(:),ick2map(:,:),Epsbar_histo(:),
     &		Rcx(:),Rcy(:),Rcz(:),
     &		icrdiff(:,:),ickdiff(:,:),
     &      ickplus(:,:),ickequ(:,:),icrequ(:,:),ipvt(:),ipvt1(:),
     &      ipvt2(:),ipvtn(:),ipvtnb(:),lizmap(:,:),ipvtliz(:),
     &      R_liz(:,:)   
		   
c**************************************************************************** 
c        Real Variables
c****************************************************************************
        real(kind) :: ed,tprime,tperp,V0,dfac,ifac,test,testl,estep,
     &          delta,wc,gvector(3,3),Epsd,g1x,g1y,g1z,g2x,g2y,g2z,g3x,
     &		g3y,g3z,origx,origy,origz,wtime,ned,aa(2),ab,tz
c**************************************************************************** 
c        Real Variables for Broyden
c****************************************************************************
        real(kind)  :: alpha,rms
c**************************************************************************** 
c        Real Allocatables
c***************************************************************************
        real(kind), allocatable  :: V(:),kt(:,:),Kc(:,:),Kcz(:),
     &		Epsbar(:),Kcx(:),Kcy(:),wn(:),dwn(:),weta(:),t(:),
     &		min_Eps(:),max_Eps(:),dEps(:),
     &		p_rho(:,:),PDOS(:,:),p_rho_GcKf(:,:)
c MB     
       real(kind), allocatable  :: V1(:),V2(:),V12(:),Vnb(:,:),Vnbc(:,:,:)
       real(kind)  :: V1a,V2a,V12a

c MPI     
       real(kind), allocatable :: logrholoc_f_nb_mpi(:,:,:), 
     &  rhof_nb_mpi(:,:,:,:), rholocf_av_nb_mpi(:,:,:),logrhotot_f_mpi(:)
       complex(kind), allocatable :: Gcf_nb_mpi(:,:,:,:),GcKf_nb_mpi(:,:,:,:)
       
        integer:: nprocs,ierr,myrank

c**************************************************************************** 
c       Complex Variables
c**************************************************************************** 
c       Complex Allocatables
c****************************************************************************
        complex(kind), allocatable :: G0t(:),G0w(:),four(:,:),work(:),
     &		chi(:,:),FTCoefs_K_to_R(:,:),FTCoefs_R_to_K(:,:)

c For odd case

       complex(kind),allocatable : :Gamma_old_nb(:,:,:,:),
     & Gamma_new_nb(:,:,:,:), sigma_nb(:,:,:,:),
     & GcKf_nb(:,:,:,:), GcKfsc_nb(:,:,:,:),GcRfsc_nb(:,:,:,:,:),
     & Gcinv_nb(:,:,:,:), Gdff_nb(:,:,:,:), GdKK_nb(:,:,:,:),
     & Gc_nb(:,:,:,:),Gc_temp(:,:,:,:),Gloc_inv_nb(:,:,:,:),
     & Gamma_rs_nb(:,:,:,:,:),Gcv_nb(:,:,:,:,:),
     & Epsb_rs_nb(:,:,:,:),Epsb_nb(:,:,:),
     & Gc_meas_nb(:,:,:,:),Gcrun_nb(:,:,:,:),Gcf_nb(:,:,:,:),
     & Gamma_nb(:,:,:,:),Gamma1_nb(:,:,:,:),GcKfsc_nb1(:,:,:,:),
     & GcKfsc_nb2(:,:,:,:),del_tpmb_nb(:,:,:,:,:),V_tpmb_nb(:,:,:,:),
     & G_tpmb_nb(:,:,:,:),G_inmb_nb(:,:,:)

       real(kind), allocatable  :: rho_meas_nb(:,:,:,:),
     & rhorun_nb(:,:,:,:), rhof_nb(:,:,:,:),
     & logrholoc_meas_nb(:,:,:),logrholoc_run_nb(:,:,:),
     & logrholoc_f_nb(:,:,:), rholoc_av_meas_nb(:,:,:),
     & rholoc_av_run_nb(:,:,:),rholocf_av_nb(:,:,:),w_tpmb_nb(:,:,:,:),
     & logrhotot_meas(:),logrhotot_run(:),logrhotot_f(:)
 
       real(kind)  :: ca

         real(kind), allocatable  :: pv(:)
          real(kind)  :: prod
          
!for Liz_embedding
       complex(kind),allocatable ::Gc_liz(:,:,:),Gcbar_liz(:,:,:),
     & Gsc_liz(:,:,:),Gc_liz_meas(:,:,:) ,Gc_liz_run(:,:,:),
     & Gc_liz_f(:,:,:),Gc_liz_mpi(:,:,:),sigma_l(:)        
!for GIPR2
       real(kind) :: GIPR2, GIPR2_meas, GIPR2_run, GIPR2_f, GIPR2_typ,
     & logGIPR2_meas, logGIPR2_run, logGIPR2_f, logGIPR2_mpi, GIPR2_mpi

c**************************************************************************** 
c        Real Variables for spectra
c****************************************************************************
        real(8), allocatable :: g1_sp(:),g2_sp(:),Sig_sp_r(:,:,:),
     &                          Sig_sp_i(:,:,:),Derivs_sp_r(:,:,:),
     &                          Derivs_sp_i(:,:,:)
     
        

c****************************************************************************
c       Define some variables.
c****************************************************************************
c       64 bin precision (kind 8) is explicitly used unless otherwise state
c	throughout this code. To use it on 32 bit system, you have to use kind = 4
c**************************************************************************** 
c       Whenever possible, the following is used to define the
c       nomeclature of subscripts
c**************************************************************************** 
c
c       Gc              cluster greens function
c       K               momentum
c       s               space
c       f               frequency
c       t               euclidean time
c       sc              script (i.e. site excluded)
cc**************************************************************************** 
c       Table of Variables
cc**************************************************************************** 
c       ndim    	Number of spatial dimensions
c       ngroup  	Number of group operations which leave the lattice ivariant
c       nwnm    	Maximum number of frequencies
c       nwn     	number of frequencies used
c       Ncm     	Maximum size of the cluster
c       Nc      	size of the cluster
c       iseed   	seed integer the random number generator
c       iter    	present iteration number
c       niter   	total number of iterations to be done
c       run     	number of runs to be executed per iteration
c       meas    	number of measurements made per run
c       ntr     	number of sweeps before G is recalculated from {s(l)}
c       nslr    	number of sweeps sinse the la	st such recalculation
c       V(p)    	Disorder potential at site p
c       isi     	switch, if isi>=1, increase the # of runs by ifac each iter.
c       isd     	switch, if isd=1 use damping (dfac) in calc. of sigma
c	isf		switch, if isf=1 initalize with an externally provided s(i)
c	isv     	switch, if isv=1 (-1) binary (continuous) disorder 
c       pi      	3.14159
c       ed      	bare orbital energy (measured relative to the chem. pot.)
c       tprime 		next-near-neighbor hopping
c	test		convergence criteria
c       testl   	convergence criteria from last iteration
c       dfac    	damping factor (for sigma, see also isd)
c       ifac    	see isi
c       Kcx     	lookup table for Kcx values
c       Kcy     	lookup table for Kcy values
c       wn      	lookup table of frequencies
c       gf      	all runs accumulator for the single-particle green function
c       gfs     	all runs accumulator for the square single-particle green
c	gft		The damped single-particle green function (only used in put)
c       G       	green function
c       GcKf    	G(k,w) on the cluster
c       GcKfsc  	cluster-excluded G_sc(k,w)
c       Gcsf    	Cluster G(r,w)
c       chi     	temperary array 
c       ii      	sqrt(-1)
c       sigma   	self energy
c       data    	temporary array. All data_ are temporary arrays. E.g., data_rho_local.
c	p_rho		partial density of states
c	G_typ		single run accumulator for the single-particle for the typical green function
c	rho_typ		The K-resolved typical density of states
c	rho_meas	The K-resolved algebraic density of states
c	Gamma		The hybridization rate
c	drhof		Error in ADOS
c	drho_log	Error in TDOS
c	d_gamma		Error in gamma
c	_temp		These are all temporary arrays
c	Gcinv		Inverse of Gscript
c	Ginv		Inverse of G 
c	contains

        
	end module Global
c************************************************************
        subroutine allocate_arrays
c       allocate some arrays to the desired size.
c************************************************************
        use Global
        implicit none
c        include 'mpif.h'
c************************************************************
        integer info,infot,icount
         

        infot=0 
        icount=0
        allocate(ickdeg(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(neighbor(6,3),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ickmap(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ick2map(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(icrequ(Nc,48),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(icrdiff(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ickdiff(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ickplus(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ickequ(Nc,48),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(Rcx(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(Rcy(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ipvt(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(ipvt2(2*Nc),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(ipvt1(2),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(Rc(ndim,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(Rcz(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(V(Nc),stat=info)
        infot=infot+info
        icount=icount+1
	    
        allocate(kt(ndim,(2*nover)**ndim),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(Kc(ndim,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (Kcx(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (Kcy(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (Kcz(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (Epsbar(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (Epsbar_histo(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(wn(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        allocate(dwn(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (weta(neta),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (t(ntm),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (min_Eps(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (max_Eps(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (dEps(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (p_rho(e_div,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (G0t(ntm),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (G0w(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (four(-nwn:nwn,ntm),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (chi(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (FTCoefs_K_to_R(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (FTCoefs_R_to_K(Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (work(lwork),stat=info)
        infot=infot+info
        icount=icount+1
        allocate (p_rho_GcKf(Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
c        write(*,*) 'mod global'


c ----        
        
        
        allocate (pv(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
c       MB        
        allocate(V1(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate(V2(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate(V12(Nc),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate(ict(-150:150,-150:150,-150:150),stat=info)
        infot=infot+info
        icount=icount+1

!MB1
        allocate (Gamma_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gamma1_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (GcKfsc_nb1(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (GcKfsc_nb2(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        
        allocate (Gamma_old_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gamma_new_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (sigma_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (GcKf_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (GcKfsc_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (GcRfsc_nb(nb,nb,Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gdff_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (GdKK_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (Gc_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gc_temp(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gloc_inv_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (rho_meas_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (rhorun_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (rhof_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (Gcinv_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (Gamma_rs_nb(nb,nb,Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
         
        allocate (Epsb_rs_nb(nb,nb,Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1 
        
        allocate (logrholoc_meas_nb(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (logrholoc_run_nb(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (logrholoc_f_nb(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (rholoc_av_meas_nb(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (rholoc_av_run_nb(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (rholocf_av_nb(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (Vnb(nb,Nc),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (Vnbc(nb,Nc,meas),stat=info)
        infot=infot+info
        icount=icount+1 

        allocate (Gcv_nb(nb,nb,Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(ipvtn(nb*Nc),stat=info)
        infot=infot+info
        icount=icount+1
         
        allocate(ipvtnb(nb),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(Epsb_nb(nb,nb,Nc),stat=info)
        infot=infot+info
        icount=icount+1
       
        allocate(Gc_meas_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(Gcrun_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(Gcf_nb(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (rhof_nb_mpi(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (logrholoc_f_nb_mpi(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1
        
        allocate (rholocf_av_nb_mpi(nb,nb,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gcf_nb_mpi(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (GcKf_nb_mpi(nb,nb,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(del_tpmb_nb(nb,nb,Nc,Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(G_tpmb_nb(nb,nb,Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1
         
        allocate(G_inmb_nb(nb*Nc,nb*Nc,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(w_tpmb_nb(nb,nb,Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(V_tpmb_nb(nb,nb,Nc,Nc),stat=info)
        infot=infot+info
        icount=icount+1 
        
        allocate (logrhotot_meas(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (logrhotot_run(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (logrhotot_f(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (logrhotot_f_mpi(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

!for Liz_embedding
        allocate (lizmap(Nc,lizsize),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gc_liz(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gc_liz_meas(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gc_liz_run(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gc_liz_f(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gc_liz_mpi(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gcbar_liz(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate (Gsc_liz(lizsize,lizsize,-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(ipvtliz(lizsize),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(R_liz(lizsize,ndim),stat=info)
        infot=infot+info
        icount=icount+1

        allocate(sigma_l(-nwn:nwn),stat=info)
        infot=infot+info
        icount=icount+1

c       write(*,*) 'mod global5'



        if(infot.ne.0.and.myrank.eq.0) then
          write(lud,*) 'allocate_arrays',infot,' failures at step',icount
          flush(lud)
          stop
        end if

        return
        end subroutine allocate_arrays

c************************************************************
        subroutine deallocate_arrays
c       deallocate the arrays allocated in subroutine allocate_arrays
c************************************************************
        use Global
        implicit none
c        include 'mpif.h'
c************************************************************
        deallocate(ickdeg)
        deallocate(ickmap)
        deallocate(ick2map)
        deallocate(icrequ)
        deallocate(icrdiff)
        deallocate(ickdiff)
        deallocate(ickplus)
        deallocate(ickequ)
        deallocate(Rcx)
        deallocate(Rcy)
        deallocate(Rcz)
        deallocate(V)
        deallocate(kt)
        deallocate(Kc)
        deallocate(Kcx)
        deallocate(Kcy)
        deallocate(Kcz)
        deallocate(Epsbar)
        deallocate(Epsbar_histo)
        deallocate(wn)
        deallocate(dwn)
        deallocate(weta)
        deallocate(ipvt)
        deallocate(work)
        deallocate(t)
        deallocate(min_Eps)
        deallocate(max_Eps)
        deallocate(dEps)
        deallocate(Rc)
        deallocate(p_rho)
        deallocate(G0t)
        deallocate(G0w)
        deallocate(four)
	    deallocate(neighbor)
        deallocate(chi)
        deallocate(FTCoefs_K_to_R)
        deallocate(FTCoefs_R_to_K)
        deallocate(p_rho_GcKf)
	    deallocate(ipvt1)
	    deallocate(ict)
        deallocate(V1)
        deallocate(V2)
        deallocate(V12)
        deallocate(Gamma_nb)
        deallocate(Gamma1_nb)
        deallocate(GcKfsc_nb1)
        deallocate(GcKfsc_nb2)
        deallocate(Gamma_old_nb)
        deallocate(Gamma_new_nb) 
        deallocate(sigma_nb)
        deallocate(GcKf_nb)
        deallocate(GcKfsc_nb) 
        deallocate(GcRfsc_nb) 
        deallocate(Gdff_nb)
        deallocate(GdKK_nb)
        deallocate(Gc_nb)
        deallocate(Gloc_inv_nb)
        deallocate(rho_meas_nb)
        deallocate(rhorun_nb)
        deallocate(rhof_nb)
        deallocate(Gcinv_nb)
        deallocate(Gamma_rs_nb)
        deallocate(Epsb_rs_nb)
        deallocate(logrholoc_meas_nb)
        deallocate(logrholoc_run_nb)
        deallocate(logrholoc_f_nb)
        deallocate(rholoc_av_meas_nb)
        deallocate(rholoc_av_run_nb)
        deallocate(rholocf_av_nb)
        deallocate(Vnb)
        deallocate(Vnbc)
        deallocate(Gcv_nb)
        deallocate(ipvtn)
        deallocate(ipvtnb)
        deallocate(Gc_temp)
        deallocate(Epsb_nb)
        deallocate(Gc_meas_nb)
        deallocate(Gcrun_nb)
        deallocate(Gcf_nb)
        deallocate(rhof_nb_mpi)
	    deallocate(logrholoc_f_nb_mpi)
	    deallocate(rholocf_av_nb_mpi)
	    deallocate(Gcf_nb_mpi)
        deallocate(GcKf_nb_mpi)
        deallocate(del_tpmb_nb)
        deallocate(V_tpmb_nb)
        deallocate(G_tpmb_nb)
        deallocate(G_inmb_nb)
        deallocate(w_tpmb_nb)
        deallocate(logrhotot_meas)
        deallocate(logrhotot_run)
        deallocate(logrhotot_f)
        deallocate(logrhotot_f_mpi)
!for LIZembedding
        deallocate(lizmap)
        deallocate(Gc_liz)
        deallocate(Gc_liz_meas)
        deallocate(Gc_liz_run)
        deallocate(Gc_liz_f)
        deallocate(Gc_liz_mpi)
        deallocate(Gcbar_liz)
        deallocate(Gsc_liz)
        deallocate(ipvtliz)
        deallocate(R_liz)
        deallocate(sigma_l)
        return
        end subroutine deallocate_arrays
c*********************************************************
	subroutine allocate_PDOS

c*********************************************************
	use Global
	implicit none
c*********************************************************
	integer info,infot

	infot=0

	allocate(PDOS(0:e_div,Nc),stat=info)
	infot = infot + info

	if(infot.ne.0.and.myrank.eq.0) then
	     write(lud,*) 'allocate_PDOS',infot,' failures'
	     flush(lud)
          stop
        end if
	      
	return
	end subroutine allocate_PDOS
	
c-------------------------------------------------------------------------

	subroutine deallocate_PDOS
c*************************************************************************
	use Global
	implicit none
c*************************************************************************
  
	deallocate(PDOS)
      
	return
	end subroutine deallocate_PDOS
c*****************************************************************************

c*************************************************************************
c	This subroutine is provided for the openmp calls and functions
c*************************************************************************
!  OpenMP runtime library to be used in conjunction with Open64 Compiler Suites.
!
!  Copyright (C) 2003 - 2009 Tsinghua University.
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
!  
!  Contact information: HPC Institute, Department of Computer Science and Technology,
!  Tsinghua University, Beijing 100084, CHINA, or:
!
!  http://hpc.cs.tsinghua.edu.cn
        module omp_lib_kinds
        integer, parameter :: omp_integer_kind   = 4
        integer, parameter :: omp_logical_kind   = 4
        integer, parameter :: omp_lock_kind      = 8
        integer, parameter :: omp_nest_lock_kind = 8
        end module omp_lib_kinds


        module omp_lib
          use omp_lib_kinds
          integer, parameter :: openmp_version = 199910
              integer (kind=omp_integer_kind) :: ID
              integer (kind=omp_lock_kind) :: lck

          interface
            subroutine omp_destroy_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(inout) :: var
            end subroutine omp_destroy_lock
          end interface

          interface
            subroutine omp_destroy_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end subroutine omp_destroy_nest_lock
          end interface

          interface
            function omp_get_dynamic ()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_get_dynamic
            end function omp_get_dynamic
          end interface

          interface
            function omp_get_max_threads ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_max_threads
            end function omp_get_max_threads
          end interface

          interface
            function omp_get_nested ()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_get_nested
            end function omp_get_nested
          end interface

          interface
            function omp_get_num_procs ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_num_procs
            end function omp_get_num_procs
          end interface

          interface
            function omp_get_num_threads ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_num_threads
            end function omp_get_num_threads
          end interface

          interface
            function omp_get_thread_num ()
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_get_thread_num
            end function omp_get_thread_num
          end interface

          interface
            function omp_get_wtick ()
              use omp_lib_kinds
              double precision :: omp_get_wtick
            end function omp_get_wtick
          end interface

          interface
            function omp_get_wtime ()
              use omp_lib_kinds
              double precision :: omp_get_wtime
            end function omp_get_wtime
          end interface

          interface
            subroutine omp_init_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(out) :: var
            end subroutine omp_init_lock
          end interface

          interface
            subroutine omp_init_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(out) :: var
            end subroutine omp_init_nest_lock
          end interface

          interface
            function omp_in_parallel ()
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_in_parallel
            end function omp_in_parallel
          end interface

          interface
            subroutine omp_set_dynamic (enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind), intent(in) :: enable
            end subroutine omp_set_dynamic
          end interface

          interface
            subroutine omp_set_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(inout) :: var
            end subroutine omp_set_lock
          end interface

          interface
            subroutine omp_set_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end subroutine omp_set_nest_lock
          end interface

          interface
            subroutine omp_set_nested (enable)
              use omp_lib_kinds
              logical (kind=omp_logical_kind), intent(in) :: enable
            end subroutine omp_set_nested
          end interface

          interface
            subroutine omp_set_num_threads (nthreads)
              use omp_lib_kinds
              integer (kind=omp_integer_kind), intent(in) :: nthreads
            end subroutine omp_set_num_threads
          end interface

          interface
            function omp_test_lock (var)
              use omp_lib_kinds
              logical (kind=omp_logical_kind) :: omp_test_lock
              integer (kind=omp_lock_kind), intent(inout) :: var
            end function omp_test_lock
          end interface

          interface
            function omp_test_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_integer_kind) :: omp_test_nest_lock
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end function omp_test_nest_lock
          end interface

          interface
            subroutine omp_unset_lock (var)
              use omp_lib_kinds
              integer (kind=omp_lock_kind), intent(inout) :: var
            end subroutine omp_unset_lock
          end interface

          interface
            subroutine omp_unset_nest_lock (var)
              use omp_lib_kinds
              integer (kind=omp_nest_lock_kind), intent(inout) :: var
            end subroutine omp_unset_nest_lock
          end interface
        end module omp_lib


