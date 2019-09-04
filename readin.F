c**************************************************************************** 
        subroutine readin
c**************************************************************************** 
	use Global
	implicit none
c**************************************************************************** 
 
        character*8 waslin
        integer ios      
#define SIMPLE_SPRNG
#define FLOAT_GEN
#define USE_MPI
#include "sprng_f.h"
!for sprng random number generator
        SPRNG_POINTER :: junkPtr
    
!****************************************************************************
!	Open the data to be outputed at the end of the simulation
!****************************************************************************
       
        if(myrank==0) call profstart('readin')
        if(myrank==0)then 

c       open the output file
cc	wtime= OMP_get_wtime()
c       open the output file
                
        open(unit=42,file='run.dat',status='unknown')
        open(unit=132,file='mb_dos.dat',status='unknown')
        open(unit=133,file='mb_dos_av.dat',status='unknown') 
        open(unit=124,file='dos_vs_iter.dat',status='unknown')
        open(unit=126,file='sigma_vs_iter.dat',status='unknown')
        open(unit=188,file='rhol_vs_iter.dat',status='unknown')
c        open(unit=101,file='gamma_k_w_aa.dat',status='unknown')
c        open(unit=102,file='gamma_k_w_bb.dat',status='unknown')
c        open(unit=103,file='gamma_k_w_ab.dat',status='unknown')
        open(unit=99,file='causality.dat',status='unknown')
        open(unit=422,file='dos.dat',status='unknown')
        open(unit=825,file='gamma_check.dat',status='unknown')
        open(unit=823,file='GckKf_rhof.dat',status='unknown')
        open(unit=821,file='sigma.dat',status='unknown')
        open(unit=822,file='gamma.dat',status='unknown')
        open(unit=824,file='sigma_K_re.dat',status='unknown')
        open(unit=826,file='sigma_K_im.dat',status='unknown')
        open(unit=827,file='gamma_K_re.dat',status='unknown')
        open(unit=828,file='gamma_K_im.dat',status='unknown')
        open(unit=321,file='eigen_cas.dat',status='unknown')
        open(unit=830,file='GIPR2.dat',status='unknown')
c	open(unit=201,file='gamma-vs-iter.dat',access='append')

c       Readin the parameters
        open(unit=5,file='in1.data',status='old')
        read(5,"(1x,a8)") waslin
        read(5,*) cluster
	    write(6,*) 'PROGRAM CTMT FOR CLUSTER ', cluster
!       Extract nc from the string cluster
        read(cluster,"(i5)",iostat=ios) Nc     ! should barf with characters
!       ibetts=0
        if (ios.gt.0) then
!         ibetts=1
          read(cluster,"(i4)",iostat=ios) Nc   ! ignore trailing charac.
          if (ios>0) then
            read(cluster,"(i3)",iostat=ios) Nc ! ignore trailing charac.
            if (ios>0) then
              read(cluster,"(i2)",iostat=ios) Nc ! ignore trailing charac.
              if (ios>0) then
        	read(cluster,"(i1)",iostat=ios) Nc ! ignore trailing charac.
        	if (ios>0) then
        	  write(lud,*) "Cannot read Nc from string ",cluster
        	endif
              endif
            endif
          endif
        endif
	write(6,*) "Nc=", Nc

        write(6,*) 'enter meas,run,niter,nover'
        read(5,13) waslin
        read(5,*) meas,run,niter,nover
        write(6,100) meas,run,niter,nover
 100    format(' meas=',i3, ' run=',i3,' niter= ',i3,' nover= ',i3)
        if(meas.lt.1000) write(6,*) 'increase meas'
        if(run.lt.10) write(6,*) 'increase run'
        write(6,*) 'enter rand #,iso,dfac,ifac'
        read(5,13) waslin
        read(5,*) iseed,iso,dfac,ifac
        write(6,101) iseed,iso,dfac,ifac
 101    format('iseed= ',i10,' iso=',i3,' dfac=',f7.4,' ifac= ',f7.4)

        isi=0
	if(ifac.gt.1.00001) isi=1
        isd=0
	if(abs(dfac-1.0).gt.0.01) isd=1

        write(6,*) 'enter ed,tprime,V'
        read(5,13) waslin
 13     format(1x,a8)
        read(5,*) ed,tprime,V0
        write(6,102) ed,tprime,V0

 102	format('ed = ',f9.6,' tprime=',f9.6,' V=',f9.6)
 
 	if(V0.ge.0) then
	  isv=1
	else
	  isv=-1
	  V0=-V0
	end if

	write(6,*) 'enter delta, nwn'
	read(5,13) waslin
	read(5,*) delta,nwn
	write(6,103) delta,nwn
 103	format('delta= ',f9.6, ' nwn=',i6)
 
 
         write(6,*) 'enter ca, tables_type'
        read(5,*) waslin
        read(5,*) ca,tables_type
        write(6,106) ca,tables_type
 106	format('ca = ',f9.6,2x,'tables_type=',i2)
        if(tables_type.lt.4) then
            write(6,*) 'enter taa, tab, Vaa, Vab'
            read(5,*) waslin
            read(5,*) aa(1),aa(2),ab,tz,V1a,V12a
            write(6,"('taa1 = ',f9.6,2x,'taa2 = ',f9.6,2x,'tab = ',f9.6,2x,'tz = ',f9.6)")aa(1),aa(2),ab,tz
            write(6,"('Vaa = ',f9.6,2x,'Vab = ',f9.6)")V1a,V12a
            write(6,*) 'enter liz linear size'
            read(5,*) waslin
            read(5,*) liz
            lizsize=liz**ndim
            write(6,"('liz linear size = ',i3,2x,'liz size = ',i3)"),liz,lizsize   
        end if
 
        close(5)
        end if !if myrank==0
        
! First, pass the variables to all processors
!       INTEGERS
        call MPI_BCAST(cluster,5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) ! cluster
        call MPI_BCAST(meas,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)      ! meas
        call MPI_BCAST(run,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       ! run
        call MPI_BCAST(niter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)     ! niter
        call MPI_BCAST(nover,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)      ! ntot
        call MPI_BCAST(iseed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)     ! iseed
        call MPI_BCAST(iso,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       ! iso
        call MPI_BCAST(isd,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       ! isd
        call MPI_BCAST(Nc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)        ! Nc
        call MPI_BCAST(nwn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       ! nwn
        call MPI_BCAST(isi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       ! isi
        call MPI_BCAST(isv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)       ! isv
        call MPI_BCAST(tables_type,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)     ! table type
        call MPI_BCAST(dfac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)        ! dfac
        call MPI_BCAST(ifac,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)        ! ifac
        call MPI_BCAST(ed,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)          ! ed
        call MPI_BCAST(delta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)       ! delta
        call MPI_BCAST(V0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)          ! V0
        call MPI_BCAST(tprime,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)      ! tprime
        call MPI_BCAST(ca,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)          ! Concentration
        
        if(tables_type.lt.4) then
        call MPI_BCAST(aa,2,MPI_REAL8,0,MPI_COMM_WORLD,ierr)          !taa
        call MPI_BCAST(ab,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)          !tab
        call MPI_BCAST(tz,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)          !tz
        call MPI_BCAST(V1a,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)        !Vaa
        call MPI_BCAST(V12a,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr)        !Vab 
        call MPI_BCAST(liz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   !liz
        call MPI_BCAST(lizsize,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)   !lizsize
        end if
!set up the random number library -- which gives an
!independent random number sequence on each processor 
        !junkPtr = init_sprng(SPRNG_CMRG,iseed,SPRNG_DEFAULT) !sprng2.0
        junkPtr = init_sprng(iseed,SPRNG_DEFAULT,3)    !sprng5.0
        if(myrank==0) call profend('readin')
        return
        end

        subroutine init_parameters
c************************************************************
c       Initialize IPT in preparation for the first iteration
c************************************************************
        use Global
        implicit none
c************************************************************
        integer n,m,i,j,ic,jc,info
        real(kind) ::  time1,time2,r1,r2,r3,e
	complex*16 :: z1
        character*72 linemc
cH:
        complex*16 :: A1,B1,C1,D1
c************************************************************
        
c	Allocate arrays
		
!	  if(iso.ge.1.) then ! readin sigma from a file	
  
      if(iso.eq.0.) then ! set sigma to a constant
	call init_coarsegraining   

	endif
	return
	end subroutine init_parameters

        subroutine init_coarsegraining
        use Global
        implicit none
c************************************************************
        integer n,m,i,j,ic,jc,info,guess
        real(kind) ::  time1,time2,r1,r2,r3,e
        complex*16 :: z1,chir(Nc)
        character*72 linemc
cH:
        complex*16 :: A1,B1,C1,D1
c************************************************************
        if(myrank==0) call profstart('init_coarsegrain')
        if(myrank==0)then
        guess=1   
	   
c ----  MB


      if(MB1) then  !inital Gamma and Fourier transofrm the bare dispersion in real space
        Gamma_new_nb(:,:,:,:)=zeroc
        Gamma_old_nb(:,:,:,:)=zeroc  
      
        if(tables_type.lt.4)then  
        do ic=1,Nc  
            chir(ic)=zeroc 
              do igroup=1,ngroup                                   
               chir(ic)=chir(ic) + 
     &		(1.d0/real(ngroup))*Epsbar(ickequ(ic,igroup))
              end do 
            !Epsbar(ic)=chir(ic)	
        end do
        do i=1,nb
           do j=1,nb
              if(i==j)then
                 do ic=1,Nc
                    Epsb_nb(i,i,ic)=aa(i)*Epsbar(ic)
                 end do
              else
                 do ic=1,Nc
                    Epsb_nb(i,j,ic)=ab*Epsbar(ic)
                 end do
              end if
           end do
        end do      
        end if!if tables_type<4
        
        if(init_sigma)then!now only works for tables_type<4
           do i=1,nb
               do j=1,nb
                   if(i.eq.j) then
                       sigma_nb(i,i,:,:)=cmplx(0.d0,-0.01d0)
                   else
                       sigma_nb(i,j,:,:)=cmplx(0.d0,-0.0d0)
                   end if
               end do
           end do
         if(liz_embedding)then
             sigma_l(:)=sigma_nb(1,1,1,:)
             call readlizlist
             call set_liz_R
         else!not liz_embedding  
!do initial coarse graining to get GcKf(K,w)
         if(tables_type<4)then
          GcKf_nb = zeroc
          do ic=1,Nc
            do n=-nwn,nwn
            do i=0,e_div
             e=estart+dfloat(i)*estep
             do jc=1,nb
                 do j=1,nb
                    if(jc.eq.j) Gc_temp(jc,jc,ic,n)=wn(n)-sigma_nb(jc,jc,ic,n)-aa(jc)*e
                    if(jc.ne.j) Gc_temp(jc,j,ic,n)=-sigma_nb(jc,j,ic,n)-ab*e
                 end do
             end do
             call zgetrf(nb,nb,Gc_temp(:,:,ic,n),nb,ipvtnb,info)
                 if (info.eq.0)then
                 call zgetri(nb,Gc_temp(:,:,ic,n),nb,ipvtnb,work,lwork,info)
             else
               print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
                 stop
             end if
             GcKf_nb(:,:,ic,n)=GcKf_nb(:,:,ic,n)+PDOS(i,ic)*estep*Gc_temp(:,:,ic,n)
            end do !eps-loop
            end do !wn-loop
          end do !ic-loop
c normalize
          GcKf_nb = GcKf_nb*real(Nc) 
!invert GcKf(K,w) to calculate Gcsript
          Gc_temp=GcKf_nb
          do ic=1,Nc
          do n=-nwn,nwn
             call zgetrf(nb,nb,Gc_temp(:,:,ic,n),nb,ipvtnb,info)
                 if (info.eq.0)then
                 call zgetri(nb,Gc_temp(:,:,ic,n),nb,ipvtnb,work,lwork,info)
             else
               print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
                 stop
             end if
          end do
          end do
          Gloc_inv_nb=Gc_temp
!calculate GcKfsc=1/(GcKf^-1+sigma_nb)
          Gc_temp=Gloc_inv_nb+sigma_nb
          do ic=1,Nc
          do n=-nwn,nwn
             call zgetrf(nb,nb,Gc_temp(:,:,ic,n),nb,ipvtnb,info)
                 if (info.eq.0)then
                 call zgetri(nb,Gc_temp(:,:,ic,n),nb,ipvtnb,work,lwork,info)
             else
               print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
                 stop
             end if
          end do
          end do
          GcKfsc_nb=Gc_temp
!invert GcKfsc_nb to calculate Gamma,temporarily restored in Gloc_inv_nb
          Gc_temp=GcKfsc_nb
          do ic=1,Nc
          do n=-nwn,nwn
             call zgetrf(nb,nb,Gc_temp(:,:,ic,n),nb,ipvtnb,info)
                 if (info.eq.0)then
                 call zgetri(nb,Gc_temp(:,:,ic,n),nb,ipvtnb,work,lwork,info)
             else
               print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
                 stop
             end if
          end do
          end do
          Gloc_inv_nb=Gc_temp
!calculate Gamma=w-Epsk-GcKfsc_nb^-1
         do i=1,nb
            do j=1,nb
            do n=-nwn,nwn
            do ic=1,Nc   
               if(i==j)Gamma_old_nb(i,i,ic,n)=wn(n)-Gloc_inv_nb(i,i,ic,n)-Epsb_nb(i,i,ic)
               if(i/=j)Gamma_old_nb(i,j,ic,n)=-Gloc_inv_nb(i,j,ic,n)-Epsb_nb(i,j,ic)
            end do
            end do
            end do
         end do
         end if! if tables_type<4
         end if! if liz_embedding
        else!if not init_sigma   
           do i=1,nb
               do j=1,nb
                   if(i.eq.j) then
                       Gamma_old_nb(i,i,:,:)=cmplx(0.d0,-0.1d0)
                   else
                       Gamma_old_nb(i,j,:,:)=cmplx(0.d0,-0.0d0)
                   end if
               end do
           end do
        end if!init_sigma
      end if!if MB1

      call set_epsb_rs_nb
!------ end MB1       
        
 	      write(6,*) '***************************************************************'					      
	      write(6,*) 'code initialized with a guess of the Hybridization Rate'
 	      write(6,*) '***************************************************************'
	      write(42,*) 'code initialized with a guess of the Hybridization Rate'	 
      end if!if myrank==0      

!	Broadcast to all the processors
        
      if(liz_embedding)then
        call MPI_BCAST(sigma_l,2*nwn+1,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(R_liz,ndim*lizsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(lizmap,Nc*lizsize,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      else
        call MPI_BCAST(Gamma_old_nb,nb*nb*Nc*(2*nwn+1),MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(Epsb_rs_nb,nb*nb*Nc*Nc,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(GcKfsc_nb,nb*nb*Nc*(2*nwn+1),MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
      end if
      if(myrank==0) call profend('init_coarsegrain')
      return
	end subroutine init_coarsegraining

      subroutine set_epsb_rs_nb
       use global
       implicit none
       complex(kind), allocatable :: Epsb_rs_nb_temp(:,:,:)
       integer :: i,j,ic,jc,info
       complex(kind) :: csum

       allocate(Epsb_rs_nb_temp(nb,nb,Nc),stat=info)
       if(info.ne.0)then
           write(*,*)'allocate failure: Epsb_rs_nb_temp'
           stop
       end if
!*************************************
!     Fourier transform to real space 
!*************************************       
       Epsb_rs_nb_temp=zeroc
       do i=1,nb
       do j=1,nb
       do ic=1,Nc
          csum=zeroc
          do jc=1,Nc
              csum=csum+Epsb_nb(i,j,jc)*FTCoefs_K_to_R(jc,ic)
          end do !jc
          Epsb_rs_nb_temp(i,j,ic)=csum
       end do !ic
       end do !j
       end do !i
!*************************************
!     Form Epsb_rs_nb(i,j,ic,jc) 
!*************************************
        do ic=1,Nc
        do jc=1,Nc
        do i=1,nb
        do j=1,nb    
          Epsb_rs_nb(i,j,ic,jc)=Epsb_rs_nb_temp(i,j,icrdiff(ic,jc))
        end do
        end do
        end do
        end do

        deallocate(Epsb_rs_nb_temp)
        end subroutine set_epsb_rs_nb
 
        subroutine readlizlist
            use global
            implicit none
            integer :: i,j,k,c,info
            integer :: d, L, lizr
            character(len=400) :: junk1, junk2
            open(unit=80,file='lizList',status='old')
            read(80,*)  junk1, d, junk1, L, junk1, lizr
            if(lizsize /= lizr**d)then
                write(*,*)'lizsize error'
                stop
            end if
            do i=1,Nc
               read(80,*) c,  junk1
               do j=1,lizsize
                  read(80,*) k , junk1
                  lizmap(c,j) = int(k)
               end do
            end do
        end subroutine        
      
      subroutine set_liz_R
         use Global
         implicit none
         integer :: i,j,k,count,n
         count=1
         n=liz/2
         R_liz(1,:)=(/0,0,0/)
         do i=-n,n
            do j=-n,n
               do k=-n,n
                  if(abs(i)+abs(j)+abs(k) /= 0)then
                      count=count+1
                      R_liz(count,:)=(/i,j,k/)
                  end if
               end do
            end do
         end do
      end subroutine set_liz_R
