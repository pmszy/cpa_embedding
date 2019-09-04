c***************************************************************************************
        subroutine geodgen				
c       This block makes the initial greens function.  This
c	function corresponds to the field configuration where all 
c 	s(i)=0.0.  It does this by summing the Fourier transform. 
c	The coarse-graining is done in parallel.  
 
c***************************************************************************************
	use Global
c	use omp_lib
	implicit none
c***************************************************************************************
c	integer,parameter:: nt_broyden = (2*5+1)*20
        integer i,j,n,ic,jc,info
        real(kind) :: time1,time2
        complex(kind) :: csum
         if(myrank==0) call profstart('geod')
c	wtime= OMP_get_wtime()
c	time1= secnds(0.0)
c***************************************************************************************
c	First, we need a map from the process rank to the cluster site
c	indices.  Each process gets its K-value, where its is determined by
c	calculate the initial reciprocal-space green function Gckf and GcKfsc.
c       Symmetrize the GcKfsc and  x <---> +/-x       
c       hence the self energy by      \ /	      
c       imposing the point group       X			    
c       symmetries.		      / \			    
c       			   y <---> +/-y 
c***************************************************************************************
c: PUT IF SWITCH for a regular CTMT

   
       
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                  Multiband 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        
      if (MB1) then
      if(liz_embedding)then
       if(free_standing)then
        call set_Gsc_liz_free
       else   
        if(myrank==0) call profstart('set_Gc_liz')
        !call set_Gc_liz
        call set_Gc_liz_cg
        if(myrank==0) call profend('set_Gc_liz')
        if(myrank==0) call profstart('Gsc')
        do n=-nwn,nwn
           call zgetrf(lizsize,lizsize,Gc_liz(:,:,n),lizsize,ipvtliz,info)
           if (info.eq.0)then
              call zgetri(lizsize,Gc_liz(:,:,n),lizsize,ipvtliz,work,lwork,info)
           else
              print*, 'info=/0 Matrix Inversion Failed in geod: Info=',info
              stop
           end if
           Gsc_liz(:,:,n)=Gc_liz(:,:,n)
           do i=1,lizsize
              !Gsc_liz store the inverse of Gscipt_liz
              Gsc_liz(i,i,n)=Gsc_liz(i,i,n)+sigma_l(n)
           end do!i   
        end do!n
        if(myrank==0) call profend('Gsc')
       end if 
      else !if not liz_embedding    
C++++++++++++++++++++++  1. FT to RS Delta  +++++++++++++++++++++++ 
        Gamma_nb(:,:,:,:)=zeroc
        Gamma_nb = Gamma_old_nb
        
        do i=1,nb
            do j=1,nb
                do ic=1,Nc  
                    chi(ic,:)=zeroc   
                    do igroup=1,ngroup                                   
                        chi(ic,:)=chi(ic,:) + 
     &		(1.d0/real(ngroup))*Gamma_nb(i,j,ickequ(ic,igroup),:)
                    end do 
                    !Gamma_nb(i,j,ic,:)=chi(ic,:)			                        
                end do
            end do
        end do

        



c 2. Fourier Transform it to cluster space
      do i=1,nb
          do j=1,nb
              do ic=1,Nc
                  do n=-nwn,nwn
                      csum=zeroc
                      do jc=1,Nc
                          csum=csum+Gamma_nb(i,j,jc,n)*FTCoefs_K_to_R(jc,ic)
                      end do
                      Gamma1_nb(i,j,ic,n)=csum
                  end do
              end do
          end do
      end do


c***************************************************************************************
c	form Gamma_ij
c***************************************************************************************
        Gamma_rs_nb=zeroc
        
        do n=-nwn,nwn
        do ic=1,Nc
        do jc=1,Nc
        do i=1,nb
        do j=1,nb    
          Gamma_rs_nb(i,j,ic,jc,n)=Gamma1_nb(i,j,icrdiff(ic,jc),n)
        end do
        end do
        end do
        end do
        end do
        
        if(Gscript.and.init_sigma)then
            GcKfsc_nb1=GcKfsc_nb                      
        do i=1,nb
            do j=1,nb
                do ic=1,Nc  
                    chi(ic,:)=zeroc   
                    do igroup=1,ngroup                                   
                        chi(ic,:)=chi(ic,:) + 
     &        (1.d0/real(ngroup))*GcKfsc_nb(i,j,ickequ(ic,igroup),:)
                    end do 
                    !GcKfsc_nb1(i,j,ic,:)=chi(ic,:)                       
                end do
            end do
        end do
        do i=1,nb
          do j=1,nb
              do ic=1,Nc
                  do n=-nwn,nwn
                      csum=zeroc
                      do jc=1,Nc
                          csum=csum+GcKfsc_nb1(i,j,jc,n)*FTCoefs_K_to_R(jc,ic)
                      end do
                      GcKfsc_nb2(i,j,ic,n)=csum
                  end do
              end do
          end do
        end do
        GcRfsc_nb=zeroc
        
        do n=-nwn,nwn
        do ic=1,Nc
        do jc=1,Nc
        do i=1,nb
        do j=1,nb    
          GcRfsc_nb(i,j,ic,jc,n)=GcKfsc_nb2(i,j,icrdiff(ic,jc),n)
        end do
        end do
        end do
        end do
        end do 
        end if
        !if(myrank==0) then
        !  write(*,*)GcKfsc_nb(:,:,:,0)
        !  stop
        !end if
        end if !if liz_embedding
        END IF !if mb1
       

c 4. Go To dategrup.f to first construct Gsc_ij(v) which needs to be inverted and construct G(v)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc        
c                  end of Multiband
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        



     
      

c***************************************************************************************
c	Initialize some variables and accumulators
c***************************************************************************************

       

!     MB1
      Gdff_nb=zeroc
      rhorun_nb=zeror
      logrholoc_run_nb=zeror
      rholoc_av_run_nb=zeror
      logrhotot_run=zeror
      rhof_nb=zeror
      logrholoc_f_nb=zeror
      rholocf_av_nb=zeror
      rho_meas_nb=zeror
      logrholoc_meas_nb=zeror
      rholoc_av_meas_nb=zeror
      logrhotot_meas=zeror

      Gc_meas_nb=zeroc
      Gcrun_nb=zeroc
      Gcf_nb=zeroc
       
      Gc_liz_meas=zeroc
      Gc_liz_run=zeroc
      Gc_liz_f=zeroc 

      GIPR2_meas=zeror
      GIPR2_run=zeror
      GIPR2_f=zeror
      logGIPR2_meas=zeror
      logGIPR2_run=zeror
      logGIPR2_f=zeror
 







c	time2= secnds(time1) 
c	wtime= OMP_get_wtime()-wtime
c	write(6,*) 'time in geod= ',time2
c	write(6,*) 'time in geod= ',wtime

         if(myrank==0) call profend('geod')
        return
        end

	subroutine HilbertTwithlogcorr(Gt,rlsigma,rhosigma)
	use Global
	implicit none
	integer :: i,j,jc,ic,n
	integer :: info
	real*8 :: rhosigma(Nc,-nwn:nwn),rlsigma(Nc,-nwn:nwn)
	real*8 :: woj,dspj,resgj,spi,dspi,r1,dwi,woji1,woji2
	real*8, allocatable  :: rlsg(:,:)
	complex(kind) :: Gt(Nc,-nwn:nwn)
	!	************************************************************

	allocate(rlsg(Nc,-nwn-1:nwn+1),stat=info)
		        
	rlsigma=zeror

	!		      oo
	!		       /
	!	rlsigma(w)= -P | dw' rhosigma(w')
	!		       /     -----------
	!	             -oo       w' - w


	!		KRAMERS-KRONIG


	    do ic=1,Nc
		do j=-nwn+1,nwn
		!	interlacing grid

		    woj=.5D0*(wn(j-1)+wn(j))
		    dspj=rhosigma(ic,j)-rhosigma(ic,j-1)
		    resgj=0.0D0

		    		 	
		    do i=-nwn,j-2

		        spi=rhosigma(ic,i)
		        dspi=rhosigma(ic,i+1)-rhosigma(ic,i)
		        dwi=wn(i+1)-wn(i)
		        woji1=wn(i)-woj
		        woji2=wn(i+1)-woj
		        r1=dlog(woji2/woji1)


		        resgj=resgj-(spi*r1 + dspi )
		        resgj=resgj-(dspi/dwi)*(woj -wn(i))*r1

		    end do

		!	 skip the interval (j-1) to j
		    		 	
		    do i=j,nwn-1

		        spi=rhosigma(ic,i)
		        dspi=rhosigma(ic,i+1)-rhosigma(ic,i)
		        dwi=wn(i+1)-wn(i)
		        woji1=wn(i)-woj
		        woji2=wn(i+1)-woj
		        r1=dlog(woji2/woji1)

		        resgj=resgj-(spi*r1 + dspi )
		        resgj=resgj-(dspi/dwi)*(woj -wn(i))*r1

		    end do


		    resgj=resgj - dspj
		    rlsg(ic,j)=resgj
		end do
	    enddo

	    rlsg(:,nwn+1)=rlsg(:,nwn)
	    rlsg(:,-nwn)=rlsg(:,-nwn+1)
	    	    	
	    do ic=1,Nc
		do i=-nwn,nwn
		    rlsigma(ic,i)=0.5d0*(rlsg(ic,i)+rlsg(ic,i+1))
		end do
	    enddo


	!c	Now load the real and imaginary part into one
	    Gt = zeroc
	    do ic=1,Nc
		do n=-nwn,nwn
		    Gt(ic,n)= rlsigma(ic,n) - ii*pi*(rhosigma(ic,n))
		end do
	    enddo

	    deallocate(rlsg)

	return
	end subroutine HilbertTwithlogcorr


!**************************************************************
	subroutine HilbertT(Gt,Gtr,Rhot)
!**********************************************************************
!c	Do Hilbert transformation to calculated the Typical Green function (TGF).
!       First the real part.  Note, this is a principal value integral
!       call profstart('HilbertT_Sumup') ! To calculate the time spent in HT.

!		      oo
!		       /
!	G_typ(k,w)= P | dw' rho_typ(k,w') -i*pi*rho_typ(k,w)
!		       /     -----------
!	             -oo       w - w'
!**********************************************************************
	use Global
	implicit none
!**************************************************************
	integer :: n,m,ic
	real(kind) :: Gtr(Nc,-nwn:nwn),Rhot(Nc,-nwn:nwn)
	complex(kind) :: Gt(Nc,-nwn:nwn),sigmacg(-nwn:nwn,Nc)
!**************************************************************

	    Gtr = zeror
	    do ic=1,Nc
		do n=-nwn,nwn
		    do m=-nwn,n-1
		        Gtr(ic,n) = Gtr(ic,n) + 
     &		        Rhot(ic,m)*dwn(m)/(wn(n)-wn(m))
		    enddo
		    do m=n+1,nwn
		        Gtr(ic,n) = Gtr(ic,n) + 
     &		        Rhot(ic,m)*dwn(m)/(wn(n)-wn(m))
		    enddo
		enddo
	    enddo

	!c	Now load the real and imaginary part into G_typ
	    Gt = zeroc
	    do ic=1,Nc
		do n=-nwn,nwn
		    Gt(ic,n)= Gtr(ic,n) - ii*pi*(Rhot(ic,n))
		end do
	    enddo
	end subroutine HilbertT

      
      subroutine set_Gc_liz
       use global
       implicit none
       integer :: i,j,idim,ic,n,info
       complex(kind) :: csum
       real(kind) :: phase

!*************************************
!     Fourier transform to real space 
!*************************************       
       Gc_liz=zeroc
       do n=-nwn,nwn
          do i=1,lizsize
          do j=1,lizsize
             csum=zeroc
             do ic=1,Nc
                phase=zeror
                do idim=1,ndim
                   phase=phase+(R_liz(i,idim)-R_liz(j,idim))*Kc(idim,ic)
                end do
                csum=csum+oner/(wn(n)-Epsbar(ic)-sigma_l(n)+eta)*exp(-ii*phase)/Nc
             end do !ic
             Gc_liz(i,j,n)=csum
          end do !i
          end do !j
       end do!n
      end subroutine set_Gc_liz
      
      subroutine set_Gc_liz_cg
       use global
       implicit none
       integer :: i,j,idim,ic,n,info
       complex(kind) :: csum
       real(kind) :: phase,e
!do initial coarse graining to get GcKf(K,w)
          GcKf_nb = zeroc
          do ic=1,Nc
            do n=-nwn,nwn
            do i=0,e_div
             e=estart+dfloat(i)*estep
             GcKf_nb(1,1,ic,n)=GcKf_nb(1,1,ic,n)+PDOS(i,ic)*estep*
     &       (oner/(wn(n)-sigma_l(n)-e+eta))
            end do !eps-loop
            end do !wn-loop
          end do !ic-loop
c normalize
          GcKf_nb = GcKf_nb*real(Nc)
       
!*************************************
!     Fourier transform to real space 
!*************************************       
       Gc_liz=zeroc
       do n=-nwn,nwn
          do i=1,lizsize
          do j=1,lizsize
             csum=zeroc
             do ic=1,Nc
                phase=zeror
                do idim=1,ndim
                   phase=phase+(R_liz(i,idim)-R_liz(j,idim))*Kc(idim,ic)
                end do
                csum=csum+GcKf_nb(1,1,ic,n)*exp(-ii*phase)/real(Nc)
             end do !ic
             Gc_liz(i,j,n)=csum
          end do !i
          end do !j
       end do!n
      end subroutine set_Gc_liz_cg
      
      subroutine set_Gsc_liz_free
       use global
       implicit none
       integer :: i,j,idim,ic,n,info
       complex(kind) :: csum,tp_liz(lizsize,lizsize)
       real(kind) :: phase

!*************************************
!     Fourier transform to real space 
!*************************************       
       tp_liz=zeroc
       do i=1,lizsize
       do j=1,lizsize
          csum=zeroc
          do ic=1,Nc
             phase=zeror
             do idim=1,ndim
                phase=phase+(R_liz(i,idim)-R_liz(j,idim))*Kc(idim,ic)
             end do
             csum=csum+Epsbar(ic)*exp(-ii*phase)/Nc
          end do !ic
          tp_liz(i,j)=csum
       end do !i
       end do !j
       Gsc_liz=zeroc
       do n=-nwn,nwn
          do i=1,lizsize
             do j=1,lizsize
                if(i.eq.j)then
                   Gsc_liz(i,j,n)=wn(n)+eta-tp_liz(i,j)
                else
                   Gsc_liz(i,j,n)=-tp_liz(i,j)
                end if 
             end do
          end do
       end do!n
      end subroutine set_Gsc_liz_free
      
