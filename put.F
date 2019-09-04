c**************************************************************************** 
        subroutine output(istate)
c**************************************************************************** 

c	This block is designed to process the binned data
c	contained in gf, etc. (where the binning was
c	done sequentially).
 
cc***************************************************************************************
	use Global
c	use omp_lib
	implicit none
c***************************************************************************************
      
       integer i,jc,istate,ic,n,rmeas,m,istart,iend
       real*8 :: r0,r1,r2,r3,r4,r5,r6,r7,r8,min_w,
     &  time1,time2,sum1(Nc),sum2(Nc),rholoc_typ_nb(nb,nb,-nwn:nwn),
     &  r(8),rhotot_typ(-nwn:nwn)  
       integer info,k,l,l1,j,check
       real*8 ::kx,ky,kz,data_nb(nb,nb,Nc,-nwn:nwn),data_rho_nb(nb,nb,-nwn:nwn),  
     &     Gc_nb_real(nb,nb,Nc,-nwn:nwn)
       real*8 :: e,rv1(3),rho_av(-nwn:nwn),rhof_av(-nwn:nwn),rhol_nb(nb,nb,-nwn:nwn)
       complex*16 :: Eps_temp(nb,nb),m_r1(Nc,-nwn:nwn),sigma_loc_temp,sigma_nb_loc(nb,nb,-nwn:nwn),
     &               Gc_nb_loc(nb,nb,-nwn:nwn),Gsc_nb_loc(nb,nb,-nwn:nwn)
       integer lworkr
       real*8 :: rho_temp_nb(nb,nb),rho_temp_nb1(nb,nb),workr(100),Eval(nb),Evect(nb,nb) 
         real test2    
         real*8 :: Gamma_ave(-nwn:nwn,3)
        lworkr=100
c***************************************************************************************
         if(myrank==0) write(*,*) 'test-put0',iter
         if(myrank==0) call profstart('output')

         if(istate.eq.0) then
          if(myrank==0) then
	      write(42,*) 'cluster type '
	      write(42,*) cluster
          if(.not. mb_tmt)then
             write(42,*) 'Results from CPA embedding for nwn,Nc='
          else
             write(42,*) 'Results from TMT embedding for nwn,Nc='
          end if   
          write(42,*) nwn,Nc
          write(42,*) 'parameters follow'
          write(42,*) 'meas,run,niter,dfac,ifac'
          write(42,2) meas,run,niter,dfac,ifac
          end if   
	     rmeas=meas
	     do n=1,niter-1
	      rmeas=int(ifac*rmeas)
 	     end do
 2        format(i5,2x,i5,2x,i5,1x,f5.2,2x,f5.2)
         if(myrank==0)then
	      write(42,*) 'random number, iso'
	      write(42,*) iseed,iso
 4        format(1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4)
          write(42,"('taa1 = ',f9.6,2x,'taa2 = ',f9.6,2x,
     &             'tab = ',f9.6,2x,'tz = ',f9.6)")aa(1),aa(2),ab,tz
          write(42,"('Vaa = ',f9.6,2x,'Vab = ',f9.6)")V1a,V12a
          write(42,"('Liz = ',i3,2x,'Lizsize = ',i4)")liz,lizsize
          write(42,*) ' '
         end if
	     return
         end if

cccccccccccccccccccccccccccccccccccccccccccc
c             MB begin
cccccccccccccccccccccccccccccccccccccccccccc
 
        IF  (MB1) THEN
         if(Liz_embedding)then
             Gc_liz_f=Gc_liz_f/real(run)
             if(mb_tmt)then
                logrholoc_f_nb=logrholoc_f_nb/real(run)
             end if
             Gc_liz_mpi=zeroc
             !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
             call MPI_REDUCE( Gc_liz_f, Gc_liz_mpi, lizsize*lizsize*(2*nwn+1),
     &            MPI_COMPLEX16,MPI_SUM, 0, MPI_COMM_WORLD, ierr)
             if(mb_tmt)then
             logrholoc_f_nb_mpi=zeror
             call MPI_REDUCE( logrholoc_f_nb, logrholoc_f_nb_mpi, nb*nb*(2*nwn+1), 
     &            MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
             end if 
             if(myrank == 0)then 
!	Normalize the runs over the no. of processors
                Gc_liz_f=Gc_liz_mpi/dfloat(nprocs)
                if(mb_tmt)then
                   logrholoc_f_nb=logrholoc_f_nb_mpi/dfloat(nprocs)
                   do n=-nwn,nwn
                      Gc_liz_f(:,:,n)=Gc_liz_f(:,:,n)*dexp(logrholoc_f_nb(1,1,n))
                      !Gc_liz_f(:,:,n)=Gc_liz_f(:,:,n)*(logrholoc_f_nb(1,1,n))
                   end do
                end if
                Gcbar_liz=Gc_liz_f 
                do n=-nwn,nwn
                call zgetrf(lizsize,lizsize,Gcbar_liz(:,:,n),lizsize,ipvtliz,info)
	            if (info.eq.0)then
	               call zgetri(lizsize,Gcbar_liz(:,:,n),lizsize,ipvtliz,work,lwork,info)
                else
                 print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
	               stop
                end if
                end do
                write(*,*)'rho(0)= ',-dimag(Gc_liz_f(1,1,0))/pi
                write(126,'(i4,20(1x,f10.4))')iter,dimag(sigma_l(0))
                write(124,'(i4,1x,f9.6,1x,f10.4)')iter,-dimag(Gc_liz_f(1,1,0))/pi
                !update sigma
                sigma_nb_loc(1,1,:)=(Gsc_liz(1,1,:)-Gcbar_liz(1,1,:))*(1.0-dfac)+dfac*sigma_l(:)
                test=zeror
                test2=zeror
                test=sum(dabs(dimag(sigma_nb_loc(1,1,:)))) 
                test2=sum(dabs(dimag(sigma_l(:)-sigma_nb_loc(1,1,:))))
                sigma_l(:)=sigma_nb_loc(1,1,:)
                if(iter==niter)then
                  do n=-nwn,nwn  
                     write(422,'(f9.4,1x,f9.6)')wn(n),-dimag(Gc_liz_f(1,1,n))/pi
                     write(821,'(f9.4,1x,f9.6)')wn(n),-dimag(sigma_l(n))
                  end do   
                end if  
                
                if(iter.ge.2) then					      
                  r1=test2/(test*dfloat(Nc))
                  write(6,*) ' test= ',r1
 	              write(6,*) '**********************************************************'		      
                  write(42,*) ' test= ',r1 
 	              write(42,*) '**********************************************************'		      
                end if
             end if !myrank==0  
         call MPI_BCAST(sigma_l,2*nwn+1,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
         if(iter==niter)then
             GIPR2_f=GIPR2_f/real(run)
             logGIPR2_f=logGIPR2_f/real(run)
             GIPR2_mpi=zeror
             !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
             call MPI_REDUCE( GIPR2_f, GIPR2_mpi, 1,
     &            MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
             logGIPR2_mpi=zeror
             call MPI_REDUCE( logGIPR2_f, logGIPR2_mpi, 1,
     &            MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
             if(myrank == 0)then 
!	Normalize the runs over the no. of processors
                GIPR2_f=GIPR2_mpi/dfloat(nprocs)
                logGIPR2_f=logGIPR2_mpi/dfloat(nprocs)
                GIPR2_typ=dexp(logGIPR2_f)
                write(830,'(f9.6,1x,f9.6)') GIPR2_f, GIPR2_typ
             end if
         end if
         else! not Liz_embedding

c***************************************************************************************
c  --- normazile DOSs
c  1.  --- run-normalize of rho_aa,rho_bb,rho_ab,rho_ba
c***************************************************************************************
       
        rhof_nb=rhof_nb/dfloat(run)    !! rhof_aa=<-1/pi(Gaa(K,w))>
        if(Gfav) Gcf_nb=Gcf_nb/dfloat(run)
         if (mb_tmt) then
         logrholoc_f_nb=logrholoc_f_nb/dfloat(run)  ! rho_loc=<Ln(-1/pi(Gaa(ii,w)))>
         
         if(iter.eq.niter) then
           rholocf_av_nb=rholocf_av_nb/dfloat(run)
           logrhotot_f=logrhotot_f/dfloat(run)
         endif
  
         end if
         

!***********************************************************************
!       DATA REDUCTION.

!       Move all data to proc 0.  These must be called from all processors
!***********************************************************************

        rhof_nb_mpi=zeror
        logrholoc_f_nb_mpi=zeror
        rholocf_av_nb_mpi=zeror
        Gcf_nb_mpi=zeror  
        logrhotot_f_mpi=zeror
        
        call MPI_REDUCE( rhof_nb, rhof_nb_mpi, nb*nb*Nc*(2*nwn+1), MPI_REAL8,
     &      MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if(mb_tmt) then
        call MPI_REDUCE( logrholoc_f_nb, logrholoc_f_nb_mpi, nb*nb*(2*nwn+1), 
     &      MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        if(iter.eq.niter) then
        call MPI_REDUCE( rholocf_av_nb, rholocf_av_nb_mpi, nb*nb*(2*nwn+1),
     &      MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE( logrhotot_f, logrhotot_f_mpi, (2*nwn+1),
     &      MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        end if
        end if
   
        if(Gfav) then
        call MPI_REDUCE( Gcf_nb, Gcf_nb_mpi, nb*nb*Nc*(2*nwn+1), MPI_COMPLEX16, 
     &     MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        end if
         
         
!	Normalize the runs over the no. of processors
        if(myrank == 0) then
        rhof_nb=rhof_nb_mpi/dfloat(nprocs)
        if(mb_tmt) then
        logrholoc_f_nb=logrholoc_f_nb_mpi/dfloat(nprocs)
        if(iter.eq.niter) then
        rholocf_av_nb=rholocf_av_nb_mpi/dfloat(nprocs)
        logrhotot_f=logrhotot_f_mpi/dfloat(nprocs)
        end if
        end if
        if(Gfav) Gcf_nb=Gcf_nb_mpi/dfloat(nprocs)
        end if

        !write(*,*)'rhof2',rhof_nb(1,:,1,0)
c  1.--- Impose point group symmetry on rho_aa, rho_bb,rho_ab, rho_ba
        if(myrank==0) then         
         data_nb=zeror
         do i=1,nb
         do j=1,nb    
         do n=-nwn,nwn
	     do ic=1,Nc
	     do igroup=1,ngroup          
              data_nb(i,j,ic,n) = data_nb(i,j,ic,n) +				      				      
     &  	       rhof_nb(i,j,ickequ(ic,igroup),n)/dfloat(ngroup)
 	     end do			      
	     end do
         end do
         end do
         end do        
        !rhof_nb=data_nb
c p-h symmetry
        if(ph_symm)then
            data_nb=zeror
            Gc_temp=zeroc
            do i=1,nb
            do j=1,nb 
            do n=-nwn,nwn
            if(Nc.eq.1)then
            data_nb(i,j,1,n)=data_nb(i,j,1,n)+halfr*(rhof_nb(i,j,1,n)+				
     &        (rhof_nb(i,j,1,-n)))
            Gc_temp(i,j,1,n)=Gc_temp(i,j,1,n)+halfr*(Gcf_nb(i,j,1,n)+                          
     &        (Gcf_nb(i,j,1,-n)))
            else
	        do ic=1,Nc
	          data_nb(i,j,ic,n)=data_nb(i,j,ic,n)+halfr*(rhof_nb(i,j,ic,n)+				
     &        (rhof_nb(i,j,ickdiff(icKpi,ic),-n)))	
	          Gc_temp(i,j,ic,n)=Gc_temp(i,j,ic,n)+halfr*(Gcf_nb(i,j,ic,n)+				
     &        (Gcf_nb(i,j,ickdiff(icKpi,ic),-n)))	
            end do
            end if            
	        end do
            end do
            end do
	  
          rhof_nb=data_nb
          Gcf_nb=Gc_temp 
        data_rho_nb=zeror
        do i=1,nb
        do j=1,nb 
        do n=-nwn,nwn
          data_rho_nb(i,j,n)=halfr*(logrholoc_f_nb(i,j,n)+logrholoc_f_nb(i,j,-n))
        end do
        end do
        end do
        logrholoc_f_nb = data_rho_nb

        endif
      
c 2.1 CHECKING Causality
        do i=1,nb 
        do ic=1,Nc
         do n=-nwn,nwn
            if(rhof_nb(i,i,ic,n).lt. 0.d0) then
c	     rhof_aa(ic,n) = -rhof_aa(ic,n)
	     write(88,*) 'CAUSALITY in orbial',i,wn(n),iter
            end if
         end do
        end do
        end do
        

       if (mb_tmt) then
        if(mb_tmt_ver8)then
          do n=-nwn,nwn
             call DSYEV('V', 'U', nb, logrholoc_f_nb(:,:,n), nb, Eval, workr, lworkr, ierr)
             Evect=logrholoc_f_nb(:,:,n)                                                                                                                           
             rho_temp_nb1=zeror                                                                                                                                        
             do i=1,nb                                                                                                                                                  
             rho_temp_nb1(i,i)=dexp(Eval(i))                                                                                                                                    
             end do                                                                                                                                                           
             call DGEMM('T','N',nb,nb,nb,oner,Evect,nb,rho_temp_nb1,nb,zeror,rho_temp_nb,nb)                                                                          
             call DGEMM('N','N',nb,nb,nb,oner,rho_temp_nb,nb,Evect,nb,zeror,rho_temp_nb1,nb)
             rholoc_typ_nb(:,:,n)=rho_temp_nb1
          end do
          do n=-nwn,nwn
          do ic=1,Nc
          call DGEMM('N','N',nb,nb,nb,oner,rholoc_typ_nb(:,:,n),nb,rhof_nb(:,:,ic,n),nb,zeror,rho_temp_nb,nb)
          rhof_nb(:,:,ic,n)=rho_temp_nb
          end do
          end do  
        else
       	rholoc_typ_nb=zeror
         do i=1,nb
         do j=1,nb   
       	 do n=-nwn,nwn
           rholoc_typ_nb(i,j,n)=dexp(logrholoc_f_nb(i,j,n))
	     end do
         end do
         end do
       if (mb_tmt_ver1) then	  
         do i=1,nb
         do ic=1,Nc
         do n=-nwn,nwn
            rhof_nb(i,i,ic,n)=rhof_nb(i,i,ic,n)*rholoc_typ_nb(i,i,n)
         end do
         end do
         end do
       
       elseif(mb_tmt_ver2) then
         do i=1,nb
         do j=1,nb  
         do ic=1,Nc
         do n=-nwn,nwn
            rhof_nb(i,j,ic,n)=rhof_nb(i,j,ic,n)*rholoc_typ_nb(1,1,n)
         end do
         end do
         end do
         end do

       elseif(mb_tmt_ver6) then
         do i=1,nb
         do j=1,nb  
         do ic=1,Nc
         do n=-nwn,nwn
            rhof_nb(i,j,ic,n)=rhof_nb(i,j,ic,n)*rholoc_typ_nb(i,j,n)
         end do
         end do
         end do
         end do
       
       elseif(mb_tmt_ver7) then
         do i=1,nb
         do ic=1,Nc
         do n=-nwn,nwn
            rhof_nb(i,i,ic,n)=rholoc_typ_nb(i,i,n)
         end do
         end do
         end do
       
       elseif(mb_tmt_ver9) then
         do i=1,nb
         do j=1,nb
         do ic=1,Nc
         do n=-nwn,nwn
            rhof_nb(i,j,ic,n)=rhof_nb(i,j,ic,n)*rholoc_typ_nb(1,1,n)/dfloat(nb)
         end do
         end do
         end do
         end do
       end if !mb_tmt_ver1 
       end if !mb_tmt_ver8
       end if !mb_tmt

c 3. --- Do Hilbert transform to get Gc_aa, Gc_bb,Gc_ab from rho_aa,rho_bb, rho_ab
        Gc_nb_real = zeror
        do i=1,nb
            do j=1,nb
            call HilbertT(Gc_nb(i,j,:,:),Gc_nb_real(i,j,:,:),rhof_nb(i,j,:,:))
            end do
        end do
        if(Gfav) Gc_nb=Gcf_nb
         
        write(*,*) 'rhof',rhof_nb(1,1,1,0),iter

c 5. --- ---	Coarse-graining steps:GcKf=Nc/N\sum_k[Gc^(-1)+\Delta-Eps(k)+Eps_bar(K)]^(-1)
c 5.1: ----     Invert cluster GF Gc(K,w)
c 5.2 now we use that ab=ba

c	invert cluster GF: Gc^(-1) 
        Gcinv_nb=zeroc
        Gc_temp=Gc_nb
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
        Gcinv_nb=Gc_temp
        end if!if (myrank==0)
        if(tables_type.ge.4)then
        call MPI_BCAST(Gcinv_nb,nb*nb*Nc*(2*nwn+1),MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
        end if
!calculate sigma_nb
        if(myrank==0)then
          if(tables_type.ge.4)then
           sigma_nb=zeroc
           do n =-nwn,nwn
	       do ic=1,Nc
             do i=1,nb
             do j=1,nb 
             sigma_nb(i,j,ic,n)=-Epsb_nb(i,j,ic)-Gamma_old_nb(i,j,ic,n)-Gcinv_nb(i,j,ic,n)
             if(i.eq.j) sigma_nb(i,i,ic,n)=sigma_nb(i,i,ic,n)+wn(n)
             end do
             end do
           end do
           end do 
          else !if tales_type <=3
           !sigma_nb=zeroc
           do n =-nwn,nwn
	       do ic=1,Nc
             do i=1,nb
             do j=1,nb 
             !sigma_nb(i,j,ic,n)=-Gamma_old_nb(i,j,ic,n)-Gcinv_nb(i,j,ic,n)
             !if(i.eq.j) sigma_nb(i,i,ic,n)=sigma_nb(i,i,ic,n)+wn(n)-aa(i)*Epsbar(ic)
             !if(i.ne.j) sigma_nb(i,j,ic,n)=sigma_nb(i,j,ic,n)-ab*Epsbar(ic)
             sigma_nb(i,j,ic,n)=(Gloc_inv_nb(i,j,ic,n)-Gcinv_nb(i,j,ic,n))*(1.d0-dfac)+dfac*sigma_nb(i,j,ic,n)
             end do
             end do
           end do
           end do 
          end if !if tables_type >=4
        end if!if myrank==0 
!Distribute the summaiton for coarse-grain into all processors
c.5.2---Construct inverted lattice GF GcKf(k,K,w)GcKf=[Gc^(-1)+\Delta-Eps(k)+Eps_bar(K)]^(-1) which is then coarse-grained
         if(tables_type.ge.4)then
          GcKf_nb=zeroc
          do ic=1,Nc
            do n=-nwn,nwn
            !do i=1,nwsc
            istart = 1+nwsc/nprocs*myrank
            if(myrank < nprocs-1) iend = nwsc/nprocs*(myrank+1)
            if(myrank == nprocs-1) iend = nwsc 
            do i=istart,iend
             Gc_temp=zeroc
             Eps_temp=zeroc
             rv1(1) = Kc(1,ic)+kt(1,i)
             rv1(2) = Kc(2,ic)+kt(2,i)
             rv1(3) = Kc(3,ic)+kt(3,i)
             Eps_temp(1,2)=cmplx(1+cos(rv1(1)),-sin(rv1(1)))
             Eps_temp(1,3)=cmplx(1+cos(rv1(2)),-sin(rv1(2)))
             Eps_temp(1,5)=cmplx(1+cos(rv1(3)),-sin(rv1(3)))
             Eps_temp(2,4)=cmplx(1+cos(rv1(2)),-sin(rv1(2)))
             Eps_temp(2,6)=cmplx(1+cos(rv1(3)),-sin(rv1(3)))  
             Eps_temp(3,4)=cmplx(1+cos(rv1(1)),-sin(rv1(1)))
             Eps_temp(3,7)=cmplx(1+cos(rv1(3)),-sin(rv1(3)))
             Eps_temp(4,8)=cmplx(1+cos(rv1(3)),-sin(rv1(3)))
             Eps_temp(5,6)=cmplx(1+cos(rv1(1)),-sin(rv1(1)))
             Eps_temp(5,7)=cmplx(1+cos(rv1(2)),-sin(rv1(2)))
             Eps_temp(6,8)=cmplx(1+cos(rv1(2)),-sin(rv1(2)))
             Eps_temp(7,8)=cmplx(1+cos(rv1(1)),-sin(rv1(1)))
             do jc=2,8
                do j=1,jc-1
                   Eps_temp(jc,j)=conjg(Eps_temp(j,jc))
                end do !j
             end do !jc
             Eps_temp=Eps_temp/4.d0  !unit:4t=1
             do jc=1,nb
                 do j=1,nb
                    Gc_temp(jc,j,ic,n)=Gcinv_nb(jc,j,ic,n)+Gamma_old_nb(jc,j,ic,n)-Eps_temp(jc,j)+Epsb_nb(jc,j,ic)
                 end do
             end do  
             call zgetrf(nb,nb,Gc_temp(:,:,ic,n),nb,ipvtnb,info)
	         if (info.eq.0)then
	         call zgetri(nb,Gc_temp(:,:,ic,n),nb,ipvtnb,work,lwork,info)
             else
               print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
	         stop
             end if
             GcKf_nb(:,:,ic,n)=GcKf_nb(:,:,ic,n)+Gc_temp(:,:,ic,n)
            end do !kt_i-loop
            end do !wn-loop
          end do !ic-loop
          call MPI_REDUCE(GcKf_nb, GcKf_nb_mpi, nb*nb*Nc*(2*nwn+1),
     &         MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
! normalize
          if(myrank==0) GcKf_nb = GcKf_nb_mpi/real(nwsc,kind)
         else!if talbes_type<=4
          if(myrank==0)then 
          write(*,*) 'using sigma to coarse grain'
          !construct sigma_nb_loc
          if(embedding1)then
          if(iter==1) sigma_nb_loc(:,:,:)=sigma_nb(:,:,1,:)    
          do i=1,nb
             do j=1,nb
                do n=-nwn,nwn
                   Gc_nb_loc(i,j,n)=sum(Gc_nb(i,j,:,n))/real(Nc)
                   !Gsc_nb_loc(i,j,n)=sum(GcKfsc_nb(i,j,:,n))/real(Nc)
                   Gsc_nb_loc(i,j,n)=sum(GcKf_nb(i,j,:,n))/real(Nc)
                end do
             end do
          end do
          do n=-nwn,nwn    
             call zgetrf(nb,nb,Gc_nb_loc(:,:,n),nb,ipvtnb,info)
             if (info.eq.0)then
                call zgetri(nb,Gc_nb_loc(:,:,n),nb,ipvtnb,work,lwork,info)
             else
                print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
                stop
             end if
             call zgetrf(nb,nb,Gsc_nb_loc(:,:,n),nb,ipvtnb,info)
             if (info.eq.0)then
                call zgetri(nb,Gsc_nb_loc(:,:,n),nb,ipvtnb,work,lwork,info)
             else
                print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
                stop
             end if
          end do!n loop
          Gsc_nb_loc=Gsc_nb_loc+sigma_nb_loc
          !sigma_nb_loc=Gsc_nb_loc-Gc_nb_loc
          sigma_nb_loc=(Gsc_nb_loc-Gc_nb_loc)*(1.d0-dfac)+dfac*sigma_nb_loc
         write(128,'(i4,20(1x,f11.4))')iter,dimag(Gc_nb_loc(:,:,36)),dimag(Gsc_nb_loc(:,:,36)),dimag(sigma_nb_loc(:,:,36))
         write(129,'(i4,20(1x,f11.4))')iter,(dimag(sigma_nb(:,:,i,36)),i=1,Ncw)
         do n=-nwn,nwn
         write(130,'(20(1x,f11.4))')wn(n),dimag(sigma_nb_loc(:,:,n)),sum(dimag(sigma_nb(1,1,:,n)))/real(Nc)
         end do
         write(130,*)
          elseif(embedding2)then
          do i=1,nb
             do j=1,nb
                do n=-nwn,nwn
                   sigma_nb_loc(i,j,n)=sum(sigma_nb(i,j,:,n))/real(Nc)
                end do
             end do
          end do
          end if!if embedding1 
          GcKf_nb = zeroc
          do ic=1,Nc
            do n=-nwn,nwn
            do i=0,e_div
             e=estart+dfloat(i)*estep
             do jc=1,nb
                 do j=1,nb
                  if(init_sigma)then
                   if(embedding1.or.embedding2)then
                     if(jc.eq.j) Gc_temp(jc,jc,ic,n)=wn(n)-sigma_nb_loc(jc,jc,n)-aa(jc)*e
                     if(jc.ne.j) Gc_temp(jc,j,ic,n)=-sigma_nb_loc(jc,j,n)-ab*e
                   else
                     if(jc.eq.j) Gc_temp(jc,jc,ic,n)=wn(n)-sigma_nb(jc,jc,ic,n)-aa(jc)*e
                     if(jc.ne.j) Gc_temp(jc,j,ic,n)=-sigma_nb(jc,j,ic,n)-ab*e
                   end if!if embedding
                  else
                    Gc_temp(jc,j,ic,n)=Gcinv_nb(jc,j,ic,n)+Gamma_old_nb(jc,j,ic,n)
                    if(jc.eq.j) Gc_temp(jc,jc,ic,n)=Gc_temp(jc,jc,ic,n)-aa(jc)*(e-Epsbar(ic))
                    if(jc.ne.j) Gc_temp(jc,j,ic,n)=Gc_temp(jc,j,ic,n)-ab*(e-Epsbar(ic))
                  endif!if init_sigma
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
          end if !if myrank==0
         end if !if tables_type>=4
         if(myrank==0)then  
             
c CAUSALITY 2
             do i=1,nb
             do ic=1,Nc
             do n=-nwn,nwn
             if(aimag( GcKf_nb(i,i,ic,n)).gt. 1.d-8) then
c	      GcKf_nb(i,i,ic,n) =cmplx(real(GcKf_nb(i,i,ic,n)),-aimag(GcKf_nb(i,i,ic,n)))
	      write(99,*) 'CAUSALITY in GLOC_',i,ic,wn(n),iter,GcKf_nb(i,i,ic,n)
!         GcKf_nb(i,i,ic,n)=cmplx(real(GcKf_nb(i,i,ic,n)),-1.d-2)
	         end if
             end do
             end do
             end do
        

c  invert local lattice GF which is needed to update Delta: Delta_new=(w-Sigma-Gc^-1)

        Gloc_inv_nb=zeroc
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

    
c 7. ---calculate new Delta(w)

   
         !Gamma_new_nb=zeroc
                   
         !Gamma_new_nb = Gamma_old_nb+(1.0D0-dfac)*(Gcinv_nb-Gloc_inv_nb)
         write(601,*)'Gc_nb'
         write(601,'(6(f9.6,2x))')Gc_nb(1,1,2,44),Gcinv_nb(1,1,2,44)
         write(601,*)'Gc_bar'
         write(601,'(6(f9.6,2x))')GcKf_nb(1,1,2,44),Gloc_inv_nb(1,1,2,44)
         write(601,*)'sigma_nb'
         write(601,'(6(f9.6,2x))')sigma_nb(1,1,2,44),sum(sigma_nb(1,1,:,44))/real(Nc)
! Calculate convergency criteria
        testl=0.d0      
        test=0.d0
        test2=0.d0        
         
         do ic=1,Nc
	     do n=-nwn,nwn
         do i=1,nb
         do j=1,nb    
	        testl = testl + dabs(dimag(Gamma_old_nb(i,j,ic,n)))
            test = test + dabs(dimag(Gamma_new_nb(i,j,ic,n))) 
            test2 = test2 + dabs(dimag(Gamma_old_nb(i,j,ic,n))
     &                    - dimag(Gamma_new_nb(i,j,ic,n)))

          	    
         end do
         end do
         end do 
         end do

         !For debug print out Gamma in the last iteration
         if(iter.eq.niter)then
        write(825,*)'Gamma_old im'
        do i=1,nb
           write(825,'(20(f7.4,1x))')-dimag(Gamma_old_nb(i,:,1,0))
        end do

        write(825,*)'Gamma_new im'
        do i=1,nb
           write(825,'(20(f7.4,1x))')-dimag(Gamma_new_nb(i,:,1,0))
        end do

        write(825,*)'Gamma_old re'
        do i=1,nb
           write(825,'(20(f7.4,1x))')-dreal(Gamma_old_nb(i,:,1,0))
        end do

        write(825,*)'Gamma_new re'
        do i=1,nb
           write(825,'(20(f7.4,1x))')-dreal(Gamma_new_nb(i,:,1,0))
        end do
        end if 
       !Gamma_old_nb = Gamma_new_nb
       end if!if myrank==0


       
       if(myrank==0)then
c 8.--- update Gscript: for this we need to invert matrixGsc=[w*1-Ebar(K)-Gamma(K)]^-1




c 2.calculate Gscript_aa,bb,ab components
c***************************************************************************************
c	form G_script to invert
c***************************************************************************************

        GcKfsc_nb=zeroc
        Gc_temp=zeroc

          do ic=1,Nc
          do n=-nwn,nwn
              do i=1,nb
              do j=1,nb
              if(tables_type.ge.4)then
              Gc_temp(i,j,ic,n)=-Epsb_nb(i,j,ic)-Gamma_old_nb(i,j,ic,n)
              if(i.eq.j) then
                 Gc_temp(i,i,ic,n)=Gc_temp(i,i,ic,n)+wn(n) 
              end if
              else     
              !Gc_temp(i,j,ic,n)=-Gamma_old_nb(i,j,ic,n)    
              !if(i.eq.j) then
                 !Gc_temp(i,i,ic,n)=Gc_temp(i,i,ic,n)+wn(n)-aa(i)*Epsbar(ic)
              !end if
              if(embedding1.or.embedding2)then
                 Gc_temp(i,j,ic,n)=Gloc_inv_nb(i,j,ic,n)+sigma_nb_loc(i,j,n)   
              else
                 Gc_temp(i,j,ic,n)=Gloc_inv_nb(i,j,ic,n)+sigma_nb(i,j,ic,n) 
              end if  
              end if !talbes_type.eq.4
              end do
              end do
              call zgetrf(nb,nb,Gc_temp(:,:,ic,n),nb,ipvtnb,info)
	          if (info.eq.0)then
	          call zgetri(nb,Gc_temp(:,:,ic,n),nb,ipvtnb,work,lwork,info)
              else
               print*, 'info=/0 Matrix Inversion Failed in put: Info=',info
	          stop
              end if
              GcKfsc_nb(:,:,ic,n)=Gc_temp(:,:,ic,n)
	      end do
          end do
!invert GcKfsc_nb to calculate Gamma_new,temporarily restored in Gloc_inv_nb
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
        do i=1,nb
           do j=1,nb
           do n=-nwn,nwn
           do ic=1,Nc   
              if(i==j)Gamma_new_nb(i,i,ic,n)=wn(n)-Gloc_inv_nb(i,i,ic,n)-Epsb_nb(i,i,ic)
              if(i/=j)Gamma_new_nb(i,j,ic,n)=-Gloc_inv_nb(i,j,ic,n)-Epsb_nb(i,j,ic)
           end do
           end do
           end do
        end do
! Calculate convergency criteria
        testl=0.d0      
        test=0.d0
        test2=0.d0        
         
         do ic=1,Nc
	     do n=-nwn,nwn
         do i=1,nb
         do j=1,nb    
	        testl = testl + dabs(dimag(Gamma_old_nb(i,j,ic,n)))
            test = test + dabs(dimag(Gamma_new_nb(i,j,ic,n))) 
            test2 = test2 + dabs(dimag(Gamma_old_nb(i,j,ic,n))
     &                    - dimag(Gamma_new_nb(i,j,ic,n)))

          	    
         end do
         end do
         end do 
         end do
       Gamma_old_nb = Gamma_new_nb

         if(Gscript_pole)then
          do ic =1,Nc
           sum1(ic) = zeror
           do n=-nwn,nwn
              sum1(ic)=sum1(ic)+dwn(n)*(-1.d0/pi)*dimag(GcKfsc_nb(1,1,ic,n))
           end do
!         write(*,*) 'Norm_GcKfsc-before taking pole=',sum1(ic), 'ic=',ic
        end do        
 
          do ic=1,Nc
            do n=-nwn,nwn
              m_r1(ic,n)=real(wn(n)-Epsbar(ic)-Gamma_new_nb(1,1,ic,n))
            end do
          end do

          do ic=1,Nc    ! ic loop starts
            min_w=abs(m_r1(ic,-nwn))
            x_pole = -nwn
            do n=-nwn,nwn
              if(min_w.gt.abs(m_r1(ic,n))) then ! finding minimum starts
                min_w=abs(m_r1(ic,n))
                x_pole=n
              end if    ! finding minimum end        
            end do ! n


          if(abs(1.-sum1(ic)) .gt. 0.01 .and. (-1.d0/pi)*aimag(Gamma_new_nb(1,1,ic,x_pole)) .lt. afact*dwn(x_pole)) then
                
            write(6,*) 'Pole of cell', ic, 'occurs at', x_pole, wn(x_pole),abs(aimag(Gamma_new_nb(1,1,ic,x_pole))),sum1(ic)
            write(42,*) 'Pole of cell', ic, 'occurs at', x_pole, wn(x_pole),abs(aimag(Gamma_new_nb(1,1,ic,x_pole))),sum1(ic)
            do n=-nwn,nwn
               if(n .eq. x_pole) then ! Calculating Gscript when pole occurs. Approximate using delta function
             GcKfsc_nb(1,1,ic,n) = cmplx(0.d0,-pi/dwn(x_pole))
                        
            else if (n .ne. x_pole) then                
              GcKfsc_nb(1,1,ic,n) = cmplx(1.d0/(wn(n)-real(wn(x_pole))),0.d0)
             !GcKfsc_nb(1,1,ic,n) = cmplx(1.d0/(wn(n)-Epsbar(ic)-real(Gamma_new_nb(1,1,ic,n))),0.d0) !                       
            end if      
           end do
          end if
         end do
        end if

c   go to geod.f to make it in real space
c *************************************************************
c 9.    write output
c  calculate local average doss \sum_K(G_aa(K)+G_bb(K)+...)    
        rhol_nb=zeror
         do i=1,nb
         do j=1,nb    
         do n=-nwn,nwn
	     do ic=1,Nc
	     rhol_nb(i,j,n) = rhol_nb(i,j,n)-dimag(GcKf_nb(i,j,ic,n))/pi
	     end do
         rhol_nb(i,j,n)=rhol_nb(i,j,n)/dfloat(Nc)
	     end do
         end do
         end do

c       
        rho_av=zeror
        rhof_av=zeror
        
	     do n=-nwn,nwn
	     do ic=1,Nc
         do i=1,nb
	        rho_av(n) = rho_av(n)-dimag(GcKf_nb(i,i,ic,n))
            rhof_av(n)=rhof_av(n)+rhof_nb(i,i,ic,n)
         end do
         end do
	      rho_av(n) = rho_av(n)/(dfloat(Nc)*pi)
          rhof_av(n)=rhof_av(n)/dfloat(Nc)       
	     end do
            
         write(124,'(i4,1x,f9.6,1x,f9.6)')iter,rho_av(0)/dfloat(nb),rhof_av(0)/dfloat(nb)
         write(188,'(i4,20(1x,f9.6))')iter,(rhol_nb(i,i,0),i=1,nb)
         write(125,'(i4,1x,4(f9.6,1x))')iter,dimag(Gamma_old_nb(1,1,2,19)),dimag(GcKf_nb(1,1,2,19)),dimag(Gcf_nb(1,1,2,19))
         write(600,'(i4,20(1x,f9.6))')iter,-dimag(Gamma_old_nb(:,:,2,44)),-dimag(Gamma_old_nb(:,:,2,-60))
         write(126,'(i4,20(1x,f9.6))')iter,dimag(sigma_nb_loc(:,:,0))


         write(*,*) 'rho_av,rhof_av',rho_av(0)/dfloat(nb),rhof_av(0)/dfloat(nb),iter  
         
         if (iter.eq. niter) then
             if(mb_tmt) rhotot_typ=dexp(logrhotot_f)
         do n=-nwn,nwn    
          write(132,'(f14.6,20(1x,f9.6))')wn(n),(rhol_nb(i,i,n),i=1,nb)
     
          write(133,'(f14.6,20(1x,f9.6))')wn(n),(rholocf_av_nb(i,i,n),i=1,nb)

          write(422,'(f9.4,1x,f9.6,1x,f9.6,1x,f9.6)')wn(n),rho_av(n)/dfloat(nb),rhof_av(n)/dfloat(nb)
!     &          ,rhotot_typ(n)/dfloat(nb)
          write(821,'(f14.6,20(1x,f9.6))')wn(n),(sigma_nb_loc(i,i,n),i=1,nb)
          write(824,'(f14.6,40(1x,f9.6))')wn(n),(dreal(sigma_nb(1,1,ic,n)),ic=1,Ncw)
          write(826,'(f14.6,40(1x,f9.6))')wn(n),(dimag(sigma_nb(1,1,ic,n)),ic=1,Ncw)
          write(822,'(f14.6,20(1x,f9.6))')wn(n),(sum(Gamma_old_nb(i,i,:,n))/real(Nc),i=1,nb)
          write(827,'(f14.6,40(1x,f9.6))')wn(n),(dreal(Gamma_old_nb(1,1,ic,n)),ic=1,Ncw)
          write(828,'(f14.6,40(1x,f9.6))')wn(n),(dimag(Gamma_old_nb(1,1,ic,n)),ic=1,Ncw)
         end do
       
       !For debug print out GcKf and rhof in the last iteration
        write(823,*)'GcKf'
        do i=1,nb
           write(823,'(20(f7.4,1x))')-dimag(GcKf_nb(i,:,1,0))/pi
        end do
        write(823,*)'rhof'
        do i=1,nb
           write(823,'(20(f7.4,1x))')rhof_nb(i,:,1,0)
        end do
        
! output Gscript_w
        open(125,file='gscript_w.dat',status='unknown')
            do ic=1,Ncw
                do n=-nwn,nwn
                    write(125,*) wn(n), -dimag(GcKfsc_nb(1,1,ic,n))
                end do
                write(125,*) '   '
            end do
            close(125)

        end if

c****check sum rule for G(i,i,w) in the last iteration
        if(iter.eq.niter)then
            if(.not. mb_tmt)then
            r=zeror
            do i=1,8
            do n=-nwn,nwn
                r(i)=r(i)+rhol_nb(i,i,n)*dwn(n)
            enddo
            end do
            else
            r=zeror
            do i=1,8
            do n=-nwn,nwn
                r(i)=r(i)+rholocf_av_nb(i,i,n)*dwn(n)
            enddo
            end do
            end if
            write(6,*)'G_11(i,i,w)_sum_rule=',r(1)
           ! write(6,*)'G_22(i,i,w)_sum_rule=',r(2)
           ! write(6,*)'G_33(i,i,w)_sum_rule=',r(3)
           ! write(6,*)'G_44(i,i,w)_sum_rule=',r(4)
           ! write(6,*)'G_55(i,i,w)_sum_rule=',r(5)
           ! write(6,*)'G_66(i,i,w)_sum_rule=',r(6)
           ! write(6,*)'G_77(i,i,w)_sum_rule=',r(7)
           ! write(6,*)'G_88(i,i,w)_sum_rule=',r(8)
            write(42,*)'G_11(i,i,w)_sum_rule=',r(1)
           ! write(42,*)'G_22(i,i,w)_sum_rule=',r(2)
           ! write(42,*)'G_33(i,i,w)_sum_rule=',r(3)
           ! write(42,*)'G_44(i,i,w)_sum_rule=',r(4)
           ! write(42,*)'G_55(i,i,w)_sum_rule=',r(5)
           ! write(42,*)'G_66(i,i,w)_sum_rule=',r(6)
           ! write(42,*)'G_77(i,i,w)_sum_rule=',r(7)
           ! write(42,*)'G_88(i,i,w)_sum_rule=',r(8)


        endif
c****check sum rule for GsKfsc
        
	     if(.true.)then
          do ic =1,Ncw
	        sum1(ic) = zeror
                do n=-nwn,nwn
                   sum1(ic)=sum1(ic)+dwn(n)*(-1.d0/pi)*dimag(GcKfsc_nb(1,1,ic,n))
                end do
	        write(*,*) 'Norm_GcKfsc=',sum1(ic), 'ic=',ic
	        write(42,*) 'Norm_GcKfsc=',sum1(ic), 'ic=',ic
          end do

          if(iter.eq.niter)then
            r3=zeror
            do n=-nwn,nwn
            do ic=1,Nc
            r3=r3+(-1.d0/pi)*dimag(GcKfsc_nb(1,1,ic,n))*dwn(n)
            enddo
            enddo
            r3=r3/dfloat(Nc)
            write(6,*)'GcKfsc_11_sum_rule=',r3
            write(42,*)'GcKfsc_11_sum_rule=',r3
          endif
 	     endif

       if(iter.ge.2) then					      
          r1=dabs(testl-test)/(test*dfloat(Nc))
          write(6,*) ' test= ',r1,'test2= ',test2/(test*dfloat(Nc)) 
 	  write(6,*) '**********************************************************'		      
          write(42,*) ' test= ',r1 ,'test2= ',test2/(test*dfloat(Nc))
 	  write(42,*) '**********************************************************'		      
        end if
        end if !if myrank==0 
       if(tables_type.le.3)then
       call MPI_BCAST(Gamma_old_nb,nb*nb*Nc*(2*nwn+1),MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(GcKfsc_nb,nb*nb*Nc*(2*nwn+1),MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
       else
       call MPI_BCAST(Gamma_old_nb,nb*nb*Nc*(2*nwn+1),MPI_COMPLEX16,0,MPI_COMM_WORLD,ierr)
       end if
       end if !Liz_embedding 
      END IF !MB1   
ccccccccccccccccccccccccccccccccccccccccccc
c      end MB
ccccccccccccccccccccccccccccccccccccccccccc      
      
      if(myrank==0) call profend('output') 
        return
        end
        
c        include  'quadpack.f90 '
