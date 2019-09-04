c***************************************************************************************
        subroutine measure
c***************************************************************************************

c***************************************************************************************
	use Global
	use omp_lib
	implicit none
c***************************************************************************************
        integer i,j,n,m,p,ic,jc,jcr,info
        real*8 :: rnb(nb,nb,-nwn:nwn),rho_temp,rho_temp_nb(nb,nb),
     &    rho_temp_nb1(nb,nb),workr(100),Eval(nb),Evect(nb,nb)
        complex*16, allocatable :: csum_nb(:,:,:,:)
        integer lworkr       
        lworkr=100
c       This program adds the present values of g, sus ......
c	to the accumulators for a single bin.  
c       Each entry in the G(i,j,w) matrix contributes to
c	the average of G(i-j,w).
c***************************************************************************************
       if(myrank==0) call profstart('measure')


c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c            I:  FT to K-space G_ij->G(K)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c cccccccccccccccccccccccccccccccccccccccccc
c            MB
cccccccccccccccccccccccccccccccccccccccccccc

      if ( MB1) then
      if(Liz_embedding)then
        if(.not. mb_tmt)then    
          Gc_liz_meas=Gc_liz_meas+Gc_liz
        else
          if(mb_tmt_ver1)then
            do n=-nwn,nwn
               Gc_liz_meas(:,:,n)=Gc_liz_meas(:,:,n)+Gc_liz(:,:,n)/(-1.d0/pi*dimag(Gc_liz(1,1,n)))
            end do
            if(iter==niter)then
               GIPR2_meas=GIPR2_meas+GIPR2
               logGIPR2_meas=logGIPR2_meas+dlog(GIPR2)
            end if    
          elseif(mb_tmt_ver2)then
            !rho_temp_nb(1,1)=zeror
            !do i=1,lizsize
            !   rho_temp_nb(1,1)=rho_temp_nb(1,1)+(-1.d0/pi*dimag(Gc_liz(1,1,n)))
            !end do
            !rho_temp_nb(1,1)=rho_temp_nb(1,1)/real(lizsize)   
          end if
        end if
        !if(myrank==0)then
        !  write(*,*)-1.d0/pi*dimag(Gc_liz(1,1,0)),logrholoc_meas_nb(1,1,0)
        !end if
      else! if not liz_embedding  
c 1:  construct nb,nb components of G(i,j,v) 
c     one can construct nb*Ncx nb*Nc matrix instead of introducing components,
      
      allocate(csum_nb(nb,nb,Nc,-nwn:nwn),stat=info)
      if(info.ne.0)then
          write(*,*) 'allocat failure : csum_nb'
          stop
      end if

        Gdff_nb = zeroc

      do n=-nwn,nwn
      do i=1,nb
      do j=1,nb    
	  do ic=1,Nc
	  do jc=1,Nc
	
	       Gdff_nb(i,j,icrdiff(ic,jc),n)=
     &     Gdff_nb(i,j,icrdiff(ic,jc),n)+Gcv_nb(i,j,ic,jc,n)/dfloat(Nc)
            
      end do
      end do  
      end do
	  end do
      end do  !end w loop
c.   2: FT to K-space to get G_aa(K),G_bb(K),G_ab(K),G_ba(K)
      
      csum_nb=zeroc
      GdKK_nb=zeroc
	  do n=-nwn,nwn
      do i=1,nb
      do j=1,nb          
	  do ick=1,Nc						      
      do icr=1,Nc				      
         csum_nb(i,j,ick,n)=csum_nb(i,j,ick,n)+FTCoefs_R_to_K(ick,icr)*Gdff_nb(i,j,icr,n)	   
	  end do
c for some Nc, it happens that at the boundary at +/- nw_max, causality is broken with 	aimag(csum).ge. 0.d0
c here we avoid this numerical problem by making it -depsi  
            if((i.eq.j .and. dimag(csum_nb(i,i,ick,n)) >= 0.d0)) then  
c	      write(*,*) n,ick,csum
                csum_nb(i,i,ick,n) = dcmplx(dreal(csum_nb(i,i,ick,n)),-depsi)
	        end if
            GdKK_nb(i,j,ick,n)=csum_nb(i,j,ick,n)
	  end do  !end K loop
      end do  !end j loop
      end do  !end i loop      
	  end do  !end w loop

         if (.not. mb_tmt) then
             do n=-nwn,nwn	
             do i=1,nb
             do j=1,nb     
	         do ick=1,Nc						      
	       rho_meas_nb(i,j,ick,n)=rho_meas_nb(i,j,ick,n)+(-(1.0/pi)*(dimag(GdKK_nb(i,j,ick,n)))) 
	         end do
	         end do
             end do
             end do
             if(Gfav) Gc_meas_nb=Gc_meas_nb+GdKK_nb
	     else
           if(iter.eq.niter)then !calculate the average density of states
               do i=1,nb
               do j=1,nb    
               do n=-nwn,nwn
                   do ic=1,Nc
                       rholoc_av_meas_nb(i,j,n)=rholoc_av_meas_nb(i,j,n)-((1./pi)*aimag(Gcv_nb(i,j,ic,ic,n)))
                   end do
               end do
               end do
               end do
             
               do n=-nwn,nwn
               do ic=1,Nc
               rho_temp=zeror
               do i=1,nb
                  rho_temp=rho_temp-((1./pi)*dimag(Gcv_nb(i,i,ic,ic,n)))
               end do
                  logrhotot_meas(n)=logrhotot_meas(n)+dlog(rho_temp)
               end do
               end do
           endif

c          if ver_1 (my)
          if (mb_tmt_ver1) then !Hanna's original ansatz	  
c     1rst construct local for each band
          rnb=zeror          
          do i=1,nb
	      do n=-nwn,nwn
             do ic=1,Nc          
             rnb(i,i,n)=rnb(i,i,n)-((1./pi)*dimag(Gcv_nb(i,i,ic,ic,n)))
             end do
             rnb(i,i,n)=rnb(i,i,n)/dfloat(Nc)
          end do 
          end do 
c   2nd construct non-local matrix
         !diagonal part 
          do i=1,nb
           do n=-nwn,nwn		  
	       do ick=1,Nc						      
        rho_meas_nb(i,i,ick,n)=rho_meas_nb(i,i,ick,n)+dabs(-(1.0/pi)*(dimag(GdKK_nb(i,i,ick,n)))/rnb(i,i,n) ) 
	       end do
	       end do
          end do

         !off-diagonal part
          do i=1,nb
            do j=1,nb
                if (i.ne.j) then
                  do n=-nwn,nwn		  
	              do ick=1,Nc						      
                  rho_meas_nb(i,j,ick,n)=rho_meas_nb(i,j,ick,n)+(-(1.0/pi)*(dimag(GdKK_nb(i,j,ick,n))) ) 
	              end do
	              end do
                end if
            end do
          end do

c   3rd construct log-prefactor for each band    
         do i=1,nb 
	     do n=-nwn,nwn        
            do ic=1,Nc
             logrholoc_meas_nb(i,i,n) = logrholoc_meas_nb(i,i,n)
     &       +dlog(dabs(-dimag(Gcv_nb(i,i,ic,ic,n))/pi))
            end do  
         end do
         end do          
         end if  !mb_tmt_ver1

          if (mb_tmt_ver2) then !rho_tot_typ basis invariance	  
c     1rst construct local for each band
          rnb=zeror          
	  do n=-nwn,nwn
          do i=1,nb
             do ic=1,Nc          
             rnb(1,1,n)=rnb(1,1,n)-((1./pi)*dimag(Gcv_nb(i,i,ic,ic,n)))
             end do
          end do 
          rnb(1,1,n)=rnb(1,1,n)/dfloat(Nc)
          end do 
c   2nd construct non-local matrix

          do n=-nwn,nwn
            do i=1,nb
              do j=1,nb
                do ick=1,Nc
                  rho_meas_nb(i,j,ick,n)=rho_meas_nb(i,j,ick,n)+(-(1.0/pi)*(dimag(GdKK_nb(i,j,ick,n)))/rnb(1,1,n) ) 
	        end do
	      end do
            end do
          end do
c   3rd construct log-prefactor for each band    
	 do n=-nwn,nwn        
         do ic=1,Nc
           rho_temp=zeror
           do i=1,nb
             rho_temp=rho_temp +dabs((1./pi)*dimag(Gcv_nb(i,i,ic,ic,n)))
           end do
           logrholoc_meas_nb(1,1,n) = logrholoc_meas_nb(1,1,n)+dlog(rho_temp)
         end do
         end do          
         end if  !mb_tmt_ver2


c          if ver_6 (new)
         if (mb_tmt_ver6) then	!Yi's ansatz  
c     1rst construct local for each band          
	      rnb=zeror          
          do i=1,nb
          do j=1,nb
             if(i.eq.j) then   
	            do n=-nwn,nwn
                do ic=1,Nc          
                rnb(i,i,n)=rnb(i,i,n)-((1./pi)*dimag(Gcv_nb(i,i,ic,ic,n)))
                end do !loop ic
                rnb(i,i,n)=rnb(i,i,n)/dfloat(Nc)
                end do !loop n
             else
                do n=-nwn,nwn
                do ic=1,Nc          
                rnb(i,j,n)=rnb(i,j,n)+dabs(((1./pi)*dimag(Gcv_nb(i,j,ic,ic,n))))
                end do !loop ic
                rnb(i,j,n)=rnb(i,j,n)/dfloat(Nc)
                end do !loop n
             end if
          end do !loop j
          end do !loop i
         
c   2nd construct non-local matrix
         !diagonal part 
          do i=1,nb
           do n=-nwn,nwn		  
	       do ick=1,Nc						      
        rho_meas_nb(i,i,ick,n)=rho_meas_nb(i,i,ick,n)+dabs(-(1.0/pi)*(dimag(GdKK_nb(i,i,ick,n)))/rnb(i,i,n) ) 
	       end do
	       end do
          end do

         !off-diagonal part
          do i=1,nb
            do j=1,nb
                if (i.ne.j) then
                  do n=-nwn,nwn		  
	              do ick=1,Nc						      
        rho_meas_nb(i,j,ick,n)=rho_meas_nb(i,j,ick,n)+(-(1.0/pi)*(dimag(GdKK_nb(i,j,ick,n)))/rnb(i,j,n) ) 
	              end do
	              end do
                end if
            end do
          end do

c   3rd construct log-prefactor for each band    
         do i=1,nb
         do j=1,nb    
	     do n=-nwn,nwn        
            do ic=1,Nc
             logrholoc_meas_nb(i,j,n) = logrholoc_meas_nb(i,j,n)
     &       +dlog(dabs(-dimag(Gcv_nb(i,j,ic,ic,n))/pi))
            end do
         end do   
         end do
         end do  
         end if !mb_tmt_ver6 

c          if ver_7 
         if (mb_tmt_ver7) then	!local ansatz  
c     1rst construct local for each band          
	      rnb=zeror          
          do i=1,nb
          do j=1,nb
                do n=-nwn,nwn
                do ic=1,Nc          
                rnb(i,j,n)=rnb(i,j,n)+(-(1./pi)*dimag(Gcv_nb(i,j,ic,ic,n)))
                end do !loop ic
                rnb(i,j,n)=rnb(i,j,n)/dfloat(Nc)
                end do !loop n
          end do !loop j
          end do !loop i
         
c   2nd construct non-local matrix

          do i=1,nb
            do j=1,nb
                  do n=-nwn,nwn		  
	              do ick=1,Nc						      
                  rho_meas_nb(i,j,ick,n)=rho_meas_nb(i,j,ick,n)+rnb(i,j,n) 
	              end do
	              end do
            end do
          end do

c   3rd construct log-prefactor for each band    
         do i=1,nb
         do j=1,nb    
	     do n=-nwn,nwn        
            do ic=1,Nc
             logrholoc_meas_nb(i,j,n) = logrholoc_meas_nb(i,j,n)
     &       +dlog(dabs(-dimag(Gcv_nb(i,j,ic,ic,n))/pi))
            end do
         end do   
         end do
         end do  
         end if !mb_tmt_ver7         
          
        if (mb_tmt_ver8) then !Matrix version local basis invariant
c     1rst calculate (rho(w))^(-1/2)
          rnb=zeror          
	      do n=-nwn,nwn
          do j=1,nb
          do i=1,nb
             do ic=1,Nc          
             rnb(i,j,n)=rnb(i,j,n)-((1./pi)*dimag(Gcv_nb(i,j,ic,ic,n)))
             end do
             rnb(i,j,n)=rnb(i,j,n)/dfloat(Nc)
          end do
          end do 
          end do
          do n=-nwn,nwn
            call DSYEV('V', 'U', nb, rnb(:,:,n), nb, Eval, workr, lworkr, ierr)
            Evect(:,:)=rnb(:,:,n)
            rho_temp_nb1=zeror
            do i=1,nb
              if(Eval(i).lt.0.d0)then
                 if(myrank==0 .and. dabs(Eval(i)).gt.1.d-8)write(321,*)
     &'rho_ave negative,i=',i,'Eval=',Eval(i),'iter=',iter,'nmeas=',nmeas,'n=',n
                 Eval(i)=1.d-8
              end if
              rho_temp_nb1(i,i)=1.d0/dsqrt(Eval(i))
            end do
            call DGEMM('T','N',nb,nb,nb,oner,Evect,nb,rho_temp_nb1,nb,zeror,rho_temp_nb,nb)
            call DGEMM('N','N',nb,nb,nb,oner,rho_temp_nb,nb,Evect,nb,zeror,rho_temp_nb1,nb)
            rnb(:,:,n)=rho_temp_nb1
          end do
c   2nd construct non-local matrix (rho(w))^(-1/2)*rho(K,w)*(rho(w))^(-1/2) 
          do n=-nwn,nwn		  
	  do ick=1,Nc
          do j=1,nb
          do i=1,nb
           rho_temp_nb(i,j)= -(1.0/pi)*(dimag(GdKK_nb(i,j,ick,n)))          
	  end do
	  end do
          call DGEMM('N','N',nb,nb,nb,oner,rnb(:,:,n),nb,rho_temp_nb,nb,zeror,rho_temp_nb1,nb)
          call DGEMM('N','N',nb,nb,nb,oner,rho_temp_nb1,nb,rnb(:,:,n),nb,zeror,rho_temp_nb,nb)
          rho_meas_nb(:,:,ick,n)=rho_meas_nb(:,:,ick,n)+rho_temp_nb(:,:)
          end do
          end do


c   3rd construct log-prefactor for each band    
	     do n=-nwn,nwn        
         do ic=1,Nc
           do j=1,nb
           do i=1,nb
             rho_temp_nb(i,j)= -((1./pi)*dimag(Gcv_nb(i,j,ic,ic,n)))
           end do
           end do
           call DSYEV('V', 'U', nb, rho_temp_nb, nb, Eval, workr, lworkr, ierr)
           Evect=rho_temp_nb
           rho_temp_nb1=zeror
           do i=1,nb
              if(Eval(i).lt.0.d0)then
                 if(myrank==0 .and. dabs(Eval(i)).gt.1.d-8) write(321,*)
     &'rho_loc negative,i=',i,'Eval=',Eval(i),'iter=',iter,'nmeas=',nmeas,'ic=',ic,'n=',n
                 Eval(i)=1.d-8
              end if
             rho_temp_nb1(i,i)=dlog(Eval(i))
           end do
           call DGEMM('T','N',nb,nb,nb,oner,Evect,nb,rho_temp_nb1,nb,zeror,rho_temp_nb,nb)                                                                                        
           call DGEMM('N','N',nb,nb,nb,oner,rho_temp_nb,nb,Evect,nb,zeror,rho_temp_nb1,nb)
           logrholoc_meas_nb(:,:,n) = logrholoc_meas_nb(:,:,n)+rho_temp_nb1
         end do
         end do          
         end if  !mb_tmt_ver8
        
         if (mb_tmt_ver9) then !Matrix version + rho_tot_typ
c     1rst calculate (rho(w))^(-1/2)
          rnb=zeror          
	      do n=-nwn,nwn
          do j=1,nb
          do i=1,nb
             do ic=1,Nc          
             rnb(i,j,n)=rnb(i,j,n)-((1./pi)*dimag(Gcv_nb(i,j,ic,ic,n)))
             end do
             rnb(i,j,n)=rnb(i,j,n)/dfloat(Nc)
          end do
          end do 
          end do
          do n=-nwn,nwn
            call DSYEV('V', 'U', nb, rnb(:,:,n), nb, Eval, workr, lworkr, ierr)
            Evect(:,:)=rnb(:,:,n)
            rho_temp_nb1=zeror
            do i=1,nb
              if(Eval(i).lt.0.d0)then
                 if(myrank==0 .and. dabs(Eval(i)).gt.1.d-8)write(321,*)
     &'rho_ave negative,i=',i,'Eval=',Eval(i),'iter=',iter,'nmeas=',nmeas,'n=',n
                 Eval(i)=1.d-8
              end if
              rho_temp_nb1(i,i)=1.d0/dsqrt(Eval(i))
            end do
            call DGEMM('T','N',nb,nb,nb,oner,Evect,nb,rho_temp_nb1,nb,zeror,rho_temp_nb,nb)
            call DGEMM('N','N',nb,nb,nb,oner,rho_temp_nb,nb,Evect,nb,zeror,rho_temp_nb1,nb)
            rnb(:,:,n)=rho_temp_nb1
          end do
c   2nd construct non-local matrix (rho(w))^(-1/2)*rho(K,w)*(rho(w))^(-1/2) 
          do n=-nwn,nwn		  
	  do ick=1,Nc
          do j=1,nb
          do i=1,nb
           rho_temp_nb(i,j)= -(1.0/pi)*(dimag(GdKK_nb(i,j,ick,n)))          
	  end do
	  end do
          call DGEMM('N','N',nb,nb,nb,oner,rnb(:,:,n),nb,rho_temp_nb,nb,zeror,rho_temp_nb1,nb)
          call DGEMM('N','N',nb,nb,nb,oner,rho_temp_nb1,nb,rnb(:,:,n),nb,zeror,rho_temp_nb,nb)
          rho_meas_nb(:,:,ick,n)=rho_meas_nb(:,:,ick,n)+rho_temp_nb(:,:)
          end do
          end do


c   3rd construct log-prefactor for each band    
	 do n=-nwn,nwn        
         do ic=1,Nc
           rho_temp=zeror
           do i=1,nb
             rho_temp=rho_temp +dabs((1./pi)*dimag(Gcv_nb(i,i,ic,ic,n)))
           end do
           logrholoc_meas_nb(1,1,n) = logrholoc_meas_nb(1,1,n)+dlog(rho_temp)
         end do
         end do          
         end if  !mb_tmt_ver9
       end if  !.not. mb_tmt
       end if  ! liz_embedding
      end if   !mb1
cccccccccccccccccccccccccccccccccccccccccccc     !end mb    
      if(myrank==0) call profend('measure')  
        return
        end
	
