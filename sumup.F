c***************************************************************************************
	subroutine sumup
c***************************************************************************************
c
c***************************************************************************************
	use Global
c	use omp_lib
	implicit none
c***************************************************************************************
	integer i,j,ic,n,m
        real(kind) :: data_nb(nb,nb,Nc,-nwn:nwn)
c***************************************************************************************
        if(myrank==0) call profstart('sumup')
c       Now normalize rho_meas and rho_Gdkk 

        iF (MB1) THEN
c  I:   Normalize rho_meas
        if(Liz_embedding)then
            Gc_liz_run=Gc_liz_meas/real(meas)
            Gc_liz_f=Gc_liz_f+Gc_liz_run
            if(mb_tmt)then
              logrholoc_run_nb=logrholoc_meas_nb/real(meas)
              logrholoc_f_nb=logrholoc_f_nb+logrholoc_run_nb
            end if
            if(iter==niter)then
              GIPR2_run=GIPR2_meas/real(meas)
              logGIPR2_run=logGIPR2_meas/real(meas)
              GIPR2_f=GIPR2_f+GIPR2_run   
              logGIPR2_f=logGIPR2_f+logGIPR2_run   
            end if    
        !if(myrank==0)then
        !  write(*,*)-1.d0/pi*dimag(Gc_liz_f(1,1,0)),logrholoc_f_nb(1,1,0)
        !end if
        else!not liz_embedding
        rhorun_nb = rho_meas_nb/real(meas)
        if(Gfav) Gcrun_nb =Gc_meas_nb/real(meas)
        
        
c. II: impose point-group symmetry

         do i=1,nb
         do j=1,nb    
         do n=-nwn,nwn
	     do ic=1,Nc
	     data_nb(i,j,ic,n) = zeror
	     do igroup=1,ngroup          
              data_nb(i,j,ic,n) = data_nb(i,j,ic,n) +				      				      
     &  	       rhorun_nb(i,j,ickequ(ic,igroup),n)/real(ngroup)
 	     end do
         end do
         end do        
	     end do
         end do        

c	     rhorun_nb=data_nb


                  
c  III: Load the measurements for this run into those for the iteration
         do i=1,nb
         do j=1,nb    
         do ic=1,Nc
	     do n =-nwn,nwn
	       rhof_nb(i,j,ic,n) = rhof_nb(i,j,ic,n) + rhorun_nb(i,j,ic,n)      
         end do
         end do
         end do
         end do
         
         if(Gfav) Gcf_nb=Gcf_nb+Gcrun_nb
         
          if (mb_tmt) then
            logrholoc_run_nb = logrholoc_meas_nb/real(meas*Nc)

            if(iter.eq.niter)then
                rholoc_av_run_nb =rholoc_av_meas_nb/real(meas*Nc)
                logrhotot_run =logrhotot_meas/real(meas*Nc)
            endif
 
           logrholoc_f_nb=logrholoc_f_nb + logrholoc_run_nb
           
           if(iter.eq.niter)then
               rholocf_av_nb=rholocf_av_nb+rholoc_av_run_nb
               logrhotot_f=logrhotot_f+logrhotot_run
           endif
 
           end if
        end if !liz_embedding
        END IF !if MB1  
c***************************************************************************************
c	Now zero the accumulator for a single run
c***************************************************************************************
        


                
        
        rho_meas_nb=zeror
        logrholoc_meas_nb=zeror
        rhorun_nb=zeror
        logrholoc_run_nb=zeror
        rholoc_av_meas_nb=zeror
        logrhotot_meas=zeror

        Gc_meas_nb=zeroc
        Gcrun_nb=zeroc

        Gc_liz_meas=zeroc
        Gc_liz_run=zeroc
          
        GIPR2_meas=zeror
        logGIPR2_meas=zeror
        GIPR2_run=zeror
        logGIPR2_run=zeror
          
        
        
        if(myrank==0) call profend('sumup')         



	return
	end

