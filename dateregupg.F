c**************************************************************************** 
        subroutine update
c**************************************************************************** 
c       Update updates the potentials by sweeping sequentially through
c       them and proposing changes.  If changes are made
c       it updates the greens functions appropriately.
 
c**************************************************************************** 
	use Global
c	use omp_lib
	implicit none
c**************************************************************************** 
 
        integer jc,check,i,j
        real(kind) :: Vt,start,end1,time,rnum
        real(kind)  :: x1,x2,y1,y2,randev,random_Cauchy
#define SIMPLE_SPRNG
#define FLOAT_GEN
#define USE_MPI
#include "sprng_f.h"
c        external randev
        
c**************************************************************************** 

        if(myrank==0) call profstart('dateregupg') 



c  New way
          if (my) then
          if(MB1) then !if MB1
           Vnb=zeror
           if(tables_type>=4)then
           if (binary) then
	         do j=1,Nc
                 do i=1,nb
                     !Vnb(i,j)=V0*dsign(1.0,ca-ran(iseed))
                     Vnb(i,j)=V0*dsign(1.0,ca-sprng())
                 end do
             end do
             if(exact)then !evaluate all disorder configurations
                 do j=1,Nc
                     do i=1,nb
                         if(BTEST(meas*myrank+nmeas-1,(nb*(j-1)+i-1)))then
                             Vnb(i,j)=V0
                         else
                             Vnb(i,j)=-V0
                         end if
                     end do
                 end do
             end if
           end if  !binary

	       if (box) then
	         do j=1,Nc
                do i=1,nb
                    !Vnb(i,j)=2.0*V0*(ran(iseed)-0.5)
                    Vnb(i,j)=2.0*V0*(sprng()-0.5)
                end do
	         end do
	       end if  !box
           else
           if (binary) then
	         do j=1,Nc
                     rnum=sprng()
                 do i=1,nb
                     !Vnb(i,j)=V1a*dsign(1.0,ca-ran(iseed))
                     Vnb(i,j)=V1a*dsign(1.0,ca-rnum)
                     !Vnb(i,j)=V1a*dsign(1.0,ca-sprng())
                 end do
             end do
           end if  !binary 

	       if (box) then
                 !if(myrank==0) write(*,*) 'BOXxxxxxxxxxxxx-dat'
	         if(iter==1)then
                   do j=1,Nc
                   do i=1,nb
                      !Vnb(i,j)=2.0*V1a*(ran(iseed)-0.5)
                      rnum=sprng() 
                      Vnb(i,j)=2.0*V1a*(rnum-0.5)
                      Vnbc(i,j,nmeas)=Vnb(i,j)
                   end do
	           end do
                 else
                   Vnb(:,:)=Vnbc(:,:,nmeas)  
                 end if
	       end if  !box    
           end if  !tables_type
	       end if  !MB1
	      end if !if my
cc  -----end original way







c!$OMP END PARALLEL

c	The Gd matrix is reconfigure for each of the runs
c	Reconfiguration is necessary to deal with roundoff error.
 
            call reg
            if(myrank==0) call profend('dateregupg')
        return
        end
!**************************************************************************** 
        subroutine reg
!       This subroutine reconfigures G given only the potential configuration
!       and the initial green's functions Gcsfsc (which corresponds to
!	all V=0).
 
c**************************************************************************** 
	use Global
c	use omp_lib
	implicit none
c**************************************************************************** 
 
        integer n,ic,jc,info,l,m,i1,j1,i,j
        real(kind) temp
c       Gcsfsc is the greens function for all V=0, hence we will
c       sweep through the lattice and reconfigure Gcsfsc to coorespond
c       to the lattice configuration.
c**************************************************************************** 

        if(MB1) then
          V_tpmb_nb=zeroc
         do i=1,nb
            do j=1,Nc
                V_tpmb_nb(i,i,j,j)=dcmplx(Vnb(i,j),0.d0)
            end do
         end do
        end if ! MB1	
	
	

c*************************************************************************************************************************
ccc	In this section, we will use lapack subroutines to invert a matrix, we will do this via LU decomposition 
ccc	using zgetri which should be called after calling the factorization subroutine zgetrf. Note, these subroutines 
ccc	are for double precision for single precision, change kind in module to 4 and the subroutines as cgetri and cgetrf, 
ccc	respectively. C. E. Ekuma: 09-02-12.
c*************************************************************************************************************************

        if (MB1) then
ccccccccccccccccccccccccccccccccccccccc
c         Multiband
ccccccccccccccccccccccccccccccccccccccc
c  construct for each wn Gtemp=w(nb,nb)-V(nb,nb,i,j)-\Delta(nb,nb,i,j)

        del_tpmb_nb(:,:,:,:,:)=zeroc
        G_tpmb_nb(:,:,:,:)=zeroc
        G_inmb_nb(:,:,:)=zeroc
        w_tpmb_nb(:,:,:,:)=zeror
        if(Gscript.and.init_sigma)then
         if(Liz_embedding)then
          Gc_liz=zeroc
          temp=zeror 
          do ic=1,Nc   
          do n=-nwn,nwn
          do i=1,lizsize
            do j=1,lizsize
               G_inmb_nb(i,j,n)=Gsc_liz(i,j,n)-V_tpmb_nb(1,1,lizmap(ic,i),lizmap(ic,j))
            end do
          end do
         
          call zgetrf(lizsize,lizsize,G_inmb_nb(1:lizsize,1:lizsize,n),lizsize,ipvtliz,info)
          if (info.eq.0)then
             call zgetri(lizsize,G_inmb_nb(1:lizsize,1:lizsize,n),lizsize,ipvtliz,work,lwork,info)
          else
             print*, 'info=/0 Matrix Inversion Failed in update: Info=',info
             write(*,*)'n=',n
             stop
          end if
          Gc_liz(:,:,n)=Gc_liz(:,:,n)+G_inmb_nb(1:lizsize,1:lizsize,n)/real(Nc)
          if(mb_tmt_ver1)then
             logrholoc_meas_nb(1,1,n)=logrholoc_meas_nb(1,1,n)+dlog(dabs(-1.d0/pi*dimag(G_inmb_nb(1,1,n))))/real(Nc)
             !logrholoc_meas_nb(1,1,n)=logrholoc_meas_nb(1,1,n)+(-1.d0/pi*dimag(G_inmb_nb(1,1,n)))/real(Nc)
          elseif(mb_tmt_ver2)then
             temp=zeror
             do i=1,lizsize
             logrholoc_meas_nb(1,1,n)=logrholoc_meas_nb(1,1,n)+dlog(dabs(-1.d0/pi*dimag(G_inmb_nb(i,i,n))))/real(Nc*lizsize)
             temp=temp+(-1.d0/pi*dimag(G_inmb_nb(i,i,n)))/real(lizsize)
             end do
             Gc_liz_meas(:,:,n)=Gc_liz_meas(:,:,n)+G_inmb_nb(1:lizsize,1:lizsize,n)/temp/real(Nc)
          end if
          if(iter==niter)then
            if(n==0) temp=temp+(dimag(G_inmb_nb(1,1,0)))**2
          end if             
          end do !n
         end do !ic
         if(iter==niter) GIPR2=temp/(Nc*dimag(Gc_liz(1,1,0)))**2 

         else !not liz_embedding 
          do n=-nwn,nwn
          do ic=1,Nc
            do jc=1,Nc
              do l=1,nb
                do m=1,nb
                  i1=(ic-1)*nb+l
                  j1=(jc-1)*nb+m
                  G_inmb_nb(i1,j1,n)=GcRfsc_nb(l,m,ic,jc,n)
                end do
              end do
            end do
          end do
          call zgetrf(nb*Nc,nb*Nc,G_inmb_nb(1,1,n),nb*Nc,ipvtn,info)
	     if (info.eq.0)then
	     call zgetri(nb*Nc,G_inmb_nb(1,1,n),nb*Nc,ipvtn,work,lwork,info)
             else
              print*, 'info=/0 Matrix Inversion Failed in update: Info=',info
	      stop
             end if
          end do 
        
         do n=-nwn,nwn
           do ic=1,Nc
             do jc=1,Nc
                 do l=1,nb
                     do m=1,nb
               i1=(ic-1)*nb+l
               j1=(jc-1)*nb+m
                 G_inmb_nb(i1,j1,n)=G_inmb_nb(i1,j1,n)-V_tpmb_nb(l,m,ic,jc)
                     end do !m
                 end do !l 
             end do !jc
           end do  !ic
         end do
         end if!if liz_embedding
        else ! if not init_sigma 
        do n=-nwn,nwn
            do i=1,nb
                do j=1,nb
                    do ic=1,Nc
                        do jc=1,Nc
            del_tpmb_nb(i,j,ic,jc,n)=Gamma_rs_nb(i,j,ic,jc,n)+Epsb_rs_nb(i,j,ic,jc)
            if ((ic .eq. jc) .and. (i .eq. j)) then
             w_tpmb_nb(i,j,ic,jc)=wn(n)
            end if
                        end do !jc
                    end do !ic
                end do !j
            end do !i
          
         do ic=1,Nc
             do jc=1,Nc
                 do l=1,nb
                     do m=1,nb
	     G_tpmb_nb(l,m,ic,jc)=w_tpmb_nb(l,m,ic,jc)-V_tpmb_nb(l,m,ic,jc)-del_tpmb_nb(l,m,ic,jc,n)
               i1=(ic-1)*nb+l
               j1=(jc-1)*nb+m
                 G_inmb_nb(i1,j1,n)=G_tpmb_nb(l,m,ic,jc)
                 if(ic.eq.jc .and. l.eq.m) G_inmb_nb(i1,j1,n)=G_inmb_nb(i1,j1,n)+eta
                     end do !m
                 end do !l 
             end do !jc
         end do  !ic
          
        end do !end wn loop
        end if !Gscipt
         end if ! if mb1
            


cccccccccccccccccccccccccccccccccccccccccc
c   end MB
cccccccccccccccccccccccccccccccccccccccccc   


      if(MB1) then !MB1 case
      if(.not. Liz_embedding)then     
	  do n=-nwn,nwn
	   call zgetrf(nb*Nc,nb*Nc,G_inmb_nb(1,1,n),nb*Nc,ipvtn,info)
	     if (info.eq.0)then
	     call zgetri(nb*Nc,G_inmb_nb(1,1,n),nb*Nc,ipvtn,work,lwork,info)
      else
           print*, 'info=/0 Matrix Inversion Failed in update: Info=',info
	  stop
         end if
c       move G_temp to Gfor each band

         do i=1,nb
             do j=1,nb
                 do i1=i,nb*Nc,nb
                    do j1=j,nb*Nc,nb
                        ic=(i1-i)/nb+1
                        jc=(j1-j)/nb+1
                        Gcv_nb(i,j,ic,jc,n)=G_inmb_nb(i1,j1,n)
                    end do            
                 end do
             end do
         end do 

      end do !wn loop
      end if
	  end if! if MB1
      
        

	

        return
        end
	

