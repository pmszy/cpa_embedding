cc**************************************************************************** 
        program fmain
	use Global
c	use omp_lib
	implicit none
 
c**************************************************************************** 
 
c       This program simulates the 3-dimensional disorded Anderson-hubbard 
c	model in the Cluster Typical Medium Theory Approximation
c
c                        H=H1+H2
c                             +
c                H1= -t SUM ((C (i,s)C(j,s)) + h.c. )
c                       <ij>sproperty of
c                               +  
c                    -t' SUM ((C (i,s)C(j,s)) + h.c. )
c                       [ij]s                           i,j n.n.n.
c
c                H2=   SUM (ed + V ) * n(i,s)
c                       i         i    
c
c	V is a random binary variable with <V>=0.  H is particle-hole 
c	symmetric whenever ed=0 and t'=0.  For a two-dimensional system, 
c	this means that 
c
c               +             +    i
c       C  --> C exp(iQ.i) = C (-1)
c        i      i             i
c                                                +              +
c       C  = SUM exp(ik.i) C  --> SUM exp(ik.i) C exp(-iQ.i) = C
c        k    i             i                    i              Q-k
c
c       or if ed=0, G (tau) = -G   (-tau) = G   (beta-tau)
c                    k          Q-k          Q-k
c       and
c                   G (iwn) = - G   (-iwn)  real parts opposite sign
c                    k           Q-k        imaginary parts equal sign
c
c	in real-frequencies, this means that 
c
c                    R          A            *
c                   G (w) = - G   (-w)  =  -G   (-w)
c                    k         Q-k           Q-k
c

c**************************************************************************** 
c	Original serial code by Mark Jarrell 1999-2000
c
c	The modification to implement Typical Medium Theory within a cluster (CTMT) is
c	by Chinedu Ekuma, Zi Yang Meng, and Hanna Terletska 2012
c
c       Checks:
c       ( ) cat *.f > onefile ; mv onefile onefile.f
c           ftnchek -usage=303 -portability=all onefile.f
c
c
c**************************************************************************** 

        integer n,m,i,j,ic,jc,info
        real(kind) ::  time1,time2,r1,r2,r3,e,time_main1,time_main2
        character*72 linemc
cH:
        complex*16 :: A1,B1,C1,D1
      
c**************************************************************************** 
!  Initialize MPI.
	    call MPI_Init(ierr)
	    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)   ! get the # processors
	    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)  ! get processor's ID
       time_main1= secnds(0.0)
c       Readin the parameters
        call readin
        
        call allocate_arrays
        if(myrank==0) write(*,*) 'allocate arrays-2'
	
c       we must now make the lookup tables
        call tables
        if(myrank == 0) then
 	    write(6,*) '***************************************************************'					      
	    write(6,*) "table completed successfully."
        end if
        isf=0
        test=0.d0	!convergence criteria
        testl=0.d0
	    call init_parameters
        do iter=1,niter					
c	  if(isi.ge.1.and.iter.gt.1) meas=int(meas*ifac)		
c	Added to enable the sampling of Large meas for histogram
        if(histo_enable) then
c	  if(isi.ge.1.and.iter.gt.1.and.iter.ne.niter) meas=int(meas*ifac)		
        if(isi.ge.1.and.iter.eq.niter) meas=int(meas*ifac)	
        end if ! End histogram 	
        if(myrank == 0) then
	     write(6,*) 'iteration',iter,' runs=',run,' meas=',meas
	     write(42,*) 'iteration',iter,' runs=',run,' meas=',meas				
        end if
        if(iter.eq.1) then
c	    Write out a data header for the binned data

	    call output(0)
  				
        end if

c         Make Gcsfsc, the greens function for all V(l)=0 and initialize all the accumulators
          call geodgen

c         Initialize the greens functions and the fields.
          call ginit

 
c         The following is the main part of the main program.  It directs
c         measurement, and updating the lattice and the d greens functions
c         between the measurements.
          do nrun=1,run
            do nmeas=1,meas
              call update
              call measure
            end do 

c           Now load the data into bins and accumulators.
            call sumup
          end do
          call output(1)
      if(myrank == 0) then
      if(iter.eq.niter) then
 	    write(6,*) "Writing Data to Output"
	  end if
        time_main2= secnds(time_main1)
 	  write(6,*) 'done iteration',iter,'at time',time_main2
 	  write(6,*) '***********************************************************'		      
 	  write(42,*) 'done iteration',iter,'at time',time_main2
 	  write(42,*) '**********************************************************'		      

      if(iter.eq.niter-1) then
 	    write(6,*) "Last Iteration. Major output will be written into dos.dat"
 	    write(6,*) "The TotalDOS and PDOS will be written into p_rho.dat"
 	    write(6,*) "To understand other outputs, see README, output section"
 	    write(6,*) '***********************************************************'		      
 	    write(42,*) "Last Iteration. Major output will be written into dos.dat"
 	    write(42,*) "The TotalDOS and PDOS will be written into p_rho.dat"
 	    write(42,*) "To understand the many outputs, see README output section"
 	    write(42,*) '**********************************************************'		      
	  end if
      end if
      end do! done iteration loop
        if(myrank==0) call profstat()
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)   ! clean up
        call deallocate_arrays
        stop

        end
