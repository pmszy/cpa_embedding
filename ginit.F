c**************************************************************************** 
        subroutine ginit
c**************************************************************************** 
c       This program makes the initial greens functions, ie. those
c       when all s=+1.  It starts with go (the greens function
c	provided by geodgen corresponding to all s=0), and then uses 
c	the greens function updating block upg to form the desired result.
 
 
c************************************************************
	use Global
	
	implicit none
c************************************************************
        integer i,j,n,ic
        real(kind)  :: x1,x2,y1,y2,randev
        external randev
#define SIMPLE_SPRNG
#define FLOAT_GEN
#define USE_MPI
#include "sprng_f.h"

c	Initialize some variables and accumulators
c       initialize the disorder potentials
	if(iter.eq.1.and.isf.eq.0) then
	
	 
c      if(MB) then
c          do i=1,Nc
c              V1(i)=V1a*sign(1.0,ca-ran(iseed))
c              V2(i)=V2a*sign(1.0,ca-ran(iseed))
c              V12(i)=V12a*sign(1.0,ca-ran(iseed))
c          end do
c      elseif(MB1) then
c          do i=1,nb
c              do j=1,Nc
c                  Vnb(i,j)=V1a*sign(1.0,ca-ran(iseed))
c              end do
c          end do
	 
c	 end if !mb 
	  
	else !if iter1
c	  use the field configuration from the last run or that
c	  was readin from sigma.dat (if isf=1).
      if(myrank==0) write(6,*) 'initialized with old potentials'
	end if



        
        return
        end
