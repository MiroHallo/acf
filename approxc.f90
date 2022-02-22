! Fortran f90 subroutines for determining the (cross-)covariance matrix by
! approximate (cross-)covariance functions (AXCF or SAXCF) (Hallo and Gallovic, 2016).
! Hallo,M., Gallovic,F. (2016): Fast and cheap approximation of Green functions uncertainty
! for waveform-based earthquake source inversions, Geophysical Journal Int., 207, 1012-1029.
!
! Authors: Miroslav Hallo and Frantisek Gallovic (8/2017)
! Charles University, Faculty of Mathematics and Physics
! Revision 1/2018: Corrected and tested for real signals
!
! Copyright (C) 2017,2018  Miroslav Hallo and František Gallovič
!
! This program is published under the GNU General Public License (GNU GPL).
!
! This program is free software: you can modify it and/or redistribute it
! or any derivative version under the terms of the GNU General Public
! License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! This code is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
! and don't remove their names from the code.
!
! You should have received copy of the GNU General Public License along
! with this program. If not, see <http://www.gnu.org/licenses/>.
!
!---------------------------------------------------------------------



SUBROUTINE axcf(N,f1,f2,L1,L12,C)
!---------------------------------------------------------------------
!    Compute AXCF of functions f1 and f2 of samples N
!    L1 and L12 are time-shift windows in samples
!    -> Output into covariance matrix C
!---------------------------------------------------------------------
    implicit none
    integer,intent(in):: N,L1,L12
	real(8),intent(in):: f1(N),f2(N)
	real(8),intent(out):: C(N,N)
    complex(8),allocatable:: smooth1(:),smooth12(:),s1(:),s2(:),s1smooth1(:),s2smooth1(:),s2smooth12(:),XC(:,:)
    integer:: i,j,FFT,Lzeros
	
! Check the length of smoothing window L1
    if(L1<2)then
	  write(*,*) ' Warning! L1 is too short'
	endif
	
! Allocation and init of variables
    FFT = 2**(int(log(dble(2*N))/log(2.d0)+0.99999999d0)+1)
    allocate(smooth1(FFT),s1(FFT),s2(FFT),s1smooth1(FFT),s2smooth1(FFT),s2smooth12(FFT),XC(FFT,FFT))
    C=0.d0
    s1=0.
    s2=0.
	s1(1:N) = f1(1:N)
	s2(1:N) = f2(1:N)
	
! Put zeros (values) before and after the signal
	Lzeros = max(L1,L12);
	s1(N+1:N+Lzeros) = f1(N)
	s2(N+1:N+Lzeros) = f2(N)
	s1(FFT-Lzeros+1:FFT) = f1(1)
	s2(FFT-Lzeros+1:FFT) = f2(1)
	
! f2 smoothing by L12 window
    if(L12>1)then
      allocate(smooth12(FFT))
      smooth12=0.
      smooth12(1:int(dble(L12)/2.+.51))=1.d0/dble(L12*FFT)
      smooth12(FFT:FFT-L12/2+1:-1)=1.d0/dble(L12*FFT)
      call four1(smooth12,FFT,1)
      call four1(s2,FFT,1)
      s2smooth12=s2*smooth12
      call four1(s2smooth12,FFT,-1)
      deallocate(smooth12)
    else
      s2smooth12=s2
    endif
	
! smoothing window L1
    smooth1=0.
    smooth1(1:int(dble(L1)/2.+.51))=1.d0/dble(L1*FFT)
    smooth1(FFT:FFT-L1/2+1:-1)=1.d0/dble(L1*FFT)
    call four1(smooth1,FFT,1)
	
! compute XCF
    do i=1,FFT
      XC(:,i)=s1(:)*cshift(s2smooth12(:),SHIFT=(FFT/2+1)-i)
      call four1(XC(:,i),FFT,1)
      XC(:,i)=XC(:,i)*smooth1(:)
      call four1(XC(:,i),FFT,-1)
    enddo

    call four1(s1,FFT,1)
    s1smooth1=s1*smooth1
    call four1(s1smooth1,FFT,-1)
    call four1(s2smooth12,FFT,1)
    s2smooth1=s2smooth12*smooth1
    call four1(s2smooth1,FFT,-1)
	
    do i=1,FFT
      XC(:,i) = XC(:,i) - s1smooth1(:)*cshift(s2smooth1(:),SHIFT=(FFT/2+1)-i)
    enddo
	
! Fill the covariance matrix by AXCF
    do i=1,N
     do j=1,N
       C(i,j)=dble(XC(i,FFT/2+1-(j-i)))
     enddo
    enddo

    deallocate(smooth1,s1,s2,s2smooth12,s1smooth1,s2smooth1,XC)
	
!---------------------------------------------------------------------
END SUBROUTINE



SUBROUTINE saxcf(N,f1,f2,L1,L12,T,C)
!---------------------------------------------------------------------
!    Compute SAXCF of functions f1 and f2 of samples N
!    L1 and L12 are time-shift windows in samples
!    T is characteristic length of the signal in samples
!    -> Output into covariance matrix C
!---------------------------------------------------------------------
    implicit none
	real(8),parameter:: PI=3.1415926535d0
    integer,intent(in):: N,L1,L12,T
    real(8),intent(in):: f1(N),f2(N)
	real(8),intent(out):: C(N,N)
    complex(8),allocatable:: s1(:),s2(:),Rfg(:),fSACF(:),smooth1(:),smooth12(:)
	real(8),allocatable:: SACF(:),tw(:)
	real(8):: taper,r
    integer:: i,j,FFT

! Check the length of smoothing window L1
    if(L1<2)then
	  write(*,*) ' Warning! L1 is too short'
	endif
	
! Allocation and init of variables
    FFT = 2**(int(log(dble(2*N))/log(2.d0)+0.99999999d0)+1)
    allocate(s1(FFT),s2(FFT),Rfg(FFT),fSACF(FFT),smooth1(FFT))
	allocate(SACF(2*N-1),tw(2*N-1))
    C=0.d0
    s1=0.
    s2=0.
	s1(1:N) = f1(1:N)
	s2(1:N) = f2(1:N)
	
! FFT of the f and g functions
    call four1(s1,FFT,1)
    call four1(s2,FFT,1)
	
! Prepare smoothing functions in freq. domain
    smooth1=0.
    smooth1(1:int(dble(L1)/2.+.51))=1.d0/dble(L1)
    smooth1(FFT:FFT-L1/2+1:-1)=1.d0/dble(L1)
    call four1(smooth1,FFT,1)
	
!  Cross-correlation of signals in freq. domain
	Rfg = conjg(s1)*s2
	
!  Convolution of triangle function (width 2*L1) and cross-correlation in freq. domain
     fSACF = Rfg - smooth1*conjg(smooth1)*Rfg

!  Smoothing by joint shift (width L12) in freq. domain
    if(L12>1)then
      allocate(smooth12(FFT))
      smooth12=0.
      smooth12(1:int(dble(L12)/2.+.51))=1.d0/dble(L12)
      smooth12(FFT:FFT-L12/2+1:-1)=1.d0/dble(L12)
      call four1(smooth12,FFT,1)
	  fSACF = smooth12*fSACF
      deallocate(smooth12)
    endif
	
!  Norm by the effective signal length
	fSACF = fSACF/dble(T)
	
!  Back to time domain
	call four1(fSACF,FFT,-1)
	fSACF = fSACF/FFT
	SACF(1:N-1) = dble(fSACF(FFT-N+2:FFT))
	SACF(N:2*N-1) = dble(fSACF(1:N))
	
!  Tapered cosine window (SACF)
	taper = 0.666d0
	tw=1.d0
	do i=1,2*N-1
	  r = dble(i-1)/dble(2*N-2)
	  if(r.le.(taper/2.d0))then
	    tw(i) = 0.5d0*( 1 + cos((2.d0*PI/taper)*(r-taper/2.d0)) )
	  elseif(r.ge.(1.d0-taper/2.d0))then
	    tw(i) = 0.5d0*( 1 + cos((2.d0*PI/taper)*(r-1+taper/2.d0)) )
	  endif
	enddo
	tw = tw**3 ! cube for bigger effect
	SACF = SACF * tw !  Taper cosine window
	
! Fill the covariance matrix by SACF
    do i=0,N-1
      C(i+1,1:N) = SACF(N-i:2*N-i-1)
    enddo
	
    deallocate(s1,s2,Rfg,fSACF,SACF,tw,smooth1)

!---------------------------------------------------------------------
END SUBROUTINE



!---------------------------------------------------------------------
! Fortran 77 subroutine four1 from Numerical Recipes used by saxcf and axcf
! Include it in separate file nr.for

      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
            tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
	  
	  
	  