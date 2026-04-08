!---------------------------------------------------------------------
! Fortran 90 example program that uses subroutines for determining
! covariance matrices by ACF, AXCF, SACF, SAXCF
!---------------------------------------------------------------------
!
! Author: Miroslav HALLO
! Charles University in Prague, Faculty of Mathematics and Physics
! E-mail: hallo@karel.troja.mff.cuni.cz
! Revision 2026/04: Initial version
! Method:
! Hallo, M., Gallovic, F. (2016): Fast and cheap approximation of Green
!         functions uncertainty for waveform-based earthquake source 
!         inversions, Geophys. J. Int., 207 1012-1029.
!         https://doi.org/10.1093/gji/ggw320
!
! Copyright (C) 2026 Miroslav Hallo
!
! This program is published under the GNU General Public License (GPL).
!
! This program is free software: you can modify it and/or redistribute
! it or any derivative version under the terms of the GNU General Public
! License as published by the Free Software Foundation, either version 3
! of the License, or (at your option) any later version.
!
! This code is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY. We would like to kindly ask you to acknowledge
! the authorsand don't remove their names from the code.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!
!---------------------------------------------------------------------

PROGRAM run_example
!---------------------------------------------------------------------
!  Main program
!---------------------------------------------------------------------
    use ApproxCovMod
    implicit none
    
    integer:: nsampl, ios, L1, L12, T, i, j, fid
    real(8):: dt, L1_sec, L12_sec, T_sec
    real(8),allocatable:: f(:), g(:), C(:,:,:)
    character(256):: line
    character(64):: infile
    
    !----------------------------------
    ! Code INPUT
	
	! Filename of the text file with input f(t) and g(t) functions
    infile = 'example_data.txt'
	! Sample interval dt [s]
    dt = 0.1d0
	! Width of joint (f and g) uniform distribution of time-shifts [s]
    L1_sec = 2.0d0
	! Width of relative (f vs g) uniform distribution of time-shifts [s]
    L12_sec = 1.0d0
	! Dominant signal length  [s]
    T_sec = 10.0d0
    
    !----------------------------------
    ! Convert seconds to samples
    L1 = ceiling(L1_sec/dt)
    L12 = ceiling(L12_sec/dt)
    T = ceiling(T_sec/dt)
    
    !----------------------------------
    ! Count number of lines of the text file with input data
    open(10,file=infile,status='old',action='read')
    nsampl = 0
    do
      read(10,'(A)',iostat=ios) line
      if(ios /= 0)then
        exit
      endif
      ! Remove whitespace
      line = adjustl(line)
      ! Skip comment or empty line
      if(line(1:1)=='#' .or. line=='')then
        cycle
      endif
      nsampl = nsampl + 1
    enddo
    
    !----------------------------------
    ! Allocate memory
    allocate(f(nsampl), g(nsampl), C(nsampl,nsampl,4))
    
    !----------------------------------
    ! Read input data
    rewind(10)
    i = 1
    do
      read(10,'(A)',iostat=ios) line
      if(ios /= 0)then
        exit
      endif
      ! Remove whitespace
      line = adjustl(line)
      ! Skip comment or empty line
      if(line(1:1)=='#' .or. line=='')then
        cycle
      endif
      read(line, *) f(i), g(i)
      i = i + 1
    enddo
    close(10)
    
    !----------------------------------
    ! Compute ACF
    call axcf(nsampl,f,f,L1,0,C(:,:,1))
    
    !----------------------------------
    ! Compute AXCF
    call axcf(nsampl,f,g,L1,L12,C(:,:,2))
    
    !----------------------------------
    ! Compute SACF
    call saxcf(nsampl,f,f,L1,0,T,C(:,:,3))
    
    !----------------------------------
    ! Compute SAXCF
    call saxcf(nsampl,f,g,L1,L12,T,C(:,:,4))
    
    !----------------------------------
    ! Save ACF to file
    fid = 101
    open(newunit=fid,form='formatted',file='example_ACF.txt',status='replace')
    do i=1,nsampl
      do j=1,nsampl
        write(fid,'(e14.6)',advance='no') C(i,j,1)
      enddo
      write(fid, *)
    enddo
    close(fid)
	
	!----------------------------------
    ! Save AXCF to file
    fid = 102
    open(newunit=fid,form='formatted',file='example_AXCF.txt',status='replace')
    do i=1,nsampl
      do j=1,nsampl
        write(fid,'(e14.6)',advance='no') C(i,j,2)
      enddo
      write(fid, *)
    enddo
    close(fid)
	
	!----------------------------------
    ! Save ACF to file
    fid = 103
    open(newunit=fid,form='formatted',file='example_SACF.txt',status='replace')
    do i=1,nsampl
      do j=1,nsampl
        write(fid,'(e14.6)',advance='no') C(i,j,3)
      enddo
      write(fid, *)
    enddo
    close(fid)
	
	!----------------------------------
    ! Save ACF to file
    fid = 104
    open(newunit=fid,form='formatted',file='example_SAXCF.txt',status='replace')
    do i=1,nsampl
      do j=1,nsampl
        write(fid,'(e14.6)',advance='no') C(i,j,4)
      enddo
      write(fid, *)
    enddo
    close(fid)
    
    !----------------------------------
    ! Deallocate memory
    deallocate(f, g, C)
    
!---------------------------------------------------------------------
END PROGRAM
