!!!!
!! File: diagonalize_zheev.f90
!! Description: Diagonalization of a HERMITEAN matrix using LAPACK's ZHEEV routine
!!              http://www.netlib.org/lapack/explore-html/df/d9a/               group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Thursday, 5th August 2021, 9:07:40 am
!! Last Modified: Thursday, 5th August 2021, 9:24:12 am
!!  
!! Copyright (c) 2021, Bruno R. de Abreu, National Center for Supercomputing Applications.
!! All rights reserved.
!! License: This program and the accompanying materials are made available to any individual
!!          under the citation condition that follows: On the event that the software is
!!          used to generate data that is used implicitly or explicitly for research
!!          purposes, proper acknowledgment must be provided in the citations section of
!!          publications. This includes both the author's name and the National Center
!!          for Supercomputing Applications. If you are uncertain about how to do
!!          so, please check this page: https://github.com/babreu-ncsa/cite-me.
!!          This software cannot be used for commercial purposes in any way whatsoever.
!!          Omitting this license when redistributing the code is strongly disencouraged.
!!          The software is provided without warranty of any kind. In no event shall the
!!          author or copyright holders be liable for any kind of claim in connection to
!!          the software and its usage.
!!!!

program diagonalize
    implicit none
    integer :: order
    parameter (order = 3)
    complex*16, allocatable :: matrix(:,:), eigvecs(:,:)
    real*8, allocatable :: eigvals(:)

    integer :: i,j

    allocate(matrix(order,order))
    allocate(eigvecs(order,order))
    allocate(eigvals(order))

    do j=1,order
        do i=1,order
            matrix(i,j) = (0.0,0.0)
        enddo
    enddo

    !! Sy for S = 1
    matrix(1,1) = (0.00, 0.00)
    matrix(1,2) = (1.00, 0.00)
    matrix(1,3) = (0.00, 0.00)
    matrix(2,1) = (-1.00, 0.00)
    matrix(2,2) = (0.00, 0.00)
    matrix(2,3) = (1.00, 0.00)
    matrix(3,1) = (0.00, 0.00)
    matrix(3,2) = (-1.00, 0.00)
    matrix(3,3) = (0.00, 0.00)
    
    call hermitean_diagonalization(matrix,order,eigvals,eigvecs)
    write(*,*) 'Eigenvalues:'
    write(*,*) eigvals
    write(*,*) 'Eigenvectors (columns):'
    write(*,*) eigvecs

end program

subroutine hermitean_diagonalization(matrix,order,eigvals,eigvecs)
    IMPLICIT NONE
    INTEGER order
    INTEGER LWMAX
    PARAMETER (LWMAX=1000)
    INTEGER INFO, LWORK
    REAL*8 W(order), RWORK(3*order-2), eigvals(order)
    COMPLEX*16 A(order,order), WORK(LWMAX)
    COMPLEX*16 matrix(order,order), eigvecs(order,order)

    A = matrix
    LWORK=-1
    call ZHEEV('V', 'U', order, A, order, W, WORK, LWORK, RWORK, INFO)
    LWORK = min(LWMAX, int(WORK(1)))
    call ZHEEV('V', 'U', order, A, order, W, WORK, LWORK, RWORK, INFO)
    if(INFO == 0) then
        eigvals = W
        eigvecs = A
    else
        write(*,*) 'Diagonalization failed.'
        return
    endif

end subroutine hermitean_diagonalization



                                            

                                            