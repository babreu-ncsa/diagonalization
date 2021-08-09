!!!!
!! File: diagonalize_zheev.f90
!! Description: Diagonalization of a HERMITEAN matrix using LAPACK's ZHEEV routine
!!              http://www.netlib.org/lapack/explore-html/df/d9a/group__complex16_h_eeigen_gaf23fb5b3ae38072ef4890ba43d5cfea2.html
!! Author: Bruno R. de Abreu  |  babreu at illinois dot edu
!! National Center for Supercomputing Applications (NCSA)
!!  
!! Creation Date: Thursday, 5th August 2021, 9:07:40 am
!! Last Modified: Monday, 9th August 2021, 12:39:12 pm
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
    parameter (order = 100)
    complex*16, allocatable :: matrix(:,:), eigvecs(:,:)
    real*8, allocatable :: eigvals(:)
    integer :: i,j

    allocate(matrix(order,order))
    allocate(eigvecs(order,order))
    allocate(eigvals(order))

    matrix = (1.0, 0.0)

    call hermitean_diagonalization(matrix,order,eigvals,eigvecs)
    write(*,*) 'Eigenvalues:'
    write(*,*) eigvals


    deallocate(matrix)
    deallocate(eigvals)
    deallocate(eigvecs)
end program diagonalize

subroutine hermitean_diagonalization(matrix,order,eigvals,eigvecs)
    IMPLICIT NONE
    COMPLEX*16, intent(in) :: matrix(order,order)
    INTEGER, intent(in) :: order
    REAL*8, intent(out) :: eigvals(order)
    COMPLEX*16, intent(out) :: eigvecs(order,order)
    INTEGER LWMAX
    PARAMETER (LWMAX=1000)
    INTEGER INFO, LWORK
    REAL*8 RWORK(3*order-2)
    COMPLEX*16, allocatable :: WORK(:)

    allocate(WORK(LWMAX))
    eigvecs = matrix
    LWORK=-1
    call ZHEEV('V', 'U', order, eigvecs, order, eigvals, WORK, LWORK, RWORK, INFO)
    LWORK = min(LWMAX, int(WORK(1)))
    deallocate(WORK)
    allocate(WORK(LWORK))
    call ZHEEV('V', 'U', order, eigvecs, order, eigvals, WORK, LWORK, RWORK, INFO)
    if(INFO == 0) then
        write(*,*) 'Diagonalization performed.'
    else
        write(*,*) 'Diagonalization failed.'
        return
    endif

end subroutine hermitean_diagonalization



                                            

                                            