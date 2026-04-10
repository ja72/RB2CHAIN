!***************************************************************************
!
! FILE:         : TRIGONOMETRY.F90
!
! PURPOSE:      : MODULE DEFININING TRIGONOMETRIC CONSTANTS AND FUNCTIONS.
!
! DEPENDENCIES  : NONE
!
! HISTORY       : Written by John Alexiou - May 22, 2007
!
!  Modified by   Date   Change(s) made
!  -----------  ------  ------------------------------------------------

    MODULE constants
    use, intrinsic :: ISO_FORTRAN_ENV, only: INT32, INT64, REAL32, REAL64, REAL128
    IMPLICIT NONE
    
        INTEGER, PARAMETER :: dimensions = 3

        ! Select integer and float precision constants
        ! NOTE: Recommended to use constants defined in ISO_C_BINDING for
        !       compatibity with modern hardware (CPU+RAM) and software
        INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(2)
        INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(4)
        INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)
        INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
        INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(32)        

        ! Constants used for various calculations (double precision)
        REAL(real64), PARAMETER :: pi    = 3.1415926535897932D0
        REAL(real64), PARAMETER :: pid2  = 1.5707963267948966D0
        REAL(real64), PARAMETER :: twopi = 6.2831853071795864D0
        REAL(real64), PARAMETER :: pi_sq = 9.8690440108935862D0
        REAL(real64), PARAMETER :: rad   = 1.7453292519943296D-2
        REAL(real64), PARAMETER :: deg   = 57.29577951308232D0
        REAL(real64), PARAMETER :: div_pi= 0.31830988618379067D0
        REAL(real64), PARAMETER :: oned3 = 3.333333333333333D-1  
        REAL(real64), PARAMETER :: tend3 = 3.3333333333333333D0
        REAL(real64), PARAMETER :: deg_per_rad  = 180.0D0/pi, rad_per_deg = pi/180.0D0
        ! limits for handling infinity or 1/infinity
        REAL(real64), PARAMETER :: tiny = 1.0D-12
        REAL(real64), PARAMETER :: huge = 1.0D+18
    CONTAINS
    
        subroutine get_block_matrix(A,x,i,j)
        real(real64), intent(in) :: A(:,:)
        real(real64), intent(inout) :: x(:,:)
        integer(int32), intent(in) :: i,j
        integer(int32) :: ni,nj
            ni = size(x,1)
            nj = size(x,2)
            
            x(1:ni,1:nj) = A(1+ni*(i-1): ni*i, 1+nj*(j-1): nj*i)
            
        end subroutine
        
        subroutine set_block_matrix(A,x,i,j)
        real(real64), intent(inout) :: A(:,:)
        real(real64), intent(in) :: x(:,:)
        integer(int32), intent(in) :: i,j
        integer(int32) :: ni,nj
            ni = size(x,1)
            nj = size(x,2)
            
            A(1+ni*(i-1): ni*i, 1+nj*(j-1): nj*i) = x(1:ni,1:nj)
            
        end subroutine

    END MODULE