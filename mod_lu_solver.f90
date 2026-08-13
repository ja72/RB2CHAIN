!********************************************************************
!*    mod_lu_solver decomposition routines used by test_lu.f90      *
!*                                                                  *
!*                 F90 version by J-P Moreau, Paris                 *
!* ---------------------------------------------------              *
!* Reference:                                                       *
!*                                                                  *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,                *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge                   *
!*  University Press, 1986" [BIBLI 08].                             *
!*                                                                  *
!********************************************************************

    !DIR$ REAL:8
    MODULE mod_lu_solver
    implicit none

    interface solve
    !module procedure solve_linear_system
    module procedure solve_linear_system_full
    !module procedure solve_linear_system_mkl
    end interface

    CONTAINS

    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the mod_lu_solver *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
    pure Subroutine LUDCMP(A,N,INDX,D,CODE)
    IMPLICIT NONE
    integer, parameter :: nmax = 100
    real, parameter :: tiny = 1.5D-16

    real(8), intent(inout), dimension(N,N) :: A
    integer, intent(in) :: N
    integer, intent(out) :: D, CODE
    integer, intent(out), dimension(N) :: INDX
    !f2py depend(N) A, indx

    real(8)  :: AMAX, DUM, SUMM, VV(NMAX)
    INTEGER :: i, j, k, imax

    D=1; CODE=0

    DO I=1,N
        AMAX=0.d0
        DO J=1,N
            IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
        END DO ! j loop
        IF(AMAX.LT.TINY) THEN
            CODE = 1
            RETURN
        END IF
        VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
        DO I=1,J-1
            SUMM = A(I,J)
            !DIR$ SIMD
            DO K=1,I-1
                SUMM = SUMM - A(I,K)*A(K,J)
            END DO ! k loop
            A(I,J) = SUMM
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
            SUMM = A(I,J)
            !DIR$ SIMD
            DO K=1,J-1
                SUMM = SUMM - A(I,K)*A(K,J)
            END DO ! k loop
            A(I,J) = SUMM
            DUM = VV(I)*DABS(SUMM)
            IF(DUM.GE.AMAX) THEN
                IMAX = I
                AMAX = DUM
            END IF
        END DO ! i loop

        IF(J.NE.IMAX) THEN
            DO K=1,N
                DUM = A(IMAX,K)
                A(IMAX,K) = A(J,K)
                A(J,K) = DUM
            END DO ! k loop
            D = -D
            VV(IMAX) = VV(J)
        END IF

        INDX(J) = IMAX
        IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

        IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            !DIR$ SIMD
            DO I=J+1,N
                A(I,J) = A(I,J)*DUM
            END DO ! i loop
        END IF
    END DO ! j loop

    RETURN
    END subroutine LUDCMP


    !  ******************************************************************
    !  * Solves the set of N linear equations A . X = B.  Here A is     *
    !  * input, not as the matrix A but rather as its mod_lu_solver decomposition, *
    !  * determined by the routine LUDCMP. INDX is input as the permuta-*
    !  * tion vector returned by LUDCMP. B is input as the right-hand   *
    !  * side vector B, and returns with the solution vector X. A, N and*
    !  * INDX are not modified by this routine and can be used for suc- *
    !  * cessive calls with different right-hand sides. This routine is *
    !  * also efficient for plain matrix inversion.                     *
    !  ******************************************************************
    pure Subroutine LUBKSB(A, N, INDX, B)
    integer, intent(in) :: N
    real(8), intent(in), dimension(N,N) :: A
    integer, intent(in), dimension(N) :: INDX
    real(8), intent(inout), dimension(N) :: B
    integer :: II, I, J, LL

    real(8)  SUMM

    II = 0

    DO I=1,N
        LL = INDX(I)
        SUMM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
            !DIR$ SIMD
            DO J=II,I-1
                SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        ELSE IF(SUMM.NE.0.d0) THEN
            II = I
        END IF
        B(I) = SUMM
    END DO ! i loop

    DO I=N,1,-1
        SUMM = B(I)
        IF(I < N) THEN
            !DIR$ SIMD
            DO J=I+1,N
                SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        END IF
        B(I) = SUMM / A(I,I)
    END DO ! i loop

    RETURN
    END subroutine LUBKSB

    !************************************************************************
    !*  Summary:                                                            *
    !*      Solves a linear system of equations of the form y=A*x+b         *
    !*                                                                      *
    !*  Description:                                                        *
    !*      The system of equations is solved for some of the values in x   *
    !*      and some of the values in y based on the index array `pivot`,   *
    !*      the system size `n` and `k` the # of known `x`'s                *
    !*                                                                      *
    !*      `pivot` contains the index of all known `y`'s first and then    *
    !*      all known `x`'s next. See example below:                        *
    !*                                                                      *
    !*      The function returns the full vector `x`. Use `y=A*x+b` for `y` *
    !*                                                                      *
    !*  Example:                                                            *
    !*      | y_1 |   | 10  -2  5  | |  1  |   | -19 |                      *
    !*      | 10  | = | -2  11  1  | | x_2 | + | -45 |                      *
    !*      | y_3 |   | -1   2  5  | |  3  |   |  -6 |                      *
    !*                                                                      *
    !*      A = [[10.0,-2.0,5.0],[-2.0,27.0,1.0],[-1.0,2.0,5.0]]            *
    !*      b = [-19.0,-45.0,-6.0]               solution:                  *
    !*      x_known = [1.0,3.0]                 } x = [1.0, 2.0, 3.0]       *
    !*      y_known = [10.0]                    } y = [2.0, 10.0, 12.0]     *
    !*      pivot = [2,1,3]        => known: y(2), x(1), x(3)               *
    !*      n = 3                                                           *
    !*      k = 2                                                           *
    !*                                                                      *
    !*  Remarks:                                                            *
    !*      The function calls mod_lu_solver decomposition `ludcmp` and solver `lubksb`*
    !*      to solve for the unknowns.                                      *
    pure subroutine solve_linear_system(A, b, x, y, pivot, n, k)
    real(8), intent(in) :: A(n,n), b(n)
    real(8), intent(inout) :: x(n), y(n)
    integer(4), intent(in) :: pivot(n), n, k
    integer(4) :: u, code, d, indx(n-k)
    real(8) :: r(n-k), A1(n-k,n-k), A3(n-k,k), b1(n-k)

    u = n-k
    if(k>0) then

        A1 = A(pivot(1:u), pivot(1:u))
        A3 = A(pivot(1:u), pivot(u+1:n))
        b1 = b(pivot(1:u))

        r = y(pivot(1:u))-matmul(A3, x(pivot(u+1:n)))-b1
        call ludcmp(A1,u,indx,d,code)
        call lubksb(A1,u,indx,r)

        x(pivot(1:u)) = r
    else if(u>0) then
        r = y - b
        A1 = A
        call ludcmp(A1,u,indx,d,code)
        call lubksb(A1,u,indx,r)
        x = r
    end if
    y = matmul(A,x)+b
    end subroutine

    pure function solve_linear_system_full(A, y) result(x)
    ! Solve system of equations y=A*x for all unknowns in x.
    ! A is n�n, x and y are n�m.
    ! integer, intent(in),value :: n, m
    real(8), intent(in) :: A(:,:), y(:,:)
    real(8) :: x(size(A,1), size(y,2))
    integer :: u, code, d, indx(size(A,1))
    real(8) ::  A1(size(A,1), size(A,1)), r(size(A,1))
    integer :: i, n, m
    n = size(A,1)
    m = size(y,2)
    if( size(A,2) /= n .or. size(y,1) /= n ) then
        error stop "Error: Incompatible dimensions of A and y in solve_linear_system_full"
    end if
    A1 = A
    call ludcmp(A1,u,indx,d,code)
    do i=1,m
        r = y(:,i)
        call lubksb(A1,u,indx,r)
        x(:,i) = r
    end do
    end function

    !! requires additional depenency on linker of
    !! mkl_blas95_lp64.lib mkl_lapack95_lp64.lib
    !!
    !pure subroutine solve_linear_system_mkl(A, b, x, y, pivot, n, k)
    !use blas95
    !use lapack95
    !real(8), intent(in) :: A(n,n), b(n)
    !real(8), intent(inout) :: x(n), y(n)
    !integer(INT32), intent(in) :: pivot(n), n, k
    !integer(INT32) :: u, code, d, indx(n-k)
    !real(8) :: r(n-k), A1(n-k,n-k), A3(n-k,k), b1(n-k)
    !
    !    u = n-k
    !    if(k>0) then
    !
    !        A1 = A(pivot(1:u), pivot(1:u))
    !        A3 = A(pivot(1:u), pivot(u+1:n))
    !        b1 = b(pivot(1:u))
    !
    !        r = y(pivot(1:u))-matmul(A3, x(pivot(u+1:n)))-b1
    !        call gesv(A1,r,indx,code)
    !        !call ludcmp(A1,u,indx,d,code)
    !        !call lubksb(A1,u,indx,r)
    !
    !        x(pivot(1:u)) = r
    !    else if(u>0) then
    !        r = y - b
    !        A1 = A
    !        call gesv(A1,r,indx,code)
    !        !call ludcmp(A1,u,indx,d,code)
    !        !call lubksb(A1,u,indx,r)
    !        x = r
    !    end if
    !    y = matmul(A,x)+b
    !end subroutine


    subroutine test_linear_system()
    integer, parameter :: n = 3
    real(8) :: A(n,n), b(n), x(n), y(n)
    integer(4) :: pivot(n), k
    real(8) :: err_x(1), err_y(2)

    A = transpose(reshape([10d0,-2d0,5d0,-2d0,11d0,1d0,-1d0,2d0,5d0],[n,n]))
    b = [-26d0,-13d0,-10d0]
    x = [1d0, 0d99, 3d0]
    y = [0d99, 10d0, 0d99]

    k = 2               ! 2 known 'x'
    pivot = [2,1,3]     ! x([1,3]) : known, y([2]) : known

    call solve_linear_system(A,b,x,y,pivot,n,k)

    if( x(2)==2d0 .and. y(1)==-5d0 .and. y(3)==8d0 ) then
    else
        ! Check failed
        stop
    end if

    end subroutine

    END MODULE