!***************************************************************************
!
! FILE:         : VECTORALGEBRA.F90
!
! PURPOSE:      : MODULE DEFININING VECTOR TYPES WITH OPERATORS.
!
! DEPENDENCIES  : NONE
!
! HISTORY       : Written by John Alexiou - May 23, 2007
!
!  Modified by   Date   Change(s) made
!  -----------  ------  ------------------------------------------------
!
    MODULE vector_algebra
    USE constants
    implicit none
    
        integer, parameter :: nvec = 3     ! Number of components in a vector

        ! Constants for fixed-axis directions
        enum, bind(c)
			enumerator:: origin  = 0
            enumerator:: x_axis = 1
            enumerator:: y_axis = 2
            enumerator:: z_axis = 3
        end enum
        
        ! Generic names for functions acting on vectors
        interface operator (.dot.)
            module procedure inner_product
        end interface
        interface operator (.outer.)
            module procedure outer_product
        end interface
        interface operator (.x.)
            module procedure cross_product_vector
            module procedure cross_operator_matrix
        end interface
        interface rot_x
            module procedure rotate_vector_x, rotation_x
        end interface
        interface rot_y
            module procedure rotate_vector_y, rotation_y
        end interface
        interface rot_z
            module procedure rotate_vector_z, rotation_z
        end interface
        interface rot
            module procedure omega_to_rotation
        end interface

        ! Unit direction vectors
        real(real64), parameter :: o_(nvec) = [ 0.0d0, 0.0d0, 0.0d0 ]    ! origin: zero vector
        real(real64), parameter :: i_(nvec) = [ 1.0d0, 0.0d0, 0.0d0 ]    ! unit x-axis, hat[i]
        real(real64), parameter :: j_(nvec) = [ 0.0d0, 1.0d0, 0.0d0 ]    ! unit y-axis, hat[j]
        real(real64), parameter :: k_(nvec) = [ 0.0d0, 0.0d0, 1.0d0 ]    ! unit z-axis, hat[k]
        
        ! 3×3 identity matrix. Used for default rotations.
        real(real64), parameter :: z3_(nvec,nvec) = reshape( [0,0,0, 0,0,0, 0,0,0], [nvec,nvec])
        real(real64), parameter :: e3_(nvec,nvec) = reshape( [1,0,0, 0,1,0, 0,0,1], [nvec,nvec])
    CONTAINS
    !-- LINEAR ALGEBRA ------------------------------

        ! Create a n×m identity matrix. 
        ! NOTE: If m is ommited then create a n×n matrix
        pure function identity(n, m) result(eye)
        integer(INT32), intent(in) :: n
        integer(INT32), intent(in), optional :: m
        real(real64), allocatable :: eye(:,:)
        integer(INT32) :: i
            if( present(m) ) then
                allocate( eye(n,m) )
                eye = 0D0
                forall(i=1:min(n,m)) eye(i,i)=1D0
            else
                allocate( eye(n,n) )
                eye = 0D0
                forall(i=1:n) eye(i,i)=1D0
            end if                
        end function
        
        ! Create a vector of n zeros with 1 in the i-th row
        pure function element_vector(n,i) result(e)
        integer(INT32), intent(in) :: n, i
        real(real64), allocatable :: e(:)    
            allocate(e(n))
            e = 0
            e(i) = 1
        end function
        
        ! Create a n×(n-1) matrix from the n×n identity matrix, by removing column i
        ! NOTE: element_antivector is mutually exclusive to element_vector
        pure function element_antivector(n,i) result(u)
        integer(INT32), intent(in) :: n, i
        real(real64), allocatable :: u(:,:)    
        integer(INT32) :: k
            allocate(u(n,n-1))
            u = 0
            do k=1,i-1
                u(k,k) = 1
            end do
            do k=i+1,n
                u(k,k-1) = 1
            end do
        end function
        
        ! Returns the sum of the component squares
        pure function sumsq(a) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec)
            real(real64) :: r
        !code start
            r = dot_product(a,a)
        end function

        ! Returns the rms value of the vector components
        pure function magnitude(a) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec)
            real(real64) :: r
        !code start
            r = dsqrt( dot_product(a,a) )
        end function


        ! Returns the magnitude of the relative vector B-A
        pure function distance(a, b) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), b(nvec)
            real(real64) :: r
        !code start
            r = magnitude(b-a)
        end function

        ! Returns the rotated vector A about the X axis
        pure function rotate_vector_x(a,theta) result(z)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), theta
            real(real64) :: z(nvec)
        !code start
            z = [   a(1), &
                    dcos(theta)*a(2)-dsin(theta)*a(nvec), &
                    dsin(theta)*a(2)+dcos(theta)*a(nvec) ]
        end function
    
        ! Returns the rotated vector A about the Y axis
        pure function rotate_vector_y(a,theta) result(z)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), theta
            real(real64) :: z(nvec)
        !code start
            z = [   dcos(theta)*a(1)+dsin(theta)*a(nvec), &
                    a(2), &
                    -dsin(theta)*a(1)+dcos(theta)*a(nvec)]
        end function
    
        ! Returns the rotated vector A about the Z axis
        pure function rotate_vector_z(a,theta) result(z)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), theta
            real(real64) :: z(nvec)
        !code start
            z = [   dcos(theta)*a(1)-dsin(theta)*a(2), &
                    dsin(theta)*a(1)+dcos(theta)*a(2), &
                    a(nvec)]
        end function
        
        ! Elementary 3×3 rotation matrix about the x-axis
        pure function rotation_x(theta) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: theta
            real(real64) :: r(nvec,nvec)
        !code start
            r(1:3,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
            r(1:3,2) = [ 0.0d0, dcos(theta), dsin(theta) ]
            r(1:3,3) = [ 0.0d0, -dsin(theta), dcos(theta) ]
        end function
        
        ! Elementary 3×3 rotation matrix about the y-axis
        pure function rotation_y(theta) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: theta
            real(real64) :: r(nvec,nvec)
        !code start
            r(1:3,1) = [ dcos(theta), 0.0d0, -dsin(theta) ]
            r(1:3,2) = [ 0.0d0, 1.0d0, 0.0d0 ]
            r(1:3,3) = [ dsin(theta), 0.0d0, dcos(theta) ]
        end function
        
        ! Elementary 3×3 rotation matrix about the z-axis
        pure function rotation_z(theta) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: theta
            real(real64) :: r(nvec,nvec)
        !code start
            r(1:3,1) = [dcos(theta), dsin(theta), 0.0d0]
            r(1:3,2) = [-dsin(theta), dcos(theta), 0.0d0]
            r(1:3,3) = [0.0d0, 0.0d0, 1.0d0 ]
        end function

        ! Normalizes a vector A such that its new magnitude is ONE
        ! May return a ZERO vector if the magnitude is less than 1e-8
        pure function direction(a) result(z)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec)
            real(real64) :: z(nvec)
        !local variables
            real(real64) :: mag 
            mag = magnitude(a)
            if (dabs(mag)>tiny) then
                z = a/mag
            else
                z = a
            endif
        end function

        ! Angle formed by two vectors. A,B are vector types
        ! Uses the arc-cosine of the dot-product for the angle
        pure function vector_angle(a,b) result(theta)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), b(nvec)
            real(real64) :: theta
        !local variables
            real(real64) :: mag_a, mag_b, mag_axb
        !code start
            mag_a = magnitude(a)
            mag_b = magnitude(b)
            mag_axb = magnitude( a .x. b) 
            if( mag_a==0.0d0 .or. mag_b==0.0d0) then
                theta = 0.0d0
            else
                theta = dasin(mag_axb/(mag_a*mag_b))
            endif
        end function
        
        ! Inner product is the sum of the products of each element in two vectors
        pure function inner_product(a,b) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(:), b(:)
            real(real64) :: r
        !code start
            r = dot_product(a, b)
        end function
        
        ! Outer product of two vectors with elements c(i,j) = a(i)*b(j)
        pure function outer_product(a,b) result(c)
            real(real64), intent(in) :: a(:), b(:)
            real(real64) :: c(size(a),size(b))
            
            c = spread(a,dim=2, ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
        end function
        
        ! Specific version of the inner product for vectors of size 3
        pure function inner_product_vector(a,b) result(r)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), b(nvec)
            real(real64) :: r
        !code start
            r = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
        end function
    
        ! Outer product is a 3×3 matrix with elements c(i,j) = a(i)*b(j)
        pure function outer_product_vector(a,b) result(c)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), b(nvec)
            real(real64) :: c(nvec,nvec)
        !code start
            c(:,1) = a(:)*b(1)
            c(:,2) = a(:)*b(2)
            c(:,3) = a(:)*b(3)
            
            ! c = spread(a,dim=2, ncopies=3)*spread(b,dim=1,ncopies=3)
        end function

        ! Returns the cross product of two vectors A, B
        pure function cross_product_vector(a,b) result(c)
        implicit none
        !argument specification
            real(real64), intent(in) :: a(nvec), b(nvec)
            real(real64) :: c(nvec)
        !code start
            c = [   a(2)*b(3) - a(3)*b(2), &
                    a(3)*b(1) - a(1)*b(3), &
                    a(1)*b(2) - a(2)*b(1) ]
        end function

        ! Returns the 3×3 matrix cross product operator
        ! NOTE: Can be used as a prefix operator (.x. a) = CROSS_OPERATOR_MATRIX(a)
        pure function cross_operator_matrix(a) result(c)
        implicit none
        !argument specifications
            real(real64), intent(in) :: a(nvec)
            real(real64) :: c(nvec,nvec)
        !code start        
            c(:,1) = [ 0.0d0, a(3), -a(2) ]
            c(:,2) = [ -a(3), 0.0d0, a(1) ]
            c(:,3) = [ a(2), -a(1), 0.0d0 ]
        end function

        ! Returns the 3×3 rotation matrix based on the axis constant (1=X, 2=Y, 3=Z) and an angle
        pure function rotation_matrix(axis, theta) result(r)
        !argument specification
            integer(INT32), intent(in) :: axis 
            real(real64), intent(in) :: theta
            real(real64) :: r(nvec,nvec)
        !code start
            select case(axis)
            case(x_axis)
                r = rotation_x(theta)
            case(y_axis)
                r = rotation_y(theta)
            case(z_axis)
                r = rotation_z(theta)
            case default
                r = e3_
            end select
            
        end function
        
        ! Returns a unit vector pointing along the x, y or z axis depending on this axis constant (1=X, 2=Y, 3=Z, other=0)
        pure function direction_vector(axis) result(z)
        integer(INT32), intent(in) :: axis 
        real(real64) :: z(nvec)            
            select case(axis)
            case(x_axis)
                z = i_
            case(y_axis)
                z = j_
            case(z_axis)
                z = k_
            case default
                z = o_
            end select
        end function
                
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
        !*      The function calls LU decomposition `ludcmp` and solver `lubksb`*
        !*      to solve for the unknowns.                                      *
        pure subroutine solve_linear_system(A, b, x, y, pivot, n, k)
        use lu
        real(real64), intent(in) :: A(n,n), b(n)
        real(real64), intent(inout) :: x(n), y(n)
        integer(INT32), intent(in) :: pivot(n), n, k
        integer(INT32) :: u, code, d, indx(n-k)
        real(real64) :: r(n-k), A1(n-k,n-k), A3(n-k,k), b1(n-k)
                        
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
        
        pure subroutine solve_linear_system_mkl(A, b, x, y, pivot, n, k)
        use blas95
        use lapack95
        real(real64), intent(in) :: A(n,n), b(n)
        real(real64), intent(inout) :: x(n), y(n)
        integer(INT32), intent(in) :: pivot(n), n, k
        integer(INT32) :: u, code, d, indx(n-k)
        real(real64) :: r(n-k), A1(n-k,n-k), A3(n-k,k), b1(n-k)
                        
            u = n-k
            if(k>0) then
                
                A1 = A(pivot(1:u), pivot(1:u))
                A3 = A(pivot(1:u), pivot(u+1:n))
                b1 = b(pivot(1:u))
                
                r = y(pivot(1:u))-matmul(A3, x(pivot(u+1:n)))-b1
                call gesv(A1,r,indx,code)
                !call ludcmp(A1,u,indx,d,code)
                !call lubksb(A1,u,indx,r)
                
                x(pivot(1:u)) = r
            else if(u>0) then
                r = y - b
                A1 = A
                call gesv(A1,r,indx,code)
                !call ludcmp(A1,u,indx,d,code)
                !call lubksb(A1,u,indx,r)
                x = r
            end if                            
            y = matmul(A,x)+b
        end subroutine
        
        
        subroutine test_linear_system()
        integer, parameter :: n = 3
        real(real64) :: A(n,n), b(n), x(n), y(n)
        integer(INT32) :: pivot(n), k
        real(real64) :: err_x(1), err_y(2)
        
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
        
        !***********************************************************************
        ! Creates a 3×3 rotation matrix based on the axis/angle formula.
        ! It uses the rotational velocty vector omg, and the scalar time t
        !
        !   angle = |omg|*t
        !   axis = omg/|omg|
        !
        !   R = 1 + SIN(angle)*[axis×] + (1-COS(angle))*[axis×]*[axis×]
        !
        function omega_to_rotation(omg, t) result(R)
            real(real64), intent(in) :: omg(3), t
            real(real64) :: R(3,3), w,s,c,v
            w = norm2(omg)
            s = sin(w*t)
            c = cos(w*t)
            v = 1d0-c
            
            R(1,1) = w*w*c+v*omg(1)**2
            R(2,2) = w*w*c+v*omg(2)**2
            R(3,3) = w*w*c+v*omg(3)**2
            R(1,2) = v*omg(1)*omg(2)-w*omg(3)*s
            R(2,1) = v*omg(1)*omg(2)+w*omg(3)*s
            R(1,3) = v*omg(1)*omg(3)+w*omg(2)*s
            R(3,1) = v*omg(1)*omg(3)-w*omg(2)*s
            R(2,3) = v*omg(2)*omg(3)-w*omg(1)*s
            R(3,2) = v*omg(2)*omg(3)+w*omg(1)*s
            
            R = R/(w*w)
            
        end function

    end module