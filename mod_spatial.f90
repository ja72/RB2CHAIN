!***************************************************************************
    
    !DIR$ REAL:8    
    module mod_spatial
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_constants
    implicit none
    
    integer, private :: iidx, jjdx
    
    ! Constants for Common directions
    enum, bind(c)
        enumerator:: origin  = 0
        enumerator:: x_axis = 1
        enumerator:: y_axis = 2
        enumerator:: z_axis = 3
    end enum

    type :: vector3
        real(r8) :: data(3)
    contains
        procedure :: to_array => vec3_to_array
        procedure :: x => vec3_get_x
        procedure :: y => vec3_get_y
        procedure :: z => vec3_get_z
    end type
    
    type :: matrix3
        !real(r8) :: a11,a21,a31, a12,a22,a32, a13,a23,a33
        real(r8) :: data(3,3)
    contains
    procedure :: to_matrix => mat3_to_matrix
    procedure :: column => mat3_get_column
    procedure :: row => mat3_get_row
    procedure :: transpose => mat3_transpose
    end type
    
    type :: quaternion
        real(r8) :: data(4)
    contains
        procedure :: to_array => q_to_array
        procedure :: to_rotation => q_rotation_matrix
        procedure :: scalar => q_get_scalar
        procedure :: vector => q_get_vector
    end type    

    interface assignment (=)
    module procedure :: asgn_vec3_array, asgn_array_vec3
    end interface

    interface vector3
    module procedure :: vec3_from_axis, vec3_from_xyz
    end interface
    
    interface norm2
        module procedure :: vec3_magnitude
    end interface
    
    interface real
        module procedure :: vec3_to_array
    end interface

    interface operator (+)
    module procedure :: vec3_add
    end interface
    interface operator (-)
    module procedure :: vec3_neg
    module procedure :: vec3_sub
    end interface
    interface operator (*)
    module procedure :: vec3_scale
    module procedure :: vec3_scale2
    end interface
    interface operator (/)
    module procedure :: vec3_div
    end interface
    interface dot
    module procedure :: vec3_dot
    end interface
    interface cross
    module procedure :: array_cross, vec3_cross
    end interface    
    interface operator (.x.)
    module procedure :: array_cross, vec3_cross
    end interface

    type(vector3), parameter :: o_ = vector3([0.0_r8, 0.0_r8, 0.0_r8])
    type(vector3), parameter :: i_ = vector3([1.0_r8, 0.0_r8, 0.0_r8])
    type(vector3), parameter :: j_ = vector3([0.0_r8, 1.0_r8, 0.0_r8])
    type(vector3), parameter :: k_ = vector3([0.0_r8, 0.0_r8, 1.0_r8])
    
    type(matrix3), parameter :: zero_mat3 = matrix3( reshape([ (0.0_r8, jjdx=1,9)], [3,3]) )
    type(matrix3), parameter :: eye_mat3 = matrix3( reshape([ &
        1.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 1.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8], [3,3]) )
    
    interface assignment (=)
    module procedure :: asgn_mat3_matrix, asgn_matrix_mat3
    end interface

    interface operator (+)
    module procedure :: mat3_add
    end interface
    interface operator (-)
    module procedure :: mat3_neg
    module procedure :: mat3_sub
    end interface
    interface operator (*)
    module procedure :: mat3_scale
    module procedure :: mat3_scale2
    module procedure :: mat3_product_vec3, vec3_product_mat3
    module procedure :: mat3_product_mat3
    end interface
    interface operator (/)
    module procedure :: mat3_div
    end interface
    interface diag
    module procedure :: mat3_from_diag
    end interface
    interface cross
    module procedure :: array_cross_matrix, vec3_cross_matrix
    end interface
    interface outer
    module procedure :: vec3_outer
    end interface
    interface mmoi
    module procedure :: vec3_mmoi
    end interface

    interface assignment (=)
    module procedure :: asgn_q_array, asgn_array_q
    end interface

    interface quaternion
    module procedure :: q_from_array_scalar, q_from_vector
    end interface

    interface operator (+)
    module procedure :: q_add
    end interface
    interface operator (-)
    module procedure :: q_neg
    module procedure :: q_sub
    end interface
    interface operator (*)
    module procedure :: q_scale
    module procedure :: q_scale2
    end interface
    interface operator (/)
    module procedure :: q_div
    end interface
    interface dot
    module procedure :: q_dot
    end interface
    interface cross
    module procedure :: q_cross
    end interface
    interface operator (.x.)
    module procedure :: q_cross
    end interface
    
    interface operator ( .o. )
        module procedure q_product
    end interface

    interface rot
        module procedure q_axis_angle, q_vector_angle
        module procedure q_rotation_matrix
        module procedure q_rotation_zyx
        module procedure q_rotate_vector
    end interface
    
    interface norm2
        module procedure :: q_magnitude
    end interface
    
    interface unit
        module procedure :: q_normalize
    end interface
    
    type(quaternion), parameter :: q_eye = quaternion([0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8])
    
    contains
    
    !-- ARRAY ALGEBRA ----------------------------------------------------------
    
    pure function outer_product(a,b) result(c)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: c(size(a),size(b))
    c = spread(source = A, dim = 2, ncopies = size(b)) * &
        spread(source = B, dim = 1, ncopies = size(a))
    end function
    
    pure function array_cross(a,b) result(c)
    real(r8), intent(in) :: a(3), b(3)
    real(r8) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    end function    
    
    pure function array_cross_matrix(a) result(c)
    real(r8), intent(in) :: a(3)
    real(r8) :: c(3,3)
    c = 0.0_r8
    c(1,2) = -a(3)
    c(1,3) = a(2)
    c(2,1) = a(3)
    c(2,3) = -a(1)
    c(3,1) = -a(2)
    c(3,2) = a(1)
    end function
    
    !-- VECTOR3 ALGEBRA ----------------------------------------------------------
    
    pure function vec3_get_x(obj) result(res)
    class(vector3), intent(in) :: obj
    real(r8) :: res
    res = obj%data(1)   
    end function

    pure function vec3_get_y(obj) result(res)
    class(vector3), intent(in) :: obj
    real(r8) :: res
    res = obj%data(2)   
    end function

    pure function vec3_get_z(obj) result(res)
    class(vector3), intent(in) :: obj
    real(r8) :: res
    res = obj%data(3)   
    end function

    pure subroutine asgn_vec3_array(v,a)
    type(vector3), intent(out) :: v
    real(8), intent(in) :: a(3)
    v = vec3_from_array(a)
    end subroutine

    pure subroutine asgn_array_vec3(a,v)
    real(8), intent(out) :: a(3)
    type(vector3), intent(in) :: v
    a = vec3_to_array(v)
    end subroutine
    
    pure function vec3_from_axis(axis) result(v)
    integer, intent(in) :: axis    
    type(vector3) :: v
    
        select case(axis)
        case (x_axis)
            v%data(1) = 1.0_r8
        case (y_axis)
            v%data(2) = 1.0_r8
        case (z_axis)
            v%data(3) = 1.0_r8
        end select
    end function
    
    pure function vec3_from_array(a) result(v)
    real(r8), intent(in) :: a(3)
    type(vector3) :: v
    v%data= a
    end function
    pure function vec3_from_xyz(x,y,z) result(v)
    real(r8), intent(in) :: x,y,z
    type(vector3) :: v
    v%data= [x,y,z]
    end function

    pure function vec3_to_array(v) result(a)
    real(r8) :: a(3)
    class(vector3), intent(in) :: v
    a = v%data
    end function

    pure function vec3_add(a,b) result(c)
    ! c = a + b
    type(vector3), intent(in) :: a,b
    type(vector3) :: c
    c%data = a%data + b%data
    end function

    pure function vec3_neg(a) result(c)
    ! c = -a
    type(vector3), intent(in) :: a
    type(vector3) :: c
    c%data = -a%data
    end function

    pure function vec3_sub(a,b) result(c)
    ! c = a - b
    type(vector3), intent(in) :: a,b
    type(vector3) :: c
    c%data = a%data - b%data
    end function

    pure function vec3_scale(a,b) result(c)
    ! c = a * b
    real(r8), intent(in) :: a
    type(vector3), intent(in) :: b
    type(vector3) :: c
    c%data = a * b%data
    end function

    pure function vec3_scale2(b,a) result(c)
    ! c = a * b
    type(vector3), intent(in) :: b
    real(r8), intent(in) :: a
    type(vector3) :: c
    c%data = a * b%data
    end function

    pure function vec3_div(a,b) result(c)
    ! c = a / b
    type(vector3), intent(in) :: a
    real(r8), intent(in) :: b
    type(vector3) :: c
    c%data = a%data / b
    end function
    
    pure function vec3_magnitude(v) result(s)
    type(vector3), intent(in) :: v
    real(r8) :: s
    s = norm2(v%data)
    end function
    
    pure function vec3_dot(a,b) result(s)
    type(vector3), intent(in) :: a,b
    real(r8) :: s
        s = dot_product(a%data, b%data)
    end function
    
    pure function vec3_cross(a, b) result(v)
    type(vector3), intent(in) :: a,b
    type(vector3) :: v
        v%data = array_cross(a%data, b%data)
    end function
    pure function vec3_cross_matrix(a) result(m)
    type(vector3), intent(in) :: a
    type(matrix3) :: m
        m%data = array_cross_matrix(a%data)
    end function

    pure function vec3_outer(a,b) result(m)
    type(vector3), intent(in) :: a,b
    type(matrix3) :: m
    m%data = reshape( &
        [a%data(1)*b%data(1), a%data(1)*b%data(2), a%data(1)*b%data(3), &
         a%data(2)*b%data(1), a%data(2)*b%data(2), a%data(2)*b%data(3), &
         a%data(3)*b%data(1), a%data(3)*b%data(2), a%data(3)*b%data(3)], [3,3] )
    end function

    pure function vec3_mmoi(v, negative) result(m)
    type(vector3), intent(in) :: v
    logical, intent(in), optional :: negative
    type(matrix3) :: m
    real(r8) :: xx,yy,zz,xy,zx,yz
    integer :: sign

    !tex: Inertia matrix
    ! $$-[\vec{v}\times][\vec{v}\times] = \begin{bmatrix}y^{2}+z^{2} & -xy & -zx\\
    !-xy & x^{2}+z^{2} & -yz\\
    !-zx & -yz & x^{2}+z^{2}
    !\end{bmatrix}$$

    if(present(negative) .and. negative) then
        sign = -1
    else
        sign  = 1
    end if
    !DIR$ FMA
    xx =  sign*v%data(1)**2
    yy =  sign*v%data(2)**2
    zz =  sign*v%data(3)**2
    xy = -sign*v%data(1)*v%data(2)
    yz = -sign*v%data(2)*v%data(3)
    zx = -sign*v%data(3)*v%data(1)

    !m%a11 = yy+zz
    !m%a22 = xx+zz
    !m%a33 = yy+xx
    !m%a12 = xy
    !m%a13 = zx
    !m%a23 = yz
    !m%a21 = xy
    !m%a31 = zx
    !m%a32 = yz
    
    m%data = reshape( [yy+zz, xy, zx, xy, xx+zz, yz, zx, yz, yy+xx], [3,3] )

    end function
    
    
    !-- MATRIX3 ALGEBRA ----------------------------------------------------------
    pure subroutine asgn_mat3_matrix(m,a)
    type(matrix3), intent(out) :: m
    real(r8), intent(in) :: a(3,3)
    m = mat3_from_matrix(a)
    end subroutine

    pure subroutine asgn_matrix_mat3(a,m)
    real(r8), intent(out) :: a(3,3)
    type(matrix3), intent(in) :: m
    a = mat3_to_matrix(m)
    end subroutine
    
    pure function mat3_from_matrix(a) result(m)
    real(r8), intent(in) :: a(3,3)
    type(matrix3) :: m
    m%data = a
    end function

    pure function mat3_to_matrix(m) result(a)
    real(r8) :: a(3,3)
    class(matrix3), intent(in) :: m
    a = m%data
    end function
    
    pure function mat3_from_diag(d) result(m)
    real(r8), intent(in) :: d(3)
    type(matrix3) :: m
        m%data = 0.0_r8
        m%data(1,1) = d(1)
        m%data(2,2) = d(2)
        m%data(3,3) = d(3)
    end function
    
    pure function mat3_add(a,b) result(c)
    type(matrix3), intent(in) :: a,b
    type(matrix3) :: c
    !DIR$ FMA
    c%data = a%data + b%data
    end function
    pure function mat3_neg(a) result(c)
    type(matrix3), intent(in) :: a
    type(matrix3) :: c
    !DIR$ FMA
    c%data = -a%data
    end function
    pure function mat3_sub(a,b) result(c)
    type(matrix3), intent(in) :: a,b
    type(matrix3) :: c
    !DIR$ FMA
    c%data = a%data - b%data
    end function
    
    pure function mat3_scale(a,b) result(c)
    real(r8), intent(in) :: a
    type(matrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%data = a * b%data
    end function

    pure function mat3_scale2(b,a) result(c)
    type(matrix3), intent(in) :: b
    real(r8), intent(in) :: a
    type(matrix3) :: c
    !DIR$ FMA
    c%data = b%data * a
    end function

    pure function mat3_div(a,b) result(c)
    type(matrix3), intent(in) :: a
    real(r8), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%data = a%data / b
    end function

    pure function mat3_product_vec3(a,b) result(c)
    type(matrix3), intent(in) :: a
    type(vector3), intent(in) :: b
    type(vector3) :: c
    !DIR$ FMA
    c%data = matmul(a%data, b%data)
    end function

    pure function vec3_product_mat3(b,a) result(c)
    type(vector3), intent(in) :: b
    type(matrix3), intent(in) :: a
    type(vector3) :: c
    !DIR$ FMA
    c%data = matmul(b%data, a%data)
    end function
    pure function mat3_product_mat3(a,b) result(c)
    type(matrix3), intent(in) :: a
    type(matrix3), intent(in) :: b
    type(matrix3) :: c
    !DIR$ FMA
    c%data = matmul(a%data, b%data)
    end function
    
    

    pure function mat3_transpose(m) result(w)
    class(matrix3), intent(in) :: m
    type(matrix3) :: w
    w%data = transpose(m%data)
    end function

    pure function mat3_get_column(m, col) result(v)
    class(matrix3), intent(in) :: m
    integer, intent(in) :: col
    type(vector3) :: v
    select case(col)
    case(1)
        v = vector3(m%data(:,1))
    case(2)
        v = vector3(m%data(:,2))
    case(3)
        v = vector3(m%data(:,3))
    case default
        v = vector3(0.0_r8, 0.0_r8, 0.0_r8)
    end select
    end function
    
    pure function mat3_get_row(m, row) result(v)
    class(matrix3), intent(in) :: m
    integer, intent(in) :: row
    type(vector3) :: v
    select case(row)
    case(1)
        v = vector3(m%data(1,:))
    case(2)
        v = vector3(m%data(2,:))
    case(3)
        v = vector3(m%data(3,:))
    case default
        v = vector3(0.0_r8, 0.0_r8, 0.0_r8)
    end select
    end function
    
    pure subroutine asgn_q_array(q,a)
    type(quaternion), intent(out) :: q
    real(r8), intent(in) :: a(4)
    q%data = a
    end subroutine

    pure subroutine asgn_array_q(a,q)
    real(r8), intent(out) :: a(4)
    type(quaternion), intent(in) :: q
    a = q%data
    end subroutine
    

    !- QUATERNION ALGEBRA -------------------------------------------------
    pure function q_from_array_scalar(a, s) result(q)
    real(r8), intent(in) :: a(3), s
    type(quaternion) :: q
        q%data = [a, s]
    end function
    
    pure function q_to_array(q) result(a)
    class(quaternion), intent(in) :: q
    real(r8) :: a(4)
        a = q%data
    end function
    
    pure function q_from_vector(v) result(q)
    type(vector3), intent(in) :: v
    type(quaternion) :: q
        q%data = [v%data, 0.0_r8]
    end function
    
    pure function q_get_scalar(q) result(s)
    class(quaternion), intent(in) :: q
    real(r8) :: s
    s = q%data(4)
    end function
    
    pure function q_get_vector(q) result(v)
    class(quaternion), intent(in) :: q
    type(vector3) :: v
    v%data = q%data(1:3)
    end function
    
    pure function q_axis_angle(axis, angle) result(q)
    integer, intent(in) :: axis    
    real(r8), intent(in) :: angle
    type(quaternion) :: q
        q = q_vector_angle(vector3(axis),angle)
    end function
    
    pure function q_vector_angle(axis, angle) result(q)
    type(vector3), intent(in) :: axis    
    real(r8), intent(in) :: angle
    type(quaternion) :: q
    real(r8) :: s, c, m
        if( abs(angle) > tiny) then
            m = norm2(axis%data)
            s = sin(angle/2)
            c = cos(angle/2)
            q%data = [s*axis%data/m, c]
        else
            q%data = [0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8]
        end if
    end function
    
    pure function q_rotation_zyx(yaw,pitch,roll) result(q)
    ! Define quaternion rotation from axis angle
    real(r8), intent(in) :: yaw,pitch,roll
    type(quaternion) :: q
    q = rot(z_axis, yaw) .o. rot(y_axis, pitch) .o. rot(x_axis, roll)
    end function
    
    
    pure function q_magnitude(a) result(s)
        type(quaternion), intent(in) :: a
        real(r8) :: s
        s = norm2(a%data)
    end function
    
    pure function q_normalize(a) result(q)
    type(quaternion), intent(in) :: a
    type(quaternion) :: q
    real(r8) :: m
        m = norm2(a%data)
        if( m >= tiny ) then
            q%data = a%data/m
        else
            q%data = a%data
        end if
    end function
    
    !-- QUATERNION ALGEBRA
    
    pure function q_add(a,b) result(c)
    type(quaternion), intent(in) :: a,b
    type(quaternion) :: c
    c%data = a%data + b%data
    end function

    pure function q_neg(a) result(c)
    type(quaternion), intent(in) :: a
    type(quaternion) :: c
    c%data = -a%data
    end function

    pure function q_sub(a,b) result(c)
    type(quaternion), intent(in) :: a,b
    type(quaternion) :: c
    c%data = a%data - b%data
    end function

    pure function q_scale(a,b) result(c)
    real(r8), intent(in) :: a
    type(quaternion), intent(in) :: b
    type(quaternion) :: c
    c%data = a*b%data
    end function

    pure function q_scale2(b,a) result(c)
    type(quaternion), intent(in) :: b
    real(r8), intent(in) :: a
    type(quaternion) :: c
    c%data = a*b%data
    end function

    pure function q_div(a,b) result(c)
    type(quaternion), intent(in) :: a
    real(r8), intent(in) :: b
    type(quaternion) :: c
    c%data = a%data/b
    end function    
        
    pure function q_dot(a,b) result(s)
        type(quaternion), intent(in) :: a,b
        real(r8) :: s
        s = dot_product(a%data, b%data)
    end function
    
    pure function q_cross(a,b) result(q)
        type(quaternion), intent(in) :: a,b
        type(quaternion) :: q
        q%data = [ array_cross(a%data(1:3), b%data(1:3)), 0.0_r8 ]
    end function
    
    pure function q_product(a,b) result(q)
        type(quaternion), intent(in) :: a,b
        type(quaternion) :: q
        real(r8) :: ax,ay,az,aw
        real(r8) :: bx,by,bz,bw
        ax = a%data(1)
        ay = a%data(2)
        az = a%data(3)
        aw = a%data(4)
        bx = b%data(1)
        by = b%data(2)
        bz = b%data(3)
        bw = b%data(4)
        !DIR$ FMA
        q%data = [ &
            aw*bx+ax*bw+ay*bz-az*by, &
            aw*by-ax*bz+ay*bw+az*bx, &
            aw*bz+ax*by-ay*bx+az*bw, &
            aw*bw-ax*bx-ay*by-az*bz ]
            
    end function
    
    pure function q_rotate_vector(q, v) result(v_rot)
    class(quaternion), intent(in) :: q
    type(vector3), intent(in) :: v
    type(vector3) :: v_rot
    real(r8) :: q_v(3), q_s, qxv(3), qxqxv(3)
    
        q_v = q%vector()
        q_s = q%scalar()
        
        qxv = cross(q_v, v%data)
        qxqxv = cross(q_v, qxv)
        
        v_rot%data = v%data + 2*q_s*qxv + 2*qxqxv
    
    end function
    
    pure function q_rotation_matrix(q) result(R)
    ! Extract rotation matrix R from a quaternion
    ! Assumes orientation is a unit quaternion
    class(quaternion), intent(in) :: q
    type(matrix3) :: R !, X, XX    
    real(r8) :: x,y,z,w,xx,yy,zz
    
    !tex: Rotation matrix from quaternion $\boldsymbol{q}=[\vec{q_v},\,q_s]$
    !$$ \mathbf{R} = \mathbf{1} + 2\, q_s\, [\vec{q_v}\times] + 2\, [\vec{q_v}\times][\vec{q_v}\times]  $$
    
   
    x = q%data(1)
    y = q%data(2)
    z = q%data(3)
    w = q%data(4)
    xx = x*x
    yy = y*y
    zz = z*z
    
    !DIR$ FMA 
    R%data = reshape( [ &
        1 - 2*(yy+zz), -2*(w*z-x*y), -2*(-x*z-w*y), &
        -2*(-w*z-x*y), 1 - 2*(xx+zz), -2*( w*x-y*z), &
        -2*(-x*z+w*y), -2*(-w*x-y*z), 1 - 2*(xx+yy) ], [3,3] )
    
    end function
    
    end module mod_spatial