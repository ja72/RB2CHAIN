!***************************************************************************
!   mod_dynamics_spatial                                                    
!                                                                           
!   structures and functions to do screw theory (3d rigid body mechanics)   
!---------------------------------------------------------------------------
    
    !DIR$ REAL:8
    module mod_dynamics_spatial
    use mod_spatial
    implicit none
    
    integer, private :: iidx, jjdx
        
    type :: vector6
        real(r8) :: data(6)
    contains
        procedure :: to_array => vec6_to_array
    end type
    
    type :: matrix6
        real(r8) :: data(6,6)
    contains
        procedure :: to_array => mat6_to_array
    end type
    
    interface vector6
    module procedure vec6_from_vectors
    end interface
    
    interface matrix6
    module procedure mat6_from_matrices
    end interface
    
    interface twist
    module procedure vec6_twist_at, vec6_pure_twist
    end interface
    interface wrench
    module procedure vec6_wrench_at, vec6_pure_wrench
    end interface
    
    interface operator(+)
    module procedure add_vec6, add_mat6, add_mat6_scalar, add_scalar_mat6
    end interface
    interface operator(-)
    module procedure neg_vec6, neg_mat6
    module procedure sub_vec6, sub_mat6, sub_mat6_scalar, sub_scalar_mat6
    end interface
    interface operator(*)
    module procedure mul_scalar_vec6, mul_scalar_mat6
    module procedure mul_vec6_scalar, mul_mat6_scalar
    module procedure mat6_mul_vec6, vec6_mul_mat6
    end interface
    interface operator(/)
    module procedure div_vec6_scalar, div_mat6_scalar
    end interface
    
    interface spi
        module procedure mat6_spi, mat6_spi_rot
    end interface
    interface spm
        module procedure mat6_spm, mat6_spm_rot
    end interface
    
    interface cross_twist
        module procedure vec6_cross_twist, vec6_cross_twist_matrix
    end interface
    interface cross_wrench
        module procedure vec6_cross_wrench, vec6_cross_wrench_matrix
    end interface
    interface operator(.xt.)
        module procedure vec6_cross_twist
    end interface
    interface operator(.xw.)
        module procedure vec6_cross_wrench
    end interface
    
    type(vector6), parameter :: zero_vec6 = vector6([(0.0_r8, jjdx=1,6)])
    type(matrix6), parameter :: zero_mat6 = matrix6(reshape([(0.0_r8, jjdx=1,36)], [6,6]))
    type(matrix6), parameter :: eye_mat6 = matrix6( reshape([ &
        1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8 ], [6,6]) )
    
    
    type :: link
        real(r8) :: mass
        type(matrix3) :: inertia
        type(vector3) :: base_cg_pos
        type(quaternion) :: base_cg_ori
        type(vector3) :: next_pos
        type(quaternion) :: next_ori
        type(vector3) :: joint_axis
        integer :: joint_type
        integer :: joint_driver
        real(r8) :: motor
    end type
    
    type :: pendulum(n)
        integer, len :: n
        type(vector3) :: gee = vector3([0.0_r8, -10.0_r8, 0.0_r8])
        type(link) :: body(n)
        integer :: parents(n)
        real(r8) :: q(n), qp(n), qpp(n), tau(n)
    contains
    procedure :: init => pendulum_init
    procedure :: solve => pendulum_solve
    end type
        
    ! Joint types
    enum, bind(c)
        enumerator:: revolute = 1       ! 1-DOF pin joint
        enumerator:: prismatic = 2      ! 1-DOF slider joint
    end enum
    
    ! Joint drivers
    enum, bind(c)
        enumerator:: known_torque = 1   ! Inverse dynamics for joint, acceleration is found from torque
        enumerator:: known_motion = 2   ! Forward dynamics for joint, torque is found from acceleration
    end enum
    
    
    
    contains
    
    pure function vec6_from_vectors(v1, v2) result(vec)
    type(vector3), intent(in) :: v1, v2
    type(vector6) :: vec
    vec%data(1:3) = v1%data
    vec%data(4:6) = v2%data
    end function
        
    pure function vec6_to_array(vec) result(arr)
    class(vector6), intent(in) :: vec
    real(r8) :: arr(6)
    arr = vec%data
    end function
    
    
    !-- VECTOR6 FUNCTIONS
    pure function vec6_twist_at(value, point, pitch) result(res)
    type(vector3), intent(in) :: value
    type(vector3), intent(in) :: point
    real(r8), intent(in), optional :: pitch
    type(vector6) :: res
    real(r8) :: top(3), bot(3)
    ! Implementation for computing twist at a specific point
    bot = value
    if( present(pitch) ) then
        top = value * pitch + cross(point, value)
    else
        top = cross(point, value)
    end if
    res%data = [ top, bot ]
    end function
    pure function vec6_pure_twist(value) result(res)
    type(vector3), intent(in) :: value
    type(vector6) :: res
    res % data = [ value%data, 0.0_r8, 0.0_r8, 0.0_r8 ]
    end function
    
    pure function vec6_wrench_at(value, point, pitch) result(res)
    type(vector3), intent(in) :: value
    type(vector3), intent(in) :: point
    real(r8), intent(in), optional :: pitch
    type(vector6) :: res
    real(r8) :: top(3), bot(3)
    ! Implementation for computing wrench at a specific point
    top = value
    if( present(pitch) ) then
        bot = value * pitch + cross(point, value)
    else
        bot = cross(point, value)
    end if
    res%data = [ top, bot ]
    end function
    pure function vec6_pure_wrench(value) result(res)
    type(vector3), intent(in) :: value
    type(vector6) :: res
    res % data = [ 0.0_r8, 0.0_r8, 0.0_r8, value%data ]
    end function
    
    pure function neg_vec6(a) result(res)
    type(vector6), intent(in) :: a
    type(vector6) :: res
    res%data = -a%data
    end function
    pure function add_vec6(a,b) result(res)
    type(vector6), intent(in) :: a, b
    type(vector6) :: res
    res%data = a%data + b%data
    end function
    pure function sub_vec6(a,b) result(res)
    type(vector6), intent(in) :: a, b
    type(vector6) :: res
    res%data = a%data - b%data
    end function
    
    pure function mul_scalar_vec6(s, a) result(res)
    real(r8), intent(in) :: s
    type(vector6), intent(in) :: a
    type(vector6) :: res
    res%data = s * a%data
    end function
    
    pure function mul_vec6_scalar(a, s) result(res)
    type(vector6), intent(in) :: a
    real(r8), intent(in) :: s
    type(vector6) :: res
    res%data = a%data * s
    end function
    
    pure function div_vec6_scalar(a, s) result(res)
    type(vector6), intent(in) :: a
    real(r8), intent(in) :: s
    type(vector6) :: res
    res%data = a%data / s
    end function
    
    !-- MATRIX6 FUNCTIONS ---------------------------------------
    pure function mat6_from_matrices(m1, m2, m3, m4) result(mat)
    type(matrix3), intent(in) :: m1, m2, m3, m4
    type(matrix6) :: mat
    mat%data(1:3,1:3) = m1%data
    mat%data(1:3,4:6) = m2%data
    mat%data(4:6,1:3) = m3%data
    mat%data(4:6,4:6) = m4%data
    end function
    pure function mat6_to_array(vec) result(arr)
    class(matrix6), intent(in) :: vec
    real(r8) :: arr(6,6)
    arr = vec%data
    end function
    pure function neg_mat6(a) result(res)
    type(matrix6), intent(in) :: a
    type(matrix6) :: res
    res%data = -a%data
    end function
    pure function add_mat6(a,b) result(res)
    type(matrix6), intent(in) :: a, b
    type(matrix6) :: res
    res%data = a%data + b%data
    end function
    pure function sub_mat6(a,b) result(res)
    type(matrix6), intent(in) :: a, b
    type(matrix6) :: res
    res%data = a%data - b%data
    end function
    
    pure function add_mat6_scalar(a, s) result(res)
    type(matrix6), intent(in) :: a
    real(r8), intent(in) :: s
    type(matrix6) :: res
    res%data = a%data + s * eye_mat6%data
    end function
    pure function add_scalar_mat6(s, a) result(res)
    real(r8), intent(in) :: s
    type(matrix6), intent(in) :: a
    type(matrix6) :: res
    res%data = s * eye_mat6%data + a%data
    end function
    pure function sub_mat6_scalar(a, s) result(res)
    type(matrix6), intent(in) :: a
    real(r8), intent(in) :: s
    type(matrix6) :: res
    res%data = a%data - s * eye_mat6%data
    end function
    pure function sub_scalar_mat6(s, a) result(res)
    real(r8), intent(in) :: s
    type(matrix6), intent(in) :: a
    type(matrix6) :: res
    res%data = s * eye_mat6%data - a%data
    end function

    pure function mul_scalar_mat6(s, a) result(res)
    real(r8), intent(in) :: s
    type(matrix6), intent(in) :: a
    type(matrix6) :: res
    res%data = s * a%data
    end function
    
    pure function mul_mat6_scalar(a, s) result(res)
    type(matrix6), intent(in) :: a
    real(r8), intent(in) :: s
    type(matrix6) :: res
    res%data = a%data * s
    end function
    
    pure function div_mat6_scalar(a, s) result(res)
    type(matrix6), intent(in) :: a
    real(r8), intent(in) :: s
    type(matrix6) :: res
    res%data = a%data / s
    end function
    
    pure function mat6_mul_vec6(a, b) result(res)
    type(matrix6), intent(in) :: a
    type(vector6), intent(in) :: b
    type(vector6) :: res
    res%data = matmul(a%data, b%data)
    end function
    
    pure function vec6_mul_mat6(a, b) result(res)
    type(vector6), intent(in) :: a
    type(matrix6), intent(in) :: b
    type(vector6) :: res
    res%data = matmul(transpose(b%data), a%data)
    end function
    

    !-- CROSS PRODUCTS
    pure function vec6_cross_twist(a, b) result(res)
    type(vector6), intent(in) :: a, b
    type(vector6) :: res
    real(r8) :: a_top(3), a_bot(3), b_top(3), b_bot(3)
        a_top = a%data(1:3) 
        a_bot = a%data(4:6)
        b_top = b%data(1:3) 
        b_bot = b%data(4:6)
        !tex: Twist Cross Product
        ! $$ \begin{bmatrix} v \\ \omega \end{bmatrix} \times \begin{bmatrix} g \\ \alpha \end{bmatrix} = 
        ! \begin{bmatrix} \omega \times g + v \times \alpha \\ \omega \times \alpha \end{bmatrix} $$
        
        res%data(1:3) = cross(a_bot, b_top) + cross(a_top, b_bot)
        res%data(4:6) = cross(a_bot, b_bot)

    end function
    
    pure function vec6_cross_wrench(a, b) result(res)
    type(vector6), intent(in) :: a, b
    type(vector6) :: res
    real(r8) :: a_top(3), a_bot(3), b_top(3), b_bot(3)
        a_top = a%data(1:3) 
        a_bot = a%data(4:6)
        b_top = b%data(1:3) 
        b_bot = b%data(4:6)
        !tex: Wrench Cross Product
        ! $$ \begin{bmatrix} v \\ \omega \end{bmatrix} \times \begin{bmatrix} f \\ \tau \end{bmatrix} = 
        ! \begin{bmatrix} \omega \times f  \\ v \times f + \omega \times \tau \end{bmatrix} $$

        
        res%data(1:3) = cross(a_bot, b_top) 
        res%data(4:6) = cross(a_top, b_top) + cross(a_bot, b_bot)

    end function
    
    pure function vec6_cross_twist_matrix(a) result(res)
    type(vector6), intent(in) :: a
    type(matrix6) :: res
    real(r8) :: a_top(3), a_bot(3)
        a_top = a%data(1:3) 
        a_bot = a%data(4:6)
        
        !tex: Twist Cross Product Matrix
        ! $$ \begin{bmatrix} v \\ \omega \end{bmatrix} \times = 
        ! \begin{bmatrix} \omega \times & v \times \\ 0 & \omega \times \end{bmatrix} $$
        
        res%data(1:3,1:3) = cross(a_bot)
        res%data(1:3,4:6) = cross(a_top)
        res%data(4:6,1:3) = 0.0_r8
        res%data(4:6,4:6) = cross(a_bot)
    end function

    pure function vec6_cross_wrench_matrix(a) result(res)
    type(vector6), intent(in) :: a
    type(matrix6) :: res
    real(r8) :: a_top(3), a_bot(3)
        a_top = a%data(1:3) 
        a_bot = a%data(4:6)
        
        !tex: Wrench Cross Product Matrix
        ! $$ \begin{bmatrix} v \\ \omega \end{bmatrix} \times = 
        ! \begin{bmatrix}  \omega \times & 0 \\ v \times  & \omega \times \end{bmatrix} $$
        
        res%data(1:3,1:3) = cross(a_bot)
        res%data(1:3,4:6) = 0.0_r8
        res%data(4:6,1:3) = cross(a_top)
        res%data(4:6,4:6) = cross(a_bot)
    end function
    
    !-- SPATIAL INERTIA FUNCTIONS ---------------------------------------
    pure function mat6_spi_rot(mass, inertia, cg, ori) result(res)
    real(r8), intent(in) :: mass
    type(matrix3), intent(in) :: inertia
    type(vector3), intent(in) :: cg
    type(quaternion), intent(in) :: ori
    type(matrix6) :: res
    type(matrix3) :: rot_inertia, R
    R = rot(ori)
    rot_inertia%data = matmul(R%data, matmul(inertia%data, transpose(R%data)))
    res = mat6_spi(mass, rot_inertia, cg)
    end function
    
    pure function mat6_spm(mass, inv_inertia, cg) result(res)
    real(r8), intent(in) :: mass
    type(matrix3), intent(in) :: inv_inertia
    type(vector3), intent(in) :: cg
    type(matrix6) :: res
    real(r8) :: cgx(3,3)
    cgx = cross(cg%data)
    
    res%data(1:3, 1:3) = (1/mass) * eye_mat3%data + matmul(transpose(cgx), matmul(inv_inertia%data, cgx))
    res%data(4:6, 1:3) = matmul(inv_inertia%data, cgx)
    res%data(1:3, 4:6) = matmul(transpose(cgx), inv_inertia%data)
    res%data(4:6, 4:6) = inv_inertia%data
    
    end function
    pure function mat6_spm_rot(mass, inv_inertia, cg, rot) result(res)
    real(r8), intent(in) :: mass
    type(matrix3), intent(in) :: inv_inertia
    type(vector3), intent(in) :: cg
    type(matrix3), intent(in) :: rot
    type(matrix6) :: res
    type(matrix3) :: rot_inertia
    rot_inertia%data = matmul(rot%data, matmul(inv_inertia%data, transpose(rot%data)))
    res = mat6_spm(mass, rot_inertia, cg)
    end function    
    pure function mat6_spi(mass, inertia, cg) result(res)
    real(r8), intent(in) :: mass
    type(matrix3), intent(in) :: inertia
    type(vector3), intent(in) :: cg
    type(matrix6) :: res
    real(r8) :: mcgx(3,3), cgx(3,3)
    cgx = cross(cg%data)
    mcgx = mass * cgx
    
    res%data(1:3, 1:3) = mass * eye_mat3%data
    res%data(4:6, 1:3) = mcgx
    res%data(1:3, 4:6) = transpose(mcgx)
    res%data(4:6, 4:6) = inertia%data - matmul(mcgx, cgx)
    
    end function
    
    !-- PENDULUM FUNCTIONS ---------------------------------------
    subroutine pendulum_init(this, rb)
    class(pendulum(*)), intent(inout) :: this
    type(link), pointer, intent(in) :: rb
    integer :: i
    this%body(:) = rb
    do i=1, this%n
        this%parents(i) = i - 1
        this%q(i) = 0.0_r8
        this%qp(i) = 0.0_r8
        this%qpp(i) = 0.0_r8
        this%tau(i) = 0.0_r8
    end do
    end subroutine
    
    pure subroutine pendulum_solve(this)
    class(pendulum(*)), intent(inout) :: this
    type(link)       :: rb
    type(vector3)    :: prev_pos
    type(quaternion) :: prev_ori
    type(vector6)    :: v_prev, a_prev
    integer          :: idx, parent_idx, child_idx
    type(vector3)    :: pos(this%n), z(this%n), cg(this%n)
    type(quaternion) :: ori(this%n)
    !type(matrix3)    :: R(this%n)
    type(vector6)    :: s(this%n), v(this%n), f(this%n), a(this%n), p(this%n), k(this%n), w(this%n)
    type(matrix6)    :: I(this%n)
    type(vector6)    :: p_A(this%n), T(this%n)
    type(matrix6)    :: I_A(this%n), RU(this%n)
    
    !n = this%n
    
    ! Forward recursion to compute kinematics_t
    do idx=1, this%n
        rb = this%body(idx)
        parent_idx = this%parents(idx)
        if( parent_idx == 0 ) then
            prev_pos = o_
            prev_ori = q_eye
            v_prev = zero_vec6
            a_prev = zero_vec6
        else
            prev_ori = ori(parent_idx)
            prev_pos = pos(parent_idx) + rot(prev_ori, this%body(parent_idx)%next_pos)
            v_prev = v(parent_idx)
            a_prev = a(parent_idx)
        end if
        z(idx) = rot(prev_ori, this%body(idx)%joint_axis)
        select case(this%body(idx)%joint_type)
        case(revolute)
            ori(idx) = prev_ori .o. rot(z(idx), this%q(idx))
            pos(idx) = prev_pos
            s(idx)   = twist(z(idx), pos(idx), 0.0_r8)
        case(prismatic)
            ori(idx) = prev_ori
            pos(idx) = prev_pos + rot(prev_ori, z(idx) * this%q(idx))
            s(idx)   = twist(z(idx))
        end select
        
        v(idx) = v_prev + s(idx) * this%qp(idx)
        k(idx) = v(idx) .xt. (s(idx) * this%qp(idx))
        
        cg(idx) = pos(idx) + rot(ori(idx) .o. this%body(idx)%base_cg_ori, this%body(idx)%base_cg_pos)
        w(idx)  = wrench( rb%mass * this%gee, cg(idx))

        select case(rb%joint_driver)
        case(known_torque)
            this%qpp(idx) = 0
            this%tau(idx) = rb%motor
        case(known_motion)
            this%qpp(idx) = rb%motor
            this%tau(idx) = 0
        end select
        
        I(idx) = spi(rb%mass, rb%inertia, cg(idx), ori(idx))
        p(idx) = v(idx) .xw. (I(idx) * v(idx))
    end do
    
    
    
    end subroutine
    
    end module mod_dynamics_spatial