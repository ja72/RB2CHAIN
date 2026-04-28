    module mod_kinematics
    use mod_rigid_body
    implicit none

    type :: kinematics
        integer  :: parent_index                          ! Array index of parent body
        integer  :: joint_driver                          ! Joint drive flag (1=known torque, 2=known motion)
        real(real64) :: q                                     ! Joint position angle/distance
        real(real64) :: qp                                    ! Joint velocity
        real(real64) :: qpp                                   ! Joint accleration
        real(real64) :: tau                                   ! Joint torque/force
        real(real64) :: pos(3)                                ! Position vector for top of joint
        real(real64) :: rot(3,3)                              ! Rotation matrix for body orientation
        real(real64) :: axis(6)                               ! Joint motion axis (twist)
        real(real64) :: cg(3)                                 ! Center of mass position vector
        real(real64) :: spi(6,6), spm(6,6)                    ! Spatial inertia and mobility
        real(real64) :: vel(6)                                ! Spatial velocity of body (twist)
        real(real64) :: applied_force(6)                      ! Applied forces on body (wrench)
        real(real64) :: weight(6)                             ! Weight wrench of body (wrench)
        real(real64) :: kappa(6)                              ! Corriolis acceleration of joint (twist)
        real(real64) :: bias(6)                               ! Centripetal forces of body (wrench)
        real(real64) :: acc(6)                                ! Spatial acceleration of the rigid body
        real(real64) :: force(6)                              ! The joint reaction force on the body(wrench)
        real(real64) :: fnet(6)                               ! The net force on the body(wrench)
        real(real64) :: facc(6)                               ! The inertial force on the body(wrench)
    contains
    procedure, pass :: velocity_at => point_velocity
    procedure, pass :: acceleration_at => point_acceleration
    procedure, pass :: calc => calc_kinematics
    end type

    type :: articulated
        type(kinematics) :: kin                               ! Body kinematics
        real(real64) :: ari(6,6)                              ! Articulated inertia matrix
        real(real64) :: arb(6)                                ! Articulated bias forces
        real(real64) :: iap(6)                                ! Articulated axis of percussion
        real(real64) :: rsp(6,6)                              ! Articulated reaction space
    contains
    procedure, pass :: init => init_articulated
    procedure, pass :: next => next_articulated
    procedure, pass :: calc => calc_articulated
    
    end type

    contains

    pure function point_velocity(kin, r) result(v)
    class(kinematics), intent(in) :: kin
    real(real64), intent(in) :: r(3)
    real(real64) :: v(3)
    v = kin%vel(1:3) + ( kin%vel(4:6) .x. r)
    end function

    pure function point_acceleration(kin, r) result(a)
    class(kinematics), intent(in) :: kin
    real(real64), intent(in) :: r(3)
    real(real64) :: a(3)
    a = kin%acc(1:3) + ( kin%acc(4:6) .x. r) + (kin%vel(4:6) .x. (kin%vel(4:6) .x. r))
    end function

    pure subroutine calc_kinematics(kin, rb, time, q, qp, parent_kin)
    ! Arguments
    class(kinematics), intent(inout) :: kin
    type(rbody), intent(in) :: rb
    real(real64), intent(in) :: time
    real(real64), intent(in) :: q, qp
    type(kinematics), intent(in), optional :: parent_kin
    ! Local variables
    real(real64) :: z(3)
    integer(int32) :: idx, n_count, jdx
    real(real64) :: prev_pos(3), prev_rot(3,3), prev_vel(6), mcx(3,3)
    real(real64) :: Ic(3,3), Mc(3,3), rot_t(3,3), I(3,3), M(3,3), cgx(3,3)

    kin%q            = q
    kin%qp           = qp
    kin%joint_driver = rb%joint_driver
    ! Note: bodies%motor should be a function of time and joint configuration
    !       right now it is just a constant.
    select case(kin%joint_driver)
    case(known_torque)
        kin%qpp = 0
        kin%tau = rb%motor
    case(known_motion)
        kin%qpp = rb%motor
        kin%tau = 0
    end select
    !kin%applied_force = screw_o_
    if( present(parent_kin)) then
        prev_pos = parent_kin%pos
        prev_rot = parent_kin%rot
        prev_vel = parent_kin%vel
    else
        prev_pos = o_
        prev_rot = E3_
        prev_vel = [o_, o_]
    end if
    kin%pos = prev_pos + matmul( prev_rot, rb%base_pos )
    kin%rot = matmul(prev_rot, rb%base_rot)
    z = matmul(kin%rot, direction_vector(rb%joint_axis))
    select case(rb%joint_type)
    case(revolute)
        kin%rot = matmul(kin%rot, rotation_matrix(rb%joint_axis, kin%q))
        kin%axis = twist(z, kin%pos)
    case(prismatic)
        kin%pos = kin%pos + z*kin%q
        kin%axis = twist(z)
        case default
        kin%axis = screw_o_
    end select
    ! Fill in the spatial inertia matrix for this body. This is a 6×6 matrix that combines the mass and inertia
    !call rb%get_spatial_mmoi_matrix(kin%pos, kin%rot, kin%cg, kin%spi, kin%spm)

    call rb%initialize_mmoi(Ic, Mc)

    rot_t = transpose(kin%rot)
    !tex: Moment of inertia in world coordinates:
    !$$ \mathcal{I}_i = {\rm R}_i\, \mathcal{I}_{\rm body} {\rm R}_i^T $$
    I = matmul(kin%rot, matmul(Ic, rot_t))
    M = matmul(kin%rot, matmul(Mc, rot_t))

    !tex: Center of mass position vector in world coordinates:
    !$$ \vec{\rm c}_i = \vec{\rm pos}_i + {\rm R}_i\, \vec{\rm cg}_i $$
    kin%cg = kin%pos + matmul(kin%rot, rb%cg)
    cgx = .x. kin%cg
    !tex: $\vec{\rm mcx} = m_i \vec{c}_i \times$
    mcx = rb%mass*cgx

    !tex: Spatial Matrix :
    !$${\bf I}_{i}=\begin{pmatrix}m_{i} & \mbox{-}m_{i}\vec{c}_{i}\times\\
    !m_{i}\vec{c}_{i}\times & \mathcal{I}_{i}-m_{i}\vec{c}_{i}\times\vec{c}_{i}\times
    !\end{pmatrix}$$
    kin%spi(1:3, 1:3) = rb%mass*E3_
    kin%spi(1:3, 4:6) = -mcx
    kin%spi(4:6, 1:3) =  mcx
    kin%spi(4:6, 4:6) = I-matmul(mcx,cgx)

    !tex: Spatial Mobility:
    !$${\bf M}_{i}=\begin{pmatrix}\tfrac{1}{m_{i}}-\vec{c}_{i}\times\mathcal{I}_{i}^{\mbox{-}1}\vec{c}_{i}\times & \vec{c}_{i}\times\mathcal{I}_{i}^{\mbox{-}1}\\
    !\mbox{-}\mathcal{I}_{i}^{\mbox{-}1}\vec{c}_{i}\times & \mathcal{I}_{i}^{\mbox{-}1}
    !\end{pmatrix}$$
    mcx = matmul(M, cgx)
    kin%spm(1:3, 1:3) = E3_/rb%mass - matmul(cgx, mcx)
    kin%spm(1:3, 4:6) = matmul(cgx,M)
    kin%spm(4:6, 1:3) = -mcx
    kin%spm(4:6, 4:6) = M

    !tex: ${\bf w}_i = \pmatrix{ m_i \vec{g} & \vec{c}_i\times m_i \vec{g} }$
    kin%weight= wrench(rb%mass*gee, kin%cg)
    !tex: ${\bf v}_i = {\bf v}_{i-1} + {\bf s}_i \dot{q}_i$
    kin%vel   = prev_vel + kin%axis*kin%qp
    !tex: $\boldsymbol{\kappa}_i = {\bf v}_i \times {\bf s}_i \dot{q}_i$
    kin%kappa = kin%vel .xt. (kin%axis*kin%qp)
    !tex: ${\bf p}_i = {\bf v}_i \times ({\bf I}_i {\bf v}_i)$
    kin%bias  = kin%vel .xw. matmul(kin%spi, kin%vel)
    kin%force = screw_o_
    kin%acc   = screw_o_
    kin%fnet  = screw_o_
    kin%facc  = screw_o_
    end subroutine

    pure subroutine init_articulated(art, kin)
    class(articulated), intent(inout) :: art
    type(kinematics), intent(in) :: kin
        art%kin = kin
        art%ari = kin%spi
        art%arb = kin%bias - (kin%weight + kin%applied_force)
    end subroutine
    
    pure subroutine next_articulated(art, child_art)
    class(articulated), intent(inout) :: art
    class(articulated), intent(in) :: child_art
    type(kinematics) :: child_kin
    
        child_kin = child_art%kin
        ! For all child objects project their inertia
        ! through the joint using the ARTICULATED INERTIA METHOD
        select case(child_kin%joint_driver)
        case(known_torque)
            ! f[k] = T[k]Q[k] + R[k]U[k]*(A[k]*(κ[k]+a[idx])+d[n])
            art%ari = art%ari + matmul( child_art%rsp, child_art%ari )
            art%arb = art%arb + child_art%iap*child_kin%tau + &
                matmul(child_art%rsp, matmul(child_art%ari, child_kin%kappa)+child_art%arb)
        case(known_motion)
            ! f[k] = A[k]*(a[idx]+s[k]*qpp[k]+κ[k])+d[k]
            art%ari = art%ari + child_art%ari
            art%arb = art%arb + child_art%arb + &
                matmul(child_art%ari, child_kin%axis*child_kin%qpp + child_kin%kappa)
        end select
    
    end subroutine
    
    pure subroutine calc_articulated(art)
    ! Arguments
    class(articulated), intent(inout) :: art   
    ! Local Variables
    real(real64) :: ars(6), axis(6)
        axis = art%kin%axis
        ars = matmul(art%ari, axis)
        art%iap = ars/dot_product(axis, ars)
        art%rsp = E6_ - outer_product(art%iap, axis)    
    end subroutine

    end module mod_kinematics