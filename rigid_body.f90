    module rigid_body
    !use constants
    !use vector_algebra
    use screw_algebra
    implicit none

        ! Gravity vector definition
        real(real64), parameter :: gee(nvec) = [0D0,-10D0,0D0]        
    
        !
        enum, bind(c)
            enumerator:: revolute = 1       ! 1-DOF pin joint
            enumerator:: prismatic = 2      ! 1-DOF slider joint
        end enum
        !
        enum, bind(c)
            enumerator:: known_torque = 1   ! Inverse dynamics for joint, acceleration is found from torque
            enumerator:: known_motion = 2   ! Forward dynamics for joint, torque is found from acceleration
        end enum

        ! Specification of rigid body type. Each body rests on top of a 1DOF joint
        type :: rbody            
            real(real64) :: mass                                  ! Mass of the body
            real(real64) :: cg(nvec) = o_                            ! Center of mass position vector relative to the top of the joint in local coordinates
            real(real64) :: I_xx, I_yy, I_zz, I_xy, I_xz, I_yz    ! Mass 3×3 symmetric moment of inertia components in local coordinates
            real(real64) :: base_pos(nvec) = o_                      ! Joint base position in previous body coordinates
            real(real64) :: base_rot(nvec,nvec) = E3_                   ! Joint base rotation in previous body coordinates
            integer  :: joint_type = revolute                 ! Joint type constant (1=revolute, 2=prismatic)
            integer  :: joint_axis = z_axis                   ! Joint direction flag (1=x-axis, 2=y-axis, 3=z-axis)
            integer  :: joint_driver = known_torque           ! Joint driver flag (1=known torque, 2=known motion)            
            real(real64) :: motor                                 ! Value to used for joint torque or motion depending on 'joint_drive' switch
            real(real64) :: Ic(nvec,nvec), Mc(nvec,nvec)                      ! Local MMOI and inverse MMOI. Calculated once only            
            logical        :: mmoi_initialized = .false.            ! Flag for local MMOI calculation
        contains
            procedure, pass :: initialize_mmoi             ! Called once to set 3×3 MMOI matrices from component values
            procedure, pass :: get_spatial_mmoi_matrix     ! Used by calc_kinematics to get 6×6 spatial inertia for each frame
            procedure, private, pass :: get_joint_properties        ! Used by calc_kinematics to get joint axis screw for each frame            
        end type
        
        ! State vector for each joint
        type :: state
        sequence            
            real(real64) :: q             ! Joint position angle/distance
            real(real64) :: qp            ! Joint velocity
        end type
                
        !
        type :: kinematics
            integer  :: parent_index                          ! Array index of parent body
            integer  :: joint_driver                          ! Joint drive flag (1=known torque, 2=known motion)
            real(real64) :: q                                     ! Joint position angle/distance
            real(real64) :: qp                                    ! Joint velocity
            real(real64) :: qpp                                   ! Joint accleration
            real(real64) :: tau                                   ! Joint torque/force
            real(real64) :: pos(nvec)                                ! Position vector for top of joint
            real(real64) :: rot(nvec,nvec)                              ! Rotation matrix for body orientation
            real(real64) :: axis(nspc)                               ! Joint motion axis (twist)
            real(real64) :: cg(nvec)                                 ! Center of mass position vector
            real(real64) :: spi(nspc,nspc), spm(nspc,nspc)                    ! Spatial inertia and mobility
            real(real64) :: vel(nspc)                                ! Spatial velocity of body (twist)
            real(real64) :: applied_force(nspc)                      ! Applied forces on body (wrench)
            real(real64) :: weight(nspc)                             ! Weight wrench of body (wrench)
            real(real64) :: kappa(nspc)                              ! Corriolis acceleration of joint (twist)
            real(real64) :: bias(nspc)                               ! Centripetal forces of body (wrench)
            real(real64) :: acc(nspc)                                ! Spatial acceleration of the rigid body
            real(real64) :: force(nspc)                              ! The joint reaction force on the body(wrench)
            real(real64) :: fnet(nspc)                               ! The net force on the body(wrench)
            real(real64) :: facc(nspc)                               ! The inertial force on the body(wrench)
        contains
            procedure, pass :: velocity_at => point_velocity
            procedure, pass :: acceleration_at => point_acceleration
        end type
        
        type :: articulated
            type(kinematics) :: kin                                 ! Body kinematics
            real(real64) :: ari(nspc,nspc)                              ! Articulated inertia matrix
            real(real64) :: arb(nspc)                                ! Articulated bias forces
            real(real64) :: iap(nspc)                                ! Articulated axis of percussion
            real(real64) :: rsp(nspc,nspc)                              ! Articulated reaction space        
        end type
        
    contains    
    
        pure subroutine initialize_mmoi(rb)
        class(rbody), intent(inout) :: rb        
        real(real64):: Ic(nvec,nvec), Mc(nvec,nvec), d
                    
            ! Set columns of local MMOI in Ic
            Ic(:,1) = [rb%I_xx, rb%I_xy, rb%I_xz]
            Ic(:,2) = [rb%I_xy, rb%I_yy, rb%I_yz]
            Ic(:,3) = [rb%I_xz, rb%I_yz, rb%I_zz]
        
            ! Find discriminate of 3×3 MMOI matrix
            d = Ic(1,1)*Ic(2,2)*Ic(3,3)+2*Ic(1,2)*Ic(2,3)*Ic(3,1)-Ic(1,1)*Ic(2,3)**2D0-Ic(2,2)*Ic(1,3)**2D0-Ic(3,3)*Ic(1,2)**2D0
            ! Set columns of inverse of MMOI in Mc
            Mc(:,1) = [Ic(2,2)*Ic(3,3)-Ic(2,3)**2D0, Ic(1,3)*Ic(2,3)-Ic(1,2)*Ic(3,3), Ic(1,2)*Ic(2,3)-Ic(1,3)*Ic(2,2)]/d
            Mc(:,2) = [Ic(1,3)*Ic(2,3)-Ic(1,2)*Ic(3,3), Ic(1,1)*Ic(3,3)-Ic(1,3)**2D0, Ic(1,2)*Ic(1,3)-Ic(1,1)*Ic(2,3)]/d
            Mc(:,3) = [Ic(1,2)*Ic(2,3)-Ic(1,3)*Ic(2,2), Ic(1,2)*Ic(1,3)-Ic(1,1)*Ic(2,3), Ic(1,1)*Ic(2,2)-Ic(1,2)**2D0]/d        
            
            rb%Ic = Ic
            rb%Mc = Mc
            rb%mmoi_initialized = .true.
        end subroutine
    
        ! Calculate top of joint position and orientation and define joint axis unit twist
        pure subroutine get_joint_properties(rb, prev_pos, prev_rot, q, pos, rot, axis)
        class(rbody), intent(in) :: rb
        real(real64), intent(in) :: prev_pos(nvec), prev_rot(nvec,nvec), q             ! previous body position & orientation. q: Joint position
        real(real64), intent(out) :: pos(nvec), rot(nvec,nvec), axis(nspc)                ! current body position & orientation. Joint axis screw
        real(real64) :: z(nvec)        
            pos = prev_pos + matmul(prev_rot, rb%base_pos)
            rot = matmul(prev_rot, rb%base_rot)        
            z = matmul(rot, direction_vector(rb%joint_axis))     
            select case(rb%joint_type)
            case(revolute)
                rot = matmul(rot, rotation_matrix(rb%joint_axis, q))
                axis = twist(z, pos)                
            case(prismatic)
                pos = pos + z*q
                axis = twist(z)
            case default                                
                axis = screw_o_
            end select
        end subroutine
        
        ! Calculate body inertial properties for current position and orientation
        pure subroutine get_spatial_mmoi_matrix(rb, pos, rot, cg, spi, spm)
        class(rbody), intent(in) :: rb
        real(real64), intent(in) :: pos(nvec), rot(nvec,nvec)
        real(real64), intent(out) :: cg(nvec), spi(nspc,nspc), spm(nspc,nspc)
        real(real64) :: rot_t(nvec,nvec), I(nvec,nvec), M(nvec,nvec), cgx(nvec,nvec), mc(nvec,nvec)
                        
            rot_t = transpose(rot)
            
            I = matmul(rot, matmul(rb%Ic, rot_t))
            M = matmul(rot, matmul(rb%Mc, rot_t))
            
            cg = pos + matmul(rot,rb%cg)
            cgx = .x. cg
            
            mc = rb%mass*cgx
            spi(1:3, 1:3) = rb%mass*E3_
            spi(1:3, 4:6) = -mc
            spi(4:6, 1:3) =  mc
            spi(4:6, 4:6) = I-matmul(mc,cgx)
            
            mc = matmul(M, cgx)
            spm(1:3, 1:3) = E3_/rb%mass - matmul(cgx, mc)
            spm(1:3, 4:6) = matmul(cgx,M)
            spm(4:6, 1:3) = -mc
            spm(4:6, 4:6) = M
        end subroutine       
        
        pure function point_velocity(kin, r) result(v)
        class(kinematics), intent(in) :: kin
        real(real64), intent(in) :: r(nvec)
        real(real64) :: v(nvec)
            v = kin%vel(1:3) + ( kin%vel(4:6) .x. r)
        end function
                                                               
        pure function point_acceleration(kin, r) result(a)
        class(kinematics), intent(in) :: kin
        real(real64), intent(in) :: r(nvec)
        real(real64) :: a(nvec)
            a = kin%acc(1:3) + ( kin%acc(4:6) .x. r) + (kin%vel(4:6) .x. (kin%vel(4:6) .x. r))
        end function
        
    end module
    