    module mod_rigid_body
    use screw_algebra
    implicit none

        ! Gravity vector definition
        real(real64), parameter :: gee(3) = [0D0,-10D0,0D0]        
    
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
            real(real64) :: mass                                ! Mass of the body
            real(real64) :: cg(3) = o_                       ! Center of mass position vector relative to the top of the joint in local coordinates
            real(real64) :: I_xx, I_yy, I_zz, I_xy, I_xz, I_yz  ! Mass 3×3 symmetric moment of inertia components in local coordinates
            real(real64) :: base_pos(3) = o_                 ! Joint base position in previous body coordinates
            real(real64) :: base_rot(3,3) = E3_           ! Joint base rotation in previous body coordinates
            integer  :: joint_type = revolute                   ! Joint type constant (1=revolute, 2=prismatic)
            integer  :: joint_axis = z_axis                     ! Joint direction flag (1=x-axis, 2=y-axis, 3=z-axis)
            integer  :: joint_driver = known_torque             ! Joint driver flag (1=known torque, 2=known motion)            
            real(real64) :: motor                               ! Value to used for joint torque or motion depending on 'joint_drive' switch
            !real(real64) :: Ic(3,3), Mc(3,3)        ! Local MMOI and inverse MMOI. Calculated once only            
            !logical        :: mmoi_initialized = .false.        ! Flag for local MMOI calculation
        contains
            procedure, pass :: initialize_mmoi                  ! Called once to set 3×3 MMOI matrices from component values
            procedure, private, pass :: get_joint_properties    ! Used by calc_kinematics to get joint axis screw for each frame
        end type
                        
    contains    
    
        pure subroutine initialize_mmoi(rb, Ic, Mc)
        class(rbody), intent(inout) :: rb        
        real(real64), intent(out) :: Ic(3,3), Mc(3,3)
        real(real64) ::  d        
                    
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
            
        end subroutine
    
        ! Calculate top of joint position and orientation and define joint axis unit twist
        pure subroutine get_joint_properties(rb, prev_pos, prev_rot, q, pos, rot, axis)
        class(rbody), intent(in) :: rb
        real(real64), intent(in) :: prev_pos(3), prev_rot(3,3), q             ! previous body position & orientation. q: Joint position
        real(real64), intent(out) :: pos(3), rot(3,3), axis(6)                ! current body position & orientation. Joint axis screw
        real(real64) :: z(3)        
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
                
        
    end module
    