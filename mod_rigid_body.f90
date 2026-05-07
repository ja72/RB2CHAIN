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
            procedure, private, pass :: rb_write
            generic :: write(formatted) => rb_write
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
        
        subroutine rb_write(rb, unit, iotype, v_list, iostat, iomsg)
        ! Custom write subroutine for rbody type. Called by 'write' statement in MAIN.f90
        class(rbody), intent(in) :: rb
        integer, intent(in) :: unit
        character(*), intent(in) :: iotype
        integer, intent(in) :: v_list(:)
        integer, intent(out) :: iostat
        character(*), intent(inout) :: iomsg
  
            !write (unit, "(A)", iostat=iostat, iomsg=iomsg) "Object says hello"
        write(unit, fmt=100, iostat=iostat, iomsg=iomsg) 'Rigid Body Properties:', new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Overall Length = ', maxval(rb%base_pos), ' m'                         , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Center of Gravity Offset = ', maxval(rb%cg), ' m'                     , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Mass = ', rb%mass, ' kg'                                              , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'I_xx = ', rb%I_xx, ' kg*m^2'                                          , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'I_yy = ', rb%I_yy, ' kg*m^2'                                          , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'I_zz = ', rb%I_zz, ' kg*m^2'                                          , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Joint Type   (1:revolute, 2:prismatic)        = ', (rb%joint_type)    , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Joint Axis   (1:x-axis, 2:y-axis, 3:z-axis)   = ', (rb%joint_axis)    , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Joint Driver (1:known_torque, 2:known_motion) = ', (rb%joint_driver)  , new_line('a')
        write(unit, fmt=101, iostat=iostat, iomsg=iomsg) 'Motor Torque/Motion = ', rb%motor, ' N*m'                             , new_line('a')
        
100     format(a,a)
101     format(1x,a,g0,a,a)
        
        end subroutine        
        
    end module
    