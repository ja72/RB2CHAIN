!***************************************************************************

    
    
    !DEC$ REAL:8    
    module mod_rigid_body
    use mod_screws
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
        type :: body_t            
            real(real64) :: mass                                ! Mass of the body
            real(real64) :: cg(3) = o_                          ! Center of mass position vector relative to the top of the joint in local coordinates
            real(real64) :: I_xx, I_yy, I_zz, I_xy, I_xz, I_yz  ! Mass 3×3 symmetric moment of inertia components in local coordinates            
            
            real(real64) :: base_pos(3) = o_                    ! Joint base position in previous body coordinates
            real(real64) :: base_rot(3,3) = E3_                 ! Joint base rotation in previous body coordinates
            integer  :: joint_type = revolute                   ! Joint type constant (1=revolute, 2=prismatic)
            integer  :: joint_axis = z_axis                     ! Joint direction flag (1=x-axis, 2=y-axis, 3=z-axis)
            integer  :: joint_driver = known_torque             ! Joint driver flag (1=known torque, 2=known motion)            
            real(real64) :: motor                               ! Value to used for joint torque or motion depending on 'joint_drive' switch
            !real(real64) :: Ic(3,3), Mc(3,3)        ! Local MMOI and inverse MMOI. Calculated once only            
            !logical        :: mmoi_initialized = .false.        ! Flag for local MMOI calculation
        contains
            procedure, pass :: get_spatial_mass_matrix          ! Used by calc_dynamics to get 6×6 spatial MMOI matrix for each frame
            procedure, private, pass :: get_joint_properties    ! Used by calc_kinematics to get joint axis screw for each frame
            procedure, private, pass :: rb_write
            generic :: write(formatted) => rb_write
        end type
                        
        
    contains    
    
        
        pure function get_spatial_mass_matrix(rb, cg, rot, inv) result(res)
        class(body_t), intent(in) :: rb
        real(real64), intent(in) :: rot(3,3), cg(3)
        logical, intent(in), optional :: inv
        real(real64) :: res(6,6)
        real(real64) :: Ic(3,3), Mc(3,3), d
        real(real64) :: mcx(3,3), cgx(3,3)
    
            cgx = cross_operator_matrix(cg)
            
            ! Set columns of local MMOI in Ic
            Ic(:,1) = [rb%I_xx, rb%I_xy, rb%I_xz]
            Ic(:,2) = [rb%I_xy, rb%I_yy, rb%I_yz]
            Ic(:,3) = [rb%I_xz, rb%I_yz, rb%I_zz]
            
            if( present(inv) ) then
                if( inv ) then
                    
                    ! Find discriminate of 3×3 MMOI matrix
                    d = Ic(1,1)*Ic(2,2)*Ic(3,3)+2*Ic(1,2)*Ic(2,3)*Ic(3,1)-Ic(1,1)*Ic(2,3)**2D0-Ic(2,2)*Ic(1,3)**2D0-Ic(3,3)*Ic(1,2)**2D0
                    ! Set columns of inverse of MMOI in Mc
                    Mc(:,1) = [Ic(2,2)*Ic(3,3)-Ic(2,3)**2D0, Ic(1,3)*Ic(2,3)-Ic(1,2)*Ic(3,3), Ic(1,2)*Ic(2,3)-Ic(1,3)*Ic(2,2)]/d
                    Mc(:,2) = [Ic(1,3)*Ic(2,3)-Ic(1,2)*Ic(3,3), Ic(1,1)*Ic(3,3)-Ic(1,3)**2D0, Ic(1,2)*Ic(1,3)-Ic(1,1)*Ic(2,3)]/d
                    Mc(:,3) = [Ic(1,2)*Ic(2,3)-Ic(1,3)*Ic(2,2), Ic(1,2)*Ic(1,3)-Ic(1,1)*Ic(2,3), Ic(1,1)*Ic(2,2)-Ic(1,2)**2D0]/d        
                    
                    Mc = matmul(rot, matmul(Mc, transpose(rot)))
                    mcx = matmul(Mc, cgx)
                    res(1:3, 1:3) = E3_/rb%mass - matmul(cgx, mcx)
                    res(1:3, 4:6) = matmul(cgx, Mc)
                    res(4:6, 1:3) = -mcx
                    res(4:6, 4:6) = Mc
                    return
                end if
            end if
            
            
            Ic = matmul(rot, matmul(Ic, transpose(rot)))
            mcx = rb%mass*cgx
            res(1:3, 1:3) = rb%mass*E3_
            res(1:3, 4:6) = -mcx
            res(4:6, 1:3) =  mcx
            res(4:6, 4:6) = Ic-matmul(mcx,cgx)
            return
        end function
    
        ! Calculate top of joint position and orientation and define joint axis unit twist
        pure subroutine get_joint_properties(rb, prev_pos, prev_rot, q, pos, rot, axis)
        class(body_t), intent(in) :: rb
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
        ! Custom write subroutine for body_t type. Called by 'write' statement in MAIN.f90
        class(body_t), intent(in) :: rb
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
        
        function rb_init_demo() result(rb)
        type(body_t) :: rb
        
            real(real64), parameter :: L = 0.09_real64, c = L/2
    
            rb%mass= 0.10_real64
            rb%I_xx = rb%mass * (0.0454304*L)**2  !1.6717776d-6
            rb%I_yy = rb%mass * (0.3209035*L)**2  !8.3413029d-5
            rb%I_zz = rb%mass * (0.3209035*L)**2  !8.3413029d-5
            rb%base_pos = L*i_
            rb%base_rot = ROT_Z(0*rad_per_deg)
            rb%cg = c*i_
            rb%joint_type = revolute
            rb%joint_axis = z_axis
            rb%joint_driver = known_torque    
            rb%motor = 0.0_real64
        
        end function
        
    end module
    