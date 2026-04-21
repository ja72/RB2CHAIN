!  RB2CHAIN.f90 
!
!  FUNCTIONS:
!  RB2CHAIN - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: RB2CHAIN
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
    program Main
    use show_matrix_mod
    use rigid_body_chain
    implicit none
    
    integer, parameter :: n = 6

    ! Variables
    type(rbody) :: rb
    type(state) :: st(n)
    type(rbchain) :: chain
    type(kinematics) :: sol(n)
    real(real64) :: t, h, qpp(n), tau(n)
    real(real64) :: f_err(n)
    integer :: i, steps, fi, sol_method, j
    integer(int64) :: tic, toc, rate
    real(real32) :: time, steps_per_sec
    !call test_linear_system()
    
    ! Set rigid body properties for ALL bodies
    rb%mass= 0.10_dp
    rb%I_xx = 1.6717776d-6            
    rb%I_yy = 8.3413029d-5
    rb%I_zz = 8.3413029d-5
    rb%base_pos = 0.09_dp*i_
    rb%base_rot = ROT_Z(0*rad_per_deg)
    rb%cg = 0.045_dp*i_
    rb%joint_type = revolute
    rb%joint_axis = z_axis
    rb%joint_driver = known_torque    
    rb%motor = 0D0
    
    ! Fill chain with 'n' bodies, based on properies of 'rb' 
    call chain%initialize_chain(n, rb)

    !Set first body joint at the origin
    chain%rb(1)%base_pos = o_
    !chain%rb(6)%joint_driver = known_motion
    !chain%rb(6)%motor = 0D0
    
    ! Set time step (not used) and initialize MMOI matrices
    call chain%prepare_simulation()    

    ! Prepare state object
    t = 0.0D0
    st%q = [(0d0,i=1,n)]   ! Set joint positions
    st%qp = [(0d0,i=1,n)]  ! Set joint velocities
    st(6)%qp = 1d0
    
    call calc_acceleration_art(chain, t, st%q, st%qp, qpp, tau, sol)
    do j=1,n
        f_err(j) = maxval(abs( sol(j)%fnet - sol(j)%facc ))
    end do
    print *, 'Force Balance Error (max abs value) for each body:'
    call show(f_err)    
    stop
    sol_method = ART_METHOD 
    !sol_method = CRB_METHOD
    
    select case (sol_method)
        case (CRB_METHOD)
             print *, 'Using Composite Rigid Body method'
        case (ART_METHOD)
             print *, 'Using Articulated Inertia method'
        case default
             print *, 'Invalid solution method selected'
             stop
        end select
    
    call SYSTEM_CLOCK(tic, rate)
            
    ! Calculate joint accelerations
    steps = 100000
    h = 5d0/steps
    do i=1,steps
        
        call chain%do_step(t, st, h, sol_method)
        
        if(mod(i,steps/20)==0) then
            write(*,*) 'step=',i,' of ',steps, ' t=', t
        end if
    end do
    
    call SYSTEM_CLOCK(toc, rate)
    time = REAL(toc - tic)/REAL(rate)
    steps_per_sec = REAL(steps)/time
    
    print *, 'total time = ', time, ' seconds'
    print *, 'speed = ', steps_per_sec, ' steps/second'
    
    open(newunit=fi, file='results.txt', action='write')
    write(fi,*) t
    write(fi,*) st%q
    write(fi,*) st%qp
    close(fi)
    
    contains
    
        
    end program
    
    


