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
    use mod_show_matrix
    use mod_body_chain
    implicit none
    
    integer, parameter :: n = 6
    real(real64), parameter :: L = 0.09_r8, c = L/2

    ! Variables
    type(body_t) :: rb
    !type(state) :: st(n)
    type(chain_t) :: chain
    type(kinematics_t) :: sol(n)
    type(state_t(n)) :: state
    real(real64) :: h
    real(real64) :: f_err(n)
    integer :: i, steps, fi, sol_method, j
    integer(int64) :: tic, toc, rate
    real(real32) :: time, steps_per_sec
    
    
    ! chain = rb_initialization(n,rb,t,q,qp,qpp,tau)
    
    rb = rb_init_demo()
    chain  = chain_init_demo(n, rb, state)
    
    ! Show rigid body properties
    print *, rb
    
    !dec$ IF DEFINED    (DEBUG)
    
    call rb_test_calculation(chain, state)
    
    !dec$ ELSE
    
    sol_method = ART_METHOD 
    call rb_do_simulation(chain, state, sol_method)
    
    !dec$ ENDIF
    
    contains
        
    subroutine rb_test_calculation(chain, state)
    use mod_show_matrix    
    type(chain_t), intent(inout) :: chain
    type(state_t(chain%n_count)), intent(inout) :: state
    type(kinematics_t) :: sol(chain%n_count)
    
    call calc_acceleration_art(chain, &
        state%t, &
        state%q, &
        state%qp, &
        state%qpp, &
        state%tau, &
        sol)
    do j=1,n
        f_err(j) = maxval(abs( sol(j)%fnet - sol(j)%facc ))
    end do
    print *, 'Force Balance Error (max abs value) for each body:'
    call show(f_err)    
    
    end subroutine
    
    subroutine rb_do_simulation(chain, state, sol_method)
    use mod_show_matrix
    type(chain_t), intent(inout) :: chain
    type(state_t(chain%n_count)), intent(inout) :: state
    integer, intent(in) :: sol_method
    
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
        
        call chain%do_step(state, h, sol_method)
        
        if(mod(i,steps/20)==0) then
            write(*,*) 'step=',i,' of ',steps, ' t=', state%t
        end if
    end do
    
    call SYSTEM_CLOCK(toc, rate)
    time = REAL(toc - tic)/REAL(rate)
    steps_per_sec = REAL(steps)/time
    
    print *, 'total time = ', time, ' seconds'
    print *, 'speed = ', steps_per_sec, ' steps/second'
    
    open(newunit=fi, file='results.txt', action='write')
    write(fi,*) state%t
    write(fi,*) state%q
    write(fi,*) state%qp
    write(fi,*) state%qpp
    write(fi,*) state%tau
    close(fi)
    end subroutine
    
        
    end program
    
    


