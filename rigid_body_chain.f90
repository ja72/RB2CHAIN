module rigid_body_chain
    use show_matrix_mod
    use mod_rigid_body
    use mod_kinematics
    implicit none
    
    enum, bind(c)
      enumerator art_method
      enumerator crb_method
    end enum

    type :: rbchain
        integer(int32) :: n_count
        real(real64) :: time, step
        type(rbody), allocatable :: bodies(:)
        integer(int32), allocatable :: parents(:)
    contains
        procedure, pass :: initialize_chain
        procedure, pass :: chain_kinematics
        procedure, pass :: chain_articulated
        procedure, pass :: chain_dynamics
        procedure, pass :: calc_acceleration_art
        procedure, pass :: calc_acceleration_crb
        procedure, pass :: do_step
    end type

    contains
    
    subroutine do_step(chain, n_count, t, q, qp, qpp, tau, h, method)
    !dec$ attributes dllexport :: do_step
    !dec$ attributes alias:'do_step' :: do_step
    !dec$ attributes value :: n_count, h, method
    class(rbchain), intent(inout) :: chain
    REAL(real64), intent(inout) :: t
    integer(int32), intent(in) :: n_count
    real(real64), intent(in) :: h
    integer, intent(in) :: method
    real(real64), intent(inout) :: q(n_count), qp(n_count), qpp(n_count), tau(n_count)
    real(real64) :: K(4,n_count), C(4,n_count)
        
        select case(method)
        case(ART_METHOD)
            C(1,:) = qp
            call calc_acceleration_art(chain, t, q, C(1,:), K(1,:), tau )
            C(2,:) = C(1,:) + h/2*K(1,:)
            call calc_acceleration_art(chain, t+h/2, q+h/2*C(1,:),C(2,:), K(2,:), tau)
            C(3,:) = C(1,:) + h/2*K(2,:)
            call calc_acceleration_art(chain, t+h/2, q+h/2*C(2,:),C(3,:), K(3,:), tau)
            C(4,:) = C(1,:) + h*K(3,:)
            call calc_acceleration_art(chain, t+h, q+h*C(3,:),C(4,:), K(4,:), tau)            
        case(CRB_METHOD)
            C(1,:) = qp
            call calc_acceleration_crb(chain, t, q, C(1,:), K(1,:), tau )
            C(2,:) = C(1,:) + h/2*K(1,:)
            call calc_acceleration_crb(chain, t+h/2, q+h/2*C(1,:),C(2,:), K(2,:), tau)
            C(3,:) = C(1,:) + h/2*K(2,:)
            call calc_acceleration_crb(chain, t+h/2, q+h/2*C(2,:),C(3,:), K(3,:), tau)
            C(4,:) = C(1,:) + h*K(3,:)
            call calc_acceleration_crb(chain, t+h, q+h*C(3,:),C(4,:), K(4,:), tau)            
        case default
        error stop 'Invalid solution method selected'
        end select
        
        t = t + h 
        q = q + h*(C(1,:)+2*C(2,:)+2*C(3,:)+C(4,:))/6
        qp = C(1,:) + h*(K(1,:)+2*K(2,:)+2*K(3,:)+K(4,:))/6
        call calc_acceleration_art(chain, t, q, qp, qpp, tau )
    end subroutine

    subroutine initialize_chain(chain, n_count, rb_prototype)
    class(rbchain), intent(out) :: chain
    integer(int32), intent(in) :: n_count
    type(rbody), intent(in), optional :: rb_prototype
    integer(int32) :: idx
        chain%n_count = n_count
        allocate(chain%bodies(n_count))
        allocate(chain%parents(n_count))
        if( present(rb_prototype) ) then
            chain%bodies = rb_prototype
        end if
        ! Default to a chain structure. Each parent is the previous body
        forall(idx=1:n_count) chain%parents(idx)=idx-1
    end subroutine

    !subroutine prepare_simulation(chain, time_step)
    !class(rbchain), intent(inout) :: chain
    !real(real64), intent(in), optional :: time_step
    !integer(int32) :: i
    !    chain%time = 0
    !    if(present(time_step)) then
    !        chain%step = time_step
    !    else
    !        chain%step = 1d-4
    !    end if
    !    do i=1,chain%n_count
    !        call chain%bodies(i)%initialize_mmoi()
    !    end do
    !end subroutine

    function chain_kinematics(chain, time, q, qp) result(kin)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: time
    real(real64), intent(in) :: q(:), qp(:)
    type(kinematics) :: kin(chain%n_count)
    type(rbody) :: rbs(chain%n_count)
    real(real64) :: z(3)
    integer(int32) :: idx, n_count, jdx
    real(real64) :: prev_pos(3), prev_rot(3,3), prev_vel(6), mcx(3,3)
    real(real64) :: Ic(3,3), Mc(3,3), rot_t(3,3), I(3,3), M(3,3), cgx(3,3)

        n_count  =chain%n_count
        rbs = chain%bodies
        chain%time = time
        do idx=1,n_count
            jdx = chain%parents(idx)
            kin(idx)%parent_index = jdx
            ! Any variable applied forces should be included here
            kin(idx)%applied_force = screw_o_

            if(jdx>0) then
                call kin(idx)%calc(rbs(idx), time, q(idx), qp(idx), kin(jdx))
            else
                call kin(idx)%calc(rbs(idx), time, q(idx), qp(idx))
            end if            
            
            !dec$ IF DEFINED    (DEBUG)
            print *, 'Calculating kinematics for body', idx
            print *, 'pos:'
            call show(kin(idx)%pos)
            print *, 'axis:'
            call show(kin(idx)%axis)
            print *, 'spi:'
            call show(kin(idx)%spi)
            print *, 'weight:'
            call show(kin(idx)%weight)
            print *, 'vel:'
            call show(kin(idx)%vel)
            print *, 'kappa:'
            call show(kin(idx)%kappa)
            print *, 'bias:'
            call show(kin(idx)%bias)
            !dec$ endif
        end do
    end function

    function chain_articulated(chain, kin) result(art)
    class(rbchain), intent(in) :: chain
    type(kinematics), intent(in) :: kin(:)
    type(articulated) :: art(chain%n_count)
    real(real64) :: ars(6)
    integer(int32) :: idx, n_count, kdx
        n_count = chain%n_count
        !allocate(art(n_count))
        do idx=n_count,1,-1
            call art(idx)%init(kin(idx))
            do kdx=idx+1,n_count
                if(kin(kdx)%parent_index==idx) then
                    call art(idx)%next(art(kdx))
                end if
            end do

            call art(idx)%calc()
            
            !dec$ IF DEFINED    (DEBUG)
            print *, 'Calculating articulated for body', idx
            print *, 'ari: (articulated inertia)'
            call show(art(idx)%ari)
            print *, 'arb: (articulated bias force)'
            call show(art(idx)%arb)
            print *, 'iap: (articulated percussiob axis)'
            call show(art(idx)%iap)
            print *, 'rsp: (artuculated reaction space)'
            call show(art(idx)%rsp)
            !dec$ endif
        end do
    end function

    function chain_dynamics(chain, art) result(kin)
    class(rbchain), intent(in) :: chain
    type(articulated),intent(in):: art(:)
    type(kinematics) :: kin(chain%n_count)
    real(real64) :: prev_acc(6), h, hh, s(6), A(6,6), f_nxt(6)
    integer(int32) :: n_count, idx, jdx, kdx
        n_count = chain%n_count
        do idx=1,n_count
            kin(idx) = art(idx)%kin
            jdx = kin(idx)%parent_index
            if(jdx==0) then
                prev_acc = screw_o_
            else
                prev_acc = kin(jdx)%acc
            end if
            ! Joint space inertia
            s = kin(idx)%axis
            A = art(idx)%ari
            h = dot_product(s, matmul(A, s))
            
            select case(kin(idx)%joint_driver)
            case(known_torque)
                kin(idx)%qpp = ( kin(idx)%tau - dot_product(s, matmul(A, prev_acc + kin(idx)%kappa) + art(idx)%arb))/h
            case(known_motion)
                kin(idx)%tau = h*kin(idx)%qpp + dot_product(s, matmul(A, prev_acc+kin(idx)%kappa)+art(idx)%arb)
            end select
            !tex: ${\bf a}_i = {\bf a}_{i-1} + {\bf s}_i \ddot{q}_i + \boldsymbol{\kappa}_i$
            kin(idx)%acc = prev_acc + s*kin(idx)%qpp + kin(idx)%kappa
            !tex: ${\bf f}_i = {\bf I}_i^A {\bf a}_i + {\bf p}_i^A$
            kin(idx)%force = matmul(A, kin(idx)%acc) + art(idx)%arb
            !tex: ${\bf f}_{\rm acc}  = {\bf I}_i {\bf a}_i + {\bf p}_i$
            kin(idx)%facc = matmul(kin(idx)%spi, kin(idx)%acc) + kin(idx)%bias
            
            !dec$ IF DEFINED    (DEBUG)
            print *, 'Calculating dynamics for body', idx
            print *, 'acc:'
            call show(kin(idx)%acc)
            print *, 'frc:'
            call show(kin(idx)%force)
            print *, 'facc:'
            call show(kin(idx)%facc)
            !dec$ ENDIF
        end do
        
        !dec$ IF DEFINED    (DEBUG)
        print *, 'qpp:'
        call show(kin(:)%qpp)
        !dec$ ENDIF
        
        do idx=n_count,1,-1
            f_nxt = screw_o_
            do kdx=idx+1,n_count
                if(kin(kdx)%parent_index==idx) then
                    f_nxt = f_nxt + kin(kdx)%force
                end if
            end do
            kin(idx)%fnet = kin(idx)%force - f_nxt + ( kin(idx)%weight + kin(idx)%applied_force )
        end do
    end function

    subroutine calc_acceleration_art(chain,t,q,qp,qpp,tau,sol)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: t, q(:), qp(:)
    real(real64), intent(inout) :: qpp(chain%n_count), tau(chain%n_count)
    type(kinematics), intent(out), optional :: sol(chain%n_count)
    type(kinematics) :: kin(chain%n_count)
    type(articulated) :: art(chain%n_count)
    integer(int32) :: n_count
    
        n_count = chain%n_count
        kin = chain%chain_kinematics(t, q, qp)
        art = chain%chain_articulated(kin)
        kin = chain%chain_dynamics(art)
        qpp = kin%qpp
        tau = kin%tau
        if( present(sol) ) then
            sol = kin
        end if
    end subroutine

    subroutine calc_acceleration_crb(chain,t,q,qp,qpp,tau)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: t, q(chain%n_count), qp(chain%n_count)
    real(real64), intent(inout) :: qpp(chain%n_count), tau(chain%n_count)
    type(kinematics) :: kin(chain%n_count)
        
    real(real64), allocatable, save :: &
        axis(:,:), spi(:,:), kappa(:), bias(:), vel(:), acc(:), force(:), &
        relative(:,:), tree(:,:), tree_t(:,:), spc(:,:), axis_t(:,:), mmoi(:,:,:), spc2(:,:)
        
    integer(int32), allocatable, save :: pivot(:), itree(:,:), eri(:)
    integer(int32) :: idx, jdx, u, kdx, row, prow, n_count, known
        
    real(real64), allocatable, save :: A(:,:), b(:), qpp_known(:), tau_known(:)
    !real(real64) :: screw_err(6), joint_err(chain%n_count), acm(3), frc(3)

        n_count = chain%n_count
        kin = chain%chain_kinematics(t, q, qp)
        known = count( kin%joint_driver == known_motion )
        
        ! Only allocate once, since this is an expensive operation. The variables
        ! have the SAVE keyword in order to re-use them on subsequent calls.
        if( .not. allocated(axis) ) then
            allocate(vel(6*n_count))
            allocate(acc(6*n_count))
            allocate(force(6*n_count))
            allocate(axis(6*n_count,n_count))
            allocate(axis_t(n_count,6*n_count))
            allocate(mmoi(6,6,n_count))
            allocate(spi(6*n_count,6*n_count))
            allocate(spc(6*n_count,6*n_count))
            allocate(spc2(6*n_count,6*n_count))
            allocate(kappa(6*n_count))
            allocate(bias(6*n_count))
            allocate(relative(6*n_count,6*n_count))
            allocate(tree(6*n_count,6*n_count))
            allocate(tree_t(6*n_count,6*n_count))
            allocate(itree(n_count,n_count))
            allocate(pivot(n_count))
            allocate(A(n_count,n_count))
            allocate(b(n_count))
            allocate(qpp_known(known))
            allocate(tau_known(n_count-known))
            axis = 0
            spi =0
            relative = identity(6*n_count,6*n_count)
            tree = identity(6*n_count,6*n_count)
            pivot = [ (idx, idx=1, n_count) ]
        end if

        kdx = n_count - known
        u = 0
        itree = 0
        do idx=1,n_count
            itree(idx,idx) = 1
            if( kin(idx)%joint_driver == known_motion ) then
                kdx = kdx + 1
                pivot(kdx) = idx
            else
                u = u + 1
                pivot(u) = idx
            end if

            row = 6*(idx-1)+1
            vel(row:row+5) = kin(idx)%vel
            acc(row:row+5) = kin(idx)%acc
            force(row:row+5) = kin(idx)%force
            axis(row:row+5,idx) = kin(idx)%axis
            spi(row:row+5, row:row+5) = kin(idx)%spi
            mmoi(:,:,idx) = kin(idx)%spi
            kappa(row:row+5) = kin(idx)%kappa
            bias(row:row+5) = kin(idx)%bias-(kin(idx)%weight + kin(idx)%applied_force)

            jdx = chain%parents(idx)
            if(jdx>0) then
                prow = 6*(jdx-1)+1
                relative(row:row+5, prow:prow+5) = -E6_
                do while(jdx>0)
                    tree(row:row+5, prow:prow+5) = E6_
                    itree(idx,jdx) = 1
                    jdx = chain%parents(jdx)
                    prow = 6*(jdx-1)+1
                end do
            end if
        end do
        tree_t = transpose(tree)
        axis_t = transpose(axis)
        ! Compose y=A*x+b and solve
        ! This takes the longest to compute.
        spc = matmul(tree_t, matmul(spi, tree))
        !spc2 = CompositeInertia(mmoi, itree, n_count)
        !allocate(eri(6*n_count))
        !eri = maxloc(abs(spc2-spc))
        !if( maxval(abs(spc2-spc))>1e-12 ) then
        !    stop 'composite error'
        !end if
        
        A = matmul(axis_t, matmul(spc, axis))
        b = matmul(axis_t, matmul(spc, kappa) + matmul(tree_t, bias) )
        
        qpp = kin%qpp ! Array of Joint Accelerations
        tau = kin%tau ! Array of Torques
        
        qpp_known = qpp(pivot(n_count-known+1:n_count))
        tau_known = tau(pivot(1:n_count-known))
        
        call solve_linear_system(A, b, qpp, tau, pivot, n_count, known)
        !call solve_linear_system_mkl(A, b, qpp, tau, pivot, n_count, known)
        
        !tau = matmul(A, qpp) + b
        acc = matmul(tree, matmul(axis, qpp) + kappa)
        force = matmul(tree_t, matmul(spi, acc) + bias)
        
    end subroutine
    
    function CompositeInertia(mmoi,itree,n_count) result(spc)
    integer, intent(in) :: n_count
    real(real64), intent(in) :: mmoi(6,6,n_count)
    integer, intent(in) :: itree(n_count,n_count)
    real(real64) :: spc(6*n_count,6*n_count)
    integer :: idx, row,col, jdx, kdx
        spc = 0
        ! Go by columns
        do jdx=1, n_count
            col = 6*(jdx-1)+1
            do idx=1, jdx-1
                row = 6*(idx-1)+1
                spc(row:row+5,col:col+5) = spc(col:col+5,row:row+5)             ! Symmetric 6×6 matrix does not need transpose
            end do
            ! scan rows
            do idx=jdx, n_count
                do kdx=idx,jdx,-1
                    if(itree(kdx,jdx)==1) then
                        row = 6*(kdx-1)+1
                        spc(row:row+5,col:col+5) = spc(row:row+5,col:col+5) + mmoi(:,:,kdx)
                    end if
                end do
            end do
        end do
        
    end function
    

end module