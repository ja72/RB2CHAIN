module rigid_body_chain
    use rigid_body
    implicit none
    
    ENUM, BIND(C)
      ENUMERATOR ART_METHOD
      ENUMERATOR CRB_METHOD
    END ENUM

    type :: rbchain
        integer(int32) :: n
        real(real64) :: time, step
        type(rbody), allocatable :: rb(:)
        integer(int32), allocatable :: parents(:)
    contains
        procedure, pass :: initialize_chain
        procedure, pass :: prepare_simulation
        procedure, pass :: calc_kinematics
        procedure, pass :: calc_articulated
        procedure, pass :: calc_dynamics
        procedure, pass :: calc_acceleration_art
        procedure, pass :: calc_acceleration_crb
        procedure, pass :: do_step
    end type

    contains
    
    subroutine do_step(chain, t, st, h, method)
    class(rbchain), intent(inout) :: chain
    REAL(real64), intent(inout) :: t
    type(state), intent(inout) :: st(:)
    real(real64), intent(in) :: h
    integer, intent(in), optional :: method
    real(real64) :: q(size(st))
    real(real64) :: K(4,size(st)), C(4,size(st))
    integer :: m
        if( present(method) ) then
            m = method
        else
            m = CRB_METHOD
        end if
        
        q = st%q
        select case(m)
        case(ART_METHOD)
            C(1,:) = st%qp
            K(1,:) = calc_acceleration_art(chain, t, q,C(1,:))
            C(2,:) = C(1,:) + h/2*K(1,:)
            K(2,:) = calc_acceleration_art(chain, t+h/2, q+h/2*C(1,:),C(1,:)+h/2*K(1,:))
            C(3,:) = C(1,:) + h/2*K(2,:)
            K(3,:) = calc_acceleration_art(chain, t+h/2, q+h/2*C(2,:),C(1,:)+h/2*K(2,:))
            C(4,:) = C(1,:) + h*K(3,:)
            K(4,:) = calc_acceleration_art(chain, t+h, q+h*C(3,:),C(1,:)+h*K(3,:))            
        case(CRB_METHOD)
            C(1,:) = st%qp
            K(1,:) = calc_acceleration_crb(chain, t, q,C(1,:))
            C(2,:) = C(1,:) + h/2*K(1,:)
            K(2,:) = calc_acceleration_crb(chain, t+h/2, q+h/2*C(1,:),C(1,:)+h/2*K(1,:))
            C(3,:) = C(1,:) + h/2*K(2,:)
            K(3,:) = calc_acceleration_crb(chain, t+h/2, q+h/2*C(2,:),C(1,:)+h/2*K(2,:))
            C(4,:) = C(1,:) + h*K(3,:)
            K(4,:) = calc_acceleration_crb(chain, t+h, q+h*C(3,:),C(1,:)+h*K(3,:))            
        case default
        end select
        
        t = t + h 
        st%q = q + h*(C(1,:)+2*C(2,:)+2*C(3,:)+C(4,:))/6
        st%qp = C(1,:) + h*(K(1,:)+2*K(2,:)+2*K(3,:)+K(4,:))/6
    end subroutine

    subroutine initialize_chain(chain, n, rb_prototype)
    class(rbchain), intent(out) :: chain
    integer(int32), intent(in) :: n
    type(rbody), intent(in), optional :: rb_prototype
    integer(int32) :: i
        chain%n = n
        allocate(chain%rb(n))
        allocate(chain%parents(n))
        if( present(rb_prototype) ) then
            chain%rb = rb_prototype
        end if
        ! Default to a chain structure. Each parent is the previous body
        forall(i=1:n) chain%parents(i)=i-1
    end subroutine

    subroutine prepare_simulation(chain, time_step)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: time_step
    integer(int32) :: i
        chain%time = 0
        chain%step = time_step
        do i=1,chain%n
            call chain%rb(i)%initialize_mmoi()
        end do
    end subroutine

    function calc_kinematics(chain, time, st) result(kin)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: time
    type(state), intent(in) :: st(:)
    type(kinematics) :: kin(chain%n)
    type(rbody) :: rbs(chain%n)
    real(real64) :: z(3)
    integer(int32) :: i, n, j
    real(real64) :: prev_pos(3), prev_rot(3,3), prev_vel(6)

        n  =chain%n
        rbs = chain%rb
        do i=1,n
            kin(i)%joint_driver = rbs(i)%joint_driver
            kin(i)%q = st(i)%q
            kin(i)%qp = st(i)%qp
            ! Note: rb%motor should be a function of time and joint configuration
            !       right now it is just a constant.
            !       Also, any variable applied forces should be included in kin%applied_force
            select case(kin(i)%joint_driver)
            case(known_torque)
                kin(i)%qpp = 0
                kin(i)%tau = rbs(i)%motor
            case(known_motion)
                kin(i)%qpp = rbs(i)%motor
                kin(i)%tau = 0
            end select
            kin(i)%applied_force = wrench_o_
            j = chain%parents(i)
            kin(i)%parent_index = j
            if(j==0) then
                prev_pos = o_
                prev_rot = E3_
                prev_vel = [o_, o_]
            else
                prev_pos = kin(j)%pos
                prev_rot = kin(j)%rot
                prev_vel = kin(j)%vel
            end if
            kin(i)%pos = prev_pos + matmul( prev_rot, rbs(i)%base_pos )
            kin(i)%rot = matmul(prev_rot, rbs(i)%base_rot)
            z = matmul(kin(i)%rot, direction_vector(rbs(i)%joint_axis))
            select case(rbs(i)%joint_type)
            case(revolute)
                kin(i)%rot = matmul(kin(i)%rot, rotation_matrix(rbs(i)%joint_axis, kin(i)%q))
                kin(i)%axis = twist(z, kin(i)%pos)
            case(prismatic)
                kin(i)%pos = kin(i)%pos + z*kin(i)%q
                kin(i)%axis = twist(z)
                case default
                kin(i)%axis = twist_o_
            end select
            call rbs(i)%get_spatial_mmoi_matrix(kin(i)%pos, kin(i)%rot, kin(i)%cg, kin(i)%spi, kin(i)%spm)
            kin(i)%weight = wrench(rbs(i)%mass*gee, kin(i)%cg)
            kin(i)%vel = prev_vel + kin(i)%axis*kin(i)%qp
            kin(i)%kappa = kin(i)%vel .xt. (kin(i)%axis*kin(i)%qp)
            kin(i)%bias = kin(i)%vel .xw. matmul(kin(i)%spi, kin(i)%vel)
            kin(i)%force = wrench_o_
            kin(i)%acc = twist_o_
        end do
    end function

    function calc_articulated(chain, kin) result(art)
    class(rbchain), intent(in) :: chain
    type(kinematics), intent(in) :: kin(:)
    type(articulated) :: art(chain%n)
    real(real64) :: ars(6)
    integer(int32) :: i, n, k
        n = chain%n
        !allocate(art(n))
        do i=n,1,-1
            art(i)%kin = kin(i)
            art(i)%ari = art(i)%kin%spi
            art(i)%arb = kin(i)%bias - kin(i)%total_force()
            ! Check all subsequent bodies if [i] is their parent
            do k=i+1,n
                if(kin(k)%parent_index==i) then
                    ! For all child objects project their inertia
                    ! through the joint using the ARTICULATED INERTIA METHOD
                    select case(kin(k)%joint_driver)
                    case(known_torque)
                        ! f[k] = T[k]Q[k] + R[k]U[k]*(A[k]*(κ[k]+a[i])+d[n])
                        art(i)%ari = art(i)%ari + matmul( art(k)%rsp, art(k)%ari )
                        art(i)%arb = art(i)%arb + art(k)%iap*kin(k)%tau + &
                            matmul(art(k)%rsp, matmul(art(k)%ari, kin(k)%kappa)+art(k)%arb)
                    case(known_motion)
                        ! f[k] = A[k]*(a[i]+s[k]*qpp[k]+κ[k])+d[k]
                        art(i)%ari = art(i)%ari + art(k)%ari
                        art(i)%arb = art(i)%arb + art(k)%arb + &
                            matmul(art(k)%ari, kin(k)%axis*kin(k)%qpp + kin(k)%kappa)
                    end select
                end if
            end do

            ars = matmul(art(i)%ari, kin(i)%axis)
            art(i)%iap = ars/dot_product(kin(i)%axis, ars)
            art(i)%rsp = E6_ - outer_product(art(i)%iap, kin(i)%axis)
        end do
    end function

    function calc_dynamics(chain, art) result(kin)
    class(rbchain), intent(in) :: chain
    type(articulated),intent(in):: art(:)
    type(kinematics) :: kin(chain%n)
    real(real64) :: prev_acc(6), h, hh, s(6), A(6,6)
    integer(int32) :: n, i, j
        n = chain%n
        do i=1,n
            kin(i) = art(i)%kin
            j = kin(i)%parent_index
            if(j==0) then
                prev_acc = twist_o_
            else
                prev_acc = kin(j)%acc
            end if
            ! Joint space inertia
            s = kin(i)%axis
            A = art(i)%ari
            h = dot_product(s, matmul(A, s))
            
            select case(kin(i)%joint_driver)
            case(known_torque)
                kin(i)%qpp = ( kin(i)%tau - dot_product(s, matmul(A, prev_acc + kin(i)%kappa) + art(i)%arb))/h
            case(known_motion)
                kin(i)%tau =h*kin(i)%qpp + dot_product(s, matmul(A, prev_acc+kin(i)%kappa)+art(i)%arb)
            end select
            kin(i)%acc = prev_acc + s*kin(i)%qpp + kin(i)%kappa
            kin(i)%force = matmul(A, kin(i)%acc) + art(i)%arb
        end do

    end function

    function calc_acceleration_art(chain,t,q,qp) result(qpp)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: t, q(chain%n), qp(chain%n)
    real(real64) :: qpp(chain%n), tau(chain%n)
    type(state) :: st(chain%n)
    type(kinematics) :: kin(chain%n)
    type(articulated) :: art(chain%n)
    integer(int32) :: n
    
        n = chain%n
        st%q = q
        st%qp = qp
        kin = chain%calc_kinematics(t, st)
        art = chain%calc_articulated(kin)
        kin = chain%calc_dynamics(art)
        qpp = kin%qpp
        tau = kin%tau
    end function

    function calc_acceleration_crb(chain,t,q,qp) result(qpp)
    class(rbchain), intent(inout) :: chain
    real(real64), intent(in) :: t, q(chain%n), qp(chain%n)
    real(real64) :: qpp(chain%n), tau(chain%n)        
    type(state) :: st(chain%n)
    type(kinematics) :: kin(chain%n)
        
    real(real64), allocatable, save :: &
        axis(:,:), spi(:,:), kappa(:), bias(:), vel(:), acc(:), force(:), &
        relative(:,:), tree(:,:), tree_t(:,:), spc(:,:), axis_t(:,:), mmoi(:,:,:), spc2(:,:)
        
    integer(int32), allocatable, save :: pivot(:), itree(:,:), eri(:)
    integer(int32) :: i, j, u, k, row, prow, n, known
        
    real(real64), allocatable, save :: A(:,:), b(:), qpp_known(:), tau_known(:)
    !real(real64) :: screw_err(6), joint_err(chain%n), acm(3), frc(3)

        n = chain%n
        st%q = q
        st%qp = qp
        kin = chain%calc_kinematics(t, st)
        known = count( kin%joint_driver == known_motion )
        
        ! Only allocate once, since this is an expensive operation. The variables
        ! have the SAVE keyword in order to re-use them on subsequent calls.
        if( .not. allocated(axis) ) then
            allocate(vel(6*n))
            allocate(acc(6*n))
            allocate(force(6*n))
            allocate(axis(6*n,n))
            allocate(axis_t(n,6*n))
            allocate(mmoi(6,6,n))
            allocate(spi(6*n,6*n))
            allocate(spc(6*n,6*n))
            allocate(spc2(6*n,6*n))
            allocate(kappa(6*n))
            allocate(bias(6*n))
            allocate(relative(6*n,6*n))
            allocate(tree(6*n,6*n))
            allocate(tree_t(6*n,6*n))
            allocate(itree(n,n))
            allocate(pivot(n))
            allocate(A(n,n))
            allocate(b(n))
            allocate(qpp_known(known))
            allocate(tau_known(n-known))
            axis = 0
            spi =0
            relative = identity(6*n,6*n)
            tree = identity(6*n,6*n)
            pivot = [ (i, i=1, n) ]
        end if

        k = n - known
        u = 0
        itree = 0
        do i=1,n
            itree(i,i) = 1
            if( kin(i)%joint_driver == known_motion ) then
                k = k + 1
                pivot(k) = i
            else
                u = u + 1
                pivot(u) = i
            end if

            row = 6*(i-1)+1
            vel(row:row+5) = kin(i)%vel
            acc(row:row+5) = kin(i)%acc
            force(row:row+5) = kin(i)%force
            axis(row:row+5,i) = kin(i)%axis
            spi(row:row+5, row:row+5) = kin(i)%spi
            mmoi(:,:,i) = kin(i)%spi
            kappa(row:row+5) = kin(i)%kappa
            bias(row:row+5) = kin(i)%bias-kin(i)%total_force()

            j = chain%parents(i)
            if(j>0) then
                prow = 6*(j-1)+1
                relative(row:row+5, prow:prow+5) = -E6_
                do while(j>0)
                    tree(row:row+5, prow:prow+5) = E6_
                    itree(i,j) = 1
                    j = chain%parents(j)
                    prow = 6*(j-1)+1
                end do
            end if
        end do
        tree_t = transpose(tree)
        axis_t = transpose(axis)
        ! Compose y=A*x+b and solve
        ! This takes the longest to compute.
        spc = matmul(tree_t, matmul(spi, tree))
        !spc2 = CompositeInertia(mmoi, itree, n)
        !allocate(eri(6*n))
        !eri = maxloc(abs(spc2-spc))
        !if( maxval(abs(spc2-spc))>1e-12 ) then
        !    stop 'composite error'
        !end if
        
        A = matmul(axis_t, matmul(spc, axis))
        b = matmul(axis_t, matmul(spc, kappa) + matmul(tree_t, bias) )
        
        qpp = kin%qpp ! Array of Joint Accelerations
        tau = kin%tau ! Array of Torques
        
        qpp_known = qpp(pivot(n-known+1:n))
        tau_known = tau(pivot(1:n-known))
        
        call solve_linear_system(A, b, qpp, tau, pivot, n, known)
        !call solve_linear_system_mkl(A, b, qpp, tau, pivot, n, known)
        
        !tau = matmul(A, qpp) + b
        acc = matmul(tree, matmul(axis, qpp) + kappa)
        force = matmul(tree_t, matmul(spi, acc) + bias)
        
    end function
    
    function CompositeInertia(mmoi,itree,n) result(spc)
    integer, intent(in) :: n
    real(real64), intent(in) :: mmoi(6,6,n)
    integer, intent(in) :: itree(n,n)
    real(real64) :: spc(6*n,6*n)
    integer :: i, row,col, j, k
        spc = 0
        ! Go by columns
        do j=1, n
            col = 6*(j-1)+1
            do i=1, j-1
                row = 6*(i-1)+1
                spc(row:row+5,col:col+5) = spc(col:col+5,row:row+5)             ! Symmetric 6×6 matrix does not need transpose
            end do
            ! scan rows
            do i=j, n
                do k=i,j,-1
                    if(itree(k,j)==1) then
                        row = 6*(k-1)+1
                        spc(row:row+5,col:col+5) = spc(row:row+5,col:col+5) + mmoi(:,:,k)
                    end if
                end do
            end do
        end do
        
    end function
    

end module