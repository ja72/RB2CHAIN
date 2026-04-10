!***************************************************************************
!
! FILE:         : SCREWTHEORY.F90
!
! PURPOSE:      : MODULE DEFININING SCREW TYPES WITH OPERATORS.
!
! DEPENDENCIES  : NONE
!
! HISTORY       : Written by John Alexiou - Mar, 2017
!
!  Modified by   Date   Change(s) made
!  -----------  ------  ------------------------------------------------
!
    MODULE screw_algebra
    !USE constants
    USE vector_algebra
    implicit none
    
        integer, parameter :: nspc = 2*nvec     ! Number of components in a screw

        ! Unit screws definition
        real(real64), parameter ::  twist_o_(nspc) = [o_, o_]
        real(real64), parameter ::  twist_i_(nspc) = [o_, i_]
        real(real64), parameter ::  twist_j_(nspc) = [o_, j_]
        real(real64), parameter ::  twist_k_(nspc) = [o_, k_]
        real(real64), parameter :: wrench_o_(nspc) = [o_, o_]
        real(real64), parameter :: wrench_i_(nspc) = [i_, o_]
        real(real64), parameter :: wrench_j_(nspc) = [j_, o_]
        real(real64), parameter :: wrench_k_(nspc) = [k_, o_]
        real(real64), parameter :: Z6_(nspc,nspc) = reshape( &
            [0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0], &
            [nspc,nspc] )
        real(real64), parameter :: E6_(nspc,nspc) = reshape( &
            [1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0, 0,0,0,0,1,0, 0,0,0,0,0,1], &
            [nspc,nspc] )
    
        ! Generic names for functions on screws
        interface twist
            module procedure twist_at, pure_twist
        end interface
        interface wrench
            module procedure wrench_at, pure_wrench
        end interface
        interface operator (.xt.)
            module procedure twist_cross_twist
            module procedure twist_cross_twist_matrix
        end interface
        interface operator (.xw.)
            module procedure twist_cross_wrench
            module procedure twist_cross_wrench_matrix
        end interface

    CONTAINS
    
        ! Defines a twist along vector z located at vector r with optional pitch h
        ! twist = (r×z+h*z, z)
        pure function twist_at(z, r, h) result(s)
            real(real64), intent(in) :: z(nvec), r(nvec)
            real(real64), intent(in), optional :: h
            real(real64) :: s(nspc)            
                if(present(h) ) then
                    s = [(r .x. z) + h*z, z]
                else
                    s = [r .x. z, z]
                end if
        end function
        
        ! Defines a wrench along vector z located at vector r with optional pitch h
        ! wrench = (z, r×z+h*z)
        pure function wrench_at(z, r, h) result(s)
            real(real64), intent(in) :: z(nvec), r(nvec)
            real(real64), intent(in), optional :: h
            real(real64) :: s(nspc)            
                if(present(h) ) then
                    s = [z, (r .x. z) + h*z]
                else
                    s = [z, r .x. z]
                end if
        end function
        
        ! Defines a pure twist along vector z
        ! twist = (z, 0)
        pure function pure_twist(z) result(s)
            real(real64), intent(in) :: z(nvec)
            real(real64) :: s(nspc)            
                s = [z, o_]
        end function
        
        ! Defines a pure wrench along vector z
        ! wrench = (0, z)
        pure function pure_wrench(z) result(s)
            real(real64), intent(in) :: z(nvec)
            real(real64) :: s(nspc)            
                s = [o_, z]
        end function
        
        ! Calculates the scalar magnitude of a twist        
        pure function twist_magnitude(s) result(m)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: m
            m =  magnitude(s(4:6))
        end function
        
        ! Calculates the unit direction of a twist
        pure function twist_direction(s) result(z)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: z(nvec)
            z =  direction(s(4:6))            
        end function
        
        ! Calculates the closest position of the twist to the origin
        pure function twist_position(s) result(r)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: r(nvec), z(nvec), m(nvec)
            z = s(4:6)
            m = s(1:3)
            r = (z .x. m)/sumsq(z)
        end function
        
        ! Calculates the scalar pitch of the twist
        pure function twist_pitch(s) result(h)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: h, z(nvec), m(nvec)
            z = s(4:6)
            m = s(1:3)
            h = (z .dot. m)/sumsq(z)            
        end function
        
        ! Calculates the scalar magnitude of a wrench
        pure function wrench_magnitude(s) result(m)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: m
            m =  magnitude(s(1:3))
        end function
        
        ! Calculates the unit direction of a wrench
        pure function wrench_direction(s) result(z)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: z(nvec)
            z =  direction(s(1:3))            
        end function
        
        ! Calculates the closest position of the wrench to the origin
        pure function wrench_position(s) result(r)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: r(nvec), z(nvec), m(nvec)
            z = s(1:3)
            m = s(4:6)
            r = (z .x. m)/sumsq(z)
        end function
        
        ! Calculates the scalar pitch of the wrench
        pure function wrench_pitch(s) result(h)
            real(real64), intent(in) :: s(nspc)
            real(real64) :: h, z(nvec), m(nvec)
            z = s(1:3)
            m = s(4:6)
            h = (z .dot. m)/sumsq(z)            
        end function
        
        ! Calculates the dot product of two screws
        pure function inner_product_screw(a, b) result(s)
            real(real64), intent(in) :: a(nspc), b(nspc)
            real(real64) :: s
            
            s = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*b(4)+a(5)*b(5)+a(6)*b(6)
            
        end function
        
        ! Calculates the outer product of two screws
        pure function outer_product_screw(a,b) result(c)        
        !argument specification
            real(real64), intent(in) :: a(nspc), b(nspc)
            real(real64) :: c(nspc,nspc)
        !code start
            c(:,1) = a(:)*b(1)
            c(:,2) = a(:)*b(2)
            c(:,3) = a(:)*b(3)
            c(:,4) = a(:)*b(4)
            c(:,5) = a(:)*b(5)
            c(:,6) = a(:)*b(6)
            
            ! c = spread(a,dim=2, ncopies=6)*spread(b,dim=1,ncopies=6)
        end function
        
        ! Returns the (twist × twist) = twist cross product
        ! result =(v, ω)×(a, b) = (ω×a + v×b, ω×b)
        pure function twist_cross_twist(t, s) result(c)
        !argument specification
            real(real64), intent(in) :: t(nspc), s(nspc)
            real(real64) :: c(nspc)
        !code start
            c = [ (t(4:6) .x. s(1:3)) + (t(1:3) .x. s(4:6)), t(4:6) .x. s(4:6) ]
        end function     
        
        ! Returns the (twist × twist) = twist cross product operator
        ! result =(v, ω)× = [ω×, v×| 0,  ω×]
        pure function twist_cross_twist_matrix(t) result(c)
        !argument specification
            real(real64), intent(in) :: t(nspc)
            real(real64) :: c(nspc, nspc), wx(nvec,nvec), vx(nvec,nvec)
        !code start
            vx = .x. t(1:3)
            wx = .x. t(4:6)
            
            c(1:3,1:3) = wx
            c(1:3,4:6) = vx            
            c(4:6,4:6) = wx
        end function        

        ! Returns the (twist × wrench) = wrench cross product
        ! result =(v, ω)×(a, b) = (ω×a,v×a + ω×b)
        pure function twist_cross_wrench(t, f) result(c)
        !argument specification
            real(real64), intent(in) :: t(nspc), f(nspc)
            real(real64) :: c(nspc)
        !code start
            c = [ t(4:6) .x. f(1:3), (t(1:3) .x. f(1:3)) + (t(4:6) .x. f(4:6)) ]
        end function     
        
        ! Returns the (twist × wrench) = wrench cross product operator
        ! result =(v, ω)× = [ω×, 0| v×,  ω×]
        pure function twist_cross_wrench_matrix(t) result(c)
        !argument specification
            real(real64), intent(in) :: t(nspc)
            real(real64) :: c(nspc, nspc), wx(nvec,nvec), vx(nvec,nvec)
        !code start
            vx = .x. t(1:3)
            wx = .x. t(4:6)
            
            c(1:3,1:3) = wx
            c(4:6,1:3) = vx            
            c(4:6,4:6) = wx
        end function        
        
        pure function material_acceleration(a,v,cg) result(acm)
        real(real64), intent(in) :: a(nspc), v(nspc)
        real(real64), optional, intent(in) :: cg(nvec)
        real(real64) :: acm(nvec), acc(nvec), alp(nvec), vee(nvec), omg(nvec)
        
            omg = v(4:6)
            alp = a(4:6)
            if(present(cg)) then
                vee = v(1:3) + ( omg .x. cg )
                acc = a(1:3) + ( alp .x. cg )
            else
                vee = v(1:3)
                acc = a(1:3)                
            end if
            acm = acc + (omg .x. vee)
            
        end function
        
    END MODULE