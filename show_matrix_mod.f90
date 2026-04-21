    module show_matrix_mod
    use, intrinsic :: iso_fortran_env
    implicit none
    
    interface show
        procedure show_vector_i, show_vector_r, show_vector_d
        procedure show_matrix_i, show_matrix_r, show_matrix_d
    end interface

    
    contains
    
        subroutine show_vector_i(v, w)
    ! Display the vector 'v' in a single column
    !   v : the array of real numbers
    !   w : the column width. default = 5
    !   s : sig. figures w-5 (calculated)
        integer, intent(in) :: v(:)
        integer, intent(in), optional :: w
        integer :: i,n,wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 5
        end if
        n = size(v)
        write( fmt, "(a,g0,a)") "(*(g",wt,".0))"
        write( * , fmt ) ( v(i), new_line("A"), i=1,n )    
    end subroutine
       
    subroutine show_vector_r(v, w)
    ! Display the vector 'v' in a single column
    !   v : the array of real numbers
    !   w : the column width. default = 12
    !   s : sig. figures w-5 (calculated)
        real(real32), intent(in) :: v(:)
        integer, intent(in), optional :: w
        integer :: i,n,dg,wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 13
        end if
        dg = wt - 7
        n = size(v)
        write( fmt, "(a,g0,a,g0,a)") "(*(g",wt,".",dg,"))"
        write( * , fmt ) ( v(i), new_line("A"), i=1,n )    
    end subroutine
    
    subroutine show_vector_d(v, w)
    ! Display the vector 'v' in a single column
    !   v : the array of real numbers
    !   w : the column width. default = 12
    !   s : sig. figures w-5 (calculated)
        real(real64), intent(in) :: v(:)
        integer, intent(in), optional :: w
        call show_vector_r(real(v),w)
    end subroutine
    
    subroutine show_matrix_i(A, w)
    ! Display the matrix 'A' in columns
    !   A : the array of integers
    !   w : the column width. default = 5
        integer, intent(in) :: A(:,:)
        integer, intent(in), optional :: w
        integer :: i,j,n,m, wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 5
        end if
        n = size(A,1)
        m = size(A,2)
        write( fmt, "(a,g0,a)") "(*(g",wt,".0))"        
        write( * , fmt ) ( (A(i,j),j=1,m), new_line("A"), i=1,n )
    end subroutine
    
    subroutine show_matrix_r(A, w)
    ! Display the matrix 'A' in columns
    !   A : the array of real numbers
    !   w : the column width. default = 12
    !   s : sig. figures w-5 (calculated)
        real(real32), intent(in) :: A(:,:)
        integer, intent(in), optional :: w
        integer :: i,j,n,m,dg,wt
        character(len=16) :: fmt
        if(present(w)) then
            wt = w
        else
            wt = 13
        end if
        dg = wt - 7
        n = size(A,1)
        m = size(A,2)
        write( fmt, "(a,g0,a,g0,a)") "(*(g",wt,".",dg,"))"
        write( * , fmt ) ( (A(i,j),j=1,m), new_line("A"), i=1,n )
    end subroutine
    
    subroutine show_matrix_d(A,w)
    ! Display the matrix 'A' in columns
    !   A : the array of dble numbers
    !   w : the column width. default = 12
    ! Converts 'A' into single precision and calls `show_matrix_r`
        real(real64), intent(in) :: A(:,:)
        integer, intent(in), optional :: w
        call show_matrix_r(real(A),w)
    end subroutine

    
    end module