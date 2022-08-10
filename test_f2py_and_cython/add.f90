! example modified from https://numpy.org/devdocs/f2py/f2py-examples.html
! compile with python -m numpy.f2py -c -m add add.f
! if it compiles, f2py is working
subroutine zadd(a,b,c,n) ! in :add:add.f
      real, dimension(n) :: a
      real, dimension(n) :: b
      real, intent(out),dimension(n) :: c
      integer, intent(in) :: n
      integer :: i
      do i=1, n 
            c(i) = a(i) + b(i)
      end do
end subroutine zadd