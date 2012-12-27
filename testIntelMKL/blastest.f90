! 
! Test Intel's MKL in Fortran
! with Fortran 
!
!   Charles O'Neill 8 April 2010
!    charles.oneill@gmail.com
!    www.caselab.okstate.edu
!
! Copyright (c) 2010 Charles O'Neill
!
! Permission is hereby granted, free of charge, to any person
! obtaining a copy of this software and associated documentation
! files (the "Software"), to deal in the Software without
! restriction, including without limitation the rights to use,
! copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the
! Software is furnished to do so, subject to the following
! conditions:
!
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
! OTHER DEALINGS IN THE SOFTWARE.
!

program blastest
  implicit none

    integer i
    integer :: k
    write(*,*) "blas testing"
    
    do i=1,10
       k = 2**i   
       call mv2(k)
    enddo

contains

  subroutine mv2(n)
    use ClockTiming
    use blas95
    integer, intent(in) :: n
    integer :: i,j,k
    real(8),allocatable :: a(:), b(:),c(:),d(:),t(:)
    real(8),allocatable :: m(:,:),mm(:,:),mmm(:,:)
    real(8) :: s, alpha
    real(8) :: e

    type (ClockTime) :: Timer
    write(*,*) 
    write(*,*) "mkl", n

    allocate(a(n))
    allocate(b(n))
    allocate(c(n))
    allocate(d(n))
    allocate(t(n))    
    allocate(m(n,n))
    allocate(mm(n,n))
    allocate(mmm(n,n))

    a = 1.0d0
    b = 0.0d0
    c = 0.0d0
    d = 1.0d0
    s = 0.0d0
    e = 0.0d0
    m = 1.0d0
    mm = 1.0d0
    alpha = 1.0d0

    ! Native Fortran 77 (non-optimal)
    do i=1,n
       t(i) = d(i)
       d(i) = a(i) 
       a(i) = t(i)
    enddo
    call set(Timer)
    do i = 1,n
       b = a
    enddo
    do i = 1,n
       b(i) = b(i) + a(i)*alpha
    enddo
    do i=1,n
       s = s + a(i)*b(i)
    enddo
    do i=1,n
       do j=1,n
          c(i) = c(i) + m(i,j) * a(j)
       enddo
    enddo
    do i=1,n
       e = e + c(i)**2
    enddo
    do i=1,n
       do j=1,n
          do k=1,n
             mmm(i,j) = mmm(i,j) +  m(i,k) * mm(k,j)
          enddo
       enddo
    enddo
    call printtime(Timer)
    write(*,*) "Native f77", s, sum(c), sum(mmm)

    !---------------------------------------------------------
    a = 1.0d0
    b = 0.0d0
    c = 0.0d0
    s = 0.0d0
    m = 1.0d0
    e = 0.0d0
    alpha = 1.0d0

    ! Fortran 90
    call set(Timer)
    t = a; d = a; a = t
    b = a
    b = b + a*alpha
    s = sum(a*b)
    c = matmul(m,a)
    e = sum(c*c)
    mmm = matmul(m,mm)
    call printtime(Timer)
    write(*,*) "F90", s, sum(c), sum(mmm)

    a = 1.0d0
    b = 0.0d0
    c = 0.0d0
    s = 0.0d0
    e = 0.0d0
    alpha = 1.0d0

    !---------------------------------------------------------
    ! MKL Blas
    call swap(a,d)
    call set(Timer)
    call copy(a,b)
    call axpy(a,b,alpha)
    s= dot(a,b)
    call gemv(m,a,c)
    e = nrm2(c)
    call gemm(m,mm,mmm)
    call printtime(Timer)
    write(*,*) "MKL BLAS",s, sum(c), sum(mmm)

  end subroutine mv2

end program blastest
