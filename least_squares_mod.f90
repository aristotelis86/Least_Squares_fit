
   module least_squares_mod
   ! Module implementing Least Squares Method 
   ! for approximating data sets.
   !
   ! Public routines:
   !              - lsq_poly
   
   implicit none

   private

   public   :: lsq_poly
   public   :: gaussian_elimination_solver
   
   contains 
   
   ! =====================================================
   subroutine lsq_poly(x, y, ord, coefs, info)
   ! Subroutine to act as interface for the rest of the 
   ! methods to calculate a polynomial fit based on
   ! Least Squares Method.
   real, intent(in)              :: x(:), y(:)
   integer, intent(in)           :: ord
   real, allocatable, intent(out)   :: coefs(:)
   integer, optional, intent(out)   :: info
   
   integer           :: n
   real, allocatable :: x_local(:)
   
   if (ord > 3 .or. ord < 0) then
      ! Handle "illegal" cases
      info = -1
      allocate(coefs(abs(ord)))
      coefs = -99999.9
      return
   end if
   
   n = size(x)
   if (ord == 0) then
      ! Handle trivial case
      info = 0
      allocate(coefs(1))
      coefs = sum(y)/(1.0*n)
      return
   end if
   
   if (n < ord+1) then
      ! Handle another illegal case
      info = 0
      allocate(coefs(abs(ord)))
      coefs = -99999.9
      return
   end if
   
   allocate(coefs(ord+1))
   allocate(x_local, source=x)
   !x_local = x_local/maxval(abs(x))
   
   if (ord == 1) then
      call linear_lsq(x_local, y, coefs, info)
      return
   else if (ord == 2) then
      call quadratic_lsq(x_local, y, coefs, info)
      return
   else
      call cubic_lsq(x_local, y, coefs, info)
      return
   end if
   
   end subroutine lsq_poly
   
   ! =====================================================
   subroutine linear_lsq(x, y, coefs, info)
   ! Solving system for linear approximation
   !        A * coefs = B
   !            |    n        sum(x)  |
   !        A = |  sum(x)    sum(x*x) |
   !
   !            |  sum(y)  |
   !        B = | sum(x*y) |
   ! 
   !        y = coefs(1) + coefs(2)*x
   
   real, intent(in)     :: x(:), y(:)
   real, intent(out)    :: coefs(2)
   integer, optional, intent(out)   :: info

   real                 :: sumx, sumy, sumxy, sumx2
   real                 :: det, detc(2)
   integer              :: n

   n = size(x)
   if (n < 2) then
      info = -1
      return
   end if
   
   sumx = sum(x)
   sumy = sum(y)
   sumx2 = dot_product(x,x)
   sumxy = dot_product(x,y)

   det = n*sumx2 - sumx*sumx
   detc(1) = sumy*sumx2 - sumx*sumxy
   detc(2) = n*sumxy - sumy*sumx

   if (abs(det) > 1.0d-8) then
      info = 0
      coefs = detc/det
   else
      info = -1
      coefs = -999999.0
   end if
   end subroutine linear_lsq

   ! =====================================================
   subroutine quadratic_lsq(x, y, coefs, info)
   ! Solving system for quadratic approximation
   !        A * coefs = B
   !            |    n       sum(x)   sum(x^2) |
   !        A = |  sum(x)   sum(x^2)  sum(x^3) |
   !            | sum(x^2)  sum(x^3)  sum(x^4) |
   !
   !            |   sum(y)   |
   !        B = |  sum(y*x)  |
   !            | sum(y*x^2) |
   ! 
   !        y = coefs(1) + coefs(2)*x + coefs(3)*x^2
   
   real, intent(in)     :: x(:), y(:)
   real, intent(out)    :: coefs(3)
   integer, optional, intent(out)   :: info
   
   integer           :: n
   real              :: sumx, sumx2, sumx3, sumx4
   real              :: sumy, sumxy, sumyx2
   real              :: A(3,3), B(3)
   
   n = size(x)
   if (n < 3) then
      info = -1
      return
   end if
   sumx = sum(x)
   sumx2 = dot_product(x,x)
   sumx3 = sum(x*x*x)
   sumx4 = sum(x*x*x*x)
   sumy = sum(y)
   sumxy = dot_product(x,y)
   sumyx2 = sum(y*x*x)
   
   A(1,1) = n
   A(1,2) = sumx
   A(1,3) = sumx2
   
   A(2,1) = sumx
   A(2,2) = sumx2
   A(2,3) = sumx3
   
   A(3,1) = sumx2
   A(3,2) = sumx3
   A(3,3) = sumx4
   
   B(1) = sumy
   B(2) = sumxy
   B(3) = sumyx2
   
   call cramer(A, B, info)
   coefs = B
   end subroutine quadratic_lsq
   
   ! ====================================================================================
   subroutine cramer(a, b, info)
   ! Solve a simple 3x3 linear system
   ! using Cramer's rule
   real, intent(in)        :: a(3,3)
   real, intent(inout)     :: b(3)
   integer, intent(out)    :: info
   
   real                    :: D, D1, D2, D3
   real                    :: a1(3,3), a2(3,3), a3(3,3)
   
   a1 = a
   a2 = a
   a3 = a
   
   D = determinant(a)
   if (abs(D) > 1.0d-8) then
      info = 0
      a1(:,1) = b
      a2(:,2) = b
      a3(:,3) = b
      
      D1 = determinant(a1)
      D2 = determinant(a2)
      D3 = determinant(a3)
      
      b(1) = D1/D
      b(2) = D2/D
      b(3) = D3/D
   else
      info = -1
   end if
   contains 
   function determinant(a) result(val)
   real              :: val
   real, intent(in)  :: a(3,3)
   
   val = 0.0
   val = val + a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2))
   val = val - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1))
   val = val + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
   end function determinant
   end subroutine cramer
   
   ! =====================================================
   subroutine cubic_lsq(x, y, coefs, info)
   ! Solving system for quadratic approximation
   !        A * coefs = B
   !            |    n       sum(x)   sum(x^2)  sum(x^3) |
   !        A = |  sum(x)   sum(x^2)  sum(x^3)  sum(x^4) |
   !            | sum(x^2)  sum(x^3)  sum(x^4)  sum(x^5) |
   !            | sum(x^3)  sum(x^4)  sum(x^5)  sum(x^6) |
   !
   !            |   sum(y)   |
   !        B = |  sum(y*x)  |
   !            | sum(y*x^2) |
   !            | sum(y*x^3) |
   ! 
   !        y = coefs(1) + coefs(2)*x + coefs(3)*x^2 + coefs(4)*x^3
   
   real, intent(in)     :: x(:), y(:)
   real, intent(out)    :: coefs(4)
   integer, optional, intent(out)   :: info
   
   real, allocatable    :: Alpha(:,:), Beta(:)
   integer              :: n

   n = size(x)
   if (n < 4) then
      info = -1
      return
   end if
   call construct_coef_matrix(x, y, 4, Alpha, Beta)
   
   call gaussian_elimination_solver(Alpha, Beta, coefs)
   end subroutine cubic_lsq
   
   ! =====================================================
   subroutine construct_coef_matrix(x, y, gen, A, B)
   ! Construct the matrices of coefficients for
   ! higher order (>3) systems.
   real, intent(in)        :: x(:), y(:)
   integer, intent(in)     :: gen
   real, allocatable, intent(out)   :: A(:,:), B(:)
   
   integer        :: i, j
   
   allocate(A(gen,gen), B(gen))
   
   do i = 1, gen
      do j = 1, gen
         A(i,j) = sum(x**(i+j-2))
      end do
      B(i) = sum(y*x**(i-1))
   end do
   end subroutine construct_coef_matrix
   
   ! =====================================================
   subroutine gaussian_elimination_solver(A, B, X)
   ! Implementing a naive gaussian elimination algorithm
   real        :: A(:,:), B(:), X(:)
   
   integer     :: h, k, nrow, ncol, i, j, i_max
   real        :: fval 
   real, allocatable :: A_aug(:,:)
   
   
   nrow = size(A,dim=1)
   ncol = size(A,dim=2)
   
   allocate(A_aug(nrow,ncol+1))
   A_aug(:,1:ncol) = A
   A_aug(:,ncol+1) = B
   
   h = 1
   k = 1
   
   endless: do
      if (h > nrow .or. k > ncol) exit endless
      
      i_max = maxloc(abs(A_aug(h:nrow,k)),1) + h - 1
      if (A_aug(i_max,k) < 1.0d-8) then
         k = k + 1
      else
         
         call swap_rows(h, i_max)
         do i = h+1, nrow
            fval = A_aug(i,k)/A_aug(h,k)
            A_aug(i,k) = 0.0
            do j = k+1, ncol
               A_aug(i,j) = A_aug(i,j) - A_aug(h,j)*fval
            end do
         end do
         
         h = h + 1
         k = k + 1
      end if
   end do endless
   
   ! Finally solve using Backward substitution
   call back_sub()
   
   contains
   ! ****************
   subroutine swap_rows(from, tow)
   integer           :: from, tow
   real, allocatable :: temp_row(:)
   
   allocate(temp_row(ncol))
   temp_row = A_aug(from,:)
   A_aug(from,:) = A_aug(tow,:)
   A_aug(tow,:) = temp_row
   deallocate(temp_row)
   end subroutine swap_rows
   
   ! ****************
   subroutine back_sub()
   real           :: sum_val
   A = A_aug(:,1:ncol)
   B = A_aug(:,ncol+1)
   X(nrow) = B(nrow)/A(nrow, ncol)
   
   do i = nrow-1, 1, -1
      sum_val = 0.0
      do j = i+1, nrow
         sum_val = sum_val + A(i, j)*X(j)
      end do
      X(i) = (1./A(i,i)) * (B(i) - sum_val)
   end do   
   end subroutine back_sub
   end subroutine gaussian_elimination_solver

   end module least_squares_mod