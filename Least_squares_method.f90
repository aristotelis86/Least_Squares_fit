
!****************************************************************************
!
!  PROGRAM: Least_squares_method
!
!  PURPOSE:  Test the module for Least Squares method using data sets.
!
!****************************************************************************

   program Least_squares_method
   use least_squares_mod
   implicit none
   
   ! Variables
   real, allocatable       :: xlin(:), ylin(:)
   real, allocatable       :: xquad(:), yquad(:)
   real, allocatable       :: xcub(:), ycub(:)
   
   integer                 :: nlin, nquad, ncub
   integer                 :: i, k, info
   real, allocatable       :: coefs_lin(:), coefs_quad(:), coefs_cub(:)
   real                    :: coefs_lin_d(2), coefs_quad_d(3), coefs_cub_d(4)
   
   real                    :: A(3,3), B(3), X(3), A_orig(3,3), B_sol(3)
   
   
   ! Body of Least_squares_method
   print *, 'Hello World'
   print *, 'We will demonstrate Least Squares method for: '
   print *, '        - Linear data + random positive error'
   print *, '        - Quadratic data + random positive error'
   print *, '        - Cubic data + random positive error'
   
   ! Open files
   open(unit=30, file='./linear_test.txt', status='old')
   open(unit=40, file='./quadratic_test.txt', status='old')
   open(unit=50, file='./cubic_test.txt', status='old')
   
   ! Read headers
   read(30,'("linear: y=a*x+b, a=",F3.1,", b=",F3.1)') coefs_lin_d(1), coefs_lin_d(2)           
   read(40,'("quadratic: y=a*x*x + b*x + c, a=",F3.1,", b=",F3.1,", c=",F3.1)') coefs_quad_d(1), coefs_quad_d(2), coefs_quad_d(3)
   read(50,'("cubic: y=a*x^3 + b*x^2 + c*x + d, a=",F3.1,", b=",F3.1,", c=",F3.1,", d=",F3.1)') coefs_cub_d(1), coefs_cub_d(2), coefs_cub_d(3), coefs_cub_d(4)
   
   ! Read sizes
   read(30,*) nlin
   read(40,*) nquad
   read(50,*) ncub
   
   ! Allocate space
   allocate(xlin(nlin), ylin(nlin))
   allocate(xquad(nquad), yquad(nquad))
   allocate(xcub(ncub), ycub(ncub))
   
   ! Read actual data and close files
   do i = 1, nlin
      read(30,*) xlin(i), ylin(i)
   end do
   close(30)
   
   do i = 1, nquad
      read(40,*) xquad(i), yquad(i)
   end do
   close(40)
   
   do i = 1, ncub
      read(50,*) xcub(i), ycub(i)
   end do
   close(50)
   
   ! Perform linear fit and display results
   call lsq_poly(xlin, ylin, 1, coefs_lin, info)
   if (info == 0) then
      write(*,'(2(A,F8.3))') 'The coefficients for linear fit should be: a=', coefs_lin_d(1), ', b=', coefs_lin_d(2)
      write(*,'(2(A,F8.3))') 'The coefficients found for linear fit are: a=', coefs_lin(2), ', b=', coefs_lin(1)
   else
      write(*,*) 'Linear fit could not be determined, please check data.'
   end if
   
   write(*,*)
   ! Perform quadratic fit and display results
   call lsq_poly(xquad, yquad, 2, coefs_quad, info)
   if (info == 0) then
      write(*,'(3(A,F8.3))') 'The coefficients for quadratic fit should be: a=', coefs_quad_d(1), ', b=', coefs_quad_d(2), ', c=', coefs_quad_d(3)
      write(*,'(3(A,F8.3))') 'The coefficients found for quadratic fit are: a=', coefs_quad(3), ', b=', coefs_quad(2), ', c=', coefs_quad(1)
   else
      write(*,*) 'Quadratic fit could not be determined, please check data.'
   end if
   
   write(*,*)
   ! Perform quadratic fit and display results
   call lsq_poly(xcub, ycub, 3, coefs_cub, info)
   if (info == 0) then
      write(*,'(4(A,F8.3))') 'The coefficients for cubic fit should be: a=', coefs_cub_d(1), ', b=', coefs_cub_d(2), ', c=', coefs_cub_d(3), ', d=', coefs_cub_d(4)
      write(*,'(4(A,F8.3))') 'The coefficients found for cubic fit are: a=', coefs_cub(4), ', b=', coefs_cub(3), ', c=', coefs_cub(2), ', d=', coefs_cub(1)
      write(*,*) coefs_cub(2)
      write(*,*) coefs_cub(1)
   else
      write(*,*) 'Quadratic fit could not be determined, please check data.'
   end if
   !
   !write(*,*)
   !write(*,*)
   !write(*,*)
   !
   !A(1,1) = 1.6
   !A(1,2) = 4.23
   !A(1,3) = 0.5
   !
   !A(2,1) = 0.1
   !A(2,2) = 6.06
   !A(2,3) = 3.87
   !
   !A(3,1) = 4.0
   !A(3,2) = 0.0
   !A(3,3) = 1.0
   !
   !B(1) = 1.
   !B(2) = 2. 
   !B(3) = 3.
   !
   !A_orig = A
   !write(*,*) 'A = '
   !do i = 1, 3
   !   write(*,'(3(F5.2),A,F5.2)') A(i,:), '  |  ', B(i)
   !end do
   !
   !call gaussian_elimination_solver(A, B, X)
   !write(*,*) 'A = '
   !do i = 1, 3
   !   write(*,'(3(F5.2),A,F5.2)') A(i,:), '  |  ', B(i)
   !end do
   !
   !B_sol = 0.0
   !do i = 1, 3
   !   do k = 1, 3
   !      B_sol(i) = B_sol(i) + A(i,k)*X(k)
   !   end do
   !end do
   !
   !write(*,*)
   !write(*,*)
   !write(*,*) 'Solution:'
   !write(*,*) ' X = ', X
   !write(*,*) ' B = ', B
   !write(*,*) ' Bsol = ', B_sol
   !
   read(*,*)
   
   end program Least_squares_method

