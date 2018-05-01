subroutine read_up(up, n, j)
  implicit none
  integer :: n, np
  integer :: i, j
  double precision :: up(10000)
  character*12 :: filename
  write(filename,'("up_",I1)')j
  print*, filename
  open(unit = 2, file = filename)
  read(2, *)
  read(2, *) n, np
  do i=1,n
     read(2, *) up(i)
  end do
  close(2)
end subroutine read_up

program test
  use libgemsrom
  implicit none
  type(gemsrom) :: gr
  double precision :: up(10000), upt(10000)
  integer :: n, i, dim
  double precision :: uhat(2000), uhat_sum(2000)
  integer :: np
  dim = 2000

  ! Create an object of type foo
  gr = gemsrom()
  uhat_sum(:) = 0.0
  ! Call bound procedures (member functions)
  call gr%initialize
  do i=0, np-1
     call read_up(up, n, i)
     print *, i, n
     call gr%get_uhat(i, up, n, uhat, dim)
     uhat_sum(:) = uhat_sum(:) +  uhat(:)
  end do
  print *, uhat_sum(2)
  call read_up(up, n, 0)
  call gr%get_u(0, upt, n, uhat_sum, dim);
  print *, upt(1)
#ifdef __GNUC__
  call gr%delete
#endif
end program test
