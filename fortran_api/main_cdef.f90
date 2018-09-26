  interface
     function create_gemsrom_c() bind(C, name="create_gemsrom")
       use iso_c_binding
       implicit none
       type(c_ptr) :: create_gemsrom_c
     end function create_gemsrom_c

     subroutine delete_gemsrom_c(gemsrom) bind(C, name="delete_gemsrom")
       use iso_c_binding
       implicit none
       type(c_ptr), value :: gemsrom
     end subroutine delete_gemsrom_c

     subroutine gemsrom_initialize_c(gemsrom, partition_id) bind(C, name="gemsrom_initialize")
       use iso_c_binding
       implicit none
       !integer(c_void) :: gemsrom_initialize_c
       ! The const qualification is translated into an intent(in)
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id
     end subroutine gemsrom_initialize_c
     
     subroutine gemsrom_get_uhat_c(gemsrom, partition_id, q, nq, qhat, n_mode) bind(C, name="gemsrom_get_uhat")
       use iso_c_binding
       implicit none
       !integer(c_void) :: gemsrom_get_uhat_c
       ! The const qualification is translated into an intent(in)
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id, nq, n_mode
       real(c_double) :: q(0:nq), qhat(0:n_mode);
     end subroutine gemsrom_get_uhat_c


     subroutine gemsrom_get_u_c(gemsrom, partition_id, q, nq, qhat, n_mode) bind(C, name="gemsrom_get_u")
       use iso_c_binding
       implicit none
       !integer(c_void) :: gemsrom_get_u_c
       ! The const qualification is translated into an intent(in)
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id, nq, n_mode
       real(c_double) :: q(0:nq), qhat(0:n_mode);
     end subroutine gemsrom_get_u_c

     subroutine gemsrom_get_deim_n_c(gemsrom, partition_id, in) bind(C, name="gemsrom_get_deim_n")
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id
       integer(c_int), intent(inout) :: in
     end subroutine gemsrom_get_deim_n_c

     subroutine gemsrom_get_deim_id_c(gemsrom, partition_id, isize, ilocal_id, ivar) bind(C, name="gemsrom_get_deim_idx")
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id, isize
       integer(c_int) :: ilocal_id(isize), ivar(isize)
     end subroutine gemsrom_get_deim_id_c

     subroutine gemsrom_get_deim_c(gemsrom, partition_id,  isize, ndof, r_s, deim_r) bind(C, name="gemsrom_calc_deim")
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id, isize, ndof
       real(c_double) :: r_s(isize), deim_r(ndof)
     end subroutine gemsrom_get_deim_c

     subroutine gemsrom_renormalize_c(gemsrom, isize, x, y) bind(C, name="gemsrom_renormalize")
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: isize
       real(c_double) :: x(isize), y(isize)
     end subroutine gemsrom_renormalize_c


     subroutine gemsrom_get_global_id_c(gemsrom, partition_id, local_id, global_id) bind(C, name="gemsrom_get_global_id")
       use iso_c_binding
       implicit none
       type(c_ptr), intent(in), value :: gemsrom
       integer(c_int), intent(in), value :: partition_id, local_id
       integer(c_int), intent(inout) :: global_id
     end subroutine gemsrom_get_global_id_c

     
  end interface
