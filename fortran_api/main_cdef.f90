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

     subroutine gemsrom_initialize_c(gemsrom) bind(C, name="gemsrom_initialize")
       use iso_c_binding
       implicit none
       !integer(c_void) :: gemsrom_initialize_c
       ! The const qualification is translated into an intent(in)
       type(c_ptr), intent(in), value :: gemsrom
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

  end interface
