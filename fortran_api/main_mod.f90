module libgemsrom
    use iso_c_binding

    private
    public :: gemsrom

    ! Yes, include is a keyword in Fortran !
    include "main_cdef.f90"

    ! We'll use a Fortan type to represent a C++ class here, in an opaque maner
    type gemsrom
        private
        type(c_ptr) :: ptr ! pointer to the Gemsrom class
    contains
        ! We can bind some functions to this type, allowing for a cleaner syntax.
#if defined(__GNUC__) || defined(__INTEL_COMPILER)
        procedure :: delete => delete_gemsrom_polymorph ! Destructor for gfortran
#else
        final :: delete_gemsrom ! Destructor
#endif
        ! Function member
        procedure :: initialize => gemsrom_initialize
        procedure :: get_u => gemsrom_get_u
        procedure :: get_uhat => gemsrom_get_uhat
        procedure :: get_deim_n => gemsrom_get_deim_n
        procedure :: get_deim_id => gemsrom_get_deim_id
        procedure :: get_deim => gemsrom_get_deim
        procedure :: renormalize => gemsrom_renormalize

     end type gemsrom

    ! This function will act as the constructor for gemsrom type
    interface gemsrom
        procedure create_gemsrom
    end interface

contains ! Implementation of the functions. We just wrap the C function here.
    function create_gemsrom()
        implicit none
        type(gemsrom) :: create_gemsrom
        create_gemsrom%ptr = create_gemsrom_c()
    end function

    subroutine delete_gemsrom(this)
        implicit none
        type(gemsrom) :: this
        call delete_gemsrom_c(this%ptr)
    end subroutine

    ! Bounds procedure needs to take a polymorphic (class) argument
    subroutine delete_gemsrom_polymorph(this)
        implicit none
        class(gemsrom) :: this
        call delete_gemsrom_c(this%ptr)
    end subroutine

    subroutine gemsrom_initialize(this, partition_id)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: partition_id
      call gemsrom_initialize_c(this%ptr, partition_id)
    end subroutine gemsrom_initialize


    subroutine gemsrom_get_uhat(this, partition_id, q, nq, qhat, n_mode)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: partition_id, nq, n_mode
      double precision :: q(nq), qhat(n_mode)
      call gemsrom_get_uhat_c(this%ptr, partition_id, q, nq, qhat, n_mode)
    end subroutine gemsrom_get_uhat

    subroutine gemsrom_get_u(this, partition_id, q, nq, qhat, n_mode)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: partition_id, nq, n_mode
      double precision :: q(nq), qhat(n_mode)
      call gemsrom_get_u_c(this%ptr, partition_id, q, nq, qhat, n_mode)
    end subroutine gemsrom_get_u

    subroutine gemsrom_get_deim_n(this, partition_id, in)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: partition_id, in
      call gemsrom_get_deim_n_c(this%ptr, partition_id, in)
    end subroutine gemsrom_get_deim_n



    subroutine gemsrom_get_deim_id(this, partition_id, isize, ilocal_id, ivar)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: partition_id, isize
      integer :: ilocal_id(:), ivar(:)
      call gemsrom_get_deim_id_c(this%ptr, partition_id, isize, ilocal_id, ivar)
    end subroutine gemsrom_get_deim_id
    
    subroutine gemsrom_get_deim(this, partition_id, isize, r_s, deim_r)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: partition_id, isize
      double precision :: r_s(:), deim_r(64000)
      call gemsrom_get_deim_c(this%ptr, partition_id, isize, r_s, deim_r)
    end subroutine gemsrom_get_deim

    subroutine gemsrom_renormalize(this, isize, x, y)
      implicit none
      class(gemsrom), intent(in) :: this
      integer :: isize
      double precision :: x(isize), y(isize)
      call gemsrom_renormalize_c(this%ptr, isize, x, y)
    end subroutine gemsrom_renormalize

end module
