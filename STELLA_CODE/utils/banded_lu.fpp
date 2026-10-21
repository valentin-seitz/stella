!###############################################################################
!                    LU factorisation of complex band matrices
!###############################################################################
! Factorises and solves the square complex band systems of stella, with n rows,
! kl sub-diagonals and ku super-diagonals, stored as in LAPACK: A(i,j) is
! ab(band_row(i, j, kl, ku), j), and the array has band_ldab(kl, ku) rows.
!
! By default, LAPACK's zgbtrf and zgbtrs are used. With USE_PLASMA_BANDED=on
! (Make) or -DSTELLA_ENABLE_PLASMA_BANDED=ON (CMake), which define the
! PLASMA_BANDED preprocessor flag, PLASMA's plasma_zgbtrf and plasma_zgbtrs are
! used instead. Then:
!  - PLASMA works on tiles of <band_nb> x <band_nb> elements, and needs more
!    rows in <ab> than LAPACK, hence <band_ldab>.
!  - The factors are not stored as in LAPACK: within every panel of <band_nb>
!    columns, the row swaps are also applied to the earlier columns of L, as in
!    a dense LU. They can only be used by <band_lu_solve>.
!  - plasma_zgbtrf opens its own OpenMP parallel region and uses the shared
!    PLASMA context, so it can't be called from inside a parallel region. The
!    callers that factorise in a loop make their loop serial when
!    <band_lu_parallel_inside> is true.
!  - plasma_zgbtrs has the same restriction, so inside a parallel region
!    <band_lu_solve> uses <serial_gbtrs_z>, which does the same operations as
!    plasma_zgbtrs with serial loops.
!###############################################################################
module banded_lu

   implicit none

   private

   public :: band_ldab, band_row
   public :: band_lu_factorise, band_lu_solve
   public :: band_lu_parallel_inside

#ifdef PLASMA_BANDED
   ! The tile size of PLASMA for the band routines. The band systems of stella
   ! are small, so the default tile size of PLASMA (256) would waste memory and time
   integer, parameter :: band_nb = 16

   ! Whether <band_lu_factorise> opens its own OpenMP parallel region
   logical, parameter :: band_lu_parallel_inside = .true.
#else
   logical, parameter :: band_lu_parallel_inside = .false.
#endif

contains

   !****************************************************************************
   !                 Leading dimension of the band storage
   !****************************************************************************
   ! The number of rows of the array <ab> that stores a band matrix with <kl>
   ! sub-diagonals and <ku> super-diagonals and its LU factors.
   !****************************************************************************
   pure integer function band_ldab(kl, ku)

      implicit none

      integer, intent(in) :: kl, ku

#ifdef PLASMA_BANDED
      ! (tiles in the upper band + tiles in the lower band + 1) * band_nb,
      ! see plasma/compute/zgbtrf.c
      band_ldab = ((kl + ku + band_nb - 1) / band_nb + (kl + band_nb - 1) / band_nb + 1) * band_nb
#else
      band_ldab = 2 * kl + ku + 1
#endif

   end function band_ldab

   !****************************************************************************
   ! The row of <ab> that stores the element A(i,j) of the band matrix, which
   ! is stored in column j of <ab>. The first <kl> rows are needed for the fill-in
   ! of the LU factorisation.
   !****************************************************************************
   pure integer function band_row(i, j, kl, ku)

      implicit none

      integer, intent(in) :: i, j, kl, ku

      band_row = kl + ku + 1 + i - j

   end function band_row

   !****************************************************************************
   !                       LU factorisation of a band matrix
   !****************************************************************************
   ! On entry, <ab> contains the n x n band matrix A, see <band_row>, and zeros
   ! everywhere else. On exit, it contains the LU factors of A, and <ipiv> the
   ! row swaps, which are used by <band_lu_solve>. <info> is the status of
   ! zgbtrf or plasma_zgbtrf, i.e., 0 on success and i > 0 if U(i,i) is zero.
   !****************************************************************************
   subroutine band_lu_factorise(n, kl, ku, ab, ipiv, info)

      use mp, only: mp_abort
#ifdef PLASMA_BANDED
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_zgbtrf, stella_plasma_usable
#endif

      implicit none

      integer, intent(in) :: n, kl, ku
      complex, dimension(:, :), intent(in out) :: ab
      integer, dimension(:), intent(out) :: ipiv
      integer, intent(out) :: info

#ifdef PLASMA_BANDED
      integer(c_int), dimension(n) :: ipiv_c
      integer :: nb_old
#endif

      !-------------------------------------------------------------------------

      if (size(ab, 1) < band_ldab(kl, ku)) call mp_abort('band_lu_factorise: ab has too few rows')

#ifdef PLASMA_BANDED
      ! The factors from plasma_zgbtrf can't be solved by LAPACK, so we can't fall back to zgbtrf
      if (.not. stella_plasma_usable()) &
         call mp_abort('band_lu_factorise: PLASMA has to be initialised, and this can not be called ' &
                       //'inside an OpenMP parallel region')
      call set_band_nb(nb_old)
      info = int(plasma_zgbtrf(int(n, c_int), int(n, c_int), int(kl, c_int), int(ku, c_int), &
                               ab, int(size(ab, 1), c_int), ipiv_c))
      call reset_nb(nb_old)
      ipiv(1:n) = ipiv_c
#else
      call zgbtrf(n, n, kl, ku, ab, size(ab, 1), ipiv, info)
#endif

   end subroutine band_lu_factorise

   !****************************************************************************
   !                     Solve with the LU factors of a band matrix
   !****************************************************************************
   ! Solves A x = b, given the factors <ab> and <ipiv> of A from
   ! <band_lu_factorise>. On exit, <b> contains x. This can be called from inside
   ! OpenMP parallel regions.
   !****************************************************************************
   subroutine band_lu_solve(n, kl, ku, ab, ipiv, b)

#ifdef PLASMA_BANDED
      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_zgbtrs, stella_plasma_usable, PlasmaNoTrans
#endif

      implicit none

      integer, intent(in) :: n, kl, ku
      complex, dimension(:, :), intent(in) :: ab
      integer, dimension(:), intent(in) :: ipiv
      complex, dimension(:), intent(in out) :: b

      integer :: info
#ifdef PLASMA_BANDED
      integer(c_int), dimension(n) :: ipiv_c
      integer :: nb_old
#endif

      !-------------------------------------------------------------------------

#ifdef PLASMA_BANDED
      if (stella_plasma_usable()) then
         ipiv_c = int(ipiv(1:n), c_int)
         call set_band_nb(nb_old)
         info = int(plasma_zgbtrs(PlasmaNoTrans, int(n, c_int), int(kl, c_int), int(ku, c_int), 1_c_int, &
                                  ab, int(size(ab, 1), c_int), ipiv_c, b, int(size(b), c_int)))
         call reset_nb(nb_old)
      else
         call serial_gbtrs_z(n, kl, ku, band_nb, ab, ipiv, b)
      end if
#else
      call zgbtrs('N', n, kl, ku, 1, ab, size(ab, 1), ipiv, b, size(b), info)
#endif

   end subroutine band_lu_solve

#ifdef PLASMA_BANDED
   !****************************************************************************
   !               Serial solve with the band factors of PLASMA
   !****************************************************************************
   ! Does the same operations as plasma_zgbtrs (see plasma/compute/pztbsm.c),
   ! with serial loops instead of tasks, so that it can be called from inside
   ! OpenMP parallel regions. The factors are those of plasma_zgbtrf with tiles
   ! of <nb> x <nb> elements:
   !  - L y = P b: for every panel of <nb> columns, first all the row swaps of
   !    the panel are applied to b, since plasma_zgbtrf also applied them to the
   !    earlier columns of L in the panel (as in a dense LU), then the columns
   !    of L in the panel are eliminated. Column j of L can have nonzeros down
   !    to row (last row of the panel) + kl.
   !  - U x = y: U has kl + ku super-diagonals, as in LAPACK.
   ! With nb = 1 this is the same as LAPACK's zgbtrs.
   !****************************************************************************
   subroutine serial_gbtrs_z(n, kl, ku, nb, ab, ipiv, b)

      implicit none

      integer, intent(in) :: n, kl, ku, nb
      complex, dimension(:, :), intent(in) :: ab
      integer, dimension(:), intent(in) :: ipiv
      complex, dimension(:), intent(in out) :: b

      integer :: kd, j0, j1, j, i, ip, ilast, ifirst
      complex :: tmp

      !-------------------------------------------------------------------------

      ! The row of <ab> that contains the diagonal
      kd = kl + ku + 1

      ! Solve L y = P b, panel by panel
      do j0 = 1, n, nb
         j1 = min(n, j0 + nb - 1)
         do i = j0, j1
            ip = ipiv(i)
            if (ip /= i) then
               tmp = b(i)
               b(i) = b(ip)
               b(ip) = tmp
            end if
         end do
         ilast = min(n, j1 + kl)
         do j = j0, j1
            if (b(j) == 0.0) cycle
            b(j + 1:ilast) = b(j + 1:ilast) - ab(kd + 1:kd + ilast - j, j) * b(j)
         end do
      end do

      ! Solve U x = y
      do j = n, 1, -1
         b(j) = b(j) / ab(kd, j)
         ifirst = max(1, j - kl - ku)
         b(ifirst:j - 1) = b(ifirst:j - 1) - ab(kd + ifirst - j:kd - 1, j) * b(j)
      end do

   end subroutine serial_gbtrs_z

   !****************************************************************************
   ! The PLASMA tile size is a global setting, so it is set to <band_nb> for every
   ! call of the band routines and reset afterwards, so that the dense routines
   ! keep the default tile size.
   !****************************************************************************
   subroutine set_band_nb(nb_old)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_get, plasma_set, PlasmaNb

      implicit none

      integer, intent(out) :: nb_old

      integer(c_int) :: info, val

      info = plasma_get(PlasmaNb, val)
      nb_old = val
      info = plasma_set(PlasmaNb, int(band_nb, c_int))

   end subroutine set_band_nb

   subroutine reset_nb(nb_old)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_set, PlasmaNb

      implicit none

      integer, intent(in) :: nb_old

      integer(c_int) :: info

      info = plasma_set(PlasmaNb, int(nb_old, c_int))

   end subroutine reset_nb
#endif

end module banded_lu
