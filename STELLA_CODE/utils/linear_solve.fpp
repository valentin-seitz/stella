module linear_solve

   implicit none

   public :: lu_decomposition
   public :: lu_back_substitution
   public :: lu_inverse
   public :: imaxloc

#ifdef PLASMA_DENSE
   private :: plasma_lu_decomposition_real
   private :: plasma_lu_decomposition_complex
   private :: plasma_lu_inverse_real
   private :: plasma_lu_inverse_complex
   private :: serial_getrs_d
   private :: serial_getrs_z
#endif

   interface lu_decomposition
      module procedure lu_decomposition_real
      module procedure lu_decomposition_complex
   end interface

   interface lu_back_substitution
      module procedure lu_back_substitution_real
      module procedure lu_back_substitution_real_complex
      module procedure lu_back_substitution_complex
      module procedure lu_back_substitution_matrix_real
      module procedure lu_back_substitution_matrix_complex
   end interface

   interface lu_inverse
      module procedure lu_inverse_real
      module procedure lu_inverse_complex
   end interface

contains

   subroutine lu_decomposition_real(lu, idx, d)

      implicit none

      real, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      real, parameter :: zero = 1.0e-20
      real, dimension(size(lu, 1)) :: vv
      real, dimension(size(lu, 2)) :: dum

      integer :: j, n, imax

#ifdef PLASMA_DENSE
      if (plasma_lu_decomposition_real(lu, idx, d)) return
#endif

      n = size(lu, 1)

      d = 1.0
      vv = maxval(abs(lu), dim=2)
      if (any(vv == 0.0)) &
         write (*, *) 'singular matrix in lu_decomposition'
      vv = 1.0 / vv
      do j = 1, n
         imax = (j - 1) + imaxloc(vv(j:n) * abs(lu(j:n, j)))
         if (j /= imax) then
            dum = lu(imax, :)
            lu(imax, :) = lu(j, :)
            lu(j, :) = dum
            d = -d
            vv(imax) = vv(j)
         end if
         idx(j) = imax
         if (lu(j, j) == 0.0) lu(j, j) = zero
         lu(j + 1:n, j) = lu(j + 1:n, j) / lu(j, j)
         lu(j + 1:n, j + 1:n) = lu(j + 1:n, j + 1:n) - spread(lu(j + 1:n, j), 2, n - j) &
                                * spread(lu(j, j + 1:n), 1, n - j)
      end do

   end subroutine lu_decomposition_real

#ifdef PLASMA_DENSE
   !****************************************************************************
   !                     LU decomposition with plasma_dgetrf
   !****************************************************************************
   ! The real version of <plasma_lu_decomposition_complex>, see there.
   !****************************************************************************
   logical function plasma_lu_decomposition_real(lu, idx, d) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_dgetrf, stella_plasma_usable

      implicit none

      real, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      integer(c_int), dimension(size(lu, 1)) :: ipiv
      integer(c_int) :: n, info
      integer :: j

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. stella_plasma_usable()) return

      n = int(size(lu, 1), c_int)
      info = plasma_dgetrf(n, n, lu, n, ipiv)
      if (info /= 0) return

      ! The sign of the permutation, i.e., -1 for every row swap
      d = 1.0
      do j = 1, n
         idx(j) = ipiv(j)
         if (ipiv(j) /= j) d = -d
      end do
      success = .true.

   end function plasma_lu_decomposition_real
#endif

   subroutine lu_decomposition_complex(lu, idx, d)

      implicit none

      complex, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      real, parameter :: zero = 1.0e-20
      real, dimension(size(lu, 1)) :: vv
      complex, dimension(size(lu, 2)) :: dum

      integer :: j, n, imax

#ifdef PLASMA_DENSE
      if (plasma_lu_decomposition_complex(lu, idx, d)) return
#endif

      n = size(lu, 1)

      d = 1.0
      vv = maxval(cabs(lu), dim=2)
      if (any(vv == 0.0)) &
         write (*, *) 'singular matrix in lu_decomposition'
      vv = 1.0 / vv
      do j = 1, n
         imax = (j - 1) + imaxloc(vv(j:n) * cabs(lu(j:n, j)))
         if (j /= imax) then
            dum = lu(imax, :)
            lu(imax, :) = lu(j, :)
            lu(j, :) = dum
            d = -d
            vv(imax) = vv(j)
         end if
         idx(j) = imax
         if (lu(j, j) == 0.0) lu(j, j) = zero
         lu(j + 1:n, j) = lu(j + 1:n, j) / lu(j, j)
         lu(j + 1:n, j + 1:n) = lu(j + 1:n, j + 1:n) - spread(lu(j + 1:n, j), 2, n - j) &
                                * spread(lu(j, j + 1:n), 1, n - j)
      end do

   end subroutine lu_decomposition_complex

#ifdef PLASMA_DENSE
   !****************************************************************************
   !                     LU decomposition with plasma_zgetrf
   !****************************************************************************
   ! PLASMA uses partial pivoting, instead of the implicit (scaled) partial
   ! pivoting of <lu_decomposition_complex>, so the factors can differ, but they
   ! are stored in the same way: L (unit diagonal) below and U on and above the
   ! diagonal, and idx(j) is the row that was swapped with row j at step j.
   ! Hence the factors can be used by <lu_back_substitution>.
   !
   ! Returns .false. if PLASMA was not used, in which case <lu> is unchanged:
   !  - if PLASMA is not initialised;
   !  - inside an OpenMP parallel region, since plasma_zgetrf uses the shared
   !    PLASMA context and opens its own parallel region;
   !  - if the matrix is singular, since then plasma_zgetrf does not return the
   !    factors. The stella routine replaces zero pivots by a small number instead.
   !****************************************************************************
   logical function plasma_lu_decomposition_complex(lu, idx, d) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_zgetrf, stella_plasma_usable

      implicit none

      complex, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      integer(c_int), dimension(size(lu, 1)) :: ipiv
      integer(c_int) :: n, info
      integer :: j

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. stella_plasma_usable()) return

      n = int(size(lu, 1), c_int)
      info = plasma_zgetrf(n, n, lu, n, ipiv)
      if (info /= 0) return

      ! The sign of the permutation, i.e., -1 for every row swap
      d = 1.0
      do j = 1, n
         idx(j) = ipiv(j)
         if (ipiv(j) /= j) d = -d
      end do
      success = .true.

   end function plasma_lu_decomposition_complex
#endif

   subroutine lu_back_substitution_real(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:), intent(in out) :: b

      integer :: i, n, ii, ll
      real :: summ

      n = size(lu, 1)
#ifdef PLASMA_DENSE
      call serial_getrs_d(lu, idx, n, 1, b, n)
#else
      ii = 0
      do i = 1, n
         ll = idx(i)
         summ = b(ll)
         b(ll) = b(i)
         if (ii /= 0) then
            summ = summ - dot_product(lu(i, ii:i - 1), b(ii:i - 1))
         else if (summ /= 0.0) then
            ii = i
         end if
         b(i) = summ
      end do
      do i = n, 1, -1
         b(i) = (b(i) - dot_product(lu(i, i + 1:n), b(i + 1:n))) / lu(i, i)
      end do
#endif

   end subroutine lu_back_substitution_real

#ifdef PLASMA_DENSE
   !****************************************************************************
   !               Back substitution with the serial PLASMA kernels
   !****************************************************************************
   ! The real version of <serial_getrs_z>, see there.
   !****************************************************************************
   subroutine serial_getrs_d(lu, idx, n, nrhs, b, ldb)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_core_dtrsm
      use plasma_interface, only: PlasmaLeft, PlasmaLower, PlasmaUpper
      use plasma_interface, only: PlasmaNoTrans, PlasmaUnit, PlasmaNonUnit

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n, nrhs, ldb
      real, dimension(ldb, nrhs), intent(in out) :: b

      real, dimension(nrhs) :: dum
      integer(c_int) :: m, k, lda, ldbc
      integer :: i

      !-------------------------------------------------------------------------

      do i = 1, n
         if (idx(i) /= i) then
            dum = b(idx(i), :)
            b(idx(i), :) = b(i, :)
            b(i, :) = dum
         end if
      end do

      m = int(n, c_int)
      k = int(nrhs, c_int)
      lda = int(size(lu, 1), c_int)
      ldbc = int(ldb, c_int)
      call plasma_core_dtrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit, &
                             m, k, 1.0, lu, lda, b, ldbc)
      call plasma_core_dtrsm(PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit, &
                             m, k, 1.0, lu, lda, b, ldbc)

   end subroutine serial_getrs_d
#endif

   subroutine lu_back_substitution_real_complex(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:), intent(in out) :: b

      integer :: i, n, ii, ll
      complex :: summ
#ifdef PLASMA_DENSE
      real, dimension(size(lu, 1), 2) :: br
#endif

    !! The dot products we use below automatically take the complex conjugate of its
    !! first argument, which we do not want here. Hence the use of conjg to undo this

      n = size(lu, 1)
#ifdef PLASMA_DENSE
      ! Solve for the real and imaginary parts as two right-hand sides
      br(:, 1) = real(b(1:n))
      br(:, 2) = aimag(b(1:n))
      call serial_getrs_d(lu, idx, n, 2, br, n)
      b(1:n) = cmplx(br(:, 1), br(:, 2))
#else
      ii = 0
      do i = 1, n
         ll = idx(i)
         summ = b(ll)
         b(ll) = b(i)
         if (ii /= 0) then
            summ = summ - dot_product(lu(i, ii:i - 1), b(ii:i - 1))
         else if (summ /= 0.0) then
            ii = i
         end if
         b(i) = summ
      end do
      do i = n, 1, -1
         b(i) = (b(i) - dot_product(lu(i, i + 1:n), b(i + 1:n))) / lu(i, i)
      end do
#endif

   end subroutine lu_back_substitution_real_complex

   subroutine lu_back_substitution_complex(lu, idx, b)

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:), intent(in out) :: b

      integer :: i, n, ii, ll
      complex :: summ

    !! The dot products we use below automatically take the complex conjugate of its
    !! first argument, which we do not want here. Hence the use of conjg to undo this

      n = size(lu, 1)
#ifdef PLASMA_DENSE
      call serial_getrs_z(lu, idx, n, 1, b, n)
#else
      ii = 0
      do i = 1, n
         ll = idx(i)
         summ = b(ll)
         b(ll) = b(i)
         if (ii /= 0) then
            summ = summ - dot_product(conjg(lu(i, ii:i - 1)), b(ii:i - 1))
         else if (summ /= 0.0) then
            ii = i
         end if
         b(i) = summ
      end do
      do i = n, 1, -1
         b(i) = (b(i) - dot_product(conjg(lu(i, i + 1:n)), b(i + 1:n))) / lu(i, i)
      end do
#endif

   end subroutine lu_back_substitution_complex

#ifdef PLASMA_DENSE
   !****************************************************************************
   !               Back substitution with the serial PLASMA kernels
   !****************************************************************************
   ! Solves A x = b for the <nrhs> right-hand sides in the first <n> rows of the
   ! columns of b (leading dimension <ldb>), given the LU factors of A from
   ! <lu_decomposition> (with or without PLASMA, the storage is the same): the
   ! row swaps in <idx>, then L y = b (unit diagonal) and U x = y.
   !
   ! plasma_core_ztrsm does not use the PLASMA context, so this can be called
   ! from inside the OpenMP parallel regions of stella.
   !****************************************************************************
   subroutine serial_getrs_z(lu, idx, n, nrhs, b, ldb)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_core_ztrsm
      use plasma_interface, only: PlasmaLeft, PlasmaLower, PlasmaUpper
      use plasma_interface, only: PlasmaNoTrans, PlasmaUnit, PlasmaNonUnit

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n, nrhs, ldb
      complex, dimension(ldb, nrhs), intent(in out) :: b

      complex, dimension(nrhs) :: dum
      integer(c_int) :: m, k, lda, ldbc
      integer :: i

      !-------------------------------------------------------------------------

      do i = 1, n
         if (idx(i) /= i) then
            dum = b(idx(i), :)
            b(idx(i), :) = b(i, :)
            b(i, :) = dum
         end if
      end do

      m = int(n, c_int)
      k = int(nrhs, c_int)
      lda = int(size(lu, 1), c_int)
      ldbc = int(ldb, c_int)
      call plasma_core_ztrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit, &
                             m, k, (1.0, 0.0), lu, lda, b, ldbc)
      call plasma_core_ztrsm(PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit, &
                             m, k, (1.0, 0.0), lu, lda, b, ldbc)

   end subroutine serial_getrs_z
#endif

   subroutine lu_back_substitution_matrix_real(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(in out) :: b

      integer :: i, j, n, ii, ll
      real :: summ

      n = size(lu, 1)
#ifdef PLASMA_DENSE
      call serial_getrs_d(lu, idx, n, size(b, 2), b, size(b, 1))
#else
      do j = 1, size(b, 2)
         ii = 0
         do i = 1, n
            ll = idx(i)
            summ = b(ll, j)
            b(ll, j) = b(i, j)
            if (ii /= 0) then
               summ = summ - dot_product(lu(i, ii:i - 1), b(ii:i - 1, j))
            else if (summ /= 0.0) then
               ii = i
            end if
            b(i, j) = summ
         end do
         do i = n, 1, -1
            b(i, j) = (b(i, j) - dot_product(lu(i, i + 1:n), b(i + 1:n, j))) / lu(i, i)
         end do
      end do
#endif

   end subroutine lu_back_substitution_matrix_real

   subroutine lu_back_substitution_matrix_complex(lu, idx, b)

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(in out) :: b

      integer :: i, j, n, ii, ll
      complex :: summ

    !! The dot products we use below automatically take the complex conjugate of its
    !! first argument, which we do not want here. Hence the use of conjg to undo this

      n = size(lu, 1)
#ifdef PLASMA_DENSE
      call serial_getrs_z(lu, idx, n, size(b, 2), b, size(b, 1))
#else
      do j = 1, size(b, 2)
         ii = 0
         do i = 1, n
            ll = idx(i)
            summ = b(ll, j)
            b(ll, j) = b(i, j)
            if (ii /= 0) then
               summ = summ - dot_product(conjg(lu(i, ii:i - 1)), b(ii:i - 1, j))
            else if (summ /= 0.0) then
               ii = i
            end if
            b(i, j) = summ
         end do
         do i = n, 1, -1
            b(i, j) = (b(i, j) - dot_product(conjg(lu(i, i + 1:n)), b(i + 1:n, j))) / lu(i, i)
         end do
      end do
#endif

   end subroutine lu_back_substitution_matrix_complex

  !!
  !! One shouldn't need to compute the inverse of a matrix if solving the linear equation A.x = y;
  !! a simple back substitution should suffice.
  !!
  !! Only compute the inverse of A when absolutely necessary!
  !!
   subroutine lu_inverse_real(lu, idx, inverse)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(out) :: inverse

      integer :: i, n

#ifdef PLASMA_DENSE
      if (plasma_lu_inverse_real(lu, idx, inverse)) return
#endif

      n = size(lu, 1)

      inverse = 0.0
      do i = 1, n
         inverse(i, i) = 1.0
      end do

      do i = 1, n
         call lu_back_substitution(lu, idx, inverse(:, i))
      end do

   end subroutine lu_inverse_real

#ifdef PLASMA_DENSE
   !****************************************************************************
   !                       Matrix inverse with plasma_dgetri
   !****************************************************************************
   ! The real version of <plasma_lu_inverse_complex>, see there.
   !****************************************************************************
   logical function plasma_lu_inverse_real(lu, idx, inverse) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_dgetri, stella_plasma_usable

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(out) :: inverse

      integer(c_int), dimension(size(lu, 1)) :: ipiv
      integer(c_int) :: n, info

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. stella_plasma_usable()) return

      n = int(size(lu, 1), c_int)
      inverse(1:n, 1:n) = lu(1:n, 1:n)
      ipiv = int(idx(1:n), c_int)
      info = plasma_dgetri(n, inverse, int(size(inverse, 1), c_int), ipiv)
      success = (info == 0)

   end function plasma_lu_inverse_real
#endif

   subroutine lu_inverse_complex(lu, idx, inverse)

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(out) :: inverse

      integer :: i, n

#ifdef PLASMA_DENSE
      if (plasma_lu_inverse_complex(lu, idx, inverse)) return
#endif

      n = size(lu, 1)

      inverse = 0.0
      do i = 1, n
         inverse(i, i) = 1.0
      end do

      do i = 1, n
         call lu_back_substitution(lu, idx, inverse(:, i))
      end do

   end subroutine lu_inverse_complex

#ifdef PLASMA_DENSE
   !****************************************************************************
   !                       Matrix inverse with plasma_zgetri
   !****************************************************************************
   ! Computes the inverse of A from its LU factors (from <lu_decomposition>, with
   ! or without PLASMA, the storage is the same as in LAPACK).
   !
   ! Returns .false. if PLASMA was not used, in which case <inverse> has to be
   ! computed by the stella routine:
   !  - if PLASMA can't be used here, see <stella_plasma_usable> in plasma_interface;
   !  - if U has a zero on the diagonal.
   !****************************************************************************
   logical function plasma_lu_inverse_complex(lu, idx, inverse) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: plasma_zgetri, stella_plasma_usable

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(out) :: inverse

      integer(c_int), dimension(size(lu, 1)) :: ipiv
      integer(c_int) :: n, info

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. stella_plasma_usable()) return

      n = int(size(lu, 1), c_int)
      inverse(1:n, 1:n) = lu(1:n, 1:n)
      ipiv = int(idx(1:n), c_int)
      info = plasma_zgetri(n, inverse, int(size(inverse, 1), c_int), ipiv)
      success = (info == 0)

   end function plasma_lu_inverse_complex
#endif

   function imaxloc(array)
      real, dimension(:), intent(in) :: array
      integer :: imaxloc
      integer, dimension(1) :: imax
      imax = maxloc(array)
      imaxloc = imax(1)
   end function imaxloc

end module linear_solve
