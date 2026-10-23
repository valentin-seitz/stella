!###############################################################################
!                      LU decomposition of dense matrices
!###############################################################################
! By default, the LU decomposition, back substitution and inverse are computed
! by the stella routines below, which use implicit (scaled) partial pivoting.
!
! When stella is compiled with PLASMA (see plasma_interface), and PLASMA is
! switched on for these routines (<use_plasma_dense>), they call PLASMA's
! plasma_?getrf, plasma_?getrs and plasma_?getri instead, or LAPACK's ?getrf,
! ?getrs and ?getri inside OpenMP parallel regions, where PLASMA can't be used
! (see <plasma_dense_usable>). PLASMA and LAPACK use partial pivoting, so the
! factors can differ from those of the stella routines, but they are stored in
! the same way: L (unit diagonal) below and U on and above the diagonal, and
! idx(j) is the row that was swapped with row j at step j. Hence the factors
! can be used by every <lu_back_substitution> and <lu_inverse>.
!
! The factors of a singular matrix are always computed by the stella routine,
! which replaces zero pivots by a small number, whereas PLASMA and LAPACK stop.
!###############################################################################
module linear_solve

   implicit none

   public :: lu_decomposition
   public :: lu_back_substitution
   public :: lu_inverse
   public :: imaxloc

   private

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

#ifdef PLASMA
      if (library_lu_decomposition_real(lu, idx, d)) return
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

   subroutine lu_decomposition_complex(lu, idx, d)

      implicit none

      complex, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      real, parameter :: zero = 1.0e-20
      real, dimension(size(lu, 1)) :: vv
      complex, dimension(size(lu, 2)) :: dum

      integer :: j, n, imax

#ifdef PLASMA
      if (library_lu_decomposition_complex(lu, idx, d)) return
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

   subroutine lu_back_substitution_real(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:), intent(in out) :: b

      integer :: i, n, ii, ll
      real :: summ

      n = size(lu, 1)

#ifdef PLASMA
      if (library_getrs_real(lu, idx, n, 1, b, n)) return
#endif

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

   end subroutine lu_back_substitution_real

   subroutine lu_back_substitution_real_complex(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:), intent(in out) :: b

      integer :: i, n, ii, ll
      complex :: summ
#ifdef PLASMA
      real, dimension(size(lu, 1), 2) :: br
#endif

    !! The dot products we use below automatically take the complex conjugate of its
    !! first argument, which we do not want here. Hence the use of conjg to undo this

      n = size(lu, 1)

#ifdef PLASMA
      ! Solve for the real and imaginary parts as two right-hand sides
      br(:, 1) = real(b(1:n))
      br(:, 2) = aimag(b(1:n))
      if (library_getrs_real(lu, idx, n, 2, br, n)) then
         b(1:n) = cmplx(br(:, 1), br(:, 2))
         return
      end if
#endif

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

#ifdef PLASMA
      if (library_getrs_complex(lu, idx, n, 1, b, n)) return
#endif

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

   end subroutine lu_back_substitution_complex

   subroutine lu_back_substitution_matrix_real(lu, idx, b)

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(in out) :: b

      integer :: i, j, n, ii, ll
      real :: summ

      n = size(lu, 1)

#ifdef PLASMA
      if (library_getrs_real(lu, idx, n, size(b, 2), b, size(b, 1))) return
#endif

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

#ifdef PLASMA
      if (library_getrs_complex(lu, idx, n, size(b, 2), b, size(b, 1))) return
#endif

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

#ifdef PLASMA
      if (library_lu_inverse_real(lu, idx, inverse)) return
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

   subroutine lu_inverse_complex(lu, idx, inverse)

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(out) :: inverse

      integer :: i, n

#ifdef PLASMA
      if (library_lu_inverse_complex(lu, idx, inverse)) return
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

   function imaxloc(array)
      real, dimension(:), intent(in) :: array
      integer :: imaxloc
      integer, dimension(1) :: imax
      imax = maxloc(array)
      imaxloc = imax(1)
   end function imaxloc

#ifdef PLASMA
   !****************************************************************************
   !                  LU decomposition with PLASMA or LAPACK
   !****************************************************************************
   ! Factorise <lu> with plasma_?getrf, or with LAPACK's ?getrf inside OpenMP
   ! parallel regions. Returns .false. if neither was used, in which case <lu>
   ! is unchanged and has to be factorised by the stella routine:
   !  - if PLASMA is switched off for the dense routines (<use_plasma_dense>);
   !  - if the matrix is singular, see the top of this module.
   !****************************************************************************
   logical function library_lu_decomposition_real(lu, idx, d) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: use_plasma_dense, plasma_dense_usable, plasma_dgetrf

      implicit none

      real, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      real, dimension(size(lu, 1), size(lu, 2)) :: a
      integer(c_int), dimension(size(lu, 1)) :: ipiv
      integer :: n, info

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. use_plasma_dense) return

      n = size(lu, 1)
      a = lu
      if (plasma_dense_usable()) then
         info = plasma_dgetrf(int(n, c_int), int(n, c_int), lu, int(n, c_int), ipiv)
      else
         call dgetrf(n, n, lu, n, ipiv, info)
      end if
      if (info /= 0) then
         lu = a
         return
      end if

      idx(1:n) = ipiv
      d = permutation_sign(idx, n)
      success = .true.

   end function library_lu_decomposition_real

   logical function library_lu_decomposition_complex(lu, idx, d) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: use_plasma_dense, plasma_dense_usable, plasma_zgetrf

      implicit none

      complex, dimension(:, :), intent(in out) :: lu
      integer, dimension(:), intent(out) :: idx
      real, intent(out) :: d

      complex, dimension(size(lu, 1), size(lu, 2)) :: a
      integer(c_int), dimension(size(lu, 1)) :: ipiv
      integer :: n, info

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. use_plasma_dense) return

      n = size(lu, 1)
      a = lu
      if (plasma_dense_usable()) then
         info = plasma_zgetrf(int(n, c_int), int(n, c_int), lu, int(n, c_int), ipiv)
      else
         call zgetrf(n, n, lu, n, ipiv, info)
      end if
      if (info /= 0) then
         lu = a
         return
      end if

      idx(1:n) = ipiv
      d = permutation_sign(idx, n)
      success = .true.

   end function library_lu_decomposition_complex

   ! The sign of the permutation of the LU decomposition, i.e., -1 for every row swap
   pure real function permutation_sign(idx, n)

      implicit none

      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n

      integer :: j

      permutation_sign = 1.0
      do j = 1, n
         if (idx(j) /= j) permutation_sign = -permutation_sign
      end do

   end function permutation_sign

   !****************************************************************************
   !                  Back substitution with PLASMA or LAPACK
   !****************************************************************************
   ! Solves A x = b for the <nrhs> right-hand sides in the first <n> rows of the
   ! columns of b (leading dimension <ldb>), given the LU factors of A from
   ! <lu_decomposition>, with plasma_?getrs, or with LAPACK's ?getrs inside
   ! OpenMP parallel regions. Returns .false. if PLASMA is switched off for the
   ! dense routines, in which case the stella routine has to solve it.
   !****************************************************************************
   logical function library_getrs_real(lu, idx, n, nrhs, b, ldb) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: use_plasma_dense, plasma_dense_usable, plasma_dgetrs, PlasmaNoTrans

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n, nrhs, ldb
      real, dimension(ldb, nrhs), intent(in out) :: b

      real, dimension(size(lu, 1), size(lu, 2)) :: a
      integer(c_int), dimension(n) :: ipiv
      integer :: info

      !-------------------------------------------------------------------------

      success = use_plasma_dense
      if (.not. success) return

      ipiv = int(idx(1:n), c_int)
      if (plasma_dense_usable()) then
         ! plasma_dgetrs does not change the factors, but has no const arguments
         a = lu
         info = plasma_dgetrs(PlasmaNoTrans, int(n, c_int), int(nrhs, c_int), a, int(size(lu, 1), c_int), &
                              ipiv, b, int(ldb, c_int))
      else
         call dgetrs('N', n, nrhs, lu, size(lu, 1), ipiv, b, ldb, info)
      end if

   end function library_getrs_real

   logical function library_getrs_complex(lu, idx, n, nrhs, b, ldb) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: use_plasma_dense, plasma_dense_usable, plasma_zgetrs, PlasmaNoTrans

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      integer, intent(in) :: n, nrhs, ldb
      complex, dimension(ldb, nrhs), intent(in out) :: b

      complex, dimension(size(lu, 1), size(lu, 2)) :: a
      integer(c_int), dimension(n) :: ipiv
      integer :: info

      !-------------------------------------------------------------------------

      success = use_plasma_dense
      if (.not. success) return

      ipiv = int(idx(1:n), c_int)
      if (plasma_dense_usable()) then
         ! plasma_zgetrs does not change the factors, but has no const arguments
         a = lu
         info = plasma_zgetrs(PlasmaNoTrans, int(n, c_int), int(nrhs, c_int), a, int(size(lu, 1), c_int), &
                              ipiv, b, int(ldb, c_int))
      else
         call zgetrs('N', n, nrhs, lu, size(lu, 1), ipiv, b, ldb, info)
      end if

   end function library_getrs_complex

   !****************************************************************************
   !                     Matrix inverse with PLASMA or LAPACK
   !****************************************************************************
   ! Computes the inverse of A from its LU factors (from <lu_decomposition>)
   ! with plasma_?getri, or with LAPACK's ?getri inside OpenMP parallel regions.
   ! Returns .false. if neither was used, in which case <inverse> has to be
   ! computed by the stella routine:
   !  - if PLASMA is switched off for the dense routines (<use_plasma_dense>);
   !  - if U has a zero on the diagonal.
   !****************************************************************************
   logical function library_lu_inverse_real(lu, idx, inverse) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: use_plasma_dense, plasma_dense_usable, plasma_dgetri

      implicit none

      real, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      real, dimension(:, :), intent(out) :: inverse

      integer(c_int), dimension(size(lu, 1)) :: ipiv
      real, dimension(size(lu, 1)) :: work
      integer :: n, info

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. use_plasma_dense) return

      n = size(lu, 1)
      inverse(1:n, 1:n) = lu(1:n, 1:n)
      ipiv = int(idx(1:n), c_int)
      if (plasma_dense_usable()) then
         info = plasma_dgetri(int(n, c_int), inverse, int(size(inverse, 1), c_int), ipiv)
      else
         call dgetri(n, inverse, size(inverse, 1), ipiv, work, n, info)
      end if
      success = (info == 0)

   end function library_lu_inverse_real

   logical function library_lu_inverse_complex(lu, idx, inverse) result(success)

      use, intrinsic :: iso_c_binding, only: c_int
      use plasma_interface, only: use_plasma_dense, plasma_dense_usable, plasma_zgetri

      implicit none

      complex, dimension(:, :), intent(in) :: lu
      integer, dimension(:), intent(in) :: idx
      complex, dimension(:, :), intent(out) :: inverse

      integer(c_int), dimension(size(lu, 1)) :: ipiv
      complex, dimension(size(lu, 1)) :: work
      integer :: n, info

      !-------------------------------------------------------------------------

      success = .false.
      if (.not. use_plasma_dense) return

      n = size(lu, 1)
      inverse(1:n, 1:n) = lu(1:n, 1:n)
      ipiv = int(idx(1:n), c_int)
      if (plasma_dense_usable()) then
         info = plasma_zgetri(int(n, c_int), inverse, int(size(inverse, 1), c_int), ipiv)
      else
         call zgetri(n, inverse, size(inverse, 1), ipiv, work, n, info)
      end if
      success = (info == 0)

   end function library_lu_inverse_complex
#endif

end module linear_solve
