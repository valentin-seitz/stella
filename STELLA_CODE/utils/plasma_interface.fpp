!###############################################################################
!                         Interface to the PLASMA library
!###############################################################################
! PLASMA (Parallel Linear Algebra Software for Multicore Architectures) is used
! for the LU operations of stella when it is compiled with USE_PLASMA=on (Make)
! or -DSTELLA_ENABLE_PLASMA=ON (CMake), which defines the PLASMA preprocessor flag.
!
! PLASMA does not ship a Fortran module, so the C routines that stella needs
! are bound here through iso_c_binding. The signatures follow
!     plasma/include/plasma_context.h
!     plasma/include/plasma_z.h, plasma_d.h
!     plasma/include/plasma_core_blas_z.h, plasma_core_blas_d.h
! and the constants follow plasma/include/plasma_types.h.
!
! The routines <stella_plasma_init> and <stella_plasma_finish> always exist, so
! that they can be called from files that are not preprocessed. Without the
! PLASMA flag they do nothing.
!###############################################################################
module plasma_interface

#ifdef PLASMA
   use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex
#endif

   implicit none

   public :: stella_plasma_init
   public :: stella_plasma_finish

#ifdef PLASMA

   private

   ! Constants from plasma/include/plasma_types.h
   integer(c_int), parameter, public :: PlasmaSuccess = 0
   integer(c_int), parameter, public :: PlasmaNoTrans = 111
   integer(c_int), parameter, public :: PlasmaTrans = 112
   integer(c_int), parameter, public :: PlasmaUpper = 121
   integer(c_int), parameter, public :: PlasmaLower = 122
   integer(c_int), parameter, public :: PlasmaNonUnit = 131
   integer(c_int), parameter, public :: PlasmaUnit = 132
   integer(c_int), parameter, public :: PlasmaLeft = 141
   integer(c_int), parameter, public :: PlasmaRight = 142
   integer(c_int), parameter, public :: PlasmaDisabled = 0
   integer(c_int), parameter, public :: PlasmaEnabled = 1
   integer(c_int), parameter, public :: PlasmaTuning = 0
   integer(c_int), parameter, public :: PlasmaNb = 1
   integer(c_int), parameter, public :: PlasmaIb = 2

   ! Whether <stella_plasma_init> has been called, the PLASMA routines that use
   ! the PLASMA context can only be called after it
   logical, public, protected :: stella_plasma_initialised = .false.

   public :: plasma_init, plasma_finalize, plasma_set, plasma_get
   public :: plasma_zgetrf, plasma_zgetrs, plasma_zgetri
   public :: plasma_dgetrf, plasma_dgetrs, plasma_dgetri
   public :: plasma_zgbtrf, plasma_zgbtrs
   public :: plasma_core_ztrsm, plasma_core_dtrsm

   interface

      !-------------------------------------------------------------------------
      !                           Context (plasma_context.h)
      !-------------------------------------------------------------------------
      integer(c_int) function plasma_init() bind(c, name='plasma_init')
         import :: c_int
      end function plasma_init

      integer(c_int) function plasma_finalize() bind(c, name='plasma_finalize')
         import :: c_int
      end function plasma_finalize

      integer(c_int) function plasma_set(param, val) bind(c, name='plasma_set')
         import :: c_int
         integer(c_int), value :: param, val
      end function plasma_set

      integer(c_int) function plasma_get(param, val) bind(c, name='plasma_get')
         import :: c_int
         integer(c_int), value :: param
         integer(c_int), intent(out) :: val
      end function plasma_get

      !-------------------------------------------------------------------------
      !                     Dense LU, complex (plasma_z.h)
      !-------------------------------------------------------------------------
      integer(c_int) function plasma_zgetrf(m, n, pA, lda, ipiv) bind(c, name='plasma_zgetrf')
         import :: c_int, c_double_complex
         integer(c_int), value :: m, n, lda
         complex(c_double_complex), dimension(*) :: pA
         integer(c_int), dimension(*) :: ipiv
      end function plasma_zgetrf

      integer(c_int) function plasma_zgetrs(trans, n, nrhs, pA, lda, ipiv, pB, ldb) bind(c, name='plasma_zgetrs')
         import :: c_int, c_double_complex
         integer(c_int), value :: trans, n, nrhs, lda, ldb
         complex(c_double_complex), dimension(*) :: pA, pB
         integer(c_int), dimension(*) :: ipiv
      end function plasma_zgetrs

      integer(c_int) function plasma_zgetri(n, pA, lda, ipiv) bind(c, name='plasma_zgetri')
         import :: c_int, c_double_complex
         integer(c_int), value :: n, lda
         complex(c_double_complex), dimension(*) :: pA
         integer(c_int), dimension(*) :: ipiv
      end function plasma_zgetri

      !-------------------------------------------------------------------------
      !                       Dense LU, real (plasma_d.h)
      !-------------------------------------------------------------------------
      integer(c_int) function plasma_dgetrf(m, n, pA, lda, ipiv) bind(c, name='plasma_dgetrf')
         import :: c_int, c_double
         integer(c_int), value :: m, n, lda
         real(c_double), dimension(*) :: pA
         integer(c_int), dimension(*) :: ipiv
      end function plasma_dgetrf

      integer(c_int) function plasma_dgetrs(trans, n, nrhs, pA, lda, ipiv, pB, ldb) bind(c, name='plasma_dgetrs')
         import :: c_int, c_double
         integer(c_int), value :: trans, n, nrhs, lda, ldb
         real(c_double), dimension(*) :: pA, pB
         integer(c_int), dimension(*) :: ipiv
      end function plasma_dgetrs

      integer(c_int) function plasma_dgetri(n, pA, lda, ipiv) bind(c, name='plasma_dgetri')
         import :: c_int, c_double
         integer(c_int), value :: n, lda
         real(c_double), dimension(*) :: pA
         integer(c_int), dimension(*) :: ipiv
      end function plasma_dgetri

      !-------------------------------------------------------------------------
      !                     Banded LU, complex (plasma_z.h)
      !-------------------------------------------------------------------------
      integer(c_int) function plasma_zgbtrf(m, n, kl, ku, pAB, ldab, ipiv) bind(c, name='plasma_zgbtrf')
         import :: c_int, c_double_complex
         integer(c_int), value :: m, n, kl, ku, ldab
         complex(c_double_complex), dimension(*) :: pAB
         integer(c_int), dimension(*) :: ipiv
      end function plasma_zgbtrf

      integer(c_int) function plasma_zgbtrs(trans, n, kl, ku, nrhs, pAB, ldab, ipiv, pB, ldb) bind(c, name='plasma_zgbtrs')
         import :: c_int, c_double_complex
         integer(c_int), value :: trans, n, kl, ku, nrhs, ldab, ldb
         complex(c_double_complex), dimension(*) :: pAB, pB
         integer(c_int), dimension(*) :: ipiv
      end function plasma_zgbtrs

      !-------------------------------------------------------------------------
      ! Serial core kernels (plasma_core_blas_{z,d}.h). These operate directly on
      ! LAPACK-layout arrays, do not use the PLASMA context and are thread safe,
      ! so they can be called from inside the OpenMP regions of stella.
      !-------------------------------------------------------------------------
      subroutine plasma_core_ztrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) &
         bind(c, name='plasma_core_ztrsm')
         import :: c_int, c_double_complex
         integer(c_int), value :: side, uplo, transa, diag, m, n, lda, ldb
         complex(c_double_complex), value :: alpha
         complex(c_double_complex), dimension(*), intent(in) :: A
         complex(c_double_complex), dimension(*) :: B
      end subroutine plasma_core_ztrsm

      subroutine plasma_core_dtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) &
         bind(c, name='plasma_core_dtrsm')
         import :: c_int, c_double
         integer(c_int), value :: side, uplo, transa, diag, m, n, lda, ldb
         real(c_double), value :: alpha
         real(c_double), dimension(*), intent(in) :: A
         real(c_double), dimension(*) :: B
      end subroutine plasma_core_dtrsm

   end interface

#endif

contains

   !****************************************************************************
   !                            Initialise PLASMA
   !****************************************************************************
   ! PLASMA has to be initialised once per MPI process, outside of any OpenMP
   ! parallel region. Its internal tuning is switched off, so that the shared
   ! PLASMA context is never modified while stella calls PLASMA concurrently.
   !****************************************************************************
   subroutine stella_plasma_init

#ifdef PLASMA
      use mp, only: mp_abort
#endif

      implicit none

#ifdef PLASMA
      integer(c_int) :: info

      !-------------------------------------------------------------------------

      info = plasma_init()
      if (info /= PlasmaSuccess) call mp_abort('stella_plasma_init: plasma_init failed')
      info = plasma_set(PlasmaTuning, PlasmaDisabled)
      if (info /= PlasmaSuccess) call mp_abort('stella_plasma_init: disabling PLASMA tuning failed')
      stella_plasma_initialised = .true.
#endif

   end subroutine stella_plasma_init

   !****************************************************************************
   !                             Finalise PLASMA
   !****************************************************************************
   subroutine stella_plasma_finish

      implicit none

#ifdef PLASMA
      integer(c_int) :: info

      !-------------------------------------------------------------------------

      info = plasma_finalize()
      stella_plasma_initialised = .false.
#endif

   end subroutine stella_plasma_finish

end module plasma_interface
