!###############################################################################
!############################### RESPONSE MATRIX ###############################
!###############################################################################
! 
! This module computes the response matrix for inverting the parallel streaming
! operator when parallel_streaming is treated implicitly.
! 
! Method:
! -------
! (Note: to see a certain part of the code, search the step it is listed under, for
! example to see where the homogeneous parallel streaming equation is solved, search 1.B)
! 
! Step 1) Get the matrix that we need to invert in implicit parallel streaming.
!         To do this we take the following steps:
! 
!         1.A) First initialise an unit impulse for each field and for each grid point in
!              the extended zed domain for a given ky and set of connected kx modes.
! 
!         1.B) Next, solve the homogeneous parallel streaming equation for the
!              distribution function, given a field with a unit impulse. This is
!              done for each field (i.e. phi, apar, bpar) that is included.
! 
!         1.C) Next, solve the field equations given the distribution function response.
! 
!              1.C.I)  Integrate over velocities to get the appropriate field operator
!                      that acts on the distibution function:
!                      For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
!                      For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
!                      For bpar:  - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g
! 
!              1.C.II) Divide by the appropriate pre-factor in the field equations.
! 
!         1.D) Construct the matrix that we need to invert by including the 
!              identity matrix where appropriate.
! 
! Step 2) LU decompose the matrix.
! 
!###############################################################################
module response_matrix

   ! Load libraries
   use netcdf

   ! Load debug flags
   use debug_flags, only: debug => response_matrix_debug

   implicit none

   ! Make routines available to other modules
   public :: init_response_matrix
   public :: finish_response_matrix
   public :: read_response_matrix
   
   ! Only initialise once
   public :: initialised_response_matrix

   private

   ! Local variables
   
   ! Only initialise once
   logical :: initialised_response_matrix = .false.

contains

!###############################################################################
!########################## INITIALISE RESPONSE MATRIX #########################
!###############################################################################

   !============================================================================
   !            Main routine called to initialise the response matrix           
   !============================================================================
   subroutine init_response_matrix

      use mp, only: proc0
      use debug_flags, only: print_extra_info_to_terminal
      use parallelisation_layouts, only: verbose
      use linear_solve, only: lu_decomposition
      use parallelisation_layouts, only: iv_idx, is_idx
      use parallelisation_layouts, only: mat_gen
      use file_units, only: unit_response_matrix
      use arrays, only: response_matrix
      use grids_kxky, only: naky

      implicit none

      ! ------------------------------------------------------------------------
      !              Set up memory and timings for response matrix              
      ! ------------------------------------------------------------------------
      
      ! Print info to command prompt
      if (verbose .and. print_extra_info_to_terminal) then
         write (*, '(A)') '############################################################'
         write (*, '(A)') '                       RESPONSE MATRIX                      '
         write (*, '(A)') '############################################################'
         write (*, '(A)') ''
      end if
      
      ! Set up response matrix utils
      call setup_response_matrix_file_io

      ! Only initialise once
      if (initialised_response_matrix) return
      initialised_response_matrix = .true.

      ! Allocate response matrix
      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))
      
      ! ------------------------------------------------------------------------
      !                      Construct the response matrix                      
      ! ------------------------------------------------------------------------
      
      ! This is the main routine for computing the response matrix, where all
      ! the calculations are done.
      call construct_response_matrix

      ! ------------------------------------------------------------------------
      !                         Close the initialisation                        
      ! ------------------------------------------------------------------------
      
      ! Write the response matrix to an output file if <mat_gen> = True
      if (proc0 .and. mat_gen) then
         close (unit=unit_response_matrix)
      end if
      
   end subroutine init_response_matrix
   
   !============================================================================
   !                     Set up .io file for response matrix                    
   !============================================================================
   subroutine setup_response_matrix_file_io

      use mp, only: proc0, job
      use parallelisation_layouts, only: mat_gen
      use system_fortran, only: systemf
      use grids_kxky, only: naky
      use file_units, only: unit_response_matrix

      implicit none

      ! Local variables
      character(len=15) :: job_str
      character(len=100) :: file_name

      ! ------------------------------------------------------------------------

      ! All matrices handled by processor i_proc and job are stored
      ! on a single file named: response_mat_job.iproc
      if (proc0 .and. mat_gen) then
         call systemf('mkdir -p mat')
         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat_'//trim(job_str)
         open (unit=unit_response_matrix, status='replace', file=file_name, &
               position='rewind', action='write', form='unformatted')
         write (unit=unit_response_matrix) naky
      end if

   end subroutine setup_response_matrix_file_io

!===============================================================================
!===============================================================================
!========================== CONSTRUCT RESPONSE MATRIX ==========================
!===============================================================================
!===============================================================================

   !============================================================================
   !               Main routine for constructing response martrix
   !============================================================================
   ! 
   ! Background:
   ! -----------
   ! Different ky modes are independent, so this creates sets of connected kx values:
   ! 
   ! After one full 2π orbit in zed, the eddy becomes sheared, shifting the mode
   ! to a higher kx value. To represent this effect, we apply the twist-and-shift
   ! boundary conditions, which connect these shifted kx values along an extended
   ! z-grid. For a given ky, the connected kx modes are spaced by,
   !                 δkx = 2π p ∂ι/∂r * (∂y/∂α) * (∂ψ/∂x) * ky
   ! where ∂ι/∂r is the radial derivative of the rotational transform.
   ! 
   ! This construction implies that not all kx modes are mutually connected. It
   ! is clear that for larger ky-values, the step size δkx becomes larger.
   ! Therefore, for each ky, we form distinct chains of connected modes, with the
   ! number of such chains determined by the spacing δkx.
   ! 
   ! The variable <neigen> gives the number of distinct chains of "eigenmodes"
   ! for a given ky, if the twist-and-shift boundary conditions are employed.
   ! If, instead, we use periodic boundary conditions, then all modes are 
   ! unconnected, and we have <neigen> = nakx for each ky-value.
   ! 
   ! What is done:
   ! ------------- 
   ! In this subroutine we construct the response matrix. This is split into two 
   ! main contributions: 
   !     1) Getting the velocity-space-integrated response. To do this we supply
   !        a unit impulse at every point in the extended zed domain for a given ky,
   !        and for a given set of connected kx values. With this, we then solve the
   !        parallel streaming equations to find the response of the distribution
   !        function to this unit impulse. This is done for each field that is
   !        included (phi, apar, bpar). Once this distribution function is obtained
   !        we integrate it with appropriate weights over (vpa, mu) and sum over
   !        species. This gives the quasineutrality component that acts on g.
   !        Then, we apply the field solve to find the matrix that we need to invert.
   !        This is essentially dividing by the correct factor from the field equations.
   !     2) LU decompose the matrix.
   !============================================================================
   subroutine construct_response_matrix

      use mp, only: proc0
      use parallelisation_layouts, only: mat_gen
      use arrays, only: response_matrix
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use file_units, only: unit_response_matrix

      implicit none

      ! Local variables
      integer :: iky, ie

      ! ------------------------------------------------------------------------
      
      ! Loop over all ky modes, as these are all independent of one another.
      do iky = 1, naky

         ! Write the ky-value to the output file
         if (proc0 .and. mat_gen) then
            write (unit=unit_response_matrix) iky, neigen(iky)
         end if

         ! ---------------------------------------------------------------------
         !              STEP 1) Find the response matrix to invert              
         ! ---------------------------------------------------------------------
         
         ! For a given ky, we need to associate an %eigen to it - to denote the
         ! connected modes. The response matrix for each ky has neigen(ky).
         if (.not. associated(response_matrix(iky)%eigen)) allocate (response_matrix(iky)%eigen(neigen(iky)))

         !> TO VALENTIN: This part could be parallelised over iky and ie as these 
         !> are all computed independently. It is not until the LU decomposition
         !> that they are needed together

         ! Loop over the independent chains of modes for a given ky value. Note
         ! that <neigen> is an integer that depends on ky, as different ky modes
         ! will have a different number of eigenmode chains. For twist-and-shift,
         ! the larger ky is, the larger δkx is, and the more independent chains
         ! there are, while all kx-modes are typically connected for the smallest ky.
         ! Hence neigen(ky_min) is typically 1 while neigen(ky_max) is typically nakx.
         do ie = 1, neigen(iky)
            call calculate_response_matrix_to_invert(iky, ie)
         end do 
         
         ! (DSO) -- Here, the parallelisation over velocity space ends. At this
         ! point every processor has the response matrix for a given ky, and the
         ! LU decomposition could be parallelised locally if need be.
         ! This is preferable to parallelisation over ky as the LU
         ! decomposition (and perhaps QN) will be dominated by the
         ! ky with the most connections if twist-and-shift is used.

         ! ---------------------------------------------------------------------
         !                   STEP 2) LU decompose the matrix 
         ! ---------------------------------------------------------------------

         ! LU decompose the matrix
         call lu_decompose_response_matrix(iky)

         ! ---------------------------------------------------------------------
         !                  Write the response matrix to a file                 
         ! ---------------------------------------------------------------------

         ! Write response matrix to file
         do ie = 1, neigen(iky)
            if (proc0 .and. mat_gen) then
               write (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%idx
               write (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%zloc
            end if
         end do

      end do

   end subroutine construct_response_matrix

   !============================================================================
   !                Step 1) Find distribution function response                 
   !============================================================================
   ! This subroutine computes the matrix we need to invert. The exact routines 
   ! called will depend on which fields are included in the simulation
   ! (phi, apar, bpar).
   !============================================================================
   subroutine calculate_response_matrix_to_invert(iky, ie)

      use mp, only: proc0
      use timers, only: region_start, region_end
      use parallelisation_layouts, only: mat_gen
      use parameters_physics, only: include_apar, include_bpar
      use grids_extended_zgrid, only: ikxmod
      use grids_extended_zgrid, only: nsegments, nzed_segment
      use grids_extended_zgrid, only: periodic
      use grids_extended_zgrid, only: iz_low, iz_up
      use parallelisation_layouts, only: vmu_lo
      use field_equations, only: nfields
      use file_units, only: unit_response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie

      ! Local variables
      integer :: idx, ikx, iseg
      integer :: iz, izl_offset, izup
      integer :: nz_ext, nresponse, nresponse_per_field
      complex, dimension(:, :), allocatable :: gext
      complex, dimension(:), allocatable :: phi_ext, apar_ext, bpar_ext
      
      !-------------------------------------------------------------------------
      
      ! Start the timer
      call region_start('response_matrix')

      ! ------------------------------------------------------------------------
      !              Set up system to compute response of the pdf               
      ! ------------------------------------------------------------------------
      ! <nz_ext> is the dimension of the extended zed domain, <nzed_segment> is
      ! the number of of unique zed values in all segments except the first. The 
      ! first segment has one extra unique zed value (all others have one grid 
      ! point in common with the previous segment due to periodicity). Finally,
      ! <nsegments> is the number of segments on the extended zed domain, in 
      ! other words, it is the number of 2π segments in z that are linked together.
      !     <nz_ext> = (number of zeds) x (number of segments on extended zed domain)
      !     <nzed_segment> = (number of zeds)
      !     <nsegments> = (Nkx -1)/neigen = (number of segments on extended zed domain)
      ! ------------------------------------------------------------------------
      
      ! Calculate the number of z-points on the extended z-domain
      nz_ext = nsegments(ie, iky) * nzed_segment + 1

      ! We need to treat the zonal mode specially to avoid double counting the end
      ! points in zed, as it is periodic so these end points are the same.
      ! This is also needed if 'periodic' is chosen for the boundary conditions.
      ! At this point <nresponse> is the number of independednt zed values for
      ! each field. This is the number of unit impulses we need to apply in order
      ! to get the response matrix.
      if (periodic(iky)) then
         nresponse_per_field = nz_ext - 1
      else
         nresponse_per_field = nz_ext
      end if

      ! If electromagnetic, we need to consider the response of phi, apar, and bpar.
      ! If we have more fields we need to apply more unit impulses -> one for each 
      ! of the fields phi, apar, and bpar.
      nresponse = nresponse_per_field * nfields

      ! Write <ie> and <nresponse> to the output file
      if (proc0 .and. mat_gen) then
         write (unit=unit_response_matrix) ie, nresponse
      end if

      ! Allocate response_matrix%eigen%zloc = size of response matrix to invert
      ! Allocate response_matrix%idx = pivot index needed for LU decomposition
      call setup_response_matrix_zloc_idx(iky, ie, nresponse)

      ! Allocate arrays on the extended zed domain
      ! Fields only have 1 index, as we are looping over ky modes, and kx and 
      ! zed are connected via the extended zed domain, while the distribution 
      ! function, gext, has 2 dimensions: the extended zed dimension, and velocity.
      allocate (phi_ext(nz_ext)); phi_ext = 0.0
      allocate (apar_ext(nz_ext)); apar_ext = 0.0
      allocate (bpar_ext(nz_ext)); bpar_ext = 0.0
      allocate (gext(nz_ext, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); gext = 0.0

      ! ------------------------------------------------------------------------
      !    Get the matrix we need to invert depending on the fields included    
      ! ------------------------------------------------------------------------
      
      ! <idx> here is the index in the extended zed domain that we are giving
      ! a unit impulse to.
      idx = 0

      ! Loop over segments, and idex this with <iseg>. 
      ! The first segment is special because it has one more unique
      ! zed value than all the other segments since the domain
      ! is [z0-pi:z0+pi], and we are including both endpoints.
      ! For iseg > 1 one endpoint is shared with the previous segment.
      izl_offset = 0

      ! Here, we apply a unit impulse at every value of zed in this segment, and find the 
      ! response of gext to this unit impulse. The impulse is provided at the <idx> value.
      ! Note - No need to obtain response to impulses at negative kx values
      ! Loop over all segments
      do iseg = 1, nsegments(ie, iky)
      
         ! Compute the index of kx that is connected in the given eigen chain in this segment.
         ! <ikxmod> gives the kx corresponding to (iseg,ie,iky)
         ! i.e. given a ky (iky), which chain are we in (ie), and within that chain
         !      which segment of 2π are we considering -> this is associated with a specific
         !      kx value. 
         ikx = ikxmod(iseg, ie, iky)

         ! Make sure the boundary points are being treated correctly depending
         ! on whether the mode is periodic or not. Here, define <izup> as the 
         ! upper zed value within a segment. If the mode is periodic, then 
         ! reduce the upper bound by one, as this is a repeated point so it is 
         ! obtained using the periodicity condition. This avoids and double-counting.
         ! Only the last segment ends on the periodic point; the upper end point of the other
         ! segments is shared with the next segment, which skips it through <izl_offset>
         if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
            izup = iz_up(iseg) - 1
         else
            izup = iz_up(iseg)
         end if
         
         ! Now apply a unit impulse at each zed location. To do this loop over the zed index, but
         ! recall that there is one less zed grid point in these connected segments as they share
         ! a grid point with the previous segment. 
         ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
         ! iz is only cycling through the zed location within a given segment. 
         do iz = iz_low(iseg) + izl_offset, izup
            idx = idx + 1
            call get_response_matrix_for_phi(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
            if (include_apar) call get_response_matrix_for_apar(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
            if (include_bpar) call get_response_matrix_for_bpar(iky, ie, idx, nz_ext, nresponse_per_field, phi_ext, apar_ext, bpar_ext, gext)
         end do

         ! Set the offset to 1 - all other connected segments need to start one point
         ! displaced as they share a point with the previous segment. 
         if (izl_offset == 0) izl_offset = 1
         
      end do

      ! Deallocate temporary arrays
      deallocate (gext, phi_ext, apar_ext, bpar_ext)

      ! Stop the timer
      call region_end('response_matrix')

   end subroutine calculate_response_matrix_to_invert

   !============================================================================
   !                      Get response matrix column for phi
   !============================================================================
   ! Here, we apply a unit impulse to phi at a given point on the extended zed
   ! domain. This location is indicated via <idx>
   ! 
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in 
   !         the same eigen chain can communicate with eachother
   ! <idx> = location of unit impulse within the segment 
   ! <nz_ext>   = length of extended zed domain for the given eigen chain 
   ! <phi_ext>  = phi on the extended zed domain
   ! <apar_ext> = apar on the extended zed domain
   ! <bpar_ext> = bpar on the extended zed domain
   ! <pdf_ext>  = distribution function on the extended zed domain
   ! 
   ! What is being done: 
   ! -------------------
   ! 1.A) First initialise (real) unit impulse phi^{n+1} (or Delta phi^{n+1})
   !      at each grid point in the extended zed domain for a given connected chain.
   !      This unit impulse is supplied at the index <idx>.
   !      If using periodic boundary conditions/or treating the zonal mode, then
   !      we need to make sure phi is also periodic in zed - so apply phase shift 
   !      to opposite mode. This is only important when the unit impulse is applied 
   !      at the boundary, otherwise both end points are automatically zero. 
   ! 1.B) Then solve the homogeneous parallel streaming equation for the response of
   !      the distribution function given a field with a unit impulse.
   ! 1.C) Then integrate over velocities to get the appropriate field operator that acts
   !      on the distibution function:
   !           For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !           For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
   !           For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g 
   !      Then we need to divide by the appropriate prefactor from the fields equations
   ! 1.D) Take into account identity matricies that appear in the 
   !      response equation. (see stella manual - to come!)
   !============================================================================
   subroutine get_response_matrix_for_phi(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)
   

      use parallelisation_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use parameters_physics, only: include_apar, include_bpar
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: periodic, phase_shift

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      ! Local variables
      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar

      ! ------------------------------------------------------------------------
      !                    1.A) Initialise a unit impulse in phi
      ! ------------------------------------------------------------------------
      ! Provide a unit impulse to phi^{n+1} (or Delta phi^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[phi], into L[f1] = R[phir] and L[f2] = i*R[phii],
      ! with f = f1 + f2; then phi = df1/dphir * phir + df2/dphii * phii.
      ! however, we see that if phir = phii = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives phi = df1/dphir * (phir + i * phii) = df1/dphir * phi
      ! ------------------------------------------------------------------------
      
      ! Initialise
      phi_ext = 0.0
      apar_ext = 0.0
      bpar_ext = 0.0
      
      ! How phi^{n+1} enters the GKE depends on whether we are solving for the
      ! non-Boltzmann pdf, h, or the guiding center pdf, 'g'
      phi_ext(idx) = time_upwind_plus
      
      ! Need to make sure that if the mode is periodic, then the boundaries match up to
      ! a phase factor. In practice this only matters if the unit impulse is at the 
      ! boundary (i.e. <idx ==1 ) otherwise phi = 0.0 at both boundary points anyway. 
      ! TOGO-GA: check division rather than multiplication -- kept division for now to 
      ! be consistent with parallel_streaming phase shift 
      if (periodic(iky) .and. idx == 1) phi_ext(nz_ext) = phi_ext(1) / phase_shift(iky)

      ! <dum> is a scratch array that takes the place of the pdf and phi
      ! at the previous time level. It also replaces the other fields (apar and bpar).
      ! This is set to zero for the response matrix approach because there is no sources
      ! coming from the previous time step in the response matrix equation.
      allocate (dum(nz_ext)); dum = 0.0

      ! ------------------------------------------------------------------------
      !      1.B) Get distribution function response to unit impulse in phi
      ! ------------------------------------------------------------------------

      ! Set the flux tube index to one - eed to check, but think this is okay as 
      ! the homogeneous equation solved here for the response matrix construction 
      ! is the same for all flux tubes in the flux tube train
      it = 1

      ! Solve for the response of the distribution function given this unit impulse in <phi>
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Calculate the RHS of the GK equation (using dum=0 as the pdf at the previous 
         ! time level, and phi_ext as the input potential) and store it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, phi_ext, dum, dum, dum, dum, pdf_ext(:, ivmu))
         
         ! Given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain; the RHS is
         ! input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
         
      end do

      deallocate (dum)

      ! ------------------------------------------------------------------------
      !          1.C) Integrate over velocity to get the matrix to invert
      ! ------------------------------------------------------------------------
      ! We now have the pdf on the extended zed domain at this ky and set of 
      ! connected kx values corresponding to a unit impulse in phi at this 
      ! location now integrate over velocities to get a square response matrix.
      ! 
      ! solve_field_equations_using_pdf_response is the operator that acts on the pdf in 
      ! quasineutrality. e.g. if using g for the pdf: 
      ! 
      !        sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]
      ! 
      ! and divide by the appropriate factors to get the field. 
      ! 
      ! This ends the parallelisation over velocity space, so every core 
      ! should have a copy of phi_ext.
      ! ------------------------------------------------------------------------
      call solve_field_equations_using_pdf_response(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie, nz_ext)

      ! ------------------------------------------------------------------------
      !                      1.D) Compute the matrix to invert                
      ! ------------------------------------------------------------------------
      ! The response matrix is the output from <get_fields_for_response_matrix>.
      ! We then just need to compute the matrix that is to be inverted in parallel streaming. 
      ! ------------------------------------------------------------------------

      offset_apar = 0
      if (include_apar) offset_apar = nresponse
      if (include_bpar) offset_bpar = offset_apar + nresponse

      ! Next need to create column in response matrix from phi_ext, apar_ext and bpar_ext
      ! The negative sign occurs because the matrix acts on the RHS of the streaming equation.
      ! For the location where the unit impulse is applied, the matrix to invert is:
      !                       (identity matrix - response matrix)
      ! For all other locations the matrix to invert is:
      !                             (- response_matrix)
      ! So add in contribution from identity matrix for the <idx> location for the field we are solving for:
      phi_ext(idx) = phi_ext(idx) - 1.0
      
      ! But everywhere else, simply add a negative sign:
      response_matrix(iky)%eigen(ie)%zloc(:nresponse, idx) = -phi_ext(:nresponse)
      if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, idx) = -apar_ext(:nresponse)
      if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, idx) = -bpar_ext(:nresponse)
      
      ! ------------------------------------------------------------------------

   end subroutine get_response_matrix_for_phi

   !============================================================================
   !                       Get response matrix for apar
   !============================================================================
   ! Here, we apply a unit impulse to apar at a given point on the extended zed
   ! domain. This location is indicated via <idx>
   ! 
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in 
   !         the same eigen chain can communicate with eachother
   ! <idx> = location of unit impulse within the segment 
   ! <nz_ext>   = length of extended zed domain for the given eigen chain 
   ! <phi_ext>  = phi on the extended zed domain
   ! <apar_ext> = apar on the extended zed domain
   ! <bpar_ext> = bpar on the extended zed domain
   ! <pdf_ext>  = distribution function on the extended zed domain
   ! 
   ! What is being done: 
   ! -------------------
   ! 1.A) First initialise (real) unit impulse apar^{n+1} at each grid point in the
   !      extended zed domain for a given connected chain. This unit impulse is 
   !      supplied at the index <idx>. If using periodic boundary conditions/or
   !      treating the zonal mode, then we need to make sure phi is also periodic
   !      in zed - so apply phase shift to opposite mode. This is only important 
   !      when the unit impulse is applied at the boundary, otherwise both end
   !      points are automatically zero. 
   ! 1.B) Then solve the homogeneous parallel streaming equation for the response of
   !      the distribution function given a field with a unit impulse.
   ! 1.C) Then integrate over velocities to get the appropriate field operator that acts
   !      on the distibution function:
   !           For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !           For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
   !           For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g
   !      Then we need to divide by the appropriate prefactor from the fields equations
   ! 1.D) Take into account identity matricies that appear in the 
   !      response equation. (see stella manual - to come!)
   !============================================================================
   subroutine get_response_matrix_for_apar(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)

      use parallelisation_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use parameters_physics, only: include_apar, include_bpar
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: periodic, phase_shift

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      ! Local variables
      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar
      
      ! ------------------------------------------------------------------------
      !                   1.A) Initialise a unit impulse in apar
      ! ------------------------------------------------------------------------
      ! Provide a unit impulse to apar^{n+1} (or Delta apar^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[apar], into L[f1] = R[aparr] and L[f2] = i*R[apari],
      ! with f = f1 + f2; then apar = df1/daparr * aparr + df2/dapari * apari.
      ! however, we see that if aparr = apari = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives apar = df1/daparr * (aparr + i * apari) = df1/daparr * apar
      ! ------------------------------------------------------------------------
      
      ! Initialise
      phi_ext = 0.0
      apar_ext = 0.0
      bpar_ext = 0.0
      apar_ext(idx) = 1.0

      ! Need to make sure that if the mode is periodic, then the boundaries match up to
      ! a phase factor. In practice this only matters if the unit impulse is at the 
      ! boundary (i.e. <idx ==1 ) otherwise phi = 0.0 at both boundary points anyway. 
      if (periodic(iky) .and. idx == 1) apar_ext(nz_ext) = apar_ext(1) / phase_shift(iky)

      ! <dum> is a scratch array that takes the place of the pdf and the fields
      ! at the previous time level. It also replaces the other fields (apar and bpar).
      ! This is set to zero for the response matrix approach because there is no sources
      ! coming from the previous time step in the response matrix equation.
      allocate (dum(nz_ext)); dum = 0.0

      ! ------------------------------------------------------------------------
      !      1.B) Get distribution function response to unit impulse in apar
      ! ------------------------------------------------------------------------

      ! Set the flux tube index to one - eed to check, but think this is okay as 
      ! the homogeneous equation solved here for the response matrix construction 
      ! is the same for all flux tubes in the flux tube train
      it = 1

      ! Solve for the response of the distribution function given this unit impulse in <apar>
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Calculate the RHS of the GK equation (using dum=0 as the pdf at the 
         ! previous time level, and phi_ext as the input potential) and store it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, dum, apar_ext * time_upwind_plus, apar_ext, dum, dum, pdf_ext(:, ivmu))
         
         ! Given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain; the RHS is
         ! input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
         
      end do

      deallocate (dum)

      ! ------------------------------------------------------------------------
      !          1.C) Integrate over velocity to get the matrix to invert
      ! ------------------------------------------------------------------------
      ! We now have the pdf on the extended zed domain at this ky and set of 
      ! connected kx values corresponding to a unit impulse in apar at this 
      ! location now integrate over velocities to get a square response matrix.
      ! 
      ! solve_field_equations_using_pdf_response is the operator that acts on the pdf in 
      ! Parallel Ampere Law. e.g. if using g for the pdf: 
      ! 
      !          β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
      ! 
      ! and divide by the appropriate factors to get the field. 
      ! 
      ! This ends the parallelisation over velocity space, so every core 
      ! should have a copy of phi_ext.
      ! ------------------------------------------------------------------------
      call solve_field_equations_using_pdf_response(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie, nz_ext)

      ! ------------------------------------------------------------------------
      !                    1.D) Compute the matrix to invert                
      ! ------------------------------------------------------------------------
      ! The response matrix is the output from <get_fields_for_response_matrix>.
      ! We then just need to compute the matrix that is to be inverted in parallel streaming. 
      ! ------------------------------------------------------------------------

      offset_apar = 0
      if (include_apar) offset_apar = nresponse
      if (include_bpar) offset_bpar = offset_apar + nresponse

      ! Next need to create column in response matrix from phi_ext, apar_ext and bpar_ext
      ! The negative sign occurs because the matrix acts on the RHS of the streaming equation.
      ! For the location where the unit impulse is applied, the matrix to invert is:
      !                       (identity matrix - response matrix)
      ! For all other locations the matrix to invert is:
      !                             (- response_matrix)
      ! So add in contribution from identity matrix for the <idx> location for the field we are solving for:
      apar_ext(idx) = apar_ext(idx) - 1.0
      
      ! But everywhere else simply add a negative sign:
      response_matrix(iky)%eigen(ie)%zloc(:nresponse, offset_apar + idx) = -phi_ext(:nresponse)
      response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, offset_apar + idx) = -apar_ext(:nresponse)
      if (include_bpar) response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, offset_apar + idx) = -bpar_ext(:nresponse) 

      
   end subroutine get_response_matrix_for_apar

   !============================================================================
   !                         Get response matrix for bpar
   !============================================================================
   ! Here, we apply a unit impulse to bpar at a given point on the extended zed
   ! domain. This location is indicated via <idx>
   ! 
   ! Variables:
   ! ----------
   ! <iky> = which ky value (these are all independent, but determines the number
   !                         of eigen chains, and the spacings within each eigen chain)
   ! <ie>  = identifies which eigen chain we are looking at. Only modes which are in 
   !         the same eigen chain can communicate with eachother
   ! <idx> = location of unit impulse within the segment 
   ! <nz_ext>   = length of extended zed domain for the given eigen chain 
   ! <phi_ext>  = phi on the extended zed domain
   ! <apar_ext> = apar on the extended zed domain
   ! <bpar_ext> = bpar on the extended zed domain
   ! <pdf_ext>  = distribution function on the extended zed domain
   ! 
   ! What is being done: 
   ! -------------------
   ! 1.A) First initialise (real) unit impulse bpar^{n+1} at each grid point in the
   !      extended zed domain for a given connected chain. This unit impulse is 
   !      supplied at the index <idx>. If using periodic boundary conditions/or
   !      treating the zonal mode, then we need to make sure phi is also periodic
   !      in zed - so apply phase shift to opposite mode. This is only important 
   !      when the unit impulse is applied at the boundary, otherwise both end
   !      points are automatically zero. 
   ! 1.B) Then solve the homogeneous parallel streaming equation for the response of
   !      the distribution function given a field with a unit impulse.
   ! 1.C) Then integrate over velocities to get the appropriate field operator that acts
   !      on the distibution function:
   !           For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !           For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
   !           For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g 
   !      Then we need to divide by the appropriate prefactor from the fields equations
   ! 1.D) Take into account identity matricies that appear in the 
   !      response equation. (see stella manual - to come!)
   !============================================================================
   subroutine get_response_matrix_for_bpar(iky, ie, idx, nz_ext, nresponse, phi_ext, apar_ext, bpar_ext, pdf_ext)
   

      use parallelisation_layouts, only: vmu_lo
      use parameters_numerical, only: time_upwind_plus
      use parameters_physics, only: include_apar, include_bpar
      use gk_implicit_terms, only: get_gke_rhs, sweep_g_zext
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: periodic, phase_shift

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, idx, nz_ext, nresponse
      complex, dimension(:), intent(out) :: phi_ext, apar_ext, bpar_ext
      complex, dimension(:, vmu_lo%llim_proc:), intent(out) :: pdf_ext

      ! Local variables
      complex, dimension(:), allocatable :: dum
      integer :: ivmu, it
      integer :: offset_apar, offset_bpar
      
      ! ------------------------------------------------------------------------
      !                  1.A) Initialise a unit impulse in bpar
      ! ------------------------------------------------------------------------
      ! Provide a unit impulse to bpar^{n+1} (or Delta bpar^{n+1}) at the location
      ! in the extended zed domain corresponding to index 'idx'
      ! note that it is sufficient to give a unit real impulse (as opposed to
      ! separately giving real and imaginary impulse) for the following reason:
      ! split homogeneous GKE, L[f] = R[bpar], into L[f1] = R[bparr] and L[f2] = i*R[bpari],
      ! with f = f1 + f2; then bpar = df1/dbparr * bparr + df2/dbpari * bpari.
      ! however, we see that if bparr = bpari = 1, L[f1] = R[1] = L[-i*f2],
      ! and thus f2 = i * f1.  This gives bpar = df1/dbparr * (bparr + i * bpari) = df1/dbparr * bpar
      ! ------------------------------------------------------------------------
      
      ! Initialise
      phi_ext = 0.0
      apar_ext = 0.0
      bpar_ext = 0.0
      
      ! how phi^{n+1} enters the GKE depends on whether we are solving for the
      ! non-Boltzmann pdf, h, or the guiding center pdf, 'g'
      bpar_ext(idx) = time_upwind_plus

      ! Need to make sure that if the mode is periodic, then the boundaries match up to
      ! a phase factor. In practice this only matters if the unit impulse is at the 
      ! boundary (i.e. <idx ==1 ) otherwise phi = 0.0 at both boundary points anyway. 
      if (periodic(iky) .and. idx == 1) bpar_ext(nz_ext) = bpar_ext(1) / phase_shift(iky)

      ! <dum> is a scratch array that takes the place of the pdf and the fields
      ! at the previous time level. It also replaces the other fields (apar and bpar).
      ! This is set to zero for the response matrix approach because there is no sources
      ! coming from the previous time step in the response matrix equation.
      allocate (dum(nz_ext)); dum = 0.0

      ! ------------------------------------------------------------------------
      !      1.B) Get distribution function response to unit impulse in bpar
      ! ------------------------------------------------------------------------

      ! Set the flux tube index to one - eed to check, but think this is okay as 
      ! the homogeneous equation solved here for the response matrix construction 
      ! is the same for all flux tubes in the flux tube train
      it = 1

      ! Solve for the response of the distribution function given this unit impulse in <bpar>
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
      
         ! Calculate the RHS of the GK equation (using dum=0 as the pdf at the 
         ! previous time level, and phi_ext as the input potential) and store it in pdf_ext
         call get_gke_rhs(ivmu, iky, ie, dum, dum, dum, dum, dum, bpar_ext, pdf_ext(:, ivmu))
         
         ! Given the RHS of the GK equation (pdf_ext), solve for the pdf at the
         ! new time level by sweeping in zed on the extended domain; the RHS is
         ! input as 'pdf_ext' and over-written with the updated solution for the pdf
         call sweep_g_zext(iky, ie, it, ivmu, pdf_ext(:, ivmu))
         
      end do

      deallocate (dum)

      ! ------------------------------------------------------------------------
      !         1.C) Integrate over velocity to get the matrix to invert
      ! ------------------------------------------------------------------------
      ! We now have the pdf on the extended zed domain at this ky and set of 
      ! connected kx values corresponding to a unit impulse in apar at this 
      ! location now integrate over velocities to get a square response matrix.
      ! 
      ! solve_field_equations_using_pdf_response is the operator that acts on the pdf in 
      ! quasineutrality and perpendicular Ampere Law. e.g. if using g for the pdf: 
      ! 
      !           sum_s Z_s n_s [ (2B/sqrt(pi)) int dvpa int dmu J_0 * g ]
      ! and
      !          - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
      ! 
      ! and divide by the appropriate factors to get the field. 
      ! 
      ! This ends the parallelisation over velocity space, so every core 
      ! should have a copy of phi_ext.
      ! ------------------------------------------------------------------------
      call solve_field_equations_using_pdf_response(pdf_ext, phi_ext, apar_ext, bpar_ext, iky, ie, nz_ext)

      ! ------------------------------------------------------------------------
      !                     1.D) Compute the matrix to invert                
      ! ------------------------------------------------------------------------
      ! The response matrix is the output from <get_fields_for_response_matrix>.
      ! We then just need to compute the matrix that is to be inverted in parallel streaming. 


      offset_apar = 0
      if (include_apar) offset_apar = nresponse
      if (include_bpar) offset_bpar = offset_apar + nresponse

      ! Next need to create column in response matrix from phi_ext, apar_ext and bpar_ext
      ! The negative sign occurs because the matrix acts on the RHS of the streaming equation.
      ! For the location where the unit impulse is applied, the matrix to invert is:
      !                       (identity matrix - response matrix)
      ! For all other locations the matrix to invert is:
      !                             (- response_matrix)
      ! So add in contribution from identity matrix for the <idx> location for the field we are solving for:
      bpar_ext(idx) = bpar_ext(idx) - 1.0
      
      ! But everywhere else simply add a negative sign:
      response_matrix(iky)%eigen(ie)%zloc(:nresponse, offset_bpar + idx) = -phi_ext(:nresponse)
      if (include_apar) response_matrix(iky)%eigen(ie)%zloc(offset_apar + 1:nresponse + offset_apar, offset_bpar + idx) = -apar_ext(:nresponse)
      response_matrix(iky)%eigen(ie)%zloc(offset_bpar + 1:nresponse + offset_bpar, offset_bpar + idx) = -bpar_ext(:nresponse) 


   end subroutine get_response_matrix_for_bpar

   !============================================================================
   !      Solve the field equations with the distribution function response
   !============================================================================
   ! We have the response of the distribution function from a unit impulse in the 
   ! fields. This is incoming as <g>. We now use this and solve the field equations.
   ! This is done in two steps:
   ! 1.C.I)  Perform the appropriate velocity integral for the fields: 
   !         For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
   !         For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
   !         For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
   ! 1.C.II) Divide by the appropriate pre-factor in the field equations
   !============================================================================
   subroutine solve_field_equations_using_pdf_response(g, phi, apar, bpar, iky, ie, nz_ext)

      use parallelisation_layouts, only: vmu_lo
      use parameters_physics, only: include_apar, include_bpar

      implicit none

      ! Arguments
      complex, dimension(:, vmu_lo%llim_proc:), intent(in) :: g
      complex, dimension(:), intent(out) :: phi, apar, bpar
      integer, intent(in) :: iky, ie, nz_ext

      ! ------------------------------------------------------------------------
      !                        1.C.I) Velocity integrals
      ! ------------------------------------------------------------------------
      ! For phi:     2B/sqrt(π) int dvpa int dmu J_0 * g
      ! For apar:    β sum_s Z_s n_s vth 2*B0/sqrt{π} \int d^2v vpar J_0 g
      ! For bpar:    - 2β sum_s n_s T_s 2*B0/sqrt{π} \int d^2v mu J_1/a_s g 
      ! ------------------------------------------------------------------------
      call integrate_over_velocity_phi
      if (include_apar) call integrate_over_velocity_apar
      if (include_bpar) call integrate_over_velocity_bpar

      ! ------------------------------------------------------------------------
      !                  1.C.II) Divide by the appropriate prefactor
      ! ------------------------------------------------------------------------
      ! Call the appropriate subroutines depending on which fields are being simulated.
      ! 
      ! In these routine we divide by the appropriate pre-factors that appear in the 
      ! field equations in front of the fields. The exact factor will depend on 
      ! which fields are being simulated (as phi and bpar are coupled), and also 
      ! on the species options (e.g. adiabatic, modified Boltzmann etc.)
      ! 
      ! These routines cannot just be called from the field_equations modules
      ! as we need to get the correct factor given the location on the extended 
      ! zed grid, so much of these routines are about mapping the prefactors onto 
      ! the extended zed domain.
      ! ------------------------------------------------------------------------
      if (include_bpar) then
         call calculate_phi_and_bpar_for_response_matrix
      else
         call calculate_phi_for_response_matrix
      end if
      if (include_apar) call get_apar_for_response_matrix
      ! ------------------------------------------------------------------------

   contains

!===============================================================================
!                              1.C.I) Velocity integrals
!===============================================================================

      !*************************************************************************
      !                       Velocity integral for phi
      !*************************************************************************
      !                   2B/sqrt(π) int dvpa int dmu J_0 * g
      !*************************************************************************
      subroutine integrate_over_velocity_phi

         use mp, only: sum_allreduce
         use parallelisation_layouts, only: vmu_lo
         use grids_species, only: nspec, spec
         use grids_extended_zgrid, only: iz_low, iz_up
         use grids_extended_zgrid, only: ikxmod
         use grids_extended_zgrid, only: nsegments
         use grids_extended_zgrid, only: periodic
         use calculations_velocity_integrals, only: integrate_species
         use calculations_gyro_averages, only: gyro_average
         use parallelisation_layouts, only: iv_idx, imu_idx, is_idx
         use parameters_numerical, only: driftkinetic_implicit
         use calculations_velocity_integrals, only: integrate_species_ffs_rm
         use parameters_physics, only: full_flux_surface
         use arrays_gyro_averages, only: j0_B_const

         implicit none

         ! Local variables
         integer :: idx, iseg, ikx, iz, ia
         integer :: izl_offset, izup
         real, dimension(nspec) :: wgt
         complex, dimension(:), allocatable :: g0
         integer :: ivmu, imu, iv, is

         ! ---------------------------------------------------------------------
         
         allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         phi = 0.
         wgt = spec%z * spec%dens_psi0

         ia = 1
         idx = 0
         izl_offset = 0

         ! ---------------------------------------------------------------------
         !                            Integrals
         ! ---------------------------------------------------------------------
         ! Do for all segments in the chain
         do iseg = 1, nsegments(ie, iky)

            ! Make sure the boundary points are being treated correctly depending
            ! on whether the mode is periodic or not. Here, define <izup> as the 
            ! upper zed value within a segment. If the mode is periodic, then 
            ! reduce the upper bound by one, as this is a repeated point so it is 
            ! obtained using the periodicity condition. This avoids and double-counting.
            ! Only the last segment ends on the periodic point; the upper end point of the other
            ! segments is shared with the next segment, which skips it through <izl_offset>
            if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
               izup = iz_up(iseg) - 1
            else
               izup = iz_up(iseg)
            end if

            ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
            ! <ikx> value on the local domain given our position on the extended domain.
            ikx = ikxmod(iseg, ie, iky)

            ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
            ! iz is only cycling through the zed location within a given segment. 
            do iz = iz_low(iseg) + izl_offset, izup
               idx = idx + 1
               
               ! ---------------------------------------------------------------
               !                            Flux tube
               ! ---------------------------------------------------------------
               if (.not. full_flux_surface .and. (.not. driftkinetic_implicit)) then
               
                  ! First get J0 * g
                  call gyro_average(g(idx, :), iky, ikx, iz, g0)
                  
                  ! Now integrate over vpa, mu and sum over species.
                  ! This returns: 2B/sqrt(π) int dvpa int dmu J_0 * g
                  call integrate_species(g0, iz, wgt, phi(idx), reduce_in=.false.)
                  
               ! ---------------------------------------------------------------
               !                       Full Flux Surface
               ! ---------------------------------------------------------------
               else
               
                  ! First multiply by B and gyroaverage, but use the piece that
                  ! is constant in alpha. This is because for FFS the response
                  ! matrix only treats the part that is constant in alpha,
                  ! and the remainder is treated using an iterative approach.
                  do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                     iv = iv_idx(vmu_lo, ivmu)
                     imu = imu_idx(vmu_lo, ivmu)
                     is = is_idx(vmu_lo, ivmu)
                     g0(ivmu) = g(idx, ivmu) * j0_B_const(iky, ikx, iz, ivmu)
                  end do
                  
                  ! Integrate over species, using the correct weights in order
                  ! to return the constant-in-alpha component of 
                  !              2B/sqrt(π) int dvpa int dmu J_0 * g
                  call integrate_species_ffs_rm(g0, wgt, phi(idx), reduce_in=.false.)
               end if
            end do
            
            ! Set the offset to 1 - all other connected segments need to start one point
            ! displaced as they share a point with the previous segment. 
            if (izl_offset == 0) izl_offset = 1
            
         end do

         ! Collect the parts of <phi> spread out across processors.
         call sum_allreduce(phi)

      end subroutine integrate_over_velocity_phi

      !*************************************************************************
      !                            Velocity integral for apar
      !*************************************************************************
      !            β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g        
      !*************************************************************************
      subroutine integrate_over_velocity_apar

         use parallelisation_layouts, only: vmu_lo, iv_idx
         use parameters_physics, only: beta
         use grids_species, only: nspec, spec
         use grids_extended_zgrid, only: iz_low, iz_up
         use grids_extended_zgrid, only: ikxmod
         use grids_extended_zgrid, only: nsegments
         use grids_extended_zgrid, only: periodic
         use calculations_velocity_integrals, only: integrate_species
         use grids_velocity, only: vpa
         use calculations_gyro_averages, only: gyro_average
         use mp, only: sum_allreduce

         implicit none

         ! Local variables
         integer :: idx, iseg, ikx, iz, ia
         integer :: ivmu, iv
         integer :: izl_offset, izup
         real, dimension(nspec) :: wgt
         complex, dimension(:), allocatable :: g0

         ! ---------------------------------------------------------------------

         allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         apar = 0.
         wgt = spec%z * spec%dens_psi0 * spec%stm_psi0 * beta

         ia = 1
         idx = 0
         izl_offset = 0

         ! ---------------------------------------------------------------------
         !                             Integrals
         ! ---------------------------------------------------------------------
         ! Do for all segments in the chain
         do iseg = 1, nsegments(ie, iky)
            
            ! Make sure the boundary points are being treated correctly depending
            ! on whether the mode is periodic or not. Here, define <izup> as the 
            ! upper zed value within a segment. If the mode is periodic, then 
            ! reduce the upper bound by one, as this is a repeated point so it is 
            ! obtained using the periodicity condition. This avoids and double-counting.
            ! Only the last segment ends on the periodic point; the upper end point of the other
            ! segments is shared with the next segment, which skips it through <izl_offset>
            if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
               izup = iz_up(iseg) - 1
            else
               izup = iz_up(iseg)
            end if

            ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
            ! <ikx> value on the local domain given our position on the extended domain.
            ikx = ikxmod(iseg, ie, iky)

            ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
            ! iz is only cycling through the zed location within a given segment. 
            do iz = iz_low(iseg) + izl_offset, izup
               idx = idx + 1
               
               ! First get J0 * g
               call gyro_average(g(idx, :), iky, ikx, iz, g0)
               
               ! Multiply by vpa to give: J0 * g * vpa
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  iv = iv_idx(vmu_lo, ivmu)
                  g0(ivmu) = g0(ivmu) * vpa(iv)
               end do
               
               ! Now integrate over vpa, mu and sum over species.
               ! This returns: β sum_s Z_s n_s vth 2*B0/sqrt{π} int d^2v vpar J_0 g
               call integrate_species(g0, iz, wgt, apar(idx), reduce_in=.false.)
               
            end do
            
            ! Set the offset to 1 - all other connected segments need to start one point
            ! displaced as they share a point with the previous segment. 
            if (izl_offset == 0) izl_offset = 1
            
         end do

         ! Collect the parts of <apar> spread out across processors. 
         call sum_allreduce(apar)

      end subroutine integrate_over_velocity_apar

      !*************************************************************************
      !                         Velocity integral for bpar                      
      !*************************************************************************
      !        - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g          
      !*************************************************************************
      subroutine integrate_over_velocity_bpar

         use parallelisation_layouts, only: vmu_lo, imu_idx
         use grids_species, only: nspec, spec
         use parameters_physics, only: beta
         use grids_extended_zgrid, only: iz_low, iz_up
         use grids_extended_zgrid, only: ikxmod
         use grids_extended_zgrid, only: nsegments
         use grids_extended_zgrid, only: periodic
         use calculations_velocity_integrals, only: integrate_species
         use grids_velocity, only: mu
         use calculations_gyro_averages, only: gyro_average_j1
         use mp, only: sum_allreduce

         implicit none

         ! Local variables
         integer :: idx, iseg, ikx, iz, ia, imu, ivmu
         integer :: izl_offset, izup
         real, dimension(nspec) :: wgt
         complex, dimension(:), allocatable :: g0

         ! ---------------------------------------------------------------------
         
         allocate (g0(vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         bpar = 0.
         wgt = -2.0 * beta * spec%temp_psi0 * spec%dens_psi0

         ia = 1
         idx = 0
         izl_offset = 0
         
         ! ---------------------------------------------------------------------
         !                          Integrals
         ! ---------------------------------------------------------------------
         ! Do for all segments in the chain
         do iseg = 1, nsegments(ie, iky)
         
            ! Make sure the boundary points are being treated correctly depending
            ! on whether the mode is periodic or not. Here, define <izup> as the 
            ! upper zed value within a segment. If the mode is periodic, then 
            ! reduce the upper bound by one, as this is a repeated point so it is 
            ! obtained using the periodicity condition. This avoids and double-counting.
            ! Only the last segment ends on the periodic point; the upper end point of the other
            ! segments is shared with the next segment, which skips it through <izl_offset>
            if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
               izup = iz_up(iseg) - 1
            else
               izup = iz_up(iseg)
            end if

            ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
            ! <ikx> value on the local domain given our position on the extended domain.
            ikx = ikxmod(iseg, ie, iky)

            ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
            ! iz is only cycling through the zed location within a given segment. 
            do iz = iz_low(iseg) + izl_offset, izup
               idx = idx + 1
               
               ! First get J1/a_s * g , where a_s is the argument of the Bessel
               ! function. Note that 'J1' in the code is actually J1/a_s
               call gyro_average_j1(g(idx, :), iky, ikx, iz, g0)
               
               ! Multiply by mu to get: J1/a_s * g * mu
               do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                  imu = imu_idx(vmu_lo, ivmu)
                  g0(ivmu) = g0(ivmu) * mu(imu)
               end do
               
               ! Now integrate over vpa, mu and sum over species.
               ! This returns: - 2β sum_s n_s T_s 2*B0/sqrt{π} int d^2v mu J_1/a_s g 
               call integrate_species(g0, iz, wgt, bpar(idx), reduce_in=.false.)
               
            end do
            
            ! Set the offset to 1 - all other connected segments need to start one point
            ! displaced as they share a point with the previous segment. 
            if (izl_offset == 0) izl_offset = 1
            
         end do

         ! Collect the parts of <bpar> spread out across processors. 
         call sum_allreduce(bpar)

      end subroutine integrate_over_velocity_bpar


!===============================================================================
!                    1.C.II) Divide by the appropriate prefactor
!===============================================================================

      !*************************************************************************
      !              Divide by the appropriate denominator for phi
      !*************************************************************************
      ! This is the case used if bpar is not included in the simulation, as the 
      ! field equations for phi and bpar are coupled. 
      !*************************************************************************
      subroutine calculate_phi_for_response_matrix

         use grids_z, only: nzgrid
         use grids_species, only: spec
         use grids_species, only: has_electron_species
         use geometry, only: dl_over_b
         use grids_extended_zgrid, only: iz_low, iz_up
         use grids_extended_zgrid, only: ikxmod
         use grids_extended_zgrid, only: nsegments
         use grids_extended_zgrid, only: periodic, phase_shift
         use grids_kxky, only: zonal_mode, akx
         use arrays, only: denominator_fields, denominator_fields_MBR
         use grids_species, only: adiabatic_option_switch
         use grids_species, only: adiabatic_option_fieldlineavg

         implicit none

         ! Local variables
         integer :: idx, iseg, ikx, iz, ia
         integer :: izl_offset, izup
         complex :: tmp
         real, dimension(:), allocatable :: denominator_seg

         ! ---------------------------------------------------------------------

         ia = 1
         idx = 0
         izl_offset = 0
         
         allocate (denominator_seg(-nzgrid:nzgrid))
         
         ! ---------------------------------------------------------------------
         !                             iky = ikx = 0 mode
         ! ---------------------------------------------------------------------
         ! Stella does not evolve the iky = ikx = 0 mode. Need to identify this
         ! mode and make sure it is set to zero.
         ! ---------------------------------------------------------------------

         ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
         ! <ikx> value on the local domain given our position on the extended domain.
         iseg = 1
         ikx = ikxmod(iseg, ie, iky)
         if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
            phi(:) = 0.0
            return
         end if

         ! ---------------------------------------------------------------------
         !                      Divide by correct field factor
         ! ---------------------------------------------------------------------
         ! Loop over all connected segments in a chain. 
         do iseg = 1, nsegments(ie, iky)
         
            ! Make sure the boundary points are being treated correctly depending
            ! on whether the mode is periodic or not. Here, define <izup> as the 
            ! upper zed value within a segment. If the mode is periodic, then 
            ! reduce the upper bound by one, as this is a repeated point so it is 
            ! obtained using the periodicity condition. This avoids and double-counting.
            ! Only the last segment ends on the periodic point; the upper end point of the other
            ! segments is shared with the next segment, which skips it through <izl_offset>
            if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
               izup = iz_up(iseg) - 1
            else
               izup = iz_up(iseg)
            end if

            ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
            ! <ikx> value on the local domain given our position on the extended domain.
            ikx = ikxmod(iseg, ie, iky)

            ! For the given value of ky, kx, store the appropriate denominator from 
            ! the Quasineutrality equation for this this segment. 
            denominator_seg = denominator_fields(iky, ikx, :)

            ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
            ! iz is only cycling through the zed location within a given segment. 
            do iz = iz_low(iseg) + izl_offset, izup
               idx = idx + 1
               phi(idx) = phi(idx) / denominator_seg(iz)
            end do

            ! Treat the periodic point correct by dividing by the phase shift
            if (periodic(iky)) phi(nz_ext) = phi (1) / phase_shift(iky)

            ! Set the offset to 1 - all other connected segments need to start one point
            ! displaced as they share a point with the previous segment. 
            if (izl_offset == 0) izl_offset = 1

         end do

         ! ---------------------------------------------------------------------
         !                           Adiabatic electrons
         ! ---------------------------------------------------------------------
         
         ! If using adiabatic electrons, or the modified Boltzmann response
         ! then we need to modify the denominator appropriately. 
         if (.not. has_electron_species(spec) .and. adiabatic_option_switch == adiabatic_option_fieldlineavg) then
         
            ! No connections for ky = 0
            if (zonal_mode(iky)) then
               iseg = 1
               tmp = sum(dl_over_b(ia, :) * phi)
               phi = phi + tmp * denominator_fields_MBR(ikxmod(1, ie, iky), :)
            end if
            
         end if

         deallocate (denominator_seg)

      end subroutine calculate_phi_for_response_matrix

      !*************************************************************************
      !           Divide by the appropriate denominator for phi + bpar
      !*************************************************************************
      ! This is the case used if bpar is included in the simulation, as the 
      ! field equations for phi and bpar are coupled. 
      !*************************************************************************
      subroutine calculate_phi_and_bpar_for_response_matrix

         use mp, only: mp_abort
         use grids_z, only: nzgrid
         use grids_kxky, only: zonal_mode, akx
         use grids_species, only: has_electron_species
         use grids_extended_zgrid, only: iz_low, iz_up
         use grids_extended_zgrid, only: ikxmod
         use grids_extended_zgrid, only: nsegments
         use grids_extended_zgrid, only: periodic, phase_shift
         use grids_species, only: spec
         use grids_species, only: adiabatic_option_switch
         use grids_species, only: adiabatic_option_fieldlineavg
         use arrays, only: denominator_fields_inv11
         use arrays, only: denominator_fields_inv13
         use arrays, only: denominator_fields_inv31
         use arrays, only: denominator_fields_inv33
         
         implicit none

         ! Local variables
         integer :: idx, iseg, ikx, iz, ia
         integer :: izl_offset, izup
         complex :: antot1, antot3
         real, dimension(:), allocatable :: gammainv11, gammainv13, gammainv31, gammainv33

         ! ---------------------------------------------------------------------

         ia = 1

         allocate (gammainv11(-nzgrid:nzgrid))
         allocate (gammainv13(-nzgrid:nzgrid))
         allocate (gammainv31(-nzgrid:nzgrid))
         allocate (gammainv33(-nzgrid:nzgrid))

         idx = 0
         izl_offset = 0

         ! ---------------------------------------------------------------------
         !                             iky = ikx = 0 mode
         ! ---------------------------------------------------------------------
         ! Stella does not evolve the iky = ikx = 0 mode. Need to identify this
         ! mode and make sure it is set to zero.
         ! ---------------------------------------------------------------------

         ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
         ! <ikx> value on the local domain given our position on the extended domain.
         iseg = 1
         ikx = ikxmod(iseg, ie, iky)
         if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
            phi(:) = 0.0
            bpar(:) = 0.0
            return
         end if

         ! ---------------------------------------------------------------------
         !                      Divide by correct field factor
         ! ---------------------------------------------------------------------
         ! Loop over all connected segments in a chain.
         do iseg = 1, nsegments(ie, iky)
         
            ! Make sure the boundary points are being treated correctly depending
            ! on whether the mode is periodic or not. Here, define <izup> as the 
            ! upper zed value within a segment. If the mode is periodic, then 
            ! reduce the upper bound by one, as this is a repeated point so it is 
            ! obtained using the periodicity condition. This avoids and double-counting.
            ! Only the last segment ends on the periodic point; the upper end point of the other
            ! segments is shared with the next segment, which skips it through <izl_offset>
            if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
               izup = iz_up(iseg) - 1
            else
               izup = iz_up(iseg)
            end if

            ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
            ! <ikx> value on the local domain given our position on the extended domain.
            ikx = ikxmod(iseg, ie, iky)

            ! For the given value of ky, kx, store the appropriate denominators from 
            ! the field equations (Quasineutrality and perpendicular Amperes law) for 
            ! this this segment. 
            gammainv11 = denominator_fields_inv11(iky, ikx, :)
            gammainv13 = denominator_fields_inv13(iky, ikx, :)
            gammainv31 = denominator_fields_inv31(iky, ikx, :)
            gammainv33 = denominator_fields_inv33(iky, ikx, :)

            ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
            ! iz is only cycling through the zed location within a given segment. 
            do iz = iz_low(iseg) + izl_offset, izup
               idx = idx + 1
               antot1 = phi(idx)
               antot3 = bpar(idx)
               phi(idx) = antot1 * gammainv11(iz) + antot3 * gammainv13(iz)
               bpar(idx) = antot1 * gammainv31(iz) + antot3 * gammainv33(iz)
            end do

            ! Treat the periodic point correct by dividing by the phase shift
            if (periodic(iky)) phi(nz_ext) = phi(1) / phase_shift(iky)
            if (periodic(iky)) bpar(nz_ext) = bpar(1) / phase_shift(iky)
            
            ! Set the offset to 1 - all other connected segments need to start one point
            ! displaced as they share a point with the previous segment. 
            if (izl_offset == 0) izl_offset = 1
            
         end do

         ! ---------------------------------------------------------------------
         !                           Adiabatic electrons
         ! ---------------------------------------------------------------------
         ! If using adiabatic electrons, or the modified Boltzmann response
         ! then we need to modify the denominator appropriately. 
         if (.not. has_electron_species(spec) .and. &
            adiabatic_option_switch == adiabatic_option_fieldlineavg) then
            call mp_abort('adiabatic electrons not yet supported for include_bpar = T. aborting.')
         end if

         deallocate (gammainv11, gammainv13, gammainv31, gammainv33)

      end subroutine calculate_phi_and_bpar_for_response_matrix

      !*************************************************************************
      !              Divide by the appropriate denominator for apar
      !*************************************************************************
      ! This is the routine used if apar is included. Note that it is separated 
      ! from the phi and bpar equations, as the field equation for apar is decoupled.
      !*************************************************************************
      subroutine get_apar_for_response_matrix

         use grids_z, only: nzgrid
         use grids_extended_zgrid, only: iz_low, iz_up
         use grids_extended_zgrid, only: ikxmod
         use grids_extended_zgrid, only: nsegments
         use grids_extended_zgrid, only: periodic, phase_shift
         use grids_kxky, only: zonal_mode, akx
         use arrays, only: kperp2

         implicit none

         ! Local variables
         integer :: idx, iseg, ikx, iz, ia
         integer :: izl_offset, izup
         real, dimension(:), allocatable :: denominator

         ! ---------------------------------------------------------------------

         ia = 1
         idx = 0
         izl_offset = 0

         allocate (denominator(-nzgrid:nzgrid))
         
         ! ---------------------------------------------------------------------
         !                             iky = ikx = 0 mode
         ! ---------------------------------------------------------------------
         ! Stella does not evolve the iky = ikx = 0 mode. Need to identify this
         ! mode and make sure it is set to zero.
         ! ---------------------------------------------------------------------

         ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
         ! <ikx> value on the local domain given our position on the extended domain.
         iseg = 1
         ikx = ikxmod(iseg, ie, iky)
         if (zonal_mode(iky) .and. abs(akx(ikx)) < epsilon(0.)) then
            apar(:) = 0.0
            return
         end if

         ! ---------------------------------------------------------------------
         !                      Divide by correct field factor
         ! ---------------------------------------------------------------------
         ! Loop over all connected segments in a chain. 
         do iseg = 1, nsegments(ie, iky)
         
            ! Make sure the boundary points are being treated correctly depending
            ! on whether the mode is periodic or not. Here, define <izup> as the 
            ! upper zed value within a segment. If the mode is periodic, then 
            ! reduce the upper bound by one, as this is a repeated point so it is 
            ! obtained using the periodicity condition. This avoids and double-counting.
            ! Only the last segment ends on the periodic point; the upper end point of the other
            ! segments is shared with the next segment, which skips it through <izl_offset>
            if (periodic(iky) .and. iseg == nsegments(ie, iky)) then
               izup = iz_up(iseg) - 1
            else
               izup = iz_up(iseg)
            end if

            ! Get the appropriate indecies. Here, the <ikxmod> routine returns the 
            ! <ikx> value on the local domain given our position on the extended domain.
            ikx = ikxmod(iseg, ie, iky)

            ! For the given value of ky, kx, store the appropriate denominator from 
            ! the Amperes equation for this this segment. 
            denominator = kperp2(iky, ikx, ia, :)

            ! The <idx> index keeps track of the location on the extended zed grid, whereas the 
            ! iz is only cycling through the zed location within a given segment. 
            do iz = iz_low(iseg) + izl_offset, izup
               idx = idx + 1
               apar(idx) = apar(idx) / denominator(iz)
            end do
            
            ! Treat the periodic point correct by dividing by the phase shift
            if (periodic(iky)) apar(nz_ext) = apar(1) / phase_shift(iky)

            ! Set the offset to 1 - all other connected segments need to start one point
            ! displaced as they share a point with the previous segment. 
            if (izl_offset == 0) izl_offset = 1
            
         end do

         deallocate (denominator)

      end subroutine get_apar_for_response_matrix

   end subroutine solve_field_equations_using_pdf_response

   !============================================================================
   !                            Allocate zloc and idx                           
   !============================================================================
   ! Sets up the response matrix storage for a given (iky, ie) eigenmode.
   ! 
   ! Allocate the following:
   ! ----------------------
   ! - response_matrix%eigen%zloc is the dimension of the response matrix for a
   !   given eigen chain 
   ! - response_matrix%idx is needed to keep track of permutations to the response
   !   matrix made during LU decomposition it will be input to LU back substitution
   !   during linear solve
   ! 
   ! Every MPI process stores its own copy of the response matrix.
   !============================================================================
   subroutine setup_response_matrix_zloc_idx(iky, ie, nresponse)

      use arrays, only: response_matrix

      implicit none

      ! Arguments
      integer, intent(in) :: iky, ie, nresponse

      ! For each ky and set of connected kx values, so we must have a response matrix
      ! that is N x N , with N = number of zeds per 2pi segment x number of 2pi segments

      ! Allocate zloc as an (nresponse x nresponse) matrix if not already allocated
      if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
         allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

      ! Allocate idx as a length-nresponse vector for LU decomposition pivots
      if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
         allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))

   end subroutine setup_response_matrix_zloc_idx
 
!===============================================================================
!                    Step 2) LU decompose the response matrix
!===============================================================================
! This section LU decomposes the response matrix for each eigen chain.
! It prepares the matrix for efficient inversion during the implicit solve.
! The decomposition uses <lu_decomposition> (PLASMA or LAPACK when available,
! see utils/linear_solve). Every MPI process factorises its own copy.
!===============================================================================

   !============================================================================
   !                       Main LU decomposition routine 
   !============================================================================
   subroutine lu_decompose_response_matrix(iky)

      use timers, only: region_start, region_end
      use arrays, only: response_matrix
      use grids_extended_zgrid, only: neigen
      use linear_solve, only: lu_decomposition

      implicit none

      ! Arguments
      integer, intent(in) :: iky

      ! Local variables
      integer :: ie
      real :: dum

      ! ------------------------------------------------------------------------
      
      ! Start the timer
      call region_start('lu_decomposition')

      ! Now that we have the reponse matrix for this ky and set of connected kx values
      ! get the LU decomposition so we are ready to solve the linear system
      do ie = 1, neigen(iky)
         call lu_decomposition(response_matrix(iky)%eigen(ie)%zloc, response_matrix(iky)%eigen(ie)%idx, dum)
      end do
      
      ! Stop the timer
      call region_end('lu_decomposition')

   end subroutine lu_decompose_response_matrix

!###############################################################################
!########################### READ RESPONSE MATRIX ##############################
!###############################################################################

   !============================================================================
   !                           Read the response matrix                         
   !============================================================================
   subroutine read_response_matrix

      use arrays, only: response_matrix
      use common_types, only: response_matrix_type
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen
      use grids_extended_zgrid, only: nsegments
      use grids_extended_zgrid, only: nzed_segment
      use grids_extended_zgrid, only: periodic
      use mp, only: proc0, job, broadcast, mp_abort
      use field_equations, only: nfields
      use file_units, only: unit_response_matrix

      implicit none

      ! Local variables
      integer :: iky, ie, nz_ext
      integer :: iky_dump, neigen_dump, naky_dump, nresponse_dump
      integer :: nresponse, nresponse_per_field
      character(len=15) :: job_str
      character(len=100) :: file_name
      integer :: ie_dump, istat
      logical, parameter :: debug = .false.

      ! ------------------------------------------------------------------------

      ! All matrices handled for the job i_job are read
      ! from a single file named: responst_mat.ijob by that
      ! jobs root process
      if (proc0) then
         write (job_str, '(I1.1)') job
         file_name = './mat/response_mat.'//trim(job_str)

         open (unit=unit_response_matrix, status='old', file=file_name, &
               action='read', form='unformatted', iostat=istat)
         if (istat /= 0) then
            print *, 'Error opening response_matrix by root processor for job ', job_str
         end if

         read (unit=unit_response_matrix) naky_dump
         if (naky /= naky_dump) call mp_abort('mismatch in naky and naky_dump')
      end if

      if (.not. allocated(response_matrix)) allocate (response_matrix(naky))

      do iky = 1, naky
         if (proc0) then
            read (unit=unit_response_matrix) iky_dump, neigen_dump
            if (iky_dump /= iky .or. neigen_dump /= neigen(iky)) &
               call mp_abort('mismatch in iky_dump/neigen_dump')
         end if

         if (.not. associated(response_matrix(iky)%eigen)) &
            allocate (response_matrix(iky)%eigen(neigen(iky)))

         ! Loop over the sets of connected kx values
         do ie = 1, neigen(iky)
         
            ! Number of zeds x number of segments
            nz_ext = nsegments(ie, iky) * nzed_segment + 1

            ! Treat zonal mode specially to avoid double counting as it is periodic
            if (periodic(iky)) then
               nresponse_per_field = nz_ext - 1
            else
               nresponse_per_field = nz_ext
            end if
            nresponse = nresponse_per_field * nfields

            if (proc0) then
               read (unit=unit_response_matrix) ie_dump, nresponse_dump
               if (ie_dump /= ie .or. nresponse /= nresponse_dump) &
                  call mp_abort('mismatch in ie/nresponse_dump')
            end if

            ! For each ky and set of connected kx values,
            ! must have a response matrix that is N x N
            ! with N = number of zeds per 2pi segment x number of 2pi segments
            if (.not. associated(response_matrix(iky)%eigen(ie)%zloc)) &
               allocate (response_matrix(iky)%eigen(ie)%zloc(nresponse, nresponse))

            ! response_matrix%idx is needed to keep track of permutations
            ! to the response matrix made during LU decomposition
            ! it will be input to LU back substitution during linear solve
            if (.not. associated(response_matrix(iky)%eigen(ie)%idx)) &
               allocate (response_matrix(iky)%eigen(ie)%idx(nresponse))
            if (proc0) then
               read (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%idx
               read (unit=unit_response_matrix) response_matrix(iky)%eigen(ie)%zloc
            end if

            call broadcast(response_matrix(iky)%eigen(ie)%idx)
            call broadcast(response_matrix(iky)%eigen(ie)%zloc)

         end do
      end do

      if (proc0) close (unit_response_matrix)

      if (debug) then
         print *, 'File', file_name, ' successfully read by root proc for job: ', job_str
      end if
   end subroutine read_response_matrix


!###############################################################################
!############################# FINISH RESPONSE MATRIX ##########################
!###############################################################################

   subroutine finish_response_matrix

      use arrays, only: response_matrix

      implicit none

      if (allocated(response_matrix)) deallocate (response_matrix)
      initialised_response_matrix = .false.

   end subroutine finish_response_matrix

end module response_matrix
