#include "cppdefs.h"
      MODULE nf_fwrite0d_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2015 Teruhisa Okada                                   !
!=======================================================================
!                                                                      !
!  This function writes out a generic floating point 2D array into an  !
!  output NetCDF file.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number.                                 !
!     model        Calling model identifier.                           !
!     ncid         NetCDF file ID.                                     !
!     ncvarid      NetCDF variable ID.                                 !
!     tindex       NetCDF time record index to write.                  !
!     A            Field to write out (real).                          !
!     SetFillVal   Logical switch to set fill value in land areas      !
!                    (optional).                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     nf_fwrite  Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      FUNCTION nf_fwrite0d (ng, model, ncid, ncvarid, tindex,           &
     &                      A, SetFillVal)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
!  Imported variable declarations.
!
      logical, intent(in), optional :: SetFillVal

      integer, intent(in) :: ng, model, ncid, ncvarid, tindex

      real(r8), intent(in) :: A
!
!  Local variable declarations.
!
      integer :: ic, Npts
      integer :: status

      integer, dimension(1) :: start, total

      integer :: nf_fwrite0d

      real(r8), dimension(1) :: Awrk
!
!  Initialize local array to avoid denormalized numbers. This
!  facilitates processing and debugging.
!
      Awrk=0.0_r8
!
!-----------------------------------------------------------------------
!  If serial or shared-memory applications and serial output, pack data
!  into a global 1D array in column-major order.
!-----------------------------------------------------------------------
!
      Awrk(1)=A
!
!-----------------------------------------------------------------------
!  Write output buffer into NetCDF file.
!-----------------------------------------------------------------------
!
      nf_fwrite0d=nf90_noerr
      IF (OutThread) THEN
        start(1)=tindex
        total(1)=1
        status=nf90_put_var(ncid, ncvarid, Awrk, start, total)
      END IF

      nf_fwrite0d=status

      RETURN
      END FUNCTION nf_fwrite0d

      END MODULE nf_fwrite0d_mod
