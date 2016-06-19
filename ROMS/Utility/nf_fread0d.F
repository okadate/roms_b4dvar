#include "cppdefs.h"
      MODULE nf_fread0d_mod
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2015 Teruhisa Okada                                   !
!=======================================================================
!                                                                      !
!  This function reads in a generic floating point 0D array from an    !
!  input NetCDF file.                                                  !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number (integer)                          !
!     model      Calling model identifier (integer)                    !
!     ncname     NetCDF file name (string)                             !
!     ncid       NetCDF file ID (integer)                              !
!     ncvname    NetCDF variable name (string)                         !
!     ncvarid    NetCDF variable ID (integer)                          !
!     tindex     NetCDF time record index to read (integer)            !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     A          Field to read in (real 0D array)                      !
!     nf_fread0d Error flag (integer)                                  !
!                                                                      !
!=======================================================================
!
      implicit none

      CONTAINS
!
!***********************************************************************
      FUNCTION nf_fread0d (ng, model, ncname, ncid,                     &
     &                     ncvname, ncvarid,                            &
     &                     tindex,                                      &
     &                     A)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf
# endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid, ncvarid, tindex

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: ncvname

      real(r8), intent(out) :: A
!
!  Local variable declarations.
!
      logical, dimension(3) :: foundit

      integer :: ic, status

      integer, dimension(1) :: start, total

      integer :: nf_fread0d

      real(r8) :: Afactor, Aoffset, Aspval

      real(r8), dimension(3) :: AttValue
      real(r8), dimension(1) :: wrk

      character (len=12), dimension(3) :: AttName
!
!-----------------------------------------------------------------------
!  Set starting and ending indices to process.
!-----------------------------------------------------------------------
!
!  Check if the following attributes: "scale_factor", "add_offset", and
!  "_FillValue" are present in the input NetCDF variable:
!
!  If the "scale_value" attribute is present, the data is multiplied by
!  this factor after reading.
!  If the "add_offset" attribute is present, this value is added to the
!  data after reading.
!  If both "scale_factor" and "add_offset" attributes are present, the
!  data are first scaled before the offset is added.
!  If the "_FillValue" attribute is present, the data having this value
!  is treated as missing and it is replaced with zero. This feature it
!  is usually related with the land/sea masking.
!
      AttName(1)='scale_factor'
      AttName(2)='add_offset  '
      AttName(3)='_FillValue  '

      CALL netcdf_get_fatt (ng, model, ncname, ncvarid, AttName,        &
     &                      AttValue, foundit,                          &
     &                      ncid = ncid)
      IF (exit_flag.ne.NoError) THEN
        nf_fread0d=ioerror
        RETURN
      END IF

      IF (.not.foundit(1)) THEN
        Afactor=1.0_r8
      ELSE
        Afactor=AttValue(1)
      END IF

      IF (.not.foundit(2)) THEN
        Aoffset=0.0_r8
      ELSE
        Aoffset=AttValue(2)
      END IF

      IF (.not.foundit(3)) THEN
        Aspval=spval_check
      ELSE
        Aspval=AttValue(3)
      END IF

      start(1) = tindex
      total(1) = 1
!
!-----------------------------------------------------------------------
!  Serial I/O: Read in requested field and scale it.
!-----------------------------------------------------------------------
!
      status=nf90_noerr
      IF (InpThread) THEN
        status=nf90_get_var(ncid, ncvarid, wrk, start, total)
        IF (status.eq.nf90_noerr) THEN
          wrk=Afactor*wrk+Aoffset
        END IF
      END IF
      IF (status.ne.nf90_noerr) THEN
        exit_flag=2
        ioerror=status
        nf_fread0d=status
        RETURN
      END IF

      A=wrk(1)
      nf_fread0d=status

#ifdef DISTRIBUTE
!
!  Broadcast new solution to other nodes.
!
      CALL mp_bcastf (ng, model, A)
#endif

      RETURN
      END FUNCTION nf_fread0d

      END MODULE nf_fread0d_mod
