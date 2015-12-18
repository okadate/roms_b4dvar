#include "cppdefs.h"
      MODULE param_adjust_mod

#if defined ADJUST_PARAM
!
!svn $Id$
!=======================================================================
!  (c) 2015 Teruhisa Okada
!=======================================================================
!                                                                      !
!  This routine time-interpolates  4DVar surface forcing increments    !
!  and then adjust nonlinear model surface forcing.  The increments    !
!  can be constant (Nfrec=1) or time interpolated between snapshots    !
!  (Nfrec>1).                                                          !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng        Nested grid number.                                    !
!     Linp      4DVar increment time index to process.                 !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: param_adjust

      CONTAINS
!
!***********************************************************************
      SUBROUTINE param_adjust (ng)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_scalars
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: itrc
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 7)
# endif
!
!  Set parameters
!
      OCEAN(ng)%p(:,iwPhy)=wPhy(ng)
      OCEAN(ng)%p(:,iwLDet)=wLDet(ng)
      OCEAN(ng)%p(:,iwSDet)=wSDet(ng)
!
# ifdef FOUR_DVAR
!
!  Adjust parameters
!
      DO itrc=1,Nparam(ng)
        IF (Lparam(itrc,ng)) THEN
          OCEAN(ng)%p(:,itrc)=OCEAN(ng)%p(:,itrc)+                      &
     &                        OCEAN(ng)%tl_p(:,itrc)
        END IF
      END DO
!
!  Set parameters standard deviation.
!
      OCEAN(ng)%e_p(:,iwPhy)=0.1_r8
      OCEAN(ng)%e_p(:,iwLDet)=1.0_r8
      OCEAN(ng)%e_p(:,iwSDet)=0.1_r8
!
!  Set parameters normalization factor.
!
      OCEAN(ng)%b_p(:,iwPhy)=1.0_r8
      OCEAN(ng)%b_p(:,iwLDet)=1.0_r8
      OCEAN(ng)%b_p(:,iwSDet)=1.0_r8
# endif
!
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 7)
# endif
      RETURN
      END SUBROUTINE param_adjust
#endif
      END MODULE param_adjust_mod
