#include "cppdefs.h"
      MODULE ini_param_mod
#if defined NONLINEAR && defined ADJUST_PARAM
!
!svn $Id$
!=======================================================================
!  (c) 2016-01-20 Teruhisa Okada                                       !
!=======================================================================
!                                                                      !
!  This routine initializes parameter to adjust.                       !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: ini_param
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE ini_param (ng)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
      USE mod_ocean
      USE mod_biology
# ifdef CHECKER
      USE mod_parallel
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      integer :: it
      integer :: Tindex
      logical, dimension(Nparam(ng)) :: Lrestart
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 2)
# endif
!
!  Initialize time stepping indices and counters.
!
      Tindex=1
      DO it=1,Nparam(ng)
        Lrestart(it)=.FALSE.
      END DO
!
!  Initialize parameters
!
# ifdef CHECKER
      if (master) print*, 'ini_param', Lparam(iwPhy,ng), Ldef_p(iwPhy,ng)
      if (master) print*, 'ini_param', OCEAN(ng)%p(:,iwPhy)
# endif
# ifdef RESTART_PARAM
      DO it=1,Nparam(ng)
        IF (.not.Ldef_p(it,ng).and.Lparam(it,ng)) THEN
          Lrestart(it)=.TRUE.
        END IF
      END DO
# endif
# ifdef EXP_PARAM
      DO it=1,Nparam(ng)
        IF (.not.Lrestart(it)) THEN
          OCEAN(ng)%p(Tindex,it)=0.0_r8
        END IF
      END DO
# else
      IF (.not.Lrestart(iAttSW)) THEN
        OCEAN(ng)%p(Tindex,iAttSW)=AttSW(ng)
      END IF
      IF (.not.Lrestart(iAttChl)) THEN
        OCEAN(ng)%p(Tindex,iAttChl)=AttChl(ng)
      END IF
      IF (.not.Lrestart(iVp0)) THEN
        OCEAN(ng)%p(Tindex,iVp0)=Vp0(ng)
      END IF
      IF (.not.Lrestart(iI_thNH4)) THEN
        OCEAN(ng)%p(Tindex,iI_thNH4)=I_thNH4(ng)
      END IF
      IF (.not.Lrestart(iD_p5NH4)) THEN
        OCEAN(ng)%p(Tindex,iD_p5NH4)=D_p5NH4(ng)
      END IF

      IF (.not.Lrestart(iK_Nitri)) THEN
        OCEAN(ng)%p(Tindex,iK_Nitri)=K_Nitri(ng)
      END IF
      IF (.not.Lrestart(iNitriR)) THEN
        OCEAN(ng)%p(Tindex,iNitriR)=NitriR(ng)
      END IF
      IF (.not.Lrestart(iK_Denit)) THEN
        OCEAN(ng)%p(Tindex,iK_Denit)=K_Denit(ng)
      END IF
      IF (.not.Lrestart(iDenitR)) THEN
        OCEAN(ng)%p(Tindex,iDenitR)=DenitR(ng)
      END IF
      IF (.not.Lrestart(iK_NO3)) THEN
        OCEAN(ng)%p(Tindex,iK_NO3)=K_NO3(ng)
      END IF

      IF (.not.Lrestart(iK_NH4)) THEN
        OCEAN(ng)%p(Tindex,iK_NH4)=K_NH4(ng)
      END IF
      IF (.not.Lrestart(iK_PO4)) THEN
        OCEAN(ng)%p(Tindex,iK_PO4)=K_PO4(ng)
      END IF
      IF (.not.Lrestart(iK_Phy)) THEN
        OCEAN(ng)%p(Tindex,iK_Phy)=K_Phy(ng)
      END IF
      IF (.not.Lrestart(iChl2C_m)) THEN
        OCEAN(ng)%p(Tindex,iChl2C_m)=Chl2C_m(ng)
      END IF
      IF (.not.Lrestart(iPhyPN)) THEN
        OCEAN(ng)%p(Tindex,iPhyPN)=PhyPN(ng)
      END IF
      IF (.not.Lrestart(iPhyCN)) THEN
        OCEAN(ng)%p(Tindex,iPhyCN)=PhyCN(ng)
      END IF

      IF (.not.Lrestart(iPhyIP)) THEN
        OCEAN(ng)%p(Tindex,iPhyIP)=PhyIP(ng)
      END IF
      IF (.not.Lrestart(iPhyIS)) THEN
        OCEAN(ng)%p(Tindex,iPhyIS)=PhyIS(ng)
      END IF
      IF (.not.Lrestart(iPhyMR)) THEN
        OCEAN(ng)%p(Tindex,iPhyMR)=PhyMR(ng)
      END IF
      IF (.not.Lrestart(iZooAE_N)) THEN
        OCEAN(ng)%p(Tindex,iZooAE_N)=ZooAE_N(ng)
      END IF
      IF (.not.Lrestart(iZooBM)) THEN
        OCEAN(ng)%p(Tindex,iZooBM)=ZooBM(ng)
      END IF

      IF (.not.Lrestart(iZooPN)) THEN
        OCEAN(ng)%p(Tindex,iZooPN)=ZooPN(ng)
      END IF
      IF (.not.Lrestart(iZooCN)) THEN
        OCEAN(ng)%p(Tindex,iZooCN)=ZooCN(ng)
      END IF
      IF (.not.Lrestart(iZooER)) THEN
        OCEAN(ng)%p(Tindex,iZooER)=ZooER(ng)
      END IF
      IF (.not.Lrestart(iZooGR)) THEN
        OCEAN(ng)%p(Tindex,iZooGR)=ZooGR(ng)
      END IF
      IF (.not.Lrestart(iZooMR)) THEN
        OCEAN(ng)%p(Tindex,iZooMR)=ZooMR(ng)
      END IF
      IF (.not.Lrestart(iK_DO)) THEN
        OCEAN(ng)%p(Tindex,iK_DO)=K_DO(ng)
      END IF

      IF (.not.Lrestart(iLDeRRN)) THEN
        OCEAN(ng)%p(Tindex,iLDeRRN)=LDeRRN(ng)
      END IF
      IF (.not.Lrestart(iLDeRRP)) THEN
        OCEAN(ng)%p(Tindex,iLDeRRP)=LDeRRP(ng)
      END IF
      IF (.not.Lrestart(iCoagR)) THEN
        OCEAN(ng)%p(Tindex,iCoagR)=CoagR(ng)
      END IF
      IF (.not.Lrestart(iSDeRRN)) THEN
        OCEAN(ng)%p(Tindex,iSDeRRN)=SDeRRN(ng)
      END IF
      IF (.not.Lrestart(iSDeRRP)) THEN
        OCEAN(ng)%p(Tindex,iSDeRRP)=SDeRRP(ng)
      END IF

      IF (.not.Lrestart(iwPhy)) THEN
        OCEAN(ng)%p(Tindex,iwPhy)=wPhy(ng)
      END IF
      IF (.not.Lrestart(iwSDet)) THEN
        OCEAN(ng)%p(Tindex,iwSDet)=wSDet(ng)
      END IF
      IF (.not.Lrestart(iwLDet)) THEN
        OCEAN(ng)%p(Tindex,iwLDet)=wLDet(ng)
      END IF
      IF (.not.Lrestart(iR_SODf)) THEN
        OCEAN(ng)%p(Tindex,iR_SODf)=R_SODf(ng)
      END IF
      IF (.not.Lrestart(iR_NH4f)) THEN
        OCEAN(ng)%p(Tindex,iR_NH4f)=R_NH4f(ng)
      END IF

      IF (.not.Lrestart(iR_PO4f)) THEN
        OCEAN(ng)%p(Tindex,iR_PO4f)=R_PO4f(ng)
      END IF
# endif
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
          OCEAN(ng)%p(:,it)=OCEAN(ng)%p(Tindex,it)
        END IF
      END DO
# ifdef CHECKER
      if (master) print*, 'ini_param', OCEAN(ng)%p(:,iwPhy)
# endif
!
!  Set parameters standard deviation.
!
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
# ifdef EXP_PARAM
          OCEAN(ng)%e_p(:,it)=pfactor(it,ng)
# else
          OCEAN(ng)%e_p(:,it)=(pmax(it,ng)-pmin(it,ng))*pfactor(it,ng)
# endif
        END IF
      END DO
!
!  Set parameters normalization factor.
!
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
          OCEAN(ng)%b_p(:,it)=1.0_r8
        END IF
      END DO
!
# ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2)
# endif
      RETURN
      END SUBROUTINE ini_param
#endif
      END MODULE ini_param_mod
