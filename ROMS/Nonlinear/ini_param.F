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
!
# ifdef PROFILE
      CALL wclock_on (ng, iNLM, 2)
# endif
!
!  Initialize time stepping indices and counters.
!
      Tindex=1
!
!  Initialize parameters
!
# ifdef CHECKER
      if (master) print*, 'ini_param', Ldef_p(iwPhy,ng), Lparam(iwPhy,ng)
      if (master) print*, 'ini_param', OCEAN(ng)%p(:,iwPhy)
# endif
# ifdef EXP_PARAM
      DO it=1,Nparam(ng)
        IF (Ldef_p(it,ng).and.Lparam(it,ng)) THEN
          OCEAN(ng)%p(Tindex,it)=0.0_r8
        END IF
      END DO
# else
#  ifdef RESTART_PARAM
      DO it=1,Nparam(ng)
        IF (Lparam(it,ng)) THEN
          Ldef_p(it,ng)=.TRUE.
        END IF
      END DO
#  endif
      IF (Ldef_p(iAttSW,ng).and.Lparam(iAttSW,ng)) THEN
        OCEAN(ng)%p(Tindex,iAttSW)=AttSW(ng)
      END IF
      IF (Ldef_p(iAttChl,ng).and.Lparam(iAttChl,ng)) THEN
        OCEAN(ng)%p(Tindex,iAttChl)=AttChl(ng)
      END IF
      IF (Ldef_p(iVp0,ng).and.Lparam(iVp0,ng)) THEN
        OCEAN(ng)%p(Tindex,iVp0)=Vp0(ng)
      END IF
      IF (Ldef_p(iI_thNH4,ng).and.Lparam(iI_thNH4,ng)) THEN
        OCEAN(ng)%p(Tindex,iI_thNH4)=I_thNH4(ng)
      END IF
      IF (Ldef_p(iD_p5NH4,ng).and.Lparam(iD_p5NH4,ng)) THEN
        OCEAN(ng)%p(Tindex,iD_p5NH4)=D_p5NH4(ng)
      END IF

      IF (Ldef_p(iK_Nitri,ng).and.Lparam(iK_Nitri,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_Nitri)=K_Nitri(ng)
      END IF
      IF (Ldef_p(iNitriR,ng).and.Lparam(iNitriR,ng)) THEN
        OCEAN(ng)%p(Tindex,iNitriR)=NitriR(ng)
      END IF
      IF (Ldef_p(iK_Denit,ng).and.Lparam(iK_Denit,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_Denit)=K_Denit(ng)
      END IF
      IF (Ldef_p(iDenitR,ng).and.Lparam(iDenitR,ng)) THEN
        OCEAN(ng)%p(Tindex,iDenitR)=DenitR(ng)
      END IF
      IF (Ldef_p(iK_NO3,ng).and.Lparam(iK_NO3,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_NO3)=K_NO3(ng)
      END IF

      IF (Ldef_p(iK_NH4,ng).and.Lparam(iK_NH4,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_NH4)=K_NH4(ng)
      END IF
      IF (Ldef_p(iK_PO4,ng).and.Lparam(iK_PO4,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_PO4)=K_PO4(ng)
      END IF
      IF (Ldef_p(iK_Phy,ng).and.Lparam(iK_Phy,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_Phy)=K_Phy(ng)
      END IF
      IF (Ldef_p(iChl2C_m,ng).and.Lparam(iChl2C_m,ng)) THEN
        OCEAN(ng)%p(Tindex,iChl2C_m)=Chl2C_m(ng)
      END IF
      IF (Ldef_p(iPhyPN,ng).and.Lparam(iPhyPN,ng)) THEN
        OCEAN(ng)%p(Tindex,iPhyPN)=PhyPN(ng)
      END IF
      IF (Ldef_p(iPhyCN,ng).and.Lparam(iPhyCN,ng)) THEN
        OCEAN(ng)%p(Tindex,iPhyCN)=PhyCN(ng)
      END IF

      IF (Ldef_p(iPhyIP,ng).and.Lparam(iPhyIP,ng)) THEN
        OCEAN(ng)%p(Tindex,iPhyIP)=PhyIP(ng)
      END IF
      IF (Ldef_p(iPhyIS,ng).and.Lparam(iPhyIS,ng)) THEN
        OCEAN(ng)%p(Tindex,iPhyIS)=PhyIS(ng)
      END IF
      IF (Ldef_p(iPhyMR,ng).and.Lparam(iPhyMR,ng)) THEN
        OCEAN(ng)%p(Tindex,iPhyMR)=PhyMR(ng)
      END IF
      IF (Ldef_p(iZooAE_N,ng).and.Lparam(iZooAE_N,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooAE_N)=ZooAE_N(ng)
      END IF
      IF (Ldef_p(iZooBM,ng).and.Lparam(iZooBM,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooBM)=ZooBM(ng)
      END IF

      IF (Ldef_p(iZooPN,ng).and.Lparam(iZooPN,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooPN)=ZooPN(ng)
      END IF
      IF (Ldef_p(iZooCN,ng).and.Lparam(iZooCN,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooCN)=ZooCN(ng)
      END IF
      IF (Ldef_p(iZooER,ng).and.Lparam(iZooER,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooER)=ZooER(ng)
      END IF
      IF (Ldef_p(iZooGR,ng).and.Lparam(iZooGR,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooGR)=ZooGR(ng)
      END IF
      IF (Ldef_p(iZooMR,ng).and.Lparam(iZooMR,ng)) THEN
        OCEAN(ng)%p(Tindex,iZooMR)=ZooMR(ng)
      END IF
      IF (Ldef_p(iK_DO,ng).and.Lparam(iK_DO,ng)) THEN
        OCEAN(ng)%p(Tindex,iK_DO)=K_DO(ng)
      END IF

      IF (Ldef_p(iLDeRRN,ng).and.Lparam(iLDeRRN,ng)) THEN
        OCEAN(ng)%p(Tindex,iLDeRRN)=LDeRRN(ng)
      END IF
      IF (Ldef_p(iLDeRRP,ng).and.Lparam(iLDeRRP,ng)) THEN
        OCEAN(ng)%p(Tindex,iLDeRRP)=LDeRRP(ng)
      END IF
      IF (Ldef_p(iCoagR,ng).and.Lparam(iCoagR,ng)) THEN
        OCEAN(ng)%p(Tindex,iCoagR)=CoagR(ng)
      END IF
      IF (Ldef_p(iSDeRRN,ng).and.Lparam(iSDeRRN,ng)) THEN
        OCEAN(ng)%p(Tindex,iSDeRRN)=SDeRRN(ng)
      END IF
      IF (Ldef_p(iSDeRRP,ng).and.Lparam(iSDeRRP,ng)) THEN
        OCEAN(ng)%p(Tindex,iSDeRRP)=SDeRRP(ng)
      END IF

      IF (Ldef_p(iwPhy,ng).and.Lparam(iwPhy,ng)) THEN
        OCEAN(ng)%p(Tindex,iwPhy)=wPhy(ng)
      END IF
      IF (Ldef_p(iwSDet,ng).and.Lparam(iwSDet,ng)) THEN
        OCEAN(ng)%p(Tindex,iwSDet)=wSDet(ng)
      END IF
      IF (Ldef_p(iwLDet,ng).and.Lparam(iwLDet,ng)) THEN
        OCEAN(ng)%p(Tindex,iwLDet)=wLDet(ng)
      END IF
      IF (Ldef_p(iR_SODf,ng).and.Lparam(iR_SODf,ng)) THEN
        OCEAN(ng)%p(Tindex,iR_SODf)=R_SODf(ng)
      END IF
      IF (Ldef_p(iR_NH4f,ng).and.Lparam(iR_NH4f,ng)) THEN
        OCEAN(ng)%p(Tindex,iR_NH4f)=R_NH4f(ng)
      END IF

      IF (Ldef_p(iR_PO4f,ng).and.Lparam(iR_PO4f,ng)) THEN
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
