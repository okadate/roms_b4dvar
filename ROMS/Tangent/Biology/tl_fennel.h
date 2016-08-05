      SUBROUTINE tl_biology (ng,tile)
!
!svn $Id$
!***********************************************************************
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                               Katja Fennel   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This routine computes the  biological sources and sinks for the     !
!  Fennel et at. (2006) ecosystem model. Then, it adds those terms     !
!  to the global biological fields.                                    !
!                                                                      !
!  This model is loosely based on the model by Fasham et al. (1990)    !
!  but it differs in many respects.  The detailed equations of the     !
!  nitrogen cycling component  are given in  Fennel et al. (2006).     !
!  Nitrogen is the  fundamental elemental  currency in this model.     !
!  This model was adapted from a code written originally  by  John     !
!  Moisan and Emanule DiLorenzo.                                       !
!                                                                      !
!  It is recommended to activate always the  "BIO_SEDIMENT" option     !
!  to ensure conservation of mass by converting the organic matter     !
!  that is sinking out of the bottom most grid cell into inorganic     !
!  nutrients (i.e.,  instantanaous remineralization  at the water-     !
!  sediment interface). Additionally, the "DENITRIFICATION" option     !
!  can be activated.  Hence, a fraction of the instantenous bottom     !
!  remineralization is  assumed to  occur  through  the  anearobic     !
!  (denitrification)  pathway  and  thus  lost  from the  pool  of     !
!  biologically availalbe fixed nitrogen. See Fennel et al. (2006)     !
!  for details.                                                        !
!                                                                      !
!  Additional  options can be  activated to  enable  simulation of     !
!  inorganic carbon and dissolved oxygen.  Accounting of inorganic     !
!  carbon is activated by the "CARBON" option,  and results in two     !
!  additional  biological  tracer  variables:  DIC and alkalinity.     !
!  See Fennel et al. (2008) for details.                               !
!                                                                      !
!  If the "pCO2_RZ" options is activated, in addition to "CARBON",     !
!  the carbonate system  routines by Zeebe and Wolf-Gladrow (2001)     !
!  are used,  while the  OCMIP  standard routines are the default.     !
!  There are two different ways of treating alkalinity.  It can be     !
!  treated diagnostically (default),  in this case alkalinity acts     !
!  like a passive tracer  that is  not affected  by changes in the     !
!  concentration of  nitrate or ammonium.  However,  if the option     !
!  "TALK_NONCONSERV" is used,  the alkalinity  will be affected by     !
!  sources and sinks in nitrate. See Fennel et al. (2008) for more     !
!  details.                                                            !
!                                                                      !
!  If the "OXYGEN" option is activated,  one additional biological     !
!  tracer variable for dissolved oxygen. "OXYGEN" can be activated     !
!  independently of the  "CARBON"  option. If "OCMIP_OXYGEN_SC" is     !
!  used, in addition to "OXYGEN",  the Schmidt number of oxygen in     !
!  seawater will be  computed  using the  formulation  proposed by     !
!  Keeling et al. (1998, Global Biogeochem. Cycles,  12, 141-163).     !
!  Otherwise, the Wanninkhof's (1992) formula will be used.            !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fennel, K., Wilkin, J., Levin, J., Moisan, J., O'Reilly, J.,      !
!      Haidvogel, D., 2006: Nitrogen cycling in the Mid Atlantic       !
!      Bight and implications for the North Atlantic nitrogen          !
!      budget: Results from a three-dimensional model.  Global         !
!      Biogeochemical Cycles 20, GB3007, doi:10.1029/2005GB002456.     !
!                                                                      !
!    Fennel, K., Wilkin, J., Previdi, M., Najjar, R. 2008:             !
!      Denitrification effects on air-sea CO2 flux in the coastal      !
!      ocean: Simulations for the Northwest North Atlantic.            !
!      Geophys. Res. Letters 35, L24608, doi:10.1029/2008GL036147.     !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef DISTRIBUTE
      IF (Lbiofile(iTLM)) THEN
#else
      IF (Lbiofile(iTLM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iTLM)=.FALSE.
        BIONAME(iTLM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iTLM, 15)
#endif
      CALL tl_biology_tile (ng, tile,                                   &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
#ifdef ADJUST_PARAM
     &                   Nparam(ng),                                    &
#endif
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
#endif
     &                   GRID(ng) % Hz,                                 &
     &                   GRID(ng) % tl_Hz,                              &
     &                   GRID(ng) % z_r,                                &
     &                   GRID(ng) % tl_z_r,                             &
     &                   GRID(ng) % z_w,                                &
     &                   GRID(ng) % tl_z_w,                             &
     &                   FORCES(ng) % srflx,                            &
     &                   FORCES(ng) % tl_srflx,                         &
#if defined OXYGEN
# ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
# else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % tl_sustr,                         &
     &                   FORCES(ng) % svstr,                            &
     &                   FORCES(ng) % tl_svstr,                         &
# endif
#endif
#ifdef ADJUST_PARAM
     &                   OCEAN(ng) % p,                                 &
     &                   OCEAN(ng) % tl_p,                              &
#endif
     &                   OCEAN(ng) % t,                                 &
     &                   OCEAN(ng) % tl_t)

#ifdef PROFILE
      CALL wclock_off (ng, iTLM, 15)
#endif

      RETURN
      END SUBROUTINE tl_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE tl_biology_tile (ng, tile,                             &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
#ifdef ADJUST_PARAM
     &                         UBp,                                     &
#endif
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
#endif
     &                         Hz, tl_Hz,                               &
     &                         z_r, tl_z_r,                             &
     &                         z_w, tl_z_w,                             &
     &                         srflx, tl_srflx,                         &
#if defined OXYGEN
# ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
# else
     &                         sustr, tl_sustr,                         &
     &                         svstr, tl_svstr,                         &
# endif
#endif
#ifdef ADJUST_PARAM
     &                         p, tl_p,                                 &
#endif
     &                         t, tl_t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel  ! (okada)
      USE mod_iounits   ! (okada)
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
#ifdef ADJUST_PARAM
      integer, intent(in) :: UBp
#endif
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: srflx(LBi:,LBj:)

      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_srflx(LBi:,LBj:)
# if defined OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)

      real(r8), intent(in) :: tl_sustr(LBi:,LBj:)
      real(r8), intent(in) :: tl_svstr(LBi:,LBj:)
#  endif
# endif
# ifdef ADJUST_PARAM
      real(r8), intent(in) :: p(:,:)

      real(r8), intent(inout) :: tl_p(:,:)
# endif
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: tl_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: tl_z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: tl_z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: tl_srflx(LBi:UBi,LBj:UBj)
# if defined OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: tl_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: tl_svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
# ifdef ADJUST_PARAM
      real(r8), intent(in) :: p(2,UBp)

      real(r8), intent(inout) :: tl_p(2,UBp)
# endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
#if defined PHOSPHORUS
      integer, parameter :: Nsink = 6
#else
      integer, parameter :: Nsink = 4
#endif

      integer :: Iter, i, ibio, isink, itrc, ivar, j, k, ks
      integer :: Iteradj

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

#if defined OXYGEN
      real(r8) :: u10squ
      real(r8) :: tl_u10squ
#endif
#ifdef OXYGEN
      real(r8), parameter :: OA0 = 2.00907_r8       ! Oxygen
      real(r8), parameter :: OA1 = 3.22014_r8       ! saturation
      real(r8), parameter :: OA2 = 4.05010_r8       ! coefficients
      real(r8), parameter :: OA3 = 4.94457_r8
      real(r8), parameter :: OA4 =-0.256847_r8
      real(r8), parameter :: OA5 = 3.88767_r8
      real(r8), parameter :: OB0 =-0.00624523_r8
      real(r8), parameter :: OB1 =-0.00737614_r8
      real(r8), parameter :: OB2 =-0.0103410_r8
      real(r8), parameter :: OB3 =-0.00817083_r8
      real(r8), parameter :: OC0 =-0.000000488682_r8
      real(r8), parameter :: rOxNO3= 8.625_r8       ! 138/16
      real(r8), parameter :: rOxNH4= 6.625_r8       ! 106/16
      real(r8) :: l2mol = 1000.0_r8/22.3916_r8      ! liter to mol (mol/l*l/m3)
      real(r8) :: molv = 22.3916_r8                 ! l/mol
      real(r8) :: rho_O2 = 1.42903_r8               ! g/l
      real(r8) :: mol2g_O2 = 1.42903_r8*22.3916_r8  ! g/mol
#endif

      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4

      real(r8) :: tl_Att, tl_ExpAtt, tl_Itop, tl_PAR
      real(r8) :: tl_Epp, tl_L_NH4, tl_L_NO3, tl_LTOT, tl_Vp
      real(r8) :: tl_Chl2C, tl_t_PPmax, tl_inhNH4

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac1, fac2, fac3, fac4, fac5, fac6
      real(r8) :: cffL, cffR, cu, dltL, dltR

      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3, tl_cff4, tl_cff5
      real(r8) :: tl_fac1, tl_fac2, tl_fac3, tl_fac4, tl_fac5, tl_fac6
      real(r8) :: tl_cffL, tl_cffR, tl_cu, tl_dltL, tl_dltR
      real(r8) :: fac, tl_fac, tlfac

      real(r8) :: total_N
      real(r8) :: tl_total_N

#ifdef OXYGEN
      real(r8) :: SchmidtN_Ox, O2satu, O2_Flux
      real(r8) :: TS, AA

      real(r8) :: tl_SchmidtN_Ox, tl_O2satu, tl_O2_Flux
      real(r8) :: tl_TS, tl_AA
#endif

      real(r8) :: N_Flux_Assim
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: N_Flux_Egest
      real(r8) :: N_Flux_NewProd, N_Flux_RegProd
      real(r8) :: N_Flux_Nitrifi
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_Remine
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo

      real(r8) :: tl_N_Flux_Assim
      real(r8) :: tl_N_Flux_CoagD, tl_N_Flux_CoagP
      real(r8) :: tl_N_Flux_Egest
      real(r8) :: tl_N_Flux_NewProd, tl_N_Flux_RegProd
      real(r8) :: tl_N_Flux_Nitrifi
      real(r8) :: tl_N_Flux_Pmortal, tl_N_Flux_Zmortal
      real(r8) :: tl_N_Flux_Remine
      real(r8) :: tl_N_Flux_Zexcret, tl_N_Flux_Zmetabo

      real(r8), dimension(Nsink) :: Wbio
      real(r8), dimension(Nsink) :: tl_Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      integer, dimension(IminS:ImaxS,N(ng)) :: tl_ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: tl_PARsur

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio1
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio2
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio3
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: tl_Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: tl_Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: tl_FC

      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bL1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: bR1
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc

      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: tl_qc

#ifdef PHOSPHORUS
      real(r8) :: L_PO4, LMIN, cff6
      real(r8) :: P_Flux

      real(r8) :: tl_L_PO4, tl_LMIN, tl_cff6
      real(r8) :: tl_P_Flux

      real(r8), parameter :: rOxPO4 = 106.0_r8   ! 106/1
#endif
#ifdef H2S
      real(r8) :: S_Flux
      real(r8) :: tl_S_Flux

      real(r8), parameter :: rOxH2S = 2.0_r8     !?
#endif

#define BLOWINGUP_CHECKER
#ifdef BLOWINGUP_CHECKER
      integer :: ii
      real(r8), dimension(NT(ng)) :: val
      character (len=8) :: valchar
#endif

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Add biological Source/Sink terms.
!-----------------------------------------------------------------------
!
!  Avoid computing source/sink terms if no biological iterations.
!
      IF (BioIter(ng).le.0) RETURN
!
!  Set time-stepping according to the number of iterations.
!
      dtdays=dt(ng)*sec2day/REAL(BioIter(ng),r8)
!
!  Set vertical sinking indentification vector.
!
      idsink(1)=iPhyt
      idsink(2)=iChlo
      idsink(3)=iSDeN
      idsink(4)=iLDeN
#if defined PHOSPHORUS
      idsink(5)=iSDeP
      idsink(6)=iLDeP
#endif

#ifdef ADJUST_PARAM
!
!  Convert parameters.
!
# ifdef EXP_PARAM
      tl_wPhy=wPhy(ng)*tl_p(nstp,iwPhy)
      tl_wSDet=wSDet(ng)*tl_p(nstp,iwSDet)
      tl_wLDet=wLDet(ng)*tl_p(nstp,iwLDet)
      tl_Chl2C_m=Chl2C_m(ng)*tl_p(nstp,iChl2C_m)
      tl_Vp0=Vp0(ng)*tl_p(nstp,iVp0)
      tl_R_SODf=R_SODf(ng)*tl_p(nstp,iR_SODf)
      tl_R_NH4f=R_NH4f(ng)*tl_p(nstp,iR_NH4f)
      tl_R_PO4f=R_PO4f(ng)*tl_p(nstp,iR_PO4f)
# else
      tl_AttSW=tl_p(nstp,iAttSW)
      tl_AttChl=tl_p(nstp,iAttChl)
      tl_Vp0=tl_p(nstp,iVp0)
      tl_I_thNH4=tl_p(nstp,iI_thNH4)
      tl_D_p5NH4=tl_p(nstp,iD_p5NH4)

      tl_K_Nitri=tl_p(nstp,iK_Nitri)
      tl_NitriR=tl_p(nstp,iNitriR)
      tl_K_Denit=tl_p(nstp,iK_Denit)
      tl_DenitR=tl_p(nstp,iDenitR)
      tl_K_NO3=tl_p(nstp,iK_NO3)

      tl_K_NH4=tl_p(nstp,iK_NH4)
      tl_K_PO4=tl_p(nstp,iK_PO4)
      tl_K_Phy=tl_p(nstp,iK_Phy)
      tl_Chl2C_m=tl_p(nstp,iChl2C_m)
      tl_PhyPN=tl_p(nstp,iPhyPN)
      tl_PhyCN=tl_p(nstp,iPhyCN)

      tl_PhyIP=tl_p(nstp,iPhyIP)
      tl_PhyIS=tl_p(nstp,iPhyIS)
      tl_PhyMR=tl_p(nstp,iPhyMR)
      tl_ZooAE_N=tl_p(nstp,iZooAE_N)
      tl_ZooBM=tl_p(nstp,iZooBM)

      tl_ZooPN=tl_p(nstp,iZooPN)
      tl_ZooCN=tl_p(nstp,iZooCN)
      tl_ZooER=tl_p(nstp,iZooER)
      tl_ZooGR=tl_p(nstp,iZooGR)
      tl_ZooMR=tl_p(nstp,iZooMR)
      tl_K_DO=tl_p(nstp,iK_DO)

      tl_LDeRRN=tl_p(nstp,iLDeRRN)
      tl_LDeRRP=tl_p(nstp,iLDeRRP)
      tl_CoagR=tl_p(nstp,iCoagR)
      tl_SDeRRN=tl_p(nstp,iSDeRRN)
      tl_SDeRRP=tl_p(nstp,iSDeRRP)

      tl_wPhy=tl_p(nstp,iwPhy)
      tl_wSDet=tl_p(nstp,iwSDet)
      tl_wLDet=tl_p(nstp,iwLDet)
      tl_R_SODf=tl_p(nstp,iR_SODf)
      tl_R_NH4f=tl_p(nstp,iR_NH4f)

      tl_R_PO4f=tl_p(nstp,iR_PO4f)
      tl_R_NH4f_max=tl_p(nstp,iR_NH4f_m)
      tl_R_PO4f_max=tl_p(nstp,iR_PO4f_m)
      tl_K_DO_npflux=tl_p(nstp,iK_DO_npf)
      tl_t_SODf=tl_p(nstp,it_SODf)
# endif
#endif
#ifdef CHECKER
!
! stdout 
!
      if (master.and.(mod(iic(ng)-1,ninfo(ng)).eq.0)) then
# ifdef ADJUST_PARAM
        write(stdout,101) 'tl 01:08', tl_p(nstp,1:8)
        write(stdout,101) 'tl 09:16', tl_p(nstp,9:16)
        write(stdout,101) 'tl 17:24', tl_p(nstp,17:24)
        write(stdout,101) 'tl 25:32', tl_p(nstp,25:32)
        write(stdout,101) 'tl 33:40', tl_p(nstp,33:40)
        write(stdout,101) 'tl 41:  ', tl_p(nstp,41:)
        write(stdout,*) ('-',i=1,78)
# endif
      end if
 101  FORMAT (a,8(1pe9.1))
#endif
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! phytoplankton
      Wbio(2)=wPhy(ng)                ! chlorophyll
      Wbio(3)=wSDet(ng)               ! small Nitrogen-detritus
      Wbio(4)=wLDet(ng)               ! large Nitrogen-detritus
#if defined PHOSPHORUS
      Wbio(5)=wSDet(ng)               ! small Phosphorus-detritus
      Wbio(6)=wLDet(ng)               ! large Phosphorus-detritus
#endif
#ifdef ADJUST_PARAM
      tl_Wbio(1)=tl_wPhy
      tl_Wbio(2)=tl_wPhy
      tl_Wbio(3)=tl_wSDet
      tl_Wbio(4)=tl_wLDet
# ifdef PHOSPHORUS
      tl_Wbio(5)=tl_wSDet
      tl_Wbio(6)=tl_wLDet
# endif
#else
      tl_Wbio(1)=0.0_r8
      tl_Wbio(2)=0.0_r8
      tl_Wbio(3)=0.0_r8
      tl_Wbio(4)=0.0_r8
# ifdef PHOSPHORUS
      tl_Wbio(5)=0.0_r8
      tl_Wbio(6)=0.0_r8
# endif
#endif
!
!  Compute inverse thickness to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
            tl_Hz_inv(i,k)=-Hz_inv(i,k)*Hz_inv(i,k)*tl_Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
            tl_Hz_inv2(i,k)=-Hz_inv2(i,k)*Hz_inv2(i,k)*                 &
     &                      (tl_Hz(i,j,k)+tl_Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
            tl_Hz_inv3(i,k)=-Hz_inv3(i,k)*Hz_inv3(i,k)*                 &
     &                      (tl_Hz(i,j,k-1)+tl_Hz(i,j,k)+               &
     &                       tl_Hz(i,j,k+1))
          END DO
        END DO
!
!  Clear tl_Bio and Bio arrays.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,ibio)=0.0_r8
              Bio1(i,k,ibio)=0.0_r8
              Bio2(i,k,ibio)=0.0_r8
              Bio3(i,k,ibio)=0.0_r8
              Bio_old(i,k,ibio)=0.0_r8
              tl_Bio(i,k,ibio)=0.0_r8
              tl_Bio_old(i,k,ibio)=0.0_r8
            END DO
          END DO
#ifdef BLOWINGUP_CHECKER
          val(ibio)=0.0_r8
#endif
        END DO
!
!  Extract biological variables from tracer arrays, place them into
!  scratch arrays, and restrict their values to be positive definite.
!  At input, all tracers (index nnew) from predictor step have
!  transport units (m Tunits) since we do not have yet the new
!  values for zeta and Hz. These are known after the 2D barotropic
!  time-stepping.
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              tlfac=0.5_r8-SIGN(0.5_r8,-t(i,j,k,nstp,ibio))
              tl_Bio_old(i,k,ibio)=tlfac*tl_t(i,j,k,nstp,ibio)

              Bio(i,k,ibio)=Bio_old(i,k,ibio)
              tl_Bio(i,k,ibio)=tl_Bio_old(i,k,ibio)
            END DO
          END DO
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
            tlfac=0.5_r8+SIGN(0.5_r8,35.0_r8-t(i,j,k,nstp,itemp))
            tl_Bio(i,k,itemp)=tlfac*tl_t(i,j,k,nstp,itemp)

            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
            tlfac=0.5_r8+SIGN(0.5_r8,t(i,j,k,nstp,isalt))
            tl_Bio(i,k,isalt)=tlfac*tl_t(i,j,k,nstp,isalt)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
          tl_PARsur(i)=PARfrac(ng)*tl_srflx(i,j)*rho0*Cp
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
!  During the iterative procedure a series of fractional time steps are
!  performed in a chained mode (splitting by different biological
!  conversion processes) in sequence of the main food chain.  In all
!  stages the concentration of the component being consumed is treated
!  in fully implicit manner, so the algorithm guarantees non-negative
!  values, no matter how strong s the concentration of active consuming
!  component (Phytoplankton or Zooplankton).  The overall algorithm,
!  as well as any stage of it, is formulated in conservative form
!  (except explicit sinking) in sense that the sum of concentration of
!  all components is conserved.
!
!
!  In the implicit algorithm, we have for example (N: nitrate,
!                                                  P: phytoplankton),
!
!     N(new) = N(old) - uptake * P(old)     uptake = mu * N / (Kn + N)
!                                                    {Michaelis-Menten}
!  below, we set
!                                           The N in the numerator of
!     cff = mu * P(old) / (Kn + N(old))     uptake is treated implicitly
!                                           as N(new)
!
!  so the time-stepping of the equations becomes:
!
!     N(new) = N(old) / (1 + cff)     (1) when substracting a sink term,
!                                         consuming, divide by (1 + cff)
!  and
!
!     P(new) = P(old) + cff * N(new)  (2) when adding a source term,
!                                         growing, add (cff * source)
!
!  Notice that if you substitute (1) in (2), you will get:
!
!     P(new) = P(old) + cff * N(old) / (1 + cff)    (3)
!
!  If you add (1) and (3), you get
!
!     N(new) + P(new) = N(old) + P(old)
!
!  implying conservation regardless how "cff" is computed. Therefore,
!  this scheme is unconditionally stable regardless of the conversion
!  rate. It does not generate negative values since the constituent
!  to be consumed is always treated implicitly. It is also biased
!  toward damping oscillations.
!
!  The iterative loop below is to iterate toward an universal Backward-
!  Euler treatment of all terms. So if there are oscillations in the
!  system, they are only physical oscillations. These iterations,
!  however, do not improve the accuaracy of the solution.
!
        ITER_LOOP: DO Iter=1,BioIter(ng)
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
#if defined OXYGEN && defined DENITRIFICATION
!-----------------------------------------------------------------------
!  Denitrification in anoxic water                      Okada 2014/02/13
!-----------------------------------------------------------------------
#endif
#ifdef ADJUST_PARAM
# include <tl_fennel_1_param.h>
!#else
!# include <tl_fennel_1.h>
#endif
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
#ifdef ADJUST_PARAM
# include <tl_fennel_2_param.h>
!#else
!# include <tl_fennel_2.h>
#endif
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
#if defined H2S && defined OXYGEN
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
#endif
#ifdef OXYGEN
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
#endif
#ifdef ADJUST_PARAM
# include <tl_fennel_3_param.h>
!#else
!# include <tl_fennel_3.h>
#endif
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
#ifdef BIO_SEDIMENT
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
#endif
#if defined BIO_SED_CONSTANT
!-----------------------------------------------------------------------
!  Elution and oxygen consumption from/by sediment. (Okada, 2014/02/13)
!-----------------------------------------------------------------------
#endif
#ifdef ADJUST_PARAM
# include <tl_fennel_4_param.h>
!#else
! include <tl_fennel_4.h>
#endif
!
        END DO ITER_LOOP
!
!-----------------------------------------------------------------------
!  Update global tracer variables: Add increment due to BGC processes
!  to tracer array in time index "nnew". Index "nnew" is solution after
!  advection and mixing and has transport units (m Tunits) hence the
!  increment is multiplied by Hz.  Notice that we need to subtract
!  original values "Bio_old" at the top of the routine to just account
!  for the concentractions affected by BGC processes. This also takes
!  into account any constraints (non-negative concentrations, carbon
!  concentration range) specified before entering BGC kernel. If "Bio"
!  were unchanged by BGC processes, the increment would be exactly
!  zero. Notice that final tracer values, t(:,:,:,nnew,:) are not
!  bounded >=0 so that we can preserve total inventory of N and
!  C even when advection causes tracer concentration to go negative.
!  (J. Wilkin and H. Arango, Apr 27, 2012)
!-----------------------------------------------------------------------
!
        DO itrc=1,NBT
          ibio=idbio(itrc)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
              tl_cff=tl_Bio(i,k,ibio)-tl_Bio_old(i,k,ibio)
              
!>            t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
              tl_t(i,j,k,nnew,ibio)=tl_t(i,j,k,nnew,ibio)+              &
     &                              tl_cff*Hz(i,j,k)+cff*tl_Hz(i,j,k)
#ifdef BLOWINGUP_CHECKER
              val(ibio)=val(ibio)+tl_t(i,j,k,nnew,ibio)
#endif
            END DO
          END DO
        END DO
!
      END DO J_LOOP

#ifdef BLOWINGUP_CHECKER
!
!  If blowing-up, set exit_flag to stop computations. (okada)
!
      DO itrc=1,NBT
        ibio=idbio(itrc)
        WRITE (valchar,'(1pe8.1)') val(ibio)
        DO ii=1,8
          IF ((valchar(ii:ii).eq.'N').or.(valchar(ii:ii).eq.'n').or.    &
     &        (valchar(ii:ii).eq.'*')) THEN
            IF (Master) WRITE (stdout,100) ibio
            exit_flag=1
          END IF
        END DO
      END DO
 100  FORMAT ('Blowing-up in tl_fennel.h, varid=',i2)
#endif

      RETURN
      END SUBROUTINE tl_biology_tile
