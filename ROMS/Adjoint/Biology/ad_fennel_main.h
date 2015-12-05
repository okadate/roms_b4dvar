      SUBROUTINE ad_biology (ng,tile)
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
      IF (Lbiofile(iADM)) THEN
#else
      IF (Lbiofile(iADM).and.(tile.eq.0)) THEN
#endif
        Lbiofile(iADM)=.FALSE.
        BIONAME(iADM)=__FILE__
      END IF
!
#ifdef PROFILE
      CALL wclock_on (ng, iADM, 15)
#endif
      CALL ad_biology_tile (ng, tile,                                   &
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % ad_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % ad_z_r,                          &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % ad_z_w,                          &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % ad_srflx,                      &
#ifdef OXYGEN
# ifdef BULK_FLUXES
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
# else
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % ad_sustr,                      &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % ad_svstr,                      &
# endif
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % ad_t)

#ifdef PROFILE
      CALL wclock_off (ng, iADM, 15)
#endif
      RETURN
      END SUBROUTINE ad_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE ad_biology_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, UBk, UBt,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
#ifdef MASKING
     &                            rmask,                                &
#endif
     &                            Hz, ad_Hz,                            &
     &                            z_r, ad_z_r,                          &
     &                            z_w, ad_z_w,                          &
     &                            srflx, ad_srflx,                      &
#ifdef OXYGEN
# ifdef BULK_FLUXES
     &                            Uwind, Vwind,                         &
# else
     &                            sustr, ad_sustr,                      &
     &                            svstr, ad_svstr,                      &
# endif
#endif
     &                            t, ad_t)
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

      real(r8), intent(inout) :: ad_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: ad_z_r(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ad_srflx(LBi:,LBj:)
# ifdef OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
#  else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)

      real(r8), intent(inout) :: ad_sustr(LBi:,LBj:)
      real(r8), intent(inout) :: ad_svstr(LBi:,LBj:)
#  endif
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(in) :: srflx(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_Hz(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(in) :: ad_z_r(LBi:UBi,LBj:UBj,UBk)
      real(r8), intent(inout) :: ad_z_w(LBi:UBi,LBj:UBj,0:UBk)
      real(r8), intent(inout) :: ad_srflx(LBi:UBi,LBj:UBj)
# ifdef OXYGEN
#  ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
#  else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)

      real(r8), intent(inout) :: ad_sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: ad_svstr(LBi:UBi,LBj:UBj)
#  endif
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
!
!  Local variable declarations.
!
#ifdef PHOSPHORUS
      integer, parameter :: Nsink = 6
#else
      integer, parameter :: Nsink = 4
#endif

      integer :: Iter, i, ibio, isink, itrc, ivar, j, k, ks
      integer :: Iteradj, kk

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

#ifdef OXYGEN
      real(r8) :: u10squ
      real(r8) :: ad_u10squ
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

      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR, PAR1
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4

      real(r8) :: ad_Att, ad_AttFac, ad_ExpAtt, ad_Itop, ad_PAR
      real(r8) :: ad_Epp, ad_L_NH4, ad_L_NO3, ad_LTOT, ad_Vp
      real(r8) :: ad_Chl2C, ad_t_PPmax, ad_inhNH4

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR

      real(r8) :: ad_cff, ad_cff1, ad_cff2, ad_cff3, ad_cff4, ad_cff5
      real(r8) :: ad_fac1, ad_fac2, ad_fac3
      real(r8) :: ad_cffL, ad_cffR, ad_cu, ad_dltL, ad_dltR
      real(r8) :: fac, adfac, adfac1, adfac2, adfac3

      real(r8) :: total_N
      real(r8) :: ad_total_N

#ifdef OXYGEN
      real(r8) :: SchmidtN_Ox, O2satu, O2_Flux
      real(r8) :: TS, AA

      real(r8) :: ad_SchmidtN_Ox, ad_O2satu, ad_O2_Flux
      real(r8) :: ad_TS, ad_AA
#endif

      real(r8) :: N_Flux_Assim
      real(r8) :: N_Flux_CoagD, N_Flux_CoagP
      real(r8) :: N_Flux_Egest
      real(r8) :: N_Flux_NewProd, N_Flux_RegProd
      real(r8) :: N_Flux_Nitrifi
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_Remine
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo

      real(r8) :: ad_N_Flux_Assim
      real(r8) :: ad_N_Flux_CoagD, ad_N_Flux_CoagP
      real(r8) :: ad_N_Flux_Egest
      real(r8) :: ad_N_Flux_NewProd, ad_N_Flux_RegProd
      real(r8) :: ad_N_Flux_Nitrifi
      real(r8) :: ad_N_Flux_Pmortal, ad_N_Flux_Zmortal
      real(r8) :: ad_N_Flux_Remine
      real(r8) :: ad_N_Flux_Zexcret, ad_N_Flux_Zmetabo

      real(r8), dimension(Nsink) :: Wbio
      real(r8), dimension(Nsink) :: ad_Wbio

      integer, dimension(IminS:ImaxS,N(ng)) :: ksource

      real(r8), dimension(IminS:ImaxS) :: PARsur
      real(r8), dimension(IminS:ImaxS) :: ad_PARsur

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio1
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio2
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: Bio_old

      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: ad_Bio
      real(r8), dimension(IminS:ImaxS,N(ng),NT(ng)) :: ad_Bio_old

      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: ad_FC

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

      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_WL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_bR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: ad_qc

#ifdef PHOSPHORUS
      real(r8) :: L_PO4, LMIN, cff6
      real(r8) :: P_Flux

      real(r8) :: ad_L_PO4, ad_LMIN, ad_cff6
      real(r8) :: ad_P_Flux

      real(r8), parameter :: rOxPO4 = 106.0_r8   ! 106/1
#endif
#ifdef H2S
      real(r8) :: S_Flux
      real(r8) :: ad_S_Flux

      real(r8), parameter :: rOxH2S = 2.0_r8     !?
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
!
!  Set vertical sinking velocity vector in the same order as the
!  identification vector, IDSINK.
!
      Wbio(1)=wPhy(ng)                ! phytoplankton
      Wbio(2)=wPhy(ng)                ! chlorophyll
      Wbio(3)=wSDet(ng)               ! small Nitrogen-detritus
      Wbio(4)=wLDet(ng)               ! large Nitrogen-detritus

!>    tl_Wbio(1)=tl_wPhy(ng)                ! phytoplankton
!>    tl_Wbio(2)=tl_wPhy(ng)                ! chlorophyll
!>    tl_Wbio(3)=tl_wSDet(ng)               ! small Nitrogen-detritus
!>    tl_Wbio(4)=tl_wLDet(ng)               ! large Nitrogen-detritus

      ad_Wbio(1)=0.0_r8
      ad_Wbio(2)=0.0_r8
      ad_Wbio(3)=0.0_r8
      ad_Wbio(4)=0.0_r8
#ifdef PHOSPHORUS
      Wbio(5)=wSDet(ng)               ! small Phosphorus-detritus
      Wbio(6)=wLDet(ng)               ! large Phosphorus-detritus

!>    tl_Wbio(5)=tl_wSDet(ng)               ! small Phosphorus-detritus
!>    tl_Wbio(6)=tl_wLDet(ng)               ! large Phosphorus-detritus

      ad_Wbio(5)=0.0_r8
      ad_Wbio(6)=0.0_r8
#endif
      J_LOOP : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Initialize adjoint private variables.
!-----------------------------------------------------------------------
!
        ad_Att=0.0_r8
        ad_AttFac=0.0_r8
        ad_ExpAtt=0.0_r8
        ad_Itop=0.0_r8
        ad_PAR=0.0_r8
        ad_Epp=0.0_r8
        ad_L_NH4=0.0_r8
        ad_L_NO3=0.0_r8
        ad_LTOT=0.0_r8
        ad_Vp=0.0_r8
        ad_Chl2C=0.0_r8
        ad_t_PPmax=0.0_r8
        ad_inhNH4=0.0_r8

        ad_cff=0.0_r8
        ad_cff1=0.0_r8
        ad_cff2=0.0_r8
        ad_cff3=0.0_r8
        ad_cff4=0.0_r8
        ad_cff5=0.0_r8
        ad_fac1=0.0_r8
        ad_fac2=0.0_r8
        ad_fac3=0.0_r8
        ad_cffL=0.0_r8
        ad_cffR=0.0_r8
        ad_cu=0.0_r8
        ad_dltL=0.0_r8
        ad_dltR=0.0_r8
        fac=0.0_r8
        adfac=0.0_r8
        adfac1=0.0_r8
        adfac2=0.0_r8
        adfac3=0.0_r8

        ad_total_N=0.0_r8
#ifdef OXYGEN
        ad_SchmidtN_Ox=0.0_r8
        ad_O2satu=0.0_r8
        ad_O2_Flux=0.0_r8
        ad_TS=0.0_r8
        ad_AA=0.0_r8
#endif
        ad_N_Flux_Assim=0.0_r8
        ad_N_Flux_CoagD=0.0_r8
        ad_N_Flux_CoagP=0.0_r8
        ad_N_Flux_Egest=0.0_r8
        ad_N_Flux_NewProd=0.0_r8
        ad_N_Flux_RegProd=0.0_r8
        ad_N_Flux_Nitrifi=0.0_r8
        ad_N_Flux_Pmortal=0.0_r8
        ad_N_Flux_Zmortal=0.0_r8
        ad_N_Flux_Remine=0.0_r8
        ad_N_Flux_Zexcret=0.0_r8
        ad_N_Flux_Zmetabo=0.0_r8
#ifdef PHOSPHORUS
        ad_L_PO4=0.0_r8
        ad_LMIN=0.0_r8
        ad_cff6=0.0_r8
        ad_P_Flux=0.0_r8
#endif
#ifdef H2S
        ad_S_Flux=0.0_r8
#endif
!
        DO k=1,N(ng)
          DO i=IminS,ImaxS
            ad_Hz_inv(i,k)=0.0_r8
            ad_Hz_inv2(i,k)=0.0_r8
            ad_Hz_inv3(i,k)=0.0_r8
            ad_WL(i,k)=0.0_r8
            ad_WR(i,k)=0.0_r8
            ad_bL(i,k)=0.0_r8
            ad_bR(i,k)=0.0_r8
            ad_qc(i,k)=0.0_r8
          END DO
        END DO
        DO i=IminS,ImaxS
          ad_PARsur(i)=0.0_r8
        END DO
        DO k=0,N(ng)
          DO i=IminS,ImaxS
            ad_FC(i,k)=0.0_r8
          END DO
        END DO
!
!  Compute inverse thickness to avoid repeated divisions.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
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
              ad_Bio(i,k,ibio)=0.0_r8
              ad_Bio_old(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
!
#include "ad_fennel_bs0.h"
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
              
!>            tl_t(i,j,k,nnew,ibio)=tl_t(i,j,k,nnew,ibio)+              &
!>   &                              tl_cff*Hz(i,j,k)+cff*tl_Hz(i,j,k)
              ad_Hz(i,j,k)=ad_Hz(i,j,k)+cff*ad_t(i,j,k,nnew,ibio)
              ad_cff=ad_cff+Hz(i,j,k)*ad_t(i,j,k,nnew,ibio)

!>            tl_cff=tl_Bio(i,k,ibio)-tl_Bio_old(i,k,ibio)
              ad_Bio_old(i,k,ibio)=ad_Bio_old(i,k,ibio)-ad_cff
              ad_Bio(i,k,ibio)=ad_Bio_old(i,k,ibio)+ad_cff
              ad_cff=0.0_r8 !okada!
            END DO
          END DO
        END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
        ITER_LOOP1: DO Iter=BioIter(ng),1,-1
!
#if defined BIO_SED_CONSTANT
!-----------------------------------------------------------------------
!  Elution and oxygen consumption from/by sediment. (Okada, 2014/02/13)
!-----------------------------------------------------------------------
#endif
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
#include "tl_fennel_4.h"
!
#ifdef OXYGEN
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
#endif
#if defined H2S && defined OXYGEN
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
#include "tl_fennel_3.h"
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
#include "ad_fennel_2.h"
!
#if defined OXYGEN && defined DENITRIFICATION
!-----------------------------------------------------------------------
!  Denitrification in anoxic water                      Okada 2014/02/13
!-----------------------------------------------------------------------
#endif
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
#include "ad_fennel_1.h"
!
        END DO ITER_LOOP1
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
!>        tl_PARsur(i)=PARfrac(ng)*tl_srflx(i,j)*rho0*Cp
          ad_srflx(i,j)=ad_srflx(i,j)+PARfrac(ng)*ad_PARsur(i)*rho0*Cp
          ad_PARsur(i)=0.0_r8  !okada!
        END DO
!
!  Extract potential temperature and salinity.
!
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tlfac=0.5_r8+SIGN(0.5_r8,t(i,j,k,nstp,isalt))
!>          tl_Bio(i,k,isalt)=tlfac*tl_t(i,j,k,nstp,isalt)
            adfac=0.5_r8+SIGN(0.5_r8,t(i,j,k,nstp,isalt))
            ad_t(i,j,k,nstp,itemp)=ad_t(i,j,k,nstp,itemp)+              &
     &                             adfac*ad_Bio(i,k,itemp)
            ad_Bio(i,k,isalt)=0.0_r8

!>          tlfac=0.5_r8+SIGN(0.5_r8,35.0_r8-t(i,j,k,nstp,itemp))
!>          tl_Bio(i,k,itemp)=tlfac*tl_t(i,j,k,nstp,itemp)
            adfac=0.5_r8+SIGN(0.5_r8,35.0_r8-t(i,j,k,nstp,itemp))
            ad_t(i,j,k,nstp,itemp)=ad_t(i,j,k,nstp,itemp)+              &
     &                             adfac*ad_Bio(i,k,itemp)
            ad_Bio(i,k,itemp)=0.0_r8
          END DO
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
!>            tl_Bio(i,k,ibio)=tl_Bio_old(i,k,ibio)
              ad_Bio_old(i,k,ibio)=ad_Bio_old(i,k,ibio)+ad_Bio(i,k,ibio)
              ad_Bio(i,k,ibio)=0.0_r8

!>            tlfac=0.5_r8-SIGN(0.5_r8,-t(i,j,k,nstp,ibio))
!>            tl_Bio_old(i,k,ibio)=tlfac*tl_t(i,j,k,nstp,ibio)
              adfac=0.5_r8-SIGN(0.5_r8,-t(i,j,k,nstp,ibio))
              ad_t(i,j,k,nstp,ibio)=ad_t(i,j,k,nstp,ibio)+              &
     &                              adfac*ad_Bio_old(i,k,ibio)
              ad_Bio_old(i,k,ibio)=0.0_r8
            END DO
          END DO
        END DO
!
!  Compute inverse thickness to avoid repeated divisions.
!
        DO k=2,N(ng)-1
          DO i=Istr,Iend
!>          tl_Hz_inv3(i,k)=-Hz_inv3(i,k)*Hz_inv3(i,k)*                 &
!>   &                      (tl_Hz(i,j,k-1)+tl_Hz(i,j,k)+               &
!>   &                       tl_Hz(i,j,k+1))
            adfac=Hz_inv3(i,k)*Hz_inv3(i,k)*ad_Hz_inv3(i,k)
            ad_Hz(i,j,k-1)=ad_Hz(i,j,k-1)-adfac
            ad_Hz(i,j,k  )=ad_Hz(i,j,k  )-adfac
            ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)-adfac
            ad_Hz_inv3(i,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
!>          tl_Hz_inv2(i,k)=-Hz_inv2(i,k)*Hz_inv2(i,k)*                 &
!>   &                      (tl_Hz(i,j,k)+tl_Hz(i,j,k+1))
            adfac=Hz_inv2(i,k)*Hz_inv2(i,k)*ad_Hz_inv2(i,k)
            ad_Hz(i,j,k  )=ad_Hz(i,j,k  )-adfac
            ad_Hz(i,j,k+1)=ad_Hz(i,j,k+1)-adfac
            ad_Hz_inv2(i,k)=0.0_r8
          END DO
        END DO
        DO k=1,N(ng)
          DO i=Istr,Iend
!>          tl_Hz_inv(i,k)=-Hz_inv(i,k)*Hz_inv(i,k)*tl_Hz(i,j,k)
            ad_Hz(i,j,k)=ad_Hz(i,j,k)-                                  &
     &                   Hz_inv(i,k)*Hz_inv(i,k)*ad_Hz_inv(i,k)
            ad_Hz_inv(i,k)=0.0_r8
          END DO
        END DO
!
      END DO J_LOOP
!
!  Set adjoint parameters
!






      RETURN
      END SUBROUTINE tl_biology_tile