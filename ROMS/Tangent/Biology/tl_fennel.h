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
     &                      LBi, UBi, LBj, UBj, N(ng), NT(ng),          &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      nstp(ng), nnew(ng),                         &
#ifdef MASKING
     &                      GRID(ng) % rmask,                           &
#endif
     &                      GRID(ng) % Hz,                              &
     &                      GRID(ng) % tl_Hz,                           &
     &                      GRID(ng) % z_r,                             &
     &                      GRID(ng) % tl_z_r,                          &
     &                      GRID(ng) % z_w,                             &
     &                      GRID(ng) % tl_z_w,                          &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % tl_srflx,                      &
#ifdef OXYGEN
# ifdef BULK_FLUXES
     &                      FORCES(ng) % Uwind,                         &
     &                      FORCES(ng) % Vwind,                         &
# else
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % tl_sustr,                      &
     &                      FORCES(ng) % svstr,                         &
     &                      FORCES(ng) % tl_svstr,                      &
# endif
#endif
     &                      OCEAN(ng) % t,                              &
     &                      OCEAN(ng) % tl_t)

#ifdef PROFILE
      CALL wclock_off (ng, iTLM, 15)
#endif
      RETURN
      END SUBROUTINE tl_biology
!
!-----------------------------------------------------------------------
      SUBROUTINE tl_biology_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, UBk, UBt,         &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            nstp, nnew,                           &
#ifdef MASKING
     &                            rmask,                                &
#endif
     &                            Hz, tl_Hz,                            &
     &                            z_r, tl_z_r,                          &
     &                            z_w, tl_z_w,                          &
     &                            srflx, tl_srflx                       &
#ifdef OXYGEN
# ifdef BULK_FLUXES
     &                            Uwind, Vwind,                         &
# else
     &                            sustr, tl_sustr,                      &
     &                            svstr, tl_svstr,                      &
# endif
#endif
     &                            t, tl_t)
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

      real(r8), intent(in) :: tl_Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: tl_z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: tl_srflx(LBi:,LBj:)
# ifdef OXYGEN
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
# ifdef OXYGEN
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
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,UBk,3,UBt)
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
      integer :: Iteradj

      integer, dimension(Nsink) :: idsink

      real(r8), parameter :: eps = 1.0e-20_r8

#ifdef OXYGEN
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
      real(r8) :: l2mol = 1000.0_r8/22.3916_r8      ! liter to mol
#endif

      real(r8) :: Att, AttFac, ExpAtt, Itop, PAR
      real(r8) :: Epp, L_NH4, L_NO3, LTOT, Vp
      real(r8) :: Chl2C, dtdays, t_PPmax, inhNH4

      real(r8) :: tl_Att, tl_ExpAtt, tl_Itop, tl_PAR
      real(r8) :: tl_Epp, tl_L_NH4, tl_L_NO3, tl_LTOT, tl_Vp
      real(r8) :: tl_Chl2C, tl_t_PPmax, tl_inhNH4

      real(r8) :: cff, cff1, cff2, cff3, cff4, cff5
      real(r8) :: fac1, fac2, fac3
      real(r8) :: cffL, cffR, cu, dltL, dltR

      real(r8) :: tl_cff, tl_cff1, tl_cff2, tl_cff3, tl_cff4, tl_cff5
      real(r8) :: tl_fac1, tl_fac2, tl_fac3
      real(r8) :: tl_cffL, tl_cffR, tl_cu, tl_dltL, tl_dltR

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
      real(r8) :: N_Flux_SumProd                  ! (okada)
      real(r8) :: N_Flux_Nitrifi
      real(r8) :: N_Flux_Pmortal, N_Flux_Zmortal
      real(r8) :: N_Flux_Remine
      real(r8) :: N_Flux_Zexcret, N_Flux_Zmetabo

      real(r8) :: tl_N_Flux_Assim
      real(r8) :: tl_N_Flux_CoagD, tl_N_Flux_CoagP
      real(r8) :: tl_N_Flux_Egest
      real(r8) :: tl_N_Flux_NewProd, tl_N_Flux_RegProd
      real(r8) :: tl_N_Flux_SumProd
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
      real(r8) :: L_PO4
      real(r8) :: P_Flux_SumProd
      real(r8) :: P_Flux_Remine

      real(r8) :: tl_L_PO4
      real(r8) :: tl_P_Flux_SumProd
      real(r8) :: tl_P_Flux_Remine

      real(r8), parameter :: rOxPO4 = 106.0_r8   ! 106/1
#endif
#ifdef H2S
      real(r8) :: S_Flux
      real(r8) :: tl_S_Flux

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
#ifdef defined PHOSPHORUS
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

      tl_Wbio(1)=tl_wPhy(ng)                ! phytoplankton
      tl_Wbio(2)=tl_wPhy(ng)                ! chlorophyll
      tl_Wbio(3)=tl_wSDet(ng)               ! small Nitrogen-detritus
      tl_Wbio(4)=tl_wLDet(ng)               ! large Nitrogen-detritus
#ifdef PHOSPHORUS
      Wbio(5)=wSDet(ng)               ! small Phosphorus-detritus
      Wbio(6)=wLDet(ng)               ! large Phosphorus-detritus

      tl_Wbio(5)=tl_wSDet(ng)               ! small Phosphorus-detritus
      tl_Wbio(6)=tl_wLDet(ng)               ! large Phosphorus-detritus
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
              tl_Bio(i,k,ibio)=0.0_r8
            END DO
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
              Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
              tl_Bio_old(i,k,ibio)=MAX(0.0_r8,tl_t(i,j,k,nstp,ibio))
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
            tl_Bio(i,k,itemp)=tl_t(i,j,k,nstp,itemp)
            Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
            tl_Bio(i,k,isalt)=tl_t(i,j,k,nstp,isalt)
          END DO
        END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
        DO i=Istr,Iend
          PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
          tl_PARsur(i)=(tl_PARfrac(ng)*srflx(i,j)+                      &
     &                  PARfrac(ng)*tl_srflx(i,j))*rho0*Cp
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
!  Compute appropriate basic state arrays I.
!
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
                Bio(i,k,ibio)=Bio_old(i,k,ibio)
              END DO
            END DO
          END DO
!
!  Extract potential temperature and salinity.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
              Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              AttFac=0.0_r8
              IF (PARsur(i).gt.0.0_r8) THEN
                DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
                  Att=(AttSW(ng)+                                       &
     &                 AttChl(ng)*Bio(i,k,iChlo)+                       &
     &                 AttFac)*                                         &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell cent
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                  cff=PhyCN(ng)*12.0_r8
                  Chl2C=MIN(Bio(i,k,iChlo)/(Bio(i,k,iPhyt)*cff+eps),    &
     &                      Chl2C_m(ng))
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                  Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                  fac1=PAR*PhyIS(ng)
                  Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                  t_PPmax=Epp*fac1
#ifdef PHOSPHORUS
!
!  Nutrient-limitation terms (Laurent et al. 2012).
!
#else
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
#endif
                  cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                  cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                  inhNH4=1.0_r8/(1.0_r8+cff1)
                  L_NH4=cff1/(1.0_r8+cff1)
                  L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                  LTOT=L_NO3+L_NH4
#ifdef PHOSPHORUS
                  cff3=Bio(i,k,iPO4_)*K_PO4(ng)
                  L_PO4=cff3/(1.0_r8+cff3)
!
!  Nitrate, ammonium and phosphate uptake by Phytoplankton.
!
#else
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
#endif
                  fac1=dtdays*t_PPmax
                  cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*             &
     &                 Bio(i,k,iPhyt)
                  cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio(i,k,iPhyt)
                  N_Flux_NewProd=Bio(i,k,iNO3_)/(1.0_r8+cff4)*cff4
                  N_Flux_RegProd=Bio(i,k,iNH4_)/(1.0_r8+cff5)*cff5
                  N_Flux_SumProd=N_Flux_NewProd+N_Flux_RegProd
                  Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
                  Bio1(i,k,iNH4_)=Bio(i,k,iNH4_)
                  Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
#ifdef PHOSPHORUS
                  Bio1(i,k,iPO4_)=Bio(i,k,iPO4_)
                  IF (LTOT.lt.L_PO4) THEN
                    Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                    Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                    Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
                    Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-                      &
     &                             PhyPN(ng)*N_Flux_SumProd
                  ELSE IF (L_PO4.lt.LTOT) THEN
                    cff4=fac1*K_PO4(ng)/(1.0_r8+cff3)*Bio(i,k,iPhyt)
                    Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff4)
                    P_Flux_SumProd=Bio(i,k,iPO4_)*cff4
                    fac1=MIN(P_Flux_SumProd/PhyPN(ng)/N_Flux_SumProd,   &
     &                       1.0_r8)
                    Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd*fac1
                    Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd*fac1
                    Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd*fac1
                    LTOT=L_PO4
                  ENDIF
#else
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
#endif
!
                  Bio(i,k,iChlo)=Bio(i,k,iChlo)+                        &
     &                           (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*     &
     &                            Chl2C_m(ng)*Bio(i,k,iChlo))/          &
     &                           (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
#ifdef OXYGEN
                  Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                        &
     &                           N_Flux_NewProd*rOxNO3+                 &
     &                           N_Flux_RegProd*rOxNH4
#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
#ifdef OXYGEN
                  fac2=MAX(Bio(i,k,iOxyg),0.0_r8)     ! O2 max
                  fac3=MAX(fac2/(K_Nitri(ng)+fac2),0.0_r8) ! MM for O2 d
                  fac1=dtdays*NitriR(ng)*fac3
#else
                  fac1=dtdays*NitriR(ng)
#endif
                  cff1=(PAR-I_thNH4(ng))/                               &
     &                 (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                  cff2=1.0_r8-MAX(0.0_r8,cff1)
                  cff3=fac1*cff2
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                  N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
              ELSE
                cff3=dtdays*NitriR(ng)
                DO k=N(ng),1,-1
                  Bio1(i,k,iNH4_)=Bio(i,k,iNH4_)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                  N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                  Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                  Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
                END DO
                DO k=N(ng),1,-1
                  Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
#ifdef PHOSPHORUS
                  Bio1(i,k,iPO4_)=Bio(i,k,iPO4_)
#endif
                END DO
              END IF
            END DO

            IF (Iteradj.ne.Iter) THEN
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
!
              fac1=dtdays*ZooGR(ng)
              cff2=dtdays*PhyMR(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
!
                  cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/              &
     &                 (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                  cff3=1.0_r8/(1.0_r8+cff1)
                  Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
                  Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
                  N_Flux_Assim=cff1*Bio(i,k,iPhyt)*ZooAE_N(ng)
                  N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)+                        &
     &                           N_Flux_Assim
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                        &
     &                           N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
                  N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),    &
     &                                    0.0_r8)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
                  Bio(i,k,iChlo)=Bio(i,k,iChlo)-                        &
     &                           cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),    &
     &                                    0.0_r8)
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                        &
     &                           N_Flux_Pmortal
#ifdef PHOSPHORUS
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                        &
     &                           PhyPN(ng)*(N_Flux_Egest+               &
     &                                      N_Flux_Pmortal)+            &
     &                           (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
#endif
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
              cff1=dtdays*ZooBM(ng)
              fac2=dtdays*ZooMR(ng)
              fac3=dtdays*ZooER(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/              &
     &                 (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                  cff2=fac2*Bio(i,k,iZoop)
                  cff3=fac1*ZooAE_N(ng)
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)/                        &
     &                           (1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
                  N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
                  N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
                  N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),    &
     &                                    0.0_r8)
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                        &
     &                           rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
#endif
#ifdef PHOSPHORUS
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                        &
     &                           ZooPN(ng)*(N_Flux_Zmetabo+             &
     &                                      N_Flux_Zexcret)
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                        &
     &                           ZooPN(ng)*N_Flux_Zmortal
#endif
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
              fac1=dtdays*CoagR(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff1=fac1*(Bio(i,k,iSDeN)+Bio(i,k,iPhyt))
                  cff2=1.0_r8/(1.0_r8+cff1)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
                  Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                  N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
                  N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
                  Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                        &
     &                           N_Flux_CoagP+N_Flux_CoagD
#ifdef PHOSPHORUS
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
                  Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                        &
     &                           PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
#endif
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
#ifdef OXYGEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
                  fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 depe
                  cff1=dtdays*SDeRRN(ng)*fac2
                  cff2=1.0_r8/(1.0_r8+cff1)
                  cff3=dtdays*LDeRRN(ng)*fac2
                  cff4=1.0_r8/(1.0_r8+cff3)
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                  Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
                  N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
# ifdef PHOSPHORUS
                  cff1=dtdays*SDeRRP(ng)*fac2
                  cff2=1.0_r8/(1.0_r8+cff1)
                  cff3=dtdays*LDeRRP(ng)*fac2
                  cff4=1.0_r8/(1.0_r8+cff3)
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
                  Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
                  P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
# endif
                END DO
              END DO
#else
              cff1=dtdays*SDeRRN(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)
              cff4=1.0_r8/(1.0_r8+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                  Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
                  N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
                END DO
              END DO
# ifdef PHOSPHORUS
              cff1=dtdays*SDeRRP(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRP(ng)
              cff4=1.0_r8/(1.0_r8+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
                  Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
                  P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
                END DO
              END DO
# endif
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
!
              DO k=1,N(ng)
                DO i=Istr,Iend
                  fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
                  fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 depe
                  cff1=dtdays*H2SOR(ng)*fac2
                  cff2=1.0_r8/(1.0_r8+cff1)
                  Bio(i,k,iH2S_)=Bio(i,k,iH2S_)*cff2
                  S_Flux=Bio(i,k,iH2S_)*cff1
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-S_Flux*rOxH2S
                END DO
              END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
              cff1=rho0*550.0_r8
              cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
              k=N(ng)
              DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
                u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
                u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+ &
     &                           (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
                SchmidtN_Ox=1638.0_r8-                                  &
     &                      Bio(i,k,itemp)*(81.83_r8-                   &
     &                                      Bio(i,k,itemp)*             &
     &                                     (1.483_r8-                   &
     &                                      Bio(i,k,itemp)*0.008004_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
                SchmidtN_Ox=1953.4_r8-                                  &
     &                      Bio(i,k,itemp)*(128.0_r8-                   &
     &                                      Bio(i,k,itemp)*             &
     &                                     (3.9918_r8-                  &
     &                                      Bio(i,k,itemp)*0.050091_r8))
# endif

                cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
                TS=LOG((298.15_r8-Bio(i,k,itemp))/                      &
     &                 (273.15_r8+Bio(i,k,itemp)))
                AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+      &
     &             Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
                O2satu=l2mol*EXP(AA)
!
!  Add in O2 gas exchange.
!
                O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
     &                         O2_Flux*Hz_inv(i,k)
              END DO
#endif
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
              DO isink=1,Nsink
                ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    qc(i,k)=Bio(i,k,ibio)
                  END DO
                END DO
!
                DO k=N(ng)-1,1,-1
                  DO i=Istr,Iend
                    FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
                  END DO
                END DO
                DO k=2,N(ng)-1
                  DO i=Istr,Iend
                    dltR=Hz(i,j,k)*FC(i,k)
                    dltL=Hz(i,j,k)*FC(i,k-1)
                    cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                    cffR=cff*FC(i,k)
                    cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                    IF ((dltR*dltL).le.0.0_r8) THEN
                      dltR=0.0_r8
                      dltL=0.0_r8
                    ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                      dltR=cffL
                    ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                      dltL=cffR
                    END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                    cff=(dltR-dltL)*Hz_inv3(i,k)
                    dltR=dltR-cff*Hz(i,j,k+1)
                    dltL=dltL+cff*Hz(i,j,k-1)
                    bR(i,k)=qc(i,k)+dltR
                    bL(i,k)=qc(i,k)-dltL
                    WR(i,k)=(2.0_r8*dltR-dltL)**2
                    WL(i,k)=(dltR-2.0_r8*dltL)**2
                  END DO
                END DO
                cff=1.0E-14_r8
                DO k=2,N(ng)-2
                  DO i=Istr,Iend
                    dltL=MAX(cff,WL(i,k  ))
                    dltR=MAX(cff,WR(i,k+1))
                    bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                    bL(i,k+1)=bR(i,k)
                  END DO
                END DO
                DO i=Istr,Iend
                  FC(i,N(ng))=0.0_r8        ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
                  bR(i,N(ng))=qc(i,N(ng))   ! default strictly monotonic
                  bL(i,N(ng))=qc(i,N(ng))   ! conditions
                  bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
                  bR(i,1)=bL(i,2)
                  bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
                  bR(i,1)=bL(i,2)
                  bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
                  bL(i,2)=qc(i,1)           ! bottom grid boxes are
                  bR(i,1)=qc(i,1)           ! re-assumed to be
                  bL(i,1)=qc(i,1)           ! piecewise constant.
#endif
                END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    dltR=bR(i,k)-qc(i,k)
                    dltL=qc(i,k)-bL(i,k)
                    cffR=2.0_r8*dltR
                    cffL=2.0_r8*dltL
                    IF ((dltR*dltL).lt.0.0_r8) THEN
                      dltR=0.0_r8
                      dltL=0.0_r8
                    ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                      dltR=cffL
                    ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                      dltL=cffR
                    END IF
                    bR(i,k)=qc(i,k)+dltR
                    bL(i,k)=qc(i,k)-dltL
                  END DO
                END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
                cff=dtdays*ABS(Wbio(isink))
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    FC(i,k-1)=0.0_r8
                    WL(i,k)=z_w(i,j,k-1)+cff
                    WR(i,k)=Hz(i,j,k)*qc(i,k)
                    ksource(i,k)=k
                  END DO
                END DO
                DO k=1,N(ng)
                  DO ks=k,N(ng)-1
                    DO i=Istr,Iend
                      IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                        ksource(i,k)=ks+1
                        FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                      END IF
                    END DO
                  END DO
                END DO
!
!  Finalize computation of flux: add fractional part.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    ks=ksource(i,k)
                    cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                    FC(i,k-1)=FC(i,k-1)+                                &
     &                        Hz(i,j,ks)*cu*                            &
     &                        (bL(i,ks)+                                &
     &                         cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-          &
     &                             (1.5_r8-cu)*                         &
     &                             (bR(i,ks)+bL(i,ks)-                  &
     &                              2.0_r8*qc(i,ks))))
                  END DO
                END DO
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bio(i,k,ibio)=qc(i,k)+                              &
     &                            (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                  END DO
                END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
                cff2=4.0_r8/16.0_r8
# ifdef OXYGEN
                cff3=115.0_r8/16.0_r8
                cff4=106.0_r8/16.0_r8
# endif
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
# ifdef DENITRIFICATION
                    Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
#  ifdef OXYGEN
                    Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3
#  endif
# else
                    Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
#  ifdef OXYGEN
                    Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
#  endif
# endif
                  END DO
                END IF
# ifdef PHOSPHORUS
                IF ((ibio.eq.iLDeP).or.                                 &
     &              (ibio.eq.iSDeP)) THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
                    Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
                  END DO
                END IF
                IF (ibio.eq.iPhyt)THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
                    Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*PhyPN(ng)
                  END DO
                END IF
# endif
# if defined H2S && defined OXYGEN
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  DO i=Istr,Iend
                    Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-                      &
     &                             0.5_r8*MIN(Bio(i,1,iOxyg),0.0_r8)
                    Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
                  END DO
                END IF
# endif
#endif
              END DO
            END IF
          END DO
!
!  End of compute basic state arrays I.
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
          DO i=Istr,Iend
            PAR=PARsur(i)
            tl_PAR=tl_PARsur(i)
            AttFac=0.0_r8
            IF (PARsur(i).gt.0.0_r8) THEN              ! day time
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
                Att=(AttSW(ng)+                                         &
     &               AttChl(ng)*Bio(i,k,iChlo)+                         &
     &               AttFac)*                                           &
     &               (z_w(i,j,k)-z_w(i,j,k-1))
                tl_Att=AttChl(ng)*tl_Bio(i,k,iPhyt)*                    &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))+                       &
     &                 (AttSW(ng)+AttChl(ng)*Bio1(i,k,iPhyt)+AttFac)*   &
     &                 (tl_z_w(i,j,k)-tl_z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                tl_ExpAtt=-ExpAtt*tl_Att
                Itop=PAR
                tl_Itop=tl_PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
                tl_PAR=(-tl_Att*PAR+tl_Itop*(1.0_r8-ExpAtt)-            &
     &                  Itop*tl_ExpAtt)/Att
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                cff=PhyCN(ng)*12.0_r8
                Chl2C=MIN(Bio(i,k,iChlo)/(Bio(i,k,iPhyt)*cff+eps),      &
     &                    Chl2C_m(ng))
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                fac1=PAR*PhyIS(ng)
                Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                t_PPmax=Epp*fac1
#ifdef PHOSPHORUS
!
!  Nutrient-limitation terms (Laurent et al. 2012).
!
#else
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
#endif
                cff1=Bio1(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio1(i,k,iNO3_)*K_NO3(ng)
                tl_cff1=tl_Bio(i,k,iNH4_)*K_NH4(ng)
                tl_cff2=tl_Bio(i,k,iNO3_)*K_NO3(ng)
                inhNH4=1.0_r8/(1.0_r8+cff1)
                tl_inhNH4=-inhNH4*inhNH4*tl_cff1
                L_NH4=cff1/(1.0_r8+cff1)
                L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                LTOT=L_NO3+L_NH4
#ifdef PHOSPHORUS
                cff3=Bio1(i,k,iPO4_)*K_PO4(ng)
                tl_cff3=tl_Bio(i,k,iPO4_)*K_PO4(ng)
                L_PO4=cff3/(1.0_r8+cff3)
!
!  Nitrate, ammonium and phosphate uptake by Phytoplankton.
!
#else
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
#endif
                fac1=dtdays*t_PPmax
                cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*Bio1(i,k,iPhyt)
                cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio1(i,k,iPhyt)
                tl_cff4=(fac1*K_NO3(ng)*(tl_inhNH4*Bio(i,k,iPhyt)+      &
     &                                   inhNH4*tl_Bio(i,k,iPhyt))-     &
     &                   cff4*tl_cff2)/(1.0_r8+cff2)
                tl_cff5=(fac1*K_NH4(ng)*tl_Bio(i,k,iPhyt)-              &
     &                   cff5*tl_cff1)/(1.0_r8+cff1)

                N_Flux_NewProd=Bio1(i,k,iNO3_)/(1.0_r8+cff4)*cff4
                N_Flux_RegProd=Bio1(i,k,iNH4_)/(1.0_r8+cff5)*cff5
                N_Flux_SumProd=N_Flux_NewProd+N_Flux_RegProd
                tl_N_Flux_NewProd=(tl_Bio(i,k,iNO3_)*cff4+              &
     &                             Bio1(i,k,iNO3_)*tl_cff4-             &
     &                             N_Flux_NewProd*tl_cff4)/(1.0_r8+cff4)
                tl_N_Flux_RegProd=(tl_Bio(i,k,iNH4_)*cff5+              &
     &                             Bio1(i,k,iNH4_)*tl_cff5-             &
     &                             N_Flux_RegProd*tl_cff5)/(1.0_r8+cff5)
                tl_N_Flux_SumProd=tl_N_Flux_NewProd+tl_N_Flux_RegProd
#ifdef PHOSPHORUS
                IF (LTOT.lt.L_PO4) THEN
!>                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
!>                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
!>                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
!>                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-                        &
!>   &                           PhyPN(ng)*N_Flux_SumProd
                  tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-                 &
     &                               Bio(i,k,iNO3_)*tl_cff4)/           &
     &                              (1.0_r8+cff4)
                  tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-                 &
     &                               Bio(i,k,iNH4_)*tl_cff5)/           &
     &                              (1.0_r8+cff5)
                  tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+tl_N_Flux_SumProd
                  tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)-                  &
     &                              PhyPN(ng)*tl_N_Flux_SumProd
                ELSE IF (L_PO4.lt.LTOT) THEN
                  cff4=fac1*K_PO4(ng)/(1.0_r8+cff3)*Bio1(i,k,iPhyt)
                  tl_cff4=(fac1*K_PO4(ng)*tl_Bio(i,k,iPhyt)-            &
     &                     cff4*tl_cff3)/(1.0_r8+cff3)
!>                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff4)
                  tl_Bio(i,k,iPO4_)=(tl_Bio(i,k,iPO4_)-                 &
     &                               Bio(i,k,iPO4_)*tl_cff4)/           &
     &                              (1.0_r8+cff4)
                  P_Flux_SumProd=Bio1(i,k,iPO4_)*cff4
                  fac1=MIN(P_Flux_SumProd/PhyPN(ng)/N_Flux_SumProd,     &
     &                     1.0_r8)
!>                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd*fac1
!>                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd*fac1
!>                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd*fac1
                  tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+                  &
     &                              tl_N_Flux_SumProd*fac1
                  tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)-                  &
     &                              tl_N_Flux_NewProd*fac1
                  tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)-                  &
     &                              tl_N_Flux_RegProd*fac1
                  LTOT=L_PO4
                ENDIF
#else
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
!>              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
                tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-                   &
     &                             Bio(i,k,iNO3_)*tl_cff4)/(1.0_r8+cff4)
                tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-                   &
     &                             Bio(i,k,iNH4_)*tl_cff5)/(1.0_r8+cff5)
                tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+tl_N_Flux_SumProd
#endif
!
!>              Bio(i,k,iChlo)=Bio(i,k,iChlo)+                          &
!>   &                         (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*       &
!>   &                          Chl2C_m(ng)*Bio(i,k,iChlo))/            &
!>   &                         (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
                tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)+                    &
     &                            (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*    &
     &                             Chl2C_m(ng)*tl_Bio(i,k,iChlo))/      &
     &                            (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
#ifdef OXYGEN
!>              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
!>   &                         N_Flux_NewProd*rOxNO3+                   &
!>   &                         N_Flux_RegProd*rOxNH4
                tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)+                    &
     &                            tl_N_Flux_NewProd*rOxNO3+             &
     &                            tl_N_Flux_RegProd*rOxNH4
#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
#ifdef OXYGEN
                fac2=MAX(Bio1(i,k,iOxyg),0.0_r8)
                fac3=MAX(fac2/(K_Nitri(ng)+fac2),0.0_r8)
                fac1=dtdays*NitriR(ng)*fac3
#else
                fac1=dtdays*NitriR(ng)
#endif
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                cff2=1.0_r8-MAX(0.0_r8,cff1)
                cff3=fac1*cff2
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)/(1.0_r8+cff3)
!>              N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                tl_N_Flux_Nitrifi=tl_Bio(i,k,iNH4_)*cff3
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
                tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+tl_N_Flux_Nitrifi
#ifdef OXYGEN
!>              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
                tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                    &
     &                            2.0_r8*tl_N_Flux_Nitrifi
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                PAR=Itop*ExpAtt
                tl_PAR=tl_Itop*ExpAtt+Itop*tl_ExpAtt
              END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
            ELSE                                       ! night time
              cff3=dtdays*NitriR(ng)
              DO k=N(ng),1,-1
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)/(1.0_r8+cff3)
!>              N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                tl_N_Flux_Nitrifi=tl_Bio(i,k,iNH4_)*cff3
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
                tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+tl_N_Flux_Nitrifi
#ifdef OXYGEN
!>              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
                tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                    &
     &                            2.0_r8*tl_N_Flux_Nitrifi
#endif
              END DO
            END IF
          END DO
!
!  Compute appropriate basic state arrays II.
!
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
                Bio(i,k,ibio)=Bio_old(i,k,ibio)
              END DO
            END DO
          END DO
!
!  Extract potential temperature and salinity.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
              Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              AttFac=0.0_r8
              IF (PARsur(i).gt.0.0_r8) THEN
                DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
                  Att=(AttSW(ng)+                                       &
     &                 AttChl(ng)*Bio(i,k,iChlo)+                       &
     &                 AttFac)*                                         &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell cent
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                  cff=PhyCN(ng)*12.0_r8
                  Chl2C=MIN(Bio(i,k,iChlo)/(Bio(i,k,iPhyt)*cff+eps),    &
     &                      Chl2C_m(ng))
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                  Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                  fac1=PAR*PhyIS(ng)
                  Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                  t_PPmax=Epp*fac1
#ifdef PHOSPHORUS
!
!  Nutrient-limitation terms (Laurent et al. 2012).
!
#else
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
#endif
                  cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                  cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                  inhNH4=1.0_r8/(1.0_r8+cff1)
                  L_NH4=cff1/(1.0_r8+cff1)
                  L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                  LTOT=L_NO3+L_NH4
#ifdef PHOSPHORUS
                  cff3=Bio(i,k,iPO4_)*K_PO4(ng)
                  L_PO4=cff3/(1.0_r8+cff3)
!
!  Nitrate, ammonium and phosphate uptake by Phytoplankton.
!
#else
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
#endif
                  fac1=dtdays*t_PPmax
                  cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*             &
     &                 Bio(i,k,iPhyt)
                  cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio(i,k,iPhyt)
                  N_Flux_NewProd=Bio(i,k,iNO3_)/(1.0_r8+cff4)*cff4
                  N_Flux_RegProd=Bio(i,k,iNH4_)/(1.0_r8+cff5)*cff5
                  N_Flux_SumProd=N_Flux_NewProd+N_Flux_RegProd
#ifdef PHOSPHORUS
                  IF (LTOT.lt.L_PO4) THEN
                    Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                    Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                    Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
                    Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-                      &
     &                             PhyPN(ng)*N_Flux_SumProd
                  ELSE IF (L_PO4.lt.LTOT) THEN
                    cff4=fac1*K_PO4(ng)/(1.0_r8+cff3)*Bio(i,k,iPhyt)
                    Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff4)
                    P_Flux_SumProd=Bio(i,k,iPO4_)*cff4
                    fac1=MIN(P_Flux_SumProd/PhyPN(ng)/N_Flux_SumProd,   &
     &                       1.0_r8)
                    Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd*fac1
                    Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd*fac1
                    Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd*fac1
                    LTOT=L_PO4
                  ENDIF
#else
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
#endif
!
                  Bio(i,k,iChlo)=Bio(i,k,iChlo)+                        &
     &                           (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*     &
     &                            Chl2C_m(ng)*Bio(i,k,iChlo))/          &
     &                           (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                        &
     &                           N_Flux_NewProd*rOxNO3+                 &
     &                           N_Flux_RegProd*rOxNH4
#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
#ifdef OXYGEN
                  fac2=MAX(Bio(i,k,iOxyg),0.0_r8)     ! O2 max
                  fac3=MAX(fac2/(K_Nitri(ng)+fac2),0.0_r8) ! MM for O2 d
                  fac1=dtdays*NitriR(ng)*fac3
#else
                  fac1=dtdays*NitriR(ng)
#endif
                  cff1=(PAR-I_thNH4(ng))/                               &
     &                 (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                  cff2=1.0_r8-MAX(0.0_r8,cff1)
                  cff3=fac1*cff2
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                  N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
              ELSE
                cff3=dtdays*NitriR(ng)
                DO k=N(ng),1,-1
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                  N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
                END DO
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
!
            fac1=dtdays*ZooGR(ng)
            cff2=dtdays*PhyMR(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
!
                cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/                &
     &               (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                cff3=1.0_r8/(1.0_r8+cff1)
                Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
                Bio1(i,k,iChlo)=Bio(i,k,iChlo)
                Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
                Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
                N_Flux_Assim=cff1*Bio(i,k,iPhyt)*ZooAE_N(ng)
                N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
                Bio1(i,k,iZoop)=Bio(i,k,iZoop)
                Bio1(i,k,iSDeN)=Bio(i,k,iSDeN)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)+                          &
     &                         N_Flux_Assim
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                          &
     &                         N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
                N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),      &
     &                                  0.0_r8)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
                Bio(i,k,iChlo)=Bio(i,k,iChlo)-                          &
     &                         cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),      &
     &                                  0.0_r8)
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                          &
     &                         N_Flux_Pmortal
#ifdef PHOSPHORUS
                Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                          &
     &                         PhyPN(ng)*(N_Flux_Egest+                 &
     &                                    N_Flux_Pmortal)+              &
     &                         (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
#endif
              END DO
            END DO
#ifdef OXYGEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
              END DO
            END DO
#endif
            IF (Iteradj.ne.Iter) THEN
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
              cff1=dtdays*ZooBM(ng)
              fac2=dtdays*ZooMR(ng)
              fac3=dtdays*ZooER(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/              &
     &                 (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                  cff2=fac2*Bio(i,k,iZoop)
                  cff3=fac1*ZooAE_N(ng)
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)/                        &
     &                           (1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
                  N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
                  N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
                  N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),    &
     &                                    0.0_r8)
                  Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                        &
     &                           rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
#endif
#ifdef PHOSPHORUS
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                        &
     &                           ZooPN(ng)*(N_Flux_Zmetabo+             &
     &                                      N_Flux_Zexcret)
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                        &
     &                           ZooPN(ng)*N_Flux_Zmortal
#endif
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
              fac1=dtdays*CoagR(ng)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  cff1=fac1*(Bio(i,k,iSDeN)+Bio(i,k,iPhyt))
                  cff2=1.0_r8/(1.0_r8+cff1)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
                  Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                  N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
                  N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
                  Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                        &
     &                           N_Flux_CoagP+N_Flux_CoagD
#ifdef PHOSPHORUS
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
                  Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                        &
     &                           PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
#endif
                END DO
              END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
#ifdef OXYGEN
              DO k=1,N(ng)
                DO i=Istr,Iend
                  fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
                  fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 depe
                  cff1=dtdays*SDeRRN(ng)*fac2
                  cff2=1.0_r8/(1.0_r8+cff1)
                  cff3=dtdays*LDeRRN(ng)*fac2
                  cff4=1.0_r8/(1.0_r8+cff3)
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                  Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
                  N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
# ifdef PHOSPHORUS
                  cff1=dtdays*SDeRRP(ng)*fac2
                  cff2=1.0_r8/(1.0_r8+cff1)
                  cff3=dtdays*LDeRRP(ng)*fac2
                  cff4=1.0_r8/(1.0_r8+cff3)
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
                  Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
                  P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
# endif
                END DO
              END DO
#else
              cff1=dtdays*SDeRRN(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)
              cff4=1.0_r8/(1.0_r8+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                  Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
                  N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
                END DO
              END DO
# ifdef PHOSPHORUS
              cff1=dtdays*SDeRRP(ng)
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRP(ng)
              cff4=1.0_r8/(1.0_r8+cff3)
              DO k=1,N(ng)
                DO i=Istr,Iend
                  Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
                  Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
                  P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
                END DO
              END DO
# endif
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
!
              DO k=1,N(ng)
                DO i=Istr,Iend
                  fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
                  fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 depe
                  cff1=dtdays*H2SOR(ng)*fac2
                  cff2=1.0_r8/(1.0_r8+cff1)
                  Bio(i,k,iH2S_)=Bio(i,k,iH2S_)*cff2
                  S_Flux=Bio(i,k,iH2S_)*cff1
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-S_Flux*rOxH2S
                END DO
              END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
              cff1=rho0*550.0_r8
              cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
              k=N(ng)
              DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
                u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
                u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+ &
     &                           (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
                SchmidtN_Ox=1638.0_r8-                                  &
     &                      Bio(i,k,itemp)*(81.83_r8-                   &
     &                                      Bio(i,k,itemp)*             &
     &                                     (1.483_r8-                   &
     &                                      Bio(i,k,itemp)*0.008004_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
                SchmidtN_Ox=1953.4_r8-                                  &
     &                      Bio(i,k,itemp)*(128.0_r8-                   &
     &                                      Bio(i,k,itemp)*             &
     &                                     (3.9918_r8-                  &
     &                                      Bio(i,k,itemp)*0.050091_r8))
# endif

                cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
                TS=LOG((298.15_r8-Bio(i,k,itemp))/                      &
     &                 (273.15_r8+Bio(i,k,itemp)))
                AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+      &
     &             Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
                O2satu=l2mol*EXP(AA)
!
!  Add in O2 gas exchange.
!
                O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
     &                         O2_Flux*Hz_inv(i,k)
              END DO
#endif
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
              DO isink=1,Nsink
                ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    qc(i,k)=Bio(i,k,ibio)
                  END DO
                END DO
!
                DO k=N(ng)-1,1,-1
                  DO i=Istr,Iend
                    FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
                  END DO
                END DO
                DO k=2,N(ng)-1
                  DO i=Istr,Iend
                    dltR=Hz(i,j,k)*FC(i,k)
                    dltL=Hz(i,j,k)*FC(i,k-1)
                    cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                    cffR=cff*FC(i,k)
                    cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                    IF ((dltR*dltL).le.0.0_r8) THEN
                      dltR=0.0_r8
                      dltL=0.0_r8
                    ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                      dltR=cffL
                    ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                      dltL=cffR
                    END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                    cff=(dltR-dltL)*Hz_inv3(i,k)
                    dltR=dltR-cff*Hz(i,j,k+1)
                    dltL=dltL+cff*Hz(i,j,k-1)
                    bR(i,k)=qc(i,k)+dltR
                    bL(i,k)=qc(i,k)-dltL
                    WR(i,k)=(2.0_r8*dltR-dltL)**2
                    WL(i,k)=(dltR-2.0_r8*dltL)**2
                  END DO
                END DO
                cff=1.0E-14_r8
                DO k=2,N(ng)-2
                  DO i=Istr,Iend
                    dltL=MAX(cff,WL(i,k  ))
                    dltR=MAX(cff,WR(i,k+1))
                    bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                    bL(i,k+1)=bR(i,k)
                  END DO
                END DO
                DO i=Istr,Iend
                  FC(i,N(ng))=0.0_r8        ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
                  bR(i,N(ng))=qc(i,N(ng))   ! default strictly monotonic
                  bL(i,N(ng))=qc(i,N(ng))   ! conditions
                  bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
                  bR(i,1)=bL(i,2)
                  bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
                  bR(i,1)=bL(i,2)
                  bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
                  bL(i,2)=qc(i,1)           ! bottom grid boxes are
                  bR(i,1)=qc(i,1)           ! re-assumed to be
                  bL(i,1)=qc(i,1)           ! piecewise constant.
#endif
                END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    dltR=bR(i,k)-qc(i,k)
                    dltL=qc(i,k)-bL(i,k)
                    cffR=2.0_r8*dltR
                    cffL=2.0_r8*dltL
                    IF ((dltR*dltL).lt.0.0_r8) THEN
                      dltR=0.0_r8
                      dltL=0.0_r8
                    ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                      dltR=cffL
                    ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                      dltL=cffR
                    END IF
                    bR(i,k)=qc(i,k)+dltR
                    bL(i,k)=qc(i,k)-dltL
                  END DO
                END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
                cff=dtdays*ABS(Wbio(isink))
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    FC(i,k-1)=0.0_r8
                    WL(i,k)=z_w(i,j,k-1)+cff
                    WR(i,k)=Hz(i,j,k)*qc(i,k)
                    ksource(i,k)=k
                  END DO
                END DO
                DO k=1,N(ng)
                  DO ks=k,N(ng)-1
                    DO i=Istr,Iend
                      IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                        ksource(i,k)=ks+1
                        FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                      END IF
                    END DO
                  END DO
                END DO
!
!  Finalize computation of flux: add fractional part.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    ks=ksource(i,k)
                    cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                    FC(i,k-1)=FC(i,k-1)+                                &
     &                        Hz(i,j,ks)*cu*                            &
     &                        (bL(i,ks)+                                &
     &                         cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-          &
     &                             (1.5_r8-cu)*                         &
     &                             (bR(i,ks)+bL(i,ks)-                  &
     &                              2.0_r8*qc(i,ks))))
                  END DO
                END DO
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bio(i,k,ibio)=qc(i,k)+                              &
     &                            (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                  END DO
                END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
                cff2=4.0_r8/16.0_r8
# ifdef OXYGEN
                cff3=115.0_r8/16.0_r8
                cff4=106.0_r8/16.0_r8
# endif
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
# ifdef DENITRIFICATION
                    Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
#  ifdef OXYGEN
                    Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3
#  endif
# else
                    Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
#  ifdef OXYGEN
                    Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
#  endif
# endif
                  END DO
                END IF
# ifdef PHOSPHORUS
                IF ((ibio.eq.iLDeP).or.                                 &
     &              (ibio.eq.iSDeP)) THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
                    Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
                  END DO
                END IF
                IF (ibio.eq.iPhyt)THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
                    Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*PhyPN(ng)
                  END DO
                END IF
# endif
# if defined H2S && defined OXYGEN
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  DO i=Istr,Iend
                    Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-                      &
     &                             0.5_r8*MIN(Bio(i,1,iOxyg),0.0_r8)
                    Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
                  END DO
                END IF
# endif
#endif
              END DO
            END IF
          END DO
!
!  End of compute basic state arrays II.
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
!
          fac1=dtdays*ZooGR(ng)
          cff2=dtdays*PhyMR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
!
              cff1=fac1*Bio1(i,k,iZoop)*Bio1(i,k,iPhyt)/                &
     &             (K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              tl_cff1=(fac1*(tl_Bio(i,k,iZoop)*Bio1(i,k,iPhyt)+         &
     &                       Bio1(i,k,iZoop)*tl_Bio(i,k,iPhyt))-        &
     &                 2.0_r8*Bio1(i,k,iPhyt)*tl_Bio(i,k,iPhyt)*cff1/   &
     &                (K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              cff3=1.0_r8/(1.0_r8+cff1)
              tl_cff3=-cff3*cff3*tl_cff1
!>            Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
!>            Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
              tl_Bio(i,k,iPhyt)=tl_cff3*Bio1(i,k,iPhyt)+                &
     &                          cff3*tl_Bio(i,k,iPhyt)
              tl_Bio(i,k,iChlo)=tl_cff3*Bio1(i,k,iChlo)+                &
     &                          cff3*tl_Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
!>            N_Flux_Assim=cff1*Bio(i,k,iPhyt)*ZooAE_N(ng)
              tl_N_Flux_Assim=(tl_cff1*Bio1(i,k,iPhyt)+                 &
     &                         cff1*tl_Bio(i,k,iPhyt))*ZooAE_N(ng)
!>            N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
              tl_N_Flux_Egest=(tl_Bio(i,k,iPhyt)*cff1+                  &
     &                         Bio1(i,k,iPhyt)*tl_cff1)*                &
     &                        (1.0_r8-ZooAE_N(ng))
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest
              tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)+tl_N_Flux_Assim
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
!>            N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              tl_N_Flux_Pmortal=tl_cff2*MAX(Bio1(i,k,iPhyt)-PhyMin(ng), &
     &                                      0.0_r8)+                    &
     &                          cff2*(0.5_r8+SIGN(0.5_r8,               &
     &                                            Bio1(i,k,iPhyt)-      &
     &                                            PhyMin(ng)))*         &
     &                               tl_Bio(i,k,iPhyt)
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)-                            &
!>   &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal
              tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)-tl_N_Flux_Pmortal
              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)-                      &
     &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Pmortal
#ifdef PHOSPHORUS
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
!>   &                       PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)+   &
!>   &                       (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+                      &
     &                          PhyPN(ng)*(tl_N_Flux_Egest+             &
     &                                     tl_N_Flux_Pmortal)+          &
     &                          (PhyPN(ng)-ZooPN(ng))*tl_N_Flux_Assim
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
          cff1=dtdays*ZooBM(ng)
          fac2=dtdays*ZooMR(ng)
          fac3=dtdays*ZooER(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=fac3*Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt)/                &
     &             (K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              tl_fac1=(fac3-fac1)*                                      &
     &                2.0_r8*Bio1(i,k,iPhyt)*tl_Bio(i,k,iPhyt)/         &
     &                (K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              cff2=fac2*Bio1(i,k,iZoop)
              tl_cff2=fac2*tl_Bio(i,k,iZoop)
              cff3=fac1*ZooAE_N(ng)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)/                            &
!>   &                       (1.0_r8+cff2+cff3)
              tl_Bio(i,k,iZoop)=(tl_Bio(i,k,iZoop)-                     &
     &                           Bio(i,k,iZoop)*tl_cff2)/               &
     &                          (1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
!>            N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
!>            N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
              tl_N_Flux_Zmortal=tl_cff2*Bio1(i,k,iZoop)+                &
     &                          cff2*tl_Bio(i,k,iZoop)
              tl_N_Flux_Zexcret=cff3*tl_Bio(i,k,iZoop)
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Zexcret
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Zmortal
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
!>            N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              tl_N_Flux_Zmetabo=cff1*(0.5_r8+SIGN(0.5_r8,               &
                                Bio1(i,k,iZoop)-ZooMin(ng)))*           &
                                tl_Bio(i,k,iZoop)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
              tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)-tl_N_Flux_Zmetabo
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Zmetabo
#ifdef OXYGEN
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
!>   &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                      &
     &                          rOxNH4*(tl_N_Flux_Zmetabo+              &
     &                                  tl_N_Flux_Zexcret)
#endif
#ifdef PHOSPHORUS
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                            &
!>   &                       ZooPN(ng)*(N_Flux_Zmetabo+N_Flux_Zexcret)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
!>   &                       ZooPN(ng)*N_Flux_Zmortal
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+                      &
     &                          ZooPN(ng)*(tl_N_Flux_Zmetabo+           &
     &                                     tl_N_Flux_Zexcret)
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+                      &
     &                          ZooPN(ng)*tl_N_Flux_Zmortal
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*(Bio1(i,k,iSDeN)+Bio1(i,k,iPhyt))
              tl_cff1=fac1*(tl_Bio(i,k,iSDeN)+tl_Bio(i,k,iPhyt))
              cff2=1.0_r8/(1.0_r8+cff1)
              tl_cff2=-tl_cff1*cff2*cff2
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)*cff2+                 &
     &                          Bio1(i,k,iPhyt)*tl_cff2
              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)*cff2+                 &
     &                          Bio1(i,k,iChlo)*tl_cff2
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2+                 &
     &                          Bio1(i,k,iSDeN)*tl_cff2
!>            N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
!>            N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
              tl_N_Flux_CoagP=tl_Bio(i,k,iPhyt)*cff1+                   &
     &                        Bio1(i,k,iPhyt)*tl_cff1
              tl_N_Flux_CoagD=tl_Bio(i,k,iSDeN)*cff1+                   &
     &                        Bio1(i,k,iSDeN)*tl_cff1
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
!>   &                       N_Flux_CoagP+N_Flux_CoagD
              tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)+                      &
     &                          tl_N_Flux_CoagP+tl_N_Flux_CoagD
#ifdef PHOSPHORUS
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                            &
!>   &                       PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)-                      &
     &                          PhyPN(ng)*tl_N_Flux_CoagD
              tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)+                      &
     &                          PhyPN(ng)*(tl_N_Flux_CoagP+             &
     &                                     tl_N_Flux_CoagD)
#endif
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
#ifdef OXYGEN
          DO k=1,N(ng)
            DO i=Istr,Iend
!
!  (okada) 
!
              fac1=MAX(Bio1(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
              fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 dependen
              cff1=dtdays*SDeRRN(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRN(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2
              tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)*cff4
!>            N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              tl_N_Flux_Remine=tl_Bio(i,k,iSDeN)*cff1+                  &
     &                         tl_Bio(i,k,iLDeN)*cff3
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Remine
              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-tl_N_Flux_Remine*     &
     &                                            rOxNH4
# ifdef PHOSPHORUS
              cff1=dtdays*SDeRRP(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
              cff3=dtdays*LDeRRP(ng)*fac2
              cff4=1.0_r8/(1.0_r8+cff3)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)*cff2
              tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)*cff4
!>            P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              tl_P_Flux_Remine=tl_Bio(i,k,iSDeP)*cff1+                  &
     &                         tl_Bio(i,k,iLDeP)*cff3
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+tl_P_Flux_Remine
# endif
            END DO
          END DO
#else
          cff1=dtdays*SDeRRN(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRN(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2
              tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)*cff4
!>            N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              tl_N_Flux_Remine=tl_Bio(i,k,iSDeN)*cff1+                  &
     &                         tl_Bio(i,k,iLDeN)*cff3
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Remine
            END DO
          END DO
# ifdef PHOSPHORUS
          cff1=dtdays*SDeRRP(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRP(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)*cff2
              tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)*cff4
!>            P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              tl_P_Flux_Remine=tl_Bio(i,k,iSDeP)*cff1+                  &
     &                         tl_Bio(i,k,iLDeP)*cff3
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+tl_P_Flux_Remine
            END DO
          END DO
# endif
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=MAX(Bio1(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
              fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 dependen
              cff1=dtdays*H2SOR(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
!>            Bio(i,k,iH2S_)=Bio(i,k,iH2S_)*cff2
              tl_Bio(i,k,iH2S_)=tl_Bio(i,k,iH2S_)*cff2
!>            S_Flux=Bio(i,k,iH2S_)*cff1
              tl_S_Flux=tl_Bio(i,k,iH2S_)*cff1
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-S_Flux*rOxH2S
              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-tl_S_Flux*rOxH2S
            END DO
          END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
          cff1=rho0*550.0_r8
          cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
          k=N(ng)
          DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
            u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
            u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+     &
     &                       (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
            SchmidtN_Ox=1638.0_r8-                                      &
     &                  Bio(i,k,itemp)*(81.83_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (1.483_r8-                      &
     &                                   Bio(i,k,itemp)*0.008004_r8))
            tl_SchmidtN_Ox=-tl_Bio(i,k,itemp)*                          &
     &                     (81.83_r8-Bio(i,k,itemp)*                    &
     &                     (2.966_r8-Bio(i,k,itemp)*0.024012_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
            SchmidtN_Ox=1953.4_r8-                                      &
     &                  Bio(i,k,itemp)*(128.0_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (3.9918_r8-                     &
     &                                   Bio(i,k,itemp)*0.050091_r8))
            tl_SchmidtN_Ox=Bio(i,k,itemp)*                              &
                           (128.0_r8-Bio(i,k,itemp)*                    &
     &                     (7.9836_r8-Bio(i,k,itemp)*0.150273_r8))
# endif

            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
            tl_cff3=-0.5_r8*cff3*tl_SchmidtN_Ox/SchmidtN_Ox
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
            TS=LOG((298.15_r8-Bio(i,k,itemp))/                          &
     &             (273.15_r8+Bio(i,k,itemp)))
            tl_TS=-tl_Bio(i,k,itemp)/(298.15_r8-Bio(i,k,itemp))-        &
     &             tl_Bio(i,k,itemp)/(273.15_r8+Bio(i,k,itemp))
            AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+          &
     &             Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+       &
     &             OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
            tl_AA=tl_TS*(OA1+TS*(2.0_r8*OA2+TS*(3.0_r8*OA3+TS*          &
     &            (4.0_r8*OA4+TS*5.0_r8*OA5))))+                        &
     &            tl_Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+     &
     &            Bio(i,k,isalt)*(tl_TS*(OB1+TS*(2.0_r8*OB2+            &
     &                                           TS*3.0_r8*OB3)))+      &
     &            OC0*2.0_r8*Bio(i,k,isalt)*tl_Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
            O2satu=l2mol*EXP(AA)
            tl_O2satu=l2mol*EXP(AA)*tl_AA
!
!  Add in O2 gas exchange.
!
!>          O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
            tl_O2_Flux=tl_cff3*(O2satu-Bio1(i,k,iOxyg))+                &
     &                 cff3*(tl_O2satu-tl_Bio(i,k,iOxyg))
!>          Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
!>   &                     O2_Flux*Hz_inv(i,k)
            tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)+                        &
     &                        tl_O2_Flux*Hz_inv(i,k)+                   &
     &                        O2_Flux*tl_Hz_inv(i,k)
          END DO
#endif
!
!  Compute appropriate basic state arrays III.
!
          DO itrc=1,NBT
            ibio=idbio(itrc)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio_old(i,k,ibio)=MAX(0.0_r8,t(i,j,k,nstp,ibio))
                Bio(i,k,ibio)=Bio_old(i,k,ibio)
              END DO
            END DO
          END DO
!
!  Extract potential temperature and salinity.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              Bio(i,k,itemp)=MIN(t(i,j,k,nstp,itemp),35.0_r8)
              Bio(i,k,isalt)=MAX(t(i,j,k,nstp,isalt), 0.0_r8)
            END DO
          END DO
!
!  Calculate surface Photosynthetically Available Radiation (PAR).  The
!  net shortwave radiation is scaled back to Watts/m2 and multiplied by
!  the fraction that is photosynthetically available, PARfrac.
!
          DO i=Istr,Iend
            PARsur(i)=PARfrac(ng)*srflx(i,j)*rho0*Cp
          END DO
!
!=======================================================================
!  Start internal iterations to achieve convergence of the nonlinear
!  backward-implicit solution.
!=======================================================================
!
          DO Iteradj=1,Iter
!
!-----------------------------------------------------------------------
!  Light-limited computations.
!-----------------------------------------------------------------------
!
!  Compute attenuation coefficient based on the concentration of
!  chlorophyll-a within each grid box.  Then, attenuate surface
!  photosynthetically available radiation (PARsur) down inot the
!  water column.  Thus, PAR at certain depth depends on the whole
!  distribution of chlorophyll-a above.
!  To compute rate of maximum primary productivity (t_PPmax), one needs
!  PAR somewhat in the middle of the gridbox, so that attenuation "Att"
!  corresponds to half of the grid box height, while PAR is multiplied
!  by it twice: once to get it in the middle of grid-box and once the
!  compute on the lower grid-box interface.
!
            DO i=Istr,Iend
              PAR=PARsur(i)
              AttFac=0.0_r8
              IF (PARsur(i).gt.0.0_r8) THEN
                DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
                  Att=(AttSW(ng)+                                       &
     &                 AttChl(ng)*Bio(i,k,iChlo)+                       &
     &                 AttFac)*                                         &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))
                  ExpAtt=EXP(-Att)
                  Itop=PAR
                  PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell cent
!
!  Compute Chlorophyll-a phytoplankton ratio, [mg Chla / (mg C)].
!
                  cff=PhyCN(ng)*12.0_r8
                  Chl2C=MIN(Bio(i,k,iChlo)/(Bio(i,k,iPhyt)*cff+eps),    &
     &                      Chl2C_m(ng))
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                  Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                  fac1=PAR*PhyIS(ng)
                  Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                  t_PPmax=Epp*fac1
#ifdef PHOSPHORUS
!
!  Nutrient-limitation terms (Laurent et al. 2012).
!
#else
!
!  Nutrient-limitation terms (Parker 1993 Ecol Mod., 66, 113-120).
!
#endif
                  cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                  cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                  inhNH4=1.0_r8/(1.0_r8+cff1)
                  L_NH4=cff1/(1.0_r8+cff1)
                  L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                  LTOT=L_NO3+L_NH4
#ifdef PHOSPHORUS
                  cff3=Bio(i,k,iPO4_)*K_PO4(ng)
                  L_PO4=cff3/(1.0_r8+cff3)
!
!  Nitrate, ammonium and phosphate uptake by Phytoplankton.
!
#else
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
#endif
                  fac1=dtdays*t_PPmax
                  cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)*             &
     &                 Bio(i,k,iPhyt)
                  cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)*Bio(i,k,iPhyt)
                  N_Flux_NewProd=Bio(i,k,iNO3_)/(1.0_r8+cff4)*cff4
                  N_Flux_RegProd=Bio(i,k,iNH4_)/(1.0_r8+cff5)*cff5
                  N_Flux_SumProd=N_Flux_NewProd+N_Flux_RegProd
#ifdef PHOSPHORUS
                  IF (LTOT.lt.L_PO4) THEN
                    Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                    Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                    Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
                    Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-                      &
     &                             PhyPN(ng)*N_Flux_SumProd
                  ELSE IF (L_PO4.lt.LTOT) THEN
                    cff4=fac1*K_PO4(ng)/(1.0_r8+cff3)*Bio(i,k,iPhyt)
                    Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff4)
                    P_Flux_SumProd=Bio(i,k,iPO4_)*cff4
                    fac1=MIN(P_Flux_SumProd/PhyPN(ng)/N_Flux_SumProd,   &
     &                       1.0_r8)
                    Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd*fac1
                    Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd*fac1
                    Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd*fac1
                    LTOT=L_PO4
                  ENDIF
#else
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+N_Flux_SumProd
#endif
!
                  Bio(i,k,iChlo)=Bio(i,k,iChlo)+                        &
     &                           (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*     &
     &                            Chl2C_m(ng)*Bio(i,k,iChlo))/          &
     &                           (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                        &
     &                           N_Flux_NewProd*rOxNO3+                 &
     &                           N_Flux_RegProd*rOxNH4
#endif
!
! The Nitrification of NH4 ==> NO3 is thought to occur only in dark and
! only in aerobic water (see Olson, R. J., 1981, JMR: (39), 227-238.).
!
!         NH4+ + 3/2 O2  ==> NO2- + H2O;  via Nitrosomonas bacteria
!         NO2-  + 1/2 O2 ==> NO3-      ;  via Nitrobacter  bacteria
!
! Note that the entire process has a total loss of two moles of O2 per
! mole of NH4. If we were to resolve NO2 profiles, this is where we
! would change the code to split out the differential effects of the
! two different bacteria types. If OXYGEN is defined, nitrification is
! inhibited at low oxygen concentrations using a Michaelis-Menten term.
!
#ifdef OXYGEN
                  fac2=MAX(Bio(i,k,iOxyg),0.0_r8)     ! O2 max
                  fac3=MAX(fac2/(K_Nitri(ng)+fac2),0.0_r8) ! MM for O2 d
                  fac1=dtdays*NitriR(ng)*fac3
#else
                  fac1=dtdays*NitriR(ng)
#endif
                  cff1=(PAR-I_thNH4(ng))/                               &
     &                 (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                  cff2=1.0_r8-MAX(0.0_r8,cff1)
                  cff3=fac1*cff2
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                  N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
!
!  Light attenuation at the bottom of the grid cell. It is the starting
!  PAR value for the next (deeper) vertical grid cell.
!
                  PAR=Itop*ExpAtt
                END DO
!
!  If PARsur=0, nitrification occurs at the maximum rate (NitriR).
!
              ELSE
                cff3=dtdays*NitriR(ng)
                DO k=N(ng),1,-1
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                  N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                  Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
#endif
                END DO
              END IF
            END DO
!
!-----------------------------------------------------------------------
!  Phytoplankton grazing by zooplankton (rate: ZooGR), phytoplankton
!  assimilated to zooplankton (fraction: ZooAE_N) and egested to small
!  detritus, and phytoplankton mortality (rate: PhyMR) to small
!  detritus. [Landry 1993 L&O 38:468-472]
!-----------------------------------------------------------------------
!
            fac1=dtdays*ZooGR(ng)
            cff2=dtdays*PhyMR(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
!
! Phytoplankton grazing by zooplankton.
!
                cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/                &
     &               (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                cff3=1.0_r8/(1.0_r8+cff1)
                Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
                Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
                N_Flux_Assim=cff1*Bio(i,k,iPhyt)*ZooAE_N(ng)
                N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
                Bio(i,k,iZoop)=Bio(i,k,iZoop)+                          &
     &                         N_Flux_Assim
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                          &
     &                         N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
                N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),      &
     &                                  0.0_r8)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
                Bio(i,k,iChlo)=Bio(i,k,iChlo)-                          &
     &                         cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),      &
     &                                  0.0_r8)
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+                          &
     &                         N_Flux_Pmortal
#ifdef PHOSPHORUS
                Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                          &
     &                         PhyPN(ng)*(N_Flux_Egest+                 &
     &                                    N_Flux_Pmortal)+              &
     &                         (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
#endif
              END DO
            END DO
!
!-----------------------------------------------------------------------
!  Zooplankton basal metabolism to NH4  (rate: ZooBM), zooplankton
!  mortality to small detritus (rate: ZooMR), zooplankton ingestion
!  related excretion (rate: ZooER).
!-----------------------------------------------------------------------
!
            cff1=dtdays*ZooBM(ng)
            fac2=dtdays*ZooMR(ng)
            fac3=dtdays*ZooER(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                &
     &               (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
                cff2=fac2*Bio(i,k,iZoop)
                cff3=fac1*ZooAE_N(ng)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)/                          &
     &                         (1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
                N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
                N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
                N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),      &
     &                                    0.0_r8)
                Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef OXYGEN
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                          &
     &                         rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
#endif
#ifdef PHOSPHORUS
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+                          &
     &                         ZooPN(ng)*(N_Flux_Zmetabo+               &
     &                                    N_Flux_Zexcret)
                Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                          &
     &                         ZooPN(ng)*N_Flux_Zmortal
#endif
              END DO
            END DO
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
            fac1=dtdays*CoagR(ng)
            DO k=1,N(ng)
              DO i=Istr,Iend
                cff1=fac1*(Bio(i,k,iSDeN)+Bio(i,k,iPhyt))
                cff2=1.0_r8/(1.0_r8+cff1)
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
                Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                N_Flux_CoagP=Bio(i,k,iPhyt)*cff1
                N_Flux_CoagD=Bio(i,k,iSDeN)*cff1
                Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                          &
     &                         N_Flux_CoagP+N_Flux_CoagD
#ifdef PHOSPHORUS
                Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
                Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                          &
     &                         PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
#endif
              END DO
            END DO
!
!-----------------------------------------------------------------------
!  Detritus recycling to NH4, remineralization.
!-----------------------------------------------------------------------
!
#ifdef OXYGEN
            DO k=1,N(ng)
              DO i=Istr,Iend
                fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
                fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 depe
                cff1=dtdays*SDeRRN(ng)*fac2
                cff2=1.0_r8/(1.0_r8+cff1)
                cff3=dtdays*LDeRRN(ng)*fac2
                cff4=1.0_r8/(1.0_r8+cff3)
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
                N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
# ifdef PHOSPHORUS
                cff1=dtdays*SDeRRP(ng)*fac2
                cff2=1.0_r8/(1.0_r8+cff1)
                cff3=dtdays*LDeRRP(ng)*fac2
                cff4=1.0_r8/(1.0_r8+cff3)
                Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
                Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
                P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
# endif
              END DO
            END DO
#else
            cff1=dtdays*SDeRRN(ng)
            cff2=1.0_r8/(1.0_r8+cff1)
            cff3=dtdays*LDeRRN(ng)
            cff4=1.0_r8/(1.0_r8+cff3)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
                Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
                N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
              END DO
            END DO
# ifdef PHOSPHORUS
            cff1=dtdays*SDeRRP(ng)
            cff2=1.0_r8/(1.0_r8+cff1)
            cff3=dtdays*LDeRRP(ng)
            cff4=1.0_r8/(1.0_r8+cff3)
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
                Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
                P_Flux_Remine=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux_Remine
              END DO
            END DO
# endif
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8) ! O2 off max
                fac2=MAX(fac1/(K_DO(ng)+fac1),0.0_r8) ! MM for O2 depe
                cff1=dtdays*H2SOR(ng)*fac2
                cff2=1.0_r8/(1.0_r8+cff1)
                Bio(i,k,iH2S_)=Bio(i,k,iH2S_)*cff2
                S_Flux=Bio(i,k,iH2S_)*cff1
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-S_Flux*rOxH2S
              END DO
            END DO
#endif
#ifdef OXYGEN
!
!-----------------------------------------------------------------------
!  Surface O2 gas exchange.
!-----------------------------------------------------------------------
!
!  Compute surface O2 gas exchange.
!
            cff1=rho0*550.0_r8
            cff2=dtdays*0.31_r8*24.0_r8/100.0_r8
            k=N(ng)
            DO i=Istr,Iend
!
!  Compute O2 transfer velocity : u10squared (u10 in m/s)
!
# ifdef BULK_FLUXES
              u10squ=Uwind(i,j)*Uwind(i,j)+Vwind(i,j)*Vwind(i,j)
# else
              u10squ=cff1*SQRT((0.5_r8*(sustr(i,j)+sustr(i+1,j)))**2+   &
     &                         (0.5_r8*(svstr(i,j)+svstr(i,j+1)))**2)
# endif
# ifdef OCMIP_OXYGEN_SC
!
!  Alternative formulation for Schmidt number (Sc will be slightly
!  smaller up to about 35 C): Compute the Schmidt number of oxygen
!  in seawater using the formulation proposed by Keeling et al.
!  (1998, Global Biogeochem. Cycles, 12, 141-163).  Input temperature
!  in Celsius.
!
              SchmidtN_Ox=1638.0_r8-                                    &
     &                    Bio(i,k,itemp)*(81.83_r8-                     &
     &                                    Bio(i,k,itemp)*               &
     &                                   (1.483_r8-                     &
     &                                    Bio(i,k,itemp)*0.008004_r8))
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
              SchmidtN_Ox=1953.4_r8-                                    &
     &                    Bio(i,k,itemp)*(128.0_r8-                     &
     &                                    Bio(i,k,itemp)*               &
     &                                   (3.9918_r8-                    &
     &                                    Bio(i,k,itemp)*0.050091_r8))
# endif

              cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
              TS=LOG((298.15_r8-Bio(i,k,itemp))/                        &
     &               (273.15_r8+Bio(i,k,itemp)))
              AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+        &
     &           Bio(i,k,isalt)*(OB0+TS*(OB1+TS*(OB2+TS*OB3)))+         &
     &           OC0*Bio(i,k,isalt)*Bio(i,k,isalt)
!
!  Convert from ml/l to mmol/m3.
!
              O2satu=l2mol*EXP(AA)
!
!  Add in O2 gas exchange.
!
              O2_Flux=cff3*(O2satu-Bio(i,k,iOxyg))
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                            &
     &                       O2_Flux*Hz_inv(i,k)
            END DO
#endif
            IF (Iteradj.ne.Iter) THEN
!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
              DO isink=1,Nsink
                ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    qc(i,k)=Bio(i,k,ibio)
                  END DO
                END DO
!
                DO k=N(ng)-1,1,-1
                  DO i=Istr,Iend
                    FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
                  END DO
                END DO
                DO k=2,N(ng)-1
                  DO i=Istr,Iend
                    dltR=Hz(i,j,k)*FC(i,k)
                    dltL=Hz(i,j,k)*FC(i,k-1)
                    cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                    cffR=cff*FC(i,k)
                    cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                    IF ((dltR*dltL).le.0.0_r8) THEN
                      dltR=0.0_r8
                      dltL=0.0_r8
                    ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                      dltR=cffL
                    ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                      dltL=cffR
                    END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                    cff=(dltR-dltL)*Hz_inv3(i,k)
                    dltR=dltR-cff*Hz(i,j,k+1)
                    dltL=dltL+cff*Hz(i,j,k-1)
                    bR(i,k)=qc(i,k)+dltR
                    bL(i,k)=qc(i,k)-dltL
                    WR(i,k)=(2.0_r8*dltR-dltL)**2
                    WL(i,k)=(dltR-2.0_r8*dltL)**2
                  END DO
                END DO
                cff=1.0E-14_r8
                DO k=2,N(ng)-2
                  DO i=Istr,Iend
                    dltL=MAX(cff,WL(i,k  ))
                    dltR=MAX(cff,WR(i,k+1))
                    bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                    bL(i,k+1)=bR(i,k)
                  END DO
                END DO
                DO i=Istr,Iend
                  FC(i,N(ng))=0.0_r8        ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
                  bL(i,N(ng))=bR(i,N(ng)-1)
                  bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
                  bR(i,N(ng))=qc(i,N(ng))   ! default strictly monotonic
                  bL(i,N(ng))=qc(i,N(ng))   ! conditions
                  bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
                  bR(i,1)=bL(i,2)
                  bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
                  bR(i,1)=bL(i,2)
                  bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
                  bL(i,2)=qc(i,1)           ! bottom grid boxes are
                  bR(i,1)=qc(i,1)           ! re-assumed to be
                  bL(i,1)=qc(i,1)           ! piecewise constant.
#endif
                END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    dltR=bR(i,k)-qc(i,k)
                    dltL=qc(i,k)-bL(i,k)
                    cffR=2.0_r8*dltR
                    cffL=2.0_r8*dltL
                    IF ((dltR*dltL).lt.0.0_r8) THEN
                      dltR=0.0_r8
                      dltL=0.0_r8
                    ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                      dltR=cffL
                    ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                      dltL=cffR
                    END IF
                    bR(i,k)=qc(i,k)+dltR
                    bL(i,k)=qc(i,k)-dltL
                  END DO
                END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
                cff=dtdays*ABS(Wbio(isink))
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    FC(i,k-1)=0.0_r8
                    WL(i,k)=z_w(i,j,k-1)+cff
                    WR(i,k)=Hz(i,j,k)*qc(i,k)
                    ksource(i,k)=k
                  END DO
                END DO
                DO k=1,N(ng)
                  DO ks=k,N(ng)-1
                    DO i=Istr,Iend
                      IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                        ksource(i,k)=ks+1
                        FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                      END IF
                    END DO
                  END DO
                END DO
!
!  Finalize computation of flux: add fractional part.
!
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    ks=ksource(i,k)
                    cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                    FC(i,k-1)=FC(i,k-1)+                                &
     &                        Hz(i,j,ks)*cu*                            &
     &                        (bL(i,ks)+                                &
     &                         cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-          &
     &                             (1.5_r8-cu)*                         &
     &                             (bR(i,ks)+bL(i,ks)-                  &
     &                              2.0_r8*qc(i,ks))))
                  END DO
                END DO
                DO k=1,N(ng)
                  DO i=Istr,Iend
                    Bio(i,k,ibio)=qc(i,k)+                              &
     &                            (FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                  END DO
                END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
                cff2=4.0_r8/16.0_r8
# ifdef OXYGEN
                cff3=115.0_r8/16.0_r8
                cff4=106.0_r8/16.0_r8
# endif
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
# ifdef DENITRIFICATION
                    Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
#  ifdef OXYGEN
                    Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3
#  endif
# else
                    Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
#  ifdef OXYGEN
                    Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
#  endif
# endif
                  END DO
                END IF
# ifdef PHOSPHORUS
                IF ((ibio.eq.iLDeP).or.                                 &
     &              (ibio.eq.iSDeP)) THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
                    Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
                  END DO
                END IF
                IF (ibio.eq.iPhyt)THEN
                  DO i=Istr,Iend
                    cff1=FC(i,0)*Hz_inv(i,1)
                    Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*PhyPN(ng)
                  END DO
                END IF
# endif
# if defined H2S && defined OXYGEN
                IF ((ibio.eq.iPhyt).or.                                 &
     &              (ibio.eq.iSDeN).or.                                 &
     &              (ibio.eq.iLDeN)) THEN
                  DO i=Istr,Iend
                    Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-                      &
     &                             0.5_r8*MIN(Bio(i,1,iOxyg),0.0_r8)
                    Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
                  END DO
                END IF
# endif
#endif
              END DO
            END IF
          END DO
!
!  End of compute basic state arrays III.
!
!-----------------------------------------------------------------------
!  Tangent linear vertical sinking terms.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of selected biological constituents
!  "Bio(:,:,isink)" in terms of a set of parabolic segments within each
!  grid box. Then, compute semi-Lagrangian flux due to sinking.
!
          SINK_LOOP: DO isink=1,Nsink
            ibio=idsink(isink)
!
!  Copy concentration of biological particulates into scratch array
!  "qc" (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for biogeochemical
!  constituent concentration.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                qc(i,k)=Bio(i,k,ibio)
                tl_qc(i,k)=tl_Bio(i,k,ibio)
              END DO
            END DO
!
            DO k=N(ng)-1,1,-1
              DO i=Istr,Iend
                FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
                tl_FC(i,k)=(tl_qc(i,k+1)-tl_qc(i,k))*Hz_inv2(i,k)+      &
     &                     (qc(i,k+1)-qc(i,k))*tl_Hz_inv2(i,k)
              END DO
            END DO
            DO k=2,N(ng)-1
              DO i=Istr,Iend
                dltR=Hz(i,j,k)*FC(i,k)
                tl_dltR=tl_Hz(i,j,k)*FC(i,k)+Hz(i,j,k)*tl_FC(i,k)
                dltL=Hz(i,j,k)*FC(i,k-1)
                tl_dltL=tl_Hz(i,j,k)*FC(i,k-1)+Hz(i,j,k)*tl_FC(i,k-1)
                cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
                tl_cff=tl_Hz(i,j,k-1)+2.0_r8*tl_Hz(i,j,k)+tl_Hz(i,j,k+1)
                cffR=cff*FC(i,k)
                tl_cffR=tl_cff*FC(i,k)+cff*tl_FC(i,k)
                cffL=cff*FC(i,k-1)
                tl_cffL=tl_cff*FC(i,k-1)+cff*tl_FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
                IF ((dltR*dltL).le.0.0_r8) THEN
                  dltR=0.0_r8
                  tl_dltR=0.0_r8
                  dltL=0.0_r8
                  tl_dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                  tl_dltR=tl_cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                  tl_dltL=tl_cffR
                END IF
!
!  Compute right and left side values (bR,bL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because bL(k+1)-bR(k) may still have different sign than
!        qc(i,k+1)-qc(i,k).  This possibility is excluded,
!        after bL and bR are reconciled using WENO procedure.
!
                cff=(dltR-dltL)*Hz_inv3(i,k)
                tl_cff=(tl_dltR-tl_dltL)*Hz_inv3(i,k)+                  &
     &                 (dltR-dltL)*tl_Hz_inv3(i,k)
                dltR=dltR-cff*Hz(i,j,k+1)
                tl_dltR=tl_dltR-tl_cff*Hz(i,j,k+1)-cff*tl_Hz(i,j,k+1)
                dltL=dltL+cff*Hz(i,j,k-1)
                tl_dltL=tl_dltL+tl_cff*Hz(i,j,k-1)+cff*tl_Hz(i,j,k-1)
                bR(i,k)=qc(i,k)+dltR
                tl_bR(i,k)=tl_qc(i,k)+tl_dltR
                bL(i,k)=qc(i,k)-dltL
                tl_bL(i,k)=tl_qc(i,k)-tl_dltL
                WR(i,k)=(2.0_r8*dltR-dltL)**2
                tl_WR(i,k)=2.0_r8*(2.0_r8*dltR-dltL)*                   &
     &                            (2.0_r8*tl_dltR-tl_dltL)
                WL(i,k)=(dltR-2.0_r8*dltL)**2
                tl_WL(i,k)=2.0_r8*(dltR-2.0_r8*dltL)*                   &
     &                            (tl_dltR-2.0_r8*tl_dltL)
              END DO
            END DO
            cff=1.0E-14_r8
            DO k=2,N(ng)-2
              DO i=Istr,Iend
                dltL=MAX(cff,WL(i,k  ))
                tl_dltL=(0.5_r8-SIGN(0.5_r8,cff-WL(i,k  )))*            &
     &                  tl_WL(i,k  )
                dltR=MAX(cff,WR(i,k+1))
                tl_dltR=(0.5_r8-SIGN(0.5_r8,cff-WR(i,k+1)))*            &
     &                  tl_WR(i,k+1)
                bR1(i,k)=bR(i,k)
                bL1(i,k+1)=bL(i,k+1)
                bR(i,k)=(dltR*bR(i,k)+dltL*bL(i,k+1))/(dltR+dltL)
                tl_bR(i,k)=(tl_dltR*bR1(i,k  )+dltR*tl_bR(i,k  )+       &
     &                      tl_dltL*bL1(i,k+1)+dltL*tl_bL(i,k+1))/      &
     &                      (dltR+dltL)-                                &
     &                      (tl_dltR+tl_dltL)*bR(i,k)/(dltR+dltL)
                bL(i,k+1)=bR(i,k)
                tl_bL(i,k+1)=tl_bR(i,k)
              END DO
            END DO
            DO i=Istr,Iend
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
              tl_FC(i,N(ng))=0.0_r8         ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              tl_bL(i,N(ng))=tl_bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
              tl_bR(i,N(ng))=2.0_r8*tl_qc(i,N(ng))-tl_bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              tl_bL(i,N(ng))=tl_bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
              tl_bR(i,N(ng))=1.5_r8*tl_qc(i,N(ng))-0.5_r8*tl_bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
              tl_bR(i,N(ng))=tl_qc(i,N(ng)) ! default strictly monotonic
              tl_bL(i,N(ng))=tl_qc(i,N(ng)) ! conditions
              tl_bR(i,N(ng)-1)=tl_qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              tl_bR(i,1)=tl_bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
              tl_bL(i,1)=2.0_r8*tl_qc(i,1)-tl_bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              tl_bR(i,1)=tl_bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
              tl_bL(i,1)=1.5_r8*tl_qc(i,1)-0.5_r8*tl_bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
              tl_bL(i,2)=tl_qc(i,1)         ! bottom grid boxes are
              tl_bR(i,1)=tl_qc(i,1)         ! re-assumed to be
              tl_bL(i,1)=tl_qc(i,1)         ! piecewise constant.
#endif
            END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                dltR=bR(i,k)-qc(i,k)
                tl_dltR=tl_bR(i,k)-tl_qc(i,k)
                dltL=qc(i,k)-bL(i,k)
                tl_dltL=tl_qc(i,k)-tl_bL(i,k)
                cffR=2.0_r8*dltR
                tl_cffR=2.0_r8*tl_dltR
                cffL=2.0_r8*dltL
                tl_cffL=2.0_r8*tl_dltL
                IF ((dltR*dltL).lt.0.0_r8) THEN
                  dltR=0.0_r8
                  tl_dltR=0.0_r8
                  dltL=0.0_r8
                  tl_dltL=0.0_r8
                ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                  dltR=cffL
                  tl_dltR=tl_cffL
                ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                  dltL=cffR
                  tl_dltL=tl_cffR
                END IF
                bR(i,k)=qc(i,k)+dltR
                tl_bR(i,k)=tl_qc(i,k)+tl_dltR
                bL(i,k)=qc(i,k)-dltL
                tl_bL(i,k)=tl_qc(i,k)-tl_dltL
              END DO
            END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
            cff=dtdays*ABS(Wbio(isink))
            tl_cff=dtdays*SIGN(1.0_r8,Wbio(isink))*tl_Wbio(isink)
            DO k=1,N(ng)
              DO i=Istr,Iend
                FC(i,k-1)=0.0_r8
                tl_FC(i,k-1)=0.0_r8
                WL(i,k)=z_w(i,j,k-1)+cff
                tl_WL(i,k)=tl_z_w(i,j,k-1)+tl_cff
                WR(i,k)=Hz(i,j,k)*qc(i,k)
                tl_WR(i,k)=tl_Hz(i,j,k)*qc(i,k)+Hz(i,j,k)*tl_qc(i,k)
                ksource(i,k)=k
              END DO
            END DO
            DO k=1,N(ng)
              DO ks=k,N(ng)-1
                DO i=Istr,Iend
                  IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                    ksource(i,k)=ks+1
                    FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                    tl_FC(i,k-1)=tl_FC(i,k-1)+tl_WR(i,ks)
                  END IF
                END DO
              END DO
            END DO
!
!  Finalize computation of flux: add fractional part.
!
            DO k=1,N(ng)
              DO i=Istr,Iend
                ks=ksource(i,k)
                cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
                tl_cu=(0.5_r8+SIGN(0.5_r8,                              &
     &                             (1.0_r8-(WL(i,k)-z_w(i,j,ks-1))*     &
     &                             Hz_inv(i,ks))))*                     &
     &                ((tl_WL(i,k)-tl_z_w(i,j,ks-1))*Hz_inv(i,ks)+      &
     &                 (WL(i,k)-z_w(i,j,ks-1))*tl_Hz_inv(i,ks))
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
                tl_FC(i,k-1)=tl_FC(i,k-1)+                              &
     &                       (tl_Hz(i,j,ks)*cu+Hz(i,j,ks)*tl_cu)*       &
     &                       (bL(i,ks)+                                 &
     &                        cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-           &
     &                            (1.5_r8-cu)*                          &
     &                            (bR(i,ks)+bL(i,ks)-                   &
     &                             2.0_r8*qc(i,ks))))+                  &
     &                       Hz(i,j,ks)*cu*                             &
     &                       (tl_bL(i,ks)+                              &
     &                        tl_cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-        &
     &                               (1.5_r8-cu)*                       &
     &                               (bR(i,ks)+bL(i,ks)-                &
     &                               2.0_r8*qc(i,ks)))+                 &
     &                        cu*(0.5_r8*(tl_bR(i,ks)-tl_bL(i,ks))+     &
     &                            tl_cu*                                &
     &                            (bR(i,ks)+bL(i,ks)-2.0_r8*qc(i,ks))-  &
     &                            (1.5_r8-cu)*                          &
     &                            (tl_bR(i,ks)+tl_bL(i,ks)-             &
     &                             2.0_r8*tl_qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
                tl_Bio(i,k,ibio)=tl_qc(i,k)+                            &
     &                           (tl_FC(i,k)-tl_FC(i,k-1))*Hz_inv(i,k)+ &
     &                           (FC(i,k)-FC(i,k-1))*tl_Hz_inv(i,k)
              END DO
            END DO

#ifdef BIO_SEDIMENT
!
!  Particulate flux reaching the seafloor is remineralized and returned
!  to the dissolved nitrate pool. Without this conversion, particulate
!  material falls out of the system. This is a temporary fix to restore
!  total nitrogen conservation. It will be replaced later by a
!  parameterization that includes the time delay of remineralization
!  and dissolved oxygen.
!
            cff2=4.0_r8/16.0_r8
# ifdef OXYGEN
            cff3=115.0_r8/16.0_r8
            cff4=106.0_r8/16.0_r8
# endif
            IF ((ibio.eq.iPhyt).or.                                     &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN)) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
# ifdef DENITRIFICATION
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
#  ifdef OXYGEN
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3
#  endif
# else
                Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
#  ifdef OXYGEN
                Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
#  endif
# endif
              END DO
            END IF
# ifdef PHOSPHORUS
            IF ((ibio.eq.iLDeP).or.                                     &
     &          (ibio.eq.iSDeP)) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
              END DO
            END IF
            IF (ibio.eq.iPhyt)THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*PhyPN(ng)
              END DO
            END IF
# endif
# if defined H2S && defined OXYGEN
            IF ((ibio.eq.iPhyt).or.                                     &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN)) THEN
              DO i=Istr,Iend
                Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-                          &
     &                         0.5_r8*MIN(Bio(i,1,iOxyg),0.0_r8)
                Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
              END DO
            END IF
# endif
#endif
          END DO SINK_LOOP
#ifdef BIO_SEDIMENT_PARAMETER
!
!  Elution and oxygen consumption parameters (okada)
!
!    SOD 0.2-4.0 g/m2/day  Sediment Oxygen Demand (WASP6)
!    NH4 15.3    mg/m2/day Average Nishimoto (2012)
!    PO4 2.0     mg/m2/day Average Nishimoto (2012)
!
          cff1=2100.0_r8/32.0_r8  !H2S elution flux from sediment
          cff2=15.3_r8/14.0_r8    !NH4 elution flux from sediment
          cff3=2.0_r8/31.0_r8     !PO4 elution flux from sediment
!
!-----------------------------------------------------------------------
!  Elution and oxygen consumption from/by sediment. (Okada, 2014/02/13)
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
            fac1=dtdays*1.05_r8**(Bio(i,1,itemp)-20.0_r8)
            cff1=fac1*Hz_inv(i,1)
            Bio(i,1,iH2S_)=Bio(i,1,iH2S_)+cff1*cff8
            Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff9  !*0.7_r8
            Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*cff10 !*0.7_r8
          END DO
#endif
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
            END DO
          END DO
        END DO

      END DO J_LOOP

      RETURN
      END SUBROUTINE tl_biology_tile
