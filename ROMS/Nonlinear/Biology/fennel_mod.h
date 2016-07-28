!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2016 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Fennel et al. (2006) model:                          !
!                                                                      !
!   AttSW    Light attenuation due to sea water [1/m].                 !
!   AttChl   Light attenuation by Chlorophyll [1/(mg_Chl m2)].         !
!   BioIter  Maximum number of iterations to achieve convergence       !
!              of the nonlinear solution.                              !
!   Chl2C_m  Maximum chlorophyll to carbon ratio [mg_Chl/mg_C].        !
!   ChlMin   Chlorophill minimum threshold value [mg_Chl/m3].          !
!   CoagR    Coagulation rate: agregation rate of SDeN + Phyt ==> LDeN !
!              [1/day].                                                !
!   D_p5NH4  Half-saturation radiation for nitrification inhibition    !
!              [Watts/m2].                                             !
!   I_thNH4  Radiation threshold for nitrification inhibition          !
!              [Watts/m2].                                             !
!   K_NH4    Inverse half-saturation for Phytoplankton NH4 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_NO3    Inverse half-saturation for Phytoplankton NO3 uptake      !
!              [m3/(mmol_N)].                                          !
!   K_Phy    Zooplankton half-saturation, squared constant for         !
!              ingestion [mmol_N/m3]^2.                                !
!   LDeRR    Large Detrital re-mineralization rate [1/day].            !
!   NitriR   Nitrification rate: oxidation of NH4 to NO3 [1/day].      !
!   PARfrac  Fraction of shortwave radiation that is available for     !
!              photosyntesis [nondimensional].                         !
!   PhyCN    Phytoplankton Carbon:Nitrogen ratio [mol_C/mol_N].        !
!   PhyIP    Phytoplankton NH4 inhibition parameter [1/(mmol_N)].      !
!   PhyIS    Phytoplankton, initial slope of the P-I curve             !
!              [1/(W m-2 day)].                                        !
!   ZooMin   Phytoplankton minimum threshold value [mmol_N/m3].        !
!   PhyMR    Phytoplankton mortality rate [1/day] to small detritus.   !
!   SDeAR    Small detritus aggregation rate into Large detritus       !
!              [1/day].                                                !
!   SDeBR    Small Detrital breakdown to NH4 rate [1/day].             !
!   SDeRR    Large Detrital re-mineralization rate [1/day].            !
!   Vp0      Eppley temperature-limited and light-limited growth       !
!              tuning parameter [nondimensional].                      !
!   wLDet    Vertical sinking velocities for Large Detritus            !
!              fraction [m/day].                                       !
!   wPhy     Vertical sinking velocity for Phytoplankton               !
!              fraction [m/day].                                       !
!   wSDet    Vertical sinking velocities for Small Detritus            !
!              fraction [m/day].                                       !
!   ZooAE_N  Zooplankton nitrogen assimilation efficiency fraction     !
!              [nondimensional].                                       !
!   ZooBM    Zooplankton basal metabolism [1/day].                     !
!   ZooCN    Zooplankton Carbon:Nitrogen ratio [mol_C/mol_N].          !
!   ZooER    Zooplankton specific excretion rate [1/day].              !
!   ZooGR    Zooplankton maximum growth rate [1/day].                  !
!   ZooMin   Zooplankton minimum threshold value [mmol_N/m3].          !
!   ZooMR    Zooplankton mortality to Detritus [1/day].                !
!   pCO2air  CO2 partial pressure in the air [ppmv].                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers
      integer :: iNO3_                  ! Nitrate concentration
      integer :: iNH4_                  ! Ammonium concentration
      integer :: iChlo                  ! Chlorophyll concentration
      integer :: iPhyt                  ! Phytoplankton concentration
      integer :: iZoop                  ! Zooplankton concentration
      integer :: iLDeN                  ! Large detritus N-concentration
      integer :: iSDeN                  ! Small detritus N-concentration
#ifdef CARBON
      integer :: iLDeC                  ! Large detritus C-concentration
      integer :: iSDeC                  ! Small detritus C-concentration
      integer :: iTIC_                  ! Total inorganic carbon
      integer :: iTAlk                  ! Total alkalinity
#endif
#ifdef OXYGEN
      integer :: iOxyg                  ! Dissolved oxygen concentration
#endif

#ifdef PHOSPHORUS
      integer :: iPO4_                  ! 
      integer :: iLDeP                  ! 
      integer :: iSDeP                  ! 
#endif
#ifdef H2S
      integer :: iH2S_                  ! 
#endif

#ifdef DIAGENESIS
!
!  Set bgc tracer identification indices.
!
      integer :: dia_count = 1
!
      integer, parameter :: NBGCPW = 11        ! Bgc variables in pore water 
      integer, parameter :: NBGCSM = 11        ! Bgc variables in sediment
      integer, parameter :: NBGCF  = 50        ! Bgc circulation variables in sediment
!
      integer :: idBpw(NBGCPW)                 ! Bgc variables in pore water
      integer :: idBsm(NBGCSM)                 ! Bgc variables in sediment
      integer :: idFpw(NBGCPW)                  ! Bgc variables in sediment
      integer :: idFsm(NBGCSM)                  ! Bgc variables in sediment
!
      integer, parameter :: iwO2_ = 1          ! O2 in pore water 
      integer, parameter :: iwNH4 = 2          ! NH4 in pore water 
      integer, parameter :: iwNO3 = 3          ! NO3 in pore water 
      integer, parameter :: iwPO4 = 4          ! PO4 in pore water 
      integer, parameter :: iwSO4 = 5          ! SO4 in pore water 
      integer, parameter :: iwH2S = 6          ! H2S in pore water 
      integer, parameter :: iwMn_ = 7          ! Mn in pore water 
      integer, parameter :: iwFe_ = 8          ! Fe in pore water 
      integer, parameter :: iwCH4 = 9          ! CH4 in pore water
      integer, parameter :: iwDOMf = 10        ! DOMfast in pore water !Nchange
      integer, parameter :: iwDOMs = 11        ! DOMslow in pore water !Nchange
!
      integer, parameter :: iPOMf = 1          ! Perticle Organic Mater fast
      integer, parameter :: iPOMs = 2          ! Perticle Organic Mater slow
      integer, parameter :: iPOMn = 3          ! Perticle Organic Mater non
      integer, parameter :: iFeOA = 4          ! FeOOHA in sediment 
      integer, parameter :: iFeOB = 5          ! FeOOHB in sediment 
      integer, parameter :: iFeOP = 6          ! FeOOH-PO4 in sediment 
      integer, parameter :: iMnOA = 7          ! MnO2A in sediment 
      integer, parameter :: iMnOB = 8          ! MnO2B in sediment 
      integer, parameter :: iS0__ = 9          ! S0 in sediment 
      integer, parameter :: iFeS_ = 10         ! FeS in sediment 
      integer, parameter :: iFes2 = 11         ! FeS2 in sediment
#endif

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Biological 2D diagnostic variable IDs.
!
      integer, allocatable :: iDbio2(:)       ! 2D biological terms

      integer  :: iNH4f = 1                   ! NH4 flux
      integer  :: iPONf = 2                   ! PON flux

      integer  :: iDNIT                       ! denitrification flux
      integer  :: iCOfx                       ! air-sea CO2 flux
      integer  :: ipCO2                       ! partial pressure of CO2
      integer  :: iO2fx                       ! air-sea O2 flux

      integer  :: iPO4f                       ! PO4 flux
      integer  :: iPOPf                       ! POP flux
      integer  :: iSODf                       ! sediment oxygen demand flux
!
!  Biological 3D diagnostic variable IDs.
!
      integer, allocatable :: iDbio3(:)       ! 3D biological terms

      integer  :: iPPro = 1                   ! primary productivity
      integer  :: iNO3u = 2                   ! NO3 uptake

      integer  :: iPmax = 3                   ! max growth rate
      integer  :: iLDIN = 4                   ! DIN limitation
      integer  :: iNitr = 5
      integer  :: iDeni = 6
      integer  :: iAssi = 7
      integer  :: iEges = 8
      integer  :: iPmor = 9
      integer  :: iZmor = 10
      integer  :: iZexc = 11
      integer  :: iZmet = 12
      integer  :: iCoaP = 13
      integer  :: iCoaD = 14
      integer  :: iSReN = 15
      integer  :: iLReN = 16
      integer  :: iLPO4                       ! PO4 limitation
      integer  :: iSReP
      integer  :: iLReP
      integer  :: iCOD_
      integer  :: iH2Sf
#endif
!
!  Biological parameters.
!
      integer, allocatable :: BioIter(:)

      real(r8), allocatable :: K_PO4(:)
      real(r8), allocatable :: PhyPN(:)
      real(r8), allocatable :: ZooPN(:)
      real(r8), allocatable :: LDeRRP(:)
      real(r8), allocatable :: SDeRRP(:)
      real(r8), allocatable :: H2SOR(:)
      real(r8), allocatable :: K_DO(:) 
      real(r8), allocatable :: K_Nitri(:)
      real(r8), allocatable :: NitriR_t(:)
      real(r8), allocatable :: g_max(:)
      real(r8), allocatable :: t_opt(:)
      real(r8), allocatable :: I_opt(:)
      real(r8), allocatable :: beta1(:)
      real(r8), allocatable :: beta2(:)
      real(r8), allocatable :: DenitR(:)
      real(r8), allocatable :: K_Denit(:)
      real(r8), allocatable :: DenitR_t(:)
      real(r8), allocatable :: PhyMR_t(:)
      real(r8), allocatable :: ZooGR_t(:)
      real(r8), allocatable :: RR_t(:)
      real(r8), allocatable :: R_SODf(:)
      real(r8), allocatable :: R_NH4f(:)
      real(r8), allocatable :: R_PO4f(:)
      real(r8), allocatable :: R_NH4f_max(:)
      real(r8), allocatable :: R_PO4f_max(:)
      real(r8), allocatable :: K_DO_npflux(:)
      real(r8), allocatable :: t_SODf(:)

      real(r8), allocatable :: AttSW(:)              ! 1/m
      real(r8), allocatable :: AttChl(:)             ! 1/(mg_Chl m2)
      real(r8), allocatable :: Chl2C_m(:)            ! mg_Chl/mg_C
      real(r8), allocatable :: ChlMin(:)             ! mg_Chl/m3
      real(r8), allocatable :: CoagR(:)              ! 1/day
      real(r8), allocatable :: D_p5NH4(:)            ! Watts/m2
      real(r8), allocatable :: I_thNH4(:)            ! Watts/m2
      real(r8), allocatable :: K_NH4(:)              ! m3/mmol_N
      real(r8), allocatable :: K_NO3(:)              ! m3/mmol_N
      real(r8), allocatable :: K_Phy(:)              ! (mmol_N/m3)^2
      real(r8), allocatable :: LDeRRN(:)             ! 1/day
      real(r8), allocatable :: LDeRRC(:)             ! 1/day
      real(r8), allocatable :: NitriR(:)             ! 1/day
      real(r8), allocatable :: PARfrac(:)            ! nondimensional
      real(r8), allocatable :: PhyCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: PhyIP(:)              ! 1/mmol_N
      real(r8), allocatable :: PhyIS(:)              ! 1/(Watts m-2 day)
      real(r8), allocatable :: PhyMin(:)             ! mmol_N/m3
      real(r8), allocatable :: PhyMR(:)              ! 1/day
      real(r8), allocatable :: SDeAR(:)              ! 1/day
      real(r8), allocatable :: SDeBR(:)              ! 1/day
      real(r8), allocatable :: SDeRRN(:)             ! 1/day
      real(r8), allocatable :: SDeRRC(:)             ! 1/day
      real(r8), allocatable :: Vp0(:)                ! nondimensional
      real(r8), allocatable :: wLDet(:)              ! m/day
      real(r8), allocatable :: wPhy(:)               ! m/day
      real(r8), allocatable :: wSDet(:)              ! m/day
      real(r8), allocatable :: ZooAE_N(:)            ! nondimensional
      real(r8), allocatable :: ZooBM(:)              ! 1/day
      real(r8), allocatable :: ZooCN(:)              ! mol_C/mol_N
      real(r8), allocatable :: ZooER(:)              ! 1/day
      real(r8), allocatable :: ZooGR(:)              ! 1/day
      real(r8), allocatable :: ZooMin(:)             ! mmol_N/m3
      real(r8), allocatable :: ZooMR(:)              ! 1/day
      real(r8), allocatable :: pCO2air(:)            ! ppmv
#ifdef TANGENT
      real(r8) :: tl_K_PO4
      real(r8) :: tl_PhyPN
      real(r8) :: tl_ZooPN
      real(r8) :: tl_LDeRRP
      real(r8) :: tl_SDeRRP
      real(r8) :: tl_H2SOR
      real(r8) :: tl_K_DO 
      real(r8) :: tl_K_Nitri
      real(r8) :: tl_NitriR_t
      real(r8) :: tl_g_max
      real(r8) :: tl_t_opt
      real(r8) :: tl_I_opt
      real(r8) :: tl_beta1
      real(r8) :: tl_beta2
      real(r8) :: tl_DenitR
      real(r8) :: tl_K_Denit
      real(r8) :: tl_DenitR_t
      real(r8) :: tl_PhyMR_t
      real(r8) :: tl_ZooGR_t
      real(r8) :: tl_RR_t
      real(r8) :: tl_R_SODf
      real(r8) :: tl_R_NH4f
      real(r8) :: tl_R_PO4f
      real(r8) :: tl_R_NH4f_max
      real(r8) :: tl_R_PO4f_max
      real(r8) :: tl_K_DO_npflux
      real(r8) :: tl_t_SODf

      real(r8) :: tl_AttSW              ! 1/m
      real(r8) :: tl_AttChl             ! 1/(mg_Chl m2)
      real(r8) :: tl_Chl2C_m            ! mg_Chl/mg_C
      real(r8) :: tl_ChlMin             ! mg_Chl/m3
      real(r8) :: tl_CoagR              ! 1/day
      real(r8) :: tl_D_p5NH4            ! Watts/m2
      real(r8) :: tl_I_thNH4            ! Watts/m2
      real(r8) :: tl_K_NH4              ! m3/mmol_N
      real(r8) :: tl_K_NO3              ! m3/mmol_N
      real(r8) :: tl_K_Phy              ! (mmol_N/m3)^2
      real(r8) :: tl_LDeRRN             ! 1/day
      real(r8) :: tl_LDeRRC             ! 1/day
      real(r8) :: tl_NitriR             ! 1/day
      real(r8) :: tl_PARfrac            ! nondimensional
      real(r8) :: tl_PhyCN              ! mol_C/mol_N
      real(r8) :: tl_PhyIP              ! 1/mmol_N
      real(r8) :: tl_PhyIS              ! 1/(Watts m-2 day)
      real(r8) :: tl_PhyMin             ! mmol_N/m3
      real(r8) :: tl_PhyMR              ! 1/day
      real(r8) :: tl_SDeAR              ! 1/day
      real(r8) :: tl_SDeBR              ! 1/day
      real(r8) :: tl_SDeRRN             ! 1/day
      real(r8) :: tl_SDeRRC             ! 1/day
      real(r8) :: tl_Vp0                ! nondimensional
      real(r8) :: tl_wLDet              ! m/day
      real(r8) :: tl_wPhy               ! m/day
      real(r8) :: tl_wSDet              ! m/day
      real(r8) :: tl_ZooAE_N            ! nondimensional
      real(r8) :: tl_ZooBM              ! 1/day
      real(r8) :: tl_ZooCN              ! mol_C/mol_N
      real(r8) :: tl_ZooER              ! 1/day
      real(r8) :: tl_ZooGR              ! 1/day
      real(r8) :: tl_ZooMin             ! mmol_N/m3
      real(r8) :: tl_ZooMR              ! 1/day
      real(r8) :: tl_pCO2air            ! ppmv
!      real(r8), allocatable :: tl_wSDet(:)
!      real(r8), allocatable :: tl_wPhy(:)
!      real(r8), allocatable :: tl_wLDet(:)
!      real(r8), allocatable :: tl_PARfrac(:)
#endif
#ifdef ADJOINT
      real(r8) :: ad_K_PO4
      real(r8) :: ad_PhyPN
      real(r8) :: ad_ZooPN
      real(r8) :: ad_LDeRRP
      real(r8) :: ad_SDeRRP
      real(r8) :: ad_H2SOR
      real(r8) :: ad_K_DO 
      real(r8) :: ad_K_Nitri
      real(r8) :: ad_NitriR_t
      real(r8) :: ad_g_max
      real(r8) :: ad_t_opt
      real(r8) :: ad_I_opt
      real(r8) :: ad_beta1
      real(r8) :: ad_beta2
      real(r8) :: ad_DenitR
      real(r8) :: ad_K_Denit
      real(r8) :: ad_DenitR_t
      real(r8) :: ad_PhyMR_t
      real(r8) :: ad_ZooGR_t
      real(r8) :: ad_RR_t
      real(r8) :: ad_R_SODf
      real(r8) :: ad_R_NH4f
      real(r8) :: ad_R_PO4f
      real(r8) :: ad_R_NH4f_max
      real(r8) :: ad_R_PO4f_max
      real(r8) :: ad_K_DO_npflux
      real(r8) :: ad_t_SODf

      real(r8) :: ad_AttSW              ! 1/m
      real(r8) :: ad_AttChl             ! 1/(mg_Chl m2)
      real(r8) :: ad_Chl2C_m            ! mg_Chl/mg_C
      real(r8) :: ad_ChlMin             ! mg_Chl/m3
      real(r8) :: ad_CoagR              ! 1/day
      real(r8) :: ad_D_p5NH4            ! Watts/m2
      real(r8) :: ad_I_thNH4            ! Watts/m2
      real(r8) :: ad_K_NH4              ! m3/mmol_N
      real(r8) :: ad_K_NO3              ! m3/mmol_N
      real(r8) :: ad_K_Phy              ! (mmol_N/m3)^2
      real(r8) :: ad_LDeRRN             ! 1/day
      real(r8) :: ad_LDeRRC             ! 1/day
      real(r8) :: ad_NitriR             ! 1/day
      real(r8) :: ad_PARfrac            ! nondimensional
      real(r8) :: ad_PhyCN              ! mol_C/mol_N
      real(r8) :: ad_PhyIP              ! 1/mmol_N
      real(r8) :: ad_PhyIS              ! 1/(Watts m-2 day)
      real(r8) :: ad_PhyMin             ! mmol_N/m3
      real(r8) :: ad_PhyMR              ! 1/day
      real(r8) :: ad_SDeAR              ! 1/day
      real(r8) :: ad_SDeBR              ! 1/day
      real(r8) :: ad_SDeRRN             ! 1/day
      real(r8) :: ad_SDeRRC             ! 1/day
      real(r8) :: ad_Vp0                ! nondimensional
      real(r8) :: ad_wLDet              ! m/day
      real(r8) :: ad_wPhy               ! m/day
      real(r8) :: ad_wSDet              ! m/day
      real(r8) :: ad_ZooAE_N            ! nondimensional
      real(r8) :: ad_ZooBM              ! 1/day
      real(r8) :: ad_ZooCN              ! mol_C/mol_N
      real(r8) :: ad_ZooER              ! 1/day
      real(r8) :: ad_ZooGR              ! 1/day
      real(r8) :: ad_ZooMin             ! mmol_N/m3
      real(r8) :: ad_ZooMR              ! 1/day
      real(r8) :: ad_pCO2air            ! ppmv
!      real(r8), allocatable :: ad_wSDet(:)
!      real(r8), allocatable :: ad_wPhy(:)
!      real(r8), allocatable :: ad_wLDet(:)
!      real(r8), allocatable :: ad_PARfrac(:)
#endif

#ifdef ADJUST_PARAM
      integer :: iAttSW   = 1
      integer :: iAttChl  = 2
      integer :: iVp0     = 3
      integer :: iI_thNH4 = 4
      integer :: iD_p5NH4 = 5
      integer :: iNitriR  = 6
      integer :: iK_NO3   = 7
      integer :: iK_NH4   = 8
      integer :: iK_Phy   = 9
      integer :: iChl2C_m = 10
      integer :: iPhyCN   = 11
      integer :: iPhyIP   = 12
      integer :: iPhyIS   = 13
      integer :: iPhyMR   = 14
      integer :: iZooAE_N = 15
      integer :: iZooBM   = 16
      integer :: iZooCN   = 17
      integer :: iZooER   = 18
      integer :: iZooGR   = 19
      integer :: iZooMR   = 20
      integer :: iLDeRRN  = 21
      integer :: iCoagR   = 22
      integer :: iSDeRRN  = 23
      integer :: iwPhy    = 24
      integer :: iwLDet   = 25
      integer :: iwSDet   = 26

      integer :: iK_Nitri = 27
      integer :: iK_Denit = 28
      integer :: iDenitR  = 29
      integer :: iK_PO4   = 30
      integer :: iPhyPN   = 31
      integer :: iZooPN   = 32
      integer :: iK_DO    = 33
      integer :: iLDeRRP  = 34
      integer :: iSDeRRP  = 35
      integer :: iR_SODf  = 36
      integer :: iR_NH4f  = 37
      integer :: iR_PO4f  = 38
      integer :: iR_NH4f_m = 39
      integer :: iR_PO4f_m = 40
      integer :: iK_DO_npf = 41
      integer :: it_SODf  = 42
#endif

      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Determine number of biological tracers.
!-----------------------------------------------------------------------
!
#ifdef CARBON
# ifdef OXYGEN
      NBT=12
# else
      NBT=11
# endif
#else
# ifdef OXYGEN
      NBT=8
# else
      NBT=7
# endif
#endif

#ifdef PHOSPHORUS
      NBT=NBT+3
#endif
#ifdef H2S
      NBT=NBT+1
#endif

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio2d=2   ! NH4f, PONf
      NDbio3d=16  ! PPro, NO3u, Pmax, LDIN, ...
# ifdef DENITRIFICATION
      NDbio2d=NDbio2d+1  ! DNIT
# endif
# ifdef CARBON
      NDbio2d=NDbio2d+2  ! COfx, pCO2
# endif
# ifdef OXYGEN
      NDbio2d=NDbio2d+2  ! O2fx ,SODf
      NDbio3d=NDbio3d+2  ! COD, H2Sf
# endif
# ifdef PHOSPHORUS
      NDbio2d=NDbio2d+2  ! PO4f, POPf
      NDbio3d=NDbio3d+3  ! LPO4, SReP, LReP
# endif
!
!  Initialize biology diagnostic indices.
!
      ic=2  ! BioDia2d
# ifdef DENITRIFICATION
      iDNIT=ic+1
      ic=ic+1
# endif
# ifdef CARBON
      iCOfx=ic+1
      ipCO2=ic+2
      ic=ic+2
# endif
# ifdef OXYGEN
      iO2fx=ic+1
      iSODf=ic+2
      ic=ic+2
# endif
# ifdef PHOSPHORUS
      iPO4f=ic+1
      iPOPf=ic+2
# endif

      ic=16  ! BioDia3d
# ifdef PHOSPHORUS
      iLPO4=ic+1
      iSReP=ic+2
      iLReP=ic+3
      ic=ic+3
# endif
# ifdef OXYGEN
      iCOD_=ic+1
      iH2Sf=ic+2
# endif
#endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(K_PO4)) THEN
        allocate ( K_PO4(Ngrids) )
      END IF
      IF (.not.allocated(PhyPN)) THEN
        allocate ( PhyPN(Ngrids) )
      END IF
      IF (.not.allocated(ZooPN)) THEN
        allocate ( ZooPN(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRP)) THEN
        allocate ( LDeRRP(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRP)) THEN
        allocate ( SDeRRP(Ngrids) )
      END IF

      IF (.not.allocated(g_max)) THEN
        allocate ( g_max(Ngrids) )
      END IF
      IF (.not.allocated(t_opt)) THEN
        allocate ( t_opt(Ngrids) )
      END IF
      IF (.not.allocated(I_opt)) THEN
        allocate ( I_opt(Ngrids) )
      END IF
      IF (.not.allocated(beta1)) THEN
        allocate ( beta1(Ngrids) )
      END IF
      IF (.not.allocated(beta2)) THEN
        allocate ( beta2(Ngrids) )
      END IF

      IF (.not.allocated(K_DO)) THEN
        allocate ( K_DO(Ngrids) )
      END IF
      IF (.not.allocated(H2SOR)) THEN
        allocate ( H2SOR(Ngrids) )
      END IF

      IF (.not.allocated(K_Nitri)) THEN
        allocate ( K_Nitri(Ngrids) )
      END IF
      IF (.not.allocated(DenitR)) THEN
        allocate ( DenitR(Ngrids) )
      END IF
      IF (.not.allocated(K_Denit)) THEN
        allocate ( K_Denit(Ngrids) )
      END IF

      IF (.not.allocated(NitriR_t)) THEN
        allocate ( NitriR_t(Ngrids) )
      END IF
      IF (.not.allocated(DenitR_t)) THEN
        allocate ( DenitR_t(Ngrids) )
      END IF
      IF (.not.allocated(PhyMR_t)) THEN
        allocate ( PhyMR_t(Ngrids) )
      END IF
      IF (.not.allocated(ZooGR_t)) THEN
        allocate ( ZooGR_t(Ngrids) )
      END IF
      IF (.not.allocated(RR_t)) THEN
        allocate ( RR_t(Ngrids) )
      END IF

      IF (.not.allocated(R_SODf)) THEN
        allocate ( R_SODf(Ngrids) )
      END IF
      IF (.not.allocated(R_NH4f)) THEN
        allocate ( R_NH4f(Ngrids) )
      END IF
      IF (.not.allocated(R_PO4f)) THEN
        allocate ( R_PO4f(Ngrids) )
      END IF
      IF (.not.allocated(R_NH4f_max)) THEN
        allocate ( R_NH4f_max(Ngrids) )
      END IF
      IF (.not.allocated(R_PO4f_max)) THEN
        allocate ( R_PO4f_max(Ngrids) )
      END IF
      IF (.not.allocated(K_DO_npflux)) THEN
        allocate ( K_DO_npflux(Ngrids) )
      END IF
      IF (.not.allocated(t_SODf)) THEN
        allocate ( t_SODf(Ngrids) )
      END IF

      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF
      IF (.not.allocated(AttSW)) THEN
        allocate ( AttSW(Ngrids) )
      END IF
      IF (.not.allocated(AttChl)) THEN
        allocate ( AttChl(Ngrids) )
      END IF
      IF (.not.allocated(Chl2C_m)) THEN
        allocate ( Chl2C_m(Ngrids) )
      END IF
      IF (.not.allocated(ChlMin)) THEN
        allocate ( ChlMin(Ngrids) )
      END IF
      IF (.not.allocated(CoagR)) THEN
        allocate ( CoagR(Ngrids) )
      END IF
      IF (.not.allocated(D_p5NH4)) THEN
        allocate ( D_p5NH4(Ngrids) )
      END IF
      IF (.not.allocated(I_thNH4)) THEN
        allocate ( I_thNH4(Ngrids) )
      END IF
      IF (.not.allocated(K_NH4)) THEN
        allocate ( K_NH4(Ngrids) )
      END IF
      IF (.not.allocated(K_NO3)) THEN
        allocate ( K_NO3(Ngrids) )
      END IF
      IF (.not.allocated(K_Phy)) THEN
        allocate ( K_Phy(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRN)) THEN
        allocate ( LDeRRN(Ngrids) )
      END IF
      IF (.not.allocated(LDeRRC)) THEN
        allocate ( LDeRRC(Ngrids) )
      END IF
      IF (.not.allocated(NitriR)) THEN
        allocate ( NitriR(Ngrids) )
      END IF
      IF (.not.allocated(PARfrac)) THEN
        allocate ( PARfrac(Ngrids) )
      END IF
      IF (.not.allocated(PhyCN)) THEN
        allocate ( PhyCN(Ngrids) )
      END IF
      IF (.not.allocated(PhyIP)) THEN
        allocate ( PhyIP(Ngrids) )
      END IF
      IF (.not.allocated(PhyIS)) THEN
        allocate ( PhyIS(Ngrids) )
      END IF
      IF (.not.allocated(PhyMin)) THEN
        allocate ( PhyMin(Ngrids) )
      END IF
      IF (.not.allocated(PhyMR)) THEN
        allocate ( PhyMR(Ngrids) )
      END IF
      IF (.not.allocated(SDeAR)) THEN
        allocate ( SDeAR(Ngrids) )
      END IF
      IF (.not.allocated(SDeBR)) THEN
        allocate ( SDeBR(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRN)) THEN
        allocate ( SDeRRN(Ngrids) )
      END IF
      IF (.not.allocated(SDeRRC)) THEN
        allocate ( SDeRRC(Ngrids) )
      END IF
      IF (.not.allocated(Vp0)) THEN
        allocate ( Vp0(Ngrids) )
      END IF
      IF (.not.allocated(wLDet)) THEN
        allocate ( wLDet(Ngrids) )
      END IF
      IF (.not.allocated(wPhy)) THEN
        allocate ( wPhy(Ngrids) )
      END IF
      IF (.not.allocated(wSDet)) THEN
        allocate ( wSDet(Ngrids) )
      END IF
      IF (.not.allocated(ZooAE_N)) THEN
        allocate ( ZooAE_N(Ngrids) )
      END IF
      IF (.not.allocated(ZooBM)) THEN
        allocate ( ZooBM(Ngrids) )
      END IF
      IF (.not.allocated(ZooCN)) THEN
        allocate ( ZooCN(Ngrids) )
      END IF
      IF (.not.allocated(ZooER)) THEN
        allocate ( ZooER(Ngrids) )
      END IF
      IF (.not.allocated(ZooGR)) THEN
        allocate ( ZooGR(Ngrids) )
      END IF
      IF (.not.allocated(ZooMin)) THEN
        allocate ( ZooMin(Ngrids) )
      END IF
      IF (.not.allocated(ZooMR)) THEN
        allocate ( ZooMR(Ngrids) )
      END IF
      IF (.not.allocated(pCO2air)) THEN
        allocate ( pCO2air(Ngrids) )
      END IF
#ifdef TANGENT
!      IF (.not.allocated(tl_wLDet)) THEN
!        allocate ( tl_wLDet(Ngrids) )
!      END IF
!      IF (.not.allocated(tl_wPhy)) THEN
!        allocate ( tl_wPhy(Ngrids) )
!      END IF
!      IF (.not.allocated(tl_wSDet)) THEN
!        allocate ( tl_wSDet(Ngrids) )
!      END IF
!      IF (.not.allocated(tl_PARfrac)) THEN
!        allocate ( tl_PARfrac(Ngrids) )
!      END IF
#endif
#ifdef ADJOINT
!      IF (.not.allocated(ad_wLDet)) THEN
!        allocate ( ad_wLDet(Ngrids) )
!      END IF
!      IF (.not.allocated(ad_wPhy)) THEN
!        allocate ( ad_wPhy(Ngrids) )
!      END IF
!      IF (.not.allocated(ad_wSDet)) THEN
!        allocate ( ad_wSDet(Ngrids) )
!      END IF
!      IF (.not.allocated(ad_PARfrac)) THEN
!        allocate ( ad_PARfrac(Ngrids) )
!      END IF
#endif
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
      IF (.not.allocated(iDbio3)) THEN
        allocate ( iDbio3(NDbio3d) )
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iNO3_=ic+1
      iNH4_=ic+2
      iChlo=ic+3
      iPhyt=ic+4
      iZoop=ic+5
      iLDeN=ic+6
      iSDeN=ic+7
      ic=ic+7
# ifdef CARBON
      iLDeC=ic+1
      iSDeC=ic+2
      iTIC_=ic+3
      iTAlk=ic+4
      ic=ic+4
# endif
# ifdef OXYGEN
      iOxyg=ic+1
      ic=ic+1
# endif

# ifdef PHOSPHORUS
      iPO4_=ic+1
      iLDeP=ic+2
      iSDeP=ic+3
      ic=ic+3
# endif
# ifdef H2S
      iH2S_=ic+1
      ic=ic+1
# endif

      RETURN
      END SUBROUTINE initialize_biology
