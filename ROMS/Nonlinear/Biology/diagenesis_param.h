!
!------------------------------------------------------------------------
!  Biogeochemical model variables and parameters.
!------------------------------------------------------------------------
!
      integer, parameter :: BgcIter = 1

      real(8), parameter :: DBL    = 0.03d0
      real(8), parameter :: dens   = 2.6d0

      real(8), parameter :: bwSO4  = 27300.0d0
      real(8), parameter :: bwMn   = 0.0d0
      real(8), parameter :: bwFe   = 0.0d0
      real(8), parameter :: bwCH4  = 0.0d0
      real(8), parameter :: bwDOMf = 118.0d0
      real(8), parameter :: bwDOMs = 0.0d0

      real(8), parameter, dimension(2) :: a_poro = (/0.78d0, 0.71d0/)
      real(8), parameter, dimension(2) :: b_poro = (/0.10d0, 0.00d0/)
      real(8), parameter, dimension(2) :: c_poro = (/0.35d0, 0.00d0/)
      real(8), parameter, dimension(2) :: a_wsm  = (/0.22d0, 0.20d0/)
      real(8), parameter, dimension(2) :: b_wsm  = (/0.26d0, 0.06d0/)
      real(8), parameter, dimension(2) :: c_wsm  = (/0.30d0, 0.30d0/)

      real(8), parameter :: D0O2  = 11.7d0
      real(8), parameter :: aO2   = 0.344d0
      real(8), parameter :: bO2   = 0.00505d0
      real(8), parameter :: D0NO3 = 9.72d0
      real(8), parameter :: aNO3  = 0.365d0
      real(8), parameter :: D0H2S = 8.74d0
      real(8), parameter :: aH2S  = 0.264d0
      real(8), parameter :: bH2S  = 0.004d0
      real(8), parameter :: D0SO4 = 4.96d0
      real(8), parameter :: aSO4  = 0.226d0
      real(8), parameter :: D0NH4 = 9.76d0
      real(8), parameter :: aNH4  = 0.398d0
      real(8), parameter :: D0Mn  = 3.04d0
      real(8), parameter :: aMn   = 0.153d0
      real(8), parameter :: D0Fe  = 3.36d0
      real(8), parameter :: aFe   = 0.148d0
      real(8), parameter :: D0PO4 = 9.76d0
      real(8), parameter :: aPO4  = 0.398d0
      real(8), parameter :: D0DOMf = 3.519d0
      real(8), parameter :: aDOMf  = 5.398d-2
      real(8), parameter :: D0DOMs = 2.069d-2
      real(8), parameter :: aDOMs  = 6.541d-3

      real(8), parameter :: z_DBw = 11.8d0
      real(8), parameter :: u_DBw = 3.51d-6
      real(8), parameter :: a_DBw = -0.378d0
      real(8), parameter :: a_DBs = 9.3d0
      real(8), parameter :: a_irr = 0.885d0
      real(8), parameter :: b_irr = 0.054d0
      real(8), parameter :: c_irr = 2.53d0
      real(8), parameter :: d_irr = 0.352d0
      real(8), parameter :: e_irr = 6.0d0
      real(8), parameter :: f_irr = 0.05d0

      real(8), parameter :: KdNH4 = 1.5d0
      real(8), parameter :: KdNO3 = 5.4d0
      real(8), parameter :: KdMn  = 13.0d0
      real(8), parameter :: KdFe  = 500.0d0
      real(8), parameter :: KdPO4 = 2.0d0

      real(8), parameter :: ratio_n    = 0.19d0
      real(8), parameter :: ratio_f    = 0.54d0
      real(8), parameter :: ratio_FA   = 0.5d0
      real(8), parameter :: ratio_MA   = 0.5d0
      real(8), parameter :: ratio_CN   = 8.00d0
      real(8), parameter :: ratio_CP   = 208.0d0
      real(8), parameter :: ratio_DOMf = 0.50d0

      real(8), parameter :: FMnO2  = 2.0d-2
      real(8), parameter :: FFeOOH = 1.0d0
      real(8), parameter :: KO2    = 20.0d0
      real(8), parameter :: KNO3   = 5.0d0
      real(8), parameter :: KMnO2  = 50000.0d0
      real(8), parameter :: KFeOOH = 100000.0d0

      real(8), parameter :: KPOMf = 2.514d-6
      real(8), parameter :: KPOMs = 1.738d-9
      real(8), parameter :: KDOMf = 0.001365d0
      real(8), parameter :: KDOMs = 6.255d-8

      real(8), parameter :: K06 = 2.5d-7
      real(8), parameter :: K07 = 5.0d-14
      real(8), parameter :: K08 = 1.7d-9
      real(8), parameter :: K09 = 1.5d-5
      real(8), parameter :: K10 = 2.0d-7
      real(8), parameter :: K11 = 5.0d-4
      real(8), parameter :: K12 = 3.0d-9
      real(8), parameter :: K13 = 3.75d-5
      real(8), parameter :: K14 = 3.0d-10
      real(8), parameter :: K15 = 7.5d-12
      real(8), parameter :: K16 = 5.0d-5
      real(8), parameter :: K17 = 6.0d-7
      real(8), parameter :: K18 = 3.0d-10
      real(8), parameter :: K19 = 7.0d-7
      real(8), parameter :: K20 = 1.3d-9
      real(8), parameter :: K21 = 9.0d-10

      real(8), parameter :: H2Sstop = 10.0d0
