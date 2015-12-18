!
!svn $Id$
!==================================================== OKADA Teruhisa ===
!  Copyright (c) 2012-2014 The Water Eng. Group Osaka Univ.            !
!    Licensed under a MIT/X style license                              !
!=======================================================================
!                                                                      !
!  Bgc Model Kernel Variables:                                    !
!                                                                      !
!  bed            Sediment properties in each bed layer:               !
!                   bed(:,:,:,ithck) => layer thickness                !
!                   bed(:,:,:,iaged) => layer age                      !
!                   bed(:,:,:,iporo) => layer porosity                 !
!                   bed(:,:,:,idiff) => layer bio-diffusivity          !
!  bottom         Exposed sediment layer properties:                   !
!                   bottom(:,:,isd50) => mean grain diameter           !
!                   bottom(:,:,idens) => mean grain density            !
!                   bottom(:,:,iwsed) => mean settling velocity        !
!                   bottom(:,:,itauc) => mean critical erosion stress  !
!                   bottom(:,:,irlen) => ripple length                 !
!                   bottom(:,:,irhgt) => ripple height                 !
!                   bottom(:,:,ibwav) => bed wave excursion amplitude  !
!                   bottom(:,:,izNik) => Nikuradse bottom roughness    !
!                   bottom(:,:,izbio) => biological bottom roughness   !
!                   bottom(:,:,izbfm) => bed form bottom roughness     !
!                   bottom(:,:,izbld) => bed load bottom roughness     !
!                   bottom(:,:,izapp) => apparent bottom roughness     !
!                   bottom(:,:,izwbl) => wave bottom roughness         !
!                   bottom(:,:,izdef) => default bottom roughness      !
!                   bottom(:,:,iactv) => active layer thickness        !
!                   bottom(:,:,ishgt) => saltation height              !
!  ero_flux       Flux from erosion.                                   !
!  settling_flux  Flux from settling.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
      TYPE T_BGCBED
!
!  Nonlinear model state.
!
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        real(r8), pointer :: avgpw(:,:,:,:)
        real(r8), pointer :: avgsm(:,:,:,:)
        real(r8), pointer :: avgpwf(:,:,:)
        real(r8), pointer :: avgsmf(:,:,:)
#endif
        real(r8), pointer :: bpw(:,:,:,:)
        real(r8), pointer :: bsm(:,:,:,:)
        real(r8), pointer :: bpwflux(:,:,:)
        real(r8), pointer :: bsmflux(:,:,:)

      END TYPE T_BGCBED

      TYPE (T_BGCBED), allocatable :: BGCBED(:)

      CONTAINS

      SUBROUTINE allocate_bgcbed (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate structure variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( BGCBED(Ngrids) )
!
!  Nonlinear model state.
!
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      IF (ANY(Aout(idBpw(:),ng))) THEN
        allocate ( BGCBED(ng) % avgpw(LBi:UBi,LBj:UBj,Nbed,NBGCPW) )
      END IF
      IF (ANY(Aout(idBsm(:),ng))) THEN
        allocate ( BGCBED(ng) % avgsm(LBi:UBi,LBj:UBj,Nbed,NBGCSM) )
      END IF
      IF (ANY(Aout(idFpw(:),ng))) THEN
        allocate ( BGCBED(ng) % avgpwf(LBi:UBi,LBj:UBj,NBGCPW) )
      END IF
      IF (ANY(Aout(idFsm(:),ng))) THEN
        allocate ( BGCBED(ng) % avgsmf(LBi:UBi,LBj:UBj,NBGCSM) )
      END IF
#endif
      allocate ( BGCBED(ng) % bpw(LBi:UBi,LBj:UBj,Nbed,NBGCPW) )
      allocate ( BGCBED(ng) % bsm(LBi:UBi,LBj:UBj,Nbed,NBGCSM) )
      allocate ( BGCBED(ng) % bpwflux(LBi:UBi,LBj:UBj,NBGCPW) )
      allocate ( BGCBED(ng) % bsmflux(LBi:UBi,LBj:UBj,NBGCSM) )

      RETURN
      END SUBROUTINE allocate_bgcbed

      SUBROUTINE initialize_bgcbed (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine initialize structure variables in the module using     !
!  first touch distribution policy. In shared-memory configuration,    !
!  this operation actually performs the propagation of the  shared     !
!  arrays  across the cluster,  unless another policy is specified     !
!  to  override the default.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
      USE mod_biology
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: i, itrc, j, k

      real(r8), parameter :: IniVal = 0.0_r8

#include "set_bounds.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=BOUNDS(ng)%LBi(tile)
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=BOUNDS(ng)%UBi(tile)
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=BOUNDS(ng)%LBj(tile)
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=BOUNDS(ng)%UBj(tile)
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=BOUNDS(ng)%LBi(tile)
      Imax=BOUNDS(ng)%UBi(tile)
      Jmin=BOUNDS(ng)%LBj(tile)
      Jmax=BOUNDS(ng)%UBj(tile)
#endif
!
!-----------------------------------------------------------------------
!  Initialize sediment structure variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN

#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
        IF (ANY(Aout(idBpw(:),ng))) THEN
          DO itrc=1,NBGCPW
            DO k=1,Nbed
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  BGCBED(ng) % avgpw(i,j,k,itrc) = IniVal
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (ANY(Aout(idBsm(:),ng))) THEN
          DO itrc=1,NBGCSM
            DO k=1,Nbed
              DO j=Jmin,Jmax
                DO i=Imin,Imax
                  BGCBED(ng) % avgsm(i,j,k,itrc) = IniVal
                END DO
              END DO
            END DO
          END DO
        END IF
        IF (ANY(Aout(idFpw(:),ng))) THEN
          DO itrc=1,NBGCPW
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                BGCBED(ng) % avgpwf(i,j,itrc) = IniVal
              END DO
            END DO
          END DO
        END IF
        IF (ANY(Aout(idFsm(:),ng))) THEN
          DO itrc=1,NBGCSM
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                BGCBED(ng) % avgsmf(i,j,itrc) = IniVal
              END DO
            END DO
          END DO
        END IF
#endif
        DO j=Jmin,Jmax

          DO itrc=1,NBGCPW
            DO k=1,Nbed
              DO i=Imin,Imax
                BGCBED(ng) % bpw(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,NBGCSM
            DO k=1,Nbed
              DO i=Imin,Imax
                BGCBED(ng) % bsm(i,j,k,itrc) = IniVal
              END DO
            END DO
          END DO
          DO itrc=1,NBGCPW
            DO i=Imin,Imax
              BGCBED(ng) % bpwflux(i,j,itrc) = IniVal
            END DO
          END DO
          DO itrc=1,NBGCSM
            DO i=Imin,Imax
              BGCBED(ng) % bsmflux(i,j,itrc) = IniVal
            END DO
          END DO
          
        END DO

      END IF

      RETURN
      END SUBROUTINE initialize_bgcbed
