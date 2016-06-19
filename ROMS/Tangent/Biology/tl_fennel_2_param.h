#include "tl_fennel_bs2.h"
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
          tl_fac1=dtdays*tl_ZooGR
          tl_cff2=dtdays*tl_PhyMR

          DO k=1,N(ng)
            DO i=Istr,Iend
#ifdef TDEPENDANCE
              fac2=ZooGR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              fac3=PhyMR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              tl_fac2=fac2*tl_Bio(i,k,itemp)*LOG(ZooGR_t(ng))
              tl_fac3=fac3*tl_Bio(i,k,itemp)*LOG(PhyMR_t(ng))

              fac1=dtdays*ZooGR(ng)*fac2
              cff2=dtdays*PhyMR(ng)*fac3
              tl_fac1=dtdays*(tl_ZooGR*fac2+ZooGR(ng)*tl_fac2)
              tl_cff2=dtdays*(tl_PhyMR*fac3+PhyMR(ng)*tl_fac3)
#endif
!
! Phytoplankton grazing by zooplankton.
!
              fac=K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt)
              tl_fac=tl_K_Phy+2.0_r8*Bio1(i,k,iPhyt)*tl_Bio(i,k,iPhyt)

!!            cff=Bio1(i,k,iPhyt)/                                      &
!!   &            (K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              cff=Bio1(i,k,iPhyt)/fac
              tl_cff=(tl_Bio(i,k,iPhyt)-cff*tl_fac)/fac

              cff1=fac1*Bio1(i,k,iZoop)*cff
              tl_cff1=tl_fac1*Bio1(i,k,iZoop)*cff+                      &
     &                fac1*tl_Bio(i,k,iZoop)*cff+                       &
     &                fac1*Bio1(i,k,iZoop)*tl_cff

              cff3=1.0_r8/(1.0_r8+cff1)
              tl_cff3=-cff3*cff3*tl_cff1

!>            Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!>            Bio1(i,k,iChlo)=Bio(i,k,iChlo)
!>            Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
!>            Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
              tl_Bio(i,k,iPhyt)=tl_cff3*Bio2(i,k,iPhyt)+                &
     &                          cff3*tl_Bio(i,k,iPhyt)
              tl_Bio(i,k,iChlo)=tl_cff3*Bio2(i,k,iChlo)+                &
     &                          cff3*tl_Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
              fac4=Bio2(i,k,iPhyt)*cff1
              tl_fac4=tl_Bio(i,k,iPhyt)*cff1+Bio(i,k,iPhyt)*tl_cff1

!!            N_Flux_Assim=Bio2(i,k,iPhyt)*cff1*ZooAE_N(ng)
              N_Flux_Assim=fac4*ZooAE_N(ng)
              tl_N_Flux_Assim=tl_fac4*ZooAE_N(ng)+fac4*tl_ZooAE_N

!!            N_Flux_Egest=Bio2(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
              N_Flux_Egest=fac4*(1.0_r8-ZooAE_N(ng))
              tl_N_Flux_Egest=tl_fac4*(1.0_r8-ZooAE_N(ng))-             &
     &                        fac4*tl_ZooAE_N

!>            Bio1(i,k,iZoop)=Bio(i,k,iZoop)
!>            Bio1(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest
              tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)+tl_N_Flux_Assim
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              cff=MAX(Bio2(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              tl_cff=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iPhyt)-PhyMin(ng)))*  &
     &               tl_Bio(i,k,iPhyt)

              N_Flux_Pmortal=cff2*cff
              tl_N_Flux_Pmortal=tl_cff2*cff+cff2*tl_cff

!>            Bio2(i,k,iPhyt)=Bio(i,k,iPhyt)
!>            Bio2(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal
              tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)-tl_N_Flux_Pmortal
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Pmortal

              cff=MAX(Bio2(i,k,iChlo)-ChlMin(ng),0.0_r8)
              tl_cff=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iChlo)-ChlMin(ng)))*  &
     &               tl_Bio(i,k,iChlo)

!>            Bio2(i,k,iChlo)=Bio(i,k,iChlo)
!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)-cff2*cff
              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)-                      &
     &                          (tl_cff2*cff+cff2*tl_cff)
#ifdef PHOSPHORUS
              cff1=PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)
              tl_cff1=tl_PhyPN*(N_Flux_Egest+N_Flux_Pmortal)+           &
     &                PhyPN(ng)*(tl_N_Flux_Egest+tl_N_Flux_Pmortal)

              cff2=(PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
              tl_cff2=(tl_PhyPN-tl_ZooPN)*N_Flux_Assim+                 &
     &                (PhyPN(ng)-ZooPN(ng))*tl_N_Flux_Assim

!>            Bio1(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+cff1+cff2
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+tl_cff1+tl_cff2
#endif
            END DO
          END DO