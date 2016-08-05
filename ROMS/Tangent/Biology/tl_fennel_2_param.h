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
!!            cff=Bio1(i,k,iPhyt)/                                      &
!!   &            (K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt))
              fac=K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt)
              tl_fac=tl_K_Phy+2.0_r8*Bio1(i,k,iPhyt)*tl_Bio(i,k,iPhyt)
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
              tl_Bio(i,k,iPhyt)=tl_cff3*Bio1(i,k,iPhyt)+                &
     &                          cff3*tl_Bio(i,k,iPhyt)
              tl_Bio(i,k,iChlo)=tl_cff3*Bio1(i,k,iChlo)+                &
     &                          cff3*tl_Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
!!            N_Flux_Assim=Bio2(i,k,iPhyt)*cff1*ZooAE_N(ng)
              fac4=Bio2(i,k,iPhyt)*cff1
              tl_fac4=tl_Bio(i,k,iPhyt)*cff1+Bio2(i,k,iPhyt)*tl_cff1
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
!!            N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
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

!!            Bio2(i,k,iChlo)=Bio(i,k,iChlo)
!!            Bio(i,k,iChlo)=Bio(i,k,iChlo)-                            &
!!   &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
              cff=MAX(Bio2(i,k,iChlo)-ChlMin(ng),0.0_r8)
              tl_cff=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iChlo)-ChlMin(ng)))*  &
     &               tl_Bio(i,k,iChlo)

!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)-cff2*cff
              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)-                      &
     &                          (tl_cff2*cff+cff2*tl_cff)
#ifdef PHOSPHORUS
!!            Bio1(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
!>   &                       PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)+   &
!>   &                       (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
              cff1=PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)
              tl_cff1=tl_PhyPN*(N_Flux_Egest+N_Flux_Pmortal)+           &
     &                PhyPN(ng)*(tl_N_Flux_Egest+tl_N_Flux_Pmortal)

              cff2=(PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
              tl_cff2=(tl_PhyPN-tl_ZooPN)*N_Flux_Assim+                 &
     &                (PhyPN(ng)-ZooPN(ng))*tl_N_Flux_Assim

!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+cff1+cff2
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+tl_cff1+tl_cff2
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
          tl_cff1=dtdays*tl_ZooBM
          tl_fac2=dtdays*tl_ZooMR
          tl_fac3=dtdays*tl_ZooER

          DO k=1,N(ng)
            DO i=Istr,Iend
!!            fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                  &
!!   &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              fac4=K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt)
              tl_fac4=tl_K_Phy+2.0_r8*Bio(i,k,iPhyt)*tl_Bio(i,k,iPhyt)
              fac5=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)
              tl_fac5=tl_fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)+            &
     &                fac3*2.0_r8*Bio(i,k,iPhyt)*tl_Bio(i,k,iPhyt)
              fac1=fac5/fac4
              tl_fac1=(tl_fac5-fac1*tl_fac4)/fac4

              cff2=fac2*Bio2(i,k,iZoop)
              tl_cff2=tl_fac2*Bio2(i,k,iZoop)+fac2*tl_Bio(i,k,iZoop)

              cff3=fac1*ZooAE_N(ng)
              tl_cff3=tl_fac1*ZooAE_N(ng)+fac1*tl_ZooAE_N

!!            Bio2(i,k,iZoop)=Bio(i,k,iZoop)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+cff2+cff3)
              tl_Bio(i,k,iZoop)=(tl_Bio(i,k,iZoop)-Bio3(i,k,iZoop)*     &
     &                          (tl_cff2+tl_cff3))/(1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
              N_Flux_Zmortal=cff2*Bio3(i,k,iZoop)
              N_Flux_Zexcret=cff3*Bio3(i,k,iZoop)
              tl_N_Flux_Zmortal=tl_cff2*Bio3(i,k,iZoop)+                &
     &                          cff2*tl_Bio(i,k,iZoop)
              tl_N_Flux_Zexcret=tl_cff3*Bio3(i,k,iZoop)+                &
     &                          cff3*tl_Bio(i,k,iZoop)

!!            Bio1(i,k,iNH4_)=Bio(i,k,iNH4_)
!!            Bio3(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Zexcret
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Zmortal
#ifdef PHOSPHORUS
!!            Bio1(i,k,iPO4_)=Bio(i,k,iPO4_)
!!            Bio2(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zexcret
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+ZooPN(ng)*N_Flux_Zmortal
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+                      &
     &                          tl_ZooPN*N_Flux_Zexcret+                &
     &                          ZooPN(ng)*tl_N_Flux_Zexcret)
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+                      &
     &                          tl_ZooPN*N_Flux_Zmortal+                &
     &                          ZooPN(ng)*tl_N_Flux_Zmortal
#endif
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
!!            N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              fac6=MAX(Bio3(i,k,iZoop)-ZooMin(ng),0.0_r8)
              tl_fac6=(0.5_r8+SIGN(0.5_r8,Bio3(i,k,iZoop)-ZooMin(ng)))* &
     &                tl_Bio(i,k,iZoop)
              N_Flux_Zmetabo=cff1*fac6
              tl_N_Flux_Zmetabo=tl_cff1*fac6+cff1*tl_fac6

!!            Bio3(i,k,iZoop)=Bio(i,k,iZoop)
!!            Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
              tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)-tl_N_Flux_Zmetabo
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Zmetabo
#ifdef PHOSPHORUS
!!            Bio2(i,k,iPO4_)=Bio(i,k,iPO4_)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zmetabo
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+                      &
     &                          ZooPN(ng)*tl_N_Flux_Zmetabo+            &
     &                          tl_ZooPN*N_Flux_Zmetabo
#endif
#ifdef OXYGEN
!!            Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
!>   &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-rOxNH4*               &
     &                          (tl_N_Flux_Zmetabo+tl_N_Flux_Zexcret)
#endif
            END DO
          END DO