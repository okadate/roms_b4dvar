#include "fennel_bs2.h"
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
              fac4=K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt)
              fac5=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)
              fac1=fac5/fac4
              cff2=fac2*Bio2(i,k,iZoop)
              cff3=fac1*ZooAE_N(ng)
!!            Bio2(i,k,iZoop)=Bio(i,k,iZoop)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
              N_Flux_Zmortal=cff2*Bio3(i,k,iZoop)
              N_Flux_Zexcret=cff3*Bio3(i,k,iZoop)
!!            Bio1(i,k,iNH4_)=Bio(i,k,iNH4_)
!!            Bio3(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
#ifdef PHOSPHORUS
!!            Bio1(i,k,iPO4_)=Bio(i,k,iPO4_)
!!            Bio2(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zexcret
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+ZooPN(ng)*N_Flux_Zmortal
#endif
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
              fac6=MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              N_Flux_Zmetabo=cff1*fac6
!!            Bio3(i,k,iZoop)=Bio(i,k,iZoop)
!!            Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef PHOSPHORUS
!!            Bio2(i,k,iPO4_)=Bio(i,k,iPO4_)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zmetabo
#endif
#ifdef OXYGEN
!!            Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
!>   &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
#endif
!             ==========================================================
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
#ifdef OXYGEN
!>            tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-rOxNH4*               &
!>   &                          (tl_N_Flux_Zmetabo+tl_N_Flux_Zexcret)
              adfac=rOxNH4*ad_Bio(i,k,iOxyg)
              ad_N_Flux_Zexcret=ad_N_Flux_Zexcret-adfac
              ad_N_Flux_Zmetabo=ad_N_Flux_Zmetabo-adfac
#endif
#ifdef PHOSPHORUS
!>            tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+                      &
!>   &                          ZooPN(ng)*tl_N_Flux_Zmetabo+            &
!>   &                          tl_ZooPN*N_Flux_Zmetabo
              ad_N_Flux_Zmetabo=ad_N_Flux_Zmetabo+                      &
     &                          ZooPN(ng)*ad_Bio(i,k,iPO4_)
              ad_ZooPN=ad_ZooPN+N_Flux_Zmetabo*ad_Bio(i,k,iPO4_)
#endif
!>            tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Zmetabo
              ad_N_Flux_Zmetabo=ad_N_Flux_Zmetabo+ad_Bio(i,k,iNH4_)

!>            tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)-tl_N_Flux_Zmetabo
              ad_N_Flux_Zmetabo=ad_N_Flux_Zmetabo-ad_Bio(i,k,iZoop)

!>            tl_N_Flux_Zmetabo=tl_cff1*fac6+cff1*tl_fac6
              ad_fac6=ad_fac6+cff1*ad_N_Flux_Zmetabo
              ad_cff1=ad_cff1+fac6*ad_N_Flux_Zmetabo
              ad_N_Flux_Zmetabo=0.0_r8

!>            tl_fac6=(0.5_r8+SIGN(0.5_r8,Bio3(i,k,iZoop)-ZooMin(ng)))* &
!>   &                tl_Bio(i,k,iZoop)
              adfac=(0.5_r8+SIGN(0.5_r8,Bio3(i,k,iZoop)-ZooMin(ng)))
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+adfac*ad_fac6
              ad_fac6=0.0_r8
!
!  Zooplankton mortality and excretion.
!
#ifdef PHOSPHORUS
!>            tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+                      &
!>   &                          tl_ZooPN*N_Flux_Zmortal+                &
!>   &                          ZooPN(ng)*tl_N_Flux_Zmortal
              ad_N_Flux_Zmortal=ad_N_Flux_Zmortal+                      &
     &                          ZooPN(ng)*ad_Bio(i,k,iSDeP)
              ad_ZooPN=ad_ZooPN+N_Flux_Zmortal*ad_Bio(i,k,iSDeP)

!>            tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+                      &
!>   &                          tl_ZooPN*N_Flux_Zexcret+                &
!>   &                          ZooPN(ng)*tl_N_Flux_Zexcret
              ad_N_Flux_Zexcret=ad_N_Flux_Zexcret+                      &
     &                          ZooPN(ng)*ad_Bio(i,k,iPO4_)
              ad_ZooPN=ad_ZooPN+N_Flux_Zexcret*ad_Bio(i,k,iPO4_)
#endif
!>            tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Zmortal
              ad_N_Flux_Zmortal=ad_N_Flux_Zmortal+ad_Bio(i,k,iSDeN)

!>            tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Zexcret
              ad_N_Flux_Zexcret=ad_N_Flux_Zexcret+ad_Bio(i,k,iNH4_)

!>            tl_N_Flux_Zexcret=tl_cff3*Bio3(i,k,iZoop)+                &
!>   &                          cff3*tl_Bio(i,k,iZoop)
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+cff3*ad_N_Flux_Zexcret
              ad_cff3=ad_cff3+Bio3(i,k,iZoop)*ad_N_Flux_Zexcret
              ad_N_Flux_Zexcret=0.0_r8

!>            tl_N_Flux_Zmortal=tl_cff2*Bio3(i,k,iZoop)+                &
!>   &                          cff2*tl_Bio(i,k,iZoop)
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+cff2*ad_N_Flux_Zmortal
              ad_cff2=ad_cff2+Bio3(i,k,iZoop)*ad_N_Flux_Zmortal
              ad_N_Flux_Zmortal=0.0_r8
!
!  Uptake zooplankton
!
!>            tl_Bio(i,k,iZoop)=(tl_Bio(i,k,iZoop)-Bio3(i,k,iZoop)*     &
!>   &                           (tl_cff2+tl_cff3))/(1.0_r8+cff2+cff3)
              adfac=ad_Bio(i,k,iZoop)/(1.0_r8+cff2+cff3)
              ad_cff3=ad_cff3-Bio3(i,k,iZoop)*adfac
              ad_cff2=ad_cff2-Bio3(i,k,iZoop)*adfac
              ad_Bio(i,k,iZoop)=adfac

!>            tl_cff3=tl_fac1*ZooAE_N(ng)+fac1*tl_ZooAE_N
              ad_fac1=ad_fac1+ZooAE_N(ng)*ad_cff3
              ad_ZooAE_N=ad_ZooAE_N+fac1*ad_cff3
              ad_cff3=0.0_r8

!>            tl_cff2=tl_fac2*Bio2(i,k,iZoop)+fac2*tl_Bio(i,k,iZoop)
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+fac2*ad_cff2
              ad_fac2=ad_fac2+Bio2(i,k,iZoop)*ad_cff2
              ad_cff2=0.0_r8

!>            tl_fac1=(tl_fac5-fac1*tl_fac4)/fac4
              ad_fac4=ad_fac4-fac1*ad_fac1/fac4
              ad_fac5=ad_fac5+ad_fac1/fac4
              ad_fac1=0.0_r8


!>            tl_fac5=tl_fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)+            &
!>   &                fac3*2.0_r8*Bio(i,k,iPhyt)*tl_Bio(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          fac3*2.0_r8*Bio(i,k,iPhyt)*ad_fac5
              ad_fac3=ad_fac3+Bio(i,k,iPhyt)*Bio(i,k,iPhyt)*ad_fac5
              ad_fac5=0.0_r8

!>            tl_fac4=tl_K_Phy+2.0_r8*Bio(i,k,iPhyt)*tl_Bio(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          2.0_r8*Bio(i,k,iPhyt)*ad_fac4
              ad_K_Phy=ad_K_Phy+ad_fac4
              ad_fac4=0.0_r8
            END DO
          END DO
!>        tl_cff1=dtdays*tl_ZooBM
!>        tl_fac2=dtdays*tl_ZooMR
!>        tl_fac3=dtdays*tl_ZooER
          ad_ZooBM=ad_ZooBM+dtdays*ad_cff1
          ad_ZooMR=ad_ZooMR+dtdays*ad_fac2
          ad_ZooER=ad_ZooER+dtdays*ad_fac3
          ad_cff1=0.0_r8
          ad_fac2=0.0_r8
          ad_fac3=0.0_r8
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
#ifdef TDEPENDANCE
              fac2=ZooGR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              fac3=PhyMR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              fac1=dtdays*ZooGR(ng)*fac2
              cff2=dtdays*PhyMR(ng)*fac3
#endif
!
! Phytoplankton grazing by zooplankton.
!
              fac=K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt)
              cff=Bio1(i,k,iPhyt)/fac
              cff1=fac1*Bio1(i,k,iZoop)*cff
              cff3=1.0_r8/(1.0_r8+cff1)
!>            Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!>            Bio1(i,k,iChlo)=Bio(i,k,iChlo)
!>            Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
!>            Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
              fac4=Bio2(i,k,iPhyt)*cff1
              N_Flux_Assim=fac4*ZooAE_N(ng)
              N_Flux_Egest=fac4*(1.0_r8-ZooAE_N(ng))
!>            Bio1(i,k,iZoop)=Bio(i,k,iZoop)
!>            Bio1(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              cff=MAX(Bio2(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              N_Flux_Pmortal=cff2*cff
!>            Bio2(i,k,iPhyt)=Bio(i,k,iPhyt)
!>            Bio2(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal
              cff=MAX(Bio2(i,k,iChlo)-ChlMin(ng),0.0_r8)
!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)-cff2*cff
#ifdef PHOSPHORUS
!!            Bio1(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
!>   &                       PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)+   &
!>   &                       (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
              cff1=PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)
              cff2=(PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+cff1+cff2
#endif
!             ==========================================================
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
#ifdef PHOSPHORUS
!>            tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)+tl_cff1+tl_cff2
              ad_cff2=ad_cff2+ad_Bio(i,k,iSDeP)
              ad_cff1=ad_cff1+ad_Bio(i,k,iSDeP)

!>            tl_cff2=(tl_PhyPN-tl_ZooPN)*N_Flux_Assim+                 &
!>   &                (PhyPN(ng)-ZooPN(ng))*tl_N_Flux_Assim
              ad_N_Flux_Assim=ad_N_Flux_Assim+                          &
     &                        (PhyPN(ng)-ZooPN(ng))*ad_cff2
              adfac=N_Flux_Assim*ad_cff2
              ad_ZooPN=ad_ZooPN-adfac
              ad_PhyPN=ad_PhyPN+adfac
              ad_cff2=0.0_r8

!>            tl_cff1=tl_PhyPN*(N_Flux_Egest+N_Flux_Pmortal)+           &
!>   &                PhyPN(ng)*(tl_N_Flux_Egest+tl_N_Flux_Pmortal)
              adfac=PhyPN(ng)*ad_cff1
              ad_N_Flux_Pmortal=ad_N_Flux_Pmortal+adfac
              ad_N_Flux_Egest=ad_N_Flux_Egest+adfac
              ad_PhyPN=ad_PhyPN+(N_Flux_Egest+N_Flux_Pmortal)*ad_cff1
              ad_cff1=0.0_r8
#endif
!>            tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)-                      &
!>   &                          (tl_cff2*cff+cff2*tl_cff)
              ad_cff=ad_cff-cff2*ad_Bio(i,k,iChlo)
              ad_cff2=ad_cff2-ad_Bio(i,k,iChlo)*cff

!>            tl_cff=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iChlo)-ChlMin(ng)))*  &
!>   &               tl_Bio(i,k,iChlo)
              adfac=0.5_r8+SIGN(0.5_r8,Bio2(i,k,iChlo)-ChlMin(ng))
              ad_Bio(i,k,iChlo)=ad_Bio(i,k,iChlo)+adfac*ad_cff
              ad_cff=0.0_r8

!>            tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Pmortal
              ad_N_Flux_Pmortal=ad_N_Flux_Pmortal+ad_Bio(i,k,iSDeN)

!>            tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)-tl_N_Flux_Pmortal
              ad_N_Flux_Pmortal=ad_N_Flux_Pmortal-ad_Bio(i,k,iPhyt)

!             ==========================================================
              cff=MAX(Bio2(i,k,iPhyt)-PhyMin(ng),0.0_r8)
!             ==========================================================

!>            tl_N_Flux_Pmortal=tl_cff2*cff+cff2*tl_cff
              ad_cff=ad_cff+cff2*ad_N_Flux_Pmortal
              ad_cff2=ad_cff2+cff*ad_N_Flux_Pmortal
              ad_N_Flux_Pmortal=0.0_r8

!>            tl_cff=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iPhyt)-PhyMin(ng)))*  &
!>   &               tl_Bio(i,k,iPhyt)
              adfac=0.5_r8+SIGN(0.5_r8,Bio2(i,k,iPhyt)-PhyMin(ng))
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+adfac*ad_cff
              ad_cff=0.0_r8
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
!>            tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)+tl_N_Flux_Egest
              ad_N_Flux_Egest=ad_N_Flux_Egest+ad_Bio(i,k,iSDeN)

!>            tl_Bio(i,k,iZoop)=tl_Bio(i,k,iZoop)+tl_N_Flux_Assim
              ad_N_Flux_Assim=ad_N_Flux_Assim+ad_Bio(i,k,iZoop)

!>            tl_N_Flux_Egest=tl_fac4*(1.0_r8-ZooAE_N(ng))-             &
!>   &                        fac4*tl_ZooAE_N
              ad_ZooAE_N=ad_ZooAE_N-fac4*ad_N_Flux_Egest
              ad_fac4=ad_fac4+(1.0_r8-ZooAE_N(ng))*ad_N_Flux_Egest
              ad_N_Flux_Egest=0.0_r8

!>            tl_N_Flux_Assim=tl_fac4*ZooAE_N(ng)+fac4*tl_ZooAE_N
              ad_ZooAE_N=ad_ZooAE_N+fac4*ad_N_Flux_Assim
              ad_fac4=ad_fac4+ZooAE_N(ng)*ad_N_Flux_Assim
              ad_N_Flux_Assim=0.0_r8

!>            tl_fac4=tl_Bio(i,k,iPhyt)*cff1+Bio2(i,k,iPhyt)*tl_cff1
              ad_cff1=ad_cff1+Bio2(i,k,iPhyt)*ad_fac4
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+cff1*ad_fac4
              ad_fac4=0.0_r8
!
! Phytoplankton grazing by zooplankton.
!
!             ==========================================================
              fac=K_Phy(ng)+Bio1(i,k,iPhyt)*Bio1(i,k,iPhyt)
              cff=Bio1(i,k,iPhyt)/fac
              cff1=fac1*Bio1(i,k,iZoop)*cff
              cff3=1.0_r8/(1.0_r8+cff1)
!             ==========================================================

!>            tl_Bio(i,k,iChlo)=tl_cff3*Bio1(i,k,iChlo)+                &
!>   &                          cff3*tl_Bio(i,k,iChlo)
              ad_cff3=ad_cff3+ad_Bio(i,k,iChlo)*Bio1(i,k,iChlo)
              ad_Bio(i,k,iChlo)=cff3*ad_Bio(i,k,iChlo)

!>            tl_Bio(i,k,iPhyt)=tl_cff3*Bio1(i,k,iPhyt)+                &
!>   &                          cff3*tl_Bio(i,k,iPhyt)
              ad_cff3=ad_cff3+ad_Bio(i,k,iPhyt)*Bio1(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=cff3*ad_Bio(i,k,iPhyt)

!>            tl_cff3=-cff3*cff3*tl_cff1
              ad_cff1=ad_cff1-cff3*cff3*ad_cff3
              ad_cff3=0.0_r8

!>            tl_cff1=tl_fac1*Bio1(i,k,iZoop)*cff+                      &
!>   &                fac1*tl_Bio(i,k,iZoop)*cff+                       &
!>   &                fac1*Bio1(i,k,iZoop)*tl_cff
              ad_cff=ad_cff+fac1*Bio1(i,k,iZoop)*ad_cff1
              ad_Bio(i,k,iZoop)=ad_Bio(i,k,iZoop)+fac1*ad_cff1*cff
              ad_fac1=ad_fac1+ad_cff1*Bio1(i,k,iZoop)*cff
              ad_cff1=0.0_r8

!>            tl_cff=(tl_Bio(i,k,iPhyt)-cff*tl_fac)/fac
              adfac=ad_cff/fac
              ad_fac=ad_fac-cff*adfac
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+adfac
              ad_cff=0.0_r8

!>            tl_fac=tl_K_Phy+2.0_r8*Bio1(i,k,iPhyt)*tl_Bio(i,k,iPhyt)
              ad_Bio(i,k,iPhyt)=ad_Bio(i,k,iPhyt)+                      &
     &                          2.0_r8*Bio1(i,k,iPhyt)*ad_fac
              ad_K_Phy=ad_K_Phy+ad_fac
              ad_fac=0.0_r8

#ifdef TDEPENDANCE
!>            tl_cff2=dtdays*(tl_PhyMR*fac3+PhyMR(ng)*tl_fac3)
              ad_fac3=ad_fac3+dtdays*PhyMR(ng)*ad_cff2
              ad_PhyMR=ad_PhyMR+dtdays*fac3*ad_cff2
              ad_cff2=0.0_r8

!>            tl_fac1=dtdays*(tl_ZooGR*fac2+ZooGR(ng)*tl_fac2)
              ad_fac2=ad_fac2+dtdays*ZooGR(ng)*ad_fac1
              ad_ZooGR=ad_ZooGR+dtdays*fac2*ad_fac1
              ad_fac1=0.0_r8

!>            tl_fac3=fac3*tl_Bio(i,k,itemp)*LOG(PhyMR_t(ng))
# ifndef UV_FIXED_TL
              ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+                      &
     &                          fac3*ad_fac3*LOG(PhyMR_t(ng))
# endif
              ad_fac3=0.0_r8

!>            tl_fac2=fac2*tl_Bio(i,k,itemp)*LOG(ZooGR_t(ng))
# ifndef UV_FIXED_TL
              ad_Bio(i,k,itemp)=ad_Bio(i,k,itemp)+                      &
     &                          fac2*ad_fac2*LOG(ZooGR_t(ng))
# endif
              ad_fac2=0.0_r8
#endif
            END DO
          END DO

!>        tl_fac1=dtdays*tl_ZooGR
          ad_ZooGR=ad_ZooGR+dtdays*ad_fac1
          ad_fac1=0.0_r8

!>        tl_cff2=dtdays*tl_PhyMR
          ad_PhyMR=ad_PhyMR+dtdays*ad_cff2
          ad_cff2=0.0_r8
