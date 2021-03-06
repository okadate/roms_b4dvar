#include "fennel_bs3.h"
!
!-----------------------------------------------------------------------
!  Coagulation of phytoplankton and small detritus to large detritus.
!-----------------------------------------------------------------------
!
          fac1=dtdays*CoagR(ng)
          tl_fac1=dtdays*tl_CoagR

          DO k=1,N(ng)
            DO i=Istr,Iend
              cff1=fac1*(Bio1(i,k,iSDeN)+Bio1(i,k,iPhyt))
              tl_cff1=tl_fac1*(Bio1(i,k,iSDeN)+Bio1(i,k,iPhyt))+          &
     &                fac1*(tl_Bio(i,k,iSDeN)+tl_Bio(i,k,iPhyt))

              cff2=1.0_r8/(1.0_r8+cff1)
              tl_cff2=-cff2*cff2*tl_cff1

!!            Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!!            Bio1(i,k,iChlo)=Bio(i,k,iChlo)
!!            Bio1(i,k,iSDeN)=Bio(i,k,iSDeN)
!>            Bio(i,k,iPhyt)=Bio(i,k,iPhyt)*cff2
!>            Bio(i,k,iChlo)=Bio(i,k,iChlo)*cff2
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
              tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)*cff2+                 &
     &                          Bio1(i,k,iPhyt)*tl_cff2
              tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)*cff2+                 &
     &                          Bio1(i,k,iChlo)*tl_cff2
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2+                 &
     &                          Bio1(i,k,iSDeN)*tl_cff2

              N_Flux_CoagP=Bio2(i,k,iPhyt)*cff1
              N_Flux_CoagD=Bio2(i,k,iSDeN)*cff1
              tl_N_Flux_CoagP=tl_Bio(i,k,iPhyt)*cff1+                   &
     &                        Bio2(i,k,iPhyt)*tl_cff1
              tl_N_Flux_CoagD=tl_Bio(i,k,iSDeN)*cff1+                   &
     &                        Bio2(i,k,iSDeN)*tl_cff1

!!            Bio1(i,k,iLDeN)=Bio(i,k,iLDeN)
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
!>   &                       N_Flux_CoagP+N_Flux_CoagD
              tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)+                      &
     &                          tl_N_Flux_CoagP+tl_N_Flux_CoagD
#ifdef PHOSPHORUS
!!            Bio1(i,k,iSDeP)=Bio(i,k,iSDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)-                      &
     &                          (tl_PhyPN*N_Flux_CoagD+                 &
     &                           PhyPN(ng)*tl_N_Flux_CoagD)

!!            Bio1(i,k,iLDeP)=Bio(i,k,iLDeP)
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                            &
!>   &                       PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
              tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)+PhyPN(ng)*            &
     &                          (tl_N_Flux_CoagP+tl_N_Flux_CoagD)+      &
     &                          tl_PhyPN*(N_Flux_CoagP+N_Flux_CoagD)
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
              fac1=MAX(Bio2(i,k,iOxyg)-6.0_r8,0.0_r8)
              tl_fac1=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)-6.0_r8))*     &
     &                tl_Bio(i,k,iOxyg)

              fac2=fac1/(K_DO(ng)+fac1)
              tl_fac2=(tl_fac1-fac2*(tl_K_DO+tl_fac1))/(K_DO(ng)+fac1)

              cff1=dtdays*SDeRRN(ng)*fac2
              tl_cff1=dtdays*(tl_SDeRRN*fac2+SDeRRN(ng)*tl_fac2)

              cff2=1.0_r8/(1.0_r8+cff1)
              tl_cff2=-cff2*cff2*tl_cff1

              cff3=dtdays*LDeRRN(ng)*fac2
              tl_cff3=dtdays*(tl_LDeRRN*fac2+LDeRRN(ng)*tl_fac2)

              cff4=1.0_r8/(1.0_r8+cff3)
              tl_cff4=-cff4*cff4*tl_cff3

!!            Bio2(i,k,iSDeN)=Bio(i,k,iSDeN)
!!            Bio2(i,k,iLDeN)=Bio(i,k,iLDeN)
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2+                 &
     &                          Bio2(i,k,iSDeN)*tl_cff2
              tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)*cff4+                 &
     &                          Bio2(i,k,iLDeN)*tl_cff4

              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              tl_N_Flux_Remine=tl_Bio(i,k,iSDeN)*cff1+                  &
     &                         Bio(i,k,iSDeN)*tl_cff1+                  &
     &                         tl_Bio(i,k,iLDeN)*cff3+                  &
     &                         Bio(i,k,iLDeN)*tl_cff3

!!            Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!!            Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
!>            Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Remine
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-N_Flux_Remine*rOxNH4
              tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)+tl_N_Flux_Remine
              tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                      &
     &                          tl_N_Flux_Remine*rOxNH4
# ifdef PHOSPHORUS
              cff1=dtdays*SDeRRP(ng)*fac2
              tl_cff1=dtdays*(tl_SDeRRP*fac2+SDeRRP(ng)*tl_fac2)

              cff2=1.0_r8/(1.0_r8+cff1)
              tl_cff2=-cff2*cff2*tl_cff1

              cff3=dtdays*LDeRRP(ng)*fac2
              tl_cff3=dtdays*(tl_LDeRRP*fac2+LDeRRP(ng)*tl_fac2)

              cff4=1.0_r8/(1.0_r8+cff3)
              tl_cff4=-cff4*cff4*tl_cff3

!!            Bio2(i,k,iSDeP)=Bio(i,k,iSDeP)
!!            Bio2(i,k,iLDeP)=Bio(i,k,iLDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)*cff2+                 &
     &                          Bio2(i,k,iSDeP)*tl_cff2
              tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)*cff4+                 &
     &                          Bio2(i,k,iLDeP)*tl_cff4

              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              tl_P_Flux=tl_Bio(i,k,iSDeP)*cff1+Bio(i,k,iSDeP)*tl_cff1+  &
     &                  tl_Bio(i,k,iLDeP)*cff3+Bio(i,k,iLDeP)*tl_cff3

!!            Bio2(i,k,iPO4_)=Bio(i,k,iPO4_)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+tl_P_Flux
# endif
            END DO
          END DO
#else
! not yet (okada)
          cff1=dtdays*SDeRRN(ng)
          cff2=1.0_r8/(1.0_r8+cff1)
          cff3=dtdays*LDeRRN(ng)
          cff4=1.0_r8/(1.0_r8+cff3)
          DO k=1,N(ng)
            DO i=Istr,Iend
!!            Bio2(i,k,iSDeN)=Bio(i,k,iSDeN)
!!            Bio2(i,k,iLDeN)=Bio(i,k,iLDeN)
!>            Bio(i,k,iSDeN)=Bio(i,k,iSDeN)*cff2
!>            Bio(i,k,iLDeN)=Bio(i,k,iLDeN)*cff4
              tl_Bio(i,k,iSDeN)=tl_Bio(i,k,iSDeN)*cff2
              tl_Bio(i,k,iLDeN)=tl_Bio(i,k,iLDeN)*cff4
              N_Flux_Remine=Bio(i,k,iSDeN)*cff1+Bio(i,k,iLDeN)*cff3
              tl_N_Flux_Remine=tl_Bio(i,k,iSDeN)*cff1+                  &
     &                         tl_Bio(i,k,iLDeN)*cff3
!!            Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
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
!!            Bio2(i,k,iSDeP)=Bio(i,k,iSDeP)
!!            Bio2(i,k,iLDeP)=Bio(i,k,iLDeP)
!>            Bio(i,k,iSDeP)=Bio(i,k,iSDeP)*cff2
!>            Bio(i,k,iLDeP)=Bio(i,k,iLDeP)*cff4
              tl_Bio(i,k,iSDeP)=tl_Bio(i,k,iSDeP)*cff2
              tl_Bio(i,k,iLDeP)=tl_Bio(i,k,iLDeP)*cff4
              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              tl_P_Flux=tl_Bio(i,k,iSDeP)*cff1+tl_Bio(i,k,iLDeP)*cff3
!!            Bio2(i,k,iPO4_)=Bio(i,k,iPO4_)
!>            Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux
              tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)+tl_P_Flux
            END DO
          END DO
# endif
#endif
#if defined H2S && defined OXYGEN
!
!-----------------------------------------------------------------------
!  H2S Oxidation. okada
!-----------------------------------------------------------------------
! not yet (okada)
          DO k=1,N(ng)
            DO i=Istr,Iend
              fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8)
              fac2=fac1/(K_DO(ng)+fac1)
# ifdef TDEPENDANCE
              fac2=fac2*(1.05_r8**(Bio(i,k,iTemp)-20.0_r8))
# endif
              cff1=dtdays*H2SOR(ng)*fac2
              cff2=1.0_r8/(1.0_r8+cff1)
!>            Bio(i,k,iH2S_)=Bio(i,k,iH2S_)*cff2
              S_Flux=Bio(i,k,iH2S_)*cff1
!>            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-S_Flux*rOxH2S
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
! pass okada
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
            tl_SchmidtN_Ox=tl_Bio(i,k,itemp)*                           &
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
            O2_Flux=cff3*(O2satu-Bio3(i,k,iOxyg))
            tl_O2_Flux=tl_cff3*(O2satu-Bio3(i,k,iOxyg))+                 &
     &                 cff3*(tl_O2satu-tl_Bio(i,k,iOxyg))

!!          Bio3(i,k,iOxyg)=Bio(i,k,iOxyg)
!>          Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
!>   &                     O2_Flux*Hz_inv(i,k)
            tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)+                        &
     &                        tl_O2_Flux*Hz_inv(i,k)+                   &
     &                        O2_Flux*tl_Hz_inv(i,k)
          END DO
#endif