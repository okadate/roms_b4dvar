#include "tl_fennel_bs1.h"
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
            IF (PARsur(i).gt.0.0_r8) THEN
              DO k=N(ng),1,-1
!
!  Compute average light attenuation for each grid cell. To include
!  other attenuation contributions like suspended sediment or CDOM
!  modify AttFac.
!
                Att=(AttSW(ng)+                                         &
     &               AttChl(ng)*Bio1(i,k,iChlo)+                        &
     &               AttFac)*                                           &
     &               (z_w(i,j,k)-z_w(i,j,k-1))
                tl_Att=AttChl(ng)*tl_Bio(i,k,iChlo)*                    &
     &                 (z_w(i,j,k)-z_w(i,j,k-1))+                       &
     &                 (AttSW(ng)+AttChl(ng)*Bio1(i,k,iChlo)+AttFac)*   &
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
                cff1=Bio1(i,k,iChlo)/(Bio1(i,k,iPhyt)*cff+eps)
                tl_cff1=(tl_Bio(i,k,iChlo)-cff1*cff*tl_Bio(i,k,iPhyt))/ &
     &                  (Bio1(i,k,iPhyt)*cff+eps)

                Chl2C=MIN(cff1,Chl2C_m(ng))
                tl_Chl2C=(0.5_r8+SIGN(0.5_r8,Chl2C_m(ng)-cff1))*tl_cff1
!
!  Temperature-limited and light-limited growth rate (Eppley, R.W.,
!  1972, Fishery Bulletin, 70: 1063-1085; here 0.59=ln(2)*0.851).
!  Check value for Vp is 2.9124317 at 19.25 degC.
!
                Vp=Vp0(ng)*0.59_r8*(1.066_r8**Bio(i,k,itemp))
                tl_Vp=Vp*tl_Bio(i,k,itemp)

                fac1=PAR*PhyIS(ng)
                tl_fac1=tl_PAR*PhyIS(ng)

                Epp=Vp/SQRT(Vp*Vp+fac1*fac1)
                tl_Epp=(tl_Vp-Epp*(Vp*tl_Vp+fac1*tl_fac1))/             &
     &                 SQRT(Vp*Vp+fac1*fac1)

                t_PPmax=Epp*fac1
                tl_t_PPmax=tl_Epp*fac1+Epp*tl_fac1
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
                tl_L_NH4=(1.0_r8-L_NH4)/(1.0_r8+cff1)*tl_cff1
                tl_L_NO3=(tl_cff2*inhNH4+cff2*tl_inhNH4-L_NO3*tl_cff2)/ &
     &                   (1.0_r8+cff2)

                LTOT=L_NO3+L_NH4
                tl_LTOT=tl_L_NO3+tl_L_NH4
#ifdef PHOSPHORUS
                cff3=Bio1(i,k,iPO4_)*K_PO4(ng)
                tl_cff3=tl_Bio(i,k,iPO4_)*K_PO4(ng)

                L_PO4=cff3/(1.0_r8+cff3)
                tl_L_PO4=(1.0_r8-L_PO4)/(1.0_r8+cff3)*tl_cff3

                LMIN=MIN(LTOT,L_PO4)
                tl_LMIN=(0.5_r8+SIGN(0.5_r8,L_PO4-LTOT))*tl_LTOT+       &
     &                  (0.5_r8-SIGN(0.5_r8,L_PO4-LTOT))*tl_L_PO4
#endif
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
                fac1=dtdays*t_PPmax*Bio1(i,k,iPhyt)
                tl_fac1=dtdays*(tl_t_PPmax*Bio1(i,k,iPhyt)+             &
     &                          t_PPmax*tl_Bio(i,k,iPhyt))

                cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)
                tl_cff4=(K_NO3(ng)*(tl_fac1*inhNH4+fac1*tl_inhNH4)-     &
     &                   cff4*tl_cff2)/(1.0_r8+cff2)

                cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)
                tl_cff5=(tl_fac1*K_NH4(ng)-cff5*tl_cff1)/(1.0_r8+cff1)

!>              Bio1(i,k,iNO3_)=Bio(i,k,iNO3_)
!>              Bio1(i,k,iNH4_)=Bio(i,k,iNH4_)
#ifdef PHOSPHORUS
                cff6=fac1*PhyPN(ng)*K_PO4(ng)/(1.0_r8+cff3)
                tl_cff6=(tl_fac1*PhyPN(ng)*K_PO4(ng)-cff6*tl_cff3)/     &
     &                  (1.0_r8+cff3)

!>              Bio1(i,k,iPO4_)=Bio(i,k,iPO4_)

                IF (LMIN.eq.L_PO4) THEN
!>                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff6)
                  tl_Bio(i,k,iPO4_)=(tl_Bio(i,k,iPO4_)-Bio(i,k,iPO4_)*  &
     &                               tl_cff6)/(1.0_r8+cff6)

                  P_Flux=Bio(i,k,iPO4_)*cff6
                  tl_P_Flux=tl_Bio(i,k,iPO4_)*cff6+                     &
     &                      Bio(i,k,iPO4_)*tl_cff6

                  fac2=L_NO3/MAX(LTOT,eps)
                  tl_fac2=(tl_L_NO3-fac2*(0.5_r8+SIGN(0.5_r8,LTOT-eps))*&
     &                              tl_LTOT)/L_NO3

                  N_Flux_NewProd=P_Flux/PhyPN(ng)*fac2
                  tl_N_Flux_NewProd=(tl_P_Flux*fac2+P_Flux*tl_fac2)/    &
     &                              PhyPN(ng)

                  fac3=L_NH4/MAX(LTOT,eps)
                  tl_fac3=(tl_L_NH4-fac3*(0.5_r8+SIGN(0.5_r8,LTOT-eps))*&
     &                              tl_LTOT)/L_NH4

                  N_Flux_RegProd=P_Flux/PhyPN(ng)*fac3
                  tl_N_Flux_RegProd=(tl_P_Flux*fac3+P_Flux*tl_fac3)/    &
     &                              PhyPN(ng)

!>                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd
                  tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)-tl_N_Flux_NewProd

!>                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd
                  tl_Bio(i,k,iNH4_)=tl_Bio(i,k,iNH4_)-tl_N_Flux_RegProd
                ELSE
!>                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                  tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio2(i,k,iNO3_)* &
     &                               tl_cff4)/(1.0_r8+cff4)

!>                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio2(i,k,iNH4_)* &
     &                               tl_cff5)/(1.0_r8+cff5)

                  N_Flux_NewProd=Bio2(i,k,iNO3_)*cff4
                  tl_N_Flux_NewProd=tl_Bio(i,k,iNO3_)*cff4+             &
     &                              Bio2(i,k,iNO3_)*tl_cff4

                  N_Flux_RegProd=Bio2(i,k,iNH4_)*cff5
                  tl_N_Flux_RegProd=tl_Bio(i,k,iNH4_)*cff5+             &
     &                              Bio2(i,k,iNH4_)*tl_cff5

                  P_Flux=(N_Flux_NewProd+N_Flux_RegProd)*PhyPN(ng)
                  tl_P_Flux=(tl_N_Flux_NewProd+tl_N_Flux_RegProd)*      &
     &                      PhyPN(ng)

!>                Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-P_Flux
                  tl_Bio(i,k,iPO4_)=tl_Bio(i,k,iPO4_)-tl_P_Flux
                END IF
#else
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio2(i,k,iNO3_)*   &
     &                             tl_cff4)/(1.0_r8+cff4)

!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio2(i,k,iNH4_)*   &
     &                             tl_cff5)/(1.0_r8+cff5)

                N_Flux_NewProd=Bio2(i,k,iNO3_)*cff4
                tl_N_Flux_NewProd=tl_Bio(i,k,iNO3_)*cff4+               &
     &                            Bio2(i,k,iNO3_)*tl_cff4

                N_Flux_RegProd=Bio2(i,k,iNH4_)*cff5
                tl_N_Flux_RegProd=tl_Bio(i,k,iNH4_)*cff5+               &
     &                            Bio2(i,k,iNH4_)*tl_cff5
#endif
!>              Bio1(i,k,iPhyt)=Bio(i,k,iPhyt)
!>              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
!>   &                         N_Flux_NewProd+N_Flux_RegProd
                tl_Bio(i,k,iPhyt)=tl_Bio(i,k,iPhyt)+                    &
     &                            tl_N_Flux_NewProd+tl_N_Flux_RegProd
!
#ifdef PHOSPHORUS
                fac1=dtdays*t_PPmax*t_PPmax*LMIN*LMIN*Chl2C_m(ng)
                tl_fac1=dtdays*2.0_r8*(t_PPmax*tl_t_PPmax*LMIN*LMIN+    &
     &                                 t_PPmax*t_PPmax*LMIN*tl_LMIN)*   &
     &                                Chl2C_m(ng)
#else
                fac1=dtdays*t_PPmax*t_PPmax*LTOT*LTOT*Chl2C_m(ng)
                tl_fac1=dtdays*2.0_r8*(t_PPmax*tl_t_PPmax*LTOT*LTOT+    &
     &                                 t_PPmax*t_PPmax*LTOT*tl_LTOT)*   &
     &                                Chl2C_m(ng)
#endif
                fac2=PhyIS(ng)*MAX(Chl2C,eps)*PAR
                tl_fac2=PhyIS(ng)*((0.5_r8+SIGN(0.5_r8,Chl2C-eps))*     &
     &                             tl_Chl2C*PAR+MAX(Chl2C,eps)*tl_PAR)

                fac3=fac1/(fac2+eps)
                tl_fac3=(tl_fac1-tl_fac2*fac3)/(fac2+eps)

!>              Bio1(i,k,iChlo)=Bio(i,k,iChlo)
!>              Bio(i,k,iChlo)=Bio(i,k,iChlo)+Bio(i,k,iChlo)*fac3
                tl_Bio(i,k,iChlo)=tl_Bio(i,k,iChlo)+                    &
     &                            tl_Bio(i,k,iChlo)*fac1+               &
     &                            Bio(i,k,iChlo)*tl_fac1
#ifdef OXYGEN
!>              Bio1(i,k,iOxyg)=Bio(i,k,iOxyg)
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
                fac=dtdays*NitriR(ng)
#ifdef OXYGEN
                fac2=MAX(Bio2(i,k,iOxyg),0.0_r8)
                tl_fac2=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)))*          &
     &                  tl_Bio(i,k,iOxyg)

                fac3=fac2/(K_Nitri(ng)+fac2)
                tl_fac3=-fac3/(K_Nitri(ng)+fac2)*tl_fac2
# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                tl_cff=cff*tl_Bio(i,k,itemp)*LOG(NitriR_t(ng))

                fac1=fac*fac3*cff
                tl_fac1=fac*tl_fac3*cff+                                &
     &                  fac*fac3*tl_cff
# else
                fac1=fac*fac3
                tl_fac1=fac*tl_fac3
# endif
#else
                fac1=fac
                tl_fac1=0.0_r8
#endif
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
                tl_cff1=(tl_PAR-cff1*tl_PAR)/                           &
     &                  (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))

                cff2=1.0_r8-MAX(0.0_r8,cff1)
                tl_cff2=-(0.5_r8-SIGN(0.5_r8,-cff1))*tl_cff1

                cff3=fac1*cff2
                tl_cff3=tl_fac1*cff2+fac1*tl_cff2

!>              Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio(i,k,iNH4_)*    &
     &                             tl_cff3)/(1.0_r8+cff3)

                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                tl_N_Flux_Nitrifi=tl_Bio(i,k,iNH4_)*cff3+               &
     &                            Bio(i,k,iNH4_)*tl_cff3

!>              Bio2(i,k,iNO3_)=Bio(i,k,iNO3_)
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
                tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+tl_N_Flux_Nitrifi
#ifdef OXYGEN
!>              Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
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
            ELSE
              fac1=dtdays*NitriR(ng)
              DO k=N(ng),1,-1
#if defined OXYGEN
                fac2=MAX(Bio2(i,k,iOxyg),0.0_r8)
                tl_fac2=(0.5_r8+SIGN(0.5_r8,Bio2(i,k,iOxyg)))*          &
     &                  tl_Bio(i,k,iOxyg)

                fac3=fac2/(K_Nitri(ng)+fac2)
                tl_fac3=-fac3/(K_Nitri(ng)+fac2)*tl_fac2

# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                tl_cff=cff*tl_Bio(i,k,itemp)*LOG(NitriR_t(ng))

                cff3=fac1*fac3*cff
                tl_cff3=fac1*tl_fac3*cff+                               &
     &                  fac1*fac3*tl_cff
# else
                cff3=fac1*fac3
                tl_cff3=fac1*tl_fac3
# endif
#else
                cff3=fac1
                tl_cff3=0.0_r8
#endif
!>              Bio2(i,k,iNH4_)=Bio(i,k,iNH4_)
!>              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                tl_Bio(i,k,iNH4_)=(tl_Bio(i,k,iNH4_)-Bio(i,k,iNH4_)*    &
     &                             tl_cff3)/(1.0_r8+cff3)

                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                tl_N_Flux_Nitrifi=tl_Bio(i,k,iNH4_)*cff3+               &
     &                            Bio(i,k,iNH4_)*tl_cff3

!>              Bio2(i,k,iNO3_)=Bio(i,k,iNO3_)
!>              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
                tl_Bio(i,k,iNO3_)=tl_Bio(i,k,iNO3_)+tl_N_Flux_Nitrifi
#ifdef OXYGEN
!>              Bio2(i,k,iOxyg)=Bio(i,k,iOxyg)
!>              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
                tl_Bio(i,k,iOxyg)=tl_Bio(i,k,iOxyg)-                    &
     &                            2.0_r8*tl_N_Flux_Nitrifi
#endif
              END DO
            END IF
          END DO
#if defined OXYGEN && defined DENITRIFICATION
!
!-----------------------------------------------------------------------
! Denitrification in anoxic water                       Okada 2014/02/13
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
            DO k=1,N(ng)
              fac1=dtdays*DenitR(ng)
              cff1=MAX(Bio(i,k,iOxyg),0.0_r8)/K_Denit(ng)
              tl_cff1=(0.5_r8+SIGN(0.5_r8,Bio(i,k,iOxyg)))*             &
     &                 tl_Bio(i,k,iOxyg)/K_Denit(ng)

              cff2=1.0_r8/(1.0_r8+cff1)
              tl_cff2=-cff2*cff2*tl_cff1
# ifdef TDEPENDANCE
              fac2=DenitR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              tl_fac2=fac2*tl_Bio(i,k,itemp)*LOG(DenitR_t(ng))

              cff3=cff2*fac1*fac2
              tl_cff3=tl_cff2*fac1*fac2+                                &
     &                cff2*tl_fac1*fac2+                                &
     &                cff2*fac1*tl_fac2
# else
              cff3=cff2*fac1
              tl_cff3=tl_cff2*fac1+cff2*tl_fac1
# endif
!>            Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff3)
              tl_Bio(i,k,iNO3_)=(tl_Bio(i,k,iNO3_)-Bio(i,k,iNO3_)*      &
     &                           tl_cff3)/(1.0_r8+cff3)
            END DO
          END DO
#endif