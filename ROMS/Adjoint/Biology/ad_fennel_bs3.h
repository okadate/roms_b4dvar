!
!  Compute appropriate basic state arrays III.
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
!>      ITER_LOOP: DO Iter=1,BioIter(ng)
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
                Att=(AttSW(ng)+                                         &
     &               AttChl(ng)*Bio(i,k,iChlo)+                         &
     &               AttFac)*                                           &
     &               (z_w(i,j,k)-z_w(i,j,k-1))
                ExpAtt=EXP(-Att)
                Itop=PAR
                PAR=Itop*(1.0_r8-ExpAtt)/Att    ! average at cell center
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
                cff1=Bio(i,k,iNH4_)*K_NH4(ng)
                cff2=Bio(i,k,iNO3_)*K_NO3(ng)
                inhNH4=1.0_r8/(1.0_r8+cff1)
                L_NH4=cff1/(1.0_r8+cff1)
                L_NO3=cff2*inhNH4/(1.0_r8+cff2)
                LTOT=L_NO3+L_NH4
#ifdef PHOSPHORUS
                cff3=Bio(i,k,iPO4_)*K_PO4(ng)
                L_PO4=cff3/(1.0_r8+cff3)
                LMIN=MIN(LTOT,L_PO4)
#endif
!
!  Nitrate and ammonium uptake by Phytoplankton.
!
                fac1=dtdays*t_PPmax*Bio(i,k,iPhyt)
                cff4=fac1*K_NO3(ng)*inhNH4/(1.0_r8+cff2)
                cff5=fac1*K_NH4(ng)/(1.0_r8+cff1)
#ifdef PHOSPHORUS
                cff6=fac1*PhyPN(ng)*K_PO4(ng)/(1.0_r8+cff3)
                IF (LMIN.eq.L_PO4) THEN
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)/(1.0_r8+cff6)
                  P_Flux=Bio(i,k,iPO4_)*cff6
                  N_Flux_NewProd=P_Flux/PhyPN(ng)*L_NO3/MAX(LTOT,eps)
                  N_Flux_RegProd=P_Flux/PhyPN(ng)*L_NH4/MAX(LTOT,eps)
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)-N_Flux_NewProd
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)-N_Flux_RegProd
                ELSE
                  Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                  Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                  N_Flux_NewProd=Bio(i,k,iNO3_)*cff4
                  N_Flux_RegProd=Bio(i,k,iNH4_)*cff5
                  P_Flux=(N_Flux_NewProd+N_Flux_RegProd)*PhyPN(ng)
                  Bio(i,k,iPO4_)=Bio(i,k,iPO4_)-P_Flux
                END IF
#else
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff4)
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff5)
                N_Flux_NewProd=Bio(i,k,iNO3_)*cff4
                N_Flux_RegProd=Bio(i,k,iNH4_)*cff5
#endif
                Bio(i,k,iPhyt)=Bio(i,k,iPhyt)+                          &
     &                         N_Flux_NewProd+N_Flux_RegProd
!
                Bio(i,k,iChlo)=Bio(i,k,iChlo)+                          &
#ifdef PHOSPHORUS
     &                         (dtdays*t_PPmax*t_PPmax*LMIN*LMIN*       &
#else
     &                         (dtdays*t_PPmax*t_PPmax*LTOT*LTOT*       &
#endif
     &                          Chl2C_m(ng)*Bio(i,k,iChlo))/            &
     &                         (PhyIS(ng)*MAX(Chl2C,eps)*PAR+eps)
#ifdef OXYGEN
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                          &
     &                         N_Flux_NewProd*rOxNO3+                   &
     &                         N_Flux_RegProd*rOxNH4
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
                fac2=MAX(Bio(i,k,iOxyg),0.0_r8)
                fac3=fac2/(K_Nitri(ng)+fac2)
# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                fac1=fac*fac3*cff
# else
                fac1=fac*fac3
# endif
#else
                fac1=fac
#endif
                cff1=(PAR-I_thNH4(ng))/                                 &
     &               (D_p5NH4(ng)+PAR-2.0_r8*I_thNH4(ng))
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
              fac1=dtdays*NitriR(ng)
              DO k=N(ng),1,-1
#if defined OXYGEN
                fac2=MAX(Bio(i,k,iOxyg),0.0_r8)
                fac3=fac2/(K_Nitri(ng)+fac2)
# ifdef TDEPENDANCE
                cff=NitriR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
                cff3=fac1*fac3*cff
# else
                cff3=fac1*fac3
# endif
#else
                cff3=fac1
#endif
                Bio(i,k,iNH4_)=Bio(i,k,iNH4_)/(1.0_r8+cff3)
                N_Flux_Nitrifi=Bio(i,k,iNH4_)*cff3
                Bio(i,k,iNO3_)=Bio(i,k,iNO3_)+N_Flux_Nitrifi
#ifdef OXYGEN
                Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-2.0_r8*N_Flux_Nitrifi
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
              cff2=1.0_r8/(1.0_r8+cff1)
# ifdef TDEPENDANCE
              fac2=DenitR_t(ng)**(Bio(i,k,itemp)-20.0_r8)
              cff3=cff2*fac1*fac2
# else
              cff3=cff2*fac1
# endif
              Bio(i,k,iNO3_)=Bio(i,k,iNO3_)/(1.0_r8+cff3)
            END DO
          END DO
#endif
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
              cff1=fac1*Bio(i,k,iZoop)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff3=1.0_r8/(1.0_r8+cff1)
              Bio(i,k,iPhyt)=cff3*Bio(i,k,iPhyt)
              Bio(i,k,iChlo)=cff3*Bio(i,k,iChlo)
!
! Phytoplankton assimilated to zooplankton and egested to small
! detritus.
!
              N_Flux_Assim=Bio(i,k,iPhyt)*cff1*ZooAE_N(ng)
              N_Flux_Egest=Bio(i,k,iPhyt)*cff1*(1.0_r8-ZooAE_N(ng))
              Bio(i,k,iZoop)=Bio(i,k,iZoop)+N_Flux_Assim
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Egest
!
! Phytoplankton mortality (limited by a phytoplankton minimum).
!
              N_Flux_Pmortal=cff2*MAX(Bio(i,k,iPhyt)-PhyMin(ng),0.0_r8)
              Bio(i,k,iPhyt)=Bio(i,k,iPhyt)-N_Flux_Pmortal
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Pmortal
              Bio(i,k,iChlo)=Bio(i,k,iChlo)-                            &
     &                       cff2*MAX(Bio(i,k,iChlo)-ChlMin(ng),0.0_r8)
#ifdef PHOSPHORUS
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+                            &
     &                       PhyPN(ng)*(N_Flux_Egest+N_Flux_Pmortal)+   &
     &                       (PhyPN(ng)-ZooPN(ng))*N_Flux_Assim
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
              fac1=fac3*Bio(i,k,iPhyt)*Bio(i,k,iPhyt)/                  &
     &             (K_Phy(ng)+Bio(i,k,iPhyt)*Bio(i,k,iPhyt))
              cff2=fac2*Bio(i,k,iZoop)
              cff3=fac1*ZooAE_N(ng)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)/(1.0_r8+cff2+cff3)
!
!  Zooplankton mortality and excretion.
!
              N_Flux_Zmortal=cff2*Bio(i,k,iZoop)
              N_Flux_Zexcret=cff3*Bio(i,k,iZoop)
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zexcret
              Bio(i,k,iSDeN)=Bio(i,k,iSDeN)+N_Flux_Zmortal
#ifdef PHOSPHORUS
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zexcret
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)+ZooPN(ng)*N_Flux_Zmortal
#endif
!
!  Zooplankton basal metabolism (limited by a zooplankton minimum).
!
              N_Flux_Zmetabo=cff1*MAX(Bio(i,k,iZoop)-ZooMin(ng),0.0_r8)
              Bio(i,k,iZoop)=Bio(i,k,iZoop)-N_Flux_Zmetabo
              Bio(i,k,iNH4_)=Bio(i,k,iNH4_)+N_Flux_Zmetabo
#ifdef PHOSPHORUS
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+ZooPN(ng)*N_Flux_Zmetabo
#endif
#ifdef OXYGEN
              Bio(i,k,iOxyg)=Bio(i,k,iOxyg)-                            &
     &                       rOxNH4*(N_Flux_Zmetabo+N_Flux_Zexcret)
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
              Bio(i,k,iLDeN)=Bio(i,k,iLDeN)+                            &
     &                       N_Flux_CoagP+N_Flux_CoagD
#ifdef PHOSPHORUS
              Bio(i,k,iSDeP)=Bio(i,k,iSDeP)-PhyPN(ng)*N_Flux_CoagD
              Bio(i,k,iLDeP)=Bio(i,k,iLDeP)+                            &
     &                       PhyPN(ng)*(N_Flux_CoagP+N_Flux_CoagD)
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
              fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8)
              fac2=fac1/(K_DO(ng)+fac1)
# ifdef TDEPENDANCE
              fac2=fac2*(RR_t(ng)**(Bio(i,k,iTemp)-20.0_r8))
# endif
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
              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux
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
              P_Flux=Bio(i,k,iSDeP)*cff1+Bio(i,k,iLDeP)*cff3
              Bio(i,k,iPO4_)=Bio(i,k,iPO4_)+P_Flux
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
              fac1=MAX(Bio(i,k,iOxyg)-6.0_r8,0.0_r8)
              fac2=fac1/(K_DO(ng)+fac1)
# ifdef TDEPENDANCE
              fac2=fac2*(1.05_r8**(Bio(i,k,iTemp)-20.0_r8))
# endif
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
# else
!
!  Calculate the Schmidt number for O2 in sea water (Wanninkhof, 1992).
!
            SchmidtN_Ox=1953.4_r8-                                      &
     &                  Bio(i,k,itemp)*(128.0_r8-                       &
     &                                  Bio(i,k,itemp)*                 &
     &                                  (3.9918_r8-                     &
     &                                   Bio(i,k,itemp)*0.050091_r8))
# endif

            cff3=cff2*u10squ*SQRT(660.0_r8/SchmidtN_Ox)
!
!  Calculate O2 saturation concentration using Garcia and Gordon
!  L&O (1992) formula, (EXP(AA) is in ml/l).
!
            TS=LOG((298.15_r8-Bio(i,k,itemp))/                          &
     &             (273.15_r8+Bio(i,k,itemp)))
            AA=OA0+TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5))))+          &
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
            Bio(i,k,iOxyg)=Bio(i,k,iOxyg)+                              &
     &                     O2_Flux*Hz_inv(i,k)
          END DO
#endif
!=======================================================================
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
!>        SINK_LOOP: DO isink=1,Nsink
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
              FC(i,N(ng))=0.0_r8            ! NO-flux boundary condition
#if defined LINEAR_CONTINUATION
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=2.0_r8*qc(i,N(ng))-bL(i,N(ng))
#elif defined NEUMANN
              bL(i,N(ng))=bR(i,N(ng)-1)
              bR(i,N(ng))=1.5_r8*qc(i,N(ng))-0.5_r8*bL(i,N(ng))
#else
              bR(i,N(ng))=qc(i,N(ng))       ! default strictly monotonic
              bL(i,N(ng))=qc(i,N(ng))       ! conditions
              bR(i,N(ng)-1)=qc(i,N(ng))
#endif
#if defined LINEAR_CONTINUATION
              bR(i,1)=bL(i,2)
              bL(i,1)=2.0_r8*qc(i,1)-bR(i,1)
#elif defined NEUMANN
              bR(i,1)=bL(i,2)
              bL(i,1)=1.5_r8*qc(i,1)-0.5_r8*bR(i,1)
#else
              bL(i,2)=qc(i,1)               ! bottom grid boxes are
              bR(i,1)=qc(i,1)               ! re-assumed to be
              bL(i,1)=qc(i,1)               ! piecewise constant.
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
                FC(i,k-1)=FC(i,k-1)+                                    &
     &                    Hz(i,j,ks)*cu*                                &
     &                    (bL(i,ks)+                                    &
     &                     cu*(0.5_r8*(bR(i,ks)-bL(i,ks))-              &
     &                         (1.5_r8-cu)*                             &
     &                         (bR(i,ks)+bL(i,ks)-                      &
     &                          2.0_r8*qc(i,ks))))
              END DO
            END DO
            DO k=1,N(ng)
              DO i=Istr,Iend
                Bio(i,k,ibio)=qc(i,k)+(FC(i,k)-FC(i,k-1))*Hz_inv(i,k)
              END DO
            END DO
#if defined BIO_SEDIMENT
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
          END DO

#if defined BIO_SED_CONSTANT
!
!  Elution and oxygen consumption parameters (okada)
!
          cff1=R_SODf(ng)/mol2g_O2    !SOD flux
          cff2=R_NH4f(ng)/14.0_r8     !NH4 elution flux from sediment
          cff3=R_PO4f(ng)/31.0_r8     !PO4 elution flux from sediment
!
!-----------------------------------------------------------------------
!  Elution and oxygen consumption from/by sediment. (Okada, 2014/02/13)
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
# ifdef TDEPENDANCE
            fac1=dtdays*(1.05_r8**(Bio(i,1,itemp)-20.0_r8))
# else
            fac1=dtdays
# endif
            cff=fac1*Hz_inv(i,1)
            Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff*cff2
# ifdef PHOSPHORUS
            Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff*cff3
# endif
# ifdef OXYGEN
            cff4=MAX(MIN(Bio(i,1,iOxyg),cff*cff1),0.0_r8)
            Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff4
#  ifdef H2S
            cff5=MIN(MIN(Bio(i,1,iOxyg),cff*cff1),0.0_r8)
            Bio(i,1,iH2S_)=Bio(i,1,iH2S_)+cff5*rOxH2S
#  endif
# endif
          END DO
#endif
          END IF
!=======================================================================
!>      END DO ITER_LOOP
        END DO
!
!  End of compute basic state arrays III.
!