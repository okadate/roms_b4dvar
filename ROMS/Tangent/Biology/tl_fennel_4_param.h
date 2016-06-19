!
!-----------------------------------------------------------------------
!  Vertical sinking terms.
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
                tl_cff1=tl_FC(i,0)*Hz_inv(i,1)+FC(i,0)*tl_Hz_inv(i,1)
# ifdef DENITRIFICATION
!>              Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1*cff2
                tl_Bio(i,1,iNH4_)=tl_Bio(i,1,iNH4_)+tl_cff1*cff2
#  ifdef OXYGEN
!>              Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff3
                tl_Bio(i,1,iOxyg)=tl_Bio(i,1,iOxyg)-tl_cff1*cff3
#  endif
# else
!>              Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff1
                tl_Bio(i,1,iNH4_)=tl_Bio(i,1,iNH4_)+tl_cff1
#  ifdef OXYGEN
!>              Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff1*cff4
                tl_Bio(i,1,iOxyg)=tl_Bio(i,1,iOxyg)-tl_cff1*cff4
#  endif
# endif
              END DO
            END IF
# ifdef PHOSPHORUS
            IF ((ibio.eq.iLDeP).or.                                     &
     &          (ibio.eq.iSDeP)) THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                tl_cff1=tl_FC(i,0)*Hz_inv(i,1)+FC(i,0)*tl_Hz_inv(i,1)

!>              Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1
                tl_Bio(i,1,iPO4_)=tl_Bio(i,1,iPO4_)+tl_cff1
              END DO
            END IF
            IF (ibio.eq.iPhyt)THEN
              DO i=Istr,Iend
                cff1=FC(i,0)*Hz_inv(i,1)
                tl_cff1=tl_FC(i,0)*Hz_inv(i,1)+FC(i,0)*tl_Hz_inv(i,1)

!>              Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff1*PhyPN(ng)
                tl_Bio(i,1,iPO4_)=tl_Bio(i,1,iPO4_)+tl_cff1*PhyPN(ng)
              END DO
            END IF
# endif
# if defined H2S && defined OXYGEN
            IF ((ibio.eq.iPhyt).or.                                     &
     &          (ibio.eq.iSDeN).or.                                     &
     &          (ibio.eq.iLDeN)) THEN
              DO i=Istr,Iend
!>              Bio(i,1,iH2S_)=Bio(i,1,iH2S_)-                          &
!>   &                         0.5_r8*MIN(Bio(i,1,iOxyg),0.0_r8)
                tl_Bio(i,1,iH2S_)=tl_Bio(i,1,iH2S_)-                    &
     &                            (0.5_r8+SIGN(0.5_r8,Bio(i,1,iOxyg)))* &
     &                            tl_Bio(i,1,iOxyg)
!>              Bio(i,1,iOxyg)=MAX(Bio(i,1,iOxyg),0.0_r8)
                tl_Bio(i,1,iOxyg)=(0.5_r8+SIGN(0.5_r8,Bio(i,1,iOxyg)))* &
     &                            tl_Bio(i,1,iOxyg)
              END DO
            END IF
# endif
#endif
          END DO SINK_LOOP
#if defined BIO_SED_CONSTANT
!
!  Elution and oxygen consumption parameters (okada)
!
# ifdef OXYGEN
          cff1=R_SODf(ng)/mol2g_O2    !SOD flux
          tl_cff1=tl_R_SODf/mol2g_O2
# endif
          cff2=R_NH4f(ng)/14.0_r8     !NH4 elution flux from sediment
          tl_cff2=tl_R_NH4f/14.0_r8
# ifdef PHOSPHORUS
          cff3=R_PO4f(ng)/31.0_r8     !PO4 elution flux from sediment
          tl_cff3=tl_R_PO4f/31.0_r8
# endif
!
!-----------------------------------------------------------------------
!  Elution and oxygen consumption from/by sediment. (Okada, 2014/02/13)
!-----------------------------------------------------------------------
!
          DO i=Istr,Iend
# ifdef TDEPENDANCE
            fac2=1.05_r8**(Bio(i,1,itemp)-20.0_r8)
            tl_fac2=fac2*tl_Bio(i,k,itemp)*LOG(1.05_r8)

            fac1=dtdays*fac2
            tl_fac1=dtdays*tl_fac2
# else
            fac1=dtdays
            tl_fac1=0.0_r8
# endif
            cff=fac1*Hz_inv(i,1)
            tl_cff=tl_fac1*Hz_inv(i,1)+fac1*tl_Hz_inv(i,1)

!>          Bio(i,1,iNH4_)=Bio(i,1,iNH4_)+cff*cff2
            tl_Bio(i,1,iNH4_)=tl_Bio(i,1,iNH4_)+tl_cff*cff2+cff*tl_cff2
# ifdef PHOSPHORUS
!>          Bio(i,1,iPO4_)=Bio(i,1,iPO4_)+cff*cff3
            tl_Bio(i,1,iPO4_)=tl_Bio(i,1,iPO4_)+tl_cff*cff3+cff*tl_cff3
# endif
# ifdef OXYGEN
            fac3=MAX(Bio(i,1,iOxyg),0.0_r8)
            tl_fac3=(0.5_r8+SIGN(0.5_r8,Bio(i,1,iOxyg)))*               &
     &              tl_Bio(i,1,iOxyg)

            cff4=MIN(fac3,cff*cff1)
            tlfac=SIGN(0.5_r8,cff*cff1-fac3)
            tl_cff4=(0.5_r8+tlfac)*tl_fac3+                             &
     &              (0.5_r8-tlfac)*(tl_cff*cff1+cff*tl_cff1)

!>          Bio(i,1,iOxyg)=Bio(i,1,iOxyg)-cff4
            tl_Bio(i,1,iOxyg)=tl_Bio(i,1,iOxyg)-tl_cff4
# endif
          END DO
#endif