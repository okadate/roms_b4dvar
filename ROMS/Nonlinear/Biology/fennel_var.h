/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2016 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the Fennel et al. (2006) ecosystem   **
**  model variables that are used in input and output NetCDF files.   **
**  The metadata information is read from "varinfo.dat".              **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iLDeN)')
                idTvar(iLDeN)=varid
              CASE ('idTvar(iSDeN)')
                idTvar(iSDeN)=varid
              CASE ('idTvar(iChlo)')
                idTvar(iChlo)=varid
# ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
              CASE ('idTvar(iLDeC)')
                idTvar(iLDeC)=varid
              CASE ('idTvar(iSDeC)')
                idTvar(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
# endif

# ifdef PHOSPHORUS
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
              CASE ('idTvar(iLDeP)')
                idTvar(iLDeP)=varid
              CASE ('idTvar(iSDeP)')
                idTvar(iSDeP)=varid
# endif
# ifdef H2S
              CASE ('idTvar(iH2S_)')
                idTvar(iH2S_)=varid
# endif
# ifdef ADJUST_PARAM
              CASE ('idPvar(iAttSW)')
                idPvar(iAttSW)=varid
              CASE ('idPvar(iAttChl)')
                idPvar(iAttChl)=varid
              CASE ('idPvar(iVp0)')
                idPvar(iVp0)=varid
              CASE ('idPvar(iI_thNH4)')
                idPvar(iI_thNH4)=varid
              CASE ('idPvar(iD_p5NH4)')
                idPvar(iD_p5NH4)=varid
              CASE ('idPvar(iNitriR)')
                idPvar(iNitriR)=varid
              CASE ('idPvar(iK_NO3)')
                idPvar(iK_NO3)=varid
              CASE ('idPvar(iK_NH4)')
                idPvar(iK_NH4)=varid
              CASE ('idPvar(iK_Phy)')
                idPvar(iK_Phy)=varid
              CASE ('idPvar(iChl2C_m)')
                idPvar(iChl2C_m)=varid
              CASE ('idPvar(iPhyCN)')
                idPvar(iPhyCN)=varid
              CASE ('idPvar(iPhyIP)')
                idPvar(iPhyIP)=varid
              CASE ('idPvar(iPhyIS)')
                idPvar(iPhyIS)=varid
              CASE ('idPvar(iPhyMR)')
                idPvar(iPhyMR)=varid
              CASE ('idPvar(iZooAE_N)')
                idPvar(iZooAE_N)=varid
              CASE ('idPvar(iZooBM)')
                idPvar(iZooBM)=varid
              CASE ('idPvar(iZooCN)')
                idPvar(iZooCN)=varid
              CASE ('idPvar(iZooER)')
                idPvar(iZooER)=varid
              CASE ('idPvar(iZooGR)')
                idPvar(iZooGR)=varid
              CASE ('idPvar(iZooMR)')
                idPvar(iZooMR)=varid
              CASE ('idPvar(iLDeRRN)')
                idPvar(iLDeRRN)=varid
              CASE ('idPvar(iCoagR)')
                idPvar(iCoagR)=varid
              CASE ('idPvar(iSDeRRN)')
                idPvar(iSDeRRN)=varid
              CASE ('idPvar(iwPhy)')
                idPvar(iwPhy)=varid
              CASE ('idPvar(iwLDet)')
                idPvar(iwLDet)=varid
              CASE ('idPvar(iwSDet)')
                idPvar(iwSDet)=varid

              CASE ('idPvar(iK_Nitri)')
                idPvar(iK_Nitri)=varid
              CASE ('idPvar(iK_Denit)')
                idPvar(iK_Denit)=varid
              CASE ('idPvar(iDenitR)')
                idPvar(iDenitR)=varid
              CASE ('idPvar(iK_PO4)')
                idPvar(iK_PO4)=varid
              CASE ('idPvar(iPhyPN)')
                idPvar(iPhyPN)=varid
              CASE ('idPvar(iZooPN)')
                idPvar(iZooPN)=varid
              CASE ('idPvar(iK_DO)')
                idPvar(iK_DO)=varid
              CASE ('idPvar(iLDeRRP)')
                idPvar(iLDeRRP)=varid
              CASE ('idPvar(iSDeRRP)')
                idPvar(iSDeRRP)=varid
              CASE ('idPvar(iR_SODf)')
                idPvar(iR_SODf)=varid
              CASE ('idPvar(iR_NH4f)')
                idPvar(iR_NH4f)=varid
              CASE ('idPvar(iR_PO4f)')
                idPvar(iR_PO4f)=varid
              CASE ('idPvar(iR_NH4f_m)')
                idPvar(iR_NH4f_m)=varid
              CASE ('idPvar(iR_PO4f_m)')
                idPvar(iR_PO4f_m)=varid
              CASE ('idPvar(iK_DO_npf)')
                idPvar(iK_DO_npf)=varid
              CASE ('idPvar(it_SODf)')
                idPvar(it_SODf)=varid
# endif
# ifdef DIAGENESIS
              CASE ('idBpw(iwO2_)')
                idBpw(iwO2_)=varid
              CASE ('idBpw(iwNH4)')
                idBpw(iwNH4)=varid
              CASE ('idBpw(iwNO3)')
                idBpw(iwNO3)=varid
              CASE ('idBpw(iwPO4)')
                idBpw(iwPO4)=varid
              CASE ('idBpw(iwSO4)')
                idBpw(iwSO4)=varid
              CASE ('idBpw(iwH2S)')
                idBpw(iwH2S)=varid
              CASE ('idBpw(iwMn_)')
                idBpw(iwMn_)=varid
              CASE ('idBpw(iwFe_)')
                idBpw(iwFe_)=varid
              CASE ('idBpw(iwCH4)')
                idBpw(iwCH4)=varid
              CASE ('idBpw(iwDOMf)')
                idBpw(iwDOMf)=varid
              CASE ('idBpw(iwDOMs)')
                idBpw(iwDOMs)=varid
                
              CASE ('idBsm(iPOMf)')
                idBsm(iPOMf)=varid
              CASE ('idBsm(iPOMs)')
                idBsm(iPOMs)=varid
              CASE ('idBsm(iPOMn)')
                idBsm(iPOMn)=varid
              CASE ('idBsm(iFeOA)')
                idBsm(iFeOA)=varid
              CASE ('idBsm(iFeOB)')
                idBsm(iFeOB)=varid
              CASE ('idBsm(iFeOP)')
                idBsm(iFeOP)=varid
              CASE ('idBsm(iMnOA)')
                idBsm(iMnOA)=varid
              CASE ('idBsm(iMnOB)')
                idBsm(iMnOB)=varid
              CASE ('idBsm(iS0__)')
                idBsm(iS0__)=varid
              CASE ('idBsm(iFeS_)')
                idBsm(iFeS_)=varid
              CASE ('idBsm(iFeS2)')
                idBsm(iFes2)=varid
                
              CASE ('idFpw(iwO2_)')
                idFpw(iwO2_)=varid
              CASE ('idFpw(iwNH4)')
                idFpw(iwNH4)=varid
              CASE ('idFpw(iwNO3)')
                idFpw(iwNO3)=varid
              CASE ('idFpw(iwPO4)')
                idFpw(iwPO4)=varid
              CASE ('idFpw(iwSO4)')
                idFpw(iwSO4)=varid
              CASE ('idFpw(iwH2S)')
                idFpw(iwH2S)=varid
              CASE ('idFpw(iwMn_)')
                idFpw(iwMn_)=varid
              CASE ('idFpw(iwFe_)')
                idFpw(iwFe_)=varid
              CASE ('idFpw(iwCH4)')
                idFpw(iwCH4)=varid
              CASE ('idFpw(iwDOMf)')
                idFpw(iwDOMf)=varid
              CASE ('idFpw(iwDOMs)')
                idFpw(iwDOMs)=varid
                
              CASE ('idFsm(iPOMf)')
                idFsm(iPOMf)=varid
              CASE ('idFsm(iPOMs)')
                idFsm(iPOMs)=varid
              CASE ('idFsm(iPOMn)')
                idFsm(iPOMn)=varid
              CASE ('idFsm(iFeOA)')
                idFsm(iFeOA)=varid
              CASE ('idFsm(iFeOB)')
                idFsm(iFeOB)=varid
              CASE ('idFsm(iFeOP)')
                idFsm(iFeOP)=varid
              CASE ('idFsm(iMnOA)')
                idFsm(iMnOA)=varid
              CASE ('idFsm(iMnOB)')
                idFsm(iMnOB)=varid
              CASE ('idFsm(iS0__)')
                idFsm(iS0__)=varid
              CASE ('idFsm(iFeS_)')
                idFsm(iFeS_)=varid
              CASE ('idFsm(iFeS2)')
                idFsm(iFes2)=varid
# endif

/*
**  Adjoint sensitivity state biological tracers.
*/

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
              CASE ('idTads(iNO3_)')
                idTads(iNO3_)=varid
              CASE ('idTads(iNH4_)')
                idTads(iNH4_)=varid
              CASE ('idTads(iPhyt)')
                idTads(iPhyt)=varid
              CASE ('idTads(iZoop)')
                idTads(iZoop)=varid
              CASE ('idTads(iLDeN)')
                idTads(iLDeN)=varid
              CASE ('idTads(iSDeN)')
                idTads(iSDeN)=varid
              CASE ('idTads(iChlo)')
                idTads(iChlo)=varid
# ifdef CARBON
              CASE ('idTads(iTIC_)')
                idTads(iTIC_)=varid
              CASE ('idTads(iTAlk)')
                idTads(iTAlk)=varid
              CASE ('idTads(iLDeC)')
                idTads(iLDeC)=varid
              CASE ('idTads(iSDeC)')
                idTads(iSDeC)=varid
# endif
# ifdef OXYGEN
              CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid
# endif

# ifdef PHOSPHORUS
              CASE ('idTads(iPO4_)')
                idTads(iPO4_)=varid
              CASE ('idTads(iLDeP)')
                idTads(iLDeP)=varid
              CASE ('idTads(iSDeP)')
                idTads(iSDeP)=varid
# endif
# ifdef H2S
              CASE ('idTads(iH2S_)')
                idTads(iH2S_)=varid
# endif
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid

              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid

              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid

              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid

              CASE ('idTbry(iwest,iSDeN)')
                idTbry(iwest,iSDeN)=varid
              CASE ('idTbry(ieast,iSDeN)')
                idTbry(ieast,iSDeN)=varid
              CASE ('idTbry(isouth,iSDeN)')
                idTbry(isouth,iSDeN)=varid
              CASE ('idTbry(inorth,iSDeN)')
                idTbry(inorth,iSDeN)=varid

              CASE ('idTbry(iwest,iLDeN)')
                idTbry(iwest,iLDeN)=varid
              CASE ('idTbry(ieast,iLDeN)')
                idTbry(ieast,iLDeN)=varid
              CASE ('idTbry(isouth,iLDeN)')
                idTbry(isouth,iLDeN)=varid
              CASE ('idTbry(inorth,iLDeN)')
                idTbry(inorth,iLDeN)=varid

              CASE ('idTbry(iwest,iChlo)')
                idTbry(iwest,iChlo)=varid
              CASE ('idTbry(ieast,iChlo)')
                idTbry(ieast,iChlo)=varid
              CASE ('idTbry(isouth,iChlo)')
                idTbry(isouth,iChlo)=varid
              CASE ('idTbry(inorth,iChlo)')
                idTbry(inorth,iChlo)=varid

#ifdef CARBON
              CASE ('idTbry(iwest,iSDeC)')
                idTbry(iwest,iSDeC)=varid
              CASE ('idTbry(ieast,iSDeC)')
                idTbry(ieast,iSDeC)=varid
              CASE ('idTbry(isouth,iSDeC)')
                idTbry(isouth,iSDeC)=varid
              CASE ('idTbry(inorth,iSDeC)')
                idTbry(inorth,iSDeC)=varid

              CASE ('idTbry(iwest,iLDeC)')
                idTbry(iwest,iLDeC)=varid
              CASE ('idTbry(ieast,iLDeC)')
                idTbry(ieast,iLDeC)=varid
              CASE ('idTbry(isouth,iLDeC)')
                idTbry(isouth,iLDeC)=varid
              CASE ('idTbry(inorth,iLDeC)')
                idTbry(inorth,iLDeC)=varid

              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid

              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid
#endif
#ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
#endif

#ifdef PHOSPHORUS
              CASE ('idTbry(iwest,iPO4_)')
                idTbry(iwest,iPO4_)=varid
              CASE ('idTbry(ieast,iPO4_)')
                idTbry(ieast,iPO4_)=varid
              CASE ('idTbry(isouth,iPO4_)')
                idTbry(isouth,iPO4_)=varid
              CASE ('idTbry(inorth,iPO4_)')
                idTbry(inorth,iPO4_)=varid
              CASE ('idTbry(iwest,iLDeP)')
                idTbry(iwest,iLDeP)=varid
              CASE ('idTbry(ieast,iLDeP)')
                idTbry(ieast,iLDeP)=varid
              CASE ('idTbry(isouth,iLDeP)')
                idTbry(isouth,iLDeP)=varid
              CASE ('idTbry(inorth,iLDeP)')
                idTbry(inorth,iLDeP)=varid
              CASE ('idTbry(iwest,iSDeP)')
                idTbry(iwest,iSDeP)=varid
              CASE ('idTbry(ieast,iSDeP)')
                idTbry(ieast,iSDeP)=varid
              CASE ('idTbry(isouth,iSDeP)')
                idTbry(isouth,iSDeP)=varid
              CASE ('idTbry(inorth,iSDeP)')
                idTbry(inorth,iSDeP)=varid
#endif
#ifdef H2S
              CASE ('idTbry(iwest,iH2S_)')
                idTbry(iwest,iH2S_)=varid
              CASE ('idTbry(ieast,iH2S_)')
                idTbry(ieast,iH2S_)=varid
              CASE ('idTbry(isouth,iH2S_)')
                idTbry(isouth,iH2S_)=varid
              CASE ('idTbry(inorth,iH2S_)')
                idTbry(inorth,iH2S_)=varid
#endif

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iPhyt)')
                idRtrc(iPhyt)=varid
              CASE ('idRtrc(iZoop)')
                idRtrc(iZoop)=varid
              CASE ('idRtrc(iLDeN)')
                idRtrc(iLDeN)=varid
              CASE ('idRtrc(iSDeN)')
                idRtrc(iSDeN)=varid
              CASE ('idRtrc(iChlo)')
                idRtrc(iChlo)=varid
#ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
              CASE ('idRtrc(iLDeC)')
                idRtrc(iLDeC)=varid
              CASE ('idRtrc(iSDeC)')
                idRtrc(iSDeC)=varid
#endif
#ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
#endif

#ifdef PHOSPHORUS
              CASE ('idRtrc(iPO4_)')
                idRtrc(iPO4_)=varid
              CASE ('idRtrc(iLDeP)')
                idRtrc(iLDeP)=varid
              CASE ('idRtrc(iSDeP)')
                idRtrc(iSDeP)=varid
#endif
#ifdef H2S
              CASE ('idRtrc(iH2S_)')
                idRtrc(iH2S_)=varid
#endif

#ifdef DIAGNOSTICS_BIO

/*
**  Biological tracers term diagnostics.
*/
              CASE ('iDbio2(iNH4f)')
                iDbio2(iNH4f)=varid
              CASE ('iDbio2(iPONf)')
                iDbio2(iPONf)=varid
# ifdef DENITRIFICATION
              CASE ('iDbio2(iDNIT)')
                iDbio2(iDNIT)=varid
# endif
# ifdef CARBON
              CASE ('iDbio2(iCOfx)')
                iDbio2(iCOfx)=varid
              CASE ('iDbio2(ipCO2)')
                iDbio2(ipCO2)=varid
# endif
# ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
              CASE ('iDbio2(iSODf)')
                iDbio2(iSODf)=varid
# endif
# ifdef PHOSPHORUS
              CASE ('iDbio2(iPO4f)')
                iDbio2(iPO4f)=varid
              CASE ('iDbio2(iPOPf)')
                iDbio2(iPOPf)=varid
# endif

              CASE ('iDbio3(iPPro)')
                iDbio3(iPPro)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
              CASE ('iDbio3(iPmax)')
                iDbio3(iPmax)=varid
              CASE ('iDbio3(iLDIN)')
                iDbio3(iLDIN)=varid
              CASE ('iDbio3(iNitr)')
                iDbio3(iNitr)=varid
              CASE ('iDbio3(iDeni)')
                iDbio3(iDeni)=varid
              CASE ('iDbio3(iAssi)')
                iDbio3(iAssi)=varid
              CASE ('iDbio3(iEges)')
                iDbio3(iEges)=varid
              CASE ('iDbio3(iPmor)')
                iDbio3(iPmor)=varid
              CASE ('iDbio3(iZmor)')
                iDbio3(iZmor)=varid
              CASE ('iDbio3(iZexc)')
                iDbio3(iZexc)=varid
              CASE ('iDbio3(iZmet)')
                iDbio3(iZmet)=varid
              CASE ('iDbio3(iCoaP)')
                iDbio3(iCoaP)=varid
              CASE ('iDbio3(iCoaD)')
                iDbio3(iCoaD)=varid
              CASE ('iDbio3(iSReN)')
                iDbio3(iSReN)=varid
              CASE ('iDbio3(iLReN)')
                iDbio3(iLReN)=varid
# ifdef PHOSPHORUS
              CASE ('iDbio3(iLPO4)')
                iDbio3(iLPO4)=varid
              CASE ('iDbio3(iSReP)')
                iDbio3(iSReP)=varid
              CASE ('iDbio3(iLReP)')
                iDbio3(iLReP)=varid
# endif
# ifdef OXYGEN
              CASE ('iDbio3(iCOD_)')
                iDbio3(iCOD_)=varid
              CASE ('iDbio3(iH2Sf)')
                iDbio3(iH2Sf)=varid
# endif
#endif
