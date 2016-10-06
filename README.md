# ROMS B4D-Var

Regional Ocean Modeling Sysytem with Biological 4D-Var for Fennel et al. (2006) biology model

Original ROMS code: https://www.myroms.org/

## Additional CPP flags

### Fennel et al.(2006) biology model OPTIONS:

```
PHOSPHORUS          use to add phosphorus variables (PO4, LdetritusP and SdetritusP)
TDEPENDANCE         use to multiply reaction rate by coefficient depending on temperature
BIO_SED_CONSTANT    use to add constant bottom fluxes of DO, NH4 and PO4
NPFLUX_BY_DO        use if bottom fluxes of NH4 and PO4 depending on DO
BIO_SED_DIAGENESIS  (not maintenance) use to activate diagenesis model for bottom fluxes
H2S                 (not maintenance) use to add hydrogen sulfide (H2S)
GROWTH1             (not maintenance) use if growth rate of optimal temperature
```

### Biogeochemical 4D-Var OPTIONS:

```
ADJUST_PARAM        use if including biogeochemical parameters in 4D-Var state
RESTART_PARAM       use to read biochemical parameters
POSITIVE_T          use to make tracers positive
POSITIVE_PARAM      use to make parameters positive
UV_FIXED_TL         use if physical model is fixed in B4D-Var
OBS_DEPTH_SSH       use if obs_depth=0 is SSH
EXP_PARAM           (not maintenance) use if exponential parameters
```

### Checker for 4D-Var OPTIONS:

```
CHECKER             use to print user's diagnostic informations
POSITIVE_CG_DELTA   use to make cg_delta positive
POSITIVE_CG_RITZ    use to make cg_ritz positive
```

### Other OPTIONS:

```
DSTARTSEC           use to change units of dstart from days to seconds
USER_ENSEMBLE       use if ensemble runs by user parameters in ocean.in
```
