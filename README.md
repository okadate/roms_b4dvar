# ROMS B4D-Var

ROMS with 4D-Var codes for Fennel's type Biogeochemical model (dev)

## CPP flags

### CPP flags for Fennel's type biogeochemical model

- PHOSPHORUS          use to add phosphorus variables (PO4, LdetritusP and SdetritusP)
- TDEPENDANCE         use to activate temperature dependance
- BIO_SED_CONSTANT    use to add constant bottom fluxes of DO, NH4 and PO4
- NPFLUX_BY_DO        use if bottom fluxes of NH4 and PO4 depending on DO
- BIO_SED_DIAGENESIS  (not maintenance)
- H2S                 (not maintenance)
- GROWTH1             (not maintenance)

### CPP flags for 4D-Var

- ADJUST_PARAM        use if including biogeochemical parameters in 4D-Var state
- RESTART_PARAM       use to read biochemical parameters
- POSITIVE_T          use to make tracers positive
- POSITIVE_PARAM      use to make parameters positive
- UV_FIXED_TL         use if physical model is fixed in B4D-Var
- OBS_DEPTH_SSH       use if obs_depth=0 is SSH
- EXP_PARAM           (not maintenance)

### CPP flags to check 4D-Var

- CHECKER             use to print user's diagnostic informations
- POSITIVE_CG_DELTA   use to make cg_delta positive
- POSITIVE_CG_RITZ    use to make cg_ritz positive

### other CPP flags

- DSTARTSEC           use to change units of dstart from days to seconds
