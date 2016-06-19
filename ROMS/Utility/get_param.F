#include "cppdefs.h"
      SUBROUTINE get_param (ng, model, msg, ncname, IniRec, Tindex)
!
!svn $Id$
!=======================================================================
!  Copyright (c) 2016 Teruhisa Okada                                   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in requested model parameter from specified      !
!  NetCDF file. It is usually used to read initial conditions.         !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     msg        Message index for StateMsg.                           !
!     ncname     Nonlinear initial conditions NetCDF file name.        !
!     IniRec     Nonlinear initial conditions time record to read.     !
!     Tindex     State variable time index to intialize.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
#ifdef FOUR_DVAR
      USE mod_fourdvar
#endif
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_strings
!
#ifdef ADJUST_PARAM
      USE nf_fread0d_mod, ONLY : nf_fread0d
#endif
      USE strings_mod, ONLY : find_string
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, msg, Tindex

      integer, intent(inout) :: IniRec

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Perfect2D, Perfect3D, foundit
      logical, dimension(NV) :: get_var, have_var

      integer :: LBi, UBi, LBj, UBj
      integer :: IDmod, InpRec, gtype, i, ifield, itrc, lstr, lend
      integer :: Nrec, mySize, ncINPid, nvatts, nvdim, status, varid
      integer :: Vsize(4), start(4), total(4)

      real(r8), parameter :: Fscl = 1.0_r8

      real(r8) :: Fmax, Fmin, INPtime, Tmax, scale, time_scale

      real(r8), allocatable :: TimeVar(:)

      character (len=6 ) :: string
      character (len=14) :: t_code
      character (len=15) :: attnam, tvarnam
      character (len=40) :: tunits
!
!-----------------------------------------------------------------------
!  Determine variables to read and their availability.
!-----------------------------------------------------------------------
!
!  Set model identification string.
!
      IF (model.eq.iNLM.or.(model.eq.0)) THEN
        string=' NLM: '                    ! nonlinear model, restart
        IDmod=iNLM
      ELSE IF (model.eq.iTLM) THEN
        string=' TLM: '                    ! tangent linear model
        IDmod=iTLM
      ELSE IF (model.eq.iRPM) THEN
        string=' RPM: '                    ! representer model
        IDmod=iRPM
      ELSE IF (model.eq.iADM) THEN
        string=' ADM: '                    ! adjoint model
        IDmod=iADM
      ELSE IF (model.eq.5) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! model or initial conditions
      ELSE IF (model.eq.6) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! model or initial conditions
      ELSE IF (model.eq.7) THEN
        string=' FRC: '                    ! impulse forcing
        IDmod=iNLM
      ELSE IF (model.eq.8) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.9) THEN
        string=' STD: '                    ! standard deviation
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.10) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! boundary conditions
      ELSE IF (model.eq.11) THEN
        string=' NRM: '                    ! normalization factor
        IDmod=iNLM                         ! surface forcing
      ELSE IF (model.eq.12) THEN
        string=' NLM: '                    ! tangent linear forcing and
        IDmod=iNLM                         ! obc increments
      END IF
!
!  Set switch to process variables for nonlinear model perfect restart.
!
      Perfect2D=.FALSE.
      Perfect3D=.FALSE.
      PerfectRST(ng)=Perfect2D.or.Perfect3D
!
!  Determine variables to read.
!
      CALL checkvars (ng, model, ncname, string, Nrec, NV, tvarnam,     &
     &                get_var, have_var)
      IF (exit_flag.ne.NoError) RETURN
!
!  Set Vsize to zero to deactivate interpolation of input data to model
!  grid in "nf_fread2d" and "nf_fread3d".
!
      DO i=1,4
        Vsize(i)=0
      END DO
      SourceFile='get_param.F'
!
!-----------------------------------------------------------------------
!  Open input NetCDF file and check time variable.
!-----------------------------------------------------------------------
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!  Open input NetCDF file.
!
      CALL netcdf_open (ng, IDmod, ncname, 0, ncINPid)
      IF (exit_flag.ne.NoError) THEN
        IF (Master) WRITE (stdout,10) TRIM(ncname)
        RETURN
      END IF
!
!  Inquire about the input time variable.
!
      CALL netcdf_inq_var (ng, IDmod, ncname,                           &
     &                     MyVarName = TRIM(tvarnam),                   &
     &                     VarID = varid,                               &
     &                     nVarDim =  nvdim,                            &
     &                     nVarAtt = nvatts)
      IF (exit_flag.ne.NoError) RETURN
!
!  Allocate input time variable and read its value(s).  Recall that
!  input time variable is a one-dimensional array with one or several
!  values.
!
      mySize=var_Dsize(1)
      IF (.not.allocated(TimeVar)) allocate (TimeVar(mySize))
      CALL netcdf_get_fvar (ng, IDmod, ncname, TRIM(tvarnam), TimeVar)
      IF (exit_flag.ne.NoError) RETURN
!
!  If using the latest time record from input NetCDF file as the
!  initialization record, assign input time.
!
      IF (LastRec(ng)) THEN
        Tmax=-1.0_r8
        DO i=1,mySize
          IF (TimeVar(i).gt.Tmax) THEN
            Tmax=TimeVar(i)
            IniRec=i
          END IF
        END DO
        INPtime=Tmax
        InpRec=IniRec
      ELSE
        IF ((IniRec.ne.0).and.(IniRec.gt.mySize)) THEN
          IF (Master)  WRITE (stdout,30) string, IniRec, TRIM(ncname),  &
     &                                   mySize
          exit_flag=2
          RETURN
        END IF
        IF (IniRec.ne.0) THEN
          InpRec=IniRec
        ELSE
          InpRec=1
        END IF
        INPtime=TimeVar(InpRec)
      END IF
      IF (allocated(TimeVar)) deallocate ( TimeVar )
!
!  Set input time scale by looking at the "units" attribute.
!
      time_scale=0.0_r8
      DO i=1,nvatts
        IF (TRIM(var_Aname(i)).eq.'units') THEN
          IF (INDEX(TRIM(var_Achar(i)),'day').ne.0) THEN
            time_scale=day2sec
          ELSE IF (INDEX(TRIM(var_Achar(i)),'second').ne.0) THEN
            time_scale=1.0_r8
          END IF
        END IF
      END DO
      IF (time_scale.gt.0.0_r8) THEN
        INPtime=INPtime*time_scale
      END IF
!
!  Set starting time index and time clock in days.  Notice that the
!  global time variables and indices are only over-written when
!  processing initial conditions (msg = 1).
!
      IF ((model.eq.0).or.(model.eq.iNLM).or.                           &
     &    (model.eq.iTLM).or.(model.eq.iRPM)) THEN
# ifdef DSTARTSEC
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(dstart))) THEN
          INPtime=dstart
        END IF
# else
        IF (((model.eq.iTLM).or.(model.eq.iRPM)).and.(msg.eq.1).and.    &
     &      (INPtime.ne.(dstart*day2sec))) THEN
          INPtime=dstart*day2sec
        END IF
# endif
        IF (msg.eq.1) THEN            ! processing initial conditions
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
# ifdef DSTARTSEC
          ntstart(ng)=NINT((time(ng)-dstart)/dt(ng))+1
# else
          ntstart(ng)=NINT((time(ng)-dstart*day2sec)/dt(ng))+1
# endif
          IF (ntstart(ng).lt.1) ntstart(ng)=1
          IF (PerfectRST(ng)) THEN
            ntfirst(ng)=1
          ELSE
            ntfirst(ng)=ntstart(ng)
          END IF
        END IF
#ifdef WEAK_CONSTRAINT
        IF (msg.eq.4) THEN
          ForceTime(ng)=time(ng)
        END IF
#endif
      ELSE IF (model.eq.iADM) THEN
        IF ((msg.eq.1).and.(INPtime.eq.0.0_r8)) THEN
          INPtime=time(ng)
        ELSE IF (msg.ne.1) THEN
          time(ng)=INPtime
          tdays(ng)=time(ng)*sec2day
        END IF
        ntstart(ng)=ntimes(ng)+1
        ntend(ng)=1
        ntfirst(ng)=ntend(ng)
      END IF
      CALL time_string (time(ng), time_code(ng))
!
!  Over-write "IniRec" to the actual initial record processed.
!
      IF (model.eq.iNLM) THEN
        IniRec=InpRec
      END IF
!
!  Set current input time, io_time .  Notice that the model time,
!  time(ng), is reset above.  This is a THREADPRIVATE variable in
!  shared-memory and this routine is only processed by the MASTER
!  thread since it is an I/O routine. Therefore, we need to update
!  time(ng) somewhere else in a parallel region. This will be done
!  with io_time variable.
!
      io_time=INPtime
!
!  Report information.
!
      lstr=SCAN(ncname,'/',BACK=.TRUE.)+1
      lend=LEN_TRIM(ncname)
      IF (Master) THEN
        CALL time_string (INPtime, t_code)
        IF (ERend.gt.ERstr) THEN
          WRITE (stdout,40) string, TRIM(StateMsg(msg)), t_code, ng,    &
     &                      Nrun, ncname(lstr:lend), InpRec, Tindex
        ELSE
          WRITE (stdout,50) string, TRIM(StateMsg(msg)), t_code, ng,    &
     &                      ncname(lstr:lend), InpRec, Tindex
        END IF
      END IF

#ifdef NONLINEAR
!
!-----------------------------------------------------------------------
!  Read in nonlinear state variables. If applicable, read in perfect
!  restart variables.
!-----------------------------------------------------------------------
!
      NLM_STATE: IF ((model.eq.iNLM).or.(model.eq.0).or.                &
     &               (model.eq.13)) THEN

# ifdef ADJUST_PARAM
!
!  Read in nonlinear parameters.
!
        DO itrc=1,Nparam(ng)
          IF (get_var(idPvar(itrc)).and.Lparam(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idPvar(itrc))), varid)
            status=nf_fread0d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idPvar(itrc)), varid,             &
     &                        InpRec,                                   &
     &                        OCEAN(ng) % p(Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idPvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idPvar(itrc))),          &
     &                            OCEAN(ng) % p(Tindex,itrc),           &
     &                            OCEAN(ng) % p(Tindex,itrc)
              END IF
            END IF
            OCEAN(ng)%p(:,itrc)=OCEAN(ng)%p(Tindex,itrc)
          END IF
        END DO
# endif
      END IF NLM_STATE
#endif

#if defined TANGENT || defined TL_IOMS
!
!-----------------------------------------------------------------------
!  Read in tangent linear state variables.
!-----------------------------------------------------------------------
!
      TLM_STATE: IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN

# ifdef ADJUST_PARAM
!
!  Read in tangent linear parameters.
!
        DO itrc=1,Nparam(ng)
          IF (get_var(idPvar(itrc)).and.Lparam(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idPvar(itrc))), varid)
            status=nf_fread0d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idPvar(itrc)), varid,             &
     &                        InpRec,                                   &
     &                        OCEAN(ng)% tl_p(Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idPvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idPvar(itrc))),          &
     &                            OCEAN(ng)%tl_p(Tindex,itrc),          &
     &                            OCEAN(ng)%tl_p(Tindex,itrc)
              END IF
            END IF
          END IF
        END DO
# endif
      END IF TLM_STATE
#endif

#ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Read in adjoint state variables.
!-----------------------------------------------------------------------
!
      ADM_STATE: IF (model.eq.iADM) THEN

# ifdef ADJUST_PARAM
!
!  Read in adjoint parameters.
!
        DO itrc=1,Nparam(ng)
          IF (get_var(idPvar(itrc)).and.Lparam(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idPvar(itrc))), varid)
            status=nf_fread0d(ng, IDmod, ncname, ncINPid,               &
     &                        Vname(1,idPvar(itrc)), varid,             &
     &                        InpRec,                                   &
     &                        OCEAN(ng) % ad_p(Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idPvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idPvar(itrc))),          &
     &                            OCEAN(ng)%ad_p(Tindex,itrc),          &
     &                            OCEAN(ng)%ad_p(Tindex,itrc)
              END IF
            END IF
          END IF
        END DO
# endif
      END IF ADM_STATE
#endif

#ifdef FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in error covariance normalization (nondimensional) factor.
!-----------------------------------------------------------------------
!
      NRM_STATE: IF ((model.eq. 5).or.                                  &
     &               (model.eq.10).or.                                  &
     &               (model.eq.11)) THEN

# ifdef ADJUST_PARAM_PASS
!
!  Read in parameter normalization factors.
!
        DO itrc=1,Nparam(ng)
          IF (get_var(idPvar(itrc)).and.(model.eq.11).and.              &
     &        Lparam(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idPvar(itrc))), varid)
            OCEAN(ng) % b_p(Tindex,itrc) = 1.0_r8 !okada! std_p(itrc)
            IF (Master) THEN
              WRITE (stdout,70) TRIM(Vname(2,idPvar(itrc))),            &
     &                          OCEAN(ng) % b_p(Tindex,itrc),           &
     &                          OCEAN(ng) % b_p(Tindex,itrc)
            END IF
          END IF
        END DO
# endif
      END IF NRM_STATE
#endif

#if defined FOUR_DVAR
!
!-----------------------------------------------------------------------
!  Read in error covariance standard deviation factors.
!-----------------------------------------------------------------------
!
      STD_STATE: IF ((model.eq.6).or.                                   &
     &               (model.eq.8).or.                                   &
     &               (model.eq.9)) THEN

# ifdef ADJUST_PARAM_PASS
!
!  Read in parameter standard deviations.
!
        DO itrc=1,Nparam(ng)
          IF (get_var(idPvar(itrc)).and.(model.eq.9).and.               &
     &        Lparam(itrc,ng)) THEN
            foundit=find_string(var_name, n_var,                        &
     &                          TRIM(Vname(1,idPvar(itrc))), varid)
            status=nf_fread0d(ng, IDmod, ncname, ncINPid,               &
                              Vname(1,idPvar(itrc)), varid,             &
     &                        InpRec,                                   &
     &                        OCEAN(ng) % e_p(Tindex,itrc))
            IF (status.ne.nf90_noerr) THEN
              IF (Master) THEN
                WRITE (stdout,60) string, TRIM(Vname(1,idPvar(itrc))),  &
     &                            InpRec, TRIM(ncname)
              END IF
              exit_flag=2
              ioerror=status
              RETURN
            ELSE
              IF (Master) THEN
                WRITE (stdout,70) TRIM(Vname(2,idPvar(itrc))),          &
     &                            OCEAN(ng) % e_p(Tindex,itrc),         &
     &                            OCEAN(ng) % e_p(Tindex,itrc)
              END IF
            END IF
          END IF
        END DO
# endif
      END IF STD_STATE
#endif

#if defined IMPULSE
!
!-----------------------------------------------------------------------
!  Read in adjoint model or tangent linear model impulse forcing terms.
!-----------------------------------------------------------------------
!
      FRC_STATE: IF (model.eq.7) THEN
!
!  pass
!
      END IF FRC_STATE
#endif

#if (defined W4DPSAS || defined TL_W4DPSAS || \
     defined W4DPSAS_SENSITIVITY)          && \
    (defined ADJUST_BOUNDARY || defined ADJUST_WSTRESS ||\
     defined ADJUST_STFLUX)
!
!-----------------------------------------------------------------------
!  Read in tangent linear forcing corrections.
!-----------------------------------------------------------------------
!
      TLM_FORCING: IF (model.eq.12) THEN
!
!  pass
!
      END IF TLM_FORCING
#endif
!
!-----------------------------------------------------------------------
!  Close input NetCDF file.
!-----------------------------------------------------------------------
!
      CALL netcdf_close (ng, IDmod, ncINPid, ncname, .FALSE.)
!
  10  FORMAT (/,a,'GET_PARAM - unable to open input NetCDF file: ',a)
  20  FORMAT (/,a,'GET_PARAM - Warning - NetCDF global attribute: ',a,  &
     &        /,18x,'for lateral boundary conditions not checked',      &
     &        /,18x,'in restart file: ',a)
  30  FORMAT (/,a,'GET_PARAM - requested input time record = ',i3,/,    &
     &        18x,'not found in input NetCDF: ',a,/,                    &
     &        18x,'number of available records = ',i3)
  40  FORMAT (/,a,'GET_PARAM - ',a,t62,'t = ',a,                        &
     &        /,19x,'(Grid ',i2.2,', Iter=',i4.4,', File: ',a,          &
     &        ', Rec=',i4.4,', Index=',i1,')')
  50  FORMAT (/,a,'GET_PARAM - ',a,t62,'t = ',a,                        &
     &        /,19x,'(Grid ',i2.2,', File: ',a,', Rec=',i4.4,           &
     &        ', Index=',i1,')')
  60  FORMAT (/,a,'GET_PARAM - error while reading variable: ',a,2x,    &
     &        'at time record = ',i3,/,18x,'in input NetCDF file: ',a)
  70  FORMAT (16x,'- ',a,/,19x,'(Min = ',1p,e15.8,                      &
     &        ' Max = ',1p,e15.8,')')

      RETURN
      END SUBROUTINE get_param
