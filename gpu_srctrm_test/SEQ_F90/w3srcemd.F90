#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRCEMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         07-Jan-2018 |
!/                  +-----------------------------------+
!/
!/    For updates see subroutine.
!/
!  1. Purpose :
!
!     Source term integration routine.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      OFFSET    R.P.  Private  Offset in time integration scheme.
!                               0.5 in original WAM, now 1.0
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. Public   Calculate and integrate source terms.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See corresponding documentation of W3SRCE.
!
!  5. Remarks :
!
!  6. Switches :
!
!       See section 9 of W3SRCE.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      REAL, PARAMETER, PRIVATE:: OFFSET = 1.
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SRCE_INIT()

      USE W3GDATMD, ONLY: NK, NSPEC, NTH, NSEAL

      IMPLICIT NONE
!
! Allocate scratch variables for GPU.
  
      WRITE(0,*) "W3SRCE_INIT: NK/NTH/NSPEC: " ,NK,NTH,NSPEC

      END SUBROUTINE W3SRCE_INIT

!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SRCE ( srce_call, IT, IMOD,          &
                          SPECOLD, SPEC,& 
                          VSIO, VDIO, SHAVEIO,         &
                          ALPHA, WN1, &
                          CG1, D_INP, U10ABS, &
                          U10DIR, AS, USTAR,&
                          USTDIR, CX, CY, &
                          ICE, ICEH, ICEF,&
                          ICEDMAX,          &
                          REFLEC, REFLED, &
                          TRNX, TRNY, BERG, &
                          FPI, DTDYN, &
                          FCUT, DTG, TAUWX, TAUWY, &
                          TAUOX, TAUOY, TAUWIX, &
                          TAUWIY, TAUWNX,&
                          TAUWNY, PHIAW, CHARN,&
                          TWS, PHIOC, WHITECAP, D50, PSIC,&
                          BEDFORMS, PHIBBL, TAUBBL,&
                          TAUICE, PHICE, COEF,&
                          SIN4T, SPR4T, SDS4T)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            F. Ardhuin             |
!/                  |            A. Roland              |
!/                  |            M. Dutour Sikiric      |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         07-Jan-2018 |
!/                  +-----------------------------------+
!/
!/    06-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    04-Feb-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    14-Feb-2000 : Exact-NL added                      ( version 2.01 )
!/    04-May-2000 : Non-central integration             ( version 2.03 )
!/    02-Feb-2001 : Xnl version 3.0                     ( version 2.07 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    27-Nov-2002 : First version of VDIA and MDIA.     ( version 3.01 )
!/    07-Oct-2003 : Output options for NN training.     ( version 3.05 )
!/    24-Dec-2004 : Multiple model version.             ( version 3.06 )
!/    23-Jun-2006 : Linear input added.                 ( version 3.09 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    04-Jul-2006 : Separation of stress computation.   ( version 3.09 )
!/    16-Apr-2007 : Miche style limiter added.          ( version 3.11 )
!/                  (J. H. Alves)
!/    25-Apr-2007 : Battjes-Janssen Sdb added.          ( version 3.11 )
!/                  (J. H. Alves)
!/    09-Oct-2007 : Adding WAM 4+ and SB1 options.      ( version 3.13 )
!/                  (F. Ardhuin)
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    19-Aug-2010 : Making treatment of 0 water depth   ( version 3.14.6 )
!/                  consistent with the rest of the model.
!/    31-Mar-2010 : Adding ice conc. and reflections    ( version 3.14.4 )
!/    15-May-2010 : Adding transparencies               ( version 3.14.4 )
!/    01-Jun-2011 : Movable bed bottom friction in BT4  ( version 4.01 )
!/    01-Jul-2011 : Energy and momentum flux, friction  ( version 4.01 )
!/    24-Aug-2011 : Uses true depth for depth-induced   ( version 4.04 )
!/    16-Sep-2011 : Initialization of TAUWAX, TAUWAY    ( version 4.04 )
!/     1-Dec-2011 : Adding BYDRZ source term package    ( version 4.04 )
!/                  ST6 and optional Hwang (2011)
!/                  stresses FLX4.
!/    14-Mar-2012 : Update of BT4, passing PSIC         ( version 4.04 )
!/    13-Jul-2012 : Move GMD (SNL3) and nonlinear filter (SNLS)
!/                  from 3.15 (HLT).                    ( version 4.08 )
!/    28-Aug-2013 : Corrected MLIM application          ( version 4.11 )
!/    10-Sep-2013 : Special treatment for IG band       ( version 4.15 )
!/    14-Nov-2013 : Make orphaned pars in par lst local ( version 4.13 )
!/    17-Nov-2013 : Coupling fraction of ice-free       ( version 4.13 )
!/                  surface to SIN and SDS. (S. Zieger)
!/    01-Avr-2014 : Adding ice thickness and floe size  ( version 4.18 )
!/    23-May-2014 : Adding ice fluxes to W3SRCE         ( version 5.01 )
!/    27-Aug-2015 : Adding inputs to function W3SIS2    ( version 5.10 )
!/    13-Dec-2015 : Implicit integration of Sice (F.A.) ( version 5.10 )
!/    30-Jul-2017 : Adds TWS in interface               ( version 6.04 )
!/    07-Jan-2018 : Allows variable ice scaling (F.A.)  ( version 6.04 )
!/    01-Jan-2018 : Add implicit source term integration ( version 6.04)
!/    01-Jan-2018 : within PDLIB (A. Roland, M. Dutour
!/    18-Aug-2018 : S_{ice} IC5 (Q. Liu)                ( version  6.06)
!/    26-Aug-2018 : UOST (Mentaschi et al. 2015, 2018)  ( version 6.06 )
!/
!/    Copyright 2009-2013 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Calculate and integrate source terms for a single grid point.
!
!  2. Method :
!
!     Physics  : see manual and corresponding subroutines.
!
!     Numerics :
!
!     Dynamic-implicit integration of the source terms based on
!     WW-II (Tolman 1992). The dynamic time step is calculated
!     given a maximum allowed change of spectral densities for
!     frequencies / wavenumbers below the usual cut-off.
!     The maximum change is given by the minimum of a parametric
!     and a relative change. The parametric change relates to a
!     PM type equilibrium range
!
!                                -1  (2pi)**4       1
!       dN(k)     =  Xp alpha  pi   ---------- ------------
!            max                       g**2     k**3 sigma
!
!                              1                                     .
!                 =  FACP ------------                              (1)
!                          k**3 sigma                                .
!
!     where
!           alpha = 0.62e-4                       (set in W3GRID)
!           Xp      fraction of PM shape          (read in W3GRID)
!           FACP    combined factor               (set in W3GRID)
!
!     The maximum relative change is given as
!
!                           /            +-                  -+ \    .
!       dN(k)     =  Xr max | N(k) , max | Nx , Xfilt N(k)    | |   (2)
!            max            \            +-               max-+ /    .
!
!     where
!           Xr      fraction of relative change   (read in W3GRID)
!           Xfilt   filter level                  (read in W3GRID)
!           Nx      Maximum parametric change (1)
!                   for largest wavenumber.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IX,IY   Int.   I   Discrete grid point counters.
!       IMOD    Int.   I   Model number.
!       SPEC    R.A.  I/O  Spectrum (action) in 1-D form.
!       ALPHA   R.A.  I/O  Nondimenional 1-D spectrum corresponding
!                          to above full spectra (Phillip's const.).
!                          Calculated separately for numerical
!                          economy on vector machine (W3SPR2).
!       WN1     R.A.   I   Discrete wavenumbers.
!       CG1     R.A.   I   Id. group velocities.
!       D_INP   Real.  I   Depth. Compared to DMIN to get DEPTH.
!       U10ABS  Real.  I   Wind speed at reference height.
!       U10DIR  Real.  I   Id. wind direction.
!       AS      Real.  I   Air-sea temp. difference.      ( !/ST3 )
!       USTAR   Real. !/O  Friction velocity.
!       USTDIR  Real!/O  Idem, direction.
!       CX-Y    Real.  I   Current velocity components.   ( !/BS1 )
!       ICE     Real   I   Sea ice concentration
!       ICEH    Real   I   Sea ice thickness
!       ICEF    Real  I/O  Sea ice maximum floe diameter  (updated)
!       ICEDMAX Real  I/O  Sea ice maximum floe diameter
!       BERG    Real   I   Iceberg damping coefficient    ( !/BS1 )
!       REFLEC  R.A.   I   reflection coefficients        ( !/BS1 )
!       REFLED  I.A.   I   reflection direction           ( !/BS1 )
!       TRNX-Y  Real   I   Grid transparency in X and Y   ( !/BS1 )
!       DELX    Real.  I   grid cell size in X direction  ( !/BS1 )
!       DELY    Real.  I   grid cell size in Y direction  ( !/BS1 )
!       DELA    Real.  I   grid cell area                 ( !/BS1 )
!       FPI     Real  I/O  Peak-input frequency.          ( !/ST2 )
!      WHITECAP R.A.   O   Whitecap statisics             ( !/ST4 )
!       DTDYN   Real   O   Average dynamic time step.
!       FCUT    Real   O   Cut-off frequency for tail.
!       DTG     Real   I   Global time step.
!       D50     Real   I   Sand grain size                ( !/BT4 )
!       BEDFORM R.A.  I/O  Bedform parameters             ( !/BT4 )
!       PSIC    Real   I   Critical Shields               ( !/BT4 )
!       PHIBBL  Real   O   Energy flux to BBL             ( !/BTx )
!       TAUBBL  R.A.   O   Momentum flux to BBL           ( !/BTx )
!       TAUICE  R.A.   O   Momentum flux to sea ice       ( !/ICx )
!       PHICE   Real   O   Energy flux to sea ice         ( !/ICx )
!     ----------------------------------------------------------------
!       Note: several pars are set to I/O to avoid compiler warnings.
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SPRn    Subr. W3SRCnMD Mean wave parameters for use in
!                               source terms.
!      W3FLXn    Subr. W3FLXnMD Flux/stress computation.
!      W3SLNn    Subr. W3SLNnMD Linear input.
!      W3SINn    Subr. W3SRCnMD Input source term.
!      W3SNLn    Subr. W3SNLnMD Nonlinear interactions.
!      W3SNLS    Subr. W3SNLSMD Nonlinear smoother.
!      W3SDSn    Subr. W3SRCnMD Whitecapping source term
!      W3SBTn    Subr. W3SBTnMD Bottom friction source term.
!      W3SDBn    Subr. W3SBTnMD Depth induced breaking source term.
!      W3STRn    Subr. W3STRnMD Triad interaction source term.
!      W3SBSn    Subr. W3SBSnMD Bottom scattering source term.
!      W3REFn    Subr. W3REFnMD Reflexions (shore, icebergs ...).
!      W3SXXn    Subr. W3SXXnMD Unclassified source term.
!      STRACE    Subr. W3SERVMD Subroutine tracing (!/S)
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - No testing is performed on the status of the grid point.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!       1.   Preparations
!         a  Set maximum change and wavenumber arrays.
!         b  Prepare dynamic time stepping.
!         c  Compute mean parameters.                       ( W3SPRn )
!         d  Compute stresses (if posible).
!         e  Prepare cut-off
!         f  Test output for !/NNT option.
!     --start-dynamic-integration-loop---------------------------------
!       2.  Calculate source terms
!         a Input.                                  ( W3SLNx, W3SINn )
!         b Nonlinear interactions.                         ( W3SNLn )
!         c Dissipation                                     ( W3SDSn )
!           1 as included in source terms                   ( W3SDSn )
!           2 optional dissipation due to different physics ( W3SWLn )
!         d Bottom friction.                                ( W3SBTn )
!       3.  Calculate cut-off frequencie(s)
!       4.  Summation of source terms and diagonal term and time step.
!       5.  Increment spectrum.
!       6.  Add tail
!         a Mean wave parameters and cut-off                ( W3SPRn )
!         b 'Seeding' of spectrum.                          ( !/SEED )
!         c Add tail
!       7.  Check if integration complete.
!     --end-dynamic-integration-loop-----------------------------------
!       8.  Save integration data.
!     -----------------------------------------------------------------
!
!  9. Switches :
!
!   !/FLX1  Wu (1980) stress computation.              ( Choose one )
!   !/FLX2  T&C (1996) stress computation.
!   !/FLX3  T&C (1996) stress computation with cap.
!   !/FLX4  Hwang (2011) stress computation (2nd order).
!
!   !/LN0   No linear input.                           ( Choose one )
!   !/LNX   User-defined bottom friction.
!
!   !/ST0   No input and dissipation.                  ( Choose one )
!   !/ST1   WAM-3 input and dissipation.
!   !/ST2   Tolman and Chalikov (1996)  input and dissipation.
!   !/ST3   WAM 4+ input and dissipation.
!   !/ST4   Ardhuin et al. (2009, 2010)
!   !/ST6   BYDB source terms after Babanin, Young, Donelan and Banner.
!   !/STX   User-defined input and dissipation.
!
!   !/NL0   No nonlinear interactions.                 ( Choose one )
!   !/NL1   Discrete interaction approximation.
!   !/NL2   Exact nonlinear interactions.
!   !/NL3   Generalized Multiple DIA.
!   !/NL4   Two Scale Approximation
!   !/NLX   User-defined nonlinear interactions.
!   !/NLS   Nonlinear HF smoother.
!
!   !/BT0   No bottom friction.                        ( Choose one )
!   !/BT1   JONSWAP bottom friction.
!   !/BT4   Bottom friction using movable bed roughness
!                  (Tolman 1994, Ardhuin & al. 2003)
!   !/BT8   Muddy bed (Dalrymple & Liu).
!   !/BT9   Muddy bed (Ng).
!   !/BTX   User-defined bottom friction.
!
!   !/IC1   Dissipation via interaction with ice according to simple
!             methods: 1) uniform in frequency or
!   !/IC2            2) Liu et al. model
!   !/IC3   Dissipation via interaction with ice according to a
!             viscoelastic sea ice model (Wang and Shen 2010).
!   !/IC4   Dissipation via interaction with ice as a function of freq.
!             (empirical/parametric methods)
!   !/IC5   Dissipation via interaction with ice according to a
!             viscoelastic sea ice model (Mosig et al. 2015).
!   !/DB0   No depth-limited breaking.                 ( Choose one )
!   !/DB1   Battjes-Janssen depth-limited breaking.
!   !/DBX   User-defined bottom friction.
!
!   !/TR0   No triad interactions.                     ( Choose one )
!/
!   !/TR1   Lumped Triad Approximation (LTA).
!   !/TRX   User-defined triad interactions.
!
!   !/BS0   No bottom scattering.                      ( Choose one )
!   !/BS1   Scattering term by Ardhuin and Magne (2007).
!   !/BSX   User-defined bottom scattering
!
!   !/XX0   No arbitrary additional source term.       ( Choose one )
!   !/XXX   User-defined bottom friction.
!
!   !/MLIM  Miche style limiter for shallow water and steepness.
!
!   !/SEED  'Seeding' of lowest frequency for suffuciently strong
!             winds.
!
!   !/NNT   Write output to file test_data_NNN.ww3 for NN training.
!
!   !/S     Enable subroutine tracing.
!   !/T     Enable general test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, NSEAL, NX, NY, SIG, TH, DMIN,&
                          DTMAX, DTMIN, FACTI1, FACTI2, FACSD, FACHFA, &
                          XFC, XFLT, XREL, XFT, FXFM, FXPM, DDEN,      &
                          FTE, FTF, FHMAX, ECOS, ESIN, IICEDISP,       &
                          ICESCALES, IICESMOOTH, MAPSF, MAPSTA, FACP,  &
                          FLAGST, NK2, MPARS
      USE W3GDATMD, ONLY: FSSOURCE, optionCall
      USE W3GDATMD, ONLY: B_JGS_NLEVEL, B_JGS_SOURCE_NONLINEAR
      USE W3WDATMD, ONLY: TIME
      USE W3ODATMD, ONLY: NDSE, NDST, NDTO, IAPROC
      USE W3IDATMD, ONLY: INFLAGS1, INFLAGS2, ICEP2
      USE W3DISPMD
      USE W3SLN1MD
      USE W3SRC4MD, ONLY : W3SPR4, W3SIN4, W3SDS4
      USE W3GDATMD, ONLY : ZZWND, FFXFM, FFXPM, FFXFA, FFXFI, FFXFD
      USE W3SNL1MD
      USE W3PARALL, ONLY : INIT_GET_ISEA
!/
      IMPLICIT NONE
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!CODENotes: Specifying the correct bounds for input variables here, and 
! on the call to the fuction are critical to avoid confusing validation
! errors.

      INTEGER, INTENT(IN)    :: srce_call, IT, IMOD
      REAL, intent(in)       :: SPECOLD(1:NSPEC)
      REAL, INTENT(OUT)      :: VSIO(1:NSPEC), VDIO(1:NSPEC)
      LOGICAL, INTENT(OUT)   :: SHAVEIO
      REAL, INTENT(IN)       :: D_INP(1:), U10ABS(1:NSEAL),            &
                                U10DIR(1:NSEAL), AS(1:), CX(1:),       &
                                CY(1:), DTG, D50, PSIC, ICE(1:),       &
                                ICEH(1:)
      INTEGER, INTENT(IN)    :: REFLED(6)
      REAL, INTENT(IN)       :: REFLEC(4), TRNX(:,:), TRNY(:,:),       &
                                BERG(:), ICEDMAX(:)
      REAL, INTENT(INOUT)    :: WN1(1:NK,1:NSEAL), CG1(1:NK,1:NSEAL),  &
                                SPEC(1:NSPEC,1:NSEAL), ALPHA(:,:),     &
                                USTAR(1:NSEAL),USTDIR(1:NSEAL),FPI(:), &
                                TAUOX(1:NSEAL), TAUOY(1:NSEAL),        &
                                TAUWX(1:NSEAL), TAUWY(1:NSEAL),        &
                                PHIAW(1:NSEAL), PHIOC(1:NSEAL),        &
                                PHICE(1:NSEAL), CHARN(1:NSEAL),        &
                                TWS(1:NSEAL), BEDFORMS(1:NSEAL,1:3),   &
                                PHIBBL(1:NSEAL), TAUBBL(1:NSEAL,1:2),  &
                                TAUICE(1:NSEAL,1:2),                   &
                                WHITECAP(1:NSEAL,1:4), TAUWIX(1:NSEAL),&
                                TAUWIY(1:NSEAL), TAUWNX(1:NSEAL),    &
                                TAUWNY(1:NSEAL),ICEF(1:NSEAL)
      REAL, INTENT(OUT)      :: DTDYN(1:NSEAL), FCUT(1:NSEAL)
      REAL, INTENT(IN)       :: COEF(:)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!CODENotes: Changes to local variables to allow managed memory to work
!properly, requires POINTERS to be defined as ALLOCATABLES and then
!allocated. 
      INTEGER                :: IK, ITH, IS, IS0, NKH, NKH1, NSTEPS,   &
                                IKS1, IS1, NSPECH, IDT, IERR, NKI, NKD,&
                                ISPEC, JSEA, ISEA, IX, IY
      REAL                   :: FHIGH, DT, AFILT, DAMAX, AFAC,         &
                                 HDT, ZWND, FP, FHIGI, KDMEAN
! Scaling factor for SIN, SDS, SNL
      REAL                   :: ICESCALELN, ICESCALEIN,          &
                                ICESCALENL, ICESCALEDS
      REAL                   :: EMEAN, FMEAN, WNMEAN, AMAX, CD, Z0,    &
                                SCAT, SMOOTH_ICEDISP,ICECOEF2, DTTOT,  &
                                PHINL, TAUWAX, TAUWAY, TAUSCX, TAUSCY

      REAL,ALLOCATABLE       :: WN_R(:,:),CG_ICE(:,:),ALPHA_LIU(:,:),  &
                                R(:,:), SPECINIT(:,:), SPEC2(:,:)
      REAL,ALLOCATABLE       :: DAM (:,:), WN2 (:,:), VSLN(:,:),       &
                                VSIN(:,:), VDIN(:,:),VSNL(:,:),        &
                                VDNL(:,:), VSDS(:,:), VDDS(:,:),       &
                                VSBT(:,:), VDBT(:,:), VS(:,:),         &
                                VD(:,:),COSI(:,:),DEPTH(:)
      LOGICAL,ALLOCATABLE    :: LLWS(:,:)
      REAL,ALLOCATABLE       :: FOUT(:,:), SOUT(:,:), DOUT(:,:)
      INTEGER                :: I,JJ
      DOUBLE PRECISION       :: ATT, ISO
      REAL                   :: FMEANS, FH1, FH2, FAGE
      REAL                   :: QCERR  = 0.   !/XNL2 and !/NNT
      REAL                   :: EBAND, DIFF, EFINISH, HSTOT,           &
                                FMEAN1, FMEANWS, MWXINIT, MWYINIT,     &
                                FACTOR, FACTOR2, DRAT, MWXFINISH,      &
                                MWYFINISH, A1BAND, B1BAND
      LOGICAL, SAVE          :: FLTEST = .FALSE., FLAGNN = .TRUE.
      LOGICAL                :: SHAVE
      LOGICAL                :: LBREAK
      LOGICAL, SAVE          :: FIRST = .TRUE.
      LOGICAL                :: PrintDeltaSmDA
      REAL                   :: eInc1, eInc2
      REAL                   :: DELX, DELY, DELA
!      REAL                   :: DeltaSRC(NSPEC), MAXDAC(NSPEC)
 
!!
!!LS  Newly added time varibles
      REAL(8)                :: sTime1, eTime1, sTime2, eTime2,        &
                                sTime3, eTime3
      REAL(8), INTENT(OUT)   :: SIN4T, SPR4T, SDS4T
     
      REAL,ALLOCATABLE       :: TMP1(:,:), TMP2(:,:), TMP3(:,:),       &
                                TMP4(:,:)
      REAL,ALLOCATABLE       :: BRLAMBDA(:,:)
!/
!/ ------------------------------------------------------------------- /
!/
  

!GPUNotes SIN4 requires these values from constants.f90, however it does
!not automatically transfer the data even if they are declared with
!copyin. This seemed to be the best place to put the updates. 
      WRITE(0,*)'TAG: W3SRCE'
      SPR4T = 0.0
      SIN4T = 0.0
      SDS4T = 0.0
!
!Allocating arrays here, this are local to routine and typically treated
!as scalars. In order to allow parallelism over sea points they have
!been converted. Deallocated at the end of this routine. 

!/
      ALLOCATE(DAM(NSPEC,NSEAL), WN2(NSPEC,NSEAL), VSLN(NSPEC,NSEAL),  &
      SPECINIT(NSPEC,NSEAL), SPEC2(NSPEC,NSEAL), VSIN(NSPEC,NSEAL),    &
      VDIN(NSPEC,NSEAL), VSNL(NSPEC,NSEAL), VDNL(NSPEC,NSEAL),         &
      VSDS(NSPEC,NSEAL), VDDS(NSPEC,NSEAL), VSBT(NSPEC,NSEAL),         &
      VDBT(NSPEC,NSEAL), VS(NSPEC,NSEAL), VD(NSPEC,NSEAL),             &
      FOUT(NK,NTH), SOUT(NK,NTH), DOUT(NK,NTH), WN_R(NK,NSEAL),        &
      CG_ICE(NK,NSEAL), ALPHA_LIU(NK,NSEAL), R(NK,NSEAL),              &
      COSI(2,NSEAL), LLWS(NSPEC,NSEAL))

      ALLOCATE(BRLAMBDA(NSPEC,NSEAL), DEPTH(NSEAL), &
      TMP1(4,NSEAL),TMP2(3,NSEAL), TMP3(2,NSEAL), TMP4(2,NSEAL))

#ifdef MM
#else
! The remaining data requirements can be copied in implicitly corretly.
! The use of IX, IY in MAPSTA stop this from working implicitly. 
!$ACC DATA  copyin(mapsta(:,:))
#endif

! Data is required in multiple places and not updated so is simple
! placed on the GPU and not moved. 
!$ACC ENTER DATA COPYIN(NK,NK2, NTH, NSPEC, NSEAL)

!$ACC KERNELS

      DEPTH(:)  = 0.
      PHIAW(:)  = 0.
      DTDYN(:)  = 0.
      CHARN(:)  = 0.
      TWS(:)    = 0.
      PHIBBL(:) = 0.
      TAUWIX(:) = 0.
      TAUWIY(:) = 0.
      TAUWNX(:) = 0.
      TAUWNY(:) = 0.
      PHICE(:) = 0.
      TAUWX(:)=0.
      TAUWY(:)=0.

      TMP3(:,:) = 0.
      TMP4(:,:) = 0.
      BRLAMBDA(:,:)=0.

      VSIN(:,:) = 0.
      VDIN(:,:) = 0.
      VDNL(:,:) = 0.
      VSDS(:,:) = 0.
      VDDS(:,:) = 0.
      VSLN(:,:) = 0.
      VSNL(:,:) = 0.
      VSBT(:,:) = 0.
      VDBT(:,:) = 0.
      ZWND   = ZZWND
      DRAT  = DAIR / DWAT
      IKS1 = 1
!GPUNotes Outer seapoint loop for source term calculations
!$ACC LOOP GANG VECTOR INDEPENDENT PRIVATE(EMEAN,FMEAN,FMEAN1,WNMEAN   )&
!$ACC PRIVATE(Z0, FMEANWS, IX, IY, TAUWAX, TAUWAY)&
!$ACC PRIVATE(AMAX, CD, ICESCALELN, ICESCALEIN, ICESCALEDS, ICESCALENL)
      DO ISEA=1, NSEAL
!GPUNotes removed INIT_GET_ISEA as its functions is simple but its
!application affected the parallelism of the kernel. 
!        CALL INIT_GET_ISEA(ISEA, JSEA)
      DTTOT=0.
      PHINL=0.
      TAUWAX=0.
      TAUWAY=0.
      TAUSCX=0.
      TAUSCY=0.
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        IF ( MAPSTA(IY,IX) .EQ. 1 .AND. FLAGST(ISEA)) THEN
          TMP1(:,ISEA)   = WHITECAP(ISEA,:)
          TMP2(:,ISEA)   = BEDFORMS(ISEA,:)

          TMP3(:,ISEA)   = TAUBBL(ISEA,:)
          TMP4(:,ISEA)   = TAUICE(ISEA,:)
          DEPTH(ISEA)  = MAX ( DMIN , D_INP(ISEA) )
          ICESCALELN = MAX(0.,MIN(1.,1.-ICE(ISEA)*ICESCALES(1)))
          ICESCALEIN = MAX(0.,MIN(1.,1.-ICE(ISEA)*ICESCALES(2)))
          ICESCALENL = MAX(0.,MIN(1.,1.-ICE(ISEA)*ICESCALES(3)))
          ICESCALEDS = MAX(0.,MIN(1.,1.-ICE(ISEA)*ICESCALES(4)))
          IS1=(IKS1-1)*NTH+1
!
!
! 1.  Preparations --------------------------------------------------- *
!
! 1.a Set maximum change and wavenumber arrays.
!
!XP     = 0.15
!FACP   = XP / PI * 0.62E-3 * TPI**4 / GRAV**2
  
!GPUNotes loop over frequencies
!$ACC LOOP 
          DO IK=1, NK
            DAM(1+(IK-1)*NTH,ISEA) = FACP / ( SIG(IK) * WN1(IK,ISEA)**3)
            WN2(1+(IK-1)*NTH,ISEA) = WN1(IK,ISEA)
          END DO
  
!GPUNotes loop over full spectrum
!$ACC LOOP COLLAPSE(2)
          DO IK=1, NK
            DO ITH=2, NTH
              IS0    = (IK-1)*NTH
              DAM(ITH+IS0,ISEA) = DAM(1+IS0,ISEA)
              WN2(ITH+IS0,ISEA) = WN2(1+IS0,ISEA)
            END DO
          END DO
!
! 1.b Prepare dynamic time stepping
!
!
  
! 1.c Set mean parameters
!
!GPUNotes calls to W3SPR4 and W3SIN4 below will contain source term specific spectral loops
!GPUNotes the sequencing is important (although maybe excessive?)
  
          IF (IT .eq. 0) THEN
             LLWS(:,ISEA) = .TRUE.
             USTAR(ISEA)=0.
             USTDIR(ISEA)=0.
          ELSE
             CALL W3SPR4(SPEC(:,ISEA), CG1(:,ISEA), WN1(:,ISEA), EMEAN,&
                         FMEAN, FMEAN1, WNMEAN, AMAX, U10ABS(ISEA),    &
                         U10DIR(ISEA), USTAR(ISEA), USTDIR(ISEA),      &
                         TAUWX(ISEA), TAUWY(ISEA), CD, Z0, CHARN(ISEA),&
                         LLWS(:,ISEA), FMEANWS,ISEA)
  
             CALL W3SIN4(SPEC(:,ISEA), CG1(:,ISEA), WN2(:,ISEA),       &
                         U10ABS(ISEA), USTAR(ISEA), DRAT, AS(ISEA),    &
                         U10DIR(ISEA), Z0, CD, TAUWX(ISEA),            &
                         TAUWY(ISEA), TAUWAX, TAUWAY,      &
                         VSIN(:,ISEA), VDIN(:,ISEA), LLWS(:,ISEA), IX, &
                         IY, BRLAMBDA(:,ISEA))
          END IF
!GPUNotes call below will contain source term specific spectral loops
!      CALL CPU_TIME(sTime1)
          CALL W3SPR4 (SPEC(:,ISEA), CG1(:,ISEA), WN1(:,ISEA), EMEAN,&
                         FMEAN, FMEAN1, WNMEAN, AMAX, U10ABS(ISEA),    &
                         U10DIR(ISEA), USTAR(ISEA), USTDIR(ISEA),      &
                         TAUWX(ISEA), TAUWY(ISEA), CD, Z0, CHARN(ISEA),&
                         LLWS(:,ISEA), FMEANWS,ISEA)
!      CALL CPU_TIME(eTime1)
!      SPR4T = SPR4T + eTime1 - sTime1
          TWS(ISEA) = 1./FMEANWS
!
! 1.c2 Stores the initial data
!
          SPECINIT(:,ISEA) = SPEC(:,ISEA)
! 1.d Stresses
!
! 1.e Prepare cut-off beyond which the tail is imposed with a power law
!
! Introduces a Long & Resio (JGR2007) type dependance on wave age
! !/ST4      FAGE   = FFXFA*TANH(0.3*U10ABS*FMEANWS*TPI/GRAV)
          FAGE   = 0.
          FHIGH  = MAX( (FFXFM + FAGE ) * MAX(FMEAN1,FMEANWS), FFXPM &
                   / USTAR(ISEA))
          FHIGI  = FFXFA * FMEAN1
       
!
! 1.f Prepare output file for !/NNT option
!
! ... Branch point dynamic integration - - - - - - - - - - - - - - - -
!GPUNotes loop for explicit time integration of source terms
!GPUNotes this might be a pain in the proverbial for tightening
!GPUNotes the seapoint and spectral loops

!GPUNotes Bodged DO WHILE LOOP, it is required to be sequential but the
!compiler panics. Setting it to loop over a larger extent and leaving
!the EXIT statements in gives us the output we require.


!$ACC LOOP SEQ
          DO I=1,120 
            NSTEPS = NSTEPS + 1
! 2.  Calculate source terms ----------------------------------------- *
! 2.a Input.
!
            CALL W3SLN1 (WN1(:,ISEA), FHIGH, USTAR(ISEA), &
            U10DIR(ISEA), VSLN(:,ISEA) )
            NKH    = MIN (NK, INT(FACTI2+FACTI1*LOG(MAX(1.E-7,FHIGH))))
            NKH1   = MIN (NK, NKH+1 )
            NSPECH = NKH1*NTH
            DT     = MIN (DTG-DTTOT, DTMAX )
            AFILT  = MAX (DAM(NSPEC,ISEA), XFLT*AMAX )
  
!        CALL CPU_TIME(sTime2) 
            CALL W3SIN4(SPEC(:,ISEA), CG1(:,ISEA), WN2(:,ISEA),        &
                        U10ABS(ISEA), USTAR(ISEA), DRAT, AS(ISEA),     &
                        U10DIR(ISEA), Z0, CD, TAUWX(ISEA), TAUWY(ISEA),&
                        TAUWAX, TAUWAY, VSIN(:,ISEA),      &
                        VDIN(:,ISEA), LLWS(:,ISEA), IX, IY,            &
                        BRLAMBDA(:,ISEA))
!        CALL CPU_TIME(eTime2)
!        SIN4T = SIN4T + eTime2 - sTime2
!
! 2.b Nonlinear interactions.
  
!GPUnotes subruoutine will contain source term specific spectral loops
  
            CALL W3SNL1 (SPEC(:,ISEA), CG1(:,ISEA), DEPTH(ISEA)*WNMEAN,&
                         VSNL(:,ISEA), VDNL(:,ISEA), ISEA )
  
! 2.c Dissipation... except for ST4
! 2.c1   as in source term package
  
!GPUNotes subroutine will contain source term specific spectral loops
  
!        CALL CPU_TIME(sTime3)
            CALL W3SDS4(SPEC(:,ISEA), WN1(:,ISEA), CG1(:,ISEA),        &
                        USTAR(ISEA), USTDIR(ISEA), DEPTH(ISEA),        &
                        VSDS(:,ISEA), VDDS(:,ISEA), IX, IY,            &
                        BRLAMBDA(:,ISEA),TMP1(1:4,ISEA), ISEA)  
          
!        CALL CPU_TIME(eTime3)
!        SDS4T = SDS4T + eTime3 - sTime3
   
!
! 2.c2   optional dissipation parameterisations
!
! 2.d Bottom interactions.
!
! 2.e Unresolved Obstacles Source Term
!
! 2.f Additional sources.
!
! 2.g Dump training data if necessary
!
! 3.  Set frequency cut-off ------------------------------------------ *
!
! 4.  Summation of source terms and diagonal term and time step ------ *
!
!
!     For input and dissipation calculate the fraction of the ice-free
!     surface. In the presence of ice, the effective water surface
!     is reduce to a fraction of the cell size free from ice, and so is
!     input :
!             SIN = (1-ICE)**ISCALEIN*SIN and SDS=(1-ICE)**ISCALEDS*SDS ------------------ *
!     INFLAGS2(4) is true if ice concentration was ever read during
!             this simulation
  
!GPUNotes The different between using the ACC WAIT and kernels is the
!knockon affect of gang and vector sizes since they must be consistant
!between kernel blocks.
            IF ( INFLAGS2(4) ) THEN
              VSNL(1:NSPECH,ISEA) = ICESCALENL*VSNL(1:NSPECH,ISEA)
              VDNL(1:NSPECH,ISEA) = ICESCALENL*VDNL(1:NSPECH,ISEA)
              VSLN(1:NSPECH,ISEA) = ICESCALELN*VSLN(1:NSPECH,ISEA)
              VSIN(1:NSPECH,ISEA) = ICESCALEIN*VSIN(1:NSPECH,ISEA)
              VDIN(1:NSPECH,ISEA) = ICESCALEIN*VDIN(1:NSPECH,ISEA)
              VSDS(1:NSPECH,ISEA) = ICESCALEDS*VSDS(1:NSPECH,ISEA)
              VDDS(1:NSPECH,ISEA) = ICESCALEDS*VDDS(1:NSPECH,ISEA)
            END IF
            NKI    = MAX( 2 , MIN ( NKH1 ,                           &
                     INT(FACTI2+FACTI1*LOG(MAX(1.E-7,FFXFI* FMEAN1)))))
!GPUNotes spectral loop up to frequency cut off
            VS(:,ISEA) = 0.
            VD(:,ISEA) = 0.
!$ACC LOOP SEQ
            DO IS=IS1, NSPECH
              VS(IS,ISEA)=VSLN(IS,ISEA) + VSIN(IS,ISEA) + VSNL(IS,ISEA)&
                     + VSDS(IS,ISEA) + VSBT(IS,ISEA)
              VD(IS,ISEA)=VDIN(IS,ISEA) + VDNL(IS,ISEA) + VDDS(IS,ISEA)&
                     + VDBT(IS,ISEA)
              DAMAX  = MIN(DAM(IS,ISEA), MAX(XREL*SPECINIT(IS,ISEA),   &
                                             AFILT))
              AFAC   = 1. / MAX(1.E-10 , ABS(VS(IS,ISEA)/DAMAX))
              DT     = MIN(DT, AFAC / (MAX ( 1.E-10, 1.+OFFSET*AFAC*   &
                                       MIN(0.,VD(IS,ISEA)) ) ) )
            END DO! end of loop on IS
!            WRITE(*,*) 'NODE_NUMBER', IX
!            IF (IX == DEBUG_NODE) WRITE(*,*) 'TIMINGS 1', DT
            DT     = MAX ( 0.5, DT ) 
! Here we have a hardlimit, which is not too usefull, at least not as a fixed con
!
            DTDYN(ISEA)  = DTDYN(ISEA) + DT
            IDT    = 1 + INT ( 0.99*(DTG-DTTOT)/DT ) ! number of iterations
            DT     = (DTG-DTTOT)/REAL(IDT)         ! actualy time step
            SHAVE  = DT.LT.DTMIN .AND. DT.LT.DTG-DTTOT ! limiter check ...
            SHAVEIO = SHAVE
            DT     = MAX ( DT , MIN (DTMIN,DTG-DTTOT) ) 
! override dt with input time step or last time step if it is bigger ... anyway
  
!CODENotes: This call is to a variable that is not defined, this causes
!issues on GPU but is just ignored on CPU. 

!        IF (srce_call .eq. srce_imp_post) DT = DTG! for implicit part
  
            HDT    = OFFSET * DT
            DTTOT  = DTTOT + DT
  
!GPUNotes calls below may be for implicit source term update
!GPUNotes would this remove the need for the NSTEPS loop?
!GPUNotes discussed in section 3.6 of manual, although maybe
!GPUNotes not it looks like Aron R code?
!GPUNotes Not used in source term test, which sets srce_direct
!GPUNotes loops below are over spectrum
  
! 5.  Increment spectrum --------------------------------------------- *
!
!GPUNotes Integrations over spectrum below active in source term test
            IF (srce_call .eq. srce_direct) THEN
!            SHAVE = .FALSE.
!           IF (IX == DEBUG_NODE) THEN
!              WRITE(*,'(A20,I20,F20.10,L20,4F20.10)') 'BEFORE', IX, DEPTH, SHAVE, SUM(VS), SUM(VD), SUM(SPEC)
!            ENDIF
              IF ( SHAVE ) THEN
!GPUNotes Running sequentially due to the recursive use of eIncX, adding
!loop independent does not force the loop to be parallel as it should.
  
!$ACC LOOP SEQ
                DO IS=IS1, NSPECH
                  eInc1 = VS(IS,ISEA)*DT/MAX(1.,(1.-HDT*VD(IS,ISEA)))
                  eInc2 = SIGN ( MIN (DAM(IS,ISEA),ABS(eInc1)) , eInc1 )
                  SPEC(IS,ISEA) = MAX ( 0. , SPEC(IS,ISEA)+eInc2 )
                END DO
              ELSE
!
!GPUNotes Running sequentially due to the recursive use of eIncX, adding
!loop independent does not force the loop to be parallel as it should.
!$ACC LOOP 
                DO IS=IS1, NSPECH
                  eInc1 = VS(IS,ISEA)*DT/MAX(1.,(1.-HDT*VD(IS,ISEA)))
                  SPEC(IS,ISEA) = MAX( 0. , SPEC(IS,ISEA)+eInc1 )
                END DO
              END IF
            END IF

! 5.b  Computes
!              atmos->wave flux PHIAW-------------------------------- *
!              wave ->BBL  flux PHIBBL------------------------------- *
!              wave ->ice  flux PHICE ------------------------------- *
!
            TMP1(3,ISEA)=0.
            HSTOT=0.
!GPUNotes Loops over spectrum - the spectrum must be properly updated first
  
!GPUNotes Loop has been refactored for ACC applications to allow the
!loops to be closely nested and collapsable.
  
!$ACC LOOP SEQ COLLAPSE(2)
            DO IK=IKS1, NK
              DO ITH=1, NTH
                FACTOR = DDEN(IK)/CG1(IK,ISEA)                  !Jacobian to get energy in band
                FACTOR2= FACTOR*GRAV*WN1(IK,ISEA)/SIG(IK)       ! coefficient to get momentum
   
            ! Wave direction is "direction to"
            ! therefore there is a PLUS sign for the stress
                IS   = (IK-1)*NTH + ITH
                COSI(1,ISEA)=ECOS(IS)
                COSI(2,ISEA)=ESIN(IS)
                PHIAW(ISEA) = PHIAW(ISEA) + VSIN(IS,ISEA) * DT * FACTOR&
                  / MAX ( 1. , (1.-HDT*VDIN(IS,ISEA))) ! semi-implict integration scheme
                PHIBBL(ISEA)= PHIBBL(ISEA)- VSBT(IS,ISEA) * DT * FACTOR&
                 / MAX ( 1. , (1.-HDT*VDBT(IS,ISEA))) ! semi-implict integration scheme
                PHINL = PHINL + VSNL(IS,ISEA) * DT * FACTOR&
                 / MAX ( 1. , (1.-HDT*VDNL(IS,ISEA))) ! semi-implict integration scheme
                IF (VSIN(IS,ISEA).GT.0.) THEN
                  TMP1(3,ISEA) = TMP1(3,ISEA) + SPEC(IS,ISEA)  * FACTOR
                END IF
                HSTOT = HSTOT + SPEC(IS,ISEA) * FACTOR
              END DO
            END DO
            TMP1(3,ISEA)=4.*SQRT(TMP1(3,ISEA))
            HSTOT=4.*SQRT(HSTOT)
            TAUWIX(ISEA)= TAUWIX(ISEA)+ TAUWX(ISEA) * DRAT *DT
            TAUWIY(ISEA)= TAUWIY(ISEA)+ TAUWY(ISEA) * DRAT *DT
            TAUWNX(ISEA)= TAUWNX(ISEA)+ TAUWAX * DRAT *DT
            TAUWNY(ISEA)= TAUWNY(ISEA)+ TAUWAY * DRAT *DT
        ! MISSING: TAIL TO BE ADDED ?
! 6.  Add tail ------------------------------------------------------- *
!   a Mean parameters
!
!GPUNotes source term specific loops over spectrum in this call
!        CALL CPU_TIME(sTime1)
            CALL W3SPR4 (SPEC(:,ISEA), CG1(:,ISEA), WN1(:,ISEA), EMEAN,&
                         FMEAN, FMEAN1, WNMEAN, AMAX, U10ABS(ISEA),    &
                         U10DIR(ISEA), USTAR(ISEA), USTDIR(ISEA),      &
                         TAUWX(ISEA), TAUWY(ISEA), CD, Z0, CHARN(ISEA),&
                         LLWS(:,ISEA), FMEANWS, ISEA)
!        CALL CPU_TIME(eTime1)
!        SPR4T = SPR4T + eTime1 - sTime1
!
! Introduces a Long & Resio (JGR2007) type dependance on wave age
            FAGE   = FFXFA*TANH(0.3*U10ABS(ISEA)*FMEANWS*TPI/GRAV)
            FH1    = (FFXFM+FAGE) * FMEAN1
   
            FH2    = FFXPM / USTAR(ISEA)
            FHIGH  = MIN ( SIG(NK) , MAX ( FH1 , FH2 ) )
            NKH    = MAX ( 2 , MIN ( NKH1 ,                           &
                   INT ( FACTI2 + FACTI1*LOG(MAX(1.E-7,FHIGH)) ) ) )
  
!GPUNotes UPDATE HOST solves issues with unnecessary last private values
!of scalars within the vector loop.
  
! 6.b Limiter for shallow water or Miche style criterion
!     Last time step ONLY !
!     uses true depth (D_INP) instead of limited depth
!
! 6.c Seeding of spectrum
!     alpha = 0.005 , 0.5 in eq., 0.25 for directional distribution
!
! 6.d Add tail
!
!GPUNotes Smaller spectral loop to add energy to tail
!GPUNotes Independence is acceptable for inner loop, outer loop
!dependence on IK within SPEC limits any collapsable option.
!$ACC LOOP SEQ COLLAPSE(2)
            DO IK=NKH+1, NK
              DO ITH=1, NTH
                SPEC(ITH+(IK-1)*NTH,ISEA) = SPEC(ITH+(IK-2)*NTH,ISEA) &
                * FACHFA         
              END DO
            END DO
!
! 6.e  Update wave-supported stress----------------------------------- *
!
!        CALL CPU_TIME(sTime2) 
            CALL W3SIN4(SPEC(:,ISEA), CG1(:,ISEA), WN2(:,ISEA),        &
                        U10ABS(ISEA), USTAR(ISEA), DRAT, AS(ISEA),     &
                        U10DIR(ISEA), Z0, CD, TAUWX(ISEA), TAUWY(ISEA),&
                        TAUWAX, TAUWAY, VSIN(:,ISEA),      &
                        VDIN(:,ISEA), LLWS(:,ISEA), IX, IY,            &
                        BRLAMBDA(:,ISEA) )
!        CALL CPU_TIME(eTime2)
!        SIN4T = SIN4T + eTime2 - sTime2
! 7.  Check if integration complete ---------------------------------- *
!
            IF (srce_call .eq. srce_imp_post) THEN
              EXIT
            ENDIF
            IF ( DTTOT .GE. 0.9999*DTG ) THEN
!            IF (IX == DEBUG_NODE) WRITE(*,*) 'DTTOT, DTG', DTTOT, DTG
              EXIT
            ENDIF
          END DO ! INTEGRATION LOOP
!
! ... End point dynamic integration - - - - - - - - - - - - - - - - - -
!
! 8.  Save integration data ------------------------------------------ *
          DTDYN(ISEA)  = DTDYN(ISEA) / REAL(MAX(1,NSTEPS))
          FCUT(ISEA)   = FHIGH * TPIINV
!
! Error escape locations
!
!
! 9.a  Computes PHIOC------------------------------------------ *
!     The wave to ocean flux is the difference between initial energy
!     and final energy, plus wind input plus the SNL flux to high freq.,
!     minus the energy lost to the bottom boundary layer (BBL)
!
!CODENotes: Brought initialisations here to avoid them being pulled and
!push to gpu.
  
!GPUNotes loop over spectrum requires spectrum to be properly updated
          EFINISH  = 0.
          MWXFINISH  = 0.
          MWYFINISH  = 0.
!$ACC LOOP REDUCTION(MWYFINISH, MWXFINISH, EFINISH)
          DO IK=1, NK
            EBAND = 0.
            A1BAND = 0.
            B1BAND = 0.
!$ACC LOOP REDUCTION(A1BAND, B1BAND, EBAND)
            DO ITH=1, NTH
              DIFF = SPECINIT(ITH+(IK-1)*NTH,ISEA)- &
                     SPEC(ITH+(IK-1)*NTH,ISEA)
              EBAND = EBAND + DIFF
              A1BAND = A1BAND + DIFF*ECOS(ITH)
              B1BAND = B1BAND + DIFF*ESIN(ITH)
            END DO
            EFINISH  = EFINISH  + EBAND * DDEN(IK) / CG1(IK,ISEA)
            MWXFINISH  = MWXFINISH  + A1BAND * DDEN(IK)                &
                         / CG1(IK,ISEA) * WN1(IK,ISEA)/SIG(IK)
            MWYFINISH  = MWYFINISH  + B1BAND * DDEN(IK)                &
                         / CG1(IK,ISEA) * WN1(IK,ISEA)/SIG(IK)
          END DO
!
! Transformation in momentum flux in m^2 / s^2
!
          TMP3(:,ISEA)=0
          TAUOX(ISEA)=(GRAV*MWXFINISH+TAUWIX(ISEA)-TMP3(1,ISEA))/DTG
          TAUOY(ISEA)=(GRAV*MWYFINISH+TAUWIY(ISEA)-TMP3(2,ISEA))/DTG
          TAUWIX(ISEA)=TAUWIX(ISEA)/DTG
          TAUWIY(ISEA)=TAUWIY(ISEA)/DTG
          TAUWNX(ISEA)=TAUWNX(ISEA)/DTG
          TAUWNY(ISEA)=TAUWNY(ISEA)/DTG
          TMP3(:,ISEA)=TMP3(:,ISEA)/DTG
!
! Transformation in wave energy flux in W/m^2=kg / s^3
!
          PHIOC(ISEA) =DWAT*GRAV*(EFINISH+PHIAW(ISEA)-PHIBBL(ISEA))/DTG
          PHIAW(ISEA) =DWAT*GRAV*PHIAW(ISEA) /DTG
          PHINL =DWAT*GRAV*PHINL /DTG
          PHIBBL(ISEA)=DWAT*GRAV*PHIBBL(ISEA)/DTG
!
! 10.1  Adds ice scattering and dissipation: implicit integration---------------- *
!     INFLAGS2(4) is true if ice concentration was ever read during
!             this simulation
!
   
          IF ( INFLAGS2(4).AND.ICE(ISEA).GT.0 ) THEN
            IF (IICEDISP) THEN
              ICECOEF2 = 1E-6
!             CALL LIU_FORWARD_DISPERSION (ICEH,ICECOEF2,DEPTH, &
!                                         SIG,WN_R,CG_ICE,ALPHA_LIU)
!
              IF (IICESMOOTH) THEN
              END IF
              WN_R(:,ISEA)=WN1(:,ISEA)
              CG_ICE(:,ISEA)=CG1(:,ISEA)
            END IF
!
            R(:,ISEA)=1 ! In case IC2 is defined but not IS2
            SPEC2(:,ISEA) = SPEC(:,ISEA)
!
            TMP4(:,ISEA) = 0.
            PHICE(ISEA) = 0.
  
!CODENotes: Previously was stored inside following IK loop, but due to
!specific values of ATT, is a null operation.
  
! First part of ice term integration: dissipation part
!
            ATT=1.
            SPEC(1+(IK-1)*NTH:NTH+(IK-1)*NTH,ISEA) =                   &
              ATT*SPEC2(1+(IK-1)*NTH:NTH+(IK-1)*NTH,ISEA)
!
  
!GPUNotes These spectral loops are for ice effects?
!GPUnotes slightly surprised they aren't enclosed in a conditional?
!GPUNotes actually maybe they are, think this is an indentation thing
!GPUNotes which I have tried to address
! Second part of ice term integration: scattering including re-distribution in directions
!
! 10.2  Fluxes of energy and momentum due to ice effects
!
!GPUNotes Moved FACTOR and FACTOR2 inside lopp for collapse.
!$ACC LOOP SEQ COLLAPSE(2)
            DO IK=1,NK
              DO ITH = 1,NTH
              FACTOR = DDEN(IK)/CG1(IK,ISEA)                  !Jacobian to get energy in band
              FACTOR2= FACTOR*GRAV*WN1(IK,ISEA)/SIG(IK)       ! coefficient to get momentum
                IS = ITH+(IK-1)*NTH
                PHICE(ISEA) = PHICE(ISEA) + (SPEC(IS,ISEA)-SPEC2(IS,   &
                              ISEA)) * FACTOR
                COSI(1,ISEA)=ECOS(IS)
                COSI(2,ISEA)=ESIN(IS)
                TMP4(:,ISEA) = TMP4(:,ISEA) -(SPEC(IS,ISEA)-SPEC2(IS,&
                ISEA))*FACTOR2*COSI(:,ISEA)
              END DO
            END DO
            PHICE(ISEA) =-1.*DWAT*GRAV*PHICE(ISEA) /DTG
            TMP4(:,ISEA)=TMP4(:,ISEA)/DTG
          ELSE
          END IF
  
! - - - - - - - - - - - - - - - - - - - - - -
! 11. Sea state dependent stress routine calls
! - - - - - - - - - - - - - - - - - - - - - -
!Note the Sea-state dependent stress calculations are primarily for high-wind
!conditions (>10 m/s).  It is not recommended to use these at lower wind
!in their current state.
!
   
! FLD1/2 requires the calculation of FPI:
!
! 12. includes shoreline reflection --------------------------------------------- *
!
!AR: this can be further simplified let's do some simple tests 1st ...
!
   
          FIRST  = .FALSE.
   
          IF (IT.EQ.0) SPEC(:,ISEA) = SPECINIT(:,ISEA)
          SPEC(:,ISEA) = MAX(0., SPEC(:,ISEA))
!LINE FROM W3WAVE
!          SIN4TOT = SIN4TOT + SIN4T
!          SPR4TOT = SPR4TOT + SPR4T
!          SDS4TOT = SDS4TOT + SDS4T

          WHITECAP(ISEA,:) = TMP1(:,ISEA)
          BEDFORMS(ISEA,:) = TMP2(:,ISEA)
          TAUBBL(ISEA,:) = TMP3(:,ISEA)
          TAUICE(ISEA,:) = TMP4(:,ISEA)
        ELSE
          USTAR   (ISEA) = UNDEF
          USTDIR(ISEA) = UNDEF
          DTDYN (ISEA) = UNDEF
          FCUT  (ISEA) = UNDEF
          SPEC(:,ISEA)  = 0.
        END IF
      END DO

!$ACC END KERNELS 
#ifdef MM
#else
!$ACC END DATA
#endif
      RETURN

! Formats
!
  101 FORMAT ('TIMESTAMP : ', A, F8.6)
!
 9006 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
              ' ------------- NEW DYNAMIC INTEGRATION LOOP',          &
              ' ------------- ')
!
 9062 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
              '               NKH        : ',I3)
!/
!/ End of W3SRCE ----------------------------------------------------- /
!/
      END SUBROUTINE W3SRCE
!/
!GPUNotes The subroutines below are for implicit options?
!GPUNotes so not used in the source term test
!/
!/ ------------------------------------------------------------------- /
      SUBROUTINE CALC_FPI( A, CG, FPI, S )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Jessica Meixner          |
!/                  |                                   |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-Jun-2018 |
!/                  +-----------------------------------+
!/
!/    06-Jul-2016 : Origination                         ( version 5.12 )
!/    06-Jul-2016 : Add SUBROUTINE SIGN_VSD_SEMI_IMPLICIT_WW3
!/                  Add optional DEBUGSRC/PDLIB           ( version 6.04 )
!/
!  1. Purpose :
!
!     Calculate equivalent peak frequency
!
!  2. Method :
!
!     Tolman and Chalikov (1996), equivalent peak frequency from source
 
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       CG      R.A.  I   Group velocities for k-axis of spectrum.
!       FPI     R.A.  O   Input 'peak' frequency.
!       S       R.A.  I   Source term (1-D version).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SRCE Subr.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S      Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, XFR, DDEN, SIG,FTE, FTTR
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), CG(NK), S(NSPEC)
      REAL, INTENT(OUT)       :: FPI
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IK
      REAL                    :: M0, M1
      REAL, ALLOCATABLE       :: SIN1A(:)
!/
!/ ------------------------------------------------------------------- /
!/
!
!     Calculate FPI: equivalent peak frequncy from wind source term
!     input
      ALLOCATE(SIN1A(NK))
!
      DO IK=1, NK
        SIN1A(IK) = 0.
        DO IS=(IK-1)*NTH+1, IK*NTH
          SIN1A(IK) = SIN1A(IK) + MAX ( 0. , S(IS) )
        END DO
      END DO
!
      M0     = 0.
      M1     = 0.
      DO IK=1, NK
        SIN1A(IK) = SIN1A(IK) * DDEN(IK) / ( CG(IK) * SIG(IK)**3 )
        M0        = M0 + SIN1A(IK)
        M1        = M1 + SIN1A(IK)/SIG(IK)
      END DO
!
      SIN1A(NK) = SIN1A(NK) / DDEN(NK)
      M0        = M0 + SIN1A(NK) * FTE
      M1        = M1 + SIN1A(NK) * FTTR
      IF ( M1 .LT. 1E-20 ) THEN
          FPI    = XFR * SIG(NK)
      ELSE
          FPI    = M0 / M1
      END IF
 
      END SUBROUTINE CALC_FPI
!/ ------------------------------------------------------------------- /!
      SUBROUTINE SIGN_VSD_SEMI_IMPLICIT_WW3(SPEC, VS, VD)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |                                   |
!/                  | Aron Roland (BGS IT&E GmbH)       |
!/                  | Mathieu Dutour-Sikiric (IRB)      |
!/                  |                                   |
!/                  |                        FORTRAN 90 |
!/                  | Last update :        01-June-2018 |
!/                  +-----------------------------------+
!/
!/    01-June-2018 : Origination.                        ( version 6.04 )
!/
!  1. Purpose : Put source term in matrix same as done always
!  2. Method :
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  6. Error messages :
!  7. Remarks
!  8. Structure :
!  9. Switches :
!
!   !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
        USE W3GDATMD, only : NTH, NK, NSPEC
        IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local PARAMETERs
!/
!/
!/ ------------------------------------------------------------------- /
!/
 
        INTEGER             :: ISP, ITH, IK, IS
        REAL, INTENT(IN)    :: SPEC(NSPEC)
        REAL, INTENT(INOUT) :: VS(NSPEC), VD(NSPEC)
        DO IS=1,NSPEC
          VD(IS) = MIN(0., VD(IS))
        END DO
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE SIGN_VSD_PATANKAR_WW3(SPEC, VS, VD)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |                                   |
!/                  | Aron Roland (BGS IT&E GmbH)       |
!/                  | Mathieu Dutour-Sikiric (IRB)      |
!/                  |                                   |
!/                  |                        FORTRAN 90 |
!/                  | Last update :        01-June-2018 |
!/                  +-----------------------------------+
!/
!/    01-June-2018 : Origination.                        ( version 6.04 )
!/
!  1. Purpose : Put source term in matrix Patankar style (experimental)
!  2. Method :
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  6. Error messages :
!  7. Remarks
!  8. Structure :
!  9. Switches :
!
!   !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!
 
        USE W3GDATMD, only : NTH, NK, NSPEC
        IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local PARAMETERs
!/
!/
!/ ------------------------------------------------------------------- /
!/
        INTEGER             :: ISP, ITH, IK, IS
        REAL, INTENT(IN)    :: SPEC(NSPEC)
        REAL, INTENT(INOUT) :: VS(NSPEC), VD(NSPEC)
        DO IS=1,NSPEC
          VD(IS) = MIN(0., VD(IS))
          VS(IS) = MAX(0., VS(IS))
        END DO
      END SUBROUTINE
!/
!/ End of module W3SRCEMD -------------------------------------------- /
!/
      END MODULE W3SRCEMD
