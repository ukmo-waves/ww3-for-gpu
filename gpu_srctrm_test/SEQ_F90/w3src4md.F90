#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRC4MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Nov-2013 |
!/                  +-----------------------------------+
!/
!/    30-Aug-2010 : Origination.                        ( version 3.14-Ifremer )
!/    02-Nov-2010 : Addding fudge factor for low freq.  ( version 4.03 )
!/    02-Sep-2011 : Clean up and time optimization      ( version 4.04 )
!/    04-Sep-2011 : Estimation of whitecap stats.       ( version 4.04 )
!/    13-Nov-2013 : Reduced frequency range with IG     ( version 4.13 )
!/
!  1. Purpose :
!
!     The 'SHOM/Ifremer' source terms based on P.A.E.M. Janssen's wind input
!     and dissipation functions by Ardhuin et al. (2009,2010)
!     and Filipot & Ardhuin (2010)
!     The wind input is converted from the original
!     WAM codes, courtesy of P.A.E.M. Janssen and J. Bidlot
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SPR4    Subr. Public   Mean parameters from spectrum.
!      W3SIN4    Subr. Public   WAM4+ input source term.
!      INSIN4    Subr. Public   Corresponding initialization routine.
!      TABU_STRESS, TABU_TAUHF, TABU_TAUHF2
!                Subr. Public   Populate various tables.
!      CALC_USTAR
!                Subr. Public   Compute stresses.
!      W3SDS4    Subr. Public   Dissipation (Ardhuin & al. / Filipot & Ardhuin)
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
!/ Public variables
!/
      !air kinematic viscosity (used in WAM)
      INTEGER, PARAMETER      :: ITAUMAX=200,JUMAX=200
      INTEGER, PARAMETER      :: IUSTAR=100,IALPHA=200, ILEVTAIL=50
      REAL                    :: DELTAUW, DELU
      ! Table for H.F. stress as a function of 2 variables
      REAL                    :: DELUST, DELALP
      ! Table for H.F. stress as a function of 3 variables
      REAL,ALLOCATABLE        :: TAUT(:,:), TAUHFT(:,:), TAUHFT2(:,:,:)
      ! Table for swell damping
      REAL                    :: DELTAIL
      REAL,    PARAMETER      :: UMAX    = 50.
      REAL,    PARAMETER      :: TAUWMAX = 2.2361 !SQRT(5.)
      INTEGER                 :: DIKCUMUL
!  Size of wave height table for integrating the PDF of wave heights
      INTEGER,    PARAMETER      :: NKHI=100
      REAL,    PARAMETER      :: FAC_KD1=1.01, FAC_KD2=1000., KHSMAX=2., KHMAX=2.
      REAL,    PARAMETER      ::KDMAX=200000.
!/
      ! Scratch variables for GPU
      INTEGER,ALLOCATABLE     :: NSMOOTH(:),IKSUP(:),IMSSMAX(:)
      REAL,ALLOCATABLE        :: S1(:), E1(:), COEF4(:)
      INTEGER,ALLOCATABLE     :: NTIMES(:)
      REAL,ALLOCATABLE        :: DK(:), HS(:), KBAR(:), DCK(:)
      REAL,ALLOCATABLE        :: EFDF(:)     ! Energy integrated over a spectral band
      REAL,ALLOCATABLE        :: BTH0(:)     !saturation spectrum
      REAL,ALLOCATABLE        :: BTH(:)   !saturation spectrum
      REAL,ALLOCATABLE        :: BTH0S(:)    !smoothed saturation spectrum
      REAL,ALLOCATABLE        :: BTHS(:)  !smoothed saturation spectrum
      REAL,ALLOCATABLE        :: SBK(:)
      REAL,ALLOCATABLE        :: SBKT(:), MSSSUM(:,:), WTHSUM(:)
      REAL,ALLOCATABLE        :: MSSSUM2(:,:)
      REAL,ALLOCATABLE        :: MSSLONG(:,:)
      REAL,ALLOCATABLE        :: QB(:), S2(:)
      REAL,ALLOCATABLE        :: PB(:),PB2(:)
      REAL, DIMENSION(:,:)   , ALLOCATABLE :: SIGTAB
      REAL, DIMENSION(:,:)   , ALLOCATABLE :: K1, K2
      REAL,ALLOCATABLE        :: EB(:), EB2(:), ALFA(:)
!/
      CONTAINS

!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SRC4_INIT()

      USE W3GDATMD, ONLY: NK, NTH, NSPEC, NKHS, NKD, NDTAB,QBI,DCKI

      IMPLICIT NONE

      WRITE(0,*) "W3SRC4_INIT: NK/NTH/NSPEC: ", NK, NTH, NSPEC
      ! allocate scratch variables for GPU (to avoid Syncing)
      ALLOCATE(NSMOOTH(NK),IKSUP(NK),S1(NK),E1(NK),COEF4(NK),NTIMES(NK)&
               ,DK(NK),HS(NK),KBAR(NK),DCK(NK),EFDF(NK),BTH0(NK),QB(NK)&
               ,S2(NK),BTH(NSPEC),BTH0S(NK),BTHS(NSPEC),SBK(NSPEC),    &
               IMSSMAX(NK),SBKT(NK),MSSSUM(NK,5),WTHSUM(NTH),PB(NSPEC),&
               MSSSUM2(NK,NTH),MSSLONG(NK,NTH),PB2(NSPEC),EB(NK),      &
               EB2(NK), ALFA(NK),K1(NK,NDTAB),K2(NK,NDTAB),            &
               SIGTAB(NK,NDTAB),DCKI(NKHS,NKD),QBI(NKHS,NKD))
!!$ACC KERNELS
!      NSMOOTH(:) = 0. 
!      IKSUP(:) = 0. 
!      S1(:) = 0. 
!      E1(:) = 0. 
!      COEF4(:) = 0. 
!      NTIMES(:) = 0. 
!      DK(:) = 0. 
!      HS(:) = 0. 
!      KBAR(:) = 0. 
!      DCK(:) = 0. 
!      EFDF(:) = 0. 
!      BTH0(:) = 0. 
!      QB(:) = 0. 
!      S2(:) = 0. 
!      BTH(:) = 0. 
!      BTH0S(:) = 0. 
!      BTHS(:) = 0. 
!      SBK(:) = 0. 
!      IMSSMAX(:) = 0. 
!      SBKT(:) = 0. 
!      MSSSUM(:,:) = 0. 
!      WTHSUM(:) = 0. 
!      PB(:) = 0. 
!      MSSSUM2(:,:) = 0. 
!      MSSLONG(:,:) = 0. 
!      PB2(:) = 0. 
!      EB(:) = 0. 
!      EB2(:) = 0. 
!      ALFA(:) = 0. 
!      K1(:,:) = 0. 
!      K2(:,:) = 0. 
!      SIGTAB(:,:) = 0. 
!      DCKI(:,:) = 0. 
!      QBI(:,:) = 0. 
!!$ACC END KERNELS      

      END SUBROUTINE W3SRC4_INIT
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SPR4 (A, CG, WN, EMEAN, FMEAN, FMEAN1, WNMEAN,     &
                    AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0,&
                    CHARN, LLWS, FMEANWS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         08-Jun-2018 |
!/                  +-----------------------------------+
!/
!/    03-Oct-2007 : Origination.                        ( version 3.13 )
!/    13-Jun-2011 : Adds f_m0,-1 as FMEAN in the outout ( version 4.04 )
!/    08-Jun-2018 : use STRACE and FLUSH                ( version 6.04 )
!/
!  1. Purpose :
!
!     Calculate mean wave parameters for the use in the source term
!     routines.
!
!  2. Method :
!
!     See source term routines.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum.
!       CG      R.A.  I   Group velocities.
!       WN      R.A.  I   Wavenumbers.
!       EMEAN   Real  O   Energy
!       FMEAN1  Real  O   Mean  frequency (fm0,-1) used for reflection
!       FMEAN   Real  O   Mean  frequency for determination of tail
!       WNMEAN  Real  O   Mean wavenumber.
!       AMAX    Real  O   Maximum of action spectrum.
!       U       Real  I   Wind speed.
!       UDIR    Real  I   Wind direction.
!       USTAR   Real I/O  Friction velocity.
!       USDIR   Real I/O  wind stress direction.
!       TAUWX-Y Real  I   Components of wave-supported stress.
!       CD      Real  O   Drag coefficient at wind level ZWND.
!       Z0      Real  O   Corresponding z0.
!       CHARN   Real  O   Corresponding Charnock coefficient
!       LLWS    L.A.  I   Wind sea true/false array for each component
!       FMEANWS Real  O   Mean frequency of wind sea, used for tail
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SRCE   Source term integration routine.
!       W3OUTP   Point output program.
!       GXEXPO   GrADS point output program.
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
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3ODATMD, ONLY: IAPROC
      USE CONSTANTS, ONLY: TPIINV
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DTH, DDEN, WWNMEANP, &
                          WWNMEANPTAIL, FTE, FTF, SSTXFTF, SSTXFTWN,&
                          SSTXFTFTAIL, SSWELLF
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK), U, UDIR
      REAL, INTENT(IN)        :: TAUWX, TAUWY
      LOGICAL, INTENT(IN)     :: LLWS(NSPEC)
      REAL, INTENT(INOUT)     :: USTAR ,USDIR
      REAL, INTENT(OUT)       :: EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX,  &
                                 CD, Z0, CHARN, FMEANWS
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IK, ITH
!/
      REAL                    :: TAUW, EBAND, EMEANWS, UNZ
!      REAL,ALLOCATABLE        :: EB(:), EB2(:), ALFA(:)
!      ALLOCATE(EB(NK), EB2(NK), ALFA(NK))
!/ ------------------------------------------------------------------- /
!/
      WRITE(0,*)'TAG: W3SPR4'
!!$ACC DATA CREATE(EB(:),EB2(:),ALFA(:))  
!!$ACC KERNELS 
      UNZ    = MAX ( 0.01 , U )
      USTAR  = MAX ( 0.0001 , USTAR )
!
      EMEAN  = 0.
      EMEANWS= 0.
      FMEANWS= 0.
      FMEAN  = 0.
      FMEAN1 = 0.
      WNMEAN = 0.
      AMAX   = 0.
!
! 1.  Integral over directions and maximum --------------------------- *
!
!GPUNotes spectral loop
!$ACC LOOP INDEPENDENT
      DO IK=1,NK
        EB(IK)  = 0.
        EB2(IK) = 0.
      END DO
!$ACC LOOP INDEPENDENT
      DO IK=1, NK
        DO ITH=1, NTH
          IS=ITH+(IK-1)*NTH
          EB(IK) = EB(IK) + A(ITH,IK)
          IF (LLWS(IS)) EB2(IK) = EB2(IK) + A(ITH,IK)
          AMAX   = MAX (AMAX, A(ITH,IK))
        END DO
      END DO
! 2.  Integrate over directions -------------------------------------- *
!
!GPUNotes directions loop only
!$ACC LOOP INDEPENDENT
      DO IK=1, NK
        ALFA(IK) = 2. * DTH * SIG(IK) * EB(IK) * WN(IK)**3
        EB(IK)   = EB(IK) * DDEN(IK) / CG(IK)
        EB2(IK)   = EB2(IK) * DDEN(IK) / CG(IK)
        EMEAN    = EMEAN  + EB(IK)
        FMEAN    = FMEAN  + EB(IK) /SIG(IK)
        FMEAN1   = FMEAN1 + EB(IK) *(SIG(IK)**(2.*WWNMEANPTAIL))
        WNMEAN   = WNMEAN + EB(IK) *(WN(IK)**WWNMEANP)
        EMEANWS  = EMEANWS+ EB2(IK)
        FMEANWS  = FMEANWS+ EB2(IK)*(SIG(IK)**(2.*WWNMEANPTAIL))
      END DO
!
! 3.  Add tail beyond discrete spectrum and get mean pars ------------ *
!     ( DTH * SIG absorbed in FTxx )
!
      EBAND  = EB(NK) / DDEN(NK)
      EMEAN  = EMEAN  + EBAND * FTE
      FMEAN  = FMEAN  + EBAND * FTF
      FMEAN1 = FMEAN1 + EBAND * SSTXFTFTAIL
      WNMEAN = WNMEAN + EBAND * SSTXFTWN
      EBAND  = EB2(NK) / DDEN(NK)
      EMEANWS = EMEANWS + EBAND * FTE
      FMEANWS = FMEANWS + EBAND * SSTXFTFTAIL
!
! 4.  Final processing
!
      FMEAN  = TPIINV * EMEAN / MAX ( 1.E-7 , FMEAN )
      IF (FMEAN1.LT.1.E-7) THEN
        FMEAN1=TPIINV * SIG(NK)
      ELSE
        FMEAN1  = TPIINV *( MAX ( 1.E-7 , FMEAN1 )                       &
                     / MAX ( 1.E-7 , EMEAN ))**(1/(2.*WWNMEANPTAIL))
      ENDIF
      WNMEAN = ( MAX ( 1.E-7 , WNMEAN )                              &
                / MAX ( 1.E-7 , EMEAN ) )**(1/WWNMEANP)
      IF (FMEANWS.LT.1.E-7.OR.EMEANWS.LT.1.E-7) THEN
        FMEANWS=TPIINV * SIG(NK)
      ELSE
        FMEANWS  = TPIINV *( MAX ( 1.E-7 , FMEANWS )                       &
                     / MAX ( 1.E-7 , EMEANWS ))**(1/(2.*WWNMEANPTAIL))
      END IF
 
! 5.  Cd and z0 ------------------------------------------------------ *
!
      TAUW = SQRT(TAUWX**2+TAUWY**2)
      Z0=0.
!!$ACC END KERNELS
      CALL CALC_USTAR(U,TAUW,USTAR,Z0,CHARN)
!!$ACC UPDATE HOST(USTAR)
!!$ACC KERNELS
      UNZ    = MAX ( 0.01 , U )
      CD     = (USTAR/UNZ)**2
      USDIR = UDIR
!!$ACC END KERNELS
!!$ACC END DATA
! 6.  Final test output ---------------------------------------------- *
!
      RETURN
!
! Formats
!
!/
!/ End of W3SPR4 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SPR4
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SIN4 (A, CG, K, U, USTAR, DRAT, AS, USDIR, Z0, CD,  &
                         TAUWX, TAUWY, TAUWNX, TAUWNY, S, D, LLWS,     &
                         IX, IY, BRLAMBDA)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Dec-2013 |
!/                  +-----------------------------------+
!/
!/    09-Oct-2007 : Origination.                        ( version 3.13 )
!/    24-Jan-2013 : Adding breaking-related input       ( version 4.16 )
!/    05-Dec-2013 : Cleaning up the ICE input           ( version 4.16 )
!/
!  1. Purpose :
!
!     Calculate diagonal and input source term for WAM4+ approach.
!
!  2. Method :
!
!       WAM-4     : Janssen et al.
!       WAM-"4.5" : gustiness effect (Cavaleri et al. )
!       SAT       : high-frequency input reduction for balance with
!                   saturation dissipation (Ardhuin et al., 2008)
!       SWELL     : negative wind input (Ardhuin et al. 2008)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       CG      R.A.  I   Group speed                              *)
!       K       R.A.  I   Wavenumber for entire spectrum.          *)
!       U       Real  I   WIND SPEED
!       USTAR   Real  I   Friction velocity.
!       DRAT    Real  I   Air/water density ratio.
!       AS      Real  I   Air-sea temperature difference
!       USDIR   Real  I   wind stress direction
!       Z0      Real  I   Air-side roughness lengh.
!       CD      Real  I   Wind drag coefficient.
!       USDIR   Real  I   Direction of friction velocity
!       TAUWX-Y Real  I   Components of the wave-supported stress.
!       TAUWNX  Real  I   Component of the negative wave-supported stress.
!       TAUWNY  Real  I   Component of the negative wave-supported stress.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative.             *)
!     ----------------------------------------------------------------
!                         *) Stored as 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!       STRACE    Subroutine tracing.                 ( !/S switch )
!       PRT2DS    Print plot of spectrum.             ( !/T0 switch )
!       OUTMAT    Print out matrix.                   ( !/T1 switch )
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output program.
!       GXEXPO   GrADS point output program.
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
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: GRAV,NU_AIR,KAPPA,TPI,FWTABLE,SIZEFWTABLE, &
                           DELAB,ABMIN
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, XFR, DDEN, SIG, SIG2, TH,   &
                          ESIN, ECOS, EC2, ZZWND, AALPHA, BBETA, ZZALP,&
                          TTAUWSHELTER, SSWELLF, DDEN2, DTH, SSINTHP,  &
                          ZZ0RAT, SSINBR
      USE W3ODATMD, ONLY: IAPROC, NDTO
      USE W3PARALL, ONLY: PRINT_MY_TIME
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), BRLAMBDA(NSPEC)
      REAL, INTENT(IN)        :: CG(NK), K(NSPEC),Z0,U, CD
      REAL, INTENT(IN)        :: USTAR, USDIR, AS, DRAT
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUWX, TAUWY, TAUWNX, TAUWNY
      LOGICAL, INTENT(OUT)    :: LLWS(NSPEC)
      INTEGER, INTENT(IN)     :: IX, IY
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS,IK,ITH, IR
      REAL                    :: FACLN1, FACLN2, LAMBDA
      REAL                    :: COSU, SINU, TAUX, TAUY, USDIRP, USTP
      REAL                    :: TAUPX, TAUPY, UST2, TAUW, TAUWB
      REAL, PARAMETER         :: EPS1 = 0.00001, EPS2 = 0.000001
      REAL                    :: Usigma           !standard deviation of U due to gustiness
      REAL                    :: USTARsigma       !standard deviation of USTAR due to gustiness
      REAL                    :: CM, UCN, ZCN, DSTAB, &
                                 Z0VISC, Z0NOZ, EB,  &
                                 EBX, EBY, AORB, AORB1, FW, UORB, TH2, &
                                 RE, FU, FUD, SWELLCOEFV, SWELLCOEFT
      REAL                    :: SMOOTH
      REAL                    :: XI,DELI1,DELI2
      REAL                    :: XJ,DELJ1,DELJ2
      REAL                    :: XK,DELK1,DELK2
      REAL                    :: CONST, CONST0, CONST2, TAU1
      REAL                    :: X,ZARG,ZLOG,UST
      REAL                    :: COSWIND, XSTRESS, YSTRESS, TAUHF
      REAL                    :: TEMP, TEMP2 
      INTEGER                 :: IND,J,I,ISTAB
      REAL                    :: DVISC, DTURB, PVISC, PTURB,ZERO 
      REAL                    :: STRESSSTABN1, STRESSSTABN2, &
                                 STRESSSTAB1, STRESSSTAB2
!/
!/ ------------------------------------------------------------------- /
!/
!
!      CALL PRINT_MY_TIME("    Calculate input source terms",NDTO)
! 1.  Preparations
!
      !JDM: Initializing values to zero, they shouldn't be used unless
      !set in another place, but seems to solve some bugs with certain
      !compilers.

! As local arrays they should not be required to be created. However
! with managed memory turned on and this statement removed the output
! fails.
!!$ACC DATA CREATE(STRESSSTAB1, STRESSSTAB2)
!!$ACC KERNELS
      PTURB = 0.
      PVISC = 0.
      D(:) = 0.
!
! 1.a  estimation of surface roughness parameters
!
      Z0VISC = 0.1*NU_AIR/MAX(USTAR,0.0001)
      Z0NOZ = MAX(Z0VISC,ZZ0RAT*Z0)
      FACLN1 = U / LOG(ZZWND/Z0NOZ)
      FACLN2 = LOG(Z0NOZ)
!
! 1.b  estimation of surface orbital velocity and displacement
!
      UORB=0.
      AORB=0.
 
!$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO IK=1, NK
        DO ITH=1, NTH
           IS=ITH+(IK-1)*NTH
           UORB = UORB + A(IS) *SIG(IK)**2 * DDEN(IK) / CG(IK)
           AORB = AORB + A(IS)             * DDEN(IK) / CG(IK) 
        END DO
      END DO
      UORB  = 2*UORB**0.5                  ! significant orbital amplitude
      AORB1 = 2*AORB**(1-0.5*SSWELLF(6))    ! half the significant wave height ... if SWELLF(6)=1
      RE = 4*UORB*AORB1 / NU_AIR           ! Reynolds number

!
! Defines the swell dissipation based on the "Reynolds number"
      
! GPUNotes Setting the loops to be more specific with ELSE IF removed
! issues surrounding PTURB and PVISC.

      IF (SSWELLF(4).GT.0) THEN
        IF (SSWELLF(7).GT.0.) THEN
          SMOOTH = 0.5*TANH((RE-SSWELLF(4))/SSWELLF(7))
          PTURB=(0.5+SMOOTH)
          PVISC=(0.5-SMOOTH)
        ELSE IF (SSWELLF(7).LE.0) THEN
          IF (RE.LE.SSWELLF(4)) THEN
            PTURB =  0.
            PVISC =  1.
          ELSE IF (RE.LE.SSWELLF(4)) THEN
            PTURB =  1.
            PVISC =  0.
          END IF
        END IF
      ELSE IF (SSWELLF(4).LE.0) THEN
        PTURB=1.
        PVISC=1.
      END IF

! GPUNotes Order of the conditions has been changed to ensure the 
! output is correct, this shouldn't be required and may be a sign
! of another issue. Can a more specific condition be used to identify
! positve entries for SSELLF(2)?

      IF (SSWELLF(2).NE.0.) THEN
        FU=ABS(SSWELLF(3))
        FUD=SSWELLF(2)
        AORB=2*AORB**0.5
        XI=(ALOG10(MAX(AORB/Z0NOZ,3.))-ABMIN)/DELAB
        IND  = MIN (SIZEFWTABLE-1, INT(XI))
        DELI1= MIN (1. ,XI-FLOAT(IND))
        DELI2= 1. - DELI1
        FW = FWTABLE(IND)*DELI2+FWTABLE(IND+1)*DELI1
      ELSE IF (SSWELLF(2).EQ.0.) THEN
        FW=ABS(SSWELLF(3))
        FU=0.
        FUD=0.
      END IF
!
! 2.  Diagonal
! Here AS is the air-sea temperature difference in degrees. Expression given by
! Abdalla & Cavaleri, JGR 2002 for Usigma. For USTARsigma ... I do not see where
! I got it from, maybe just made up from drag law ...
!
      UST=USTAR
      ISTAB=3
      TAUX = UST**2* COS(USDIR)
      TAUY = UST**2* SIN(USDIR)
!
! Loop over the resolved part of the spectrum
!
      ! ChrisB: STRESSSTAB[N] is a 2D array and does not reduce
      ! properly in an ACC loop. Only element 3 is ever using in
      ! the first dimension and the second dimension is size 2.
      ! Splitting into two seperate scalars simplifies the handling
      ! in ACC reduce operations.

      !STRESSSTAB(ISTAB,:)=0.
      !STRESSSTABN(ISTAB,:)=0.
      STRESSSTAB1 = 0.
      STRESSSTAB2 = 0.
      STRESSSTABN1 = 0.
      STRESSSTABN1 = 0.
!
! Coupling coefficient times density ratio DRAT
!
      CONST0=BBETA*DRAT/(KAPPA**2)
!
!GPUNotes loops over full spectrum
!         Private is not required, applied implicitly if not written
!         here. Is it possible to convert multiple variables into arrays
!         to run this loop in parallel?

!$ACC LOOP SEQ PRIVATE(STRESSSTAB1, STRESSSTAB2) 
      DO IK=1, NK
        TAUPX=TAUX-ABS(TTAUWSHELTER)*STRESSSTAB1
        TAUPY=TAUY-ABS(TTAUWSHELTER)*STRESSSTAB2
        USTP=MIN((TAUPX**2+TAUPY**2)**0.25,MAX(UST,0.3))
        USDIRP=ATAN2(TAUPY,TAUPX)
        COSU   = COS(USDIRP) ! CB - these lines are problematic as the LAST 
        SINU   = SIN(USDIRP) ! CB - value is used later outside the loop
!$ACC LOOP INDEPENDENT & 
!$ACC      REDUCTION(+:STRESSSTABN1,STRESSSTABN2,STRESSSTAB1,STRESSSTAB2)
        DO ITH=1,NTH
          IS=1+(IK-1)*NTH
          CM=K(IS)/SIG2(IS) !inverse of phase speed
          UCN=USTP*CM+ZZALP  !this is the inverse wave age
          ! the stress is the real stress (N/m^2) divided by
          CONST2=DDEN2(IS)/CG(IK) &        !Jacobian to get energy in band
              *GRAV/(SIG(IK)/K(IS)*DRAT) ! coefficient to get momentum
          CONST=SIG2(IS)*CONST0
          ZCN=ALOG(K(IS)*Z0)
          SWELLCOEFV=-SSWELLF(5)*DRAT*2*K(IS)*SQRT(2*NU_AIR*SIG2(IS))
          SWELLCOEFT=-DRAT*SSWELLF(1)*16*SIG2(IS)**2/GRAV

          IR=ITH+(IK-1)*NTH
          !WRITE(0,*)'IR: ', IR
          COSWIND=(ECOS(IR)*COSU+ESIN(IR)*SINU)
          IF (COSWIND.GT.0.01) THEN
            X=COSWIND*UCN
            ! this ZARG term is the argument of the exponential
            ! in Janssen 1991 eq. 16.
            ZARG=KAPPA/X
            ! ZLOG is ALOG(MU) where MU is defined by Janssen 1991 eq. 15
            ! MU=
            ZLOG=ZCN+ZARG
            IF (ZLOG.LT.0.) THEN
              ! The source term Sp is beta * omega * X**2
              ! as given by Janssen 1991 eq. 19
              DSTAB = CONST*EXP(ZLOG)*ZLOG**4*UCN*UCN*COSWIND**SSINTHP
              LLWS(IR)=.TRUE.
            ELSE
              DSTAB = 0.
              LLWS(IR)=.FALSE.
            END IF
!
!  Added for consistency with ECWAM implsch.F
!
            IF (28.*CM*USTAR*COSWIND.GE.1) THEN
              LLWS(IR)=.TRUE.
            END IF
          ELSE  ! (COSWIND.LE.0.01)
            DSTAB = 0.
            LLWS(IR)=.FALSE.
          END IF
!
          IF ((SSWELLF(1).NE.0.AND.DSTAB.LT.1E-7*SIG2(IR)) &
              .OR.SSWELLF(3).GT.0) THEN
            DVISC=SWELLCOEFV
            DTURB=SWELLCOEFT*(FW*UORB+(FU+FUD*COSWIND)*USTP)
            DSTAB = DSTAB + PTURB*DTURB +  PVISC*DVISC
! GPUNotes Only update D(:) at the end of the time step, partially removes 
! loop dependency. 

          END IF
! Sums up the wave-supported stress
!
          ! Wave direction is "direction to"
          ! therefore there is a PLUS sign for the stress
          !TEMP2=CONST2*DSTAB(ISTAB,IS)*A(IS)
          TEMP2=CONST2*DSTAB*A(IR)
          !IF (DSTAB(ISTAB,IS).LT.0) THEN
          IF (DSTAB.LT.0) THEN
            !STRESSSTABN(ISTAB,1)=STRESSSTABN(ISTAB,1)+TEMP2*ECOS(IS)
            !STRESSSTABN(ISTAB,2)=STRESSSTABN(ISTAB,2)+TEMP2*ESIN(IS)
            STRESSSTABN1=STRESSSTABN1+TEMP2*ECOS(IR)
            STRESSSTABN2=STRESSSTABN2+TEMP2*ESIN(IR)
          ELSE
            !STRESSSTAB(ISTAB,1)=STRESSSTAB(ISTAB,1)+TEMP2*ECOS(IS)
            !STRESSSTAB(ISTAB,2)=STRESSSTAB(ISTAB,2)+TEMP2*ESIN(IS)
            STRESSSTAB1=STRESSSTAB1+TEMP2*ECOS(IR)
            STRESSSTAB2=STRESSSTAB2+TEMP2*ESIN(IR)
          END IF
          D(IR) = DSTAB
        END DO

! GPUNotes Calls to STRESSSTAB1/2 inside the outer loops means that
! STRESSSTAB1/2 variable is local to the loop and not required outside.
! This prevents issues arrising with TAUPX as the variable is not used
! later on in the code, hence STRESSSTABN1/2 is unaffected.
      XSTRESS=STRESSSTAB1
      YSTRESS=STRESSSTAB2
      END DO

      TAUWNX =STRESSSTABN1
      TAUWNY =STRESSSTABN2
      !------------
      ! By adding the calls to the outer loop that is run sequentially
      ! we can avoid the need of redoing the calls as the last value is
      ! already being used. P.s. the outer loop must be seq anyway.
      !------------
      ! ChrisB: Need to repeat code from lines 548 - 554 here
      ! as COSU and SINU need to be the last calculated 
      ! values from the IK loop
!      TAUPX=TAUX-ABS(TTAUWSHELTER)*STRESSSTAB1
!      TAUPY=TAUY-ABS(TTAUWSHELTER)*STRESSSTAB2
!      USDIRP=ATAN2(TAUPY,TAUPX)
!      COSU   = COS(USDIRP) ! CB - these lines are problematic as the LAST 
!      SINU   = SIN(USDIRP) ! CB - value is used later outside the loop
      !------------

      S(:) = D(:) * A(:)
!
! ... Test output of arrays
!
      ! Computes the high-frequency contribution
      ! the difference in spectal density (kx,ky) to (f,theta)
      ! is integrated in this modified CONST0
      CONST0=DTH*SIG(NK)**5/((GRAV**2)*tpi) &
         *TPI*SIG(NK) / CG(NK)  !conversion WAM (E(f,theta) to WW3 A(k,theta)
      TEMP=0.
!GPUNotes loop over directions
!$ACC LOOP INDEPENDENT
      DO ITH=1,NTH
        IS=ITH+(NK-1)*NTH
        COSWIND=(ECOS(IS)*COSU+ESIN(IS)*SINU)
        TEMP=TEMP+A(IS)*(MAX(COSWIND,0.))**3
      END DO
 
      TAUPX=TAUX-ABS(TTAUWSHELTER)*XSTRESS
      TAUPY=TAUY-ABS(TTAUWSHELTER)*YSTRESS
      USTP=(TAUPX**2+TAUPY**2)**0.25
      USDIRP=ATAN2(TAUPY,TAUPX)
 
      UST=USTP
      ! finds the values in the tabulated stress TAUHFT
      XI=UST/DELUST
      IND  = MAX(1,MIN (IUSTAR-1, INT(XI)))
      DELI1= MAX(MIN (1. ,XI-FLOAT(IND)),0.)
      DELI2= 1. - DELI1
      XJ=MAX(0.,(GRAV*Z0/MAX(UST,0.00001)**2-AALPHA) / DELALP)
      J    = MAX(1 ,MIN (IALPHA-1, INT(XJ)))
      DELJ1= MAX(0.,MIN (1.      , XJ-FLOAT(J)))
      DELJ2=1. - DELJ1
      IF (TTAUWSHELTER.GT.0) THEN
        XK = CONST0*TEMP / DELTAIL
        I = MIN (ILEVTAIL-1, INT(XK))
        DELK1= MIN (1. ,XK-FLOAT(I))
        DELK2=1. - DELK1
        TAU1 =((TAUHFT2(IND,J,I)*DELI2+TAUHFT2(IND+1,J,I)*DELI1 )*DELJ2 &
               +(TAUHFT2(IND,J+1,I)*DELI2+TAUHFT2(IND+1,J+1,I)*DELI1)*DELJ1)*DELK2 &
              +((TAUHFT2(IND,J,I+1)*DELI2+TAUHFT2(IND+1,J,I+1)*DELI1 )*DELJ2 &
               +(TAUHFT2(IND,J+1,I+1)*DELI2+TAUHFT2(IND+1,J+1,I+1)*DELI1)*DELJ1)*DELK1
      ELSE
        TAU1 =(TAUHFT(IND,J)*DELI2+TAUHFT(IND+1,J)*DELI1 )*DELJ2 &
         +(TAUHFT(IND,J+1)*DELI2+TAUHFT(IND+1,J+1)*DELI1)*DELJ1
      END IF
      TAUHF = CONST0*TEMP*UST**2*TAU1
      TAUWX = XSTRESS+TAUHF*COS(USDIRP)
      TAUWY = YSTRESS+TAUHF*SIN(USDIRP)
!    
! Reduces tail effect to make sure that wave-supported stress
! is less than total stress, this is borrowed from ECWAM Stresso.F
!
      TAUW = (TAUWX**2+TAUWY**2)**0.5
      UST2   = MAX(USTAR,EPS2)**2
      TAUWB = MIN(TAUW,MAX(UST2-EPS1,EPS2**2))
      IF (TAUWB.LT.TAUW) THEN
        TAUWX=TAUWX*TAUWB/TAUW
        TAUWY=TAUWY*TAUWB/TAUW
      END IF
!!$ACC END KERNELS
      RETURN
!
! Formats
!
!/
!/ End of W3SIN4 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SIN4
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSIN4(FLTABS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |                         SHOM      |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Aug-2010 |
!/                  +-----------------------------------+
!/
!/    30-Aug-2010 : Origination.                        ( version 3.14-Ifremer )
!
!  1. Purpose :
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     ----------------------------------------------------------------
!      FLTABS    Logical
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
!      W3SIN4    Subr. W3SRC3MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: TPIINV, RADE, GRAV
      USE W3ODATMD,  ONLY: NDSE
      USE W3SERVMD,  ONLY: EXTCDE
      USE W3DISPMD,  ONLY: WAVNU2
      USE W3GDATMD,  ONLY: SIG, DSIP, NK, NTH, TTAUWSHELTER,             &
                           SSDSDTH, SSDSCOS, TH, DTH, XFR, ECOS, ESIN,   &
                           SSDSC,  SSDSBRF1, SSDSBCK, SSDSBINT, SSDSPBK, &
                           SSDSABK, SSDSHCK, IKTAB, DCKI, SATINDICES,    &
                           SATWEIGHTS, CUMULW, NKHS, NKD, NDTAB, QBI
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      LOGICAL, INTENT(IN)     :: FLTABS
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER  SDSNTH, ITH, I_INT, J_INT, IK, IK2, ITH2 , IS, IS2
      INTEGER  IKL, ID, ICON, IKD, IKHS, IKH, TOTO
      REAL     C, C2
      REAL     DIFF1, DIFF2, BINF, BSUP, CGG, PROF
      REAL     KIK, DHS, KD, KHS, KH, XT, GAM, DKH, PR, W, EPS
      REAL     DKD
!      ! Now declared in W3SRC4_INIT
!      REAL, DIMENSION(:,:)   , ALLOCATABLE :: SIGTAB
!      REAL, DIMENSION(:,:)   , ALLOCATABLE :: K1, K2
!        ALLOCATE(K1(NK,NDTAB))
!        ALLOCATE(K2(NK,NDTAB))
!        ALLOCATE(SIGTAB(NK,NDTAB))
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Initializations ------------------------------------------------ *
!
! These precomputed tables are written in mod_def.ww3
!
      IF (FLTABS) THEN
        CALL TABU_STRESS
        CALL TABU_TAUHF(SIG(NK) )      !tabulate high-frequency stress: 2D table
        IF (TTAUWSHELTER.GT.0) THEN
          CALL TABU_TAUHF2(SIG(NK) )   !tabulate high-frequency stress: 3D table
        END IF
      END IF
!
! 2.  SPONTANEOUS BREAKING
! 2.a Precomputes the indices for integrating the spectrum to get saturation (TEST 4xx )
!
      IF (SSDSDTH.LT.180) THEN
        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADE)),NTH/2-1)
        SATINDICES(:,:)=1
        SATWEIGHTS(:,:)=0.
!GPUNotes loops over full spectrum
        DO ITH=1,NTH
          DO I_INT=ITH-SDSNTH, ITH+SDSNTH
            J_INT=I_INT
            IF (I_INT.LT.1)  J_INT=I_INT+NTH
            IF (I_INT.GT.NTH) J_INT=I_INT-NTH
            SATINDICES(I_INT-(ITH-SDSNTH)+1,ITH)=J_INT
            SATWEIGHTS(I_INT-(ITH-SDSNTH)+1,ITH)=          &
                   COS(TH(ITH)-TH(J_INT))**SSDSCOS
          END DO
        END DO
      ELSE
        SATINDICES(:,:)=1
        SATWEIGHTS(:,:)=1.
      END IF
!/ ------------------------------------------------------------------- /
!
! Precomputes QBI and DCKI (TEST 500)
!
      IF (SSDSBCK.GT.0) THEN
!
! Precomputes the indices for integrating the spectrum over frequency bandwidth
!
        BINF=(1-SSDSBINT) ! Banner et al 2002: Hp=4*sqrt(int_0.7^1.3fp E df), SSDSBINT=0.3
        BSUP=(1+SSDSBINT)
        KIK=0.
!
! High frequency tail for convolution calculation
!
 
        SIGTAB=0. !contains frequency for upper windows boundaries
        IKTAB=0  ! contains indices for upper windows boundaries
 
!GPUNotes NDTAB is dimension for dissipation table, size=2000
!GPUNotes nested loops over frequency
        DO ID=1,NDTAB
          TOTO=0
          PROF=REAL(ID)
          DO IKL=1,NK ! last window starts at IK=NK
            CALL WAVNU2(SIG(IKL), PROF, KIK, CGG, 1E-7, 15, ICON)
            K1(IKL,ID)=KIK  ! wavenumber lower boundary (is directly related to the frequency indices, IK)
            K2(IKL,ID)=((BSUP/BINF)**2.)*K1(IKL,ID)! wavenumber upper boundary
            SIGTAB(IKL,ID)=SQRT(GRAV*K2(IKL,ID)*TANH(K2(IKL,ID)*ID)) ! corresponding frequency upper boundary
            IF(SIGTAB(IKL,ID) .LE. SIG(1)) THEN
              IKTAB(IKL,ID)=1
            END IF
            IF(SIGTAB(IKL,ID) .GT. SIG(NK)) THEN
              IKTAB(IKL,ID)=NK+TOTO       ! in w3sds4 only windows with IKSUP<=NK will be kept
              TOTO=1
            END IF
            DO IK=1,NK-1
              DIFF1=0.
              DIFF2=0.
              IF(SIG(IK)<SIGTAB(IKL,ID) .AND. SIG(IK+1)>=SIGTAB(IKL,ID)) THEN
                DIFF1=SIGTAB(IKL,ID)-SIG(IK)   ! seeks the indices of the upper boundary
                DIFF2=SIG(IK+1)-SIGTAB(IKL,ID)! the indices of lower boudary = IK
                IF (DIFF1<DIFF2) THEN
                  IKTAB(IKL,ID)=IK
                ELSE
                  IKTAB(IKL,ID)=IK+1
                END IF
              END IF
            END DO
          END DO
        END DO
!
! Tabulates DCKI and QBI
!
        DHS=KHSMAX/NKHS ! max value of KHS=KHSMAX
        DKH=KHMAX/NKHI  ! max value of KH=KHMAX
        DKD=KDMAX/NKD
        DCKI=0.
        QBI =0.
!GPUnotes Another nested loop over dissipation table and frequency dimensions
        DO IKD=1,NKD
          KHS=0.
          KD=(FAC_KD1**(IKD-FAC_KD2))
          XT=TANH(KD)
          GAM=1.0314*(XT**3)-1.9958*(XT**2)+1.5522*XT+0.1885
          GAM=GAM/2.15
          DO IKHS=1,NKHS  ! max value of KHS=1.
            KH=0.
            KHS=KHS+DHS
            DO IKH=1,NKHI
              KH=KH+DKH
              PR=(4.*KH/(KHS**2.))*exp(-(2*((KH/KHS)**2.)))
!              W=1.5*(((KHS)/(SQRT(2.)*GAM*XT))**2.)*(1-exp(-(((KH)/(GAM*XT))**4.))) !CK2002 parameterization
              W=SSDSABK*(((KHS)/(SQRT(2.)*GAM*XT))**2.)*(1-exp(-(((KH)/(GAM*XT))**SSDSPBK)))
              EPS=-((((SSDSBCK/(XT**SSDSHCK))*KH)**3.)/4)*SQRT(GRAV/XT)
              DCKI(IKHS, IKD)= DCKI(IKHS, IKD)+PR*W*EPS*DKH
              QBI(IKHS, IKD) = QBI(IKHS, IKD) +PR*W*    DKH
            END DO
          END DO
        END DO
 
        WHERE ( QBI .GT. 1. )
          QBI = 1.
        END WHERE
 
      ELSE
        IKTAB(:,:)=1
        DCKI(:,:) =0.
        QBI(:,:)  =0.
      END IF
!
!/ ------------------------------------------------------------------- /
!                        CUMULATIVE EFFECT
!/ ------------------------------------------------------------------- /
!
! Precomputes the weights for the cumulative effect (TEST 441 and 500)
!
      DIKCUMUL = 0
      IF (SSDSC(3).NE.0) THEN
!       DIKCUMUL is the integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
        DIKCUMUL = NINT(SSDSBRF1/(XFR-1.))
!        WRITE(6,*) 'INSIN4b:',DIKCUMUL
        CUMULW(:,:)=0.
!GPUNotes loops over full spectrum twice
        DO IK=1,NK
          C = GRAV/SIG(IK)   ! Valid in deep water only
          !C = SIG(IK)/K(IK) ! Valid in all water depth ???
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            DO IK2=1,IK-DIKCUMUL
              C2 = GRAV/SIG(IK2) ! Valid in deep water only
              !C2 = SIG(IK2)/K(IK2) ! Valid in all water depth ???
              DO ITH2=1,NTH
                IS2=ITH2+(IK2-1)*NTH
                CUMULW(IS2,IS)=SQRT(C**2+C2**2-2*C*C2*ECOS(1+ABS(ITH2-ITH))) & ! = deltaC
                                   *DSIP(IK2)/(0.5*C2) * DTH                   ! = dk*dtheta (Valid in deep water only)
              END DO
            END DO
          END DO
        END DO
      ELSE
        CUMULW(:,:)=0.
      END IF
!/
!/ End of INSIN4 ----------------------------------------------------- /
!/
      END SUBROUTINE INSIN4
! ----------------------------------------------------------------------
      SUBROUTINE TABU_STRESS
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    23-Jun-2006 : Origination.                        ( version 3.13 )
!/     adapted from WAM, original:P.A.E.M. JANSSEN    KNMI AUGUST 1990
!/     adapted version (subr. STRESS): J. BIDLOT    ECMWF OCTOBER 2004
!/     Table values were checkes against the original f90 result and found to
!/     be identical (at least at 0.001 m/s accuracy)
!/
!  1. Purpose :
!     TO GENERATE friction velocity table TAUT(TAUW,U10)=SQRT(TAU).
!     METHOD.
!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH
!                  Z1=Z0/SQRT(1-TAUW/TAU)
!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!
!     Initialization for source term routine.
!
!  2. Method :
!
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
!      W3SIN3    Subr. W3SRC3MD Corresponding source term.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: KAPPA, GRAV
      USE W3GDATMD, ONLY: ZZWND, AALPHA, ZZ0MAX
      IMPLICIT NONE
      INTEGER, PARAMETER      :: NITER=10
      REAL   , PARAMETER      :: XM=0.50, EPS1=0.00001
!     VARIABLE.   TYPE.     PURPOSE.
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *XNU*       REAL      KINEMATIC VISCOSITY OF AIR.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
!      *EPS1*      REAL      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
!                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
! ----------------------------------------------------------------------
      INTEGER I,J,ITER
      REAL ZTAUW,UTOP,CDRAG,WCD,USTOLD,TAUOLD
      REAL X,UST,ZZ0,ZNU,F,DELF,ZZ00
      ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX))
      DELU    = UMAX/FLOAT(JUMAX)
      DELTAUW = TAUWMAX/FLOAT(ITAUMAX)
!GPUNotes loops over stress table dimensions and calculation iteration
      DO I=0,ITAUMAX
        ZTAUW   = (REAL(I)*DELTAUW)**2
        DO J=0,JUMAX
           UTOP    = FLOAT(J)*DELU
           CDRAG   = 0.0012875
           WCD     = SQRT(CDRAG)
           USTOLD  = UTOP*WCD
           TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
           DO ITER=1,NITER
             X   = ZTAUW/TAUOLD
             UST = SQRT(TAUOLD)
             ZZ00=AALPHA*TAUOLD/GRAV
             IF (ZZ0MAX.NE.0) ZZ00=MIN(ZZ00,ZZ0MAX)
             ! Corrects roughness ZZ00 for quasi-linear effect
             ZZ0 = ZZ00/(1.-X)**XM
             !ZNU = 0.1*nu_air/UST  ! This was removed by Bidlot in 1996
             !ZZ0 = MAX(ZNU,ZZ0)
             F   = UST-KAPPA*UTOP/(ALOG(ZZWND/ZZ0))
             DELF= 1.-KAPPA*UTOP/(ALOG(ZZWND/ZZ0))**2*2./UST &
                      *(1.-(XM+1)*X)/(1.-X)
             UST = UST-F/DELF
             TAUOLD= MAX(UST**2., ZTAUW+EPS1)
          END DO
          TAUT(I,J)  = SQRT(TAUOLD)
        END DO
      END DO
      I=ITAUMAX
      J=JUMAX
!
!  Force zero wind to have zero stress (Bidlot 1996)
!
!GPUNotes loops over single dimension of stress table
      DO I=0,ITAUMAX
        TAUT(I,0)=0.0
      END DO
      RETURN
      END SUBROUTINE TABU_STRESS
!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF(SIGMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    27-Feb-2004 : Origination in WW3                  ( version 2.22.SHOM )
!/     the resulting table was checked to be identical to the original f77 result
!/    14-Aug-2006 : Modified following Bidlot           ( version 2.22.SHOM )
!/    18-Aug-2006 : Ported to version 3.09
!
!  1. Purpose :
!
!     Tabulation of the high-frequency wave-supported stress
!
!  2. Method :
!
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!     See tech. Memo ECMWF 03 december 2003 by Bidlot & Janssen
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SIGMAX   Real  I   maximum frequency * TPI
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
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
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: KAPPA, GRAV
      USE W3GDATMD, ONLY: AALPHA, BBETA, ZZALP, XFR, FACHFE, ZZ0MAX
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, intent(in) :: SIGMAX  !  maximum frequency
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!       USTARM  R.A.  Maximum friction velocity
!       ALPHAM  R.A.  Maximum Charnock Coefficient
!       WLV     R.A.  Water levels.
!       UA      R.A.  Absolute wind speeds.
!       UD      R.A.  Absolute wind direction.
!       U10     R.A.  Wind speed used.
!       U10D    R.A.  Wind direction used.
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      REAL                    :: USTARM, ALPHAM
      REAL                    :: CONST1, OMEGA, OMEGAC
      REAL                    :: UST, ZZ0,OMEGACC, CM
      INTEGER, PARAMETER      :: JTOT=250
      REAL, ALLOCATABLE       :: W(:)
      REAL                    :: ZX,ZARG,ZMU,ZLOG,ZZ00,ZBETA
      REAL                    :: Y,YC,DELY
      INTEGER                 :: J,K,L
      REAL                    :: X0
!
      ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA))
      USTARM = 5.
      ALPHAM = 20.*AALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)
      CONST1 = BBETA/KAPPA**2
      OMEGAC = SIGMAX
!
      TAUHFT(0:IUSTAR,0:IALPHA)=0. !table initialization
!
      ALLOCATE(W(JTOT))
      W(2:JTOT-1)=1.
      W(1)=0.5
      W(JTOT)=0.5
      X0 = 0.05
!
!GPUNotes loops over stress table dimensions
      DO L=0,IALPHA
        DO K=0,IUSTAR
          UST      = MAX(REAL(K)*DELUST,0.000001)
          ZZ00       = UST**2*AALPHA/GRAV
          IF (ZZ0MAX.NE.0) ZZ00=MIN(ZZ00,ZZ0MAX)
          ZZ0       = ZZ00*(1+FLOAT(L)*DELALP/AALPHA)
          OMEGACC  = MAX(OMEGAC,X0*GRAV/UST)
          YC       = OMEGACC*SQRT(ZZ0/GRAV)
          DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          ! For a given value of UST and ALPHA,
          ! the wave-supported stress is integrated all the way
          ! to 0.05*g/UST
          DO J=1,JTOT
            Y        = YC+REAL(J-1)*DELY
            OMEGA    = Y*SQRT(GRAV/ZZ0)
            ! This is the deep water phase speed
            CM       = GRAV/OMEGA
            !this is the inverse wave age, shifted by ZZALP (tuning)
            ZX       = UST/CM +ZZALP
            ZARG     = MIN(KAPPA/ZX,20.)
            ZMU      = MIN(GRAV*ZZ0/CM**2*EXP(ZARG),1.)
            ZLOG     = MIN(ALOG(ZMU),0.)
            ZBETA        = CONST1*ZMU*ZLOG**4
            ! Power of Y in denominator should be FACHFE-4
            TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
          END DO
        END DO
      END DO
      DEALLOCATE(W)
      RETURN
      END SUBROUTINE TABU_TAUHF
 
!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF2(SIGMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  | Last update 2013/01/24            |
!/                  +-----------------------------------+
!/
!/    15-May-2007 : Origination in WW3                  ( version 3.10.SHOM )
!/    24-Jan-2013 : Allows to read in table             ( version 4.08 )
!
!  1. Purpose :
!
!     Tabulation of the high-frequency wave-supported stress as a function of
!     ustar, alpha (modified Charnock), and tail energy level
!
!  2. Method :
!
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!     See tech. Memo ECMWF 03 december 2003 by Bidlot & Janssen
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SIGMAX   Real  I   maximum frequency*TPI
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
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
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: KAPPA, GRAV
      USE W3GDATMD, ONLY: AALPHA, BBETA, ZZALP, XFR, FACHFE,  &
                          TTAUWSHELTER, ZZ0MAX
      USE W3ODATMD, ONLY: NDSE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, intent(in) :: SIGMAX  !  maximum frequency * TPI
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!       USTARM  R.A.  Maximum friction velocity
!       ALPHAM  R.A.  Maximum Charnock Coefficient
!       WLV     R.A.  Water levels.
!       UA      R.A.  Absolute wind speeds.
!       UD      R.A.  Absolute wind direction.
!       U10     R.A.  Wind speed used.
!       U10D    R.A.  Wind direction used.
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      REAL                    :: USTARM, ALPHAM, LEVTAILM
      REAL                    :: CONST1, OMEGA, OMEGAC, LEVTAIL
      REAL                    :: UST, UST0, ZZ0,OMEGACC, CM
      REAL                    :: TAUW, TAUW0
      INTEGER, PARAMETER      :: JTOT=250
      REAL, ALLOCATABLE       :: W(:)
      REAL                    :: ZX,ZARG,ZMU,ZLOG,ZBETA
      REAL                    :: Y,YC,DELY
      INTEGER                 :: I, J, K, L
      REAL                    :: X0, INSIGMAX, INAALPHA, INBBETA, INZZALP, INKAPPA, INGRAV
      INTEGER                 :: INIUSTAR, INIALPHA, INILEVTAIL, IERR
      CHARACTER(160)          :: FNAMETAB
      LOGICAL                 :: NOFILE
      CHARACTER(LEN=10), PARAMETER :: VERGRD = '2018-06-08'
      CHARACTER(LEN=35), PARAMETER :: IDSTR = 'WAVEWATCH III ST4 TABLE FOR STRESS '
      CHARACTER(LEN=10)       :: VERTST=' '
      CHARACTER(LEN=35)       :: IDTST=' '
!
      ALLOCATE(TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL))
      FNAMETAB='ST4TABUHF2.bin'
      NOFILE=.TRUE.
      OPEN (993,FILE=FNAMETAB,FORM='UNFORMATTED',IOSTAT=IERR,STATUS='OLD')
      IF (IERR.EQ.0) THEN
        READ(993,IOSTAT=IERR) IDTST, VERTST, INSIGMAX, INAALPHA, INBBETA, INIUSTAR,  &
                              INIALPHA, INILEVTAIL, INZZALP, INKAPPA, INGRAV
        IF (VERTST.EQ.VERGRD.AND.IDTST.EQ.IDSTR.AND.IERR.EQ.0             &
            .AND.INSIGMAX.EQ.SIGMAX.AND.INAALPHA.EQ.AALPHA.AND.INBBETA.EQ.BBETA) THEN
          IF (INIUSTAR.EQ.IUSTAR.AND.INIALPHA.EQ.IALPHA.AND.INILEVTAIL.EQ.ILEVTAIL.AND. &
              INZZALP.EQ.ZZALP.AND.INGRAV.EQ.GRAV.AND.INKAPPA.EQ.KAPPA) THEN
            NOFILE=.FALSE.
          ELSE
            CLOSE(993)
          END IF
        END IF
      END IF
!
      USTARM = 5.
      ALPHAM = 20.*AALPHA
      LEVTAILM = 0.05
      DELUST  = USTARM/REAL(IUSTAR)
      DELALP  = ALPHAM/REAL(IALPHA)
      DELTAIL = ALPHAM/REAL(ILEVTAIL)
      CONST1  = BBETA/KAPPA**2
      OMEGAC  = SIGMAX
800   CONTINUE
      IF ( NOFILE ) THEN
        WRITE(NDSE,*) 'Filling 3D look-up table for SIN4. please wait'
        WRITE(NDSE,*)  IDSTR, VERGRD, SIGMAX, AALPHA, BBETA, IUSTAR, IALPHA,  &
                       ILEVTAIL, ZZALP, KAPPA, GRAV
!
        TAUHFT(0:IUSTAR,0:IALPHA)=0.  !table initialization
!
        ALLOCATE(W(JTOT))
        W(2:JTOT-1)=1.
        W(1)=0.5
        W(JTOT)=0.5
        X0 = 0.05
!
!GPUNotes Loops over stress table dimensions
        DO K=0,IUSTAR
          UST0      = MAX(REAL(K)*DELUST,0.000001)
          DO L=0,IALPHA
            UST=UST0
            ZZ0       = UST0**2*(AALPHA+FLOAT(L)*DELALP)/GRAV
            OMEGACC  = MAX(OMEGAC,X0*GRAV/UST)
            YC       = OMEGACC*SQRT(ZZ0/GRAV)
            DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          ! For a given value of UST and ALPHA,
          ! the wave-supported stress is integrated all the way
          ! to 0.05*g/UST
            DO I=0,ILEVTAIL
              LEVTAIL=REAL(I)*DELTAIL
              TAUHFT(K,L)=0.
              TAUHFT2(K,L,I)=0.
              TAUW0=UST0**2
              TAUW=TAUW0
              DO J=1,JTOT
                Y        = YC+REAL(J-1)*DELY
                OMEGA    = Y*SQRT(GRAV/ZZ0)
                ! This is the deep water phase speed
                CM       = GRAV/OMEGA
                !this is the inverse wave age, shifted by ZZALP (tuning)
                ZX       = UST0/CM +ZZALP
                ZARG     = MIN(KAPPA/ZX,20.)
                ZMU      = MIN(GRAV*ZZ0/CM**2*EXP(ZARG),1.)
                ZLOG     = MIN(ALOG(ZMU),0.)
                ZBETA        = CONST1*ZMU*ZLOG**4
                ! Power of Y in denominator should be FACHFE-4
                TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
                ZX       = UST/CM +ZZALP
                ZARG     = MIN(KAPPA/ZX,20.)
                ZMU      = MIN(GRAV*ZZ0/CM**2*EXP(ZARG),1.)
                ZLOG     = MIN(ALOG(ZMU),0.)
                ZBETA        = CONST1*ZMU*ZLOG**4
                ! Power of Y in denominator should be FACHFE-4
                TAUHFT2(K,L,I)  = TAUHFT2(K,L,I)+W(J)*ZBETA*(UST/UST0)**2/Y*DELY
                TAUW=TAUW-W(J)*UST**2*ZBETA*LEVTAIL/Y*DELY
                UST=SQRT(MAX(TAUW,0.))
              END DO
            END DO
          END DO
        END DO
        DEALLOCATE(W)
        OPEN (993,FILE=FNAMETAB,FORM='UNFORMATTED',IOSTAT=IERR,STATUS='UNKNOWN')
        WRITE(993) IDSTR, VERGRD, SIGMAX, AALPHA, BBETA, IUSTAR, IALPHA, ILEVTAIL, ZZALP, KAPPA, GRAV
        WRITE(993) TAUHFT(0:IUSTAR,0:IALPHA)
        WRITE(993) TAUHFT2
        CLOSE(993)
        !DO K=0,IUSTAR
        !  DO L=0,IALPHA
        !    DO I=0,ILEVTAIL
        !      WRITE(995,*) K,L,I,MAX(REAL(K)*DELUST,0.000001),AALPHA+FLOAT(L)*DELALP,REAL(I)*DELTAIL,TAUHFT(K,L),TAUHFT2(K,L,I)
        !      END DO
        !    END DO
        !  END DO
!
      ELSE
        WRITE(NDSE,*) 'Reading 3D look-up table for SIN4 from file.'
        READ(993,ERR=2000,IOSTAT=IERR ) TAUHFT(0:IUSTAR,0:IALPHA)
        READ(993,ERR=2000,IOSTAT=IERR ) TAUHFT2
        CLOSE(993)
        END IF
!
      GOTO 2001
2000  NOFILE=.TRUE.
      GOTO 800
2001  CONTINUE
      RETURN
      END SUBROUTINE TABU_TAUHF2
 
!/ ------------------------------------------------------------------- /
      SUBROUTINE CALC_USTAR(WINDSPEED,TAUW,USTAR,Z0,CHARN)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    27-Feb-2004 : Origination in WW3                  ( version 2.22-SHOM )
!/     the resulting table was checked to be identical to the original f77 result
!/    14-Aug-2006 : Modified following Bidlot           ( version 2.22-SHOM )
!/    18-Aug-2006 : Ported to version 3.09
!/    03-Apr-2010 : Adding output of Charnock parameter ( version 3.14-IFREMER )
!
!  1. Purpose :
!
!     Compute friction velocity based on wind speed U10
!
!  2. Method :
!
!     Computation of u* based on Quasi-linear theory
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       U10,TAUW,USTAR,Z0
!     ----------------------------------------------------------------
!       WINDSPEED Real  I   10-m wind speed ... should be NEUTRAL
!       TAUW      Real  I   Wave-supported stress
!       USTAR     Real  O   Friction velocity.
!       Z0        Real  O   air-side roughness length
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
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
!       !/S      Enable subroutine tracing.
!       !/T      Enable test output.
!
! 10. Source code :
!-----------------------------------------------------------------------------!
      USE CONSTANTS, ONLY: GRAV, KAPPA
      USE W3GDATMD,  ONLY: AALPHA, ZZWND
      IMPLICIT NONE
      REAL, intent(in) :: WINDSPEED,TAUW
      REAL, intent(out) :: USTAR, Z0, CHARN
      ! local variables
      REAL SQRTCDM1
      REAL XI,DELI1,DELI2,XJ,delj1,delj2
      REAL TAUW_LOCAL

      INTEGER IND,J
!!$ACC DATA COPYOUT(USTAR)
!!$ACC KERNELS
      TAUW_LOCAL=MAX(MIN(TAUW,TAUWMAX),0.)
      XI      = SQRT(TAUW_LOCAL)/DELTAUW
      IND     = MIN ( ITAUMAX-1, INT(XI)) ! index for stress table
      DELI1   = MIN(1.,XI - REAL(IND))  !interpolation coefficient for stress table
      DELI2   = 1. - DELI1
      XJ      = WINDSPEED/DELU
      J       = MIN ( JUMAX-1, INT(XJ) )
      DELJ1   = MIN(1.,XJ - REAL(J))
      DELJ2   = 1. - DELJ1
      USTAR=(TAUT(IND,J)*DELI2+TAUT(IND+1,J  )*DELI1)*DELJ2 &
       + (TAUT(IND,J+1)*DELI2+TAUT(IND+1,J+1)*DELI1)*DELJ1

! Determines roughness length
!
      SQRTCDM1  = MIN(WINDSPEED/USTAR,100.0)
      Z0  = ZZWND*EXP(-KAPPA*SQRTCDM1)
      IF (USTAR.GT.0.001) THEN
        CHARN = GRAV*Z0/USTAR**2
      ELSE
        CHARN = AALPHA
      END IF
!
!!$ACC END KERNELS
!!$ACC END DATA
      RETURN
      END SUBROUTINE CALC_USTAR
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SDS4 (A, K, CG, USTAR, USDIR, DEPTH, SRHS,      &
                         DDIAG, IX, IY, BRLAMBDA, WHITECAP )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  !            F. Ardhuin             !
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-Jun-2018 |
!/                  +-----------------------------------+
!/
!/    30-Aug-2010 : Clean up from common ST3-ST4 routine( version 3.14-Ifremer )
!/    23-Jan-2012 : Add output of lambdas to be used in SIN
!/    13-Nov-2013 : Reduced frequency range with IG1 switch
!/    06-Jun-2018 : Add optional DEBUGSRC              ( version 6.04 )
!/
!  1. Purpose :
!
!     Calculate whitecapping source term and diagonal term of derivative.
!
!  2. Method :
!
!       This codes does either one or the other of
!       Ardhuin et al. (JPO 2010)
!       Filipot & Ardhuin (JGR 2012)
!       the choice depends on SDSBCK
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IX, IY    Int   I   Grid Index
!       A         R.A.  I   Action density spectrum (1-D).
!       K         R.A.  I   Wavenumber for entire spectrum.          *)
!       USTAR     Real  I   Friction velocity.
!       USDIR     Real  I   wind stress direction.
!       DEPTH     Real  I   Water depth.
!       S         R.A.  O   Source term (1-D version).
!       D         R.A.  O   Diagonal term of derivative.             *)
!       BRLAMBDA  R.A.  O   Phillips' Lambdas
!     ----------------------------------------------------------------
!                         *) Stored in 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!       STRACE    Subroutine tracing.                 ( !/S switch )
!       PRT2DS    Print plot of spectrum.             ( !/T0 switch )
!       OUTMAT    Print out matrix.                   ( !/T1 switch )
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output program.
!       GXEXPO   GrADS point output program.
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
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS,ONLY: GRAV, DAIR, DWAT, PI, TPI, RADE, DEBUG_NODE
      USE W3GDATMD, ONLY: NSPEC, NTH, NK, SSDSBR, DDEN,              &
                          SSDSC, EC2, ES2, ESC,                      &
                          SIG, SSDSP, ECOS, ESIN, DTH, DSIP,         &
                          SSDSISO, SSDSDTH, SSDSBR2, SSDSBM,         &
                          SSDSBRFDF, SSDSBCK, SSDSBINT, IKTAB, DCKI, &
                          SATINDICES, SATWEIGHTS, CUMULW, NKHS, NKD, &
                          NDTAB, QBI
      USE W3ODATMD, ONLY: UNDEF, FLOGRD
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, OPTIONAL, INTENT(IN) :: IX, IY
      REAL, INTENT(INOUT)     :: A(NSPEC), K(NK), CG(NK)
      REAL, INTENT(IN)        :: DEPTH, USTAR, USDIR
      REAL, INTENT(OUT)       :: SRHS(NSPEC), DDIAG(NSPEC), BRLAMBDA(NSPEC)
      REAL, INTENT(OUT)       :: WHITECAP(1:4)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IS2, IS0, IKL, IKC, ID, NKL, ISPEC
      INTEGER                 :: IK, IK1, ITH, IK2, JTH, ITH2,         &
                                 IKHS, IKD, SDSNTH, IT, IKM, NKM
      REAL                    :: COSWIND, ASUM, SDIAGISO
      REAL                    :: COEF1, COEF2, COEF3
      REAL                    :: FACTURB, DTURB, BREAKFRACTION
      REAL                    :: RENEWALFREQ, EPSR
      REAL                    :: GAM, XT
      REAL                    :: FACSAT, DKHS, FACSTRAIN
      INTEGER                 :: NTHSUM
      REAL                    :: MSSPCS,MSSPC2,MSSPS2,MSSP,MSSD,MSSTH
      REAL                    :: MICHE, X, FACHF
      REAL                    :: TSTR, TMAX, DT, T, MFT

!      ! Now declared in W3SRC4_INIT
!      INTEGER,ALLOCATABLE     :: NSMOOTH(:),IKSUP(:),IMSSMAX(:)
!      REAL,ALLOCATABLE        :: S1(:), E1(:), COEF4(:)
!      INTEGER,ALLOCATABLE     :: NTIMES(:)
!      REAL,ALLOCATABLE        :: DK(:), HS(:), KBAR(:), DCK(:)
!      REAL,ALLOCATABLE        :: EFDF(:)     ! Energy integrated over a spectral band
!      REAL,ALLOCATABLE        :: BTH0(:)     !saturation spectrum
!      REAL,ALLOCATABLE        :: BTH(:)   !saturation spectrum
!      REAL,ALLOCATABLE        :: BTH0S(:)    !smoothed saturation spectrum
!      REAL,ALLOCATABLE        :: BTHS(:)  !smoothed saturation spectrum
!      REAL,ALLOCATABLE        :: SBK(:)
!      REAL,ALLOCATABLE        :: SBKT(:), MSSSUM(:,:), WTHSUM(:)
!      REAL,ALLOCATABLE        :: MSSSUM2(:,:)
!      REAL,ALLOCATABLE        :: MSSLONG(:,:)
!      REAL,ALLOCATABLE        :: QB(:), S2(:)
!      REAL,ALLOCATABLE        :: PB(:),PB2(:)
!/
!/ ------------------------------------------------------------------- /
!/
!      ! Now allocated in W3SRC4_INIT
!      ALLOCATE(NSMOOTH(NK),IKSUP(NK),S1(NK),E1(NK),COEF4(NK),NTIMES(NK)&
!               ,DK(NK),HS(NK),KBAR(NK),DCK(NK),EFDF(NK),BTH0(NK),QB(NK)&
!               ,S2(NK),BTH(NSPEC),BTH0S(NK),BTHS(NSPEC),SBK(NSPEC),    &
!               !IMSSMAX(NK),SBKT(NK),MSSSUM(NK,5),WTHSUM(NTH),PB(NSPEC),&
!               IMSSMAX(NK),SBKT(NK),MSSSUM(NK,5),PB(NSPEC),&
!               MSSSUM2(NK,NTH),MSSLONG(NK,NTH),PB2(NSPEC))

!!$ACC DATA CREATE(NSMOOTH(:),IKSUP(:),S1(:),E1(:),COEF4(:),NTIMES(:)  )&
!!$ACC      CREATE(DK(:),HS(:),KBAR(:),DCK(:),EFDF(:),BTH0(:),QB(:)    )&
!!$ACC      CREATE(S2(:),BTH0S(:),BTHS(:),SBK(:),IMSSMAX(:),WTHSUM(:)  )&
!!$ACC      CREATE(S2(:),BTH0S(:),BTHS(:),SBK(:),IMSSMAX(:) )&
!!$ACC      CREATE(SBKT(:),MSSSUM(:,:),PB(:),PB2(:),MSSLONG(:,:)       )&
!!$ACC      CREATE(MSSSUM2(:,:), BTH(:))
!
!----------------------------------------------------------------------
!
! 0.  Pre-Initialization to zero out arrays. All arrays should be reset
!     within the computation, but these are helping with some bugs
!     found in certain compilers

!!$ACC KERNELS
!CODENotes: Removed pre-initialization, this creates additional 
!data transfers and are not needed for the mini-app to run.
       
! 1.  Initialization and numerical factors
!
      FACTURB=SSDSC(5)*USTAR**2/GRAV*DAIR/DWAT
      BREAKFRACTION=0.
      RENEWALFREQ=0.
      IK1=1
 
      !IF (IX == DEBUG_NODE) WRITE(*,'(A20,4F20.10)') 'ST4 DISSIP ANFANG', SUM(SRHS), SUM(DDIAG)
      NTHSUM=MIN(FLOOR(SSDSC(10)+0.5),NTH-1)  ! number of angular bins for enhanced modulation
      IF (NTHSUM.GT.0) THEN
        WTHSUM(1:NTHSUM)=1
        WTHSUM(NTHSUM+1)=SSDSC(10)+0.5-NTHSUM
      ELSE
        WTHSUM(1)=2*SSDSC(10)
      END IF
!!$ACC END KERNELS
!
! 2.   Estimation of spontaneous breaking
!GPUNotes Attempted to use ACC wait but the code continues to break
!using this, requires further research as to why. 

!!$ACC KERNELS
!
      IF ( (SSDSBCK-SSDSC(1)).LE.0 ) THEN
!
! 2.a  Case of a direction-dependent breaking term (TEST441)
!
        EPSR = SQRT(SSDSBR)
!
! 2.a.1 Computes saturation
!
        SDSNTH = MIN(NINT(SSDSDTH/(DTH*RADE)),NTH/2-1)
!$ACC LOOP INDEPENDENT
        DO IK=1,NK
          DO ITH=1,NTH
            MSSLONG(IK,ITH) = 0.
            MSSSUM2(IK,ITH) = 0.
            MSSSUM(IK,ITH) = 0. 
!       SSDSDIK is the integer difference in frequency bands
!       between the "large breakers" and short "wiped-out waves"
          END DO
          BTH(IK) = 0.
        END DO

!GPUNotes: Loops over full spectrum
!$ACC LOOP INDEPENDENT
        DO  IK=IK1, NK
 
          FACSAT=SIG(IK)*K(IK)**3*DTH
          IS0=(IK-1)*NTH
          BTH(IS0+1)=0.
          ASUM = SUM(A(IS0+1:IS0+NTH))
          BTH0(IK)=ASUM*FACSAT
 
          IF (SSDSDTH.GE.180) THEN  ! integrates around full circle
            BTH(IS0+1:IS0+NTH)=BTH0(IK)
          ELSE
!
! straining effect: first finds the mean direction of mss, then applies cos^2
!                   straining
!
            IF (SSDSC(8).GT.0.OR.SSDSC(11).GT.0) THEN
              IKC = MAX(1,IK-DIKCUMUL)
              IMSSMAX (IK) = 1
              MSSP   = 0.
              MSSPC2 = 0.
              MSSPS2 = 0.
              MSSPCS = 0.
!
! Sums the contributions to the directional MSS for all ITH
!
!GPUNotes: loops over directions and sum iteration within IK loop above
!CODENotes: Internalise as much code as possible for collapsable loops
!$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO ITH=1,NTH
                DO JTH=-NTHSUM,NTHSUM
                IS=ITH+(IK-1)*NTH
                MSSLONG(IK,ITH) = K(IK)**2 * A(IS) * DDEN(IK) / CG(IK) ! contribution to MSS
                   ITH2 = 1+MOD(ITH-1+JTH+NTH,NTH)
                   MSSSUM2(IK:NK,ITH2) = MSSSUM2(IK:NK,ITH2)+MSSLONG(IK,ITH)*WTHSUM(ABS(JTH)+1)
                END DO
                MSSPC2 = MSSPC2 +MSSLONG(IK,ITH)*EC2(ITH)
                MSSPS2 = MSSPS2 +MSSLONG(IK,ITH)*ES2(ITH)
                MSSPCS = MSSPCS +MSSLONG(IK,ITH)*ESC(ITH)
                MSSP   = MSSP +MSSLONG(IK,ITH)
              END DO
!
! Now sums over IK
!
              MSSSUM  (IK:NK,1) = MSSSUM (IK:NK,1) +MSSP
              MSSSUM  (IK:NK,3) = MSSSUM (IK:NK,3) +MSSPC2
              MSSSUM  (IK:NK,4) = MSSSUM (IK:NK,4) +MSSPS2
              MSSSUM  (IK:NK,5) = MSSSUM (IK:NK,5) +MSSPCS
!
! Direction of long wave mss summed up to IK
!
              MSSD=0.5*(ATAN2(2*MSSSUM(IK,5),MSSSUM(IK,3)-MSSSUM(IK,4)))
              IF (MSSD.LT.0) MSSD = MSSD + PI
              IMSSMAX (IK)=1+NINT(MSSD *NTH/TPI)
!
! mss along perpendicular direction
!
              MSSSUM  (IK,2)  = MAX(0.,MSSSUM(IK,4)*COS(MSSD)**2   &
                       -2*MSSSUM(IK,5)*SIN(MSSD)*COS(MSSD)+MSSSUM(IK,3)*SIN(MSSD)**2)
 
            END IF ! SSDSC(8).GT.0) THEN
 
!GPUNotes: loop over directions
!$ACC LOOP INDEPENDENT PRIVATE(BTH)
            DO ITH=1,NTH            ! partial integration
              IS=ITH+(IK-1)*NTH
!
! Testing straining effect of long waves on short waves
! extended from Longuet-Higgins and Stewart (JFM 1960, eq. 2.27) the amplitude modulation
! in deep water is equal to the long wave slope k*a cos(theta1-theta2)
! Here we assume that the saturation is modulated as (1 + 2*SSDSC(8) *  sqrt(mss_theta) )
! where mss_theta is the mss in direction ITH.
!
! Note: SSDSC(8) is 2 times the MTF: close to 4x2=8 for very small ka but can be 20 or larger for ak~0.1
!
              IF (SSDSC(8).GT.0) THEN
!
                MSSTH=(MSSSUM(IKC,1)-MSSSUM(IKC,2))*EC2(1+ABS(ITH-IMSSMAX (IKC))) &
                        +MSSSUM(IKC,2)*ES2(1+ABS(ITH-IMSSMAX (IKC)))
!
                FACSTRAIN=1+SSDSC(8)*SQRT(MSSTH)+SSDSC(11)*SQRT(MSSSUM2(IKC,ITH))
                FACSAT=SIG(IK)*K(IK)**3*DTH*FACSTRAIN
              END IF
 
              BTH(IS)=DOT_PRODUCT(SATWEIGHTS(:,ITH),         &
                       A(IS0+SATINDICES(:,ITH)) )*FACSAT
!              BTH(IS)=SUM(  SATWEIGHTS(:,ITH)*         &
!                     A(IS0+SATINDICES(:,ITH)) )*FACSAT
            END DO
 
            IF (SSDSISO.EQ.1) THEN
              BTH0(IK)=SUM(A(IS0+1:IS0+NTH))*FACSAT
            ELSE
              BTH0(IK)=MAXVAL(BTH(IS0+1:IS0+NTH))
            END IF
          END IF
 
        END DO !NK END
!
!   Optional smoothing of B and B0 over frequencies
!
        IF (SSDSBRFDF.GT.0.AND.SSDSBRFDF.LT.NK/2) THEN
!$ACC LOOP INDEPENDENT
          DO IK=1,NK
            BTH0S(IK)=BTH0(IK)
            BTHS(IK)=BTH(IK)
            NSMOOTH(IK)=1
          END DO
!GPUNotes: loops over full spectrum - assume need to be sequential
!CODENotes: Internalise as much code as possible for collapsable loops
!original indentation has been left for code taken from outer loops
!$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO IK=1, SSDSBRFDF
            DO ITH=1,NTH
            BTH0S(1+SSDSBRFDF)=BTH0S(1+SSDSBRFDF)+BTH0(IK)
            NSMOOTH(1+SSDSBRFDF)=NSMOOTH(1+SSDSBRFDF)+1
              IS=ITH+(IK-1)*NTH
              BTHS(ITH+SSDSBRFDF*NTH)=BTHS(ITH+SSDSBRFDF*NTH)+BTH(IS)
            END DO
          END DO
!$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO IK=IK1+1+SSDSBRFDF,1+2*SSDSBRFDF
            DO ITH=1,NTH
            BTH0S(1+SSDSBRFDF)=BTH0S(1+SSDSBRFDF)+BTH0(IK)
            NSMOOTH(1+SSDSBRFDF)=NSMOOTH(1+SSDSBRFDF)+1
              IS=ITH+(IK-1)*NTH
              BTHS(ITH+SSDSBRFDF*NTH)=BTHS(ITH+SSDSBRFDF*NTH)+BTH(IS)
            END DO
          END DO
!$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO IK=SSDSBRFDF,IK1,-1
            DO ITH=1,NTH
            BTH0S(IK)=BTH0S(IK+1)-BTH0(IK+SSDSBRFDF+1)
            NSMOOTH(IK)=NSMOOTH(IK+1)-1
              IS=ITH+(IK-1)*NTH
              BTHS(IS)=BTHS(IS+NTH)-BTH(IS+(SSDSBRFDF+1)*NTH)
            END DO
          END DO
!$ACC LOOP INDEPENDENT COLLAPSE(2)  
          DO IK=IK1+1+SSDSBRFDF,NK-SSDSBRFDF
            DO ITH=1,NTH
            BTH0S(IK)=BTH0S(IK-1)-BTH0(IK-SSDSBRFDF-1)+BTH0(IK+SSDSBRFDF)
            NSMOOTH(IK)=NSMOOTH(IK-1)
              IS=ITH+(IK-1)*NTH
              BTHS(IS)=BTHS(IS-NTH)-BTH(IS-(SSDSBRFDF+1)*NTH)+BTH(IS+(SSDSBRFDF)*NTH)
            END DO
          END DO
!$ACC LOOP INDEPENDENT COLLAPSE(2)
          DO IK=NK-SSDSBRFDF+1,NK
            DO ITH=1,NTH
            BTH0S(IK)=BTH0S(IK-1)-BTH0(IK-SSDSBRFDF)
            NSMOOTH(IK)=NSMOOTH(IK-1)-1
              IS=ITH+(IK-1)*NTH
              BTHS(IS)=BTHS(IS-NTH)-BTH(IS-(SSDSBRFDF+1)*NTH)
            END DO
          END DO
!
!    final division by NSMOOTH
!
!GPUNotes created loop instead of array syntax.
!$ACC LOOP INDEPENDENT
          DO IK=1,NK
            BTH0(IK)=MAX(0.,BTH0S(IK)/NSMOOTH(IK))
          END DO

!$ACC LOOP INDEPENDENT
          DO IK=IK1,NK
            IS0=(IK-1)*NTH
            BTH(IS0+1:IS0+NTH)=MAX(0.,BTHS(IS0+1:IS0+NTH)/NSMOOTH(IK))
          END DO
        END IF
!
!  2.a.2  Computes spontaneous breaking dissipation rate
!
!GPUNotes Loop over ferquencies
!$ACC LOOP INDEPENDENT
        DO  IK=IK1, NK
!
!  Correction of saturation level for shallow-water kinematics
!
          IF (SSDSBM(0).EQ.1) THEN
            MICHE=1.
          ELSE
            X=TANH(MIN(K(IK)*DEPTH,10.))
            MICHE=(X*(SSDSBM(1)+X*(SSDSBM(2)+X*(SSDSBM(3)+X*SSDSBM(4)))))**2 ! Correction of saturation level for shallow-water kine
          END IF
          COEF1=(SSDSBR*MICHE)
!
!  Computes isotropic part
!
          SDIAGISO = SSDSC(2) * SIG(IK)*SSDSC(6)*(MAX(0.,BTH0(IK)/COEF1-1.))**2
!
!  Computes anisotropic part and sums isotropic part
!
          COEF2=SSDSC(2) * SIG(IK)*(1-SSDSC(6))/(COEF1*COEF1)
          COEF3=-2.*SIG(IK)*K(IK)*FACTURB
          DDIAG((IK-1)*NTH+1:IK*NTH) = SDIAGISO + &
                                   COEF2*((MAX(0.,BTH((IK-1)*NTH+1:IK*NTH)-COEF1))**SSDSP)
!            IF (IX == DEBUG_NODE) THEN
!              WRITE(*,'(A10,I10,10F15.6)') 'ST4 D3',IK,BTH0(IK),SUM(BTH((IK-1)*NTH+1:IK*NTH)),COEF1,COEF2,COEF3,SSDSP,SDIAGISO
!            ENDIF
        END DO

          !IF(IX == DEBUG_NODE) WRITE(*,'(A20,4F15.6)') 'ST4 DISSIP 1', SUM(SRHS), SUM(DDIAG), SUM(BTH)
 
!
! Computes Breaking probability
!
!GPUNotes created loop instead of array syntax.
!$ACC LOOP INDEPENDENT
        DO ISPEC=1,NSPEC
          PB(ISPEC) = (MAX(SQRT(BTH(ISPEC))-EPSR,0.))**2
!
! Multiplies by 28.16 = 22.0 * 1.6² * 1/2 with
!  22.0 (Banner & al. 2000, figure 6)
!  1.6  the coefficient that transforms  SQRT(B) to Banner et al. (2000)'s epsilon
!  1/2  factor to correct overestimation of Banner et al. (2000)'s breaking probability due to zero-crossing analysis
!
          PB(ISPEC) = PB(ISPEC) * 28.16
        END DO
!/
      END IF ! End of test for (Ardhuin et al. 2010)'s spontaneous dissipation source term
!
! 2.b             Computes spontaneous breaking for T500 //////////////
!
!GPUNotes Not sure that the loops in 2b are activated for the source
!GPUNotes term test
      IF (SSDSBCK.GT.0) THEN ! test for (Filipot et al. 2010)'s disspation source term
        E1 = 0.
        HS = 0.
        SRHS  = 0.
        DDIAG = 0.
        PB2  = 0.
!
! Computes Wavenumber spectrum E1 integrated over direction and computes dk
!
!GPUNotes loops over full spectrum
        E1(:)=0.
!$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO IK=IK1, NK
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            E1(IK)=E1(IK)+(A(IS)*SIG(IK))*DTH
          END DO
          DK(IK)=DDEN(IK)/(DTH*SIG(IK)*CG(IK))
        END DO
!
! Gets windows indices of IKTAB
!
        ID=MIN(NINT(DEPTH),NDTAB)
        IF (ID < 1) THEN
          ID = 1
        ELSE IF(ID > NDTAB) THEN
          ID = NDTAB
        END IF
!
! loop over wave scales
!
!GPUNotes loop over frequencies
        HS=0.
        KBAR=0.
        EFDF=0.
        NKL=0. !number of windows
!$ACC LOOP INDEPENDENT
        DO IKL=1,NK
          IKSUP(IKL)=IKTAB(IKL,ID)
          IF (IKSUP(IKL) .LE. NK) THEN
            EFDF(IKL) = DOT_PRODUCT(E1(IKL:IKSUP(IKL)-1),DK(IKL:IKSUP(IKL)-1))
            IF (EFDF(IKL) .NE. 0) THEN
              KBAR(IKL) = DOT_PRODUCT(K(IKL:IKSUP(IKL)-1)*E1(IKL:IKSUP(IKL)-1), &
                                      DK(IKL:IKSUP(IKL)-1)) / EFDF(IKL)
            ELSE
              KBAR(IKL)=0.
            END IF
! estimation of Significant wave height of a given scale
            HS(IKL) = 4*SQRT(EFDF(IKL))
            NKL = NKL+1
          END IF
        END DO
!
! Computes Dissipation and breaking probability in each scale
!
        DCK=0.
        QB =0.
        DKHS = KHSMAX/NKHS
!GPUNotes loop over ferquencies
!$ACC LOOP INDEPENDENT
        DO IKL=1, NKL
          IF (HS(IKL) .NE. 0. .AND. KBAR(IKL) .NE. 0.)  THEN
!
! gets indices for tabulated dissipation DCKI and breaking probability QBI
!
            IKD = FAC_KD2+ANINT(LOG(KBAR(IKL)*DEPTH)/LOG(FAC_KD1))
            IKHS= 1+ANINT(KBAR(IKL)*HS(IKL)/DKHS)
            IF (IKD > NKD) THEN    ! Deep water
              IKD = NKD
            ELSE IF (IKD < 1) THEN ! Shallow water
              IKD = 1
            END IF
            IF (IKHS > NKHS) THEN
              IKHS = NKHS
            ELSE IF (IKHS < 1) THEN
              IKHS = 1
            END IF
            XT = TANH(KBAR(IKL)*DEPTH)
!
!  Gamma corrected for water depth
!
            GAM=1.0314*(XT**3)-1.9958*(XT**2)+1.5522*XT+0.1885
!
! Computes the energy dissipated for the scale IKL
! using DCKI which is tabulated in INSIN4
!
            DCK(IKL)=((KBAR(IKL)**(-2.5))*(KBAR(IKL)/(2*PI)))*DCKI(IKHS,IKD)
!
! Get the breaking probability for the scale IKL
!
            QB(IKL) = QBI(IKHS,IKD) ! QBI is tabulated in INSIN4
          ELSE
            DCK(IKL)=0.
            QB(IKL) =0.
          END IF
        END DO
!
! Distributes scale dissipation over the frequency spectrum
!
!$ACC LOOP INDEPENDENT
        DO IK=1,NK
          S1(IK) = 0.
          S2(IK) = 0.
          NTIMES(IK) = 0
        END DO
!GPUNotes loop over freqeuncies
!$ACC LOOP INDEPENDENT
        DO IKL=1, NKL
          IF (EFDF(IKL) .GT. 0.) THEN
            S1(IKL:IKSUP(IKL))    = S1(IKL:IKSUP(IKL)) + &
                                     DCK(IKL)*E1(IKL:IKSUP(IKL)) / EFDF(IKL)
            S2(IKL:IKSUP(IKL))    = S2(IKL:IKSUP(IKL)) + &
                                     QB(IKL) *E1(IKL:IKSUP(IKL)) / EFDF(IKL)
            NTIMES(IKL:IKSUP(IKL)) = NTIMES(IKL:IKSUP(IKL)) + 1
          END IF
        END DO
!
! Finish the average
!
!$ACC LOOP INDEPENDENT
        DO IK=1,NK
          IF (NTIMES(IK) .GT. 0) THEN
            S1(IK) = S1(IK) / NTIMES(IK)
            S2(IK) = S2(IK) / NTIMES(IK)
          ELSE
            S1(IK) = 0.
            S2(IK) = 0.
          END IF
        END DO
! goes back to action for dissipation source term

!GPUNotes created loop instead of array syntax division.
!$ACC LOOP INDEPENDENT
        DO IK=1,NK
          S1(IK) = S1(IK) / SIG(IK)
        END DO
!
! Makes Isotropic distribution
!
        ASUM = 0.
!GPUNotes loop over frequencies
!$ACC LOOP INDEPENDENT
        DO IK = 1, NK
          ASUM = (SUM(A(((IK-1)*NTH+1):(IK*NTH)))*DTH)
          IF (ASUM.GT.1.E-8) THEN
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) DDIAG(IS)  = S1(IK)/ASUM
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) PB2(IS) = S2(IK)/ASUM
          ELSE
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) DDIAG(IS)  = 0.
            FORALL (IS=1+(IK-1)*NTH:IK*NTH) PB2(IS) = 0.
          END IF
          IF (PB2(1+(IK-1)*NTH).GT.0.001) THEN
            BTH0(IK) = 2.*SSDSBR
          ELSE
            BTH0(IK) = 0.
          END IF
        END DO
!
!GPUNotes created loop instead of array syntax.
!$ACC LOOP INDEPENDENT
        DO ISPEC=1,NSPEC
          PB(ISPEC) = (1-SSDSC(1))*PB2(ISPEC)*A(ISPEC) + &
                      SSDSC(1)*PB(ISPEC)
        END DO
!
      END IF   ! END OF TEST ON SSDSBCK
!
! 3.   Computes Lambda from breaking probability
!
! Compute Lambda = PB* l(k,th)
! with l(k,th)=1/(2*pi²)= the breaking crest density
!
!GPUNotes created loop instead of array syntax.
!$ACC LOOP INDEPENDENT
      DO ISPEC=1,NSPEC
        BRLAMBDA(ISPEC) = PB(ISPEC) / (2.*PI**2.)
      END DO
!     
!/ ------------------------------------------------------------------- /
!             WAVE-TURBULENCE INTERACTION AND CUMULATIVE EFFECT
!/ ------------------------------------------------------------------- /
!
! loop over spectrum

!
!GPUNotes Loops over the full spectrum
      SBKT(:)=0.
      SBK(:)=0.
!$ACC LOOP INDEPENDENT GANG COLLAPSE(2)
      DO  IK=IK1, NK
        DO ITH=1,NTH
          IS=ITH+(IK-1)*NTH
!
! Computes cumulative effect from Breaking probability
!
          RENEWALFREQ = 0.
          IF (SSDSC(3).NE.0 .AND. IK.GT.DIKCUMUL) THEN
!GPUNotes loop over frequencies
!$ACC LOOP INDEPENDENT VECTOR
            DO IK2=IK1,IK-DIKCUMUL
              IF (BTH0(IK2).GT.SSDSBR) THEN
                IS2=(IK2-1)*NTH
                RENEWALFREQ=RENEWALFREQ+DOT_PRODUCT(CUMULW(IS2+1:IS2+NTH,IS),BRLAMBDA(IS2+1:IS2+NTH))
              END IF
            END DO
          END IF
!
! Computes wave turbulence interaction
!
          COSWIND=(ECOS(IS)*COS(USDIR)+ESIN(IS)*SIN(USDIR))
          DTURB=-2.*SIG(IK)*K(IK)*FACTURB*COSWIND  ! Theory -> stress direction
!
! Add effects
!
          SBK(IS)   = DDIAG(IS)*A(IS)
          SBKT(IK)  = SBKT(IK) + SBK(IS)
          DDIAG(IS) = DDIAG(IS) + (SSDSC(3)*RENEWALFREQ+DTURB)
        END DO
      END DO
!
! COMPUTES SOURCES TERM from diagonal term
!
!GPUNotes created loop instead of array syntax.
!$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO  IK=IK1, NK
        DO ITH=1,NTH
          IS=ITH+(IK-1)*NTH
          SRHS(IS) = DDIAG(IS) * A(IS)
        END DO
      END DO
      !IF(IX == DEBUG_NODE) WRITE(*,'(A10,4F20.10)') 'ST4 DISSIP 2', SUM(SRHS), SUM(DDIAG), SSDSC(3)*RENEWALFREQ, DTURB
!
! Adds non-diagonal part: high and low frequency generation
!
      IF (SSDSC(1).GT.0) THEN
!GPUNotes loops over partial spectrum
!$ACC LOOP INDEPENDENT COLLAPSE(3)
        DO IK2 = IK1+DIKCUMUL, NK
          DO IK = IK2-DIKCUMUL, IK2-1
            DO ITH=1,NTH
              IS2=ITH+(IK2-1)*NTH
              IS=ITH+(IK-1)*NTH
              SRHS(IS) = SRHS(IS) + SSDSC(1)*ABS(SBK(IS2))/DIKCUMUL
            END DO
          END DO
        END DO
      END IF
!
      IF (SSDSC(9).GT.0) THEN
!GPUNotes loops over frequencies x2
!$ACC LOOP INDEPENDENT PRIVATE(SRHS)
        DO IK2 = IK1, NK-DIKCUMUL
          IKM= MIN(NK,IK2+2*DIKCUMUL)
          NKM=IKM-(IK2+DIKCUMUL)+1
          FACHF=SSDSC(9)/FLOAT(NKM*NTH)
!$ACC LOOP INDEPENDENT
          DO IK = IK2+DIKCUMUL, IKM
            SRHS((IK-1)*NTH+1:IK*NTH) = SRHS((IK-1)*NTH+1:IK*NTH) + ABS(SBKT(IK2))*FACHF
          END DO
        END DO
      END IF
 
        !IF(IX == DEBUG_NODE) WRITE(*,'(A10,4F20.10)') 'ST4 DISSIP 3', SUM(SRHS), SUM(DDIAG), SUM(A)
!
!  COMPUTES WHITECAP PARAMETERS
!
      WHITECAP(1:2) = 0.
      IF ( .NOT. (FLOGRD(5,7).OR.FLOGRD(5,8) ) ) THEN
        GOTO 1000 
      END IF
!
!
! precomputes integration of Lambda over direction
! times wavelength times a (a=5 in Reul&Chapron JGR 2003) times dk
!
!CODENotes: It for all cases in the miniapp the following is not used
!this means that GPU optimisation has no affect for speed or validation.
!GPUNotes loops over frequencies
!$ACC LOOP INDEPENDENT
      DO IK=1,NK
        COEF4(IK) = SUM(BRLAMBDA((IK-1)*NTH+1:IK*NTH) * DTH) *(2*PI/K(IK)) *  &
                    SSDSC(7) * DDEN(IK)/(DTH*SIG(IK)*CG(IK))
!                   NB: SSDSC(7) is WHITECAPWIDTH
      END DO
!/
!CODENotes: If for all cases in the miniapp the following is not used
!then any GPU optimisation will have no affect for speed or validation.
      IF ( FLOGRD(5,7) ) THEN
!
! Computes the Total WhiteCap Coverage (a=5. ; Reul and Chapron, 2003)
!
!GPUNotes loops over frequencies
!$ACC LOOP INDEPENDENT
        DO IK=IK1,NK
          WHITECAP(1) = WHITECAP(1) + COEF4(IK) * (1-WHITECAP(1))
        END DO
      END IF
!/
!CODENotes: If for all cases in the miniapp the following is not used
!then any GPU optimisation will have no affect for speed or validation.
      IF ( FLOGRD(5,8) ) THEN
!
! Calculates the Mean Foam Thickness for component K(IK) => Fig.3, Reul and Chapron, 2003
!
!GPUNotes loops over frequencies and iterations
!$ACC LOOP INDEPENDENT
        DO IK=IK1,NK
!    Duration of active breaking (TAU*)
          TSTR = 0.8 * 2*PI/SIG(IK)
!    Time persistence of foam (a=5.)
          TMAX = 5.  * 2*PI/SIG(IK)
          DT   = TMAX / 50
          MFT  = 0.
          DO IT = 1, 50
! integration over time of foam persistance
            T = FLOAT(IT) * DT
! Eq. 5 and 6 of Reul and Chapron, 2003
            IF ( T .LT. TSTR ) THEN
              MFT = MFT + 0.4 / (K(IK)*TSTR) * T * DT
            ELSE
              MFT = MFT + 0.4 / K(IK) * EXP(-1*(T-TSTR)/3.8) * DT
            END IF
          END DO
          MFT = MFT / TMAX
!
! Computes foam-layer thickness (Reul and Chapron, 2003)
!
          WHITECAP(2) = WHITECAP(2) + COEF4(IK) * MFT
        END DO
      END IF
!
! End of output computing
1000  CONTINUE
!!$ACC END KERNELS
!!$ACC END DATA





      RETURN
!
! Formats
!
!/
!/ End of W3SDS4 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SDS4
 
 
      END MODULE W3SRC4MD

