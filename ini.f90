!#include "symbol.inc"

      MODULE PREC
      INTEGER, PARAMETER :: q =SELECTED_REAL_KIND(10)
      INTEGER, PARAMETER :: qs=SELECTED_REAL_KIND(5)
!
! this parameter controls the step width in for numerical differentiation
! in some VASP routines
! 1E-5 is very reliable yielding at least 7 digits in all contributions
! to the forces
! 1E-4, however, is better suited for second derivatives
! for reasons of consistentcy with previous versions 1E-5
      REAL(q), PARAMETER :: fd_displacement=1E-5

!  Some important Parameters, to convert to a.u.
!  - AUTOA  = 1. a.u. in Angstroem
!  - RYTOEV = 1 Ry in Ev
!  - EVTOJ  = 1 eV in Joule
!  - AMTOKG = 1 atomic mass unit ("proton mass") in kg
!  - BOLKEV = Boltzmanns constant in eV/K
!  - BOLK   = Boltzmanns constant in Joule/K

      REAL(q), PARAMETER :: AUTOA=0.529177249_q,RYTOEV=13.605826_q
      REAL(q), PARAMETER :: CLIGHT = 137.037  ! speed of light in a.u.
      REAL(q), PARAMETER :: EVTOJ=1.60217733E-19_q,AMTOKG=1.6605402E-27_q, &
     &           BOLKEV=8.6173857E-5_q,BOLK=BOLKEV*EVTOJ

      REAL(q), PARAMETER :: EVTOKCAL=23.06
! FELECT = (the electronic charge)/(4*pi*the permittivity of free space)
!         in atomic units this is just e^2
! EDEPS = electron charge divided by the permittivity of free space
!         in atomic units this is just 4 pi e^2
! HSQDTM = (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)
!
      REAL(q),PARAMETER  :: PI0 =3.141592653589793238_q,TPI=2*PI0
      REAL(q),PARAMETER  :: PI = 3.14159265358979323846264338327950_q
      COMPLEX(q),PARAMETER  :: CITPI = (0._q,1._q)*TPI
      REAL(q),PARAMETER  :: FELECT = 2*AUTOA*RYTOEV,EDEPS=4*PI0*2*RYTOEV*AUTOA,&
     &                   HSQDTM = RYTOEV*AUTOA*AUTOA
! vector field A times momentum times e/ (2 m_e c) is an energy
! magnetic moments are supplied in Bohr magnetons
! e / (2 m_e c) A(r) p(r) = energy
! e / (2 m_e c) m_s x ( r - r_s) / (r-r_s)^3 hbar nabla =
! e^2 hbar^2 / (2 m_e^2 c^2) 1/ lenght^3 = energy
! conversion factor from magnetic moment to energy
! checked independently in SI by Gilles de Wijs
      REAL(q),PARAMETER :: MAGMOMTOENERGY=1/CLIGHT**2*AUTOA**3*RYTOEV

! dimensionless number connecting input and output magnetic moments
! AUTOA e^2 (2 m_e c^2)
      REAL(q),PARAMETER :: MOMTOMOM=AUTOA/CLIGHT/CLIGHT/2

      REAL(q),PARAMETER :: AUTOA2=AUTOA *AUTOA
      REAL(q),PARAMETER :: AUTOA3=AUTOA2*AUTOA
      REAL(q),PARAMETER :: AUTOA4=AUTOA2*AUTOA2
      REAL(q),PARAMETER :: AUTOA5=AUTOA3*AUTOA2
      END MODULE PREC

      MODULE ini
      USE prec
!**********************************************************************
!
!  this module implements a small timer utility to time
!  subroutine
!  (see START_TIMING)
!**********************************************************************

      ! this allows to set a maxmimum of 10 internested timers
      ! that should suffice for VASP
      INTEGER, PARAMETER, PRIVATE :: maxtimer=10

      INTEGER,SAVE :: used_timers=0
      CHARACTER (LEN=6), PRIVATE  :: timer_tag(maxtimer)
      REAL(q)      :: timer_vpu(maxtimer),timer_cpu(maxtimer)

      ! this allows to set a maxmimum of registered allocates
      INTEGER, PARAMETER, PRIVATE :: maxalloc=20
      INTEGER,SAVE :: used_allocs=0
      CHARACTER (LEN=10), PRIVATE :: alloc_tag(maxalloc)
      REAL(q),PRIVATE      :: alloc_event(maxalloc)=0
      REAL(q),PRIVATE      :: alloc_total=0


      INTEGER, PRIVATE ::  MINPGF,MAJPGF,ISWPS,IOOPS,IVCSW
      REAL(q), PRIVATE ::  UTIME,STIME,ETIME,RSIZM,AVSIZ,DAYTIM


      CONTAINS

!***********************************************************************
!
! timing routines
! START_TIMING(TAG) 
! registers a new timing routine with a specific name and 
! initialises the timer
!
!***********************************************************************

      SUBROUTINE START_TIMING(TAG)
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      REAL(q) TV,TC

      CALL SEARCH_TIMING(TAG,ENTRY)
      IF (ENTRY==0) RETURN
    
!      CALL VTIME(TV,TC)
      timer_vpu(ENTRY)=TV
      timer_cpu(ENTRY)=TC
      timer_tag(ENTRY)=TAG 

      END SUBROUTINE


      SUBROUTINE SEPERATOR_TIMING(IU)

      IF (IU>0) WRITE(IU,100)
100   FORMAT(2X,'  --------------------------------------------')
      END SUBROUTINE SEPERATOR_TIMING


      SUBROUTINE SEARCH_TIMING(TAG,ENTRY)
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      
      ! search for entry
      DO ENTRY=1,used_timers
         IF (timer_tag(ENTRY)==TAG) THEN
            RETURN
         ENDIF
      END DO

      IF (ENTRY>maxtimer) THEN
         ! no more entry available
         ENTRY=0
         WRITE(0,*) 'internal ERROR in SEARCH_TIMING: no more timing slot available'
      ELSE
         used_timers=used_timers+1
      ENDIF

      
      END SUBROUTINE SEARCH_TIMING

!***********************************************************************
!
! the following routines can be used to keep track of allocate
! and deallocate commands
! registered allocate calles also supply a tag
!
!***********************************************************************

      SUBROUTINE REGISTER_ALLOCATE(NALLOC, TAG)
      IMPLICIT NONE
      REAL(q) NALLOC
      CHARACTER (LEN=*), OPTIONAL :: TAG
      INTEGER ENTRY

      IF (PRESENT(TAG)) THEN
         CALL SEARCH_ALLOC(TAG, ENTRY)
         alloc_tag(ENTRY)  =TAG
         alloc_event(ENTRY)=alloc_event(ENTRY)+AINT(NALLOC/1000)
      END IF


      alloc_total=alloc_total+AINT(NALLOC/1000)

      END SUBROUTINE

      SUBROUTINE DEREGISTER_ALLOCATE(NALLOC, TAG)
      IMPLICIT NONE
      REAL(q) NALLOC
      CHARACTER (LEN=*), OPTIONAL :: TAG
      INTEGER ENTRY

      IF (PRESENT(TAG)) THEN
         CALL SEARCH_ALLOC(TAG, ENTRY)
         alloc_event(ENTRY)=alloc_event(ENTRY)-AINT(NALLOC/1000)
      END IF

      alloc_total=alloc_total-AINT(NALLOC/1000)

      END SUBROUTINE

      FUNCTION QUERRY_ALLOCATE()
      IMPLICIT NONE
      INTEGER QUERRY_ALLOCATE

      QUERRY_ALLOCATE=alloc_total
      END FUNCTION QUERRY_ALLOCATE


      SUBROUTINE DUMP_ALLOCATE(IU)
      IMPLICIT NONE
      INTEGER IU
      INTEGER ENTRY

      IF (IU>=0) THEN
      WRITE(IU,'(/1X,A,F10.0,A/A/)') 'total amount of memory used by VASP on root node',alloc_total,' kBytes', &
                                '========================================================================'

      DO ENTRY=1,used_allocs
         WRITE(IU,'(3X,A,A,F10.0,A)') alloc_tag(ENTRY),':  ',alloc_event(ENTRY),' kBytes'
      ENDDO
      WRITE(IU,*)
      ENDIF
      END SUBROUTINE DUMP_ALLOCATE

      SUBROUTINE DUMP_ALLOCATE_TAG(IU,TAG)
      IMPLICIT NONE
      INTEGER IU
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG

      IF (IU>=0) THEN
      WRITE(IU,'(/1X,A,A,F10.0,A/A/)') 'memory high mark on root node inside ',TAG, alloc_total,' kBytes', &
                                '========================================================================'

      DO ENTRY=1,used_allocs
         WRITE(IU,'(3X,A,A,F10.0,A)') alloc_tag(ENTRY),':  ',alloc_event(ENTRY),' kBytes'
      ENDDO
      WRITE(IU,*)
      ENDIF
      END SUBROUTINE DUMP_ALLOCATE_TAG


      SUBROUTINE SEARCH_ALLOC(TAG,ENTRY)
      IMPLICIT NONE
      INTEGER ENTRY
      CHARACTER (LEN=*) :: TAG
      
      ! search for entry
      DO ENTRY=1,used_allocs
         IF (alloc_tag(ENTRY)==TAG) THEN
            RETURN
         ENDIF
      END DO

      IF (ENTRY>maxalloc) THEN
         ! no more entry available
         ENTRY=0
         WRITE(0,*) 'internal ERROR in SEARCH_ALLOC: no more registered allocation slots available'
      ELSE
         used_allocs=used_allocs+1
      ENDIF

      
      END SUBROUTINE SEARCH_ALLOC

      FUNCTION SEARCH_ALLOC_MEMORY(TAG)
      IMPLICIT NONE
      REAL(q)  SEARCH_ALLOC_MEMORY
      CHARACTER (LEN=*) :: TAG
      INTEGER ENTRY
      
      ! search for entry
      DO ENTRY=1,used_allocs
         IF (alloc_tag(ENTRY)==TAG) THEN
            SEARCH_ALLOC_MEMORY=alloc_event(ENTRY)
            RETURN
         ENDIF
      END DO

      SEARCH_ALLOC_MEMORY=0
      END FUNCTION SEARCH_ALLOC_MEMORY


!***********************************************************************
!
! dump some information on paging memory etc.
!
!***********************************************************************

      SUBROUTINE INIT_FINAL_TIMING()
      INTEGER IERR

!      CALL TIMING(0,UTIME,STIME,ETIME,MINPGF,MAJPGF, &
!     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)
      IF (IERR/=0) ETIME=0._q

      END SUBROUTINE INIT_FINAL_TIMING

      SUBROUTINE DUMP_FINAL_TIMING(TIU6)

      INTEGER TIU6
      ! local
      INTEGER IERR
      INTEGER NODE_ME, IONODE


!      CALL TIMING(0,UTIME,STIME,DAYTIM,MINPGF,MAJPGF, &
!     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)

      IF (TIU6>=0) THEN

      ETIME=DAYTIM-ETIME

      TOTTIM=UTIME+STIME
      WRITE(TIU6,*) ' '
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(A)') &
     &   ' General timing and accounting informations for this job:'
      WRITE(TIU6,'(A)') &
     &   ' ========================================================'
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(17X,A,F12.3)') ' Total CPU time used (sec): ',TOTTIM
      WRITE(TIU6,'(17X,A,F12.3)') '           User time (sec): ',UTIME
      WRITE(TIU6,'(17X,A,F12.3)') '         System time (sec): ',STIME
      WRITE(TIU6,'(17X,A,F12.3)') '        Elapsed time (sec): ',ETIME
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(17X,A,F12.0)') '  Maximum memory used (kb): ',RSIZM
      WRITE(TIU6,'(17X,A,F12.0)') '  Average memory used (kb): ',AVSIZ
      WRITE(TIU6,*) ' '
      WRITE(TIU6,'(17X,A,I12)')   '         Minor page faults: ',MINPGF
      WRITE(TIU6,'(17X,A,I12)')   '         Major page faults: ',MAJPGF
      WRITE(TIU6,'(17X,A,I12)')   'Voluntary context switches: ',IVCSW
      ENDIF

      END SUBROUTINE

      END MODULE

      MODULE BASE
      USE prec
!
! header type information
!
        TYPE info_struct
!only  INFO
! mode information
        LOGICAL LREAL                ! real space projection/ reciprocal space proj.
        LOGICAL LOVERL               ! vanderbilt type PP read in ?
        LOGICAL LCORE                ! any partial core read in ?
        LOGICAL LCHCON               ! charge density constant during run
        LOGICAL LCHCOS               ! allways same as above in cur. impl.
        LOGICAL LONESW               ! use all bands simultaneous
        LOGICAL LONESW_AUTO          ! switch automatically between LONESW and DIIS
        LOGICAL LPRECONDH            ! precondition the subspace rotation matrix
        LOGICAL LDAVID               ! use block davidson
        LOGICAL LEXACT_DIAG          ! use exact diagonalization
        LOGICAL LRMM                 ! use RMM-DIIS algorithm
        LOGICAL LORTHO               ! orthogonalize
        LOGICAL LABORT,LSOFT         ! soft / hard stop
        LOGICAL LSTOP                ! stop with this iteration
        LOGICAL LPOTOK               ! local potential ok ?
        LOGICAL LMIX                 ! mixing done
        LOGICAL LMETAGGA             ! metagga 
        LOGICAL LASPH                ! aspherical radial PAW 
! exchange correlation spin
        LOGICAL LXLDA                ! calculate LDA exchange only in POTLOK
        INTEGER ISPIN                ! spin 1=no 2 =yes
        REAL(q)    RSPIN             ! 2 for spinpolar. 1 for non spinpol.
! electronic relaxation
        INTEGER ISTART               ! how to start up
        INTEGER ICHARG               ! initial charge density
        INTEGER INIWAV               ! how to initialize wavefunctions
        INTEGER INICHG               ! how to initialize charge
        INTEGER NELM                 ! maximal number of el-steps in SCC
        INTEGER NELMALL              ! maximal number of el-steps in direct optimization
        INTEGER NELMIN               ! minimal number of el-steps
        INTEGER NELMDL               ! number of delay el-steps
        INTEGER IALGO                ! algorithm for el-relax
        INTEGER NDAV                 ! number of steps in RMM-DIIS
        REAL(q)    TIME                 ! timestep for el (IALGO>50)
        REAL(q)    WEIMIN,EBREAK,DEPER  ! control tags for elm
        REAL(q)    EDIFF                ! accuracy for electronic relaxation
        LOGICAL LDIAG                ! level reordering allowed or not
        LOGICAL LSUBROT              ! sub space rotation to optimize rotation matrix
        LOGICAL LPDIAG               ! sub space rotation before iterat. diag.
        LOGICAL LCDIAG               ! recalculate eigenvalues after iterat. diag.
! cutoff information
        REAL(q)    ENMAX             ! cutoff for calculations
        REAL(q)    ENINI             ! cutoff during delay
        REAL(q)    ENAUG             ! cutoff for augmentation charges
! some important things
        REAL(q)    EALLAT            ! total energy of all atoms
        REAL(q) NELECT               ! number of electrons
        REAL(q) NUP_DOWN             ! spin multiplicity
        INTEGER NBANDTOT             ! total number of bands
        INTEGER MCPU                 ! max number of proc (dimensioned)
        INTEGER NCPU                 ! actual number of proc
        LOGICAL LCORR                ! correction to forces
        INTEGER IQUASI               ! eigenvalue/occupation number corrections
        INTEGER TURBO                ! turbo mode
        INTEGER IFOLD                ! folded eigenproblem: 1=LFOLD,2=LFOLDHpsi
        INTEGER IRESTART             ! whether to restart: 2=restart with 2 optimized vectors
        INTEGER IHARMONIC            ! harmonic Ritz values
        INTEGER NREBOOT              ! number of reboots
        INTEGER NMIN                 ! reboot dimension
        REAL(q) EREF                 ! reference energy to select bands
        LOGICAL NLSPLINE             ! use spline interpolation to construct projection operator
! characters allways last
        CHARACTER*40 SZNAM1          ! header of INCAR
        CHARACTER*6  SZPREC          ! precision information
      END TYPE

      TYPE in_struct
!only  IO
        LOGICAL LOPEN                ! files open at startup
        INTEGER IU0                  ! unit for error
        INTEGER IU6                  ! unit for stdout
        INTEGER IU5                  ! unit for stdin
        INTEGER NWRITE               ! how much information is written out
        INTEGER IDIOT                ! how much information is written out
        INTEGER ICMPLX               ! size of a complex item upon IO
        INTEGER MRECL                ! maximal size of record length
        LOGICAL LREALD               ! no LREAL read in
        LOGICAL LMUSIC               ! jF (just a joke)
        LOGICAL LFOUND               ! WAVECAR exists ?
        LOGICAL LWAVE                ! write WAVECAR
        LOGICAL LCHARG               ! write CHGCAR
        LOGICAL LVTOT                ! write total local potential to LOCPOT
        LOGICAL LVHAR                ! write Hartree potential to LOCPOT
        LOGICAL LPDENS               ! write partial density (charge density for one band)
        INTEGER LORBIT               ! write orbit/dos
        LOGICAL LELF                 ! write elf
        LOGICAL LOPTICS              ! calculate/write optical matrix elements
        LOGICAL LPETIM               ! timing information
        INTEGER IUVTOT               ! unit for local potential
        LOGICAL INTERACTIVE          ! vasp runs interactive 
        INTEGER IRECLW               ! record lenght for WAVECAR
      END TYPE

      TYPE mixing
!only MIX
        INTEGER IUBROY               ! unit for broyden mixer
        REAL(q)     AMIX             ! mixing parameter A
        REAL(q)     BMIX             ! mixing parameter B
        REAL(q)     AMIX_MAG         ! mixing parameter A for magnetization
        REAL(q)     BMIX_MAG         ! mixing parameter B for magnetization
      REAL(q)     AMIN             ! minimal mixing parameter A
        REAL(q)     WC               ! weight factor for Johnsons method
        INTEGER  IMIX                ! type of mixing
        INTEGER  INIMIX              ! initial mixing matrix
        INTEGER  MIXPRE              ! form of metric for mixing
        LOGICAL  LRESET              ! reset mixer on next call (set when ions move)
        LOGICAL  HARD_RESET          ! force hard reset of mixer (force full reset regardless of MAXMIX)
        INTEGER  MAXMIX              ! maximum number of mixing steps (if positive LRESET does not apply)
        INTEGER  NEIG                ! number of eigenvalues
        INTEGER  MREMOVE             ! how many vectors are removed once
          ! iteration depth is reached
        REAL(q) :: EIGENVAL(512)     ! eigenvalues of dielectric matrix
        REAL(q) :: AMEAN             ! mean eigenvalue
        LOGICAL :: MIXFIRST          ! mix before diagonalization (or after)
      END TYPE

      TYPE symmetry
!only  SYMM
        INTEGER,POINTER:: ROTMAP(:,:,:) ! 
        REAL(q),POINTER:: TAU(:,:)      ! jF
        REAL(q),POINTER:: TAUROT(:,:)   ! jF
        REAL(q),POINTER:: WRKROT(:)     ! jF
        REAL(q),POINTER:: PTRANS(:,:)   ! jF
        REAL(q),POINTER:: MAGROT(:,:)   ! jF
        INTEGER,POINTER:: INDROT(:)     ! jF
        INTEGER ISYM                    ! symmetry on/of
        INTEGER NROT                    ! number of rotations
        INTEGER NPTRANS                 ! number of primitive translations
      END TYPE

      TYPE prediction
!only PRED
        INTEGER IWAVPR               ! prediction of wavefunctions
        INTEGER INIPRE               ! initialized yes/no
        INTEGER IPRE                 ! what was done in wavefunction predic.
        INTEGER IUDIR                ! unit for prediction of wavefunction
        INTEGER ICMPLX               ! size of complex word
        REAL(q)  ALPHA,BETA
      END TYPE

      TYPE dipol
!only DIP
      INTEGER IDIPCO                 ! direction (0 no dipol corrections)
        LOGICAL LCOR_DIP             ! correct potential
        REAL(q) :: POSCEN(3)         ! position of center
        REAL(q) :: DIPOLC(3)         ! calculated dipol
        REAL(q) :: QUAD              ! trace of quadrupol
        INTEGER :: INDMIN(3)         ! position of minimum
        REAL(q) :: EDIPOL,EMONO,E_ION_EXTERN
      REAL(q),POINTER :: FORCE(:,:)
       REAL(q) :: VACUUM(2)          ! vacuum level
      END TYPE

      TYPE smear_struct
!only SMEAR_LOOP
        INTEGER ISMCNT               !
        REAL(q)    SMEARS(200)       !
      END TYPE


      TYPE paco_struct
!only PACO
        INTEGER NPACO                ! number of grid points for pair corr.
        REAL(q)    APACO             ! cutoff
        REAL(q),POINTER :: SIPACO(:) 
      END TYPE                       
                             
      TYPE energy
        REAL(q)    TOTENMGGA         ! total energy for METAGGA calculation
        REAL(q)    TOTENASPH         ! total energy for aspherical GGA
        REAL(q)    EBANDSTR          ! bandstructure energy
        REAL(q)    DENC              ! -1/2 hartree (d.c.)
        REAL(q)    XCENC             ! -V(xc)+E(xc) (d.c.)
        REAL(q)    EXCG              ! E(xc) (LDA+GGA)
        REAL(q)    EXCM              ! E(xc) (metaGGA)
        REAL(q)    EXLDA             ! LDA excchange energy
        REAL(q)    ECLDA             ! LDA correlation energy
        REAL(q)    EXGGA             ! GGA exchange energy
        REAL(q)    ECGGA             ! GGA correlation energy
        REAL(q)    EXHF              ! Hartree-Fock exchange energy
        REAL(q)    EXHF_ACFDT        ! difference between HF energy, and exchange energy in ACFDT
        REAL(q)    EDOTP             ! Electric field \dot Polarization
        REAL(q)    TEWEN             ! Ewald energy
        REAL(q)    PSCENC            ! alpha Z (V(q->0) Z)
        REAL(q)    EENTROPY          ! Entropy term
        REAL(q)    PAWPS,PAWAE       ! paw double counting corrections
        REAL(q)    PAWPSG,PAWAEG     ! paw xc energies (LDA+GGA)
        REAL(q)    PAWCORE           ! exchange correlation energy of core (LDA+GGA)
        REAL(q)    PAWPSM,PAWAEM     ! paw xc energies (metaGGA)
        REAL(q)    PAWCOREM          ! exchange correlation energy of core (metaGGA)
        REAL(q)    PAWPSAS,PAWAEAS   ! paw xc energies (aspherical)
        COMPLEX(q) CVZERO            ! average local potential
      END TYPE

      TYPE type_info
!only T_INFO
        CHARACTER*40 SZNAM2           ! name of poscar file
        INTEGER NTYPD                 ! dimension for types
        INTEGER NTYP                  ! number of types
        INTEGER NTYPPD                ! dimension for types inc. empty spheres
        INTEGER NTYPP                 ! number of types empty spheres
        INTEGER NIOND                 ! dimension for ions
        INTEGER NIONPD                ! dimension for ions inc. empty spheres
        INTEGER NIONS                 ! actual number of ions
        INTEGER NIONP                 ! actual number of ions inc. empty spheres
        LOGICAL LSDYN                 ! selective dynamics (yes/ no)
        LOGICAL LDIRCO                ! positions in direct/recproc. lattice
        REAL(q), POINTER :: POSION(:,:)  ! positions usually same as DYN%POSION
        LOGICAL,POINTER ::  LSFOR(:,:) ! selective dynamics
        INTEGER, POINTER :: ITYP(:)   ! type for each ion
        INTEGER, POINTER :: NITYP(:)  ! number of ions for each type
        REAL(q), POINTER :: POMASS(:) ! mass for each ion type
        REAL(q), POINTER :: RWIGS(:)  ! wigner seitz radius for each ion type
        REAL(q), POINTER :: ROPT(:)   ! optimization radius for each type
        REAL(q), POINTER :: ATOMOM(:) ! initial local spin density for each ion
        REAL(q), POINTER :: DARWIN_R(:) ! parameter for darwin like mass term at each ion
        REAL(q), POINTER :: DARWIN_V(:) ! parameter for darwin like mass term at each ion 
        REAL(q), POINTER :: VCA(:)    ! weight of each species for virtual crystal approximation 
        REAL(q), POINTER :: ZCT(:)    ! "charge transfer" charges for non-scf calculations
        REAL(q), POINTER :: RGAUS(:)  ! widths for Gaussian CT charge distributions
        CHARACTER (LEN=2), POINTER :: TYPE(:)  ! type information for each ion
      END TYPE


!      CHARACTER (LEN=10) :: INCAR='INCAR'

      CONTAINS
!
! small subroutine which tries to give good dimensions for 1 dimension
!
      SUBROUTINE MAKE_STRIDE (N)
      INTEGER N,NEW
      INTEGER, PARAMETER :: NGOOD=16

      NEW=(N+NGOOD)/NGOOD
      NEW=NEW*NGOOD+1
      N=NEW

      END SUBROUTINE

      END MODULE BASE

!**************** SUBROUTINE SPLCOF, SPLCOF_N0 *************************
! RCS:  $Id: ini.F,v 1.3 2002/08/14 13:59:39 kresse Exp $
!
!  Subroutine for calculating spline-coefficients
!  using the routines of the book 'numerical  recipes'
!  on input P(1,N) must contain x-values
!           P(2,N) must contain function-values
!  YP is the first derivatives at the first point
!  if >= 10^30 natural boundary-contitions (y''=0) are used
!
!  for point N always natural boundary-conditions are used in
!  SPLCOF, whereas SPLCOF_N0 assume 0 derivative at N
!  SPLCOF_NDER allows to specify a boundary condition
!  at both end points
!
!***********************************************************************

      SUBROUTINE SPLCOF(P,N,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO

      P(N,4)=0.0_q
      P(N,3)=0.0_q
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE



      SUBROUTINE SPLCOF_N0(P,N,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO
      YNP=0
      IF (YNP> .99E30_q) THEN
        QN=0
        UN=0
      ELSE
        QN=0.5_q
        UN=(3._q/(P(N,1)-P(N-1,1)))*(YNP-(P(N,2)-P(N-1,2))/ &
     &             (P(N,1)-P(N-1,1)))
      ENDIF
      P(N,4)=(UN-QN*P(N-1,3))/(QN*P(N-1,4)+1.)
      P(N,3)=0  ! never used
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE


      SUBROUTINE SPLCOF_NDER(P,N,NDIM,Y1P,YNP)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!
!     determination of spline coefficients
!     ------------------------------------
!     f = ((d*dx+c)*dx+b)*dx+a
!         between adjacent x - values
!
!     result
!     P-ARRAY
!     P(I,1) = X(I)
!     P(I,2) = A(I) = F(I)
!     P(I,3) = B(I)
!     P(I,4) = C(I)
!     P(I,5) = D(I)
!
      IF (Y1P> .99E30_q) THEN
        P(1,4)=0.0_q
        P(1,3)=0.0_q
      ELSE
        P(1,4)=-.5_q
        P(1,3)=(3._q/(P(2,1)-P(1,1)))*((P(2,2)-P(1,2))/ &
     &             (P(2,1)-P(1,1))-Y1P)
      ENDIF

      DO 20 I=2,N-1
        S=(P(I,1)-P(I-1,1))/(P(I+1,1)-P(I-1,1))
        R=S*P(I-1,4)+2._q
        P(I,4)=(S-1._q)/R
        P(I,3)=(6*((P(I+1,2)-P(I,2))/(P(I+1,1)-P(I,1))- &
     &          (P(I,2)-P(I-1,2))/(P(I,1)-P(I-1,1)))/ &
     &          (P(I+1,1)-P(I-1,1))-S*P(I-1,3))/R
   20 ENDDO
      IF (YNP> .99E30_q) THEN
        QN=0
        UN=0
      ELSE
        QN=0.5_q
        UN=(3._q/(P(N,1)-P(N-1,1)))*(YNP-(P(N,2)-P(N-1,2))/ &
     &             (P(N,1)-P(N-1,1)))
      ENDIF
      P(N,4)=(UN-QN*P(N-1,3))/(QN*P(N-1,4)+1.)
      P(N,3)=0  ! never used
!
      DO 30 I=N-1,1,-1
        P(I,4)=P(I,4)*P(I+1,4)+P(I,3)
  30  ENDDO
!
      DO 50 I=1,N-1
        S= P(I+1,1)-P(I,1)
        R=(P(I+1,4)-P(I,4))/6
        P(I,5)=R/S
        P(I,4)=P(I,4)/2.0_q
        P(I,3)=(P(I+1,2)-P(I,2))/S-(P(I,4)+R)*S
   50 ENDDO
      RETURN
      END SUBROUTINE

!
!  helper routine, which copies X and Y arrays to P
!  and than performes the fit on the array Y
!
      SUBROUTINE SPLCPY(X,Y,P,NAC,NDIM,Y1P)
      USE prec

      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
      DIMENSION X(NAC)
      DIMENSION Y(NAC)
      DO 100 N=1,NAC
        P(N,1)=X(N)
        P(N,2)=Y(N)
  100 CONTINUE
      CALL SPLCOF(P,NAC,NDIM,Y1P)
      RETURN
      END SUBROUTINE
!
!  helper routine, which evaluates the spline fit at a specific
!  position
!
      SUBROUTINE SPLVAL(X,F,FDER,P,NAC,NDIM)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION P(NDIM,5)
!  interval bisectioning
      I=1
      J=NAC
      IF (X   <P(I,1)) GO TO 60
      IF (X   <P(J,1)) GO TO 70
      K=J-1
      GOTO 90
   60 K=1
      GOTO 90
   70 K=(I+J)/2
      IF(I==K) GOTO 90
      IF (X   <P(K,1)) GO TO 80
      I=K
      GOTO 70
   80 J=K
      GOTO 70
!
   90 DX=X   -P(K,1)
      F   =((P(K,5)*DX+P(K,4))*DX+P(K,3))*DX+P(K,2)
      FDER=(3.0_q*P(K,5)*DX+2.0_q*P(K,4))*DX+P(K,3)
      END SUBROUTINE


!***********************************************************************
!
! system name date and time
!
!***********************************************************************


      SUBROUTINE MEMORY_CHECK(LOOP,STR)
      USE prec
      IMPLICIT REAL(q) (A-B,D-H,O-Z)

      CHARACTER (LEN=*) STR
      REAL(q) SUM

!      CALL TIMING(0,UTIME,STIME,DAYTIM,MINPGF,MAJPGF, &
!     &            RSIZM,AVSIZ,ISWPS,IOOPS,IVCSW,IERR)
      IF (IERR/=0) WRITE(*,*) 'WARNING main: call to TIMING failed.'
      IF (IERR/=0) WRITE(*,*) 'WARNING main: call to TIMING failed.'
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,'(A)') &
     &   ' General timing and accounting informations for this job:'
      WRITE(*,'(A)') &
     &   ' ========================================================'
      WRITE(*,*) ' '
      WRITE(*,'(17X,A,F12.0)') '  Maximum memory used (kb): ',RSIZM
      WRITE(*,'(17X,A,F12.0)') '  Average memory used (kb): ',AVSIZ
      WRITE(*,*) ' '
      WRITE(*,'(17X,A,I12)')   '         Minor page faults: ',MINPGF
      WRITE(*,'(17X,A,I12)')   '         Major page faults: ',MAJPGF
      WRITE(*,'(17X,A,I12)')   'Voluntary context switches: ',IVCSW

      END SUBROUTINE


!         FUNCTION ERRF(x)
!           USE PREC
!           INTEGER :: SIG
!           REAL(q) :: x,errf,t,y
!           REAL(q),PARAMETER :: a1 =  0.254829592
!           REAL(q),PARAMETER :: a2 = -0.284496736
!           REAL(q),PARAMETER :: a3 =  1.421413741
!           REAL(q),PARAMETER :: a4 = -1.453152027
!           REAL(q),PARAMETER :: a5 =  1.061405429
!           REAL(q),PARAMETER :: p  =  0.3275911
! 
!           SIG = 1
!           IF (x < 0.) SIG = -1
!           x = ABS(x)
!           t = 1._q/(1._q + p*x)
!           y = 1._q - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*EXP(-x*x)
!           errf=SIG*y
!         END FUNCTION ERRF

