      MODULE VASP_POTCAR
         USE pseudo         ! VASP_Pseudo
         USE base           ! VASP_base
         USE ini            ! VASP_ini/PREC
         USE setexm
!         USE core_rel
         USE PSEUDO_struct  ! VASP_Pseudo_struct
         USE rhfatm         ! VASP_rhfatm
         USE GlobalMath
         USE atomdata
         USE aeatom
         USE excor_atom
         USE exx_pseudo
         USE hf_pseudo
         USE numerov_mod
         USE paw_sub
         USE pseudodata
         USE pseudo_sub
         USE radialsr
         USE pseudo_atom  ! ATOMPAW_Pseudo
         USE atompaw_report
!
     CONTAINS
      SUBROUTINE vasp_pseudo(ifinput, ifen, Grid, Orbit, Pot_AE, Pot_FC, success)!(ifinput,ifen,Grid)
         IMPLICIT COMPLEX(q)   (C)
         IMPLICIT REAL(q)  (A-B,D-H,O-Z)

!         TYPE(GridInfo), INTENT(IN) :: Grid
         TYPE(GridInfo) :: Grid
!         REAL(8), INTENT(IN) :: coreden(:)
         Type(OrbitInfo), INTENT(IN) :: Orbit
         TYPE(PotentialInfo), INTENT(IN) :: Pot_AE
         TYPE(PotentialInfo), INTENT(IN) :: Pot_FC
         REAL(q), ALLOCATABLE :: PotHr(:), PotXCr(:), PotAEr(:), PotAECr(:), PotATr(:)
         REAL(q), ALLOCATABLE :: Pot_eff(:), Pot_teff(:)
         REAL(q), ALLOCATABLE :: PotAE(:), PotAE00(:), PotPS(:), PotPSC(:), PotPSCr(:)
         REAL(q), ALLOCATABLE :: POTAE_EFF(:), DPOTAE_EFF(:), POTPS_EFF(:)
         REAL(q), ALLOCATABLE :: pdensity(:), den(:) ,cpdensity(:)
         INTEGER :: nbase, N, irc, irc_core, irc_shap
         REAL(q) :: Q_00 , Q_00c, tq, alpha, beta
         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: PP
         TYPE(INFO_STRUCT) :: INFO
         TYPE(in_struct) IO
!         TYPE(PseudoInfo) :: PAW

         INTEGER :: ifinput,ifen,Z
         INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS, LMAX, LMMAX, LI, MI, LMI
         INTEGER :: LMAX_TABLE, LYMAX 
         REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
         REAL(q) ROOT(2), QR, Qloc
         REAL(q) :: DHARTREE, QCORE,SCALE, DOUBLEAE, EXCG
         REAL(q), ALLOCATABLE :: RHO(:,:,:), POT(:,:,:), V(:,:,:), RHOAE00(:), RHOPS00(:)
         REAL(q), ALLOCATABLE :: POTAEC(:), POTPSC_TEST(:), POT_TEST(:), POTAE_TEST(:), POTPS_TEST(:)
         REAL(q), ALLOCATABLE :: CRHODE(:,:)
         REAL(q), ALLOCATABLE :: RHOLM(:), DLM(:)!, GAUSSIAN(:)
         CHARACTER(LEN=2) :: TYPE(1)
         LOGICAL ::   LPAW,  success, UNSCREN, LXCADD
!         INTEGER, EXTERNAL :: MAXL1
         REAL(q), EXTERNAL :: ERRF

         ALLOCATE (P(1))
         NTYP = 1
         NTYPD = 1
         LDIM = 8
         LDIM2 = (LDIM*(LDIM+1))/2
         LMDIM = 64
         POMASS(1) = 14.0
         RWIGS(1) = 1.5
         TYPE(1) = 'Si'
!         VCA(1) = 1.0
         LPAW = .TRUE.
         IU0 = 6
         IU6 = 7
         IU7 = 8
         IU9 = 13
        IU11 = 15
        IU15 = 19
        IU17 = 21
        IU19 = 23

         IU8 = 12
        IU10 = 14
        IU12 = 16
        IU13 = 17
        IU14 = 18
        IU16 = 20
        IU18 = 22
        IU20 = 24

         Call RD_PSEUDO(INFO, P, NTYP, NTYPD, LDIM, LDIM2, LMDIM, POMASS,     &
     &                  RWIGS, TYPE,                       &
     &                  CHANNELS, IU0, IU6, -1, LPAW)
       PP => P(1)
       SCALE = 2.0*sqrt(PI0)

       CALL SET_SIMP(PP%R)
       ALLOCATE(RHOPS00(PP%R%NMAX), RHOAE00(PP%R%NMAX))

      OPEN(UNIT=7,FILE='VASP_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
      OPEN(UNIT=8,FILE='VASP_WPS',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=7,FILE='VASP_WAE',STATUS='OLD')
         OPEN(UNIT=8,FILE='VASP_WPS',STATUS='OLD')
      ENDIF
      DO i=1, CHANNELS
        DO j=1, PP%R%NMAX
           WRITE(IU6,*) PP%R%R(j), PP%WAE(j,i)
!       if (mod(i,2) == 0)WRITE(IU6,*) PP%R%R(j), PP%WAE(j,i)
           WRITE(IU7,*) PP%R%R(j), PP%WPS(j,i)
        ENDDO
        WRITE(IU6,*)
        WRITE(IU7,*)
      ENDDO

      OPEN(UNIT=13,FILE='VASP_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=13,FILE='VASP_CORE',STATUS='OLD')
      ENDIF
      DO j=1, PP%R%NMAX
         WRITE(IU9,*) PP%R%R(j), PP%RHOAE(j)
      ENDDO
!
      OPEN(UNIT=15,FILE='VASP_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=15,FILE='VASP_PCORE',STATUS='OLD')
      ENDIF
      DO j=1, PP%R%NMAX
         WRITE(IU11,*) PP%R%R(j), PP%RHOPS(j)
      ENDDO


!      CALL SIMPI(PP%R,PP%RHOAE, QCORE)
!      WRITE(6,*) 'QCORE=', QCORE, QCORE*SCALE
!      CALL SIMPI(PP%R,PP%RHOPS, QTEST)
!      WRITE(6,*) 'QCORE=', QTEST, QTEST*SCALE

      LMMAX = (PP%LMAX+1)**2
      LYMAX =MAXL1(PP)*2


      PP%LMAX_CALC=LMAX
      
      CALL SET_PAW_AUG(NTYPD, PP, IU6, PP%LMAX_CALC, .TRUE.)
!      WRITE(6,*) 'LYMAX=', LYMAX
      ALLOCATE(RHO(PP%R%NMAX, LMMAX,1), V(PP%R%NMAX, LMMAX,1))
      ALLOCATE(POT(PP%R%NMAX, LMMAX,1), POTAEC(PP%R%NMAX))
      ALLOCATE(POTAE_EFF(PP%R%NMAX),DPOTAE_EFF(PP%R%NMAX), POTPS_EFF(PP%R%NMAX) )
      ALLOCATE(POT_TEST(PP%R%NMAX),POTAE_TEST(PP%R%NMAX), POTPSC_TEST(PP%R%NMAX) )
      ALLOCATE(POTPS_TEST(PP%R%NMAX))
!      ALLOCATE(V1(PP%R%NMAX, LMMAX,1), V2(PP%R%NMAX, LMMAX,1))
      ALLOCATE(CRHODE(LDIM,LDIM))
!      ALLOCATE(CRHODE(LDIM,LDIM), POT-AE(PP%R%NMAX))
      ALLOCATE(RHOLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
      ALLOCATE(DLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
!

!
      LMAX_TABLE=6; CALL YLM3ST_(LMAX_TABLE)
!      DO j=1, PP%R%NMAX
!         GAUSSIAN(j) = exp(-PP%R%R(j)**2/PP%R%R(1)**2/10000)
!         WRITE(6,*) 'GAUSSIAN=', GAUSSIAN(j)
!      ENDDO 
! CRHODE is n_ij (n_occ)    # 价电子轨道的占据数
      CALL SET_CRHODE_ATOM(CRHODE,PP)
! RHOLM is f_LM            #  价电子占据数构造 f_LM
      CALL TRANS_RHOLM(CRHODE, RHOLM, PP)
! RHO is density of valence # 价电子密度 n_lm^v = f_LM\langle\ps_i|\psi_j\rangle
      RHO = 0
      CALL RAD_CHARGE(RHO(:,:,1), PP%R, RHOLM(:),PP%LMAX, PP%LPS, PP%WAE)
      RHOAE00 = 0
      DO i =1,CHANNELS
         LI = PP%LPS(i)
         DO MI = 1, 2*LI+1
            LMI = LI*LI + MI
!           WRITE(6,*) 'TRHOV=', RHO(PP%R%NMAX, LI+1,1)
            RHOAE00(:) = RHOAE00(:)+RHO(:,LMI,1)
         ENDDO
      ENDDO
!      CALL SIMPI(PP%R,RHOAE00+PP%RHOAE*SCALE, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST

!!!!!!!!! POTAE = V_H[n_v]+V_XC[n_v+n_c] != POTAEC  !!!!!!!!!     
!      POT = 0
!      POTAEC = 0
!!      RHO = 0
!      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
!      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
!                   RHO, PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
!      POTAE_TEST(:) = POT(:,1,1)

!!!!!!!!!!!!!!!!!!!!!!!! POTAEC = V_H[n_Zc] = V_H[n_c]+ Z/r !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Z=INT(PP%ZVALF_ORIG+QCORE)
!      Z=INT(QCORE)

      CALL RAD_POT_HAR(0, PP%R, POTAEC, PP%RHOAE, DHARTEE)
      DO j=1, PP%R%NMAX
         POTAEC(j) = POTAEC(j)/SCALE -FELECT/PP%R%R(j)*Z
      ENDDO
!      DO j =1, PP%R%NMAX
!         WRITE(6,*) POTAEC(j)
!      ENDDO
!      WRITE(6,*)

!!!!!!!!!!!!!!!!!!!!!!!! POTAE_TEST = V_H[n_v] +V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!
      POT = 0
      POT_TEST = 0
      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
                   RHO, PP%RHOAE, POT_TEST, POT, DOUBLEAE, EXCG)
      POTAE_TEST(:) =  -POT(:,1,1)/SCALE !+ 1.0/PP%R%R(:)

!!!!!!!!!!!!!!!!!!!!!!!! POTAE_EFF = V_H[n_v+n_Zc] +V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!
      POT = 0
      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
                   RHO, PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
      POTAE_EFF(:) =  -POT(:,1,1)/SCALE

      OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='OLD')
      ENDIF
!
      DO j=1, PP%R%NMAX
         WRITE(IU15,'(8f20.8)') PP%R%R(j), PP%POTAE(j), POTAE_TEST(j),  &
     &                          (PP%POTAE(j)-POTAE_TEST(j))*PP%R%R(j),  &
     &                          POTAEC(j), POTAE_EFF(j)
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!! POTPS_EFF = A*sin(qloc*r)/r !!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 1, PP%R%NMAX
         POTPS_EFF(j) = POTAE_EFF(j)
      ENDDO
      CALL GRAD(PP%R, POTAE_EFF, DPOTAE_EFF)
!      CALL GRAD(PP%R, POTAEC, DPOTAE_EFF)

!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX)/POTAEC(PP%R%NMAX)
      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX-31)/POTAE_EFF(PP%R%NMAX-31)
!      WRITE(6,*) 'alpha= ', DPOTAE_EFF(PP%R%NMAX), POTAE_EFF(PP%R%NMAX), alpha
      beta = 1.0

      CALL SOLVEBESL_Q(ROOT, alpha, beta, 0, 1)
!      WRITE(6,*) 'ROOT=' , ROOT(1), ROOT(2)
!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX)
      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX-31)
!      alpha = POTAE_EFF(PP%R%NMAX)*PP%R%R(PP%R%NMAX)/sin(QR)
      alpha = POTAE_EFF(PP%R%NMAX-31)*PP%R%R(PP%R%NMAX-31)/sin(QR)
!      WRITE(6,*) 'alpha=' , alpha
      DO j = 1, PP%R%NMAX-31
         QR = Qloc*PP%R%R(j)
         POTPS_EFF(j) = alpha*sin(QR)/PP%R%R(j)
!         WRITE(6,*) PP%R%R(j), POTPS_EFF(j), POTAE_EFF(j)
!         WRITE(6,'(5f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j), PP%POTAE(j), PP%POTPS(j)
      ENDDO

      RHO = 0
      POTPSC_TEST = 0
      POT_TEST = 0
!!!!!!!!!!!!!  POTPS = V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c]  !!!!!!!!!!!!!!!!!!!!!!!!
! RHO is \tilde RHO + \hat n
      CALL RAD_CHARGE(RHO(:,:,1), PP%R, RHOLM(:),PP%LMAX, PP%LPS, PP%WPS)
      CALL RAD_AUG_CHARGE(  RHO(:,:,1), PP%R, RHOLM, PP%LMAX, PP%LPS,  &
              LYMAX, PP%AUG, PP%QPAW )

      RHOPS00 = 0
      DO i =1,CHANNELS
         LI = PP%LPS(i)
         DO MI = 1, 2*LI+1
            LMI = LI*LI + MI
!           WRITE(6,*) 'TRHOV=', RHO(PP%R%NMAX, LI+1,1)
            RHOPS00(:) = RHOPS00(:)+RHO(:,LMI,1)
         ENDDO
      ENDDO
!      CALL SIMPI(PP%R,RHOPS00, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST
      POT = 0
      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
                   RHO, PP%RHOPS, POTPSC_TEST, POT, DOUBLEAE, EXCG)
!      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
!                   RHO, PP%RHOPS, PP%POTPSC, POT, DOUBLEAE, EXCG)
      POTPS_TEST(:) =  -POT(:,1,1)/SCALE

      OPEN(UNIT=25,FILE='VASP_POTPS',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=25,FILE='VASP_POTPS',STATUS='OLD')
      ENDIF

      DO j=1, PP%R%NMAX
         WRITE(25,'(6f20.8)') PP%R%R(j), PP%POTPS(j), POTPS_TEST(j), POTAE_EFF(j), POTPS_EFF(j)
      ENDDO

!!!!!!!!!!!!!  POTPSC =! POTPS_EFF - (V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c])  !!!!!!!!!!!!!!!!!!!!!!!!
      POTPSC_TEST(:) = POTPS_EFF(:) - POTPS_TEST(:)

      OPEN(UNIT=21,FILE='VASP_POTPSC',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=21,FILE='VASP_POTPSC',STATUS='OLD')
      ENDIF

      DO j=1, PP%R%NMAX
         WRITE(IU17,'(6f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j)- POTPS_TEST(j),  &
     &                         (PP%POTPSC(j) - POTPS_EFF(j)+POTPS_TEST(j))*PP%R%R(j)    
      ENDDO


      WRITE(6,*) 'VALUE=', PP%ZVALF_ORIG
      WRITE(6,*) 'NPSPTS=', NPSPTS
      WRITE(6,*) 'PSGMAX=', PP%PSGMAX
      WRITE(6,*) 'NMAX=', PP%R%NMAX

!             WRITE(IU19,*)
!             WRITE(IU19,*) PP%PSDMAX
!             WRITE(IU19,*)

!   ---------------- !!!!!! POT IN RECIPROCAL SPACE FROM POT_V[n_Zc] !!!!!  ----------------------    !     

!!      CALL FOURPOT_TO_Q( PP%RDEP, POT, PP%PSP(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)
!      PP%PSP_TEST(:,2) = 0.0
!      CALL FOURPOT_TO_Q( PP%R%R(PP%R%NMAX), PP%POTPSC, PP%PSP_TEST(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)

!          DO j=1, SIZE(PP%PSP,1)
!             WRITE(IU19,'(4f20.8)') PP%PSP(j,2), -PP%PSP_TEST(j,2)
!          ENDDO

!!!! -------------- UNIT IN VASP     E: Hartree   r: Angstrom ------------------ !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  --------------  UNIT IN AtomPAW  E: Ry        r: A.u.  ------------------ !!!!

!      WRITE(6,*) 'irc=', PAW%irc
!      Grid%n = PAW%irc
!      WRITE(6,*) 'N=', Grid%n
      ALLOCATE(PotHr(Grid%n), PotXCr(Grid%n), PotAEr(Grid%n), PotAECr(Grid%n), PotATr(Grid%n))
      ALLOCATE(Pot_eff(Grid%n), Pot_teff(Grid%n))
      ALLOCATE(PotAE00(Grid%n), PotPS(Grid%n))
      ALLOCATE(PotAE(Grid%n), PotPSC(Grid%n))
      ALLOCATE(den(Grid%n))

!      Q_00 = integrator(Grid, FC%coreden, 1, Grid%n)
!      WRITE(6,*) 'Q_00=', Q_00

      OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='OLD')
      ENDIF
!!!   FC%coreden = n(r)*r^2  ==> consider the UNIT_EXCHANGE
!!!   PP%RHOAE = AUTOAE * FC%coreden/(SCALE *AUTOA*AUTOA)
      DO j=1,Grid%n 
         WRITE(IU10,*) Grid%r(j)*AUTOA, FC%coreden(j)/(SCALE*AUTOA)
      ENDDO
!
      OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='OLD')
      ENDIF
!!!   PAW%phi = phi(r)*r  ==> consider the UNIT_EXCHANGE
!!!   PP%WAE = SQRT(AUTOAE) * PAW%phi/AUTOA  !!! SO THAT the DENSITY IS SAME !!!
      DO i=1, PAW%nbase
         DO j=1, Grid%n
             WRITE(IU8,*) Grid%r(j)*AUTOA, PAW%phi(j,i)/SQRT(AUTOA)
         ENDDO
         WRITE(IU8,*)
      ENDDO

      OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='OLD')
      ENDIF
          WRITE(6,*) 'AEPOT calcualted'

!!!!!!!!! POTAEC = V_Z + V_H[n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          den = 0.d0
          call SetPOT(Grid, FC%coreden, den, PotHr, PotXCr, .TRUE.)
!!!!      VASP_V_Z = -Pot_AE%rvn(j)/Grid%r(j)*RYTOEV/SCALE
!!!!      VASP_V_H[n_c] = -PotHr(j)/Grid%r(j)*RYTOEV/(AUTOA^2*SCALE)   !!!!  FOR the factor FELECT
          PotAECr(:) = PotHr(:) + Pot_AE%rvn(:)/SCALE

!!!!!!!!! POTAE = V_H[n_v] + V_XC[n_c+n_v] !!!!!!!!!!!!!!!!!!!!!!!!!!!
          call SetPOT(Grid, FC%coreden, FC%valeden, PotHr, PotXCr, .FALSE., .TRUE.)
!!!!      VASP_POTAE = -POTAEr(j)/Grid%r(j)*RYTOEV/(AUTOA^2*SCALE)
          PotAEr(:) = PotHr(:)+PotXCr(:)

!!!!!!!!! POTAE_EFF = V_Z+V_H[n_c+n_v]+V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!     
!!!!!!!!! POTAE_EFF ===> Pot_AE%rv
          call SetPOT(Grid, FC%coreden, FC%valeden, PotHr, PotXCr, .TRUE., .TRUE.)
          PotATr(:) = PotHr(:)+PotXCr(:) + Pot_AE%rvn(:) 
!!!!      POT_AE%rv = PotAECr + PotAEr

          DO j=1,Grid%n 
             WRITE(IU16,'(8f20.8)') Grid%r(j)*AUTOA, -PotAEr(j)/Grid%r(j)*RYTOEV, &    !!!  VASP_POT_AE
     &                                               -PotAECr(j)/Grid%r(j)*RYTOEV !, & !!!  VASP_POT_AEC 
!    &                      -PotATr(j)/Grid%r(j)*RYTOEV, Pot_AE%rv(j)/Grid%r(j)*RYTOEV !!!  POTAE_EFF
          ENDDO

!!!!!!!!! tcore_den = sum_i B_i*sin(q_i r)/r !!!!!!!!!!!!!!!!!     
         Call setcoretail2(Grid, FC%coreden)
         OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
            OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='OLD')
         ENDIF
!!!   FC%tcore = \tilde{n}(r)*r^2  ==> consider the UNIT_EXCHANGE
!!!   PP%RHOPS = AUTOAE * FC%tcore/(SCALE *AUTOA*AUTOA)
         DO j=1,Grid%n 
            WRITE(IU12,*) Grid%r(j)*AUTOA, PAW%tcore(j)*AUTOA
         ENDDO
!            irc_core= FindGridIndex(Grid, PAW%rc_core)
!            Q_00 = integrator(Grid, PAW%tcore, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!
         Call SetPAWOptions2(ifinput,ifen,Grid, Orbit,Pot_FC,success)

         OPEN(UNIT=17,FILE='ATOM_WPS',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
            OPEN(UNIT=17,FILE='ATOM_WPS',STATUS='OLD')
         ENDIF
!!!   PAW%tphi = tphi(r)*r  ==> consider the UNIT_EXCHANGE
!!!   PP%WPS = SQRT(AUTOA) * PAW%tphi/AUTOA  !!! SO THAT the DENSITY IS SAME !!!
         DO i=1, PAW%nbase
            DO j=1, Grid%n
               WRITE(IU13,*) Grid%r(j)*AUTOA, PAW%tphi(j,i)/sqrt(AUTOA)
            ENDDO
            WRITE(IU13,*)
         ENDDO

!      ALLOCATE(RHOV(Grid%n), d(Grid%n))
!      RHOV = 0
!      DO io =1, PAW%nbase
!!         WRITE(6,*) 'OCC2=', Orbit%occ(io)
!         if (mod(io,2) .ne. 0) RHOV(:) = RHOV(:)+PAW%phi(:,io)**2*Orbit%occ(io)
!      ENDDO
!      DO j=1, Grid%n
!          WRITE(26,*) Grid%r(j)*AUTOA, RHOV(j)/AUTOA**2/SCALE
!!           if (mod(i,2) == 0) WRITE(IU8,*) Grid%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
!!          WRITE(IU13,*) Grid%r(j)*AUTOA, PAW%tphi(j,i)/sqrt(AUTOA)
!      ENDDO

         Call Report_Pseudobasis(Grid,PAW,ifen)
         Call Set_PAW_MatrixElements(Grid,PAW)
         CALL logderiv(Grid,FCPot,PAW)
         CALL ftprod(Grid)
!!   ---------------- Method 1 ----------------------
         CALL FindVlocfromVeff(Grid,FCOrbit,PAW)

!!   ---------------- Method 2 ----------------------
         DO i = 1, Grid%n
            Pot_eff(i) = (PotAEr(i)+PotAECr(i))/Grid%r(i)
         ENDDO

!!!!!!!!!!!!!!!!!!!! POT_tVeff FROM POT_EFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
         call SetPOT_TEFF(Grid, Pot_eff, Pot_teff)
         
         OPEN(UNIT=22,FILE='ATOM_POTPS',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
             OPEN(UNIT=22,FILE='ATOM_POTPS',STATUS='OLD')
         ENDIF
         DO j=1,Grid%n 
            WRITE(IU18,*) Grid%r(j)*AUTOA, -Pot_teff(j)*RYTOEV
         ENDDO
!         STOP   !! POTPS_EFF TEST CORRECT

!!!!!!!!!! PSEUDO Calculations !!!!!!!!!!!!!!!!!!!!         
!         nbase=PAW%nbase
!
!!         N=FindGridIndex(Grid, PAW%rc_shap)
!!         PAW%irc= N
         
!      d=PAW%den-PAW%tden
!      tq=integrator(Grid,d,1,irc)
!      write(6,*) ' abinit tq = ', tq

!         Q_00= 0.d0
!         do ib=1,nbase
!            do ic=1,nbase
!               Q_00=Q_00+PAW%wij(ib,ic)*PAW%oij(ib,ic)
!            enddo
!         enddo
!         write(6,*) 'TEST: Complete q00 ', Q_00
!
!
!!         PAW%rc_shape
!!         N=FindGridIndex(Grid, PAW%rc_shape)
!!         integrator(Grid,FC%valeden)
!         DO i =1, Grid%n
!            pdensity(i) = PAW%hatshape(i)*Q_00
!            pdensity(i) = FC%valeden(i)+FC%coreden(i) 
!         ENDDO

         Q_00 = 0.d0; Q_00_c = 0.d0
         ALLOCATE(pdensity(Grid%n),cpdensity(Grid%n), PotPSCr(Grid%n))
         pdensity = 0.d0
         cpdensity = 0.d0
         PotPSCr = 0.d0

!         WRITE(25, *) PAW%hatden
!            WRITE(6,*) 'unscreened POTPS calcualted'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        irc=max(PAW%irc,PAW%irc_shap,PAW%irc_vloc,PAW%irc_core)
         irc= FindGridIndex(Grid, PAW%rc_shap)
!         N = Grid%n
!            Q_00 = integrator(Grid, PAW%den, 1, N) 
!           write(6,*) 'Q_00 for atom ', N, Q_00, PAW%den
!            Q_00 = integrator(Grid, PAW%den, 1, irc) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid, PAW%tden, 1, N) 
!           write(6,*) 'Q_00 for atom ', N, Q_00
!            Q_00 = integrator(Grid, PAW%tden, 1, irc) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = (integrator(Grid, PAW%den, 1, irc)-integrator(Grid,PAW%tden,1, irc)) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid, FC%coreden, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid, FC%coreden, 1, N) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid, PAW%tcore, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid, PAW%tcore, 1, N) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!             Q_00= (integrator(Grid, FC%coreden, 1, irc)-integrator(Grid, PAW%tcore,1, irc))
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = (integrator(Grid, PAW%den, 1, irc)-integrator(Grid,PAW%tden,1, irc))  &
!             + (integrator(Grid, FC%coreden, 1, irc)-integrator(Grid, PAW%tcore,1, irc))
         Q_00 = (integrator(Grid, PAW%den, 1, irc)-integrator(Grid,PAW%tden,1, irc)) 
         Q_00c = (integrator(Grid, FC%coreden, 1, irc)-integrator(Grid, PAW%tcore,1, irc))

!        Q_00 = 2.4
            
!        Q_00= integrator(Grid, Q_00*PAW%hatden, 1, irc)
!        write(6,*) 'Q_00 for atom ', irc, Q_00
!        Q_00= integrator(Grid, (10-Q_00)*PAW%hatden, 1, irc)
!        write(6,*) 'Q_00c for atom ', irc, Q_00c
!        pdensity=PAW%tden
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! tn = tn_v + that n = tn_v + \sum Qg
         pdensity=(PAW%tden+(Q_00)*PAW%hatden)
!        pdensity=PAW%tden+Q_00*PAW%hatden
!        cpdensity=PAW%tcore+Q_00c*PAW%hatden

!!!!!!!!! POTPS = V_H[tn_v+that_n]+V_XC[tn_v+that_n+tn_c] !!!!!!!!!!!!!!!!!     
         call SetPOT(Grid, PAW%tcore, pdensity, PotHr, PotXCr, .FALSE., LXCADD=.TRUE.)

         DO j=1,Grid%n 
            PotPS(j) = (PotHr(j)+PotXCr(j))/Grid%r(j)
         ENDDO

!!!!!!!!!!!!!  POTPSC =! POTPS_EFF - (V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c])  !!!!!!!!!!!!!!!!!!!!!!!!
         DO j =1, Grid%n
            PotPSC(j)= -Pot_teff(j)+PotPS(j)
         ENDDO

         OPEN(UNIT=24,FILE='ATOM_POTPSC',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
             OPEN(UNIT=24,FILE='ATOM_POTPSC',STATUS='OLD')
         ENDIF

         DO j=1, Grid%n
            WRITE(IU20,*) Grid%r(j)*AUTOA, PotPSC(j)*RYTOEV   !, -PotAECr(j)/Grid%r(j)*RYTOEV
         ENDDO

!!   ---------------- !!!!!! POT_V[\tilde{n}_Zc] FROM POT_V[n_Zc] !!!!!  ----------------------    !!     
         call SetPOT_TEFF(Grid, PotAECr/Grid%r, PotPSC)
!         DO j=1, Grid%n
!            WRITE(IU20,*) Grid%r(j)*AUTOA, -PotPSC(j)*RYTOEV   !, -PotAECr(j)/Grid%r(j)*RYTOEV
!         ENDDO

!   ---------------- !!!!!! POT IN RECIPROCAL SPACE FROM POT_V[n_Zc] !!!!!  ----------------------    !     
!        CALL FOURPOT_TO_Q()

!        deallocate(coreden)
         DEALLOCATE(RHO, V, RHOAE00, CRHODE, RHOLM)
         DEALLOCATE(POT, POTAE_EFF, DPOTAE_EFF, POTAEC)
         DEALLOCATE(PotHr, PotXCr, PotAEr, PotATr,PotAECr, PotPS, PotAE, PotAE00,  PotPSC)
         DEALLOCATE(pdensity, PotPSCr)

         CLOSE(IU6)
         CLOSE(IU8)
         CLOSE(IU9)
         CLOSE(IU10)
         CLOSE(IU11)
         CLOSE(IU12)
         CLOSE(IU13)
!         CLOSE(IU14)
         CLOSE(IU15)
         CLOSE(IU16)
         CLOSE(IU17)
         CLOSE(IU18)
         CLOSE(IU19)
         CLOSE(IU20)

        CALL Report_Pseudopotential(Grid,PAW)
        CALL SPMatrixElements(Grid,FCPot,FC,PAW)

        RETURN
! 
        END SUBROUTINE vasp_pseudo

        SUBROUTINE SetPOT(Grid2, coreden, valeden, POTHR, POTXCR, LCORADD, LXCADD)
         USE atomdata
         USE aeatom
         USE atomdata
!        TYPE(GridInfo) :: Grid2
        TYPE(GridInfo), INTENT(IN) :: Grid2
        INTEGER :: N
        REAL(q), INTENT(IN) :: coreden(:)
        REAL(q), INTENT(IN) :: valeden(:)
        REAL(q) :: density(Grid2%n)
        REAL(q) :: POTXCR(Grid2%n), POTHR(Grid2%n)
        REAL(q) qc, ecoul, v0, etxc, eex, SCALE
        LOGICAL :: LCORADD
        LOGICAL, OPTIONAL :: LXCADD

!        SCALE = 2.0*sqrt(PI0)

        POTHR =0.d0
        POTXCR =0.d0
!        N= Grid2%n
!        call InitPot(POT, Grid2%n)

!        call Get_Nuclearpotential(Grid2, Pot_AE)
!        WRITE(6,*) 'RVN=', Pot_AE%rvn(Grid2%n)

!        allocate (density(Grid2%n), Vrxc_tmp(Grid2%n))
!        allocate (POTHR(Grid2%n), POTXCR(Grid2%n))
        density = 0.d0
        DO i = 1, Grid2%n 
!           density(i) = coreden(i)+FCOrbit%den(i)
           density(i) = coreden(i)+valeden(i)
        ENDDO

!        Grid2%n = N
!        Grid2%n=FindGridIndex(Grid, PAW%rc)
!        call poisson(Grid2, qc, coreden, POT%rvh, ecoul, v0)
        IF (.NOT. LCORADD) THEN
            call poisson(Grid2, qc, valeden, POTHR, ecoul, v0)
        ELSE
            call poisson(Grid2, qc, density, POTHR, ecoul, v0)
        ENDIF
!        call poisson(Grid2, qc, density, POT%rvh, ecoul, v0)
!        DO i = 1, Grid2%n 
!           density(i) = coreden(i)+FCOrbit%den(i)
!           WRITE(25,*) Grid2%r(Grid2%n), qc, POT%rvh(Grid2%n)
!            POTHR(i) = POTHR(i)*FELECT*SCALE/2.0*AUTOA
!        ENDDO
        IF (PRESENT(LXCADD)) &
     &     call exch(Grid2, density, POTXCR, etxc,eex)

!        deallocate(density, Vrxc_tmp)
!        deallocate(POTHR, POTXCR)
         RETURN
        END SUBROUTINE SetPOT


        SUBROUTINE SetPOT_TEFF(Grid2, POTAE, POTPS, RC_CUT)
!!!!!!!!!!!!!!!!!!!!!!!! POTPS_EFF = A*sin(qloc*r)/r !!!!!!!!!!!!!!!!!!!!!!!!!!!
         TYPE(GridInfo) :: Grid2
         REAL(q) :: PotPS(Grid2%n), PotAE(Grid2%n)
         REAL(q) :: ql(2), qr, xx(2), bb(2), alpha, beta
         REAL(q) :: RC, g, gp, gpp, gg
         REAL(q) :: jbes, jbesp, jbespp, amat(2,2), al(2)
         REAL(q), OPTIONAL :: RC_CUT
         INTEGER :: irc, i

         IF (PRESENT(RC_CUT)) THEN
            irc= FindGridIndex(Grid2, RC_CUT)
            rc = RC_CUT
         ELSE
            irc = PAW%irc_vloc
            rc = Grid2%r(irc)
         ENDIF

         alpha=1-rc*Gfirstderiv(Grid2,irc,POTAE)/POTAE(irc)
         beta=1.0d0
         call solvbes(xx,alpha,beta,0,1)
         ql(1)=xx(1)/rc

            qr=ql(1)*rc
            call jbessel(jbes,jbesp,jbespp, 0, 1, qr)

         DO i=irc,Grid2%n
            POTPS(i) = POTAE(i)
         ENDDO

         al(1) = POTAE(irc)*rc/SIN(xx(1))
         do i=1,irc-1
            qr=ql(1)*Grid2%r(i)
            PotPS(i)=al(1)*SIN(qr)/Grid2%r(i)
         enddo

         RETURN
         END SUBROUTINE SetPOT_TEFF

!*******************************************************************
!
! set the array CRHODE such that it corresponds
! to the atomic reference occupancies
! 
!*******************************************************************

      SUBROUTINE SET_CRHODE_ATOM(CRHODE, PP)
      USE prec
      USE pseudo
      IMPLICIT NONE
      REAL(q) :: CRHODE(:,:)

      TYPE (potcar) :: PP
      INTEGER LM, LL, LHI, LOW, MMAX, L, LP, M

      CRHODE=0

      LOW=1
      LM =1
      block: DO
         LL=PP%LPS(LOW)
         ! search block with same L
         DO LHI=LOW,PP%LMAX
            IF (LL/=PP%LPS(LHI)) EXIT
         ENDDO
         LHI=LHI-1
         MMAX=2*LL+1

         DO L =LOW,LHI
         DO LP=LOW,LHI
            DO M =0,MMAX-1
               CRHODE(LM+(L-LOW)*MMAX+M,LM+(LP-LOW)*MMAX+M)=PP%QATO(L,LP)
            ENDDO
         ENDDO
         ENDDO
      
         ! set new LOW value and LM value and go on
         LM=LM+(LHI-LOW+1)*MMAX
         LOW=LHI+1
         IF (LOW > PP%LMAX) EXIT block
      ENDDO block
      END SUBROUTINE SET_CRHODE_ATOM



!*******************************************************************
!
!  transform the real part of the occupancies RHO(lm,l'm') 
!  to RHO(ll',L,M) using Clebsch Gordan coefficients'
!
!' RHO(ll',LM) =  sum C(LM,ll',mm')  RHO(lm,l'm')
!  where C(LM,ll',mm') = \int Y_LM Y_lm Y_l'm' d Omega
!
!  the storage layout of RHO(llp,LM) is somewhat akward
!  for each l lp pair, (2l+1) (2lp+1) elements must be stored
!  they are stored in the order
!    Lmin,M=0 ... Lmin,M=2*Lmin+1
!     ...
!    Lmax,M=0 ... Lmax,M=2*Lmax+1
!  where Lmin and Lmax are given by the triangular rule
!  Lmin = | l-l' |  and Lmax = | l+l'|
!  certain elements in this array will be always zero, because
!  the sum rule allows only L=Lmin,Lmin+2,...,Lmax
!
!*******************************************************************


    SUBROUTINE TRANS_RHOLM( RHOLLMM, RHOLM, P)
      USE pseudo
      USE asa
      IMPLICIT NONE
      TYPE (potcar) P
      REAL(q) ::RHOLLMM(:,:)   ! net augmentation charge
      REAL(q) RHOLM(:)       ! local charge for each L,M
    ! local varible
      INTEGER CH1,CH2,LL,LLP,LM,LMP,LMINDX,ISTART,IEND,IC,M,MP
      INTEGER IBASE,JBASE,LMIN,LMAX
      REAL(q) FAKT
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,P%LMAX
      LMP=LM
      DO CH2=CH1,P%LMAX

         ! quantum numbers l and lp of these two channels
         LL =P%LPS(CH1)
         LLP=P%LPS(CH2)

         CALL YLM3LOOKUP(LL,LLP,LMINDX)
         ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         
         ! transform coefficients
         FAKT=1
         IF (CH1 /= CH2) THEN
            FAKT=FAKT*2
         ENDIF
         ! JS(IC) is pointing to L*L+M+1, we must subtract LMIN*LMIN
         JBASE=IBASE-LMIN*LMIN
         
         RHOLM(IBASE+1:IBASE+(2*LL+1)*(2*LLP+1))=0
!#ifdef debug
!         WRITE(0,*) 'CALC_RHOLM: RHOLLMM is',CH1,CH2
!         DO MP=1,2*LLP+1
!            WRITE(0,'(10F10.7)') (RHOLLMM(LM+M-1,LMP+MP-1),M=1,2*LL+1)
!         ENDDO
!#endif
         DO M =1,2*LL+1
            DO MP=1,2*LLP+1
               LMINDX=LMINDX+1
               
               ISTART=INDCG(LMINDX)
               IEND  =INDCG(LMINDX+1)
               DO  IC=ISTART,IEND-1
                  RHOLM(JS(IC)+JBASE)= RHOLM(JS(IC)+JBASE)+ YLM3(IC)*FAKT* &
                       &         RHOLLMM(LM+M-1,LMP+MP-1)
               ENDDO
            ENDDO
         ENDDO

!#ifdef debug
!      WRITE(0,*) 'CALC_RHOLM: augmentation charges are',CH1,CH2
!      DO MP=LMIN,LMAX
!         WRITE(0,'(I3,10F10.7)') MP,(RHOLM(JBASE+MP**2+M),M=1,MP*2+1)
!      ENDDO
!#endif

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO


      END SUBROUTINE TRANS_RHOLM

      FUNCTION MAXL1(P)
      USE prec
      USE pseudo
      IMPLICIT NONE
      INTEGER MAXL1
      TYPE (potcar) P
! local varibale
      INTEGER I,LTMP,CHANNELS

      MAXL1=0

      CHANNELS=P%LMAX
      LTMP=0
      DO I=1,CHANNELS
         LTMP=MAX( P%LPS(I),LTMP )
      ENDDO
      MAXL1=LTMP

      END FUNCTION MAXL1


!*******************************************************************
!
!  start up procedure for PAW
!  checks internal consistency of the QPAW
!  and sets up the compensation charges
!  on the radial grid, and spline coefficients which are used
!  to interpolate compensation charges in us.F
!
!*******************************************************************

      SUBROUTINE SET_PAW_AUG(NTYPD, P, IU6, LMAX_CALC, LCOMPAT)
        USE radial
!        USE LDAPLUSU_MODULE
        IMPLICIT NONE
        INTEGER NTYPD
        INTEGER IU6
        TYPE (potcar),TARGET :: P
        TYPE (rgrid),POINTER :: R
        INTEGER LMAX_CALC  ! if -1 mimic US PP (see below)
      ! local variables
        INTEGER NTYP,CHANNELS,LMAX,L,LP,N,I, LMAX_MIX
        INTEGER, PARAMETER :: NQ=2
        REAL(q)  QQ(NQ)
        REAL(q)  A(NQ),B,ALPHA
        REAL(q)  SUM,QR,BJ,STEP,X
        REAL(q),PARAMETER ::  TH=1E-6_q
        LOGICAL LCOMPAT


      ! if LMAX_CALC == -1 the PAW will run in a special mode
      ! resulting in essentially US-PP like behaviour
      ! i.e. all terms are linearized around the atomic reference configuration
      ! and one center terms are not evaluated

!        IF (LMAX_CALC ==-1) THEN
!           MIMIC_US=.TRUE.
!        ELSE
!           MIMIC_US=.FALSE.
!        ENDIF

        LMAX_MIX=0

        typ: DO NTYP=1,NTYPD
           IF (.NOT. ASSOCIATED(P%QPAW)) CYCLE

           ! maximal L mixed in Broyden mixer
           LMAX_MIX=MAX(LMAX_MIX,P%LMAX_CALC)

           R => P%R
           CALL  RAD_ALIGN(R,R%RMAX) ! reallign RMAX with grid

           CALL RAD_CHECK_QPAW( R, P%LMAX, &
                P%WAE, P%WPS, P%QPAW, P%QTOT , P%LPS, P%RDEP )
!           IF ((USELDApU().OR.LCALC_ORBITAL_MOMENT()).AND.INTEGRALS_LDApU()) THEN
!              CALL OVERLAP_AE(R,P%RDEP,NTYP,NTYPD,P%LMAX,P%WAE,P%LPS)
!           ENDIF

           CHANNELS=P%LMAX
           DO L=1, CHANNELS
              DO LP=1, P%LMAX
                 IF (P%LPS(L)==P%LPS(LP)) THEN
                    P%QION(L,LP)=P%QPAW(L,LP,0)
                 ENDIF
              ENDDO
           ENDDO

           LMAX=0
           DO I=1,CHANNELS
              LMAX=MAX( P%LPS(I),LMAX )
           ENDDO

           LMAX=LMAX*2+1              ! maximum l in augmentation charges
                                      ! to allow use of one-center dipole operators increase by one
           ALLOCATE(P%QDEP(NPSRNL,5,0:LMAX), &
                    P%AUG (R%NMAX,0:LMAX) )

!           IF (IU6>=0) WRITE(IU6,"(' L augmenation charges for type=',I4)") NTYP

           ll: DO L=0,LMAX

        ! find q values
              CALL AUG_SETQ(L,R,R%RMAX,QQ,A,LCOMPAT)
!              IF (IU6>=0) WRITE(IU6,2) L,QQ,A

        ! setup augmentation charge on radial grid rho(r) r^2

              DO N=1,R%NMAX
                 SUM=0
                 IF (R%R(N) <= R%RMAX) THEN
                    DO I=1,NQ
                       QR=QQ(I)*R%R(N)
                       CALL SBESSEL( QR, BJ, L)
                       SUM=SUM+BJ*A(I)*R%R(N)*R%R(N)
                    ENDDO
                 ENDIF
                 P%AUG(N,L)=SUM
              ENDDO

!        ! setup spline for augmentation charge
!              ! the spline ends at PSDMAX*(NPSRNL-1)/NPSRNL see SETDEP
!              STEP= R%RMAX/(NPSRNL-1)
!        ! this resets PSDMAX which is from now on no longer identical to R%RMAX
!              P(NTYP)%PSDMAX=NPSRNL*STEP
!
!              DO N=1,NPSRNL
!                 X=STEP*(N-1)
!                 SUM=0
!                 DO I=1,NQ
!                    QR=QQ(I)*X
!                    CALL SBESSEL( QR, BJ, L)
!                    SUM=SUM+BJ*A(I)
!                 ENDDO
!                 P(NTYP)%QDEP(N,1,L) = X
!                 P(NTYP)%QDEP(N,2,L) = SUM
!              ENDDO
!              ! derivative at startpoint
!              X=STEP/1000
!              SUM=0
!              DO I=1,NQ
!                 QR=QQ(I)*X
!                 CALL SBESSEL( QR, BJ, L)
!                 SUM=SUM+BJ*A(I)
!              ENDDO
!              SUM=(SUM-P(NTYP)%QDEP(1,2,L))/X
!              ! derivative is zero for all but L=1
!              IF (L/=1) THEN
!                 SUM=0
!              ENDIF
!              CALL SPLCOF(P(NTYP)%QDEP(1,1,L),NPSRNL,NPSRNL,SUM)
           ENDDO ll
!           P(NTYP)%AUG_SOFT =>P(NTYP)%AUG
!
        ENDDO typ
!
!! mix only L=SET_LMAX_MIX_TO components in Broyden mixer
!!        LMAX_MIX=MIN(LMAX_MIX,SET_LMAX_MIX_TO)
!
!! if US-PP are mimicked no need to mix any onsite components
!!        IF (MIMIC_US) LMAX_MIX=-1
!
      END SUBROUTINE SET_PAW_AUG

      SUBROUTINE SOLVEBESL_Q(ROOT, alpha, beta, L, NQ)
        USE radial

        IMPLICIT REAL(q) (A-H,O-Z)
        TYPE (rgrid) R
        INTEGER,INTENT(IN) :: L, NQ
        REAL(q),INTENT(IN) :: alpha, beta
        REAL(q),INTENT(OUT) :: ROOT(NQ)

        REAL(q),PARAMETER :: Dh=1.D-2, tol=1.D-14

        INTEGER :: NROOT
        REAL(8) :: dum,y1,y2,jbes,jbesp,QR,QX,hh

        QR = Dh; NROOT = 0
        DO WHILE (NROOT < NQ)
           CALL SBESSE3(QR, BJ, BJP, BJPP, L)
           y1 = alpha*BJ + beta*QR*BJP
           QR = QR + Dh
           CALL SBESSE3(QR, BJ, BJP, BJPP, L)
           y2 = alpha*BJ + beta*QR*BJP

           DO WHILE( y1*y2 >= 0.D0)
              QR = QR + Dh
              CALL SBESSE3(QR, BJ, BJP, BJPP, L)
              y2 = alpha*BJ + beta*QR*BJP
          ENDDO

          hh = Dh; QX = QR
          DO WHILE (hh > tol)
             hh = 0.5D0*hh
             IF(y1*y2 < 0) THEN
                     QX = QX - hh
             ELSE
                     QX = QX + hh
             ENDIF
             CALL SBESSE3(QX, BJ, BJP, BJPP, L)
              y2 = alpha*BJ + beta*QR*BJP
          ENDDO
          NROOT = NROOT+1
          ROOT(NROOT) = QX
        ENDDO
        
      END SUBROUTINE SOLVEBESL_Q

      END MODULE VASP_POTCAR
