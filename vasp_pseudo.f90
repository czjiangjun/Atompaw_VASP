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
!      SUBROUTINE vasp_pseudo(ifinput, ifen, Grid, Orbit, Pot_AE, Pot_FC, IMESH)  !(ifinput,ifen,Grid)
      SUBROUTINE vasp_pseudo(ifinput, ifen, Grid, Orbit, Pot_AE, Pot_FC)  !(ifinput,ifen,Grid)
         IMPLICIT COMPLEX(q)   (C)
         IMPLICIT REAL(q)  (A-B,D-H,O-Z)

!         TYPE(GridInfo), INTENT(IN) :: Grid
         TYPE (rgrid)    :: PR        ! radial grid for Projector
         TYPE (rgrid)    :: FR        ! radial grid for Projector
         TYPE(GridInfo) :: Grid
!         REAL(8), INTENT(IN) :: coreden(:)
         Type(OrbitInfo), INTENT(IN) :: Orbit
         TYPE(PotentialInfo), INTENT(IN) :: Pot_AE
         TYPE(PotentialInfo), INTENT(IN) :: Pot_FC
         REAL(q), ALLOCATABLE :: PotHr(:), PotXCr(:), PotAEr(:), PotAECr(:), PotATr(:)
         REAL(q), ALLOCATABLE :: Pot_eff(:), Pot_teff(:)
         REAL(q), ALLOCATABLE :: PotAE(:), PotAE00(:), PotPS(:), PotPSC(:), PotPSCr(:)
         REAL(q), ALLOCATABLE :: POTAE_EFF(:), DPOTAE_EFF(:), POTPS_EFF(:)
         REAL(q), ALLOCATABLE :: POTPS_G_EFF(:), POTPS_G(:), CORPS_G(:), RHOPS_G(:)
         REAL(q), ALLOCATABLE :: pdensity(:), den(:) ,cpdensity(:)
         REAL(q), ALLOCATABLE :: TMP(:),DTMP(:),VTMP(:,:,:), DIJ(:,:), DION(:,:)
         REAL(q), ALLOCATABLE :: PARWKINAE(:,:), PARWKINPS(:,:)
         INTEGER :: nbase, N, irc, irc_core, irc_shap
         REAL(q) :: Q_v, Q_00 , Q_00c, tq, alpha, beta, random_x
         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: PP
         TYPE(INFO_STRUCT) :: INFO
         TYPE(in_struct) IO
!         TYPE(PseudoInfo) :: PAW

         INTEGER :: ifinput,ifen,Z
         INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS, LMAX, LMMAX, LI, MI, LMI
         INTEGER :: LMAX_TABLE, LYMAX 
         INTEGER :: CH1, CH2, LL, LLP, LMIN, LMAIN
         REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
         REAL(q) ROOT(2), QR, Qloc, R_CUT, QQ
         REAL(q) :: DHARTREE, QCORE,SCALE, DOUBLEAE, EXCG
         REAL(q), ALLOCATABLE :: RHO(:,:,:), POT(:,:,:), V(:,:,:), RHOAE00(:), RHOPS00(:)
         REAL(q), ALLOCATABLE :: POTAEC(:), POTPSC_TEST(:), POT_TEST(:), POTAE_TEST(:), POTPS_TEST(:)
         REAL(q), ALLOCATABLE :: PSPNL(:), PSPRNL(:)
         REAL(q), ALLOCATABLE :: POTPSC_CHECK(:)
         REAL(q), ALLOCATABLE :: CRHODE(:,:)
         REAL(q), ALLOCATABLE :: RHOLM(:), DLM(:), DLLMM(:,:,:)!, GAUSSIAN(:)
         REAL(q), ALLOCATABLE :: DHXC(:,:), QPAW(:,:,:)
         CHARACTER(LEN=2) :: TYPE(1)
         LOGICAL ::   LPAW, UNSCREN, LXCADD
!         INTEGER, EXTERNAL :: MAXL1
!         LOGICAL, OPTIONAL :: success
         LOGICAL :: success
!         INTEGER, OPTIONAL :: IMESH
         INTEGER :: IMESH
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
         IU5 = 7
         IU6 = -1
         IU7 = 8
         IU9 = 13
        IU11 = 15
        IU15 = 19
        IU17 = 21
        IU19 = 23
        IU21 = 25
        IU23 = 27
        IU25 = 29
        IU27 = 31

         IU8 = 12
        IU10 = 14
        IU12 = 16
        IU13 = 17
        IU14 = 18
        IU16 = 20
        IU18 = 22
        IU20 = 24
        IU22 = 26

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
           WRITE(IU5,*) PP%R%R(j), PP%WAE(j,i)
!       if (mod(i,2) == 0)WRITE(IU6,*) PP%R%R(j), PP%WAE(j,i)
           WRITE(IU7,*) PP%R%R(j), PP%WPS(j,i)
        ENDDO
        WRITE(IU5,*)
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


      CALL SIMPI(PP%R,PP%RHOAE, QCORE)
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
      ALLOCATE(POTPSC_CHECK(PP%R%NMAX))
      ALLOCATE(POTPS_G_EFF(PP%R%NMAX))
      ALLOCATE(POTPS_G(NPSPTS), CORPS_G(NPSPTS), RHOPS_G(NPSPTS))
      ALLOCATE(POTPS_TEST(PP%R%NMAX))
!      ALLOCATE(V1(PP%R%NMAX, LMMAX,1), V2(PP%R%NMAX, LMMAX,1))
      ALLOCATE(CRHODE(LMDIM,LMDIM))
!      ALLOCATE(CRHODE(LDIM,LDIM), POT-AE(PP%R%NMAX))
      ALLOCATE(RHOLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
!      ALLOCATE(DLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
!

!
      LMAX_TABLE=6; CALL YLM3ST_(LMAX_TABLE)
!      DO j=1, PP%R%NMAX
!         GAUSSIAN(j) = exp(-PP%R%R(j)**2/PP%R%R(1)**2/10000)
!         WRITE(6,*) 'GAUSSIAN=', GAUSSIAN(j)
!      ENDDO 
! CRHODE is n_ij (n_occ)    # 价电子轨道的占据数 考虑 lm,l'm'
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
!     &            RHO, PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
!        CALL POP_XC_TYPE
!      POTAE_TEST(:) = POT(:,1,1)

!!!!!!!!!!!!!!!!!!!!!!!! POTAEC = V_H[n_Zc] = V_H[n_c]+ Z/r !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Z=INT(PP%ZVALF_ORIG+QCORE*SCALE)
!      Z=INT(QCORE)

      CALL RAD_POT_HAR(0, PP%R, POTAEC, PP%RHOAE, DHARTEE)
      DO j=1, PP%R%NMAX
!         WRITE(6,*) 'TEST =', Z/PP%R%R(j)*FELECT
         POTAEC(j) =   POTAEC(j)/SCALE - Z/PP%R%R(j)*FELECT  !!!   V_H[n_Zc] =  V_H[n_c] - Z/r 
!         POTAEC(j) =   POTAEC(j)/SCALE - (Z-PP%ZVALF_ORIG)/PP%R%R(j)*FELECT  !!!   V_H[n_Zc] =  V_H[n_c] - Z/r 
      ENDDO
!      DO j =1, PP%R%NMAX
!         WRITE(91,*) PP%R%R(j), POTAEC(j)
!      ENDDO
!      WRITE(6,*)

!!!!!!!!!!!!!!!!!!!!!!!! POTAE_TEST =  V_H[n_v] +V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!
      POT = 0
      POT_TEST = 0
      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
     &            RHO, PP%RHOAE, POT_TEST, POT, DOUBLEAE, EXCG)
      CALL POP_XC_TYPE
      POTAE_TEST(:) =  POT(:,1,1)/SCALE                       !!!    POTAE = V_H[n_v] +V_XC[n_v+n_c] 

!      DO i =1,CHANNELS
!         LI = PP%LPS(i)
!         DO MI = 1, 2*LI+1
!            LMI = LI*LI + MI
!            POTAE_EFF(:) =  POTAE_EFF(:) + POT(:,LMI,1)
!        ENDDO
!      ENDDO

!      POTAE_TEST(:) =  -POT(:,1,1)/SCALE + PP%ZVALF_ORIG/PP%R%R(:)/SCALE/SQRT(2.0)
!      DO j=1, PP%R%NMAX
!         POTAE_TEST(j) = - POT(j,1,1)/SCALE + FELECT*Z*(1.0-ERRF(PP%R%R(j)/2.0/AUTOA))/PP%R%R(j)
!         WRITE(6,'(8f20.8)') PP%R%R(j), POT(j,1,1)/SCALE, ERRF(PP%R%R(j)/AUTOA), &
!     &             PP%POTAE(j)+POT(j,1,1)/SCALE, FELECT*Z*ERRF(PP%R%R(j)/PP%R%R(200))/PP%R%R(j) , &
!     &             FELECT*Z/PP%R%R(j)/(PP%POTAE(j)+POT(j,1,1)/SCALE )
!      ENDDO

!!!!!!!!!!!!!!!!!   POTAE_EFF = V_H[n_Zc] + V_H[n_v] +V_XC[n_v+n_c]  !!!!!!!!!!!!!!!!!!
      POT = 0
      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
     &            RHO, PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
      CALL POP_XC_TYPE
      POTAE_EFF(:) =  - POT(:,1,1)/SCALE

      OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='OLD')
      ENDIF
!
      DO j=1, PP%R%NMAX
         WRITE(IU15,'(8f20.8)') PP%R%R(j), PP%POTAE(j), - POTAE_TEST(j),  &
     &                           -POTAEC(j), &
     &                           POTAE_EFF(j)
      ENDDO

      PP%POTAE(:) = - POTAE_TEST(:)

!!!!!!!!!!!!!!!!!!!!!!!! POTPS_EFF = A*sin(qloc*r)/r !!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 1, PP%R%NMAX
         POTPS_EFF(j) = POTAE_EFF(j)
      ENDDO
      CALL GRAD(PP%R, POTAE_EFF, DPOTAE_EFF)
      IMESH = 41
!      IMESH = 0

!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX)/POTAEC(PP%R%NMAX)
      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX-IMESH)/POTAE_EFF(PP%R%NMAX-IMESH)
!      WRITE(6,*) 'alpha= ', DPOTAE_EFF(PP%R%NMAX), POTAE_EFF(PP%R%NMAX), alpha
      beta = 1.0

      CALL SOLVEBESL_Q(ROOT, alpha, beta, 0, 1)
!      WRITE(6,*) 'ROOT=' , ROOT(1), ROOT(2)
!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX)
      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX-IMESH)
      alpha = POTAE_EFF(PP%R%NMAX-IMESH)*PP%R%R(PP%R%NMAX-IMESH)/sin(QR)
!      alpha = POTAE_EFF(PP%R%NMAX-41)*PP%R%R(PP%R%NMAX-41)/sin(QR)
!      WRITE(6,*) 'alpha=' , alpha
!      DO j = 1, PP%R%NMAX
      DO j = 1, PP%R%NMAX-IMESH
         QR = Qloc*PP%R%R(j)
         POTPS_EFF(j) = alpha*sin(QR)/PP%R%R(j)
!         WRITE(6,*) PP%R%R(j), POTPS_EFF(j), POTAE_EFF(j)
!         WRITE(6,'(5f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j), PP%POTAE(j), PP%POTPS(j)
      ENDDO

      RHO = 0
      POTPSC_TEST = 0
      POTPSC_CHECK = 0
      POT_TEST = 0
!!!!!!!!!!!!!  POTPS = V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c]  !!!!!!!!!!!!!!!!!!!!!!!!
! RHO is \tilde RHO_v + \hat n
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
     &            RHO, PP%RHOPS, POTPSC_TEST, POT, DOUBLEAE, EXCG)
!      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
!                   RHO, PP%RHOPS, PP%POTPSC, POT, DOUBLEAE, EXCG)
      CALL POP_XC_TYPE
      POTPS_TEST(:) =  POT(:,1,1)/SCALE

      OPEN(UNIT=25,FILE='VASP_POTPS',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=25,FILE='VASP_POTPS',STATUS='OLD')
      ENDIF

      DO j=1, PP%R%NMAX
         WRITE(25,'(6f20.8)') PP%R%R(j), PP%POTPS(j), -POTPS_TEST(j), POTAE_EFF(j), POTPS_EFF(j)
      ENDDO

      PP%POTPS(:) = -POTPS_TEST(:)

!!!!!!!!!!!!!  Method_1:   POTPSC = POTPS_EFF - (V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c])  !!!!!!!!!!!!!!!!!!!!!!!!
      POTPSC_TEST(:) =  - POTPS_EFF(:) - POTPS_TEST(:)
!      DO j=1, PP%R%NMAX
!         WRITE(94,'(6f20.8)') PP%R%R(j), PP%POTPSC(j), POTPSC_TEST(j)
!      ENDDO

!!!!!!!!!!!!   Method_2:   POTPSC = A*sin(qloc*r)/r FROM POTAEC DIRACTLY     !!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 1, PP%R%NMAX
         POTPSC_CHECK(j) = POTAEC(j)
      ENDDO
      CALL GRAD(PP%R, POTAEC, DPOTAE_EFF)

      IMESH = 38
!      IMESH = 0
!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX)/POTAEC(PP%R%NMAX)
      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX-IMESH)/POTAEC(PP%R%NMAX-IMESH)
!      WRITE(6,*) 'alpha= ', DPOTAE_EFF(PP%R%NMAX), POTAE_EFF(PP%R%NMAX), alpha
      beta = 1.0

      CALL SOLVEBESL_Q(ROOT, alpha, beta, 0, 1)
!      WRITE(6,*) 'ROOT=' , ROOT(1), ROOT(2)
!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX)
      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX-IMESH)
!      alpha = POTAE_EFF(PP%R%NMAX)*PP%R%R(PP%R%NMAX)/sin(QR)
      alpha = POTAEC(PP%R%NMAX-IMESH)*PP%R%R(PP%R%NMAX-IMESH)/sin(QR)
!      WRITE(6,*) 'alpha=' , alpha
      DO j = 1, PP%R%NMAX-IMESH
         QR = Qloc*PP%R%R(j)
         POTPSC_CHECK(j) = alpha*sin(QR)/PP%R%R(j)
!         WRITE(6,*) PP%R%R(j), POTPS_EFF(j), POTAE_EFF(j)
!         WRITE(6,'(5f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j), PP%POTAE(j), PP%POTPS(j)
      ENDDO
!      DO j = 1, PP%R%NMAX
!         WRITE(91,*) PP%R%R(j), POTPSC_TEST(j), POTPSC_CHECK(j)
!         WRITE(6,'(5f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j), PP%POTAE(j), PP%POTPS(j)
!      ENDDO



      OPEN(UNIT=21,FILE='VASP_POTPSC',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=21,FILE='VASP_POTPSC',STATUS='OLD')
      ENDIF
!      PP%POTPSC(:) = POTPSC_TEST(:)

!   ---------------- !!!!!! POT IN RECIPROCAL SPACE FROM PSPCORE !   (TESTED)!!!!!  ----------------------    !     
!     FR%NMAX = PP%R%NMAX
!     ALLOCATE(FR%R(FR%NMAX))
!     POTPS_G(:) = 0.0
!     POTPS_G_EFF(:) = 0.0
!
!     DO j=1, PP%R%NMAX
!        FR%R(j) = PP%R%R(j)
!     ENDDO
!!!!!!!!!!!!     POT_G_EFF-part.1-core = Z_VAL*4\pi/q^2*(1-EXP(-q^2/4)) DIRACTLY     !!!!!!!!!!!!!!!!!!!!!!!!
!     DO j=1, SIZE(PP%PSP,1)
!        QQ=(j-1)*PP%PSGMAX/ SIZE(PP%PSP,1)
!        POTPS_G(j) = (PP%ZVALF_ORIG)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA))
!     ENDDO
!
!!!!!!!!!!!!     POT_G-part.1-core = A*sin(q*rloc)/q FROM POT_G_EFF-part.1 DIRACTLY     !!!!!!!!!!!!!!!!!!!!!!!!
!      DO j=1, PP%R%NMAX
!         QQ=(j-1)*PP%PSGMAX/ SIZE(PP%PSP,1)
!         PP%R%R(j) = QQ
!         POTPS_G_EFF(j) = (PP%ZVALF_ORIG)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA))
!      ENDDO
!      CALL GRAD(PP%R, POTPS_G_EFF, DPOTAE_EFF)
!!      IMESH = 0
!      IMESH = 277
!
!!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX)/POTAEC(PP%R%NMAX)
!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX-IMESH)/POTPS_G_EFF(PP%R%NMAX-IMESH)
!!      WRITE(6,*) 'alpha= ', DPOTAE_EFF(PP%R%NMAX), POTAE_EFF(PP%R%NMAX), alpha
!      beta = 1.0
!
!      CALL SOLVEBESL_Q(ROOT, alpha, beta, 0, 1)
!!      WRITE(6,*) 'ROOT=' , ROOT(1), ROOT(2)
!!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX)
!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX-IMESH)
!      alpha = POTPS_G_EFF(PP%R%NMAX-IMESH)*PP%R%R(PP%R%NMAX-IMESH)/sin(QR)
!!      alpha = POTAE_EFF(PP%R%NMAX-41)*PP%R%R(PP%R%NMAX-41)/sin(QR)
!!      WRITE(6,*) 'alpha=' , alpha
!!      DO j = 1, PP%R%NMAX
!      DO j = 1, PP%R%NMAX-IMESH
!         QR = Qloc*PP%R%R(j)
!         POTPS_G(j) = alpha*sin(QR)/PP%R%R(j)
!!         WRITE(6,*) PP%R%R(j), POTPS_EFF(j), POTAE_EFF(j)
!      ENDDO
!
!!!!!!!!!!!!     POT_G-part.2-valence = - 4\pi\rho_c(q)/q^2(1-exp(-q^2/4)) DIRACTLY     !!!!!!!!!!!!!!!!!!!!!!!!
!
!      DO j=1, SIZE(PP%PSP,1)
!         QQ=(j-1)*PP%PSGMAX/ SIZE(PP%PSP,1)
!         POTPS_G(j)= POTPS_G(j)-PP%PSPCOR(j)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA))
!!          POTPS_G(j)= -PP%PSPCOR(j)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/1.2*AUTOA*AUTOA)) &
!!     &     +        (PP%ZVALF_ORIG)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA))
!         WRITE(78,'(8f20.8)') PP%PSP(j,1),  POTPS_G(j)   !,  & 
!!     & PP%PSP(j,2)+ PP%PSPCOR(j)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA)), &
!!     & PP%PSP(j,2)-PP%ZVALF_ORIG*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA)), &
!!     &      (PP%ZVALF_ORIG)*4*PI0*FELECT/QQ/QQ*(1-EXP(-QQ*QQ/4*AUTOA*AUTOA))
!                 
!      ENDDO
!
!   ---------------- !!!!!! FOR CHECK !!!!! -------------------
!      DO j=1, PP%R%NMAX
!         PP%R%R(j) = FR%R(j)
!      ENDDO
!
!      POTPSC_CHECK = 0.0
!      CALL POTTORHO( PP%ZVALF_ORIG, NPSPTS, POTPS_G, PP%PSGMAX/NPSPTS, &
!     &            .TRUE. , PP%R%NMAX, PP%R%R ,  POTPSC_CHECK )                        
!
!      DO j=1, PP%R%NMAX
!         WRITE(22,'(6f20.8)') PP%R%R(j), PP%POTPSC(j), POTPSC_CHECK(j)
!      ENDDO
!   -----------------------------------------------------------------------------------------    !

!   ---------------- !!!!!! POT IN RECIPROCAL SPACE FROM POTPSC !!!!!  ----------------------    !     

!!      CALL FOURPOT_TO_Q( PP%RDEP, POT, PP%PSP(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)
      POTPS_G(:) = 0.0
!      POTPS_TEST(:) = POTPS_TEST(:)
!      DO i = 1, PP%R%NMAX
!          CALL RANDOM_SEED()
!          CALL RANDOM_NUMBER(random_x)
!          POTPSC_TEST(I) = -random_x
!          WRITE(6,*) 'TEST=', POTPS_TEST(I)
!      ENDDO

!      POTPSC_CHECK(:) = PP%POTPSC(:)
!      POTPSC_CHECK(:) = - POTPS_EFF(:)  !PP%POTPS(:)
!      POTPSC_CHECK(:) = POTPSC_TEST(:)
!      CALL FOURPOT_TO_Q_CHECK( PP%R%R(PP%R%NMAX), PP%ZVALF_ORIG, PP%POTPSC,   &
      CALL FOURPOT_TO_Q_CHECK( PP%R%R(PP%R%NMAX), PP%ZVALF_ORIG, POTPSC_TEST,   &
!      CALL FOURPOT_TO_Q_CHECK( PP%R%R(PP%R%NMAX), PP%ZVALF_ORIG, POTPSC_CHECK,   &
     &             POTPS_G, SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)
!      CALL FOURPOT_TO_Q( PP%R%R(PP%R%NMAX), POTPS_TEST, POTPS_G, SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)

      OPEN(UNIT=23,FILE='VASP_G_POTLOC',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=23,FILE='VASP_G_POTLOC',STATUS='OLD')
      ENDIF
      DO j=1, SIZE(PP%PSP,1)
         WRITE(IU19,'(6f20.8)') PP%PSP(j,1), PP%PSP(j,2), POTPS_G(j), PP%PSPCOR(j),  &
     &                          PP%PSPRHO(j), PP%PSPTAU(j)
      ENDDO
!      PP%PSP(:,2) = POTPS_G(:)
!   ---------------- !!!!!! FOR CHECK !!!!! -------------------
      
!      POTPSC_CHECK(:) = 0.0
!!       CALL POTTORHO( PP%ZVALF_ORIG, NPSPTS, PP%PSP(:,2), PP%PSGMAX/NPSPTS, &
!!     &            .TRUE. , PP%R%NMAX, PP%R%R ,  POTPSC_CHECK )                        
!       CALL POTTORHO( PP%ZVALF_ORIG, NPSPTS, POTPS_G, PP%PSGMAX/NPSPTS, &
!     &            .TRUE. , PP%R%NMAX, PP%R%R ,  POTPSC_CHECK )                        

      DO j=1, PP%R%NMAX

         WRITE(IU17,'(6f20.8)') PP%R%R(j), PP%POTPSC(j), POTPSC_TEST(j),  &
     &                         PP%POTPSC(j) + PP%POTPS(j), POTPSC_CHECK(j)
      ENDDO

!   ---------------- !!!!!! CORE-CHG IN RECIPROCAL SPACE FROM RHOPS !!!!!  ----------------------    !     
      CORPS_G(:) = 0.0
      DO j=1, PP%R%NMAX
         POTPSC_CHECK(j) = PP%RHOPS(j)/SCALE/PP%R%R(j)/PP%R%R(j)
      ENDDO
      CALL FOURPOT_TO_Q( PP%R%R(PP%R%NMAX), POTPSC_CHECK,   &
     &             CORPS_G, SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)

      OPEN(UNIT=25,FILE='VASP_G_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=25,FILE='VASP_G_PCORE',STATUS='OLD')
      ENDIF
      DO j=1, SIZE(PP%PSP,1)
         WRITE(IU21,'(6f20.8)') PP%PSP(j,1), PP%PSPCOR(j), -CORPS_G(j)
      ENDDO

      PP%PSPCOR(:) = -CORPS_G(:)

!   ---------------- !!!!!! FOR CHECK !!!!! -------------------
!      RHOPS00 = 0.d0
!       CALL POTTORHO( PP%ZVALF_ORIG+10, NPSPTS, PP%PSPCOR, PP%PSGMAX/NPSPTS, &
!     &            .FALSE. , PP%R%NMAX, PP%R%R ,  RHOPS00 )                        

!   ---------------- !!!!!! PSEUDO-CHG IN RECIPROCAL SPACE FROM  PSPRHO !!!!!  ----------------------    !     
      RHOPS_G(:) = 0.0
      RHOPS00(:) = RHOPS00/SCALE
      CALL FOURPOT_TO_Q( PP%R%R(PP%R%NMAX), RHOPS00,   &
     &             RHOPS_G, SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)

      OPEN(UNIT=27,FILE='VASP_G_PSPRHO',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=27,FILE='VASP_G_PSPRHO',STATUS='OLD')
      ENDIF
      DO j=1, SIZE(PP%PSP,1)
         QQ = (j-1)*PP%PSGMAX/SIZE(PP%PSP,1)
         WRITE(IU23,'(6f20.8)') PP%PSP(j,1), PP%PSPRHO(j), -RHOPS_G(j)/2.0
      ENDDO
      PP%PSPRHO(j) = -RHOPS_G(j)/2.0

!   ---------------- !!!!!! FOR CHECK !!!!! -------------------
!      DO j=1, PP%R%NMAX
!         WRITE(77,'(6f20.8)') PP%R%R(j), PP%RHOPS(j), RHOPS00(j)/SCALE
!      ENDDO
!      RHOPS00 = 0.d0
!       CALL POTTORHO( PP%ZVALF_ORIG, NPSPTS, PP%PSPRHO, PP%PSGMAX/NPSPTS, &
!     &            .FALSE. , PP%R%NMAX, PP%R%R ,  RHOPS00 )                        

      WRITE(6,*) 'VALUE=', PP%ZVALF_ORIG
      WRITE(6,*) 'NPSPTS=', NPSPTS
      WRITE(6,*) 'PSGMAX=', PP%PSGMAX
      WRITE(6,*) 'NMAX=', PP%R%NMAX

!   ---------------- !!!!!! DION !!!!!  ----------------------    !     
 
      ALLOCATE(TMP(PP%R%NMAX),DTMP(PP%R%NMAX))
      ALLOCATE(PARWKINAE(PP%R%NMAX, PP%LMAX), PARWKINPS(PP%R%NMAX, PP%LMAX))
      ALLOCATE(QPAW(PP%LMAX, PP%LMAX, 0:2*PP%LMAX))
      ALLOCATE(DIJ(PP%LMAX, PP%LMAX), DION(PP%LMAX, PP%LMAX))

!!    T|\psi_i\rangle = (-\nabla^2-\dfrac{l(l+1)}2) \psi_i     
      DO I=1,PP%LMAX
         CALL GRAD(PP%R,PP%WAE(:,I),TMP)
         CALL GRAD(PP%R,TMP,DTMP)
         DO J =1, PP%R%NMAX
            PARWKINAE(J,I) = -(DTMP(J)-PP%LPS(I)*(PP%LPS(I)+1)/PP%R%R(J)/PP%R%R(J)*PP%WAE(J,I))*HSQDTM
         ENDDO
!         WRITE(6,'(5F20.10)') (PP%R%R(J),PP%WAE(J,I),PP%WPS(J,I), & !PP%PARWKINAE(J,I),PP%PARWKINPS(J,I), &
!        &   -(DTMP(J)-PP%LPS(I)*(PP%LPS(I)+1)/PP%R%R(J)/PP%R%R(J)*PP%WAE(J,I))*HSQDTM, &
!        &   TMP(J),J=1,PP%R%NMAX)
      ENDDO
      
      DO I=1,PP%LMAX
         CALL GRAD(PP%R,PP%WPS(:,I),TMP)
         CALL GRAD(PP%R,TMP,DTMP)
         DO J =1, PP%R%NMAX
            PARWKINPS(J,I) = -(DTMP(J)-PP%LPS(I)*(PP%LPS(I)+1)/PP%R%R(J)/PP%R%R(J)*PP%WPS(J,I))*HSQDTM
         ENDDO
!         WRITE(6,'(5F20.10)') (PP%R%R(J),PP%WAE(J,I),PP%WPS(J,I), & !PP%PARWKINAE(J,I),PP%PARWKINPS(J,I), &
!        &   -(DTMP(J)-PP%LPS(I)*(PP%LPS(I)+1)/PP%R%R(J)/PP%R%R(J)*PP%WPS(J,I))*HSQDTM, &
!        &   TMP(J),J=1,PP%R%NMAX)
      ENDDO

!!  Dion_{ij}=\langle\phi_i|T+V_eff|\phi_j\rangle-\langle\tilde\phi_i|T+\tilde{V}_eff|\tilde\phi_j\rangle
!!      -\int\hat{Q}^{00}_{ij}(r)\tilde{V}(r)\mathrm{d}r

      !! DIJ=\langle\phi_i|T+V_eff|\phi_j\rangle-\langle\tilde\phi_i|T+\tilde{V}_eff|\tilde\phi_j\rangle
      DIJ = 0.d0
      DO CH1=1, PP%LMAX
      DO CH2=1, PP%LMAX
         TMP = 0.d0
         IF(PP%LPS(CH1) == PP%LPS(CH2)) THEN
                 DO I=1, PP%R%NMAX
                    TMP(I)=PP%WAE(I,CH1)*PARWKINAE(I,CH2)-PP%WAE(I,CH1)*PARWKINPS(I,CH2) &
     &                    +PP%WAE(I,CH1)*POTAE_EFF(I)*PP%WAE(I,CH2)-PP%WPS(I,CH1)*POTPS_EFF(I)*PP%WPS(I,CH2)
                 ENDDO
                 CALL SIMPI(PP%R, TMP, DIJ(CH1, CH2))
         ENDIF
      ENDDO
      ENDDO

!!   QPAW =\int(\phi_i\phi_i-\tilde\phi_i\tilde\phi_j)r^L\mathrm{d}r 
      DO CH1=1, PP%LMAX
      DO CH2=1, PP%LMAX
         LL=PP%LPS(CH1)
         LLP=PP%LPS(CH2)
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         DO LMAIN=LMIN, LMAX, 2
            TMP(:)=(PP%WAE(:,CH1)*PP%WAE(:,CH2)-      &
     &         PP%WPS(:,CH1)*PP%WPS(:,CH2))*PP%R%R(:)**LMAIN
            CALL SIMPI(PP%R, TMP, SUM)
            QPAW(CH1,CH2, LMAIN)=SUM
!            WRITE(6,*) 'QPAW:', PP%QPAW(CH1,CH2, LMAIN), SUM
         ENDDO
      ENDDO
      ENDDO

!!!! " unscreen "
!!      -\int\hat{Q}^{00}_{ij}(r)\tilde{V}(r)\mathrm{d}r

      LYMAX=MAXVAL(PP%LPS(1:PP%LMAX))
      IF (ASSOCIATED(PP%QPAW)) LYMAX=LYMAX*2
      LMMAX=(LYMAX+1)**2
      ALLOCATE(VTMP(PP%R%NMAX,LMMAX,1))
      ALLOCATE(DLM(PP%LMDIM*PP%LMDIM), DLLMM(PP%LMDIM,PP%LMDIM,1), DHXC(PP%LMDIM,PP%LMDIM))
!      SCALE=1/(2*SQRT(PI))      
      VTMP=0; VTMP(:,1,1)=POTPS_EFF(:)*SCALE
!      ! Reconstruct the PAW strength parameters for the reference system
      CALL RAD_POT_WEIGHT(PP%R,1,LYMAX,VTMP)
      DLM=0; DLLMM=0; DHXC=0
      CALL RAD_AUG_PROJ(VTMP(:,:,1),PP%R,DLM,PP%LMAX,PP%LPS,LYMAX,PP%AUG,QPAW)
      CALL TRANS_DLM(DLLMM(:,:,1),DLM,PP)

      DHXC=-DLLMM(:,:,1)
!      ! Compute D_{ij} and Q_{ij}
      LM=1
      DO CH1=1,PP%LMAX
      LMP=1
      DO CH2=1,PP%LMAX
         DIJ(CH1,CH2)=DIJ(CH1,CH2)+DHXC(LM,LMP)
!         QIJ(I,J)=PP%QION(I,J)
         LMP=LMP+2*PP%LPS(CH2)+1
      ENDDO
      LM=LM+2*PP%LPS(CH1)+1
      ENDDO      

!! Make DION hermitian
      DO CH1=1,PP%LMAX
      DO CH2=1,PP%LMAX
         DION(CH1,CH2)=(DIJ(CH1,CH2)+DIJ(CH2,CH1))/2.0/5.0/(floor((CH1-1)/2.0)+1)
         DION(CH2,CH1)=DION(CH1,CH2)
!         WRITE(6,*) 'DION:', DION(CH1,CH2), PP%DION(CH1,CH2)!, DION(CH1,CH2)/PP%DION(CH1,CH2)
         PP%DION(CH1, CH2) = DION(CH1, CH2)
      ENDDO      
      ENDDO      

!      CALL WRITE_POTCAR(INFO, PP)
!      STOP

      DEALLOCATE(VTMP, DLLMM, DHXC)
      DEALLOCATE(TMP, DTMP, PARWKINAE, PARWKINPS, QPAW)
      DEALLOCATE(DIJ, DION)

!   ---------------- !!!!!! PROJECTOR IN REAL SPACE FROM  PROJECTOR !!!!!  ----------------------    !     
      ALLOCATE(PSPNL(NPSNL), PSPRNL(NPSRNL))
!      OPEN(UNIT=29,FILE='VASP_PSPRNL',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=29,FILE='VASP_PSPRNL',STATUS='OLD')
!      ENDIF
      PP%R%NMAX = NPSRNL
      PR%NMAX = NPSRNL
      ALLOCATE(PR%R(PR%NMAX))
      DO i=1, 4
      RHOPS_G = 0.0
         DO j=1, NPSRNL
            PP%R%R(j)=PP%PSPRNL(j,1,i)
            PR%R(j)=PP%PSPRNL(j,1,i)
            PSPRNL(j)=PP%PSPRNL(j,2,i)
!            WRITE(IU25,'(6f20.8)') PP%PSPRNL(j,1,i), PP%PSPRNL(j,2,i)
            WRITE(75,'(6f20.8)') PP%R%R(j), PSPRNL(j)
         ENDDO

         CALL FOURPOT_TO_Q( PP%R%R(NPSRNL), PSPRNL,   &
!         CALL FOURPOT_TO_Q( PR%R(NPSRNL), PSPRNL,   &
     &             RHOPS_G, NPSRNL, PP%PSMAXN/NPSRNL, PP%R, IU6)
!     &             RHOPS_G, NPSRNL, PP%PSMAXN/NPSRNL, PR, IU6)
         WRITE(IU25,*)
!         DO j=1, NPSNL
!            QQ=PP%PSMAXN/NPSNL*(j-1)
!            WRITE(78,'(6f20.8)')QQ, RHOPS_G(j)
!         ENDDO
!         WRITE(78,*)
      ENDDO

!             WRITE(IU19,*)
!             WRITE(IU19,*) PP%PSDMAX
!             WRITE(IU19,*)

!      OPEN(UNIT=31,FILE='VASP_PSPNL',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=31,FILE='VASP_PSPNL',STATUS='OLD')
!      ENDIF
      DO i=1, 4
         DO j=1, NPSNL
            QQ=PP%PSMAXN/NPSNL*(j-1)
            WRITE(76,*) QQ, j, PP%PSPNL(j,i)
         ENDDO
         WRITE(76,*)
      ENDDO

!             WRITE(IU19,*)
      DEALLOCATE(PR%R, PSPNL, PSPRNL)
!      STOP

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

      Q_00 = integrator(Grid, FC%coreden, 1, Grid%n)
      Q_v = Pot_AE%q - Q_00
      paw%chag_val = Q_v
!      WRITE(6,*) 'Q_00=', Q_00

      OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='OLD')
      ENDIF
!!!   FC%coreden = n(r)*r^2  ==> consider the UNIT_EXCHANGE
!!!   PP%RHOAE = AUTOA * FC%coreden/(SCALE *AUTOA*AUTOA)
      DO j=1,Grid%n 
!         WRITE(IU10,'(6f20.8)') Grid%r(j)*AUTOA, FC%coreden(j)/(SCALE*AUTOA)/(SCALE*AUTOA)
         WRITE(IU10,'(6f20.8)') Grid%r(j)*AUTOA, FC%coreden(j)/(SCALE*AUTOA)
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
             WRITE(IU8,'(6f20.8)') Grid%r(j)*AUTOA, PAW%phi(j,i)/SQRT(AUTOA)*PAW%wfrate(i), PAW%phi(j,i)*PAW%wfrate(i)   &
     &                  , PAW%ophi(j,i)        
         ENDDO
         WRITE(IU8,*)
      ENDDO

      OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='OLD')
      ENDIF
          WRITE(6,*) 'AEPOT calcualted'

!!!!!!!!! POTAEC = V_H[n_Zc] = V_Z + V_H[n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          den = 0.d0
          call SetPOT(Grid, FC%coreden, den, PotHr, PotXCr, .TRUE.)
!!!!      VASP_V_Z = - Pot_AE%rvn(j)/Grid%r(j)*RYTOEV      !!!!  FOR the factor FELECT
!!!!      VASP_V_H[n_c] = PotHr(j)/Grid%r(j)*RYTOEV

          PotAECr(:) =  PotHr(:) + Pot_AE%rvn(:)                     !!!!  V_Z + V_H[n_c] 
!           PotAECr(:) = PotHr(:)                                    !!!!   V_H[n_c]

!!!!!!!!! POTAE = \Delta V_Z* + V_H[n_v] + V_XC[n_c+n_v] !!!!!!!!!!!!!!!!!!!!!!!!!!!
          call SetPOT(Grid, FC%coreden, FC%valeden, PotHr, PotXCr, .FALSE., .TRUE.)
!!!!      VASP_POTAE = -POTAEr(j)/Grid%r(j)*RYTOEV
          PotAEr(:) = PotHr(:)+PotXCr(:)

!!!!!!!!!!!!!!!!!   POTAE_EFF = V_H[n_Zc] + V_H[n_v]+ V_XC[n_v+n_c]  !!!!!!!!!!!!!!!!!!
!!!!!!!!! POTAE_EFF ===> Pot_AE%rv
!          call SetPOT(Grid, FC%coreden, FC%valeden, PotHr, PotXCr, .TRUE., .TRUE.)
!          PotATr(:) = PotHr(:)+PotXCr(:)
          PotATr(:) = -( PotAECr(:) + PotAEr(:))

          DO j=1,Grid%n 
             WRITE(IU16,'(8f20.8)') Grid%r(j)*AUTOA, - PotAEr(j)/Grid%r(j)*RYTOEV,   &  !!!  VASP_POT_AE
     &                                               - PotAECr(j)/Grid%r(j)*RYTOEV , &  !!!  VASP_POT_AEC 
     &                       PotATr(j)/Grid%r(j)*RYTOEV !, Pot_AE%rv(j)/Grid%r(j)*RYTOEV !!!  POTAE_EFF
          ENDDO

!!!!!!!!! tcore_den = sum_i B_i*sin(q_i r)/r !!!!!!!!!!!!!!!!!     
!         Grid%r(:)=Grid%r(:)*AUTOA
!         FC%coreden(:)=FC%coreden(:)/(SCALE*AUTOA)
         Call setcoretail2(Grid, FC%coreden)
         OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
            OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='OLD')
         ENDIF
!!!   FC%tcore = \tilde{n}(r)*r^2  ==> consider the UNIT_EXCHANGE
!!!   PP%RHOPS = AUTOAE * FC%tcore/(SCALE *AUTOA*AUTOA)
         DO j=1,Grid%n 
            WRITE(IU12,'(6f20.8)') Grid%r(j)*AUTOA, PAW%tcore(j)/(SCALE*AUTOA)
         ENDDO
!            irc_core= FindGridIndex(Grid, PAW%rc_core)
!            Q_00 = integrator(Grid, PAW%tcore, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!         Grid%r(:)=Grid%r(:)/AUTOA
!         FC%coreden(:) = FC%coreden(:)*(SCALE*AUTOA)
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
               WRITE(IU13,'(6f20.8)') Grid%r(j)*AUTOA, PAW%tphi(j,i)/sqrt(AUTOA)*PAW%wfrate(i), PAW%tphi(j,i)*PAW%wfrate(i)   &
     &                        , PAW%otphi(j,i)
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
!!!!!!!!!!!!!!!!!   POTAE_EFF = V_H[n_v] + V_H[n_Zc] +V_XC[n_v+n_c]  !!!!!!!!!!!!!!!!!!
         DO i = 1, Grid%n
            Pot_eff(i) =  PotATr(i)/Grid%r(i)
         ENDDO

!!!!!!!!!!!!!!!!!!!! POT_tVeff FROM POT_EFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
         R_CUT = PP%R%R(PP%R%NMAX-44)/AUTOA
         call SetPOT_TEFF(Grid, Pot_eff, Pot_teff, R_CUT)
         
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

         OPEN(UNIT=22,FILE='ATOM_POTPS',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
             OPEN(UNIT=22,FILE='ATOM_POTPS',STATUS='OLD')
         ENDIF
         DO j=1,Grid%n 
            WRITE(IU18,'(6f20.8)') Grid%r(j)*AUTOA, -PotPS(j)*RYTOEV, -Pot_eff(j)*RYTOEV, -Pot_teff(j)*RYTOEV
         ENDDO



!         STOP   !! POTPS_EFF TEST CORRECT

!!!!!!!!!!!!!  POTPSC = POTPS_EFF - (V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c])  !!!!!!!!!!!!!!!!!!!!!!!!
         DO j =1, Grid%n
            PotPSC(j)= -Pot_teff(j) - PotPS(j)
         ENDDO

         OPEN(UNIT=24,FILE='ATOM_POTPSC',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
             OPEN(UNIT=24,FILE='ATOM_POTPSC',STATUS='OLD')
         ENDIF

!         DO j=1, Grid%n
!            WRITE(IU20,'(6f20.8)') Grid%r(j)*AUTOA, PotPSC(j)*RYTOEV   
!         ENDDO

!!   ---------------- !!!!!! POT_V[\tilde{n}_Zc] FROM POT_V[n_Zc] !!!!!  ----------------------    !!     
         R_CUT = PP%R%R(PP%R%NMAX-40)/AUTOA
         call SetPOT_TEFF(Grid, PotAECr/Grid%r, PotPSC, R_CUT)
         DO j=1, Grid%n
            WRITE(IU20,'(6f20.8)') Grid%r(j)*AUTOA, PotPSC(j)*RYTOEV
         ENDDO

!   ---------------- !!!!!! POT IN RECIPROCAL SPACE FROM POT_V[n_Zc] !!!!!  ----------------------    !     
!        CALL FOURPOT_TO_Q()

!   ---------------- !!!!!! PROJECTOR IN RECIPROCAL SPACE ] !!!!!  ----------------------    !     
!        
         OPEN(UNIT=26,FILE='ATOM_PROJECT',STATUS='UNKNOWN',IOSTAT=IERR)
         DO i=1, PAW%nbase
           DO j=1, Grid%n
            WRITE(IU22,'(6f20.8)') Grid%r(j)*AUTOA, PAW%otp(j,i)!*SQRT(RYTOEV)
           ENDDO
         ENDDO


!        deallocate(coreden)

         CLOSE(IU5)
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

        CALL GENERATE_POTCARLIB(INFO, PP, CHANNELS)

        CALL GENERATE_POTCARDATA (INFO, PP, CHANNELS, Grid, FC%coreden, PAW, PotPS)

        WRITE(6, 300)
        CALL WRITE_POTCAR(INFO, PP, CHANNELS)

300     FORMAT(' !!! ----  GENERATE POTCAR  ---- !!! ',/)

         DEALLOCATE(RHO, V, RHOAE00, CRHODE, RHOLM)
         DEALLOCATE(POT, POTAE_EFF, DPOTAE_EFF, POTAEC)
         DEALLOCATE(PotHr, PotXCr, PotAEr, PotATr,PotAECr, PotPS, PotAE, PotAE00,  PotPSC)
         DEALLOCATE(pdensity, PotPSCr)
!        STOP
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
         ELSE
            irc = PAW%irc_vloc
         ENDIF
            rc = Grid2%r(irc)

         WRITE(6,*) 'For VASP: irc = ', irc, 'r_vloc=', rc

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

      SUBROUTINE GENERATE_POTCARDATA (INFO, FROM_PP, CHANNELS, Grid, COREDEN, PAW, POTPS)
         IMPLICIT COMPLEX(q)   (C)
         IMPLICIT REAL(q)  (A-B,D-H,O-Z)
!         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: FROM_PP
         TYPE(INFO_STRUCT) :: INFO
         TYPE(in_struct) IO
         TYPE(PseudoInfo) :: PAW
         Type(GridInfo) :: Grid

!      CHARACTER*40 SZNAMP         ! header
      CHARACTER*80 :: CSEL
      CHARACTER(LEN=4) :: PARAM_CHARACTER
      REAL(q)  ::   POTCAR_PARAM, COREDEN(Grid%n), POTAE(Grid%n), POTPS(Grid%n)
      REAL(q)  ::   POTCAR_G_PARAM, POTCAR_R_PARAM
!      REAL(q)  ::   WAE(Grid%n,PAW%nbase), WPS(Grid%n,PAW%nbase)
      REAL(q)  ::   WAE(Grid%n), WPS(Grid%n)
      INTEGER  ::   XC_TYPE, L1, L2, NMAX, NRANGE
      REAL(q), ALLOCATABLE :: POTCAR_DATA(:)
      REAL(q), ALLOCATABLE :: SPLINE_COEF(:,:), SPLINE_VALUE(:), SPLINE_DATA(:)
      LOGICAL  ::   PARAM_LOG

      INTEGER :: ifinput,ifen,Z, IU6, Kcut
      INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS, LMAX, LMMAX, LI, MI, LMI
      INTEGER :: LMAX_TABLE, LYMAX, NCUT 
      INTEGER :: CH1, CH2, LL, LLP, LMIN, LMAIN
      REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
      REAL(q) ROOT(2), QR, Qloc, R_CUT
      REAL(q) :: DHARTREE, QCORE,SCALE, DOUBLEAE, EXCG
!      REAL(q), ALLOCATABLE :: PotAE(:), PotAE00(:), PotPS(:), PotPSC(:), PotPSCr(:)
      REAL(q), ALLOCATABLE :: POTAE_EFF(:), DPOTAE_EFF(:), POTPS_EFF(:)
      REAL(q), ALLOCATABLE :: POTPS_G(:), CORPS_G(:), RHOPS_G(:)
      REAL(q), ALLOCATABLE :: RHO(:,:,:), POT(:,:,:), V(:,:,:), RHOAE00(:), RHOPS00(:)
      REAL(q), ALLOCATABLE :: POTAEC(:), POTPSC_TEST(:), POT_TEST(:), POTAE_TEST(:), POTPS_TEST(:)
      REAL(q), ALLOCATABLE :: POTPSC_CHECK(:)
      REAL(q), ALLOCATABLE :: CRHODE(:,:)
      REAL(q), ALLOCATABLE :: RHOLM(:), DLM(:), DLLMM(:,:,:)!, GAUSSIAN(:)
      REAL(q), ALLOCATABLE :: DHXC(:,:), QPAW(:,:,:)
      REAL(q), ALLOCATABLE :: TMP(:),DTMP(:),VTMP(:,:,:), DIJ(:,:), DION(:,:)
      REAL(q), ALLOCATABLE :: PARWKINAE(:,:), PARWKINPS(:,:)
      CHARACTER(LEN=2) :: TYPE(1)
      LOGICAL ::   LPAW, UNSCREN, LXCADD
!      INTEGER, EXTERNAL :: MAXL1
!      LOGICAL, OPTIONAL :: success
      LOGICAL :: success
!      INTEGER, OPTIONAL :: IMESH
      INTEGER :: IMESH
      REAL(q), EXTERNAL :: ERRF

      IU6 = -1
      LDIM = 8
      LDIM2 = (LDIM*(LDIM+1))/2
      LMDIM = 64
      ALLOCATE(POTCAR_DATA(NPSPTS))
      ALLOCATE(SPLINE_VALUE(FROM_PP%R%NMAX))
      ALLOCATE(SPLINE_DATA(Grid%n))
      ALLOCATE(SPLINE_COEF(1:3, 1:Grid%n))
      ALLOCATE(QPAW(FROM_PP%LMAX, FROM_PP%LMAX, 0:2*FROM_PP%LMAX))
      ALLOCATE(TMP(FROM_PP%R%NMAX))

!!!!!!!!_GENERATE_DATA_VASP_RHOAE_!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      DO I = 1, Grid%n
         COREDEN(I) = PAW%core(I)*AUTOA
      ENDDO
      
!      DO I=1, FROM_PP%R%NMAX
!         IF (FROM_PP%R%R(I) .GT. Grid%r(PAW%irc)*AUTOA) EXIT
!      ENDDO
!      NCUT = I+1

      CALL SPLINE(Grid%r*AUTOA, COREDEN, SPLINE_COEF(1,1:Grid%n),  &
     &       SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
      DO I=1, FROM_PP%R%NMAX
          SPLINE_VALUE(I) = ISPLINE(FROM_PP%R%R(I), Grid%r*AUTOA, COREDEN,  &
     &      SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
          FROM_PP%RHOAE(I) = SPLINE_VALUE(I)
      ENDDO
!      DO I=1, FROM_PP%R%NMAX
!         WRITE(96,*) FROM_PP%R%R(I), FROM_PP%RHOAE(I)
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

!!!!!!!!_GENERATE_DATA_VASP_RHOPS_!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      DO I = 1, Grid%n
         COREDEN(I) = PAW%tcore(I)*AUTOA
      ENDDO
      CALL SPLINE(Grid%r*AUTOA, COREDEN, SPLINE_COEF(1,1:Grid%n),  &
     &       SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
      DO I=1, FROM_PP%R%NMAX
          SPLINE_VALUE(I) = ISPLINE(FROM_PP%R%R(I), Grid%r*AUTOA, COREDEN,  &
     &      SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
          FROM_PP%RHOPS(I) = SPLINE_VALUE(I)
      ENDDO
!      DO I=1, FROM_PP%R%NMAX
!         WRITE(97,*) FROM_PP%R%R(I), FROM_PP%RHOPS(I)
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

!!!!!!!!_GENERATE_DATA_WAVE-Function_!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      DO I = 1, PAW%nbase
         !!!!!   pseudo  wave-function   !!!!!
           DO K = 1, Grid%n
              WPS(K) = PAW%tphi(K,I)/SQRT(AUTOA)*PAW%wfrate(I)
           ENDDO
           
           DO K = Grid%n, 1, -1
              IF (PAW%tphi(K,I) .GT. 1.0E-6) EXIT
           ENDDO
           Kcut = K

           DO K = FROM_PP%R%NMAX, 1, -1
              IF (FROM_PP%WPS(K,I) .GT. 1.0E-6) EXIT
           ENDDO
           NCUT = K
           
           SPLINE_VALUE = 0.d0
           CALL SPLINE(Grid%r*AUTOA, WPS, SPLINE_COEF(1,1:Kcut),  &
     &             SPLINE_COEF(2,1:Kcut), SPLINE_COEF(3,1:Kcut), Kcut)
           DO K=1, NCUT
                  SPLINE_VALUE(K) = ISPLINE(FROM_PP%R%R(K), Grid%r*AUTOA, WPS,  &
     &           SPLINE_COEF(1,1:Kcut),  SPLINE_COEF(2,1:Kcut), SPLINE_COEF(3,1:Kcut), Kcut)
                FROM_PP%WPS(K,I) = SPLINE_VALUE(K)
           ENDDO

!           DO K=1, FROM_PP%R%NMAX
!                          WRITE(98,*) FROM_PP%R%R(K), FROM_PP%WPS(K,I), SPLINE_VALUE(K)
!           ENDDO
!                          WRITE(98,*) 
        !!!!!!   ae  wave-function   !!!!!
           DO K = 1, Grid%n
              WAE(K) = PAW%phi(K,I)/SQRT(AUTOA)*PAW%wfrate(I)
           ENDDO
           CALL SPLINE(Grid%r*AUTOA, WAE, SPLINE_COEF(1,1:Kcut),  &
     &             SPLINE_COEF(2,1:Kcut), SPLINE_COEF(3,1:Kcut), Kcut)
           DO K=1, NCUT
                  SPLINE_VALUE(K) = ISPLINE(FROM_PP%R%R(K), Grid%r*AUTOA, WAE,  &
     &           SPLINE_COEF(1,1:Kcut),  SPLINE_COEF(2,1:Kcut), SPLINE_COEF(3,1:Kcut), Kcut)
                 FROM_PP%WAE(K,I) = SPLINE_VALUE(K)
!                          WRITE(99,*) FROM_PP%R%R(K), FROM_PP%WAE(K,I), SPLINE_VALUE(K)
           ENDDO

!           DO K=1, FROM_PP%R%NMAX
!                          WRITE(99,*) FROM_PP%R%R(K), FROM_PP%WAE(K,I), SPLINE_VALUE(K)
!           ENDDO
!                          WRITE(99,*) 
      ENDDO
!      STOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
       FROM_PP%ZVALF_ORIG = paw%chag_val

       SCALE = 2.0*sqrt(PI0)

       CALL SET_SIMP(FROM_PP%R)

       QPAW = 0.d0
!!   QPAW =\int(\phi_i\phi_i-\tilde\phi_i\tilde\phi_j)r^L\mathrm{d}r 
      DO CH1=1, FROM_PP%LMAX
      DO CH2=1, FROM_PP%LMAX
         LL= FROM_PP%LPS(CH1)
         LLP= FROM_PP%LPS(CH2)
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         DO LMAIN=LMIN, LMAX, 2
            TMP(:)=(FROM_PP%WAE(:,CH1)*FROM_PP%WAE(:,CH2)-      &
     &         FROM_PP%WPS(:,CH1)*FROM_PP%WPS(:,CH2))*FROM_PP%R%R(:)**LMAIN
            CALL SIMPI(FROM_PP%R, TMP, SUM)
            QPAW(CH1,CH2, LMAIN)=SUM
            FROM_PP%QPAW(CH1,CH2, LMAIN)=SUM
!           WRITE(6,*) 'TEST_QPAW:', FROM_PP%QPAW(CH1,CH2, LMAIN)
         ENDDO
      ENDDO
      ENDDO


       ALLOCATE(RHOPS00(FROM_PP%R%NMAX), RHOAE00(FROM_PP%R%NMAX))

      CALL SIMPI(FROM_PP%R,FROM_PP%RHOAE, QCORE)
!      WRITE(6,*) 'QCORE=', QCORE, QCORE*SCALE
!      CALL SIMPI(FROM_PP%R,FROM_PP%RHOPS, QTEST)
!      WRITE(6,*) 'QCORE=', QTEST, QTEST*SCALE

      LMMAX = (FROM_PP%LMAX+1)**2
      LYMAX =MAXL1(FROM_PP)*2


      FROM_PP%LMAX_CALC=LMAX
      
      CALL SET_PAW_AUG(NTYPD, FROM_PP, IU6, FROM_PP%LMAX_CALC, .TRUE.)
!        WRITE(*,*) 'ZVAL=', FROM_ZVALF_ORIG
!      WRITE(6,*) 'LYMAX=', LYMAX
      ALLOCATE(RHO(FROM_PP%R%NMAX, LMMAX,1), V(FROM_PP%R%NMAX, LMMAX,1))
      ALLOCATE(POT(FROM_PP%R%NMAX, LMMAX,1), POTAEC(FROM_PP%R%NMAX))
      ALLOCATE(POTAE_EFF(FROM_PP%R%NMAX),DPOTAE_EFF(FROM_PP%R%NMAX), POTPS_EFF(FROM_PP%R%NMAX) )
      ALLOCATE(POT_TEST(FROM_PP%R%NMAX),POTAE_TEST(FROM_PP%R%NMAX), POTPSC_TEST(FROM_PP%R%NMAX) )
      ALLOCATE(POTPSC_CHECK(FROM_PP%R%NMAX))
      ALLOCATE(POTPS_G(NPSPTS), CORPS_G(NPSPTS), RHOPS_G(NPSPTS))
      ALLOCATE(POTPS_TEST(FROM_PP%R%NMAX))
!      ALLOCATE(V1(PP%R%NMAX, LMMAX,1), V2(PP%R%NMAX, LMMAX,1))
      ALLOCATE(CRHODE(LMDIM,LMDIM))
!      ALLOCATE(CRHODE(LDIM,LDIM), POT-AE(PP%R%NMAX))
      ALLOCATE(RHOLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
!      ALLOCATE(DLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
!
!
      LMAX_TABLE=6; CALL YLM3ST_(LMAX_TABLE)
!      DO j=1, PP%R%NMAX
!         GAUSSIAN(j) = exp(-PP%R%R(j)**2/PP%R%R(1)**2/10000)
!         WRITE(6,*) 'GAUSSIAN=', GAUSSIAN(j)
!      ENDDO 
! CRHODE is n_ij (n_occ)    # 价电子轨道的占据数(展开) n_lm,l'm'
      CALL SET_CRHODE_ATOM(CRHODE,FROM_PP)
! RHOLM is f_LM            #  价电子占据数构造 f_LM
      CALL TRANS_RHOLM(CRHODE, RHOLM, FROM_PP)
! RHO is density of valence # 价电子密度 n_lm^v = f_LM\langle\psi_i|\psi_j\rangle
      RHO = 0
      CALL RAD_CHARGE(RHO(:,:,1), FROM_PP%R, RHOLM(:), FROM_PP%LMAX, FROM_PP%LPS, FROM_PP%WAE)
      RHOAE00 = 0
      DO i =1,CHANNELS
         LI = FROM_PP%LPS(i)
         DO MI = 1, 2*LI+1
            LMI = LI*LI + MI
!            WRITE(6,*) 'TRHOV=', RHO(FROM_PP%R%NMAX, LI+1,1)
            DO j =1, FROM_PP%R%NMAX
             RHOAE00(j) = RHOAE00(j)+RHO(j,LMI,1)
            ENDDO
         ENDDO
      ENDDO
!      CALL SIMPI(PP%R,RHOAE00+PP%RHOAE*SCALE, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST, FROM_PP%ZVALF_ORIG

!!!!!!!!! POTAE = V_H[n_v]+V_XC[n_v+n_c] != POTAEC  !!!!!!!!!     
!      POT = 0
!      POTAEC = 0
!      RHO = 0
!      CALL PUSH_XC_TYPE(FROM_PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
!      CALL RAD_POT(FROM_PP%R, 1, 1, 1, .FALSE., &
!     &            RHO, FROM_PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
!      CALL POP_XC_TYPE
!      POTAE_TEST(:) = POT(:,1,1)
!!!!!!!!!!!!!!!!!!!!!!!! POTAEC = V_H[n_Zc] = V_H[n_c]+ Z/r !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Z=ANINT(FROM_PP%ZVALF_ORIG+QCORE*SCALE)
!      Z=INT(QCORE)

      CALL RAD_POT_HAR(0, FROM_PP%R, POTAEC, FROM_PP%RHOAE, DHARTEE)
      DO j=1, FROM_PP%R%NMAX
!         WRITE(6,*) 'TEST =', Z/PP%R%R(j)*FELECT
         POTAEC(j) =   POTAEC(j)/SCALE - Z/FROM_PP%R%R(j)*FELECT  !!!   V_H[n_Zc] =  V_H[n_c] - Z/r 
!         POTAEC(j) =   POTAEC(j)/SCALE - (Z-PP%ZVALF_ORIG)/PP%R%R(j)*FELECT  !!!   V_H[n_Zc] =  V_H[n_c] - Z/r 
      ENDDO
!      DO j =1, FROM_PP%R%NMAX
!         WRITE(95,*) FROM_PP%R%R(j), POTAEC(j)
!      ENDDO
!      WRITE(6,*) !'THE_NUCLEAR_CHARGE: ',Z

!!!!!!!!!!!!!!!!!!!!!!!! POTAE_TEST =  V_H[n_v] +V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!
      POT = 0
      POT_TEST = 0
      CALL PUSH_XC_TYPE(FROM_PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(FROM_PP%R, 1, 1, 1, .FALSE., &
     &            RHO, FROM_PP%RHOAE, POT_TEST, POT, DOUBLEAE, EXCG)
      CALL POP_XC_TYPE
      FROM_PP%POTAE(:) = - POT(:,1,1)/SCALE                       !!!    POTAE = V_H[n_v] +V_XC[n_v+n_c] 
!      POTAE_TEST(:) =  -POT(:,1,1)/SCALE + PP%ZVALF_ORIG/PP%R%R(:)/SCALE/SQRT(2.0)
!      DO j=1, PP%R%NMAX
!         POTAE_TEST(j) = - POT(j,1,1)/SCALE + FELECT*Z*(1.0-ERRF(PP%R%R(j)/2.0/AUTOA))/PP%R%R(j)
!         WRITE(6,'(8f20.8)') PP%R%R(j), POT(j,1,1)/SCALE, ERRF(PP%R%R(j)/AUTOA), &
!     &             PP%POTAE(j)+POT(j,1,1)/SCALE, FELECT*Z*ERRF(PP%R%R(j)/PP%R%R(200))/PP%R%R(j) , &
!     &             FELECT*Z/PP%R%R(j)/(PP%POTAE(j)+POT(j,1,1)/SCALE )
!      ENDDO
!      DO j =1, FROM_PP%R%NMAX
!         WRITE(95,*) FROM_PP%R%R(j), -POTAE_TEST(j)
!      ENDDO
!      STOP

!!!!!!!!!!!!!!!!!   POTAE_EFF = V_H[n_Zc] + V_H[n_v] +V_XC[n_v+n_c]  !!!!!!!!!!!!!!!!!!
      POT = 0
      CALL PUSH_XC_TYPE(FROM_PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(FROM_PP%R, 1, 1, 1, .FALSE., &
     &            RHO, FROM_PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
      CALL POP_XC_TYPE
      POTAE_EFF(:) =  - POT(:,1,1)/SCALE
!      DO i =1,CHANNELS
!         LI = FROM_PP%LPS(i)
!         DO MI = 1, 2*LI+1
!            LMI = LI*LI + MI
!            POTAE_EFF(:) =  POTAE_EFF(:) - POT(:,LMI,1)
!        ENDDO
!      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!! POTPS_EFF = A*sin(qloc*r)/r !!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 1, FROM_PP%R%NMAX
         POTPS_EFF(j) = POTAE_EFF(j)
      ENDDO
      CALL GRAD(FROM_PP%R, POTAE_EFF, DPOTAE_EFF)

      IMESH = 41
!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX)/POTAEC(PP%R%NMAX)
      alpha = 1.0-DPOTAE_EFF(FROM_PP%R%NMAX-41)/POTAE_EFF(FROM_PP%R%NMAX-41)
!      WRITE(6,*) 'alpha= ', DPOTAE_EFF(PP%R%NMAX), POTAE_EFF(PP%R%NMAX), alpha
      beta = 1.0

      CALL SOLVEBESL_Q(ROOT, alpha, beta, 0, 1)
!      WRITE(6,*) 'ROOT=' , ROOT(1), ROOT(2)
!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX)
      QR = ROOT(1); Qloc = QR/FROM_PP%R%R(FROM_PP%R%NMAX-IMESH)
!      alpha = POTAE_EFF(PP%R%NMAX)*PP%R%R(PP%R%NMAX)/sin(QR)
      alpha = POTAE_EFF(FROM_PP%R%NMAX-IMESH)*FROM_PP%R%R(FROM_PP%R%NMAX-IMESH)/sin(QR)
!      WRITE(6,*) 'alpha=' , alpha
      DO j = 1, FROM_PP%R%NMAX-IMESH
         QR = Qloc*FROM_PP%R%R(j)
         POTPS_EFF(j) = alpha*sin(QR)/FROM_PP%R%R(j)
!         WRITE(6,*) PP%R%R(j), POTPS_EFF(j), POTAE_EFF(j)
!         WRITE(6,'(5f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j), PP%POTAE(j), PP%POTPS(j)
      ENDDO

      RHO = 0
      POTPSC_TEST = 0
      POTPSC_CHECK = 0
      POT_TEST = 0
!!!!!!!!!!!!!  POTPS = V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c]  !!!!!!!!!!!!!!!!!!!!!!!!
! RHO is \tilde RHO + \hat n
      CALL RAD_CHARGE(RHO(:,:,1), FROM_PP%R, RHOLM(:), FROM_PP%LMAX, FROM_PP%LPS, FROM_PP%WPS)
      CALL RAD_AUG_CHARGE(  RHO(:,:,1), FROM_PP%R, RHOLM, FROM_PP%LMAX, FROM_PP%LPS,  &
              LYMAX, FROM_PP%AUG, FROM_PP%QPAW )

      RHOPS00 = 0
      DO i =1,CHANNELS
         LI = FROM_PP%LPS(i)
         DO MI = 1, 2*LI+1
            LMI = LI*LI + MI
!           WRITE(6,*) 'TRHOV=', RHO(PP%R%NMAX, LI+1,1)
            RHOPS00(:) = RHOPS00(:)+RHO(:,LMI,1)
         ENDDO
      ENDDO
!      CALL SIMPI(PP%R,RHOPS00, QTEST)
      WRITE(6,*) 'CHECK_POTCAR_DATA_GENERATE'
      POT = 0
      CALL PUSH_XC_TYPE(FROM_PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(FROM_PP%R, 1, 1, 1, .FALSE., &
     &            RHO, FROM_PP%RHOPS, POTPSC_TEST, POT, DOUBLEAE, EXCG)
!      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
!                   RHO, PP%RHOPS, PP%POTPSC, POT, DOUBLEAE, EXCG)
      CALL POP_XC_TYPE
      POTPS_TEST(:) =  POT(:,1,1)/SCALE

      DO I = 1, Grid%n
         POTPS(I) = -POTPS(I)*RYTOEV
      ENDDO
      CALL SPLINE(Grid%r*AUTOA, POTPS, SPLINE_COEF(1,1:Grid%n),  &
     &       SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
      DO I=1, FROM_PP%R%NMAX
          SPLINE_VALUE(I) = ISPLINE(FROM_PP%R%R(I), Grid%r*AUTOA, POTPS,  &
     &      SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
          FROM_PP%POTPS(I) = SPLINE_VALUE(I)
      ENDDO
!      DO j =1, FROM_PP%R%NMAX
!         WRITE(95,*) FROM_PP%R%R(j), FROM_PP%POTPS(j)
!      ENDDO
!      STOP
!!!!!!!!!!!!!  Method_1:   POTPSC = POTPS_EFF - (V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c])  !!!!!!!!!!!!!!!!!!!!!!!!
      POTPSC_TEST(:) =  - POTPS_EFF(:) - POTPS_TEST(:)
      FROM_PP%POTPSC(:) = POTPSC_TEST(:)

!      DO j =1, FROM_PP%R%NMAX
!         WRITE(95,*) FROM_PP%R%R(j), POTPSC_TEST(j)
!      ENDDO
!      STOP
!!!!!!!!!!!!   Method_2:   POTPSC = A*sin(qloc*r)/r FROM POTAEC DIRACTLY     !!!!!!!!!!!!!!!!!!!!!!!!
      DO j = 1, FROM_PP%R%NMAX
         POTPSC_CHECK(j) = POTAEC(j)
      ENDDO
      CALL GRAD(FROM_PP%R, POTAEC, DPOTAE_EFF)

      IMESH = 38
!      alpha = 1.0-DPOTAE_EFF(PP%R%NMAX)/POTAEC(PP%R%NMAX)
      alpha = 1.0-DPOTAE_EFF(FROM_PP%R%NMAX-IMESH)/POTAEC(FROM_PP%R%NMAX-IMESH)
!      WRITE(6,*) 'alpha= ', DPOTAE_EFF(PP%R%NMAX), POTAE_EFF(PP%R%NMAX), alpha
      beta = 1.0

      CALL SOLVEBESL_Q(ROOT, alpha, beta, 0, 1)
!      WRITE(6,*) 'ROOT=' , ROOT(1), ROOT(2)
!      QR = ROOT(1); Qloc = QR/PP%R%R(PP%R%NMAX)
      QR = ROOT(1); Qloc = QR/FROM_PP%R%R(FROM_PP%R%NMAX-IMESH)
!      alpha = POTAE_EFF(PP%R%NMAX)*PP%R%R(PP%R%NMAX)/sin(QR)
      alpha = POTAEC(FROM_PP%R%NMAX-IMESH)*FROM_PP%R%R(FROM_PP%R%NMAX-IMESH)/sin(QR)
!      WRITE(6,*) 'alpha=' , alpha
      DO j = 1, FROM_PP%R%NMAX-IMESH
         QR = Qloc*FROM_PP%R%R(j)
         POTPSC_CHECK(j) = alpha*sin(QR)/FROM_PP%R%R(j)
!         WRITE(6,*) PP%R%R(j), POTPS_EFF(j), POTAE_EFF(j)
!         WRITE(6,'(5f20.8)') PP%R%R(j), PP%POTPSC(j), POTPS_EFF(j), PP%POTAE(j), PP%POTPS(j)
      ENDDO

!      DO j =1, FROM_PP%R%NMAX
!         WRITE(95,*) FROM_PP%R%R(j), POTPSC_TEST(j), POTPSC_CHECK(j)
!      ENDDO
!      STOP
!   ---------------- !!!!!! POT IN RECIPROCAL SPACE FROM POTPS !!!!!  ----------------------    !     

!!      CALL FOURPOT_TO_Q( PP%RDEP, POT, PP%PSP(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)
      POTPS_G(:) = 0.0
!      POTPS_TEST(:) = POTPS_TEST(:)
!      DO i = 1, PP%R%NMAX
!          CALL RANDOM_SEED()
!          CALL RANDOM_NUMBER(random_x)
!          POTPSC_TEST(I) = -random_x
!          WRITE(6,*) 'TEST=', POTPS_TEST(I)
!      ENDDO

!      POTPSC_CHECK(:) = FROM_PP%POTPSC(:)
!      POTPSC_CHECK(:) =  POTPSC_TEST(:)
!      CALL FOURPOT_TO_Q_CHECK( FROM_PP%R%R(FROM_PP%R%NMAX), FROM_PP%ZVALF_ORIG, POTPSC_TEST,   &
      CALL FOURPOT_TO_Q_CHECK( FROM_PP%R%R(FROM_PP%R%NMAX), FROM_PP%ZVALF_ORIG, POTPSC_CHECK,   &
     &             POTPS_G, SIZE(FROM_PP%PSP,1), FROM_PP%PSGMAX/SIZE(FROM_PP%PSP,1), FROM_PP%R, IU6)


!      CALL FOURPOT_TO_Q( PP%R%R(PP%R%NMAX), POTPS_TEST, POTPS_G, SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)

      OPEN(UNIT=33,FILE='ATOM_G_POTLOC',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=33,FILE='ATOM_G_POTLOC',STATUS='OLD')
      ENDIF
      DO j =1, SIZE(FROM_PP%PSP,1)
         WRITE(33,*) FROM_PP%PSP(j,1), POTPS_G(j), FROM_PP%PSP(j,2)
      ENDDO
      CLOSE(33)

      DO j=1, SIZE(FROM_PP%PSP,1)
         FROM_PP%PSP(j,2) = POTPS_G(j)!, FROM_PP%PSPCOR(j) ,  &
!     &                          FROM_PP%PSPRHO(j) = FROM_PP%PSPTAU(j)
      ENDDO

!      STOP

!   ---------------- !!!!!! FOR CHECK !!!!! -------------------
      
!      POTPSC_CHECK(:) = 0.0
!!       CALL POTTORHO( PP%ZVALF_ORIG, NPSPTS, PP%PSP(:,2), PP%PSGMAX/NPSPTS, &
!!     &            .TRUE. , PP%R%NMAX, PP%R%R ,  POTPSC_CHECK )                        
!       CALL POTTORHO( FROM_PP%ZVALF_ORIG, NPSPTS, POTPS_G, FROM_PP%PSGMAX/NPSPTS, &
!     &            .TRUE. , FROM_PP%R%NMAX, FROM_PP%R%R ,  POTPSC_CHECK )                        


!   ---------------- !!!!!! CORE-CHG IN RECIPROCAL SPACE FROM RHOPS !!!!!  ----------------------    !     
      CORPS_G(:) = 0.0
      CALL FOURPOT_TO_Q( FROM_PP%R%R(FROM_PP%R%NMAX), FROM_PP%RHOPS,   &
     &             CORPS_G, SIZE(FROM_PP%PSP,1), FROM_PP%PSGMAX/SIZE(FROM_PP%PSP,1), FROM_PP%R, IU6)
!      DO j =1, SIZE(FROM_PP%PSP,1)
!         WRITE(93,*) FROM_PP%PSP(j,1),  FROM_PP%PSPCOR(j), -CORPS_G(j)
!      ENDDO

      DO j=1, SIZE(FROM_PP%PSP,1)
         FROM_PP%PSPCOR(j) = -CORPS_G(j)
      ENDDO

!      STOP
!   ---------------- !!!!!! PSEUDO-CHG IN RECIPROCAL SPACE FROM  PSPRHO !!!!!  ----------------------    !     
      RHOPS_G(:) = 0.0
      CALL FOURPOT_TO_Q( FROM_PP%R%R(FROM_PP%R%NMAX), RHOPS00+FROM_PP%RHOPS,   &
     &             RHOPS_G, SIZE(FROM_PP%PSP,1),FROM_PP%PSGMAX/SIZE(FROM_PP%PSP,1), FROM_PP%R, IU6)

      OPEN(UNIT=37,FILE='ATOM_G_PSPRHO',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=37,FILE='ATOM_G_PSPRHO',STATUS='OLD')
      ENDIF
      DO j =1, SIZE(FROM_PP%PSP,1)
         WRITE(37,*) FROM_PP%PSP(j,1), -RHOPS_G(j)/4.0/SQRT(PI0), FROM_PP%PSPRHO(j)
      ENDDO
      CLOSE(37)

      DO j=1, SIZE(FROM_PP%PSP,1)
         FROM_PP%PSPRHO(j) = -RHOPS_G(j)/4.0/SQRT(PI0)
      ENDDO

!      STOP

!      WRITE(6,*) 'VALUE=', FROM_PP%ZVALF_ORIG
!      WRITE(6,*) 'NPSPTS=', NPSPTS
!      WRITE(6,*) 'PSGMAX=', FROM_PP%PSGMAX
!      WRITE(6,*) 'NMAX=', FROM_PP%R%NMAX
!   ---------------- !!!!!! DION !!!!!  ----------------------    !     
 
      ALLOCATE(DTMP(FROM_PP%R%NMAX))
      ALLOCATE(PARWKINAE(FROM_PP%R%NMAX, FROM_PP%LMAX), PARWKINPS(FROM_PP%R%NMAX, FROM_PP%LMAX))
      ALLOCATE(DIJ(FROM_PP%LMAX, FROM_PP%LMAX), DION(FROM_PP%LMAX, FROM_PP%LMAX))

!!    T|\psi_i\rangle = (-\nabla^2-\dfrac{l(l+1)}2) \psi_i     
      DO I=1,FROM_PP%LMAX
         CALL GRAD(FROM_PP%R,FROM_PP%WAE(:,I),TMP)
         CALL GRAD(FROM_PP%R,TMP,DTMP)
         DO J =1, FROM_PP%R%NMAX
            PARWKINAE(J,I) = -(DTMP(J)-FROM_PP%LPS(I)*(FROM_PP%LPS(I)+1)/FROM_PP%R%R(J)/FROM_PP%R%R(J)*FROM_PP%WAE(J,I))*HSQDTM
         ENDDO
!         WRITE(6,'(5F20.10)') (PP%R%R(J),PP%WAE(J,I),PP%WPS(J,I), & !PP%PARWKINAE(J,I),PP%PARWKINPS(J,I), &
!        &   -(DTMP(J)-PP%LPS(I)*(PP%LPS(I)+1)/PP%R%R(J)/PP%R%R(J)*PP%WAE(J,I))*HSQDTM, &
!        &   TMP(J),J=1,PP%R%NMAX)
      ENDDO
      
      DO I=1,FROM_PP%LMAX
         CALL GRAD(FROM_PP%R,FROM_PP%WPS(:,I),TMP)
         CALL GRAD(FROM_PP%R,TMP,DTMP)
         DO J =1, FROM_PP%R%NMAX
            PARWKINPS(J,I) = -(DTMP(J)-FROM_PP%LPS(I)*(FROM_PP%LPS(I)+1)/FROM_PP%R%R(J)/FROM_PP%R%R(J)*FROM_PP%WPS(J,I))*HSQDTM
         ENDDO
!         WRITE(6,'(5F20.10)') (PP%R%R(J),PP%WAE(J,I),PP%WPS(J,I), & !PP%PARWKINAE(J,I),PP%PARWKINPS(J,I), &
!        &   -(DTMP(J)-PP%LPS(I)*(PP%LPS(I)+1)/PP%R%R(J)/PP%R%R(J)*PP%WPS(J,I))*HSQDTM, &
!        &   TMP(J),J=1,PP%R%NMAX)
      ENDDO

!!  Dion_{ij}=\langle\phi_i|T+V_eff|\phi_j\rangle-\langle\tilde\phi_i|T+\tilde{V}_eff|\tilde\phi_j\rangle
!!      -\int\hat{Q}^{00}_{ij}(r)\tilde{V}(r)\mathrm{d}r

      !! DIJ=\langle\phi_i|T+V_eff|\phi_j\rangle-\langle\tilde\phi_i|T+\tilde{V}_eff|\tilde\phi_j\rangle
      DIJ = 0.d0
      DO CH1=1, FROM_PP%LMAX
      DO CH2=1, FROM_PP%LMAX
         TMP = 0.d0
         IF (FROM_PP%LPS(CH1) == FROM_PP%LPS(CH2)) THEN
                 DO I=1, FROM_PP%R%NMAX
                    TMP(I)= FROM_PP%WAE(I,CH1)*PARWKINAE(I,CH2)-FROM_PP%WAE(I,CH1)*PARWKINPS(I,CH2) &
     & +FROM_PP%WAE(I,CH1)*POTAE_EFF(I)*FROM_PP%WAE(I,CH2)-FROM_PP%WPS(I,CH1)*POTPS_EFF(I)*FROM_PP%WPS(I,CH2)
                 ENDDO
                 CALL SIMPI(FROM_PP%R, TMP, DIJ(CH1, CH2))
         ENDIF
      ENDDO
      ENDDO

!!!! " unscreen "
!!      -\int\hat{Q}^{00}_{ij}(r)\tilde{V}(r)\mathrm{d}r

      LYMAX=MAXVAL(FROM_PP%LPS(1:FROM_PP%LMAX))
      IF (ASSOCIATED(FROM_PP%QPAW)) LYMAX=LYMAX*2
      LMMAX=(LYMAX+1)**2
      ALLOCATE(VTMP(FROM_PP%R%NMAX,LMMAX,1))
      ALLOCATE(DLM(FROM_PP%LMDIM*FROM_PP%LMDIM), DLLMM(FROM_PP%LMDIM,FROM_PP%LMDIM,1), DHXC(FROM_PP%LMDIM,FROM_PP%LMDIM))
!      SCALE=1/(2*SQRT(PI))      
      VTMP=0; VTMP(:,1,1)=POTPS_EFF(:)*SCALE
!      ! Reconstruct the PAW strength parameters for the reference system
      CALL RAD_POT_WEIGHT(FROM_PP%R,1,LYMAX,VTMP)
      DLM=0; DLLMM=0; DHXC=0
      CALL RAD_AUG_PROJ(VTMP(:,:,1),FROM_PP%R,DLM,FROM_PP%LMAX,FROM_PP%LPS,LYMAX,FROM_PP%AUG,QPAW)
      CALL TRANS_DLM(DLLMM(:,:,1),DLM,FROM_PP)

      DHXC=-DLLMM(:,:,1)
!      ! Compute D_{ij} and Q_{ij}
      LM=1
      DO CH1=1,FROM_PP%LMAX
      LMP=1
      DO CH2=1,FROM_PP%LMAX
         DIJ(CH1,CH2)=DIJ(CH1,CH2)+DHXC(LM,LMP)
!         QIJ(I,J)=PP%QION(I,J)
         LMP=LMP+2*FROM_PP%LPS(CH2)+1
      ENDDO
      LM=LM+2*FROM_PP%LPS(CH1)+1
      ENDDO      

!! Make DION hermitian
      DO CH1=1,FROM_PP%LMAX
      DO CH2=1,FROM_PP%LMAX
         DION(CH1,CH2)=(DIJ(CH1,CH2)+DIJ(CH2,CH1))/2.0/5.0/(floor((CH1-1)/2.0)+1)
         DION(CH2,CH1)=DION(CH1,CH2)
!         WRITE(6,*) 'DION:', DION(CH1,CH2), FROM_PP%DION(CH1,CH2)!, DION(CH1,CH2)/PP%DION(CH1,CH2)
         FROM_PP%DION(CH1,CH2) = DION(CH1,CH2) !, DION(CH1,CH2)/PP%DION(CH1,CH2)
      ENDDO      
      ENDDO      

      DEALLOCATE(POTCAR_DATA, SPLINE_VALUE, SPLINE_DATA, SPLINE_COEF)

      DEALLOCATE(RHOPS00, RHOAE00)
      DEALLOCATE(RHO, V, POT, POTAEC, POTAE_EFF, DPOTAE_EFF, POTPS_EFF)
      DEALLOCATE(POT_TEST,POTAE_TEST, POTPSC_TEST, POTPSC_CHECK, POTPS_G, CORPS_G, RHOPS_G)
      DEALLOCATE(POTPS_TEST, CRHODE, RHOLM)
!      DEALLOCATE(V1, V2, CRHODE, POT-AE, DLM)

      DEALLOCATE(VTMP, DLLMM, DHXC)
      DEALLOCATE(TMP, DTMP, PARWKINAE, PARWKINPS, QPAW)
      DEALLOCATE(DIJ, DION)

!      STOP

      RETURN

       END SUBROUTINE GENERATE_POTCARDATA

      SUBROUTINE GENERATE_POTCARLIB (INFO, FROM_PP, CHANNELS)
!         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: FROM_PP
         TYPE(INFO_STRUCT) :: INFO
!         TYPE(in_struct) IO

!      CHARACTER*40 SZNAMP         ! header
      CHARACTER*80 :: CSEL
      CHARACTER(LEN=4) :: PARAM_CHARACTER
      REAL(q)  ::   POTCAR_PARAM
      REAL(q)  ::   POTCAR_G_PARAM, POTCAR_R_PARAM
!      REAL(q)  ::   WAE(Grid%n,PAW%nbase), WPS(Grid%n,PAW%nbase)
      INTEGER  ::   LSTATE, NWRITE, XC_TYPE, L1, L2, NMAX, CHANNELS
      REAL(q), ALLOCATABLE :: POTCAR_DATA(:), VASP_PSNL_CHECK(:), VASP_PSRNL_CHECK(:), VASP_PROJ_R(:)
      REAL(q), ALLOCATABLE :: SPLINE_COEF(:,:), SPLINE_VALUE(:), SPLINE_DATA(:)
      LOGICAL  ::   PARAM_LOG

      ALLOCATE(POTCAR_DATA(NPSPTS))
!     ALLOCATE(SPLINE_VALUE(FROM_PP%R%NMAX))

      OPEN(UNIT=30,FILE='POTCAR',STATUS='OLD',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=30,FILE='POTCAR',STATUS='OLD')
      ENDIF

!      OPEN(UNIT=88,FILE='POTCAR_NEW',STATUS='UNKNOWN',IOSTAT=IERR)
      OPEN(UNIT=89,FILE='POTCAR.lib',STATUS='UNKNOWN',FORM='FORMATTED', IOSTAT=IERR)
      IF (IERR/=0) THEN
!         OPEN(UNIT=88,FILE='POTCAR_NEW',STATUS='UNKNOWN')
         OPEN(UNIT=89,FILE='POTCAR.lib',STATUS='UNKNOWN')
      ENDIF

      LPAW = .FALSE.
      INFO%LOVERL=.FALSE.
      INFO%LCORE =.FALSE.
      
!      READ(30,'(A40)',END=100,ERR=100) SZNAMP
!      WRITE(88,'(A40)') SZNAMP
!      READ(30,*) 
!      WRITE(88,*) FROM_PP%ZVALF
!!!!!   HEAD, ELEMENT & PSCTR PARAMETER   !!!!!
!      DO I=1,26
!        READ(30,'(A80)') CSEL
!!        WRITE(88,'(A80)') CSEL
!        WRITE(89,'(A80)') CSEL
!      ENDDO
!      
!      READ(30, '(4X, I1,A8)') LSTATE, CSEL
!      WRITE(89, '(4X, I1,A8)') LSTATE, CSEL
!      DO I=1, LSTATE+1
!        READ(30,'(A80)') CSEL
!!        WRITE(88,'(A80)') CSEL
!        WRITE(89,'(A80)') CSEL
!      ENDDO

!      DO I=1, 8
!        READ(30,'(A80)') CSEL
!!        WRITE(88,'(A80)') CSEL
!        WRITE(89,'(A80)') CSEL
!      ENDDO
200   READ(30,'(A80)') CSEL
      IF (CSEL(1:32) .EQ. 'END of PSCTR-controll parameters') THEN
         GOTO 305
      ELSE
         WRITE(89,'(A80)') CSEL
         GOTO 200
      ENDIF
305      WRITE(89, '(A80)') CSEL

!      DO I=1,18
        READ(30,'(A80)') CSEL
!!        WRITE(88,'(A80)') CSEL
        WRITE(89,'(A80)') CSEL
!      ENDDO
      
!!!!!   LOCAL PART POTENTIAL   !!!!!
      READ(30,*) POTCAR_PARAM
      READ(30,*) (POTCAR_DATA(I),I=1,NPSPTS)
!      WRITE(88,'(F19.13)') POTCAR_PARAM
      WRITE(89,'(F19.13)') POTCAR_PARAM
!      WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
      WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)

      READ(30,'(1X,A1)') CSEL
      IF (CSEL=='g') THEN
!!!!!          XC_TYPE         !!!!!
!         WRITE(88,*) ' gradient corrections used for XC '
         WRITE(89,*) 'gradient corrections used for XC '
         READ(30,*) XC_TYPE
!         WRITE(88,'(I12)') XC_TYPE
         WRITE(89,*) XC_TYPE

        READ(30,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
!        WRITE(89,'(A80)') CSEL
      ENDIF
!!!!!      PSPCOR             !!!!!
      IF (CSEL=='c') THEN
!         WRITE(88,*) ' core charge-density (partial) '
         WRITE(89,*) 'core charge-density (partial) '
         READ(30,*) (POTCAR_DATA (I),I=1,NPSPTS)
!         WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         READ(30,'(1X, A1)') CSEL
      ELSE
!         NULLIFY(P(NTYP)%PSPCOR)
      ENDIF

!!!!!  KINETIC ENERGY PARTIAL     !!!!!
      IF (CSEL=='k') THEN
!        WRITE(88,'(A80)') CSEL
         WRITE(89,*) 'kinetic energy density (partial)'
         READ(30,*) (POTCAR_DATA (I),I=1,NPSPTS)
!        WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         READ(30,'(1X, A1)') CSEL
!        WRITE(88,'(A80)') CSEL
      ELSE
         NULLIFY(FROM_PP%PSPTAU)
      ENDIF

!!!!!  KINETIC ENERGY DENSITY VALENCE    !!!!!
      IF (CSEL=='K') THEN
         WRITE(89,*) 'Kinetic energy density valence'
         READ(30,*) (POTCAR_DATA (I),I=1,NPSPTS)
!        WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         READ(30,'(1X,A1)') CSEL
      ELSE
         NULLIFY(FROM_PP%PSPTAUVAL)
      ENDIF

!!!!!   ATOMIC PSEUDO CHARGE-DENSITY     !!!!!
!      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      IF(NWRITE >=0) WRITE(89,*) 'atomic pseudo charge-density'
      READ(30,*) (POTCAR_DATA (I),I=1,NPSPTS)
!      WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
      WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)

!!!!!   NON LOCAL PART PROJECTOR   !!!!!
!      WRITE(6,*) 'nbase=',PAW%nbase, PAW%l(PAW%nbase)
      READ(30,*) POTCAR_G_PARAM, PARAM_CHARACTER
!      WRITE(88,500) POTCAR_G_PARAM, PARAM_CHARACTER
      WRITE(89,500) POTCAR_G_PARAM, PARAM_CHARACTER

      OPEN(UNIT=29,FILE='VASP_PSPNL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=29,FILE='VASP_PSPNL',STATUS='UNKNOWN')
      ENDIF

      OPEN(UNIT=31,FILE='VASP_PSRNL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=31,FILE='VASP_PSRNL',STATUS='UNKNOWN')
      ENDIF


      DO L = 0,  PAW%l(PAW%nbase)
        READ(30,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
        WRITE(89,'(A80)') CSEL
        READ(30,*) L1, NL1, POTCAR_R_PARAM
!        WRITE(88,'(2I12, F19.14)') L1, NL1, POTCAR_R_PARAM
        WRITE(89,'(2I12, F19.14)') L1, NL1, POTCAR_R_PARAM
        READ(30,*) (POTCAR_DATA (I),I=1,4)
!        WRITE(88,'(F19.14, 2F24.13)') (POTCAR_DATA (I),I=1,4)
        WRITE(89,'(F19.14, 2F24.13)') (POTCAR_DATA (I),I=1,4)

        DO LI = 1, NL1
!        ALLOCATE(VASP_PROJ_R(NPSNL), VASP_PSNL_CHECK(NPSNL), VASP_PSRNL_CHECK(NSPNL))
        !!!!!   RECIPROCAL SPACE PART   !!!!!
           READ(30,'(A80)') CSEL
!           WRITE(88,'(A80)') CSEL
           WRITE(89,'(A80)') CSEL
!           WRITE(6,*) 'NPSNL=', NPSNL
           READ(30,*) (POTCAR_DATA (I),I=1,NPSNL)
!           WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSNL)
           WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSNL)

           DO I=1, NPSNL
              WRITE(29,*) POTCAR_G_PARAM/NPSNL*(I-1), POTCAR_DATA(I)
           ENDDO


       !!!!!   REAL SPACE PART   !!!!!
           READ(30,'(A80)') CSEL
!           WRITE(88,'(A80)') CSEL
           WRITE(89,'(A80)') CSEL
           READ(30,*) (POTCAR_DATA (I),I=1,NPSNL)
!           WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSNL)
           WRITE(89,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSNL)


           DO I=1, NPSNL
              WRITE(31,*) POTCAR_R_PARAM/NPSNL*(I-1), POTCAR_DATA(I)
!              FROM_PP%R%R(I)=POTCAR_R_PARAM/NPSNL*(I-1)
!             VASP_PSRNL_CHECK(I)=POTCAR_DATA(I)
           ENDDO
!           FROM_PP%R%NMAX = NPSNL

!     !!!!!   FOR CHECK !!!!!
!           DO I=1, NPSNL
!              WRITE(35,*)VASP_PROJ_R(I), POTCAR_DATA(I)
!           ENDDO
!           VASP_PSNL_CHECK = 0.d0
!           CALL FOURPOT_TO_Q( VASP_PROJ_R(NPSNL), VASP_PSRNL_CHECK,   &
!     &                   VASP_PSNL_CHECK, NPSNL, POTCAR_G_PARAM/NPSNL, FROM_PP%R, IU6)
!           DO I=1, NPSNL
!              WRITE(35,*)POTCAR_G_PARAM/NPSNL*(I-1), VASP_PSNL_CHECK(I)
!           ENDDO
!           STOP
!           DEALLOCATE(VASP_PROJ_R, VASP_PSNL_CHECK, VASP_PSRNL_CHECK)

        ENDDO
        WRITE(29,*)
        WRITE(31,*)
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!        PAW PART        !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) NMAX, POTCAR_PARAM
!      WRITE(88,600) FROM_PP%R%NMAX, POTCAR_PARAM
      FROM_PP%R%NMAX = NMAX
      WRITE(89,600) FROM_PP%R%NMAX, POTCAR_PARAM
      READ(30,*) 
!      WRITE(88,*)'(5E20.12)' 
      WRITE(89,'(A9,32X)')'(5E20.12)' 

      DO I =1, 2
!!!!!   augmentation charge (non spherical)   !!!!!
         READ(30,'(A80)') CSEL
!         WRITE(88,'(A80)') CSEL
         WRITE(89,'(A80)') CSEL
         READ(30,*) (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS)
!         WRITE(88,'(5E20.12)') (POTCAR_DATA (J),J=1,NRANGE)
         WRITE(89,'(5E20.12)') (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS)
      ENDDO
!!!!!            grid              !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      FROM_PP%R%R(1:FROM_PP%R%NMAX) = POTCAR_DATA(1:FROM_PP%R%NMAX)
!      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (FROM_PP%R%R(I),I=1,FROM_PP%R%NMAX)

!!!!!            aepotential              !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      DO I = 1, Grid%n
!         POTAE(I) = -POTAE(I)/Grid%r(I)*RYTOEV
!      ENDDO
!      CALL SPLINE(Grid%r*AUTOA, POTAE, SPLINE_COEF(1,1:Grid%n),  &
!     &       SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!      DO I=1, FROM_PP%R%NMAX
!          SPLINE_VALUE(I) = ISPLINE(FROM_PP%R%R(I), Grid%r*AUTOA, POTAE,  &
!     &      SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!                  WRITE(99,'(4F22.13)') FROM_PP%R%R(I), FROM_PP%POTAE(I), &
!     &                 SPLINE_VALUE(I), (FROM_PP%POTAE(I)-SPLINE_VALUE(I))*FROM_PP%R%R(I)
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      WRITE(88,'(5E20.12)') (SPLINE_VALUE (I),I=1,FROM_PP%R%NMAX)

!!!!!        core charge-density          !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      DO I = 1, Grid%n
!!         COREDEN(I) = COREDEN(I)/(SCALE*AUTOA)/(SCALE*AUTOA) 
!         COREDEN(I) = COREDEN(I)/(SCALE*AUTOA)
!      ENDDO
!      CALL SPLINE(Grid%r*AUTOA, COREDEN, SPLINE_COEF(1,1:Grid%n),  &
!     &       SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!      DO I=1, FROM_PP%R%NMAX
!          SPLINE_VALUE(I) = ISPLINE(FROM_PP%R%R(I), Grid%r*AUTOA, COREDEN,  &
!     &      SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!!                  WRITE(99,*) FROM_PP%R%R(I), FROM_PP%RHOAE(I), -SPLINE_VALUE(I)
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      WRITE(88,'(5E20.12)') (-SPLINE_VALUE (I),I=1,FROM_PP%R%NMAX)

!!!!!       kinetic energy-density        !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!    mkentic energy-density pseudized    !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!      local pseudo-potential core      !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!     pseudo-potential valence only      !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!     core-charge density (pseudized)     !!!!!
      READ(30,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
      WRITE(89,'(A80)') CSEL
      READ(30,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(89,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DO I = 1, Grid%n
!         COREDEN(I) = PAW%tcore(I)*AUTOA
!      ENDDO
!      CALL SPLINE(Grid%r*AUTOA, COREDEN, SPLINE_COEF(1,1:Grid%n),  &
!     &       SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!      DO I=1, FROM_PP%R%NMAX
!          SPLINE_VALUE(I) = ISPLINE(FROM_PP%R%R(I), Grid%r*AUTOA, COREDEN,  &
!     &      SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!!                  WRITE(99,*) FROM_PP%R%R(I), FROM_PP%RHOPS(I), SPLINE_VALUE(I)
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!      WRITE(88,'(5E20.12)') (SPLINE_VALUE (I),I=1,FROM_PP%R%NMAX)

!!!!!     wave-function     !!!!!
!      WRITE(6,*) 'nbase=', PAW%nbase
      DO I = 1, PAW%nbase
         !!!!!   pseudo  wave-function   !!!!!
           READ(30,'(A80)') CSEL
!           WRITE(88,'(A80)') CSEL
           WRITE(89,'(A80)') CSEL
           READ(30,*) (POTCAR_DATA (K),K=1,FROM_PP%R%NMAX)
           WRITE(89,'(5E20.12)') (POTCAR_DATA (K),K=1,FROM_PP%R%NMAX)
!           DO K = 1, Grid%n
!              WPS(K) = PAW%tphi(K,I)/SQRT(AUTOA)
!           ENDDO
!           CALL SPLINE(Grid%r*AUTOA, WPS, SPLINE_COEF(1,1:Grid%n),  &
!     &             SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!           DO K=1, FROM_PP%R%NMAX
!                  SPLINE_VALUE(K) = ISPLINE(FROM_PP%R%R(K), Grid%r*AUTOA, WPS,  &
!     &           SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!!                  WRITE(99,*) FROM_PP%R%R(K), FROM_PP%WPS(K,I), SPLINE_VALUE(K)
!           ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!           WRITE(88,'(5E20.12)') (SPLINE_VALUE (K),K=1,FROM_PP%R%NMAX)

         !!!!!   ae  wave-function   !!!!!
           READ(30,'(A80)') CSEL
!           WRITE(88,'(A80)') CSEL
           WRITE(89,'(A80)') CSEL
           READ(30,*) (POTCAR_DATA (K),K=1,FROM_PP%R%NMAX)
           WRITE(89,'(5E20.12)') (POTCAR_DATA (K),K=1,FROM_PP%R%NMAX)
!           DO K = 1, Grid%n
!              WAE(K) = PAW%phi(K,I)/SQRT(AUTOA)
!           ENDDO
!           CALL SPLINE(Grid%r*AUTOA, WAE, SPLINE_COEF(1,1:Grid%n),  &
!     &             SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!           DO K=1, FROM_PP%R%NMAX
!                  SPLINE_VALUE(K) = ISPLINE(FROM_PP%R%R(K), Grid%r*AUTOA, WAE,  &
!     &           SPLINE_COEF(1,1:Grid%n),  SPLINE_COEF(2,1:Grid%n), SPLINE_COEF(3,1:Grid%n), Grid%n)
!!                  WRITE(99,*) FROM_PP%R%R(K), FROM_PP%WAE(K,I), SPLINE_VALUE(K)
!           ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!           WRITE(88,'(5E20.12)') (SPLINE_VALUE (K),K=1,FROM_PP%R%NMAX)
      ENDDO
!      WRITE(88,*) 'End of Dataset'
      WRITE(89,*) 'End of Dataset'
100   CONTINUE

500      FORMAT(F19.13, 6X, A1)
600      FORMAT(I12, F19.15)
      DEALLOCATE(POTCAR_DATA)
!      DEALLOCATE(SPLINE_COEF, SPLINE_VALUE, SPLINE_DATA)
      CLOSE(30)
!      CLOSE(88)
      CLOSE(89)
      RETURN

       END SUBROUTINE GENERATE_POTCARLIB


      SUBROUTINE WRITE_POTCAR (INFO, FROM_PP, CHANNELS)
!         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: FROM_PP
         TYPE(INFO_STRUCT) :: INFO

!      CHARACTER*40 SZNAMP         ! header
      CHARACTER*80 :: CSEL
      CHARACTER(LEN=4) :: PARAM_CHARACTER
      REAL(q)  ::   POTCAR_PARAM, COREDEN(Grid%n), POTAE(Grid%n)
      REAL(q)  ::   POTCAR_G_PARAM, POTCAR_R_PARAM
!      REAL(q)  ::   WAE(Grid%n,PAW%nbase), WPS(Grid%n,PAW%nbase)
!      REAL(q)  ::   WAE(Grid%n), WPS(Grid%n)
      INTEGER  ::   LSTATE, NWRITE, XC_TYPE, L1, L2, NMAX, CHANNELS
      REAL(q), ALLOCATABLE :: POTCAR_DATA(:)
!      REAL(q), ALLOCATABLE :: SPLINE_COEF(:,:), SPLINE_VALUE(:), SPLINE_DATA(:)
      LOGICAL  ::   PARAM_LOG

      ALLOCATE(POTCAR_DATA(NPSPTS))

      OPEN(UNIT=89,FILE='POTCAR.lib',STATUS='OLD',FORM='FORMATTED', IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=89,FILE='POTCAR.lib',STATUS='OLD', FORM='FORMATTED')
      ENDIF

      OPEN(UNIT=88,FILE='POTCAR_NEW',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=88,FILE='POTCAR_NEW',STATUS='UNKNOWN')
      ENDIF

      LPAW = .FALSE.
      INFO%LOVERL=.FALSE.
      INFO%LCORE =.FALSE.
      
!      READ(89,'(A40)',END=100,ERR=100) SZNAMP
!      WRITE(88,'(A40)') SZNAMP
!      READ(389*) 
!      WRITE(88,*) FROM_PP%ZVALF
!!!!!   HEAD, ELEMENT & PSCTR PARAMETER   !!!!!
!      DO I=1,26
!        READ(89,'(A80)') CSEL
!!        WRITE(88,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
!      ENDDO
!      
!      READ(89, '(4X, I1,A8)') LSTATE, CSEL
!      WRITE(88, '(4X, I1,A8)') LSTATE, CSEL
!      DO I=1, LSTATE+1
!        READ(89,'(A80)') CSEL
!!        WRITE(88,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
!      ENDDO

!      DO I=1,8
!        READ(89,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
!      ENDDO
200   READ(89,'(A80)') CSEL
      IF (CSEL(1:32) .EQ. 'END of PSCTR-controll parameters') THEN
         GOTO 304
      ELSE
         WRITE(88,'(A80)') CSEL
         GOTO 200
      ENDIF
304      WRITE(88, '(A80)') CSEL

!      DO I=1,18
        READ(89,'(A80)') CSEL
!        WRITE(88,'(A80)') CSEL
        WRITE(88,'(A80)') CSEL
!      ENDDO
      
!!!!!   LOCAL PART POTENTIAL   !!!!!
      READ(89,*) POTCAR_PARAM
      READ(89,*) (POTCAR_DATA(I),I=1,NPSPTS)
      WRITE(88,'(F19.13)') POTCAR_PARAM
      WRITE(88,'(5E16.8)') (FROM_PP%PSP(I,2),I=1,SIZE(FROM_PP%PSP,1))

      READ(89,'(1X,A1)') CSEL
      IF (CSEL=='g') THEN
!!!!!          XC_TYPE         !!!!!
         WRITE(88,*) 'gradient corrections used for XC '
         READ(89,*) XC_TYPE
         WRITE(88,'(I12)') XC_TYPE

         READ(89,'(1X,A1)') CSEL
!        WRITE(88,'(A80)') CSEL
      ENDIF
!!!!!      PSPCOR             !!!!!
      IF (CSEL=='c') THEN
!         WRITE(88,*) ' core charge-density (partial) '
         WRITE(88,*) 'core charge-density (partial) '
         READ(89,*) (POTCAR_DATA (I),I=1,NPSPTS)
         WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         READ(89,'(1X, A1)') CSEL
!         WRITE(88,'(A80)') CSEL
      ELSE
!         NULLIFY(P(NTYP)%PSPCOR)
      ENDIF

!!!!!  KINETIC ENERGY PARTIAL     !!!!!
      IF (CSEL=='k') THEN
         WRITE(88,*) 'kinetic energy density (partial)'
         READ(89,*) (POTCAR_DATA (I),I=1,NPSPTS)
         WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         READ(89,'(1X,A1)') CSEL
      ELSE
         NULLIFY(FROM_PP%PSPTAU)
      ENDIF

!!!!!  KINETIC ENERGY DENSITY VALENCE    !!!!!
      IF (CSEL=='K') THEN
         WRITE(88,*) 'Kinetic energy density valence'
         READ(89,*) (POTCAR_DATA (I),I=1,NPSPTS)
         WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSPTS)
         READ(89,'(1X,A1)') CSEL
      ELSE
         NULLIFY(FROM_PP%PSPTAUVAL)
      ENDIF

!!!!!   ATOMIC PSEUDO CHARGE-DENSITY     !!!!!
!      READ(89,'(A80)') CSEL
      IF(NWRITE >=0) WRITE(88,*) 'atomic pseudo charge-density'
      READ(89,*) (POTCAR_DATA (I),I=1,NPSPTS)
      WRITE(88,'(5E16.8)') (FROM_PP%PSPRHO(I),I=1,SIZE(FROM_PP%PSP,1))

!!!!!   NON LOCAL PART PROJECTOR   !!!!!
!      WRITE(6,*) 'nbase=',PAW%nbase, PAW%l(PAW%nbase)
      READ(89,*) POTCAR_G_PARAM, PARAM_CHARACTER
      WRITE(88,500) POTCAR_G_PARAM, PARAM_CHARACTER

      OPEN(UNIT=29,FILE='VASP_PSNL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=29,FILE='VASP_PSNL',STATUS='UNKNOWN')
      ENDIF

      OPEN(UNIT=31,FILE='VASP_PSRNL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=31,FILE='VASP_PSRNL',STATUS='UNKNOWN')
      ENDIF


      DO L = 0,  PAW%l(PAW%nbase)
        READ(89,'(A80)') CSEL
        WRITE(88,'(A80)') CSEL
        READ(89,'(2I12, F19.14)') L1, NL1, POTCAR_R_PARAM
        WRITE(88,'(2I12, F19.14)') L1, NL1, POTCAR_R_PARAM
        READ(89,*) (POTCAR_DATA (I),I=1,4)
        WRITE(88,'(F19.14, 2F24.13)') ((FROM_PP%DION(I,J),I=2*L+1, 2*(L+1)), J=2*L+1,2*(L+1))
!        WRITE(88,'(F19.14, 2F24.13)') (DION(I),I=1,4)

        DO LI = 1, NL1
        !!!!!   RECIPROCAL SPACE PART   !!!!!
           READ(89,'(A80)') CSEL
           WRITE(88,'(A80)') CSEL
!           WRITE(6,*) 'NPSNL=', NPSNL
           READ(89,*) (POTCAR_DATA (I),I=1,NPSNL)
           WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSNL)

!           DO I=1, NPSNL
!              WRITE(29,*) POTCAR_G_PARAM/NPSNL*(I-1), POTCAR_DATA(I)
!           ENDDO


       !!!!!   REAL SPACE PART   !!!!!
           READ(89,'(A80)') CSEL
           WRITE(88,'(A80)') CSEL
           READ(89,*) (POTCAR_DATA (I),I=1,NPSNL)
           WRITE(88,'(5E16.8)') (POTCAR_DATA (I),I=1,NPSNL)


!           DO I=1, NPSNL
!              WRITE(31,*) POTCAR_R_PARAM/NPSNL*(I-1), POTCAR_DATA(I)
!           ENDDO


        ENDDO
!        WRITE(29,*)
!        WRITE(31,*)
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!        PAW PART        !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) NMAX, POTCAR_PARAM
      WRITE(88,600) FROM_PP%R%NMAX, POTCAR_PARAM
      READ(89,*) 
      WRITE(88,'(A9,32X)')'(5E20.12)' 

       READ(89,'(A80)') CSEL
       WRITE(88,'(A80)') CSEL
       WRITE(6,*) 'CSEL=', CSEL
       READ(89,*) (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS)
       WRITE(6,*) 'CHANNELS=', CHANNELS
       WRITE(6,*) 'TEST=', (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS) 
       WRITE(88,'(5E20.12)') ((FROM_PP%QPAW (I,J,0),J=1,CHANNELS),I=1,CHANNELS)
!!!!!   uccopancies in atom   !!!!!
!      NRANGE = (PAW%l(PAW%nbase)+1)**4
!!!!!   augmentation charge (non spherical)   !!!!!
       READ(89,'(A80)') CSEL
       WRITE(88,'(A80)') CSEL
       WRITE(6,*) 'CSEL=', CSEL
       READ(89,*) (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS)
       WRITE(6,*) 'CHANNELS=', CHANNELS
       WRITE(6,*) 'TEST=', (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS) 
       WRITE(88,'(5E20.12)') (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS)
!       WRITE(88,'(5E20.12)') (FROM_PP%QPAW (J),J=1,CHANNELS*CHANNELS)
!!!!!   uccopancies in atom   !!!!!
!      READ(89,'(A80)') CSEL
!      WRITE(88,'(A80)') CSEL
!      READ(89,*) (POTCAR_DATA (J),J=1,CHANNELS*CHANNELS)
!      WRITE(88,'(5E20.12)') (POTCAR_DATA (J),J=1, CHANNELS*CHANNELS)
!!!!!            grid              !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(88,'(5E20.12)') (FROM_PP%R%R(I),I=1,FROM_PP%R%NMAX)

!!!!!            aepotential              !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!      WRITE(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      WRITE(88,'(5E20.12)') (FROM_PP%POTAE(I), I=1,FROM_PP%R%NMAX)

!!!!!        core charge-density          !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      WRITE(88,'(5E20.12)') (FROM_PP%RHOAE(I), I=1,FROM_PP%R%NMAX)

!!!!!       kinetic energy-density        !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!    mkentic energy-density pseudized    !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!      local pseudo-potential core      !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(88,'(5E20.12)') (FROM_PP%POTPSC(I),I=1,FROM_PP%R%NMAX)
!!!!!     pseudo-potential valence only      !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
      WRITE(88,'(5E20.12)') (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!     core-charge density (pseudized)     !!!!!
      READ(89,'(A80)') CSEL
      WRITE(88,'(A80)') CSEL
      READ(89,*) (POTCAR_DATA (I),I=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(88,'(5E20.12)') (FROM_PP%RHOPS(I),I=1,FROM_PP%R%NMAX)

!!!!!     wave-function     !!!!!
!      WRITE(6,*) 'nbase=', PAW%nbase
      DO I = 1, PAW%nbase
         !!!!!   pseudo  wave-function   !!!!!
           READ(89,'(A80)') CSEL
           WRITE(88,'(A80)') CSEL
           READ(89,*) (POTCAR_DATA (K),K=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
           WRITE(88,'(5E20.12)') (FROM_PP%WPS(K,I),K=1,FROM_PP%R%NMAX)

         !!!!!   ae  wave-function   !!!!!
           READ(89,'(A80)') CSEL
           WRITE(88,'(A80)') CSEL
           READ(89,*) (POTCAR_DATA (K),K=1,FROM_PP%R%NMAX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
           WRITE(88,'(5E20.12)') (FROM_PP%WAE (K,I),K=1,FROM_PP%R%NMAX)
      ENDDO
      WRITE(88,*) 'End of Dataset'
100   CONTINUE

500      FORMAT(F19.13, 6X, A1)
600      FORMAT(I12, F19.15)
      DEALLOCATE(POTCAR_DATA)
      CLOSE(88)
      CLOSE(89)

       END SUBROUTINE WRITE_POTCAR

      subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
      implicit none
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer i, j, gap
      double precision h

         gap = n-1
         b(1:n) = 0.d0
         c(1:n) = 0.d0
         d(1:n) = 0.d0
! check input
         if ( n < 2 ) return
         if ( n < 3 ) then
             b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
             c(1) = 0.
             d(1) = 0.
             b(2) = b(1)
             c(2) = 0.
             d(2) = 0.
             return
         end if
!
! step 1: preparation
!
        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, gap
           d(i) = x(i+1) - x(i)
           b(i) = 2.0*(d(i-1) + d(i))
           c(i+1) = (y(i+1) - y(i))/d(i)
           c(i) = c(i+1) - c(i)
        end do
!
! step 2: end conditions 
!
        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = 0.0
        c(n) = 0.0
        if (n /= 3) then
           c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
           c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
           c(1) = c(1)*d(1)**2/(x(4)-x(1))
           c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
        end if
!
! step 3: forward elimination 
!
        do i = 2, n
           h = d(i-1)/b(i-1)
           b(i) = b(i) - h*d(i-1)
           c(i) = c(i) - h*c(i-1)
        end do
!
! step 4: back substitution
!
        c(n) = c(n)/b(n)
        do j = 1, gap
           i = n-j
           c(i) = (c(i) - d(i)*c(i+1))/b(i)
        end do
!
! step 5: compute spline coefficients
!
        b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
        do i = 1, gap
           b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
           d(i) = (c(i+1) - c(i))/d(i)
           c(i) = 3.*c(i)
        end do
        c(n) = 3.0*c(n)
        d(n) = d(n-1)
      end subroutine spline

      function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
        implicit none
        double precision ispline
        integer n
        double precision  u, x(n), y(n), b(n), c(n), d(n)
        integer i, j, k
        double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
        if (u <= x(1)) then
           ispline = y(1)
           return
        end if
        if (u >= x(n)) then
           ispline = y(n)
           return
        end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
        i = 1
        j = n+1
        do while (j > i+1)
           k = (i+j)/2
           if (u < x(k)) then
              j=k
           else
              i=k
           end if
        end do
!*
!  evaluate spline interpolation
!*
        dx = u - x(i)
        ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      end function ispline


      END MODULE VASP_POTCAR
