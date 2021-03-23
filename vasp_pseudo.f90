      MODULE VASP_POTCAR
         USE pseudo         ! VASP_Pseudo
         USE base           ! VASP_base
         USE ini            ! VASP_ini/PREC
         USE setexm
!         USE core_rel
         USE PSEUDO_struct  ! VASP_Pseudo_struct
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
      SUBROUTINE vasp_pseudo(ifinput, ifen, Grid0, Orbit, Pot_AE, Pot_FC, success)!(ifinput,ifen,Grid)
         IMPLICIT COMPLEX(q)   (C)
         IMPLICIT REAL(q)  (A-B,D-H,O-Z)

         TYPE(GridInfo), INTENT(IN) :: Grid0
!         REAL(8), INTENT(IN) :: coreden(:)
         Type(OrbitInfo), INTENT(IN) :: Orbit
         TYPE(PotentialInfo), INTENT(IN) :: Pot_AE
         TYPE(PotentialInfo), INTENT(IN) :: Pot_FC
         REAL(q), ALLOCATABLE :: PotHr(:), PotXCr(:), PotAEr(:), PotAECr(:)
         REAL(q), ALLOCATABLE :: Pot_eff(:), Pot_teff(:)
         REAL(q), ALLOCATABLE :: PotAE(:), PotAE00(:), PotPS(:), PotPSC(:), PotPSCr(:)
         REAL(q), ALLOCATABLE :: POTAE_TEST(:)
         REAL(q), ALLOCATABLE :: pdensity(:), d(:) ,cpdensity(:)
         INTEGER :: nbase, N, irc, irc_core, irc_shap
         REAL(q) :: Q_00 , Q_00c,tq
         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: PP
         TYPE(INFO_STRUCT) :: INFO
!         TYPE(PseudoInfo) :: PAW

         INTEGER :: ifinput,ifen,Z
         INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS, LMAX, LMMAX, LI
         INTEGER :: LMAX_TABLE
         REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
         REAL(q) :: DHARTREE, QTEST,SCALE, DOUBLEAE, EXCG
         REAL(q), ALLOCATABLE :: RHO(:,:,:), POT(:,:,:), V(:,:,:), RHOAE00(:), RHOV(:)
         REAL(q), ALLOCATABLE :: POTAEC(:)
         REAL(q), ALLOCATABLE :: CRHODE(:,:)
         REAL(q), ALLOCATABLE :: RHOLM(:), DLM(:)
         CHARACTER(LEN=2) :: TYPE(1)
         LOGICAL ::   LPAW,  success, UNSCREN

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
      OPEN(UNIT=7,FILE='VASP_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
      OPEN(UNIT=8,FILE='VASP_WPS',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=7,FILE='VASP_WAE',STATUS='OLD')
         OPEN(UNIT=8,FILE='VASP_WPS',STATUS='OLD')
      ENDIF
         DO i=1, CHANNELS
          DO j=1, PP%R%NMAX
          WRITE(IU6,*) PP%R%R(j), PP%WAE(j,i)
!          if (mod(i,2) == 0)WRITE(IU6,*) PP%R%R(j), PP%WAE(j,i)
          WRITE(IU7,*) PP%R%R(j), PP%WPS(j,i)
          ENDDO
          WRITE(IU6,*)
          WRITE(IU7,*)
         ENDDO
!
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

      CALL SET_SIMP(PP%R)

!      CALL SIMPI(PP%R,PP%RHOAE, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST, QTEST*SCALE
!      CALL SIMPI(PP%R,PP%RHOPS, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST, QTEST*SCALE

      LMMAX = (PP%LMAX+1)**2
      ALLOCATE(RHO(PP%R%NMAX, LMMAX,1), V(PP%R%NMAX, LMMAX,1), RHOAE00(PP%R%NMAX))
      ALLOCATE(POT(PP%R%NMAX, LMMAX,1), POTAEC(PP%R%NMAX))
      ALLOCATE(POTAE_TEST(PP%R%NMAX))
!      ALLOCATE(V1(PP%R%NMAX, LMMAX,1), V2(PP%R%NMAX, LMMAX,1))
      ALLOCATE(CRHODE(LDIM,LDIM))
!      ALLOCATE(CRHODE(LDIM,LDIM), POT-AE(PP%R%NMAX))
      ALLOCATE(RHOLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
      ALLOCATE(DLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
      RHO = 0
!
!
!!!!!!!!!!!!!    RELATION BETWEEN PP%WAE and RHO_V:   \int PP%WAE**2 dR = \int RHO_V*SCALE dR   !!!!!!
!TEST FOR \int PP%WAE**2 dR
!      DO io =1,CHANNELS
!         LI = PP%LPS(io)
!         WRITE(6,*) 'OCC1=', PP%QATO(io,io)
!         RHOAE00(:) = RHOAE00(:)+PP%WAE(:,io)**2*PP%QATO(io,io)*(2*LI+1)
!      ENDDO
!      CALL SIMPI(PP%R,RHOAE00, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST
!
      LMAX_TABLE=6; CALL YLM3ST_(LMAX_TABLE)
      CALL SET_CRHODE_ATOM(CRHODE,PP)
      CALL TRANS_RHOLM(CRHODE, RHOLM, PP)
      CALL RAD_CHARGE(RHO(:,:,1), PP%R, RHOLM(:),PP%LMAX, PP%LPS, PP%WAE)

      OPEN(UNIT=25,FILE='VASP_VAL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=25,FILE='VASP_VAL',STATUS='OLD')
      ENDIF
!TEST FOR \int PP%RHO_V*SCALE dR
!         DO j=1, PP%R%NMAX
!             WRITE(25,*) PP%R%R(j), RHOAE00(j), RHO(j,1,1)*SCALE
!          ENDDO

!      RHOAE00(:) = RHO(:,1,1)
!      CALL SIMPI(PP%R,RHOAE00, QTEST)
!      WRITE(6,*) 'QTEST=', QTEST*SCALE
!      WRITE(6,*) 'QTEST=', QTEST
!      STOP
!!!!!!!!! RHO(:) = RHOAE(:) + RHO_V(:) !!!!!!!!!
!      RHO(:,1,1)=RHOAE00(:)+PP%RHOAE(:)

!      Z=INT(PP%ZVAL_ORIG+PP%ZCORE)
!      CALL POT(RHO, Z, PP%R, V)
!      RHO(:,1,1)=RHOAE00(:)


!!!!!!!!! POTAE = V_H[n_v]+V_XC[n_v+n_c] != POTAEC  !!!!!!!!!     
      POT = 0
      POTAEC = 0
!      RHO = 0
      CALL PUSH_XC_TYPE(PP%LEXCH, 1.0_q, 1.0_q, 1.0_q, 1.0_q, 0.0_q)
      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
       RHO, PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
      POTAE_TEST(:) = POT(:,1,1)
!!!!!!!!!! VASP_POT is THE SAME AS RAD_POT !!!!!!!!!!!!!!!!!!!!!!!!!!!
!      POTAEC = 0
!      POT = 0
!      CALL VASP_POT(RHO,0, PP%R, POT, PP%RHOAE, POTAEC) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!! POTAEC = V_H[n_Zc] = V_H[n_c]+ Z/r !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! POTAEC calected undirected !!!!!!!!!!!!!!!!!!!!
!      RHO(:,1,1) = RHO(:,1,1)+PP%RHOAE(:) 
!      POT = 0
!      CALL VASP_POT(RHO,14, PP%R, POT) !
!!!!!! V_H[n_Zc] = V_H[n_c+n_v]+Z/r+V_XC[n_c+n_v] - V_H[n_v]-V_XC[n_c+n_v]
!      POTAEC(:) = POT(:,1,1)-POTAE_TEST(:)

!!!!!!!!! POTAEC calculated Directed !!!!!!!!!
      CALL RAD_POT_HAR(0, PP%R, POTAEC, PP%RHOAE, DHARTEE)
      DO j=1, PP%R%NMAX
         POTAEC(j) = POTAEC(j)-FELECT*SCALE/PP%R%R(j)*14.00
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!! POT_EFF = V_H[n_v+n_Zc] +V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!!!!!!!!!!!
!      POTAE_TEST = 0
!      POT = 0
!      POTAEC(:) = POTAEC(:)/SCALE
!      CALL RAD_POT(PP%R, 1, 1, 1, .FALSE., &
!       RHO, PP%RHOAE, POTAEC, POT, DOUBLEAE, EXCG)
!      POTAE_TEST(:) =  POT(:,1,1)

      RHO(:,1,1) = RHO(:,1,1)+PP%RHOAE(:) 
      POT = 0
      CALL VASP_POT(RHO,14, PP%R, POT) !

!      CALL VASP_POT(RHO, 0, PP%R, V2, PP%RHOAE)
!      CALL VASP_POT(RHO, 0, PP%R, V, PP%RHOAE)

      OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='OLD')
      ENDIF


          DO j=1, PP%R%NMAX
!             WRITE(IU15,'(4f20.8)') PP%R%R(j), PP%POTAE(j), POT(j,1,1)-POTAE_TEST(j), POTAEC(j)
             WRITE(IU15,'(4f20.8)') PP%R%R(j), PP%POTAE(j) !,  POT(j,1,1)
!             WRITE(IU15,'(4f20.8)') PP%R%R(j), POT(j,1,1) !,POTAE_TEST(j)
          ENDDO

!!!!!!!!!!!!! POTPS = V_H[tn_v+tn_aug]+V_XC[tn_v+tn_aug+tn_c] !!!!!!!!!!!!!!!!!!!!!!!!
      RHO = 0
      POT=0

      DLM = 0
      CALL RAD_CHARGE(RHO(:,:,1), PP%R, RHOLM(:),PP%LMAX, PP%LPS, PP%WPS)

      CALL RAD_POT_WEIGHT( PP%R, 1, 0, POT)

      CALL RAD_PROJ(  POT(:,:,1)  , PP%R,-1._q, DLM, PP%LMAX, PP%LPS, PP%WPS )
!      CALL RAD_AUG_PROJ( POT(:,:,1), PP%R, DLM, PP%LMAX, PP%LPS, &
!                  0, PP%AUG, PP%QPAW )

      CALL VASP_POT(RHO, 0, PP%R, POT, PP%RHOPS)

      OPEN(UNIT=21,FILE='VASP_POTPS',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=21,FILE='VASP_POTPS',STATUS='OLD')
      ENDIF
          DO j=1, PP%R%NMAX
             WRITE(IU17,'(4f20.8)') PP%R%R(j),  PP%POTPS(j)    !,POT(j,1,1)/SCALE 
!             WRITE(IU17,'(4f20.8)') PP%R%R(j),  POT(j,1,1)/SCALE, PP%POTPS(j)
          ENDDO
      
!!!!!!!!!!!!! POTPSC = POT_tEFF - POTPS !!!!!!!!!!!!!!!!!!!!!!!!
!      CALL RAD_CHARGE(RHO(:,:,1), PP%R, RHOLM(:),PP%LMAX, PP%LPS, PP%WPS)
!      POT=0
!      CALL VASP_POT(RHO, 0, PP%R, POT, PP%RHOPS)

      WRITE(6,*) 'VALUE=', PP%ZVALF
      WRITE(6,*) 'NPSPTS=', NPSPTS
      WRITE(6,*) 'PSGMAX=', PP%PSGMAX
      WRITE(6,*) 'NMAX=', PP%R%NMAX

      OPEN(UNIT=23,FILE='VASP_POTPSC',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=23,FILE='VASP_POTPSC',STATUS='OLD')
      ENDIF

      CALL POTTORHO( 0.0, NPSPTS, PP%PSP(:,2), PP%PSGMAX/NPSPTS, &
                 .TRUE. , PP%R%NMAX, PP%R%R ,  PP%POTPSC )                        

          DO j=1, PP%R%NMAX
!             WRITE(IU19,'(4f20.8)') PP%R%R(j), PP%R%R(j)*PP%POTPSC(j)/RYTOEV/AUTOA
!             WRITE(IU19,'(4f20.8)') PP%R%R(j), PP%POTPSC(j), PP%POTAE(j), PP%POTPS(j)
             WRITE(IU19,'(4f20.8)') PP%R%R(j), PP%POTPSC(j)+PP%POTPS(j)
          ENDDO

!             WRITE(IU19,*)
!             WRITE(IU19,*) PP%PSDMAX
!             WRITE(IU19,*)

!!      CALL FOURPOT_TO_Q( PP%RDEP, POT, PP%PSP(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)
!      PP%PSP_TEST(:,2) = 0.0
!      CALL FOURPOT_TO_Q( PP%R%R(PP%R%NMAX), PP%POTPSC, PP%PSP_TEST(:,2), SIZE(PP%PSP,1), PP%PSGMAX/ SIZE(PP%PSP,1), PP%R, IU6)

!          DO j=1, SIZE(PP%PSP,1)
!             WRITE(IU19,'(4f20.8)') PP%PSP(j,2), -PP%PSP_TEST(j,2)
!          ENDDO
!      STOP

!!        WRITE(6,*) 'N=', Grid0%n
      ALLOCATE(PotHr(Grid0%n), PotXCr(Grid0%n), PotAEr(Grid0%n), PotAECr(Grid0%n))
      ALLOCATE(Pot_eff(Grid0%n), Pot_teff(Grid0%n))
      ALLOCATE(PotAE00(Grid0%n), PotPS(Grid0%n))
      ALLOCATE(PotAE(Grid0%n), PotPSC(Grid0%n))

      OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='OLD')
      ENDIF
          WRITE(6,*) 'AEPOT calcualted'

!!!!!!!!! POTAE = V_H[n_v]+V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!     
          call SetPOT(Grid0, FC%coreden, FC%valeden, PotHr, PotXCr, .FALSE.)

!!!!!!!!! POT_EFF = V_H[n_v+n_Zc]+V_XC[n_v+n_c] !!!!!!!!!!!!!!!!!     
          call SetPOT(Grid0, FC%coreden, FC%valeden, PotHr, PotXCr, .TRUE.)
!          PotAEr(:) = PotHr(:)+PotXCr(:)+Pot_AE%rvn(:)*8.0/28.0
          PotAEr(:) = PotHr(:)+PotXCr(:)+Pot_AE%rvn(:)*8.0/28.0
!          DO j = 1, Grid0%n
!             WRITE(6,*) Grid0%r(j)*AUTOA, PotHr(j), PotXCr(j), Pot_AE%rvn(j)
!          ENDDO

!!!!!!!!! POTAEC = V_H[n_c]+V_Z !!!!!!!!!!!!!!!!!     
          FC%valeden = 0.0
          call SetPOT(Grid0, FC%valeden, FC%coreden, PotHr, PotXCr, .FALSE.)
             PotAECr(:) = PotHr(:)+Pot_AE%rvn(:)

!!!!!!!!! tcore_den = sum_i B_isin(q_i r)/r !!!!!!!!!!!!!!!!!     
         Call setcoretail2(Grid0, FC%coreden)
!!         Call setcoretail2(Grid0, coreden)
      OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='OLD')
      ENDIF
          DO j=1,Grid0%n 
             WRITE(IU10,*) Grid0%r(j)*AUTOA, FC%coreden(j)*AUTOA
          ENDDO
!
      OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='OLD')
      ENDIF
          DO j=1,Grid0%n 
!             WRITE(IU12,*) Grid0%r(j)*AUTOA, PAW%tcore(j)/SCALE/AUTOA
             WRITE(IU12,*) Grid0%r(j)*AUTOA, PAW%tcore(j)*AUTOA
          ENDDO
!            irc_core= FindGridIndex(Grid0, PAW%rc_core)
!            Q_00 = integrator(Grid0, PAW%tcore, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!
         Call SetPAWOptions2(ifinput,ifen,Grid0, Orbit,Pot_FC,success)

      OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
      OPEN(UNIT=17,FILE='ATOM_WPS',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='OLD')
         OPEN(UNIT=17,FILE='ATOM_WPS',STATUS='OLD')
      ENDIF
         DO i=1, PAW%nbase
          DO j=1, Grid0%n
             WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
!           if (mod(i,2) == 0) WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
             WRITE(IU13,*) Grid0%r(j)*AUTOA, PAW%tphi(j,i)/sqrt(AUTOA)
          ENDDO
          WRITE(IU8,*)
          WRITE(IU13,*)
         ENDDO

!      ALLOCATE(RHOV(Grid0%n), d(Grid0%n))
!      RHOV = 0
!      DO io =1, PAW%nbase
!!         WRITE(6,*) 'OCC2=', Orbit%occ(io)
!         if (mod(io,2) .ne. 0) RHOV(:) = RHOV(:)+PAW%phi(:,io)**2*Orbit%occ(io)
!      ENDDO
!      DO j=1, Grid0%n
!          WRITE(26,*) Grid0%r(j)*AUTOA, RHOV(j)/AUTOA**2/SCALE
!!           if (mod(i,2) == 0) WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
!!          WRITE(IU13,*) Grid0%r(j)*AUTOA, PAW%tphi(j,i)/sqrt(AUTOA)
!      ENDDO

!!!!!!!!!!!!!!!!!!!! POT_tVeff FROM POT_EFF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
!!   ---------------- Method 1 ----------------------
!         Call Report_Pseudobasis(Grid0,PAW,ifen)
!         Call Set_PAW_MatrixElements(Grid0,PAW)
!         CALL logderiv(Grid,FCPot,PAW)
!         CALL ftprod(Grid)
!         CALL FindVlocfromVeff(Grid,FCOrbit,PAW)

!!   ---------------- Method 2 ----------------------
         call SetPOT_TEFF(Grid0, PotAEr, Pot_teff)
         
          DO j=1,Grid0%n 

             WRITE(IU16,*) Grid0%r(j)*AUTOA, -PotAEr(j)/Grid0%r(j)!*SCALE*RYTOEV , Pot_teff(j)!*SCALE*RYTOEV   !!  POT_EFF
!             WRITE(IU16,*) Grid0%r(j)*AUTOA, PotAEr(j), Pot_teff(j)*Grid0%r(j)!*SCALE*RYTOEV   !!  POT_EFF
          ENDDO

!!!!!!!!!! PSEUDO Calculations !!!!!!!!!!!!!!!!!!!!         
!         nbase=PAW%nbase
!
!!         N=FindGridIndex(Grid0, PAW%rc_shap)
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
!!         N=FindGridIndex(Grid0, PAW%rc_shape)
!!         integrator(Grid,FC%valeden)
!         DO i =1, Grid0%n
!            pdensity(i) = PAW%hatshape(i)*Q_00
!            pdensity(i) = FC%valeden(i)+FC%coreden(i) 
!         ENDDO

         Q_00 = 0.d0; Q_00_c = 0.d0
         ALLOCATE(pdensity(Grid0%n),cpdensity(Grid0%n), PotPSCr(Grid0%n))
         pdensity = 0.d0
         cpdensity = 0.d0
         PotPSCr = 0.d0

!         WRITE(25, *) PAW%hatden
!            WRITE(6,*) 'unscreened POTPS calcualted'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            irc=max(PAW%irc,PAW%irc_shap,PAW%irc_vloc,PAW%irc_core)
            irc= FindGridIndex(Grid0, PAW%rc_shap)
!            N = Grid0%n
!            Q_00 = integrator(Grid0, PAW%den, 1, N) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, PAW%den, 1, irc) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, PAW%tden, 1, N) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, PAW%tden, 1, irc) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = (integrator(Grid0, PAW%den, 1, irc)-integrator(Grid0,PAW%tden,1, irc)) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, FC%coreden, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, FC%coreden, 1, N) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, PAW%tcore, 1, irc_core) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = integrator(Grid0, PAW%tcore, 1, N) 
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!             Q_00= (integrator(Grid0, FC%coreden, 1, irc)-integrator(Grid0, PAW%tcore,1, irc))
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            Q_00 = (integrator(Grid0, PAW%den, 1, irc)-integrator(Grid0,PAW%tden,1, irc))  &
!             + (integrator(Grid0, FC%coreden, 1, irc)-integrator(Grid0, PAW%tcore,1, irc))
            Q_00 = (integrator(Grid0, PAW%den, 1, irc)-integrator(Grid0,PAW%tden,1, irc)) 
            Q_00c = (integrator(Grid0, FC%coreden, 1, irc)-integrator(Grid0, PAW%tcore,1, irc))

!            Q_00 = 2.4
            
!             Q_00= integrator(Grid0, Q_00*PAW%hatden, 1, irc)
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!             Q_00= integrator(Grid0, (10-Q_00)*PAW%hatden, 1, irc)
!           write(6,*) 'Q_00 for atom ', irc, Q_00
!            pdensity=PAW%tden
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! tn = tn_v + that n = tn_v + \sum Qg
              pdensity=PAW%tden!+(Q_00+Q_00c)*PAW%hatden
!              pdensity=PAW%tden+Q_00*PAW%hatden
!              cpdensity=PAW%tcore+Q_00c*PAW%hatden

!!!!!!!!! POTPS = V_H[tn_v+that_n]+V_XC[tn_v+that_n+tn_c] !!!!!!!!!!!!!!!!!     
            call SetPOT(Grid0, PAW%tcore, pdensity, PotHr, PotXCr, .FALSE.)
!            call SetPOT(Grid0, cpdensity, pdensity, PotHr, PotXCr, .FALSE.)

            DO j=1,Grid0%n 
              PotPS(j) = -(PotHr(j)+PotXCr(j))/Grid0%r(j)
            ENDDO

!!!!!!!!!!!!!!!!!!!  POT_PSP = POT_tVEFF - Pot_PS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!           call SetPOT(Grid0,PAW%tcore, pdensity , PotHr, PotXCr, .FALSE.)
            DO j =1, Grid0%n
              PotPSC(j)= Pot_teff(j)-PotPS(j)*SCALE
!              PotPSC(j) = PotPSC(j)/SCALE
            ENDDO

!!   ---------------- POTAEC IN REAL ----------------------
!            call SetPOT_TEFF(Grid0, PotAECr, PotPSC)

         OPEN(UNIT=22,FILE='ATOM_POTPS',STATUS='UNKNOWN',IOSTAT=IERR)
         OPEN(UNIT=24,FILE='ATOM_POTPSC',STATUS='UNKNOWN',IOSTAT=IERR)
         IF (IERR/=0) THEN
             OPEN(UNIT=22,FILE='ATOM_POTPS',STATUS='OLD')
             OPEN(UNIT=24,FILE='ATOM_POTPSC',STATUS='OLD')
         ENDIF

          DO j=1, Grid0%n
             WRITE(IU18,*) Grid0%r(j)*AUTOA, PotPS(j)*RYTOEV!*SCALE
             WRITE(IU20,*) Grid0%r(j)*AUTOA, -PotPSC(j)*RYTOEV!*SCALE
          ENDDO



         STOP

        CALL Report_Pseudopotential(Grid,PAW)
        CALL SPMatrixElements(Grid,FCPot,FC,PAW)

!        CALL FOURPOT_TO_Q()

!        deallocate(coreden)
      DEALLOCATE(RHO, V, RHOAE00, CRHODE, RHOLM)
      DEALLOCATE(POT, POTAE_TEST, POTAEC)
      DEALLOCATE(PotHr, PotXCr, PotAEr,PotAECr, PotPS, PotAE, PotAE00,  PotPSC)
      DEALLOCATE(pdensity, PotPSCr)
        CLOSE(IU6)
        CLOSE(IU8)
        CLOSE(IU9)
        CLOSE(IU10)
        CLOSE(IU11)
        CLOSE(IU12)
        CLOSE(IU13)
!        CLOSE(IU14)
        CLOSE(IU15)
        CLOSE(IU16)
        CLOSE(IU17)
        CLOSE(IU18)
        CLOSE(IU19)
        CLOSE(IU20)
        RETURN
! 
        END SUBROUTINE vasp_pseudo

        SUBROUTINE SetPOT(Grid2, coreden, valeden, POTHR, POTXCR, LADD)
         USE atomdata
         USE aeatom
         USE atomdata
        TYPE(GridInfo) :: Grid2
        INTEGER :: N
        REAL(q), INTENT(IN) :: coreden(:)
        REAL(q), INTENT(IN) :: valeden(:)
        REAL(q) :: density(Grid2%n)
        REAL(q) :: POTXCR(Grid2%n), POTHR(Grid2%n)
        REAL(q) qc, ecoul, v0, etxc, eex, SCALE
        LOGICAL :: LADD

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
!        WRITE(6,*) 'Grid_N=',Grid2%n
!        call poisson(Grid2, qc, coreden, POT%rvh, ecoul, v0)
        IF (.NOT. LADD) THEN
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

        call exch(Grid2, density, POTXCR, etxc,eex)
!        DO i = 1, Grid2%n 
!            POTXCR(i) = POTXCR(i)*FELECT*SCALE/2.0*AUTOA
!        ENDDO
!        Grid2%n = N



!        DO i = 1, Grid2%n
!            WRITE(6, '(3f15.8)') Grid2%r(i), POTHR(i), POTXCR(i)
!        ENDDO
!        call exch(Grid2, coreden, POT%rvx, etxc,eex)
!        call exch(Grid2, valeden, POT%rvx, etxc,eex)
!        DO i = 1, Grid2%n
!           POT%rvx(i) =  Vrxc_tmp(i)-POT%rvx(i)
!        ENDDO

!         POT%rv=(POT%rvh+POT%rvx)*SCALE
!         POT%rv=(POT%rvh+POT%rvx)*FELECT*SCALE/2.0*AUTOA

!        deallocate(density, Vrxc_tmp)
!        deallocate(POTHR, POTXCR)
        END SUBROUTINE SetPOT


        SUBROUTINE SetPOT_TEFF(Grid2, POTAE, POTPS)
         TYPE(GridInfo) :: Grid2
         REAL(q) :: PotPS(Grid2%n), PotAE(Grid2%n), PotAEr(Grid2%n)
         REAL(q) :: ql(2), qr, rc, xx(2), bb(2), alpha, beta
         REAL(q) :: g, gp, gpp, gg
         REAL(q) :: jbes, jbesp, jbespp, amat(2,2), al(2)
         INTEGER :: irc, i

      irc=PAW%irc_vloc
      rc=PAW%rc_vloc

      PotAEr = 0.d0
      POTPS = 0.d0
      PotAEr= POTAE
      alpha=1-rc*Gfirstderiv(Grid2,irc,PotAEr)/PotAEr(irc)
      beta=1.0d0
      call solvbes(xx,alpha,beta,0,2)
      ql(1:2)=xx(1:2)/rc

      DO i=1,2
        qr=ql(i)*rc
        call jbessel(jbes,jbesp,jbespp, 0, 2, qr)
        jbespp=2.d0*ql(i)*jbesp+jbespp*ql(i)*ql(i)*rc
        jbesp=jbes+jbesp*ql(i)*rc
        jbes=jbes*rc
        amat(1,i)=jbes
        amat(2,i)=jbespp
      ENDDO

      bb(1)=PotAEr(irc)
      bb(2)=Gsecondderiv(Grid2, irc,PotAEr)

      det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
      al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
      al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det

      PotPS=PotAEr/Grid%r

        do i=1,irc-1
           qr=ql(1)*Grid2%r(i)
           call jbessel(g,gp,gpp,0,2,qr)
           PotPS(i)=al(1)*g/Grid2%r(i)
           qr=ql(2)*Grid2%r(i)
           call jbessel(g,gp,gpp,0,2,qr)
           PotPS(i)=PotPS(i)+al(2)*g/Grid2%r(i)
        enddo

        END SUBROUTINE SetPOT_TEFF

!*********************** SUBROUTINE POT *********************************      
! Setup the onsite potential V (on a radial grid)
!
! Input
! RHO(:,:,:) charge density, in (charge,magnetization) format 
! Z          nuclear charge
! R          grid
! RHOC(:)    partial core density (optional)
! POTC(:)    frozen core potential (optional)
!
! Output
! V(:,:,:)   potential, in (charge,magnetization) format
!************************************************************************
      SUBROUTINE VASP_POT(RHO,Z,R,V,RHOC,VC,EXC,ADD_GGA)
      
      USE ini
      USE radial
      USE setexm
      USE gridmod
      
      IMPLICIT NONE
      
      TYPE(rgrid) R
      TYPE(GridInfo) :: Grid0
      
      INTEGER, INTENT(IN) :: Z
      REAL(q), DIMENSION(:,:,:), INTENT(IN) :: RHO
      REAL(q), DIMENSION(:,:,:), INTENT(OUT) :: V
      REAL(q), DIMENSION(:), OPTIONAL, INTENT(IN) :: RHOC
      REAL(q), DIMENSION(:), OPTIONAL, INTENT(IN) :: VC
      REAL(q), OPTIONAL, INTENT(OUT):: EXC
      LOGICAL, OPTIONAL, INTENT(IN) :: ADD_GGA
! local variables
      INTEGER ISPIN,LMAX,NMAX
      INTEGER K,L,M,LM,ISP
      REAL(q) SCALE,SUM,QINT,QINF
      REAL(q), DIMENSION(:,:,:), ALLOCATABLE :: WORK1
      REAL(q), DIMENSION(:,:),ALLOCATABLE :: WORK2

      LOGICAL,PARAMETER :: TREL=.TRUE. ! use relativistic corrections to exchange
      LOGICAL,PARAMETER :: TLDA=.TRUE. ! calculate LDA contribution seperately
      ! TLDA=.FALSE. works only for Perdew Burke Ernzerhof
      ! in this case non spherical contributions are missing
      REAL(q) :: EXCG,DHARTREE,DEXC,DVXC,DEXC_GGA,DVXC_GGA,DOUBLEC
      LOGICAL :: ADDGGA

      ADDGGA=ISGGA()
      IF (PRESENT(ADD_GGA)) ADDGGA=ADD_GGA
      
! get dimensions and perform some checks
      LMAX=INT(SQRT(REAL(SIZE(RHO,2)))-1)      
      IF (((LMAX+1)*(LMAX+1))/=SIZE(RHO,2)) THEN
         WRITE(*,*) 'POT: LMAX and the 2nd dimension of RHO do not match:',((LMAX+1)*(LMAX+1)),SIZE(RHO,2)
         STOP
      ENDIF
      IF (SIZE(V,2)<SIZE(RHO,2)) THEN
         WRITE(*,*) 'POT: 2nd dimension of V too small',SIZE(V,2),SIZE(RHO,2)
      ENDIF
      ISPIN=SIZE(RHO,3)
      IF (ISPIN/=1.AND.ISPIN/=2.AND.ISPIN/=4) THEN
         WRITE(*,*) 'POT: ISPIN /= 1,2, or 4:',ISPIN
         STOP
      ENDIF
      NMAX=SIZE(RHO,1)
      IF (NMAX<R%NMAX) THEN
         WRITE(*,*) 'POT: Grid inconsistency (1):',R%NMAX,NMAX
         STOP
      ENDIF
      IF (PRESENT(RHOC)) THEN
         IF (SIZE(RHOC)<R%NMAX) THEN
            WRITE(*,*) 'POT: Grid inconsistency (2):',R%NMAX,NMAX,SIZE(RHOC)
            STOP
         ENDIF
      ENDIF
      IF (PRESENT(VC)) THEN
         IF (SIZE(VC)<R%NMAX) THEN
            WRITE(*,*) 'POT: Grid inconsistency (3):',R%NMAX,NMAX,SIZE(VC)
            STOP
         ENDIF
      ENDIF

      SCALE=2*SQRT(PI0)

      V=0
      ISPIN = 1

!      IF (ISPIN==1.OR.ISPIN==2) THEN
         ALLOCATE(WORK1(NMAX,(LMAX+1)*(LMAX+1),ISPIN),WORK2(NMAX,ISPIN))
         WORK1=RHO
!      ELSEIF (ISPIN==4) THEN
!         ALLOCATE(WORK1(NMAX,(LMAX+1)*(LMAX+1),2),WORK2(NMAX,2))
!         CALL RAD_MAG_DENSITY(RHO,WORK1,LMAX,R)
!      ENDIF

!========================================================================
! Hartree potential
!========================================================================
      DHARTREE=0

      DO L=0,LMAX
      DO M=0,2*L
         LM=L*L+M+1
         CALL RAD_POT_HAR(L,R,V(:,LM,1),WORK1(:,LM,1),SUM)
         DHARTREE=DHARTREE+SUM
      ENDDO
      ENDDO
!      WRITE(6,*) 'DHARTREE=', DHARTREE
!      WRITE(6,*) 'V=', V(:,1,1)

!#ifdef testout
!      ! integrate charge
      CALL SIMPI(R,WORK1(:,1,1),QINT)
!      ! infer integrated charge from Hartree potential
      QINF=V(R%NMAX,1,1)*R%R(R%NMAX)/FELECT/2/SQRT(PI0)
!      ! write
      WRITE(*,'(A,F14.7,A,F14.7)') 'POT: Q_int=',QINT*SCALE,' Q_inf=',QINF
!#endif      

      IF (PRESENT(VC)) THEN
         DO K=1,R%NMAX
            V(K,1,1)=V(K,1,1)+VC(K)*SCALE
         ENDDO
      ENDIF
!      IF (ISPIN==2.OR.ISPIN==4) V(:,:,2)=V(:,:,1)

!========================================================================
! add nuclear potential
!========================================================================

      DO K=1,R%NMAX
         V(K,1,1)=V(K,1,1)-FELECT*SCALE*Z/R%R(K)
         WRITE(77, *) R%R(K), V(K,1,1)
      ENDDO
!      IF (ISPIN==2.OR.ISPIN==4) V(:,1,2)=V(:,1,1)

!========================================================================
! LDA exchange correlation energy, potential
! and double counting corrections
!========================================================================
      DEXC=0
      DVXC=0
      
      DO ISP=1,MIN(ISPIN,2)
         DO K=1,R%NMAX
            WORK2(K,ISP)=WORK1(K,1,ISP)/(SCALE*R%R(K)*R%R(K))
         ENDDO
      ENDDO
      ! add partial core charge if present
      IF (PRESENT(RHOC)) THEN
         DO K=1,R%NMAX
            WORK2(K,1)=WORK2(K,1)+RHOC(K)/(SCALE*R%R(K)*R%R(K))
         ENDDO
      ENDIF
      
 lda: IF (TLDA) THEN
            CALL RAD_LDA_XC(R,TREL,LMAX,WORK2(:,1),WORK1(:,:,1),V(:,:,1),DEXC,DVXC,.TRUE.)
      ENDIF lda
!========================================================================
! GGA if required
!========================================================================
      DEXC_GGA=0
      DVXC_GGA=0
      
! gga: IF (ISGGA()) THEN
  gga: IF (ADDGGA) THEN
            CALL RAD_GGA_XC(R,TLDA,WORK2(:,1),WORK1(:,1,1),V(:,1,1),DEXC_GGA,DVXC_GGA)
      ENDIF gga
!========================================================================
! done
!========================================================================
      ! Classical double counting correction:
      ! E_dc = -1/2 \int rho(r) V_H[rho(r)] dr + E_xc[rho+rhoc]
      !          -  \int rho(r) V_xc[rho(r)+rhoc(r)] dr
!      DOUBLEC= -DHARTREE/2+DEXC-DVXC+DEXC_GGA-DVXC_GGA

      EXCG= DEXC+DEXC_GGA
      IF (PRESENT(EXC)) EXC=EXCG
      
      IF (ISPIN==2) THEN
         WORK1(1:R%NMAX,:,1)=(V(1:R%NMAX,1:SIZE(RHO,2),1)+V(1:R%NMAX,1:SIZE(RHO,2),2))/2
         WORK1(1:R%NMAX,:,2)=(V(1:R%NMAX,1:SIZE(RHO,2),1)-V(1:R%NMAX,1:SIZE(RHO,2),2))/2
         V=WORK1
      ENDIF
!      IF (ISPIN==4) CALL RAD_MAG_DIRECTION(RHO,WORK1,V,LMAX,R)
      
      DEALLOCATE(WORK1,WORK2)
      RETURN
      END SUBROUTINE VASP_POT


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


      END MODULE VASP_POTCAR
