      MODULE VASP_POTCAR
         USE pseudo         ! VASP_Pseudo
         USE base           ! VASP_base
         USE ini            ! VASP_ini/PREC
         USE PSEUDO_struct  ! VASP_Pseudo_struct
         USE GlobalMath
         USE atomdata
         USE aeatom
         USE excor
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
      SUBROUTINE vasp_pseudo(ifinput, ifen, Grid0, Orbit, PotAE, Pot, success)!(ifinput,ifen,Grid)
         IMPLICIT COMPLEX(q)   (C)
         IMPLICIT REAL(q)  (A-B,D-H,O-Z)

         TYPE(GridInfo), INTENT(IN) :: Grid0
!         REAL(8), INTENT(IN) :: coreden(:)
         Type(OrbitInfo), INTENT(IN) :: Orbit
         TYPE(PotentialInfo), INTENT(IN) :: PotAE
         TYPE(PotentialInfo), INTENT(IN) :: Pot

         TYPE(potcar), TARGET, ALLOCATABLE :: P(:)
         TYPE(potcar), POINTER :: PP
         TYPE(INFO_STRUCT) :: INFO
!         TYPE(PseudoInfo) :: PAW

         INTEGER :: ifinput,ifen,Z
         INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS, LMAX, LMMAX, LI
         INTEGER :: LMAX_TABLE
         REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
         REAL(q) :: QTEST,SCALE
         REAL(q), ALLOCATABLE :: RHO(:,:,:), V(:,:,:), RHOAE00(:), RHOV(:)
         REAL(q), ALLOCATABLE :: CRHODE(:,:)
         REAL(q), ALLOCATABLE :: RHOLM(:)
         CHARACTER(LEN=2) :: TYPE(1)
         LOGICAL ::   LPAW,  success

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

         IU8 = 12
        IU10 = 14
        IU12 = 16
        IU13 = 17
        IU14 = 18
        IU16 = 20

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
      LMMAX = (PP%LMAX+1)**2
      ALLOCATE(RHO(PP%R%NMAX, LMMAX,1), V(PP%R%NMAX, LMMAX,1), RHOAE00(PP%R%NMAX))
      ALLOCATE(CRHODE(LDIM,LDIM))
      ALLOCATE(RHOLM(SIZE(CRHODE,1)*SIZE(CRHODE,1)))
      RHO = 0
!      DO io =1,CHANNELS
!         LI = PP%LPS(io)
!!         WRITE(6,*) 'OCC1=', PP%QATO(io,io)
!         RHOAE00(:) = RHOAE00(:)+PP%WAE(:,io)**2*PP%QATO(io,io)*(2*LI+1)
!      ENDDO

      LMAX_TABLE=6; CALL YLM3ST_(LMAX_TABLE)
      CALL SET_CRHODE_ATOM(CRHODE,PP)
      CALL TRANS_RHOLM(CRHODE, RHOLM, PP)
      CALL RAD_CHARGE(RHO(:,:,1), PP%R, RHOLM(:),PP%LMAX, PP%LPS, PP%WAE)
      RHOAE00(:) = RHO(:,1,1)

      OPEN(UNIT=15,FILE='VASP_VAL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=15,FILE='VASP_VAL',STATUS='OLD')
      ENDIF
          DO j=1, PP%R%NMAX
             WRITE(IU11,*) PP%R%R(j), RHOAE00(j) 
          ENDDO
!      RHO(:,1,1)=RHOAE00(:)+PP%RHOAE(:)

!      Z=INT(PP%ZVAL_ORIG+PP%ZCORE)
!      CALL POT(RHO, Z, PP%R, V)
!      RHO(:,1,1)=RHOAE00(:)

      CALL SIMPI(PP%R,RHOAE00, QTEST)
      WRITE(6,*) 'QTEST=', QTEST

      V=0
      CALL VASP_POT(RHO, 1, PP%R, V, PP%RHOAE)

      OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='OLD')
      ENDIF
          DO j=1, PP%R%NMAX
             WRITE(IU15,*) PP%R%R(j), -V(j,1,1)*SCALE/RYTOEV*2.0, PP%POTAE(j) 
          ENDDO


!!        WRITE(6,*) 'N=', Grid0%n
          call SetPOTAE(Grid0, FC%coreden, FC%valeden, Pot)
      OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='OLD')
      ENDIF
          DO j=1,Grid0%n 
             WRITE(IU16,*) Grid0%r(j)*AUTOA,(-Pot%rv(j))/(Grid0%r(j))
          ENDDO

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
             WRITE(IU12,*) Grid0%r(j)*AUTOA, PAW%tcore(j)*AUTOA
          ENDDO
!
         Call SetPAWOptions2(ifinput,ifen,Grid0, Orbit,Pot,success)

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
!
      ALLOCATE(RHOV(Grid0%n))
      RHOV = 0
      DO io =1, PAW%nbase
!         WRITE(6,*) 'OCC2=', Orbit%occ(io)
         if (mod(io,2) .ne. 0) RHOV(:) = RHOV(:)+PAW%phi(:,io)**2*Orbit%occ(io)
      ENDDO
      OPEN(UNIT=16,FILE='ATOM_VAL',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=16,FILE='ATOM_VAL',STATUS='OLD')
      ENDIF
          DO j=1,Grid0%n 
             WRITE(IU12,*) Grid0%r(j)*AUTOA, RHOV(j)*AUTOA, FC%valeden(j)
          ENDDO

         Call Report_Pseudobasis(Grid0,PAW,ifen)
         Call Set_PAW_MatrixElements(Grid,PAW)
         CALL logderiv(Grid,FCPot,PAW)
         CALL ftprod(Grid)

        CALL FindVlocfromVeff(Grid,FCOrbit,PAW)

        CALL Report_Pseudopotential(Grid,PAW)

        CALL SPMatrixElements(Grid,FCPot,FC,PAW)

!        deallocate(coreden)
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
        RETURN
! 
        END SUBROUTINE vasp_pseudo

        SUBROUTINE SetPOTAE(Grid2, coreden, valeden, POT)
         USE atomdata
         USE aeatom
         USE atomdata
        TYPE(GridInfo), INTENT(IN) :: Grid2
        TYPE(PotentialInfo) :: POT
        REAL(8), INTENT(IN) :: coreden(:)
        REAL(8), INTENT(IN) :: valeden(:)
        REAL(8), ALLOCATABLE :: density(:), Vrxc_tmp(:)
        REAL(8) qc, ecoul, v0, etxc, eex, SCALE

        SCALE = 2.0*sqrt(PI0)

        call InitPot(POT, Grid2%n)

        call Get_Nuclearpotential(Grid2, AEPot)
!        WRITE(6,*) 'RVN=', AEPot%rvn(Grid2%n)

        allocate (density(Grid2%n), Vrxc_tmp(Grid2%n))
        density = 0.d0
        DO i = 1, Grid2%n
!           density(i) = coreden(i)+FCOrbit%den(i)
           density(i) = coreden(i)+valeden(i)
        ENDDO

!        call poisson(Grid2, qc, coreden, POT%rvh, ecoul, v0)
        call poisson(Grid2, qc, valeden, POT%rvh, ecoul, v0)
!        call poisson(Grid2, qc, density, POT%rvh, ecoul, v0)
        WRITE(6,*) 'qc, ecoul, v0=', qc,ecoul, v0, POT%rvh(Grid2%n)

        call exch(Grid2, density, POT%rvx, etxc,eex)




!        DO i = 1, Grid2%n
!           Vrxc_tmp(i) =  POT%rvx(i)
!        ENDDO
!        call exch(Grid2, coreden, POT%rvx, etxc,eex)
!        call exch(Grid2, valeden, POT%rvx, etxc,eex)
!        DO i = 1, Grid2%n
!           POT%rvx(i) =  Vrxc_tmp(i)-POT%rvx(i)
!        ENDDO

!         POT%rv=(POT%rvh+POT%rvx)*SCALE
         POT%rv=AEPot%rvn+POT%rvh+POT%rvx

        deallocate(density, Vrxc_tmp)
        END SUBROUTINE SetPOTAE

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
      WRITE(6,*) 'DHARTREE=', DHARTREE
      WRITE(6,*) 'V=', V(:,1,1)

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
      IF (ISPIN==2.OR.ISPIN==4) V(:,:,2)=V(:,:,1)

!========================================================================
! add nuclear potential
!========================================================================

      DO K=1,R%NMAX
         V(K,1,1)=V(K,1,1)-FELECT*SCALE*Z/R%R(K)
      ENDDO
      IF (ISPIN==2.OR.ISPIN==4) V(:,1,2)=V(:,1,1)

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
