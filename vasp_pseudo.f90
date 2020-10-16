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

         TYPE(potcar), ALLOCATABLE :: P(:)
         TYPE(INFO_STRUCT) :: INFO
!         TYPE(PseudoInfo) :: PAW

         INTEGER :: ifinput,ifen
         INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS
         REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
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

!      OPEN(UNIT=7,FILE='VASP_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
!      OPEN(UNIT=8,FILE='VASP_WPS',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=7,FILE='VASP_WAE',STATUS='OLD')
!         OPEN(UNIT=8,FILE='VASP_WPS',STATUS='OLD')
!      ENDIF
!         DO i=1, CHANNELS
!          DO j=1, P(1)%R%NMAX
!          WRITE(IU6,*) P(1)%R%R(j), P(1)%WAE(j,i)
!!          if (mod(i,2) == 0)WRITE(IU6,*) P(1)%R%R(j), P(1)%WAE(j,i)
!          WRITE(IU7,*) P(1)%R%R(j), P(1)%WPS(j,i)
!          ENDDO
!          WRITE(IU6,*)
!          WRITE(IU7,*)
!         ENDDO
!
!      OPEN(UNIT=13,FILE='VASP_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=13,FILE='VASP_CORE',STATUS='OLD')
!      ENDIF
!          DO j=1, P(1)%R%NMAX
!             WRITE(IU9,*) P(1)%R%R(j), P(1)%RHOAE(j)
!          ENDDO
!
!      OPEN(UNIT=15,FILE='VASP_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=15,FILE='VASP_PCORE',STATUS='OLD')
!      ENDIF
!          DO j=1, P(1)%R%NMAX
!             WRITE(IU11,*) P(1)%R%R(j), P(1)%RHOPS(j)
!          ENDDO

      OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=19,FILE='VASP_POTAE',STATUS='OLD')
      ENDIF
          DO j=1, P(1)%R%NMAX
             WRITE(IU15,*) P(1)%R%R(j), P(1)%POTAE(j)!, P(1)%POTPS(j)
          ENDDO

!!        WRITE(6,*) 'N=', Grid0%n
          call SetPOTAE(Grid0, FC%coreden, FC%valeden, Pot)
      OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=20,FILE='ATOM_POTAE',STATUS='OLD')
      ENDIF
          DO j=1,Grid0%n 
             WRITE(IU16,*) Grid0%r(j)*AUTOA, (-Pot%rv(j))/(Grid0%r(j))
          ENDDO

         Call setcoretail2(Grid0, FC%coreden)
!!         Call setcoretail2(Grid0, coreden)
!      OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='OLD')
!      ENDIF
!          DO j=1,Grid0%n 
!             WRITE(IU10,*) Grid0%r(j)*AUTOA, FC%coreden(j)*AUTOA
!          ENDDO
!
!      OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=16,FILE='ATOM_PCORE',STATUS='OLD')
!      ENDIF
!          DO j=1,Grid0%n 
!             WRITE(IU12,*) Grid0%r(j)*AUTOA, PAW%tcore(j)*AUTOA
!          ENDDO
!
         Call SetPAWOptions2(ifinput,ifen,Grid0, Orbit,Pot,success)

!      OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
!      OPEN(UNIT=17,FILE='ATOM_WPS',STATUS='UNKNOWN',IOSTAT=IERR)
!      IF (IERR/=0) THEN
!         OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='OLD')
!         OPEN(UNIT=17,FILE='ATOM_WPS',STATUS='OLD')
!      ENDIF
!         DO i=1, PAW%nbase
!          DO j=1, Grid0%n
!             WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
!!           if (mod(i,2) == 0) WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
!             WRITE(IU13,*) Grid0%r(j)*AUTOA, PAW%tphi(j,i)/sqrt(AUTOA)
!          ENDDO
!          WRITE(IU8,*)
!          WRITE(IU13,*)
!         ENDDO
!

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
        TYPE(GridInfo), INTENT(IN) :: Grid2
        TYPE(PotentialInfo) :: POT
        REAL(8), INTENT(IN) :: coreden(:)
        REAL(8), INTENT(IN) :: valeden(:)
        REAL(8), ALLOCATABLE :: density(:), Vrxc_tmp(:)
        REAL(8) qc, ecoul, v0, etxc, eex, SCALE

        SCALE = 2.0*sqrt(PI0)

        call InitPot(POT, Grid2%n)

        call Get_Nuclearpotential(Grid2, AEPot)
!        WRITE(6,*) AEPot%rvn(1:Grid2%n)

        allocate (density(Grid2%n), Vrxc_tmp(Grid2%n))
        density = 0.d0
        DO i = 1, Grid2%n
!           density(i) = coreden(i)+FCOrbit%den(i)
           density(i) = coreden(i)+valeden(i)
        ENDDO

!        call poisson(Grid2, qc, coreden, POT%rvh, ecoul, v0)
        call poisson(Grid2, qc, valeden, POT%rvh, ecoul, v0)
!        call poisson(Grid2, qc, density, POT%rvh, ecoul, v0)
!        WRITE(6,*) qc, v0

        call exch(Grid2, density, POT%rvx, etxc,eex)
!        DO i = 1, Grid2%n
!           Vrxc_tmp(i) =  POT%rvx(i)
!        ENDDO
!        call exch(Grid2, coreden, POT%rvx, etxc,eex)
!        call exch(Grid2, valeden, POT%rvx, etxc,eex)
!        DO i = 1, Grid2%n
!           POT%rvx(i) =  Vrxc_tmp(i)-POT%rvx(i)
!        ENDDO

         POT%rv=AEPot%rvn+(POT%rvh+POT%rvx)
!         POT%rv=AEPot%rvn+POT%rvh+POT%rvx)

        deallocate(density, Vrxc_tmp)
        END SUBROUTINE SetPOTAE
       END MODULE VASP_POTCAR
