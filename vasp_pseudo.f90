      SUBROUTINE vasp_pseudo(ifinput, ifen, Grid0, coreden, Orbit, Pot, success)!(ifinput,ifen,Grid)
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
         IMPLICIT COMPLEX(q)   (C)
         IMPLICIT REAL(q)  (A-B,D-H,O-Z)

         TYPE(potcar), ALLOCATABLE :: P(:)
         TYPE(INFO_STRUCT) :: INFO
         TYPE(GridInfo) :: Grid0
         TYPE(PotentialInfo) :: Pot
         Type(OrbitInfo):: Orbit
!         TYPE(PseudoInfo) :: PAW

         INTEGER :: ifinput,ifen
         INTEGER :: NTYP, NTYPD, LDIM, LDIM2, LMDIM,CHANNELS
         REAL(q) ZVALF(1),POMASS(1),RWIGS(1), VCA(1)  ! valence, mass, wigner seitz radius
         REAL(8),allocatable :: coreden(:)
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
         IU0 = 6
         IU6 = 8
         IU8 = 12
         IU9 = 13
        IU10 = 14
         LPAW = .TRUE.

         Call RD_PSEUDO(INFO, P, NTYP, NTYPD, LDIM, LDIM2, LMDIM, POMASS,     &
     &                  RWIGS, TYPE,                       &
     &                  CHANNELS, IU0, IU6, -1, LPAW)

      OPEN(UNIT=8,FILE='VASP_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=8,FILE='VASP_WAE',STATUS='OLD')
      ENDIF
         DO i=1, CHANNELS
          DO j=1, P(1)%R%NMAX
!          WRITE(IU6,*) P(1)%R%R(j), P(1)%WAE(j,i)
           if (mod(i,2) == 0)WRITE(IU6,*) P(1)%R%R(j), P(1)%WAE(j,i)
          ENDDO
          WRITE(IU6,*)
         ENDDO

      OPEN(UNIT=13,FILE='VASP_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=13,FILE='VASP_CORE',STATUS='OLD')
      ENDIF
          DO j=1, P(1)%R%NMAX
             WRITE(IU9,*) P(1)%R%R(j), P(1)%RHOAE(j)
          ENDDO

      OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=12,FILE='ATOM_WAE',STATUS='OLD')
      ENDIF
         DO i=1, PAW%nbase
          DO j=1, Grid0%n
!             WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
           if (mod(i,2) == 0) WRITE(IU8,*) Grid0%r(j)*AUTOA, PAW%phi(j,i)/sqrt(AUTOA)
          ENDDO
          WRITE(IU8,*)
         ENDDO
!
      OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='UNKNOWN',IOSTAT=IERR)
      IF (IERR/=0) THEN
         OPEN(UNIT=14,FILE='ATOM_CORE',STATUS='OLD')
      ENDIF
          DO j=1,Grid0%n 
             WRITE(IU10,*) Grid0%r(j)*AUTOA, FC%coreden(j)*AUTOA
          ENDDO

         allocate(coreden(Grid0%n))
!        WRITE(6,*) 'N=', Grid0%n
         Call setcoretail2(Grid0, coreden)

         Call SetPAWOptions2(ifinput,ifen,Grid0,Orbit,Pot,success)

         Call Report_Pseudobasis(Grid0,PAW,ifen)
         WRITE(6, *) "TEST"
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
        RETURN
! 
        END SUBROUTINE vasp_pseudo
