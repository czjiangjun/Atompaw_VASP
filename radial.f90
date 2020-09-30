!#include "symbol.inc"
!#define vector
!*******************************************************************
! RCS:  $Id: radial.F,v 1.11 2003/06/27 13:22:22 kresse Exp kresse $
!
!  MODULE which supports operations on radial grid 
!  all routines written by gK with the exception of the
!  routines required for metaGGA, which were written by
!  Robin Hirsch, Dec 2000
!
!*******************************************************************
  MODULE radial
    USE prec

    ! structure which is used for the logarithmic grid
    ! the grid points are given by R(i) = RSTART * exp [H (i-1)]
    TYPE rgrid
       REAL(q)  :: RSTART          ! starting point
       REAL(q)  :: REND            ! endpoint
       REAL(q)  :: RMAX            ! radius of augmentation sphere
       REAL(q)  :: D               ! R(N+1)/R(N) = exp(H)
       REAL(q)  :: H               !
       REAL(q),POINTER :: R(:)     ! radial grid (r-grid)
       REAL(q),POINTER :: SI(:)    ! integration prefactors on r-grid
       INTEGER  :: NMAX            ! number of grid points
    END TYPE rgrid

    ! This parameter determines at which magnetization the aspherical contributions
    ! to the one center magnetization are truncated in the non collinear case
    !   Without any truncation the aspherical terms for non magnetic atoms
    ! tend to yield spurious but meaningless contributions to the potential
    ! so that convergence to the groundstate can not be achieved
    ! for details see the routines RAD_MAG_DIRECTION and RAD_MAG_DENSITY
    REAL(q), PARAMETER :: MAGMIN=1E-2

    ! for non collinear calculations, setting
    ! the parameter USE_AVERAGE_MAGNETISATION  means that the aspherical 
    ! contributions to the one center magnetisation are projected onto the
    ! average magnetization direction in the PAW sphere instead of the
    ! local moment of the spherical magnetization density at
    ! each grid-point
    ! USE_AVERAGE_MAGNETISATION improves the numerical stability significantly
    ! and must be set
    LOGICAL :: USE_AVERAGE_MAGNETISATION=.TRUE.

    !
    ! parameters that determine the spherical integration mesh 
    ! for potentials in the GGA/LDA case
    !
    INTEGER, PARAMETER, PRIVATE :: PHPTS_FACT=3, THPTS_FACT=3

!#ifndef noPAW
  CONTAINS

!*******************************************************************
!
! RAD_ALIGN
! alligns R%RMAX to the radial grid
! this should be used to avoid that R%RMAX is slightly off
! from one of the grid points (due to rounding errors on reading)
!
!*******************************************************************

    SUBROUTINE RAD_ALIGN(R,RMAX)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (rgrid) R
      REAL(q) RMAX

      RMAX=RMAX-1E-5_q

      DO N=1,R%NMAX
         IF (R%R(N).GT. RMAX) THEN
            RMAX=R%R(N)
            EXIT
         ENDIF
      ENDDO

      CALL SET_SIMP(R)
    END SUBROUTINE RAD_ALIGN

!*******************************************************************
!
!  RAD_CHECK_QPAW
!  subroutine checks the consistency of the array QPAW
!  with the all electron and pseustored wavefunctions WAE and WPS
!    Q(PAW,ll') = \int (phi^AE(l,n,r)^2  -  phi^PS(l',n,r)^2) dr
!  if there are large errors the program stops and reports
!  the error
!  in all cases QPAW is corrected, so that it is exactly equal to
!  the integral
!
!  remember the wavefunction psi(r) can be obtained from the stored
!  coefficients using:
!      psi(r) = \sum_lmn Y_lm(r) w_ln(r) / r
!
!*******************************************************************

    SUBROUTINE RAD_CHECK_QPAW( R, CHANNELS, WAE, WPS , QPAW, QTOT, L, RDEP )
      IMPLICIT NONE

      TYPE (rgrid) R
      INTEGER CHANNELS
      REAL(q) :: WAE(:,:),WPS(:,:)   ! AE and soft wavefunctions
      REAL(q) :: QPAW(:,:,0:)        ! moments of compensation charge
      REAL(q) :: QTOT(:,:)           ! L=0 moments of AE charge
      INTEGER :: L(:)
      REAL(q) :: RDEP                ! outermost radius of augmentation charge 
! local variables
      REAL(q) :: RHOT(R%NMAX)
      INTEGER CH1,CH2,I
      REAL(q) :: RES
      INTEGER :: LL,LLP,LMIN,LMAX,LMAIN
      INTEGER :: IRMAX,NMAX_STORE

      DO IRMAX=1,R%NMAX-1
         ! because of rounding errors we need to be somewhat sloppy for the comparison
         IF (RDEP>0 .AND. R%R(IRMAX)-RDEP > -5E-3) EXIT
      ENDDO

      ! set the simpson weights in accordance with IRMAX
      NMAX_STORE=R%NMAX; R%NMAX=IRMAX
!      WRITE(*,*) NMAX_STORE, IRMAX
      CALL SET_SIMP(R)

      DO CH1=1,CHANNELS
      DO CH2=1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
         IF (LL==LLP) THEN
            RHOT=0
            DO I=1,IRMAX
               RHOT(I)=WAE(I,CH1)*WAE(I,CH2)
            ENDDO
            CALL SIMPI(R, RHOT , RES)
            QTOT(CH1,CH2)=RES
         ENDIF
      ENDDO
      ENDDO

      ! reset the simpson weights
      R%NMAX=NMAX_STORE
      CALL SET_SIMP(R)

      DO CH1=1,CHANNELS
      DO CH2=1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         DO LMAIN=LMIN,LMAX,2
            DO I=1,R%NMAX
               RHOT(I)=(WAE(I,CH1)*WAE(I,CH2)-WPS(I,CH1)*WPS(I,CH2))*R%R(I)**LMAIN
            ENDDO
            CALL SIMPI(R, RHOT , RES)
         ! screwed if we do not have the correct integrated charge
            IF ( ABS(RES-QPAW(CH1,CH2,0)) > 1E-4 .AND. LMAIN==0 ) THEN
               WRITE(0,1) CH1,CH2,RES,QPAW(CH1,CH2,0),RHOT(R%NMAX-4:R%NMAX)
 1             FORMAT('internal error RAD_CHECK_QPAW: QPAW is incorrect',/ &
               '      channels',2I3,' QPAW=',E20.10,' int=',10E20.10)
               STOP
            ENDIF
            QPAW(CH1,CH2,LMAIN)=RES
         ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE RAD_CHECK_QPAW


!*******************************************************************
!
!  RAD_AUG_CHARGE
!  add the compensation charge density on the radial grid
!  to RHO
!  RHO(L,M,r) = RHO(L,M,r) +  Q(r,L) Q(PAW,ll' L) RHOLM(ll',LM)
!
!  Q(r,L) are the L-dependent 1-normalized compensation charges
!
!*******************************************************************

!*******************************************************************
!
!  RAD_CHARGE
!  calculate the soft pseudo/AE charge density on the radial grid
!  from a set of frozen wavefunctions and the ll LM dependent
!  occupancies
!  the normalised wavefunction psi(r) can be obtained from the stored 
!  coefficients w(l,n,r) = W(:,:) in the following way
!      psi(r) = \sum_lmn Y_lm(r) w_ln(r) / r
!  here 
!'      rho_LM(r)= \sum_n w_ln(r) w_l'n'(r)  RHOLM(ll',LM)
!  is calculated and stored in RHO(2l+1+m,..) l=0,..,LMAX, m=0,..,2*l
!  RHOLM is the occupancy of each channel (see TRANS_RHOLM)
!  thus the charge density can be obtained from the stored
!  coefficients rho_lm(r) in the following way
!     rho(r) =  \sum_lm rho_lm(r) * Y_lm(r)  / r^2
!
!*******************************************************************

    SUBROUTINE RAD_CHARGE( RHOAE, R, RHOLM, CHANNELS, L, W )
      IMPLICIT NONE

      REAL(q) :: RHOAE(:,:)
      TYPE (rgrid) R
      INTEGER CHANNELS, L(:)
      REAL(q) :: RHOLM(:)
      REAL(q) :: W(:,:)
! local variables
      REAL(q) :: RHOT(R%NMAX)
      INTEGER CH1,CH2,LL,LLP,I
      INTEGER IBASE,JBASE,LMIN,LMAX,LMAIN,MMAIN,LMMAIN
   ! loop over all channels (l,epsilon)
      IBASE=0

      DO CH1=1,CHANNELS
      DO CH2=CH1,CHANNELS
         DO I=1,R%NMAX
            RHOT(I)=W(I,CH1)*W(I,CH2)
         ENDDO
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependent charge
         ! loop could be replaced by matrix vector DGEMV
         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DO I=1,R%NMAX
               RHOAE(I,LMMAIN)=RHOAE(I,LMMAIN)+RHOT(I)*RHOLM(LMMAIN+JBASE)
            ENDDO
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      ENDDO
      ENDDO

    END SUBROUTINE RAD_CHARGE



!*******************************************************************
!
!  FLIP_RAD
!  Flips an arbitrary real(!) array from total, magnetization to
!  up,down spin
!   
!  Robin Hirschl 20010109
!*******************************************************************
    SUBROUTINE FLIP_RAD(WORKIN,WORKOUT,N)
      IMPLICIT NONE

      REAL(q) :: WORKIN(:,:),WORKOUT(:,:)
      INTEGER I,N
      REAL(q) TEMP

      DO I=1,N
         TEMP=WORKIN(I,1)
         WORKOUT(I,1)=(WORKIN(I,1)+WORKIN(I,2))/2._q
         WORKOUT(I,2)=(TEMP-WORKIN(I,2))/2._q
      ENDDO
    END SUBROUTINE FLIP_RAD

!*******************************************************************
!
!  RAD_PROJ
!  calculate
!'  D(ll',LM)= D(ll',LM)+ A *\int dr phi_l'(r) pot_lm(r) phi_l(r) dr
!  on the radial grid
!  the potential is given by
!       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
!*******************************************************************

    SUBROUTINE RAD_PROJ( POT, R, A, DLM, CHANNELS, L, W )
      IMPLICIT NONE

      REAL(q) :: POT(:,:)     ! radial potential V(r,L,M)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q) :: W(:,:)       ! wavefunctions phi(r,l)
      REAL(q) :: A            ! scaling factor
      INTEGER CHANNELS, L(:)
! local variables
      REAL(q) :: RHOT(R%NMAX),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I
      INTEGER IBASE,JBASE,LMIN,LMAX,LMAIN,MMAIN,LMMAIN
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
         DO I=1,R%NMAX
            RHOT(I)=W(I,CH1)*W(I,CH2)
         ENDDO
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN

         DO LMAIN=LMIN,LMAX,2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            SUM=0

            ! integrate RHOT POT(L,M) (the potentials are already weighted)
            DO I=1,R%NMAX
               SUM=SUM+RHOT(I)*POT(I,LMMAIN)
            ENDDO
            DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)+SUM*A
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE RAD_PROJ

!*******************************************************************
!
!  RAD_PROJ_KINPOT
!  calculate the integral
!'   D(ll',LM) =  D(ll',LM) + A * < nabla phi_l' | V(L,M) | nabla phi_l >
!  on the radial grid
!  only spherical contributions are used
!
!*******************************************************************

    SUBROUTINE RAD_PROJ_KINPOT( POT, R, A, DLM, CHANNELS, L, W )
      IMPLICIT NONE

      REAL(q) :: POT(:)       ! radial kinetic energy potential V(r,L,M)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q) :: W(:,:)       ! wavefunctions phi(r,l)
      REAL(q) :: A            ! scaling factor
      INTEGER CHANNELS, L(:)
! local variables
      REAL(q) :: K(R%NMAX),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I
      INTEGER IBASE,JBASE,LMIN,LMAX,LMAIN,MMAIN,LMMAIN
   ! loop over all channels (l,epsilon)
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
         
         IF (LL == LLP) THEN

            ! here one has to calculate the spherical contribution
            ! to the kinetic energy density

            DO I=1,R%NMAX
               K(I)=W(I,CH1)*W(I,CH2)
            ENDDO

            ! Lmin and Lmax
            LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
            JBASE=IBASE-LMIN*LMIN

            LMAIN=0
            MMAIN=1

              ! currently LMMAIN is allways 1, and JBASE=IBASE
              LMMAIN=LMAIN*LMAIN+MMAIN
              ! integrate K *  POT (the potential is already weighted)
              SUM=0
              DO I=1,R%NMAX
                 SUM=SUM+K(I)*POT(I)
              ENDDO
              DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)+SUM*A
         ENDIF

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO

    END SUBROUTINE RAD_PROJ_KINPOT

!*******************************************************************
!
! calculate the gradient of a function on the radial grid
! nabla rho(r) = d/dr rho(r)
! (adopted from some routine of a program called atoms)
!
!*******************************************************************

      SUBROUTINE GRAD_(R,RH,DRH)
      IMPLICIT NONE
      TYPE (rgrid) R      ! defines the radial grid
      REAL(q) RH(R%NMAX)  ! charge density
      REAL(q) DRH(R%NMAX) ! gradient of the function
      REAL(q) H
      INTEGER NR,NM2,I

      H=R%H
      NR=R%NMAX
      NM2=NR-2
! For the first and second point the forward differences are used
      DRH(1)=((6._q*RH(2)+20._q/3._q*RH(4)+1.2_q*RH(6)) &
     &      -(2.45_q*RH(1)+7.5_q*RH(3)+3.75_q*RH(5)+1._q/6._q*RH(7)))/H
      DRH(2)=((6._q*RH(3)+20._q/3._q*RH(5)+1.2_q*RH(7)) &
     &        -(2.45_q*RH(2)+7.5_q*RH(4)+3.75_q*RH(6)+1._q/6._q*RH(8)))/H
! Five points formula
      DO  I=3,NM2
         DRH(I)=((RH(I-2)+8._q*RH(I+1))-(8._q*RH(I-1)+RH(I+2)))/12._q/H
      ENDDO
! Five points formula for the last two points ('backward differences')
      DRH(NR-1)=(-1._q/12._q*RH(NR-4)+0.5_q*RH(NR-3)-1.5_q*RH(NR-2) &
     &           +5._q/6._q*RH(NR-1)+0.25_q*RH(NR))/H
      DRH(NR)=  (0.25_q*RH(NR-4)-4._q/3._q*RH(NR-3)+3._q*RH(NR-2) &
     &           -4._q*RH(NR-1)+25._q/12._q*RH(NR))/H
! account for logarithmic mesh
      DO  I=1,NR
         DRH(I)=DRH(I)/R%R(I)
      ENDDO

      RETURN
      END SUBROUTINE

!*******************************************************************
!
! calculate the gradient of a function on the radial grid
! using 2.order differentiation (this is good enough 
! and less suseptible to noise than the previous routine)
! nabla rho(r) = d/dr rho(r)
!
!*******************************************************************

      SUBROUTINE GRAD(R,RH,DRH)
      IMPLICIT NONE
      TYPE (rgrid) R      ! defines the radial grid
      REAL(q) RH(R%NMAX)  ! charge density
      REAL(q) DRH(R%NMAX) ! gradient of the function
      REAL(q) H
      INTEGER NR,NM1,I

      H=R%H
      NR=R%NMAX
      NM1=NR-1
! 1. point use first order differantion
      DRH(1)=(RH(2)-RH(1))/H
! three point formula
      DO  I=2,NM1
         DRH(I)=(RH(I+1)-RH(I-1))/2/H
      ENDDO
! last point
      DRH(NR)=(RH(NR)-RH(NR-1))/H

      DO  I=1,NR
         DRH(I)=DRH(I)/R%R(I)
      ENDDO

      RETURN
      END SUBROUTINE

!*******************************************************************
!
!  AUG_SETE
!  calculate coefficients for generalized gaussian
!   exp ( r^2 / alpha)
!  not used
!*******************************************************************


    SUBROUTINE AUG_SETE(R, ALPHA, A, TH)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (rgrid) R

      REAL(q)  ALPHA, A, TH
      REAL(q)  TMP(R%NMAX)

      RC= R%RMAX
      ALPHA= LOG( TH) / (RC*RC)

      DO N=1,R%NMAX
         TMP(N)=0
         IF (R%R(N) <= RC) THEN
            VAL=EXP(ALPHA*R%R(N)*R%R(N))
            TMP(N)=R%R(N)*R%R(N)*VAL
         ENDIF
      ENDDO
      CALL SIMPI(R,TMP,SUM)
      A = 1/SUM

    END SUBROUTINE AUG_SETE

!*******************************************************************
!
!  AUG_SETQ
!  find a set of 2 q_i and coefficients A_i such that
!     sum_i  A_i j_l(q_i,Rc)  =0
!     sum_i  A_i j_l(q_i,Rc)' =0'
!     sum_i  A_i j_l(q_i,Rc)''=0
!  and
!     sum_i int_0^Rc A_i j_l(q_i r) r^(l+2) dr = 1
!  the augmentation charge is then represented by
!     sum_i  A_i j_l(q_i, r)
!  PS: I confess I have an obsession with spherical Besselfunctions
!*******************************************************************


    SUBROUTINE AUG_SETQ(L,R,RMAX,QQ,A,LCOMPAT)
      IMPLICIT REAL(q) (A-H,O-Z)
      TYPE (rgrid) R

      PARAMETER (NQ=2)
      REAL(q)  QQ(NQ)
      REAL(q)  A(NQ)
      LOGICAL LCOMPAT
      REAL(q)  AMAT(NQ,NQ)
      REAL(q)  B(3)
      INTEGER  IPIV(NQ)

      REAL(q)  TMP(R%NMAX)
!-----gaussian integration
      PARAMETER (M=32)
      REAL(q)  WR(M),RA(M)
      INTEGER IFAIL
      EXTERNAL GAUSSI2
!-----------------------------------------------------------------------
! search for q so that j(q Rc)=0
!-----------------------------------------------------------------------
      CALL AUG_BEZERO(QQ,L,NQ)
!-----------------------------------------------------------------------
! second set the matrix
!-----------------------------------------------------------------------
      ITYPE=0
      CALL GAUSSI(GAUSSI2,0.0_q,RMAX, ITYPE,M,WR,RA,IFAIL)

      DO I=1,NQ
         QQ(I)=QQ(I)/RMAX
         QR  =QQ(I)*RMAX
         CALL SBESSE3( QR, BJ, BJP, BJPP, L)
         BJP =BJP *QQ(I)
         BJPP=BJPP*QQ(I)*QQ(I)
         AMAT(1,I) = BJP
!  I could not find an analytical expression for \int j_l(qr) r^(2+l) dr
!  I will check again probably it is simple enough
         IF (LCOMPAT) THEN
           DO N=1,R%NMAX
              IF (R%R(N) < RMAX) THEN
                 QR=QQ(I)*R%R(N)
                 CALL SBESSEL(QR,BJ,L)
                 TMP(N)=BJ*R%R(N)**(2+L)
              ELSE
                 TMP(N)=0
              ENDIF
           ENDDO
           CALL SIMPI(R,TMP,SUM)
           AMAT(2,I) = SUM
         ELSE

! Gauss integration, more accurate !
           SUM2=0
           DO N=1,M
              QR=QQ(I)*RA(N)
              CALL SBESSEL(QR,BJ,L)
              SUM2=SUM2+BJ*RA(N)**(2+L)*WR(N)
           ENDDO

           AMAT(2,I) = SUM2
         ENDIF
      ENDDO
!-----------------------------------------------------------------------
!  solve the linear equations  B_i = AMAT_ii' A_i'
!-----------------------------------------------------------------------
      IFAIL=0
      B(1)=0
      B(2)=1

      A(1)=0
      A(2)=0

      CALL DGETRF( NQ, NQ, AMAT, NQ, IPIV, IFAIL )
      CALL DGETRS('N', NQ, 1, AMAT, NQ, IPIV, B, NQ, IFAIL)
      A=B(1:2)

      B(1)=0
      B(2)=0
      B(3)=0

      DO I=1,NQ
         SUM=0
         QR  =QQ(I)*RMAX
         CALL SBESSE3( QR, BJ, BJP, BJPP, L)
         B(1)=B(1)+BJ*A(I)
         B(2)=B(2)+BJP*A(I)*QQ(I)
         B(3)=B(3)+BJPP*A(I)*QQ(I)*QQ(I)
      ENDDO
    END SUBROUTINE

!*******************************************************************
!  SUBROUTINE BEZERO
!  searches for NQ zeros j(qr)
!  i/o:
!         XNULL(NQ) result
!         L           quantum number l
!  great-full spaghetti code (written by gK)
!********************************************************************

    SUBROUTINE AUG_BEZERO(XNULL,L,NQ)
      IMPLICIT REAL(q) (A-H,O-Z)
      PARAMETER (STEP=.1_q, BREAK= 1E-10_q)
      DIMENSION XNULL(NQ)
! initialization
      X=STEP
      N=0
! entry point for next q_n
  30  CALL SBESSE2(X, BJ1, DUMMY,  L)
! coarse search
  10  X=X+STEP
      CALL SBESSE2(X, BJ2, DUMMY,  L)
! found one point
      IF(BJ1*BJ2 < 0) THEN
        ETA=0.0_q
! interval bisectioning
        SSTEP=STEP
        XX   =X
  20    SSTEP=SSTEP/2
        IF (BJ1*BJ2 < 0) THEN
          XX=XX-SSTEP
        ELSE
          XX=XX+SSTEP
        ENDIF
        CALL SBESSE2(XX, BJ2, DUMMY,  L)
        IF (SSTEP > BREAK) GOTO 20

        N=N+1
        XNULL(N)=XX
        IF (N == NQ) RETURN
        GOTO 30
      ENDIF
      GOTO 10

    END SUBROUTINE

!********************************************************************
!
!  SUBROUTINE SIMPI
!  integrate a function on the logarithmic grid
!  uses the previously setup weights
!
!********************************************************************

    SUBROUTINE SIMPI(R,F,FI)
      IMPLICIT NONE
      TYPE (rgrid) R
      REAL(q)  F(:),FI,SUM
      INTEGER   K

      SUM=0
!OCL SCALAR
      DO K=1,R%NMAX
         SUM=SUM+F(K)*R%SI(K)
      ENDDO

      FI=SUM

    END SUBROUTINE

!********************************************************************
!
!  SUBROUTINE SET_SIMP
!  setup weights for simpson integration on radial grid
!  any radial integral can then be evaluated by just summing all
!  radial grid points with the weights SI
!
!  int dr = sum_i si(i) * f(i)
!  the factors  R%R(K)  *R%H stem from the logarithmic grid
!********************************************************************

    SUBROUTINE SET_SIMP(R)
      IMPLICIT NONE
      TYPE (rgrid) R
      INTEGER   K

      ALLOCATE(R%SI(R%NMAX))
      R%SI=0
      DO K=R%NMAX,3,-2
         R%SI(K)=    R%R(K)  *R%H/3._q+R%SI(K)
         R%SI(K-1)=4*R%R(K-1)*R%H/3._q
         R%SI(K-2)=  R%R(K-2)*R%H/3._q
      ENDDO
    END SUBROUTINE
!#endif
  END MODULE radial
