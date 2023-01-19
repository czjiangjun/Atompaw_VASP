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

    SUBROUTINE RAD_AUG_CHARGE( RHO, R, RHOLM, CHANNELS, L, &
         LYMAX, AUG, QPAW )
      USE ini
      IMPLICIT NONE

      REAL(q) :: RHO(:,:)    ! charge on radial grid
      TYPE (rgrid) R
      REAL(q) :: RHOLM(:)    ! occupancy of each llLM channel
      REAL(q) :: AUG(:,0:)   ! 1-normalized L-dep compensation charge
      REAL(q) :: QPAW(:,:,0:)! moments of compensation charge Q(PAW,ll L)
      INTEGER :: LYMAX       ! maximum L
      INTEGER CHANNELS, L(:)
   ! local variables
      REAL(q) :: RHOLMT((LYMAX+1)*(LYMAX+1))
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I,LMAX
      INTEGER IBASE,JBASE,LMIN,LMAIN,MMAIN,LMMAIN
!-----------------------------------------------------------------------
! first contract to L and M dependent Q(LM)
!-----------------------------------------------------------------------
     IBASE=0
      RHOLMT=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependent charge
         DO LMAIN=LMIN,MIN(LMAX,LYMAX),2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            RHOLMT(LMMAIN)=RHOLMT(LMMAIN)+QPAW(CH1,CH2,LMAIN)*RHOLM(LMMAIN+JBASE)
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
!-----------------------------------------------------------------------
! then add to charge on radial grid
!-----------------------------------------------------------------------
      ! could be replaced by matrix vecor DGEMV
      DO LMAIN=0,LYMAX
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DO I=1,R%NMAX
               RHO(I,LMMAIN)=RHO(I,LMMAIN)+AUG(I,LMAIN)*RHOLMT(LMMAIN)
            ENDDO
         ENDDO
      ENDDO
!#ifdef debug
!      WRITE(0,*) 'RAD_AUG_CHARGE: L compensation charges are'
!      DO LL=0,LYMAX
!         WRITE(0,'(I3,10F10.5)') LL,(RHOLMT(LL**2+I)*2*SQRT(PI),I=1,LL*2+1)
!      ENDDO
!#endif
    END SUBROUTINE RAD_AUG_CHARGE


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
!  RAD_AUG_PROJ
!  calculate the integral
!   D(ll'LM) =  \int V(r,L,M) Q(r,L) Q(PAW,ll' L) dr
!  on a radial grid
!  Q(r,L) are the L-dependent 1-normalized compensation charges
!  the potential is given by
!       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
!*******************************************************************

    SUBROUTINE RAD_AUG_PROJ( POT, R, DLM, CHANNELS, L, &
         LYMAX, AUG, QPAW )
      USE ini
      IMPLICIT NONE

      REAL(q) :: POT(:,:)
      TYPE (rgrid) R
      REAL(q) :: DLM(:)
      REAL(q) :: AUG(:,0:)    ! 1-normalized L-dep compensation charge
      REAL(q) :: QPAW(:,:,0:) ! moments of compensation charges Q(PAW,ll L)
      INTEGER :: LYMAX        ! maximum L
      INTEGER CHANNELS,L(:)
   ! local variables
      REAL(q) :: RHOLMT((LYMAX+1)*(LYMAX+1)),SUM
      INTEGER CH1,CH2,LL,LLP,LM,LMP,I,LMAX
      INTEGER IBASE,JBASE,LMIN,LMAIN,MMAIN,LMMAIN
!-----------------------------------------------------------------------
! first calculate \int V(L,M) Q(L,M)
!-----------------------------------------------------------------------
      RHOLMT=0
      ! could be replaced by matrix vecor DGEMV
      DO LMAIN=0,LYMAX
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            SUM=0
            DO I=1,R%NMAX
               SUM=SUM+POT(I,LMMAIN)*AUG(I,LMAIN)
            ENDDO
            RHOLMT(LMMAIN)=SUM
         ENDDO
      ENDDO
!-----------------------------------------------------------------------
! than multiply with QPAW(llp, L) and add the DLM
!-----------------------------------------------------------------------
      IBASE=0

      LM=1
      DO CH1=1,CHANNELS
      LMP=LM
      DO CH2=CH1,CHANNELS
      ! quantum numbers l and lp of these two channels
         LL =L(CH1)
         LLP=L(CH2)
      ! Lmin and Lmax
         LMIN=ABS(LL-LLP) ; LMAX=ABS(LL+LLP)
         JBASE=IBASE-LMIN*LMIN
      ! add to LM dependet charge
         DO LMAIN=LMIN,MIN(LYMAX,LMAX),2
         DO MMAIN=1,LMAIN*2+1
            LMMAIN=LMAIN*LMAIN+MMAIN
            DLM(LMMAIN+JBASE)=DLM(LMMAIN+JBASE)-RHOLMT(LMMAIN)*QPAW(CH1,CH2,LMAIN)
         ENDDO
         ENDDO

      IBASE=IBASE+(2*LL+1)*(2*LLP+1)
      LMP=LMP+2*LLP+1
      ENDDO
      LM =LM +2*LL +1
      ENDDO
!#ifdef debug
!      WRITE(0,*) 'RAD_AUG_PROJ: int V Q(LM) is'
!      DO LL=0,LYMAX
!         WRITE(0,'(I3,10F10.5)') LL,(RHOLMT(LL**2+I)*2*SQRT(PI),I=1,LL*2+1)
!      ENDDO
!#endif
    END SUBROUTINE RAD_AUG_PROJ

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

!*******************************************************************
!
!  RAD_POT
!  calculate the radial potential from the radial chargedensity
!  the charge density rho(r) is given by
!     rho(r) =  \sum_lm rho_lm(r) * Y_lm(r)  / r^2
!  the potential is given by
!       V(r) =  \sum_lm pot_lm(r) * Y_lm(r)
!
!  where rho_lm(r) is stored in RHO(2l+1+m,..) l=0,..,LMAX, m=0,..,2*l
!  and   pot_lm(r) is stored in POT(2l+1+m,..)
!
! in many places we use a scaling factor 2 sqrt(pi)
! the "real" charge density for L=0 angular quantum number is
!   n_0= rho_00 Y_00 = RHO(r,0) / (2 sqrt(pi))
! for other channels it is
!   n_lm= rho_lm Y_lm
! in comments we will always distinct between n_L and rho_L
!
! the charge density is supplied as 
!  total=RHO(:,:,1), magnetization=RHO(:,:,2)
! the potential however is returned for spin up and spin down 
!  spin up=POT(:,:,1), spin down=POT(:,:,2)
!
!*******************************************************************


    SUBROUTINE RAD_POT( R, ISPIN, LMAX, LMAX_CALC, LASPH, &
         RHO, RHOC, POTC, POT, DOUBLEC, EXCG)
      USE ini
      USE setexm

      IMPLICIT NONE
      INTEGER LMAX, ISPIN, LMAX_CALC
      LOGICAL :: LASPH          ! non spherical contributions
      REAL(q) :: RHOC(:)        ! core charge for exchange correlation
      REAL(q) :: POTC(:)        ! froze core potential
      REAL(q) :: RHO(:,:,:)     ! charge distribution see above
      ! RHO(:,:,1) contains total charge distribution
      ! RHO(:,:,2) magnetization charge distribution
      REAL(q) :: POT(:,:,:)     ! potential
      ! POT(:,:,1) up   component
      ! POT(:,:,2) down component
      TYPE (rgrid) :: R
      REAL(q) :: EXCG           ! exchange energy only
      REAL(q) :: DOUBLEC        ! double counting corrections
    ! local variables
      REAL(q) RHOT(R%NMAX,ISPIN)
      INTEGER K,N,I,L,M,LM
      REAL(q) SCALE,SUM
      
      LOGICAL,PARAMETER :: TREL=.TRUE. ! use relativistic corrections to exchange
      LOGICAL,PARAMETER :: TLDA=.TRUE. ! calculate LDA contribution seperately
      ! TLDA=.FALSE. works only for Perdew Burke Ernzerhof
      ! in this case non spherical contributions are missing
      REAL(q) :: DHARTREE,DEXC_LDA,DVXC_LDA,DEXC_GGA,DVXC_GGA
      REAL(q) :: TMP((LMAX+1)*(LMAX+1),ISPIN)
!
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DEXC1, DVXC1
      REAL(q) T1(R%NMAX),T2(R%NMAX),V1(R%NMAX)

      POT=0

      SCALE=2*SQRT(PI)
      N=R%NMAX

      DHARTREE=0
      DO L=0,LMAX_CALC
      DO M=0,2*L
         LM=L*L+M+1
         CALL RAD_POT_HAR(L,R,POT(:,LM,1),RHO(:,LM,1),SUM)
!         IF (ISPIN==2) POT(:,LM,2)=POT(:,LM,1)
         DHARTREE=DHARTREE+SUM
         TMP(LM,1)=SUM
      ENDDO
      ! WRITE(0,'(I2,10F12.7)') L, (TMP(LM,1),LM=L*L+1,(L+1)*(L+1))
      ENDDO
      DO K=1,N
         POT(K,1,1)=POT(K,1,1)+POTC(K)*SCALE
      ENDDO
!      IF (ISPIN==2) POT(:,1,2)=POT(:,1,1)
!========================================================================
! exchange correlation energy, potential
! and double counting corrections
!========================================================================
      DEXC_LDA=0
      DVXC_LDA=0
      DEXC_GGA=0
      DVXC_GGA=0

!      IF (ISPIN==1) THEN
        DO K=1,N
          RHOT(K,1)=(RHO(K,1,1)+RHOC(K))/ (SCALE*R%R(K)*R%R(K)) ! charge density rho_0 Y(0)
        ENDDO

!      ELSE
!        DO K=1,N
!         RHOT(K,1)=(RHO(K,1,1)+RHOC(K))/ (SCALE*R%R(K)*R%R(K)) ! charge density n_0=rho_0 Y(0)
!         RHOT(K,2)= RHO(K,1,2)/ (SCALE*R%R(K)*R%R(K))          ! magnetization
!        ENDDO

!      ENDIF
      ! WRITE(*,'(9F14.7)') POT(R%NMAX-10,:,:)
!========================================================================
! LDA+ GGA if required
!========================================================================
!      IF (.NOT. LASPH) THEN
       IF (TLDA) THEN
            
!            IF (ISPIN==1) THEN
               CALL RAD_LDA_XC( R, TREL, LMAX_CALC, RHOT(:,1), RHO(:,:,1), POT(:,:,1), DEXC_LDA, DVXC_LDA, .TRUE.)
!            ELSE
!               CALL RAD_LDA_XC_SPIN( R, TREL, LMAX_CALC, RHOT, RHO, POT, DEXC_LDA, DVXC_LDA, .TRUE.)
!            ENDIF
            
       ELSE
            
!            IF (ISPIN==1) THEN
               CALL RAD_GGA_XC( R, TLDA, RHOT(:,1), RHO(:,1,1), POT(:,1,1), DEXC_GGA, DVXC_GGA)
!            ELSE
!               DO K=1,N
!                  RHOT(K,1)=(RHO(K,1,1)+RHOC(K)+RHO(K,1,2))/(2*SCALE*R%R(K)*R%R(K)) ! up
!                  RHOT(K,1)=MAX(RHOT(K,1), 1E-7_q)
!                  RHOT(K,2)=(RHO(K,1,1)+RHOC(K)-RHO(K,1,2))/(2*SCALE*R%R(K)*R%R(K)) ! down
!                  RHOT(K,2)=MAX(RHOT(K,2), 1E-7_q)
!               ENDDO
!               
!               CALL RAD_GGA_XC_SPIN( R, TLDA, RHOT, RHO(:,1,:), POT(:,1,:), DEXC_GGA, DVXC_GGA)
!            ENDIF
         ENDIF
!      ELSE
!         CALL RAD_GGA_ASPH( R, ISPIN, LMAX_CALC, &
!              RHO, RHOC, POTC, DEXC_LDA, DEXC_GGA, DVXC_LDA, DVXC_GGA, POT )
!      ENDIF
!========================================================================
! thats it
!========================================================================
      DOUBLEC= -DHARTREE/2+DEXC_LDA-DVXC_LDA+DEXC_GGA-DVXC_GGA
      EXCG= DEXC_LDA+DEXC_GGA
!#ifdef debug
!      WRITE(*,1)  N,-DHARTREE/2, -DVXC_LDA-DVXC_GGA, DEXC_LDA+DEXC_GGA
!1     FORMAT(' -Hartree, -vxc, exc:',I4,3F14.7)
!#endif

    END SUBROUTINE RAD_POT


!*******************************************************************
!
!  RAD_POT_WEIGHT
!  because POT will be required only for integration we
!  multiply now with the weights
!
!*******************************************************************

    SUBROUTINE RAD_POT_WEIGHT( R, ISPIN, LMAX, POT)
      IMPLICIT NONE
      REAL(q) :: POT(:,:,:)     ! radial potential
      TYPE (rgrid) :: R
      INTEGER LMAX,ISPIN,I,L,M,LM,K

      DO I=1,ISPIN
      DO L=0,LMAX
      DO M=0,2*L
         LM=L*L+M+1
         DO K=1,R%NMAX
           POT(K,LM,I)=POT(K,LM,I)*R%SI(K)
         ENDDO      
      ENDDO
      ENDDO
      ENDDO
    END SUBROUTINE RAD_POT_WEIGHT

!*******************************************************************
!
! Coloumb potential
! we use the equation
! POT_lm(R) = 4 pi/(2l+1)[ 1/R^(l+1) \int_0^R RHO_lm(r) r^(l) dr
!                       + R^l \int_R^Infinity RHO_lm(r) r^(-l-1) dr ]
!           = 4 pi \int_0^Infinity dr RHO_lm(r) (r,R)<^l / (r,R)>^(l+1)
!
!  rho(r) = \sum_lm  RHO_lm(r) Y_lm(r) / r^2
!    V(r) = \sum_lm  POT_lm(r) Y_lm(r)
!
! which can be obtained by partial integration of
! V  = \int dR 1/R^2  \int  rho(r) dr
!
!*******************************************************************


      SUBROUTINE RAD_POT_HAR(LL,R,POT,RHO,DHARTREE)
      USE prec
      IMPLICIT NONE

      TYPE (rgrid) :: R
      INTEGER LL
      REAL(q) :: RHO(:),POT(:)
      REAL(q) T1(R%NMAX),T2(R%NMAX),V1(R%NMAX),V2(R%NMAX),RL(R%NMAX)
      REAL(q) DHARTREE,H3,EXT
      INTEGER N,I,K,L

      N=R%NMAX

      I=0
    ! additional factor r stems from logarithmic grid
      DO K=N,1,-1
         RL(K)=R%R(K)**LL
         I=I+1
         T2(I)=RHO(K)/RL(K)
         T1(I)=RL(K)*R%R(K)*RHO(K)
      ENDDO
      H3=R%H/ 3.0_q
    ! integrate inward (assuming zero moment for grid point NMAX)
    ! V1 = \int_R^Inf RHO(r) r^l dr
      V1(1)=0
      DO L=3,N,2
        V1(L)  =V1(L-2)+H3*(T1(L-2)+4.0_q*T1(L-1)+T1(L))
        V1(L-1)=V1(L-2)+H3*(1.25_q*T1(L-2)+2.0_q*T1(L-1)-0.25_q*T1(L))
      ENDDO
      IF (MOD(N,2)==0) V1(N)=V1(N-2)+H3*(T1(N-2)+4.0_q*T1(N-1)+T1(N))
    ! V2 = \int_R^Inf RHO(r) r^(-l-1) dr
      V2(1)=0
      DO L=3,N,2
         V2(L)  =V2(L-2)+H3*(T2(L-2)+4.0_q*T2(L-1)+T2(L))
         V2(L-1)=V2(L-2)+H3*(1.25_q*T2(L-2)+2.0_q*T2(L-1)-0.25_q*T2(L))
      ENDDO
      IF (MOD(N,2)==0) V2(N)=V2(N-2)+H3*(T2(N-2)+4.0_q*T2(N-1)+T2(N))

     ! EXT total moment (i.e. charge)
     ! EXT-V1(K) = \int_0^R RHO(r) r^l dr 
      EXT=V1(N)
      I=0
      DO K=N,1,-1
         I=I+1
         POT(I)=(V2(K)*RL(I)+(EXT-V1(K))/(R%R(I)*RL(I)))*(FELECT*4*PI/(2*LL+1))
      ENDDO
    ! double counting corrections are simply obtained by summing over grid
    ! since int Y_lm(Omega_r) Y_lm(Omega_r) d_Omega_r= 1
      DHARTREE=0
      DO K=1,N
         DHARTREE=DHARTREE+POT(K)*(RHO(K)*R%SI(K))
      ENDDO
      END SUBROUTINE




!*******************************************************************
!
! calculate the LDA contribution to the exchange correlation
! energy and to the potentials
! non spin polarised case and below the spinpolarised case
!
! following Bloechl we use a Taylor expansion of the
! exchange correlation energy around the spherical density
!  E_xc = \int d Omega r^2 dr
!'        eps_xc( n_0(r)) + 1/2 v_xc'(n_0) Y_L rho_L(r) Y_L rho_L(r)
! (tick denotes derivative w.r.t. n_0)
! the potential is defined as
!  v_L(r) = var E_xc / var n_L(r) / Y_L
! (since  the factor  Y_L is applied later in the program (see above)
!  we have to divide by Y_L)
!
!
!*******************************************************************

   SUBROUTINE RAD_LDA_XC( R, TREL, LMAX_CALC, RHOT, RHO, POT, DEXC, DVXC, LASPH)
      USE prec
      USE setexm
      IMPLICIT NONE
      TYPE (rgrid) :: R
      LOGICAL TREL      ! relativistic corrections to exchange
      LOGICAL LASPH     ! wether to calculate aspherical corrections
      INTEGER LMAX_CALC ! maximum L
      REAL(q) RHOT(:)   ! spherical charge + core charge density
      REAL(q) RHO(:,:)  ! charge distribution see above
      REAL(q) POT(:,:)  ! potential
      REAL(q) DEXC      ! exchange correlation energy
      REAL(q) DVXC      ! V_xc rho
! local
      REAL(q) EXCA(4)
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DVXC1, RHOPA, SCALE
      INTEGER K,LM

      SCALE=2*SQRT(PI)
      DO K=1,R%NMAX
! MM
!        IF (RHOT(K) <= 0 ) CYCLE
         IF (RHOT(K) <= 1E-99_q ) CYCLE
! MM
         CALL EXCOR_DER_PARA(RHOT(K),LMAX_CALC,EXCA,TREL)

         SIM_FAKT= R%SI(K)*R%R(K)*R%R(K)
         RHOP    = RHO(K,1)/(R%R(K)*R%R(K))

         ! store v_xc(r) / Y_0
         POT(K,1)=POT(K,1)+EXCA(2)*SCALE
         !  \int v_xc(r)  rho_0(r) Y_0 4 Pi r^2 dr
         DVXC    =DVXC    +EXCA(2) *  RHOP*SCALE*SIM_FAKT
         !  \int eps_xc(r) 4 pi r^2 dr
         DEXC    =DEXC    +EXCA(1) *4*PI*SIM_FAKT
         EXT=0
         VXT=0
         IF (LASPH) THEN
         DO LM=2,(LMAX_CALC+1)*(LMAX_CALC+1)
            RHOPA=RHO(K,LM)/ (R%R(K)*R%R(K)) ! rho_L
            ! corrections to energy
            ! 1/2 \int v_xc p(n_0) Y_L rho_L(r) Y_L rho_L(r) d Omega r^2 dr
            !  = 1/2 \int v_xc p(n_0) rho_L(r) rho_L(r) r^2 dr
            EXT =EXT +(EXCA(3)*RHOPA*RHOPA)/2
            ! correction to L=0 potential
            ! 1/2  \int v_xc''(n_0) Y_L rho_L(r) Y_L rho_L(r) d Omega =
            ! \int Y_0 (v_xc''(n_0) rho_L(r) rho_L(r) Y_0/2) d Omega
            DVXC1 = (EXCA(4)*RHOPA*RHOPA)/2/SCALE
            ! \int 1/2 v_xc''(n_0) Y(L) rho_L(r) Y(L) rho_L(r) Y(0) rho_0 d Omega r^2 dr
            ! =\int  Y_0/2 v_xc''(n_0) rho_L(r) rho_L(r)  r^2 dr
            VXT=VXT + DVXC1*RHOP
            POT(K,1)=POT(K,1)+DVXC1

            ! L/=0 potentials
            ! \int Y_L (v_xc p(n_0) \rho_L)
            POT(K,LM)=POT(K,LM)+ (EXCA(3)*RHOPA)
            ! \int v_xc p(n_0) \rho_L Y_L \rho_L Y_L d Omega r^2 dr
            ! =  \int v_xc p(n_0)  \rho_L \rho_L r^2 dr
            VXT=VXT+EXCA(3)*RHOPA*RHOPA
         ENDDO
         ENDIF
         DEXC=DEXC + EXT*SIM_FAKT
         DVXC=DVXC + VXT*SIM_FAKT
      ENDDO

   END SUBROUTINE RAD_LDA_XC


!*******************************************************************
!
! calculate the GGA contribution to the exchange correlation
! energy and to the potentials
! only spherical contributions are accounted for
!
!*******************************************************************

    SUBROUTINE RAD_GGA_XC( R, TLDA, RHOT, RHO, POT, DEXC_GGA, DVXC_GGA)
      USE prec
      IMPLICIT NONE
      TYPE (rgrid) :: R
      LOGICAL TLDA      ! include LDA contributions  (usually .FALSE.)
      REAL(q) RHOT(:)   ! spherical charge + core charge density
      REAL(q) RHO(:)    ! spherical charge
      REAL(q) POT(:)    ! potential
      REAL(q) DEXC_GGA  ! exchange correlation energy
      REAL(q) DVXC_GGA  ! V_xc rho
! local
      REAL(q) SIM_FAKT, RHOP, EXT, VXT, DEXC1, DVXC1, SCALE
      REAL(q) T1(R%NMAX),T2(R%NMAX),V1(R%NMAX)
      INTEGER K

      SCALE=2*SQRT(PI)


      CALL GRAD(R,RHOT,T1)
      DO K=1,R%NMAX
!#ifdef vector
!         CALL GGA91_WB(RHOT(K)*AUTOA3,T1(K)*AUTOA4,EXT,DEXC1,DVXC1)
!#else
         CALL GGAALL(RHOT(K)*AUTOA3,T1(K)*AUTOA4,EXT,DEXC1,DVXC1,.NOT.TLDA)
!#endif
         SIM_FAKT=R%SI(K)*SCALE
         DEXC_GGA=DEXC_GGA+(EXT*RYTOEV)*RHOT(K)*(SCALE*R%R(K)*R%R(K))*SIM_FAKT
         
         !  store d f/ d (d rho )   in T2
         T2(K) = DVXC1*RYTOEV*AUTOA
         !  store d f/ d rho  in T1
         T1(K) = DEXC1*RYTOEV
      ENDDO
      CALL GRAD(R,T2,V1)
      
      DO K=1,R%NMAX
         VXT     = T1(K) - V1(K) - 2*T2(K)/ R%R(K)
         SIM_FAKT=R%SI(K)*SCALE
         DVXC_GGA=DVXC_GGA+VXT*RHO(K)*SIM_FAKT
         POT(K)=POT(K)+VXT*SCALE
      ENDDO

    END SUBROUTINE RAD_GGA_XC

  END MODULE radial
