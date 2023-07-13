      MODULE sed_settling_mod
!
!svn $Id: sed_settling.F 1054 2021-03-06 19:47:12Z arango $
!=======================================================================
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine computes the vertical settling (sinking) of suspended  !
!  sediment via a semi-Lagrangian advective flux algorithm. It uses a  !
!  parabolic,  vertical reconstructuion of the suspended  sediment in  !
!  the water column with PPT/WENO constraints to avoid oscillations.   !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Colella, P. and P. Woodward, 1984: The piecewise parabolic method   !
!    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201. !
!                                                                      !
!  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially        !
!    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.              !
!                                                                      !
!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.    !
!    Arango, 2008:  Development of a three-dimensional,  regional,     !
!    coupled wave, current, and sediment-transport model, Computers    !
!    & Geosciences, 34, 1284-1306.                                     !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: sed_settling
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_settling (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/Sediment/sed_settling.F"
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, iNLM, 16, 66, MyFile)
      CALL sed_settling_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nstp(ng), nnew(ng),                       &
     &                        GRID(ng) % Hz,                            &
     &                        GRID(ng) % z_w,                           &
     &                        SEDBED(ng) % settling_flux,               &
     &                        OCEAN(ng) % t)
      CALL wclock_off (ng, iNLM, 16, 77, MyFile)
!
      RETURN
      END SUBROUTINE sed_settling
!
!***********************************************************************
      SUBROUTINE sed_settling_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS,         &
     &                              nstp, nnew,                         &
     &                              Hz, z_w,                            &
     &                              settling_flux,                      &
     &                              t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(inout) :: settling_flux(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
!
!  Local variable declarations.
!
      integer :: i, indx, ised, j, k, ks
      real(r8) :: cff, cu, cffL, cffR, dltL, dltR
      integer, dimension(IminS:ImaxS,N(ng)) :: ksource
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv2
      real(r8), dimension(IminS:ImaxS,N(ng)) :: Hz_inv3
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qc
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: qL
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WR
      real(r8), dimension(IminS:ImaxS,N(ng)) :: WL
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Add sediment vertical sinking (settling) term.
!-----------------------------------------------------------------------
!
!  Compute inverse thicknesses to avoid repeated divisions.
!
      J_LOOP : DO j=Jstr,Jend
        DO k=1,N(ng)
          DO i=Istr,Iend
            Hz_inv(i,k)=1.0_r8/Hz(i,j,k)
          END DO
        END DO
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            Hz_inv2(i,k)=1.0_r8/(Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            Hz_inv3(i,k)=1.0_r8/(Hz(i,j,k-1)+Hz(i,j,k)+Hz(i,j,k+1))
          END DO
        END DO
!
!  Copy concentration of suspended sediment into scratch array "qc"
!  (q-central, restrict it to be positive) which is hereafter
!  interpreted as a set of grid-box averaged values for sediment
!  concentration.
!
        SED_LOOP: DO ised=1,NST
          indx=idsed(ised)
          DO k=1,N(ng)
            DO i=Istr,Iend
              qc(i,k)=t(i,j,k,nstp,indx)
            END DO
          END DO
!
!-----------------------------------------------------------------------
!  Vertical sinking of suspended sediment.
!-----------------------------------------------------------------------
!
!  Reconstruct vertical profile of suspended sediment "qc" in terms
!  of a set of parabolic segments within each grid box. Then, compute
!  semi-Lagrangian flux due to sinking.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              FC(i,k)=(qc(i,k+1)-qc(i,k))*Hz_inv2(i,k)
            END DO
          END DO
          DO k=2,N(ng)-1
            DO i=Istr,Iend
              dltR=Hz(i,j,k)*FC(i,k)
              dltL=Hz(i,j,k)*FC(i,k-1)
              cff=Hz(i,j,k-1)+2.0_r8*Hz(i,j,k)+Hz(i,j,k+1)
              cffR=cff*FC(i,k)
              cffL=cff*FC(i,k-1)
!
!  Apply PPM monotonicity constraint to prevent oscillations within the
!  grid box.
!
              IF ((dltR*dltL).le.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
!
!  Compute right and left side values (qR,qL) of parabolic segments
!  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
!
!  NOTE: Although each parabolic segment is monotonic within its grid
!        box, monotonicity of the whole profile is not guaranteed,
!        because qL(k+1)-qR(k) may still have different sign than
!        qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!        are reconciled using WENO procedure.
!
              cff=(dltR-dltL)*Hz_inv3(i,k)
              dltR=dltR-cff*Hz(i,j,k+1)
              dltL=dltL+cff*Hz(i,j,k-1)
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
              WR(i,k)=(2.0_r8*dltR-dltL)**2
              WL(i,k)=(dltR-2.0_r8*dltL)**2
            END DO
          END DO
          cff=1.0E-14_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              dltL=MAX(cff,WL(i,k  ))
              dltR=MAX(cff,WR(i,k+1))
              qR(i,k)=(dltR*qR(i,k)+dltL*qL(i,k+1))/(dltR+dltL)
              qL(i,k+1)=qR(i,k)
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8              ! no-flux boundary condition
            qR(i,N(ng))=qc(i,N(ng))         ! default strictly monotonic
            qL(i,N(ng))=qc(i,N(ng))         ! conditions
            qR(i,N(ng)-1)=qc(i,N(ng))
            qL(i,2)=qc(i,1)                 ! bottom grid boxes are
            qR(i,1)=qc(i,1)                 ! re-assumed to be
            qL(i,1)=qc(i,1)                 ! piecewise constant.
          END DO
!
!  Apply monotonicity constraint again, since the reconciled interfacial
!  values may cause a non-monotonic behavior of the parabolic segments
!  inside the grid box.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              dltR=qR(i,k)-qc(i,k)
              dltL=qc(i,k)-qL(i,k)
              cffR=2.0_r8*dltR
              cffL=2.0_r8*dltL
              IF ((dltR*dltL).lt.0.0_r8) THEN
                dltR=0.0_r8
                dltL=0.0_r8
              ELSE IF (ABS(dltR).gt.ABS(cffL)) THEN
                dltR=cffL
              ELSE IF (ABS(dltL).gt.ABS(cffR)) THEN
                dltL=cffR
              END IF
              qR(i,k)=qc(i,k)+dltR
              qL(i,k)=qc(i,k)-dltL
            END DO
          END DO
!
!  After this moment reconstruction is considered complete. The next
!  stage is to compute vertical advective fluxes, FC. It is expected
!  that sinking may occurs relatively fast, the algorithm is designed
!  to be free of CFL criterion, which is achieved by allowing
!  integration bounds for semi-Lagrangian advective flux to use as
!  many grid boxes in upstream direction as necessary.
!
!  In the two code segments below, WL is the z-coordinate of the
!  departure point for grid box interface z_w with the same indices;
!  FC is the finite volume flux; ksource(:,k) is index of vertical
!  grid box which contains the departure point (restricted by N(ng)).
!  During the search: also add in content of whole grid boxes
!  participating in FC.
!
          cff=dt(ng)*ABS(Wsed(ised,ng))
          DO k=1,N(ng)
            DO i=Istr,Iend
              FC(i,k-1)=0.0_r8
              WL(i,k)=z_w(i,j,k-1)+cff
              WR(i,k)=Hz(i,j,k)*qc(i,k)
              ksource(i,k)=k
            END DO
          END DO
          DO k=1,N(ng)
            DO ks=k,N(ng)-1
              DO i=Istr,Iend
                IF (WL(i,k).gt.z_w(i,j,ks)) THEN
                  ksource(i,k)=ks+1
                  FC(i,k-1)=FC(i,k-1)+WR(i,ks)
                END IF
              END DO
            END DO
          END DO
!
!  Finalize computation of flux: add fractional part.
!
          DO k=1,N(ng)
            DO i=Istr,Iend
              ks=ksource(i,k)
              cu=MIN(1.0_r8,(WL(i,k)-z_w(i,j,ks-1))*Hz_inv(i,ks))
              FC(i,k-1)=FC(i,k-1)+                                      &
     &                  Hz(i,j,ks)*cu*                                  &
     &                  (qL(i,ks)+                                      &
     &                   cu*(0.5_r8*(qR(i,ks)-qL(i,ks))-                &
     &                       (1.5_r8-cu)*                               &
     &                       (qR(i,ks)+qL(i,ks)-2.0_r8*qc(i,ks))))
            END DO
          END DO
          DO i=Istr,Iend
            DO k=1,N(ng)
              t(i,j,k,nnew,indx)=t(i,j,k,nnew,indx)+(FC(i,k)-FC(i,k-1))
            END DO
            settling_flux(i,j,ised)=FC(i,0)
          END DO
        END DO SED_LOOP
      END DO J_LOOP
!
      RETURN
      END SUBROUTINE sed_settling_tile
      END MODULE sed_settling_mod
