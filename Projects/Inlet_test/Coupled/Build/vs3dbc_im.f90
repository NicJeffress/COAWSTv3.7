      MODULE vs3dbc_mod
!
!svn $Id: vs3dbc_im.F 732 2008-09-07 01:55:51Z jcwarner $
!=======================================================================
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine sets lateral boundary conditions for total 3D       !
!  Vstokes-velocity.                                                   !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: vs3dbc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE vs3dbc (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
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
      CALL vs3dbc_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj, N(ng),                      &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  OCEAN(ng) % v_stokes)
      RETURN
      END SUBROUTINE vs3dbc
!
!***********************************************************************
      SUBROUTINE vs3dbc_tile (ng, tile,                                 &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       v_stokes)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_boundary
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(inout) :: v_stokes(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k, Jmin, Jmax
      real(r8), parameter :: eps = 1.0E-20_r8
      real(r8) :: Ce, Cx, cff, dVde, dVdt, dVdx, tau
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
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
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=v_stokes(i  ,Jstr  ,k)-                    &
     &                       v_stokes(i-1,Jstr  ,k)
              grad(i,Jstr+1)=v_stokes(i  ,Jstr+1,k)-                    &
     &                       v_stokes(i-1,Jstr+1,k)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dVdt=v_stokes(i,Jstr+1,k)-v_stokes(i,Jstr+1,k)
                dVde=v_stokes(i,Jstr+1,k)-v_stokes(i,Jstr+2,k)
                IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(i,Jstr+1)+grad(i+1,Jstr+1))).gt.        &
     &              0.0_r8) THEN
                  dVdx=grad(i  ,Jstr+1)
                ELSE
                  dVdx=grad(i+1,Jstr+1)
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
                Cx=0.0_r8
                Ce=dVdt*dVde
                v_stokes(i,Jstr,k)=(cff*v_stokes(i,Jstr  ,k)+           &
     &                             Ce *v_stokes(i,Jstr+1,k)-            &
     &                             MAX(Cx,0.0_r8)*grad(i  ,Jstr)-       &
     &                             MIN(Cx,0.0_r8)*grad(i+1,Jstr))/      &
     &                             (cff+Ce)
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask(i,Jstr)
              END IF
            END DO
          END DO
!
!  Southern edge, clamped boundary condition.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v_stokes(i,Jstr,k)=BOUNDARY(ng)%vstokes_south(i,k)
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
     &                             GRID(ng)%vmask(i,Jstr)
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr+1,k)
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
       &                           GRID(ng)%vmask(i,Jstr)
              END IF
            END DO
          END DO
!
!  Southern edge, nested.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%nested) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr+1,k)
                v_stokes(i,Jstr,k)=v_stokes(i,Jstr,k)*                  &
       &                           GRID(ng)%vmask(i,Jstr)
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isV3Sd,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                v_stokes(i,Jstr,k)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=v_stokes(i  ,Jend  ,k)-                    &
     &                       v_stokes(i-1,Jend  ,k)
              grad(i,Jend+1)=v_stokes(i  ,Jend+1,k)-                    &
     &                       v_stokes(i-1,Jend+1,k)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dVdt=v_stokes(i,Jend,k)-v_stokes(i,Jend  ,k)
                dVde=v_stokes(i,Jend,k)-v_stokes(i,Jend-1,k)
                IF ((dVdt*dVde).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(i,Jend)+grad(i+1,Jend))).gt.0.0_r8) THEN
                  dVdx=grad(i  ,Jend)
                ELSE
                  dVdx=grad(i+1,Jend)
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
                Cx=0.0_r8
                Ce=dVdt*dVde
                v_stokes(i,Jend+1,k)=(cff*v_stokes(i,Jend+1,k)+         &
     &                               Ce *v_stokes(i,Jend  ,k)-          &
     &                               MAX(Cx,0.0_r8)*grad(i  ,Jend+1)-   &
     &                               MIN(Cx,0.0_r8)*grad(i+1,Jend+1))/  &
     &                               (cff+Ce)
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask(i,Jend+1)
              END IF
            END DO
          END DO
!
!  Northern edge, clamped boundary condition.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v_stokes(i,Jend+1,k)=BOUNDARY(ng)%v_north(i,k)
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                               GRID(ng)%vmask(i,Jend+1)
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend,k)
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                                GRID(ng)%vmask(i,Jend+1)
            END IF
          END DO
        END DO
!
!  Northern edge, nested.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%nested) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend,k)
                v_stokes(i,Jend+1,k)=v_stokes(i,Jend+1,k)*              &
     &                                GRID(ng)%vmask(i,Jend+1)
          END DO
        END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isV3Sd,ng)%closed) THEN
          DO k=1,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                v_stokes(i,Jend+1,k)=0.0_r8
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Istr-1,j)=v_stokes(Istr-1,j+1,k)-                    &
     &                       v_stokes(Istr-1,j  ,k)
              grad(Istr  ,j)=v_stokes(Istr  ,j+1,k)-                    &
     &                       v_stokes(Istr  ,j  ,k)
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dVdt=v_stokes(Istr,j,k)-v_stokes(Istr  ,j,k)
                dVdx=v_stokes(Istr,j,k)-v_stokes(Istr+1,j,k)
                IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(Istr,j-1)+grad(Istr,j))).gt.0.0_r8) THEN
                  dVde=grad(Istr,j-1)
                ELSE
                  dVde=grad(Istr,j  )
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
                Cx=dVdt*dVdx
                Ce=0.0_r8
                v_stokes(Istr-1,j,k)=(cff*v_stokes(Istr-1,j,k)+         &
     &                                Cx *v_stokes(Istr  ,j,k)-         &
     &                                MAX(Ce,0.0_r8)*grad(Istr-1,j-1)-  &
     &                                MIN(Ce,0.0_r8)*grad(Istr-1,j  ))/ &
     &                                (cff+Cx)
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
              END IF
            END DO
          END DO
!
!  Western edge, clamped boundary condition.
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                v_stokes(Istr-1,j,k)=BOUNDARY(ng)%vstokes_west(j,k)
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                v_stokes(Istr-1,j,k)=v_stokes(Istr,j,k)
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
              END IF
            END DO
          END DO
!
!  Western edge, nested.
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%nested) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
                v_stokes(Istr-1,j,k)=v_stokes(Istr,j,k)
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
            END DO
          END DO
!
!  Western edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(iwest,isV3Sd,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%west(j)) THEN
                v_stokes(Istr-1,j,k)=gamma2(ng)*v_stokes(Istr,j,k)
                v_stokes(Istr-1,j,k)=v_stokes(Istr-1,j,k)*              &
     &                               GRID(ng)%vmask(Istr-1,j)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isV3Sd,ng)%radiation) THEN
          DO k=1,N(ng)
            DO j=JstrV-1,Jend
              grad(Iend  ,j)=v_stokes(Iend  ,j+1,k)-                    &
     &                       v_stokes(Iend  ,j  ,k)
              grad(Iend+1,j)=v_stokes(Iend+1,j+1,k)-                    &
     &                       v_stokes(Iend+1,j  ,k)
            END DO
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dVdt=v_stokes(Iend,j,k)-v_stokes(Iend  ,j,k)
                dVdx=v_stokes(Iend,j,k)-v_stokes(Iend-1,j,k)
                IF ((dVdt*dVdx).lt.0.0_r8) dVdt=0.0_r8
                IF ((dVdt*(grad(Iend,j-1)+grad(Iend,j))).gt.0.0_r8) THEN
                  dVde=grad(Iend,j-1)
                ELSE
                  dVde=grad(Iend,j  )
                END IF
                cff=MAX(dVdx*dVdx+dVde*dVde,eps)
                Cx=dVdt*dVdx
                Ce=0.0_r8
                v_stokes(Iend+1,j,k)=(cff*v_stokes(Iend+1,j,k)+         &
     &                               Cx *v_stokes(Iend  ,j,k)-          &
     &                               MAX(Ce,0.0_r8)*grad(Iend+1,j-1)-   &
     &                               MIN(Ce,0.0_r8)*grad(Iend+1,j  ))/  &
     &                               (cff+Cx)
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
              END IF
            END DO
          END DO
!
!  Eastern edge, clamped boundary condition.
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%clamped) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                v_stokes(Iend+1,j,k)=BOUNDARY(ng)%vstokes_east(j,k)
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%gradient) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                v_stokes(Iend+1,j,k)=v_stokes(Iend,j,k)
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
              END IF
            END DO
          END DO
!
!  Eastern edge, nested.
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%nested) THEN
          DO k=1,N(ng)
            DO j=JstrV,Jend
                v_stokes(Iend+1,j,k)=v_stokes(Iend,j,k)
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
            END DO
          END DO
!
!  Eastern edge, closed boundary condition: free slip (gamma2=1)  or
!                                           no   slip (gamma2=-1).
!
        ELSE IF (LBC(ieast,isV3Sd,ng)%closed) THEN
          IF (NSperiodic(ng)) THEN
            Jmin=JstrV
            Jmax=Jend
          ELSE
            Jmin=Jstr
            Jmax=JendR
          END IF
          DO k=1,N(ng)
            DO j=Jmin,Jmax
              IF (LBC_apply(ng)%east(j)) THEN
                v_stokes(Iend+1,j,k)=gamma2(ng)*v_stokes(Iend,j,k)
                v_stokes(Iend+1,j,k)=v_stokes(Iend+1,j,k)*              &
     &                               GRID(ng)%vmask(Iend+1,j)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF ((LBC_apply(ng)%south(Istr-1).and.                         &
     &        LBC_apply(ng)%west (Jstr  )).or.                          &
     &        (LBC(iwest,isV3Sd,ng)%nested.and.                         &
     &         LBC(isouth,isV3Sd,ng)%nested)) THEN
            DO k=1,N(ng)
              v_stokes(Istr-1,Jstr,k)=0.5_r8*                           &
     &                                (v_stokes(Istr  ,Jstr  ,k)+       &
     &                                 v_stokes(Istr-1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF ((LBC_apply(ng)%south(Iend+1).and.                         &
     &        LBC_apply(ng)%east (Jstr  )).or.                          &
     &        (LBC(ieast,isV3Sd,ng)%nested.and.                         &
     &         LBC(isouth,isV3Sd,ng)%nested)) THEN
            DO k=1,N(ng)
              v_stokes(Iend+1,Jstr,k)=0.5_r8*                           &
     &                                (v_stokes(Iend  ,Jstr  ,k)+       &
     &                                 v_stokes(Iend+1,Jstr+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF ((LBC_apply(ng)%north(Istr-1).and.                         &
     &        LBC_apply(ng)%west (Jend+1)).or.                          &
     &        (LBC(iwest,isV3Sd,ng)%nested.and.                         &
     &         LBC(inorth,isV3Sd,ng)%nested)) THEN
            DO k=1,N(ng)
              v_stokes(Istr-1,Jend+1,k)=0.5_r8*                         &
     &                                  (v_stokes(Istr-1,Jend  ,k)+     &
     &                                   v_stokes(Istr  ,Jend+1,k))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF ((LBC_apply(ng)%north(Iend+1).and.                         &
     &        LBC_apply(ng)%east (Jend+1)).or.                          &
     &        (LBC(ieast,isV3Sd,ng)%nested.and.                         &
     &         LBC(inorth,isV3Sd,ng)%nested)) THEN
            DO k=1,N(ng)
              v_stokes(Iend+1,Jend+1,k)=0.5_r8*                         &
     &                                  (v_stokes(Iend+1,Jend  ,k)+     &
     &                                   v_stokes(Iend  ,Jend+1,k))
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE vs3dbc_tile
      END MODULE vs3dbc_mod
