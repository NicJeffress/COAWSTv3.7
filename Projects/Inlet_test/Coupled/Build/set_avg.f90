      MODULE set_avg_mod
!
!svn $Id: set_avg.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine accumulates and computes output time-averaged       !
!  fields.  Due to synchronization, the time-averaged fields are       !
!  computed in delayed mode. All averages are accumulated at the       !
!  beggining of the next time-step.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: set_avg
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_avg (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/set_avg.F"
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
      CALL wclock_on (ng, iNLM, 5, 62, MyFile)
      CALL set_avg_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nrhs(ng),                                          &
     &                   kstp(ng))
      CALL wclock_off (ng, iNLM, 5, 102, MyFile)
!
      RETURN
      END SUBROUTINE set_avg
!
!***********************************************************************
      SUBROUTINE set_avg_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         Nout,                                    &
     &                         Kout)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_average
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_scalars
      USE mod_bbl
      USE mod_sedbed
      USE mod_sediment
!
      USE exchange_2d_mod
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE uv_rotate_mod, ONLY : uv_rotate2d
      USE uv_rotate_mod, ONLY : uv_rotate3d
      USE vorticity_mod, ONLY : vorticity_tile
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Kout
      integer, intent(in) :: Nout
!
!
!  Local variable declarations.
!
      integer :: i, it, itrc, j, k
      real(r8) :: fac
      real(r8) :: pfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: rfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: ufac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: vfac(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: potvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: relvor(LBi:UBi,LBj:UBj,N(ng))
      real(r8) :: potvor_bar(LBi:UBi,LBj:UBj)
      real(r8) :: relvor_bar(LBi:UBi,LBj:UBj)
      real(r8), allocatable :: wrk(:,:)
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
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  Compute vorticity fields.
!-----------------------------------------------------------------------
!
      IF (Aout(id2dPV,ng).or.Aout(id2dRV,ng).or.                        &
     &    Aout(id3dPV,ng).or.Aout(id3dRV,ng)) THEN
        CALL vorticity_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Kout, Nout,                                &
     &                       GRID(ng) % pmask,                          &
     &                       GRID(ng) % umask,                          &
     &                       GRID(ng) % vmask,                          &
     &                       GRID(ng) % fomn,                           &
     &                       GRID(ng) % h,                              &
     &                       GRID(ng) % om_u,                           &
     &                       GRID(ng) % on_v,                           &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % z_r,                            &
     &                       OCEAN(ng) % pden,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % ubar,                          &
     &                       OCEAN(ng) % vbar,                          &
     &                       OCEAN(ng) % zeta,                          &
     &                       potvor, relvor,                            &
     &                       potvor_bar, relvor_bar)
      END IF
!
!-----------------------------------------------------------------------
!  Initialize time-averaged arrays when appropriate.  Notice that
!  fields are initilized twice during re-start.  However, the time-
!  averaged fields are computed correctly.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.1)).or.                          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1)).or.             &
     &    ((nrrec(ng).gt.0).and.(iic(ng).eq.ntstart(ng)))) THEN
!
!  Initialize state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_full,                        &
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .FALSE., .FALSE.,                 &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_full,                        &
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%bedldu(i,j,it)
              END DO
            END DO
          END IF
          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%bedldv(i,j,it)
              END DO
            END DO
          END IF
        END DO
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
!
!  Initialize surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbrs(i,j)=BBL(ng)%bustrc(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbrs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbrs(i,j)=BBL(ng)%bvstrc(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbws(i,j)=BBL(ng)%bustrw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbws,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbws(i,j)=BBL(ng)%bvstrw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbcs(i,j)=BBL(ng)%bustrcwmax(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbcs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbcs(i,j)=BBL(ng)%bvstrcwmax(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUVwc,ng)) THEN
          allocate (wrk(LBi:UBi,LBj:UBj))
          wrk(LBi:UBi,LBj:UBj)=0.0_r8
          wrk=sqrt(BBL(ng)%bustrcwmax*BBL(ng)%bustrcwmax+               &
     &             BBL(ng)%bvstrcwmax*BBL(ng)%bvstrcwmax+1.0E-10_r8)
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUVwc(i,j)=wrk(i,j)
            END DO
          END DO
          deallocate(wrk)
        END IF
        IF (Aout(idUbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbot(i,j)=BBL(ng)%Ubot(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbot,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbot(i,j)=BBL(ng)%Vbot(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbur,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgUbur(i,j)=BBL(ng)%Ur(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbvr,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgVbvr(i,j)=BBL(ng)%Vr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
!
!  Initialize  values.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=OCEAN(ng)%ubar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=OCEAN(ng)%vbar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2rs(i,j)=MIXING(ng)%rustr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2rs(i,j)=MIXING(ng)%rvstr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=OCEAN(ng)%u_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=OCEAN(ng)%v_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3Sd(i,j,k)=OCEAN(ng)%W_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3St,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3St(i,j,k)=OCEAN(ng)%wstvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3rs(i,j,k)=MIXING(ng)%rustr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3rs(i,j,k)=MIXING(ng)%rvstr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWztw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWztw(i,j)=OCEAN(ng)%zetaw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWqsp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgwqsp(i,j)=OCEAN(ng)%qsp(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWbeh,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgwbeh(i,j)=OCEAN(ng)%bh(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWamp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWamp(i,j)=FORCES(ng)%Hwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWam2,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWam2(i,j)=FORCES(ng)%Hwave(i,j)*           &
     &                                 FORCES(ng)%Hwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWlen,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlen(i,j)=FORCES(ng)%Lwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdir,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdir(i,j)=FORCES(ng)%Dwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdip,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdip(i,j)=FORCES(ng)%Dwavep(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWptp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWptp(i,j)=FORCES(ng)%Pwave_top(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWpbt,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWpbt(i,j)=FORCES(ng)%Pwave_bot(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWorb,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWorb(i,j)=FORCES(ng)%Uwave_rms(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdif,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdif(i,j)=FORCES(ng)%Dissip_fric(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdib,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdib(i,j)=FORCES(ng)%Dissip_break(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiw(i,j)=FORCES(ng)%Dissip_wcap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwav(i,j)=OCEAN(ng)%uWave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwav(i,j)=OCEAN(ng)%vWave(i,j)
            END DO
          END DO
        END IF
!
!  Initialize vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=relvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
!
!  Initialize quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Accumulate time-averaged fields.
!-----------------------------------------------------------------------
!
      ELSE IF (iic(ng).gt.ntsAVG(ng)) THEN
!
!  Accumulate state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=AVERAGE(ng)%avgzeta(i,j)+        &
     &                                 OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=AVERAGE(ng)%avgu2d(i,j)+          &
     &                                OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=AVERAGE(ng)%avgv2d(i,j)+          &
     &                                OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          CALL uv_rotate2d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_full,                        &
     &                      OCEAN(ng) % ubar(:,:,Kout),                 &
     &                      OCEAN(ng) % vbar(:,:,Kout),                 &
     &                      AVERAGE(ng)%avgu2dE,                        &
     &                      AVERAGE(ng)%avgv2dN)
        END IF
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=AVERAGE(ng)%avgu3d(i,j,k)+    &
     &                                    OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=AVERAGE(ng)%avgv3d(i,j,k)+    &
     &                                    OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          CALL uv_rotate3d (ng, tile, .TRUE., .FALSE.,                  &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      GRID(ng) % CosAngler,                       &
     &                      GRID(ng) % SinAngler,                       &
     &                      GRID(ng)%rmask_full,                        &
     &                      OCEAN(ng) % u(:,:,:,Nout),                  &
     &                      OCEAN(ng) % v(:,:,:,Nout),                  &
     &                      AVERAGE(ng)%avgu3dE,                        &
     &                      AVERAGE(ng)%avgv3dN)
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=AVERAGE(ng)%avgw3d(i,j,k)+    &
     &                                    OCEAN(ng)%W(i,j,k)*           &
     &                                    GRID(ng)%pm(i,j)*             &
     &                                    GRID(ng)%pn(i,j)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=AVERAGE(ng)%avgwvel(i,j,k)+  &
     &                                     OCEAN(ng)%wvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=AVERAGE(ng)%avgrho(i,j,k)+    &
     &                                    OCEAN(ng)%rho(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=AVERAGE(ng)%avgt(i,j,k,it)+&
     &                                       OCEAN(ng)%t(i,j,k,Nout,it)
                END DO
              END DO
            END DO
          END IF
        END DO
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=SEDBED(ng)%avgbedldu(i,j,  &
     &                                                            it)+  &
     &                                       SEDBED(ng)%bedldu(i,j,it)
              END DO
            END DO
          END IF
          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=SEDBED(ng)%avgbedldv(i,j,  &
     &                                                            it)+  &
     &                                       SEDBED(ng)%bedldv(i,j,it)
              END DO
            END DO
          END IF
        END DO
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=AVERAGE(ng)%avgAKv(i,j,k)+    &
     &                                    MIXING(ng)%Akv(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=AVERAGE(ng)%avgAKt(i,j,k)+    &
     &                                    MIXING(ng)%Akt(i,j,k,itemp)
              END DO
            END DO
          END DO
        END IF
!
!  Accumulate surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=AVERAGE(ng)%avgsus(i,j)+          &
     &                                FORCES(ng)%sustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=AVERAGE(ng)%avgsvs(i,j)+          &
     &                                FORCES(ng)%svstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=AVERAGE(ng)%avgbus(i,j)+          &
     &                                FORCES(ng)%bustr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=AVERAGE(ng)%avgbvs(i,j)+          &
     &                                FORCES(ng)%bvstr(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=AVERAGE(ng)%avgstf(i,j)+          &
     &                                FORCES(ng)%stflx(i,j,itemp)
            END DO
          END DO
        END IF
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=AVERAGE(ng)%avgu2Sd(i,j)+        &
     &                                 OCEAN(ng)%ubar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=AVERAGE(ng)%avgv2Sd(i,j)+        &
     &                                 OCEAN(ng)%vbar_stokes(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2rs(i,j)=AVERAGE(ng)%avgu2rs(i,j)+        &
     &                                 MIXING(ng)%rustr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2rs(i,j)=AVERAGE(ng)%avgv2rs(i,j)+        &
     &                                 MIXING(ng)%rvstr2d(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=AVERAGE(ng)%avgu3Sd(i,j,k)+  &
     &                                     OCEAN(ng)%u_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=AVERAGE(ng)%avgv3Sd(i,j,k)+  &
     &                                     OCEAN(ng)%v_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3Sd(i,j,k)=AVERAGE(ng)%avgw3Sd(i,j,k)+  &
     &                                     OCEAN(ng)%W_stokes(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idW3St,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3St(i,j,k)=AVERAGE(ng)%avgw3St(i,j,k)+  &
     &                                     OCEAN(ng)%wstvel(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3rs(i,j,k)=AVERAGE(ng)%avgu3rs(i,j,k)+  &
     &                                     MIXING(ng)%rustr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3rs(i,j,k)=AVERAGE(ng)%avgv3rs(i,j,k)+  &
     &                                     MIXING(ng)%rvstr3d(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idWztw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWztw(i,j)=AVERAGE(ng)%avgWztw(i,j)+        &
     &                                 OCEAN(ng)%zetaw(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWqsp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWqsp(i,j)=AVERAGE(ng)%avgWqsp(i,j)+        &
     &                                 OCEAN(ng)%qsp(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWbeh,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbeh(i,j)=AVERAGE(ng)%avgWbeh(i,j)+        &
     &                                 OCEAN(ng)%bh(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWamp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWamp(i,j)=AVERAGE(ng)%avgWamp(i,j)+        &
     &                                 FORCES(ng)%Hwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWam2,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWam2(i,j)=AVERAGE(ng)%avgWam2(i,j)+        &
     &                                 FORCES(ng)%Hwave(i,j)*           &
     &                                 FORCES(ng)%Hwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWlen,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlen(i,j)=AVERAGE(ng)%avgWlen(i,j)+        &
     &                                 FORCES(ng)%Lwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdir,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdir(i,j)=AVERAGE(ng)%avgWdir(i,j)+        &
     &                                 FORCES(ng)%Dwave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdip,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdip(i,j)=AVERAGE(ng)%avgWdip(i,j)+        &
     &                                 FORCES(ng)%Dwavep(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWptp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWptp(i,j)=AVERAGE(ng)%avgWptp(i,j)+        &
     &                                 FORCES(ng)%Pwave_top(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWpbt,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWpbt(i,j)=AVERAGE(ng)%avgWpbt(i,j)+        &
     &                                 FORCES(ng)%Pwave_bot(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWorb,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWorb(i,j)=AVERAGE(ng)%avgWorb(i,j)+        &
     &                                 FORCES(ng)%Uwave_rms(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdif,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdif(i,j)=AVERAGE(ng)%avgWdif(i,j)+        &
     &                                 FORCES(ng)%Dissip_fric(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdib,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdib(i,j)=AVERAGE(ng)%avgWdib(i,j)+        &
     &                                 FORCES(ng)%Dissip_break(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idWdiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiw(i,j)=AVERAGE(ng)%avgWdiw(i,j)+        &
     &                                 FORCES(ng)%Dissip_wcap(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idUwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwav(i,j)=AVERAGE(ng)%avgUwav(i,j)+        &
     &                                 OCEAN(ng)%uWave(i,j)
            END DO
          END DO
        END IF
        IF (Aout(idVwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwav(i,j)=AVERAGE(ng)%avgVwav(i,j)+        &
     &                                 OCEAN(ng)%vWave(i,j)
            END DO
          END DO
        END IF
!
!  Accumulate vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=AVERAGE(ng)%avgpvor2d(i,j)+    &
     &                                   potvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=AVERAGE(ng)%avgrvor2d(i,j)+    &
     &                                   relvor_bar(i,j)
            END DO
          END DO
        END IF
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=AVERAGE(ng)%avgpvor3d(i,j, &
     &                                                             k)+  &
     &                                       potvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=AVERAGE(ng)%avgrvor3d(i,j, &
     &                                                             k)+  &
     &                                       relvor(i,j,k)
              END DO
            END DO
          END DO
        END IF
!
!  Accumulate quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=AVERAGE(ng)%avgZZ(i,j)+            &
     &                               OCEAN(ng)%zeta(i,j,Kout)*          &
     &                               OCEAN(ng)%zeta(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=AVERAGE(ng)%avgU2(i,j)+            &
     &                               OCEAN(ng)%ubar(i,j,Kout)*          &
     &                               OCEAN(ng)%ubar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=AVERAGE(ng)%avgV2(i,j)+            &
     &                               OCEAN(ng)%vbar(i,j,Kout)*          &
     &                               OCEAN(ng)%vbar(i,j,Kout)
            END DO
          END DO
        END IF
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=AVERAGE(ng)%avgUU(i,j,k)+      &
     &                                   OCEAN(ng)%u(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%u(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=AVERAGE(ng)%avgVV(i,j,k)+      &
     &                                   OCEAN(ng)%v(i,j,k,Nout)*       &
     &                                   OCEAN(ng)%v(i,j,k,Nout)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=AVERAGE(ng)%avgUV(i,j,k)+      &
     &                                   0.25_r8*                       &
     &                                   (OCEAN(ng)%u(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%u(i+1,j  ,k,Nout))* &
     &                                   (OCEAN(ng)%v(i  ,j  ,k,Nout)+  &
     &                                    OCEAN(ng)%v(i  ,j+1,k,Nout))
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=AVERAGE(ng)%avgHuon(i,j,k)+  &
     &                                     GRID(ng)%Huon(i,j,k)
              END DO
            END DO
          END DO
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=AVERAGE(ng)%avgHvom(i,j,k)+  &
     &                                     GRID(ng)%Hvom(i,j,k)
              END DO
            END DO
          END DO
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=AVERAGE(ng)%avgTT(i,j,k,  &
     &                                                          it)+    &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)*     &
     &                                        OCEAN(ng)%t(i,j,k,        &
     &                                                    Nout,it)
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=AVERAGE(ng)%avgUT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%u(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i-1,j,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i  ,j,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=AVERAGE(ng)%avgVT(i,j,k,  &
     &                                                          it)+    &
     &                                        0.5_r8*                   &
     &                                        OCEAN(ng)%v(i,j,k,Nout)*  &
     &                                        (OCEAN(ng)%t(i,j-1,k,     &
     &                                                     Nout,it)+    &
     &                                         OCEAN(ng)%t(i,j  ,k,     &
     &                                                     Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=AVERAGE(ng)%avgHuonT(i,&
     &                                                       j,k,it)+   &
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Huon(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i-1,j,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i  ,j,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=AVERAGE(ng)%avgHvomT(i,&
     &                                                       j,k,it)+   &
     &                                           0.5_r8*                &
     &                                           GRID(ng)%Hvom(i,j,k)*  &
     &                                           (OCEAN(ng)%t(i,j-1,k,  &
     &                                                        Nout,it)+ &
     &                                            OCEAN(ng)%t(i,j  ,k,  &
     &                                                        Nout,it))
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Convert accumulated sums into time-averages, if appropriate.
!  Notice that we need to apply periodic conditions, if any, since
!  the full I- and J-ranges are different.
!-----------------------------------------------------------------------
!
      IF (((iic(ng).gt.ntsAVG(ng)).and.                                 &
     &     (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                          &
     &     ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))).or.          &
     &    ((iic(ng).ge.ntsAVG(ng)).and.(nAVG(ng).eq.1))) THEN
        IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
          IF (nAVG(ng).eq.1) THEN
            AVGtime(ng)=time(ng)
          ELSE
            AVGtime(ng)=AVGtime(ng)+REAL(nAVG(ng),r8)*dt(ng)
          END IF
        END IF
!
!  Set time-averaged factors for each C-grid variable type. Notice that
!  the I- and J-ranges are all grid types are the same for convinience.
!
        fac=1.0_r8/REAL(nAVG(ng),r8)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            pfac(i,j)=fac
            rfac(i,j)=fac
            ufac(i,j)=fac
            vfac(i,j)=fac
          END DO
        END DO
!
!  Process state variables.
!
        IF (Aout(idFsur,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgzeta(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgzeta(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgzeta)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgzeta)
          END IF
        END IF
        IF (Aout(idUbar,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2d(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgu2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2d)
          END IF
        END IF
        IF (Aout(idVbar,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2d(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgv2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2d)
          END IF
        END IF
        IF (Aout(idu2dE,ng).and.Aout(idv2dN,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgu2dE(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2dE(i,j)
              AVERAGE(ng)%avgv2dN(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2dN(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2dE)
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2dN)
            CALL mp_exchange2d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2dE,                    &
     &                          AVERAGE(ng)%avgv2dN)
          END IF
        END IF
        IF (Aout(idUvel,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3d(i,j,k)=ufac(i,j)*                    &
     &                                    AVERAGE(ng)%avgu3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3d)
          END IF
        END IF
        IF (Aout(idVvel,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3d(i,j,k)=vfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgv3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3d)
          END IF
        END IF
        IF (Aout(idu3dE,ng).and.Aout(idv3dN,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgu3dE(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3dE(i,j,k)
                AVERAGE(ng)%avgv3dN(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3dN(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3dE)
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3dN)
            CALL mp_exchange3d (ng, tile, iNLM, 2,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3dE,                    &
     &                          AVERAGE(ng)%avgv3dN)
          END IF
        END IF
        IF (Aout(idOvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgw3d(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgw3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgw3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgw3d)
          END IF
        END IF
        IF (Aout(idWvel,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgwvel(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgwvel(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgwvel)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgwvel)
          END IF
        END IF
        IF (Aout(idDano,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgrho(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgrho(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrho)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrho)
          END IF
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTvar(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgt(i,j,k,it)=rfac(i,j)*                 &
     &                                       AVERAGE(ng)%avgt(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgt(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgt(:,:,:,it))
            END IF
          END IF
        END DO
        DO it=1,NST
          IF (Aout(idUbld(it),ng)) THEN
            DO j=JstrR,JendR
              DO i=Istr,IendR
                SEDBED(ng)%avgbedldu(i,j,it)=ufac(i,j)*                 &
     &                                       SEDBED(ng)%avgbedldu(i,j,  &
     &                                                            it)
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                SEDBED(ng)%avgbedldu(:,:,it))
              CALL mp_exchange2d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            SEDBED(ng)%avgbedldu(:,:,it))
            END IF
          END IF
          IF (Aout(idVbld(it),ng)) THEN
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                SEDBED(ng)%avgbedldv(i,j,it)=vfac(i,j)*                 &
     &                                       SEDBED(ng)%avgbedldv(i,j,  &
     &                                                            it)
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v2d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj,               &
     &                                SEDBED(ng)%avgbedldv(:,:,it))
              CALL mp_exchange2d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            SEDBED(ng)%avgbedldv(:,:,it))
            END IF
          END IF
        END DO
        IF (Aout(idVvis,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKv(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKv(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKv)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKv)
          END IF
        END IF
        IF (Aout(idTdif,ng)) THEN
          DO k=0,N(ng)
            DO j=JstrR,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgAKt(i,j,k)=rfac(i,j)*                    &
     &                                    AVERAGE(ng)%avgAKt(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_w3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 0, N(ng),       &
     &                              AVERAGE(ng)%avgAKt)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgAKt)
          END IF
        END IF
!
!  Process surface and bottom fluxes.
!
        IF (Aout(idUsms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgsus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgsus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsus)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsus)
          END IF
        END IF
        IF (Aout(idVsms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgsvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgsvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgsvs)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgsvs)
          END IF
        END IF
        IF (Aout(idUbms,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgbus(i,j)=ufac(i,j)*                        &
     &                                AVERAGE(ng)%avgbus(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbus)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbus)
          END IF
        END IF
        IF (Aout(idVbms,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgbvs(i,j)=vfac(i,j)*                        &
     &                                AVERAGE(ng)%avgbvs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgbvs)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgbvs)
          END IF
        END IF
        IF (Aout(idTsur(itemp),ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgstf(i,j)=rfac(i,j)*                        &
     &                                AVERAGE(ng)%avgstf(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgstf)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgstf)
          END IF
        END IF
!
!  Process radiation stresses.
!
        IF (Aout(idU2Sd,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2Sd(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2Sd(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2Sd)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2Sd)
          END IF
        END IF
        IF (Aout(idV2Sd,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2Sd(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2Sd(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2Sd)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2Sd)
          END IF
        END IF
        IF (Aout(idU2rs,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgu2rs(i,j)=ufac(i,j)*                       &
     &                                 AVERAGE(ng)%avgu2rs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgu2rs)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu2rs)
          END IF
        END IF
        IF (Aout(idV2rs,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgv2rs(i,j)=vfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgv2rs(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgv2rs)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv2rs)
          END IF
        END IF
        IF (Aout(idU3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3Sd(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3Sd)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3Sd)
          END IF
        END IF
        IF (Aout(idV3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3Sd(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3Sd)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3Sd)
          END IF
        END IF
        IF (Aout(idU3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgu3rs(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgu3rs(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgu3rs)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgu3rs)
          END IF
        END IF
        IF (Aout(idV3rs,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgv3RS(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgv3RS(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgv3RS)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgv3RS)
          END IF
        END IF
        IF (Aout(idW3Sd,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgW3Sd(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgW3Sd(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgW3Sd)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgW3Sd)
          END IF
        END IF
        IF (Aout(idW3St,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgW3St(i,j,k)=rfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgW3St(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgW3St)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgW3St)
          END IF
        END IF
        IF (Aout(idWztw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWztw(i,j)=rfac(i,j)*                       &
     &                                 AVERAGE(ng)%avgWztw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWztw)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWztw)
          END IF
        END IF
        IF (Aout(idWqsp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWqsp(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWqsp(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWqsp)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWqsp)
          END IF
        END IF
        IF (Aout(idWbeh,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWbeh(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWbeh(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWbeh)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWbeh)
          END IF
        END IF
        IF (Aout(idWamp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWamp(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWamp(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWamp)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWamp)
          END IF
        END IF
        IF (Aout(idWam2,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWam2(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWam2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWam2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWam2)
          END IF
        END IF
        IF (Aout(idWlen,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWlen(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWlen(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWlen)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWlen)
          END IF
        END IF
        IF (Aout(idWdir,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdir(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdir(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdir)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdir)
          END IF
        END IF
        IF (Aout(idWdip,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdip(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdip(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdip)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdip)
          END IF
        END IF
        IF (Aout(idWptp,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWptp(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWptp(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWptp)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWptp)
          END IF
        END IF
        IF (Aout(idWpbt,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWpbt(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWpbt(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWpbt)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWpbt)
          END IF
        END IF
        IF (Aout(idWorb,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWorb(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWorb(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWorb)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWorb)
          END IF
        END IF
        IF (Aout(idWdif,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdif(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdif(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdif)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdif)
          END IF
        END IF
        IF (Aout(idWdib,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdib(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdib(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdib)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdib)
          END IF
        END IF
        IF (Aout(idWdiw,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgWdiw(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgWdiw(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgWdiw)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgWdiw)
          END IF
        END IF
        IF (Aout(idUwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgUwav(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgUwav(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgUwav)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUwav)
          END IF
        END IF
        IF (Aout(idVwav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgVwav(i,j)=rfac(i,j)*                       &
     &                                  AVERAGE(ng)%avgVwav(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgVwav)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVwav)
          END IF
        END IF
!
!  Process vorticity fields.
!
        IF (Aout(id2dPV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgpvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgpvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgpvor2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor2d)
          END IF
        END IF
        IF (Aout(id2dRV,ng)) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              AVERAGE(ng)%avgrvor2d(i,j)=pfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgrvor2d(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgrvor2d)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor2d)
          END IF
        END IF
        IF (Aout(id3dPV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgpvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgpvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgpvor3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgpvor3d)
          END IF
        END IF
        IF (Aout(id3dRV,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgrvor3d(i,j,k)=pfac(i,j)*                 &
     &                                      AVERAGE(ng)%avgrvor3d(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_p3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgrvor3d)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgrvor3d)
          END IF
        END IF
!
!  Process quadratic fields.
!
        IF (Aout(idZZav,ng)) THEN
          DO j=JstrR,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgZZ(i,j)=rfac(i,j)*                         &
     &                               AVERAGE(ng)%avgZZ(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgZZ)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgZZ)
          END IF
        END IF
        IF (Aout(idU2av,ng)) THEN
          DO j=JstrR,JendR
            DO i=Istr,IendR
              AVERAGE(ng)%avgU2(i,j)=ufac(i,j)*                         &
     &                               AVERAGE(ng)%avgU2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgU2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgU2)
          END IF
        END IF
        IF (Aout(idV2av,ng)) THEN
          DO j=Jstr,JendR
            DO i=IstrR,IendR
              AVERAGE(ng)%avgV2(i,j)=vfac(i,j)*                         &
     &                               AVERAGE(ng)%avgV2(i,j)
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v2d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              AVERAGE(ng)%avgV2)
            CALL mp_exchange2d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgV2)
          END IF
        END IF
        IF (Aout(idUUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgUU(i,j,k)=ufac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUU(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUU)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUU)
          END IF
        END IF
        IF (Aout(idVVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgVV(i,j,k)=vfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgVV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgVV)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgVV)
          END IF
        END IF
        IF (Aout(idUVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                AVERAGE(ng)%avgUV(i,j,k)=rfac(i,j)*                     &
     &                                   AVERAGE(ng)%avgUV(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_r3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgUV)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgUV)
          END IF
        END IF
        IF (Aout(idHUav,ng)) THEN
          DO k=1,N(ng)
            DO j=JstrR,JendR
              DO i=Istr,IendR
                AVERAGE(ng)%avgHuon(i,j,k)=ufac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHuon(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_u3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHuon)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHuon)
          END IF
        END IF
        IF (Aout(idHVav,ng)) THEN
          DO k=1,N(ng)
            DO j=Jstr,JendR
              DO i=IstrR,IendR
                AVERAGE(ng)%avgHvom(i,j,k)=vfac(i,j)*                   &
     &                                     AVERAGE(ng)%avgHvom(i,j,k)
              END DO
            END DO
          END DO
          IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
            CALL exchange_v3d_tile (ng, tile,                           &
     &                              LBi, UBi, LBj, UBj, 1, N(ng),       &
     &                              AVERAGE(ng)%avgHvom)
            CALL mp_exchange3d (ng, tile, iNLM, 1,                      &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          NghostPoints,                           &
     &                          EWperiodic(ng), NSperiodic(ng),         &
     &                          AVERAGE(ng)%avgHvom)
          END IF
        END IF
        DO it=1,NT(ng)
          IF (Aout(idTTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgTT(i,j,k,it)=rfac(i,j)*                &
     &                                       AVERAGE(ng)%avgTT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_r3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgTT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgTT(:,:,:,it))
            END IF
          END IF
          IF (Aout(idUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgUT(i,j,k,it)=ufac(i,j)*                &
     &                                       AVERAGE(ng)%avgUT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgUT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgUT(:,:,:,it))
            END IF
          END IF
          IF (Aout(idVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgVT(i,j,k,it)=vfac(i,j)*                &
     &                                       AVERAGE(ng)%avgVT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgVT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgVT(:,:,:,it))
            END IF
          END IF
          IF (Aout(iHUTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=JstrR,JendR
                DO i=Istr,Iend
                  AVERAGE(ng)%avgHuonT(i,j,k,it)=ufac(i,j)*             &
     &                                    AVERAGE(ng)%avgHuonT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_u3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHuonT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHuonT(:,:,:,it))
            END IF
          END IF
          IF (Aout(iHVTav(it),ng)) THEN
            DO k=1,N(ng)
              DO j=Jstr,Jend
                DO i=IstrR,IendR
                  AVERAGE(ng)%avgHvomT(i,j,k,it)=vfac(i,j)*             &
     &                                    AVERAGE(ng)%avgHvomT(i,j,k,it)
                END DO
              END DO
            END DO
            IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
              CALL exchange_v3d_tile (ng, tile,                         &
     &                                LBi, UBi, LBj, UBj, 1, N(ng),     &
     &                                AVERAGE(ng)%avgHvomT(:,:,:,it))
              CALL mp_exchange3d (ng, tile, iNLM, 1,                    &
     &                            LBi, UBi, LBj, UBj, 1, N(ng),         &
     &                            NghostPoints,                         &
     &                            EWperiodic(ng), NSperiodic(ng),       &
     &                            AVERAGE(ng)%avgHvomT(:,:,:,it))
            END IF
          END IF
        END DO
      END IF
      RETURN
      END SUBROUTINE set_avg_tile
      END MODULE set_avg_mod
