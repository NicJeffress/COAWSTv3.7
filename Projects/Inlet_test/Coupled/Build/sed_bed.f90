      MODULE sed_bed_mod
!
!svn $Id: sed_bed.F 1054 2021-03-06 19:47:12Z arango $
!==================================================== John C. Warner ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group      Hernan G. Arango   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes sediment bed layer stratigraphy.              !
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
      PUBLIC  :: sed_bed
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sed_bed (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_sedbed
      USE mod_stepping
      USE mod_bbl
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/Sediment/sed_bed.F"
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
      CALL wclock_on (ng, iNLM, 16, 57, MyFile)
      CALL sed_bed_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
     &                   BBL(ng) % bustrc,                              &
     &                   BBL(ng) % bvstrc,                              &
     &                   BBL(ng) % bustrw,                              &
     &                   BBL(ng) % bvstrw,                              &
     &                   BBL(ng) % bustrcwmax,                          &
     &                   BBL(ng) % bvstrcwmax,                          &
     &                   OCEAN(ng) % t,                                 &
     &                   SEDBED(ng) % ero_flux,                         &
     &                   SEDBED(ng) % settling_flux,                    &
     &                   SEDBED(ng) % bed_thick,                        &
     &                   SEDBED(ng) % bed,                              &
     &                   SEDBED(ng) % bed_frac,                         &
     &                   SEDBED(ng) % bed_mass,                         &
     &                   SEDBED(ng) % bottom)
      CALL wclock_off (ng, iNLM, 16, 90, MyFile)
!
      RETURN
      END SUBROUTINE sed_bed
!
!***********************************************************************
      SUBROUTINE sed_bed_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
     &                         bustrc, bvstrc,                          &
     &                         bustrw, bvstrw,                          &
     &                         bustrcwmax, bvstrcwmax,                  &
     &                         t,                                       &
     &                         ero_flux, settling_flux,                 &
     &                         bed_thick,                               &
     &                         bed, bed_frac, bed_mass,                 &
     &                         bottom)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sediment
!
      USE bc_3d_mod, ONLY : bc_r3d_tile
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
      real(r8), intent(in) :: bustrc(LBi:,LBj:)
      real(r8), intent(in) :: bvstrc(LBi:,LBj:)
      real(r8), intent(in) :: bustrw(LBi:,LBj:)
      real(r8), intent(in) :: bvstrw(LBi:,LBj:)
      real(r8), intent(in) :: bustrcwmax(LBi:,LBj:)
      real(r8), intent(in) :: bvstrcwmax(LBi:,LBj:)
      real(r8), intent(inout):: bed_thick(LBi:,LBj:,:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ero_flux(LBi:,LBj:,:)
      real(r8), intent(inout) :: settling_flux(LBi:,LBj:,:)
      real(r8), intent(inout) :: bed(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_frac(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: bed_mass(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: bottom(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Ksed, i, ised, j, k, ks
      integer :: bnew
      real(r8), parameter :: eps = 1.0E-14_r8
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: thck_avail, thck_to_add
      real(r8), dimension(IminS:ImaxS,NST) :: dep_mass
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tau_w
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
      bnew=nnew
!
!-----------------------------------------------------------------------
! Compute sediment bed layer stratigraphy.
!-----------------------------------------------------------------------
!
      DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
          tau_w(i,j)=SQRT(bustrcwmax(i,j)*bustrcwmax(i,j)+              &
     &                    bvstrcwmax(i,j)*bvstrcwmax(i,j))
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Update bed properties according to ero_flux and dep_flux.
!-----------------------------------------------------------------------
!
      J_LOOP : DO j=Jstr,Jend
        SED_LOOP: DO ised=1,NST
!
!  The deposition and resuspension of sediment on the bottom "bed"
!  is due to precepitation flux FC(:,0), already computed, and the
!  resuspension (erosion, hence called ero_flux). The resuspension is
!  applied to the bottom-most grid box value qc(:,1) so the total mass
!  is conserved. Restrict "ero_flux" so that "bed" cannot go negative
!  after both fluxes are applied.
!
          DO i=Istr,Iend
            dep_mass(i,ised)=0.0_r8
!
! Apply morphology factor.
!
            ero_flux(i,j,ised)=ero_flux(i,j,ised)*morph_fac(ised,ng)
            settling_flux(i,j,ised)=settling_flux(i,j,ised)*            &
     &                              morph_fac(ised,ng)
!
            IF ((ero_flux(i,j,ised)-settling_flux(i,j,ised)).lt.        &
     &           0.0_r8) THEN
!
!  If first time step of deposit, then store deposit material in
!  temporary array, dep_mass.
!
              IF ((time(ng).gt.(bed(i,j,1,iaged)+1.1_r8*dt(ng))).and.   &
     &            (bed(i,j,1,ithck).gt.newlayer_thick(ng))) THEN
                dep_mass(i,ised)=settling_flux(i,j,ised)-               &
     &                           ero_flux(i,j,ised)
              END IF
                bed(i,j,1,iaged)=time(ng)
            END IF
!
!  Update bed mass arrays.
!
            bed_mass(i,j,1,nnew,ised)=MAX(bed_mass(i,j,1,bnew,ised)-    &
     &                                   (ero_flux(i,j,ised)-           &
     &                                    settling_flux(i,j,ised)),     &
     &                                    0.0_r8)
            DO k=2,Nbed
              bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k,nstp,ised)
            END DO
          END DO
        END DO SED_LOOP
!
!  If first time step of deposit, create new layer and combine bottom
!  two bed layers.
!
        DO i=Istr,Iend
          cff=0.0_r8
!
!  Determine if deposition ocurred here.
!
          IF (Nbed.gt.1) THEN
            DO ised=1,NST
              cff=cff+dep_mass(i,ised)
            END DO
            IF (cff.gt.0.0_r8) THEN
!
!  Combine bottom layers.
!
              bed(i,j,Nbed,iporo)=0.5_r8*(bed(i,j,Nbed-1,iporo)+        &
     &                                    bed(i,j,Nbed,iporo))
              bed(i,j,Nbed,iaged)=0.5_r8*(bed(i,j,Nbed-1,iaged)+        &
     &                                    bed(i,j,Nbed,iaged))
              DO ised=1,NST
                bed_mass(i,j,Nbed,nnew,ised)=                           &
     &                             bed_mass(i,j,Nbed-1,nnew,ised)+      &
     &                             bed_mass(i,j,Nbed  ,nnew,ised)
              END DO
!
!  Push layers down.
!
              DO k=Nbed-1,2,-1
                bed(i,j,k,iporo)=bed(i,j,k-1,iporo)
                bed(i,j,k,iaged)=bed(i,j,k-1,iaged)
                DO ised =1,NST
                  bed_mass(i,j,k,nnew,ised)=bed_mass(i,j,k-1,nnew,ised)
                END DO
              END DO
!
!  Set new top layer parameters.
!
              DO ised=1,NST
                bed_mass(i,j,2,nnew,ised)=MAX(bed_mass(i,j,2,nnew,ised)-&
     &                                    dep_mass(i,ised),0.0_r8)
                bed_mass(i,j,1,nnew,ised)=dep_mass(i,ised)
              END DO
            END IF
          END IF !NBED=1
!
! Recalculate thickness and fractions for all layers.
!
          DO k=1,Nbed
            cff3=0.0_r8
            DO ised=1,NST
              cff3=cff3+bed_mass(i,j,k,nnew,ised)
            END DO
            IF (cff3.eq.0.0_r8) THEN
              cff3=eps
            END IF
            bed(i,j,k,ithck)=0.0_r8
            DO ised=1,NST
              bed_frac(i,j,k,ised)=bed_mass(i,j,k,nnew,ised)/cff3
              bed(i,j,k,ithck)=MAX(bed(i,j,k,ithck)+                    &
     &                         bed_mass(i,j,k,nnew,ised)/               &
     &                         (Srho(ised,ng)*                          &
     &                          (1.0_r8-bed(i,j,k,iporo))),0.0_r8)
            END DO
          END DO
        END DO
      END DO J_LOOP
!
!  End of Suspended Sediment only section.
!
!
!  Ensure top bed layer thickness is greater or equal than active layer
!  thickness. If need to add sed to top layer, then entrain from lower
!  levels. Create new layers at bottom to maintain Nbed.
!
      J_LOOP2 : DO j=Jstr,Jend
        DO i=Istr,Iend
!
!  Calculate active layer thickness, bottom(i,j,iactv).
!
          bottom(i,j,iactv)=MAX(0.0_r8,                                 &
     &                          0.007_r8*                               &
     &                          (tau_w(i,j)-bottom(i,j,itauc))*rho0)+   &
     &                          6.0_r8*bottom(i,j,isd50)
!
!
! Apply morphology factor.
!
          bottom(i,j,iactv)=MAX(bottom(i,j,iactv)*morph_fac(1,ng),      &
     &                          bottom(i,j,iactv))
!
          IF (bottom(i,j,iactv).gt.bed(i,j,1,ithck)) THEN
            IF (Nbed.eq.1) THEN
              bottom(i,j,iactv)=bed(i,j,1,ithck)
            ELSE
              thck_to_add=bottom(i,j,iactv)-bed(i,j,1,ithck)
              thck_avail=0.0_r8
              Ksed=1                                        ! initialize
              DO k=2,Nbed
                IF (thck_avail.lt.thck_to_add) THEN
                  thck_avail=thck_avail+bed(i,j,k,ithck)
                  Ksed=k
                END IF
              END DO
!
!  Catch here if there was not enough bed material.
!
              IF (thck_avail.lt.thck_to_add) THEN
                bottom(i,j,iactv)=bed(i,j,1,ithck)+thck_avail
                thck_to_add=thck_avail
              END IF
!
!  Update bed mass of top layer and fractional layer.
!
              cff2=MAX(thck_avail-thck_to_add,0.0_r8)/                  &
     &             MAX(bed(i,j,Ksed,ithck),eps)
              DO ised=1,NST
                cff1=0.0_r8
                DO k=1,Ksed
                  cff1=cff1+bed_mass(i,j,k,nnew,ised)
                END DO
                cff3=cff2*bed_mass(i,j,Ksed,nnew,ised)
                bed_mass(i,j,1   ,nnew,ised)=cff1-cff3
                bed_mass(i,j,Ksed,nnew,ised)=cff3
              END DO
!
!  Update thickness of fractional layer ksource_sed.
!
              bed(i,j,Ksed,ithck)=MAX(thck_avail-thck_to_add,0.0_r8)
!
!  Upate bed fraction of top layer.
!
              cff3=0.0_r8
              DO ised=1,NST
                cff3=cff3+bed_mass(i,j,1,nnew,ised)
              END DO
              IF (cff3.eq.0.0_r8) THEN
                cff3=eps
              END IF
              DO ised=1,NST
                bed_frac(i,j,1,ised)=bed_mass(i,j,1,nnew,ised)/cff3
              END DO
!
!  Upate bed thickness of top layer.
!
              bed(i,j,1,ithck)=bottom(i,j,iactv)
!
!  Pull all layers closer to the surface.
!
              DO k=Ksed,Nbed
                ks=Ksed-2
                bed(i,j,k-ks,ithck)=bed(i,j,k,ithck)
                bed(i,j,k-ks,iporo)=bed(i,j,k,iporo)
                bed(i,j,k-ks,iaged)=bed(i,j,k,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k-ks,ised)=bed_frac(i,j,k,ised)
                  bed_mass(i,j,k-ks,nnew,ised)=bed_mass(i,j,k,nnew,ised)
                END DO
              END DO
!
!  Add new layers onto the bottom. Split what was in the bottom layer to
!  fill these new empty cells. ("ks" is the number of new layers).
!
              ks=Ksed-2
              cff=1.0_r8/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(i,j,k,ithck)=bed(i,j,Nbed-ks,ithck)*cff
                bed(i,j,k,iaged)=bed(i,j,Nbed-ks,iaged)
                DO ised=1,NST
                  bed_frac(i,j,k,ised)=bed_frac(i,j,Nbed-ks,ised)
                  bed_mass(i,j,k,nnew,ised)=                            &
     &                             bed_mass(i,j,Nbed-ks,nnew,ised)*cff
                END DO
              END DO
            END IF  ! Nbed > 1
          END IF  ! increase top bed layer
        END DO
      END DO J_LOOP2
!
!-----------------------------------------------------------------------
! Store old bed thickness.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
            bed_thick(i,j,3)=bed_thick(i,j,nnew)
            bed_thick(i,j,nnew)=0.0_r8
            DO k=1,Nbed
              bed_thick(i,j,nnew)=bed_thick(i,j,nnew)+                  &
     &                            bed(i,j,k,ithck)
            END DO
          END DO
        END DO
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            bed_thick(:,:,nnew))
        END IF
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            bed_thick(:,:,3))
        END IF
!
!-----------------------------------------------------------------------
!  Apply periodic or gradient boundary conditions to property arrays.
!-----------------------------------------------------------------------
!
      DO ised=1,NST
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_frac(:,:,:,ised))
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed_mass(:,:,:,nnew,ised))
      END DO
      CALL mp_exchange4d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, NST,          &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed_frac,                                     &
     &                    bed_mass(:,:,:,nnew,:))
      DO i=1,MBEDP
        CALL bc_r3d_tile (ng, tile,                                     &
     &                    LBi, UBi, LBj, UBj, 1, Nbed,                  &
     &                    bed(:,:,:,i))
      END DO
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 1, Nbed, 1, MBEDP,        &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    bed)
!
      RETURN
      END SUBROUTINE sed_bed_tile
      END MODULE sed_bed_mod
