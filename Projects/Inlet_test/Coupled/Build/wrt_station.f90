      SUBROUTINE wrt_station (ng)
!
!svn $Id: wrt_station.F 1054 2021-03-06 19:47:12Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2021 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes out data into stations NetCDF file.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_bbl
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_mixing
      USE mod_ncparam
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_sedbed
      USE mod_sediment
      USE mod_stepping
!
      USE extract_sta_mod, ONLY : extract_sta2d, extract_sta3d
      USE uv_rotate_mod,   ONLY : uv_rotate2d
      USE uv_rotate_mod,   ONLY : uv_rotate3d
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB, NposR, NposW, LBi, UBi, LBj, UBj
      integer :: Fcount, i, ifield, k, np, status, tile
!
      real(dp) :: scale
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
      real(r8), dimension(Nstation(ng)*Nbed) :: XposB, YposB, ZposB
      real(r8), dimension(Nstation(ng)*Nbed) :: bsta
      real(r8), dimension(Nstation(ng)*(N(ng))) :: XposR, YposR, ZposR
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: XposW, YposW, ZposW
      real(r8), dimension(Nstation(ng)*(N(ng)+1)) :: rsta
      real(r8), allocatable :: Ur2d(:,:)
      real(r8), allocatable :: Vr2d(:,:)
      real(r8), allocatable :: Ur3d(:,:,:)
      real(r8), allocatable :: Vr3d(:,:,:)
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Utility/wrt_station.F"
!
      SourceFile=MyFile
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out station data at RHO-points.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, 93, MyFile)) RETURN
!
!  Set time record index.
!
      STA(ng)%Rindex=STA(ng)%Rindex+1
      Fcount=STA(ng)%Fcount
      STA(ng)%Nrec(Fcount)=STA(ng)%Nrec(Fcount)+1
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
      Cgrid=.FALSE.
!
!  Set positions for generic extraction routine.
!
      NposB=Nstation(ng)*Nbed
      NposR=Nstation(ng)*N(ng)
      NposW=Nstation(ng)*(N(ng)+1)
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
        DO k=1,N(ng)
          np=k+(i-1)*N(ng)
          XposR(np)=SCALARS(ng)%SposX(i)
          YposR(np)=SCALARS(ng)%SposY(i)
          ZposR(np)=REAL(k,r8)
        END DO
        DO k=0,N(ng)
          np=k+1+(i-1)*(N(ng)+1)
          XposW(np)=SCALARS(ng)%SposX(i)
          YposW(np)=SCALARS(ng)%SposY(i)
          ZposW(np)=REAL(k,r8)
        END DO
        DO k=1,Nbed
          np=k+(i-1)*Nbed
          XposB(np)=SCALARS(ng)%SposX(i)
          YposB(np)=SCALARS(ng)%SposY(i)
          ZposB(np)=REAL(k,r8)
        END DO
      END DO
!
!  Write out model time (s).
!
      CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                     &
     &                      TRIM(Vname(1,idtime)), time(ng:),           &
     &                      (/STA(ng)%Rindex/), (/1/),                  &
     &                      ncid = STA(ng)%ncid,                        &
     &                      varid = STA(ng)%Vid(idtime))
      IF (FoundError(exit_flag, NoError, 150, MyFile)) RETURN
!
!  Write out free-surface (m).
!
      IF (Sout(idFsur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idFsur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zeta(:,:,kstp(ng)),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idFsur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idFsur))
        IF (FoundError(exit_flag, NoError, 165, MyFile)) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idWztw, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%zetaw(:,:),                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWztw)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWztw))
        IF (FoundError(exit_flag, NoError, 176, MyFile)) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idWqsp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%qsp(:,:),                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWqsp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWqsp))
        IF (FoundError(exit_flag, NoError, 187, MyFile)) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idWbeh, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%bh(:,:),                   &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWbeh)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWbeh))
        IF (FoundError(exit_flag, NoError, 198, MyFile)) RETURN
      END IF
!
!  Define time-varying bathymetry.
!
      IF (Sout(idbath,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idbath, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%h,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idbath)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idbath))
        IF (FoundError(exit_flag, NoError, 218, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbar, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%ubar(:,:,kstp(ng)),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbar))
        IF (FoundError(exit_flag, NoError, 235, MyFile)) RETURN
      END IF
!
!  Write out 2D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVbar,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbar, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng)%vbar(:,:,kstp(ng)),            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbar)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbar))
        IF (FoundError(exit_flag, NoError, 251, MyFile)) RETURN
      END IF
!
!  Write out 2D Eastward and Northward momentum components (m/s) at
!  RHO-points
!
      IF (Sout(idu2dE,ng).and.Sout(idv2dN,ng)) THEN
        IF (.not.allocated(Ur2d)) THEN
          allocate (Ur2d(LBi:UBi,LBj:UBj))
            Ur2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        IF (.not.allocated(Vr2d)) THEN
          allocate (Vr2d(LBi:UBi,LBj:UBj))
            Vr2d(LBi:UBi,LBj:UBj)=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate2d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % ubar(:,:,kstp(ng)),                   &
     &                    OCEAN(ng) % vbar(:,:,kstp(ng)),                   &
     &                    Ur2d, Vr2d)
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idu2dE, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Ur2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idu2dE)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu2dE))
        IF (FoundError(exit_flag, NoError, 292, MyFile)) RETURN
        CALL extract_sta2d (ng, iNLM, Cgrid, idv2dN, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, Vr2d,                                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idv2dN)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv2dN))
        IF (FoundError(exit_flag, NoError, 303, MyFile)) RETURN
        deallocate (Ur2d)
        deallocate (Vr2d)
      END IF
!
!  Write out 3D momentum component (m/s) in the XI-direction.
!
      IF (Sout(idUvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idUvel, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%u(:,:,:,nrhs(ng)),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUvel))
        IF (FoundError(exit_flag, NoError, 325, MyFile)) RETURN
      END IF
!
!  Write out 3D momentum component (m/s) in the ETA-direction.
!
      IF (Sout(idVvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvel, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%v(:,:,:,nrhs(ng)),             &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvel))
        IF (FoundError(exit_flag, NoError, 342, MyFile)) RETURN
      END IF
!
!  Write out 3D Eastward and Northward momentum components (m/s) at
!  RHO-points.
!
      IF (Sout(idu3dE,ng).and.Sout(idv3dN,ng)) THEN
        IF (.not.allocated(Ur3d)) THEN
          allocate (Ur3d(LBi:UBi,LBj:UBj,N(ng)))
          Ur3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        IF (.not.allocated(Vr3d)) THEN
          allocate (Vr3d(LBi:UBi,LBj:UBj,N(ng)))
          Vr3d(LBi:UBi,LBj:UBj,1:N(ng))=0.0_r8
        END IF
        tile=MyRank
        CALL uv_rotate3d (ng, tile, .FALSE., .TRUE.,                    &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    GRID(ng) % CosAngler,                         &
     &                    GRID(ng) % SinAngler,                         &
     &                    GRID(ng) % rmask_full,                        &
     &                    OCEAN(ng) % u(:,:,:,nrhs(ng)),                    &
     &                    OCEAN(ng) % v(:,:,:,nrhs(ng)),                    &
     &                    Ur3d, Vr3d)
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idu3dE, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Ur3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idu3dE)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idu3dE))
        IF (FoundError(exit_flag, NoError, 384, MyFile)) RETURN
        CALL extract_sta3d (ng, iNLM, Cgrid, idv3dN, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, Vr3d,                                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idv3dN)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idv3dN))
        IF (FoundError(exit_flag, NoError, 396, MyFile)) RETURN
        deallocate (Ur3d)
        deallocate (Vr3d)
      END IF
!
!  Write out vertical velocity (m/s).
!
      IF (Sout(idWvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idWvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%wvel,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWvel))
        IF (FoundError(exit_flag, NoError, 416, MyFile)) RETURN
      END IF
!
!  Write S-coordinate "omega" vertical velocity (m3/s).
!
      IF (Sout(idOvel,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idOvel, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng)%W,                         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idOvel)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idOvel))
        IF (FoundError(exit_flag, NoError, 433, MyFile)) RETURN
      END IF
!
!  Write out tracer type variables.
!
      DO i=1,NT(ng)
        ifield=idTvar(i)
        IF (Sout(ifield,ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, ifield, r3dvar,          &
     &                        LBi, UBi, LBj, UBj, 1, N(ng),             &
     &                        scale, OCEAN(ng)%t(:,:,:,nrhs(ng),i),         &
     &                        NposR, XposR, YposR, ZposR, rsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idTvar(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/N(ng),Nstation(ng),1/),               &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Tid(i))
          IF (FoundError(exit_flag, NoError, 452, MyFile)) RETURN
        END IF
      END DO
!
!  Write out density anomaly.
!
      IF (Sout(idDano,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idDano, r3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng)%rho,                       &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idDano)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idDano))
        IF (FoundError(exit_flag, NoError, 470, MyFile)) RETURN
      END IF
!
!  Write out vertical viscosity coefficient.
!
      IF (Sout(idVvis,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idVvis, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akv,                      &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVvis)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVvis))
        IF (FoundError(exit_flag, NoError, 624, MyFile)) RETURN
      END IF
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      IF (Sout(idTdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idTdif, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%Akt(:,:,:,itemp),         &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idTdif)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idTdif))
        IF (FoundError(exit_flag, NoError, 641, MyFile)) RETURN
      END IF
!
!  Write out turbulent kinetic energy.
!
      IF (Sout(idMtke,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtke, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%tke(:,:,:,nrhs(ng)),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idMtke)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idMtke))
        IF (FoundError(exit_flag, NoError, 679, MyFile)) RETURN
      END IF
!
!  Write out turbulent kinetic energy times length scale.
!
      IF (Sout(idMtls,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idMtls, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, MIXING(ng)%gls(:,:,:,nrhs(ng)),          &
     &                      NposW, XposW, YposW, ZposW, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idMtls)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng)+1,Nstation(ng),1/),               &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idMtls))
        IF (FoundError(exit_flag, NoError, 696, MyFile)) RETURN
      END IF
!
!  Write out surface net heat flux.
!
      IF (Sout(idTsur(itemp),ng)) THEN
        ifield=idTsur(itemp)
        scale=rho0*Cp
        CALL extract_sta2d (ng, iNLM, Cgrid, ifield, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%stflx(:,:,itemp),         &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,ifield)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(ifield))
        IF (FoundError(exit_flag, NoError, 819, MyFile)) RETURN
      END IF
!
!  Write out surface U-momentum stress.
!
      IF (Sout(idUsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUsms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%sustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUsms))
        IF (FoundError(exit_flag, NoError, 976, MyFile)) RETURN
      END IF
!
!  Write out surface V-momentum stress.
!
      IF (Sout(idVsms,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVsms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%svstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVsms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVsms))
        IF (FoundError(exit_flag, NoError, 992, MyFile)) RETURN
      END IF
!
!  Write out bottom U-momentum stress.
!
      IF (Sout(idUbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbms, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bustr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbms))
        IF (FoundError(exit_flag, NoError, 1008, MyFile)) RETURN
      END IF
!
!  Write out bottom V-momentum stress.
!
      IF (Sout(idVbms,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbms, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng)%bvstr,                    &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbms)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbms))
        IF (FoundError(exit_flag, NoError, 1024, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom U-stress.
!
      IF (Sout(idUbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbrs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbrs))
        IF (FoundError(exit_flag, NoError, 1043, MyFile)) RETURN
      END IF
!
!  Write out current-induced, bottom V-stress.
!
      IF (Sout(idVbrs,ng)) THEN
        scale=-rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbrs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrc,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbrs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbrs))
        IF (FoundError(exit_flag, NoError, 1059, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom U-stress.
!
      IF (Sout(idUbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbws, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbws)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbws))
        IF (FoundError(exit_flag, NoError, 1075, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bottom V-wave stress.
!
      IF (Sout(idVbws,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbws, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrw,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbws)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbws))
        IF (FoundError(exit_flag, NoError, 1091, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom U-stress.
!
      IF (Sout(idUbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbcs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bustrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbcs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbcs))
        IF (FoundError(exit_flag, NoError, 1107, MyFile)) RETURN
      END IF
!
!  Write out maximum wind and current, bottom V-stress.
!
      IF (Sout(idVbcs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbcs, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%bvstrcwmax,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbcs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbcs))
        IF (FoundError(exit_flag, NoError, 1123, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital U-velocity.
!
      IF (Sout(idUbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ubot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbot)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbot))
        IF (FoundError(exit_flag, NoError, 1139, MyFile)) RETURN
      END IF
!
!  Write out wind-induced, bed wave orbital V-velocity.
!
      IF (Sout(idVbot,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbot, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vbot,                        &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbot)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbot))
        IF (FoundError(exit_flag, NoError, 1155, MyFile)) RETURN
      END IF
!
!  Write out bottom U-velocity above bed.
!
      IF (Sout(idUbur,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idUbur, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Ur,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idUbur)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idUbur))
        IF (FoundError(exit_flag, NoError, 1171, MyFile)) RETURN
      END IF
!
!  Write out bottom V-velocity above bed.
!
      IF (Sout(idVbvr,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idVbvr, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, BBL(ng)%Vr,                          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idVbvr)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idVbvr))
        IF (FoundError(exit_flag, NoError, 1187, MyFile)) RETURN
      END IF
!
!  Write out sediment fraction of each size class in each bed layer.
!
      DO i=1,NST
        IF (Sout(idfrac(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, idfrac(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed_frac(:,:,:,i),      &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idfrac(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idfrac(i)))
          IF (FoundError(exit_flag, NoError, 1208, MyFile)) RETURN
        END IF
!
!  Write out sediment mass of each size class in each bed layer.
!
        IF (Sout(idBmas(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, idBmas(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale,                                    &
     &                        SEDBED(ng)%bed_mass(:,:,:,nrhs(ng),i),        &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idBmas(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idBmas(i)))
          IF (FoundError(exit_flag, NoError, 1226, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment properties in each bed layer.
!
      DO i=1,MBEDP
        IF (Sout(idSbed(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, iNLM, Cgrid, idSbed(i), b3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, Nbed,              &
     &                        scale, SEDBED(ng)%bed(:,:,:,i),           &
     &                        NposB, XposB, YposB, ZposB, bsta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idSbed(i))), rsta,         &
     &                          (/1,1,STA(ng)%Rindex/),                 &
     &                          (/Nbed,Nstation(ng),1/),                &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idSbed(i)))
          IF (FoundError(exit_flag, NoError, 1245, MyFile)) RETURN
        END IF
      END DO
!
!  Write out exposed sediment layer properties.
!
      DO i=1,MBEDP
        IF (Sout(idBott(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, iNLM, Cgrid, idBott(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBED(ng)%bottom(:,:,i),          &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                 &
     &                          TRIM(Vname(1,idBott(i))), rsta,         &
     &                          (/1,STA(ng)%Rindex/),                   &
     &                          (/Nstation(ng),1/),                     &
     &                          ncid = STA(ng)%ncid,                    &
     &                          varid = STA(ng)%Vid(idBott(i)))
          IF (FoundError(exit_flag, NoError, 1267, MyFile)) RETURN
        END IF
      END DO
!
!  Write out 2D U-momentum Stokes drift velocity.
!
      IF (Sout(idU2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idU2Sd, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % ubar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU2Sd)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU2Sd))
        IF (FoundError(exit_flag, NoError, 1456, MyFile)) RETURN
      END IF
!
!  Write out 2D V-momentum Stokes drift velocity.
!
      IF (Sout(idV2Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idV2Sd, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % vbar_stokes,             &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV2Sd)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV2Sd))
        IF (FoundError(exit_flag, NoError, 1472, MyFile)) RETURN
      END IF
!
!  Write out total 2D  U-stress.
!
      IF (Sout(idU2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idU2rs, u2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rustr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU2rs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU2rs))
        IF (FoundError(exit_flag, NoError, 1488, MyFile)) RETURN
      END IF
!
!  Write out total 2D  V-stress.
!
      IF (Sout(idV2rs,ng)) THEN
        scale=rho0
        CALL extract_sta2d (ng, iNLM, Cgrid, idV2rs, v2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, MIXING(ng) % rvstr2d,                &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV2rs)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV2rs))
        IF (FoundError(exit_flag, NoError, 1504, MyFile)) RETURN
      END IF
!
!  Write out 3D U-momentum Stokes drift velocity.
!
      IF (Sout(idU3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idU3Sd, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % u_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU3Sd)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU3Sd))
        IF (FoundError(exit_flag, NoError, 1522, MyFile)) RETURN
      END IF
!
!  Write out 3D V-momentum stokes velocity.
!
      IF (Sout(idV3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idV3Sd, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, OCEAN(ng) % v_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV3Sd)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV3Sd))
        IF (FoundError(exit_flag, NoError, 1539, MyFile)) RETURN
      END IF
!
!  Write out 3D Omega-momentum stokes velocity.
!
      IF (Sout(idW3Sd,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3Sd, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % W_stokes,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3Sd)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3Sd))
        IF (FoundError(exit_flag, NoError, 1556, MyFile)) RETURN
      END IF
!
!  Write out 3D W-momentum stokes velocity
!
      IF (Sout(idW3St,ng)) THEN
        scale=1.0_dp
        CALL extract_sta3d (ng, iNLM, Cgrid, idW3St, w3dvar,            &
     &                      LBi, UBi, LBj, UBj, 0, N(ng),               &
     &                      scale, OCEAN(ng) % wstvel,                  &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idW3St)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idW3St))
        IF (FoundError(exit_flag, NoError, 1573, MyFile)) RETURN
      END IF
!
!  Write out 3D total  U-stress.
!
      IF (Sout(idU3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, iNLM, Cgrid, idU3rs, u3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rustr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idU3rs)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idU3rs))
        IF (FoundError(exit_flag, NoError, 1590, MyFile)) RETURN
      END IF
!
!  Write out 3D total  V-stress.
!
      IF (Sout(idV3rs,ng)) THEN
        scale=rho0
        CALL extract_sta3d (ng, iNLM, Cgrid, idV3rs, v3dvar,            &
     &                      LBi, UBi, LBj, UBj, 1, N(ng),               &
     &                      scale, MIXING(ng) % rvstr3d,                &
     &                      NposR, XposR, YposR, ZposR, rsta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idV3rs)), rsta,              &
     &                        (/1,1,STA(ng)%Rindex/),                   &
     &                        (/N(ng),Nstation(ng),1/),                 &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idV3rs))
        IF (FoundError(exit_flag, NoError, 1607, MyFile)) RETURN
      END IF
!
!  Write out wind-induced wave height.
!
      IF (Sout(idWamp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWamp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Hwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWamp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWamp))
        IF (FoundError(exit_flag, NoError, 1626, MyFile)) RETURN
      END IF
!
!  Write out wind-induced wave length.
!
      IF (Sout(idWlen,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWlen, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Lwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWlen)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWlen))
        IF (FoundError(exit_flag, NoError, 1644, MyFile)) RETURN
      END IF
!
!  Write out wind-induced wave direction.
!
      IF (Sout(idWdir,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdir, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwave,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdir)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdir))
        IF (FoundError(exit_flag, NoError, 1680, MyFile)) RETURN
      END IF
!
!  Write out wind-induced peak wave direction.
!
      IF (Sout(idWdip,ng)) THEN
        scale=rad2deg
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdip, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dwavep,                  &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdip)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdip))
        IF (FoundError(exit_flag, NoError, 1698, MyFile)) RETURN
      END IF
!
!  Write out wind-induced surface wave period.
!
      IF (Sout(idWptp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWptp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_top,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWptp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWptp))
        IF (FoundError(exit_flag, NoError, 1716, MyFile)) RETURN
      END IF
!
!  Write out wind-induced bottom wave period.
!
      IF (Sout(idWpbt,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWpbt, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Pwave_bot,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWpbt)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWpbt))
        IF (FoundError(exit_flag, NoError, 1734, MyFile)) RETURN
      END IF
!
!  Write out wind-induced wave bottom orbital velocity.
!
      IF (Sout(idWorb,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWorb, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Uwave_rms,              &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWorb)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWorb))
        IF (FoundError(exit_flag, NoError, 1752, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to bottom friction.
!
      IF (Sout(idWdif,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdif, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_fric,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdif)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdif))
        IF (FoundError(exit_flag, NoError, 1770, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to breaking.
!
      IF (Sout(idWdib,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdib, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_break,           &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdib)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdib))
        IF (FoundError(exit_flag, NoError, 1789, MyFile)) RETURN
      END IF
!
!  Write out wave dissipation due to white capping.
!
      IF (Sout(idWdiw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWdiw, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, FORCES(ng) % Dissip_wcap,            &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWdiw)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWdiw))
        IF (FoundError(exit_flag, NoError, 1805, MyFile)) RETURN
      END IF
!
!  Write out  quasi-static sea level adjustment.
!
      IF (Sout(idWztw,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWztw, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % zeta(:,:,kstp(ng)),          &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWztw)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWztw))
        IF (FoundError(exit_flag, NoError, 1908, MyFile)) RETURN
      END IF
!
!  Write out  quasi-static pressure.
!
      IF (Sout(idWqsp,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWqsp, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % qsp,                     &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWqsp)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWqsp))
        IF (FoundError(exit_flag, NoError, 1924, MyFile)) RETURN
      END IF
!
!  Write out  Bernoulli head.
!
      IF (Sout(idWbeh,ng)) THEN
        scale=1.0_dp
        CALL extract_sta2d (ng, iNLM, Cgrid, idWbeh, r2dvar,            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, OCEAN(ng) % bh,                      &
     &                      Nstation(ng), Xpos, Ypos, psta)
        CALL netcdf_put_fvar (ng, iNLM, STA(ng)%name,                   &
     &                        TRIM(Vname(1,idWbeh)), psta,              &
     &                        (/1,STA(ng)%Rindex/), (/Nstation(ng),1/), &
     &                        ncid = STA(ng)%ncid,                      &
     &                        varid = STA(ng)%Vid(idWbeh))
        IF (FoundError(exit_flag, NoError, 1940, MyFile)) RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize stations NetCDF file to disk.
!-----------------------------------------------------------------------
!
      CALL netcdf_sync (ng, iNLM, STA(ng)%name, STA(ng)%ncid)
      RETURN
      END SUBROUTINE wrt_station
