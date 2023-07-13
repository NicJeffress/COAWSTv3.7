      SUBROUTINE read_SedPar (model, inp, out, Lwrite)
!
!svn $Id: sediment_inp.h 1054 2021-03-06 19:47:12Z arango $
!=======================================================================
!                                                                      !
!  This routine reads in cohesive and non-cohesive sediment model      !
!  parameters.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
      USE mod_sedflocs
      USE inp_decode_mod
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status
      logical, dimension(Ngrids) :: Lbed
      logical, dimension(Ngrids) :: Lbottom
      logical, dimension(NCS,Ngrids) :: Lmud
      logical, dimension(NNS,Ngrids) :: Lsand
      real(r8), dimension(Ngrids) :: Rbed
      real(r8), dimension(NCS,Ngrids) :: Rmud
      real(r8), dimension(NNS,Ngrids) :: Rsand
      real(dp), dimension(nRval) :: Rval
      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(nCval) :: Cval
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NST                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in cohesive and non-cohesive model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lsediment')
              Npts=load_l(Nval, Cval, Ngrids, Lsediment)
            CASE ('Hadvection')
              IF (itracer.lt.NST) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idsed(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idsed(iTrcStr), idsed(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Hadvection)
            CASE ('Vadvection')
              IF (itracer.lt.NST) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              itrc=idsed(itracer)
              Npts=load_tadv(Nval, Cval, line, nline, itrc, igrid,      &
     &                       itracer, idsed(iTrcStr), idsed(iTrcEnd),   &
     &                       Vname(1,idTvar(itrc)),                     &
     &                       Vadvection)
            CASE ('LBC(isTvar)')
              IF (itracer.lt.NST) THEN
                itracer=itracer+1
              ELSE
                itracer=1                      ! next nested grid
              END IF
              ifield=isTvar(idsed(itracer))
              Npts=load_lbc(Nval, Cval, line, nline, ifield, igrid,     &
     &                      idsed(iTrcStr), idsed(iTrcEnd),             &
     &                      Vname(1,idTvar(idsed(itracer))), LBC)
            CASE ('NEWLAYER_THICK')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                newlayer_thick(ng)=Rbed(ng)
              END DO
            CASE ('MINLAYER_THICK')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                minlayer_thick(ng)=Rbed(ng)
              END DO
            CASE ('BEDLOAD_COEFF')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
              DO ng=1,Ngrids
                bedload_coeff(ng)=Rbed(ng)
              END DO
            CASE ('SG_ZWBL')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  sg_zwbl(ng)=Rbed(ng)
                END DO
            CASE ('SEDSLOPE_CRIT_WET')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  sedslope_crit_wet(ng)=Rbed(ng)
                END DO
            CASE ('SEDSLOPE_CRIT_DRY')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  sedslope_crit_dry(ng)=Rbed(ng)
                END DO
            CASE ('SLOPEFAC_WET')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  slopefac_wet(ng)=Rbed(ng)
                END DO
            CASE ('SLOPEFAC_DRY')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  slopefac_dry(ng)=Rbed(ng)
                END DO
            CASE ('BEDLOAD_VANDERA_ALPHAW')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  bedload_vandera_alphaw(ng)=Rbed(ng)
                END DO
            CASE ('BEDLOAD_VANDERA_ALPHAC')
              Npts=load_r(Nval, Rval, Ngrids, Rbed)
                DO ng=1,Ngrids
                  bedload_vandera_alphac(ng)=Rbed(ng)
                END DO
            CASE ('Hout(ithck)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(ithck)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(iaged)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iaged)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(iporo)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iporo)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(idiff)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(idiff)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbed(ng)
              END DO
            CASE ('Hout(idsurs)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsurs
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsrrw)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsrrw
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsbtw)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsbtw
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsucr)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsucr
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsutr)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsutr
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idstcr)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idstcr
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsttr)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsttr
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(isd50)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(isd50)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idens)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idens)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
             END DO
            CASE ('Hout(iwsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iwsed)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(itauc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(itauc)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(irlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irlen)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(irhgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irhgt)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(ibwav)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ibwav)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izdef)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izdef)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izapp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izapp)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izNik)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izNik)
             DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbio)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbio)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbfm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbfm)
              DO ng=1,Ngrids
               Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izbld)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbld)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(izwbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izwbl)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(iactv)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iactv)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(ishgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ishgt)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(imaxD)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(imaxD)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idnet)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idnet)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idtbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idtbl)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idubl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idubl)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idfdw)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idfdw)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idzrw)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idzrw)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idksd)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idksd)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idusc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idusc)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idpcx)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idpcx)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Hout(idpwc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idpwc)
              DO ng=1,Ngrids
                Hout(i,ng)=Lbottom(ng)
              END DO
            CASE ('SAND_SD50')
              IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Sd50(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_CSED')
              IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand )
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Csed(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_SRHO')
              IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Srho(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_WSED')
              IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Wsed(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_ERATE')
              IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  Erate(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TAU_CE')
              IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  tau_ce(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TAU_CD')
              IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  tau_cd(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_POROS')
              IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  poros(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_RXN')
              IF (.not.allocated(sed_rxn)) allocate (sed_rxn(NST,Ngrids))
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  sed_rxn(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNU2')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  nl_tnu2(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNU4')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  nl_tnu4(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('ad_SAND_TNU2')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_tnu2(i,ng)=Rsand(itrc,ng)
                  tl_tnu2(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('ad_SAND_TNU4')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_tnu4(i,ng)=Rsand(itrc,ng)
                  tl_tnu4(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Sponge')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerSponge(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_AKT_BAK')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  Akt_bak(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_AKT_fac')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  ad_Akt_fac(i,ng)=Rsand(itrc,ng)
                  tl_Akt_fac(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_TNUDG')
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  Tnudg(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_MORPH_FAC')
              IF (.not.allocated(morph_fac)) THEN
                allocate (morph_fac(NST,Ngrids))
              END IF
              Npts=load_r(Nval, Rval, NNS, Ngrids, Rsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=NCS+itrc
                  morph_fac(i,ng)=Rsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Ltsrc', 'SAND_Ltracer')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerSrc(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Ltclm')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LtracerCLM(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('SAND_Tnudge')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsed(NCS+itrc)
                  LnudgeTCLM(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idsand)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSfrac)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idfrac(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSmass)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idBmas(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSUbld)')
              DO ng=1,Ngrids
                DO itrc=NCS+1,NST
                  IF (idUbld(itrc).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idUbld'
                    exit_flag=5
                    RETURN
                  END IF
                END DO
              END DO
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idUbld(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iSVbld)')
              DO ng=1,Ngrids
                DO itrc=NCS+1,NST
                  IF (idVbld(itrc).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idVbld'
                    exit_flag=5
                    RETURN
                  END IF
                END DO
              END DO
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idVbld(NCS+itrc)
                  Hout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idsand)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSsand)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idsurT(idsed(NCS+itrc))
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSfrac)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idfrac(NCS+itrc)
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSmass)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idBmas(NCS+itrc)
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSUbld)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idUbld(NCS+itrc)
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSVbld)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idVbld(NCS+itrc)
                  Qout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idsand)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTvar(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSTTav)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idTTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSUTav)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idUTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSVTav)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idVTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(SHUTav)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=iHUTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(SHVTav)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=iHVTav(idsed(NCS+itrc))
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSUbld)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idUbld(NCS+itrc)
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iSVbld)')
              Npts=load_l(Nval, Cval, NNS, Ngrids, Lsand)
              DO ng=1,Ngrids
                DO itrc=1,NNS
                  i=idVbld(NCS+itrc)
                  Aout(i,ng)=Lsand(itrc,ng)
                END DO
              END DO
            CASE ('MUD_SD50')
              IF (.not.allocated(Sd50)) allocate (Sd50(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Sd50(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_CSED')
              IF (.not.allocated(Csed)) allocate (Csed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud )
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Csed(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_SRHO')
              IF (.not.allocated(Srho)) allocate (Srho(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Srho(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_WSED')
              IF (.not.allocated(Wsed)) allocate (Wsed(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Wsed(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_ERATE')
              IF (.not.allocated(Erate)) allocate (Erate(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  Erate(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TAU_CE')
              IF (.not.allocated(tau_ce)) allocate (tau_ce(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  tau_ce(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TAU_CD')
              IF (.not.allocated(tau_cd)) allocate (tau_cd(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  tau_cd(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_POROS')
              IF (.not.allocated(poros)) allocate (poros(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  poros(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_RXN')
              IF (.not.allocated(sed_rxn)) allocate (sed_rxn(NST,Ngrids))
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  sed_rxn(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNU2')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  nl_tnu2(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNU4')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  nl_tnu4(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('ad_MUD_TNU2')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_tnu2(i,ng)=Rmud(itrc,ng)
                  tl_tnu2(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('ad_MUD_TNU4')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_tnu4(i,ng)=Rmud(itrc,ng)
                  tl_tnu4(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Sponge')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerSponge(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_AKT_BAK')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  Akt_bak(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_AKT_fac')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  ad_Akt_fac(i,ng)=Rmud(itrc,ng)
                  tl_Akt_fac(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_TNUDG')
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  Tnudg(i,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_MORPH_FAC')
              IF (.not.allocated(morph_fac)) THEN
                allocate (morph_fac(NST,Ngrids))
              END IF
              Npts=load_r(Nval, Rval, NCS, Ngrids, Rmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  morph_fac(itrc,ng)=Rmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Ltsrc', 'MUD_Ltracer')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerSrc(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Ltclm')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LtracerCLM(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('MUD_Tnudge')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsed(itrc)
                  LnudgeTCLM(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(idmud)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMfrac)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idfrac(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMmass)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idBmas(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMUbld)')
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  IF (idUbld(itrc).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idUbld'
                    exit_flag=5
                    RETURN
                  END IF
                END DO
              END DO
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idUbld(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Hout(iMVbld)')
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  IF (idVbld(itrc).eq.0) THEN
                    IF (Master) WRITE (out,30) 'idVbld'
                    exit_flag=5
                    RETURN
                  END IF
                END DO
              END DO
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idVbld(itrc)
                  Hout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(idmud)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iSmud)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idsurT(idsed(itrc))
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iMfrac)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idfrac(itrc)
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iMmass)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idBmas(itrc)
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iMUbld)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idUbld(itrc)
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(iMVbld)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idVbld(itrc)
                  Qout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(idmud)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTvar(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMTTav)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idTTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMUTav)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idUTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMVTav)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idVTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(MHUTav)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=iHUTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(MHVTav)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=iHVTav(idsed(itrc))
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMUbld)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idUbld(itrc)
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Aout(iMVbld)')
              Npts=load_l(Nval, Cval, NCS, Ngrids, Lmud)
              DO ng=1,Ngrids
                DO itrc=1,NCS
                  i=idVbld(itrc)
                  Aout(i,ng)=Lmud(itrc,ng)
                END DO
              END DO
            CASE ('Qout(ithck)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(ithck)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(iaged)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iaged)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(iporo)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(iporo)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(idiff)')
              Npts=load_l(Nval, Cval, Ngrids, Lbed)
              i=idSbed(idiff)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbed(ng)
              END DO
            CASE ('Qout(isd50)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(isd50)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(idens)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(idens)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
             END DO
            CASE ('Qout(iwsed)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iwsed)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(itauc)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(itauc)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(irlen)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irlen)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(irhgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(irhgt)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(ibwav)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ibwav)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izdef)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izdef)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izapp)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izapp)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izNik)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izNik)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izbio)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbio)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izbfm)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbfm)
              DO ng=1,Ngrids
               Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izbld)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izbld)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(izwbl)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(izwbl)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(iactv)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(iactv)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
            CASE ('Qout(ishgt)')
              Npts=load_l(Nval, Cval, Ngrids, Lbottom)
              i=idBott(ishgt)
              DO ng=1,Ngrids
                Qout(i,ng)=Lbottom(ng)
              END DO
          END SELECT
        END IF
      END DO
  10  IF (Master) WRITE (out,40) line
      exit_flag=4
      RETURN
  20  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Master.and.Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lsediment(ng)) THEN
            WRITE (out,50) ng
            WRITE (out,60)
            DO itrc=1,NST
              WRITE (out,70) itrc, Sd50(itrc,ng), Csed(itrc,ng),        &
     &                       Srho(itrc,ng), Wsed(itrc,ng),              &
     &                       Erate(itrc,ng), poros(itrc,ng)
            END DO
            WRITE (out,80)
            DO itrc=1,NST
              i=idsed(itrc)
              WRITE (out,70) itrc, tau_ce(itrc,ng), tau_cd(itrc,ng),    &
     &                       nl_tnu2(i,ng), nl_tnu4(i,ng),              &
     &                       Akt_bak(i,ng), Tnudg(i,ng)
            END DO
            WRITE (out,90)
            DO itrc=1,NST
              WRITE (out,70) itrc,  morph_fac(itrc,ng)
            END DO
            WRITE (out,100) newlayer_thick(ng)
            WRITE (out,110) minlayer_thick(ng)
            WRITE (out,120) bedload_coeff(ng)
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerSponge(i,ng)) THEN
                WRITE (out,150) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning ON  sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSponge(i,ng), 'LtracerSponge',   &
     &              i, 'Turning OFF sponge on tracer ', i,              &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerSrc(i,ng)) THEN
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc', i,      &
     &              'Turning ON  point sources/Sink on tracer ', i,     &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerSrc(i,ng), 'LtracerSrc', i,      &
     &              'Turning OFF point sources/Sink on tracer ', i,     &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LtracerCLM(i,ng)) THEN
                WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning ON  processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LtracerCLM(i,ng), 'LtracerCLM', i,      &
     &              'Turning OFF processing of climatology tracer ', i, &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            DO itrc=1,NST
              i=idsed(itrc)
              IF (LnudgeTCLM(i,ng)) THEN
                WRITE (out,150) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning ON  nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              ELSE
                WRITE (out,150) LnudgeTCLM(i,ng), 'LnudgeTCLM', i,      &
     &              'Turning OFF nudging of climatology tracer ', i,    &
     &              TRIM(Vname(1,idTvar(i)))
              END IF
            END DO
            IF ((nHIS(ng).gt.0).and.ANY(Hout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NST
                i=idTvar(idsed(itrc))
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idTvar)',                                     &
     &              'Write out sediment', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idfrac(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idfrac)',                                     &
     &              'Write out bed fraction, sediment ', itrc,          &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idBmas(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idmass)',                                     &
     &              'Write out mass, sediment ', itrc,                  &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idUbld(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idUbld)',                                     &
     &              'Write out bed load at U-points, sediment ', itrc,  &
     &              TRIM(Vname(1,i))
                i=idVbld(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idVbld)',                                     &
     &              'Write out bed load at V-points, sediment ', itrc,  &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBEDP
                i=idSbed(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &              'Hout(idSbed)',                                     &
     &              'Write out BED property ', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBOTP
                i=idBott(itrc)
                IF (Hout(i,ng)) WRITE (out,160) Hout(i,ng),             &
     &             'Hout(idBott)',                                      &
     &             'Write out BOTTOM property ', itrc, TRIM(Vname(1,i))
              END DO
            END IF
            IF ((nQCK(ng).gt.0).and.ANY(Qout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NST
                i=idTvar(idsed(itrc))
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idTvar)',                                     &
     &              'Write out sediment', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idsurT(idsed(itrc))
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idTvar)',                                     &
     &              'Write out surface sediment', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idfrac(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idfrac)',                                     &
     &              'Write out bed fraction, sediment ', itrc,          &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idBmas(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idfrac)',                                     &
     &              'Write out mass, sediment ', itrc,                  &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idUbld(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idUbld)',                                     &
     &              'Write out bed load at U-points, sediment ', itrc,  &
     &              TRIM(Vname(1,i))
                i=idVbld(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idVbld)',                                     &
     &              'Write out bed load at V-points, sediment ', itrc,  &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBEDP
                i=idSbed(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &              'Qout(idSbed)',                                     &
     &              'Write out BED property ', itrc, TRIM(Vname(1,i))
              END DO
              DO itrc=1,MBOTP
                i=idBott(itrc)
                IF (Qout(i,ng)) WRITE (out,160) Qout(i,ng),             &
     &             'Qout(idBott)',                                      &
     &             'Write out BOTTOM property ', itrc, TRIM(Vname(1,i))
              END DO
            END IF
            IF ((nAVG(ng).gt.0).and.ANY(Aout(:,ng))) THEN
              WRITE (out,'(1x)')
              DO itrc=1,NST
                i=idTvar(idsed(itrc))
                IF (Aout(i,ng)) WRITE (out,160) Aout(i,ng),             &
     &              'Aout(idTvar)',                                     &
     &              'Write out averaged sediment', itrc,                &
     &              TRIM(Vname(1,i))
              END DO
              DO itrc=1,NST
                i=idsed(itrc)
                IF (Aout(idTTav(i),ng)) WRITE (out,160)                 &
     &              Aout(idTTav(i),ng), 'Aout(idTTav)',                 &
     &              'Write out averaged <t*t> for tracer ', i,          &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NST
                i=idsed(itrc)
                IF (Aout(idUTav(i),ng)) WRITE (out,160)                 &
     &              Aout(idUTav(i),ng), 'Aout(idUTav)',                 &
     &              'Write out averaged <u*t> for tracer ', i,          &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NST
                i=idsed(itrc)
                IF (Aout(idVTav(i),ng)) WRITE (out,160)                 &
     &              Aout(idVTav(i),ng), 'Aout(idVTav)',                 &
     &              'Write out averaged <v*t> for tracer ', i,          &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NST
                i=idsed(itrc)
                IF (Aout(iHUTav(i),ng)) WRITE (out,160)                 &
     &              Aout(iHUTav(i),ng), 'Aout(iHUTav)',                 &
     &              'Write out averaged <Huon*t> for tracer ', i,       &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NST
                i=idsed(itrc)
                IF (Aout(iHVTav(i),ng)) WRITE (out,160)                 &
     &              Aout(iHVTav(i),ng), 'Aout(iHVTav)',                 &
     &              'Write out averaged <Hvom*t> for tracer ', i,       &
     &              TRIM(Vname(1,idTvar(i)))
              END DO
              DO itrc=1,NST
                i=idUbld(itrc)
                IF (Aout(i,ng)) WRITE (out,160) Aout(i,ng),             &
     &              'Aout(idUbld)',                                     &
     &              'Write out U-bedload, sediment ', itrc,             &
     &              TRIM(Vname(1,i))
                i=idVbld(itrc)
                IF (Aout(i,ng)) WRITE (out,160) Aout(i,ng),             &
     &              'Aout(idVbld)',                                     &
     &              'Write out V-bedload, sediment ', itrc,             &
     &              TRIM(Vname(1,i))
              END DO
            END IF
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant input parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NST
          Sd50(i,ng)=Sd50(i,ng)*0.001_r8
          Wsed(i,ng)=Wsed(i,ng)*0.001_r8
          tau_ce(i,ng)=tau_ce(i,ng)/rho0
          tau_cd(i,ng)=tau_cd(i,ng)/rho0
          nl_tnu4(idsed(i),ng)=SQRT(ABS(nl_tnu4(idsed(i),ng)))
          IF (Tnudg(idsed(i),ng).gt.0.0_r8) THEN
            Tnudg(idsed(i),ng)=1.0_r8/(Tnudg(idsed(i),ng)*86400.0_r8)
          ELSE
            Tnudg(idsed(i),ng)=0.0_r8
          END IF
        END DO
      END DO
  30  FORMAT (/,' READ_SedPar - variable info not yet loaded, ', a)
  40  FORMAT (/,' READ_SedPar - Error while processing line: ',/,a)
  50  FORMAT (/,/,' Sediment Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/)
  60  FORMAT (/,1x,'Size',5x,'Sd50',8x,'Csed',8x,'Srho',8x,'Wsed',      &
     &        8x,'Erate',7x,'poros',/,1x,'Class',4x,'(mm)',7x,          &
     &        '(kg/m3)',5x,'(kg/m3)',5x,'(mm/s)',5x,'(kg/m2/s)',4x,     &
     &        '(nondim)',/)
  70  FORMAT (2x,i2,2x,6(1x,1p,e11.4))
  80  FORMAT (/,9x,'tau_ce',6x,'tau_cd',6x,'nl_tnu2',5x,'nl_tnu4',5x,   &
     &        'Akt_bak',6x,'Tnudg',/,9x,'(N/m2)',6x,'(N/m2)',6x,        &
     &        '(m2/s)',6x,'(m4/s)',7x,'(m2/s)',6x,'(day)',/)
  90  FORMAT (/,9x,'morph_fac',/,9x,'(nondim)',/)
 100  FORMAT (/,' New bed layer formed when deposition exceeds ',e12.5, &
     &        ' (m).')
 110  FORMAT (' Two first layers are combined when 2nd layer smaller ', &
     &         'than ',e12.5,' (m).')
 120  FORMAT (' Rate coefficient for bed load transport = ',e12.5,/)
 130  FORMAT (' Transition for mixed sediment =',e12.5,/)
 140  FORMAT (' Transition for cohesive sediment =',e12.5,/)
 150  FORMAT (10x,l1,2x,a,'(',i2.2,')',t30,a,i2.2,':',1x,a)
 160  FORMAT (10x,l1,2x,a,t29,a,i2.2,':',1x,a)
 170  FORMAT (/,9x,'sed_rxn',/,9x,'(1/d)',/)
      RETURN
      END SUBROUTINE read_SedPar
