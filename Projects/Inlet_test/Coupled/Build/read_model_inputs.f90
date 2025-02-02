      SUBROUTINE read_model_inputs
!
!=======================================================================
!                                                                      !
!  This routine reads in model input parameters of dt and              !
!  number of grids for each model.                                     !
!                                                                      !
!=======================================================================
!
      USE mct_coupler_params
      USE mod_iounits
      USE swan_iounits
      implicit none
!
      include 'mpif.h'
!
!  Imported variable declarations.
!
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, isval, iw, ia, inp, out, status
      integer :: MyRank, MyError, MyMaster, DT, num, den
      integer :: cdecode_line, cload_i, cload_r, indx, indx2, test
      integer :: Ivalue
      real(m8), dimension(100) :: Rval
      real(m8) :: FAC
      character (len=1 ), parameter :: blank = ' '
      character (len=1 ) :: KEY
      character (len=40) :: KeyWord
      character (len=160) :: line
      character (len=160) :: aline
      character (len=160) :: saveline1, saveline2, saveline3
      character (len=160), dimension(100) :: Cval
!
      inp=1
      out=stdout
      MyMaster=0
      CALL mpi_comm_rank (MPI_COMM_WORLD, MyRank, MyError)
!
!     Read ROMS input file
!
      OPEN (inp, FILE=TRIM(Iname), FORM='formatted', STATUS='old',      &
     &      ERR=10)
      GO TO 30
 10   WRITE (out,20) Iname
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Iname)
!     exit_flag=4
      RETURN
 20   FORMAT (/,' READ MODEL INPUTS - Unable to open roms input file.', &
     &        /,a80)
 30   CONTINUE
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=15,END=40) line
        status=cdecode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Ngrids')
              Npts=cload_i(Nval, Rval, 1, Ivalue)
              Nocn_grids=Ivalue
              IF (Nocn_grids.le.0) THEN
                IF (MyRank.eq.MyMaster) WRITE (out,290)'Ngrids', Ngrids,&
     &            'Ngrids must be greater than zero.'
!                exit_flag=5
                RETURN
              END IF
              allocate (dtocn(Nocn_grids))
            CASE ('DT')
              Npts=cload_r(Nval, Rval, Nocn_grids, dtocn)
          END SELECT
!         IF (exit_flag.ne.NoError) RETURN
        END IF
      END DO
 15   IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 40   CLOSE (inp)
 290  FORMAT (/,'read model inputs - Invalid dimension parameter,',a,i4,&
     &        /,15x,a)
!
!     Read SWAN input file
!
      iw=1
      OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',  &
     &      ERR=110)
      GO TO 130
 110  WRITE (out,120) Wname(iw)
      IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,         &
     &                        TRIM(Wname(iw))
!     exit_flag=4
      RETURN
 120  FORMAT (/,' READ MODEL INPUTS - Unable to open swan input file.', &
     &        /,a80)
 130  CONTINUE
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=116,END=140) line
        aline=ADJUSTL(line)
!
!  dont read comment lines
!
        isval=INDEX(aline,'&')
        IF(isval.eq.0) THEN
!
!  here we look for number of swan grids
!
          IF(aline(1:7).eq.'NSGRIDS') THEN
            read(aline(8:12),'(i5)') Nwav_grids
            allocate (dtwav(Nwav_grids))
          END IF
!
!  here we look for model time step
!
          IF(aline(1:4).eq.'COMP') THEN
            DO i=1,3
              indx=INDEX(aline,blank)
              aline=aline(indx+1:LEN(aline))
            END DO
            DO i=1,1
              indx=INDEX(aline,blank)
              read(aline(1:indx),'(i10)') DT
              dtwav(iw)=REAL(DT,m8)
              aline=aline(indx+1:LEN(aline))
            END DO
            READ(aline,'(a1)') KEY
            IF (KEY.eq.'D') THEN
              FAC = 24.0_m8*3600.0_m8
            ELSE IF (KEY.eq.'H') THEN
              FAC = 3600.0_m8
            ELSE IF (KEY.eq.'M') THEN
              FAC = 60.0_m8
            ELSE
              FAC = 1.0_m8
            ENDIF
            dtwav(iw)=dtwav(iw)*FAC
          END IF
        END IF
      END DO
 116  IF (MyRank.eq.MyMaster) WRITE (out,60) line
!     exit_flag=4
      RETURN
 140  CLOSE (inp)
!
      DO iw=2,Nwav_grids
        OPEN (inp, FILE=TRIM(Wname(iw)), FORM='formatted', STATUS='old',&
     &        ERR=110)
        GO TO 135
        IF (MyRank.eq.MyMaster) WRITE(out,*) 'MyRank = ', MyRank,       &
     &                          TRIM(Wname(iw))
!       exit_flag=4
        RETURN
 135    CONTINUE
        DO WHILE (.TRUE.)
          READ (inp,'(a)',ERR=115,END=160) line
          aline=ADJUSTL(line)
          IF(aline(1:4).eq.'COMP') THEN
            DO i=1,3
              indx=INDEX(aline,blank)
              aline=aline(indx+1:LEN(aline))
            END DO
            DO i=1,1
              indx=INDEX(aline,blank)
              read(aline(1:indx),'(i10)') DT
              dtwav(iw)=REAL(DT,m8)
              aline=aline(indx+1:LEN(aline))
            END DO
            READ(aline,'(a1)') KEY
            IF (KEY.eq.'D') THEN
              FAC = 24.0_m8*3600.0_m8
            ELSE IF (KEY.eq.'H') THEN
              FAC = 3600.0_m8
            ELSE IF (KEY.eq.'M') THEN
              FAC = 60.0_m8
            ELSE
              FAC = 1.0_m8
            ENDIF
            dtwav(iw)=dtwav(iw)*FAC
          END IF
        END DO
 115    IF (MyRank.eq.MyMaster) WRITE (out,60) line
!       exit_flag=4
        RETURN
 160    CLOSE (inp)
      END DO
  60  FORMAT (/,'read model inputs - Error while processing line: ',/,a)
      RETURN
      END SUBROUTINE read_model_inputs
