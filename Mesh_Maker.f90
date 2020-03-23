! This file contains subroutines for mesh maker

!************************************************************************
      SUBROUTINE MESHM
!************************************************************************
!
!          .....................................................
!          .                                                   .
!          .  TOUGH2, MODULE MESHM, VERSION 1.0, MARCH 1991    .
!          .....................................................
!
!-----THIS IS THE EXECUTIVE ROUTINE FOR THE *MESHMAKER* MODULE.
!
!     IT IS CALLED WHEN DATA BLOCK "MESHM" IS PRESENT IN THE INPUT
!     FILE, AND WILL THEN PROCEED TO PERFORM INTERNAL MESH GENERATION.
!
      USE SVZ
      USE MMC
      USE MMD
      USE MESH_ARRAY   !数组
      USE FILEUNITS   !文件的编号

      IMPLICIT NONE
      
      INTEGER :: IDUM1,K,IOS
      CHARACTER VER*5,WORD*5
      LOGICAL EX
      DIMENSION VER(6)
      SAVE VER
      DATA VER/'RZ2D ','RZ2DL','XYZ  ','PATTY','MINC ','     '/
!
      ALLOCATE(e(6000000))      !可分配数组
      
      BACKSPACE(NUINPUT)
      READ(NUINPUT,'(5X,I1)',IOSTAT=IOS) IDUM1
      IF (IOS.NE.0) IDUM1=0
      IF (IDUM1.EQ.0.AND.MOP2(2).EQ.0) IELL=5
      IF (IDUM1.EQ.0.AND.MOP2(2).NE.0) IELL=MOP2(2)
      IF (IDUM1.GE.5.AND.IDUM1.LE.9.AND.MOP2(2).EQ.0) THEN
        IELL=IDUM1
        MOP2(2)=IELL
      ENDIF
!
  100 CONTINUE
      READ(NUINPUT,1) WORD
    1 FORMAT(A5)
      DO 2 K=1,6
      IF(WORD.EQ.VER(K)) GOTO(11,11,12,13,14,101),K
    2 CONTINUE
!
      write(36,4) WORD
    4 FORMAT(' MESHMAKER HAS READ UNKNOWN BLOCK LABEL ',A5,'  ---  ', &
      'STOP EXECUTION')
      STOP
!
   11 CONTINUE
      write(36,21) IELL
   21 format(/,' ',131('*')/' *',20X,'MESHMAKER - RZ2D: GENERATE 2-D R-Z', &
      ' MESH WITH',I2,'-CHARACTER ELEMENT NAMES',38X,'*'/' ',131('*'))
      CALL RZ2D
      CALL WRZ2D(K)
      GOTO 100
!
   12 CONTINUE
      write(36,22) IELL
   22 format(/,' ',131('*')/' *',26X,'MESHMAKER - XYZ: GENERATE', &
      ' CARTESIAN MESH WITH ',I1,'-CHARACTER ELEMENT NAMES',32X,'*'/' ',131('*'))
      CALL GXYZ
      GOTO 100
!
   13 CONTINUE
      GOTO 100
!
   14 CONTINUE
      MESHMC=MESHMC+1
      write(36,24) 
   24 format(/,' ',131('*')/' *',18X,'MESHMAKER - MINC: GENERATE MULTIPLE' &
      ,' INTERACTING CONTINUA MESH FOR FRACTURED MEDIUM',29X,'*'/ &
      ' ',131('*')/)
!
      INQUIRE(FILE='MINC',EXIST=EX)
      IF(EX) GOTO 112
      write(36, 113)
  113 FORMAT(' FILE *MINC* DOES NOT EXIST --- OPEN AS A NEW FILE')
      OPEN(10,FILE='MINC',STATUS='NEW')
      GOTO 120
!
  112 write(36, 114)
  114 FORMAT(' FILE *MINC* EXISTS --- OPEN AS AN OLD FILE')
      OPEN(10,FILE='MINC',STATUS='OLD')
!
  120 CONTINUE
!
      CALL MINC
      GOTO 100
!
  101 CONTINUE
      write(36, 102)
  102 FORMAT(/' ',131('*')//' MESH GENERATION COMPLETE --- EXIT FROM ', &
      'MODULE *MESHMAKER*'/)

      DEALLOCATE(e)

      RETURN
      END

!************************************************************************      
      SUBROUTINE RZ2D
!************************************************************************      
      USE SVZ
      USE MMD
      USE FILEUNITS
      IMPLICIT NONE
      
      REAL(8) :: S,X,DR,RLOG,AM,AN,XM,XL,XR,SMID,RAT,ZN,ZTOT,REL
      INTEGER :: NRAD,NEQU,NLOG,NRADC,NEQUC,NLOGC,NINT,K,NINT1,I,IX
      CHARACTER SUB*5,WORD*5,CBLANK*4
!
      DIMENSION SUB(4)
      SAVE SUB
!
      DATA SUB/'RADII','EQUID','LOGAR','LAYER'/
      DATA CBLANK/'    '/
!
      S(X)=X*(X**NLOG-1.D0)/(X-1.D0)
!
      NRAD=0
      NEQU=0
      NLOG=0
      NRADC=0
      NEQUC=0
      NLOGC=0
      NINT=0
!
  204 READ(NUINPUT,200) WORD
  200 FORMAT(A5)
!
      DO 201 K=1,4
      IF(WORD.EQ.SUB(K)) GOTO 202
  201 CONTINUE
!
      write(36, 203) WORD
  203 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL ',A5,'  ---  STOP', &
      ' EXECUTION')
      STOP
!
  202 GOTO(211,212,213,214),K
!
  211 CONTINUE
!-----READ A SET OF INTERFACE RADII.
      READ(NUINPUT,220) NRAD  !NRAD要读取的半径的个数
  220 FORMAT(I5,5X,7E10.4)
!
      NINT1=NINT+1    !内存段的起点
      NINT=NINT+NRAD   !内存段的重点
      NRADC=NRADC+NRAD   !首次读入的半径数
      CALL ALLOCATE_RC(NINT1-1,NINT) !动态分配内存
      READ(NUINPUT,221) (RC(I),I=NINT1,NINT)
  221 FORMAT(8E10.4)
      GOTO 204
!
  212 CONTINUE
!-----READ DATA FOR GENERATING SOME EQUIDISTANT INTERFACES.
!
      READ(NUINPUT,220) NEQU,DR
      NINT1=NINT+1
      NINT=NINT+NEQU
      NEQUC=NEQUC+NEQU
      CALL ALLOCATE_RC(NINT1-1,NINT)
      DO 222 I=NINT1,NINT
  222 RC(I)=RC(I-1)+DR  !半径
      GOTO 204
!
  213 CONTINUE
!-----READ DATA FOR GENERATING INTERFACES WITH LOGARITHMICALLY
!     INCREASING DISTANCE.
!
      READ(NUINPUT,220) NLOG,RLOG,DR
      NLOGC=NLOGC+NLOG
!      IF(DR.EQ.0.D0) DR=RC(NINT)-RC(NINT-1)
       IF(abs(DR).le.0.0) DR=RC(NINT)-RC(NINT-1)
!
      write(36, 109)
  109 FORMAT(/' ***** FIND APPROPRIATE FACTOR FOR', &
      ' INCREASING RADIAL DISTANCES *****')
!
      IX=0
      AM=(RLOG-RC(NINT))/DR
      AN=FLOAT(NLOG)
      XM=AM/AN
      IF(abs(XM-1.0).le.0.0) GOTO 40
      IF(AM.GT.AN) THEN  !半径增大
      XL=1.D0
      XR=(AM+1.D0)**(1.D0/AN)
      XR=MIN(XR,XM)
      ELSE   !半径减小
      XL=AM/(AM+1.D0)
      XL=MAX(XL,XM)
      XR=1.D0
      ENDIF
!     FIND PROPER INCREMENT FACTOR THROUGH ITERATIVE BISECTING求得增量因子
   41 XM=0.5D0*(XL+XR)
      IF(XR.GT.2.D0*XL) XM=SQRT(XL*XR)
      SMID=S(XM)
      IX=IX+1
      IF(IX.GE.1000) GOTO 43
      IF(XR-XL.LE.XL*1.D-10) GOTO 40
      IF(SMID.LE.AM) XL=XM
      IF(SMID.GE.AM) XR=XM
      GOTO 41
!
   40 RAT=XM
!
      write(36, 30) IX,RAT
   30 FORMAT(6X,' AFTER',I4,' ITERATIONS FOUND INCREMENT FACTOR', &
      ' RAT =',E14.8)
!
      CALL ALLOCATE_RC(NINT,NINT+NLOG)   !记录半径
      DO 223 I=1,NLOG
      DR=DR*RAT
  223 RC(NINT+I)=RC(NINT+I-1)+DR
!
      NINT=NINT+NLOG  
      GOTO 204  !半径总数
!
  214 CONTINUE
!-----COME HERE AFTER ALL NINT INTERFACES ARE ASSIGNED, AND READ
!     DATA ON LAYER THICKNESSES.
!
      READ(NUINPUT,220) NLAY
      ALLOCATE(HHHH(NLAY))        !分配层数及读取层厚度
      READ(NUINPUT,221) (HHHH(I),I=1,NLAY)
      DO 219 I=2,NLAY
      IF(abs(HHHH(I)).le.0.0) HHHH(I)=HHHH(I-1)
  219 IF(HHHH(I).LT.0.D0) HHHH(I)=0.D0
      IF(HHHH(1).LT.0.D0) HHHH(1)=0.D0
!
      WRITE(36,216)
  216 FORMAT(//' ***** DATA ON GRID LAYERS *****'//                     &
     &       '    NO.     THICKNESS     NODAL DEPTH    TOTAL DEPTH'/)   
      ZN=0.D0 
      DO 218 I=1,NLAY 
      ZN=ZN+HHHH(I)/2.D0
      IF(I.GT.1) ZN=ZN+HHHH(I-1)/2.D0 
      ZTOT=ZN+HHHH(I)/2.D0 
      WRITE(36,217) I,HHHH(I),ZN,ZTOT 
  218 END DO 
  217 FORMAT(3X,I4,3(3X,E12.6)) 
!  
!-----NOW ASSIGN ALL GEOMETRY DATA FOR A LAYER OF UNIT THICKNESS.
!
      NELEMT=NINT-1
      ALLOCATE(A(NELEMT))         !可分配数组
      ALLOCATE(DDDD(NELEMT))
      ALLOCATE(V(NELEMT))
      DO 15 I=1,NELEMT
      A(I)=2.D0*PI*RC(I+1)
      DDDD(I)=(RC(I+1)-RC(I))/2.D0
   15 v(i)=pi*(rc(i+1)+rc(i))*(rc(i+1)-rc(i))
!
!=====PRINT ALL GEOMETRIC QUANTITIES=====
!
      write(36, 19) CBLANK(:IELL-5)
   19 FORMAT(//' * * * * * M E S H  G E O M E T R Y * * * * *',10X, &
      'VOLUME AND AREA ARE GIVEN FOR HEIGHT = 1 METER'// &
      ' ELEM',A,'     REL           RCON           D              V      ', &
      '       A'/)
!
      DO 17 I=1,NELEMT
      REL=RC(I+1)-DDDD(I)
   17 write(36,18) I,REL,RC(I+1),DDDD(I),V(I),A(I)
   18 FORMAT(' ',I4,5(2X,E12.5))
!
      REWIND 4
      write(36, 20)
   20 FORMAT(/' WRITE FILE *MESH* FROM INPUT DATA')
      WRITE(4,3) MOP2(2),NELEMT,NRADC,NEQUC,NLOGC,RC(1),RC(NELEMT+1)
    3 FORMAT('ELEME',I1,' -- ',4I5,F10.5,F10.3)
!
      RETURN
!
   43 CONTINUE
      write(36,42) IX,XL,XR,XM,SMID,AM
   42 FORMAT(' AT ITERATION',I3,'  XL=',E20.14,'  XR=', &
      E20.14,'  XM=',E12.6,'  SMID=',E12.6,'  AM=',E12.6)
      write(36, 44)
   44 FORMAT(' CONVERGENCE FAILURE IN SUBROUTINE RZ2D'/ &
      ' STOP EXECUTION')
      STOP
      END

!************************************************************************
      SUBROUTINE WRZ2D(K)
!************************************************************************
      USE MADIM
      USE MESH_ARRAY
      USE MMD
      USE SVZ
      IMPLICIT NONE
      LOGICAL     DECIMAL
      INTEGER     INDEL,NSKIP1,NSKIP2,IELLM9,I,IDIG,LDIG,IDIGTEN,K,MI,II, &
                  J,MJ,IOS,MI1,II1,J1,MJ1,L
      REAL*8      REL,ZJ,AIJ,VIJ,AI,D1,D2,AII
      CHARACTER   NII*1,NII1*1,NOV*1,NOV1*1
      CHARACTER*4 CBLANK,CFORMAT
      CHARACTER*9 CI,CI1
      CHARACTER   CFELEM1*44,CFELEM2*45,CFCONN1*28,CFCONN2*29, &
                  CFCONN3*34,CFCONN4*33
!
      DATA CBLANK/'    '/
      DATA CFORMAT/'(I?)'/
      DATA CFELEM1/'(Ix,Iy,6X,A,A5,2E10.4,10X,E10.4 ,10X,E10.4 )'/
      DATA CFELEM2/'(3A1,Ix,6X,A,A5,2E10.4,10X,E10.4 ,10X,E10.4 )'/
      DATA CFCONN1/'(2(Ix,Iy),11X,2A,1H1,3E10.4)'/
      DATA CFCONN2/'(2(3A1,Ix),11X,2A,1H1,3E10.4)'/
      DATA CFCONN3/'(2(3A1,Ix),11X,2A,1H3,3E10.4,2H1.)'/
      DATA CFCONN4/'(2(Ix,Iy),11X,2A,1H3,3E10.4,2H1.)'/
!
      NSKIP1=1
      NSKIP2=2
      IELLM9=9-IELL      
      IF(MOP2(2).NE.0) THEN
      DO 1000 I=6,1,-1
        IF (10**I.GT.NELEMT) IDIG=I
        IF (10**I.GT.NLAY) LDIG=I
 1000 CONTINUE
      IF (IDIG+LDIG.LE.IELL) THEN
         DECIMAL=.TRUE.
         LDIG=IELL-IDIG
         IDIGTEN=1
      ELSE
         DECIMAL=.FALSE.
         IDIG=IELL-3
         IDIGTEN=10**IDIG
      ENDIF
      WRITE(CFELEM1(3:3),'(I1)') LDIG
      WRITE(CFELEM1(6:6),'(I1)') IDIG
      WRITE(CFELEM2(7:7),'(I1)') IDIG
      WRITE(CFCONN1(5:5),'(I1)') LDIG
      WRITE(CFCONN1(8:8),'(I1)') IDIG
      WRITE(CFCONN2(9:9),'(I1)') IDIG
      WRITE(CFCONN3(9:9),'(I1)') IDIG
      WRITE(CFCONN4(5:5),'(I1)') LDIG
      WRITE(CFCONN4(8:8),'(I1)') IDIG
      ENDIF
!
      IF (MOP2(22).GT.0) THEN
        CFELEM1(28:32)='20.14'
        CFELEM1(34:35)='20'
        CFELEM1(39:43)='20.14'      
        CFELEM2(29:33)='20.14'
        CFELEM2(35:36)='20'
        CFELEM2(40:44)='20.14'      
      ENDIF
!
      NII=' '
      IF(K.EQ.2) GOTO 1
!
!-----COME HERE TO GENERATE MESH BY VERTICAL COLUMNS
      IF(MOP2(2).EQ.0.AND.MOP2(27).EQ.0) THEN
!     TOUGH2 V2.1      
      DO 4 I=NSKIP1,NELEMT
        MI=MOD(I,100)
        II=I/100
        IF(II.GE.1) NII=NA(II)      !NA是什么？？？
        REL=RC(I+1)-DDDD(I)
        ZJ=-HHHH(1)/2.D0
      DO 4 J=1,NLAY
        MJ=MOD(J-1,35)+1
!     MJ RUNS FROM 1 TO 35 FOR CONSECUTIVE LAYERS.
        NOV=NA((J-1)/35+10)
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 35 LAYERS.
        IF(J.GT.1) ZJ=ZJ-HHHH(J)/2.D0-HHHH(J-1)/2.D0
        AIJ=0.D0
        IF(J.EQ.1) AIJ=AIJ+V(I)
        IF(J.EQ.NLAY) AIJ=AIJ+V(I)
        VIJ=V(I)*HHHH(J)
      
        IF(MOP2(22).EQ.0) THEN  
            WRITE(4,5) NOV,NA(MJ),NII,MI,VIJ,AIJ,REL,ZJ
        ELSE IF(MOP2(22).GT.0) THEN
            WRITE(4,55) NOV,NA(MJ),NII,MI,VIJ,AIJ,REL,ZJ
        ENDIF
!-----CONVENTION FOR ELEMENT NAMES.
!     NOV = A, B, C, .... LABELS 1ST, 2ND, 3RD, ... SET OF 35 LAYERS.
!     NA = 1, 2, 3, ... 9, A, B, C, ... LABELS LAYERS.
!     NII = 1, 2, 3, ... AFTER 100, 200, 300, ... RADIAL ELEMENTS.
!     MI = 1, 2, 3, ... 99; NUMBERS RADIAL ELEMENTS IN EXCESS OF
!          FULL HUNDREDS.
    4 CONTINUE
    5 FORMAT(3A1,I2,10X,'    1',2E10.4,10X,E10.4,10X,E10.4)
   55 FORMAT(3A1,I2,10X,'    1',2E10.4,10X,E20.14,20X,E20.14) 
!
      ELSE
!     FROM iTOUGH2  
      INDEL=0     
      WRITE(CFORMAT(3:3),'(I1)',IOSTAT=IOS) IELL
      DO I=NSKIP1,NELEMT
         IF (MOP2(27).EQ.1) THEN
            CONTINUE
         ELSE IF (.NOT.DECIMAL) THEN
            MI=MOD(I,IDIGTEN)
            II=I/IDIGTEN
            IF (II.GE.1) NII=NA(II)
         ENDIF
         REL=RC(I+1)-DDDD(I)
         ZJ=-HHHH(1)/2.D0
         DO 444 J=1,NLAY
            IF(J.GT.1) ZJ=ZJ-HHHH(J)/2.0D0-HHHH(J-1)/2.0D0
            AIJ=0.0D0
            IF(J.EQ.1) AIJ=AIJ+V(I)
            IF(J.EQ.NLAY) AIJ=AIJ+V(I)
            VIJ=V(I)*HHHH(J)
            IF (MOP2(27).EQ.1) THEN
               INDEL=INDEL+1
               WRITE(CI(1:IELL),CFORMAT) INDEL
               WRITE(4,7000) CI(:IELL),CBLANK(:IELLM9),'    1', &
                             VIJ,AIJ,REL,ZJ
 7000 FORMAT(A,6X,A,A5,2E10.4,2(10X,E10.4))          
            ELSE IF (DECIMAL) THEN
               WRITE(4,CFELEM1) J,I,CBLANK(:IELLM9),'    1', &
                                VIJ,AIJ,REL,ZJ
            ELSE
               MJ=MOD(J-1,61)+1
!     MJ RUNS FROM 1 TO 61 FOR CONSECUTIVE LAYERS.
               NOV=NA((J-1)/61+10)
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 61 LAYERS.
               WRITE(4,CFELEM2) NOV,NA(MJ),NII,MI,CBLANK(:IELLM9), &
                                 '    1',VIJ,AIJ,REL,ZJ
            ENDIF
  444    CONTINUE
      ENDDO
      ENDIF
!
      GOTO 2
!
    1 CONTINUE
!-----COME HERE TO GENERATE MESH BY HORIZONTAL LAYERS
      ZJ=-HHHH(1)/2.D0
      IF(MOP2(2).EQ.0.AND.MOP2(27).EQ.0) THEN
!     TOUGH2 V2.1      
      DO 3 J=1,NLAY
      MJ=MOD(J-1,35)+1
!     MJ RUNS FROM 1 TO 35 FOR CONSECUTIVE LAYERS.
      NOV=NA((J-1)/35+10)
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 35 LAYERS.
      IF(J.GT.1) ZJ=ZJ-HHHH(J)/2.D0-HHHH(J-1)/2.D0
!.....1-12-95
      nii=' '
      DO 3 I=NSKIP1,NELEMT
      MI=MOD(I,100)
      II=I/100
      IF(II.GE.1) NII=NA(II)
      REL=RC(I+1)-DDDD(I)
      AIJ=0.D0
      IF(J.EQ.1) AIJ=AIJ+V(I)
      IF(J.EQ.NLAY) AIJ=AIJ+V(I)
      VIJ=V(I)*HHHH(J)
      IF(MOP2(22).EQ.0) THEN
        WRITE(4,5) NOV,NA(MJ),NII,MI,VIJ,AIJ,REL,ZJ
      ELSE IF(MOP2(22).GT.0) THEN
        WRITE(4,55) NOV,NA(MJ),NII,MI,VIJ,AIJ,REL,ZJ
      ENDIF
    3 CONTINUE
!      
      ELSE
!     FROM iTOUGH2       
      DO J=1,NLAY
         IF (.NOT.DECIMAL) THEN
            MJ=MOD(J-1,61)+1
!     MJ RUNS FROM 1 TO 61 FOR CONSECUTIVE LAYERS.
            NOV=NA((J-1)/61+10)
         ENDIF
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 61 LAYERS.
         IF (J.GT.1) ZJ=ZJ-HHHH(J)/2.0D0-HHHH(J-1)/2.0D0
         DO 33 I=NSKIP1,NELEMT
            REL=RC(I+1)-DDDD(I)
            AIJ=0.0D0
            IF(J.EQ.1) AIJ=AIJ+V(I)
            IF(J.EQ.NLAY) AIJ=AIJ+V(I)
            VIJ=V(I)*HHHH(J)
            IF (MOP2(27).EQ.1) THEN
               INDEL=INDEL+1
               WRITE(CI(1:IELL),CFORMAT) INDEL
               WRITE(4,7000) CI(:IELL),CBLANK(:IELLM9), &
                             '    1',VIJ,AIJ,REL,ZJ
            ELSE IF (DECIMAL) THEN
               WRITE(4,CFELEM1) J,I,CBLANK(:IELLM9), &
                                '    1',VIJ,AIJ,REL,ZJ
            ELSE
               MI=MOD(I,IDIGTEN)
               II=I/IDIGTEN
               IF(II.GE.1) NII=NA(II)
               WRITE(4,CFELEM2) NOV,NA(MJ),NII,MI,CBLANK(:IELLM9), &
                                '    1',VIJ,AIJ,REL,ZJ
            ENDIF
   33    CONTINUE
      ENDDO
      ENDIF
!
    2 CONTINUE
!
      WRITE(4,6)
    6 FORMAT('     ')
      WRITE(4,7)
    7 FORMAT('CONNE')
!
!-----ASSIGN HORIZONTAL CONNECTIONS.
!
      NII=' '
      NII1=' '
      IF(MOP2(2).EQ.0.AND.MOP2(27).EQ.0) THEN
!     TOUGH2 V2.1      
      DO 8 I=NSKIP2,NELEMT
      MI=MOD(I,100)
      II=I/100
      IF(II.GE.1) NII=NA(II)
      MI1=MOD(I-1,100)
      II1=(I-1)/100
      IF(II1.GE.1) NII1=NA(II1)
      DO 8 J=1,NLAY
      MJ=MOD(J-1,35)+1
!     MJ RUNS FROM 1 TO 35 FOR CONSECUTIVE LAYERS.
      NOV=NA((J-1)/35+10)
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 35 LAYERS.
!     DO 8 J=1,1
      AIJ=A(I-1)*HHHH(J)
    8 WRITE(4,9) NOV,NA(MJ),NII1,MI1,NOV,NA(MJ),NII,MI,DDDD(I-1),DDDD(I), &
      AIJ
    9 FORMAT(3A1,I2,3A1,I2,19X,1H1,3E10.4)
      ELSE
!     FROM iTOUGH2       
      DO I=NSKIP2,NELEMT      
         IF (.NOT.DECIMAL) THEN
            MI=MOD(I,IDIGTEN)
            II=I/IDIGTEN
            IF(II.GE.1) NII=NA(II)
            MI1=MOD(I-1,IDIGTEN)
            II1=(I-1)/IDIGTEN
            IF(II1.GE.1) NII1=NA(II1)
         ENDIF
         DO 88 J=1,NLAY
            AIJ=A(I-1)*HHHH(J)
            IF (MOP2(27).EQ.1) THEN
               WRITE(CI(1:IELL), CFORMAT) (I-1)*NLAY+J
               WRITE(CI1(1:IELL),CFORMAT) (I-2)*NLAY+J
               WRITE(4,914) CI1(:IELL),CI(:IELL),CBLANK(:IELLM9),   &
                            CBLANK(:IELLM9),DDDD(I-1),DDDD(I),AIJ
  914 FORMAT(2A,11X,2A,1H1,4E10.4)      
            ELSE IF (DECIMAL) THEN
               WRITE(4,CFCONN1) J,I-1,J,I,CBLANK(:IELLM9),          &
                                CBLANK(:IELLM9),DDDD(I-1),DDDD(I),AIJ
            ELSE
               MJ=MOD(J-1,61)+1
!     MJ RUNS FROM 1 TO 61 FOR CONSECUTIVE LAYERS.
               NOV=NA((J-1)/61+10)
!     NOV RUNS THROUGH A, B, C, ... FOR THE 1ST, 2ND, 3RD, ... SET
!     OF 61 LAYERS.
               WRITE(4,CFCONN2) NOV,NA(MJ),NII1,MI1,NOV,NA(MJ),NII,     &
                                MI,CBLANK(:IELLM9),CBLANK(:IELLM9),     &
                                DDDD(I-1),DDDD(I),AIJ
            ENDIF
   88    CONTINUE
      ENDDO
      ENDIF
!
      IF(NLAY.LE.1) GOTO 42
!
!-----NOW FOR THE VERTICAL CONNECTIONS.
!
      NII=' '
      IF(MOP2(2).EQ.0.AND.MOP2(27).EQ.0) THEN
!     TOUGH2 V2.1            
      DO 40 I=NSKIP1,NELEMT
      MI=MOD(I,100)
      II=I/100
      IF(II.GE.1) NII=NA(II)
      AI=V(I)
!
      DO 40 J=2,NLAY
      J1=J-1
      MJ=MOD(J-1,35)+1
      NOV=NA((J-1)/35+10)
      MJ1=MOD(J1-1,35)+1
      NOV1=NA((J1-1)/35+10)
      D1=HHHH(J1)/2.D0
      D2=HHHH(J)/2.D0
   40 WRITE(4,41) NOV1,NA(MJ1),NII,MI,NOV,NA(MJ),NII,MI,D1,D2,AI
   41 FORMAT(3A1,I2,3A1,I2,19X,1H3,3E10.4,2H1.)
      ELSE
!     FROM iTOUGH2             
      DO I=NSKIP1,NELEMT
         IF (.NOT.DECIMAL) THEN
            MI=MOD(I,IDIGTEN)
            II=I/IDIGTEN
            IF(II.GE.1) NII=NA(II)
         ENDIF
         AII=V(I)
!
         DO 44 J=2,NLAY
            J1=J-1
            MJ=MOD(J-1,61)+1
            NOV=NA((J-1)/61+10)
            MJ1=MOD(J1-1,61)+1
            NOV1=NA((J1-1)/61+10)
            D1=HHHH(J1)/2.0D0
            D2=HHHH(J)/2.0D0
            IF (MOP2(27).EQ.1) THEN
               WRITE(CI(1:IELL), CFORMAT) (I-1)*NLAY+J
               WRITE(CI1(1:IELL),CFORMAT) (I-1)*NLAY+J-1
               WRITE(4,915) CI1(:IELL),CI(:IELL),CBLANK(:IELLM9),   &
                            CBLANK(:IELLM9),D1,D2,AII
  915 FORMAT(2A,11X,2A,1H3,3E10.4,2H1.)
            ELSE IF (DECIMAL) THEN
               WRITE(4,CFCONN4) J1,I,J,I,CBLANK(:IELLM9),           &
                                CBLANK(:IELLM9),D1,D2,AII
            ELSE
               WRITE(4,CFCONN3) NOV1,NA(MJ1),NII,MI,NOV,NA(MJ),NII, &
                     MI,CBLANK(:IELLM9),CBLANK(:IELLM9),D1,D2,AII
            ENDIF
   44    CONTINUE
      ENDDO
      ENDIF
!
   42 CONTINUE
!
      WRITE(4,6)
!
      REWIND 4
      READ(4,*) 
   43 FORMAT(A)
      L=0
      DO 441 I=NSKIP1,NELEMT
      DO 441 J=1,NLAY
      L=L+1
  441 READ(4,43) E(L)
  443 FORMAT(A9)  
!
      CALL PRZ2D(K,NSKIP1,NELEMT,NLAY)
!
      RETURN
      END

!************************************************************************
      SUBROUTINE PRZ2D(KK,NSKIP1,NELEMT,NLAY)
!************************************************************************      
!
!-----MAKE STRUCTURED PRINTOUT OF R-Z GRID.
!
      USE MESH_ARRAY
      USE MADIM
      USE SVZ, ONLY:IELL
      IMPLICIT NONE
!
      INTEGER :: KK,NSKIP1,NELEMT,NLAY
      INTEGER :: NIK1,NIK2
      INTEGER :: I,K,NRZ
!     FIND LOCATION OF ELEMENT NAME WITH R-INDEX I, Z-INDEX K.
!     MESH ORGANIZED "BY COLUMNS" (KK = 1)
      NIK1(I,K)=(I-1)*NLAY+K
!     MESH ORGANIZED "BY LAYERS" (KK = 2)
      NIK2(I,K)=(K-1)*(NELEMT-NSKIP1+1)+I
!
      NRZ=NELEMT*NLAY
      write(36,31) NELEMT,NLAY,NRZ
   31 FORMAT(/' ',131('*')/' *',20X,'2-D R-Z MESH WITH NR*NLAY = ', &
      I4,' *',I4,'  =',I4,' GRID BLOCKS',52X,'*'/' ',131('*'))
!
      write(36,18) NLAY,NELEMT
   18 FORMAT(' *',120X,9X,'*'/ &
      ' *',20X,'THE MESH WILL BE PRINTED AS VERTICAL SLICES', &
      66X,'*'/' *',120X,9X,'*'/ &
      ' *',20X,'LAYERS GO FROM K = 1 TO K = NLAY =',I4,71X,'*'/ &
      ' *',120X,9X,'*'/ &
      ' *',20X,'RADIAL GRID BLOCKS GO IN COLUMNS FROM I = 1 TO', &
      ' I = NR =',I4,50X,'*'/ &
      ' *',120X,9X,'*'/' ',131('*')/)
!
      WRITE(36, 4) (I,I=1,NELEMT)
    4 FORMAT('   COLUMN I =',I9,9I10/(12X,10I10))
      WRITE(36,'(A)') ' LAYERS'
!
      DO 1 K=1,NLAY
      IF(KK.EQ.1) write(36,6) K,(E(NIK1(I,K))(:IELL),I=NSKIP1,NELEMT)
      IF(KK.EQ.2) write(36,6) K,(E(NIK2(I,K))(:IELL),I=NSKIP1,NELEMT)
    1 CONTINUE
!
    6 FORMAT(/,' K =',I6,2X,10(1X,A)/(12X,10(1X,A)))
!
      RETURN
      END

!************************************************************************
      SUBROUTINE GXYZ
!************************************************************************      
!
!-----PROGRAM FOR GENERATING 3-D X-Y-Z MESHES.
!
!     GRID BLOCK DIMENSIONS ARE DX, DY, AND DZ.
!
      USE MADIM
      USE MESH_ARRAY
      USE MMD
      USE SVZ
      USE FILEUNITS
      USE OUTPUT, ONLY:NNXYZ
      IMPLICIT NONE

      CHARACTER NW*2,NXYZ*2   
      CHARACTER DOM*5
      CHARACTER CFORMAT*4,CBLANK*10,CFELEM1*34,CFELEM2*33
      CHARACTER*9 CI,CI1,CJ,CJ1,CK,CK1

      INTEGER :: NX,NY,NZ,NO,N,NX1,NY1,NZ1,I,NN,MESHDIM,I1,ILOOP,   &
                 IDIG,JDIG,KDIG,IM,NOVI,IM1,NOVI1,J,J1,JM,NOVJ,     &
                 JM1,JM11,NOVJ1,NOVJ11,K,K1,KM,NOVK,KM1,NOVK1,NOV,  &
                 NOV1,NOV2,NOV3,NOV4,NOV5,L,NLOC,IC,N1,N2,N3,INDEL, &
                 IOS
      REAL(8) :: D1,D2,D3
      REAL(8) :: DEG,DEL,RAD,BET,BETA,XI,YJ,ZK,AIJ,AYZ,AXZ,AXY,VV
      REAL(8),ALLOCATABLE :: DX(:),DY(:),DZ(:),DXTEMP(:),DYTEMP(:),DZTEMP(:)
      DIMENSION D1(2),D2(2),D3(2)
      DIMENSION NW(3)
      SAVE NW,NX,NY,NZ
      DATA NW/'NX','NY','NZ'/
      DATA NX,NY,NZ/0,0,0/
      DATA CBLANK/'          '/
      DATA CFORMAT/'(I?)'/      
      DATA CFELEM1/'(3A1,I2,10X,A5,2E10.4,10X,3E10.4 )'/
      DATA CFELEM2/'(4A,A5,2E10.4,10X,2E10.4 ,E10.4 )'/
!
      DOM='    1'
      WRITE(CFORMAT(3:3),'(I1)',IOSTAT=IOS) IELL
!
      IF(MOP2(22).GT.0) THEN
        CFELEM1(29:33)='20.14'
        CFELEM2(21:25)='20.14'
        CFELEM2(28:32)='20.14'
        IF(MOP2(27).GT.0) CFELEM2(2:2)='2'
      ENDIF
      
      REWIND 4
      READ(NUINPUT,101) DEG
!
  103 READ(NUINPUT,100) NXYZ,NO,DEL
  100 FORMAT(A2,3X,I5,E10.4)
      DO 102 N=1,3
      IF(NXYZ.EQ.NW(N)) GOTO(111,112,113),N
  102 CONTINUE
      GOTO 120
!
  111 NX1=NX+1
      NX=NX+NO
      IF(.NOT.ALLOCATED(DX)) THEN
        ALLOCATE(DX(NX))
      ELSE
        ALLOCATE(DXTEMP(NX1-1))
        DXTEMP=DX       
        DEALLOCATE(DX)
        ALLOCATE(DX(NX))
        DX(1:SIZE(DXTEMP))=DXTEMP
        DEALLOCATE(DXTEMP)
      ENDIF
!     IF(DEL.NE.0.D0) GOTO 104
      IF(abs(DEL).gt.0.0) GOTO 104
      READ(NUINPUT,101) (DX(I),I=NX1,NX)
  101 FORMAT(8E10.4)
      GOTO 103
  104 DO 105 I=NX1,NX
  105 DX(I)=DEL
      GOTO 103
!
  112 NY1=NY+1
      NY=NY+NO
      IF(.NOT.ALLOCATED(DY)) THEN
        ALLOCATE(DY(NY))
      ELSE
        ALLOCATE(DYTEMP(NY1-1))
        DYTEMP=DY       
        DEALLOCATE(DY)
        ALLOCATE(DY(NY))
        DY(1:SIZE(DYTEMP))=DYTEMP
        DEALLOCATE(DYTEMP)
      ENDIF
!     IF(DEL.NE.0.D0) GOTO 106
      IF(abs(DEL).gt.0.0) GOTO 106
      READ(NUINPUT,101) (DY(I),I=NY1,NY)
      GOTO 103
  106 DO 107 I=NY1,NY
  107 DY(I)=DEL
      GOTO 103
!
  113 NZ1=NZ+1
      NZ=NZ+NO
      IF(.NOT.ALLOCATED(DZ)) THEN
        ALLOCATE(DZ(NZ))
      ELSE
        ALLOCATE(DZTEMP(NZ1-1))
        DZTEMP=DZ       
        DEALLOCATE(DZ)
        ALLOCATE(DZ(NZ))
        DZ(1:SIZE(DZTEMP))=DZTEMP
        DEALLOCATE(DZTEMP)
      ENDIF

!     IF(DEL.NE.0.D0) GOTO 108
      IF(abs(DEL).gt.0.0) GOTO 108
      READ(NUINPUT,101) (DZ(I),I=NZ1,NZ)
      GOTO 103
  108 DO 109 I=NZ1,NZ
  109 DZ(I)=DEL
      GOTO 103
!
  120 CONTINUE
      RAD=PI*DEG/1.8D2
      BET=SIN(RAD)
      BETA=COS(RAD)
!
!-----MAKE PRINTOUT OF MESH SPECIFICATIONS.
      NN=NX*NY*NZ
      MESHDIM=0
      IF (NX.GT.1) THEN
         MESHDIM=MESHDIM+1
!         IJKXYZ(MESHDIM)=1
      ENDIF
      IF (NY.GT.1) THEN
         MESHDIM=MESHDIM+1
!         IJKXYZ(MESHDIM)=2
      ENDIF
      IF (NZ.GT.1) THEN
         MESHDIM=MESHDIM+1
!         IJKXYZ(MESHDIM)=3
      ENDIF
      write(36,40) NX,NY,NZ,NN,DEG
   40 FORMAT(/' GENERATE CARTESIAN MESH WITH NX*NY*NZ = ', &
      I4,' *',I4,' *',I4,'  =  ',I5,'  GRID BLOCKS'/1X,131('*')// &
      '   THE X-AXIS IS HORIZONTAL; THE Y-AXIS IS ROTATED BY ',E10.4, &
      ' DEGREES AGAINST THE HORIZONTAL'// &
      '   THE GRID INCREMENTS ARE AS FOLLOWS')
      write(36, 41)NX,(NA(I),I=1,9),(DX(I),I=1,NX)
   41 FORMAT(//'   DELTA-X  (',I4,' INCREMENTS)'//13X,9(A1,10X),'10'// &
      (10X,10(1X,E10.4)))
      write(36,42) NY,(NA(I),I=1,9),(DY(I),I=1,NY)
   42 FORMAT(//'   DELTA-Y  (',I4,' INCREMENTS)'//13X,9(A1,10X),'10'// &
      (10X,10(1X,E10.4)))
      write(36,43) NZ,(NA(I),I=1,9),(DZ(I),I=1,NZ)
   43 FORMAT(//'   DELTA-Z  (',I4,' INCREMENTS)'//13X,9(A1,10X),'10'// &
      (10X,10(1X,E10.4)))
      write(36, 44)
   44 FORMAT(/' ',131('*'))
!
      IF(MOP2(2).NE.0) THEN
      IDIG=0
      JDIG=0
      KDIG=0
      IF (NX.GT.1) IDIG=IELL/MESHDIM
      IF (NY.GT.1) JDIG=IELL/MESHDIM
      IF (NZ.GT.1) KDIG=IELL/MESHDIM
      IF (NX.GE.NY.AND.NX.GE.NZ) THEN
          IDIG=IDIG+MOD(IELL,MESHDIM)
      ELSE IF (NY.GE.NZ) THEN
          JDIG=JDIG+MOD(IELL,MESHDIM)
      ELSE
          KDIG=KDIG+MOD(IELL,MESHDIM)
      ENDIF
      ENDIF
!      
      ILOOP=0
   30 CONTINUE
      ILOOP=ILOOP+1
      INDEL=0
!
         IF(ILOOP.EQ.1) WRITE(4,20) MOP2(2),NX,NY,NZ
   20    FORMAT('ELEME',I1,'  NX=',I4,'  NY=',I4,'  NZ=',I4)
         IF(ILOOP.EQ.2) WRITE(4,21)
   21    FORMAT('CONNE')
!
      XI=DX(1)/2.D0
      DO 1 I=1,NX
      I1=I+1
      IF (MOP2(27).EQ.1) THEN
         CONTINUE
      ELSE IF(MOP2(2).EQ.0) THEN
         IM=MOD(I-1,99)+1
         NOVI=(I-1)/99
         IM1=MOD(I1-1,99)+1
         NOVI1=(I1-1)/99
      ELSE
         WRITE(CFORMAT(3:3),'(I1)') IDIG
         IF (IDIG.NE.0) THEN
            WRITE(CI(:IDIG),CFORMAT) I
            WRITE(CI1(:IDIG),CFORMAT) I+1
         ENDIF
      ENDIF
      IF(I.GT.1) XI=XI+DX(I)/2.D0+DX(I-1)/2.D0
      D1(1)=DX(I)/2.D0
      IF(I.LT.NX) D1(2)=DX(I1)/2.D0
!
      YJ=DY(1)/2.D0
      DO 1 J=1,NY
      J1=J+1
!      IF (IELL.EQ.0) THEN
!         JM=MOD(J-1,62)+1
!         NOVJ=(J-1)/62
!         JM1=MOD(J1-1,62)+1
!         NOVJ1=(J1-1)/62
!      ELSEIF(IELL.EQ.5) THEN
      IF (MOP2(27).EQ.1) THEN
         CONTINUE
      ELSE IF (MOP2(2).EQ.0) THEN
         JM=MOD(J-1,61)+1
         NOVJ=(J-1)/61
         JM1=MOD(J1-1,61)+1
         JM11=MOD(j-2,61)+1
         NOVJ1=(J1-1)/61
         NOVJ11=(j-2)/61
      ELSE
         WRITE(CFORMAT(3:3),'(I1)') JDIG
         IF (JDIG.NE.0) THEN
            WRITE(CJ(:JDIG),CFORMAT) J
            WRITE(CJ1(:JDIG),CFORMAT) J+1
         ENDIF
      ENDIF
      IF(J.GT.1) YJ=YJ+DY(J)/2.D0+DY(J-1)/2.D0
      D2(1)=DY(J)/2.D0
      IF(J.LT.NY) D2(2)=DY(J1)/2.D0
!
      ZK=-DZ(1)/2.D0
      DO 1 K=1,NZ
      K1=K+1
!      IF (IELL.EQ.0) THEN
!         KM=MOD(K-1,62)+1
!         NOVK=(K-1)/62
!         KM1=MOD(K1-1,62)+1
!         NOVK1=(K1-1)/62
!!     SET CUMULATIVE "OVERFLOW" PARAMETERS, TO START AT NA(10) = 'A'
!      NOV=NOVK+((NZ-1)/62+1)*NOVJ+((NZ-1)/62+1)*((NY-1)/62+1)*NOVI+10
!      NOV1=NOVK+((NZ-1)/62+1)*NOVJ+((NZ-1)/62+1)*((NY-1)/62+1)*NOVI1+10
!      NOV2=NOVK+((NZ-1)/62+1)*NOVJ1+((NZ-1)/62+1)*((NY-1)/62+1)*NOVI+10
!      NOV3=NOVK1+((NZ-1)/62+1)*NOVJ+((NZ-1)/62+1)*((NY-1)/62+1)*NOVI+10
!      ELSEIF(IELL.EQ.5) THEN
      IF (MOP2(27).EQ.1) THEN
         CONTINUE
      ELSE IF (MOP2(2).EQ.0) THEN
         KM=MOD(K-1,61)+1
         NOVK=(K-1)/61
         KM1=MOD(K1-1,61)+1
         NOVK1=(K1-1)/61
!     SET CUMULATIVE "OVERFLOW" PARAMETERS, TO START AT NA(10) = 'A'
      NOV=NOVK+((NZ-1)/61+1)*NOVJ+((NZ-1)/61+1)*((NY-1)/61+1)*NOVI+10
      NOV1=NOVK+((NZ-1)/61+1)*NOVJ+((NZ-1)/61+1)*((NY-1)/61+1)*NOVI1+10
      NOV2=NOVK+((NZ-1)/61+1)*NOVJ1+((NZ-1)/61+1)*((NY-1)/61+1)*NOVI+10
      NOV3=NOVK1+((NZ-1)/61+1)*NOVJ+((NZ-1)/61+1)*((NY-1)/61+1)*NOVI+10
      NOV4=NOVJ1+((NY-1)/61+1)*NOVI1+10
      NOV5=NOVJ11+((NY-1)/61+1)*NOVI1+10
      ELSE
         WRITE(CFORMAT(3:3),'(I1)') KDIG
         IF (KDIG.NE.0) THEN
            WRITE(CK(:KDIG),CFORMAT) K
            WRITE(CK1(:KDIG),CFORMAT) K+1
         ENDIF
      ENDIF
      IF(K.GT.1) ZK=ZK-DZ(K)/2.D0-DZ(K-1)/2.D0
      D3(1)=DZ(K)/2.D0
      IF(K.LT.NZ) D3(2)=DZ(K1)/2.D0
!
      IF(ILOOP.EQ.2) GOTO 31
      VV=DX(I)*DY(J)*DZ(K)
!
      AIJ=0.D0
      IF(K.EQ.1 ) AIJ=AIJ+DX(I)*DY(J)
      IF(K.EQ.NZ) AIJ=AIJ+DX(I)*DY(J)
      IF (MOP2(27).EQ.1) THEN
         INDEL=INDEL+1
         WRITE(CI(1:IELL),CFORMAT) INDEL
         WRITE(4,CFELEM2) CI(:IELL),CBLANK(:15-IELL),DOM,VV,AIJ,    &
                      XI,YJ,ZK
  912 FORMAT(2A,A5,2E10.4,10X,3G10.4)     
      ELSE IF (MOP2(2).EQ.0) THEN
         WRITE(4,CFELEM1)  NA(NOV),NA(KM),NA(JM),IM,DOM,VV,AIJ, &
                           XI,YJ,ZK
      ELSE
         WRITE(4,CFELEM2) CI(:IDIG),CJ(:JDIG),CK(:KDIG),        &
                       CBLANK(:15-IELL),DOM,VV,AIJ,XI,YJ,ZK
      ENDIF
      GOTO 1
   31 CONTINUE
!
      AYZ=DY(J)*DZ(K)
      AXZ=DX(I)*DZ(K)
      AXY=DX(I)*DY(J)
      IF (MOP2(27).EQ. 1) THEN
         WRITE(CI (1:IELL),CFORMAT) (I-1)*NY*NZ+(J-1)*NZ+K
         WRITE(CI1(1:IELL),CFORMAT) (I  )*NY*NZ+(J-1)*NZ+K
         WRITE(CJ1(1:IELL),CFORMAT) (I-1)*NY*NZ+(J  )*NZ+K
         WRITE(CK1(1:IELL),CFORMAT) (I-1)*NY*NZ+(J-1)*NZ+K+1
         IF(I.LT.NX) WRITE(4,913) CI(:IELL),CI1(:IELL),             &
                     CBLANK(:18-2*IELL),NA(1),D1(1),D1(2),AYZ
         IF(J.LT.NY) WRITE(4,913) CI(:IELL),CJ1(:IELL),             &
                     CBLANK(:18-2*IELL),NA(2),D2(1),D2(2),AXZ,BET
         IF(K.LT.NZ) WRITE(4,913) CI(:IELL),CK1(:IELL),             &
                     CBLANK(:18-2*IELL),NA(3),D3(1),D3(2),AXY,BETA
  913 FORMAT(2A,11X,A,A1,4E10.4)
      ELSE IF (MOP2(2).EQ.0) THEN
         IF(I.LT.NX) WRITE(4,1111) &
         NA(NOV),NA(KM),NA(JM),IM,NA(NOV1),NA(KM),NA(JM),IM1,NA(1), &
         D1(1),D1(2),AYZ
         IF(J.LT.NY) WRITE(4,1111) &
         NA(NOV),NA(KM),NA(JM),IM,NA(NOV2),NA(KM),NA(JM1),IM,NA(2), &
         D2(1),D2(2),AXZ,BET
         IF(K.LT.NZ) WRITE(4,1111) &
         NA(NOV),NA(KM),NA(JM),IM,NA(NOV3),NA(KM1),NA(JM),IM,NA(3), &
         D3(1),D3(2),AXY,BETA
 1111 FORMAT(2(3A1,I2),19X,A1,4E10.4)
      ELSE
         IF(I.LT.NX) WRITE(4,911) CI(:IDIG),CJ(:JDIG),CK(:KDIG), &
            CI1(:IDIG),CJ(:JDIG),CK(:KDIG), &
            CBLANK(:18-2*IELL),NA(1),D1(1),D1(2),AYZ
         IF(J.LT.NY) WRITE(4,911) CI(:IDIG),CJ(:JDIG),CK(:KDIG), &
            CI(:IDIG),CJ1(:JDIG),CK(:KDIG), &
            CBLANK(:18-2*IELL),NA(2),D2(1),D2(2),AXZ,BET
         IF(K.LT.NZ) WRITE(4,911) CI(:IDIG),CJ(:JDIG),CK(:KDIG), &
            CI(:IDIG),CJ(:JDIG),CK1(:KDIG), &
            CBLANK(:18-2*IELL),NA(3),D3(1),D3(2),AXY,BETA
  911 FORMAT(6A,11X,A,A1,4E10.4)
      ENDIF
!
    1 CONTINUE
      WRITE(4,22)
   22 FORMAT('     ')
      IF(ILOOP.LT.2) GOTO 30
!
      REWIND 4
      READ(4,50) DOM
   50 FORMAT(A5)
      L=0
      DO 51 I=1,NX
      DO 51 J=1,NY
      DO 51 K=1,NZ
      L=L+1
      READ(4,50) E(L)
   51 CONTINUE
      CALL PCAR(NX,NY,NZ)
!
      RETURN
      END

!************************************************************************
      SUBROUTINE PCAR(NX,NY,NZ)
!************************************************************************      
!
!-----MAKE STRUCTURED PRINTOUT OF X-Y-Z GRID.
!
      USE MADIM
      USE MESH_ARRAY
      USE SVZ
      IMPLICIT NONE
      
      INTEGER :: NLOC
      INTEGER :: NX,NY,NZ,I,J,K,IC,N1,N2,N3
      CHARACTER*1 L1,L2,L3
      CHARACTER*2 M1,M2,M3
      DIMENSION L1(6),L2(6),L3(6),M1(6),M2(6),M3(6)
      SAVE L1,L2,L3,M1,M2,M3
      DATA L1/'I','I','J','J','K','K'/
      DATA L2/'J','K','I','K','I','J'/
      DATA L3/'K','J','K','I','J','I'/
      DATA M1/'NX','NX','NY','NY','NZ','NZ'/
      DATA M2/'NY','NZ','NX','NZ','NX','NY'/
      DATA M3/'NZ','NY','NZ','NX','NY','NX'/
!
!     FIND LOCATION OF ELEMENT NAME WITH X-INDEX I, Y-INDEX J,
!     AND Z-INDEX K
      NLOC(I,J,K)=(I-1)*NY*NZ+(J-1)*NZ+K
!
      IF(NX.LE.NY.AND.NY.LE.NZ) IC=1
      IF((NX.LE.NZ.AND.NZ.LT.NY).OR. &
         (NX.LT.NZ.AND.NZ.LE.NY)) IC=2
      IF((NY.LE.NX.AND.NX.LT.NZ).OR. &
         (NY.LT.NX.AND.NX.LE.NZ)) IC=3
      IF((NY.LE.NZ.AND.NZ.LT.NX).OR. &
         (NY.LT.NZ.AND.NZ.LE.NX)) IC=4
      IF((NZ.LE.NX.AND.NX.LT.NY).OR. &
         (NZ.LT.NX.AND.NX.LE.NY)) IC=5
      IF((NZ.LE.NY.AND.NY.LT.NX).OR. &
         (NZ.LT.NY.AND.NY.LE.NX)) IC=6
!
      write(36,31) NX,NY,NZ
   31 format(/,/, &
      ' ',131('*')/' *',20X,'CARTESIAN MESH WITH NX*NY*NZ = ', &
      I4,' *',I4,' *',I4,'  GRID BLOCKS',49X,'*'/' ',131('*'))
!
      GOTO(11,12,13,14,15,16),IC
   11 CONTINUE
      N1=NX
      N2=NY
      N3=NZ
      GOTO 17
!
   12 CONTINUE
      N1=NX
      N2=NZ
      N3=NY
      GOTO 17
!
   13 CONTINUE
      N1=NY
      N2=NX
      N3=NZ
      GOTO 17
!
   14 CONTINUE
      N1=NY
      N2=NZ
      N3=NX
      GOTO 17
!
   15 CONTINUE
      N1=NZ
      N2=NX
      N3=NY
      GOTO 17
!
   16 CONTINUE
      N1=NZ
      N2=NY
      N3=NX
      GOTO 17
!
   17 CONTINUE
      write(36,18) L1(IC),L1(IC),M1(IC),N1
   18 FORMAT(' *',120X,9X,'*'/ &
      ' *',20X,'THE MESH WILL BE PRINTED AS SLICES FOR ',A1, &
      ' = 1 TO ',A1,' = ',A2,' =',I4,49X,'*')
      write(36,19) L2(IC),L2(IC),M2(IC),N2
   19 FORMAT(' *',120X,9X,'*'/ &
      ' *',20X,'IN EACH MESH SLICE, ROWS WILL GO FROM  ',A1, &
      ' = 1 TO ',A1,' = ',A2,' =',I4,49X,'*')
      write(36,20) L3(IC),L3(IC),M3(IC),N3
   20 FORMAT(' *',120X,9X,'*'/ &
      ' *',20X,'IN EACH ROW, COLUMNS WILL GO FROM      ',A1, &
      ' = 1 TO ',A1,' = ',A2,' =',I4,49X,'*'/' *',120X,9X,'*')
!
      DO 2 I=1,N1
      write(36,3) L1(IC),I
    3 FORMAT(' ',131('*')//' SLICE WITH ',A1,' =',I4/)
      write(36,4) L3(IC),(K,K=1,N3)
    4 FORMAT('   COLUMN ',A1,' =',I9,9I10/(12X,10I10))
      WRITE(36,'(A)') ' ROWS'
      DO 5 J=1,N2
!
      GOTO(21,22,23,24,25,26),IC
   21 write(36,6) L2(IC),J,(E(NLOC(I,J,K))(:IELL),K=1,N3)
      GOTO 5
   22 write(36,6) L2(IC),J,(E(NLOC(I,K,J))(:IELL),K=1,N3)
      GOTO 5
   23 write(36,6) L2(IC),J,(E(NLOC(J,I,K))(:IELL),K=1,N3)
      GOTO 5
   24 write(36,6) L2(IC),J,(E(NLOC(K,I,J))(:IELL),K=1,N3)
      GOTO 5
   25 write(36,6) L2(IC),J,(E(NLOC(J,K,I))(:IELL),K=1,N3)
      GOTO 5
   26 write(36,6) L2(IC),J,(E(NLOC(K,J,I))(:IELL),K=1,N3)
      GOTO 5
    6 FORMAT(/,'  ',A1,' = ',I4,2X,10(1X,A)/(12X,10(1X,A)))
    5 CONTINUE
    2 CONTINUE
!
      RETURN
      END

!************************************************************************
      SUBROUTINE MINC
!************************************************************************      
!
!-----THIS IS A MODIFIED AND ENHANCED VERSION OF PROGRAM GMINC ---
!     EXECUTIVE ROUTINE FOR MAKING A "SECONDARY" FRACTURED-POROUS MEDIUM MESH
!
!-----THIS VERSION CAN ASSIGN VERTICAL CONNECTIONS BETWEEN MATRIX
!     CONTINUA (11-19-84).
!     FURTHER MODIFICATIONS (5-26-88):
!     (1) INCLUDE ELEMENT X,Y,Z-COORDINATES
!     (2) ELIMINATE "BOUNDARY ELEMENTS" FROM MINC-PROCEDURE
!
!***** GMINC WAS DEVELOPED BY KARSTEN PRUESS
!                          AT LAWRENCE BERKELEY LABORATORY. ************
!
!     THE PROGRAM GENERATES ONE-, TWO-, OR THREE-DIMENSIONAL MESHES
!     FOR FLOW SIMULATIONS IN FRACTURED POROUS MEDIA.
!
!     GMINC IMPLEMENTS THE METHOD OF
!                  MULTIPLE INTERACTING CONTINUA (MINC)
!     AS DEVELOPED BY PRUESS AND NARASIMHAN.
!
!     REFERENCES:
!
!     (1) K. PRUESS AND T.N. NARASIMHAN, A PRACTICAL METHOD FOR
!         MODELING FLUID AND HEAT FLOW IN FRACTURED POROUS MEDIA,
!         PAPER SPE-10509, PRESENTED AT THE SIXTH SPE-SYMPOSIUM ON
!         RESERVOIR SIMULATION, NEW ORLEANS, LA; (FEBRUARY 1982);
!         ALSO: SOCIETY OF PETROLEUM ENGINEERS JOURNAL, VOL. 25, NO.1,
!               PP. 14-26, 1985.
!
!     (2) K. PRUESS AND T.N. NARASIMHAN, ON FLUID RESERVES AND THE
!         PRODUCTION OF SUPERHEATED STEAM FROM FRACTURED, VAPOR-
!         DOMINATED GEOTHERMAL RESERVOIRS, J. GEOPHYS. RES. 87 (B11),
!         9329-9339, 1982.
!
!     (3) K. PRUESS AND K. KARASAKI, PROXIMITY FUNCTIONS FOR MODELING
!         FLUID AND HEAT FLOW IN RESERVOIRS WITH STOCHASTIC FRACTURE
!         DISTRIBUTIONS, PAPER PRESENTED AT EIGTH STANFORD WORKSHOP
!         ON GEOTHERMAL RESERVOIR ENGINEERING, STANFORD, CA.
!         (DECEMBER 1982).
!
!     (4) K. PRUESS, GMINC - A MESH GENERATOR FOR FLOW SIMULATIONS IN
!         FRACTURED RESERVOIRS, LAWRENCE BERKELEY LABORATORY REPORT
!         LBL-15227, 1983.
!
!***********************************************************************
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
!-----READ DATA ON MINC-PARTITIONING.
!
      CALL PART
!
!!!!! CALL TO GEOM INSERTED 11-19-84, TO OBTAIN ALL VOL(M) BEFORE
!     PROCESSING CONNECTIONS.
      CALL GEOM
!
!-----READ ELEMENT DATA FROM FILE *MESH* AND PROCESS SEQUENTIALLY.
      CALL MINCME
!
      RETURN
      END

!************************************************************************
      SUBROUTINE PART
!************************************************************************
!
!     READ SPECIFICATIONS OF MINC-PARTITIONING FROM FILE *INPUT*
!
      USE MMC
      USE MMD
      USE MINCMODULE
      USE FILEUNITS
      IMPLICIT NONE
      
      INTEGER :: I,L,M
      CHARACTER WORD*5
      DIMENSION WORD(16)
!
      READ(NUINPUT,5020) (WORD(I),I=1,16)
 5020 FORMAT(16A5)
      IF(WORD(1).EQ.'PART ') GOTO 920
!
      write(36,901) WORD(1)
  901 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL "',A5,'" IN MINC MODULE', &
      ' --- STOP EXECUTION ---')
      STOP
!
  920 CONTINUE
!
!****READ DATA FOR MULTIPLE INTERACTING CONTINUA.***********************
!
!***** *J* IS THE NUMBER OF MULTIPLE INTERACTING CONTINUA.
!*
!***** *NVOL* (.LE.J) IS THE NUMBER OF EXPLICITLY SPECIFIED VOLUME
!*     FRACTIONS.
!*
!***** *WHERE* SPECIFIES WHETHER EXPLICITLY PROVIDED VOLUME FRACTIONS
!     ARE GIVEN STARTING AT THE OUTSIDE (FRACTURE) OR INSIDE (MATRIX).
!*
!***** *PAR* IS AN ARRAY WITH PARAMETERS FOR SPECIFYING FRACTURE
!*     DISTRIBUTIONS.
!
      DO 902 LPROXI=1,6
      IF(WORD(2).EQ.FTYPE(LPROXI)) GOTO 903
  902 CONTINUE
      write(36,904) WORD(2)
  904 FORMAT(51H HAVE READ UNKNOWN PROXIMITY FUNCTION IDENTIFIER  *,A5,2 &
      2H* ---   STOP EXECUTION)
      STOP
!
  903 CONTINUE
!-----INDEX *LPROXI* LABELS THE TYPE OF PROXIMITY FUNCTION SELECTED.
!
!
      READ(NUINPUT,1) JMINCD,NVOL,WHERE,(PARMINC(I),I=1,7)
    1 FORMAT(2I3,A4,7E10.4)
      IF(.NOT.ALLOCATED(VOL)) ALLOCATE(VOL(JMINCD))
!
!-----READ A SET OF VOLUME FRACTIONS-----
      IF(WHERE.EQ.'OUT ') READ(NUINPUT,2) (VOL(M),M=1,NVOL)
      IF(WHERE.EQ.'IN  ') READ(NUINPUT,2) (VOL(JMINCD+1-M),M=1,NVOL)
    2 FORMAT(8E10.4)
!
!----- END OF MINC-DATA ------------------------------------------------
!
      IF((LPROXI.EQ.2 .OR. LPROXI.EQ.3) .AND.abs(PARMINC(2)).le.0.0) PARMINC(2)=PARMINC(1)
      IF(LPROXI.EQ.3 .AND.abs(PARMINC(3)).le.0.0) PARMINC(3)=PARMINC(2)
!
!-----IDENTIFY CHOICE MADE FOR GLOBAL MATRIX-MATRIX FLOW.
!
      DO 802 MM=1,2
      IF(WORD(4).EQ.DUAL(MM)) GOTO 803
  802 CONTINUE
      MM=0
!
  803 CONTINUE
      write(36,801) WORD(4)
  801 FORMAT(/' CHOICE OF MATRIX-MATRIX FLOW HANDLING: "',A5,'"')
      write(36, 804)
  804 FORMAT(/'                       THE OPTIONS ARE: "     "', &
      ' (DEFAULT), NO GLOBAL MATRIX-MATRIX FLOW; GLOBAL FLOW ONLY', &
      ' THROUGH FRACTURES'/ &
      40X,'"MMVER", GLOBAL MATRIX-MATRIX FLOW IN VERTICAL DIRECTION', &
      ' ONLY'/ &
      40X,'"MMALL", GLOBAL MATRIX-MATRIX FLOW IN ALL DIRECTIONS'/)
      IF(MM.EQ.2.AND.JMINCD.NE.2) write(36, 806)
  806 FORMAT( &
      ' !!!!! WARNING !!!!! THE "MMALL" OPTION SHOULD ONLY', &
      ' BE USED WITH TWO CONTINUA, WHERE IT AMOUNTS TO A DUAL-', &
      'PERMEABILITY TREATMENT')
!
! 805 CONTINUE
      RETURN
      END

!************************************************************************
      SUBROUTINE GEOM
!************************************************************************      
!
!     CALCULATE GEOMETRY PARAMETERS OF SECONDARY (FRACTURED-POROUS) MESH
!
      USE MMD
      USE MINCMODULE
      IMPLICIT NONE

!     FUNCTIONS
      REAL(8) :: PROX
!      
      INTEGER :: M,NVOL1
      REAL(8) :: DELTA,VEX,VF,TVOL,XL,XR,XMID,XMD,U,VV,W
      REAL(8), ALLOCATABLE :: X(:)
      SAVE DELTA
      DATA DELTA/1.D-8/
!
      IF (.NOT.ALLOCATED(X)) ALLOCATE(X(JMINCD))
      IF (.NOT.ALLOCATED(AMINCD)) ALLOCATE(AMINCD(JMINCD))
      IF (.NOT.ALLOCATED(DMINCD)) ALLOCATE(DMINCD(JMINCD))
      X=0.0D0
      IF(NVOL.GE.JMINCD) GOTO 3
!
!-----COME HERE TO ASSIGN EQUAL VOLUMINA TO SUBDIVISIONS WHICH HAVE NOT
!     BEEN EXPLICITLY SPECIFIED-----
!
      VEX=0.D0
      DO 4 M=1,NVOL
      IF(WHERE.EQ.'OUT ') VEX=VEX+VOL(M)
    4 IF(WHERE.EQ.'IN  ') VEX=VEX+VOL(JMINCD+1-M)
!     VEX IS THE TOTAL EXPLICITLY ASSIGNED VOLUME FRACTION.
!
      IF(VEX.GT.1.D0) GOTO 10
      IF(abs(VEX-1.).le.0.0) GOTO 3
!
      VF=(1.D0-VEX)/FLOAT(JMINCD-NVOL)
!-----VF IS THE VOLUME FRACTION FOR PARTITIONS WHICH ARE NOT
!     EXPLICITLY ASSIGNED.
      NVOL1=NVOL+1
      DO 5 M=NVOL1,JMINCD
      IF(WHERE.EQ.'OUT ') VOL(M)=VF
    5 IF(WHERE.EQ.'IN  ') VOL(JMINCD+1-M)=VF
      GOTO 3
!
   10 CONTINUE
!-----COME HERE IF EXPLICITLY ASSIGNED VOLUMINA EXCEED 100%-----
      write(36,11) VEX
   11 FORMAT(42H PROGRAM STOPS BECAUSE TOTAL VOLUME VEX = ,E12.5,40H > 1 &
      00%  ---  NEED TO CORRECT INPUT DATA)
      STOP
!
    3 CONTINUE
!
!-----NOW FIND DISTANCES FROM FRACTURES WHICH CORRESPOND TO
!     DESIRED VOLUME FRACTIONS.
!     INDEXING STARTS AT THE OUTSIDE; I.E. *1* IS THE OUTERMOST
!     VOLUME ELEMENT, AND *J* IS THE INNERMOST ONE.
!
!-----INITIALIZE TOTAL VOLUME FRACTION.
      TVOL=0.D0
!
!-----FIRST INTERFACE WILL BE AT FRACTURE FACE.
      X(1)=0.D0
      DMINCD(1)=0.D0
      AMINCD(1)=(1.D0-VOL(1))*PROX(1.D-10)/1.D-10
!
!-----INITIALIZE SEARCH INTERVAL.
!
      XL=0.D0
      XR=VOL(2)/AMINCD(1)
!
      DO 30 M=2,JMINCD
!
!-----COMPUTE TOTAL FRACTION OF MATRIX VOLUME.
      TVOL=TVOL+VOL(M)/(1.D0-VOL(1))
      IF(M.EQ.JMINCD) TVOL=1.D0-1.D-9
!
      CALL INVER(TVOL,XMID,XL,XR)
!
      X(M)=XMID
!
      XMD=XMID*DELTA
      AMINCD(M)=(1.D0-VOL(1))*(PROX(XMID+XMD)-PROX(XMID-XMD))/(2.D0*XMD)
!
      DMINCD(M)=(X(M)-X(M-1))/2.D0
!
!-----PUT LEFT END OF NEXT ITERATION INTERVAL AT PRESENT X.
      XL=XMID
!
   30 CONTINUE
!
!
!-----COME HERE TO COMPUTE A QUASI-STEADY VALUE FOR INNERMOST
!     NODAL DISTANCE.
!
      GOTO (41,42,43,44,45,46,47,48,49,50),LPROXI
!
   41 CONTINUE
!----- ONE-D CASE.
      DMINCD(JMINCD)=(PARMINC(1)-2.D0*X(JMINCD-1))/6.D0
      GOTO 40
!
   42 CONTINUE
!----- TWO-D CASE.
      U=PARMINC(1)-2.D0*X(JMINCD-1)
      VV=PARMINC(2)-2.D0*X(JMINCD-1)
      DMINCD(JMINCD)=U*VV/(4.D0*(U+VV))
      GOTO 40
!
   43 CONTINUE
!----- THRED CASE.
      U=PARMINC(1)-2.D0*X(JMINCD-1)
      VV=PARMINC(2)-2.D0*X(JMINCD-1)
      W=PARMINC(3)-2.D0*X(JMINCD-1)
      DMINCD(JMINCD)=3.D0*U*VV*W/(10.D0*(U*VV+VV*W+U*W))
      GOTO 40
!
   44 CONTINUE
   45 CONTINUE
   46 CONTINUE
   47 CONTINUE
   48 CONTINUE
   49 CONTINUE
   50 CONTINUE
      DMINCD(JMINCD)=(X(JMINCD)-X(JMINCD-1))/5.D0
!
   40 CONTINUE
!
!-----PRINT OUT GEOMETRY DATA.
      write(36, 27)
   27 FORMAT(//' ==================== GEOMETRY DATA, NORMALIZED TO', &
      ' A DOMAIN OF UNIT VOLUME ========================='//)
!
      write(36, 23)
   23 FORMAT(100H  CONTINUUM     IDENTIFIER       VOLUME      NODAL DIST&
      ANCE      INTERFACE AREA   INTERFACE DISTANCE)
      write(36, 24)
   24 FORMAT(84X,14HFROM FRACTURES/)
      write(36,25) VOL(1),DMINCD(1)
   25 FORMAT(26H  1-FRACTURES      * *    ,2(4X,E12.5))
      write(36,26) AMINCD(1),X(1)
   26 FORMAT(66X,E12.5,7X,E12.5)
!
      DO 100 M=2,JMINCD
      write(36,101) M,NA(M),VOL(M),DMINCD(M)
  101 FORMAT(1H ,I2,1H-,6HMATRIX,9X,1H*,A1,1H*,8X,E12.5,4X,E12.5)
      IF(M.NE.JMINCD) write(36,102) AMINCD(M),X(M)
  102 FORMAT(66X,E12.5,7X,E12.5)
  100 CONTINUE
!
      write(36, 103)
  103 FORMAT(/1X,99('='))
!
      RETURN
      END

!************************************************************************
      SUBROUTINE MINCME
!************************************************************************      
!
!=====THIS ROUTINE WORKS SEQUENTIALLY THROUGH THE ELEMENTS OF THE
!     PRIMARY MESH, ASSIGNING ALL SECONDARY ELEMENTS AND INTRA-BLOCK
!     CONNECTIONS.
!     PROCESS PRIMARY MESH FROM FILE *MESH*, AND WRITE
!     SECONDARY MESH ON FILE *MINC*
!
      USE MADIM
      USE NN
      USE BC
!      USE E1
      USE SVZ
      USE MMC
      USE MMD
      USE MINCMODULE
      USE OUTPUT, ONLY: NNXYZ
      IMPLICIT NONE

      INTEGER :: I,IDUM,ISIZEF,II,IOLDSIZE,IELLM1,IELLM15,IELLM18,  &
                 IOS,MA,MA1,INEWSIZE,M,MA2,NBC,NELP,NELAP,NCONP,N,  &
                 ISOT,IAC,N1,N2,IFCHAR,K,IFCHAR1,IFCHAR2
      REAL(8) :: VOLX,VV,X,Y,Z,AH,AHTX,D1,D2,AREAX,BETAX,AB,AX,PMX
      CHARACTER*1 HB,EL1,EL10,EM1
      CHARACTER*5 DENT(16),CMA
      CHARACTER*8 EL2,EM2
      CHARACTER*9 EC1,EC2,ELE
      CHARACTER   CBLANK*10
      CHARACTER   CFELIN*33,CFELFOUT*30,CFELMOUT*30,CPMX*30
      CHARACTER*9,ALLOCATABLE,DIMENSION(:)::ELEMA,CTEMP
!
      SAVE HB
      DATA HB/' '/,CBLANK/'          '/
!
      REWIND 4
      REWIND 10
      IF (NNXYZ(1).NE.0) THEN
        WRITE(10,26) MOP2(2),NNXYZ(1),NNXYZ(2),NNXYZ(3)
      ELSE
        WRITE(10,6) MOP2(2)
      ENDIF
    6 FORMAT('ELEME',I1)
   26 FORMAT('ELEME',I1,'  NX=',I4,'  NY=',I4,'  NZ=',I4)
!
      write(36, 2)
    2 FORMAT(/' READ PRIMARY MESH FROM FILE *MESH*')
!
!*****READ ELEMENT DATA FROM FILE *MESH*.*******************************
      READ(4,1) (DENT(I),I=1,16)
    1 FORMAT(16A5)
      IF(DENT(1).EQ.'ELEME') GOTO 3
      write(36,4) DENT(1)
    4 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL "',A5,'" ON FILE *MESH*', &
      ' --- STOP EXECUTION ---')
      STOP
!
    3 CONTINUE
!
      REWIND(4)
      READ(4,'(A5,I1)') DENT(1),IDUM
      IF (IDUM.EQ.0.AND.MOP2(2).EQ.0) IELL=5
      IF (IDUM.EQ.0.AND.MOP2(2).NE.0) IELL=MOP2(2)
      IF (IDUM.GE.5.AND.MOP2(2).EQ.0) THEN
        IELL=IDUM
        MOP2(2)=IELL
      ENDIF
!
      NEL=0
      NELA=0
      ISIZEF=1
      ALLOCATE(ELEMA(MNEL),STAT=iI)      
      IOLDSIZE=MNEL
!
      IELLM1=IELL-1
      IELLM15=15-IELL
      IELLM18=18-2*IELL
      CFELIN='(A1,Ax,Ayy,A5,3E10.4,3E10.4 ,A30)'
      CFELFOUT='(Ax,Ayy,I5,3E10.4,3E10.4 ,A30)'
      CFELMOUT='(Ax,Ayy,I5,3E10.4,3E10.4 )'
      WRITE(CFELIN(6:6),'(I1)')    IELLM1
      WRITE(CFELIN(9:10),'(I2)')    IELLM15
      WRITE(CFELFOUT(3:3),'(I1)')    IELL
      WRITE(CFELFOUT(6:7),'(I2)')    IELLM15
      WRITE(CFELMOUT(3:3),'(I1)')    IELL
      WRITE(CFELMOUT(6:7),'(I2)')    IELLM15
      IF(MOP2(22).GT.0) THEN
        CFELIN(24:28)='20.14'
        CFELFOUT(21:25)='20.14'
        CFELMOUT(21:25)='20.14'
      ENDIF
    9 READ(4,CFELIN) EL1,EL2(:IELLM1),CBLANK(:IELLM15),CMA,VOLX,AHTX,PMX,X,Y,Z,CPMX      
!    9 READ(4,10) EL1,EL2(:IELLM1),CBLANK(:IELLM15),CMA,VOLX,AHTX,X,Y,Z
!   10 FORMAT(A1,2A,A5,E10.4,E10.4,10X,3E10.4)
      IF(EL1.EQ.' '.AND.EL2(:IELLM1).EQ.CBLANK(:IELLM1)) GOTO 40
      READ(CMA,'(I5)',IOSTAT=IOS) MA
      IF (IOS.NE.0) CALL GETNMAT(CMA,MA)
      NEL=NEL+1
      IF(VOLX.LE.0..AND.NELA.EQ.0) NELA=NEL-1
      IF(NEL.GT.NELA.AND.NELA.NE.0) GOTO 8
!-----COME HERE FOR ACTIVE ELEMENTS
!      IF (ICHAR(EL1).GE.65) EL1=CHAR(MOD(ICHAR(EL1)-64,32)+48)
!      IF (EL1.EQ.'1') EL1=' '
      ELE(:IELL)=HB//EL2(:IELLM1)
      VV=VOL(1)*VOLX
      MA1=MA+1
      AH=AHTX
      IF(MM.EQ.1) AH=VOL(1)*AHTX
      WRITE(10,CFELFOUT) ELE(:IELL),CBLANK(:IELLM15),MA1,VV,AH,PMX,X,Y,Z,CPMX
!     STORE FRACTURE ELEMENT NAME FOR LATER CROSS REFERENCING
      IF(NEL.GT.IOLDSIZE) THEN
        ISIZEF=ISIZEF+1
        INEWSIZE=IOLDSIZE*ISIZEF
        ALLOCATE(CTEMP(IOLDSIZE),STAT=iI)
        CTEMP=ELEMA
        DEALLOCATE(ELEMA,STAT=iI)
        ALLOCATE(ELEMA(INEWSIZE),STAT=iI)
        ELEMA(1:IOLDSIZE)=CTEMP(1:IOLDSIZE)
        DEALLOCATE(CTEMP,STAT=iI)
        IOLDSIZE=INEWSIZE
      ENDIF
      ELEMA(NEL)(:IELL)=ELE(:IELL)
!      
!-----FOR EACH PRIMARY ELEMENT, ASSIGN *J* SECONDARY ELEMENTS.
      DO 11 M=2,JMINCD
      ELE=NA(M)//EL2(:IELLM1)     
      VV=VOL(M)*VOLX
      MA2=MA+2
!      MA2=MA
!      MA2(1:1)='*'
      AH=0.D0
      IF(MM.EQ.1) AH=VOL(M)*AHTX
      IF(MOP2(8).GT.0) THEN
        IFCHAR=0
        DO 6666 K=1,IELL-1
         IF (ELE(K:K).NE.' ') THEN
            IFCHAR=IFCHAR+1
         ELSE IF (ELE(K:K).EQ.' '.AND.IFCHAR.NE.0) THEN
            ELE(K:K)='0'
         ENDIF   
 6666   CONTINUE
      ENDIF
      WRITE(10,CFELMOUT) ELE(:IELL),CBLANK(:IELLM15),MA2,VV,AH,1.0D0,X,Y,Z
!
   11 CONTINUE
!
      GOTO 9
    8 CONTINUE
!-----COME HERE FOR INACTIVE ELEMENTS. THEY WILL NOT BE SUBJECTED TO
!     THE MINC PROCEDURE, BUT WILL BE HANDLED AS A SINGLE CONTINUUM.
      ELE=EL1//EL2(:IELLM1)
      WRITE(10,CFELFOUT) ELE(:IELL),CBLANK(:IELLM15),MA,VOLX,AHTX,PMX,X,Y,Z,CPMX
!
      GOTO 9
!
   40 CONTINUE
!
      WRITE(10,103)
  103 FORMAT('     ')
      WRITE(10,104)
  104 FORMAT('CONNE')
!
!-----NOW REDEFINE ELEMENT COUNTERS.
      NBC=0
      IF(NELA.NE.0) NBC=NEL-NELA
      IF(NELA.EQ.0) NELA=NEL
!     PARAMETERS FOR PRIMARY MESH
      NELP=NEL
      NELAP=NELA
!     PARAMETERS FOR SECONDARY MESH
      NELA=JMINCD*NELA
      NEL=NELA+NBC
!
!*****READ CONNECTION DATA.*********************************************
!
      N=0
      NCONP=0
      READ(4,1) DENT(1)
      IF(DENT(1).EQ.'CONNE') GOTO 1200
      write(36,1201) DENT(1)
 1201 FORMAT(' HAVE READ UNKNOWN BLOCK LABEL "',A5,'" ON FILE *MESH*', &
      ' --- STOP EXECUTION ---')
      STOP
!
 1200 CONTINUE
      READ(4,20) EL1,EL2(:IELLM1),EM1,EM2(:IELLM1), &
                     CBLANK(:IELLM18),ISOT,D1,D2,AREAX,BETAX
   20 FORMAT(2(A1,A),7X,A,I5,4E10.4)
      IF(EL1.EQ.' '.AND.EL2(:IELLM1).EQ.CBLANK(:IELLM1)) GOTO 1400
      IF(EL1.EQ.'+') GOTO 1400
      IF(MOP2(2).NE.0) THEN
        IF (ICHAR(EL1).GE.65) EL1=CHAR(MOD(ICHAR(EL1)-64,32)+48)
        IF (EL1.EQ.'1') EL1=' '
        IF (ICHAR(EM1).GE.65) EM1=CHAR(MOD(ICHAR(EM1)-64,32)+48)
        IF (EM1.EQ.'1') EM1=' '
      ENDIF
      N=N+1
      NCONP=NCONP+1
!
!-----NOW DETERMINE ACTIVITY STATUS OF ELEMENTS AT THIS CONNECTION.
!     ASSIGN THE "WOULD-BE" FRACTURE ELEMENTS, AND SEE WHETHER THEY
!     APPEAR IN THE LIST OF ACTIVE ELEMENTS
      EC1=HB//EL2(:IELLM1)
      EC2=HB//EM2(:IELLM1)
!     INITIALIZE ACTIVITY INDEX
      IAC=0
      N1=NELAP+1
      N2=NELAP+1
!
      DO 15 I=1,NELAP
      IF(ELEMA(I)(:IELL).NE.EC1(:IELL)) GOTO 16
      IAC=IAC+1
      N1=I
   16 IF(ELEMA(I)(:IELL).NE.EC2(:IELL)) GOTO 17
      IAC=IAC+1
      N2=I
   17 CONTINUE
      IF(IAC.EQ.2) GOTO 18
   15 CONTINUE
!
      IF(N1.GT.NELAP) EC1(:IELL)=EL1//EL2(:IELL-1)
      IF(N2.GT.NELAP) EC2(:IELL)=EM1//EM2(:IELL-1)
   18 CONTINUE
      IF(MOP2(8).GT.0) THEN
        IFCHAR1=0
        IFCHAR2=0
        DO 6667 K=1,IELL-1
         IF (EC1(K:K).NE.' ') THEN
            IFCHAR1=IFCHAR1+1
         ELSE IF (EC1(K:K).EQ.' '.AND.IFCHAR1.NE.0) THEN
            EC1(K:K)='0'
         ENDIF
         IF (EC2(K:K).NE.' ') THEN
            IFCHAR2=IFCHAR2+1
         ELSE IF (EC2(K:K).EQ.' '.AND.IFCHAR2.NE.0) THEN
            EC2(K:K)='0'
         ENDIF
 6667   CONTINUE   
      ENDIF
!
!-----GENERATE GLOBAL FRACTURE CONNECTION DATA.
      WRITE(10,105) EC1(:IELL),EC2(:IELL),CBLANK(:IELLM18), &
                    ISOT,D1,D2,AREAX,BETAX
  105 FORMAT(2A,7X,A,I5,4E10.4)
!
!
      IF(MM.EQ.0) GOTO 1200
!-----COME HERE FOR GLOBAL MATRIX-MATRIX FLOW CONNECTIONS.
      AB=ABS(BETAX)
      IF((AB.EQ.1.D0.AND.MM.EQ.1).OR.MM.EQ.2) GOTO 33
      GOTO 1200
!
   33 CONTINUE
!
      IF(N1.GT.NELAP.AND.N2.GT.NELAP) GOTO 1200
!!!!! DO 223 - LOOP INSERTED 11-19-84 TO ASSIGN GLOBAL CONNECTIONS
!     BETWEEN MATRIX CONTINUA.
      DO 223 M=2,JMINCD
      N=N+1
      IF(MM.EQ.1) AX=AREAX*VOL(M)
      IF(MM.EQ.2) AX=AREAX
      IF(N1.LE.NELAP) EC1=NA(M)//EL2(:IELLM1)
      IF(N2.LE.NELAP) EC2=NA(M)//EM2(:IELLM1)
      IF(MOP2(8).GT.0) THEN
        IFCHAR1=0
        IFCHAR2=0
        DO 6668 K=1,IELL-1
         IF (EC1(K:K).NE.' ') THEN
            IFCHAR1=IFCHAR1+1
         ELSE IF (EC1(K:K).EQ.' '.AND.IFCHAR1.NE.0) THEN
            EC1(K:K)='0'
         ENDIF
         IF (EC2(K:K).NE.' ') THEN
            IFCHAR2=IFCHAR2+1
         ELSE IF (EC2(K:K).EQ.' '.AND.IFCHAR2.NE.0) THEN
            EC2(K:K)='0'
         ENDIF
 6668   CONTINUE   
      ENDIF

      WRITE(10,105) EC1(:IELL),EC2(:IELL),CBLANK(:IELLM18), &
                    ISOT,D1,D2,AX,BETAX
!
  223 CONTINUE
!
      GOTO 1200
!
!-----END OF GLOBAL CONNECTION DATA.-------------------------------------------
!
 1400 CONTINUE
      write(36,5) NELP,NELAP,NCONP
    5 FORMAT('            THE PRIMARY MESH HAS',I7,' ELEMENTS (',I7, &
      ' ACTIVE) AND',I7,' CONNECTIONS (INTERFACES) BETWEEN THEM')
!
!-----NOW ASSIGN INTRA-BLOCK CONNECTION DATA.
!     LOOP AGAIN OVER ELEMENTS
      REWIND 4
      READ(4,1) (DENT(I),I=1,16)
      I=0
  110 READ(4,CFELIN) EL1,EL2(:IELLM1),CBLANK(:IELLM15),MA,VOLX,AHTX  
      IF(EL1.EQ.' '.AND.EL2(:IELLM1).EQ.CBLANK(:IELLM1)) GOTO 111
      IF(MOP2(2).NE.0) THEN
      IF (ICHAR(EL1).GE.65) EL1=CHAR(MOD(ICHAR(EL1)-64,32)+48)
      IF (EL1.EQ.'1') EL1=' '      
      ENDIF
      I=I+1
!-----FOR INACTIVE ELEMENTS, DO NOT ASSIGN INTRABLOCK CONNECTIONS
      IF(I.GT.NELAP) GOTO 111
!-----COME HERE TO ASSIGN INTRABLOCK CONNECTIONS FOR ACTIVE ELEMENTS
!     ONLY.
      EL10=HB
      BETAX=0.D0
      ISOT=1      
      DO 112 M=2,JMINCD
         N=N+1
         AREAX=VOLX*AMINCD(M-1)
         EC1=EL10//EL2(:IELLM1)
         EC2=NA(M)//EL2(:IELLM1)
         IF(MOP2(8).GT.0) THEN
            IFCHAR1=0
            IFCHAR2=0
            DO 6669 K=1,IELL-1
                IF (EC1(K:K).NE.' ') THEN
                    IFCHAR1=IFCHAR1+1
                ELSE IF (EC1(K:K).EQ.' '.AND.IFCHAR1.NE.0) THEN
                    EC1(K:K)='0'
                ENDIF
                IF (EC2(K:K).NE.' ') THEN
                    IFCHAR2=IFCHAR2+1
                ELSE IF (EC2(K:K).EQ.' '.AND.IFCHAR2.NE.0) THEN
                    EC2(K:K)='0'
                ENDIF
 6669       CONTINUE            
         ENDIF
         WRITE(10,105) EC1(:IELL),EC2(:IELL),CBLANK(:IELLM18), &
                       ISOT,DMINCD(M-1),DMINCD(M),AREAX
         EL10=NA(M)         
  112 CONTINUE
      GOTO 110
!
  111 CONTINUE
      NCON=N
      write(36,113) NEL,NELA,NCON
  113 FORMAT(/' WRITE SECONDARY MESH ON FILE *MINC*'/ &
      '          THE SECONDARY MESH HAS',I7,' ELEMENTS (',I7,' ACTIVE)', &
      ' AND',I7,' CONNECTIONS (INTERFACES) BETWEEN THEM')
!
      WRITE(10,103)
      ENDFILE 10
      DEALLOCATE(ELEMA,STAT=iI)      
!
      RETURN
      END

!************************************************************************
      SUBROUTINE INVER(F,X,XL,XR)
!************************************************************************      
!
!===== THIS ROUTINE INVERTS THE PROXIMITY FUNCTION, TO GIVE A
!     DISTANCE *X* FROM FRACTURE FACES FOR A DESIRED FRACTION *F* OF
!     MATRIX VOLUME.
!
      IMPLICIT NONE
      
      REAL(8) :: PROX
      REAL(8) :: F,X,XL,XR,TOL,FR,XMID,FMID
      SAVE TOL
      DATA TOL/1.D-10/
!
!-----CHECK AND ADJUST UPPER LIMIT OF SEARCH INTERVAL.
   22 FR=PROX(XR)
      IF(FR.GT.F) GOTO 20
      XR=2.D0*XR
      GOTO 22
!
!-----PERFORM ITERATIVE BISECTING, TO OBTAIN A SEQUENCE OF NESTED
!     INTERVALS CONTAINING THE DESIRED POINT, X.
   20 XMID=(XR+XL)/2.D0
      IF(XR-XL.LE.TOL*XR) GOTO 21
      FMID=PROX(XMID)
      IF(FMID.LE.F) XL=XMID
      IF(FMID.GE.F) XR=XMID
      GOTO 20
!
   21 CONTINUE
!-----COME HERE FOR CONVERGENCE.
      X=XMID
      RETURN
      END

!************************************************************************
      REAL*8 FUNCTION PROX(X)
!************************************************************************      
!
!-----THE PROXIMITY FUNCTION PROX(X) REPRESENTS THE FRACTION OF
!     MATRIX VOLUME [VM=(1.-VOL(1))*V0 WITHIN A DOMAIN V0] WHICH
!     IS WITHIN A DISTANCE X FROM THE FRACTURES.
!
!     CALCULATE PROXIMITY FUNCTIONS FOR DIFFERENT MATRIX BLOCK SHAPES
!
      USE MINCMODULE
      IMPLICIT NONE
      
      REAL(8) :: X,A,B,C,D,U,V,W,VR,ARGT,VT,VRT2,VTT2,VTOT
      SAVE A,B,C,D
!     NOW ASSIGN DATA FOR STANFORD LARGE RESERVOIR MODEL.
      DATA A,B,C,D/.263398D0,.190754D0,.2032D0,.191262D0/
!
      GOTO(1,2,3,4,4,4,1,1,1,1),LPROXI
!
    1 CONTINUE
!----- ONE-D CASE.
      PROX=2.D0*X/PARMINC(1)
      IF(X.GE.PARMINC(1)/2.D0) PROX=1.D0
      RETURN
!
    2 CONTINUE
!----- TWO-D CASE.
!     THE MATRIX BLOCKS HAVE THICKNESS OF PARMINC(1) AND PARMINC(2),
!     RESPECTIVELY, MEASURED PERPENDICULAR TO THE FRACTURES.
!     THE PROXIMITY FUNCTION IS VALID FOR ARBITRARY ANGLE
!     BETWEEN THE FRACTURE SETS.
      PROX=2.D0*(PARMINC(1)+PARMINC(2))*X/(PARMINC(1)*PARMINC(2)) &
      -4.D0*X*X/(PARMINC(1)*PARMINC(2))
      IF(X.GE.PARMINC(1)/2.D0 .OR. X.GE.PARMINC(2)/2.D0) PROX=1.D0
      RETURN
!
    3 CONTINUE
!----- THREE DIMENSIONAL CASE.
      U=2.D0*X/PARMINC(1)
      V=2.D0*X/PARMINC(2)
      W=2.D0*X/PARMINC(3)
      PROX=U*V*W-(U*V+U*W+V*W)+U+V+W
      IF(U.GE.1.D0 .OR. V.GE.1.D0 .OR. W.GE.1.D0) PROX=1.D0
      RETURN
!
    4 CONTINUE
!***** MATRIX OF STANFORD LARGE RESERVOIR MODEL *****
!
!
!     RECTANGULAR BLOCKS IN LAYERS B1,B2,M1,M2,T1.
      VR=8.D0*X**3-(8.D0*B+4.D0*A)*X**2+(4.D0*A*B+2.D0*B**2)*X
      IF(X.GE.B/2.D0) VR=A*B*B
!
!     TRIANGULAR BLOCKS IN LAYERS B1,B2,M1,M2,T1.
      ARGT = 2.D0
      VT=(6.D0+4.D0*SQRT(ARGT))*X**3 &
      -(A*(6.D0+4.D0*SQRT(ARGT))/2.D0+2.D0*B*(2.D0+SQRT(ARGT)))*X**2 &
      +(A*B*(2.D0+SQRT(ARGT))+B*B)*X
      IF(X.GE.B/(2.D0+SQRT(ARGT))) VT=A*B*B/2.D0
!
!     RECTANGULAR BLOCKS IN LAYER T2.
      VRT2=8.D0*X**3-(8.D0*D+4.D0*C)*X**2+(4.D0*C*D+2.D0*D**2)*X
      IF(X.GE.D/2.D0) VRT2=C*D*D
!
!     TRIANGULAR BLOCKS IN LAYER T2.
      VTT2=(6.D0+4.D0*SQRT(ARGT))*X**3 &
      -(C*(6.D0+4.D0*SQRT(ARGT))/2.D0+2.D0*D*(2.D0+SQRT(ARGT)))*X**2 &
      +(C*D*(2.D0+SQRT(ARGT))+D*D)*X
      IF(X.GE.D/(2.D0+SQRT(ARGT))) VTT2=C*D*D/2.D0
!
      IF(LPROXI.EQ.4) GOTO 14
      IF(LPROXI.EQ.5) GOTO 15
      IF(LPROXI.EQ.6) GOTO 16
!
!***** NOW COMPUTE TOTAL MATRIX VOLUME WITHIN DISTANCE X.
   14 V=5.D0*(5.D0*VR+4.D0*VT)+5.D0*VRT2+4.D0*VTT2
!-----AVERAGE PROXIMITY FUNCTION FOR ENTIRE ROCK LOADING.
!
      VTOT=35.D0*A*B**2+7.D0*C*D**2
!     VOLUME FRACTION.
      PROX=V/VTOT
      RETURN
!
   15 PROX=(5.D0*VR+4.D0*VT)/(7.D0*A*B*B)
!-----PROXIMITY FUNCTION FOR FIVE BOTTOM LAYERS.
      RETURN
!
   16 PROX=(5.D0*VRT2+4.D0*VTT2)/(7.D0*C*D*D)
!-----PROXIMITY FUNCTION FOR TOP LAYER.
      RETURN
!
!
      END

