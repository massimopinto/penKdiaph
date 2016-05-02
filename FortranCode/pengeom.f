CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C        PPPPP   EEEEEE  N    N   GGGG   EEEEEE   OOOO   M    M        C
C        P    P  E       NN   N  G    G  E       O    O  MM  MM        C
C        P    P  E       N N  N  G       E       O    O  M MM M        C
C        PPPPP   EEEE    N  N N  G  GGG  EEEE    O    O  M    M        C
C        P       E       N   NN  G    G  E       O    O  M    M        C
C        P       EEEEEE  N    N   GGGG   EEEEEE   OOOO   M    M        C
C                                                                      C
C                                                   (version 2014).    C
C                                                                      C
C  Constructive quadric geometry subroutines for Monte Carlo simula-   C
C  tion of radiation transport in material systems, using PENELOPE or  C
C  detailed simulation procedures.                                     C
C                                                                      C
C  ****  All lengths are assumed to be given in the same _arbitrary_   C
C        unit, which in PENELOPE is the cm.                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENELOPE/PENGEOM (version 2014)                                     C
C  Copyright (c) 2001-2014                                             C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Universitat de   C
C  Barcelona makes no representations about the suitability of this    C
C  software for any purpose. It is provided "as is" without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C  The module TRACK_mod contains particle track variables (position,
C  direction, current body and material), which are shared with the
C  steering main program. When PENGEOM is used with PENELOPE, this
C  module is commented off, because it is defined in the source listing
C  'penelope.f'. To use PENGEOM with other Monte Carlo codes, uncomment
C  the four lines starting with 'c>' below, include the statement 'USE
C  TRACK_mod' in all subprograms that define or modify the position or
C  direction of the particle, and compile the file 'pengeom.f' before
C  any other source files that use the present subroutines.
C
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c>    MODULE TRACK_mod   ! Particle TRACK variables.
c>    SAVE   ! Saves all items in the module.
c>    DOUBLE PRECISION :: X,Y,Z,U,V,W   ! Position and direction.
c>    INTEGER*4 :: IBODY,MAT   ! Current body and material.
c>    END MODULE TRACK_mod
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE PENGEOM_mod
C
C  ****  Geometry definition parameters and I/O quantities.
C
      SAVE  ! Saves all items in the module.
C  ----  Geometry array sizes.
C        Maximum numbers of surfaces, bodies, and limiting elements.
      INTEGER*4, PARAMETER :: NS=10000, NB=5000, NXG=250
C        Number of bodies in the material system (given by PENGEOM).
      INTEGER*4 :: NBODY
C
C  ----  Body aliases (user labels).
      CHARACTER*4 :: BALIAS(NB)='    '
C
C  ----  Body materials. MATER(KB) is the material in body KB.
      INTEGER*4 :: MATER(NB)=0
C
C  ----  Detector definition.
C        KDET(KB)=ID if body KB is part of detector ID.
      INTEGER*4 :: KDET(NB)=0
C
C  ****  Warning messages for accidental undershots or round-off errors
C        are issued when LVERB=.TRUE.
C
      LOGICAL :: LVERB=.FALSE.
C
C  ****  Last step features (output from subroutine STEP).
C
C  ----  Travelled path length, including segments in void volumes.
      DOUBLE PRECISION :: DSTOT
C  ----  Label of the last surface crossed by the particle before
C        entering a material body (defined only when NCROSS /= 0).
      INTEGER*4 :: KSLAST
C
      END MODULE PENGEOM_mod
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


C  *********************************************************************
C                       SUBROUTINE GEOMIN
C  *********************************************************************
      SUBROUTINE GEOMIN(PARINP,NPINP,NMATG,NBOD,IRD,IWR)
C
C     Reads the geometry-definition file and sets up the arrays used
C  for tracking particles through the system.
C
C  Input arguments:
C     PARINP .... array containing optional parameters, which may
C                 replace the ones entered from the input file.
C     NPINP ..... index of the highest component of PARINP that carries
C                 an effective parameter value.
C                 When NPINP is negative or zero, the numerical values
C                 of the geometry file are kept unaltered.
C     IRD ....... input unit of the geometry definition file (opened in
C                 the main program).
C     IWR ....... output unit (opened in the main program). At output
C                 this file contains a duplicate of the definition file
C                 with the effective parameter values and with elements
C                 of the various kinds labeled in strictly increasing
C                 order. This part of the output file describes the
C                 actual geometry used in the simulation. After the END
C                 line of the geometry definition block, subroutine
C                 GEOMIN writes a detailed report with the structure of
C                 the tree of modules and with information on redundant
C                 (duplicated) surfaces.
C
C  Output arguments:
C     NMATG ..... number of different materials in full bodies (exclud-
C                 ing void regions).
C     NBOD ...... Number of defined bodies.
C
      USE PENGEOM_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), CHARACTER*8 (L),
     1    INTEGER*4 (I-K,M-N)
      LOGICAL LOPCLO
C ----  Improvements suggested by F. Tola and M. Moreau  ---------------
      CHARACTER*1 CHR(6)  ! Used to check the format of surface indices.
      CHARACTER*1 CCR,CNL  ! Used to identify 'C\r' and 'C\n' comments.
C ----------------------------------------------------------------------
      CHARACTER ALIAS*5,ALIAB*5,C5*5,C5C*5,C4*4,C1*1,CA*1
      DIMENSION CA(36)
      DATA CA/ '0','1','2','3','4','5','6','7','8','9','A','B','C','D',
     1     'E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S',
     1     'T','U','V','W','X','Y','Z'/
C
      CHARACTER*12 GFILE
      CHARACTER*72 BLINE,DEFB,DEFS,DEF
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (LNUL='00000000',LSUR='SURFACE ',LIND='INDICES=',
     1           LBOD='BODY    ',LMAT='MATERIAL',LMOD='MODULE  ',
     2          LOPEN=')       ',LEND='END     ',LONE='11111111',
     3           LXSC='X-SCALE=',LYSC='Y-SCALE=',LZSC='Z-SCALE=',
     4           LTHE='  THETA=',LPHI='    PHI=',LOME='  OMEGA=',
     5           LXSH='X-SHIFT=',LYSH='Y-SHIFT=',LZSH='Z-SHIFT=',
     6           LAXX='    AXX=',LAXY='    AXY=',LAXZ='    AXZ=',
     7           LAYY='    AYY=',LAYZ='    AYZ=',LAZZ='    AZZ=',
     8            LAX='     AX=', LAY='     AY=', LAZ='     AZ=',
     9            LA0='     A0=',LDEG=') DEG   ',LRAD=') RAD   ',
     A          LCLON='CLONE   ',LSUA='SURFACE*',LINC='INCLUDE ',
     B           LINA='INCLUDE*',LFIL='   FILE=')
      DIMENSION PARINP(NPINP),LARRAY(8),KQ(5),KM(NS)
      DIMENSION IDESC(NB),IDONE(NB),ISCL(NS),IBCL(NB),IBOR(NB)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QBODY/KBODY(NB,NXG),KBOMO(NB)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
      DIMENSION DEFB(NB),DEFS(NS),ALIAS(NS),ALIAB(NB),KSFF(NS)
C
C  ************  Initialise parameters.
C
C  ****  Character constants used to identify comment lines of the forms
C  'C\r' and 'C\n', where '\r'=CHAR(13) is the carriage return character
C  and '\n'=CHAR(10) is the new line character.
      CCR=CHAR(13)
      CNL=CHAR(10)
C
C  ****  Surface coefficients, alias and KPLANE.
C
      NSFF=0  ! Number of fixed (starred) surfaces.
      DO KS=1,NS
        ALIAS(KS)='    0' ! Surface aliases (user labels).
        KSFF(KS)=0  ! Clonable-free surfaces.
        AXX(KS)=0.0D0
        AXY(KS)=0.0D0
        AXZ(KS)=0.0D0
        AYY(KS)=0.0D0
        AYZ(KS)=0.0D0
        AZZ(KS)=0.0D0
        AX(KS)=0.0D0
        AY(KS)=0.0D0
        AZ(KS)=0.0D0
        A0(KS)=0.0D0
        KPLANE(KS)=0
C  KPLANE=1 if the surface is a plane.
      ENDDO
C  **** Alias, material, mother, daughters, surfaces and side pointers
C  of bodies. The last second component of the double arrays is the
C  number of used components. For example, KSURF(KB,NXG) is the number
C  of surfaces that limit the body or module KB. The number of values
C  stored in array KFLAG(KB,--), however, is given by KSURF(KB,NXG).
      DO KB=1,NB
C  **** The array KDET(KB) is used to label bodies that are parts of
C  impact detectors. Detectors must be defined in the main program,
C  after the call to subroutine GEOMIN.
        KDET(KB)=0  ! KDET(KB).ne.0 if body KB is part of a detector.
        ALIAB(KB)='    0'  ! Body aliases (user labels).
        MATER(KB)=0  ! Material in body KB.
        KBOMO(KB)=0  ! 0 for bodies, 1 for modules.
        KMOTH(KB)=0  ! Mother of body KB (must be unique).
        DO K2=1,NXG
          KBODY(KB,K2)=0  ! Limiting bodies of body KB (keep it small).
C  KBODY(KB,K2), is the K2-th limiting body of body KB.
          KDGHT(KB,K2)=0  ! Daughters of module KB (should be small).
C  KDGHT(KB,K2), is the K2-th daughter (body or module) of module KB.
          KSURF(KB,K2)=0  ! Limiting surfaces of body KB.
C  KSURF(KB,K2), is the K2-th limiting surface of body KB.
          KFLAG(KB,K2)=5  ! Surface side pointers of material bodies.
C  KFLAG(KB,KS)=1, if KS is a limiting surface of KB and KB is inside KS
C                  (i.e. side pointer =-1).
C  KFLAG(KB,KS)=2, if KS is a limiting surface of KB and KB is outside
C                  KS (i.e. side pointer =+1).
C  KFLAG(KB,KS)=3, if KB is a body and KS does not directly limit KB,
C                  but KS appears in the definition of a body that
C                  limits KB.
C  KFLAG(KB,KS)=4, if KB is a module and KS limits one of its daughters
C                  (bodies and submodules), but does not appear in the
C                  definition of KB.
C  KFLAG(KB,KS)=5, otherwise.
        ENDDO
      ENDDO
      NMATG=0
      NSURF=0
      NBODYS=0
      ICLONE=0
      KEEPL=0
      NWARN=0
C
C  ************  Reading the geometry input file.
C
      IR=IRD
      IW=IWR
      NINCL=0
      C1=CA(NINCL+1)
      IF(IW.EQ.IR) THEN
        WRITE(IW,'(''SUBROUTINE GEOMIN. Input arguments.'')')
        WRITE(IW,'(''IRD ='',I4,'',  IWR ='',I4)') IRD,IWR
        WRITE(IW,'(''*** The input and output units must be '',
     1     ''different.'')')
        STOP
      ENDIF
C
C  ****  Look for an I/O unit that is not opened, to be connected to
C        possible INCLUDEd files.
      IRI=MAX(IR,IW)+1
 9    CONTINUE
      INQUIRE(UNIT=IRI,OPENED=LOPCLO)
      IF(LOPCLO) THEN
        IRI=IRI+1
        GOTO 9
      ENDIF
C
C
 1    CONTINUE
      READ(IR,'(A72)') BLINE
      READ(BLINE,'(A8)') LKEYW
      IF(LKEYW.NE.LNUL) THEN
        IF(IR.EQ.IRD) THEN
          WRITE(IW,'(A72)') BLINE
        ELSE
          WRITE(IW,'(A2,A72)') 'C ',BLINE
        ENDIF
        GOTO 1
      ELSE
        IF(IR.EQ.IRD) THEN
          WRITE(IW,'(64(''0''))')
        ENDIF
      ENDIF
C
 2    CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 2
      ENDIF
      READ(BLINE,'(A8,1X,A4,7A8)',ERR=900) LKEYW,C4,(LARRAY(I),I=1,7)
      IF(LKEYW.EQ.LSUR.OR.LKEYW.EQ.LSUA) THEN
        GOTO 100
      ELSE IF(LKEYW.EQ.LBOD) THEN
        GOTO 200
      ELSE IF(LKEYW.EQ.LMOD) THEN
        GOTO 300
      ELSE IF(LKEYW.EQ.LCLON) THEN
        GOTO 400
      ELSE IF(LKEYW.EQ.LINC.OR.LKEYW.EQ.LINA) THEN
        GOTO 500
      ELSE IF(LKEYW.EQ.LEND) THEN
        GOTO 600
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
C
C  ************  Surfaces.
C
 100  CONTINUE
      IF(BLINE(9:9).NE.'('.OR.BLINE(14:14).NE.')') THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect label format.'')')
        STOP
      ENDIF
      READ(BLINE,'(9X,A4,7A8)',ERR=900) C4,(LARRAY(I),I=1,7)
      IF(IR.EQ.IRD) THEN
        C5=C4//'0'
      ELSE IF(KEEPL.EQ.1) THEN
        C5=C4//'0'
      ELSE
        C5=C4//C1
      ENDIF
      IF(NSURF.GT.0) THEN
        DO KS0=1,NSURF
          IF(C5.EQ.ALIAS(KS0)) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,'(''*** Same label for two surfaces.'')')
            STOP
          ENDIF
        ENDDO
      ENDIF
      NSURF=NSURF+1
      KS=NSURF
      WRITE(IW,'(A8,''('',I4,7A8)') LKEYW,KS,(LARRAY(I),I=1,7)
      IF(KS.GT.NS) THEN
        WRITE(IW,'(''*** The parameter NS must be increased.'')')
        STOP
      ENDIF
      WRITE(DEF,'(A8,''('',I4,7A8)') LKEYW,KS,(LARRAY(I),I=1,7)
      ALIAS(KS)=C5
      IF(LKEYW.EQ.LSUA) THEN
        KSFF(KS)=1
        NSFF=NSFF+1
      ENDIF
C  ****  Indices.
 190  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 190
      ENDIF
      READ(BLINE,'(A8,5(A1,I2),A1)',ERR=900)
     1  LKEYW,(CHR(I),KQ(I),I=1,5),CHR(6)
      WRITE(IW,'(A8,5(A1,I2),A1)') LKEYW,(CHR(I),KQ(I),I=1,5),CHR(6)
C  ****  Check parentheses and commas.
      IF(CHR(1).NE.'('.OR.CHR(2).NE.','.OR.CHR(3).NE.','.
     1   OR.CHR(4).NE.','.OR.CHR(5).NE.','.OR.CHR(6).NE.')') THEN
        WRITE(IW,'(''*** Incorrect format of surface indices.'')')
        STOP
      ENDIF
C  ****  Now check the values of the indices:
      IMODE=MAX(ABS(KQ(1)),ABS(KQ(2)),ABS(KQ(3)),ABS(KQ(4)),ABS(KQ(5)))
      IF(IMODE.GT.1.OR.LKEYW.NE.LIND) THEN
        WRITE(IW,'(''*** Incorrect surface indices.'')')
        STOP
      ENDIF
      IF(IMODE.EQ.0) GOTO 103
C
C  ****  Reduced form.
C
      XSCALE=1.0D0
      YSCALE=1.0D0
      ZSCALE=1.0D0
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
C  ****  Parameters of the quadric.
 101  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 101
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I4,A8)',ERR=900)
     1  LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 102
      IF(ICHPAR.GT.0) THEN
        IF(ICHPAR.LE.NPINP) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LANGLE
          WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
          STOP
        ENDIF
      ENDIF
C  ****  Transformation parameters.
      IF(LKEYW.EQ.LXSC) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        IF(VALUE.LT.1.0D-15) THEN
          WRITE(IW,'(''*** Scale factor less than 1.0D-15.'')')
          STOP
        ENDIF
        XSCALE=VALUE
      ELSE IF(LKEYW.EQ.LYSC) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        IF(VALUE.LT.1.0D-15) THEN
          WRITE(IW,'(''*** Scale factor less than 1.0D-15.'')')
          STOP
        ENDIF
        YSCALE=VALUE
      ELSE IF(LKEYW.EQ.LZSC) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        IF(VALUE.LT.1.0D-15) THEN
          WRITE(IW,'(''*** Scale factor less than 1.0D-15.'')')
          STOP
        ENDIF
        ZSCALE=VALUE
      ELSE IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 101
C  ****  Expanded quadric.
 102  CONTINUE
      QXX=KQ(1)/XSCALE**2
      QXY=0.0D0
      QXZ=0.0D0
      QYY=KQ(2)/YSCALE**2
      QYZ=0.0D0
      QZZ=KQ(3)/ZSCALE**2
      QX=0.0D0
      QY=0.0D0
      QZ=KQ(4)/ZSCALE
      Q0=KQ(5)
C  ****  Rotated and shifted quadric.
      CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1            QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
      IF(MAX(ABS(MIN(QXX,QXY,QXZ,QYY,QYZ,QZZ)),
     1   ABS(MAX(QXX,QXY,QXZ,QYY,QYZ,QZZ))).LT.1.0D-30) KPLANE(KS)=1
      AXX(KS)=QXX
      AXY(KS)=QXY
      AXZ(KS)=QXZ
      AYY(KS)=QYY
      AYZ(KS)=QYZ
      AZZ(KS)=QZZ
      AX(KS)=QX
      AY(KS)=QY
      AZ(KS)=QZ
      A0(KS)=Q0
      WRITE(IW,'(64(''0''))')
      DEFS(KS)=DEF
      GOTO 2
C
C  ****  Implicit form.
C
 103  CONTINUE
      QXX=0.0D0
      QXY=0.0D0
      QXZ=0.0D0
      QYY=0.0D0
      QYZ=0.0D0
      QZZ=0.0D0
      QX=0.0D0
      QY=0.0D0
      QZ=0.0D0
      Q0=0.0D0
C
 193  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 193
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I4,A8)',ERR=900)
     1  LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 107
      IF(LKEYW.EQ.LONE) GOTO 104
      IF(ICHPAR.GT.0) THEN
        IF(ICHPAR.LE.NPINP) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LANGLE
          WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
          STOP
        ENDIF
      ENDIF
      IF(LKEYW.EQ.LAXX) THEN
        QXX=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAXX,QXX,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAXY) THEN
        QXY=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAXY,QXY,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAXZ) THEN
        QXZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAXZ,QXZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAYY) THEN
        QYY=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAYY,QYY,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAYZ) THEN
        QYZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAYZ,QYZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAZZ) THEN
        QZZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAZZ,QZZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAX) THEN
        QX=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAX,QX,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAY) THEN
        QY=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAY,QY,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LAZ) THEN
        QZ=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LAZ,QZ,ICHPAR,LOPEN
      ELSE IF(LKEYW.EQ.LA0) THEN
        Q0=VALUE
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LA0,Q0,ICHPAR,LOPEN
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 193
C  ****  Transformation parameters.
 104  CONTINUE
      WRITE(IW,'(64(''1''))')
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
C
 105  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 105
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I4,A8)',ERR=900)
     1  LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 106
      IF(ICHPAR.GT.0) THEN
        IF(ICHPAR.LE.NPINP) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LANGLE
          WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
          STOP
        ENDIF
      ENDIF
C
      IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 105
C  ****  Rotation and translation of the surface.
 106  CONTINUE
      CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1            QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
C
 107  CONTINUE
      TSTL=MIN(QXX,QXY,QXZ,QYY,QYZ,QZZ)
      TSTU=MAX(QXX,QXY,QXZ,QYY,QYZ,QZZ)
      IF(MAX(ABS(TSTL),ABS(TSTU)).LT.1.0D-30) KPLANE(KS)=1
      AXX(KS)=QXX
      AXY(KS)=QXY
      AXZ(KS)=QXZ
      AYY(KS)=QYY
      AYZ(KS)=QYZ
      AZZ(KS)=QZZ
      AX(KS)=QX
      AY(KS)=QY
      AZ(KS)=QZ
      A0(KS)=Q0
      WRITE(IW,'(64(''0''))')
      DEFS(KS)=DEF
      GOTO 2
C
C  ************  Bodies.
C
 200  CONTINUE
      IF(BLINE(9:9).NE.'('.OR.BLINE(14:14).NE.')') THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect label format.'')')
        STOP
      ENDIF
      READ(BLINE,'(9X,A4,7A8)',ERR=900) C4,(LARRAY(I),I=1,7)
      IF(IR.EQ.IRD) THEN
        C5=C4//'0'
      ELSE IF(KEEPL.EQ.1) THEN
        C5=C4//'0'
      ELSE
        C5=C4//C1
      ENDIF
      IF(NBODYS.GT.0) THEN
        DO KB0=1,NBODYS
          IF(C5.EQ.ALIAB(KB0)) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,
     1        '(''*** Same label for two bodies (or modules).'')')
            STOP
          ENDIF
        ENDDO
      ENDIF
      NBODYS=NBODYS+1
      WRITE(IW,'(A8,''('',I4,7A8)') LKEYW,NBODYS,(LARRAY(I),I=1,7)
      IF(NBODYS.GT.NB) THEN
        WRITE(IW,'(''*** The parameter NB must be increased.'')')
        STOP
      ENDIF
      WRITE(DEF,'(A8,''('',I4,7A8)') LKEYW,NBODYS,(LARRAY(I),I=1,7)
      ALIAB(NBODYS)=C5
C
 295  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 295
      ENDIF
      READ(BLINE,'(A8,1X,I4)',ERR=900) LKEYW,IMAT
      IF(LKEYW.NE.LMAT) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect material definition line.'')')
        STOP
      ENDIF
      WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,IMAT
      IF(IMAT.LT.0) IMAT=0
      MATER(NBODYS)=IMAT
      NMATG=MAX(NMATG,IMAT)
C
 201  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 201
      ENDIF
      READ(BLINE,'(A8)') LKEYW
      IF(LKEYW.EQ.LNUL) GOTO 208
      READ(BLINE,'(9X,A4)',ERR=900) C4
      IF(IR.EQ.IRD) THEN
        C5=C4//'0'
      ELSE IF(KEEPL.EQ.1) THEN
        C5=C4//'0'
      ELSE
        C5=C4//C1
      ENDIF
      IF(LKEYW.EQ.LSUR.OR.LKEYW.EQ.LSUA) THEN
        READ(BLINE,'(30X,I2)',ERR=900) INDS
C  ****  Surface.
        DO KS0=1,NSURF
          IF(C5.EQ.ALIAS(KS0)) THEN
            KS=KS0
            GOTO 202
          ENDIF
        ENDDO
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Undefined surface label.'')')
        STOP
 202    CONTINUE
        WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KS,INDS
C
        KST=KSURF(NBODYS,NXG)
        IF(KST.GT.0) THEN
          KST0=KST
          DO K=1,KST0
            IF(KS.EQ.KSURF(NBODYS,K)) THEN
              IF(KFLAG(NBODYS,K).LT.3) THEN
                WRITE(IW,'(''*** The last limiting surface has '',
     1            ''been defined twice.'')')
                STOP
              ELSE
                KST=K  ! KS limits a limiting body.
                GOTO 203
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        KST=KST+1
        IF(KST.GE.NXG) THEN
          WRITE(IW,
     1      '(''*** The number of limiting surfaces is too large.'')')
          STOP
        ENDIF
        KSURF(NBODYS,NXG)=KST
        KSURF(NBODYS,KST)=KS
 203    CONTINUE
        IF(INDS.EQ.-1) THEN
          KFLAG(NBODYS,KST)=1
        ELSE IF(INDS.EQ.1) THEN
          KFLAG(NBODYS,KST)=2
        ELSE
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Check side pointer value.'')')
          STOP
        ENDIF
      ELSE IF(LKEYW.EQ.LBOD) THEN
C  ****  Body.
        IF(NBODYS.GT.1) THEN
          DO KB0=1,NBODYS-1
            IF(C5.EQ.ALIAB(KB0)) THEN
              KB=KB0
              GOTO 204
            ENDIF
          ENDDO
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 204    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.0) THEN
          WRITE(IW,'(''*** This body is a module.'')')
          STOP
        ENDIF
C
        KN1=KSURF(KB,NXG)
        DO 205 KS1=1,KN1
          KSURF1=KSURF(KB,KS1)
          KN2=KSURF(NBODYS,NXG)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODYS,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 205
          ENDDO
          KST=KN2+1
          IF(KST.GE.NXG) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,
     1        '(''*** The number of limiting surfaces is too large.'')')
            STOP
          ENDIF
          KSURF(NBODYS,NXG)=KST
          KSURF(NBODYS,KST)=KSURF1
          KFLAG(NBODYS,KST)=3
 205    CONTINUE
        KBODY(NBODYS,NXG)=KBODY(NBODYS,NXG)+1
        KBODY(NBODYS,KBODY(NBODYS,NXG))=KB
      ELSE IF(LKEYW.EQ.LMOD) THEN
C  ****  Module.
        IF(NBODYS.GT.1) THEN
          DO KB0=1,NBODYS-1
            IF(C5.EQ.ALIAB(KB0)) THEN
              KB=KB0
              GOTO 206
            ENDIF
          ENDDO
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 206    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.1) THEN
          WRITE(IW,'(''*** This module is a body.'')')
          STOP
        ENDIF
C
        KN1=KSURF(KB,NXG)
        DO 207 KS1=1,KN1
          KSURF1=KSURF(KB,KS1)
          IF(KFLAG(KB,KS1).GT.2) GOTO 207
          KN2=KSURF(NBODYS,NXG)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODYS,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 207
          ENDDO
          KST=KN2+1
          IF(KST.GE.NXG) THEN
            WRITE(IW,'(A72)') BLINE
            WRITE(IW,
     1        '(''*** The number of limiting surfaces is too large.'')')
            STOP
          ENDIF
          KSURF(NBODYS,NXG)=KST
          KSURF(NBODYS,KST)=KSURF1
          KFLAG(NBODYS,KST)=3
 207    CONTINUE
        KBODY(NBODYS,NXG)=KBODY(NBODYS,NXG)+1
        KBODY(NBODYS,KBODY(NBODYS,NXG))=KB
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 201
 208  CONTINUE
      WRITE(IW,'(64(''0''))')
C
      DEFB(NBODYS)=DEF
      GOTO 2
C
C  ************  Modules.
C
 300  CONTINUE
      IF(BLINE(9:9).NE.'('.OR.BLINE(14:14).NE.')') THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect label format.'')')
        STOP
      ENDIF
      READ(BLINE,'(9X,A4,7A8)',ERR=900) C4,(LARRAY(I),I=1,7)
      IF(IR.EQ.IRD) THEN
        C5=C4//'0'
      ELSE IF(KEEPL.EQ.1) THEN
        C5=C4//'0'
      ELSE
        C5=C4//C1
      ENDIF
      IF(NBODYS.GT.0) THEN
        DO KB0=1,NBODYS
        IF(C5.EQ.ALIAB(KB0)) THEN
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,
     1      '(''*** Same label for two bodies (or modules).'')')
          STOP
        ENDIF
        ENDDO
      ENDIF
      NBODYS=NBODYS+1
      WRITE(IW,'(A8,''('',I4,7A8)') LKEYW,NBODYS,(LARRAY(I),I=1,7)
      IF(NBODYS.GT.NB) THEN
        WRITE(IW,'(''*** The parameter NB must be increased.'')')
        STOP
      ENDIF
      WRITE(DEF,'(A8,''('',I4,7A8)') LKEYW,NBODYS,(LARRAY(I),I=1,7)
      ALIAB(NBODYS)=C5
C
 391  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 391
      ENDIF
      READ(BLINE,'(A8,1X,I4)',ERR=900) LKEYW,IMAT
      WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,IMAT
      IF(LKEYW.NE.LMAT) THEN
        WRITE(IW,'(''*** Incorrect material definition line.'')')
        STOP
      ENDIF
      IF(IMAT.LT.0) IMAT=0
      MATER(NBODYS)=IMAT
      NMATG=MAX(NMATG,IMAT)
C
      KDGHT(NBODYS,NXG)=1
      KDGHT(NBODYS,1)=NBODYS
C
C  ****  Limiting surfaces and components.
C
 301  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 301
      ENDIF
      READ(BLINE,'(A8)') LKEYW
      IF(LKEYW.EQ.LNUL.OR.LKEYW.EQ.LONE) THEN
        KDT=KDGHT(NBODYS,NXG)
C  ****  Sort daughters in increasing order.
        IF(KDT.GT.1) THEN
          DO KI=1,KDT-1
            KBMIN=KDGHT(NBODYS,KI)
            KMIN=KI
            DO KJ=KI+1,KDT
              IF(KDGHT(NBODYS,KJ).LT.KBMIN) THEN
                KBMIN=KDGHT(NBODYS,KJ)
                KMIN=KJ
              ENDIF
            ENDDO
            IF(KMIN.NE.KI) THEN
              KSAVE=KDGHT(NBODYS,KI)
              KDGHT(NBODYS,KI)=KDGHT(NBODYS,KMIN)
              KDGHT(NBODYS,KMIN)=KSAVE
            ENDIF
          ENDDO
        ENDIF
        IF(LKEYW.EQ.LONE) GOTO 308
        IF(LKEYW.EQ.LNUL) GOTO 312
      ENDIF
C
C  ****  Limiting surface.
C
      IF(LKEYW.EQ.LSUR.OR.LKEYW.EQ.LSUA) THEN
        READ(BLINE,'(9X,A4,17X,I2)',ERR=900) C4,INDS
        IF(IR.EQ.IRD) THEN
          C5=C4//'0'
        ELSE IF(KEEPL.EQ.1) THEN
          C5=C4//'0'
        ELSE
          C5=C4//C1
        ENDIF
        DO KS0=1,NSURF
          IF(C5.EQ.ALIAS(KS0)) THEN
            KS=KS0
            GOTO 302
          ENDIF
        ENDDO
        WRITE(IW,'(A8,''('',A4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,C5,INDS
        WRITE(IW,'(''*** Undefined surface label.'')')
        STOP
 302    CONTINUE
        WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1    LKEYW,KS,INDS
C
        KST=KSURF(NBODYS,NXG)
        IF(KST.GT.0) THEN  ! Check whether KS is in the list.
          KST0=KST
          DO K=1,KST0
            IF(KS.EQ.KSURF(NBODYS,K)) THEN
              IF(KFLAG(NBODYS,K).LT.3) THEN
                WRITE(IW,'(''*** The last limiting surface has '',
     1            ''been defined twice.'')')
                STOP
              ELSE
                KST=K
                GOTO 303
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        KST=KST+1
        IF(KST.GE.NXG) THEN
          WRITE(IW,
     1      '(''*** The number of limiting surfaces is too large.'')')
          STOP
        ENDIF
        KSURF(NBODYS,NXG)=KST
        KSURF(NBODYS,KST)=KS
 303    CONTINUE
        IF(INDS.EQ.-1) THEN
          KFLAG(NBODYS,KST)=1
        ELSE IF(INDS.EQ.1) THEN
          KFLAG(NBODYS,KST)=2
        ELSE
          WRITE(IW,'(''*** Check side pointer value.'')')
          STOP
        ENDIF
C
C  ****  Body.
C
      ELSE IF(LKEYW.EQ.LBOD) THEN
        READ(BLINE,'(9X,A4)',ERR=900) C4
        IF(IR.EQ.IRD) THEN
          C5=C4//'0'
        ELSE IF(KEEPL.EQ.1) THEN
          C5=C4//'0'
        ELSE
          C5=C4//C1
        ENDIF
        IF(NBODYS.GT.1) THEN  ! Check whether KB is in the list.
          DO KB0=1,NBODYS-1
            IF(C5.EQ.ALIAB(KB0)) THEN
              KB=KB0
              GOTO 304
            ENDIF
          ENDDO
          WRITE(IW,'(A8,''('',A4,'')'')') LKEYW,C4
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 304    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.0) THEN
          WRITE(IW,'(''*** This body is a module.'')')
          STOP
        ENDIF
        IF(KMOTH(KB).GT.0.AND.KMOTH(KB).NE.NBODYS) THEN
          WRITE(IW,'(''*** You are trying to assign two mothers to '',
     1      ''the last body.'')')
          STOP
        ENDIF
        KMOTH(KB)=NBODYS
        KDT=KDGHT(NBODYS,NXG)+1
        KDGHT(NBODYS,NXG)=KDT
        KDGHT(NBODYS,KDT)=KB
        DO K2=1,KBODY(KB,NXG)
          IF(KMOTH(KBODY(KB,K2)).EQ.0) KMOTH(KBODY(KB,K2))=NBODYS
        ENDDO
C  ****  Assign genealogical attributes to the sisters of body KB.
        K2M=KBODY(KB,NXG)
        DO K2=1,K2M
          KBD=KBODY(KB,K2)
          IF(KMOTH(KBD).EQ.0) KMOTH(KBD)=NBODYS
          IDGHT=0
          DO K=1,KDT
            IF(KDGHT(NBODYS,K).EQ.KBD) IDGHT=K
          ENDDO
          IF(IDGHT.EQ.0) THEN
            KDGHT(NBODYS,NXG)=KDGHT(NBODYS,NXG)+1
            KDGHT(NBODYS,KDGHT(NBODYS,NXG))=KB
          ENDIF
        ENDDO
C  ****  Surfaces of the sister bodies.
        KN1=KSURF(KB,NXG)
        DO 317 KS1=1,KN1
          IF(KFLAG(KB,KS1).GT.3) GOTO 317
          KSURF1=KSURF(KB,KS1)
          KN2=KSURF(NBODYS,NXG)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODYS,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 317
          ENDDO
          KN2=KN2+1
          IF(KN2.GE.NXG) THEN
            WRITE(IW,
     1        '(''*** The number of limiting surfaces is too large.'')')
            STOP
          ENDIF
          KSURF(NBODYS,NXG)=KN2
          KSURF(NBODYS,KN2)=KSURF1
          KFLAG(NBODYS,KN2)=4
 317    CONTINUE
C
C  ****  Module.
C
      ELSE IF(LKEYW.EQ.LMOD) THEN  ! Check whether KB is in the list.
        READ(BLINE,'(9X,A4)',ERR=900) C4
        IF(IR.EQ.IRD) THEN
          C5=C4//'0'
        ELSE IF(KEEPL.EQ.1) THEN
          C5=C4//'0'
        ELSE
          C5=C4//C1
        ENDIF
        IF(NBODYS.GT.0) THEN
          DO KB0=1,NBODYS-1
            IF(C5.EQ.ALIAB(KB0)) THEN
              KB=KB0
              GOTO 305
            ENDIF
          ENDDO
          WRITE(IW,'(A8,''('',A4,'')'')') LKEYW,C4
          WRITE(IW,'(''*** Undefined body label.'')')
          STOP
        ENDIF
 305    CONTINUE
        WRITE(IW,'(A8,''('',I4,'')'')') LKEYW,KB
        IF(KBOMO(KB).NE.1) THEN
          WRITE(IW,'(''*** This module is a body.'')')
          STOP
        ENDIF
        IF(KMOTH(KB).GT.0.AND.KMOTH(KB).NE.NBODYS) THEN
          WRITE(IW,'(''*** You are trying to assign two mothers to '',
     1      ''the last module.'')')
          STOP
        ENDIF
        KMOTH(KB)=NBODYS
        KDT=KDGHT(NBODYS,NXG)+1
        KDGHT(NBODYS,NXG)=KDT
        KDGHT(NBODYS,KDT)=KB
C  ****  Limiting surfaces.
        KN1=KSURF(KB,NXG)
        DO 307 KS1=1,KN1
          IF(KFLAG(KB,KS1).GT.2) GOTO 307
          KSURF1=KSURF(KB,KS1)
          KN2=KSURF(NBODYS,NXG)
          DO KS2=1,KN2
            KSURF2=KSURF(NBODYS,KS2)
            IF(KSURF2.EQ.KSURF1) GOTO 307
          ENDDO
          KN2=KN2+1
          IF(KN2.GE.NXG) THEN
            WRITE(IW,
     1        '(''*** The number of limiting surfaces is too large.'')')
            STOP
          ENDIF
          KSURF(NBODYS,NXG)=KN2
          KSURF(NBODYS,KN2)=KSURF1
          KFLAG(NBODYS,KN2)=4
 307    CONTINUE
C
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 301
C
C  ****  Transformation parameters.
C
 308  CONTINUE
      WRITE(IW,'(64(''1''))')
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
C
 309  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 309
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I4,A8)',ERR=900)
     1  LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 310
      IF(ICHPAR.GT.0) THEN
        IF(ICHPAR.LE.NPINP) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LANGLE
          WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
          STOP
        ENDIF
      ENDIF
      IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 309
C
C  ****  Rotation and translation of the module contents (its surfaces).
C
 310  CONTINUE
      DO KS=1,NSURF
        KM(KS)=0
      ENDDO
      DO KB=1,NBODYS
        KBMOTH=KB  ! We need to transform all the descendants.
 311    CONTINUE
        IF(KBMOTH.EQ.NBODYS) THEN
          KNS=KSURF(KB,NXG)
          DO KSS=1,KNS
            KS=KSURF(KB,KSS)
            IF(KFLAG(KB,KSS).LT.5.AND.KM(KS).EQ.0
     1        .AND.KSFF(KS).EQ.0) THEN
              QXX=AXX(KS)
              QXY=AXY(KS)
              QXZ=AXZ(KS)
              QYY=AYY(KS)
              QYZ=AYZ(KS)
              QZZ=AZZ(KS)
              QX=AX(KS)
              QY=AY(KS)
              QZ=AZ(KS)
              Q0=A0(KS)
              CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1                    QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
              AXX(KS)=QXX
              AXY(KS)=QXY
              AXZ(KS)=QXZ
              AYY(KS)=QYY
              AYZ(KS)=QYZ
              AZZ(KS)=QZZ
              AX(KS)=QX
              AY(KS)=QY
              AZ(KS)=QZ
              A0(KS)=Q0
              KM(KS)=1
            ENDIF
          ENDDO
        ELSE
         KBMOTH=KMOTH(KBMOTH)  ! Moves a generation up (grandmother).
         IF(KBMOTH.GT.0) GOTO 311
        ENDIF
      ENDDO
C
 312  CONTINUE
      WRITE(IW,'(64(''0''))')
C
      KBOMO(NBODYS)=1
      DEFB(NBODYS)=DEF
      GOTO 2
C
C  ************  Clone a module.
C
 400  CONTINUE
      IF(BLINE(9:9).NE.'('.OR.BLINE(14:14).NE.')') THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect label format.'')')
        STOP
      ENDIF
      WRITE(IW,'(''C '')')
      WRITE(IW,'(''C '',62(''>''))')
      WRITE(IW,'(''C ************  Cloned module: '')')
      WRITE(IW,'(''C '')')
      ICLONE=ICLONE+1
      READ(BLINE,'(9X,A4,7A8)',ERR=900) C4,(LARRAY(I),I=1,7)
      IF(IR.EQ.IRD) THEN
        C5=C4//'0'
      ELSE IF(KEEPL.EQ.1) THEN
        C5=C4//'0'
      ELSE
        C5=C4//C1
      ENDIF
      C5C=C5
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
         WRITE(IW,'(A72)') BLINE
      ELSE
         WRITE(IW,'(''# '',A72)') BLINE
      ENDIF
      IF(NBODYS.GT.0) THEN
        DO KB0=1,NBODYS
        IF(C5.EQ.ALIAB(KB0)) THEN
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,
     1      '(''*** Same label for two bodies or modules.'')')
          STOP
        ENDIF
        ENDDO
      ENDIF
C
 401  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 401
      ENDIF
      READ(BLINE,'(A8,1X,A4)',ERR=900) LKEYW,C4
      IF(IR.EQ.IRD) THEN
        C5=C4//'0'
      ELSE IF(KEEPL.EQ.1) THEN
        C5=C4//'0'
      ELSE
        C5=C4//C1
      ENDIF
      WRITE(IW,'(''# '',A8,''('',A4,'')'')') LKEYW,C4
      IF(LKEYW.NE.LMOD) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** The cloned object must be a module.'')')
        STOP
      ENDIF
      IF(NBODYS.EQ.0) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** This module is not defined.'')')
        STOP
      ENDIF
      DO KB0=1,NBODYS
        IF(C5.EQ.ALIAB(KB0)) THEN
          KORIG=KB0
          IF(KBOMO(KORIG).NE.1) THEN
            WRITE(IW,'(''*** The cloned object must be a module.'')')
            WRITE(IW,'(''*** The selected object is a body.'')')
            STOP
          ENDIF
          GOTO 402
        ENDIF
      ENDDO
      WRITE(IW,'(''*** The label does not correspond to a module.'')')
      STOP
C
 402  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 402
      ENDIF
      READ(BLINE,'(A8)') LKEYW
      IF(LKEYW.EQ.LONE.OR.LKEYW.EQ.LNUL) THEN
        GOTO 403
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
C
C  ****  Transformation parameters.
C
 403  CONTINUE
      OMEGA=0.0D0
      THETA=0.0D0
      PHI=0.0D0
      XSHIFT=0.0D0
      YSHIFT=0.0D0
      ZSHIFT=0.0D0
      IF(LKEYW.EQ.LNUL) GOTO 405
C
 404  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 404
      ENDIF
      READ(BLINE,'(A8,1X,E22.15,1X,I4,A8)',ERR=900)
     1  LKEYW,VALUE,ICHPAR,LANGLE
      IF(LKEYW.EQ.LNUL) GOTO 405
      IF(ICHPAR.GT.0) THEN
        IF(ICHPAR.LE.NPINP) THEN
          VALUE=PARINP(ICHPAR)
          ICHPAR=-ICHPAR  ! Switches off the parameter changing option.
        ELSE
          WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LANGLE
          WRITE(IW,'(''*** NPINP is too small (check PARINP).'')')
          STOP
        ENDIF
      ENDIF
      IF(LKEYW.EQ.LOME) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          OMEGA=VALUE
        ELSE
          WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          OMEGA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LTHE) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          THETA=VALUE
        ELSE
          WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          THETA=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LPHI) THEN
        IF(LANGLE.EQ.LRAD) THEN
          WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LRAD
          PHI=VALUE
        ELSE
          WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1      LKEYW,VALUE,ICHPAR,LDEG
          PHI=VALUE*PI/180.0D0
        ENDIF
      ELSE IF(LKEYW.EQ.LXSH) THEN
        WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        XSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LYSH) THEN
        WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        YSHIFT=VALUE
      ELSE IF(LKEYW.EQ.LZSH) THEN
        WRITE(IW,'(''# '',A8,''('',1P,E22.15,'','',I4,A8)')
     1    LKEYW,VALUE,ICHPAR,LOPEN
        ZSHIFT=VALUE
      ELSE
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ENDIF
      GOTO 404
 405  CONTINUE
C  ****  Determine all the descendants of module KORIG.
      ND=1
      IDESC(1)=KORIG  ! Descendants of the cloned module.
      IDONE(1)=0      ! The descendants have not yet been identified.
 406  CONTINUE
      NDC=ND
      KDG=0
      DO I=1,NDC
        IF(IDONE(I).EQ.0) THEN
          KB=IDESC(I)
          IF(KBODY(KB,NXG).GT.0) THEN
            DO J=1,KBODY(KB,NXG)
              ND=ND+1
              IDESC(ND)=KBODY(KB,J)  ! New descendant.
              IDONE(ND)=0
              KDG=KDG+1
            ENDDO
          ELSE IF(KDGHT(KB,NXG).GT.0) THEN
            DO J=1,KDGHT(KB,NXG)
              IF(KDGHT(KB,J).NE.KB) THEN
                ND=ND+1
                IDESC(ND)=KDGHT(KB,J)  ! New descendant.
                IDONE(ND)=0
                KDG=KDG+1
              ENDIF
            ENDDO
          ENDIF
          IDONE(I)=1  ! The descendants of KB=IDESC(I) have been listed.
        ENDIF
      ENDDO
      IF(KDG.GT.0) GOTO 406
C
      IN=0
      DO I=1,NB
        IBCL(I)=0  ! Label of a cloned body or module.
        IBOR(I)=0  ! Label of the original body or module.
      ENDDO
      KSD=NSURF
      DO I=1,NS
        ISCL(I)=0  ! Label of a cloned surface.
      ENDDO
      KBD=NBODYS
      DO 407 KBB=1,NBODYS
        DO ID=ND,1,-1
          KB=IDESC(ID)
          IF(KBB.EQ.KB) THEN
            KBD=KBD+1
            IBCL(KB)=KBD
            IBOR(KBD)=KB
            MATER(KBD)=MATER(KB)
C  ****  Clone the surfaces of the original module and its descendants.
            KSURF(KBD,NXG)=KSURF(KB,NXG)
            DO KSS=KSURF(KB,NXG),1,-1
              KS=KSURF(KB,KSS)
              IF(KFLAG(KB,KSS).LT.3.AND.ISCL(KS).EQ.0) THEN
                IF(KSFF(KS).EQ.0) THEN
                  KSD=KSD+1
                  IF(KSD.GT.NS) THEN
                    WRITE(IW,
     1                '(''*** The parameter NS must be increased.'')')
                    STOP
                  ENDIF
                  DEFS(KSD)=DEFS(KS)
                  ISCL(KS)=KSD
                  BLINE=DEFS(KSD)
                  QXX=AXX(KS)
                  QXY=AXY(KS)
                  QXZ=AXZ(KS)
                  QYY=AYY(KS)
                  QYZ=AYZ(KS)
                  QZZ=AZZ(KS)
                  QX=AX(KS)
                  QY=AY(KS)
                  QZ=AZ(KS)
                  Q0=A0(KS)
                  CALL ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,
     1                        QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0)
                  AXX(KSD)=QXX
                  AXY(KSD)=QXY
                  AXZ(KSD)=QXZ
                  AYY(KSD)=QYY
                  AYZ(KSD)=QYZ
                  AZZ(KSD)=QZZ
                  AX(KSD)=QX
                  AY(KSD)=QY
                  AZ(KSD)=QZ
                  A0(KSD)=Q0
                  KSFF(KSD)=0
                  WRITE(IW,'(A8,''('',I4,A51)') LSUR,KSD,BLINE(14:64)
                  WRITE(IW,'(A8,''('',4(I2,'',''),I2,'')'')')
     1              LIND,IN,IN,IN,IN,IN
                  IF(ABS(AXX(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAXX,AXX(KSD),IN,LOPEN
                  IF(ABS(AXY(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAXY,AXY(KSD),IN,LOPEN
                  IF(ABS(AXZ(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAXZ,AXZ(KSD),IN,LOPEN
                  IF(ABS(AYY(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAYY,AYY(KSD),IN,LOPEN
                  IF(ABS(AYZ(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAYZ,AYZ(KSD),IN,LOPEN
                  IF(ABS(AZZ(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAZZ,AZZ(KSD),IN,LOPEN
                  IF(ABS(AX(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAX,AX(KSD),IN,LOPEN
                  IF(ABS(AY(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAY,AY(KSD),IN,LOPEN
                  IF(ABS(AZ(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LAZ,AZ(KSD),IN,LOPEN
                  IF(ABS(A0(KSD)).GT.1.0D-20)
     1              WRITE(IW,'(A8,''('',1P,E22.15,'','',I4,A8)')
     2              LA0,A0(KSD),IN,LOPEN
                  WRITE(IW,'(64(''0''))')
                ELSE
                  ISCL(KS)=KS
                ENDIF
              ENDIF
            ENDDO
            GOTO 407
          ENDIF
        ENDDO
 407  CONTINUE
C  ****  Clone the original module and its descendants.
      DO KB=NBODYS+1,KBD
        IF(KB.GT.NB) THEN
          WRITE(IW,'(''*** The parameter NB must be increased.'')')
          STOP
        ENDIF
        KBO=IBOR(KB)
        IF(KMOTH(KBO).GT.0) THEN
          KMOTH(KB)=IBCL(KMOTH(KBO))
        ELSE
          KMOTH(KB)=0
        ENDIF
        KBOMO(KB)=KBOMO(KBO)
        DEFB(KB)=DEFB(KBO)
        BLINE=DEFB(KBO)
        IF(KBOMO(KB).EQ.0) THEN
          LKEYW=LBOD
          WRITE(IW,'(A8,''('',I4,A51)') LBOD,KB,BLINE(14:64)
        ELSE IF(KBOMO(KB).EQ.1) THEN
          LKEYW=LMOD
          WRITE(IW,'(A8,''('',I4,A51)') LMOD,KB,BLINE(14:64)
        ELSE
          WRITE(IW,'(''KBOMO('',I4,'') ='',I4)') KB,KBOMO(KB)
          WRITE(IW,'(''*** Something wrong...'')')
          STOP
        ENDIF
        WRITE(IW,'(A8,''('',I4,'')'')') LMAT,MATER(KB)
        DO KS=1,KSURF(KB,NXG)
          KSURF(KB,KS)=ISCL(KSURF(KBO,KS))
          KFLAG(KB,KS)=KFLAG(KBO,KS)
          IF(KFLAG(KB,KS).LT.3) THEN
            IF(KFLAG(KB,KS).EQ.1) THEN
              INDS=-1
            ELSE
              INDS=+1
            ENDIF
            WRITE(IW,'(A8,''('',I4,''), SIDE POINTER=('',I2,'')'')')
     1        LSUR,KSURF(KB,KS),INDS
          ENDIF
        ENDDO
        IF(KBOMO(KB).EQ.0) THEN
          KBODY(KB,NXG)=KBODY(KBO,NXG)
          DO I=1,KBODY(KB,NXG)
            KBODY(KB,I)=IBCL(KBODY(KBO,I))
            KBB=KBODY(KB,I)
            IF(KBOMO(KBB).EQ.0) THEN
              LKEYW=LBOD
              WRITE(IW,'(A8,''('',I4,'')'')') LBOD,KBB
            ELSE IF(KBB.NE.KB) THEN
              LKEYW=LMOD
              WRITE(IW,'(A8,''('',I4,'')'')') LMOD,KBB
            ENDIF
            IF(KBB.GE.KB) THEN
              WRITE(IW,'(''*** The limiting body or module is not '',
     1          ''yet defined'')')
              STOP
            ENDIF
          ENDDO
        ELSE
          KDGHT(KB,NXG)=KDGHT(KBO,NXG)
          DO I=1,KDGHT(KB,NXG)
            KDGHT(KB,I)=IBCL(KDGHT(KBO,I))
            KBB=KDGHT(KB,I)
            IF(KBB.NE.KB) THEN
              IF(KBOMO(KBB).EQ.0) THEN
                LKEYW=LBOD
                WRITE(IW,'(A8,''('',I4,'')'')') LBOD,KBB
              ELSE
                LKEYW=LMOD
                WRITE(IW,'(A8,''('',I4,'')'')') LMOD,KBB
              ENDIF
              IF(KBB.GE.KB) THEN
                WRITE(IW,'(''*** The limiting body or module is not '',
     1            ''yet defined'')')
                STOP
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        WRITE(IW,'(64(''0''))')
        IF(KB.EQ.KBD) THEN
          WRITE(IW,'(''C '')')
          WRITE(IW,'(''C ************  End of cloned module.'')')
          WRITE(IW,'(''C '',62(''<''))')
          WRITE(IW,'(''C '')')
        ENDIF
      ENDDO
      ALIAB(KBD)=C5C
      KMOTH(KBD)=0
      NBODYS=KBD
      NSURF=KSD
      GOTO 2
C
C  ************  Included geometry file.
C
 500  CONTINUE
      IF(LKEYW.EQ.LINA) THEN
        KEEPL=1  ! Keep the user labels of the included file elements.
      ELSE
        KEEPL=0
      ENDIF
      WRITE(IW,'(''C '')')
      WRITE(IW,'(''C '',62(''>''))')
      WRITE(IW,'(''C ************  Included file:  '')')
      IF(KEEPL.EQ.1) WRITE(IW,'(
     1''C The included elements keep their user labels'')')
      WRITE(IW,'(''C '')')
      WRITE(IW,'(''# '',A72)') BLINE
 501  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 501
      ENDIF
      IF(BLINE(9:9).NE.'('.OR.BLINE(22:22).NE.')') THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** Incorrect label format.'')')
        STOP
      ENDIF
      READ(BLINE,'(A8,1X,A12,6A8)',ERR=900)
     1  LKEYW,GFILE,(LARRAY(I),I=1,6)
      IF(LKEYW.NE.LFIL) THEN
        WRITE(IW,'(A72)') BLINE
        WRITE(IW,'(''*** What do you mean?'')')
        STOP
      ELSE
         WRITE(IW,'(''# '',A72)') BLINE
      ENDIF
 502  CONTINUE
      READ(IR,'(A72)') BLINE
      IF(((BLINE(1:1).EQ.'C'.OR.BLINE(1:1).EQ.'c').AND.
     1  (BLINE(2:2).EQ.CCR.OR.BLINE(2:2).EQ.CNL.OR.BLINE(2:2).EQ.' '))
     2  .OR.BLINE(1:1).EQ.'#') THEN
        WRITE(IW,'(A72)') BLINE
        GOTO 502
      ENDIF
      READ(BLINE,'(8A8)') LKEYW,(LARRAY(I),I=1,6)
      IF(LKEYW.EQ.LNUL) THEN
        IF(IR.EQ.IRI) THEN
          WRITE(IW,'(A72)') BLINE
          WRITE(IW,'(''*** Too many include levels.'')')
          STOP
        ELSE
          IR=IRI
C  ****  The alias of elements in the included geometry are altered by
C  appending an extra character, to prevent having duplicated labels.
          IF(KEEPL.EQ.0) THEN
            NINCL=NINCL+1
            IF(NINCL.GT.35) THEN
              WRITE(IW,'(A,I3)') 'No. of included files =',NINCL
              WRITE(IW,'(''*** Too many included files.'')')
              STOP
            ENDIF
            C1=CA(NINCL+1)
          ENDIF
C
C  ****  Remove leading blanks in filename.
 503      CONTINUE
          IF(GFILE(1:1).EQ.' ') THEN
            GFILE=GFILE(2:)
            GOTO 503
          ENDIF
C
          OPEN(IR,FILE=GFILE)
        ENDIF
        GOTO 1
      ENDIF
C
C  ************  End-line in the input file.
C
 600  CONTINUE
      IF(IR.EQ.IRI) THEN
        CLOSE(IR)
        WRITE(IW,'(''C '')')
        WRITE(IW,'(''C ************  End of included file.'')')
          WRITE(IW,'(''C '',62(''<''))')
        WRITE(IW,'(''C '')')
        IR=IRD
        C1='0'
        GOTO 2
      ENDIF
      IF(NBODYS.EQ.1) THEN
        KDGHT(1,NXG)=1
        KDGHT(1,1)=1
      ENDIF
C
      DO KB=1,NBODYS
        C5=ALIAB(KB)
        IF(C5(5:5).EQ.'0') BALIAS(KB)=C5(1:4)
      ENDDO
C
C  ************  Check for motherless bodies or modules.
C
      MLESS=0
      KBENC=0
      DO KBB=1,NBODYS
        IF(KMOTH(KBB).EQ.0) THEN
          MLESS=MLESS+1
          KBENC=KBB
        ENDIF
      ENDDO
      IF(MLESS.EQ.1) THEN
        IF(KBOMO(KBENC).EQ.1) GOTO 602  ! There is a root module.
      ENDIF
C
C  ************  Define the enclosure.
C
      IF(NBODYS.GT.NB-1) THEN
        WRITE(IW,'(''*** The parameter NB must be increased.'')')
        STOP
      ENDIF
      IF(NSURF.GT.NS-1) THEN
        WRITE(IW,'(''*** The parameter NS must be increased.'')')
        STOP
      ENDIF
C  ****  The next line serves only to avoid a warning issued by
C        certain compilers.
      NT=0
      KS=NSURF+1
C  ****  The enclosure is a sphere centred at the origin and with a
C        radius of 1.0E7 length units.
      AXX(KS)=1.0D0
      AYY(KS)=1.0D0
      AZZ(KS)=1.0D0
      A0(KS)=-1.0D14
      KB=NBODYS+1
      KBOMO(KB)=1
      KSURF(KB,NXG)=1
      KSURF(KB,1)=KS
      KFLAG(KB,1)=1
      KDGHT(KB,NXG)=1
      KDGHT(KB,1)=KB
      DO KBB=1,NBODYS
        IF(KMOTH(KBB).EQ.0) THEN
          NT=KDGHT(KB,NXG)+1
          KDGHT(KB,NXG)=NT
          KDGHT(KB,NT)=KBB
          KN1=KSURF(KBB,NXG)
          DO 601 KS1=1,KN1
            IF(KFLAG(KBB,KS1).GT.3) GOTO 601
            KSURF1=KSURF(KBB,KS1)
            KN2=KSURF(KB,NXG)
            DO KS2=1,KN2
              IF(KSURF(KB,KS2).EQ.KSURF1) GOTO 601
            ENDDO
            KN2=KN2+1
            IF(KN2.GE.NXG) THEN
              WRITE(IW,'(''*** The parameter NXG is too small.'')')
              STOP
            ENDIF
            KSURF(KB,NXG)=KN2
            KSURF(KB,KN2)=KSURF1
            KFLAG(KB,KN2)=4
 601      CONTINUE
          KMOTH(KBB)=KB
          KN3=KN1+1
          KSURF(KBB,NXG)=KN3
          KSURF(KBB,KN3)=KS
          KFLAG(KBB,KN3)=1
        ENDIF
      ENDDO
      NSURF=KS
      NBODYS=KB
C  ****  Sort daughters in increasing order.
      IF(NT.GT.1) THEN
        DO KI=1,NT-1
          KBMIN=KDGHT(NBODYS,KI)
          KMIN=KI
          DO KJ=KI+1,NT
            IF(KDGHT(NBODYS,KJ).LT.KBMIN) THEN
              KBMIN=KDGHT(NBODYS,KJ)
              KMIN=KJ
            ENDIF
          ENDDO
          IF(KMIN.NE.KI) THEN
            KSAVE=KDGHT(NBODYS,KI)
            KDGHT(NBODYS,KI)=KDGHT(NBODYS,KMIN)
            KDGHT(NBODYS,KMIN)=KSAVE
          ENDIF
        ENDDO
      ENDIF
C
 602  CONTINUE
      WRITE(IW,'(A8,1X,55(''0''))') LKEYW
      NBOD=NBODYS
      NBODY=NBODYS
C
C  ****  Duplicated surfaces (within round-off tolerance) are removed.
C
      WRITE(IW,'(///,''*****************************************'')')
      WRITE(IW,'(''****     PENGEOM (version 2014)      ****'')')
      WRITE(IW,'(''****  Constructive Quadric Geometry  ****'')')
      WRITE(IW,'(''*****************************************'',/)')
C
      IF(NSFF.GT.0) THEN
        WRITE(IW,'(/,''WARNING: The system contains fixed (starred) '',
     1    ''surfaces, which are'',/9X,''not affected by translations'',
     2    '' and rotations. Hence, any'',/9X,''translation or '',
     3    ''rotation that does not leave these'',/9X,''surfaces '',
     4    ''invariant will distort the system.'',/)')
      ENDIF
C
      IF(NSURF.LT.2) GOTO 704
      IWRITE=0
      TOL=1.0D-14
      DO KS=1,NSURF
        KM(KS)=0
      ENDDO
      DO 703 KS=1,NSURF-1
        IF(KM(KS).NE.0) GOTO 703
        F=MAX(AXX(KS),AXY(KS),AXZ(KS),AYY(KS),AYZ(KS),AZZ(KS),
     1        AX(KS),AY(KS),AZ(KS),A0(KS))
        FM=MIN(AXX(KS),AXY(KS),AXZ(KS),AYY(KS),AYZ(KS),AZZ(KS),
     1         AX(KS),AY(KS),AZ(KS),A0(KS))
        IF(ABS(FM).GT.ABS(F)) F=FM
        DO 702 KST=KS+1,NSURF
          IF(KM(KST).NE.0) GOTO 702
          FP=MAX(AXX(KST),AXY(KST),AXZ(KST),AYY(KST),AYZ(KST),AZZ(KST),
     1           AX(KST),AY(KST),AZ(KST),A0(KST))
          FM=MIN(AXX(KST),AXY(KST),AXZ(KST),AYY(KST),AYZ(KST),AZZ(KST),
     1           AX(KST),AY(KST),AZ(KST),A0(KST))
          IF(ABS(FM).GT.ABS(FP)) FP=FM
          FFP=F/FP
          RFFP=1.0D0/FFP
C
          TST=0.0D0
          IF(ABS(AXX(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AXX(KS)-AXX(KST)*FFP)/AXX(KS)))
          ELSE IF(ABS(AXX(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AXX(KS)*RFFP-AXX(KST))/AXX(KST)))
          ENDIF
          IF(ABS(AXY(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AXY(KS)-AXY(KST)*FFP)/AXY(KS)))
          ELSE IF(ABS(AXY(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AXY(KS)*RFFP-AXY(KST))/AXY(KST)))
          ENDIF
          IF(ABS(AXZ(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AXZ(KS)-AXZ(KST)*FFP)/AXZ(KS)))
          ELSE IF(ABS(AXZ(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AXZ(KS)*RFFP-AXZ(KST))/AXZ(KST)))
          ENDIF
          IF(ABS(AYY(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AYY(KS)-AYY(KST)*FFP)/AYY(KS)))
          ELSE IF(ABS(AYY(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AYY(KS)*RFFP-AYY(KST))/AYY(KST)))
          ENDIF
          IF(ABS(AYZ(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AYZ(KS)-AYZ(KST)*FFP)/AYZ(KS)))
          ELSE IF(ABS(AYZ(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AYZ(KS)*RFFP-AYZ(KST))/AYZ(KST)))
          ENDIF
          IF(ABS(AZZ(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AZZ(KS)-AZZ(KST)*FFP)/AZZ(KS)))
          ELSE IF(ABS(AZZ(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AZZ(KS)*RFFP-AZZ(KST))/AZZ(KST)))
          ENDIF
          IF(ABS(AX(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AX(KS)-AX(KST)*FFP)/AX(KS)))
          ELSE IF(ABS(AX(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AX(KS)*RFFP-AX(KST))/AX(KST)))
          ENDIF
          IF(ABS(AY(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AY(KS)-AY(KST)*FFP)/AY(KS)))
          ELSE IF(ABS(AY(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AY(KS)*RFFP-AY(KST))/AY(KST)))
          ENDIF
          IF(ABS(AZ(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AZ(KS)-AZ(KST)*FFP)/AZ(KS)))
          ELSE IF(ABS(AZ(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((AZ(KS)*RFFP-AZ(KST))/AZ(KST)))
          ENDIF
          IF(ABS(A0(KS)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((A0(KS)-A0(KST)*FFP)/A0(KS)))
          ELSE IF(ABS(A0(KST)).GT.1.0D-16) THEN
            TST=MAX(TST,ABS((A0(KS)*RFFP-A0(KST))/A0(KST)))
          ENDIF
C
          IF(TST.LT.TOL) THEN
            IF(IWRITE.EQ.0) THEN
              WRITE(IW,'(/,''************  Removal of duplicated '',
     1          ''user-defined surfaces.'',/)')
              IWRITE=1
            ENDIF
            WRITE(IW,'(''SURFACE ('',I4,'') is replaced by SURFACE ('',
     1        I4,'')'' )') KST,KS
C  ****  Check whether the two surface functions have the same global
C        sign (ISPF=0) or not (ISPF=1).
            ISPF=0
            TF=TOL*F
            IF(ABS(AXX(KS)).GT.TF.AND.AXX(KS)*AXX(KST).LT.0.0D0) ISPF=1
            IF(ABS(AXY(KS)).GT.TF.AND.AXY(KS)*AXY(KST).LT.0.0D0) ISPF=1
            IF(ABS(AXZ(KS)).GT.TF.AND.AXZ(KS)*AXZ(KST).LT.0.0D0) ISPF=1
            IF(ABS(AYY(KS)).GT.TF.AND.AYY(KS)*AYY(KST).LT.0.0D0) ISPF=1
            IF(ABS(AYZ(KS)).GT.TF.AND.AYZ(KS)*AYZ(KST).LT.0.0D0) ISPF=1
            IF(ABS(AZZ(KS)).GT.TF.AND.AZZ(KS)*AZZ(KST).LT.0.0D0) ISPF=1
            IF(ABS( AX(KS)).GT.TF.AND. AX(KS)* AX(KST).LT.0.0D0) ISPF=1
            IF(ABS( AY(KS)).GT.TF.AND. AY(KS)* AY(KST).LT.0.0D0) ISPF=1
            IF(ABS( AZ(KS)).GT.TF.AND. AZ(KS)* AZ(KST).LT.0.0D0) ISPF=1
            IF(ABS( A0(KS)).GT.TF.AND. A0(KS)* A0(KST).LT.0.0D0) ISPF=1
C
            IF(TST.GT.1.0D-15) THEN
              WRITE(IW,'(A,1P,4E15.7)') 'F,FP,F/FP,TST =',F,FP,FFP,TST
              WRITE(IW,'(A,1P,3E15.7)') 'AXX,AXXP,DIFF =',
     1          AXX(KS),AXX(KST),AXX(KS)-AXX(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AXY,AXYP,DIFF =',
     1          AXY(KS),AXY(KST),AXY(KS)-AXY(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AXZ,AXZP,DIFF =',
     1          AXZ(KS),AXZ(KST),AXZ(KS)-AXZ(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AYY,AYYP,DIFF =',
     1          AYY(KS),AYY(KST),AYY(KS)-AYY(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AYZ,AYZP,DIFF =',
     1          AYZ(KS),AYZ(KST),AYZ(KS)-AYZ(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AZZ,AZZP,DIFF =',
     1          AZZ(KS),AZZ(KST),AZZ(KS)-AZZ(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AX,AXP,DIFF   =',
     1          AX(KS),AX(KST),AX(KS)-AX(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AY,AYP,DIFF   =',
     1          AY(KS),AY(KST),AY(KS)-AY(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'AZ,AZP,DIFF   =',
     1          AZ(KS),AZ(KST),AZ(KS)-AZ(KST)*FFP
              WRITE(IW,'(A,1P,3E15.7)') 'A0,A0P,DIFF   =',
     1          A0(KS),A0(KST),A0(KS)-A0(KST)*FFP
              WRITE(IW,'(A)')  ' '
            ENDIF
C
            AXX(KST)=AXX(KS)
            AXY(KST)=AXY(KS)
            AXZ(KST)=AXZ(KS)
            AYY(KST)=AYY(KS)
            AYZ(KST)=AYZ(KS)
            AZZ(KST)=AZZ(KS)
            AX(KST)=AX(KS)
            AY(KST)=AY(KS)
            AZ(KST)=AZ(KS)
            A0(KST)=A0(KS)
            KPLANE(KST)=KPLANE(KS)
            KM(KST)=KS
C
            KBI=1
 701        CONTINUE
            DO KB=KBI,NBODYS
              KSL=0
              KSLP=0
              NSB=KSURF(KB,NXG)
              DO K=1,NSB
                IF(KSURF(KB,K).EQ.KS) KSL=K
                IF(KSURF(KB,K).EQ.KST) KSLP=K
              ENDDO
              IF(KSLP.GT.0) THEN
                KSURF(KB,KSLP)=KS
C  ****  If the implicit equations of surfaces KS and KST differ by a
C        global sign, the side pointer of KST must be reversed.
                IF(ISPF.EQ.1) THEN
                  IF(KFLAG(KB,KSLP).EQ.1) THEN
                    KFLAG(KB,KSLP)=2
                  ELSE IF(KFLAG(KB,KSLP).EQ.2) THEN
                    KFLAG(KB,KSLP)=1
                  ENDIF
                ENDIF
                IF(KSL.GT.0) THEN
                  KFL=KFLAG(KB,KSL)
                  KFLP=KFLAG(KB,KSLP)
                  IF(MIN(KFL,KFLP).LT.3) THEN
                    KFLAG(KB,KSL)=MIN(KFL,KFLP)
C
                    IF(KFL.NE.KFLP.AND.MAX(KFL,KFLP).LT.3) THEN
                      IF(KBOMO(KB).EQ.0) THEN
                      WRITE(IW,'(''*** ERROR: BODY('',I4,'') is li'',
     1  ''mited by two equivalent surfaces. Probably, '',/11X,''this '',
     2  ''body cannot be resolved because it is small and located'',
     3  /11X,''far from the origin.'')') KB
                      ELSE
                      WRITE(IW,'(''*** ERROR: MODULE('',I4,'') is li'',
     1  ''mited by two equivalent surfaces. Probably, '',/11X,''this '',
     2  ''module cannot be resolved because it is small and located'',
     3  /11X,''far from the origin.'')') KB
                      ENDIF
                      STOP
                    ENDIF
C
                  ELSE IF(KBOMO(KB).EQ.0) THEN
                     KFLAG(KB,KSL)=3
                  ELSE IF(KBOMO(KB).EQ.1) THEN
                     KFLAG(KB,KSL)=4
                  ELSE
                     KFLAG(KB,KSL)=5
                  ENDIF
                  DO K=KSLP,NSB-1
                    KSURF(KB,K)=KSURF(KB,K+1)
                    KFLAG(KB,K)=KFLAG(KB,K+1)
                  ENDDO
                  KSURF(KB,NSB)=0
                  KFLAG(KB,NSB)=5
                  KSURF(KB,NXG)=NSB-1
                  IF(KB.LT.NBODYS) THEN
                    KBI=KB
                    GOTO 701
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
 702    CONTINUE
 703  CONTINUE
 704  CONTINUE
C
      WRITE(IW,'(//,''************  Genealogical tree. '',/)')
      DO KB=1,NBODYS
        IF(KBOMO(KB).EQ.0) THEN
          WRITE(IW,'(/,''*** BODY   ='',I5,'',  KMOTH ='',
     1      I5,'',  MAT ='',I3)') KB,KMOTH(KB),MATER(KB)
          IF(KBODY(KB,NXG).GT.0)
     1    WRITE(IW,'(''KBODY ='',15I5,100(/,7X,15I5))')
     2      (KBODY(KB,K2),K2=1,KBODY(KB,NXG))
        ELSE IF(KBOMO(KB).EQ.1) THEN
          WRITE(IW,'(/,''*** MODULE ='',I5,'',  KMOTH ='',
     1      I5,'',  MAT ='',I3)') KB,KMOTH(KB),MATER(KB)
          WRITE(IW,'(''KDGHT ='',15I5,100(/,7X,15I5))')
     1      (KDGHT(KB,K2),K2=1,KDGHT(KB,NXG))
        ELSE
          WRITE(IW,'(//,''*** ERROR: the label '',I5,
     1      '' does not correspond to a body.'')') KB
          STOP
        ENDIF
        WRITE(IW,'(''KSURF ='',15I5,100(/,7X,15I5))')
     1    (KSURF(KB,KS),KS=1,KSURF(KB,NXG))
        WRITE(IW,'(''KFLAG ='',15I5,100(/,7X,15I5))')
     1    (KFLAG(KB,KS),KS=1,KSURF(KB,NXG))
      ENDDO
      IF(MLESS.EQ.1) THEN
        IF(KBOMO(KBENC).EQ.1) THEN
          WRITE(IW,'(/,''The module'',I5,'' is the enclosure.'')')
     1      KBENC
        ENDIF
      ENDIF
C
C  ****  Surface consistency test (F. Tola).
C
      DO KB=1,NBODYS-1
        KB1=KMOTH(KB)
        DO I=1,KSURF(KB,NXG)
          KS=KSURF(KB,I)
          KF=KFLAG(KB,I)
          DO J=1,KSURF(KB1,NXG)
            IF(KSURF(KB1,J).EQ.KS) THEN
C ---- Surface has daughter and mother at opposite sides.
              IF((KF.EQ.1.AND.KFLAG(KB1,J).EQ.2).OR.
     1          (KF.EQ.2.AND.KFLAG(KB1,J).EQ.1)) THEN
                IF(KBOMO(KB).EQ.0) THEN
                  WRITE(IW,'(//,''*** ERROR: the SURFACE ('',I4,''),'',
     1              '' which limits BODY ('', I4,'') and MODULE ('',
     2              I4,''),'',/11X,''has inconsistent side pointers.''
     3              )') KS,KB,KB1
                ELSE
                  WRITE(IW,'(//,''*** ERROR: the SURFACE ('',I4,''),'',
     1              '' which limits MODULE ('', I4,'') and MODULE ('',
     2              I4,''),'',/11X,''has inconsistent side pointers.''
     3              )') KS,KB,KB1
                ENDIF
                STOP
              ENDIF
              GOTO 801
            ENDIF
          ENDDO
 801      CONTINUE
        ENDDO
      ENDDO
C
C  ****  Easiness test.
C
      NBU=0
      NSU=0
      DO KB=1,NBODYS
        IF(KBOMO(KB).EQ.0) THEN
          NBU=MAX(NBU,KBODY(KB,NXG))
        ELSE IF(KBOMO(KB).EQ.1) THEN
          NBU=MAX(NBU,KDGHT(KB,NXG))
        ENDIF
        NSE=0
        DO K=1,KSURF(KB,NXG)
          IF(KFLAG(KB,K).LT.5) NSE=NSE+1
        ENDDO
        NSU=MAX(NSU,NSE)
      ENDDO
      WRITE(IW,'(//,''************  Adequacy of the geometry defin'',
     1  ''ition.'')')
      WRITE(IW,'(/,''The largest number of bodies in a module or'')')
      WRITE(IW,'(''     bodies limiting a single body is ........'',
     1  ''.... '',I4)') NBU
      WRITE(IW,'(/,''The largest number of limiting surfaces for'')')
      WRITE(IW,'(''     a single body or module is ...............'',
     1  ''... '',I4,/)') NSU
      IF(NBODYS.LT.15.AND.NSURF.LT.15) THEN
        WRITE(IW,'(/,''The simulation of this geometry will be rela'',
     1    ''tively fast,'')')
        WRITE(IW,'(''     no further optimization seems to be re'',
     1    ''quired.'')')
      ELSE IF(NBU.LT.10.AND.NSU.LT.10) THEN
        WRITE(IW,'(/,''The simulation of this geometry will be rela'',
     1    ''tively fast,'')')
        WRITE(IW,'(''      no further optimization seems to be re'',
     1    ''quired.'')')
      ELSE IF(NBU.LT.15.AND.NSU.LT.20) THEN
        WRITE(IW,'(/,''The simulation of this geometry is expected '',
     1    ''to be slow,'')')
        WRITE(IW,'(''     try to split complex bodies into seve'',
     1    ''ral modules.'')')
      ELSE IF(NBU.LT.25.AND.NSU.LT.30) THEN
        WRITE(IW,'(/,''The simulation of this geometry will be ve'',
     1    ''ry slow, you should'')')
        WRITE(IW,'(''     try to optimize the structure of the '',
     1    ''tree of modules.'')')
      ELSE
        WRITE(IW,'(/,''Simulating this geometry will be extremely '',
     1    ''slow.'')')
      ENDIF
      WRITE(IW,'(/''************  The end.'')')
      RETURN
C
  900 CONTINUE
      WRITE(IW,'(A72)') BLINE
      WRITE(IW,'(''*** Wrong input format.'')')
      STOP
      END
C  *********************************************************************
C                       SUBROUTINE ROTSHF
C  *********************************************************************
      SUBROUTINE ROTSHF(OMEGA,THETA,PHI,DX,DY,DZ,
     1                  AXX,AXY,AXZ,AYY,AYZ,AZZ,AX,AY,AZ,A0)
C
C     This subroutine rotates and shifts a quadric surface.
C
C  Input parameters:
C     OMEGA, THETA, PHI ... Euler rotation angles,
C     DX, DY, DZ .......... components of the displacement vector,
C     AXX, ..., A0 ........ coefficients of the initial quadric.
C
C  Output parameters:
C     AXX, ..., A0 ........ coefficients of the transformed quadric.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(3,3),A2(3,3),B2(3,3),A1(3),B1(3),D1(3)
C
C  ****  Initial quadric.
C
      B2(1,1)=AXX
      B2(1,2)=0.5D0*AXY
      B2(1,3)=0.5D0*AXZ
      B2(2,1)=B2(1,2)
      B2(2,2)=AYY
      B2(2,3)=0.5D0*AYZ
      B2(3,1)=B2(1,3)
      B2(3,2)=B2(2,3)
      B2(3,3)=AZZ
      B1(1)=AX
      B1(2)=AY
      B1(3)=AZ
      B0=A0
      D1(1)=DX
      D1(2)=DY
      D1(3)=DZ
C
C  ****  Rotation matrix.
C
      STHETA=SIN(THETA)
      CTHETA=COS(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      SOMEGA=SIN(OMEGA)
      COMEGA=COS(OMEGA)
C
      R(1,1)=CPHI*CTHETA*COMEGA-SPHI*SOMEGA
      R(1,2)=-CPHI*CTHETA*SOMEGA-SPHI*COMEGA
      R(1,3)=CPHI*STHETA
      R(2,1)=SPHI*CTHETA*COMEGA+CPHI*SOMEGA
      R(2,2)=-SPHI*CTHETA*SOMEGA+CPHI*COMEGA
      R(2,3)=SPHI*STHETA
      R(3,1)=-STHETA*COMEGA
      R(3,2)=STHETA*SOMEGA
      R(3,3)=CTHETA
C
C  ****  Rotated quadric.
C
      DO I=1,3
        A1(I)=0.0D0
        DO J=1,3
          A1(I)=A1(I)+R(I,J)*B1(J)
          A2(I,J)=0.0D0
          DO M=1,3
            DO K=1,3
              A2(I,J)=A2(I,J)+R(I,K)*B2(K,M)*R(J,M)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C  ****  Shifted-rotated quadric.
C
      DO I=1,3
        A2D=0.0D0
        DO J=1,3
          A2D=A2D+A2(I,J)*D1(J)
        ENDDO
        B1(I)=A1(I)-2.0D0*A2D
        B0=B0+D1(I)*(A2D-A1(I))
      ENDDO
C
      AXX=A2(1,1)
      AXY=A2(1,2)+A2(2,1)
      AXZ=A2(1,3)+A2(3,1)
      AYY=A2(2,2)
      AYZ=A2(2,3)+A2(3,2)
      AZZ=A2(3,3)
      AX=B1(1)
      AY=B1(2)
      AZ=B1(3)
      A0=B0
      IF(ABS(AXX).LT.1.0D-16) AXX=0.0D0
      IF(ABS(AXY).LT.1.0D-16) AXY=0.0D0
      IF(ABS(AXZ).LT.1.0D-16) AXZ=0.0D0
      IF(ABS(AYY).LT.1.0D-16) AYY=0.0D0
      IF(ABS(AYZ).LT.1.0D-16) AYZ=0.0D0
      IF(ABS(AZZ).LT.1.0D-16) AZZ=0.0D0
      IF(ABS(AX).LT.1.0D-16) AX=0.0D0
      IF(ABS(AY).LT.1.0D-16) AY=0.0D0
      IF(ABS(AZ).LT.1.0D-16) AZ=0.0D0
      IF(ABS(A0).LT.1.0D-16) A0=0.0D0
      RETURN
      END
C  *********************************************************************
C                       FUNCTION FSURF
C  *********************************************************************
      SUBROUTINE FSURF(KS,A,B,C)
C
C     Calculates the parameters of the master function of the surface KS
C  and the ray (X,Y,Z)+S*(U,V,W).
C
      USE PENGEOM_mod
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
C
      IF(KPLANE(KS).EQ.0) THEN
        A=U*(AXX(KS)*U+AXY(KS)*V+AXZ(KS)*W)
     1   +V*(AYY(KS)*V+AYZ(KS)*W)+W*AZZ(KS)*W
        XXX=AXX(KS)*X+AXY(KS)*Y+AXZ(KS)*Z+AX(KS)
        YYY=AYY(KS)*Y+AYZ(KS)*Z+AY(KS)
        ZZZ=AZZ(KS)*Z+AZ(KS)
        B=U*(AXX(KS)*X+XXX)+V*(AXY(KS)*X+AYY(KS)*Y+YYY)
     1   +W*(AXZ(KS)*X+AYZ(KS)*Y+AZZ(KS)*Z+ZZZ)
        C=X*XXX+Y*YYY+Z*ZZZ+A0(KS)
      ELSE
        A=0.0D0
        B=U*AX(KS)+V*AY(KS)+W*AZ(KS)
        C=X*AX(KS)+Y*AY(KS)+Z*AZ(KS)+A0(KS)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE LOCATE
C  *********************************************************************
      SUBROUTINE LOCATE
C
C     This subroutine determines the body that contains the point with
C  coordinates (X,Y,Z). The effects of numerical round-off errors are
C  avoided by considering fuzzy surfaces, which swell or shrink slightly
C  when the particle crosses them.
C
C  Input values (module TRACK_mod):
C     X, Y, Z ... coordinates of the particle,
C     U, V, W ... direction of movement.
C
C  Output values (module TRACK_mod):
C     IBODY ..... body where the particle moves,
C     MAT  ...... material in IBODY,
C                    MAT=0, indicates a void region.
C
      USE PENGEOM_mod
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (FUZZL=1.0D-12)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
C
      DO I=1,NSURF
        KSP(I)=0
      ENDDO
      KB0=NBODYS
 100  CONTINUE
      DO 101 KSS=1,KSURF(KB0,NXG)
        KS=KSURF(KB0,KSS)
        IF(KSP(KS).NE.0.OR.KFLAG(KB0,KSS).GT.4) GOTO 101
C
        CALL FSURF(KS,A,B,C)
        ABSA=ABS(A)
        IF(ABSA.GT.1.0D-36) THEN
          FUZZ=FUZZL*(B*B-4.0D0*A*C)/ABSA
        ELSE
          FUZZ=FUZZL*ABS(B)
        ENDIF
C
        IF(C.LT.-FUZZ) THEN
          KSP(KS)=1
        ELSE IF(C.GT.FUZZ) THEN
          KSP(KS)=2
        ELSE
C  ****  Point close to the surface.
          IF(B.LT.0.0D0) THEN
            KSP(KS)=1  ! Particle moving 'inwards'.
          ELSE
            KSP(KS)=2  ! Particle moving 'outwards'.
          ENDIF
        ENDIF
 101  CONTINUE
C
C  ****  Determine the module or body that contains the point.
C
      DO 102 KBB=1,KDGHT(KB0,NXG)
        KB=KDGHT(KB0,KBB)
        DO KSS=1,KSURF(KB,NXG)
          KS=KSURF(KB,KSS)
          KF=KFLAG(KB,KSS)
          IF(KF.LT.3.AND.KSP(KS).NE.KF) GOTO 102
        ENDDO
        IF(KB.EQ.KB0) THEN
          IBODY=KB  ! The particle is inside the body or module KB.
          MAT=MATER(KB)
          RETURN
        ELSE IF(KDGHT(KB,NXG).GT.1) THEN
          KB0=KB    ! The point is inside a submodule.
          GOTO 100
        ELSE
          IBODY=KB  ! The particle is inside a sister body or module.
          MAT=MATER(KB)
          RETURN
        ENDIF
 102  CONTINUE
      IBODY=NBODYS+1  ! The particle is outside the enclosure.
      MAT=0
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEPSI
C  *********************************************************************
      SUBROUTINE STEPSI(KB,S,IS,NSC)
C
C     Calculates the intersections of the trajectory with the limiting
C  surfaces of body KB. The intersections are added to the list and
C  sorted in decreasing order. This subroutine works only when called
C  from within subroutine STEP.
C
      USE PENGEOM_mod
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (FUZZL=1.0D-12)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
      PARAMETER (NS2M=2*NS)
      DIMENSION S(NS2M),IS(NS2M)
C
C  ************  Determine surface crossings.
C
      DO 100 KSS=1,KSURF(KB,NXG)
C  ****  Intersections with a given surface are calculated only once.
C        The side pointer of a surface must be changed each time the
C        surface is crossed.
        KFL=KFLAG(KB,KSS)
        IF(KFL.GT.4) GOTO 100
        KS=KSURF(KB,KSS)
        IF(KSP(KS).NE.0) GOTO 100
        CALL FSURF(KS,A,B,C)
        ABSA=ABS(A)
        ABSB=ABS(B)
C
C  ****  Plane, a single root.
C
        IF(ABSA.LT.1.0D-36) THEN
          IF(ABSB.GT.0.0D0) THEN
            IF(C.LT.-FUZZL) THEN  ! SP=-1
              KSP(KS)=1
            ELSE IF(C.GT.FUZZL) THEN  ! SP=+1
              KSP(KS)=2
            ELSE  ! Point close to the surface.
              IF(B.LT.0.0D0) THEN
                KSP(KS)=1  ! Particle moving 'inwards'.
              ELSE
                KSP(KS)=2  ! Particle moving 'outwards'.
              ENDIF
              GO TO 100
            ENDIF
            T1=-C/B
            IF(T1.GT.0.0D0) THEN
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=T1
            ENDIF
          ELSE  ! Ray parallel to the plane.
            IF(C.LT.0.0D0) THEN
              KSP(KS)=1
            ELSE
              KSP(KS)=2
            ENDIF
          ENDIF
C
C  ****  Non-planar surface, two roots.
C
        ELSE
          DISCR=B*B-4.0D0*A*C
          FUZZ=FUZZL*DISCR/ABSA
          IF(C.LT.-FUZZ) THEN  ! SP=-1
            IAMBIG=0
            KSP(KS)=1
          ELSE IF(C.GT.FUZZ) THEN  ! SP=+1
            IAMBIG=0
            KSP(KS)=2
          ELSE
            IAMBIG=1  ! Point close to the surface.
            IF(B.LT.0.0D0) THEN
              KSP(KS)=1  ! Particle moving 'inwards'.
            ELSE
              KSP(KS)=2  ! Particle moving 'outwards'.
            ENDIF
          ENDIF
C
          IF(DISCR.LT.1.0D-36) GO TO 100  ! No true intersections.
C
          IF(IAMBIG.EQ.0) THEN
            R2A=0.5D0/A
            DELTA=SQRT(DISCR)*ABS(R2A)
            SH=-B*R2A
            T1=SH-DELTA
            IF(T1.GT.0.0D0) THEN
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=T1
            ENDIF
            T2=SH+DELTA
            IF(T2.GT.0.0D0) THEN
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=T2
            ENDIF
          ELSE
            IF(B*A.LT.0.0D0) THEN
              R2A=0.5D0/A
              DELTA=SQRT(DISCR)*ABS(R2A)
              SH=-B*R2A
              T2=SH+DELTA
              NSC=NSC+1
              IS(NSC)=KS
              S(NSC)=MAX(T2,0.0D0)
            ENDIF
          ENDIF
        ENDIF
 100  CONTINUE
C
C  ****  Sort surface distances in decreasing order.
C
      IF(NSC.GT.1) THEN
        DO KI=1,NSC-1
          SMAX=S(KI)
          KMAX=KI
          DO KJ=KI+1,NSC
            IF(S(KJ).GT.SMAX) THEN
              SMAX=S(KJ)
              KMAX=KJ
            ENDIF
          ENDDO
          IF(KMAX.NE.KI) THEN
            SMAX=S(KI)
            S(KI)=S(KMAX)
            S(KMAX)=SMAX
            KKMAX=IS(KI)
            IS(KI)=IS(KMAX)
            IS(KMAX)=KKMAX
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEPLB
C  *********************************************************************
      SUBROUTINE STEPLB(KB,IERR)
C
C     Helps finding the body or module that has the given side pointers
C  for the analyzed surfaces. Subroutine STEPLB works only when invoked
C  from within subroutine STEP. It moves through the tree of modules a
C  single step.
C
      USE PENGEOM_mod
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QBODY/KBODY(NB,NXG),KBOMO(NB)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
C
C  ****  Analyze the current body or module.
C
      IF(KBOMO(KB).EQ.0) THEN
C  ****  Body.
        NLBOD=KBODY(KB,NXG)
        IF(NLBOD.GT.0) THEN
          DO 100 KBB=1,NLBOD
            KBS=KBODY(KB,KBB)
            DO KSS=1,KSURF(KBS,NXG)
              KS=KSURF(KBS,KSS)
              KF=KFLAG(KBS,KSS)
              IF(KF.LT.3.AND.KSP(KS).NE.KF) GOTO 100
            ENDDO
            IBODY=KBS
            IF(KDGHT(IBODY,NXG).GT.1) THEN
              IERR=-1  ! The particle is inside a sister module.
            ELSE
              IERR=0   ! The particle is inside a sister body.
              MAT=MATER(IBODY)
            ENDIF
            RETURN
 100      CONTINUE
        ENDIF
        DO KSS=1,KSURF(KB,NXG)
          KS=KSURF(KB,KSS)
          KF=KFLAG(KB,KSS)
          IF(KF.LT.3.AND.KSP(KS).NE.KF) GOTO 300
        ENDDO
        IBODY=KB
        IERR=0   ! The particle remains in the same body.
        MAT=MATER(IBODY)
        RETURN
      ELSE
C  ****  Module.
        DO 200 KBB=1,KDGHT(KB,NXG)
          KBD=KDGHT(KB,KBB)
          DO KSS=1,KSURF(KBD,NXG)
            KS=KSURF(KBD,KSS)
            KF=KFLAG(KBD,KSS)
            IF(KF.LT.3.AND.KSP(KS).NE.KF) GOTO 200
          ENDDO
          IBODY=KBD
          IF(KBD.EQ.KB) THEN
            IERR=0  !  The particle remains within the current module.
            MAT=MATER(IBODY)
          ELSE
            IF(KDGHT(KBD,NXG).GT.1) THEN
              IERR=-1  ! The particle is inside a submodule.
            ELSE
              IERR=0   ! The particle is inside a single body.
              MAT=MATER(IBODY)
            ENDIF
          ENDIF
          RETURN
 200    CONTINUE
      ENDIF
C
C  ****  The particle is outside the current body or module.
C
 300  CONTINUE
      IERR=1
      IBODY=KMOTH(KB)
      IF(IBODY.EQ.0) THEN
        IBODY=NBODYS+1
        MAT=0
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE STEP
C  *********************************************************************
      SUBROUTINE STEP(DS,DSEF,NCROSS)
C
C     This subroutine handles the geometrical part of the track simula-
C  tion. The particle starts from the point (X,Y,Z) and travels a length
C  DS in the direction (U,V,W) within the material where it moves. When
C  the track leaves the initial material, the particle is stopped just
C  after entering the next material body (void regions with MAT=0 are
C  crossed automatically). Furthermore, when the particle arrives from
C  a void region, it is stopped just after entering the first material
C  body.
C
C  Input values (module TRACK_mod):
C     X, Y, Z ... coordinates of the initial point,
C     U, V, W ... direction cosines of the displacement,
C     IBODY ..... body where the initial point is located,
C     MAT ....... material in body IBODY.
C  NB: When a particle track is started, the variables IBODY and MAT
C  must be set by calling subroutine LOCATE.
C
C  Input argument:
C     DS ........ path length to travel.
C
C  Output arguments:
C     DSEF....... travelled path length before leaving the initial
C                 material or completing the jump (less than DS if the
C                 track crosses an interface),
C     NCROSS .... = 0 if the whole step is contained in the initial
C                   material,
C                 .gt.0 if the particle has crossed an interface, i.e.
C                   if it has entered a new material.
C
C  Output values (module TRACK_mod):
C     X, Y, Z ... coordinates of the final position,
C     IBODY ..... body where the final point is located,
C     MAT ....... material in IBODY. The value MAT=0 indicates that the
C                 particle has escaped from the system.
C
C  Output values (module PENGEOM_mod):
C     DSTOT ..... travelled path length, including path segments in void
C                 volumes.
C     KSLAST .... when NCROSS.ne.0, the output value of KSLAST is the
C        label of the last surface crossed by the particle before
C        entering a material body. KSLAST is used for rendering in 3D.
C
      USE TRACK_mod
      USE PENGEOM_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/QSURF/AXX(NS),AXY(NS),AXZ(NS),AYY(NS),AYZ(NS),AZZ(NS),
     1    AX(NS),AY(NS),AZ(NS),A0(NS),NSURF,KPLANE(NS)
      COMMON/QTREE/NBODYS,KMOTH(NB),KDGHT(NB,NXG),KSURF(NB,NXG),
     1  KFLAG(NB,NXG),KSP(NS),NWARN
      PARAMETER (NS2M=2*NS)
      DIMENSION S(NS2M),IS(NS2M)
C
      DSEF=0.0D0
      DSTOT=0.0D0
      NCROSS=0
      KSLAST=0
C
      NSC=0        ! Number of surface crossings ahead of the particle.
      DO I=1,NSURF
        KSP(I)=0   ! Side pointers of the evaluated surfaces.
      ENDDO
      MAT0=MAT     ! Initial material.
      IF(MAT.EQ.0) THEN
        DSRES=1.0D35  ! In vacuum particles fly freely.
      ELSE
        DSRES=DS   ! Residual path length.
      ENDIF
C
C  ************  The particle enters from outside the enclosure.
C
      IF(IBODY.GT.NBODYS) THEN
        KB1=NBODYS
        CALL STEPSI(KB1,S,IS,NSC)
        IF(NSC.EQ.0) GOTO 300
        NSCT=NSC
        NST=KSURF(KB1,NXG)
        DO 101 KI=NSCT,1,-1
C  ****  The particle crosses a surface.
          KSLAST=IS(KI)
          IF(KSP(KSLAST).EQ.1) THEN
            KSP(KSLAST)=2
          ELSE
            KSP(KSLAST)=1
          ENDIF
          DSP=S(KI)
          DSEF=DSEF+DSP
          DSTOT=DSTOT+DSP
          X=X+DSP*U
          Y=Y+DSP*V
          Z=Z+DSP*W
          NSC=NSC-1
          IF(NSC.GT.0) THEN
            DO I=1,NSC
              S(I)=S(I)-DSP
            ENDDO
          ENDIF
C
          DO KSS=1,NST
            KS1=KSURF(KB1,KSS)
            KF=KFLAG(KB1,KSS)
            IF(KF.LT.3.AND.KSP(KS1).NE.KF) GOTO 101
          ENDDO
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C  ****  The particle enters the enclosure.
 100      CONTINUE
          CALL STEPLB(KB1,IERR)
C  ****  The particle enters a submodule.
          IF(IERR.EQ.-1) THEN
            KB1=IBODY
            CALL STEPSI(KB1,S,IS,NSC)
            GOTO 100
          ELSE
C  ****  The particle enters a material body.
            IF(MAT.NE.0) THEN
              NCROSS=1
              RETURN
            ELSE
              KB1=IBODY
              CALL STEPSI(KB1,S,IS,NSC)
              GOTO 200
            ENDIF
          ENDIF
C  ****  At this point the program has left the DO loop.
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 101    CONTINUE
        GOTO 300
      ENDIF
C
C  ************  Surface crossings.
C
      IBODYL=IBODY
      IF(LVERB) NERR=0
 102  CONTINUE
      KB1=IBODY
      CALL STEPSI(KB1,S,IS,NSC)
      CALL STEPLB(KB1,IERR)
C
C  ****  Evidence of round-off errors.
C
      IF(IERR.NE.0) THEN
        IF(NSC.GT.0) THEN
C  ****  When a surface is very close, we move the particle beyond it.
          IF(S(NSC).LT.1.0D-10) THEN
            KSLAST=IS(NSC)
            IF(KSP(KSLAST).EQ.1) THEN
              KSP(KSLAST)=2
            ELSE
              KSP(KSLAST)=1
            ENDIF
            DSP=S(NSC)
            X=X+DSP*U
            Y=Y+DSP*V
            Z=Z+DSP*W
            IF(MAT.EQ.MAT0) THEN
              DSEF=DSEF+DSP
              DSRES=DSRES-DSP
            ENDIF
            DSTOT=DSTOT+DSP
            NSC=NSC-1
            IF(IBODY.LE.NBODYS) GOTO 102
          ENDIF
        ENDIF
        IF(LVERB) THEN
          NERR=NERR+1
          IF(NWARN.LT.100.AND.MAT.NE.0) THEN
          WRITE(26,'(''*** WARNING, STEP: Accidental undershot or r'',
     1      ''ound-off error?'',/4X,''IBODY0 ='',I5,'',  IBODY ='',I5,
     2      '',  IERR ='',I3,'',  NERR ='',I4)') KB1,IBODY,IERR,NERR
          WRITE(26,'(''(X,Y,Z) = '',1P,E15.7,'', '',E15.7,'', '',
     1      E15.7)') X,Y,Z
          WRITE(26,'(''(U,V,W) = '',1P,E15.7,'', '',E15.7,'', '',
     1      E15.7)') U,V,W
          DO KSS=1,KSURF(KB1,NXG)
            KS=KSURF(KB1,KSS)
            KFLO=KFLAG(KB1,KSS)
            IF(KFLO.LT.3) THEN
              DO KI=NSC,1,-1
                IF(KS.EQ.IS(KI)) THEN
                  SW=S(KI)
                  GO TO 103
                ENDIF
              ENDDO
              SW=0.0D0
 103          CONTINUE
              CALL FSURF(KS,A,B,C)
              IF(KFLO.EQ.KSP(KS)) THEN
                WRITE(26,'(A,I5,2X,2I3,1P,2E15.7)')
     1            'Surface,FL0,FL,C,S =',KS,KFLO,KSP(KS),C,SW
              ELSE
                WRITE(26,'(A,I5,2X,2I3,1P,2E15.7,A)')
     1            'Surface,FL0,FL,C,S =',KS,KFLO,KSP(KS),C,SW,' *'
                KSLAST=KS
              ENDIF
            ENDIF
          ENDDO
          NWARN=NWARN+1
          ENDIF
        ENDIF
        IF(IBODY.LE.NBODYS) GOTO 102
      ENDIF
C
      IF(KDET(IBODY).NE.KDET(IBODYL).OR.MAT.NE.MAT0) THEN
        NCROSS=1
        DSEF=0.0D0
        RETURN
      ENDIF
C
C  ****  The particle remains in the same material.
C
      IF(MAT.NE.0.AND.DSRES.LT.S(NSC)) THEN
        IF(MAT.EQ.MAT0) DSEF=DSEF+DSRES
        DSTOT=DSTOT+DSRES
        X=X+DSRES*U
        Y=Y+DSRES*V
        Z=Z+DSRES*W
        RETURN
      ENDIF
C
C  ************  New position.
C
 200  CONTINUE
      IF(NSC.EQ.0) THEN
        IF(MAT.EQ.MAT0) DSEF=DSEF+DSRES
        DSTOT=DSTOT+DSRES
        X=X+DSRES*U
        Y=Y+DSRES*V
        Z=Z+DSRES*W
        RETURN
      ENDIF
      NSCT=NSC
      MATL=MAT
      IBODYL=IBODY
      DO 203 KI=NSCT,1,-1
C  ****  The step ends within the body.
        IF(DSRES.LT.S(KI)) THEN
          IF(MAT.EQ.MAT0) DSEF=DSEF+DSRES
          DSTOT=DSTOT+DSRES
          X=X+DSRES*U
          Y=Y+DSRES*V
          Z=Z+DSRES*W
          RETURN
        ENDIF
C  ****  The particle crosses a surface.
        KSLAST=IS(KI)
        IF(KSP(KSLAST).EQ.1) THEN
          KSP(KSLAST)=2
        ELSE
          KSP(KSLAST)=1
        ENDIF
        DSP=S(KI)
        X=X+DSP*U
        Y=Y+DSP*V
        Z=Z+DSP*W
        IF(MAT.EQ.MAT0) THEN
          DSEF=DSEF+DSP
          DSRES=DSRES-DSP
        ENDIF
        DSTOT=DSTOT+DSP
        NSC=NSC-1
        IF(NSC.GT.0) THEN
          DO I=1,NSC
            S(I)=S(I)-DSP
          ENDDO
        ENDIF
        CALL STEPLB(KB1,IERR)
C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 201    CONTINUE
        KB1=IBODY
        IF(IERR.EQ.-1) THEN
C  ****  The particle enters a submodule.
          CALL STEPSI(KB1,S,IS,NSC)
          CALL STEPLB(KB1,IERR)
          GOTO 201
        ELSE IF(IERR.EQ.1) THEN
C  ****  The particle leaves the body or module.
          IF(IBODY.LE.NBODYS) THEN
            CALL STEPSI(KB1,S,IS,NSC)
            CALL STEPLB(KB1,IERR)
            GOTO 201
          ELSE
C  ****  The particle leaves the enclosure.
            IF(MAT.NE.MATL) NCROSS=NCROSS+1
            GOTO 300
          ENDIF
        ENDIF
C  ****  The particle continues flying when it enters a void region...
        IF(MAT.EQ.0) THEN
          IF(MATL.EQ.MAT0) NCROSS=NCROSS+1
          MATL=0
          DSRES=1.0D35
          GOTO 202
C  ****  The particle continues flying when it enters a new body of the
C        same material which is not part of a different detector...
        ELSE IF(MAT.EQ.MATL) THEN
          IF(KDET(IBODY).EQ.KDET(IBODYL)) THEN
            GOTO 202
          ELSE
            NCROSS=NCROSS+1
            RETURN
          ENDIF
C  ****  ... and stops when it penetrates a new material body or a
C        detector.
        ELSE   ! IF(MAT.NE.MATL) THEN
          NCROSS=NCROSS+1
          RETURN
        ENDIF
 202    CONTINUE
        CALL STEPSI(KB1,S,IS,NSC)
        GOTO 200
C  ****  At this point the program has left the DO loop.
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 203  CONTINUE
C
C  ************  The particle leaves the enclosure.
C
 300  CONTINUE
      DSP=1.0D36
      IBODY=NBODYS+1
      MAT=0
      IF(MAT.EQ.MAT0) DSEF=DSEF+DSP
      DSTOT=DSTOT+DSP
      X=X+DSP*U
      Y=Y+DSP*V
      Z=Z+DSP*W
      RETURN
      END
