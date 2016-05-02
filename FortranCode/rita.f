CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  Subroutine package RITA            (Francesc Salvat. 8 June, 2014)  C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     This package contains subroutines for random sampling from single-
C  variate discrete and continuous probability distributions. Discrete
C  distributions are sampled by using the aliasing method of Walker. The
C  sampling from continuous distributions is performed by means of the
C  RITA (Rational Inverse Transform with Aliasing) algorithm. These
C  methods are described in Chapter 1 of the PENELOPE manual; they are
C  among the fastest sampling techniques available for arbitrary
C  numerical distributions.
C
C  *********************************************************************
C                       FUNCTION IRND
C  *********************************************************************
      FUNCTION IRND(FA,IA,N)
C
C  Random sampling from a discrete probability distribution using
C  Walker's aliasing algorithm.
C
C  The arrays F and IA are determined by the initialisation routine
C  IRND0, which must be invoked before using function IRND.
C
C  Input arguments:
C    FA(1:N) ... cutoff values.
C    IA(1:N) ... alias values.
C    N ......... number of different values of the random variable.
C
C  Output argument:
C    IRND ...... sampled value.
C
C  Other subprograms needed: function RAND and subroutine IRND0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION FA(N),IA(N)
      EXTERNAL RAND
C
      RN=RAND(1.0D0)*N+1.0D0
      IRND=INT(RN)
      TST=RN-IRND
      IF(TST.GT.FA(IRND)) IRND=IA(IRND)
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RITA
C  *********************************************************************
      FUNCTION RITA()
C
C  Random sampling of a continuous variable using the RITA (Rational
C  Inverse Transform with Aliasing) method.
C
C  The needed numerical parameters are determined by the initialisa-
C  tion subroutine RITA0, which must be invoked before using function
C  RITA; these parameters are stored in common block /CRITAA/.
C
C  Other subprograms needed: EXTERNAL functions PDF and RAND,
C                            subroutines RITA0, IRND0 and RITAI0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NM=512)
      COMMON/CRITAA/XA(NM),AA(NM),BA(NM),FA(NM),IA(NM),NPM1A
      EXTERNAL RAND
C  ****  Selection of the interval (Walker's aliasing).
      RN=RAND(1.0D0)*NPM1A+1.0D0
      K=INT(RN)
      TST=RN-K
      IF(TST.LT.FA(K)) THEN
        I=K
        RR=TST
        D=FA(K)
      ELSE
        I=IA(K)
        RR=TST-FA(K)
        D=1.0D0-FA(K)
      ENDIF
C  ****  Sampling from the rational inverse cumulative distribution.
      IF(RR.GT.1.0D-12) THEN
        RITA=XA(I)+((1.0D0+AA(I)+BA(I))*D*RR
     1    /(D*D+(AA(I)*D+BA(I)*RR)*RR))*(XA(I+1)-XA(I))
      ELSE
        RITA=XA(I)+RAND(2.0D0)*(XA(I+1)-XA(I))
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       FUNCTION RITAI
C  *********************************************************************
      FUNCTION RITAI()
C
C  Random sampling of a continuous variable using the RITA (Rational
C  Inverse Transform with Aliasing) method, with binary search within
C  pre-calculated index intervals.
C
C  The needed numerical parameters are determined by the initialisation
C  subroutine RITAI0, which must be invoked before using function RITAI.
C  These parameters are stored in common block /CRITA/.
C
C  Other subprograms needed: EXTERNAL functions PDF and RAND,
C                            subroutine RITAI0.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NM=512)
      COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A(NM),B(NM),IL(NM),IU(NM),
     1  NPM1
      EXTERNAL RAND
C
C  ****  Selection of the interval
C        (binary search within pre-calculated limits).
C
      RU=RAND(1.0D0)
      ITN=RU*NPM1+1.0D0
      I=IL(ITN)
      J=IU(ITN)
      IF(J-I.LT.2) GO TO 2
    1 K=(I+J)/2
      IF(RU.GT.PAC(K)) THEN
        I=K
      ELSE
        J=K
      ENDIF
      IF(J-I.GT.1) GO TO 1
C
C  ****  Sampling from the rational inverse cumulative distribution.
C
    2 CONTINUE
      RR=RU-PAC(I)
      D=DPAC(I)  ! DPAC(I)=PAC(I+1)-PAC(I)
      IF(D.GT.1.0D-12) THEN
        RITAI=XT(I)+((1.0D0+A(I)+B(I))*D*RR/(D*D+(A(I)*D+B(I)*RR)*RR))
     1       *(XT(I+1)-XT(I))
      ELSE
        RITAI=XT(I)+RAND(2.0D0)*(XT(I+1)-XT(I))
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE IRND0
C  *********************************************************************
      SUBROUTINE IRND0(W,F,K,N)
C
C  Initialisation of Walker's aliasing algorithm for random sampling
C  from discrete probability distributions.
C
C  Input arguments:
C    N ........ number of different values of the random variable.
C    W(1:N) ... corresponding point probabilities (not necessarily
C               normalised to unity).
C  Output arguments:
C    F(1:N) ... cutoff values.
C    K(1:N) ... alias values.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION W(N),F(N),K(N)
      COMMON/CRITAN/CNORM  ! Normalising constant, output.
C  ****  Renormalisation.
      CNORM=0.0D0
      DO I=1,N
        IF(W(I).LT.0.0D0) STOP 'IRND0. Negative point probability.'
        CNORM=CNORM+W(I)
      ENDDO
      CNORM=1.0D0/CNORM
      FACT=DBLE(N)*CNORM
      DO I=1,N
        K(I)=I
        F(I)=W(I)*FACT
      ENDDO
      IF(N.EQ.1) RETURN
C  ****  Cutoff and alias values.
      DO I=1,N-1
        HLOW=1.0D0
        HIGH=1.0D0
        ILOW=0
        IHIGH=0
        DO J=1,N
          IF(K(J).EQ.J) THEN
            IF(F(J).LT.HLOW) THEN
              HLOW=F(J)
              ILOW=J
            ELSE IF(F(J).GT.HIGH) THEN
              HIGH=F(J)
              IHIGH=J
            ENDIF
          ENDIF
        ENDDO
        IF(ILOW.EQ.0.OR.IHIGH.EQ.0) RETURN
        K(ILOW)=IHIGH
        F(IHIGH)=HIGH+HLOW-1.0D0
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RITA0
C  *********************************************************************
      SUBROUTINE RITA0(PDF,XLOW,XHIGH,N,NU,ERRM,IWR)
C
C  Initialisation of the RITA algorithm for random sampling of a
C  continuous random variable X from a probability distribution function
C  PDF(X) defined in the interval (XLOW,XHIGH). The external function
C  PDF(X) --not necessarily normalised-- must be provided by the user.
C  N is the number of points in the sampling grid. These points are
C  determined by means of an adaptive strategy that minimizes local
C  interpolation errors. The first NU grid points are uniformly spaced
C  in (XLOW,XHIGH); when NU is negative, the initial grid consists of
C  -NU points logarithmically spaced (in this case, XLOW must be
C  nonnegative).
C
C  ERRM is a measure of the interpolation error (the largest value of
C  the absolute error of the rational interpolation integrated over each
C  grid interval).
C
C  ****  Interpolation coefficients and PDF tables are printed on
C        separate files (UNIT=IWR) if IWR is greater than zero.
C
C  Other subprograms needed: EXTERNAL function PDF,
C                            subroutines RITAI0 and IRND0.
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (NM=512)
      COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A(NM),B(NM),IL(NM),IU(NM),
     1  NPM1
      COMMON/CRITAA/XA(NM),AA(NM),BA(NM),FA(NM),IA(NM),NPM1A
      COMMON/CRITAN/CNORM  ! Normalising constant, output.
      EXTERNAL PDF
C  ****  Initialisation of the RITA algorithm.
      CALL RITAI0(PDF,XLOW,XHIGH,N,NU,ERRM,IWR)
C  ****  Walker's aliasing; cutoff and alias values.
      NPM1A=NPM1
      DO I=1,NPM1
        XA(I)=XT(I)
        AA(I)=A(I)
        BA(I)=B(I)
      ENDDO
      XA(NPM1+1)=XT(NPM1+1)
      SAVE=CNORM
      CALL IRND0(DPAC,FA,IA,NPM1A)
      CNORM=SAVE
      FA(NPM1+1)=1.0D0
      IA(NPM1+1)=NPM1A
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RITAI0
C  *********************************************************************
      SUBROUTINE RITAI0(PDF,XLOW,XHIGH,N,NU,ERRM,IWR)
C
C  Initialisation of the RITA algorithm for random sampling of a
C  continuous random variable X from a probability distribution function
C  PDF(X) defined in the interval (XLOW,XHIGH). The external function
C  PDF(X) --not necessarily normalised-- must be provided by the user.
C  N is the number of points in the sampling grid. These points are
C  determined by means of an adaptive strategy that minimizes local
C  interpolation errors. The first NU grid points are uniformly spaced
C  in (XLOW,XHIGH); when NU is negative, the initial grid consists of
C  -NU points logarithmically spaced (in this case, XLOW must be
C  nonnegative).
C
C  ERRM is a measure of the interpolation error (the largest value of
C  the absolute error of the rational interpolation integrated over each
C  grid interval).
C
C  ****  Interpolation coefficients and PDF tables are printed on
C        separate files (UNIT=IWR) if IWR is greater than zero.
C
C  Other subprograms needed: EXTERNAL function PDF.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (EPS=1.0D-10, ZERO=1.0D-75, ZEROT=0.1D0*ZERO)
      PARAMETER (NM=512)
C
C  The information used by the sampling function RITAI is exported
C  through the following common block,
      COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A(NM),B(NM),IL(NM),IU(NM),
     1  NPM1
C  where
C    XT(I) ..... grid points, in increasing order.
C    PAC(I) .... value of the cumulative pdf at XT(I).
C    DPAC(I) ... probability of the I-th interval.
C    A(I), B(I) ... rational inverse cumulative distribution parameters.
C    IL(I) .... largest J for which PAC(J) < (I-1)/(NP-1).
C    IU(I) .... smallest J for which PAC(J) > I/(NP-1).
C    NPM1 ..... numner of grid points minus one, NP-1 (8.LE.NP.LE.NM).
C
      COMMON/CRITAN/CNORM
C    CNORM .... normalising constant of the external PDF, output.
C
      DIMENSION ERR(NM),C(NM)
      PARAMETER (NIP=51)
      DIMENSION XS(NIP),YS(NIP),SUMI(NIP)
      EXTERNAL PDF
C
      IF(N.LT.9) THEN
        WRITE(6,'('' Error in RITAI0: N must be larger than 8.'',
     1    /,'' N='',I11)') N
        STOP 'RITAI0: N must be larger than 8.'
      ENDIF
      IF(N.GT.NM) THEN
        WRITE(6,'('' Error in RITAI0: N must be less than NM=512.'',
     1    /,'' N='',I11)') N
        STOP 'RITAI0: N must be less than NM=512.'
      ENDIF
      IF(XLOW.GT.XHIGH-EPS) THEN
        WRITE(6,'('' Error in RITAI0: XLOW must be larger than XHIGH.'',
     1    /,'' XLOW='',1P,E13.6,'', XHIGH ='',E13.6)') XLOW,XHIGH
        STOP 'RITAI0: XLOW must be larger than XHIGH.'
      ENDIF
C
C  ****  We start with a grid of NUNIF points uniformly spaced in the
C        interval (XLOW,XHIGH).
C
      IF(NU.GE.0) THEN
        NUNIF=MIN(MAX(8,NU),N/2)
        NP=NUNIF
        DX=(XHIGH-XLOW)/DBLE(NP-1)
        XT(1)=XLOW
        DO I=1,NP-1
          XT(I+1)=XLOW+I*DX
        ENDDO
        XT(NP)=XHIGH
      ELSE
C  ****  If NU.LT.0,the NUNIF points are logarithmically spaced.
C        XLOW must be greater than or equal to zero.
        NUNIF=MIN(MAX(8,-NU),N/2)
        NP=NUNIF
        IF(XLOW.LT.0.0D0) THEN
          WRITE(6,'('' Error in RITAI0: XLOW and NU are negative.'',
     1    /,'' XLOW='',1P,E14.7, '',  NU='',I11)') XLOW,NU
          STOP 'RITAI0: XLOW and NU are negative.'
        ENDIF
        XT(1)=XLOW
        IF(XLOW.LT.1.0D-16) THEN
          XT(2)=XLOW+1.0D-6*(XHIGH-XLOW)
          I1=2
        ELSE
          I1=1
        ENDIF
        FX=EXP(LOG(XHIGH/XT(I1))/DBLE(NP-I1))
        DO I=I1,NP-1
          XT(I+1)=XT(I)*FX
        ENDDO
        XT(NP)=XHIGH
      ENDIF
C
      DO I=1,NP-1
        DX=XT(I+1)-XT(I)
        DXI=DX/DBLE(NIP-1)
        PDFMAX=0.0D0
        DO K=1,NIP
          XS(K)=XT(I)+DBLE(K-1)*DXI
          YS(K)=MAX(PDF(XS(K)),ZEROT)
          PDFMAX=MAX(PDFMAX,YS(K))
        ENDDO
C  ****  Simpson's integration.
        CONS=DXI*3.3333333333333333D-1*0.5D0
        SUMI(1)=0.0D0
        DO K=2,NIP
          XIH=XS(K)-0.5D0*DXI
          YSH=MAX(PDF(XIH),ZEROT)
          PDFMAX=MAX(PDFMAX,YSH)
          SUMI(K)=SUMI(K-1)+CONS*(YS(K-1)+4.0D0*YSH+YS(K))
        ENDDO
C
        DPAC(I)=SUMI(NIP)
        FACT=1.0D0/DPAC(I)
        DO K=1,NIP
          SUMI(K)=FACT*SUMI(K)
        ENDDO
C  ****  When the PDF vanishes at one of the interval end points, its
C        value is modified.
        IF(YS(1).LT.ZERO) YS(1)=1.0D-5*PDFMAX
        IF(YS(NIP).LT.ZERO) YS(NIP)=1.0D-5*PDFMAX
C
        PLI=YS(1)*FACT
        PUI=YS(NIP)*FACT
        B(I)=1.0D0-1.0D0/(PLI*PUI*DX*DX)
        A(I)=(1.0D0/(PLI*DX))-1.0D0-B(I)
        C(I)=1.0D0+A(I)+B(I)
        IF(C(I).LT.ZERO) THEN
          A(I)=0.0D0
          B(I)=0.0D0
          C(I)=1.0D0
        ENDIF
C
C  ****  ERR(I) is the integral of the absolute difference between the
C  rational interpolation and the true PDF, extended over the interval
C  (XT(I),XT(I+1)). Calculated using the trapezoidal rule.
C
        ICASE=1
  100   CONTINUE
        ERR(I)=0.0D0
        DO K=1,NIP
          RR=SUMI(K)
          PAP=DPAC(I)*(1.0D0+(A(I)+B(I)*RR)*RR)**2/
     1       ((1.0D0-B(I)*RR*RR)*C(I)*(XT(I+1)-XT(I)))
          IF(K.EQ.1.OR.K.EQ.NIP) THEN
            ERR(I)=ERR(I)+0.5D0*ABS(PAP-YS(K))
          ELSE
            ERR(I)=ERR(I)+ABS(PAP-YS(K))
          ENDIF
        ENDDO
        ERR(I)=ERR(I)*DXI
C  ****  If ERR(I) is too large, the PDF is approximated by a uniform
C        distribution.
        IF(ERR(I).GT.0.10D0*DPAC(I).AND.ICASE.EQ.1) THEN
          B(I)=0.0D0
          A(I)=0.0D0
          C(I)=1.0D0
          ICASE=2
          GO TO 100
        ENDIF
      ENDDO
      XT(NP)=XHIGH
      A(NP)=0.0D0
      B(NP)=0.0D0
      C(NP)=0.0D0
      ERR(NP)=0.0D0
      DPAC(NP)=0.0D0
C
C  ****  New grid points are added by halving the subinterval with the
C        largest absolute error.
C
  200 CONTINUE
      ERRM=0.0D0
      LMAX=1
      DO I=1,NP-1
C  ****  ERRM is the largest of the interval errors ERR(I).
        IF(ERR(I).GT.ERRM) THEN
          ERRM=ERR(I)
          LMAX=I
        ENDIF
      ENDDO
C
      NP=NP+1
      DO I=NP,LMAX+1,-1
        XT(I)=XT(I-1)
        A(I)=A(I-1)
        B(I)=B(I-1)
        C(I)=C(I-1)
        ERR(I)=ERR(I-1)
        DPAC(I)=DPAC(I-1)
      ENDDO
      XT(LMAX+1)=0.5D0*(XT(LMAX)+XT(LMAX+2))
      DO I=LMAX,LMAX+1
        DX=XT(I+1)-XT(I)
        DXI=(XT(I+1)-XT(I))/DBLE(NIP-1)
        PDFMAX=0.0D0
        DO K=1,NIP
          XS(K)=XT(I)+DBLE(K-1)*DXI
          YS(K)=MAX(PDF(XS(K)),ZEROT)
          PDFMAX=MAX(PDFMAX,YS(K))
        ENDDO
C  ****  Simpson's integration.
        CONS=DXI*3.3333333333333333D-1*0.5D0
        SUMI(1)=0.0D0
        DO K=2,NIP
          XIH=XS(K)-0.5D0*DXI
          YSH=MAX(PDF(XIH),ZEROT)
          SUMI(K)=SUMI(K-1)+CONS*(YS(K-1)+4.0D0*YSH+YS(K))
        ENDDO
C
        DPAC(I)=SUMI(NIP)
        FACT=1.0D0/DPAC(I)
        DO K=1,NIP
          SUMI(K)=FACT*SUMI(K)
        ENDDO
C
        IF(YS(1).LT.ZERO) YS(1)=1.0D-5*PDFMAX
        IF(YS(NIP).LT.ZERO) YS(NIP)=1.0D-5*PDFMAX
        PLI=YS(1)*FACT
        PUI=YS(NIP)*FACT
        B(I)=1.0D0-1.0D0/(PLI*PUI*DX*DX)
        A(I)=(1.0D0/(PLI*DX))-1.0D0-B(I)
        C(I)=1.0D0+A(I)+B(I)
        IF(C(I).LT.ZERO) THEN
          A(I)=0.0D0
          B(I)=0.0D0
          C(I)=1.0D0
        ENDIF
C
        ICASE=1
  300   CONTINUE
        ERR(I)=0.0D0
        DO K=1,NIP
          RR=SUMI(K)
          PAP=DPAC(I)*(1.0D0+(A(I)+B(I)*RR)*RR)**2/
     1       ((1.0D0-B(I)*RR*RR)*C(I)*(XT(I+1)-XT(I)))
          IF(K.EQ.1.OR.K.EQ.NIP) THEN
            ERR(I)=ERR(I)+0.5D0*ABS(PAP-YS(K))
          ELSE
            ERR(I)=ERR(I)+ABS(PAP-YS(K))
          ENDIF
        ENDDO
        ERR(I)=ERR(I)*DXI
C
        IF(ERR(I).GT.0.10D0*DPAC(I).AND.ICASE.EQ.1) THEN
          B(I)=0.0D0
          A(I)=0.0D0
          C(I)=1.0D0
          ICASE=2
          GO TO 300
        ENDIF
      ENDDO
C
      IF(NP.LT.N) GO TO 200
      NPM1=NP-1
C
C  ****  Renormalisation.
C
      CNORM=0.0D0
      DO I=1,NPM1
        CNORM=CNORM+DPAC(I)
      ENDDO
      CNORM=1.0D0/CNORM
      ERRM=0.0D0
      DO I=1,NPM1
        DPAC(I)=DPAC(I)*CNORM
        ERR(I)=ERR(I)*CNORM
        ERRM=MAX(ERRM,ERR(I))
      ENDDO
C
      PAC(1)=0.0D0
      DO I=1,NPM1
        PAC(I+1)=PAC(I)+DPAC(I)
      ENDDO
      PAC(NP)=1.0D0
C
C  ****  Pre-calculated limits for the initial binary search in
C        subroutine RITAI.
C
      BIN=1.0D0/DBLE(NPM1)
      IL(1)=1
      DO I=2,NPM1
        PTST=(I-1)*BIN
        DO J=IL(I-1),NP
          IF(PAC(J).GT.PTST) THEN
            IL(I)=J-1
            IU(I-1)=J
            GO TO 400
          ENDIF
        ENDDO
  400   CONTINUE
      ENDDO
      IU(NPM1)=NP
      IL(NP)=NP-1
      IU(NP)=NP
C
C  ****  Print interpolation tables (only when IWR.GT.0).
C
      IF(IWR.GT.0) THEN
        OPEN(IWR,FILE='param.dat')
        WRITE(IWR,1000) ERRM
 1000   FORMAT(1X,'#  Interpolation error =',1P,E15.7)
        WRITE(IWR,1001) CNORM
 1001   FORMAT(1X,'# Normalising constant =',1P,E15.7)
        WRITE(IWR,1002)
 1002   FORMAT(1X,'#',5X,'X',11X,'PDF(X)',10X,'A',13X,'B',13X,'C',
     1    11X,'error')
        DO I=1,NPM1
          PDFE=MAX(PDF(XT(I)),ZEROT)*CNORM
          WRITE(IWR,'(1P,7E14.6)') XT(I),PDFE,A(I),B(I),C(I),ERR(I)
        ENDDO
        CLOSE(IWR)
C
        OPEN(IWR,FILE='table.dat')
        WRITE(IWR,1000) ERRM
        WRITE(IWR,1001) CNORM
        WRITE(IWR,2000)
 2000   FORMAT(1X,'#',6X,'X',13X,'PDF_ex',10X,'PDF_ap',11X,'err')
        DO I=1,NPM1
          DX=(XT(I+1)-XT(I))/DBLE(NIP-1)
          DO K=1,NIP,5
            XTAU=XT(I)+(K-1)*DX
            P1=MAX(PDF(XTAU)*CNORM,ZEROT)
C  ****  Rational interpolation.
            TAU=(XTAU-XT(I))/(XT(I+1)-XT(I))
            CON1=2.0D0*B(I)*TAU
            CON2=C(I)-A(I)*TAU
            IF(ABS(CON1).GT.1.0D-10*ABS(CON2)) THEN
              ETA=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))/CON1
            ELSE
              ETA=TAU/CON2
            ENDIF
            P2=DPAC(I)*(1.0D0+(A(I)+B(I)*ETA)*ETA)**2
     1        /((1.0D0-B(I)*ETA*ETA)*C(I)*(XT(I+1)-XT(I)))
            WRITE(IWR,'(1P,5E16.8)') XTAU,P1,P2,(P1-P2)/P1
          ENDDO
        ENDDO
        CLOSE(IWR)
C
        OPEN(IWR,FILE='limits.dat')
        WRITE(IWR,3000)
 3000   FORMAT(1X,'#  I',6X,'PAC(ITL)',9X,'(I-1)/NPM1',11X,'I/NPM1',
     1    12X,'PAC(ITU)',12X,'PAC(I)')
        DO I=1,NPM1
          WRITE(IWR,'(I5,1P,5E19.11)') I,PAC(IL(I)),(I-1)*BIN,I*BIN,
     1      PAC(IU(I)),PAC(I)
          IF(PAC(IL(I)).GT.(I-1)*BIN+EPS.OR.
     1      PAC(IU(I)).LT.I*BIN-EPS) THEN
            WRITE(IWR,3001)
 3001       FORMAT(' #  WARNING: The first four values should be in in',
     1        'creasing order.')
          ENDIF
        ENDDO
        CLOSE(IWR)
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RITAV
C  *********************************************************************
      SUBROUTINE RITAV(X,PDF,CDF)
C
C  This subroutine gives the values of the (normalised) pdf, PDF, and of
C  its cumulative distribution function, CDF, at the point X. These
C  values are calculated from the RITA approximation.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NM=512)
      COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A(NM),B(NM),IL(NM),IU(NM),
     1  NPM1
C
      IF(X.GT.XT(NPM1+1)) THEN
        PDF=0.0D0
        CDF=1.0D0
      ELSE IF(X.LT.XT(1)) THEN
        PDF=0.0D0
        CDF=0.0D0
      ELSE
        I=1
        I1=NPM1+1
    1   IT=(I+I1)/2
        IF(X.GT.XT(IT)) THEN
          I=IT
        ELSE
          I1=IT
        ENDIF
        IF(I1-I.GT.1) GO TO 1
        TAU=(X-XT(I))/(XT(I+1)-XT(I))
        CON1=2.0D0*B(I)*TAU
        CON2=1.0D0+B(I)+A(I)*(1.0D0-TAU)
        IF(ABS(CON1).GT.1.0D-10*ABS(CON2)) THEN
          ETA=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))/CON1
        ELSE
          ETA=TAU/CON2
        ENDIF
        PDF=DPAC(I)*(1.0D0+(A(I)+B(I)*ETA)*ETA)**2
     1    /((1.0D0-B(I)*ETA*ETA)*(1.0D0+A(I)+B(I))*(XT(I+1)-XT(I)))
        CDF=PAC(I)+ETA*DPAC(I)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE RITAM
C  *********************************************************************
      SUBROUTINE RITAM(XD,XU,XM0,XM1,XM2)
C
C  Calculation of (restricted) momenta of a pdf, PDF(X), obtained from
C  its RITA approximation.
C
C     XD, XU ... limits of the integration interval.
C     XM0 ...... 0th order moment.
C     XM1 ...... 1st order moment.
C     XM2 ...... 2nd order moment.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      PARAMETER (NM=512)
      COMMON/CRITA/XT(NM),PAC(NM),DPAC(NM),A(NM),B(NM),IL(NM),IU(NM),
     1  NPM1
      PARAMETER (NIP=51)
      DIMENSION XS(NIP),YS(NIP)
C
      XM0=0.0D0
      XM1=0.0D0
      XM2=0.0D0
      DO I=1,NPM1
        IF(XT(I+1).GE.XD.AND.XT(I).LE.XU) THEN
          X1=MAX(XT(I),XD)
          X2=MIN(XT(I+1),XU)
          DX=(X2-X1)/DBLE(NIP-1)
C
          DO K=1,NIP
            XS(K)=X1+DBLE(K-1)*DX
C  ****  Value of the RITA rational pdf at the point XS(K).
            TAU=(XS(K)-XT(I))/(XT(I+1)-XT(I))
            CON1=2.0D0*B(I)*TAU
            CON2=1.0D0+B(I)+A(I)*(1.0D0-TAU)
            IF(ABS(CON1).GT.1.0D-10*ABS(CON2)) THEN
              ETA=CON2*(1.0D0-SQRT(1.0D0-2.0D0*TAU*CON1/CON2**2))/CON1
            ELSE
              ETA=TAU/CON2
            ENDIF
            YS(K)=DPAC(I)*(1.0D0+(A(I)+B(I)*ETA)*ETA)**2
     1        /((1.0D0-B(I)*ETA*ETA)*(1.0D0+A(I)+B(I))*(XT(I+1)-XT(I)))
          ENDDO
C
C  ****  Simpson's integration.
C
          CONS=DX*3.3333333333333333D-1
          SUM=0.0D0
          DO L=3,NIP,2
            SUM=SUM+YS(L-2)+4.0D0*YS(L-1)+YS(L)
          ENDDO
          XM0=XM0+SUM*CONS
C
          DO K=1,NIP
            YS(K)=YS(K)*XS(K)
          ENDDO
          SUM=0.0D0
          DO L=3,NIP,2
            SUM=SUM+YS(L-2)+4.0D0*YS(L-1)+YS(L)
          ENDDO
          XM1=XM1+SUM*CONS
C
          DO K=1,NIP
            YS(K)=YS(K)*XS(K)
          ENDDO
          SUM=0.0D0
          DO L=3,NIP,2
            SUM=SUM+YS(L-2)+4.0D0*YS(L-1)+YS(L)
          ENDDO
          XM2=XM2+SUM*CONS
        ENDIF
      ENDDO
      RETURN
      END


C  *********************************************************************
C                         SUBROUTINE RAND0
C  *********************************************************************
      SUBROUTINE RAND0(N)
C
C  In parallel calculations, we need to initialise RAND() in such a way
C  that different processors produce truly independent sequences of
C  random numbers. This can be accomplished by feeding each processor
C  with initial seeds that belong to a single long sequence generated by
C  RAND(), but are far enough from each other to avoid overlap of the
C  sub-sequences generated by the different processors. The list below
C  was obtained by running a program written by Andreu Badal and Josep
C  Sempau [Comp. Phys. Commun. 175 (2006) 440–450]. It contains pairs of
C  seeds that belong to a long sequence and whose separation is 10**14
C  calls. That is, if we start with the N-th pair of seeds, after 10**14
C  calls to RAND() we obtain the (N+1)-th pair.
C
C  A call to the present subroutine RAND0(N) loads the N-th pair of
C  seeds in the list. Thus, in parallel simulations, we can initialise
C  RAND() in the different processors by calling RAND0 with different
C  values of the argument N.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      COMMON/RSEED/ISEED1,ISEED2
      DIMENSION IS1(1001),IS2(1001)
C
      IS1(1   )=         1; IS2(1   )=         1
      IS1(2   )=1088794366; IS2(2   )= 722792456
      IS1(3   )=1993751964; IS2(3   )= 753089694
      IS1(4   )= 610005387; IS2(4   )=1748134360
      IS1(5   )=  27944595; IS2(5   )= 774572312
      IS1(6   )=1108394934; IS2(6   )= 620441713
      IS1(7   )= 580582953; IS2(7   )=  73753031
      IS1(8   )= 479095664; IS2(8   )= 801971873
      IS1(9   )= 204282520; IS2(9   )=1624172377
      IS1(10  )= 327796803; IS2(10  )= 915478773
      IS1(11  )= 918882992; IS2(11  )= 858672133
      IS1(12  )= 210516912; IS2(12  )= 158496513
      IS1(13  )=1155963390; IS2(13  )= 704005330
      IS1(14  )=1929580647; IS2(14  )= 324681933
      IS1(15  )=2061216471; IS2(15  )= 209323758
      IS1(16  )=1963942692; IS2(16  )=1686825906
      IS1(17  )=2124896501; IS2(17  )= 460463303
      IS1(18  )= 233299547; IS2(18  )=2134948096
      IS1(19  )=1675549170; IS2(19  )=1144329343
      IS1(20  )=1421385813; IS2(20  )= 672532778
      IS1(21  )=2069007070; IS2(21  )=1309916099
      IS1(22  )=1905029285; IS2(22  )=1817420725
      IS1(23  )=1056478052; IS2(23  )=1943244114
      IS1(24  )=1628262065; IS2(24  )= 782316895
      IS1(25  )=1368245324; IS2(25  )=2070604605
      IS1(26  )= 414218397; IS2(26  )= 519652740
      IS1(27  )=1032135553; IS2(27  )=1729805542
      IS1(28  )=1709285619; IS2(28  )=1859376254
      IS1(29  )=1363836736; IS2(29  )= 433738184
      IS1(30  )=1307886537; IS2(30  )=1841770585
      IS1(31  )= 944675654; IS2(31  )=1438406465
      IS1(32  )= 733913518; IS2(32  )= 219173860
      IS1(33  )=1994697060; IS2(33  )=1473267531
      IS1(34  )=1816284687; IS2(34  )=1806014371
      IS1(35  )=1204077982; IS2(35  )= 313513580
      IS1(36  )=1322304106; IS2(36  )=1896720184
      IS1(37  )=1152681299; IS2(37  )=1347076753
      IS1(38  )=1910575054; IS2(38  )=1881544152
      IS1(39  )=1158165315; IS2(39  )=1769280665
      IS1(40  )=1256652686; IS2(40  )=1678458141
      IS1(41  )= 149156960; IS2(41  )= 257442270
      IS1(42  )=1531563057; IS2(42  )= 769379539
      IS1(43  )=  33227674; IS2(43  )=1940657859
      IS1(44  )= 563899375; IS2(44  )=1495802342
      IS1(45  )= 534312669; IS2(45  )=1627605675
      IS1(46  )=  20883767; IS2(46  )=1860146428
      IS1(47  )= 654784964; IS2(47  )=1449089950
      IS1(48  )=1622565819; IS2(48  )= 300004830
      IS1(49  )= 720217438; IS2(49  )=  59015257
      IS1(50  )=1162709863; IS2(50  )=1966621140
      IS1(51  )= 360537627; IS2(51  )= 133123709
      IS1(52  )= 345935498; IS2(52  )=1389524816
      IS1(53  )= 876665055; IS2(53  )=1171092501
      IS1(54  )=1998564166; IS2(54  )=1037986676
      IS1(55  )=2066661466; IS2(55  )= 366166383
      IS1(56  )=1279875030; IS2(56  )= 582766249
      IS1(57  )=1955021266; IS2(57  )= 116351171
      IS1(58  )= 195816083; IS2(58  )= 933489046
      IS1(59  )=2055228393; IS2(59  )= 283595456
      IS1(60  )=  93961337; IS2(60  )= 611643703
      IS1(61  )=1446789139; IS2(61  )=1248992867
      IS1(62  )= 414524648; IS2(62  )= 803732805
      IS1(63  )=1327723283; IS2(63  )= 201031458
      IS1(64  )= 535511919; IS2(64  )= 938162313
      IS1(65  )=  86014814; IS2(65  )=   7981320
      IS1(66  )=1646051327; IS2(66  )=1690586644
      IS1(67  )=1717255367; IS2(67  )= 569003207
      IS1(68  )=1794770573; IS2(68  )=1043545502
      IS1(69  )= 853934312; IS2(69  )= 848156009
      IS1(70  )= 903676328; IS2(70  )=1753511888
      IS1(71  )= 888673974; IS2(71  )=2014364429
      IS1(72  )=1763908015; IS2(72  )= 835386162
      IS1(73  )= 415969825; IS2(73  )= 751041626
      IS1(74  )=2041881831; IS2(74  )= 527275421
      IS1(75  )=1223957624; IS2(75  )= 282891294
      IS1(76  )= 216252620; IS2(76  )=1933221826
      IS1(77  )= 956133864; IS2(77  )=1506952100
      IS1(78  )=1265794918; IS2(78  )= 713084270
      IS1(79  )=1474654961; IS2(79  )=  88797131
      IS1(80  )= 122865758; IS2(80  )=2073278963
      IS1(81  )=    258943; IS2(81  )= 664687714
      IS1(82  )=1152463120; IS2(82  )= 839536051
      IS1(83  )= 181367474; IS2(83  )=1038508765
      IS1(84  )= 618933039; IS2(84  )= 131404490
      IS1(85  )=1185139338; IS2(85  )=1987460903
      IS1(86  )=2078683375; IS2(86  )=1391654211
      IS1(87  )=1157287301; IS2(87  )=1870939725
      IS1(88  )= 490572205; IS2(88  )=1828672152
      IS1(89  )= 713452544; IS2(89  )=1880501392
      IS1(90  )=1399731654; IS2(90  )= 661442337
      IS1(91  )=1434784182; IS2(91  )=1598489021
      IS1(92  )= 157980824; IS2(92  )= 512999623
      IS1(93  )=1631668015; IS2(93  )=2021483107
      IS1(94  )= 695840037; IS2(94  )= 316143310
      IS1(95  )=2011902133; IS2(95  )= 460681770
      IS1(96  )=1482456963; IS2(96  )=1832621179
      IS1(97  )= 983630146; IS2(97  )=1244860854
      IS1(98  )= 424488068; IS2(98  )=1894080738
      IS1(99  )= 592109479; IS2(99  )= 254844002
      IS1(100 )= 698713089; IS2(100 )= 852473215
      IS1(101 )= 698429770; IS2(101 )=1978724894
      IS1(102 )=1991339714; IS2(102 )= 590904772
      IS1(103 )=1812612029; IS2(103 )=1455273894
      IS1(104 )= 878555690; IS2(104 )= 849354664
      IS1(105 )=1415741666; IS2(105 )=1989849407
      IS1(106 )= 740336773; IS2(106 )=2094756349
      IS1(107 )=1357150877; IS2(107 )=1752614313
      IS1(108 )= 446288602; IS2(108 )= 605474927
      IS1(109 )= 156511135; IS2(109 )=1112863746
      IS1(110 )=1315731039; IS2(110 )=1420661365
      IS1(111 )=1590179518; IS2(111 )= 797341276
      IS1(112 )= 210183789; IS2(112 )=1297307595
      IS1(113 )=  94234820; IS2(113 )=1379679575
      IS1(114 )= 273023900; IS2(114 )=1595938728
      IS1(115 )= 868831745; IS2(115 )=1149726246
      IS1(116 )= 170384549; IS2(116 )= 129431617
      IS1(117 )=  13185587; IS2(117 )= 487985593
      IS1(118 )=1362898234; IS2(118 )= 342281893
      IS1(119 )= 574105532; IS2(119 )= 245617051
      IS1(120 )=1822103760; IS2(120 )=1034640984
      IS1(121 )= 658812725; IS2(121 )=1829379549
      IS1(122 )=1005671726; IS2(122 )= 183096918
      IS1(123 )=1254435204; IS2(123 )=1402468728
      IS1(124 )=1798478003; IS2(124 )= 432519190
      IS1(125 )= 601935466; IS2(125 )= 253536637
      IS1(126 )= 341666647; IS2(126 )= 118329947
      IS1(127 )=1713183513; IS2(127 )=2027318311
      IS1(128 )=1196735493; IS2(128 )=1530795526
      IS1(129 )=1348986425; IS2(129 )= 668877863
      IS1(130 )= 532005772; IS2(130 )= 203741901
      IS1(131 )=1141813962; IS2(131 )=1863091192
      IS1(132 )=1741278758; IS2(132 )= 104388873
      IS1(133 )=  86859261; IS2(133 )=2120737529
      IS1(134 )= 726923420; IS2(134 )=1212110048
      IS1(135 )=1058541527; IS2(135 )= 166531451
      IS1(136 )=2131549752; IS2(136 )=1447596029
      IS1(137 )=2033958033; IS2(137 )= 926546635
      IS1(138 )= 969574276; IS2(138 )=1513067587
      IS1(139 )=1924153065; IS2(139 )=1224270071
      IS1(140 )=1789201564; IS2(140 )= 244925517
      IS1(141 )= 970004038; IS2(141 )= 827424326
      IS1(142 )= 584997635; IS2(142 )= 900020443
      IS1(143 )=1121567821; IS2(143 )=1551288142
      IS1(144 )=1901695144; IS2(144 )= 802507892
      IS1(145 )=1381949729; IS2(145 )= 338664653
      IS1(146 )= 931183121; IS2(146 )= 255723333
      IS1(147 )= 425040923; IS2(147 )=1183865313
      IS1(148 )=2063648383; IS2(148 )=1931656776
      IS1(149 )= 372520795; IS2(149 )=1381463141
      IS1(150 )= 112959009; IS2(150 )= 328713331
      IS1(151 )= 156172654; IS2(151 )=2028805919
      IS1(152 )=1206630112; IS2(152 )=1317701868
      IS1(153 )=1141325584; IS2(153 )=1635316096
      IS1(154 )=1226401966; IS2(154 )= 891065751
      IS1(155 )=1044869640; IS2(155 )=1695983251
      IS1(156 )=1468275435; IS2(156 )= 827674970
      IS1(157 )= 791167482; IS2(157 )= 645339068
      IS1(158 )= 464714547; IS2(158 )=1616229175
      IS1(159 )=1084778504; IS2(159 )= 563204166
      IS1(160 )= 258987949; IS2(160 )=1152030385
      IS1(161 )= 479486796; IS2(161 )= 238097920
      IS1(162 )=1499316991; IS2(162 )=1401468685
      IS1(163 )= 519780735; IS2(163 )= 481196391
      IS1(164 )=1679611087; IS2(164 )= 604709095
      IS1(165 )=1691996345; IS2(165 )= 989109993
      IS1(166 )= 196557364; IS2(166 )= 113178392
      IS1(167 )= 851784008; IS2(167 )=1222308139
      IS1(168 )= 143277176; IS2(168 )=1931319584
      IS1(169 )=1973626791; IS2(169 )=1586075498
      IS1(170 )= 585402145; IS2(170 )= 125579035
      IS1(171 )=1926622811; IS2(171 )=1383430112
      IS1(172 )= 614117445; IS2(172 )= 809143743
      IS1(173 )= 209722147; IS2(173 )=1485956780
      IS1(174 )= 445830939; IS2(174 )=1964280618
      IS1(175 )=1346542997; IS2(175 )= 748234912
      IS1(176 )= 576549607; IS2(176 )= 554392763
      IS1(177 )=1852906063; IS2(177 )= 675586306
      IS1(178 )=1415103532; IS2(178 )= 839112213
      IS1(179 )= 868356749; IS2(179 )=1226343583
      IS1(180 )=1373222177; IS2(180 )=1149205884
      IS1(181 )=1024191502; IS2(181 )= 938910203
      IS1(182 )= 463832557; IS2(182 )= 441736082
      IS1(183 )= 599161815; IS2(183 )= 880842978
      IS1(184 )=  27305702; IS2(184 )=1259875018
      IS1(185 )=1622662871; IS2(185 )=1545669937
      IS1(186 )=1314825492; IS2(186 )=1230485856
      IS1(187 )=1919359439; IS2(187 )= 649279764
      IS1(188 )= 776491967; IS2(188 )= 940088497
      IS1(189 )= 153029369; IS2(189 )= 604610332
      IS1(190 )=1075161827; IS2(190 )= 333444224
      IS1(191 )= 914109165; IS2(191 )= 803254235
      IS1(192 )=1944937918; IS2(192 )=1451340862
      IS1(193 )=1738752454; IS2(193 )= 499708706
      IS1(194 )= 321360177; IS2(194 )=1192481545
      IS1(195 )=1078034966; IS2(195 )= 318667790
      IS1(196 )=1983953435; IS2(196 )=1865542330
      IS1(197 )=1964072934; IS2(197 )=1092645796
      IS1(198 )=1951113931; IS2(198 )=1549447078
      IS1(199 )=1059124001; IS2(199 )=1047835649
      IS1(200 )=1695266076; IS2(200 )= 163616804
      IS1(201 )= 996688518; IS2(201 )=2093478795
      IS1(202 )=1354445549; IS2(202 )= 335951295
      IS1(203 )= 998773555; IS2(203 )=1219483032
      IS1(204 )=   4058649; IS2(204 )= 559033528
      IS1(205 )= 470886335; IS2(205 )=1667952318
      IS1(206 )= 158745647; IS2(206 )=2012602568
      IS1(207 )= 461380034; IS2(207 )= 799155945
      IS1(208 )= 105655873; IS2(208 )=1187620434
      IS1(209 )= 915491195; IS2(209 )= 858595438
      IS1(210 )= 247536109; IS2(210 )= 727545379
      IS1(211 )=1343291889; IS2(211 )=1203423504
      IS1(212 )=2001578188; IS2(212 )= 990991500
      IS1(213 )=1002265824; IS2(213 )= 219121455
      IS1(214 )= 501880984; IS2(214 )= 846802413
      IS1(215 )=1411611767; IS2(215 )= 443023724
      IS1(216 )=1514586313; IS2(216 )= 974179120
      IS1(217 )= 261179117; IS2(217 )=  18434323
      IS1(218 )=  20916375; IS2(218 )=2145251247
      IS1(219 )=1863207976; IS2(219 )= 176589398
      IS1(220 )= 969221669; IS2(220 )= 683394530
      IS1(221 )=1783116228; IS2(221 )= 282077422
      IS1(222 )=1847032911; IS2(222 )=1518960264
      IS1(223 )= 754180209; IS2(223 )= 423392320
      IS1(224 )=1927259077; IS2(224 )= 249198027
      IS1(225 )= 285155942; IS2(225 )=1861278512
      IS1(226 )=1748963900; IS2(226 )=1181877087
      IS1(227 )= 691637076; IS2(227 )=1593810530
      IS1(228 )= 378323627; IS2(228 )= 196639057
      IS1(229 )= 985536851; IS2(229 )= 991565678
      IS1(230 )=1218539427; IS2(230 )=1989132277
      IS1(231 )= 552569869; IS2(231 )=1861318300
      IS1(232 )= 300092931; IS2(232 )= 550437007
      IS1(233 )=1582757652; IS2(233 )=1690689155
      IS1(234 )= 482469061; IS2(234 )= 126744526
      IS1(235 )= 225088150; IS2(235 )=1140456485
      IS1(236 )= 566037433; IS2(236 )= 654854217
      IS1(237 )=1644348452; IS2(237 )= 737109295
      IS1(238 )= 337491116; IS2(238 )=2098977589
      IS1(239 )= 655153746; IS2(239 )=  43299124
      IS1(240 )=1499772543; IS2(240 )=1730234211
      IS1(241 )= 839738220; IS2(241 )=1673889598
      IS1(242 )= 528696432; IS2(242 )=1902853997
      IS1(243 )=  50552080; IS2(243 )=2075206178
      IS1(244 )= 130137340; IS2(244 )=1283599409
      IS1(245 )=1013040075; IS2(245 )= 676361106
      IS1(246 )=1793063152; IS2(246 )=1860713192
      IS1(247 )=1912078503; IS2(247 )= 259429094
      IS1(248 )= 695750580; IS2(248 )=1364366801
      IS1(249 )= 851302736; IS2(249 )=1708079864
      IS1(250 )=1365371254; IS2(250 )= 227135534
      IS1(251 )=1436544680; IS2(251 )= 196014388
      IS1(252 )= 255720485; IS2(252 )=1188024965
      IS1(253 )= 517569477; IS2(253 )=  63939330
      IS1(254 )= 379457723; IS2(254 )= 417116556
      IS1(255 )=1714565676; IS2(255 )= 870984368
      IS1(256 )= 427585641; IS2(256 )= 200189082
      IS1(257 )=1982058823; IS2(257 )=1003464933
      IS1(258 )=1738564860; IS2(258 )=1247323567
      IS1(259 )= 487708829; IS2(259 )= 462209958
      IS1(260 )=1646850245; IS2(260 )=  61645060
      IS1(261 )=1590006138; IS2(261 )=1754830438
      IS1(262 )=1733095787; IS2(262 )=1907130822
      IS1(263 )=1660369047; IS2(263 )=1665129257
      IS1(264 )=1104503590; IS2(264 )=2005102976
      IS1(265 )= 557868773; IS2(265 )=1956949606
      IS1(266 )= 665813373; IS2(266 )=2050201793
      IS1(267 )=1832810072; IS2(267 )= 206523161
      IS1(268 )= 420890754; IS2(268 )=1367078059
      IS1(269 )= 517192491; IS2(269 )=1150925658
      IS1(270 )=1564894415; IS2(270 )=2058874982
      IS1(271 )=  88748046; IS2(271 )=1574175136
      IS1(272 )= 214105914; IS2(272 )=1747446780
      IS1(273 )= 667844668; IS2(273 )= 188322609
      IS1(274 )=1127730224; IS2(274 )= 995741669
      IS1(275 )= 367330131; IS2(275 )= 816142314
      IS1(276 )=1246444093; IS2(276 )= 549605711
      IS1(277 )= 847974161; IS2(277 )= 183325985
      IS1(278 )=1721669301; IS2(278 )= 479407779
      IS1(279 )=1915327339; IS2(279 )=1088013018
      IS1(280 )=1578470586; IS2(280 )=1065050626
      IS1(281 )= 940315134; IS2(281 )=1848214072
      IS1(282 )=1473252673; IS2(282 )=   2877464
      IS1(283 )=1827676712; IS2(283 )=1664447670
      IS1(284 )= 839285700; IS2(284 )=1640026298
      IS1(285 )= 751020328; IS2(285 )=1588681006
      IS1(286 )=1091191738; IS2(286 )=1790306835
      IS1(287 )= 177083683; IS2(287 )=1254397576
      IS1(288 )=2055613182; IS2(288 )= 939654007
      IS1(289 )=1473470878; IS2(289 )= 335189253
      IS1(290 )=1800767926; IS2(290 )= 290320395
      IS1(291 )=1676463528; IS2(291 )=1603536943
      IS1(292 )=1650288797; IS2(292 )= 388374267
      IS1(293 )=2035708356; IS2(293 )= 518387205
      IS1(294 )=1452341404; IS2(294 )= 985322233
      IS1(295 )=1661040488; IS2(295 )=  68406361
      IS1(296 )= 895503595; IS2(296 )=1876688389
      IS1(297 )=1951727335; IS2(297 )= 185661402
      IS1(298 )= 195345739; IS2(298 )=1532771577
      IS1(299 )= 268247973; IS2(299 )=1395541411
      IS1(300 )=1057595183; IS2(300 )= 128171866
      IS1(301 )=1466997063; IS2(301 )=1372373334
      IS1(302 )=1588823091; IS2(302 )=1473252323
      IS1(303 )=1825606993; IS2(303 )= 398379605
      IS1(304 )= 836937167; IS2(304 )=2001897494
      IS1(305 )=  60259114; IS2(305 )=2113852924
      IS1(306 )=1298040193; IS2(306 )= 597799372
      IS1(307 )= 892667157; IS2(307 )=  98544655
      IS1(308 )= 349769327; IS2(308 )=1119519495
      IS1(309 )=1659599388; IS2(309 )= 848436478
      IS1(310 )= 347450508; IS2(310 )= 197988152
      IS1(311 )=1564752903; IS2(311 )= 303075472
      IS1(312 )= 271104778; IS2(312 )=2101408514
      IS1(313 )=1071650412; IS2(313 )= 557206316
      IS1(314 )=2129784998; IS2(314 )= 319291050
      IS1(315 )= 149442067; IS2(315 )=1161643665
      IS1(316 )=1382871443; IS2(316 )= 623775604
      IS1(317 )= 217752411; IS2(317 )= 888310836
      IS1(318 )=1265937666; IS2(318 )= 276467372
      IS1(319 )= 569940604; IS2(319 )=1477486653
      IS1(320 )= 379102003; IS2(320 )= 377940861
      IS1(321 )= 885729895; IS2(321 )=1974466429
      IS1(322 )= 163050126; IS2(322 )=  23068033
      IS1(323 )=1194934955; IS2(323 )=1962429405
      IS1(324 )=1988644587; IS2(324 )=1176217709
      IS1(325 )= 223953149; IS2(325 )= 165457549
      IS1(326 )= 173062795; IS2(326 )=1058081267
      IS1(327 )=1975967738; IS2(327 )= 793527446
      IS1(328 )= 903886181; IS2(328 )=  72550470
      IS1(329 )=1844109661; IS2(329 )=1278970903
      IS1(330 )= 224746454; IS2(330 )=2027692323
      IS1(331 )= 381257506; IS2(331 )= 782649282
      IS1(332 )=1521365813; IS2(332 )=1328897351
      IS1(333 )= 838883107; IS2(333 )= 636334681
      IS1(334 )=1958225287; IS2(334 )=1397155038
      IS1(335 )=1880542285; IS2(335 )=2136705686
      IS1(336 )= 351096309; IS2(336 )=  67624347
      IS1(337 )= 276181341; IS2(337 )= 720002598
      IS1(338 )=1483813475; IS2(338 )=1496522845
      IS1(339 )=1721418406; IS2(339 )= 298856548
      IS1(340 )=1646984747; IS2(340 )= 613517180
      IS1(341 )=1115726648; IS2(341 )=1979141747
      IS1(342 )= 497311653; IS2(342 )= 431235843
      IS1(343 )=1847923343; IS2(343 )=1440287460
      IS1(344 )=1612185030; IS2(344 )=1770007478
      IS1(345 )=1436854618; IS2(345 )=2062850297
      IS1(346 )=1289356410; IS2(346 )= 773503574
      IS1(347 )=1618534459; IS2(347 )= 405022273
      IS1(348 )=1923690772; IS2(348 )=1950961468
      IS1(349 )=2005241207; IS2(349 )= 854563799
      IS1(350 )=1191066623; IS2(350 )= 808720040
      IS1(351 )= 463667715; IS2(351 )= 466536804
      IS1(352 )=1792026494; IS2(352 )=1569987550
      IS1(353 )=1145504660; IS2(353 )=1622132321
      IS1(354 )= 885012929; IS2(354 )=1092533602
      IS1(355 )=1224139137; IS2(355 )=1840751652
      IS1(356 )=1930672614; IS2(356 )=1637040668
      IS1(357 )=1073822199; IS2(357 )=1783545033
      IS1(358 )= 921097169; IS2(358 )=1938213600
      IS1(359 )=1425751390; IS2(359 )=1172781558
      IS1(360 )=1214921245; IS2(360 )= 825459365
      IS1(361 )=2040403296; IS2(361 )=1533648867
      IS1(362 )= 896839067; IS2(362 )=1828625926
      IS1(363 )=1663221476; IS2(363 )= 623151978
      IS1(364 )=1820439500; IS2(364 )=1190628682
      IS1(365 )=1787778713; IS2(365 )=2139868031
      IS1(366 )=1312411209; IS2(366 )= 124016638
      IS1(367 )=1511089893; IS2(367 )=1123651614
      IS1(368 )= 630843631; IS2(368 )=1212169916
      IS1(369 )=1564923744; IS2(369 )= 514797409
      IS1(370 )= 258126650; IS2(370 )=1007434416
      IS1(371 )=1135521143; IS2(371 )= 630960850
      IS1(372 )= 830276638; IS2(372 )= 851724397
      IS1(373 )= 376035913; IS2(373 )= 808391452
      IS1(374 )=2023656286; IS2(374 )= 465517081
      IS1(375 )=1766701603; IS2(375 )=1993165647
      IS1(376 )=1259750908; IS2(376 )= 951063358
      IS1(377 )=1644557611; IS2(377 )=1363160828
      IS1(378 )=1583850975; IS2(378 )=1328161074
      IS1(379 )= 579027304; IS2(379 )=1626248155
      IS1(380 )=1175583557; IS2(380 )=1137630999
      IS1(381 )= 253025339; IS2(381 )= 222102409
      IS1(382 )=1864136836; IS2(382 )=1013284156
      IS1(383 )= 447381646; IS2(383 )= 720482175
      IS1(384 )=1422107210; IS2(384 )= 101344372
      IS1(385 )=  60187744; IS2(385 )=1036992166
      IS1(386 )= 736865928; IS2(386 )=1011413694
      IS1(387 )=1223364420; IS2(387 )=1661421549
      IS1(388 )= 199571836; IS2(388 )=1149316001
      IS1(389 )= 118600702; IS2(389 )= 498570418
      IS1(390 )=1053021159; IS2(390 )=1200634496
      IS1(391 )=1054745378; IS2(391 )=1927760137
      IS1(392 )=1594338747; IS2(392 )=1093514040
      IS1(393 )= 694107219; IS2(393 )= 541441173
      IS1(394 )=1957890210; IS2(394 )= 910785465
      IS1(395 )=1470833484; IS2(395 )= 426201828
      IS1(396 )= 585868751; IS2(396 )=1069807689
      IS1(397 )=1134131445; IS2(397 )=1067457516
      IS1(398 )= 921089340; IS2(398 )=1038655815
      IS1(399 )= 616921523; IS2(399 )=1366192583
      IS1(400 )=  41330973; IS2(400 )= 255560572
      IS1(401 )= 164189022; IS2(401 )=  49014916
      IS1(402 )=1161928238; IS2(402 )= 209936365
      IS1(403 )=2020361273; IS2(403 )=1950362287
      IS1(404 )= 511719771; IS2(404 )= 325061593
      IS1(405 )= 454126255; IS2(405 )=1574510902
      IS1(406 )=1859816564; IS2(406 )=1632823687
      IS1(407 )=1432776991; IS2(407 )= 227252761
      IS1(408 )= 765835313; IS2(408 )=2029746355
      IS1(409 )= 547912517; IS2(409 )= 591050613
      IS1(410 )=1863964385; IS2(410 )= 712242677
      IS1(411 )=1995501485; IS2(411 )=1312458862
      IS1(412 )=1905501124; IS2(412 )= 276433488
      IS1(413 )= 651011325; IS2(413 )= 278589745
      IS1(414 )= 616628955; IS2(414 )=1231979682
      IS1(415 )= 625576690; IS2(415 )=  76923407
      IS1(416 )= 214901782; IS2(416 )= 956950803
      IS1(417 )=1874851100; IS2(417 )= 630418924
      IS1(418 )=1796015654; IS2(418 )=1799191741
      IS1(419 )= 378100074; IS2(419 )= 405443763
      IS1(420 )=2115599125; IS2(420 )=1158325172
      IS1(421 )=  51748422; IS2(421 )=1108481669
      IS1(422 )=1519507484; IS2(422 )=1078568949
      IS1(423 )= 231731663; IS2(423 )=1639105310
      IS1(424 )=1813510342; IS2(424 )=1756686695
      IS1(425 )=1180641209; IS2(425 )= 817484587
      IS1(426 )= 547335020; IS2(426 )=1740348176
      IS1(427 )= 981294631; IS2(427 )=1645483561
      IS1(428 )= 267267642; IS2(428 )= 525489921
      IS1(429 )= 975682868; IS2(429 )=1681199536
      IS1(430 )=1070208555; IS2(430 )= 572418479
      IS1(431 )=1998085041; IS2(431 )=1713183034
      IS1(432 )= 173175676; IS2(432 )= 486314861
      IS1(433 )=1393518568; IS2(433 )=  48960372
      IS1(434 )=1341932065; IS2(434 )=1865938542
      IS1(435 )= 801752013; IS2(435 )= 341065424
      IS1(436 )=2132704165; IS2(436 )= 494928752
      IS1(437 )= 422454854; IS2(437 )= 743233701
      IS1(438 )=1402307772; IS2(438 )= 303431257
      IS1(439 )=1557964062; IS2(439 )=1825819623
      IS1(440 )=1894865913; IS2(440 )=1365259674
      IS1(441 )= 622187675; IS2(441 )=1865578472
      IS1(442 )=1326023361; IS2(442 )= 122041713
      IS1(443 )=1730100218; IS2(443 )= 365085301
      IS1(444 )= 243501496; IS2(444 )=1010792790
      IS1(445 )=1899088012; IS2(445 )= 162527744
      IS1(446 )= 421847337; IS2(446 )=2076217683
      IS1(447 )=1760083121; IS2(447 )= 891099538
      IS1(448 )=2057638875; IS2(448 )=1503480695
      IS1(449 )= 195671618; IS2(449 )=2085854796
      IS1(450 )=1810716138; IS2(450 )= 477554292
      IS1(451 )=1801885889; IS2(451 )= 710751106
      IS1(452 )=1154047452; IS2(452 )=1841903658
      IS1(453 )= 667706413; IS2(453 )=1964945942
      IS1(454 )=2017974505; IS2(454 )=1663725788
      IS1(455 )= 884726865; IS2(455 )= 522792338
      IS1(456 )= 926150544; IS2(456 )=1559103945
      IS1(457 )= 690922273; IS2(457 )=  59346276
      IS1(458 )= 712157685; IS2(458 )= 287197618
      IS1(459 )=2059778138; IS2(459 )= 815135778
      IS1(460 )=1190461438; IS2(460 )=1905576318
      IS1(461 )=1564333110; IS2(461 )= 542814824
      IS1(462 )= 271883897; IS2(462 )= 163552060
      IS1(463 )= 344754143; IS2(463 )=1429455140
      IS1(464 )= 668346479; IS2(464 )= 907377254
      IS1(465 )= 906777901; IS2(465 )= 973540088
      IS1(466 )= 211937991; IS2(466 )=1829082247
      IS1(467 )=2057662804; IS2(467 )= 466664141
      IS1(468 )= 685469316; IS2(468 )= 801959481
      IS1(469 )=1806627086; IS2(469 )=1933314854
      IS1(470 )=  27530138; IS2(470 )=1590842779
      IS1(471 )=1972873114; IS2(471 )=1790212125
      IS1(472 )=2028079049; IS2(472 )= 909199739
      IS1(473 )= 568804232; IS2(473 )= 180866254
      IS1(474 )=1897611427; IS2(474 )=1205049755
      IS1(475 )=1529982236; IS2(475 )= 727044114
      IS1(476 )=1930649184; IS2(476 )=1157145550
      IS1(477 )= 579071696; IS2(477 )=1002305204
      IS1(478 )=1522526588; IS2(478 )= 102634606
      IS1(479 )= 502645745; IS2(479 )=2002850332
      IS1(480 )=1455547110; IS2(480 )=1863533555
      IS1(481 )=1865827035; IS2(481 )=  87808690
      IS1(482 )=2146621882; IS2(482 )=2037357581
      IS1(483 )=1348397757; IS2(483 )= 252539429
      IS1(484 )=  74374264; IS2(484 )= 584457262
      IS1(485 )= 909905711; IS2(485 )= 641791254
      IS1(486 )= 406627724; IS2(486 )=1970801280
      IS1(487 )=1869585944; IS2(487 )=2086692085
      IS1(488 )= 471312487; IS2(488 )=  58607088
      IS1(489 )= 314332810; IS2(489 )=1762002696
      IS1(490 )= 714059427; IS2(490 )= 892062884
      IS1(491 )=  88551984; IS2(491 )= 707506711
      IS1(492 )=1764182800; IS2(492 )=1566043351
      IS1(493 )=1660801101; IS2(493 )=  91743653
      IS1(494 )=2053618389; IS2(494 )=  87543326
      IS1(495 )=1872198017; IS2(495 )=1419845282
      IS1(496 )=1718051970; IS2(496 )=  16608354
      IS1(497 )= 106783453; IS2(497 )=1579552005
      IS1(498 )= 198639753; IS2(498 )=1260841572
      IS1(499 )= 444341049; IS2(499 )= 185823881
      IS1(500 )= 822270701; IS2(500 )= 703588888
      IS1(501 )=1741688389; IS2(501 )=1199341216
      IS1(502 )=1740532989; IS2(502 )= 506564970
      IS1(503 )= 306294863; IS2(503 )=1989939786
      IS1(504 )=1282800170; IS2(504 )= 909139593
      IS1(505 )= 869069844; IS2(505 )= 759737034
      IS1(506 )=1099376549; IS2(506 )= 592166139
      IS1(507 )= 530526621; IS2(507 )= 372525993
      IS1(508 )= 132237605; IS2(508 )=1235725904
      IS1(509 )=1250152263; IS2(509 )= 734059529
      IS1(510 )= 997470030; IS2(510 )= 853534414
      IS1(511 )=1214905199; IS2(511 )=1227201813
      IS1(512 )=1024791465; IS2(512 )=1264083624
      IS1(513 )= 215462734; IS2(513 )=1323361591
      IS1(514 )=1746861828; IS2(514 )=1171508714
      IS1(515 )=  98049234; IS2(515 )=     76692
      IS1(516 )=1941536515; IS2(516 )=1557540564
      IS1(517 )=1062043665; IS2(517 )=1536279542
      IS1(518 )=1933220889; IS2(518 )= 531737550
      IS1(519 )=1920928412; IS2(519 )=1961452165
      IS1(520 )= 924668593; IS2(520 )=1126181753
      IS1(521 )= 483780576; IS2(521 )=1943663885
      IS1(522 )=1172795790; IS2(522 )= 902336756
      IS1(523 )=  57254329; IS2(523 )= 548344687
      IS1(524 )= 742985485; IS2(524 )=  75812010
      IS1(525 )=1138516383; IS2(525 )= 704793701
      IS1(526 )= 516069233; IS2(526 )= 658536656
      IS1(527 )= 767025613; IS2(527 )=1696634101
      IS1(528 )=1183876758; IS2(528 )= 436794231
      IS1(529 )=  68594352; IS2(529 )=1019241011
      IS1(530 )=1445194529; IS2(530 )=1652427192
      IS1(531 )=1984022317; IS2(531 )= 634620520
      IS1(532 )=1581637534; IS2(532 )= 715105076
      IS1(533 )=1852560766; IS2(533 )=1849389822
      IS1(534 )=1786797677; IS2(534 )=1775016593
      IS1(535 )= 282102855; IS2(535 )= 812059288
      IS1(536 )=2120900530; IS2(536 )=1367713004
      IS1(537 )= 483020346; IS2(537 )= 224667086
      IS1(538 )=1686388582; IS2(538 )=1056110078
      IS1(539 )=1080296513; IS2(539 )=1000944206
      IS1(540 )=1128089199; IS2(540 )= 211017438
      IS1(541 )=1873948292; IS2(541 )=1459653839
      IS1(542 )=1478295042; IS2(542 )=2091011370
      IS1(543 )= 780336939; IS2(543 )=1878440217
      IS1(544 )=  89453090; IS2(544 )= 496618994
      IS1(545 )=  84723786; IS2(545 )= 193890549
      IS1(546 )= 916751048; IS2(546 )= 529107147
      IS1(547 )= 937673116; IS2(547 )= 141932466
      IS1(548 )=1522160349; IS2(548 )= 617355429
      IS1(549 )=1323199052; IS2(549 )= 759342156
      IS1(550 )=1238578537; IS2(550 )=1328836664
      IS1(551 )=1112923558; IS2(551 )=1026465383
      IS1(552 )= 505271372; IS2(552 )=1384592778
      IS1(553 )= 566594858; IS2(553 )= 934194365
      IS1(554 )=1039240942; IS2(554 )=2009330113
      IS1(555 )=1832500056; IS2(555 )=1947683833
      IS1(556 )= 360636801; IS2(556 )=1001734064
      IS1(557 )= 669874416; IS2(557 )=1595554733
      IS1(558 )=1832363129; IS2(558 )=1805004882
      IS1(559 )=1913361231; IS2(559 )=1861860225
      IS1(560 )=1042646163; IS2(560 )=1027660606
      IS1(561 )=1997920265; IS2(561 )=1340547150
      IS1(562 )= 212027369; IS2(562 )=1980768873
      IS1(563 )= 955366244; IS2(563 )=1250568151
      IS1(564 )=  62449968; IS2(564 )=2054390312
      IS1(565 )= 407925365; IS2(565 )= 387411782
      IS1(566 )= 136158279; IS2(566 )= 868480095
      IS1(567 )= 383555850; IS2(567 )=1375373714
      IS1(568 )=1578435951; IS2(568 )= 212038261
      IS1(569 )= 358815004; IS2(569 )=1539319712
      IS1(570 )= 207947798; IS2(570 )=  84668320
      IS1(571 )=1234390761; IS2(571 )= 410483487
      IS1(572 )=1493324788; IS2(572 )= 487843369
      IS1(573 )=2056136989; IS2(573 )=1938329879
      IS1(574 )=1636698515; IS2(574 )= 698986119
      IS1(575 )= 304252344; IS2(575 )=1411426968
      IS1(576 )= 541389504; IS2(576 )= 660279563
      IS1(577 )=1976796554; IS2(577 )= 721255515
      IS1(578 )=1994743666; IS2(578 )= 260308500
      IS1(579 )=1129912793; IS2(579 )=  70624725
      IS1(580 )= 394341156; IS2(580 )= 160687023
      IS1(581 )= 698480071; IS2(581 )=1131283564
      IS1(582 )=  15953128; IS2(582 )=1327808851
      IS1(583 )=2047761093; IS2(583 )=1655969917
      IS1(584 )=1514379033; IS2(584 )= 628302318
      IS1(585 )= 455080996; IS2(585 )=2123515005
      IS1(586 )=1145649301; IS2(586 )=1563293737
      IS1(587 )=1631296493; IS2(587 )=1232136613
      IS1(588 )= 376739480; IS2(588 )= 694184162
      IS1(589 )=  62700700; IS2(589 )=1621822750
      IS1(590 )=1443922028; IS2(590 )= 394412632
      IS1(591 )=2034833661; IS2(591 )=  56369217
      IS1(592 )= 179925148; IS2(592 )=  46595906
      IS1(593 )= 943763422; IS2(593 )= 359712423
      IS1(594 )= 294551736; IS2(594 )=2091885037
      IS1(595 )=1612117876; IS2(595 )=1416720025
      IS1(596 )=2060353278; IS2(596 )=1999037873
      IS1(597 )=2027668315; IS2(597 )=1597386667
      IS1(598 )= 565740086; IS2(598 )=1681038972
      IS1(599 )=  85508641; IS2(599 )= 422362053
      IS1(600 )= 632136951; IS2(600 )= 967303111
      IS1(601 )= 287527126; IS2(601 )=2006085515
      IS1(602 )= 872188325; IS2(602 )= 973826090
      IS1(603 )=1528678858; IS2(603 )= 870128621
      IS1(604 )=2133000311; IS2(604 )=1630341425
      IS1(605 )=   9267403; IS2(605 )=  82683218
      IS1(606 )=1055725918; IS2(606 )=1071217436
      IS1(607 )=1601015878; IS2(607 )= 460901436
      IS1(608 )=1173244981; IS2(608 )= 575155810
      IS1(609 )= 190740363; IS2(609 )= 209424492
      IS1(610 )= 373750904; IS2(610 )=1037445515
      IS1(611 )=1903087315; IS2(611 )= 199145625
      IS1(612 )= 935216143; IS2(612 )= 767883538
      IS1(613 )= 616643835; IS2(613 )= 249787085
      IS1(614 )=1269743498; IS2(614 )=2038689023
      IS1(615 )=1906186820; IS2(615 )= 455180313
      IS1(616 )=1341988859; IS2(616 )=1474088526
      IS1(617 )= 999778032; IS2(617 )=1005749219
      IS1(618 )= 620954154; IS2(618 )=1241433020
      IS1(619 )=1863565816; IS2(619 )=1065421906
      IS1(620 )=1497283145; IS2(620 )=2115805116
      IS1(621 )= 406679876; IS2(621 )= 510949186
      IS1(622 )= 912988730; IS2(622 )=  49564907
      IS1(623 )=1764585681; IS2(623 )= 313681775
      IS1(624 )= 493773352; IS2(624 )= 791156315
      IS1(625 )=2086722153; IS2(625 )= 107352467
      IS1(626 )=1286084681; IS2(626 )= 208200070
      IS1(627 )= 808109356; IS2(627 )=1731413771
      IS1(628 )=1486614524; IS2(628 )=2064085170
      IS1(629 )=2129645386; IS2(629 )=1284378691
      IS1(630 )=1014423030; IS2(630 )= 701300786
      IS1(631 )=  60658846; IS2(631 )= 397903983
      IS1(632 )=1048804021; IS2(632 )=1593351979
      IS1(633 )=1842774497; IS2(633 )=1285982663
      IS1(634 )= 854481213; IS2(634 )=1322235277
      IS1(635 )= 798953202; IS2(635 )=1538431839
      IS1(636 )=1052424288; IS2(636 )=1822372594
      IS1(637 )= 601068089; IS2(637 )=1425648861
      IS1(638 )=1849885612; IS2(638 )= 733138526
      IS1(639 )=1269451977; IS2(639 )= 917070258
      IS1(640 )=2144361786; IS2(640 )=1842934549
      IS1(641 )=1496984691; IS2(641 )=1297819612
      IS1(642 )= 438239309; IS2(642 )=1148023460
      IS1(643 )=1737313181; IS2(643 )= 645117283
      IS1(644 )= 430834434; IS2(644 )= 284660368
      IS1(645 )=1143775914; IS2(645 )=1381744399
      IS1(646 )=1995585326; IS2(646 )=2120820043
      IS1(647 )= 135161363; IS2(647 )=1799867404
      IS1(648 )= 836364692; IS2(648 )=1029908703
      IS1(649 )=1757227577; IS2(649 )= 546266239
      IS1(650 )= 623798852; IS2(650 )= 585377965
      IS1(651 )=1098341577; IS2(651 )= 636341909
      IS1(652 )=1520291052; IS2(652 )= 913917239
      IS1(653 )= 975535163; IS2(653 )=1598597453
      IS1(654 )=1447444469; IS2(654 )=1937942110
      IS1(655 )=1432720174; IS2(655 )= 839452541
      IS1(656 )=1295341632; IS2(656 )=2103887297
      IS1(657 )=1300522179; IS2(657 )= 809881664
      IS1(658 )=1425114463; IS2(658 )=1208521323
      IS1(659 )=1257800427; IS2(659 )=1205340370
      IS1(660 )=1186021021; IS2(660 )=1115454768
      IS1(661 )=1164381967; IS2(661 )=1695441674
      IS1(662 )=1339932055; IS2(662 )= 626668376
      IS1(663 )=1054295865; IS2(663 )=1077719907
      IS1(664 )= 908887630; IS2(664 )= 375160191
      IS1(665 )=1172794729; IS2(665 )=1568832201
      IS1(666 )= 192588897; IS2(666 )= 917870514
      IS1(667 )= 643461236; IS2(667 )= 484049433
      IS1(668 )= 315353239; IS2(668 )= 397852714
      IS1(669 )= 931623820; IS2(669 )=1720458459
      IS1(670 )=1580993163; IS2(670 )= 734491598
      IS1(671 )=1513146206; IS2(671 )=1441389702
      IS1(672 )= 701699279; IS2(672 )=1253582219
      IS1(673 )=1597280672; IS2(673 )= 920294785
      IS1(674 )=2083951540; IS2(674 )=1899768161
      IS1(675 )=1604320624; IS2(675 )= 154338943
      IS1(676 )=1832407686; IS2(676 )=2039753469
      IS1(677 )=1522755360; IS2(677 )=1558890156
      IS1(678 )=1696355490; IS2(678 )=1445912335
      IS1(679 )=1308491933; IS2(679 )=1267480279
      IS1(680 )=1946363807; IS2(680 )=1896823905
      IS1(681 )=2028942171; IS2(681 )=1457946439
      IS1(682 )= 510086891; IS2(682 )= 536540300
      IS1(683 )=1486809993; IS2(683 )=1547406458
      IS1(684 )=1317061925; IS2(684 )=1591791104
      IS1(685 )= 989995158; IS2(685 )=1604821909
      IS1(686 )= 326298337; IS2(686 )=1171640652
      IS1(687 )= 152659804; IS2(687 )= 495837027
      IS1(688 )=1222941036; IS2(688 )= 796199205
      IS1(689 )= 907350872; IS2(689 )= 727966425
      IS1(690 )=1686084314; IS2(690 )=1613446594
      IS1(691 )=1642955746; IS2(691 )=1880495478
      IS1(692 )= 247284465; IS2(692 )=1706304962
      IS1(693 )=1611723103; IS2(693 )=1706608231
      IS1(694 )= 719689499; IS2(694 )=  31477369
      IS1(695 )=1079226399; IS2(695 )=1782920006
      IS1(696 )=1420885629; IS2(696 )=1072981519
      IS1(697 )=2084453400; IS2(697 )= 224386433
      IS1(698 )=1047203160; IS2(698 )= 614309249
      IS1(699 )= 586617884; IS2(699 )=1582667003
      IS1(700 )=  11984589; IS2(700 )=1558266095
      IS1(701 )= 755652237; IS2(701 )= 866088075
      IS1(702 )=2017712409; IS2(702 )=1954340696
      IS1(703 )= 589844984; IS2(703 )= 314030536
      IS1(704 )= 129229479; IS2(704 )= 272111716
      IS1(705 )= 870276471; IS2(705 )= 465006906
      IS1(706 )=1612019958; IS2(706 )=  21334935
      IS1(707 )=1315108425; IS2(707 )=1818825397
      IS1(708 )=1258602567; IS2(708 )=1066619326
      IS1(709 )= 166077102; IS2(709 )= 909070060
      IS1(710 )=1848198256; IS2(710 )= 385880783
      IS1(711 )= 488779996; IS2(711 )= 241674453
      IS1(712 )=1098298571; IS2(712 )=1348905710
      IS1(713 )= 561394508; IS2(713 )=1782802528
      IS1(714 )= 254791280; IS2(714 )= 354432011
      IS1(715 )=1214976755; IS2(715 )= 774386703
      IS1(716 )=  90779321; IS2(716 )=1418378337
      IS1(717 )=1355457939; IS2(717 )=1690663694
      IS1(718 )=1541432462; IS2(718 )=1040751740
      IS1(719 )= 394291901; IS2(719 )=1603027222
      IS1(720 )=1239001540; IS2(720 )= 146841931
      IS1(721 )=  38818735; IS2(721 )=1839009879
      IS1(722 )=1674954341; IS2(722 )=1920405940
      IS1(723 )=1640316191; IS2(723 )=1265350958
      IS1(724 )=1937228975; IS2(724 )= 996533450
      IS1(725 )= 423424194; IS2(725 )= 561329945
      IS1(726 )= 323756417; IS2(726 )=1682627390
      IS1(727 )=1549193098; IS2(727 )= 805029685
      IS1(728 )=1360079132; IS2(728 )= 925188637
      IS1(729 )=  78241293; IS2(729 )=1040585429
      IS1(730 )=1074975265; IS2(730 )= 797619830
      IS1(731 )=1582340080; IS2(731 )= 721022974
      IS1(732 )= 629043128; IS2(732 )= 610470736
      IS1(733 )=1657680582; IS2(733 )= 963797322
      IS1(734 )=1011918751; IS2(734 )=1560942165
      IS1(735 )=  42122891; IS2(735 )=1481369897
      IS1(736 )=2043026443; IS2(736 )=1863587933
      IS1(737 )=1332181389; IS2(737 )= 854812560
      IS1(738 )=1986410252; IS2(738 )=1040318983
      IS1(739 )= 521386266; IS2(739 )=1950110774
      IS1(740 )= 216896173; IS2(740 )=1685734611
      IS1(741 )=1867435681; IS2(741 )=2057062478
      IS1(742 )=1401811081; IS2(742 )=1014131666
      IS1(743 )=1538740757; IS2(743 )= 620335187
      IS1(744 )= 923455231; IS2(744 )=1901856320
      IS1(745 )= 597449802; IS2(745 )= 706565272
      IS1(746 )= 683455885; IS2(746 )=1033766701
      IS1(747 )= 828147504; IS2(747 )=1580010438
      IS1(748 )= 451283302; IS2(748 )= 781324118
      IS1(749 )=2048256218; IS2(749 )=1184873148
      IS1(750 )=1583574204; IS2(750 )=1032841150
      IS1(751 )=1652016656; IS2(751 )=  33056864
      IS1(752 )=1532136667; IS2(752 )=1581149947
      IS1(753 )=1462299459; IS2(753 )= 687082954
      IS1(754 )= 698944454; IS2(754 )=2122885802
      IS1(755 )=1195045208; IS2(755 )=1678424394
      IS1(756 )= 502707485; IS2(756 )=1444358879
      IS1(757 )= 941731361; IS2(757 )=1717503286
      IS1(758 )= 206769672; IS2(758 )=1166864868
      IS1(759 )=1312487368; IS2(759 )=1171921906
      IS1(760 )=  70908405; IS2(760 )=1544257314
      IS1(761 )= 666380866; IS2(761 )=1849868712
      IS1(762 )= 305790335; IS2(762 )=2037569416
      IS1(763 )= 334326322; IS2(763 )=1721074287
      IS1(764 )=2101405245; IS2(764 )=1240524239
      IS1(765 )= 760831995; IS2(765 )= 351651496
      IS1(766 )=2076975119; IS2(766 )= 732446407
      IS1(767 )=1457683031; IS2(767 )= 930496443
      IS1(768 )=1947001405; IS2(768 )= 590201248
      IS1(769 )= 275187592; IS2(769 )=  62024761
      IS1(770 )=1289777261; IS2(770 )= 542257293
      IS1(771 )= 363521237; IS2(771 )= 517555072
      IS1(772 )=1063822524; IS2(772 )=  90991909
      IS1(773 )= 312343146; IS2(773 )=1445115042
      IS1(774 )=2081825980; IS2(774 )=1071980321
      IS1(775 )=1790415941; IS2(775 )= 818819165
      IS1(776 )=  28581357; IS2(776 )= 877661732
      IS1(777 )=1803064654; IS2(777 )=1992656309
      IS1(778 )=1030876307; IS2(778 )= 513683199
      IS1(779 )= 492391370; IS2(779 )=1206877439
      IS1(780 )= 211680892; IS2(780 )=1783884173
      IS1(781 )= 744879183; IS2(781 )= 984195787
      IS1(782 )= 490863602; IS2(782 )=1663478249
      IS1(783 )= 756240663; IS2(783 )= 927897638
      IS1(784 )= 980102031; IS2(784 )=1228241271
      IS1(785 )=1517579622; IS2(785 )=1565288529
      IS1(786 )=2018942706; IS2(786 )=1486937165
      IS1(787 )= 914892050; IS2(787 )=1011671153
      IS1(788 )=2071502238; IS2(788 )= 910410508
      IS1(789 )=1144674111; IS2(789 )=1035198034
      IS1(790 )= 551571043; IS2(790 )=1704515839
      IS1(791 )=2067235260; IS2(791 )=  62628367
      IS1(792 )=1183655237; IS2(792 )=1825596188
      IS1(793 )=1720738448; IS2(793 )=1426908311
      IS1(794 )=1428394554; IS2(794 )=1331528227
      IS1(795 )=1236406902; IS2(795 )= 739342228
      IS1(796 )= 314637105; IS2(796 )=1635162390
      IS1(797 )=1631561757; IS2(797 )=1259987681
      IS1(798 )=1113570671; IS2(798 )= 941650185
      IS1(799 )=1316684934; IS2(799 )=1653069887
      IS1(800 )=1026427146; IS2(800 )= 713191356
      IS1(801 )=1616622265; IS2(801 )=1445073589
      IS1(802 )=1184120268; IS2(802 )= 844684601
      IS1(803 )=1916487469; IS2(803 )= 957905046
      IS1(804 )=1196992255; IS2(804 )= 143852508
      IS1(805 )= 810274414; IS2(805 )= 670905422
      IS1(806 )= 994554837; IS2(806 )=1752024033
      IS1(807 )= 631323451; IS2(807 )=1938843572
      IS1(808 )=2108796165; IS2(808 )= 686849224
      IS1(809 )= 675199957; IS2(809 )=2066177454
      IS1(810 )= 500608835; IS2(810 )= 501044809
      IS1(811 )= 318482166; IS2(811 )= 929306814
      IS1(812 )=1833923717; IS2(812 )= 224005423
      IS1(813 )=1473405260; IS2(813 )= 584253050
      IS1(814 )=1922717185; IS2(814 )= 725230049
      IS1(815 )= 358747736; IS2(815 )=1894346138
      IS1(816 )=1262935388; IS2(816 )=1123083929
      IS1(817 )= 170752885; IS2(817 )= 282349087
      IS1(818 )=1766873876; IS2(818 )=1639448540
      IS1(819 )=1327238154; IS2(819 )=2086656898
      IS1(820 )= 608102000; IS2(820 )=1953835572
      IS1(821 )=1933730221; IS2(821 )= 592600179
      IS1(822 )=  20093493; IS2(822 )=1802818520
      IS1(823 )=2047560831; IS2(823 )=1690864475
      IS1(824 )=2120624346; IS2(824 )=1399722254
      IS1(825 )= 728200766; IS2(825 )=1251271247
      IS1(826 )=1996802725; IS2(826 )=1182594334
      IS1(827 )=1732977781; IS2(827 )= 494629971
      IS1(828 )=1333989141; IS2(828 )=1463491202
      IS1(829 )= 414434960; IS2(829 )= 427148665
      IS1(830 )=2058685774; IS2(830 )=1258425844
      IS1(831 )=1570688571; IS2(831 )=1718768035
      IS1(832 )=1361600109; IS2(832 )= 732095097
      IS1(833 )=1464532988; IS2(833 )=1592327040
      IS1(834 )=1154029608; IS2(834 )= 457075509
      IS1(835 )= 504833970; IS2(835 )= 282459181
      IS1(836 )= 720633547; IS2(836 )=1754749459
      IS1(837 )=1542782084; IS2(837 )= 556876143
      IS1(838 )= 498983980; IS2(838 )=1448363633
      IS1(839 )= 815335207; IS2(839 )=1938426616
      IS1(840 )= 250652265; IS2(840 )=1560814150
      IS1(841 )=  61330462; IS2(841 )=1822327770
      IS1(842 )= 252780969; IS2(842 )=2058641830
      IS1(843 )=2103405990; IS2(843 )= 532243551
      IS1(844 )=1908261717; IS2(844 )=2064263729
      IS1(845 )=1977211665; IS2(845 )= 777733094
      IS1(846 )= 237237934; IS2(846 )= 528493150
      IS1(847 )=  73678412; IS2(847 )= 219112977
      IS1(848 )=1963696531; IS2(848 )=1929981191
      IS1(849 )= 101047790; IS2(849 )= 216743219
      IS1(850 )=   2093800; IS2(850 )=  74011539
      IS1(851 )= 632655512; IS2(851 )= 198607329
      IS1(852 )= 941788307; IS2(852 )=1227252584
      IS1(853 )= 545304972; IS2(853 )=1963545088
      IS1(854 )=  55853502; IS2(854 )= 498296470
      IS1(855 )=1891044982; IS2(855 )= 212219604
      IS1(856 )= 713967553; IS2(856 )=1094926756
      IS1(857 )= 126818203; IS2(857 )= 569771356
      IS1(858 )=1701715475; IS2(858 )= 840368587
      IS1(859 )= 904419474; IS2(859 )=1160456593
      IS1(860 )= 914755144; IS2(860 )=1602166631
      IS1(861 )= 719105598; IS2(861 )= 768277780
      IS1(862 )= 436565842; IS2(862 )=1526045329
      IS1(863 )= 642339928; IS2(863 )=1308502226
      IS1(864 )= 952737893; IS2(864 )= 889261161
      IS1(865 )= 750349739; IS2(865 )= 419661629
      IS1(866 )=1240092349; IS2(866 )= 385969468
      IS1(867 )= 947883679; IS2(867 )= 858658062
      IS1(868 )= 412164731; IS2(868 )= 227225801
      IS1(869 )=1077025161; IS2(869 )=1809495121
      IS1(870 )= 623779545; IS2(870 )=1748437918
      IS1(871 )=1462115422; IS2(871 )=1828054930
      IS1(872 )=1793988879; IS2(872 )= 971499218
      IS1(873 )= 401387846; IS2(873 )=2063529218
      IS1(874 )= 939523117; IS2(874 )= 466501459
      IS1(875 )=1325434731; IS2(875 )= 933144734
      IS1(876 )=1844466921; IS2(876 )=1723628496
      IS1(877 )= 612243172; IS2(877 )= 675935430
      IS1(878 )= 122322491; IS2(878 )= 700754464
      IS1(879 )=1118985067; IS2(879 )=  77681872
      IS1(880 )= 673808420; IS2(880 )=2027510724
      IS1(881 )=1900746742; IS2(881 )= 875745816
      IS1(882 )=1951431584; IS2(882 )=1326338561
      IS1(883 )=1185595160; IS2(883 )=1508083812
      IS1(884 )=  20387986; IS2(884 )=   2502650
      IS1(885 )=1448512176; IS2(885 )= 406078533
      IS1(886 )= 341831642; IS2(886 )= 444907341
      IS1(887 )=1749622481; IS2(887 )=1458991053
      IS1(888 )= 179921081; IS2(888 )=1671951076
      IS1(889 )= 928183806; IS2(889 )=1991458904
      IS1(890 )= 603248036; IS2(890 )=1824888100
      IS1(891 )=1063724212; IS2(891 )=1890874257
      IS1(892 )=  78830689; IS2(892 )=2043989239
      IS1(893 )=1753470474; IS2(893 )= 830139537
      IS1(894 )=1566789535; IS2(894 )=1388524135
      IS1(895 )=1518518357; IS2(895 )=1787113559
      IS1(896 )= 891266992; IS2(896 )=1954215738
      IS1(897 )= 725221948; IS2(897 )= 471108830
      IS1(898 )= 396459481; IS2(898 )= 560156443
      IS1(899 )=  54855828; IS2(899 )=1250479705
      IS1(900 )=1671412588; IS2(900 )= 238648368
      IS1(901 )=1457198819; IS2(901 )=1108923041
      IS1(902 )=1864168313; IS2(902 )=2034120136
      IS1(903 )= 737458311; IS2(903 )= 420578061
      IS1(904 )= 693032926; IS2(904 )=1415068309
      IS1(905 )= 549217560; IS2(905 )= 285822649
      IS1(906 )=2028105476; IS2(906 )= 257057932
      IS1(907 )=   5253877; IS2(907 )= 467436652
      IS1(908 )= 124062287; IS2(908 )= 913845906
      IS1(909 )= 983237044; IS2(909 )=1573260196
      IS1(910 )=1688115577; IS2(910 )=1088025993
      IS1(911 )= 456193194; IS2(911 )=1237163793
      IS1(912 )= 818167884; IS2(912 )=2137528272
      IS1(913 )=    670825; IS2(913 )= 326538226
      IS1(914 )=1256025768; IS2(914 )=1517888550
      IS1(915 )=1603246774; IS2(915 )= 955989622
      IS1(916 )=1723321062; IS2(916 )= 756662277
      IS1(917 )= 548919448; IS2(917 )= 670807456
      IS1(918 )= 765198119; IS2(918 )=1636394764
      IS1(919 )=1940460545; IS2(919 )=1546053813
      IS1(920 )=1253731346; IS2(920 )= 462240916
      IS1(921 )= 452873281; IS2(921 )=1640963727
      IS1(922 )= 563455527; IS2(922 )= 678314347
      IS1(923 )= 296152006; IS2(923 )=1921509503
      IS1(924 )=1488389520; IS2(924 )=2075130919
      IS1(925 )=1416431702; IS2(925 )= 400649975
      IS1(926 )= 187540584; IS2(926 )= 280946768
      IS1(927 )= 265064824; IS2(927 )= 292386889
      IS1(928 )=  17812467; IS2(928 )= 158666141
      IS1(929 )= 678154378; IS2(929 )= 273032591
      IS1(930 )=1008995524; IS2(930 )= 225822851
      IS1(931 )= 932640361; IS2(931 )=1794357721
      IS1(932 )=1256165221; IS2(932 )=1362164825
      IS1(933 )=1566130220; IS2(933 )= 880937875
      IS1(934 )=1719055144; IS2(934 )=1475807990
      IS1(935 )=1511255203; IS2(935 )=1398326134
      IS1(936 )=  40137809; IS2(936 )= 844274026
      IS1(937 )=1514777596; IS2(937 )=1176190656
      IS1(938 )= 862983829; IS2(938 )=1444976675
      IS1(939 )=1419376580; IS2(939 )= 897595798
      IS1(940 )=2004253892; IS2(940 )= 394099144
      IS1(941 )=1208598747; IS2(941 )= 710801376
      IS1(942 )=2134068686; IS2(942 )=1199555698
      IS1(943 )= 383571065; IS2(943 )=1799022351
      IS1(944 )=1896509659; IS2(944 )=1062349110
      IS1(945 )=2089619889; IS2(945 )= 103925910
      IS1(946 )=  15951717; IS2(946 )=1115129379
      IS1(947 )=1209658212; IS2(947 )=1067867972
      IS1(948 )=1533825359; IS2(948 )= 707403901
      IS1(949 )= 878388702; IS2(949 )= 641697588
      IS1(950 )= 370529890; IS2(950 )= 306771059
      IS1(951 )=1054078483; IS2(951 )=1574620834
      IS1(952 )=1512913863; IS2(952 )= 619850280
      IS1(953 )= 516041141; IS2(953 )=1250978720
      IS1(954 )= 964083750; IS2(954 )=1794316764
      IS1(955 )=1222510673; IS2(955 )=1010199648
      IS1(956 )= 716095488; IS2(956 )= 890400554
      IS1(957 )= 257132284; IS2(957 )= 233348130
      IS1(958 )=1317716326; IS2(958 )=2012396977
      IS1(959 )=1982982843; IS2(959 )= 584095052
      IS1(960 )=1376656999; IS2(960 )=1982442382
      IS1(961 )=1158326173; IS2(961 )= 841561919
      IS1(962 )= 223831123; IS2(962 )= 115318433
      IS1(963 )=1612316526; IS2(963 )= 190139923
      IS1(964 )= 811660944; IS2(964 )=1411875518
      IS1(965 )= 461792825; IS2(965 )=1500187934
      IS1(966 )=2080858235; IS2(966 )=2112319914
      IS1(967 )=  34295036; IS2(967 )= 281625837
      IS1(968 )=1796061661; IS2(968 )=1731981711
      IS1(969 )= 338906098; IS2(969 )= 474926566
      IS1(970 )=1049296891; IS2(970 )=1071577019
      IS1(971 )=1254386847; IS2(971 )= 867276511
      IS1(972 )= 823835412; IS2(972 )=1169979512
      IS1(973 )=  12773203; IS2(973 )=1151341486
      IS1(974 )= 145388856; IS2(974 )=1918711308
      IS1(975 )=1163564225; IS2(975 )= 588792593
      IS1(976 )=1585160972; IS2(976 )=1831535360
      IS1(977 )= 151747303; IS2(977 )= 893561329
      IS1(978 )=2103142931; IS2(978 )=1524302572
      IS1(979 )=1076384122; IS2(979 )=1181095863
      IS1(980 )=1047454590; IS2(980 )= 212580845
      IS1(981 )=1391900713; IS2(981 )=1595457237
      IS1(982 )=1980182019; IS2(982 )=2099452891
      IS1(983 )=1134036187; IS2(983 )=1086465811
      IS1(984 )=1561015123; IS2(984 )=2087543080
      IS1(985 )=1091162968; IS2(985 )=2039625632
      IS1(986 )= 905907344; IS2(986 )=1707901257
      IS1(987 )=1097576632; IS2(987 )= 399477627
      IS1(988 )= 854163724; IS2(988 )= 727734495
      IS1(989 )= 993622338; IS2(989 )=1208219252
      IS1(990 )=1137639335; IS2(990 )=1202540324
      IS1(991 )=1405884001; IS2(991 )=1769921362
      IS1(992 )=1002962410; IS2(992 )= 726546017
      IS1(993 )=1906752935; IS2(993 )=1700064870
      IS1(994 )= 547343311; IS2(994 )= 947713211
      IS1(995 )= 154484285; IS2(995 )=1863597615
      IS1(996 )=1283816755; IS2(996 )= 282974211
      IS1(997 )=1206555620; IS2(997 )=1748596686
      IS1(998 )= 830620896; IS2(998 )= 922833376
      IS1(999 )= 470833195; IS2(999 )=2082830087
      IS1(1000)=1275774316; IS2(1000)= 515177342
      IS1(1001)= 731077242; IS2(1001)=2126489340
C
      IF(N.GT.0.AND.N.LT.1002) THEN
        ISEED1=IS1(N)
        ISEED2=IS2(N)
      ELSE
        ISEED1=1
        ISEED2=1
      ENDIF
C
      RETURN
      END
C  *********************************************************************
C                         FUNCTION RAND (Random number generator)
C  *********************************************************************
      FUNCTION RAND(DUMMY)
C
C  This is an adapted version of subroutine RANECU written by F. James
C  (Comput. Phys. Commun. 60 (1990) 329-344), which has been modified to
C  give a single random number at each call.
C
C  The 'seeds' ISEED1 and ISEED2 must be initialised in the main program
C  and transferred through the named common block /RSEED/.
C
C  Some compilers incorporate an intrinsic random number generator with
C  the same name (but with different argument lists). To avoid conflict,
C  it is advisable to declare RAND as an external function in all sub-
C  programs that call it.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      PARAMETER (USCALE=1.0D0/2.147483563D9)
      COMMON/RSEED/ISEED1,ISEED2
C
      I1=ISEED1/53668
      ISEED1=40014*(ISEED1-I1*53668)-I1*12211
      IF(ISEED1.LT.0) ISEED1=ISEED1+2147483563
C
      I2=ISEED2/52774
      ISEED2=40692*(ISEED2-I2*52774)-I2*3791
      IF(ISEED2.LT.0) ISEED2=ISEED2+2147483399
C
      IZ=ISEED1-ISEED2
      IF(IZ.LT.1) IZ=IZ+2147483562
      RAND=IZ*USCALE
C
      RETURN
      END
