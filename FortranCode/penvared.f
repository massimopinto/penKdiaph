CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    PPPPP   EEEEEE  N    N  V    V    AA    RRRRR   EEEEEE  DDDDD     C
C    P    P  E       NN   N  V    V   A  A   R    R  E       D    D    C
C    P    P  E       N N  N  V    V  A    A  R    R  E       D    D    C
C    PPPPP   EEEE    N  N N  V    V  AAAAAA  RRRRR   EEEE    D    D    C
C    P       E       N   NN   V  V   A    A  R  R    E       D    D    C
C    P       EEEEEE  N    N    VV    A    A  R   R   EEEEEE  DDDDD     C
C                                                                      C
C                                                   (version 2014).    C
C                                                                      C
C     The present routines permit to apply basic variance-reduction    C
C  methods with PENELOPE.                                              C
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
C  software for any purpose. It is provided 'as is' without express    C
C  or implied warranty.                                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The variance-reduction methods considered here are designed to speed
C  up the calculation of local quantities (e.g., the average deposited
C  energy or the average particle fluence) that are obtained by scoring
C  contributions of individual interaction events. Other quantities,
C  such as the distribution of energy deposited in a body (which, in the
C  case of a scintillation detector, may be transformed into the pulse-
C  height spectrum) are determined by contributions of complete showers,
C  and require strict control of the shower evolution through the whole
C  geometry.
C
C
C  ---->   SUBROUTINE VSPLIT(NSPLIT)
C  This subroutine splits the current particle into NSPLIT identical
C  particles, defines their weights appropriately, and stores NSPLIT-1
C  of them into the secondary stack. The current particle continues with
C  a reduced weight.
C  Note: NSPLIT must be larger than one (NSPLIT>1).
C
C
C  ---->   SUBROUTINE VRR(PSURV)
C          SUBROUTINE VKILL(PKILL)
C  These subroutines apply the Russian roulette technique. The particle
C  is killed with probability PKILL=1-PSURV; if it survives, its weight
C  is increased by a factor 1/PSURV = 1/(1-PKILL).
C  Note: PSURV and PKILL must be larger than zero and less than one.
C
C
C  ---->   SUBROUTINE JUMPF(DSMAX,DS)
C          SUBROUTINE KNOCKF(DE,ICOL)
C  These two subroutines perform interaction forcing. Their action is to
C  artificially insert 'forced' interactions of selected kinds randomly
C  along the particle trajectory. This is accomplished by replacing the
C  inverse mean free path (IMFP) by a larger value, FORCE(.)*IMFP, where
C  FORCE(IBODY,KPAR,ICOL) is the forcing factor specified by the user.
C  Notice that the forcing factor must be larger than unity. To keep the
C  simulation unbiased, interactions are allowed to affect the state of
C  the projectile only with probability WFORCE=1/FORCE(.), which is less
C  than unity, and, at the same time, secondary particles generated in
C  forced interactions are assigned a weight smaller than that of the
C  projectile by a factor =WFORCE.
C
C  To apply simulation forcing, the MAIN program must call subroutines
C  'JUMPF' and 'KNOCKF' instead of the usual subroutines 'JUMP' and
C  'KNOCK'. Moreover, subroutine START _must_ be called before starting
C  a track and after each interface crossing, even for photons.
C
C  Different kinds of interactions are identified by the integer label
C  ICOL, as indicated in the following table:
C
C     +----+-----------------+-----------------+-----------------+
C     |ICOL|electron (KPAR=1)|photon (KPAR=2)  |positron (KPAR=3)|
C     +----+-----------------+-----------------+-----------------+
C     | 1  |hinge            |Rayleigh         |hinge            |
C     +----+-----------------+-----------------+-----------------+
C     | 2  |elastic          |Compton          |elastic          |
C     +----+-----------------+-----------------+-----------------+
C     | 3  |inelastic        |photoabsorption  |inelastic        |
C     +----+-----------------+-----------------+-----------------+
C     | 4  |bremsstrahlung   |pair production  |bremsstrahlung   |
C     +----+-----------------+-----------------+-----------------+
C     | 5  |inner-shell ion. |not defined      |inner-shell ion. |
C     +----+-----------------+-----------------+-----------------+
C     | 6  |not defined      |not defined      |annihilation     |
C     +----+-----------------+-----------------+-----------------+
C     | 7  |delta scattering |delta scattering |delta scattering |
C     +----+-----------------+-----------------+-----------------+
C     | 8  |not defined      |not defined      |not defined      |
C     +----+-----------------+-----------------+-----------------+
C
C  The forcing factors FORCE(.) have to be specified by the user in the
C  MAIN program; they are transferred through the module PENVARED_mod
C  (see below). Forcing factors must be larger than, or equal to unity;
C  obviously, the value FORCE(.)=1.0D0 means 'no forcing'.
C
C  For the sake of simplicity, the forcing factors are considered to be
C  independent of the particle energy. Although this scheme is flexible
C  enough for many practical uses, the FORCE(.) values can also be
C  varied during the simulation. To ensure consistency, any modification
C  of the forcing parameters should be defined immediately after a call
C  to subroutine START.
C
C  Combined with interaction forcing, the subroutine KNOCKF can also
C  apply bremsstrahlung splitting. This technique consists of sampling a
C  'normal' photon in each radiative event, and then splitting it into a
C  number IBRSPL (>1) of photons, all them with the energy and polar
C  emission angle of the 'normal' photon, weight equal to 1/IBRSPL times
C  the 'normal' weight, and azimuthal emission angles sampled uniformly
C  in the interval from zero to 2*PI. This method is computationally
C  simple (a single DO loop) and very effective to increase the number
C  of photons emitted by electrons and positrons, because each 'split'
C  photon is obtained by only applying a rotation to the direction of
C  the normal photon.
C
C  The bremsstrahlung splitting numbers IBRSPL(IBODY) must be specified
C  by the user in the MAIN program, they are transferred through the
C  module PENVARED_mod. Splitting numbers must be equal to, or greater
C  than 1.
C
C
C  ---->   SUBROUTINE JUMPW(DS)
C          SUBROUTINE KNOCKW(DE,ICOL)
C  These two subroutines implement Woodcock's delta-scattering method
C  for photons. The method takes advantage of the high penetration of
C  photons to simplify the tracking of these particles through material
C  systems with complex geometries. Photons are transported freely
C  across the system using an augmented IMFP (= linear attenuation
C  coefficient), STMAX, which is larger than the actual IMFPs in all the
C  materials crossed by a trajectory ray. The event at the end of each
C  free flight may be either a real interaction (ICOL=1 to 4) or a delta
C  interaction (ICOL=7). Delta interactions occur with probability
C  1-(ST/STMAX), where ST is the actual IMFP in the current material.
C  Evidently, the simulation remains unbiased. This procedure avoids the
C  need of computing intersections of particle rays with interfaces, at
C  the expense of having to determine which material is at the end of
C  each free flight. Notice that the method will improve the efficiency
C  only for those geometries where locating a particle (i.e., finding
C  the material at its current position) is faster than normal tracking.
C
C  To apply the delta-scattering method, the MAIN program must call
C  subroutines 'JUMPW' and 'KNOCKW' instead of the usual subroutines
C  'JUMP' and 'KNOCK'. Moreover, subroutine START _must_ be called
C  before starting a photon track and after each interface crossing.
C
C  The augmented IMFP of photons with energy E is obtained by interpol-
C  ation of a table of values at the grid energies, which is stored in
C  the common block
C     COMMON/CWOODC/STMV(NEGP),STMC(NEGP),STMAX
C  We use the same interpolation formula as in subroutine GIMFP (see
C  the source file PENELOPE.F),
C     STMAX=EXP(STMV(KE)+(STMV(KE+1)-STMV(KE))*XEK)+STMC(KE)
C  The arrays STMV and STMC are calculated by subroutine JUMPW the
C  first time it is invoked. The value of STMAX is updated at each
C  call to subroutine JUMPW.
C
C  Woodcock's method can be applied together with interaction forcing.
C  In this case, the interaction forcing parameters must be loaded
C  before the first call to subroutine JUMPW. Of course, the augmented
C  IMFP has to be larger than the IMFP for the forced interactions.
C  Notice that, when the geometry is described with PENGEOM, the delta-
C  scattering method is incompatible with the use of impact detectors,
C  which require strict control of interface crossings.
C
C  Note: The parameter NMS in module PENELOPE_mod defines the maximum
C  number of secondary particles in the stack. The default value,
C  NMS=1000, which is large enough for most applications, may be
C  insufficient with heavy splitting and interaction forcing. To prevent
C  saturation of the secondary stack, it is advisable to apply these
C  techniques only when the particle weight is in a limited range. Wild
C  variations in particle weights should be avoided at any cost, because
C  they usually cause an increase of variance.
C


C  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      MODULE PENVARED_mod
C
C  ****  Variance-reduction parameters.
C
      SAVE  ! Saves all items in the module.
C  ----  Parameter values defined for each body. NBV is the maximum
C        number of bodies. When using PENGEOM, NBV must have the same
C        value as the parameter NB in module PENGEOM_mod.
      INTEGER*4, PARAMETER :: NBV=5000
C  ----  Forcing factors, FORCE(IBODY,KPAR,ICOL).
      DOUBLE PRECISION, DIMENSION (NBV,3,8) :: FORCE=1.0D0
      DOUBLE PRECISION :: WFORCE
C  ----  Bremsstrahlung splitting numbers, IBRSPL(IBODY).
      INTEGER*4, DIMENSION (NBV) :: IBRSPL=1
C  ----  Energy deposited in the last event (analogue simulation).
      DOUBLE PRECISION :: DEA
C
      END MODULE PENVARED_mod
C  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


C  *********************************************************************
C                       SUBROUTINE VSPLIT
C  *********************************************************************
      SUBROUTINE VSPLIT(NSPLIT)
C
C  This subroutine splits the current particle into NSPLIT identical
C  particles, defines their weights appropriately, and stores NSPLIT-1
C  of them into the secondary stack. The current particle continues with
C  the reduced weight.
C
C  NOTE: NSPLIT must be larger than one.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION ILBA(5)
C
      WGHT=WGHT/NSPLIT
      ILBA(1)=ILB(1)+1  ! Split particles are treated as secondaries.
      ILBA(2)=ILB(2)
      ILBA(3)=ILB(3)
      ILBA(4)=0
      ILBA(5)=ILB(5)
C  ****  Particles 2, ..., NSPLIT are stored in the secondary stack.
      DO I=2,NSPLIT
        CALL STORES(E,X,Y,Z,U,V,W,WGHT,KPAR,ILBA,IPOL)
      ENDDO
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE VRR
C  *********************************************************************
      SUBROUTINE VRR(PSURV)
C
C  This subroutine applies the Russian roulette technique. PSURV is the
C  survival probability; when the particle survives, its weight is
C  increased by a factor 1/PSURV.
C
C  NOTE: PSURV must be larger than zero and less than one.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      EXTERNAL RAND
C
      IF(RAND(1.0D0).GT.PSURV) THEN
        E=0.0D0  ! The particle is killed.
        WGHT=0.0D0
      ELSE
        WGHT=WGHT/PSURV
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE VKILL
C  *********************************************************************
      SUBROUTINE VKILL(PKILL)
C
C  This subroutine applies the Russian roulette technique. The particle
C  is killed with probability PKILL; if it survives, its weight is in-
C  creased by a factor 1/(1-PKILL).
C
C  NOTE: PKILL must be larger than zero and less than one.
C
      USE TRACK_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
C
      EXTERNAL RAND
C
      IF(RAND(1.0D0).LT.PKILL) THEN
        E=0.0D0
        WGHT=0.0D0
      ELSE
        WGHT=WGHT/(1.0D0-PKILL)
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE JUMPF
C  *********************************************************************
      SUBROUTINE JUMPF(DSMAX,DS)
C
C  Modified subroutine 'JUMP' for interaction forcing.
C
C  Calculation of the free path from the starting point to the position
C  of the next event and of the probabilities of occurrence of different
C  events.
C
C  Arguments:
C    DSMAX .... maximum allowed step length (input),
C    DS ....... segment length (output).
C
C  Output, through module PENELOPE_mod:
C    E0STEP ... energy at the beginning of the segment,
C    DESOFT ... energy loss due to soft interactions along the step,
C    SSOFT .... stopping power due to soft interactions,
C               = DESOFT/step_length.
C
      USE TRACK_mod
      USE PENELOPE_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LFORC
C  ****  Energy grid and interpolation constants for the current energy.
      COMMON/CEGRID/EMIN,EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  DW1EL(MAXMAT,NEGP),DW2EL(MAXMAT,NEGP),
     5  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     6  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  DW1PL(MAXMAT,NEGP),DW2PL(MAXMAT,NEGP),
     5  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     6  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C  ****  Current state and IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DSR,W1,W2,T1,T2
      COMMON/CJUMP1/ELAST1,ELAST2,MHINGE,KSOFTE,KSOFTI,KDELTA
C  ****  Interaction forcing parameters.
      COMMON/CFORCF/POR(8),P0(8),IBR,LFORC(8)
C
      EXTERNAL RAND
C
      IF(KPAR.EQ.1) THEN
C
C  ************  Electrons (KPAR=1).
C
        IF(MHINGE.EQ.1) THEN
          IF(E.LT.ELAST1) THEN
            XEL=LOG(E)
            XE=1.0D0+(XEL-DLEMP1)*DLFC
            KE=XE
            XEK=XE-KE
            CALL EIMFP(1)
            DO KCOL=2,5
              POR(KCOL)=P(KCOL)  ! Int. forcing modifies the IMFPs.
            ENDDO
            ELAST1=E
          ENDIF
          DS=DSR
          RETURN
        ENDIF
C
        E0STEP=E
        IF(E.LT.ELAST2) THEN
          XEL=LOG(E)
          XE=1.0D0+(XEL-DLEMP1)*DLFC
          KE=XE
          XEK=XE-KE
          CALL EIMFP(2)
          DO KCOL=2,5
            POR(KCOL)=P(KCOL)  ! Int. forcing modifies the IMFPs.
          ENDDO
          ELAST2=E
          ELAST1=E
        ENDIF
C  ****  Interaction forcing.
        TFP=0.0D0
        P0(2)=POR(2)
        P(2)=POR(2)
        DO KCOL=3,5
          IF(FORCE(IBODY,1,KCOL).GT.1.0D0.AND.P(KCOL).GT.1.0D-16) THEN
            P0(KCOL)=POR(KCOL)
            P(KCOL)=POR(KCOL)*FORCE(IBODY,1,KCOL)
            TFP=TFP+(POR(KCOL)-P0(KCOL))
            LFORC(KCOL)=.TRUE.
          ELSE
            LFORC(KCOL)=.FALSE.
          ENDIF
        ENDDO
C
C  ****  Inverse hard mean free path (interaction probability per unit
C        path length).
C
        ST=P(2)+P(3)+P(4)+P(5)+P(8)
        DSMAXP=DSMAX
C
C  ****  Soft stopping interactions.
C        KSOFTI=1, soft stopping is active,
C        KSOFTI=0, soft stopping is not active.
C
        IF(W1.GT.1.0D-20) THEN
          KSOFTI=1
C  ****  The maximum step length, DSMAXP, is determined in terms of the
C        input DSMAX value (which is specified by the user) and the mean
C        free path for hard interactions (1/ST).
          DSMC=4.0D0/ST
          IF(DSMAXP.GT.DSMC) THEN
            DSMAXP=DSMC
          ELSE IF(DSMAXP.LT.1.0D-8) THEN
            DSMAXP=DSMC
          ENDIF
C  ****  The value of DSMAXP is randomized to eliminate dose artifacts
C        at the end of the first step.
          DSMAXP=(0.5D0+RAND(1.0D0)*0.5D0)*DSMAXP
C
C  ****  Upper bound for the interaction probability along the step
C        (including soft energy straggling).
C
          EDE0=W1*DSMAXP
          VDE0=W2*DSMAXP
          FSEDE=MAX(1.0D0-DW1EL(MAT,KE)*EDE0,0.75D0)
          FSVDE=MAX(1.0D0-DW2EL(MAT,KE)*EDE0,0.75D0)
          EDEM=EDE0*FSEDE
          VDEM=VDE0*FSVDE
          W21=VDEM/EDEM
          IF(EDEM.GT.9.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+3.0D0*SQRT(VDEM)),EMIN)
          ELSE IF(EDEM.GT.3.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+SQRT(3.0D0*VDEM)),EMIN)
          ELSE
            ELOWER=MAX(E-1.5D0*(EDEM+W21),EMIN)
          ENDIF
          XE1=1.0D0+(LOG(ELOWER)-DLEMP1)*DLFC
          KE1=XE1
          XEK1=XE1-KE1
          STLWR=EXP(SETOT(MAT,KE1)+(SETOT(MAT,KE1+1)
     1      -SETOT(MAT,KE1))*XEK1)
          ST=MAX(ST,STLWR+TFP)
        ELSE
          KSOFTI=0
          DESOFT=0.0D0
          SSOFT=0.0D0
        ENDIF
C
C  ****  Soft elastic scattering.
C        KSOFTE=1, soft scattering is active,
C        KSOFTE=0, soft scattering is not active.
C
        IF(T1.GT.1.0D-20) THEN
          KSOFTE=1
        ELSE
          KSOFTE=0
        ENDIF
C
C  ****  Delta interactions.
C        KDELTA=0, a hard interaction follows,
C        KDELTA=1, a delta interaction follows.
C
        DST=-LOG(RAND(2.0D0))/ST
        IF(DST.LT.DSMAXP) THEN
          KDELTA=0
        ELSE
          DST=DSMAXP
          KDELTA=1
        ENDIF
C
        IF(KSOFTE+KSOFTI.EQ.0) THEN
          MHINGE=1
          DS=DST
        ELSE
          DS=DST*RAND(3.0D0)
          DSR=DST-DS
          IF(KSOFTI.EQ.1) THEN
            IF(DST.LT.1.0D-8) THEN
              SSOFT=W1
              DESOFT=SSOFT*DST
            ELSE
              EDE0=W1*DST
              VDE0=W2*DST
              FSEDE=MAX(1.0D0-DW1EL(MAT,KE)*EDE0,0.75D0)
              FSVDE=MAX(1.0D0-DW2EL(MAT,KE)*EDE0,0.75D0)
              EDE=EDE0*FSEDE
              VDE=VDE0*FSVDE
C  ****  Generation of random values DE with mean EDE and variance VDE.
              SIGMA=SQRT(VDE)
              IF(SIGMA.LT.0.333333333D0*EDE) THEN
C  ****  Truncated Gaussian distribution.
                DESOFT=EDE+RNDG3()*SIGMA
              ELSE
                RU=RAND(4.0D0)
                EDE2=EDE*EDE
                VDE3=3.0D0*VDE
                IF(EDE2.LT.VDE3) THEN
                  PNULL=(VDE3-EDE2)/(VDE3+3.0D0*EDE2)
                  IF(RU.LT.PNULL) THEN
                    DESOFT=0.0D0
                    SSOFT=0.0D0
                    IF(KSOFTE.EQ.0) THEN
                      MHINGE=1
                      DS=DST
                    ELSE
                      KSOFTI=0
                    ENDIF
                    RETURN
                  ELSE
C  ****  Uniform distribution.
                    DESOFT=1.5D0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0D0-PNULL)
                  ENDIF
                ELSE
                  DESOFT=EDE+(2.0D0*RU-1.0D0)*SQRT(VDE3)
                ENDIF
              ENDIF
              SSOFT=DESOFT/DST
            ENDIF
          ENDIF
        ENDIF
        RETURN
      ELSE IF(KPAR.EQ.3) THEN
C
C  ************  Positrons (KPAR=3).
C
        IF(MHINGE.EQ.1) THEN
          IF(E.LT.ELAST1) THEN
            XEL=LOG(E)
            XE=1.0D0+(XEL-DLEMP1)*DLFC
            KE=XE
            XEK=XE-KE
            CALL PIMFP(1)
            DO KCOL=2,6
              POR(KCOL)=P(KCOL)  ! Int. forcing modifies the IMFPs.
            ENDDO
            ELAST1=E
          ENDIF
          DS=DSR
          RETURN
        ENDIF
C
        E0STEP=E
        IF(E.LT.ELAST2) THEN
          XEL=LOG(E)
          XE=1.0D0+(XEL-DLEMP1)*DLFC
          KE=XE
          XEK=XE-KE
          CALL PIMFP(2)
          DO KCOL=2,6
            POR(KCOL)=P(KCOL)  ! Int. forcing modifies the IMFPs.
          ENDDO
          ELAST2=E
          ELAST1=E
        ENDIF
C  ****  Interaction forcing.
        TFP=0.0D0
        P0(2)=POR(2)
        P(2)=POR(2)
        DO KCOL=3,6
          IF(FORCE(IBODY,3,KCOL).GT.1.0D0.AND.P(KCOL).GT.1.0D-16) THEN
            P0(KCOL)=POR(KCOL)
            P(KCOL)=POR(KCOL)*FORCE(IBODY,3,KCOL)
            TFP=TFP+(POR(KCOL)-P0(KCOL))
            LFORC(KCOL)=.TRUE.
          ELSE
            LFORC(KCOL)=.FALSE.
          ENDIF
        ENDDO
C
C  ****  Inverse hard mean free path (interaction probability per unit
C        path length).
C
        ST=P(2)+P(3)+P(4)+P(5)+P(6)+P(8)
        DSMAXP=DSMAX
C
C  ****  Soft stopping interactions.
C        KSOFTI=1, soft stopping is active,
C        KSOFTI=0, soft stopping is not active.
C
        IF(W1.GT.1.0D-20) THEN
          KSOFTI=1
C  ****  The maximum step length, DSMAXP, is determined in terms of the
C        input DSMAX value (which is specified by the user) and the mean
C        free path for hard interactions (1/ST).
          DSMC=4.0D0/ST
          IF(DSMAXP.GT.DSMC) THEN
            DSMAXP=DSMC
          ELSE IF(DSMAXP.LT.1.0D-8) THEN
            DSMAXP=DSMC
          ENDIF
C  ****  The value of DSMAXP is randomized to eliminate dose artifacts
C        at the end of the first step.
          DSMAXP=(0.5D0+RAND(1.0D0)*0.5D0)*DSMAXP
C
C  ****  Upper bound for the interaction probability along the step
C        (including soft energy straggling).
C
          EDE0=W1*DSMAXP
          VDE0=W2*DSMAXP
          FSEDE=MAX(1.0D0-DW1PL(MAT,KE)*EDE0,0.75D0)
          FSVDE=MAX(1.0D0-DW2PL(MAT,KE)*EDE0,0.75D0)
          EDEM=EDE0*FSEDE
          VDEM=VDE0*FSVDE
          W21=VDEM/EDEM
          IF(EDEM.GT.9.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+3.0D0*SQRT(VDEM)),EMIN)
          ELSE IF(EDEM.GT.3.0D0*W21) THEN
            ELOWER=MAX(E-(EDEM+SQRT(3.0D0*VDEM)),EMIN)
          ELSE
            ELOWER=MAX(E-1.5D0*(EDEM+W21),EMIN)
          ENDIF
          XE1=1.0D0+(LOG(ELOWER)-DLEMP1)*DLFC
          KE1=XE1
          XEK1=XE1-KE1
          STLWR=EXP(SPTOT(MAT,KE1)+(SPTOT(MAT,KE1+1)
     1      -SPTOT(MAT,KE1))*XEK1)
          ST=MAX(ST,STLWR+TFP)
        ELSE
          KSOFTI=0
          DESOFT=0.0D0
          SSOFT=0.0D0
        ENDIF
C
C  ****  Soft elastic scattering.
C        KSOFTE=1, soft scattering is active,
C        KSOFTE=0, soft scattering is not active.
C
        IF(T1.GT.1.0D-20) THEN
          KSOFTE=1
        ELSE
          KSOFTE=0
        ENDIF
C
C  ****  Delta interactions.
C        KDELTA=0, a hard interaction follows,
C        KDELTA=1, a delta interaction follows.
C
        DST=-LOG(RAND(2.0D0))/ST
        IF(DST.LT.DSMAXP) THEN
          KDELTA=0
        ELSE
          DST=DSMAXP
          KDELTA=1
        ENDIF
C
        IF(KSOFTE+KSOFTI.EQ.0) THEN
          MHINGE=1
          DS=DST
        ELSE
          DS=DST*RAND(3.0D0)
          DSR=DST-DS
          IF(KSOFTI.EQ.1) THEN
            IF(DST.LT.1.0D-8) THEN
              SSOFT=W1
              DESOFT=SSOFT*DST
            ELSE
              EDE0=W1*DST
              VDE0=W2*DST
              FSEDE=MAX(1.0D0-DW1PL(MAT,KE)*EDE0,0.75D0)
              FSVDE=MAX(1.0D0-DW2PL(MAT,KE)*EDE0,0.75D0)
              EDE=EDE0*FSEDE
              VDE=VDE0*FSVDE
C  ****  Generation of random values DE with mean EDE and variance VDE.
              SIGMA=SQRT(VDE)
              IF(SIGMA.LT.0.333333333D0*EDE) THEN
C  ****  Truncated Gaussian distribution.
                DESOFT=EDE+RNDG3()*SIGMA
              ELSE
                RU=RAND(4.0D0)
                EDE2=EDE*EDE
                VDE3=3.0D0*VDE
                IF(EDE2.LT.VDE3) THEN
                  PNULL=(VDE3-EDE2)/(VDE3+3.0D0*EDE2)
                  IF(RU.LT.PNULL) THEN
                    DESOFT=0.0D0
                    SSOFT=0.0D0
                    IF(KSOFTE.EQ.0) THEN
                      MHINGE=1
                      DS=DST
                    ELSE
                      KSOFTI=0
                    ENDIF
                    RETURN
                  ELSE
C  ****  Uniform distribution.
                    DESOFT=1.5D0*(EDE+VDE/EDE)*(RU-PNULL)/(1.0D0-PNULL)
                  ENDIF
                ELSE
                  DESOFT=EDE+(2.0D0*RU-1.0D0)*SQRT(VDE3)
                ENDIF
              ENDIF
              SSOFT=DESOFT/DST
            ENDIF
          ENDIF
        ENDIF
        RETURN
      ELSE
C
C  ************  Photons (KPAR=2).
C
        IF(E.LT.ELAST1) THEN
          XEL=LOG(E)
          XE=1.0D0+(XEL-DLEMP1)*DLFC
          KE=XE
          XEK=XE-KE
          CALL GIMFP
          ELAST1=E
C  ****  Interaction forcing.
          P0(1)=P(1)
          DO KCOL=2,4
            IF(FORCE(IBODY,2,KCOL).GT.1.0D0.AND.P(KCOL).GT.1.0D-16) THEN
              P0(KCOL)=P(KCOL)
              P(KCOL)=P(KCOL)*FORCE(IBODY,2,KCOL)
              LFORC(KCOL)=.TRUE.
            ELSE
              LFORC(KCOL)=.FALSE.
            ENDIF
          ENDDO
          ST=P(1)+P(2)+P(3)+P(4)+P(8)
        ENDIF
C
        DS=-LOG(RAND(2.0D0))/ST
      ENDIF
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE KNOCKF
C  *********************************************************************
      SUBROUTINE KNOCKF(DE,ICOL)
C
C  Modified subroutine 'KNOCK' for interaction forcing.
C
C  Simulation of random hinges and hard interaction events.
C
C  Output arguments:
C    DE ....... net energy deposited in the event, = DEA*WFORCE.
C               The weight WGHT of the running particle is not modified.
C               The effective (weighted) energy deposited in the event
C               is DE*WGHT.
C               --> Use with care: simulated energy deposition spectra
C               are biased (because particle energy losses and deposited
C               energies are not balanced).
C    ICOL ..... kind of interaction suffered by the particle.
C
C  Output, through module PENVARED_mod:
C    WFORCE ... relative weight of the secondary particles released in
C               the interaction, equal to 1/FORCE(.).
C    DEA ...... energy deposited in the interaction of a particle with
C               unit weight (analogue simulation).
C
      USE TRACK_mod
      USE PENELOPE_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LFORC,LCOL
      CHARACTER*2 LASYMB
      PARAMETER (PI=3.1415926535897932D0, TWOPI=PI+PI)
      PARAMETER (REV=5.10998928D5)  ! Electron rest energy (eV)
      PARAMETER (RREV=1.0D0/REV, TREV=2.0D0*REV)
C
      COMMON/CHIST/ILBA(5)
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZ(MAXMAT,30),NELEM(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      COMMON/CEGRID/EMIN,EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Element data.
      COMMON/CADATA/ATW(99),EPX(99),RSCR(99),ETA(99),EB(99,30),
     1  ALW(99,30),CP0(99,30),IFI(99,30),IKS(99,30),NSHT(99),LASYMB(99)
C  ****  E/P inelastic collisions.
      PARAMETER (NO=512)
      COMMON/CEIN/EXPOT(MAXMAT),OP2(MAXMAT),F(MAXMAT,NO),UI(MAXMAT,NO),
     1  WRI(MAXMAT,NO),KZ(MAXMAT,NO),KS(MAXMAT,NO),NOSC(MAXMAT)
C  ****  Compton scattering.
      PARAMETER (NOCO=512)
      COMMON/CGCO/FCO(MAXMAT,NOCO),UICO(MAXMAT,NOCO),FJ0(MAXMAT,NOCO),
     1  PTRSH(MAXMAT,NOCO),KZCO(MAXMAT,NOCO),KSCO(MAXMAT,NOCO),
     2  NOSCCO(MAXMAT)
C  ****  Bremsstrahlung emission.
      PARAMETER (NBW=32)
      COMMON/CEBR/WB(NBW),PBCUT(MAXMAT,NEGP),WBCUT(MAXMAT,NEGP),
     1  PDFB(MAXMAT,NEGP,NBW),DPDFB(MAXMAT,NEGP,NBW),
     2  PACB(MAXMAT,NEGP,NBW),ZBR2(MAXMAT)
C  ****  Electron simulation tables.
      COMMON/CEIMFP/SEHEL(MAXMAT,NEGP),SEHIN(MAXMAT,NEGP),
     1  SEISI(MAXMAT,NEGP),SEHBR(MAXMAT,NEGP),SEAUX(MAXMAT,NEGP),
     2  SETOT(MAXMAT,NEGP),CSTPE(MAXMAT,NEGP),RSTPE(MAXMAT,NEGP),
     3  DEL(MAXMAT,NEGP),W1E(MAXMAT,NEGP),W2E(MAXMAT,NEGP),
     4  DW1EL(MAXMAT,NEGP),DW2EL(MAXMAT,NEGP),
     5  RNDCE(MAXMAT,NEGP),AE(MAXMAT,NEGP),BE(MAXMAT,NEGP),
     6  T1E(MAXMAT,NEGP),T2E(MAXMAT,NEGP)
C  ****  Positron simulation tables.
      COMMON/CPIMFP/SPHEL(MAXMAT,NEGP),SPHIN(MAXMAT,NEGP),
     1  SPISI(MAXMAT,NEGP),SPHBR(MAXMAT,NEGP),SPAN(MAXMAT,NEGP),
     2  SPAUX(MAXMAT,NEGP),SPTOT(MAXMAT,NEGP),CSTPP(MAXMAT,NEGP),
     3  RSTPP(MAXMAT,NEGP),W1P(MAXMAT,NEGP),W2P(MAXMAT,NEGP),
     4  DW1PL(MAXMAT,NEGP),DW2PL(MAXMAT,NEGP),
     5  RNDCP(MAXMAT,NEGP),AP(MAXMAT,NEGP),BP(MAXMAT,NEGP),
     6  T1P(MAXMAT,NEGP),T2P(MAXMAT,NEGP)
C
      COMMON/CELSEP/EELMAX(MAXMAT),PELMAX(MAXMAT),
     1              RNDCEd(MAXMAT,NEGP),RNDCPd(MAXMAT,NEGP)
C  ****  Secondary stack.
      COMMON/CERSEC/IERSEC
C  ****  Current state and IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DSR,W1,W2,T1,T2
      COMMON/CJUMP1/ELAST1,ELAST2,MHINGE,KSOFTE,KSOFTI,KDELTA
      COMMON/CFORCF/POR(8),P0(8),IBR,LFORC(8)
C
      EXTERNAL RAND
C
C  ****  Allow writing multiple stack overflow warnings.
      IF(IERSEC.NE.0) IERSEC=0
C
      WFORCE=1.0D0
      IF(KPAR.EQ.1) THEN
        GO TO 1000
      ELSE IF(KPAR.EQ.2) THEN
        GO TO 2000
      ELSE IF(KPAR.EQ.3) THEN
        GO TO 3000
      ELSE
        STOP 'KNOCKF: Incorrect particle type.'
      ENDIF
C
C  ************  Electrons (KPAR=1).
C
 1000 CONTINUE
      IF(MHINGE.EQ.1) GO TO 1100
C
C  ****  Hinge, artificial soft event (ICOL=1).
C
      ICOL=1
      MHINGE=1
C
C  ****  Energy loss.
C
      IF(KSOFTI.EQ.1) THEN
        DE=DESOFT
        E=E-DE
        IF(E.LT.EABS(1,MAT)) THEN
          DE=E0STEP
          DEA=DE
          E=0.0D0
          RETURN
        ENDIF
        DEA=DE
        E0STEP=E0STEP-SSOFT*(DST-DSR)
        IF(KSOFTE.EQ.0) RETURN
        XEL=LOG(E0STEP)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
      ELSE
        DEA=0.0D0
        DE=0.0D0
      ENDIF
C
C  ****  Angular deflection.
C
      IF(T1E(MAT,KE+1).GT.-78.3D0) THEN
        T1=EXP(T1E(MAT,KE)+(T1E(MAT,KE+1)-T1E(MAT,KE))*XEK)
        T2=EXP(T2E(MAT,KE)+(T2E(MAT,KE+1)-T2E(MAT,KE))*XEK)
      ELSE
        T1=0.0D0
        T2=0.0D0
      ENDIF
      IF(T1.LT.1.0D-20) RETURN
C  ****  1st and 2nd moments of the angular distribution.
      EMU1=0.5D0*(1.0D0-EXP(-DST*T1))
      EMU2=EMU1-(1.0D0-EXP(-DST*T2))/6.0D0
C  ****  Sampling from a two-bar histogram with these moments.
      PNUM=2.0D0*EMU1-3.0D0*EMU2
      PDEN=1.0D0-2.0D0*EMU1
      PMU0=PNUM/PDEN
      PA=PDEN+PMU0
      RND=RAND(2.0D0)
      IF(RND.LT.PA) THEN
        CDT=1.0D0-2.0D0*PMU0*(RND/PA)
      ELSE
        CDT=1.0D0-2.0D0*(PMU0+(1.0D0-PMU0)*((RND-PA)/(1.0D0-PA)))
      ENDIF
      DF=TWOPI*RAND(3.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ************  Hard event.
C
 1100 CONTINUE
      MHINGE=0
C  ****  A delta interaction (ICOL=7) occurs when the maximum
C        allowed step length is exceeded.
      IF(KDELTA.EQ.1) THEN
        ICOL=7
        DEA=0.0D0
        DE=0.0D0
        RETURN
      ENDIF
C  ****  Interaction forcing.
      DO KCOL=3,5
        IF(FORCE(IBODY,1,KCOL).GT.1.0D0.AND.POR(KCOL).GT.1.0D-16) THEN
          P0(KCOL)=POR(KCOL)
          P(KCOL)=POR(KCOL)*FORCE(IBODY,1,KCOL)
          LFORC(KCOL)=.TRUE.
        ELSE
          LFORC(KCOL)=.FALSE.
        ENDIF
      ENDDO
      IBR=MAX(IBRSPL(IBODY),1)
C  ****  Random sampling of the interaction type.
      STNOW=P(2)+P(3)+P(4)+P(5)+P(8)
      STS=MAX(STNOW,ST)*RAND(4.0D0)
      SS=P(2)
      IF(SS.GT.STS) GO TO 1200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 1300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 1400
      SS=SS+P(5)
      IF(SS.GT.STS) GO TO 1500
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 1800
C  ****  A delta interaction (ICOL=7) may occur when the total
C        interaction probability per unit path length, ST, is
C        larger than STNOW.
      ICOL=7
      DEA=0.0D0
      DE=0.0D0
      RETURN
C
C  ****  Hard elastic collision (ICOL=2).
C
 1200 ICOL=2
      IF(E.GE.EELMAX(MAT)) THEN
        TRNDC=RNDCE(MAT,KE)+(RNDCE(MAT,KE+1)-RNDCE(MAT,KE))*XEK
        TA=EXP(AE(MAT,KE)+(AE(MAT,KE+1)-AE(MAT,KE))*XEK)
        TB=BE(MAT,KE)+(BE(MAT,KE+1)-BE(MAT,KE))*XEK
        CALL EELa(TA,TB,TRNDC,RMU)
      ELSE
        TRNDC=RNDCEd(MAT,KE)+(RNDCEd(MAT,KE+1)-RNDCEd(MAT,KE))*XEK
        CALL EELd(TRNDC,RMU)  ! Uses the ELSEPA database.
      ENDIF
      CDT=1.0D0-(RMU+RMU)
      DF=TWOPI*RAND(5.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      DEA=0.0D0
      DE=0.0D0
      RETURN
C
C  ****  Hard inelastic collision (ICOL=3).
C
 1300 ICOL=3
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      DELTA=DEL(MAT,KE)+(DEL(MAT,KE+1)-DEL(MAT,KE))*XEK
      CALL EINa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IOSC)
C  ****  Scattering angles (primary electron).
      DF=TWOPI*RAND(6.0D0)
C  ****  Delta ray.
      IF(ES.GT.EABS(1,MAT)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
C  ****  New energy and direction.
      DEA=DE
      DE=DEA*WFORCE
      IF(LCOL) THEN
        IF(EP.GT.EABS(1,MAT)) THEN
          E=EP
          CALL DIRECT(CDT,DF,U,V,W)
        ELSE
          DEA=DEA+EP
          DE=DE+EP
          E=0.0D0
        ENDIF
      ENDIF
      RETURN
C
C  ****  Hard bremsstrahlung emission (ICOL=4).
C
 1400 ICOL=4
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      CALL EBRa(E,DE,MAT)
C  ****  Bremsstrahlung photon.
      IF(DE.GT.EABS(2,MAT)) THEN
        CALL EBRaA(E,DE,CDTS,MAT)
        WSPLIT=WGHT*WFORCE/IBR
        DO I=1,IBR
          DFS=TWOPI*RAND(7.0D0)
          US=U
          VS=V
          WS=W
          CALL DIRECT(CDTS,DFS,US,VS,WS)
          ILBA(1)=ILB(1)+1
          ILBA(2)=KPAR
          ILBA(3)=ICOL
          ILBA(4)=0
          ILBA(5)=ILB(5)
          CALL STORES(DE,X,Y,Z,US,VS,WS,WSPLIT,2,ILBA,0)
        ENDDO
      ENDIF
C  ****  New energy.
      DEA=DE
      DE=DE*WFORCE
      IF(LCOL) THEN
        E=E-DEA
        IF(E.LT.EABS(1,MAT)) THEN
          DEA=DEA+E
          DE=DE+E
          E=0.0D0
        ENDIF
      ENDIF
      RETURN
C
C  ****  Ionization of an inner shell (ICOL=5).
C
 1500 ICOL=5
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      DELTA=DEL(MAT,KE)+(DEL(MAT,KE+1)-DEL(MAT,KE))*XEK
      CALL ESIa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IZA,ISA)
C  ****  Scattering angles (primary electron).
      DF=TWOPI*RAND(8.0D0)
C  ****  Delta ray.
      IF(ES.GT.EABS(1,MAT)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
C  ****  Atomic relaxation.
      IF(IZA.GT.2) THEN
        ILBA(3)=ICOL
        WGHTA=WGHT
        WGHT=WGHT*WFORCE
        CALL RELAX(IZA,ISA)
        WGHT=WGHTA
      ENDIF
C  ****  New energy and direction.
      DEA=DE
      DE=DEA*WFORCE
      IF(LCOL) THEN
        IF(EP.GT.EABS(1,MAT)) THEN
          E=EP
          CALL DIRECT(CDT,DF,U,V,W)
        ELSE
          DEA=DEA+EP
          DE=DE+EP
          E=0.0D0
        ENDIF
      ENDIF
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 1800 ICOL=8
      DEA=0.0D0
      DE=0.0D0
      CALL EAUX
      RETURN
C
C  ************  Photons (KPAR=2).
C
 2000 CONTINUE
C
      STS=ST*RAND(1.0D0)
      SS=P(1)
      IF(SS.GT.STS) GO TO 2100
      SS=SS+P(2)
      IF(SS.GT.STS) GO TO 2200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 2300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 2400
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 2800
C
C  ****  Rayleigh scattering (ICOL=1).
C
 2100 CONTINUE
      DEA=0.0D0
      DE=0.0D0
      CALL GRAa(E,CDT,IEFF,MAT)
C  ****  Delta interaction. Introduced to correct for the use of an
C        upper bound of the Rayleigh attenuation coefficient.
      IF(IEFF.EQ.0) THEN
        ICOL=7
        RETURN
      ENDIF
      ICOL=1
C
      IF(IPOL.EQ.1) THEN
        CALL DIRPOL(CDT,DF,0.0D0,SP1,SP2,SP3,U,V,W)
      ELSE
        DF=TWOPI*RAND(2.0D0)
        CALL DIRECT(CDT,DF,U,V,W)
      ENDIF
      ILB(1)=ILB(1)+1
      ILB(2)=KPAR
      ILB(3)=ICOL
      RETURN
C
C  ****  Compton scattering (ICOL=2).
C
 2200 ICOL=2
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      CALL GCOa(E,DE,EP,CDT,ES,CDTS,MAT,IZA,ISA)
      US=U
      VS=V
      WS=W
      DF=-1.0D0
      IF(IZA.GT.0.AND.ISA.LT.17) THEN
        ILBA(3)=ICOL
        WGHTA=WGHT
        WGHT=WGHT*WFORCE
        CALL RELAX(IZA,ISA)
        WGHT=WGHTA
      ENDIF
C  ****  Compton electron.
      IF(ES.GT.EABS(1,MAT)) THEN
        IF(DF.LT.-0.5D0) DF=TWOPI*RAND(4.0D0)
        DFS=DF+PI
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=IZA*1000000+ISA
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
      DEA=DE
      DE=DEA*WFORCE
C  ****  New direction and energy.
      IF(LCOL) THEN
        IF(EP.GT.EABS(2,MAT)) THEN
          IF(IPOL.EQ.1) THEN
            ECDT=E*RREV*(1.0D0-CDT)
            CONS=ECDT*ECDT/(1.0D0+ECDT)
            CALL DIRPOL(CDT,DF,CONS,SP1,SP2,SP3,U,V,W)
          ELSE
            DF=TWOPI*RAND(3.0D0)
            CALL DIRECT(CDT,DF,U,V,W)
          ENDIF
          E=EP
        ELSE
          DEA=DEA+EP
          DE=DE+EP
          E=0.0D0
        ENDIF
        ILB(1)=ILB(1)+1
        ILB(2)=KPAR
        ILB(3)=ICOL
      ENDIF
      RETURN
C
C  ****  Photoelectric absorption (ICOL=3).
C
 2300 ICOL=3
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      CALL GPHa(ES,IZA,ISA)
C  ****  Delta interaction. Introduced to correct for the use of an
C        upper bound of the photoelectric attenuation coefficient.
      IF(IZA.EQ.0) THEN
        ICOL=7
        DEA=0.0D0
        DE=0.0D0
        RETURN
      ENDIF
C
      IF(ES.GT.EABS(1,MAT)) THEN
        CALL SAUTER(ES,CDTS)
        DFS=TWOPI*RAND(5.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=IZA*1000000+ISA
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
      IF(ISA.LT.17) THEN
        ILBA(3)=ICOL
        WGHTA=WGHT
        WGHT=WGHT*WFORCE
        CALL RELAX(IZA,ISA)
        WGHT=WGHTA
      ENDIF
      DEA=E
      DE=E*WFORCE
      IF(LCOL) E=0.0D0
      RETURN
C
C  ****  Electron-positron pair production (ICOL=4).
C
 2400 ICOL=4
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      CALL GPPa(EE,CDTE,EP,CDTP,IZA,ISA)
      DE=E
C  ****  Electron.
      IF(EE.GT.EABS(1,MAT)) THEN
        DF=TWOPI*RAND(6.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTE,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(EE,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
C  ****  Positron.
      IF(EP.GT.EABS(3,MAT)) THEN
        DF=TWOPI*RAND(7.0D0)
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTP,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(EP,X,Y,Z,US,VS,WS,WGHT*WFORCE,3,ILBA,0)
C  ****  The positron carries a 'latent' energy of 1022 keV.
        DE=DE-TREV
      ELSE
        WGHTA=WGHT
        WGHT=WGHT*WFORCE
        CALL PANaR(EABS(2,MAT))
        WGHT=WGHTA
      ENDIF
C  ****  Atomic relaxation after triplet production.
      IF(ISA.LT.17) THEN
        ILBA(3)=ICOL
        WGHTA=WGHT
        WGHT=WGHT*WFORCE
        CALL RELAX(IZA,ISA)
        WGHT=WGHTA
      ENDIF
C
      DEA=DE
      DE=DE*WFORCE
      IF(LCOL) E=0.0D0
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 2800 ICOL=8
      DEA=0.0D0
      DE=0.0D0
      CALL GAUX
      RETURN
C
C  ************  Positrons (KPAR=3).
C
 3000 CONTINUE
      IF(MHINGE.EQ.1) GO TO 3100
C
C  ****  Hinge, artificial soft event (ICOL=1).
C
      ICOL=1
      MHINGE=1
C
C  ****  Energy loss.
C
      IF(KSOFTI.EQ.1) THEN
        DE=DESOFT
        E=E-DE
        IF(E.LT.EABS(3,MAT)) THEN
          CALL PANaR(EABS(2,MAT))  ! Annihilation at rest.
          DE=E0STEP+TREV
          DEA=DE
          E=0.0D0
          RETURN
        ENDIF
        DEA=DE
        E0STEP=E0STEP-SSOFT*(DST-DSR)
        IF(KSOFTE.EQ.0) RETURN
        XEL=LOG(E0STEP)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
      ELSE
        DEA=0.0D0
        DE=0.0D0
      ENDIF
C
C  ****  Angular deflection.
C
      IF(T1E(MAT,KE+1).GT.-78.3D0) THEN
        T1=EXP(T1P(MAT,KE)+(T1P(MAT,KE+1)-T1P(MAT,KE))*XEK)
        T2=EXP(T2P(MAT,KE)+(T2P(MAT,KE+1)-T2P(MAT,KE))*XEK)
      ELSE
        T1=0.0D0
        T2=0.0D0
      ENDIF
      IF(T1.LT.1.0D-20) RETURN
C  ****  1st and 2nd moments of the angular distribution.
      EMU1=0.5D0*(1.0D0-EXP(-DST*T1))
      EMU2=EMU1-(1.0D0-EXP(-DST*T2))/6.0D0
C  ****  Sampling from a two-bar histogram with these moments.
      PNUM=2.0D0*EMU1-3.0D0*EMU2
      PDEN=1.0D0-2.0D0*EMU1
      PMU0=PNUM/PDEN
      PA=PDEN+PMU0
      RND=RAND(2.0D0)
      IF(RND.LT.PA) THEN
        CDT=1.0D0-2.0D0*PMU0*(RND/PA)
      ELSE
        CDT=1.0D0-2.0D0*(PMU0+(1.0D0-PMU0)*((RND-PA)/(1.0D0-PA)))
      ENDIF
      DF=TWOPI*RAND(3.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      RETURN
C
C  ************  Hard event.
C
 3100 CONTINUE
      MHINGE=0
C  ****  A delta interaction (ICOL=7) occurs when the maximum
C        allowed step length is exceeded.
      IF(KDELTA.EQ.1) THEN
        ICOL=7
        DEA=0.0D0
        DE=0.0D0
        RETURN
      ENDIF
C  ****  Interaction forcing.
      DO KCOL=3,6
        IF(FORCE(IBODY,3,KCOL).GT.1.0D0.AND.POR(KCOL).GT.1.0D-16) THEN
          P0(KCOL)=POR(KCOL)
          P(KCOL)=POR(KCOL)*FORCE(IBODY,3,KCOL)
          LFORC(KCOL)=.TRUE.
        ELSE
          LFORC(KCOL)=.FALSE.
        ENDIF
      ENDDO
      IBR=MAX(IBRSPL(IBODY),1)
C  ****  Random sampling of the interaction type.
      STNOW=P(2)+P(3)+P(4)+P(5)+P(6)+P(8)
      STS=MAX(STNOW,ST)*RAND(4.0D0)
      SS=P(2)
      IF(SS.GT.STS) GO TO 3200
      SS=SS+P(3)
      IF(SS.GT.STS) GO TO 3300
      SS=SS+P(4)
      IF(SS.GT.STS) GO TO 3400
      SS=SS+P(5)
      IF(SS.GT.STS) GO TO 3500
      SS=SS+P(6)
      IF(SS.GT.STS) GO TO 3600
      SS=SS+P(8)
      IF(SS.GT.STS) GO TO 3800
C  ****  A delta interaction (ICOL=7) may occur when the total
C        interaction probability per unit path length, ST, is
C        larger than STNOW.
      ICOL=7
      DEA=0.0D0
      DE=0.0D0
      RETURN
C
C  ****  Hard elastic collision (ICOL=2).
C
 3200 ICOL=2
      IF(E.GE.PELMAX(MAT)) THEN
        TRNDC=RNDCP(MAT,KE)+(RNDCP(MAT,KE+1)-RNDCP(MAT,KE))*XEK
        TA=EXP(AP(MAT,KE)+(AP(MAT,KE+1)-AP(MAT,KE))*XEK)
        TB=BP(MAT,KE)+(BP(MAT,KE+1)-BP(MAT,KE))*XEK
        CALL EELa(TA,TB,TRNDC,RMU)
      ELSE
        TRNDC=RNDCPd(MAT,KE)+(RNDCPd(MAT,KE+1)-RNDCPd(MAT,KE))*XEK
        CALL PELd(TRNDC,RMU)  ! Uses the ELSEPA database.
      ENDIF
      CDT=1.0D0-(RMU+RMU)
      DF=TWOPI*RAND(5.0D0)
      CALL DIRECT(CDT,DF,U,V,W)
      DEA=0.0D0
      DE=0.0D0
      RETURN
C
C  ****  Hard inelastic collision (ICOL=3).
C
 3300 ICOL=3
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      DELTA=DEL(MAT,KE)+(DEL(MAT,KE+1)-DEL(MAT,KE))*XEK
      CALL PINa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IOSC)
C  ****  Scattering angles (primary positron).
      DF=TWOPI*RAND(6.0D0)
C  ****  Delta ray.
      IF(ES.GT.EABS(1,MAT)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
C  ****  New energy and direction.
      DEA=DE
      DE=DEA*WFORCE
      IF(LCOL) THEN
        IF(EP.GT.EABS(3,MAT)) THEN
          E=EP
          CALL DIRECT(CDT,DF,U,V,W)
        ELSE
          CALL PANaR(EABS(2,MAT))  ! Annihilation at rest.
          DEA=DEA+EP+TREV
          DE=DE+EP+TREV
          E=0.0D0
        ENDIF
      ENDIF
      RETURN
C
C  ****  Hard bremsstrahlung emission (ICOL=4).
C
 3400 ICOL=4
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      CALL EBRa(E,DE,MAT)
C  ****  Bremsstrahlung photon.
      IF(DE.GT.EABS(2,MAT)) THEN
        CALL EBRaA(E,DE,CDTS,MAT)
        WSPLIT=WGHT*WFORCE/IBR
        DO I=1,IBR
          DFS=TWOPI*RAND(7.0D0)
          US=U
          VS=V
          WS=W
          CALL DIRECT(CDTS,DFS,US,VS,WS)
          ILBA(1)=ILB(1)+1
          ILBA(2)=KPAR
          ILBA(3)=ICOL
          ILBA(4)=0
          ILBA(5)=ILB(5)
          CALL STORES(DE,X,Y,Z,US,VS,WS,WSPLIT,2,ILBA,0)
        ENDDO
      ENDIF
C  ****  New energy.
      DEA=DE
      DE=DE*WFORCE
      IF(LCOL) THEN
        E=E-DEA
        IF(E.LT.EABS(3,MAT)) THEN
          CALL PANaR(EABS(2,MAT))  ! Annihilation at rest.
          DEA=DEA+E+TREV
          DE=DE+E+TREV
          E=0.0D0
        ENDIF
      ENDIF
      RETURN
C
C  ****  Ionization of an inner shell (ICOL=5).
C
 3500 ICOL=5
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      DELTA=DEL(MAT,KE)+(DEL(MAT,KE+1)-DEL(MAT,KE))*XEK
      CALL PSIa(E,DELTA,DE,EP,CDT,ES,CDTS,MAT,IZA,ISA)
C  ****  Scattering angles (primary electron).
      DF=TWOPI*RAND(8.0D0)
C  ****  Delta ray.
      IF(ES.GT.EABS(1,MAT)) THEN
        DFS=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDTS,DFS,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(ES,X,Y,Z,US,VS,WS,WGHT*WFORCE,1,ILBA,0)
      ENDIF
C  ****  Atomic relaxation.
      IF(IZA.GT.2) THEN
        ILBA(3)=ICOL
        WGHTA=WGHT
        WGHT=WGHT*WFORCE
        CALL RELAX(IZA,ISA)
        WGHT=WGHTA
      ENDIF
C  ****  New energy and direction.
      DEA=DE
      DE=DEA*WFORCE
      IF(LCOL) THEN
        IF(EP.GT.EABS(3,MAT)) THEN
          E=EP
          CALL DIRECT(CDT,DF,U,V,W)
        ELSE
          CALL PANaR(EABS(2,MAT))  ! Annihilation at rest.
          DEA=DEA+EP+TREV
          DE=DE+EP+TREV
          E=0.0D0
        ENDIF
      ENDIF
      RETURN
C
C  ****  Positron annihilation in flight (ICOL=6).
C
 3600 ICOL=6
      IF(LFORC(ICOL)) THEN  ! Forced interaction.
        WFORCE=P0(ICOL)/P(ICOL)
        IF(RAND(10.0D0).LT.WFORCE) THEN
          LCOL=.TRUE.
        ELSE
          LCOL=.FALSE.
        ENDIF
      ELSE  ! Unforced interaction.
        LCOL=.TRUE.
      ENDIF
C
      CALL PANa(E,E1,CDT1,E2,CDT2,MAT)
      DF=TWOPI*RAND(9.0D0)
      IF(E1.GT.EABS(2,MAT)) THEN
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDT1,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(E1,X,Y,Z,US,VS,WS,WGHT*WFORCE,2,ILBA,0)
      ENDIF
      IF(E2.GT.EABS(2,MAT)) THEN
        DF=DF+PI
        US=U
        VS=V
        WS=W
        CALL DIRECT(CDT2,DF,US,VS,WS)
        ILBA(1)=ILB(1)+1
        ILBA(2)=KPAR
        ILBA(3)=ICOL
        ILBA(4)=0
        ILBA(5)=ILB(5)
        CALL STORES(E2,X,Y,Z,US,VS,WS,WGHT*WFORCE,2,ILBA,0)
      ENDIF
      DEA=E+TREV
      DE=DEA*WFORCE
      IF(LCOL) E=0.0D0
      RETURN
C
C  ****  Auxiliary fictitious mechanism (ICOL=8).
C
 3800 ICOL=8
      DEA=0.0D0
      DE=0.0D0
      CALL PAUX
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE JUMPW
C  *********************************************************************
      SUBROUTINE JUMPW(DS)
C
C  Modified subroutine 'JUMP' for the simulation of photons with
C  Woodcock's method (delta scattering).
C  This subroutine does not modify the weight of the running photon.
C
C  Sampling of the free path from the starting point to the position
C  of the next tentative event.
C
C  Output argument:
C    DS ...... step length (output).
C
      USE TRACK_mod
      USE PENELOPE_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LINIT,LFORC
C  ****  Composition data.
      COMMON/COMPOS/STF(MAXMAT,30),ZT(MAXMAT),AT(MAXMAT),RHO(MAXMAT),
     1  VMOL(MAXMAT),IZC(MAXMAT,30),NELEM(MAXMAT)
      COMMON/CECUTR/ECUTR(MAXMAT)
C  ****  Energy grid and interpolation constants for the current energy.
      COMMON/CEGRID/EMIN,EL,EU,ET(NEGP),DLEMP(NEGP),DLEMP1,DLFC,
     1  XEL,XE,XEK,KE
C  ****  Current state and IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DSR,W1,W2,T1,T2
      COMMON/CJUMP1/ELAST1,ELAST2,MHINGE,KSOFTE,KSOFTI,KDELTA
C  ****  Woodcock method. Largest attenuation coefficient.
      COMMON/CWOODC/STMV(NEGP),STMC(NEGP),STMAX
      DIMENSION FF(8)
C  ****  Interaction forcing parameters.
      COMMON/CFORCF/POR(8),P0(8),IBR,LFORC(8)
C
      DATA LINIT/.TRUE./
      SAVE LINIT
C
      EXTERNAL RAND
C
C  ****  Initialisation (performed only at the first call to JUMPW).
C
      IF(LINIT) THEN
        WRITE(26,1000)
 1000   FORMAT(/3X,'Woodcock''s delta-scattering method is active')
C  **** Largest forcing factors.
        DO ICOL=1,8
          FF(ICOL)=1.0D0
          DO KB=1,NBV
            FF(ICOL)=MAX(FORCE(KB,2,ICOL),FF(ICOL))
          ENDDO
        ENDDO
C  **** Largest attenuation coefficient.
        XEK=0.0D0
        DO IE=1,NEGP
          IF(IE.EQ.NEGP) THEN
            KE=NEGP-1
            XEK=1.0D0
          ELSE
            KE=IE
          ENDIF
          STMVI=0.0D0
          STMCI=0.0D0
          DO MM=1,NMAT
            MAT=MM
            CALL GIMFP
            SUMV=P(2)*FF(2)+P(4)*FF(4)+P(8)*FF(8)
            STMVI=MAX(STMVI,SUMV)
            SUMC=P(1)*FF(1)+P(3)*FF(3) ! Rayleigh, photoeffect.
            STMCI=MAX(STMCI,SUMC)
          ENDDO
          STMV(IE)=LOG(MAX(STMVI,1.0D-80))
          STMC(IE)=STMCI
        ENDDO
        DS=0.0D0
        LINIT=.FALSE.
      ENDIF
C
C  ****  Verify that the particle is a photon.
C
      IF(KPAR.NE.2) THEN
        WRITE(26,'(/3X,A,I3)') 'KPAR =',KPAR
        WRITE(26,'(/3X,A)') 'JUMPW: The particle is not a photon.'
        STOP 'JUMPW: The particle is not a photon.'
      ENDIF
C
C  ****  Length of the next free flight.
C
      IF(E.LT.ELAST1) THEN
        XEL=LOG(E)
        XE=1.0D0+(XEL-DLEMP1)*DLFC
        KE=XE
        XEK=XE-KE
        ELAST1=E
        STMAX=EXP(STMV(KE)+(STMV(KE+1)-STMV(KE))*XEK)+STMC(KE)
      ENDIF
      DS=-LOG(RAND(2.0D0))/STMAX
      RETURN
      END
C  *********************************************************************
C                       SUBROUTINE KNOCKW
C  *********************************************************************
      SUBROUTINE KNOCKW(DE,ICOL)
C
C  Modified subroutine 'KNOCK' for the simulation of photon interactions
C  with Woodcock's method (delta scattering).
C  This subroutine does not modify the weight of the running photon.
C
C  Output arguments:
C    DE ....... net energy deposited in the event, = DEA*WFORCE.
C               The weight WGHT of the running particle is not modified.
C               The effective (weighted) energy deposited in the event
C               is DE*WGHT.
C    ICOL ..... kind of interaction experienced by the particle.
C
      USE TRACK_mod
      USE PENELOPE_mod
      USE PENVARED_mod
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      LOGICAL LFORC
C  ****  Woodcock method. Largest attenuation coefficient.
      COMMON/CWOODC/STMV(NEGP),STMC(NEGP),STMAX
C  ****  Current IMFPs.
      COMMON/CJUMP0/P(8),ST,DST,DSR,W1,W2,T1,T2
      COMMON/CJUMP1/ELAST1,ELAST2,MHINGE,KSOFTE,KSOFTI,KDELTA
C  ****  Interaction forcing parameters.
      COMMON/CFORCF/POR(8),P0(8),IBR,LFORC(8)
C
      EXTERNAL RAND
C
      IF(KPAR.NE.2) THEN
        WRITE(26,'(/3X,A,I3)') 'KPAR =',KPAR
        WRITE(26,'(/3X,A)') 'KNOCKW: The particle is not a photon.'
        STOP 'KNOCKW: The particle is not a photon.'
      ENDIF
C
      CALL GIMFP
C  ****  Interaction forcing.
      P0(1)=P(1)
      DO KI=2,4
        IF(FORCE(IBODY,2,KI).GT.1.0D0.AND.P(KI).GT.1.0D-16) THEN
          P0(KI)=P(KI)
          P(KI)=P(KI)*FORCE(IBODY,2,KI)
          LFORC(KI)=.TRUE.
        ELSE
          LFORC(KI)=.FALSE.
        ENDIF
      ENDDO
C
      ST=P(1)+P(2)+P(3)+P(4)+P(8)
      IF(RAND(1.0D0)*STMAX.GT.ST) THEN  ! Delta interaction.
        DE=0.0D0
        WFORCE=1.0D0
        ICOL=7
      ELSE
        CALL KNOCKF(DE,ICOL)  ! Real interaction.
      ENDIF
C
      RETURN
      END
