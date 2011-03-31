*
* with gcc4 (for RANDOM_NUMBER) or gcc3 (simple uniform distribution)
*
* simple track model: - initial direction (COSL,0.,SINL) 
*                     - solenoidal magnectic field B=(0,0,Bz)
*                     - simplified jacobian 
*                       (parabolic in arclength)
*                     - distortions in curvilinear system
*                     - multiple scattering
*
      PROGRAM TEST 
*
      PARAMETER (NLOC=2) ! number of 'local' parameters
      PARAMETER (MP=5+NLOC)
      PARAMETER (MP2=(MP+1)*MP/2)
      DOUBLE PRECISION DPAR(MP), DCOV(MP2), DPSEED(MP2), 
     &   DJAC(5,5),DAUX(15),
     &   AJACP,AJACN,PM2L,PL2M
*
      DIMENSION CLPAR(5),DIRT(3),DIRM(3,2),TMP(5),SARC(100),
     &   DIF(2),SIG(2),PREC(2),AJACP(5,5),AJACN(5,5),
     &   DIRZ(3),DIRU(3),DIRV(3),PM2L(2,2),PL2M(2,2),BETA0(2),
     &   CLERR(5),LDER(NLOC),DERLC(2,NLOC),THSIG(2),THPREC(2)
*
      DATA DIRZ / 0.0, 0.0, 1.0 / ! Z direction
*
* (type of) local system
*
C     ITYPE=0 ! curvilinear track parameter (q/p,lambda,phi,x_t,y_t)
      ITYPE=1 ! curvilinear local system    (q/p,v',w',v,w), (v,w)=(x_t,y_t)
C     ITYPE=2 ! 'magnetic'  local system    (q/p,v',w',v,w), (v,w)=(x_t,z)
*               (perpendicular, parallel to B)
      print *, ' GBLTST ITYPE: ', ITYPE
*
      BZ=1.0 ! magnetic field 
      BFAC=BZ*0.2998
*     initial track direction (from prefit)
      SINL=0.3               ! sin(lambda)
      COSL=SQRT(1.0-SINL**2) ! cos(lambla)
      PHI=0.                 ! move in X direction
*     preset to unit matrix
      DO I=1,5
        DO J=1,5
          AJACP(I,J)=0.0D0
          AJACN(I,J)=0.0D0
        ENDDO
        AJACP(I,I)=1.0D0
        AJACN(I,I)=1.0D0
      ENDDO
      BETA0(1)=0.0
      BETA0(2)=0.0
*     distortions to track parameters
      CLERR(1)= 0.001
      CLERR(2)=-0.1
      CLERR(3)= 0.2
      CLERR(4)=-0.15
      CLERR(5)= 0.25
*
      IP=0
      DO I=1,MP
        DO J=1,I
          DPSEED(IP+J)=0.0D0
        ENDDO
        IP=IP+I
        IF (I.LE.5) DPSEED(IP)=1.0D0/DBLE(CLERR(I)**2)
      ENDDO
C      print *, ' DPSEED '
C      CALL DBPRV(6,DPSEED,MP)
*
      NTR=1000!0 ! number of tracks
      SCHI2=0.0
      DO ITR=1,NTR
*
        N=10
        S=0.0            ! at vertex
        S0=10.           ! first measurement
        STEP1=1.0/COSL   ! constant steps
        STEP2=2.0/COSL   !    in RPhi
        DERT0=1.0
*       curvilinear distortions d(q/p,theta,phi,xt,yt)
        DO K=1,5
C          CLPAR(K)=0.0
C          CLPAR(K)=CLERR(K)
          CLPAR(K)=CLERR(K)*UNRM()
        ENDDO
*
        IPNT1=0        ! first measurement, distortions defined here
*
        THETA0=1.0E-3  ! 1 mrad multiple scatt.
        IF (ITYPE.EQ.1) THEN
          THSIG(1)=THETA0         ! MS error in v' (dx_t/dz_t)
          THSIG(2)=THETA0         ! MS error in w' (dy_t/dz_t)
        ELSEIF (ITYPE.EQ.2) THEN
          THSIG(1)=THETA0/COSL    ! MS error in v' (dx_t/ds_2D)
          THSIG(2)=THETA0/COSL**2 ! MS error in w' (dz/ds_2D)
        ELSE
          THSIG(1)=THETA0         ! MS error in lambda
          THSIG(2)=THETA0/COSL    ! MS error in phi
        ENDIF
        THPREC(1)=1.0/THSIG(1)**2 ! diagonal of inverse
        THPREC(2)=1.0/THSIG(2)**2 ! covariance matrix
*
        SIG(1)=1.0E-3         ! 10 mu resolution
        SIG(2)=1.0E-3         ! 10 mu resolution
        PREC(1)=1.0/SIG(1)**2 ! diagonal of inverse
        PREC(2)=1.0/SIG(2)**2 ! covariance matrix
*
        CALL GBLINI(1)
C        CALL GBLINP(1,1)
        CALL GBLADP(IPNT)
        SARC(IPNT)=S
* create track
        S=S0
*
        DO I=1,N
C         print *, ' measurement ', I, S, OFF, SLP
* track direction 
          CPHI=COS(PHI)
          SPHI=SIN(PHI)
          DIRT(1)= COSL*CPHI
          DIRT(2)= COSL*SPHI
          DIRT(3)= SINL
* U = Z x T / |Z x T|
          DIRU(1)=-SPHI
          DIRU(2)= CPHI
          DIRU(3)= 0.
* V = T x U
          DIRV(1)=-SINL*CPHI
          DIRV(2)=-SINL*SPHI
          DIRV(3)= COSL
* V = Z
          IF (ITYPE.EQ.2) THEN
            DIRV(1)=0.
            DIRV(2)=0.
            DIRV(3)=1.
          ENDIF
* RPhi meas. ( = Y as track moves in X direction)
          EPS=0.
          IF (MOD(I,2).EQ.0) EPS=0.5
          DIRM(1,1)=EPS ! direction of measurement
          DIRM(2,1)=SQRT(1.0-EPS**2)
          DIRM(3,1)=0.0
* Z meas.
          EPS=0.
C          IF (MOD(I,2).EQ.0) EPS=0.3
          DIRM(1,2)=0.0 ! direction of measurement
          DIRM(2,2)=EPS
          DIRM(3,2)=SQRT(1.0-EPS**2)
* projection (du/dm: measurement to local)
          DO L=1,2
            PM2L(L,1)=
     &      DBLE(DIRM(1,L)*DIRU(1)+DIRM(2,L)*DIRU(2)+DIRM(3,L)*DIRU(3))
            PM2L(L,2)=
     &      DBLE(DIRM(1,L)*DIRV(1)+DIRM(2,L)*DIRV(2)+DIRM(3,L)*DIRV(3))
          ENDDO
* projection (dm/du: local to measurement)
          DET=PM2L(1,1)*PM2L(2,2)-PM2L(1,2)*PM2L(2,1)
          PL2M(1,1)= PM2L(2,2)/DET
          PL2M(1,2)=-PM2L(1,2)/DET
          PL2M(2,1)=-PM2L(2,1)/DET
          PL2M(2,2)= PM2L(1,1)/DET
* measurement - prediction in measurement system with error
          DO M=1,2
            DIF(M)=CLPAR(4)*SNGL(PL2M(M,1))
     &            +CLPAR(5)*SNGL(PL2M(M,2))+SIG(M)*UNRM()
          ENDDO
* distort with local parameters
C          DIF(1)=DIF(1)+DERT0*0.0075
C          DIF(2)=DIF(2)-DERT0*0.0025
*
          CALL GBLADP(IPNT)
          SARC(IPNT)=S
          CALL GBLADM(PL2M,DIF,PREC)
          IF (IPNT1.EQ.0) IPNT1=IPNT
* local or global parameters
          LDER(1)=4711
          LDER(2)=4712
          DERLC(1,1)= DERT0
          DERLC(2,1)=0.0
          DERLC(1,2)=0.0
          DERLC(2,2)=-DERT0
C          CALL GBLADL(2,DERLC,IRET)
C          CALL GBLADG(2,LDER,DERLC,IRET)
          DERT0=-DERT0
* propagate by STEP1
          DS=STEP1
          CALL GBLTJC(ITYPE,STEP1,BFAC,COSL,AJACN)
*
          DO K=1,5
            TMP(K)=CLPAR(K)
          ENDDO
          DO K=1,5
            CLPAR(K)=0.0
            DO L=1,5
              CLPAR(K)=CLPAR(K)+SNGL(AJACN(K,L))*TMP(L)
            ENDDO
          ENDDO
          S=S+STEP1*0.5
          CALL GBLADP(IPNT)
          SARC(IPNT)=S
          S=S+STEP1*0.5
          IF (I.LT.N) THEN
             CALL GBLADP(IPNT)
             SARC(IPNT)=S
             CALL GBLADS(BETA0,THPREC)
          ENDIF
* scatter a little
          CLPAR(2)=CLPAR(2)+THSIG(1)*UNRM()
          CLPAR(3)=CLPAR(3)+THSIG(2)*UNRM()
* propagate by STEP1
          CALL GBLTJC(ITYPE,STEP2,BFAC,COSL,AJACN)
*
          DO K=1,5
            TMP(K)=CLPAR(K)
          ENDDO
          DO K=1,5
            CLPAR(K)=0.0
            DO L=1,5
              CLPAR(K)=CLPAR(K)+SNGL(AJACN(K,L))*TMP(L)
            ENDDO
          ENDDO
          S=S+STEP2
        ENDDO
*
        NPNT=IPNT
        DO I=1,NPNT
          CALL GBLQRJ(I,IPRV,INXT)
C          print *, ' GBLQRY ', I, IPRV, INXT, SARC(I)
          IF (IPRV.GT.0) THEN
*           simple curvilinear jacobians ('parabola')
*              need only part for offsets d(xt',yt')/d(q/p,lambda,phi,xt,yt)
            CALL GBLTJC(ITYPE,SARC(IPRV)-SARC(I),BFAC,COSL,AJACP)
            CALL GBLTJC(ITYPE,SARC(INXT)-SARC(I),BFAC,COSL,AJACN)
*
            CALL GBLADJ(I,AJACP,AJACN)
          ENDIF
        ENDDO
*
C        CALL GBLADX(IPNT1,DPSEED)
*
C        CALL GBLDMP
        CALL GBLFIT('',MRANK,NP,NDF,CHI2,WLS)
C        print *, ' GBLFIT ', MRANK, NP, NDF, CHI2, WLS
        SCHI2=SCHI2+CHI2/FLOAT(NDF)
C        CALL GBLFIT('HT',MRANK,NDF,CHI2,WLS)
        IF (ITR.LE.1) THEN
           print *, ' GBLFIT ', MRANK, NP, NDF, CHI2, WLS
           CALL GBLRES(IPNT1,DPAR,DCOV)
           print *, ' DPAR ', (DPAR(K),K=1,NP)
           print *, ' DCOV(loc) '
           CALL DBPRV(6,DCOV,NP)
           CALL GBLL2C(ITYPE,COSL,DJAC)
           CALL DBAVAT(DCOV,DJAC,DAUX,5,5)
           print *, ' DCOV(cvl) '
           CALL DBPRV(6,DAUX,5)
        ENDIF
*
C        CALL GBLMP2(IRET) ! prepare input for Millepede-II
*
      ENDDO
*
      SCHI2=SCHI2/FLOAT(NTR)
      print *, ' === <Chi2/ndf> === ', NTR, SCHI2
*
      END
*
*     simplified jacobian (quadratic in DS)
*
      SUBROUTINE GBLTJC(ITYPE,DS,BFAC,COSL,AJAC)
      DOUBLE PRECISION AJAC(5,5)
*
      IF (ITYPE.EQ.1) THEN
*       curvilinear local system (q/p,v',w',v,w), (v,w)=(x_t,y_t)
        AJAC(1,1)= 1.0D0
        AJAC(2,1)=-DBLE(BFAC*DS*COSL)
        AJAC(2,2)= 1.0D0
        AJAC(3,3)= 1.0D0
        AJAC(4,1)=-DBLE(0.5*BFAC*DS*DS*COSL)
        AJAC(4,2)= DBLE(DS)
        AJAC(4,4)= 1.0D0
        AJAC(5,3)= DBLE(DS)
        AJAC(5,5)= 1.0D0
      ELSEIF (ITYPE.EQ.2) THEN
*       local system (q/p,v',w',v,w), (v,w)=(x_t,z)
        AJAC(1,1)= 1.0D0
        AJAC(2,1)=-DBLE(BFAC*DS)
        AJAC(2,2)= 1.0D0
        AJAC(3,3)= 1.0D0
        AJAC(4,1)=-DBLE(0.5*BFAC*DS*DS*COSL)
        AJAC(4,2)= DBLE(DS*COSL)
        AJAC(4,4)= 1.0D0
        AJAC(5,3)= DBLE(DS*COSL)
        AJAC(5,5)= 1.0D0
      ELSE
*       curvilinear track parameter (q/p,lambda,phi,x_t,y_t)
        AJAC(1,1)= 1.0D0
        AJAC(2,2)= 1.0D0
        AJAC(3,1)=-DBLE(BFAC*DS)
        AJAC(3,3)= 1.0D0
        AJAC(4,1)=-DBLE(0.5*BFAC*DS*DS*COSL)
        AJAC(4,3)= DBLE(DS*COSL)
        AJAC(4,4)= 1.0D0
        AJAC(5,2)= DBLE(DS)
        AJAC(5,5)= 1.0D0
      ENDIF
*
      RETURN
      END
*
*     simplified transformation local to curvilinear
*
      SUBROUTINE GBLL2C(ITYPE,COSL,DJAC)
      DOUBLE PRECISION DJAC(5,5)
*     preset with unit matrix
      DO K=1,5
        DO L=1,5
          DJAC(K,L)=0.0D0
        ENDDO
        DJAC(K,K)=1.0D0
      ENDDO
*
      IF (ITYPE.EQ.1) THEN
*       curvilinear local system (q/p,v',w',v,w), (v,w)=(x_t,y_t)
        DJAC(2,2)=0.0D0
        DJAC(2,3)=1.0D0/DBLE(COSL)
        DJAC(3,2)=1.0D0
        DJAC(3,3)=0.0D0
      ELSEIF (ITYPE.EQ.2) THEN
*       local system (q/p,v',w',v,w), (v,w)=(x_t,z)
        DJAC(2,2)=0.0D0
        DJAC(2,3)=1.0D0
        DJAC(3,2)=DBLE(COSL)**2
        DJAC(3,3)=0.0D0
        DJAC(5,5)=DBLE(COSL)
      ENDIF
*
      RETURN
      END
*
      FUNCTION UNRM() ! simple unit normal
      DIMENSION R(12)
C      DATA IX / 4711 /      ! use for gcc3
      CALL RANDOM_NUMBER(R) ! comment out for gcc3
      UNRM=-6.0
      DO K=1,12
C      R(K)=UNIF(IX)         ! use for gcc3
      UNRM=UNRM+R(K)
      ENDDO
      END

      FUNCTION UNIF(IX)
* portable random number generator using the recursion:
* IX = 16807 * IX MOD (2**31-1)
      K1 = IX/127773
      IX = 16807* (IX - K1*127773) - K1 * 2836
      IF (IX.LT.0) IX = IX + 2147483647
      UNIF = IX*4.656612875E-10
      RETURN
      END

