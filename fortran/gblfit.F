*
* C. Kleinwort,  DESY-FH1, www.terascale.de
* March 2011, Analysis Centre: Statistics Tools Group
*
      SUBROUTINE GBLFIT(CDW,MRANK,NP,NDF,CHI2,WLS)
*
* Input
*       CDW   : string defining iterations for outlier down weighting,
*               one char per iteration (C: Cauchy, H: Huber, T: Tukey)
* Output
*       MRANK : rank of measurements
*       NP    : number of track parameter at given point
*       NDF   : degrees of freedom
*       CHI2  : Chi2
*       WLS   : lost measurements: N-sum(weight)
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      CHARACTER*(*) CDW
*
      MRANK=0
      NP=0
      NDF=0
      CHI2=0.
      WLS=0.
      LCDW=LEN_TRIM(CDW)
      ITER=0
*
      CALL GBLPRP
      IF (NLVL.LT.3) RETURN
*
  100 CONTINUE
*
      CALL GBLMAT
      MRANK=NCRV+2*NDIM ! rank of trajectory
      NP   =5+NLOC      ! number of track parameters
*
      IF (LPR.GT.1) 
     &   PRINT *, ' GBLFIT: NPAR, NBND, NBDR ', NPAR, NBND, NCRV+NLOC
      INV=1 ! calc band part of covariance matrix
      CALL SQMIBB(CMAT,BVEC,NPAR,NCRV+NLOC,NBND,INV,NRANK)
      IF (NRANK.NE.NPAR) THEN ! rank deficit
         IF (LPR.GT.1) PRINT *, ' GBLFIT: rank deficit ', NRANK, NPAR
         MRANK=-MRANK
         CHI2=0.
         NDF=0.
         ALS=0.
         RETURN 
      ENDIF
*
      ITER=ITER+1
      IDWM=0 ! down weighting method (Tukey, Huber, Cauchy)
      IF (ITER.LE.LCDW) IDWM=(INDEX('TtHhCc',CDW(ITER:ITER))+1)/2
*
      CALL GBLCH2(IDWM,CHI2,SWGT)
      NDF=NDAT-NPAR+NDFSD
* iterate ?
      IF (IDWM.GT.0) GOTO 100
*
      WLS=FLOAT(NDAT)-SWGT
      NLVL=4 ! fitted
*
      RETURN
      END
*
      SUBROUTINE GBLPRP
*
* Prepare broken lines trajectory
*
#include "gbldef.inc"
#include "gblcom.inc"
* index in symmetric matrix     
      IJSYM(I,J)=MIN(I,J)+(MAX(I,J)*MAX(I,J)-MAX(I,J))/2
*
      IF (NLVL.NE.2) RETURN     ! no points with jacobians added
      IF (NOFF.LT.2) RETURN     ! too few parameters
* number of parameters
      NLOC=MIN(NLOC,MXBDR-NCRV) ! limit number of (external) local par.
      NPAR=NOFF*NDIM+NCRV+NLOC  ! number of parameters
      IF (NPAR.GT.MXPAR) RETURN ! too many parameters
*
      CALL GBLDAT
*
      NBNDSD=0
      NDFSD=0
      IF (IXSD.NE.0) THEN     ! contribution of seed to linear equation system
*
         CALL GBLJAC(IXSD,1)  ! get transposed jacobian broken lines -> local
         CALL DBAVAT(DPSD,DJAC,CMATSD,5+NLOC,NPAR)
*
         DO I=1,NPAR
            DO J=1,I
               IJ=IJSYM(I,J)
               IF (CMATSD(IJ).NE.0.0D0.AND.MIN(I,J).GT.NCRV+NLOC) 
     &            NBNDSD=MAX(NBNDSD,IABS(I-J)) ! update band width
            ENDDO
            IF (CMATSD(IJ).GT.0.0D0) NDFSD=NDFSD+1
         ENDDO
*
      ELSE                    ! no seed
         NPAR2=(NPAR*NPAR+NPAR)/2
         DO I=1,NPAR2
            CMATSD(I)=0.0D0
         ENDDO
      ENDIF
*
      NLVL=3 ! prepared
*
      RETURN
      END
*
* initialize
*
      SUBROUTINE GBLINI(LPRNT)
*
* Input
*       LPRNT  : print level
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DATA IFIRST / 1 /
*
      NDIM=2  ! 2D offsets for track in 3D
      ICRD0=0
      GOTO 10
*
      ENTRY GBLINP(LPRNT,ICOORD)
*
* Input
*       LPRNT    : print level
*       ICOORD   : coordinate (1: u_1, 2: u_2) to use
*
      NDIM=1  ! 1D offsets for track in 2D
      ICRD0=0
      IF (ICOORD.GT.1) ICRD0=1
*
   10 CONTINUE
*
      NLVL=0
      LPR=LPRNT
      NCRV=1  ! with Q/P
      NLOC=0
      NPNT=0
      LPNT=0
      NMEAS=0
      NSCAT=0
      LSCAT=0
      NOFF=0
      NKNK=0
      MDAT=MXDAT
      IXSD=0  ! external seed parameter offset
*
      IF (IFIRST.EQ.1) THEN
        IFIRST=0
        IF (LPR.GT.0) PRINT *, ' GENBRL $Rev$'
      ENDIF
*
      RETURN
      END
*
* add point to trajectory
*
      SUBROUTINE GBLADP(IPOINT)
*
* Output
*       IPOINT : identifier 
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      IPOINT=0
      IF (NPNT.GE.MXPNT) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADP: too many points, ignoring new point '
         RETURN
      ENDIF
      NPNT=NPNT+1
      IMEAS(NPNT)=0
      ISCAT(NPNT)=0
      IOFF(NPNT)=0
      IKNK(NPNT)=0
      IDLC(1,NPNT)=0
      IDLC(2,NPNT)=0
      IDGL(1,NPNT)=0
      IDGL(2,NPNT)=0
      IPOINT=NPNT
      DO K=1,8
        JMAT(K,NPNT)=0.0D0
        SMAT(K,NPNT)=0.0D0
      ENDDO
      DO K=1,4
        DVEC(K,NPNT)=0.0D0
      ENDDO
*
      IF (LPNT.GT.0) THEN
         IF (ISCAT(LPNT).LE.0.AND.NOFF.GT.1) THEN
            NOFF=NOFF-1
            IOFF(LPNT)=-IOFF(LPNT)
         ENDIF
      ENDIF
      NOFF=NOFF+1
      IOFF(NPNT)=NOFF
      IO2P(NOFF)=NPNT
      LPNT=NPNT
      IF (NOFF.GT.2.AND.LSCAT.GT.0) THEN
         NKNK=NKNK+1
         IKNK(LSCAT)=NKNK
         LSCAT=0
      ENDIF
*
      NLVL=1 ! adding points
      RETURN
      END
*
* add measurement to current point
*
      SUBROUTINE GBLADM(PROJ,RES,PREC)
*
* Input
*       PROJ(2,2)  : projection matrix of measurement directions into local system (dm/du)
*       RES(2)     : residuals (m)
*       PREC(2)    : diagonal of inverse covariance matrix
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION PROJ
      DIMENSION PROJ(2,2),RES(2),PREC(2)
*
      IF (NPNT.LE.0) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADM: no point defined'
         RETURN
      ENDIF
      IF (NMEAS+NSCAT.GE.MXPNT) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADM: too many measurement+scatterers ', NMEAS+NSCAT
         RETURN
      ENDIF
      IF (IMEAS(NPNT).LE.0) THEN 
         NMEAS=NMEAS+1
         IMEAS(NPNT)=NMEAS
         PMAT(1,NMEAS)=PROJ(1,1)
         PMAT(2,NMEAS)=PROJ(1,2)
         PMAT(3,NMEAS)=PROJ(2,1)
         PMAT(4,NMEAS)=PROJ(2,2)
         YVEC(1,NMEAS)=RES(1)
         YVEC(2,NMEAS)=RES(2)
         WVEC(1,NMEAS)=PREC(1)
         WVEC(2,NMEAS)=PREC(2)
      ELSE
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADM: measurement already defined for point ', NPNT
      ENDIF
*
      RETURN
      END
*
* add local derivatives to measurement
*
      SUBROUTINE GBLADL(NDER,DER,IRET)
*
* Input
*       NDER     : number of local derivatives
*       DER(2,*) : local derivatives
* Output
*       IRET     : number of non zero derivatives added
*
      DIMENSION DER(2,NDER)
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      IRET=0
      IF (NPNT.LE.0) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADL: no point defined'
         RETURN
      ENDIF
*
      DO IM=1,2
*
         IF (IDLC(IM,NPNT).EQ.0) THEN      ! local derivs not yet defined
*
            IF (LDAT.GT.MDAT-NDER) RETURN  ! enough space left
            MDAT0=MDAT
            DO K=NDER,1,-1
               IF (DER(IM,K).NE.0.0) THEN
                  NLOC=MAX(NLOC,K)
                  IDAT(MDAT)=K
                  ADAT(MDAT)=DER(IM,K)
                  MDAT=MDAT-1
               ENDIF
            ENDDO
*
            IDLC(IM,NPNT)=MDAT
            IDAT(MDAT)=MDAT0-MDAT
            ADAT(MDAT)=0.0
            MDAT=MDAT-1
            IRET=IRET+MDAT0-MDAT
*
         ELSE
            IF (LPR.GT.0) PRINT *, 
     &        ' GBLADL: local derivatives already defined for point ',
     &        NPNT
         ENDIF
*
      ENDDO
*
      RETURN
      END
*
* add global derivatives to measurement
*
      SUBROUTINE GBLADG(NDER,LDER,DER,IRET)
*
* Input
*       NDER     : number of derivatives
*       LDER(*)  : labels for global derivatives
*       DER(2,*) : global derivatives
* Output
*       IRET     : number of non zero derivatives added
*
      DIMENSION LDER(NDER),DER(2,NDER)
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      IRET=0
      IF (NPNT.LE.0) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADG: no point defined'
         RETURN
      ENDIF
*
      DO IM=1,2
*
         IF (IDGL(IM,NPNT).EQ.0) THEN      ! local derivs not yet defined
*
            IF (LDAT.GT.MDAT-NDER) RETURN  ! enough space left
            MDAT0=MDAT
            DO K=NDER,1,-1
               IF (DER(IM,K).NE.0.0) THEN
                  IDAT(MDAT)=LDER(K)
                  ADAT(MDAT)=DER(IM,K)
                  MDAT=MDAT-1
               ENDIF
            ENDDO
*
            IDGL(IM,NPNT)=MDAT
            IDAT(MDAT)=MDAT0-MDAT
            ADAT(MDAT)=0.0
            MDAT=MDAT-1
            IRET=IRET+MDAT0-MDAT
*
         ELSE
            IF (LPR.GT.0) PRINT *, 
     &        ' GBLADG: global derivatives already defined for point ',
     &        NPNT
         ENDIF
*
      ENDDO
*
      RETURN
      END
*
* add (thin) scatterer to current point
*
      SUBROUTINE GBLADS(RES,PREC)
*
* Input
*       RES(2)     : values for initial kinks (in case of iterating)
*       PREC(2)    : diagonal of inverse (multiple scattering) covariance matrix
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DIMENSION RES(2),PREC(2)
*
      IF (NPNT.LE.0) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADS: no point defined'
         RETURN
      ENDIF
      IF (NMEAS+NSCAT.GE.MXPNT) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADS: too many measurement+scatterers ', NMEAS+NSCAT
         RETURN
      ENDIF
      IF (PREC(1).LE.0.0.OR.PREC(2).LE.0.0) THEN
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADS: invalid scattering precision ', PREC
         RETURN
      ENDIF
*
      IF (ISCAT(NPNT).LE.0) THEN
         JSCAT=MXPNT-NSCAT
         NSCAT=NSCAT+1
         ISCAT(NPNT)=JSCAT
         LSCAT=NPNT
         YVEC(1,JSCAT)=RES(1)
         YVEC(2,JSCAT)=RES(2)
         WVEC(1,JSCAT)=PREC(1)
         WVEC(2,JSCAT)=PREC(2)
      ELSE
         IF (LPR.GT.0) PRINT *, 
     &     ' GBLADM: scatterer already defined for point ', NPNT
      ENDIF
*
      RETURN
      END
*
* dump trajectory definition
*
      SUBROUTINE GBLDMP
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      IF (NPNT.LE.0) THEN
         PRINT *, ' GBLDMP no trajectory defined '
         RETURN
      ENDIF
*
      PRINT *
      PRINT *, '    GBLDMP trajectory definition '
      PRINT *, '-------------------------------------'
      PRINT *, ' number of points       ', NPNT
      PRINT *, ' number of offsets      ', NOFF
      PRINT *, ' number of scatterers   ', NSCAT
      PRINT *, ' number of measurements ', NMEAS
      PRINT *, ' number of kinks        ', NKNK
      IF (NLOC.GT.0)
     &PRINT *, ' number of local par.   ', NLOC
      PRINT *
*
      DO I=1,NPNT
         PRINT *, ' Point          ', I
         IF (IOFF(I).GT.0)  PRINT *, '    Offset      ', IOFF(I)
         IF (ISCAT(I).GT.0) PRINT *, '    Scatterer   ', ISCAT(I)
         IF (IMEAS(I).GT.0) THEN
            IF (IOFF(I).GT.0) THEN
               PRINT *, '    Measurement ', IMEAS(I),
     &            ' using offset  ',IOFF(I)
            ELSE
               PRINT *, '    Measurement ', IMEAS(I),
     &            ' using offsets ',-IOFF(I)-1,' ,',-IOFF(I)
            ENDIF
         ENDIF
         IF (IKNK(I).GT.0)  PRINT *, '    Kink        ', IKNK(I),
     &       ' using offsets ',IOFF(I)-1,'..',IOFF(I)+1
         IF (IMEAS(I).LE.0.AND.ISCAT(I).LE.0) THEN
            IF (IOFF(I).GT.0) THEN
               PRINT *, '    Prediction               ',
     &            ' using offset  ',IOFF(I)

            ELSE
               PRINT *, '    Prediction               ',
     &            ' using offsets ',-IOFF(I)-1,' ,',-IOFF(I)
            ENDIF
         ENDIF
      ENDDO
*
      RETURN
      END
*
* query for jacobians
*
      SUBROUTINE GBLQRJ(IPOINT,IPREV,INEXT)
*
* Input
*       IPOINT  : point to query for jacobians
* Output
*       IPREV   : previous point, need jacobian IPOINT -> IPREV
*       INEXT   : next     point, need jacobian IPOINT -> INEXT
*            -1 : illegal point
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      IPREV=-1
      INEXT=-1
      IF (IPOINT.LT.1.OR.IPOINT.GT.NPNT) THEN
        IF (LPR.GT.0) PRINT *, ' GBLQRJ invalid point ', IPOINT
        RETURN
      ENDIF
      JOFF=IOFF(IPOINT)
      IF (JOFF.GT.0) THEN
         JOFFP=MAX(JOFF-1,1   )
         JOFFN=MIN(JOFF+1,NOFF)
         IPREV=IO2P(JOFFP)
         INEXT=IO2P(JOFFN)
      ELSE
         IPREV=IO2P(-JOFF-1)
         INEXT=IO2P(-JOFF)
      ENDIF
*
      RETURN
      END
*
* add jacobians
*
      SUBROUTINE GBLADJ(IPOINT,AJACP,AJACN)
*
* Input
*       IPOINT     : point to add jacobians
*       AJACP(5,5) : jacobian to previous point
*       AJACN(5,5) : jacobian to next point
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION AJACP, AJACN
      DIMENSION AJACP(5,5),AJACN(5,5)
*
      IF (IPOINT.LT.1.OR.IPOINT.GT.NPNT) THEN
         IF (LPR.GT.0) PRINT *, ' GBLADJ invalid point ', IPOINT
         RETURN
      ENDIF
*     IPOINT -> IPREV
                                         ! J-
      JMAT(1,IPOINT)=AJACP(4,4)
      JMAT(2,IPOINT)=AJACP(4,5)
      JMAT(3,IPOINT)=AJACP(5,4)
      JMAT(4,IPOINT)=AJACP(5,5)
                                         ! S-
      SMAT(1,IPOINT)=AJACP(4,2)
      SMAT(2,IPOINT)=AJACP(4,3)
      SMAT(3,IPOINT)=AJACP(5,2)
      SMAT(4,IPOINT)=AJACP(5,3)
                                         ! d-
      DVEC(1,IPOINT)=AJACP(4,1)
      DVEC(2,IPOINT)=AJACP(5,1)
*     IPOINT -> INEXT
                                         ! J+
      JMAT(5,IPOINT)=AJACN(4,4)
      JMAT(6,IPOINT)=AJACN(4,5)
      JMAT(7,IPOINT)=AJACN(5,4)
      JMAT(8,IPOINT)=AJACN(5,5)
                                         ! S+
      SMAT(5,IPOINT)=AJACN(4,2)
      SMAT(6,IPOINT)=AJACN(4,3)
      SMAT(7,IPOINT)=AJACN(5,2)
      SMAT(8,IPOINT)=AJACN(5,3)
                                         ! d+
      DVEC(3,IPOINT)=AJACN(4,1)
      DVEC(4,IPOINT)=AJACN(5,1)
*
      NLVL=2 ! adding jacobians
      RETURN
      END
*
      SUBROUTINE GBLDAT
*
*     generate DATa record
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION TC2L(2,2),
     &    WP(2,2),WN(2,2),WJP(2,2),WJN(2,2),
     &    DWP(2),DWN(2),WJS(2,2),WJI(2,2),
     &    WIP(2,2),WIN(2,2),DIP(2),DIN(2),
     &    WPP(2,2),WPN(2,2),DPP(2),DPN(2)
*
      L0=ICRD0
*
      LDWM=0 ! start without downweighting
* loop over points
      NDAT=0 ! number of data items
      LDAT=0 ! length of data
      IPMN=NPAR+1 ! minimum parameter with non zero derivatives
      DO JPNT=1,NPNT
*
         JOFF=IOFF(JPNT)               ! offset at point
         IPAR0=NDIM*(JOFF-1)+NCRV+NLOC ! offset in parameter number
* kink at point ?
         IF (IKNK(JPNT).GT.0) THEN
            JSCAT=ISCAT(JPNT)
*
            CALL GBLDER(JPNT,-1,WP,WJP,DWP)
            CALL GBLDER(JPNT, 2,WN,WJN,DWN)
            DO J=1,2
               DO K=1,2
                  WJS(J,K)=WJP(J,K)+WJN(J,K) !  W- * J- + W+ * J+
               ENDDO
            ENDDO
*
            DO M=1,NDIM ! pseudo measurements
               KDAT=MDAT-(NCRV+3*NDIM)
               IF (LDAT+3.GT.KDAT) GOTO 9999
               KDAT0=KDAT
*
               LDAT0=LDAT
               NDAT=NDAT +1
               IDAT(LDAT+1)=3             ! length of item
               IDAT(LDAT+2)=JPNT          ! point
               IDAT(LDAT+3)=-M            ! pseudo measurement
               ADAT(LDAT+1)=YVEC(M,JSCAT) ! (start) value
               ADAT(LDAT+2)=WVEC(M,JSCAT) ! precision (1/error**2)
               ADAT(LDAT+3)=1.0           ! (additional) weight
               LDAT=LDAT+3
*
               IDAT(KDAT+1)= 1
               ADAT(KDAT+1)= SNGL(-DWP(M)-DWN(M))
               KDAT=KDAT+NCRV
               DO L=1,NDIM
                  IDAT(KDAT+1)= IPAR0+L-NDIM
                  ADAT(KDAT+1)= SNGL(  WP(M,L+L0))
                  IDAT(KDAT+2)= IPAR0+L
                  ADAT(KDAT+2)= SNGL(-WJS(M,L+L0))
                  IDAT(KDAT+3)= IPAR0+L+NDIM
                  ADAT(KDAT+3)= SNGL(  WN(M,L+L0))
                  KDAT=KDAT+3
               ENDDO
*
               DO K=KDAT0+1,KDAT
                  IF (ADAT(K).NE.0.0) THEN ! copy non zero derivatives
                     LDAT=LDAT+1
                     IDAT(LDAT)=IDAT(K)
                     ADAT(LDAT)=ADAT(K)
                     IPMN=MIN(IPMN,IDAT(LDAT))
                  ENDIF
               ENDDO
               IDAT(LDAT0+1)=LDAT-LDAT0
            ENDDO
*
         ENDIF
* measurement at point ?
         JMEAS=IMEAS(JPNT)        ! measurement at point
         IF (JMEAS.GT.0) THEN
            IP=0
            DO J=1,2
               DO K=1,2
                  IP=IP+1
                  TC2L(J,K)=PMAT(IP,JMEAS)      !  P
               ENDDO
            ENDDO
*
            IPAR0=NDIM*(JOFF-1)+NCRV+NLOC       ! offset in parameter number
            IF (JOFF.LT.0) THEN                 ! need interpolation
               IPAR0=NDIM*(-JOFF-2)+NCRV+NLOC   ! offset in parameter number
               CALL GBLDER(JPNT,-1,WP,WJP,DWP)
               CALL GBLDER(JPNT, 2,WN,WJN,DWN)
               DO J=1,2
                  DO K=1,2
                     WJS(J,K)=WJP(J,K)+WJN(J,K) !  W- * J- + W+ * J+
                  ENDDO
               ENDDO
               CALL GBLINV(2,WJS,WJI)           ! (W- * J- + W+ * J+)^-1 (=N)
               CALL GBLMLT(2,2,WJI,WP,WIP)      !  N * W-
               CALL GBLMLT(2,2,WJI,WN,WIN)      !  N * W+
               CALL GBLMLT(2,1,WJI,DWP,DIP)     !  N * W- * d-
               CALL GBLMLT(2,1,WJI,DWN,DIN)     !  N * W+ * d+
               CALL GBLMLT(2,2,TC2L,WIP,WPP)    !  P * N * W-
               CALL GBLMLT(2,2,TC2L,WIN,WPN)    !  P * N * W+
               CALL GBLMLT(2,1,TC2L,DIP,DPP)    !  P * N * W- * d-
               CALL GBLMLT(2,1,TC2L,DIN,DPN)    !  P * N * W+ * d+
            ENDIF
*
            DO M=1,2
               IF (WVEC(M,JMEAS).LE.0.0) GOTO 100 ! no precision ?
               KDAT=MDAT-(NCRV+2*NDIM)
               IF (LDAT+3.GT.KDAT) GOTO 9999
               KDAT0=KDAT
*
               LDAT0=LDAT
               NDAT=NDAT+1
               IDAT(LDAT+1)=3               ! length of item
               IDAT(LDAT+2)=JPNT            ! point
               IDAT(LDAT+3)=M               ! measurement
               ADAT(LDAT+1)=YVEC(M,JMEAS)   ! value
               ADAT(LDAT+2)=WVEC(M,JMEAS)   ! precision (1/error**2)
               ADAT(LDAT+3)=1.0             ! (additional) weight
               LDAT=LDAT+3
*
               IF (JOFF.GT.0) THEN         ! measurement at offset
                  DO L=1,NDIM
                     IDAT(KDAT+1)= IPAR0+L
                     ADAT(KDAT+1)= SNGL(TC2L(M,L+L0))
                     KDAT=KDAT+1
                  ENDDO
               ELSE                         ! interpolation between offsets
*
                  IDAT(KDAT+1)= 1
                  ADAT(KDAT+1)= SNGL(-DPP(M)-DPN(M))
                  KDAT=KDAT+NCRV
                  DO L=1,NDIM
                     IDAT(KDAT+1)= IPAR0+L
                     IDAT(KDAT+2)= IPAR0+L+NDIM
                     ADAT(KDAT+1)= SNGL(WPP(M,L+L0))
                     ADAT(KDAT+2)= SNGL(WPN(M,L+L0))
                     KDAT=KDAT+2
                  ENDDO
               ENDIF
               DO K=KDAT0+1,KDAT
                  IF (ADAT(K).NE.0.0) THEN ! copy non zero derivatives
                     LDAT=LDAT+1
                     IDAT(LDAT)=IDAT(K)
                     ADAT(LDAT)=ADAT(K)
                     IPMN=MIN(IPMN,IDAT(LDAT))
                  ENDIF
               ENDDO
               IDAT(LDAT0+1)=LDAT-LDAT0
* check for local derivatives
               ILCL=IDLC(M,JPNT)
               IF (ILCL.LE.0) GOTO 99
               NLCL=IDAT(ILCL)
               IF (LDAT+NLCL.GT.MDAT) GOTO 99
               DO K=1,NLCL
                  IF (IDAT(ILCL+K).LE.NLOC) THEN
                     LDAT=LDAT+1
                     IDAT(LDAT)=IDAT(ILCL+K)+NCRV
                     ADAT(LDAT)=ADAT(ILCL+K)
                     IPMN=MIN(IPMN,IDAT(LDAT))
                  ENDIF
               ENDDO
               IDAT(LDAT0+1)=LDAT-LDAT0
   99          CONTINUE
*
  100          CONTINUE
            ENDDO
*
         ENDIF
*
      ENDDO
*
      IF (IXSD.NE.0) THEN
         IF (DPSD(1).GT.0.0D0) IPMN=1 ! external 'seed' for curvature
      ENDIF
*     check minimum parameter with non zero derivatives
      IF (IPMN.GT.NCRV) THEN          ! curvature undefined
         IPOFF=NCRV
         NCRV=0
*        correct parameter indices 
         IF (IPOFF.GT.0) THEN
            NPAR=NPAR-IPOFF
            LDAT=0
            DO I=1,NDAT
               LTEM=IDAT(LDAT+1)
               DO J=LDAT+4,LDAT+LTEM 
                  IDAT(J)=IDAT(J)-IPOFF
               ENDDO
               LDAT=LDAT+LTEM
            ENDDO
*
         ENDIF
      ENDIF
*
 9999 CONTINUE
*
      RETURN
      END
*
      SUBROUTINE GBLDER(IPOINT,IDIR,W,WJ,DW)
*
*     get (matrices and vectors for) derivatives
*          W ((+/-)S^-1), W*J, W*d
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION JM(2,2),SM(2,2),DV(2),W(2,2),WJ(2,2),DW(2),SGN
* 
* (parts of) jacobians to previous/next offset
      SGN=1.0D0
      IF (IDIR.LT.0) SGN=-1.0D0
      JP=(IABS(IDIR)-1)*2
      IP=JP*2
      DO J=1,2
         DO K=1,2
            IP=IP+1
            JM(J,K)=JMAT(IP  ,IPOINT)     ! J
            SM(J,K)=SMAT(IP  ,IPOINT)*SGN ! S
         ENDDO
         JP=JP+1
         DV(J)=DVEC(JP,IPOINT)            ! d
      ENDDO
*
      CALL GBLINV(2,SM,W)                 ! W
      CALL GBLMLT(2,2,W,JM,WJ)            ! W*J
      CALL GBLMLT(2,1,W,DV,DW)            ! W*d
*
      RETURN
      END
*
      SUBROUTINE GBLMAT
*
*     build MATrix and rhs vector 
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION DVAL,DWGH
* index in symmetric matrix     
      IJSYM(I,J)=MIN(I,J)+(MAX(I,J)*MAX(I,J)-MAX(I,J))/2
*
      DO I=1,NPAR
        BVEC(I)=0.0D0
      ENDDO
      NPAR2=(NPAR*NPAR+NPAR)/2
      DO I=1,NPAR2
        CMAT(I)=CMATSD(I)
      ENDDO
      NBND=NBNDSD ! band size
*
      LDAT=0
      NS=0
      NZ=0
      NT=0
      DO I=1,NDAT
        LTEM=IDAT(LDAT+1)
        DVAL=DBLE(ADAT(LDAT+1))              ! value
        DWGH=DBLE(ADAT(LDAT+3)*ADAT(LDAT+2)) ! (total) weight
        DO J=LDAT+4,LDAT+LTEM ! update matrix
          IJ=IDAT(J)
          BVEC(IJ)=BVEC(IJ)+DBLE(ADAT(J))*DVAL*DWGH
          NT=NT+1
          IF (ADAT(J).NE.0.0) NZ=NZ+1
          DO K=LDAT+4,J
            IK=IDAT(K)
            JK=IJSYM(IJ,IK)
            CMAT(JK)=CMAT(JK)+DBLE(ADAT(J))*DBLE(ADAT(K))*DWGH
            IF (MIN(IJ,IK).GT.NCRV+NLOC) NBND=MAX(NBND,IABS(IJ-IK)) ! update band width
            NS=NS+1
          ENDDO
        ENDDO
        LDAT=LDAT+LTEM
      ENDDO
*
      RETURN
      END
*
      SUBROUTINE GBLCH2(IDWM,CHI2,SWGT)
*
*     calculate Chi2
*
      DOUBLE PRECISION DCHI2
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DIMENSION DWINT(0:3) ! Integral(weight*normal_distribution)
      DATA DWINT / 1.0, 0.8737, 0.9326, 0.8228 /
*
      CHI2=0.0
      SWGT=0.0
*
      LDAT=0
      DO I=1,NDAT
        LTEM=IDAT(LDAT+1)
        VAL=ADAT(LDAT+1)      ! value
        WGH=ADAT(LDAT+3)      ! down weighting
        BRL=0.0               ! predction from broken lines
        DO J=LDAT+4,LDAT+LTEM 
          IJ=IDAT(J)
          BRL=BRL+BVEC(IJ)*ADAT(J)
        ENDDO
        SRES=ABS(BRL-VAL)*SQRT(ADAT(LDAT+2))
        CHI2=CHI2+SRES*SRES*WGH
        SWGT=SWGT+WGH
*       down weighting
        IF (IDWM.EQ.1) THEN
          IF(SRES.LT.4.6851) THEN          ! Tukey 
            WGH=(1.0-0.045558*SRES*SRES)**2  ! 1/4.6851**2
          ELSE
            WGH=0.0
          ENDIF
        ELSEIF (IDWM.EQ.2) THEN
          IF (SRES.LT.1.345) THEN          ! Huber
            WGH=1.0
          ELSE
            WGH=1.345/SRES
          ENDIF
        ELSEIF (IDWM.EQ.3) THEN
          WGH=1.0/(1.0+(SRES/2.3849)**2)   ! Cauchy
        ENDIF
        ADAT(LDAT+3)=WGH
        LDAT=LDAT+LTEM
      ENDDO
*
      CHI2=CHI2/DWINT(LDWM) ! renormalize Chi2
      LDWM=IDWM
*
      IF (IXSD.NE.0) THEN
         CALL DBAVAT(CMATSD,BVEC,DCHI2,NPAR,1)
         CHI2=CHI2+SNGL(DCHI2)
      ENDIF
*
      RETURN
      END
*
* generate Millepede-II record
*
      SUBROUTINE GBLMP2(IRET)
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DIMENSION DERLC(MXPAR),DERSD(MXAUX2),LABSD(MXAUX2)
*
* Output
*   IRET      : number MillePede measurements in record
*
      IRET=0
      CALL GBLPRP
      IF (NLVL.LT.3) RETURN
*
      NSD=0
      IF (IXSD.NE.0) THEN
         NPAR2=(NPAR*NPAR+NPAR)/2
         DO I=1,NPAR2
            IF (CMATSD(I).NE.0.0D0) THEN
               NSD=NSD+1
               LABSD(NSD)=I
               DERSD(NSD)=SNGL(CMATSD(I))
            ENDIF
         ENDDO
      ENDIF
*
      LDAT=0
      DO I=1,NDAT
        LTEM=IDAT(LDAT+1)
        IPNT=IDAT(LDAT+2)
        IM  =IDAT(LDAT+3)
* local derivatives
        DO K=1,NPAR
          DERLC(K)=0.0
        ENDDO
        DO J=LDAT+4,LDAT+LTEM
          DERLC(IDAT(J))=ADAT(J)
        ENDDO
        IGBL=0
        NGBL=0
        IF (IM.GT.0) IGBL=IDGL(IM,IPNT)
        IF (IGBL.GT.0) NGBL=IDAT(IGBL)
        SIG=1.0/SQRT(ADAT(LDAT+2))
*
        CALL MILLE(NPAR,DERLC,NGBL,ADAT(IGBL+1),IDAT(IGBL+1),
     &             ADAT(LDAT+1),SIG) ! add data
*
        LDAT=LDAT+LTEM
      ENDDO
*
      IF (NSD.GT.0) CALL MILLSD(NSD,DERSD,LABSD)
*
      CALL ENDLE ! complete, write record (many sets)
      IRET=NDAT
*
      RETURN
      END
*
      SUBROUTINE GBLJAC(IPOINT,ITRANS)
*
* get jacobian (transposed): broken lines parameter (q/p,..,u_i,..)
*     to local parameter (q/p,alpha,u)
*
* Input
*    IPOINT    : (signed) point 
*    ITRANS   =0 not transposed
*             =1 transpoded
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION 
     &   W(2,2),WJ(2,2),DW(2),
     &   WP(2,2),WJP(2,2),WN(2,2),WJN(2,2),
     &   DWP(2),DWN(2),WJS(2,2),WJI(2,2),
     &   WIP(2,2),WIN(2,2),DIP(2),DIN(2),
     &   WPP(2,2),WPN(2,2),DPP(2),DPN(2)
*
      INDX(I,J)=((I-1)*NPAR+J)*(1-ITRANS)+((J-1)*MP+I)*ITRANS
*
      MP=5+NLOC
      JPNT=IABS(IPOINT)
*
      JOFF=IOFF(JPNT)
*
      DO K=1,MP*NPAR
         DJAC(K)=0.0D0 ! reset Jacobi matrix
      ENDDO
*
      L0=ICRD0
      IF (JOFF.GT.0) THEN
* at offset
         IF (IPOINT.GT.0) THEN ! forward
            IDIR=2
            IF (JOFF.EQ.NOFF) IDIR=1
         ELSE                  ! backward
            IDIR=1
            IF (JOFF.EQ.1) IDIR=2
         ENDIF
         KOFF=JOFF+2*IDIR-3 ! 2nd offset for slope measurement
         CALL GBLDER(JPNT,IDIR,W,WJ,DW)
*
         IP1=NCRV+NLOC+NDIM*(JOFF-1)
         IP2=NCRV+NLOC+NDIM*(KOFF-1)
*
         DO L=1,NDIM
            DJAC(INDX(2  ,IP1+L))=-WJ(1,L+L0)
            DJAC(INDX(3  ,IP1+L))=-WJ(2,L+L0)
            DJAC(INDX(2  ,IP2+L))= W(1,L+L0)
            DJAC(INDX(3  ,IP2+L))= W(2,L+L0)
            DJAC(INDX(3+L,IP1+L))= 1.0D0
         ENDDO
         IF (NCRV.GT.0) THEN                        ! correct for curvature
            DJAC(INDX(2,1)) = -DW(1)
            DJAC(INDX(3,1)) = -DW(2)
         ENDIF
*
      ELSE
* between offsets
         CALL GBLDER(JPNT,-1,WP,WJP,DWP)
         CALL GBLDER(JPNT, 2,WN,WJN,DWN)
         IP1=NDIM*(-JOFF-2)+NCRV+NLOC
         IP2=IP1+NDIM
         DO J=1,2
            DO K=1,2
               WJS(J,K)=WJP(J,K)+WJN(J,K) !  W- * J- + W+ * J+
            ENDDO
         ENDDO
         CALL GBLINV(2,WJS,WJI)          ! (W- * J- + W+ * J+)^-1 (=N)
*        derivatives for u_int
         CALL GBLMLT(2,2,WJI,WP,WIP)     !  N * W-
         CALL GBLMLT(2,2,WJI,WN,WIN)     !  N * W+
         CALL GBLMLT(2,1,WJI,DWP,DIP)    !  N * W- * d-
         CALL GBLMLT(2,1,WJI,DWN,DIN)    !  N * W+ * d+
*        derivatives for alpha_int
         CALL GBLMLT(2,2,WJN,WIP,WPP)    !  W+ * J+ * N * W-
         CALL GBLMLT(2,2,WJP,WIN,WPN)    !  W- * J- * N * W+
         CALL GBLMLT(2,1,WJN,DIP,DPP)    !  W+ * J+ * N * W- * d-
         CALL GBLMLT(2,1,WJP,DIN,DPN)    !  W- * J- * N * W+ * d+
*
         DO L=1,NDIM
*           du_int/du-
            DJAC(INDX(4,IP1+L))= WIP(1,L+L0)
            DJAC(INDX(5,IP1+L))= WIP(2,L+L0)
*           du_int/du+
            DJAC(INDX(4,IP2+L))= WIN(1,L+L0)
            DJAC(INDX(5,IP2+L))= WIN(2,L+L0)
*           dalpha_int/du-
            DJAC(INDX(2,IP1+L))=-WPP(1,L+L0) 
            DJAC(INDX(3,IP1+L))=-WPP(2,L+L0)
*           dalpha_int/du+
            DJAC(INDX(2,IP2+L))= WPN(1,L+L0) 
            DJAC(INDX(3,IP2+L))= WPN(2,L+L0)
         ENDDO
         IF (NCRV.GT.0) THEN              ! correct for curvature
*           du_int/dQbyP
            DJAC(INDX(4,1)) =-DIP(1)-DIN(1)
            DJAC(INDX(5,1)) =-DIP(2)-DIN(2)
*           dalpha_int/dQbyP
            DJAC(INDX(2,1)) = DPP(1)-DPN(1)
            DJAC(INDX(3,1)) = DPP(2)-DPN(2)
         ENDIF
*
      ENDIF
*
      IF (NCRV.GT.0) THEN   ! curvature
         DJAC(INDX(1,1))=1.0D0
      ENDIF
*
      DO K=1,NLOC           ! local parameters
         DJAC(INDX(K+5,NCRV+K))= 1.0D0
      ENDDO
*
      RETURN
      END
*
*
      SUBROUTINE GBLRES(IPOINT,DPAR,DCOV)
*
* get parameters and covariance matrix at point
*
* Input
*    IPOINT    : (signed) point 
*                (<0: side towards previous point)
*                (>0: side towards next point)
* Output
*    DPAR(NP)  : corrections (double precision)
*                (NP is number of track parameters: 
*                 5 + #local par.)
*    DCOV(NP2) : covariance matrix (double precision,
*                symmetric storage mode, NP2=(NP+1)*NP/2)

*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION DPAR(*),DCOV(*)
*
      MP=5+NLOC
      MP2=(MP*MP+MP)/2
      DO K=1,MP
         DPAR(K)=0.0D0
      ENDDO
      DO K=1,MP2
         DCOV(K)=0.0D0
      ENDDO
*
      JPNT=IABS(IPOINT)
      IF (JPNT.LT.1.OR.JPNT.GT.NPNT) THEN
        IF (LPR.GT.0) PRINT *, ' GBLRES invalid point ', IPOINT
        RETURN
      ENDIF
*
      IF (NLVL.LT.4) RETURN ! fit not yet performed or failed
*
      CALL GBLJAC(IPOINT,0) ! get jacobian broken lines -> local
*
      CALL DBGAX(DJAC,BVEC,DPAR,MP,NPAR)
      CALL DBAVAT(CMAT,DJAC,DCOV,NPAR,MP)
*
      RETURN
      END
*
      SUBROUTINE GBLADX(IPOINT,DPRC)
*
* add (inverse covariance matrix from) external seed
*
* Input
*    IPOINT    : (signed) point 
*                (<0: side towards previous point)
*                (>0: side towards next point)
*    DPRC(NP2) :  precision matrix (inverse covariance) from external seed
*                (NP is number of track parameters: 5 + #local par.,
*                 double precision, symmetric storage mode, NP2=(NP+1)*NP/2)
*
#include "gbldef.inc"
#include "gblcom.inc"
*
      DOUBLE PRECISION DPRC(*)
*
      JPNT=IABS(IPOINT)
      IF (JPNT.LT.1.OR.JPNT.GT.NPNT) THEN
        IF (LPR.GT.0) PRINT *, ' GBLADX invalid point ', IPOINT
        RETURN
      ENDIF
*
      IF (NLVL.GE.3) RETURN ! fit already prepared or performed 
*
      MP=5+NLOC
      MP2=(MP*MP+MP)/2
      IXSD=IPOINT
      DO K=1,MP2
         DPSD(K)=DPRC(K)
      ENDDO
*
      RETURN
      END
*
      SUBROUTINE GBLINV(ND,A,B) ! invert matrix A (dimension ND), B=A^-1
      DOUBLE PRECISION A(*),B(*)
      IF (ND.LE.1) THEN
        B(1)=1.0/A(1)
      ELSE
        DET=A(1)*A(4)-A(2)*A(3)
        B(1)= A(4)/DET
        B(2)=-A(2)/DET
        B(3)=-A(3)/DET
        B(4)= A(1)/DET
      ENDIF
      RETURN
      END
*
      SUBROUTINE GBLMLT(N1,N2,A,B,C) ! matrix multiplication C=A*B 
*     A is N1*N1, B N1*N2, C=N1*N2
      DOUBLE PRECISION A(*),B(*),C(*),SUM
      DO I=1,N1
        DO J=1,N2
          IOFF=I
          JOFF=(J-1)*N1
          SUM=0.0D0
          DO K=1,N1
            SUM=SUM+A(IOFF)*B(JOFF+K)
            IOFF=IOFF+N1
          ENDDO
          C((J-1)*N1+I)=SUM
        ENDDO
      ENDDO
      RETURN 
      END

