!

! Code converted using TO_F90 by Alan Miller
! Date: 2012-03-03  Time: 16:59:41

! C. Kleinwort,  DESY-FH1, www.terascale.de
! March 2011, Analysis Centre: Statistics Tools Group

!> \mainpage General information
!!
!!  \section intro_sec Introduction
!!
!!  For a track with an initial trajectory from a prefit of the
!!  measurements (internal seed) or an external prediction
!!  (external seed) the description of multiple scattering is
!!  added by offsets in a local system. Along the initial
!!  trajectory points are defined with can describe a measurement
!!  or a (thin) scatterer or both. Measurements are arbitrary
!!  functions of the local track parameters at a point (e.g. 2D:
!!  position, 4D: slope+position). The refit provides corrections
!!  to the local track parameters (in the local system) and the
!!  corresponding covariance matrix at any of those points.
!!  Outliers can be down-weighted by use of M-estimators.
!!
!!  The broken lines trajectory is defined by (2D) offsets at the
!!  first and last point and all points with a scatterer. The
!!  prediction for a measurement is obtained by interpolation of
!!  the enclosing offsets and for triplets of adjacent offsets
!!  kink angles are determined. This requires for all points the
!!  jacobians for propagation to the previous and next offset.
!!  These are calculated from the point-to-point jacobians along
!!  the initial trajectory.
!!
!!  Additional local or global parameters can be added and the
!!  trajectories can be written to special binary files for
!!  calibration and alignment with Millepede-II.
!!  (V. Blobel, NIM A, 566 (2006), pp. 5-13).
!!
!!  The conventions for the coordinate systems follow:
!!  Derivation of Jacobians for the propagation of covariance
!!  matrices of track parameters in homogeneous magnetic fields
!!  A. Strandlie, W. Wittek, NIM A, 566 (2006) 687-698.
!!
!!  \section call_sec Calling sequence
!!
!!    -# Initialize trajectory:\n
!!            <tt>CALL gblini(..)</tt>
!!    -# For all points on initial trajectory:
!!        - Create points and add appropriate attributes:\n
!!           - <tt>CALL gbladp(..)</tt>
!!           - <tt>CALL gbladm(..)</tt>
!!           - Add additional local or global parameters to measurement:\n
!!             - <tt>CALL gbladl(..)</tt>
!!             - <tt>CALL gbladg(..)</tt>
!!           - <tt>CALL gblads(..)</tt>
!!    -# Optionally add external seed:\n
!!            <tt>CALL gbladx(..)</tt>
!!    -# Construct and fit trajectory,
!!       get Chi2, Ndf (and weight lost by M-estimators):\n
!!            <tt>CALL gblfit(..)</tt>
!!    -# For any point on initial trajectory:
!!        - Get corrections and covariance matrix for track parameters:\n
!!            <tt>CALL gblres(..)</tt>
!!    -# Optionally write trajectory to MP binary file:\n
!!            <tt>CALL gblmp2(..)</tt>
!!
!!  \section impl_sec Implementation
!!
!!  Linear algebra routines are taken from Millepede-II
!!  (by V. Blobel, University Hamburg).
!!  Only 2D (or 1D) measurements are implemented.
!!
!!  \section ref_sec References
!!    - V. Blobel, C. Kleinwort, F. Meier,
!!      Fast alignment of a complex tracking detector using advanced track models,
!!      Computer Phys. Communications (2011), doi:10.1016/j.cpc.2011.03.017
!!    - C. Kleinwort, General Broken Lines as advanced track fitting method,
!!      NIM A, 673 (2012), 107-110, doi:10.1016/j.nima.2012.01.024
!!
!! \file
!! Trajectory construction and fitting.

!> Perform fit of trajectory.
!!
!! Optionally iterate for outlier down-weighting.
!! \param [in] CDW    string defining iterations for outlier down weighting,
!!                    one char per iteration (C: Cauchy, H: Huber, T: Tukey)
!! \param [out] MRANK rank of measurements
!! \param [out] NP    number of track parameter at given point
!! \param [out] NDF   degrees of freedom
!! \param [out] CHI2  Chi2
!! \param [out] WLS   lost measurements: N-sum(weight)
!!

SUBROUTINE gblfit(cdw,mrank,np,ndf,chi2,wls)

    use gblmod


    CHARACTER (LEN=*), INTENT(IN)            :: cdw
    INTEGER, INTENT(OUT)                     :: mrank
    INTEGER, INTENT(OUT)                     :: np
    INTEGER, INTENT(OUT)                     :: ndf
    REAL, INTENT(OUT)                        :: chi2
    REAL, INTENT(OUT)                        :: wls




    mrank=0
    np=0
    ndf=-1
    chi2=0.
    wls=0.
    lcdw=len_trim(cdw)
    iter=0

    CALL gblprp
    IF (nlvl < 3) RETURN

    iterate: DO
        CALL gblmat
        mrank=ncrv+2*ndim ! rank of trajectory
        np   =5+nloc      ! number of track parameters

        IF (lpr > 1) PRINT *, ' GBLFIT: NPAR, NBND, NBDR ', npar, nbnd, ncrv+nloc
        inv=1 ! calc band part of covariance matrix
        CALL sqmibb(cmat,bvec,npar,ncrv+nloc,nbnd,inv,nrank)
        IF (nrank /= npar) THEN ! rank deficit
            IF (lpr > 1) PRINT *, ' GBLFIT: rank deficit ', nrank, npar
            mrank=-mrank
            chi2=0.
            ndf=0
            als=0.
            RETURN
        END IF

        iter=iter+1
        idwm=0 ! down weighting method (Tukey, Huber, Cauchy)
        IF (iter <= lcdw) idwm=(INDEX('TtHhCc',cdw(iter:iter))+1)/2

        CALL gblch2(idwm,chi2,swgt)
        ndf=ndat-npar
        ! iterate ?
        IF (idwm <= 0) EXIT
    END DO iterate

    wls=FLOAT(ndat)-swgt
    nlvl=4 ! fitted

    RETURN
END SUBROUTINE gblfit

!> Prepare broken lines trajectory.
!!

SUBROUTINE gblprp

    use gblmod

    ! index in symmetric matrix
    ijsym(i,j)=MIN(i,j)+(MAX(i,j)*MAX(i,j)-MAX(i,j))/2

    IF (nlvl == 1) call gblcjc ! calc jacobians to prev/next offset
    IF (nlvl /= 2) RETURN      ! no points with jacobians added
    IF (noff < 2) RETURN       ! too few parameters
    ! number of parameters
    nloc=MIN(nloc,mxbdr-ncrv) ! limit number of (external) local par.
    npar=noff*ndim+ncrv+nloc  ! number of parameters
    IF (npar > mxpar) RETURN ! too many parameters

    CALL gbldat

    nlvl=3 ! prepared

    RETURN
END SUBROUTINE gblprp

!> Calculate broken lines from point to point jacobians.
!!
SUBROUTINE gblcjc

    use gblmod

    DOUBLE PRECISION :: ajaci(5,5),ajacn(5,5),ajact(5,5),dsum

    !     forward propagation (all)

    nstep=0
    lpnt=1                  ! last point with offset

    DO ipnt=2,npnt
        !        full jacobian (to next point)
        k=0
        DO i=1,5
            DO j=1,5
                k=k+1
                ajacn(i,j)=djp2p(k,ipnt)
            END DO
        END DO

        IF (nstep == 0) THEN  ! start integrated jacobian
            DO i=1,5
                DO j=1,5
                    ajaci(i,j)=ajacn(i,j)
                END DO
            END DO
        ELSE                  ! update integrated jacobian
            DO i=1,5           ! multiply
                DO j=1,5
                    dsum=0.0D0
                    DO l=1,5
                        dsum=dsum+ajacn(i,l)*ajaci(l,j)
                    END DO
                    ajact(i,j)=dsum
                END DO
            END DO
            DO i=1,5           ! copy
                DO j=1,5
                    ajaci(i,j)=ajact(i,j)
                END DO
            END DO
        END IF
        nstep=nstep+1
        CALL gblinu(ajaci,ajact)
        !        IPNT -> PREV
        ! J-
        jmat(1,ipnt)=ajact(4,4)
        jmat(2,ipnt)=ajact(4,5)
        jmat(3,ipnt)=ajact(5,4)
        jmat(4,ipnt)=ajact(5,5)
        ! S-
        smat(1,ipnt)=ajact(4,2)
        smat(2,ipnt)=ajact(4,3)
        smat(3,ipnt)=ajact(5,2)
        smat(4,ipnt)=ajact(5,3)
        ! d-
        dvec(1,ipnt)=ajact(4,1)
        dvec(2,ipnt)=ajact(5,1)

        IF (ioff(ipnt) > 0) THEN
            !           LPNT -> NEXT
            ! J+
            jmat(5,lpnt)=ajaci(4,4)
            jmat(6,lpnt)=ajaci(4,5)
            jmat(7,lpnt)=ajaci(5,4)
            jmat(8,lpnt)=ajaci(5,5)
            ! S+
            smat(5,lpnt)=ajaci(4,2)
            smat(6,lpnt)=ajaci(4,3)
            smat(7,lpnt)=ajaci(5,2)
            smat(8,lpnt)=ajaci(5,3)
            ! d+
            dvec(3,lpnt)=ajaci(4,1)
            dvec(4,lpnt)=ajaci(5,1)

            nstep=0    ! restart integration
            lpnt=ipnt
        END IF
    END DO

    !     backward propagation (without scatterers)

    nstep=0

    DO ipnt=npnt-1,1,-1

        IF (ioff(ipnt) > 0) THEN ! skip offsets
            nstep=0
            GO TO 999
        END IF
        !        full jacobian (to next point)
        k=0
        DO i=1,5
            DO j=1,5
                k=k+1
                ajacn(i,j)=djp2p(k,ipnt+1)
            END DO
        END DO

        IF (nstep == 0) THEN  ! start integrated jacobian
            DO i=1,5
                DO j=1,5
                    ajaci(i,j)=ajacn(i,j)
                END DO
            END DO
        ELSE                  ! update integrated jacobian
            DO i=1,5           ! multiply
                DO j=1,5
                    dsum=0.0D0
                    DO l=1,5
                        dsum=dsum+ajaci(i,l)*ajacn(l,j)
                    END DO
                    ajact(i,j)=dsum
                END DO
            END DO
            DO i=1,5            ! copy
                DO j=1,5
                    ajaci(i,j)=ajact(i,j)
                END DO
            END DO
        END IF
        nstep=nstep+1
        !        IPNT -> NEXT
        ! J+
        jmat(5,ipnt)=ajaci(4,4)
        jmat(6,ipnt)=ajaci(4,5)
        jmat(7,ipnt)=ajaci(5,4)
        jmat(8,ipnt)=ajaci(5,5)
        ! S+
        smat(5,ipnt)=ajaci(4,2)
        smat(6,ipnt)=ajaci(4,3)
        smat(7,ipnt)=ajaci(5,2)
        smat(8,ipnt)=ajaci(5,3)
        ! d+
        dvec(3,ipnt)=ajaci(4,1)
        dvec(4,ipnt)=ajaci(5,1)

999 CONTINUE
    END DO

    nlvl=2
    RETURN
END SUBROUTINE gblcjc

!> Initialize.
!!
!! \param [in] LPRNT  print level
!! \param [in] ICOORD coordinate (1: u_1, 2: u_2) to use
!!

SUBROUTINE gblini(lprnt,icoord)

    use gblmod

    INTEGER, INTENT(IN)                      :: lprnt
    INTEGER, OPTIONAL, INTENT(IN)            :: icoord

    DATA ifirst / 1 /

    ndim=2  ! 2D offsets for track in 3D
    icrd0=0
    !    print *, " test1 ", Present(icoord)

    IF (PRESENT(icoord)) THEN
        IF (icoord > 0) THEN
            ndim=1  ! 1D offsets for track in 2D
            icrd0=0
            IF (icoord > 1) icrd0=1
        END IF
    END IF

    nlvl=0
    lpr=lprnt
    ncrv=1  ! with Q/P
    nloc=0
    npnt=0
    lpnt=0
    nmeas=0
    nscat=0
    lscat=0
    noff=0
    nknk=0
    mdat=mxdat
    ixsd=0  ! external seed parameter offset

    IF (ifirst == 1) THEN
        ifirst=0
        IF (lpr > 0) PRINT *, ' GENBRL $Rev: 23 $'
    END IF

    RETURN
END SUBROUTINE gblini

!> Add point to trajectory.
!!
!! \param [in]  AJAC   jacobian from previous point
!! \param [out] IPOINT identifier
!!

SUBROUTINE gbladp(ajac,ipoint)

    use gblmod

    DOUBLE PRECISION, INTENT(IN)             :: ajac(5,5)
    INTEGER, INTENT(OUT)                     :: ipoint

    ipoint=0
    IF (npnt >= mxpnt) THEN
        IF (lpr > 0) PRINT *, ' GBLADP: too many points, ignoring new point '
        RETURN
    END IF
    npnt=npnt+1
    imeas(npnt)=0
    iscat(npnt)=0
    ioff(npnt)=0
    iknk(npnt)=0
    idlc(1,npnt)=0
    idlc(2,npnt)=0
    idgl(1,npnt)=0
    idgl(2,npnt)=0
    ipoint=npnt
    k=0
    DO i=1,5
        DO j=1,5
            k=k+1
            djp2p(k,ipoint)=ajac(i,j)
        END DO
    END DO
    DO k=1,8
        jmat(k,npnt)=0.0D0
        smat(k,npnt)=0.0D0
    END DO
    DO k=1,4
        dvec(k,npnt)=0.0D0
    END DO

    IF (lpnt > 0) THEN
        IF (iscat(lpnt) <= 0.AND.noff > 1) THEN
            noff=noff-1
            ioff(lpnt)=-ioff(lpnt)
        END IF
    END IF
    noff=noff+1
    ioff(npnt)=noff
    lpnt=npnt
    IF (noff > 2.AND.lscat > 0) THEN
        nknk=nknk+1
        iknk(lscat)=nknk
        lscat=0
    END IF

    nlvl=1 ! adding points
    RETURN
END SUBROUTINE gbladp

!> Add 2D measurement to current point.
!!
!! \param [in] PROJ   projection matrix of measurement directions
!!                    into local system (dm/du)
!! \param [in] RES    residuals (m)
!! \param [in] PREC   diagonal of inverse covariance matrix
!!

SUBROUTINE gbladm(proj,res,prec)

    use gblmod


    DOUBLE PRECISION, INTENT(IN)             :: proj(2,2)
    REAL, INTENT(IN)                         :: res(2)
    REAL, INTENT(IN)                         :: prec(2)



    IF (npnt <= 0) THEN
        IF (lpr > 0) PRINT *, ' GBLADM: no point defined'
        RETURN
    END IF
    IF (nmeas+nscat >= mxpnt) THEN
        IF (lpr > 0) PRINT *,  &
        ' GBLADM: too many measurement+scatterers ', nmeas+nscat
        RETURN
    END IF
    IF (imeas(npnt) <= 0) THEN
        nmeas=nmeas+1
        imeas(npnt)=nmeas
        pmat(1,nmeas)=proj(1,1)
        pmat(2,nmeas)=proj(1,2)
        pmat(3,nmeas)=proj(2,1)
        pmat(4,nmeas)=proj(2,2)
        yvec(1,nmeas)=res(1)
        yvec(2,nmeas)=res(2)
        wvec(1,nmeas)=prec(1)
        wvec(2,nmeas)=prec(2)
    ELSE
        IF (lpr > 0) PRINT *,  &
        ' GBLADM: measurement already defined for point ', npnt
    END IF

    RETURN
END SUBROUTINE gbladm

!> Add local derivatives to measurement.
!!
!! \param [in]  NDER    number of local derivatives
!! \param [in]  DER     local derivatives
!! \param [out] IRET    number of non zero derivatives added

SUBROUTINE gbladl(nder,der,iret)

    use gblmod

    INTEGER, INTENT(IN)                      :: nder
    REAL, INTENT(IN)                         :: der(2,nder)
    INTEGER, INTENT(OUT)                     :: iret


    iret=0
    IF (npnt <= 0) THEN
        IF (lpr > 0) PRINT *, ' GBLADL: no point defined'
        RETURN
    END IF

    DO im=1,2
  
        IF (idlc(im,npnt) == 0) THEN      ! local derivs not yet defined
    
            IF (ldat > mdat-nder) RETURN  ! enough space left
            mdat0=mdat
            DO k=nder,1,-1
                IF (der(im,k) /= 0.0) THEN
                    nloc=MAX(nloc,k)
                    idat(mdat)=k
                    adat(mdat)=der(im,k)
                    mdat=mdat-1
                END IF
            END DO
    
            idlc(im,npnt)=mdat
            idat(mdat)=mdat0-mdat
            adat(mdat)=0.0
            mdat=mdat-1
            iret=iret+mdat0-mdat
    
        ELSE
            IF (lpr > 0) PRINT *,  &
            ' GBLADL: local derivatives already defined for point ', npnt
        END IF
  
    END DO

    RETURN
END SUBROUTINE gbladl

!> Add global derivatives to measurement.
!!
!! \param [in]  NDER    number of local derivatives
!! \param [in]  LDER    labels for global derivatives
!! \param [in]  DER     local derivatives
!! \param [out] IRET    number of non zero derivatives added

SUBROUTINE gbladg(nder,lder,der,iret)

    use gblmod

    INTEGER, INTENT(IN)                      :: nder
    INTEGER, INTENT(IN)                         :: lder(nder)
    REAL, INTENT(IN)                         :: der(2,nder)
    INTEGER, INTENT(OUT)                     :: iret


    iret=0
    IF (npnt <= 0) THEN
        IF (lpr > 0) PRINT *, ' GBLADG: no point defined'
        RETURN
    END IF

    DO im=1,2
  
        IF (idgl(im,npnt) == 0) THEN      ! local derivs not yet defined
    
            IF (ldat > mdat-nder) RETURN  ! enough space left
            mdat0=mdat
            DO k=nder,1,-1
                IF (der(im,k) /= 0.0) THEN
                    idat(mdat)=lder(k)
                    adat(mdat)=der(im,k)
                    mdat=mdat-1
                END IF
            END DO
    
            idgl(im,npnt)=mdat
            idat(mdat)=mdat0-mdat
            adat(mdat)=0.0
            mdat=mdat-1
            iret=iret+mdat0-mdat
    
        ELSE
            IF (lpr > 0) PRINT *,  &
            ' GBLADG: global derivatives already defined for point ', npnt
        END IF
  
    END DO

    RETURN
END SUBROUTINE gbladg

!> Add (thin) scatterer to current point.
!!
!! \param [in] RES    values for initial kinks (in case of iterating)
!! \param [in] PREC   diagonal of inverse (multiple scattering)
!!                    covariance matrix

SUBROUTINE gblads(res,prec)

    use gblmod

    REAL, INTENT(IN)                         :: res(2)
    REAL, INTENT(IN)                         :: prec(2)


    IF (npnt <= 0) THEN
        IF (lpr > 0) PRINT *, ' GBLADS: no point defined'
        RETURN
    END IF
    IF (nmeas+nscat >= mxpnt) THEN
        IF (lpr > 0) PRINT *,  &
        ' GBLADS: too many measurement+scatterers ', nmeas+nscat
        RETURN
    END IF
    IF (prec(1) <= 0.0.OR.prec(2) <= 0.0) THEN
        IF (lpr > 0) PRINT *, ' GBLADS: invalid scattering precision ', prec
        RETURN
    END IF

    IF (iscat(npnt) <= 0) THEN
        jscat=mxpnt-nscat
        nscat=nscat+1
        iscat(npnt)=jscat
        lscat=npnt
        yvec(1,jscat)=res(1)
        yvec(2,jscat)=res(2)
        wvec(1,jscat)=prec(1)
        wvec(2,jscat)=prec(2)
    ELSE
        IF (lpr > 0) PRINT *, ' GBLADM: scatterer already defined for point ', npnt
    END IF

    RETURN
END SUBROUTINE gblads

!> Dump trajectory definition.

SUBROUTINE gbldmp

    use gblmod

    IF (npnt <= 0) THEN
        PRINT *, ' GBLDMP no trajectory defined '
        RETURN
    END IF

    PRINT *
    PRINT *, '    GBLDMP trajectory definition '
    PRINT *, '-------------------------------------'
    PRINT *, ' number of points       ', npnt
    PRINT *, ' number of offsets      ', noff
    PRINT *, ' number of scatterers   ', nscat
    PRINT *, ' number of measurements ', nmeas
    PRINT *, ' number of kinks        ', nknk
    IF (nloc > 0) PRINT *, ' number of local par.   ', nloc
    PRINT *

    DO i=1,npnt
        PRINT *, ' Point          ', i
        IF (ioff(i) > 0)  PRINT *, '    Offset      ', ioff(i)
        IF (iscat(i) > 0) PRINT *, '    Scatterer   ', iscat(i)
        IF (imeas(i) > 0) THEN
            IF (ioff(i) > 0) THEN
                PRINT *, '    Measurement ', imeas(i), ' using offset  ',ioff(i)
            ELSE
                PRINT *, '    Measurement ', imeas(i),  &
                ' using offsets ',-ioff(i)-1,' ,',-ioff(i)
            END IF
        END IF
        IF (iknk(i) > 0)  PRINT *, '    Kink        ', iknk(i),  &
        ' using offsets ',ioff(i)-1,'..',ioff(i)+1
        IF (imeas(i) <= 0.AND.iscat(i) <= 0) THEN
            IF (ioff(i) > 0) THEN
                PRINT *, '    Prediction               ', ' using offset  ',ioff(i)
      
            ELSE
                PRINT *, '    Prediction               ',  &
                ' using offsets ',-ioff(i)-1,' ,',-ioff(i)
            END IF
        END IF
    END DO

    RETURN
END SUBROUTINE gbldmp

!> Generate DATa record.

SUBROUTINE gbldat

    use gblmod

    DOUBLE PRECISION :: tc2l(2,2), wp(2,2),wn(2,2),wjp(2,2),wjn(2,2),  &
    dwp(2),dwn(2),wjs(2,2),wji(2,2),pn(2,2), wpp(2,2),wpn(2,2),dws(2),dps(2),det
    DOUBLE PRECISION :: diag(mxaux), eigen(mxaux*mxaux), work(mxaux), &
    eigent(mxaux*mxaux), seed(mxaux*mxaux)
    INTEGER          :: iwork(mxaux)

    l0=icrd0

    ldwm=0 ! start without downweighting
    ! loop over points
    ndat=0 ! number of data items
    ldat=0 ! length of data
    ipmn=npar+1 ! minimum parameter with non zero derivatives
    DO jpnt=1,npnt
  
        joff=ioff(jpnt)               ! offset at point
        ipar0=ndim*(joff-1)+ncrv+nloc ! offset in parameter number
        ! kink at point ?
        IF (iknk(jpnt) > 0) THEN
            jscat=iscat(jpnt)
    
            CALL gblder(jpnt,-1,wp,wjp,dwp)
            CALL gblder(jpnt, 2,wn,wjn,dwn)
            DO j=1,2
                dws(j)=dwp(j)+dwn(j)          !  W- * d- + W+ * d+
                DO k=1,2
                    wjs(j,k)=wjp(j,k)+wjn(j,k) !  W- * J- + W+ * J+
                END DO
            END DO
    
            DO m=1,ndim ! pseudo measurements
                kdat=mdat-(ncrv+3*ndim)
                IF (ldat+3 > kdat) GO TO 9999
                kdat0=kdat
      
                ldat0=ldat
                ndat=ndat +1
                idat(ldat+1)=3             ! length of item
                idat(ldat+2)=jpnt          ! point
                idat(ldat+3)=-m            ! pseudo measurement
                adat(ldat+1)=yvec(m,jscat) ! (start) value
                adat(ldat+2)=wvec(m,jscat) ! precision (1/error**2)
                adat(ldat+3)=1.0           ! (additional) weight
                ldat=ldat+3
      
                idat(kdat+1)= 1
                adat(kdat+1)= -SNGL(dws(m))
                kdat=kdat+ncrv
                DO l=1,ndim
                    idx=kdat+l
                    idat(idx)= ipar0+l-ndim
                    adat(idx)= SNGL(  wp(m,l+l0))
                    idx=idx+ndim
                    idat(idx)= ipar0+l
                    adat(idx)= SNGL(-wjs(m,l+l0))
                    idx=idx+ndim
                    idat(idx)= ipar0+l+ndim
                    adat(idx)= SNGL(  wn(m,l+l0))
                END DO
                kdat=kdat+3*ndim
      
                DO k=kdat0+1,kdat
                    IF (adat(k) /= 0.0) THEN ! copy non zero derivatives
                        ldat=ldat+1
                        idat(ldat)=idat(k)
                        adat(ldat)=adat(k)
                        ipmn=MIN(ipmn,idat(ldat))
                    END IF
                END DO
                idat(ldat0+1)=ldat-ldat0
            END DO
    
        END IF
        ! measurement at point ?
        jmeas=imeas(jpnt)        ! measurement at point
        IF (jmeas > 0) THEN
            ip=0
            DO j=1,2
                DO k=1,2
                    ip=ip+1
                    tc2l(j,k)=pmat(ip,jmeas)      !  P
                END DO
            END DO
    
            ipar0=ndim*(joff-1)+ncrv+nloc       ! offset in parameter number
            IF (joff < 0) THEN                 ! need interpolation
                ipar0=ndim*(-joff-2)+ncrv+nloc   ! offset in parameter number
                CALL gblder(jpnt,-1,wp,wjp,dwp)
                CALL gblder(jpnt, 2,wn,wjn,dwn)
                DO j=1,2
                    dws(j)=dwp(j)+dwn(j)          !  W- * d- + W+ * d+
                    DO k=1,2
                        wjs(j,k)=wjp(j,k)+wjn(j,k) !  W- * J- + W+ * J+
                    END DO
                END DO
                !                                               ! (W- * J- + W+ * J+)^-1 (=N)
                det=wjs(1,1)*wjs(2,2)-wjs(1,2)*wjs(2,1)
                wji(1,1)= wjs(2,2)/det
                wji(1,2)=-wjs(1,2)/det
                wji(2,1)=-wjs(2,1)/det
                wji(2,2)= wjs(1,1)/det
                !                                               !  P * N
                pn(1,1)=tc2l(1,1)*wji(1,1)+tc2l(1,2)*wji(2,1)
                pn(1,2)=tc2l(1,1)*wji(1,2)+tc2l(1,2)*wji(2,2)
                pn(2,1)=tc2l(2,1)*wji(1,1)+tc2l(2,2)*wji(2,1)
                pn(2,2)=tc2l(2,1)*wji(1,2)+tc2l(2,2)*wji(2,2)
                !                                               !  P * N * W-
                wpp(1,1)=pn(1,1)*wp(1,1)+pn(1,2)*wp(2,1)
                wpp(1,2)=pn(1,1)*wp(1,2)+pn(1,2)*wp(2,2)
                wpp(2,1)=pn(2,1)*wp(1,1)+pn(2,2)*wp(2,1)
                wpp(2,2)=pn(2,1)*wp(1,2)+pn(2,2)*wp(2,2)
                !                                               !  P * N * W+
                wpn(1,1)=pn(1,1)*wn(1,1)+pn(1,2)*wn(2,1)
                wpn(1,2)=pn(1,1)*wn(1,2)+pn(1,2)*wn(2,2)
                wpn(2,1)=pn(2,1)*wn(1,1)+pn(2,2)*wn(2,1)
                wpn(2,2)=pn(2,1)*wn(1,2)+pn(2,2)*wn(2,2)
                !                                               !  P * N * (W+ * d+ + W- * d-)
                dps(1)=dws(1)*pn(1,1)+dws(2)*pn(1,2)
                dps(2)=dws(1)*pn(2,1)+dws(2)*pn(2,2)
      
            END IF
    
            DO m=1,2
                IF (wvec(m,jmeas) <= 0.0) CYCLE ! no precision ?
                kdat=mdat-(ncrv+2*ndim)
                IF (ldat+3 > kdat) GO TO 9999
                kdat0=kdat
      
                ldat0=ldat
                ndat=ndat+1
                idat(ldat+1)=3               ! length of item
                idat(ldat+2)=jpnt            ! point
                idat(ldat+3)=m               ! measurement
                adat(ldat+1)=yvec(m,jmeas)   ! value
                adat(ldat+2)=wvec(m,jmeas)   ! precision (1/error**2)
                adat(ldat+3)=1.0             ! (additional) weight
                ldat=ldat+3
      
                IF (joff > 0) THEN         ! measurement at offset
                    DO l=1,ndim
                        idat(kdat+1)= ipar0+l
                        adat(kdat+1)= SNGL(tc2l(m,l+l0))
                        kdat=kdat+1
                    END DO
                ELSE                         ! interpolation between offsets
        
                    idat(kdat+1)= 1
                    adat(kdat+1)= -SNGL(dps(m))
                    kdat=kdat+ncrv
                    DO l=1,ndim
                        idx=kdat+l
                        idat(idx)= ipar0+l
                        adat(idx)= SNGL(wpp(m,l+l0))
                        idx=idx+ndim
                        idat(idx)= ipar0+l+ndim
                        adat(idx)= SNGL(wpn(m,l+l0))
                    END DO
                    kdat=kdat+2*ndim
                END IF
                DO k=kdat0+1,kdat
                    IF (adat(k) /= 0.0) THEN ! copy non zero derivatives
                        ldat=ldat+1
                        idat(ldat)=idat(k)
                        adat(ldat)=adat(k)
                        ipmn=MIN(ipmn,idat(ldat))
                    END IF
                END DO
                idat(ldat0+1)=ldat-ldat0
                ! check for local derivatives
                ilcl=idlc(m,jpnt)
                IF (ilcl <= 0) CYCLE
                nlcl=idat(ilcl)
                IF (ldat+nlcl > mdat) CYCLE
                DO k=1,nlcl
                    IF (idat(ilcl+k) <= nloc) THEN
                        ldat=ldat+1
                        idat(ldat)=idat(ilcl+k)+ncrv
                        adat(ldat)=adat(ilcl+k)
                        ipmn=MIN(ipmn,idat(ldat))
                    END IF
                END DO
                idat(ldat0+1)=ldat-ldat0
            END DO
    
        END IF
  
    END DO

    IF (ixsd /= 0) THEN
        IF (dpsd(1) > 0.0D0) ipmn=1 ! external 'seed' for curvature

        nb=ncrv+nloc
        npsd=nb+2*ndim      ! fit parameters
        nplc=5+nloc         ! local track parameters

        CALL gbljac(ixsd,1,ioff1)  ! get transposed jacobian broken lines -> local
        CALL devrot(nplc,diag,eigen,dpsd,work,iwork)
        DO k=1,nplc  ! transpose
            DO l=1,nplc
                eigent((k-1)*nplc+l)=eigen((l-1)*nplc+k)
            END DO
        END DO
        CALL gblmlt(nplc,npsd,eigent,djac,seed)
        ip0=(ioff1-1)*ndim

        DO i=1, nplc
            IF (diag(i) > 0.0) THEN
                IF (ldat+3+npsd > mdat) GO TO 9999
                ldat0=ldat
                ndat=ndat+1                  ! 'virtual' measurement
                idat(ldat+1)=3               ! length of item
                idat(ldat+2)=ixsd            ! point
                idat(ldat+3)=i               ! measurement
                adat(ldat+1)=0.0             ! value
                adat(ldat+2)=diag(i)         ! precision (1/error**2)
                adat(ldat+3)=1.0             ! (additional) weight
                ldat=ldat+3
                DO j=1,npsd
                    IF (seed((j-1)*nplc+i) /= 0.0) THEN
                        ldat=ldat+1
                        IF (j > nb) THEN
                            idat(ldat)=j+ip0
                        ELSE
                            idat(ldat)=j
                        END IF
                        adat(ldat)=seed((j-1)*nplc+i)
                    END IF
                END DO
                idat(ldat0+1)=ldat-ldat0
            END IF
        END DO
    END IF

    !     check minimum parameter with non zero derivatives
    IF (ipmn > ncrv) THEN          ! curvature undefined
        ipoff=ncrv
        ncrv=0
        !        correct parameter indices
        IF (ipoff > 0) THEN
            npar=npar-ipoff
            ldat=0
            DO i=1,ndat
                ltem=idat(ldat+1)
                DO j=ldat+4,ldat+ltem
                    idat(j)=idat(j)-ipoff
                END DO
                ldat=ldat+ltem
            END DO
    
        END IF
    END IF

9999 CONTINUE

     RETURN
 END SUBROUTINE gbldat

 !> get (matrices and vectors for) derivatives.
 !!
 !! \param [in]  IPOINT   point
 !! \param [in]  IDIR     direction
 !! \param [out] W        W (=(+/-)S^-1)
 !! \param [out] WJ       W*J
 !! \param [out] DW       W*d

 SUBROUTINE gblder(ipoint,idir,w,wj,dw)

     use gblmod


     INTEGER, INTENT(IN OUT)                  :: ipoint
     INTEGER, INTENT(IN)                      :: idir
     DOUBLE PRECISION, INTENT(OUT)            :: w(2,2)
     DOUBLE PRECISION, INTENT(OUT)            :: wj(2,2)
     DOUBLE PRECISION, INTENT(OUT)            :: dw(2)

     DOUBLE PRECISION :: jm(2,2),sm(2,2),dv(2), det

     ! (parts of) jacobians to previous/next offset
     jp=(IABS(idir)-1)*2
     ip=jp*2
     DO j=1,2
         DO k=1,2
             ip=ip+1
             jm(j,k)=jmat(ip  ,ipoint)     ! J
             sm(j,k)=smat(ip  ,ipoint)     ! S
         END DO
         jp=jp+1
         dv(j)=dvec(jp,ipoint)            ! d
     END DO

     det=sm(1,1)*sm(2,2)-sm(1,2)*sm(2,1)
     IF (idir < 0) det=-det
     !                                         ! W
     w(1,1)= sm(2,2)/det
     w(1,2)=-sm(1,2)/det
     w(2,1)=-sm(2,1)/det
     w(2,2)= sm(1,1)/det
     !                                         ! W*J
     wj(1,1)=w(1,1)*jm(1,1)+w(1,2)*jm(2,1)
     wj(1,2)=w(1,1)*jm(1,2)+w(1,2)*jm(2,2)
     wj(2,1)=w(2,1)*jm(1,1)+w(2,2)*jm(2,1)
     wj(2,2)=w(2,1)*jm(1,2)+w(2,2)*jm(2,2)
     !                                         ! W*d
     dw(1)=dv(1)*w(1,1)+dv(2)*w(1,2)
     dw(2)=dv(1)*w(2,1)+dv(2)*w(2,2)

     RETURN
 END SUBROUTINE gblder

 !> Build MATrix and rhs vector.

 SUBROUTINE gblmat

     use gblmod

     DOUBLE PRECISION :: dval,dwgh

     DATA mpar / 0 /, lpar / 0 /, lbdr / 0 /, lbnd / 0 /
     ! index in symmetric matrix
     ijsym(i,j)=(i*i-i)/2+j

     DO i=1,npar
         bvec(i)=0.0D0
     END DO
     !     'smart' clear
     IF (npar > mpar) THEN
         mpar2=(mpar*mpar+mpar)/2
         DO i=mpar2+1,npar2
             cmat(i)=0.0D0
         END DO
         mpar=npar
     END IF

     IF (lpar > 0) THEN
         DO i=1,lpar
             ij0=(i*i-i)/2
             DO j=1,MIN(lbdr,i) ! clear border
                 cmat(ij0+j)=0.0D0
             END DO
             DO j=1,MAX(i-lbnd,1),i ! clear band
                 cmat(ij0+j)=0.0D0
             END DO
         END DO
     END IF

     ldat=0
     nbdr=ncrv+nloc
     DO i=1,ndat
         ltem=idat(ldat+1)
         dval=DBLE(adat(ldat+1))              ! value
         dwgh=DBLE(adat(ldat+3)*adat(ldat+2)) ! (total) weight
         DO j=ldat+4,ldat+ltem ! update matrix
             ij=idat(j)
             bvec(ij)=bvec(ij)+DBLE(adat(j))*dval*dwgh
             DO k=ldat+4,j
                 ik=idat(k)
                 jk=ijsym(ij,ik) ! parameter labels must be sorted: IK<=IJ !
                 cmat(jk)=cmat(jk)+DBLE(adat(j))*DBLE(adat(k))*dwgh
                 IF (ik > nbdr) nbnd=MAX(nbnd,ij-ik) ! update band width
             END DO
         END DO
         ldat=ldat+ltem
     END DO

     lpar=npar
     lbdr=nbdr
     lbnd=nbnd

     RETURN
 END SUBROUTINE gblmat

 !> Calculate Chi2.
 !!
 !! \param [in]  IDWM  down-weighting method (0-3)
 !! \param [out] CHI2  Chi2
 !! \param [out] SWGT  sum of weights

 SUBROUTINE gblch2(idwm,chi2,swgt)

     use gblmod


     INTEGER, INTENT(IN)                      :: idwm
     REAL, INTENT(OUT)                        :: chi2
     REAL, INTENT(OUT)                        :: swgt


     DIMENSION dwint(0:3) ! Integral(weight*normal_distribution)
     DATA dwint / 1.0, 0.8737, 0.9326, 0.8228 /

     chi2=0.0
     swgt=0.0

     ldat=0
     DO i=1,ndat
         ltem=idat(ldat+1)
         val=adat(ldat+1)      ! value
         wgh=adat(ldat+3)      ! down weighting
         brl=0.0               ! predction from broken lines
         DO j=ldat+4,ldat+ltem
             ij=idat(j)
             brl=brl+SNGL(bvec(ij)*adat(j))
         END DO
         sres=ABS(brl-val)*SQRT(adat(ldat+2))
         chi2=chi2+sres*sres*wgh
         swgt=swgt+wgh
         !       down weighting
         IF (idwm == 1) THEN
             IF(sres < 4.6851) THEN          ! Tukey
                 wgh=(1.0-0.045558*sres*sres)**2  ! 1/4.6851**2
             ELSE
                 wgh=0.0
             END IF
         ELSE IF (idwm == 2) THEN
             IF (sres < 1.345) THEN          ! Huber
                 wgh=1.0
             ELSE
                 wgh=1.345/sres
             END IF
         ELSE IF (idwm == 3) THEN
             wgh=1.0/(1.0+(sres/2.3849)**2)   ! Cauchy
         END IF
         adat(ldat+3)=wgh
         ldat=ldat+ltem
     END DO

     chi2=chi2/dwint(ldwm) ! renormalize Chi2
     ldwm=idwm

     RETURN
 END SUBROUTINE gblch2

 !> generate Millepede-II record
 !!
 !! \param [out] IRET  number MillePede measurements in record

 SUBROUTINE gblmp2(iret)

     use gblmod


     INTEGER, INTENT(OUT)                     :: iret

     DIMENSION derlc(mxpar)

     iret=0
     CALL gblprp
     IF (nlvl < 3) RETURN

     ldat=0
     DO i=1,ndat
         ltem=idat(ldat+1)
         ipnt=idat(ldat+2)
         im  =idat(ldat+3)
         ! local derivatives
         DO k=1,npar
             derlc(k)=0.0
         END DO
         DO j=ldat+4,ldat+ltem
             derlc(idat(j))=adat(j)
         END DO
         igbl=0
         ngbl=0
         IF (im > 0) igbl=idgl(im,ipnt)
         IF (igbl > 0) ngbl=idat(igbl)
         sig=1.0/SQRT(adat(ldat+2))
  
         CALL mille(npar,derlc,ngbl,adat(igbl+1),idat(igbl+1),  &
         adat(ldat+1),sig) ! add data
  
         ldat=ldat+ltem
     END DO

     CALL endle ! complete, write record (many sets)
     iret=ndat

     RETURN
 END SUBROUTINE gblmp2

 !> Get jacobian (transposed).
 !! broken lines parameter (q/p,..,u_i,..) to local parameter (q/p,alpha,u)
 !!
 !! \param [in]  IPOINT  (signed) point
 !! \param [in]  ITRANS  =0 not transposed, =1 transposed
 !! \param [out] IOFF1   offsets IOFF1,IOFF1+1 needed
 !!                      to define offset and slope at IPOINT

 SUBROUTINE gbljac(ipoint,itrans,ioff1)

     use gblmod


     INTEGER, INTENT(IN)                      :: ipoint
     INTEGER, INTENT(IN)                      :: itrans
     INTEGER, INTENT(OUT)                     :: ioff1

     DOUBLE PRECISION :: w(2,2),wj(2,2),dw(2),  &
     wp(2,2),wjp(2,2),wn(2,2),wjn(2,2), dwp(2),dwn(2),wjs(2,2),wji(2,2),  &
     wip(2,2),win(2,2),dip(2),din(2), wpp(2,2),wpn(2,2),dpp(2),dpn(2)

     indx(i,j)=((i-1)*np+j)*(1-itrans)+((j-1)*mp+i)*itrans

     np=ncrv+nloc+2*ndim
     mp=5+nloc
     jpnt=IABS(ipoint)

     joff=ioff(jpnt)

     DO k=1,mp*np
         djac(k)=0.0D0 ! reset Jacobi matrix
     END DO

     l0=icrd0
     IF (joff > 0) THEN
         ! at offset
         IF (ipoint > 0) THEN ! forward
             idir=2
             IF (joff == noff) idir=-1
         ELSE                  ! backward
             idir=-1
             IF (joff == 1) idir=2
         END IF
         koff=joff+2*IABS(idir)-3 ! 2nd offset for slope measurement
         CALL gblder(jpnt,idir,w,wj,dw)
  
         IF (idir > 0) THEN
             ioff1=joff
             ioff2=koff
             ip1=ncrv+nloc
             ip2=ncrv+nloc+ndim
         ELSE
             ioff2=joff
             ioff1=koff
             ip2=ncrv+nloc
             ip1=ncrv+nloc+ndim
         END IF
  
         DO l=1,ndim
             djac(indx(2  ,ip1+l))=-wj(1,l+l0)
             djac(indx(3  ,ip1+l))=-wj(2,l+l0)
             djac(indx(2  ,ip2+l))= w(1,l+l0)
             djac(indx(3  ,ip2+l))= w(2,l+l0)
             djac(indx(3+l,ip1+l))= 1.0D0
         END DO
         IF (ncrv > 0) THEN                        ! correct for curvature
             djac(indx(2,1)) = -dw(1)
             djac(indx(3,1)) = -dw(2)
         END IF
  
     ELSE
         ! between offsets
         CALL gblder(jpnt,-1,wp,wjp,dwp)
         CALL gblder(jpnt, 2,wn,wjn,dwn)
         ioff1=-joff-1
         ioff2=-joff
         ip1=ncrv+nloc
         ip2=ncrv+nloc+ndim
  
         DO j=1,2
             DO k=1,2
                 wjs(j,k)=wjp(j,k)+wjn(j,k) !  W- * J- + W+ * J+
             END DO
         END DO
         CALL gblinv(2,wjs,wji)          ! (W- * J- + W+ * J+)^-1 (=N)
         !        derivatives for u_int
         CALL gblmlt(2,2,wji,wp,wip)     !  N * W-
         CALL gblmlt(2,2,wji,wn,win)     !  N * W+
         CALL gblmlt(2,1,wji,dwp,dip)    !  N * W- * d-
         CALL gblmlt(2,1,wji,dwn,din)    !  N * W+ * d+
         !        derivatives for alpha_int
         CALL gblmlt(2,2,wjn,wip,wpp)    !  W+ * J+ * N * W-
         CALL gblmlt(2,2,wjp,win,wpn)    !  W- * J- * N * W+
         CALL gblmlt(2,1,wjn,dip,dpp)    !  W+ * J+ * N * W- * d-
         CALL gblmlt(2,1,wjp,din,dpn)    !  W- * J- * N * W+ * d+
  
         DO l=1,ndim
             !           du_int/du-
             djac(indx(4,ip1+l))= wip(1,l+l0)
             djac(indx(5,ip1+l))= wip(2,l+l0)
             !           du_int/du+
             djac(indx(4,ip2+l))= win(1,l+l0)
             djac(indx(5,ip2+l))= win(2,l+l0)
             !           dalpha_int/du-
             djac(indx(2,ip1+l))=-wpp(1,l+l0)
             djac(indx(3,ip1+l))=-wpp(2,l+l0)
             !           dalpha_int/du+
             djac(indx(2,ip2+l))= wpn(1,l+l0)
             djac(indx(3,ip2+l))= wpn(2,l+l0)
         END DO
         IF (ncrv > 0) THEN              ! correct for curvature
             !           du_int/dQbyP
             djac(indx(4,1)) =-dip(1)-din(1)
             djac(indx(5,1)) =-dip(2)-din(2)
             !           dalpha_int/dQbyP
             djac(indx(2,1)) = dpp(1)-dpn(1)
             djac(indx(3,1)) = dpp(2)-dpn(2)
         END IF
  
     END IF

     IF (ncrv > 0) THEN   ! curvature
         djac(indx(1,1))=1.0D0
     END IF

     DO k=1,nloc           ! local parameters
         djac(indx(k+5,ncrv+k))= 1.0D0
     END DO

     RETURN
 END SUBROUTINE gbljac

 !> Get parameters and covariance matrix at point.
 !!
 !! \param [in]  IPOINT  (signed) point
 !!                      (<0: side towards previous point,
 !!                       >0: side towards next point)
 !! \param [out] DPAR    corrections (NP double precision values)
 !!                      (NP is number of track parameters: 5 + local par.)
 !! \param [out] DCOV    covariance matrix (NP2 double precision values,
 !!                      symmetric storage mode, NP2=(NP+1)*NP/2)

 SUBROUTINE gblres(ipoint,dpar,dcov)

     use gblmod


     INTEGER, INTENT(IN OUT)                  :: ipoint
     DOUBLE PRECISION, INTENT(OUT)            :: dpar(*)
     DOUBLE PRECISION, INTENT(OUT)            :: dcov(*)

     DOUBLE PRECISION :: caux(mxaux2),baux(mxaux)

     nb=ncrv+nloc
     nb2=(nb*nb+nb)/2
     np=nb+2*ndim
     mp=5+nloc
     mp2=(mp*mp+mp)/2

     jpnt=IABS(ipoint)
     IF (jpnt < 1.OR.jpnt > npnt) THEN
         IF (lpr > 0) PRINT *, ' GBLRES invalid point ', ipoint
         DO k=1,mp
             dpar(k)=0.0D0
         END DO
         DO k=1,mp2
             dcov(k)=0.0D0
         END DO
         RETURN
     END IF

     IF (nlvl < 4) RETURN ! fit not yet performed or failed

     CALL gbljac(ipoint,0,ioff1) ! get jacobian broken lines -> local
     !     get compressed covariance matrix, result vector
     DO i=1,nb  ! border
         baux(i)=bvec(i)
     END DO
     DO i=1,nb2 ! border
         caux(i)=cmat(i)
     END DO
     ip0=(ioff1-1)*ndim
     k=nb2
     DO i=nb+1,np
         io=i+ip0
         ko=(io*io-io)/2
         baux(i)=bvec(io) ! band part
         DO j=1,nb        ! mixed part
             k=k+1
             caux(k)=cmat(ko+j)
         END DO
         ko=ko+ip0
         DO j=nb+1,i      ! band part
             k=k+1
             caux(k)=cmat(ko+j)
         END DO
     END DO

     CALL dbgax(djac,baux,dpar,mp,np)
     CALL dbavat(caux,djac,dcov,np,mp)

     RETURN
 END SUBROUTINE gblres

 !> Add (inverse covariance matrix from) external seed.
 !!
 !! \param [in]  IPOINT  (signed) point
 !!                      (<0: side towards previous point,
 !!                       >0: side towards next point)
 !! \param [in]  DPRC    precision matrix (inverse covariance) from
 !!                      external seed (NP2 double precision values,
 !!                      symmetric storage mode, NP2=(NP+1)*NP/2,
 !!                      NP is number of track parameters: 5 + local par.)

 SUBROUTINE gbladx(ipoint,dprc)

     use gblmod


     INTEGER, INTENT(IN)                      :: ipoint
     DOUBLE PRECISION, INTENT(IN)             :: dprc(*)


     jpnt=IABS(ipoint)
     IF (jpnt < 1.OR.jpnt > npnt) THEN
         IF (lpr > 0) PRINT *, ' GBLADX invalid point ', ipoint
         RETURN
     END IF

     IF (nlvl >= 3) RETURN ! fit already prepared or performed

     mp=5+nloc
     mp2=(mp*mp+mp)/2
     ixsd=ipoint
     DO k=1,mp2
         dpsd(k)=dprc(k)
     END DO

     RETURN
 END SUBROUTINE gbladx

 !> Invert 1*1 or 2*2 matrix.
 !! Invert matrix A (dimension ND), B=A^-1

 SUBROUTINE gblinv(nd,a,b)

     INTEGER, INTENT(IN)                      :: nd
     DOUBLE PRECISION, INTENT(IN)             :: a(*)
     DOUBLE PRECISION, INTENT(OUT)            :: b(*)

     DOUBLE PRECISION :: DET

     IF (nd <= 1) THEN
         b(1)=1.0/a(1)
     ELSE
         det=a(1)*a(4)-a(2)*a(3)
         b(1)= a(4)/det
         b(2)=-a(2)/det
         b(3)=-a(3)/det
         b(4)= a(1)/det
     END IF
     RETURN
 END SUBROUTINE gblinv

 !> Matrix multiplication.
 !! Matrix multiplication C=A*B, A is N1*N1, B N1*N2, C=N1*N2

 SUBROUTINE gblmlt(n1,n2,a,b,c)


     INTEGER, INTENT(IN)                      :: n1
     INTEGER, INTENT(IN)                      :: n2
     DOUBLE PRECISION, INTENT(IN)             :: a(*)
     DOUBLE PRECISION, INTENT(IN)             :: b(*)
     DOUBLE PRECISION, INTENT(OUT)            :: c(*)

     DOUBLE PRECISION :: sum

     DO i=1,n1
         DO j=1,n2
             ioff=i
             joff=(j-1)*n1
             sum=0.0D0
             DO k=1,n1
                 sum=sum+a(ioff)*b(joff+k)
                 ioff=ioff+n1
             END DO
             c((j-1)*n1+i)=sum
         END DO
     END DO
     RETURN
 END SUBROUTINE gblmlt


 !> Calculate offset part of inverse jacobian.
 !!
 !! \param [in]  A (5*5) matrix A
 !! \param [out] B (5*5) matrix B with last 2 rows of inverse(A)
 !!
 SUBROUTINE gblinu(a,b)


     DOUBLE PRECISION, INTENT(IN)             :: a(5,5)
     DOUBLE PRECISION, INTENT(OUT)            :: b(5,5)
     DOUBLE PRECISION :: a33(3,3),a23(2,3),a22(2,2),det

     !     A33 = A(1..3,1..3)^-1 * det(A(1..3,1..3)
     a33(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
     a33(2,1) = a(2,3)*a(3,1)-a(2,1)*a(3,3)
     a33(3,1) = a(2,1)*a(3,2)-a(2,2)*a(3,1)
     a33(1,2) = a(1,3)*a(3,2)-a(1,2)*a(3,3)
     a33(2,2) = a(1,1)*a(3,3)-a(1,3)*a(3,1)
     a33(3,2) = a(1,2)*a(3,1)-a(1,1)*a(3,2)
     a33(1,3) = a(1,2)*a(2,3)-a(1,3)*a(2,2)
     a33(2,3) = a(1,3)*a(2,1)-a(1,1)*a(2,3)
     a33(3,3) = a(1,1)*a(2,2)-a(1,2)*a(2,1)
     det=a(1,1)*a33(1,1)+a(1,2)*a33(2,1)+a(1,3)*a33(3,1)
     !     A23 = A(4..5,1..3) * A33 / det(A(1..3,1..3))
     a23(1,1) = (a(4,1)*a33(1,1)+a(4,2)*a33(2,1)+a(4,3)*a33(3,1))/det
     a23(1,2) = (a(4,1)*a33(1,2)+a(4,2)*a33(2,2)+a(4,3)*a33(3,2))/det
     a23(1,3) = (a(4,1)*a33(1,3)+a(4,2)*a33(2,3)+a(4,3)*a33(3,3))/det
     a23(2,1) = (a(5,1)*a33(1,1)+a(5,2)*a33(2,1)+a(5,3)*a33(3,1))/det
     a23(2,2) = (a(5,1)*a33(1,2)+a(5,2)*a33(2,2)+a(5,3)*a33(3,2))/det
     a23(2,3) = (a(5,1)*a33(1,3)+a(5,2)*a33(2,3)+a(5,3)*a33(3,3))/det
     !     A22 = A(4..5,4..5) - A23 * A((1..3,4..5)
     a22(1,1) = a(4,4)-a23(1,1)*a(1,4)-a23(1,2)*a(2,4)-a23(1,3)*a(3,4)
     a22(1,2) = a(4,5)-a23(1,1)*a(1,5)-a23(1,2)*a(2,5)-a23(1,3)*a(3,5)
     a22(2,1) = a(5,4)-a23(2,1)*a(1,4)-a23(2,2)*a(2,4)-a23(2,3)*a(3,4)
     a22(2,2) = a(5,5)-a23(2,1)*a(1,5)-a23(2,2)*a(2,5)-a23(2,3)*a(3,5)
     !     (4..5,4..5) of A^-1 = A22^-1
     det = a22(1,1)*a22(2,2)-a22(1,2)*a22(2,1)
     b(4,4) = a22(2,2)/det
     b(4,5) =-a22(1,2)/det
     b(5,4) =-a22(2,1)/det
     b(5,5) = a22(1,1)/det
     !     (4..5,1..3) of A^-1 = -A22^-1 * A23
     b(4,1) = -b(4,4)*a23(1,1)-b(4,5)*a23(2,1)
     b(4,2) = -b(4,4)*a23(1,2)-b(4,5)*a23(2,2)
     b(4,3) = -b(4,4)*a23(1,3)-b(4,5)*a23(2,3)
     b(5,1) = -b(5,4)*a23(1,1)-b(5,5)*a23(2,1)
     b(5,2) = -b(5,4)*a23(1,2)-b(5,5)*a23(2,2)
     b(5,3) = -b(5,4)*a23(1,3)-b(5,5)*a23(2,3)

     RETURN
 END SUBROUTINE gblinu

