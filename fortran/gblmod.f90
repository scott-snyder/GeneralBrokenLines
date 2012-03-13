!> \file
!! Trajectory data.

!> Definitions and data fields.
module gblmod
    IMPLICIT NONE
    SAVE
    INTEGER, PARAMETER  :: mxpar=250  !< max. number of fit parameters
    INTEGER, PARAMETER  :: mxpar2=(mxpar*mxpar+mxpar)/2
    INTEGER, PARAMETER  :: mxbnd=10   !< max. band size
    INTEGER, PARAMETER  :: mxbdr=10   !< max. border size
    INTEGER, PARAMETER  :: mxaux=mxbdr+5  !< max. number of local track parameters
    INTEGER, PARAMETER  :: mxaux2=(mxaux*mxaux+mxaux)/2
    INTEGER, PARAMETER  :: mxpnt=1000 !< max. points
    INTEGER, PARAMETER  :: mxdat=10*mxpnt !< max. integrated size of data blocks
!
    INTEGER :: nlvl   !< level of (preparation of) trajectory (fit)
    INTEGER :: lpr    !< print level
    INTEGER :: ndim   !< dimension of trajectory
                      !! (2: 2D offsets for track in 3D, 1: 1D offsets for track in 2D)
    INTEGER :: npnt   !< number of points
    INTEGER :: nmeas  !< number (of points) with measurement
    INTEGER :: nscat  !< number (of points) with scatterer
    INTEGER :: noff   !< number (of points) with offsets (as fit parameter)
    INTEGER :: nknk   !< number of (multiple scattering) kinks
    INTEGER :: lpnt   !< last point with offset
    INTEGER :: lscat  !< last point with scatterer
    INTEGER :: npar   !< number of fit parameters
    INTEGER :: ncrv   !< flag (0/1) for usage of Q/P as fit parameter
    INTEGER :: nbnd   !< band size
    INTEGER :: nloc   !< number of additional local track parameters
    INTEGER :: ixsd   !< point with external seed (or 0)
    INTEGER :: icrd0  !< first active offset coordinate (0: u_1, 1: u_2)
    INTEGER :: ldwm   !< last used down-weighting method
    DOUBLE PRECISION, DIMENSION(mxpar)  :: bvec        !< right hand side of linear equation system (A*x=b)
    DOUBLE PRECISION, DIMENSION(mxpar2) :: cmat        !< (sym.) matrix of linear equation system (A*x=b)
    DOUBLE PRECISION, DIMENSION(mxaux2) :: dpsd        !< external seed (precision matrix)
    DOUBLE PRECISION, DIMENSION(mxaux*mxaux) :: djac   !< jacobian for transformation fit to track parameter
    DOUBLE PRECISION, DIMENSION(25,mxpnt) :: djp2p     !< point-to-point jacobian (from previous point)
    DOUBLE PRECISION, DIMENSION(8,mxpnt)  :: jmat      !< matrix 'J', offset part of jacobian to next/previous point with offset
    DOUBLE PRECISION, DIMENSION(8,mxpnt)  :: smat      !< matrix 'S', slope part of jacobian to next/previous point with offset
    DOUBLE PRECISION, DIMENSION(4,mxpnt)  :: dvec      !< vector 'd', Q/P part of jacobian to next/previous point with offset
    DOUBLE PRECISION, DIMENSION(4,mxpnt)  :: pmat      !< projection matrix from measurement to local coordinate system
    INTEGER, DIMENSION(mxpnt) :: imeas    !< measurement index (for point)
    INTEGER, DIMENSION(mxpnt) :: iscat    !< scatterer index (for point)
    INTEGER, DIMENSION(mxpnt) :: ioff     !< offset at point (or from previous point)
    INTEGER, DIMENSION(mxpnt) :: iknk     !< kink index (for point)
    INTEGER, DIMENSION(2,mxpnt) :: idlc   !< index of local derivatives (for point)
    INTEGER, DIMENSION(2,mxpnt) :: idgl   !< index of global derivatives (for point)
    REAL, DIMENSION(2,mxpnt) :: yvec      !< residual (for measurement or kink)
    REAL, DIMENSION(2,mxpnt) :: wvec      !< precision (for measurement or kink, 1/sigma**2)
    INTEGER :: ndat  !< number of data blocks with measurements or kinks
    INTEGER :: ldat  !< integrated size of data blocks
    INTEGER :: mdat  !< max. integrated size
    INTEGER, DIMENSION(mxdat) :: idat  !< integer part of data blocks (lower part, 1..mdat)
                                       !! or local/global derivatives (upper part, mdat+1..mxdat)
    REAL, DIMENSION(mxdat)    :: adat  !< float part of data blocks or local/global derivatives
!
end module gblmod
