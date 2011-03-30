*
* C. Kleinwort , DESY-FH1, www.terascale.de
* February 2011, Analysis Centre: Statistics Tools Group
*
***** Collection of utility routines from V. Blobel ***** 
*
*======================================================================
*
*     V. Blobel, Univ. Hamburg
*     Numerical subprograms used in MP-II: matrix equations, 
*        and matrix products, double precision
*
*     Solution by inversion
*        SQMINV
*
*     Solution by Cholesky decomposition of bordered band matrix
CHK      SQMIBB
*
*     Matrix/vector products
*-        DBDOT     dot vector product
*-        DBAXPY    multiplication and addition
*-        DBSVX     symmetric matrix vector
*        DBGAX     general matrix vector
*        DBAVAT    AVAT product
*-        DBMPRV    print parameter and matrix 
CHK      DBPRV     print matrix 
*
*     matrix inversion -------------------------------------------------

      SUBROUTINE SQMINV(V,B,N,NRANK,DIAG,NEXT)   ! matrix inversion  
*     obtain solution of a system of linear equations with symmetric
*     matrix and the inverse.
*
*                    - - -
*        CALL SQMINV(V,B,N,NRANK,...,...)      solve  V * X = B
*                    - -   -----
*
*           V = symmetric N-by-N matrix in symmetric storage mode
*               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
*               replaced by inverse matrix
*           B = N-vector, replaced by solution vector
*
*     DIAG(N) =  double precision scratch array
*     NEXT(N) =  integer aux array
*
*     Method of solution is by elimination selecting the  pivot  on  the
*     diagonal each stage. The rank of the matrix is returned in  NRANK.
*     For NRANK ne N, all remaining  rows  and  cols  of  the  resulting
*     matrix V and the corresponding elements of  B  are  set  to  zero.
*
      DOUBLE PRECISION V(*),B(N),DIAG(N),VKK,VJK,EPS
      INTEGER        NEXT(N)
      PARAMETER     (EPS=1.0D-10)
*     ...
      NEXT0=1
      DO I=1,N
       NEXT(I)=I+1                ! set "next" pointer
       DIAG(I)=ABS(V((I*I+I)/2))  ! save abs of diagonal elements
      END DO
      NEXT(N)=-1                  ! end flag

      NRANK=0
      DO I=1,N                    ! start of loop
       K  =0
       VKK=0.0D0

       J=NEXT0
       LAST=0 
 05    IF(J.GT.0) THEN
          JJ=(J*J+J)/2
          IF(ABS(V(JJ)).GT.MAX(ABS(VKK),EPS*DIAG(J))) THEN
             VKK=V(JJ)              
             K=J
             L=LAST
          END IF
          LAST=J
          J=NEXT(LAST)
          GOTO 05   
       END IF 

       IF(K.NE.0) THEN            ! pivot found
          KK=(K*K+K)/2
          IF(L.EQ.0) THEN
             NEXT0=NEXT(K)
          ELSE 
             NEXT(L)=NEXT(K)
          END IF
          NEXT(K)=0               ! index is used, reset
          NRANK=NRANK+1           ! increase rank and ...
          VKK    =1.0/VKK
          V(KK)  =-VKK
          B(K)   =B(K)*VKK
          JK     =KK-K
          JL     =0
          DO J=1,N                ! elimination
           IF(J.EQ.K) THEN
              JK=KK
              JL=JL+J
           ELSE
              IF(J.LT.K) THEN
                 JK=JK+1
              ELSE
                 JK=JK+J-1
              END IF
              VJK  =V(JK)
              V(JK)=VKK*VJK
              B(J) =B(J)-B(K)*VJK
              LK   =KK-K
              DO L=1,J
               JL=JL+1
               IF(L.EQ.K) THEN
                  LK=KK
               ELSE
                  IF(L.LT.K) THEN
                     LK=LK+1
                  ELSE
                     LK=LK+L-1
                  END IF
                  V(JL)=V(JL)-V(LK)*VJK
               END IF
              END DO
           END IF
          END DO
       ELSE
          DO K=1,N
           IF(NEXT(K).NE.0) THEN
              B(K)=0.0D0       ! clear vector element
              DO J=1,K
               IF(NEXT(J).NE.0) V((K*K-K)/2+J)=0.0D0  ! clear matrix row/col
              END DO
           END IF
          END DO
          GOTO 10
       END IF
      END DO             ! end of loop
 10   DO IJ=1,(N*N+N)/2
       V(IJ)=-V(IJ)      ! finally reverse sign of all matrix elements
      END DO
      END

      SUBROUTINE DBGAX(A,X,Y,M,N)                ! gen. matrix, vector
*     multiply general M-by-N matrix A and N-vector X 
*
*                    - -   - -
*        CALL  DBGAX(A,X,Y,M,N)          Y   :=   A   *    X
*                        -               M       M*N       N
*
*        where A = general M-by-N matrix (A11 A12 ... A1N  A21 A22 ...) 
*              X = N vector 
*              Y = M vector
*
      DOUBLE PRECISION A(*),X(*),Y(*)
*     ...
      IJ=0
      DO I=1,M
       Y(I)=0.0D0
       DO J=1,N
        IJ=IJ+1
        Y(I)=Y(I)+A(IJ)*X(J)
       END DO  
      END DO 
      END

      SUBROUTINE DBAVAT(V,A,W,N,M)               ! A V AT product    
*     multiply symmetric N-by-N matrix from the left with general M-by-N
*     matrix and from the right with the transposed of the same  general
*     matrix  to  form  symmetric  M-by-M   matrix   (used   for   error
*     propagation).
*
*                    - -   - -                                   T
*        CALL DBAVAT(V,A,W,N,M)         W   =   A   *   V   *   A
*                        -             M*M     M*N     N*N     N*M
*
*        where V = symmetric N-by-N matrix
*              A = general N-by-M matrix
*              W = symmetric M-by-M matrix
*
      DOUBLE PRECISION V,A,W,CIK 
      DIMENSION V(*),A(*),W(*)
*     ...
      DO I=1,(M*M+M)/2
       W(I)=0.0                ! reset output matrix
      END DO
      IL=-N
      IJS=0
      DO I=1,M                 ! do I
       IJS=IJS+I-1             !
       IL=IL+N                 !
       LKL=0                   !
       DO K=1,N                !   do K
        CIK=0.0D0              !  
        LKL=LKL+K-1            ! 
        LK=LKL                 !
        DO L=1,K               !     do L 
         LK=LK+1               !     .
         CIK=CIK+A(IL+L)*V(LK) !     .
        END DO                 !     end do L
        DO L=K+1,N             !     do L 
         LK=LK+L-1             !     . 
         CIK=CIK+A(IL+L)*V(LK) !     .
        END DO                 !     end do L 
        JK=K                   !
        IJ=IJS                 !
        DO J=1,I               !     do J 
         IJ=IJ+1               !     .
         W(IJ)=W(IJ)+CIK*A(JK) !     .
         JK=JK+N               !     . 
        END DO                 !     end do J
       END DO                  !   end do K 
      END DO                   ! end do I 
      END

*                                                 090817 C. Kleinwort, DESY-FH1
*     --------- bordered band matrix -----------

      SUBROUTINE SQMIBB(V,B,N,NBDR,NBND,INV,NRANK)   
*     obtain solution of a system of linear equations with symmetric
*     bordered band matrix 
*     on request inverse is calculated  
*
*                    - - - ---- ---- ---
*        CALL SQMIBB(V,B,N,NBND,NBDR,INV,NRANK,...,...)      solve  V * X = B
*                    - -                 -----
*
*         V    = symmetric N-by-N matrix in symmetric storage mode
*                V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
*                optionally replaced by inverse matrix
*         B    = N-vector, replaced by solution vector
*
*         NBRD = border size
*         NBND = band width
*         INV =1 calculate band part of inverse (for pulls)
*             >1 calculate complete inverse
*
*
* Double precision scratch arrays:
*     VBND(N*(NBND+1)) = storage of band   part
*     VBDR(N* NBDR)    = storage of border part 
*     AUX (N* NBDR)    = intermediate results
*
* cost[dot ops] ~= (N-NBDR)*(NBDR+NBND+1)**2 + NBDR**3/3 (leading term, solution only)
*
      DOUBLE PRECISION V(*),B(N)
*
#include "gbldef.inc"
      DOUBLE PRECISION VBND(MXPAR*(MXBND+1)),VBDR(MXPAR*MXBDR),
     +   AUX(MXPAR*MXBDR),VBK((MXBDR*MXBDR+MXBDR)/2),
     +   VZRU(MXBDR),SCDIAG(MXBDR)
      INTEGER SCFLAG(MXBDR)
*          
      SAVE NPRI
      DATA NPRI / 100 /
*           ...
      NRANK=0
      NB1=NBDR+1      
      MP1=NBND+1
      NMB=N-NBDR
*     copy band part
      DO I=NB1,N
         IP=(I*(I+1))/2
         IS=0
         DO J=I,MIN(N,I+NBND)
            IP=IP+IS
            IS=J
            IB=J-I+1
            VBND(IB+(I-NB1)*MP1)=V(IP)
         ENDDO
      ENDDO
*     copy border part
      IF (NBDR.GT.0) THEN
         IOFF=0
         DO I=1,NBDR
            IP=(I*(I+1))/2
            IS=0
            DO J=I,N
               IP=IP+IS
               IS=J
               VBDR(IOFF+J)=V(IP)
            ENDDO
            IOFF=IOFF+N
         ENDDO
      ENDIF 
*
      CALL DBCDEC(VBND,MP1,NMB,AUX)
* use? CALL DBFDEC(VBND,MP1,NMB) ! modified decomp., numerically more stable
C      CALL DBCPRB(VBND,MP1,NMB)
      IP=1
      DO I=1, NMB
         IF (VBND(IP).LE.0.0D0) THEN
            NPRI=NPRI-1
            IF (NPRI.GE.0) THEN
               IF (VBND(IP).EQ.0.0D0) THEN
                  print *, ' SQMIBB matrix singular', N, NBDR, NBND   
               ELSE
                  print *, ' SQMIBB matrix not positive definite', 
     +                  N, NBDR, NBND
               ENDIF 
            ENDIF     
*           return zeros        
            DO IP=1,N
               B(IP)=0.0D0
            ENDDO   
            DO IP=1,(N*N+N)/2
               V(IP)=0.0D0
            ENDDO   
            RETURN
         ENDIF
         IP=IP+MP1
      ENDDO
      NRANK=NMB              
*
      IF (NBDR.EQ.0) THEN ! special case NBDR=0
*      
         CALL DBCSLV(VBND,MP1,NMB,B,B) 
         IF (INV.GT.0) THEN
            IF (INV.GT.1) THEN
               CALL DBCINV(VBND,MP1,NMB,V)
            ELSE
               CALL DBCINB(VBND,MP1,NMB,V)
            ENDIF
         ENDIF        
*       
      ELSE ! general case NBDR>0
*        
         IOFF=NB1 
         DO IB=1,NBDR
*           solve for aux. vectors
            CALL DBCSLV(VBND,MP1,NMB,VBDR(IOFF),AUX(IOFF)) 
*           zT ru
            VZRU(IB)=B(IB)
            DO I=0,NMB-1
               VZRU(IB)=VZRU(IB)-B(NB1+I)*AUX(IOFF+I)
            END DO
            IOFF=IOFF+N
         ENDDO
*        solve for band part only   
         CALL DBCSLV(VBND,MP1,NMB,B(NB1),B(NB1))
*        Ck - cT z
         IP=0
         IOFF=NB1 
         DO IB=1,NBDR
            JOFF=NB1 
            DO JB=1,IB
               IP=IP+1
               VBK(IP)=V(IP)
               DO I=0,NMB-1
                  VBK(IP)=VBK(IP)-VBDR(IOFF+I)*AUX(JOFF+I)
               ENDDO
               JOFF=JOFF+N
            ENDDO
            IOFF=IOFF+N                      
         ENDDO
*        solve border part
         CALL SQMINV(VBK,VZRU,NBDR,NRANKB,SCDIAG,SCFLAG)        
         IF (NRANKB.EQ.NBDR) THEN
            NRANK=NRANK+NBDR
         ELSE
            NPRI=NPRI-1
            IF (NPRI.GE.0) 
     +      print *, ' SQMIBB undef border ', N, NBDR, NBND, NRANKB   
            DO IB=1,NBDR
               VZRU(IB)=0.0D0
            ENDDO   
            DO IP=(NBDR*NBDR+NBDR)/2,1,-1
               VBK(IP)=0.0D0
            ENDDO   
         ENDIF   
*        smoothed data points
         IOFF=NB1
         DO IB=1, NBDR
            B(IB) = VZRU(IB)
            DO I=0,NMB-1
               B(NB1+I)=B(NB1+I)-B(IB)*AUX(IOFF+I)
            ENDDO
            IOFF=IOFF+N
         ENDDO         
*        inverse requested ?        
         IF (INV.GT.0) THEN
            IF (INV.GT.1) THEN
               CALL DBCINV(VBND,MP1,NMB,V)
            ELSE
               CALL DBCINB(VBND,MP1,NMB,V)     
            ENDIF        
*           expand/correct from NMB to N      
            IP1=(NMB*NMB+NMB)/2
            IP2=(N*N+N)/2
            DO I=NMB-1,0,-1
               J0=0
               IF (INV.EQ.1) J0=MAX(0,I-NBND)
               DO J=I,J0,-1
                  V(IP2)=V(IP1)
                  IOFF=NB1
                  DO IB=1,NBDR
                     JOFF=NB1
                     DO JB=1,NBDR
                        IJ=MAX(IB,JB)
                        IJ=(IJ*IJ-IJ)/2+MIN(IB,JB)
                        V(IP2)=V(IP2)+VBK(IJ)*AUX(IOFF+I)*AUX(JOFF+J)
                        JOFF=JOFF+N
                     ENDDO
                     IOFF=IOFF+N
                  ENDDO
                  IP1=IP1-1
                  IP2=IP2-1
               ENDDO
               IP1=IP1-J0
               IP2=IP2-J0
*
               DO IB=NBDR,1,-1
                  V(IP2)=0.0D0
                  JOFF=NB1
                  DO JB=1,NBDR
                     IJ=MAX(IB,JB)
                     IJ=(IJ*IJ-IJ)/2+MIN(IB,JB)
                     V(IP2)=V(IP2)-VBK(IJ)*AUX(I+JOFF)
                     JOFF=JOFF+N
                  ENDDO
                  IP2=IP2-1
               ENDDO
            ENDDO
*          
            DO IP=(NBDR*NBDR+NBDR)/2,1,-1
               V(IP2)=VBK(IP)
               IP2=IP2-1 
            ENDDO
*
         ENDIF                     
      ENDIF     
*            
      END
      SUBROUTINE DBPRV(LUN,V,N)
*     Prints the symmetric N-by-N matrix V
*
      DOUBLE PRECISION V
      DIMENSION V(*)
      PARAMETER (ISTP=6) 
*
      WRITE(LUN,101)
*
      DO I=1,N
         IPS=(I*I-I)/2
        IPE=IPS+I
        IP =IPS
  100   CONTINUE
        IPN=IP+ISTP
        WRITE(LUN,102), I, IP+1-IPS, (V(K),K=IP+1,MIN(IPN,IPE))
        IF (IPN.LT.IPE) THEN
           IP=IPN
           GOTO 100
        ENDIF   
      ENDDO    
      RETURN
  101 FORMAT(1X,'--- DBPRV -----------------------------------')
  102 FORMAT(1X,2I3,6G12.4)
*
      END
