
    SUBROUTINE PPFBET (PR,P,Q,TOL,IFLAG,X)
    !
    !-----------------------------------------------------------------------
    !   PPFBET   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
    !            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
    !            MARYLAND 20899
    !
    !   FOR: EVALUATING THE INVERSE CUMULATIVE DISTRIBUTION FUNCTION
    !        (PERCENT POINT FUNCTION) OF THE BETA(P,Q) DISTRIBUTION.
    !        FOR A GIVEN PROBABILITY PR.  THE PERCENT POINT X IS COMPUTED
    !        TO A SPECIFIED ABSOLUTE ACCURACY TOL WHEN POSSIBLE. THE
    !        METHOD OF BRENT, AS DESCRIBED IN THE REFERENCES BELOW, IS
    !        USED TO COMPUTE THE APPROXIMATE ZERO OF I(X,P,Q)-PR WHERE
    !        I(X,P,Q) IS THE CUMULATIVE DISTRIBUTION FUNCTION OF THE
    !        BETA(P,Q) DISTRIBUTION EVALUATED AT X.  THIS METHOD DOES
    !        NOT REQUIRE DERIVATIVES.
    !
    !   NOTE: THE CONSTANT EPS FOR MACHINE FLOATING POINT PRECISION IS
    !         MACHINE DEPENDENT.
    !
    !   SUBPROGRAMS CALLED: CDFBET (BETA CUMULATIVE DISTRIBUTION FUNCTION)
    !
    !   CURRENT VERSION COMPLETED OCTOBER 13, 1987
    !
    !   REFERENCES:
    !
    !   1) PRESS, WILLIAM H., FLANNERY, BRIAN P., TEUKOLSKY, SAUL A.,
    !      AND VETTERLING, WILLIAM T., 'NUMERICAL RECIPES - THE ART OF
    !      SCIENTIFIC COMPUTING', CAMBRIDGE UNIVERSITY PRESS, 1986,
    !      PP. 251-254.
    !
    !   2) BRENT, RICHARD P., 'ALGORITHMS FOR MINIMIZATION WITHOUT
    !      DERIVATIVES', PRENTICE-HALL, 1973, CH. 3-4.
    !-----------------------------------------------------------------------
    !   DEFINITION OF PASSED PARAMETERS:
    !
    !    * PR = A PROBABILITY VALUE IN THE INTERVAL [0,1] (REAL)
    !
    !     * P = THE FIRST PARAMETER (>0) OF THE BETA(P,Q) DISTRIBUTION
    !           (REAL)
    !
    !     * Q = THE SECOND PARAMETER (>0) OF THE BETA(P,Q) DISTRIBUTION
    !           (REAL)
    !
    !   * TOL = THE REQUIRED ACCURACY (>=1.0E-8) OF THE PERCENT POINT X
    !           (REAL)
    !
    !   IFLAG = ERROR INDICATOR ON OUTPUT (INTEGER)   INTERPRETATION:
    !             0 -> NO ERRORS DETECTED
    !           1,2 -> ERROR FLAGS FROM SUBROUTINE CDFBET
    !             3 -> PR<0 OR PR>1
    !             4 -> P<=0 OR Q<=0
    !             5 -> TOL<1.0E-8
    !             6 -> THE CDF'S AT THE ENDPOINTS HAVE THE SAME SIGN - NO
    !                  VALUE OF X IS DEFINED (THIS SHOULD NEVER OCCUR)
    !             7 -> MAXIMUM ITERATIONS EXCEEDED - CURRENT VALUE OF X
    !                  RETURNED
    !
    !       X = THE COMPUTED PERCENT POINT (REAL)
    !
    !   * INDICATES PARAMETERS REQUIRING INPUT VALUES
    !-----------------------------------------------------------------------
    !
    !--- DEFINE MAXIMUM ITERATIONS AND MACHINE FLOATING POINT PRECISION
    !
    DATA ITMAX,EPS / 50,1.0E-12 /
    !
    !--- CHECK VALIDITY OF INPUT ARGUMENTS
    !
    IFLAG = 0
    IF (PR.LT.0.0.OR.PR.GT.1.0) THEN
        IFLAG = 3
        RETURN
        !
    ENDIF
    IF (AMIN1(P,Q).LE.0.0) THEN
        IFLAG = 4
        RETURN
        !
    ENDIF
    IF (TOL.LT.1.0E-8) THEN
        IFLAG = 5
        RETURN
        !
    ENDIF
    A = 0.0
    B = 1.0
    FA = -PR
    FB = 1.0-PR
    !
    !--- CHECK FOR ZERO BEING BRACKETED
    !
    IF (FB*FA.GT.0.0) THEN
        IFLAG = 6
        RETURN
        !
    ENDIF
    FC = FB
    DO 10 ITER = 1, ITMAX
        IF (FB*FC.GT.0.0) THEN
            !
            !--- RENAME A,B,C AND ADJUST BOUNDING INTERVAL D
            !
            C = A
            FC = FA
            D = B-A
            E = D
        ENDIF
        IF (ABS(FC).LT.ABS(FB)) THEN
            A = B
            B = C
            C = A
            FA = FB
            FB = FC
            FC = FA
        ENDIF
        !
        !--- CONVERGENCE CHECK
        !
        TOL1 = 2.0*EPS*ABS(B)+0.5*TOL
        XM = 0.5*(C-B)
        IF (ABS(XM).LE.TOL1.OR.FB.EQ.0.0) THEN
            X = B
            RETURN
            !
        ENDIF
        IF (ABS(E).GE.TOL1.AND.ABS(FA).GT.ABS(FB)) THEN
            !
            !--- ATTEMPT INVERSE QUADRATIC INTERPOLATION
            !
            S = FB/FA
            IF (A.EQ.C) THEN
                U = 2.0*XM*S
                V = 1.0-S
            ELSE
                V = FA/FC
                R = FB/FC
                U = S*(2.0*XM*V*(V-R)-(B-A)*(R-1.0))
                V = (V-1.0)*(R-1.0)*(S-1.0)
            ENDIF
            !
            !--- CHECK WHETHER IN BOUNDS
            !
            IF (U.GT.0.0) V = -V
            U = ABS(U)
            IF (2.0*U.LT.AMIN1(3.0*XM*V-ABS(TOL1*V),ABS(E*V))) THEN
                !
                !--- ACCEPT INTERPLATION
                !
                E = D
                D = U/V
            ELSE
                !
                !--- INTERPOLATION FAILED, USE BISECTION
                !
                D = XM
                E = D
            ENDIF
        ELSE
            !
            !--- BOUNDS DECREASING TOO SLOWLY, USE BISECTION
            !
            D = XM
            E = D
        ENDIF
        !
        !--- MOVE LAST BEST GUESS TO A
        !
        A = B
        FA = FB
        !
        !--- EVALUATE NEW TRIAL ZERO
        !
        IF (ABS(D).GT.TOL1) THEN
            B = B+D
        ELSE
            B = B+SIGN(TOL1,XM)
        ENDIF
        CALL CDFBET (B,P,Q,EPS,IFLAG,CDF)
        IF (IFLAG.NE.0) RETURN
        FB = CDF-PR
10  CONTINUE
    IFLAG = 7
    X = B
    RETURN
    !
    END


    SUBROUTINE CDFBET (X,P,Q,EPS,IFLAG,CDFX)
    !
    !-----------------------------------------------------------------------
    !   CDFBET   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
    !            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
    !            MARYLAND  20899
    !
    !   FOR: COMPUTING THE CUMULATIVE DISTRIBUTION FUNCTION OF THE BETA
    !        DISTRIBUTION (ALSO KNOWN AS THE INCOMPLETE BETA RATIO) TO A
    !        SPECIFIED ACCURACY (TRUNCATION ERROR IN THE INFINITE SERIES).
    !        THE ALGORITHM, DESCRIBED IN REFERENCE 2, IS A MODIFICATION OF
    !        THE ALGORITHM OF REFERENCE 1.  THREE FEATURES HAVE BEEN ADDED:
    !
    !        1) A PRECISE METHOD OF MEETING THE TRUNCATION ACCURACY,
    !        2) A CONSTANT W USED IN DETERMINING FOR WHICH X VALUES THE
    !           RELATION I(X,P,Q) = 1 - I(1-X,Q,P) IS TO BE USED, AND
    !        3) A CONSTANT UFLO >= THE UNDERFLOW LIMIT ON THE COMPUTER.
    !
    !   SUBPROGRAMS CALLED: DGAMLN (LOG OF GAMMA FUNCTION)
    !
    !   CURRENT VERSION COMPLETED OCTOBER 24, 1986
    !
    !   REFERENCES:
    !
    !   1) MAJUMDER, K.L. AND BHATTACHARJEE, G.P., 'THE INCOMPLETE BETA
    !      INTEGRAL', ALGORITHM AS 63, APPLIED STATISTICS, VOL. 22, NO. 3,
    !      1973, PP. 409-411.
    !
    !   2) REEVE, CHARLES P., 'AN ALGORITHM FOR COMPUTING THE BETA !.D.F.
    !      TO A SPECIFIED ACCURACY', STATISTICAL ENGINEERING DIVISION
    !      NOTE 86-3, OCTOBER 1986.
    !-----------------------------------------------------------------------
    !   DEFINITION OF PASSED PARAMETERS:
    !
    !      * X = VALUE AT WHICH THE !.D.F. IS TO BE COMPUTED (REAL)
    !
    !      * P = FIRST PARAMETER OF THE BETA FUNCTION (>0) (REAL)
    !
    !      * Q = SECOND PARAMETER OF THE BETA FUNCTION (>0) (REAL)
    !
    !    * EPS =  THE DESIRED ABSOLUTE ACCURACY OF THE !.D.F. (>0) (REAL)
    !
    !    IFLAG = ERROR INDICATOR ON OUTPUT (INTEGER)   INTERPRETATION:
    !            0 -> NO ERRORS DETECTED
    !            1 -> EITHER P OR Q OR EPS IS <= UFLO
    !            2 -> NUMBER OF TERMS EVALUATED IN THE INFINITE SERIES
    !                 EXCEEDS JMAX
    !
    !     CDFX = THE !.D.F. EVALUATED AT X (REAL)
    !
    !   * INDICATES PARAMETERS REQUIRING INPUT VALUES
    !-----------------------------------------------------------------------
    !
    LOGICAL LL
    DOUBLE PRECISION DP,DQ,DGAMLN
    !CCC DATA JMAX,W,UFLO / 5000,20.0,1E-100 /
    DATA JMAX,W,UFLO / 5000,20.0,1E-30 /
    CDFX = 0.0
    !
    !--- CHECK FOR VALIDITY OF ARGUMENTS P, Q, AND EPS
    !
    IF (P.LE.UFLO.OR.Q.LE.UFLO.OR.EPS.LE.UFLO) THEN
        IFLAG = 1
        RETURN
    ENDIF
    IFLAG = 0
    !
    !--- CHECK FOR SPECIAL CASES OF X
    !
    IF (X.LE.0.0) RETURN
    IF (X.GE.1.0) THEN
        CDFX = 1.0
    ELSE
        !
        !--- SWITCH ARGUMENTS IF NECESSARY
        !
        LL = P+W.GE.(P+Q+2.0*W)*X
        IF (LL) THEN
            XY = X
            YX = 1.0-XY
            PQ = P
            QP = Q
        ELSE
            YX = X
            XY = 1.0-YX
            QP = P
            PQ = Q
        ENDIF
        !
        !--- EVALUATE THE BETA P.D.F. AND CHECK FOR UNDERFLOW
        !
        DP = DBLE(PQ-1.0)*DLOG(DBLE(XY))-DGAMLN(PQ)
        DQ = DBLE(QP-1.0)*DLOG(DBLE(YX))-DGAMLN(QP)
        PDFL = SNGL(DGAMLN(PQ+QP)+DP+DQ)
        IF (PDFL.LT.ALOG(UFLO)) THEN
        ELSE
            U = EXP(PDFL)*XY/PQ
            R = XY/YX
10          IF (QP.LE.1.0) GO TO 20
            !
            !--- INCREMENT PQ AND DECREMENT QP
            !
            IF (U.LE.EPS*(1.0-(PQ+QP)*XY/(PQ+1.0))) GO TO 40
            CDFX = CDFX+U
            PQ = PQ+1.0
            QP = QP-1.0
            U = QP*R*U/PQ
            GO TO 10
20          V = YX*U
            YXEPS = YX*EPS
            !
            !--- INCREMENT PQ
            !
            DO 30 J = 0, JMAX
                IF (V.LE.YXEPS) GO TO 40
                CDFX = CDFX+V
                PQ = PQ+1.0
                V = (PQ+QP-1.0)*XY*V/PQ
30          CONTINUE
            IFLAG = 2
        ENDIF
40      IF (.NOT.LL) CDFX = 1.0-CDFX
    ENDIF
    RETURN
    END

    DOUBLE PRECISION FUNCTION DGAMLN (X)
    !
    !-----------------------------------------------------------------------
    !   DGAMLN   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
    !            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
    !            MARYLAND  20899
    !
    !   FOR: COMPUTING THE DOUBLE PRECISION LOG OF THE GAMMA FUNCTION WITH
    !        SINGLE PRECISION PARAMETER X>0.  THE MAXIMUM TRUNCATION ERROR
    !        IN THE INFINITE SERIES (SEE REFERENCE 1) IS DETERMINED BY THE
    !        CONSTANT XMIN.  WHEN X<XMIN A RECURRENCE RELATION IS USED IN
    !        ORDER TO ACHIEVE THE REQUIRED ABSOLUTE ACCURACY.  THE TABLE
    !        BELOW GIVES THE MINIMUM VALUE OF X WHICH YIELDS THE CORRES-
    !        PONDING ABSOLUTE ACCURACY IN DGAMLN(X) ASSUMING THE MACHINE
    !        CARRIES ENOUGH DIGITS WHEN THOSE TO THE LEFT OF THE DECIMAL
    !        ARE CONSIDERED (SEE REFERENCE 2 FOR FURTHER DISCUSSION).  IF
    !        THE LATTER CONDITION IS NOT MET, AN ERROR MESSAGE IS PRINTED.
    !
    !        THE CYBER 180/855 AT NBS CARRIES ABOUT 15 DIGITS IN SINGLE
    !        PRECISION, THEREFORE THE PRE-SET VALUE OF ABSACC IS 10**(-15)
    !        AND THE CORRESPONDING VALUE OF XMIN IS 6.894.  ON A DIFFERENT
    !        MACHINE THESE CONSTANTS SHOULD BE CHANGED ACCORDINGLY.
    !
    !           XMIN    ACCURACY    XMIN    ACCURACY    XMIN    ACCURACY
    !          ------   --------   ------   --------   ------   --------
    !           1.357     1E-3      4.592     1E-12    15.539     1E-21
    !           2.037     1E-6      6.894     1E-15    23.330     1E-24
    !           3.059     1E-9     10.351     1E-18    35.025     1E-27
    !
    !   NOTE: THIS IS EXACTLY THE SAME SOFTWARE AS SUBROUTINE DGAMLN
    !
    !   SUBPROGRAMS CALLED: -NONE-
    !
    !   CURRENT VERSION COMPLETED MAY 1, 1989
    !
    !   REFERENCES:
    !
    !   1) ABRAMOWITZ, MILTON AND STEGUN, IRENE, 'HANDBOOK OF MATHEMATICAL
    !      FUNCTIONS', NBS APPLIED MATHEMATICS SERIES 55, NOV. 1970,
    !      EQ. 6.1.40, P 257.
    !
    !   2) REEVE, CHARLES P., 'ACCURATE COMPUTATION OF THE LOG OF THE GAMMA
    !      FUNCTION', STATISTICAL ENGINEERING DIVISION NOTE 86-1, OCTOBER
    !      1986.
    !-----------------------------------------------------------------------
    !
    DOUBLE PRECISION ABSACC,B1,B2,B3,B4,B5,B6,B7,B8,C,DX,Q,R,XMIN,XN
    DATA XMIN,ABSACC / 6.894D0,1.0D-15 /
    DATA C / 0.918938533204672741780329736D0 /
    DATA B1 / 0.833333333333333333333333333D-1 /
    DATA B2 / -0.277777777777777777777777778D-2 /
    DATA B3 / 0.793650793650793650793650794D-3 /
    DATA B4 / -0.595238095238095238095238095D-3 /
    DATA B5 / 0.841750841750841750841750842D-3 /
    DATA B6 / -0.191752691752691752691752692D-2 /
    DATA B7 / 0.641025641025641025641025641D-2 /
    DATA B8 / -0.295506535947712418300653595D-1 /
    !
    !--- TERMINATE EXECUTION IF X<=0.0
    !
    IF (X.LE.0.0) STOP  '*** X<=0.0 IN FUNCTION DGAMLN ***'
    DX = DBLE(X)
    N = MAX0(0,INT(XMIN-DX+1.0D0))
    XN = DX+DBLE(N)
    R = 1.0D0/XN
    Q = R*R
    DGAMLN = R*(B1+Q*(B2+Q*(B3+Q*(B4+Q*(B5+Q*(B6+Q*(B7+Q*B8)))))))+C &
        +(XN-0.5D0)*DLOG(XN)-XN
    !
    !--- USE RECURRENCE RELATION WHEN N>0 (X<XMIN)
    !
    IF (N.GT.0) THEN
        Q = 1.0D0
        DO 20 I  =  0, N-1
            Q = Q*(DX+DBLE(I))
20      CONTINUE
        DGAMLN=DGAMLN-DLOG(Q)
    ENDIF
    !
    !--- PRINT WARNING IF ABSOLUTE ACCURACY HAS NOT BEEN ATTAINED
    !
    IF (DGAMLN+ABSACC.EQ.DGAMLN) THEN
        !CCC    PRINT *,' ********* WARNING FROM FUNCTION DGAMLN *********'
        !CCC    PRINT *,' REQUIRED ABSOLUTE ACCURACY NOT ATTAINED FOR X = ',X
    ENDIF
    RETURN
    END