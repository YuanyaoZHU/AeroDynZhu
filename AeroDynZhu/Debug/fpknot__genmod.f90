        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 31 01:12:44 2020
        MODULE FPKNOT__genmod
          INTERFACE 
            SUBROUTINE FPKNOT(X,M,T,N,FPINT,NRDATA,NRINT,NEST,ISTART)
              INTEGER(KIND=4) :: NEST
              INTEGER(KIND=4) :: M
              REAL(KIND=4) :: X(M)
              REAL(KIND=4) :: T(NEST)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: FPINT(NEST)
              INTEGER(KIND=4) :: NRDATA(NEST)
              INTEGER(KIND=4) :: NRINT
              INTEGER(KIND=4) :: ISTART
            END SUBROUTINE FPKNOT
          END INTERFACE 
        END MODULE FPKNOT__genmod
