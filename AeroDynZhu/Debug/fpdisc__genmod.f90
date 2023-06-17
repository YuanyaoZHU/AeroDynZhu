        !COMPILER-GENERATED INTERFACE MODULE: Fri Jul 31 01:12:44 2020
        MODULE FPDISC__genmod
          INTERFACE 
            SUBROUTINE FPDISC(T,N,K2,B,NEST)
              INTEGER(KIND=4) :: NEST
              INTEGER(KIND=4) :: K2
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: T(N)
              REAL(KIND=4) :: B(NEST,K2)
            END SUBROUTINE FPDISC
          END INTERFACE 
        END MODULE FPDISC__genmod
