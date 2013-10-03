      MODULE CONSTANTS2
!
!     Determine double precision and set constants.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: DBLE_PREC = KIND(0.0D0)
      REAL(KIND=DBLE_PREC), PARAMETER :: ZERO  = 0.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: ONE   = 1.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: TWO   = 2.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: THREE = 3.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: FOUR  = 4.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: FIVE  = 5.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: SIX   = 6.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: SEVEN = 7.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: EIGHT = 8.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: NINE  = 9.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: TEN   = 10.0_DBLE_PREC
      REAL(KIND=DBLE_PREC), PARAMETER :: HALF  = ONE/TWO
      REAL(KIND=DBLE_PREC), PARAMETER :: PI    = FOUR*ATAN(ONE)
      REAL(KIND=DBLE_PREC), PARAMETER :: DBLE127 = TWO**7-ONE
      REAL(KIND=DBLE_PREC), PARAMETER :: DBLE254 = TWO*DBLE127
      INTEGER, PARAMETER :: OUTPUT_UNIT          = 50
      INTEGER, PARAMETER :: SUMMARY_UNIT          = 51
      LOGICAL :: STOP_PROGRAM                        = .FALSE.
      LOGICAL ::                   VERBOSE                     = .TRUE.
      INTEGER(KIND=1), DIMENSION(:,:), ALLOCATABLE, SAVE :: SNP_BIT



!
!     Length constants.
!
      INTEGER, PARAMETER :: CHROMOSOME_LEN = 8
      INTEGER, PARAMETER :: LOCUS_LONG_LEN = 16
      INTEGER, PARAMETER :: NAME_LEN = 8
      INTEGER, PARAMETER :: PHENOTYPE_LEN = 10
      INTEGER, PARAMETER :: REAL_LEN_OUTPUT = 16
      INTEGER(KIND=1),      PARAMETER ::  ZERO_SINGLE_BYTE     = 0
      INTEGER, PARAMETER :: BITS_PER_ELEMENT            = BIT_SIZE(ZERO_SINGLE_BYTE)
      INTEGER, PARAMETER :: BITS_PER_GENOTYPE           = 2
      INTEGER, PARAMETER :: GENOTYPE_BIT_LAST_POSITION  = BITS_PER_ELEMENT-BITS_PER_GENOTYPE
      INTEGER(KIND=1),      PARAMETER :: THREE_SINGLE_BYTE     = 3
      REAL(KIND=DBLE_PREC)       :: ABSENT = HUGE(ZERO)
!
!
!     Other constants 
      CHARACTER(LEN=NAME_LEN) :: FEMALE  = "F"
      CHARACTER(LEN=NAME_LEN) :: MALE = "M"
      CHARACTER(LEN=32) :: IMPUTATION_METHOD  = "HIGHEST_POSTERIOR_PROBABILITY"
      INTEGER :: MAXIMUM_ITERATIONS = 200


! Begin hard coded constants, that should be in a Control file      
      CHARACTER(LEN=32) :: INPUTFILE = "input"
      INTEGER :: PEOPLE 
      INTEGER :: BINARY_SNPS
      INTEGER :: FLANKING_SNPS
      INTEGER :: MODEL
!      CHARACTER(LEN=9) :: KERNEL_PATH = "kernels.c"
! End hard coded

      LOGICAL :: ALL_SNPS_XLINKED
      REAL(KIND=DBLE_PREC) :: CONVERGENCE_CRITERION = TEN**(-4)
      REAL(KIND=DBLE_PREC), DIMENSION(2) :: TUNING_CONSTANT
! = 0
      !REAL(KIND=DBLE_PREC), DIMENSION(2) :: TUNING_CONSTANT = (/TEN**(-2),TEN**3/)
!     
!     DEFINE A DATA STRUCTURE TO STORE SNP DATA.         
!     TYPE CAN HAVE THE FOLLOWING VALUES:
!     A => AUTOSOMAL
!     M => MITOCHONDRIAL
!     X => X-LINKED
!     Y => Y-LINKED
!     Z => PSEUDO-AUTOSOMAL REGION OF X AND Y
!     LINKAGE_BREAK INDICATES IF THERE IS A LINKAGE BREAK
!     BETWEEN THE PREVIOUS SNP AND THE CURRENT SNP.
!     
      TYPE SNP_DATA
         CHARACTER(LEN=LOCUS_LONG_LEN) :: LONG_NAME
         CHARACTER(LEN=CHROMOSOME_LEN) :: CHROMOSOME
         INTEGER :: BASE_PAIR = 0
         INTEGER :: MISSING_GENOTYPES = 0
         INTEGER :: GROUP_ID = 0
         INTEGER :: NEXT_IN_GROUP = 0
         LOGICAL :: LINKAGE_BREAK
         LOGICAL :: OMIT
         LOGICAL :: RETAINED_IN_MODEL = .FALSE.
         LOGICAL :: ASSIGNED_WEIGHT = .FALSE.
         REAL(KIND=DBLE_PREC) :: WEIGHT = ONE
         CHARACTER(LEN=1) :: TYPE
      END TYPE SNP_DATA
      TYPE(SNP_DATA), DIMENSION(:), ALLOCATABLE, SAVE :: SNP_LOCUS
!
!     DEFINE A DATA STRUCTURE TO STORE PERSON DATA.  THE COMPONENT NUMBER_TYPE
!     REPRESENTS A PHENOTYPE BY ITS POSITION IN AN ORDERED LIST ON PHENOTYPES.
!     NOTE THAT THE COMPONENT VARIABLE_CHAR CONTAINS THE CHARACTER STRING OF
!     THE UNTRANSFORMED VARIABLE VALUE.
!
      TYPE PERSON_DATA
         CHARACTER(LEN=NAME_LEN) :: PEDIGREE_NAME
         CHARACTER(LEN=NAME_LEN) :: NAME
         CHARACTER(LEN=NAME_LEN) :: FATHER
         CHARACTER(LEN=NAME_LEN) :: MOTHER
         CHARACTER(LEN=NAME_LEN) :: SEX
         CHARACTER(LEN=NAME_LEN) :: TWIN
         INTEGER :: GROUP_NUMBER
         INTEGER :: PEDIGREE_NUMBER
         INTEGER :: POPULATION_NUMBER = 0
         INTEGER :: MISSING_GENOTYPES = 0
         INTEGER :: XLINKED_HET = 0
         LOGICAL :: TYPED_AT_SNPS = .FALSE.
         CHARACTER(LEN=PHENOTYPE_LEN), DIMENSION(:), POINTER :: PHENOTYPE => NULL()
         CHARACTER(LEN=REAL_LEN_OUTPUT), DIMENSION(:), POINTER :: VARIABLE_CHAR => NULL()
         INTEGER, DIMENSION(:), POINTER :: NUMBER_TYPE => NULL()
         REAL(KIND=DBLE_PREC), DIMENSION(:), POINTER :: VARIABLE => NULL()
      END TYPE PERSON_DATA
      TYPE(PERSON_DATA), DIMENSION(:), ALLOCATABLE, SAVE :: PERSON
      END MODULE CONSTANTS2

      MODULE IMPUTATION

      USE CONSTANTS
      USE LOCUS_AND_PEDIGREE_STRUCTURES

      !REAL, DIMENSION(:,:,:), ALLOCATABLE :: SNP_PENETRANCE
      !INTEGER :: I
      CHARACTER(LEN=32) :: ARG
      
! BEGIN GPU VARIABLES
      INTEGER :: SCALER,CHUNKSIZE 
      INTEGER, DIMENSION(256) :: RETURN_VEC
      LOGICAL :: TEST_GPU
! END GPU SECTION

      REAL :: DELTA
      REAL :: LAMBDA
      INTEGER :: PLATFORM_ID
      INTEGER :: DEVICE_ID 
      INTEGER :: GENO_DIM
      INTEGER :: MAX_HAPLOTYPES,MAX_EXTENDED_HAPLOTYPES
      INTEGER :: HAPLOTYPE_MODE
      INTEGER :: TOTAL_REGIONS

      INTEGER :: CSNP,DSNP,HALF_MAX,HAPLOTYPES,H,I,ITERATION,J,K,L,M,MARKERS
      INTEGER :: MAX_WINDOW,N,QUARTER_MAX,REGION,REGIONS,SNP,START
      LOGICAL :: GENOTYPE_IMPUTATION,BEST_HAPLOTYPE_PAIR,DIPLOID,LOW_FREQUENCY
      INTEGER :: INT_GENOTYPE_IMPUTATION
! GPU CONSTANTS
      INTEGER :: BLOCK_WIDTH = 256
      INTEGER :: SMALL_BLOCK_WIDTH = 64
!
      INTEGER, DIMENSION(2) :: BEST_PAIR,IMPUTED_GENOTYPE,WINDOW
      INTEGER, DIMENSION(3) :: WINDOW_TRANSITION
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PERMUTATION
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: BOUNDARY,HAPLOTYPE

      REAL :: BEST_PROB,D,E,FREQ
      REAL :: LIKELIHOOD,P,PENETRANCE,Q,T,T1,T2,MAX_PROB
      REAL, PARAMETER :: EPS = TEN**(-10),DROP_CRITERION = TEN**(-8)
      REAL, DIMENSION(0:3) :: POSTERIOR_PROB
      REAL, ALLOCATABLE, DIMENSION(:) :: FREQUENCY,OLD_FREQUENCY
      REAL, ALLOCATABLE, DIMENSION(:) :: CURRENT_WEIGHT,WEIGHT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ACTIVE_HAPLOTYPE
      LOGICAL :: ENABLE_FORTRAN
      INTEGER :: MAX_REGION_SIZE
      INTEGER :: HAPMODE_DENOVO = 0
      INTEGER :: HAPMODE_GUIDE = 1
      INTEGER :: LAST_REGION_SNP,CSNP_START,CSNP_END
      INTEGER :: PREVIOUS_LEFT_MARKER,PREVIOUS_RIGHT_MARKER
     
!
      CONTAINS
!
      SUBROUTINE NEW_GENOTYPE_HAPLOTYPE_IMPUTATION()
      !SUBROUTINE NEW_GENOTYPE_HAPLOTYPE_IMPUTATION(SNP_PENETRANCE)
!
!     THIS SUBROUTINE ESTIMATES HAPLOTYPE FREQUENCIES BY A PENALIZED MM
!     ALGORITHM AND IMPUTES MISSING SNP GENOTYPES OR HAPLOTYPES. MODEL
!     1 IS INTENDED FOR GENOTYPE IMPUTATION AND MODEL 2 FOR HAPLOTYPING.
!
      IMPLICIT NONE
!
      !REAL, DIMENSION(BINARY_SNPS,0:2,PEOPLE) :: SNP_PENETRANCE

!
!     INITIALIZE CONSTANTS.
!
      MAX_REGION_SIZE = BINARY_SNPS/TOTAL_REGIONS + 1
      GENOTYPE_IMPUTATION = MODEL==1
      !OPEN(200,file="FORTRAN_GENOTYPES")
      HALF_MAX = MAX_HAPLOTYPES/2
      QUARTER_MAX = MAX_HAPLOTYPES/4
      IF (HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
         MAX_WINDOW = 3*FLANKING_SNPS
      ELSE
         MAX_WINDOW = 2*FLANKING_SNPS+1
      END IF
      BEST_HAPLOTYPE_PAIR = IMPUTATION_METHOD=="BEST_HAPLOTYPE_PAIR"
      IF(GENOTYPE_IMPUTATION) THEN
         INT_GENOTYPE_IMPUTATION=1
         print *,"Running genotype imputation"
      ELSE
         INT_GENOTYPE_IMPUTATION=0
         print *,"Running haplotype phasing"
      END IF
      REGIONS = COUNT(SNP_LOCUS%LINKAGE_BREAK)
!
!     ALLOCATE ARRAYS.
!
print *,"regions fortran: ",REGIONS
print *,"max_window fortran: ",MAX_WINDOW
print *,"max_hap fortran: ",MAX_HAPLOTYPES
      ALLOCATE(BOUNDARY(2,REGIONS))
      ALLOCATE(HAPLOTYPE(MAX_WINDOW,MAX_HAPLOTYPES))
      ALLOCATE(OLD_FREQUENCY(MAX_HAPLOTYPES))
      ALLOCATE(CURRENT_WEIGHT(MAX_HAPLOTYPES))
      ALLOCATE(WEIGHT(MAX_HAPLOTYPES))
      ALLOCATE(FREQUENCY(MAX_HAPLOTYPES))
      ALLOCATE(ACTIVE_HAPLOTYPE(MAX_HAPLOTYPES))

      CALL init_buffers(ACTIVE_HAPLOTYPE,MAX_WINDOW,MAX_HAPLOTYPES,&
      MAX_EXTENDED_HAPLOTYPES,MAX_REGION_SIZE, PEOPLE,BINARY_SNPS,&
      INT_GENOTYPE_IMPUTATION, HAPLOTYPES,MARKERS,WINDOW(1),&
      PREVIOUS_LEFT_MARKER,WINDOW(2),PREVIOUS_RIGHT_MARKER,FREQUENCY,&
      WEIGHT,HAPLOTYPE,FLANKING_SNPS)
      !CALL init_buffers(FLANKING_SNPS,ACTIVE_HAPLOTYPE,MAX_WINDOW,&
      !MAX_HAPLOTYPES,&
      !MAX_EXTENDED_HAPLOTYPES,MAX_REGION_SIZE, PEOPLE,BINARY_SNPS,&
      !SNP_PENETRANCE, INT_GENOTYPE_IMPUTATION, HAPLOTYPES,MARKERS,WINDOW(1),&
      !PREVIOUS_LEFT_MARKER,FREQUENCY,WEIGHT,HAPLOTYPE)

!
!     COMPUTE THE NUMBER OF CHROMOSOME REGIONS AND THE BOUNDARIES OF EACH REGION.
!
      K = 0
      print*,"BINARY_SNPS = ",BINARY_SNPS
      print*,"FLANKING_SNPS = ",FLANKING_SNPS
      print*,"REGIONS",REGIONS
      DO L = 1,BINARY_SNPS
         IF (SNP_LOCUS(L)%LINKAGE_BREAK) THEN
            K = K+1
            BOUNDARY(1,K) = L
            print*,"start of REGION ",K," is ",L
            IF (K>1) THEN
               BOUNDARY(2,K-1) = L-1
               print*,"end of REGION ",K-1," is ",L-1
            END IF
         END IF
      END DO
      BOUNDARY(2,REGIONS) = BINARY_SNPS
      print*,"end of REGION ",REGIONS," is ",BINARY_SNPS
!
!     PROCESS EACH CHROMOSOME REGION.
!


      DO REGION = 1,REGIONS
         print*,"REGION:",REGION
         IF (BOUNDARY(2,REGION)-BOUNDARY(1,REGION)+1 > MAX_REGION_SIZE) THEN
            print *,"Region",REGION,"of length",BOUNDARY(2,REGION)-BOUNDARY(1,REGION)+1,"exceeded maximum size",MAX_REGION_SIZE
            STOP 1
         END IF
         CALL init_region_buffers(BOUNDARY(1,REGION),BOUNDARY(2,REGION));
!
!     LOOP OVER ALL SNPS IN THE CURRENT REGION.  CSNP IS THE CENTRAL SNP.
!
         HAPLOTYPES = 1
         FREQUENCY(1:MAX_HAPLOTYPES) = ONE/REAL(HAPLOTYPES,KIND=DBLE_PREC)
         MARKERS = 0
         WINDOW(1) = BOUNDARY(1,REGION)
         WINDOW(2) = WINDOW(1)-1
         CSNP_START = 0
         WINDOW_TRANSITION(1) = BOUNDARY(1,REGION)+FLANKING_SNPS
         WINDOW_TRANSITION(2) = WINDOW_TRANSITION(1)+FLANKING_SNPS
         WINDOW_TRANSITION(3) = BOUNDARY(2,REGION)
         print*,"REGION start: ",BOUNDARY(1,REGION),", end:",BOUNDARY(2,REGION)+FLANKING_SNPS
         IF (HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
            LAST_REGION_SNP = BOUNDARY(2,REGION)
         ELSE
            LAST_REGION_SNP = BOUNDARY(2,REGION)+FLANKING_SNPS
         END IF
         SNP = BOUNDARY(1,REGION)
         PREVIOUS_LEFT_MARKER = WINDOW(1)
         PREVIOUS_RIGHT_MARKER = WINDOW(2)
! SNP = 381
         DO WHILE (SNP<=LAST_REGION_SNP)
         !DO SNP = BOUNDARY(1,REGION),BOUNDARY(2,REGION)+FLANKING_SNPS
            print*," SNP = ",SNP
            IF (HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
               IF(SNP-FLANKING_SNPS<BOUNDARY(1,REGION)) THEN
                  ! THE FIRST WINDOW
                  WINDOW(1) = BOUNDARY(1,REGION)
                  WINDOW(2) = WINDOW(1)+FLANKING_SNPS*2-1
                  CSNP_START =  1
                  CSNP_END = CSNP_START+FLANKING_SNPS-1
               ELSE IF(SNP-FLANKING_SNPS == BOUNDARY(1,REGION)) THEN
                  ! THE SECOND WINDOW
                  WINDOW(1) = BOUNDARY(1,REGION)
                  WINDOW(2) = WINDOW(1)+MAX_WINDOW-1 
                  CSNP_START = FLANKING_SNPS + 1
                  CSNP_END = CSNP_START+FLANKING_SNPS-1
               ELSE IF (SNP+2*FLANKING_SNPS>BOUNDARY(2,REGION)) THEN
                  ! THE LAST WINDOW
                  WINDOW(1) = WINDOW(1)+FLANKING_SNPS
                  WINDOW(2) = BOUNDARY(2,REGION)
                  CSNP_START = FLANKING_SNPS + 1
                  CSNP_END = CSNP_START+FLANKING_SNPS-1
               ELSE
                  ! ALL OTHER WINDOWS
                  WINDOW(1) = WINDOW(1)+FLANKING_SNPS
                  WINDOW(2) = WINDOW(1)+MAX_WINDOW-1
                  CSNP_START = FLANKING_SNPS + 1
                  CSNP_END = CSNP_START+FLANKING_SNPS-1
               END IF
               ! CHECK TO SEE THAT WE DON'T OVERFLOW OFF THE REGION
               IF(WINDOW(1)+CSNP_END-1>BOUNDARY(2,REGION)) THEN
                  CSNP_END=CSNP_END-(WINDOW(1)+CSNP_END-BOUNDARY(2,REGION)-1)
                  print *,"CSNP_END overflowed the region"
               END IF
            ELSE
!
!     IMPUTATION IS PERFORMED IN A WINDOW EXTENDING FROM WINDOW(1) TO WINDOW(2).
!     MOVE THE WINDOW OVER A SPACE. CSNP IS THE CENTRAL SNP WHERE IMPUTATION OCCURS.
!
               IF (SNP<WINDOW_TRANSITION(1)) THEN
                  WINDOW(2) = WINDOW(2)+1
               ELSE IF (SNP<=WINDOW_TRANSITION(2)) THEN
                  CSNP_START = CSNP_START+1
                  WINDOW(2) = WINDOW(2)+1
               ELSE IF (SNP<=WINDOW_TRANSITION(3)) THEN
                  WINDOW(1) = WINDOW(1)+1
                  WINDOW(2) = WINDOW(2)+1
               ELSE
                  WINDOW(1) = WINDOW(1)+1
               END IF
               IF (WINDOW(2)>=BOUNDARY(2,REGION)) THEN
                  WINDOW(2) = BOUNDARY(2,REGION)
               END IF
               CSNP_END = CSNP_START
            END IF
            MARKERS = WINDOW(2)-WINDOW(1)+1
            print*,"WINDOW start: ",WINDOW(1),", end: ",WINDOW(2)
            print*,"CSNP_START = ",CSNP_START
            print*,"CSNP_END = ",CSNP_END
            print*,"DSNP RANGE: ",WINDOW(1)+CSNP_START-1,WINDOW(1)+CSNP_END-1
            print*,"MARKERS = ",MARKERS
if (1==1) THEN

            IF (HAPLOTYPE_MODE==HAPMODE_GUIDE .OR. CSNP_START>=0) THEN

               IF (HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
                  print *,"HAPLOTYPES BEFORE REFERENCE COPY: ",HAPLOTYPES
                  CALL copy_ref_haplotypes
                  print *,"HAPLOTYPES AFTER REFERENCE COPY: ",HAPLOTYPES
               ELSE IF (HAPLOTYPE_MODE==HAPMODE_DENOVO) THEN
   !
   !     DUPLICATE EACH HAPLOTYPE AND APPEND ONE OF THE TWO ALLELES TO EACH
   !     EXTENSION. IN THE PRESENT DUPLICATION SCHEME, HAPLOTYPES APPEAR IN
   !     REVERSE DICTIONARY ORDER. THIS CHOICE AVOIDS SORTING OF HAPLOTYPES
   !     WHEN THEY ARE CROPPED AND CONSOLIDATED.
   !
                  CALL double_haplotypes
               END IF
   !
   !     INITIALIZE HAPLOTYPE FREQUENCIES.
   !
               !FREQUENCY(1:MAX_HAPLOTYPES) = ONE/REAL(HAPLOTYPES,KIND=DBLE_PREC)
               !OLD_FREQUENCY(1:MAX_HAPLOTYPES) = ZERO
               OLD_FREQUENCY(1:MAX_HAPLOTYPES) = FREQUENCY(1:MAX_HAPLOTYPES)
   
   !
   !     ESTIMATE HAPLOTYPE FREQUENCIES BY THE PENALIZED MM ALGORITHM.
   !
               CALL init_window_buffers
               DO ITERATION = 1,MAXIMUM_ITERATIONS
                  CALL init_iteration_buffers
   !
   !     INITIALIZE THE WEIGHTED HAPLOTYPE COUNTS OVER THE ENTIRE SAMPLE.
   !
                  IF(ITERATION==1) THEN
                      IF (HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
                         CALL precompute_penetrance
                      ELSE IF (HAPLOTYPE_MODE==HAPMODE_DENOVO) THEN
                         CALL impute_penetrance_matrix
                         CALL precompute_penetrance_fast(CSNP_START)
                      END IF
                  END IF
                  CALL compute_haplotype_weights (ITERATION,WEIGHT,HAPLOTYPES,MARKERS)
                  !print *,"compute hapweight at iteration",ITERATION
   !
   !     COMPUTE THE PENALTY INGREDIENTS OF THE MM UPDATE.
   !
                  D = ZERO
                  E = ZERO
                  DO I = 1,MAX_HAPLOTYPES
                     IF(ACTIVE_HAPLOTYPE(I)==1) THEN
                        IF (WEIGHT(I)<EPS) THEN
                           WEIGHT(I) = EPS
                           !print *,"WEIGHT",(I-1)," bounded to ",EPS
                        END IF
                        IF (FREQUENCY(I)>DELTA) THEN
                           D = D+WEIGHT(I)
                        ELSE
                           E = E+WEIGHT(I)
                        END IF
                     END IF 
                  END DO
                  T = D+E
   !
   !     PERFORM THE MM UPDATES OF THE HAPLOTYPE FREQUENCIES. THESE REDUCE TO
   !     TO THE EM UPDATES WHENEVER LAMBDA, D, OR E IS 0.
   !
                  IF (LAMBDA*D*E>ZERO) THEN
                     Q = TWO*E/(D+E+LAMBDA+SQRT((D+E+LAMBDA)**2-FOUR*LAMBDA*E))
                     DO I = 1,MAX_HAPLOTYPES
                        IF (ACTIVE_HAPLOTYPE(I)==1) THEN
                           IF (FREQUENCY(I)>DELTA) THEN
                              FREQUENCY(I) = WEIGHT(I)/(T-LAMBDA*Q)
                           ELSE
                              FREQUENCY(I) = WEIGHT(I)/(T-LAMBDA*Q+LAMBDA)
                           END IF
                           IF (FREQUENCY(I)<EPS) THEN
                             FREQUENCY(I) = EPS
                             print *,"HAP",(I-1)," bounded to freq ",FREQUENCY(I)
                           ELSE
                             print *,"HAP",(I-1)," estimated at freq ",FREQUENCY(I)
                           END IF
                        END IF
                     END DO
                     print *,"Total haplotypes: ",HAPLOTYPES
                  ELSE IF (T>ZERO) THEN
                     FREQUENCY(1:MAX_HAPLOTYPES) = WEIGHT(1:MAX_HAPLOTYPES)/T
                     DO I = 1,MAX_HAPLOTYPES
                        IF (ACTIVE_HAPLOTYPE(I)==1) THEN
                          print *,"HAP",(I-1)," (T>0) estimated at freq ",FREQUENCY(I)
                        END IF
                     END DO
                  ELSE
                     FREQUENCY(1:MAX_HAPLOTYPES) = ONE/REAL(MAX_HAPLOTYPES,KIND=DBLE_PREC)
                     DO I = 1,MAX_HAPLOTYPES
                        IF (ACTIVE_HAPLOTYPE(I)==1) THEN
                          print *,"HAP",(I-1)," (FLAT) estimated at freq ",FREQUENCY(I)
                        END IF
                     END DO
                  END IF
   !
   !     CHECK FOR CONVERGENCE. UPDATE THE OLD HAPLOTYPE FREQUENCIES AND PROCEED
   !     TO THE NEXT ITERATION WHEN THE CONVERGENCE TEST FAILS.
   !
                  T = ZERO
                  DO I = 1,MAX_HAPLOTYPES
                     IF (ACTIVE_HAPLOTYPE(I)==1) THEN
                        T = T+ABS(FREQUENCY(I) - OLD_FREQUENCY(I))
                     END IF
                  END DO
                  !print *,"L1 raw norm is ",T
                  T = T/HAPLOTYPES
                  !print *,"L1 adjusted norm is ",T
                  IF (T<CONVERGENCE_CRITERION) THEN
                     EXIT
                  ELSE
                     OLD_FREQUENCY(1:MAX_HAPLOTYPES) = FREQUENCY(1:MAX_HAPLOTYPES)
                  END IF
               END DO ! END OF CONVERGENCE LOOP
               print *,"Converged at iteration",ITERATION,"of",MAXIMUM_ITERATIONS

            END IF ! END IF WINDOW HAS CHANGED
!
!     IMPUTE AT THE CENTRAL SNP. DSNP IS ITS POSITION IN THE SNP DATA FILE.
!
            IF(HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
               CALL prep_impute_genotypes_guide
            END IF
            IF (CSNP_START>=1) THEN
               IF(HAPLOTYPE_MODE==HAPMODE_DENOVO) THEN
                  DSNP = WINDOW(1)+CSNP_START-1
                  print*,"passing csnp_start,dsnp:",CSNP_START,DSNP
                  IF(GENO_DIM==4) THEN
                     CALL impute_haploid_genotypes_denovo(CSNP_START,DSNP)
                  ELSE
                     CALL impute_diploid_genotypes_denovo(CSNP_START,DSNP)
                  END IF
               ELSE IF(HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
                  IF(GENO_DIM==4) THEN
                     CALL impute_haploid_genotypes_guide(CSNP_START,CSNP_END,WINDOW(1))
                  ELSE
                     CALL impute_diploid_genotypes_guide(CSNP_START,CSNP_END,WINDOW(1))
                  END IF
               END IF
            END IF
!
!     DROP LOW FREQUENCY HAPLOTYPES.

            IF(HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
               print *,"Skipping haplotype pruning"
            ELSE
               IF (CSNP_START>=0) THEN
                  CALL prune_haplotypes
               END IF ! END IF WINDOW HAS CHANGED
            END IF ! END HAPLOTYPE_MODE
end if
            IF(HAPLOTYPE_MODE==HAPMODE_GUIDE) THEN
               print *,"Skipping haplotype pruning"
               SNP = SNP + FLANKING_SNPS
            ELSE
               SNP = SNP + 1
            END IF
            PREVIOUS_LEFT_MARKER = WINDOW(1)
            PREVIOUS_RIGHT_MARKER = WINDOW(2)
         END DO !LOOP OVER ALL SNPS
      END DO !LOOP REGIONS
      CALL cleanup()
      !CLOSE(200)
!print*,"EXITING NEW GENOTYPE IMPUTATION"
      END SUBROUTINE NEW_GENOTYPE_HAPLOTYPE_IMPUTATION

      SUBROUTINE SORT_REALS(LIST)
!
!     THIS SUBROUTINE PERFORMS A HEAP SORT ON A LIST OF REALS.
!     SEE: NIJENHUIS A AND WILF HS (1978) "COMBINATORIAL ALGORITHMS
!     FOR COMPUTERS AND CALCULATORS, 2ND ED", CHAPTER 15, ACADEMIC PRESS.
!
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(INOUT) :: LIST
      REAL :: TEMP

      
      !REAL(KIND=DBLE_PREC) :: TEMP
      !REAL(KIND=DBLE_PREC), DIMENSION(:), INTENT(INOUT) :: LIST

      INTEGER :: I,J,K,L,N

      N = SIZE(LIST)
      IF (N<=1) RETURN
      L = 1+N/2
      K = N
      DO
         IF (L>1) THEN
            L = L-1
            TEMP = LIST(L)
         ELSE
            TEMP = LIST(K)
            LIST(K) = LIST(1)
            K = K-1
            IF (K<=1) THEN
               LIST(1) = TEMP
               RETURN
            END IF
         END IF
         I = L
         J = L+L
         DO WHILE (J<=K)
            IF (J<K) THEN
               IF (LIST(J+1)>LIST(J)) J = J+1
            END IF
            IF (LIST(J)>TEMP) THEN
               LIST(I) = LIST(J)
               I = J
               J = J+J
            ELSE
               J = K+1
            END IF
         END DO
         LIST(I) = TEMP
      END DO
      END SUBROUTINE SORT_REALS      
      END MODULE IMPUTATION

      PROGRAM IMPUTE

      USE IMPUTATION
      USE MATRIX_COMPLETION_ROUTINES
      
      CALL MAIN

      CONTAINS

      SUBROUTINE MAIN

      IMPLICIT NONE

      !REAL(8) :: VALIDATION_FRACTION

      CALL load_constants(MODEL,PEOPLE,BINARY_SNPS,TOTAL_REGIONS,&
      HAPLOTYPE_MODE,FLANKING_SNPS,MAX_HAPLOTYPES,MAX_EXTENDED_HAPLOTYPES, &
      PLATFORM_ID,DEVICE_ID,DELTA,LAMBDA,GENO_DIM)
     print *,"geno dim ",GENO_DIM

      ALLOCATE(SNP_LOCUS(BINARY_SNPS))
      DO I=1,BINARY_SNPS
        SNP_LOCUS(I)%LINKAGE_BREAK = .FALSE.
        SNP_LOCUS(I)%OMIT = .FALSE.
      END DO
      CHUNKSIZE = BINARY_SNPS/TOTAL_REGIONS
      I = 1
      DO WHILE(I<=BINARY_SNPS)
        SNP_LOCUS(I)%LINKAGE_BREAK = .TRUE.
        !print *,"Installing a linkage break at SNP",I
        I = I + CHUNKSIZE
      END DO

! Just to test regions:      

      ALLOCATE(PERSON(PEOPLE))
      DO I=1,PEOPLE
        PERSON(I)%TYPED_AT_SNPS = .TRUE.
        PERSON(I)%SEX = MALE
      END DO
! END This would probably get set up somewhere else in Mendel

      !IF (GENO_DIM==3) THEN
        !ALLOCATE(SNP_PENETRANCE(0:2,BINARY_SNPS,PEOPLE))
      !ELSE IF (GENO_DIM==4) THEN
        !ALLOCATE(SNP_PENETRANCE(0:3,BINARY_SNPS,PEOPLE))
       
      !END IF
!!      ALLOCATE(FRACTIONAL_COUNT(PEOPLE,BINARY_SNPS))

!      TEST_GPU = .FALSE.
      TEST_GPU = .TRUE.
      IF (TEST_GPU) THEN
      print *,"calling initgpu"
        CALL init_gpu(PLATFORM_ID,DEVICE_ID)
!     RUN THE PROGRAM AS NORMAL 
        !CALL read_stream(PEOPLE,BINARY_SNPS)
        !CALL read_stream(PEOPLE,BINARY_SNPS,SNP_PENETRANCE)
        IF (MODEL==1 .OR. MODEL==2) THEN
          CALL NEW_GENOTYPE_HAPLOTYPE_IMPUTATION()
          !CALL NEW_GENOTYPE_HAPLOTYPE_IMPUTATION(SNP_PENETRANCE)
        ELSE IF (MODEL==3) THEN
          CALL MATRIX_COMPLETION_IMPUTATION
        END IF
      END IF

      END SUBROUTINE MAIN

      END PROGRAM IMPUTE
