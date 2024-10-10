
!! create a module called "IH_INNERS" with a number of subroutines
! **** Global variables

MODULE IH_INNERS

  IMPLICIT NONE

  ! makes all variables and functions in the module accessible to other program units

  PUBLIC

  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)

  INTEGER, SAVE :: nodeSize ! minimum number of cases in a node
  INTEGER, SAVE :: mTry ! maximum number of covariates to try
  INTEGER, SAVE :: minEvent ! minimum number of events in a node
  INTEGER, SAVE :: ERT ! 0/1 1 = use extremely randomized tree
  INTEGER, SAVE :: uniformSplit ! 0/1 1 = random cutoff comes from values
  INTEGER, SAVE :: nt ! number of time points
  INTEGER, SAVE :: rule ! 1 if logrank 0 if truncated mean
  INTEGER, SAVE :: nLevs ! maximum number of levels

  ! used in setupInners:
  INTEGER, SAVE :: np ! number of covariates
  INTEGER, SAVE :: nAll ! number of cases
  INTEGER, SAVE :: sampleSize
  INTEGER, SAVE :: nTree ! number of trees
  INTEGER, SAVE :: nrNodes ! maximum number of nodes in a tree

  ! used in tsurvTree
  INTEGER, SAVE :: replace
  INTEGER, SAVE :: n  ! number of sampled cases

  ! used in setUpBasics & tsurvTree
  ! for survival probability, the index of nearest time point
  INTEGER, SAVE :: sIndex

  ! -------------------------------------------------------!


  ! Allocatable Arrays to store info for survival analysis !
  ! -------------------------------------------------------!
   ! number of categories in each np covariate
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nCat

    ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: delta


      ! used in setupInners:
    ! censoring indicator 1 = not censored
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: deltaAll

        ! propensity score for sampled cases
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: propensityALL

  ! -------------------------------------------------------!


    ! the probability for a random split: this is called in setUpBasics
  REAL(dp), SAVE :: rs
    ! probability mass vector of survival function for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: pr
  REAL(dp), DIMENSION(:), ALLOCATABLE, SAVE :: dt
    ! covariates to be considered for split for sampled cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: x

    ! used in tsurvTree
    ! the stratified random split coefficient
  REAL(dp), SAVE :: stratifiedSplit
    ! the fraction above sIndex time that survival time lies; also used in setUpBasics
  REAL(dp), SAVE :: sFraction

      ! used in setUpInners
    ! covariates to be considered for split for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: xAll
    ! probability mass vector of survival function for all cases
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: prAll
    ! vector of propensity scores for all patients/stages
  !REAL(dp), DIMENSION(:,:), ALLOCATABLE, SAVE :: propensityAll

    ! -------------------------------------------------------!

      ! used in setUpBasics
      ! TRUE = time differences vector has been allocated
  LOGICAL, SAVE :: dtAllocated = .FALSE.

      ! used in setUpInners
    ! TRUE = all other allocatables have been allocated
  LOGICAL, SAVE :: isAllocated = .FALSE.

      ! used in tsurvTree

    ! TRUE = using survival probability
  LOGICAL, SAVE :: isSurvival = .FALSE.

    ! **************************************************************** !

  ! define a type called "NODE"
  ! node in the decision tree with various arrays for survival function, mean, survival probability, and matrix

  TYPE NODE
    INTEGER :: nNode
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: survFunc
    REAL(dp), DIMENSION(:), ALLOCATABLE :: mean
    REAL(dp), DIMENSION(:), ALLOCATABLE :: survProb
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: matrix
  END TYPE

  ! create allocatable array of type "NODE" that represents a collection of decision trees. This is called "trees"

  TYPE(NODE), DIMENSION(:), ALLOCATABLE, SAVE :: trees

    ! **************************************************************** !

  ! define a type called "ForestValues"
  ! store values related to the entire forest of decision trees.
  ! used in setUpInners

  TYPE ForestValues
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: survFunc
    REAL(dp), DIMENSION(:), ALLOCATABLE :: mean
    REAL(dp), DIMENSION(:), ALLOCATABLE :: survProb
  END TYPE

  ! instance of type "ForestValues" to store values related to the entire forest of decision trees

  TYPE(ForestValues), SAVE :: forest

   ! **************************************************************** !

  ! this defines a function called "sampleWithReplace" within the module "IH_INNERS"
  ! this is used in tsurvTREE

  CONTAINS

! sample an array of indices allowing for duplicates
!   nCases: integer, the total number of indices to sample
!   n: integer, the number of indices to draw
! returns an array (1:n) of the indices sampled

! define a function called "sampleWithReplace" that takes two input arguments "nCases" and "n"
! returns an array called "RESULT" of integers


FUNCTION sampleWithReplace(nCases, n) RESULT(array)

! all variables must be explicitly declared
  IMPLICIT NONE

  ! input arguments should be integers

  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: n

  ! named array with a size of n, indexed from 1...n

  INTEGER, DIMENSION(1:n) :: array

  INTEGER :: i

  ! declares a real variable "rnd" to store a random number

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  ! executes a loop from 1: n

  DO i = 1, n
  ! in the loop, samples a random integer within the range (0, 1) and scaled it by nCases, then floor this and increase by 1
    array(i) = 1 + floor(rnd(0.d0, 1.d0)*nCases)
  END DO

END FUNCTION sampleWithReplace

  ! **************************************************************** !

  ! this is used in tsurvTREE


! define a function called "sampleWithOutReplace"
! takes two input arguments "nCases", "n", and returns an array called "array"
FUNCTION sampleWithOutReplace(nCases, n) RESULT(array)

! sample an array of indices without allowing duplicates
!   nCases: integer, the total number of indices to sample
!   n: integer, the number of indices to draw
! returns an array (1:n) of the indices sampled

! all variables must be explicitly defined
  IMPLICIT NONE

  ! the input arguments should be integers

  INTEGER, INTENT(IN) :: nCases
  INTEGER, INTENT(IN) :: n

  ! declare an array called "array" that indexes 1:n

  INTEGER, DIMENSION(1:n) :: array

  ! declare variables used in the function: j = variable sampled, used = logical array tracking whether an index has been used
  ! rnd = real number, EXTERNAL :: rnd = rnd is an external function

  INTEGER :: j, cnt

  LOGICAL, DIMENSION(1:nCases) :: used

  REAL(dp) :: rnd

  EXTERNAL :: rnd

  ! initializes all the "used" array to FALSE to track if index has been used

  used(:) = .FALSE.


  ! initialize the loop counter to 1
  cnt = 1

  ! only continue while loop counter <= number of indices to be sampled

  DO WHILE (cnt <= n)
    ! draw a number between 0 and 1, scale it by nCases, round down using "floor", add 1
    j = 1 + floor(rnd(0.d0, 1.d0)*nCases)

    ! if the number has already been drawn cycle throught the next sampling iteration
    if (used(j)) CYCLE

    ! if the number has not previously been drawn, store in the array at the current position of the counter (cnt) and set used
    array(cnt) = j
    used(j) = .TRUE.

    ! update counter
    cnt = cnt + 1
  END DO

END FUNCTION sampleWithOutReplace

! ******************************* subroutine: tfindSplit ******************************** !





  ! Identify the optimal split

!! input parameters:
  !   nCases : integer, the number of elements in input casesIn
!   casesIn : integer(:), the indices of the cases in this node
!   nv : integer, the number of covariates to include in search
!   varsIn : integer(:), covariates to include in search

!! outputs
!   splitVar : integer, the index of the selected variable for splitting
!   cutoffBest : real(:), the cutoff (<= go to 'left' node)
!   splitFound : integer, 0 = no split; 1 = found a split
!   casesOut : integer(:), elements of casesIn that go left; ind if yes, 0
!     otherwise
!   nCuts : integer, the number of cutoff values returned
!   lft : integer, the number of cases in the left node


! define a subroutine called "tfindSplit"
! subroutine: program unit that performs a specific task or a set of task but doesn't return a value
! delete input "propensity"

SUBROUTINE tfindSplit(nCases, casesIn, nv, varsIn, &
                    & splitVar, cutoffBest, splitFound, casesOut, nCuts, lft)

  ! all variables must be explicitly defined


  IMPLICIT NONE


  ! variables inputs

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  INTEGER, INTENT(IN) :: nv
  INTEGER, DIMENSION(1:nv), INTENT(IN) :: varsIn

  ! variable outputs
  INTEGER, INTENT(OUT) :: splitVar
      ! nLevs is global variable
  REAL(dp), DIMENSION(1:nLevs), INTENT(OUT) :: cutoffBest
  INTEGER, INTENT(OUT) :: splitFound
  INTEGER, DIMENSION(1:nCases), INTENT(OUT) :: casesOut
  INTEGER, INTENT(OUT) :: nCuts
  INTEGER, INTENT(OUT) :: lft

  ! integer variables are declared for loop counters, indices, and tracking
  ! also declare allocatable arrays to store data

  INTEGER :: cnt, i, ikv, j, jj, k, kv, l, nUncensored, ptr, rightNode
  INTEGER :: rUnifSet, set, splitLeft, splitLeftFinal, tieCovariate
  INTEGER :: tieValue, variablesTried
  INTEGER, DIMENSION(1:nCases) :: cases, dSorted, tcases
  INTEGER, DIMENSION(1:nv) :: variables
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, indSingles, leftCases, rightCases
  INTEGER, DIMENSION(:), ALLOCATABLE :: uncensoredIndices

  ! declares real variables and arrays for handling floating-point numbers

  REAL(dp) :: cutoff, maxValueSplit, maxValueXm, rUnif, valuej
  REAL(dp), DIMENSION(1:nt) :: atRiskLeft, atRiskRight, D, denJ
  REAL(dp), DIMENSION(1:nt) :: eventsLeft, eventsRight, numJ, pd1, pd2, Rcum
  REAL(dp), DIMENSION(1:nCases) :: xSorted
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: prl, prr
  !REAL(dp), DIMENSION(:,:), ALLOCATABLE :: propensityLeft, propensityRight

  ! logical variables are declared to handle boolean values, create an allocatable array called "singles"

  LOGICAL :: randomSplit
  LOGICAL, DIMENSION(:), ALLOCATABLE :: singles

  ! declare a number, "rnd", which is an external function

  REAL(dp) :: rnd

  EXTERNAL :: rnd


  ! determine if this is to be a random split
  ! randomSplit is set to true if a random number is less than or equal to rs (global var prob of random split)

  randomSplit = rnd(0.d0, 1.d0) <= rs

  ! splitFound is flag indicating if a split was found (OUTPUT)
  ! 0 = no split found; 1 = split found
  !! initialize splitFound to 0

  splitFound = 0

  !! initialize tracking variables

  ! index of variable with largest critical value (OUTPUT)
  ! negative number indicates no split found
  ! initialize to negative (no split found)
  splitVar = -1

  ! largest critical value
  ! initialize as 0
  maxValueSplit = 0.0

  ! cutoff that gives largest critical value or is randomly selected (OUTPUT)
  cutoffBest = 0.0

  ! set casesOut to 0 (OUTPUT)
  casesOut = 0

  ! set number of cutoffs to 0 (OUTPUT)
  nCuts = 0

  ! indices of variables to be explored, set equal to input: covariates to include in search
  variables = varsIn

  ! location of last available parameter in variables: the number of covariates to include in search
  ptr = nv

  ! tracks the number of variables tried
  variablesTried = 0

  ! creates an array called "tcases"" and initializes it with consecutive integers from 1 to nCases

  tcases = (/(i,i=1,nCases)/)

  ! if this condition is true, execute the RETURN statement which exits the subroutine bc splitting not possible

  IF (nCases < 2*nodeSize) RETURN


    !! iterates from 1 to nv
    !! process is repeated for each covariate in the loop, and at each step, the code explores the possibility of splitting the data based on the selected covariate to find the optimal split for building a decision tree

  DO i = 1, nv

    ! if number of variables tried is equal to maximum allowed, EXIT the loop
    ! if mTry successful splits already explored, exit
    IF (variablesTried .EQ. mTry) EXIT

    ! randomly select an index, then subset covariate for that index on which to split
    ikv = 1 + floor(rnd(0.d0, 1.d0)*ptr)
    kv = variables(ikv)

    ! move this index to the end of the list; shift pointer down

    ! stores the value at the current pointer position (ptr) in the variable j
    j = variables(ptr)

    ! write kv to the current pointer position in the variables array
    variables(ptr) = kv

    ! Moves the temporarily stored value j to the position indicated by the randomly selected index ikv
    variables(ikv) = j

    ! decrease pointer position by 1

    ptr = ptr - 1

    ! pull appropriate covariate. If un-ordered factor use mean survival time
    ! for each factor level

    ! covariate kv is an unordered factor (categorical variable with more than one category)
    IF (nCat(kv) .GT. 1) THEN

    ! if true, calls the "getCovariate"" subroutine to get the appropriate covariate-level survival values (xSorted)
      CALL getCovariate(nCases, casesIn, kv, xSorted)
    ELSE

    ! If false, it directly assigns the covariate values from the full list of caviarate in the samples to be considered to xSorted
    ! kv = covariate under consideration
      xSorted = x(casesIn,kv)
    END IF

    ! initialize cases array with the original set of indices (casesIn)

    cases = casesIn


    ! Calls a built in function qsort4 and hpsort_eps_epw to sort the covariate values (xSorted) and track the corresponding indices (cases)

    ! sort the covariate and track the indices
    CALL qsort4 (xSorted, cases, 1, nCases)
!    CALL hpsort_eps_epw(nCases, xSorted, cases, 1d-8)


    ! sort event indicator data accordingly to the sorted indices in cases array
    dSorted = delta(cases)

    ! ******************** splitBoundaries ********************

    ! each parent node is splitting into two daughter nodes: a left node and a right node
    ! decision to assign a case to the left or right node is based on a splitting criterion

    ! Left Node: The left node contains cases that satisfy the condition specified by the splitting criterion
    ! Right Node: cases that DO NOT satisfy the condition specified by the splitting criterion
    ! splitting criteria based on identifying the minimum cases for left and right splits,



    ! identify minimum cases for left and right splits based on
    ! minimum uncensored cases, minimum node size, and assurance
    ! that all equal valued cases are included in the minimum nodes

    ! initialize variables rUnif and rUnifSet. rUnif will later represent a randomly chosen split point.
    ! Set rUnifSet to -1 to indicate split point hasn't been chosen yet

rUnif = 0.d0
rUnifSet = -1

! cases that are not-censored

! selecting indices of cases that are not censored (where dSorted .EQ. 1)
uncensoredIndices = pack(tcases, dSorted .EQ. 1)

! stores the number of uncensored cases
nUncensored = size(uncensoredIndices)



! if too few cases to meet minimum number of uncensored cases, CYCLE to the next covariate
! this condition checks if there are enough cases to proceed with a split

IF (nUncensored .LT. (minEvent * 2)) CYCLE

!! able to split and satisfy minimum number of events in each node

! cases to left include all indices up to and including minEvent case
! must have at least nodeSize cases

! maximum of the indices corresponding to the minEvent-th uncensored case or nodeSize
! ensures that the left split includes the minimum number of events
splitLeft = max(uncensoredIndices(minEvent), nodeSize)

! move splitLeft up to include cases with equivalent values of x
! ensures that splitLeft is moved up to include all cases with the same x value, with a tolerance for equivalence
splitLeft = count(xSorted .LE. (xSorted(splitLeft) + 1e-8))

! cases to right
! include all indices down to and including nUncensored - minEvent + 1 case
! must have at least nodeSize cases

! rightNode is set to the minimum of indices corresponding to the (nUncensored - minEvent + 1)-th uncensored case or nCases - nodeSize + 1
! ensures that the right split has at least nodeSize cases and includes the minimum required number of events

rightNode = min(uncensoredIndices(nUncensored - minEvent + 1), &
                  & nCases - nodeSize + 1)

! move rightNode down to include cases with equivalent values of x
! splitLeftFinal is the last possible case for the left node
! splitLeftFinal is adjusted to include all cases with covariate values less than those of the selected right split point

splitLeftFinal = count(xSorted .LT. xSorted(rightNode))

! if the splitLeft index is above the splitLeftFinal index cycle,
! split is not possible, so cycle to the next covariate
IF (splitLeft .GT. splitLeftFinal) CYCLE

! Reset rUnifSet to 0
rUnifSet = 0

! if randomSplit (input) is FALSE and ERT is used instead, proceed with further computation
! ERT is global variable

IF ((.NOT. randomSplit) .AND. (ERT .EQ. 1)) THEN

!************* getUniformSplit *****************

  !! split based on extremely randomized trees and uniform split inputs
! ERT: ensemble learning method of decision trees where splits are chosen randomly
! random split points add additional level of randomness
! regular RF splits chosen from range of values in subset of features

! if ERT splitLeft = splitLeftFinal, which is the splitting point
! and cutoff is set to random value or mid-point depending on uniERT

!! extremely randomized tree methods used

!! if uniformSplit is 0, the cutoff for split isn't determined from uniform distribution

      IF (uniformSplit .EQ. 0) THEN

        ! if the cutoff is not determined from a uniform distribution
        ! randomly sample available indices to identify the last case
        ! of the left split to define the minimum; cutoff is the
        ! mid-point {x[r] + x[r+1]} / 2

        ! only indices for which x[r] < x[r+1] can be sampled

        ! create integer array to generate indices within the range from splitLeft to splitLeftFinal
        ind = (/ (i, i = splitLeft, splitLeftFinal) /)

        ! create logical array by comparing each element of xSorted at indiced by "ind" with corresponding element at next index must be < next value
        singles = xSorted(ind) .LT. xSorted(ind + 1)

        ! select only indices from "ind" where element in "singles" array is TRUE where x[r] < x[r + 1]
        indSingles = pack(ind, singles)

        ! count of cases where x[r] < x[r + 1]
        cnt = size(indSingles)

        splitLeftFinal = indSingles(1 + floor(rnd(0.d0, 1.d0)*cnt))
        ! the last required case in the left node is now splitLeftFinal
        splitLeft = splitLeftFinal

        ! update the rUnif variable with the midpoint of the selected interval
        rUnif = (xSorted(splitLeftFinal) + xSorted(splitLeftFinal+1))/2.0

        ! set this equal to 1 to indicate the cutoff has been determined
        rUnifSet = 1

        ! if uniformSplit = 1, means cutoff is randomly selected within the range of values satisfying the allowed cases in the left and right nodes

      ELSE IF (uniformSplit .EQ. 1) THEN
        ! randomly select a value in the range of values that satisfy the
        ! allowed cases in the left/right nodes

        ! set rUnif cutoff randomly within the range of the left and right split points
        rUnif = rnd(0.d0, 1.d0) * (xSorted(splitLeftFinal+1) - &
              & xSorted(splitLeft)) + xSorted(splitLeft)

              ! set rUnifSet to 1 to indicate a cutoff has been determined
        rUnifSet = 1

        ! identify the first case that splits to the right of this value
        ! the preceding case is the last case to the left node
        splitLeftFinal = nCases - count(xSorted > rUnif)

        ! the last required case in the left node is now splitLeftFinal
        splitLeft = splitLeftFinal

        ! end of (uniformSplit .EQ. 1)

      END IF

      ! end of IF (uniformSplit .EQ. 0)

    END IF

    ! -1 is returned if cannot satisfy minimum requirements for nodes
    ! cycle to next covariate
    ! NOTE: it's initialized earlier as -1

IF (rUnifSet .EQ. -1) CYCLE

! increment the number of covariates that have been explored
variablesTried = variablesTried + 1

!***************** maxValue ***************

! when ERT is used, this section is responsible for calculating the test statistic for each potential split point in covariate "kv"

! set initial values for outputs

! initialize variables to keep track of current maximum value of test stat
set = 0
maxValueXm = 0.d0
cutOff = 0.d0

! separate cases into left and right nodes based on the current split point

leftCases = cases(1:(splitLeft-1))
rightCases = cases(splitLeft:nCases)

! get probabilities of cases in left node ("prl") and right node ("prr")
! also get propensity of cases in left node and right node

prl = pr(leftCases,:)
!propensityLeft = propensity(leftCases,:)
prr = pr(rightCases,:)
!propensityRight = propensity(rightCases,:)



! calculate number of events in left and right nodes based on separated cases
! thisis done by summing the product of probabilities * delta
! we include an adjustment by the propensity score

eventsLeft = sum(prl * &
                 & spread(dSorted(1:(splitLeft-1)), 2, nt), DIM = 1)

                 !! merge(prl/propensityLeft, 0.d0, propensityLeft == 0.d0) * &
                 !! & spread(dSorted(1:(splitLeft-1)), 2, nt)

eventsRight = sum(prr * &
                 & spread(dSorted(splitLeft:nCases), 2, nt), DIM = 1)

                 !! sum(merge(prr/propensityRight, 0.d0, propensityRight == 0.d0) * &
                 !! & spread(dSorted(splitLeft:nCases), 2, nt), DIM = 1)


! the sum of probabilities for the left and right cases
! we adjust these by the propensity score

pd1 = sum( prl, DIM = 1)

! sum( merge(prl/propensityLeft, 0.d0, propensityLeft == 0.d0), DIM = 1)

pd2 = sum( prr, DIM = 1)
! sum( merge(prr/propensityRight, 0.d0, propensityRight == 0.d0), DIM = 1)

! at risk initialized to number of cases for the first time point (j = 1)

atRiskLeft(1) = splitLeft - 1
atRiskRight(1) = nCases - splitLeft + 1

! calculates remaining elements of arRiskLeft and atRiskRight baseed on previous values
! for each time point ("j"), the number at risk in the left and right nodes updated based on cum number of events

! note: this is similar to the "Rb" portion in calcvaluesingle that calculates # at risk cases at each time point
! also note, that we use pd1 and pd2 which are already aggregated
! therefore, we need to adjust the values of prl and prr themselves

DO j = 2, nt
atRiskLeft(j) = atRiskLeft(j-1) - pd1(j-1)
atRiskRight(j) = atRiskRight(j-1) - pd2(j-1)
END DO

! if logrank, do calculations that do not depend on node occupancy
! if the rule == 1, we're using the logrank test for splitting rule
    !! therefore we would call a subroutine called "logrankSetUp"

    IF (rule == 1) THEN
      CALL logrankSetUp(atRiskLeft, atRiskRight, eventsLeft, eventsRight, &
                      & numJ, denJ)
    END IF

    ! loop over each potential split point within range splitLeft to splitLeftFinal
    ! calculate test statistic ("valuej") for each potential split point ("j") based on rule
    ! rule is either logrank or mean split

    ! initialize counter variable to 1
    cnt = 1

    ! start a loop over each potential split point ("j") within the range
    ! splitLeft: minimum # of uncensored cases required in left node
    ! splitLeftFinal: last case that can be included in left node while satisfying min node siqe requirement

    DO j = splitLeft, splitLeftFinal

      ! at risk indicators for jth case
      ! number of events for jth case

      ! at risk indicator for j-th case
      ! note, that we also need to adjust these splits

      pd1 = prr(cnt,:)

      ! merge(prr(cnt,:)/propensityRight(cnt,:), 0.d0, propensityRight(cnt,:) == 0.d0)

      cnt = cnt + 1

      ! number of events for the j-th case
      ! this is adjusted by the propensity score
     D = prr(cnt,:) * delta(cases(j))

     ! D = merge(prr(cnt,:) / propensityRight(cnt,:), 0.d0, propensityRight(cnt,:) == 0.d0 ) * delta(cases(j))

      ! initializes cumulative risk at the first time point to 0

      Rcum(1) = 0.0

      ! set cumulative risk at second time point to prob of event happening at first point
      Rcum(2) = pd1(1)

      ! initialize and compute cumulative sum of "pd1" with conditional check of threshold

      ! starting at time point 3 until the total number of points
      DO k = 3, nt

      ! if probability of event at previous time point is greater than a small threshold
        IF (pd1(k-1) .GT. 1d-8) THEN

        ! then update the cumulative risk at the current time point by adding prob at previous time point
        ! note: these are already adjusted for propensity score, so we don't need to adjust them again
          Rcum(k) = Rcum(k-1) + pd1(k-1)

          ! otherwise the prob isn't greater than threshold so no change in cumulative risk at current time point
ELSE

! update cumulative risk at current time point without adding additional risk
Rcum(k) = Rcum(k-1)
END IF
END DO

! calculate the complement of the cumulative risk t each time point
! gives number at risk at each time point (the proportion)
! note, that each rcum is inherently adjusted by propensity score since the probabilities are propensity adjusted

Rcum = 1.d0 - Rcum
! number at risk

! add the jth case to the left node
! add the cumulative risk to number at risk in the left node to adjust counts of indiv at risk based on prob of survival
atRiskLeft = atRiskLeft + Rcum

! remove the jth case from the right node
! remove the cumulative risk from the number at risk in the right node
atRiskRight = atRiskRight - Rcum

! number of events

! add the jth case to the left node
! add the number of events for the j-th case ("D") to total number of events in left node
! note that D is already propensity adjusted, as are the eventsLeft

eventsLeft = eventsLeft + D

! remove the jth case from the right node
! subtract number of events for j-th case ("D") from number of events in right node
eventsRight = eventsRight - D

! if the potential split point ("j") is not the last case with this covariate value, cycle
! by cycling, we examine the next case with the same covariate value
IF (xSorted(j) .GE. (xSorted(j+1) - 1d-8)) CYCLE

! calculate test statistic
! if the splitting rule is the log-rank test,
IF (rule == 1) THEN

! call a subroutine to calculate the log-rank test statistic based on updated number at risk and events in left and right nodes
CALL logrank(atRiskLeft, atRiskRight, eventsLeft, numJ, &
               & denJ, valuej)
ELSE

! otherwise, call a subroutine called "meanSplit" to calculate TS based on mean split rule
! we note that these inputs are all propensity adjusted

CALL meanSplit(atRiskLeft, atRiskRight, eventsLeft, eventsRight, valuej)
END IF

! if "set" = 0 AKA it's the first value being processed in the loop
      ! or if the current test statistic ("valuej") is greater than current maximum TS, enter the block

      IF ((set .EQ. 0) .OR. (valuej .GT. maxValueXm)) THEN

      ! if a uniform split is used AKA rUnifSet == 1


        ! if first value or value > current max, save
        IF (rUnifSet .EQ. 1) THEN

        ! set the cutoff to the predetermined random value "rUnif"
          cutoff = rUnif
        ELSE

        ! otherwise, set the cutoff equal to the midpoint between the current and next covar values
          cutoff = (xSorted(j) + xSorted(j+1))/2.d0

          ! end if for if rUnifSet == 1
        END IF

        ! update the maximum value and corresponding cutoff if current TS greater than current max

         ! update the maxValueXm to the current test statistic "valuej"
        maxValueXm = valuej

        ! reset the tie value to 1
        tieValue = 1

        ! set == 1 to indicate a valid split has been found
        set = 1

        ! else, if "set" isn't 0 or valuej < maxValueXm), and the current stat is close to the max (with tolerance), enter the block

ELSE IF (valuej > (maxValueXm - 1d-8)) THEN
! if value is a tie, randomly determine if cutoff should be taken
! tieValue is incremented
tieValue = tieValue + 1

! generate a random number between 0, 1
! if this is < the reciprocal of the current tie count, choose new cutoff as midpoint between current and next covariate values


IF (rnd(0.d0, 1.d0) < (1.d0 / REAL(tieValue))) THEN
cutoff = (xSorted(j) + xSorted(j+1))/2.d0
END IF

END IF

! end cycling across each time point

END DO

! if the loop was not able to find a split point, cycle to the next covariate

! if not successful, cycle to next covariate
! this condition should never be true
IF (set .EQ. 0) CYCLE

! if successful, determine if it yields the maximum value of the
! covariates considered

! if current split point has higher maximum value than previously considered splits,
! or if it's the first non-0 maximum (splitVar == -1), enter the code block

    IF ((splitVar .EQ. -1) .OR. (maxValueXm .GT. maxValueSplit)) THEN

      ! if first non-zero or largest value, keep cutoff and value and
      ! reset tie counter to 1


        ! if current split point yields higher maximum val or is the first non-0 max, update best split info

        ! update the best split information by settinv variable index to current variable index
        ! maximum value to current maximum value ("maxValueXm")
        ! reset tie counter "tieCovariate" to 1

      splitVar = kv
      maxValueSplit = maxValueXm
      tieCovariate = 1



      ! count the number of cases in the left node,assigns variable "cases" to "casesOut"
      lft = count(xSorted .LE. cutoff)
      casesOut = cases

      ! initializes "cutoffBest" variable to 0

      cutoffBest = 0

      ! if number of categories for current variable "kv" is <= 1, set first element to current cutoff
      ! also set nCuts to 1

      IF (nCat(kv) .LE. 1) THEN
        cutoffBest(1) = cutoff
        nCuts = 1

        ! otherwise, proceed to next block of code
      ELSE


        ! for factors with more than 1 category, identify factor values contained in left cases
        l = 1

        ! iterates over the categories "j"
        DO j = 1, nCat(kv)

        ! iterates over the left cases "jj"
          DO jj = 1, lft

          ! identifies factor values contained in left cases and updates the "cutoffBest" array accordingly

            IF (nint(x(casesOut(jj),kv)) .EQ. j) THEN
              cutoffBest(l) = j
              l = l + 1
              EXIT
            END IF
          END DO
        END DO

        ! nCuts variable set to number of identified factor values

        nCuts = l - 1
      END IF

      ! if a random split was triggered, break out of loop over covariates
      IF (randomSplit) EXIT

      ! if this condition isn't met, ((splitVar .EQ. -1) .OR. (maxValueXm .GT. maxValueSplit))
! and the value of maxValueXm is larger than the previous minus a small tolerance then these values are basically the same

ELSE IF (maxValueXm .GT. (maxValueSplit-1d-8)) THEN

! if equal to current maximum value, increment tie counter and randomly
! select the cutoff with equal probability for each tie
tieCovariate = tieCovariate + 1

! update the split variable, maximum split value

IF (rnd(0.d0, 1.d0) < (1.0 / tieCovariate)) THEN
splitVar = kv
maxValueSplit = maxValueXm
! count the number of cases in the left node
lft = count(xSorted .LE. cutoff)
casesOut = cases

cutoffBest = 0

IF (nCat(kv) .LE. 1) THEN
cutoffBest(1) = cutoff
nCuts = 1
ELSE
! for factors, identify factor values contained in left cases

! for categorical variables with more than one category, identify factor values contained in left cases and updates "cutoffBest" array

l = 1
DO j = 1, nCat(kv)
DO jj = 1, lft
IF (nint(x(casesOut(jj),kv)) .EQ. j) THEN
cutoffBest(l) = j
l = l + 1
EXIT
END IF
END DO
END DO

! set nCuts based on the identified factor values

nCuts = l - 1
END IF
END IF
END IF

! end the loop over covariates
END DO

! if no split was possible,exit the subroutine by using "RETURN"
if (splitVar .EQ. -1) RETURN

! if successful at finding a split set flag and return
splitFound = 1

RETURN

END SUBROUTINE tfindSplit

  ! ******************************* subroutine: kaplan ******************************** !
  !** used in subroutine meanSplit, calcValueSingle


!! create a new subroutine called "kaplan" for the KM estimator

! Calculate the Kaplan Meier estimator
! ns integer, the number of time points
! nj real(:), at risk
! oj real(:), events
! z real(:), estimator


!! takes 3 inputs: number of time points, number at risk at each time point, number of events (failures) at each point,

! propensity at the end was deleted
SUBROUTINE kaplan(ns, nj, oj, z)
  IMPLICIT NONE

  ! array representing # of time points
  INTEGER, INTENT(IN) :: ns

  ! array representing number at risk at each time point (already weighted by propensity score)
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: nj

  ! number of events at each time point (already weighted by propensity score)
  REAL(dp), DIMENSION(1:ns), INTENT(IN) :: oj

  ! array where the Kaplan-Meier estimator will be stored
  REAL(dp), DIMENSION(1:ns), INTENT(OUT) :: z


  ! local variable used as a loop index
  INTEGER :: i

  ! local arrays for weighted at-risk and events
  REAL(dp), DIMENSION(1:ns) :: Y_w, d_w

  ! weight the number at risk and events by the propensity score (delta * propensity)
  Y_w = nj
  d_w = oj

  ! calculate the Kaplan-Meier estimator for the first time point
  z(1) = (Y_w(1) - d_w(1)) / Y_w(1)

  ! if there is only one time point, exit the subroutine
  IF (ns .LT. 2) RETURN

  ! loop through each of the remaining time points
  DO i = 2, ns
    IF (Y_w(i) > 1d-8) THEN
      z(i) = ((Y_w(i) - d_w(i)) / Y_w(i)) * z(i-1)
    ELSE
      z(i) = z(i-1)
    END IF

    ! ensure no negative survival probabilities
    IF (z(i) < 0) THEN
      z(i) = 0
    END IF
  END DO

END SUBROUTINE

  ! ************************** subroutine: meanSplit ************************************** !
  ! ** this is used in subroutine "tfindSplit"
  ! ** there, if log-rank test isn't used (rule !=1), use the meanSplit rule

  !!!!!!! create a new subroutine called "meanSplit"

! Truncated Mean test
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! Z: real, truncated mean

! takes as iput: N1j, N2j, O1j, O2j representing at-risk and event counts for G1 and G2


SUBROUTINE meanSplit(N1j, N2j, O1j, O2j, Z)
  ! enforces explicit declaration of all variables, preventing use of implicitly declared vars

  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j

  !!!!! output variable: is truncated mean test statistic

  REAL(dp), INTENT(OUT) :: Z

  !! initiate arrays to store KM estimates

  REAL(dp), DIMENSION(1:nt) :: E1, E2

  !! call "kaplan" subroutine to calculate KM estimates for group 1
  ! deleted "propensity1" input
  CALL kaplan(nt, N1j, O1j, E1)

  !! call "kaplan" subroutine to call KM estimates for group 2

  ! deleted "propensity2" input
  CALL kaplan(nt, N2j, O2j, E2)

  !! calculate the truncated mean test-statistic:
  ! sum of (difference between KM estimates *time_intervals)
  ! note: dt is global variable

  Z = sum((E1 - E2) * dt)

  ! squares the TS
  Z = Z*Z

END SUBROUTINE

  ! ******************************** subroutine: logRankSetUp ******************************** !
   !** used in subroutine tfindSplit

! Log rank test set up
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! numJ: real(:), numerator
! denJ: real(:), denominator

! takes inputs N1j, N2j, O1j, O2j
! provides outputs numJ (numerator), denJ (denominator)

SUBROUTINE logRankSetUp(N1j, N2j, O1j, O2j, numJ, denJ)

! all variables must be declared

  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O2j

  !! output variables: numerator and denominator

  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: numJ
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: denJ

  ! local variable used as a loop index

  INTEGER :: i

  REAL(dp) :: Nj, Oj

    ! initialize the numerator and denominator to 0

  numJ = 0.d0
  denJ = 0.d0

  ! loop through each of the time points (nt) is total number of time points

  DO i = 1, nt

  ! if there are 0 individuals at risk in group 1 or group2 (based on a certain tolerance),
  ! CYCLE to skip iterations of the do loop

    IF (N1j(i) .LT. 1d-8) CYCLE
    IF (N2j(i) .LT. 1d-8) CYCLE

    ! time points for which both events have individuals at risk

    ! calculate the sum of individuals at risk in the two groups
    ! number of individuals at risk for type 1 or 2 events
    Nj = N1j(i) + N2j(i)

    ! calculate the number of events in both groups
    ! number of events of type 1 or type 2
    Oj = O1j(i) + O2j(i)

    ! calculate the numerator at the i-th time point-- ratio of events / total at risk at ith time point

    numJ(i) = Oj / Nj

    ! contribution to the denominator of the log-rank test statistic at ith time point

    denJ(i) = numJ(i) * (Nj - Oj) / (Nj * Nj)

  END DO

END SUBROUTINE logrankSetUp


 ! ******************************** subroutine: logRank ******************************** !
 !** used in subroutine tfindSplit

! Log rank test
! N1j: real(:), at risk in group 1
! N2j: real(:), at risk in group 2
! O1j: real(:), events in group 1
! O2j: real(:), events in group 2
! numJ: real(:), numerator
! denJ: real(:), denominator
! Z: real, test value

! takes several inputs, and outputs a scalar called "Z" which is the test statistic

SUBROUTINE logRank(N1j, N2j, O1j, numJ, denJ, Z)

  ! ensures all variables must be explicitly declared

  IMPLICIT NONE

  ! declares input and output variables and specifies types, dimensions, intent

  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: N2j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: O1j
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: numJ
  REAL(dp), DIMENSION(1:nt), INTENT(IN) :: denJ
  REAL(dp), INTENT(OUT) :: Z

  ! declare an internal loop variable to be used

  INTEGER :: i

  ! declare local variables "den" and "num" to store numerator and denominator of TS

  REAL(dp) :: den, num

  ! initialize the numerator and denom to 0

  num = 0.d0
  den = 0.d0

  ! loop over each of the time points
  DO i = 1, nt

  ! if either N1j(i) or N2j(i) is effectively 0 (no one at risk in either group), CYCLE to next iteration of time
    IF (N1j(i) .LT. 1d-8) CYCLE
    IF (N2j(i) .LT. 1d-8) CYCLE
    ! time points for which both events have individuals at risk

    ! accumulate numerator and denom values for log-rank TS

    num = num + O1j(i) - numJ(i) * N1j(i)
    den = den + denJ(i) * N1j(i) * N2j(i)

  END DO

  ! if the denom is greater than 0 (with some tolerance), set

  IF (den .GT. 1d-8) THEN

  ! calculate LR statistic and set Z to num^2 / denom
    Z = num * num / den
  ELSE

  ! otherwise, set Z to 0 to avoid division by a small number
    Z = 0.d0
  END IF

END SUBROUTINE

  ! *********************************** subroutine: calcValueSingle ***************************** !
  !** used in subroutine getCovariate

! estimate the survival function and mean survival time
!   nCases: integer, the number of elements in casesIn
!   casesIn: integer, the indices of the subset for which the value is
!     calculated

!! outputs: estimated survival function and estimated mean survival time
!   survFunc: real(:), the estimated survival function
!   mean: real, the estimated mean survival
! deleted input "propensity"
SUBROUTINE calcValueSingle(nCases, casesIn,survFunc, mean)

  ! ensures all variables must be explicitly declared

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn
  REAL(dp), DIMENSION(1:nt), INTENT(OUT) :: survFunc
  REAL(dp), INTENT(OUT) :: mean

  ! declare loop counter "i" and arrays of Nj, Ok, Rb

  INTEGER :: i

  REAL(dp), DIMENSION(1:nt) :: Nj, Oj, Rb
  ! initialize propensity score input

  ! initialize values of "survFunc" and "mean" as 0

  survFunc = 0.d0
  mean = 0.d0

  ! number of at risk cases at each time point
  ! {nt}

  ! sums the elements of the "pr" array along the first dimension AKA cumulative number of at-risk cases
          !! this is now adjusted by propensity score (each case is a probability, and each adjusted by prop score)
  ! select all the time points of the pr array for the rows specified by "casesIn"
  ! extracts the subset of probabilities corresponding to the chosen cases
  ! NOTE: "pr" is not defined in this subroutine
  ! NOTE: "propensity" is also not defined in this subroutine, but pr and propensity are the same size

  ! Calculate Rb while accounting for zeros in the propensity
  Rb = sum(pr(casesIn, :), DIM = 1)

  ! sum(merge(pr(casesIn, :) / propensity(casesIn, :), 0.d0, propensity(casesIn, :) == 0.d0), DIM = 1)


    !!!!!!!
    !!!!!!! Adding print statement
    !!!!!!!

    !! Print values of variables
    !print *, "Number of at risk cases at each time point:"
    !print *, "Rb =", Rb

    !!!!!!!
    !!!!!!!
    !!!!!!!

  ! initialize Nj(1) with the number of cases at the first time point

  Nj(1) = nCases

  ! loop from the next time point until the last
  DO i = 2, nt

  ! calculate number of cases at risk at each subsequent time point
  ! subtract cumulative at-risk cases form the previous total
  ! this doesn't need to be weighted by the propensity score, since the Rb part is weighted

    Nj(i) = Nj(i-1) - Rb(i-1)

  END DO

  ! number of events at each time point
  ! {nt}

  ! looping through each of the time points,
  DO i = 1, nt

  ! calculate the number of events by summing (delta * probability for each indiv at that time point)
  ! we need to adjust the probability of the number of events by the corresponding propensity
    Oj(i) = sum((pr(casesIn, i)*delta(casesIn)))

    ! sum(merge(((pr(casesIn, i)*delta(casesIn))/propensity(casesIn, i)), 0.d0, propensity(casesIn, i) == 0.d0))


  END DO

    !!!!!!!
    !!!!!!! Adding print statement
    !!!!!!!

    ! Print values of variables
    !print *, "Number of at events at each time point calcvaluesingle:"
    !print *, "Oj =", Oj
    !print *, "Number of cases at risk at each time point calcvaluesingle:"
    !print *, "Nj =", Nj


    !!!!!!!
    !!!!!!!
    !!!!!!!

  ! Kaplan-Meier estimate survival function
  ! {nt}

  !! call the "kaplan" subroutine to calculate the KM estimate of survival func
  ! these are stored in the "survFunc" array
  ! nt: number of time points (ns in "kaplan")-- dimension 1:ns
  ! Nj: number at risk-- dimension 1:ns
  ! Oj: number of events-- dimension 1:ns
  ! survFunc: array where the KM estimator is stored ("z")
  ! we rely on "kaplan" to handle the propensity score calculations individually, so we don't need to weight Nj or Oj
! now, the input values into the KM function are already weighted by the propensity score
  CALL kaplan(nt, Nj, Oj, survFunc)

  ! mean survival time

  ! calculate the mean survival time by summing (estimated survFunc * dt)
  ! where dt are the time intervals (global variable)

  mean = sum(survFunc * dt)

  RETURN
END SUBROUTINE

  ! ********************************** subroutine: getCovariate ****************************** !
  !** used in subroutine tfindSplit

! For factor covariates, calculated the mean survival time and
! use as covariate values for splitting
! nCases integer, number of cases under consideration
! casesIn, integer(:), cases under consideration
! propensity
! kv, integer, covariate under consideration
! array, real(:), covariate vector as mean survival times

! delete input "propensity"
SUBROUTINE getCovariate(nCases, casesIn, kv, array)

! ensures all variables must be explicitly declared
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nCases
  INTEGER, DIMENSION(1:nCases), INTENT(IN) :: casesIn

  ! the covariate being considered
  INTEGER, INTENT(IN) :: kv
  REAL(dp), DIMENSION(1:nCases), INTENT(OUT) :: array

  INTEGER :: i
  INTEGER, DIMENSION(1:nCases) :: ix
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind

  REAL(dp) :: mean
  REAL(dp), DIMENSION(1:nt) :: survFunc

  LOGICAL, DIMENSION(1:nCases) :: inSubset

  ! initialize array to 0

  array = 0.d0

  ! convert covariate to integers to ensure correct equality tests
  ix = nint(x(casesIn,kv))

  ! loop over all the levels of the covariates
  !! nCat: number of categories in the covariate... I think this is global param

  DO i = 1, nCat(kv)

    ! identify cases that have the current covariate level
    inSubset = ix .EQ. i

    ! ensure that there are individuals in the lth level

    ! if there are no individuals in the current covariate level, skip to the next level
    IF (count(inSubset) .EQ. 0) CYCLE

    ! pack indices and estimate survival

    ! packs indices of cases in the current subset
    ! select only elements from "casesIn" where the  element "inSubet is true"
    ! AKA select only individuals who have a covariate level matching the index

    ind = pack(casesIn, inSubset)


    ! call the "calcValueSingle" subroutine to calculate the mean survival time for the current subset

    ! calculate the mean survival time for each individual in this subset
    ! AKA used to calculate mean survival tiem for cases sharing the same level of a particular covariate

    ! use size(ind) = nCases: as number of elements in indizes of subset
    ! ind = casesIn: indices of subset for which value is calculated
    ! outputs are "survFunc" and "mean"
  ! delete input "subsetPropensity"
    CALL calcValueSingle(size(ind), ind, survFunc, mean)

    ! update the "array" with the mean survival time for individuals in current subset

    WHERE (inSubset) array = mean

    ! end the loop over covariate levels

  END DO

END SUBROUTINE



  ! ******************************* subroutine: tsurvTree ********************************* !

!! grow each tree
!! output arrays for forest-level survival functions

! forestSurvFunc, real(:), survival function averaged over forest
! forestMean, real, mean survival averaged over forest
! forestSurvProb, real, survival probability averaged over forest


SUBROUTINE tsurvTree(forestSurvFunc, forestMean, forestSurvProb)
  IMPLICIT NONE

  REAL(dp), DIMENSION(1:nAll*nt), INTENT(OUT) :: forestSurvFunc
  REAL(dp), DIMENSION(1:nAll), INTENT(OUT) :: forestMean
  REAL(dp), DIMENSION(1:nAll), INTENT(OUT) :: forestSurvProb

  ! declare loop counters "i, iTree, j, k lft, m, nc, ncur, splitFound, splitVar"

  INTEGER :: i, iTree, j, k, lft, m, nc, ncur, splitFound, splitVar

  ! declare arrays to store indices

  INTEGER, DIMENSION(1:sampleSize) :: indices, jdex, xrand

  ! declare other arrays to store information

  INTEGER, DIMENSION(1:np) :: newstat, pindices
  INTEGER, DIMENSION(1:np, 1:nrNodes) :: cstat
  INTEGER, DIMENSION(1:nrNodes, 1:2) :: stm
  INTEGER, DIMENSION(1:nAll) :: allStatus
  INTEGER, DIMENSION(:), ALLOCATABLE :: ind, indOut, leftCases, rightCases, pind

  ! declre real variables to store values

  REAL(dp) :: srs
  REAL(dp), DIMENSION(1:nLevs) :: cutoffBest
  REAL(dp), DIMENSION(1:nrNodes) :: mean, survProb
  REAL(dp), DIMENSION(1:nAll) :: xm
  REAL(dp), DIMENSION(1:nt, 1:nrNodes) :: survFunc
  REAL(dp), DIMENSION(1:nrNodes, 1:(5+nLevs)) :: nMatrix
  REAL(dp), DIMENSION(1:nt, 1:nAll) :: tforestSurvFunc

  ! declare logical arrays to handle conditions for tree building

  LOGICAL, DIMENSION(1:np) :: cand
  LOGICAL, DIMENSION(1:nAll) :: tst

  ! initialize the output arrays to 0

  tforestSurvFunc = 0.d0
  forestMean = 0.d0
  forestSurvProb = 0.d0

  ! iterate over each tree in the forest
  ! nTree isn't an input parameter but it's probably called later with this subroutine

  DO iTree = 1, nTree

  ! initialize these variable values to a default value of 0 or 1 (they were declared above)

    survFunc = 0.d0
    mean = 0.d0
    survProb = 0.d0
    nMatrix = 0.0
    allStatus = 1

    ! sample data and set local variables x, pr, and delta to the selected
    ! subset

    !! if replacement is allowed == 1, (replace is not an input)

    IF (replace .EQ. 1) THEN

    !! then sample the data (nAll = number of cases) with a specified sampleSize

      xrand = sampleWithReplace(nAll, sampleSize)

      ! assign sample size to "n"
      n = sampleSize

      ! select the covariates corresponding to the random pt index
      x = xAll(xrand,:)

      ! select the probabilities corresponding to the random pt index
      pr = prAll(xrand,:)

      ! select the delta corresponding to the random pt index
      delta = deltaAll(xrand)

      ! select the propensity score corresponding to the random pt index
      propensity = propensityAll(xrand)



      !! NOTE: xAll, prAll, deltaAll are not inputs here, but they will be defined when this subroutine is called later

      ! if replacement is not allowed AND total number of cases isn't equal to samp size,

    ELSE IF (nAll .NE. sampleSize) THEN

    ! sample witout replacement

      xrand = sampleWithoutReplace(nAll, sampleSize)
      n = sampleSize
      x = xAll(xrand,:)
      pr = prAll(xrand,:)
      delta = deltaAll(xrand)
      propensity = propensityAll(xrand)

      ! if replacement not allowed AND total cases = samp size,

    ELSE

    ! use the entire dataset

      n = sampleSize
      x = xAll
      pr = prAll
      delta = deltaAll
      propensity = propensityAll

      ! end if statement about sampling with/without relacement
    END IF

    ! cutoff for identifying covariates to be explored

    ! calculates a cutoffvalue "srs" by dividing "stratifiedSplit" parameter by number of covars
    ! stratifiedSplit is a global parameter which will be input into dtrSurv


    ! stratifiedSplit = stratified random split coefficient

    srs = stratifiedSplit / REAL(np)

    ! indices for all cases

    ! initialize arrays to store indices for cases
    indices = (/(i,i=1,n)/)
    jdex = indices

    ! indices for all covariates

    ! initialize array to store indices for covariates
    pindices = (/(i,i=1,np)/)

    !! initialize first node

    ! calculate survival function and mean survival of the node


      ! call the "calcValueSingle" subroutine which returns an estimated survival function and mean survival time as output
    ! n = nCases: as number of elements in indices of subset (number of cases)
    ! indices = casesIn: indices of subset for which value is calculated (number of events)


    ! ---------------- to do: input the propensity score (same dimensions as delta) ------------------ !
    CALL calcValueSingle(n, indices, survFunc(:,1), mean(1))

    !!!!!!!
    !!!!!!! Adding print statement
    !!!!!!!

    ! Print values of variables
    !print *, "calcValueSingle for first node:"
    !print *, "surv function and mean surv time =", survFunc(:,1)

    !!!!!!!
    !!!!!!!
    !!!!!!!


    ! if the "isSurvival" flag is TRUE,

    IF (isSurvival) THEN

    ! calculate estimated survival probability using linear interpolation
    ! after getting KM estimate from "calcValueSingle" can get values between observed times
    ! store this result in survProb(1)
            ! sIndex is an input from setUpBasics
            ! this is input into dtrSurv.R for mean survival and mean survival probability
            ! the index of the last time point <= survival time
        ! sFraction (setUpBasics): the fractional location of survival time between the sIndex and the next index time point


      ! estimate survival probability at SurvivalTime
      survProb(1) = survFunc(sIndex,1) * (1.d0 - sFraction) + &
                  & survFunc(sIndex+1,1) * sFraction

      ! if the calculated survival probability is less than a small threshold, set to 0
      IF (survProb(1) .LT. 1d-8) survProb(1) = 0.d0

      ! end IF statement using "isSurvival" flag
    END IF

    ! determine if the node can split based on basic minimum requirements

    !if the number of cases ("n") is <= "nodeSize": minimum number of cases in a node (global) OR
    ! if there are no uncensored events, mark node as terminal

    if (n .LE. nodeSize .OR. sum(delta(indices)) .LE. 1) THEN

    ! mark node as terminal AKA it doesn't split further
      nMatrix(1,1) = -1
    ELSE

    ! otherwise, requirements are met, mark node as having potential to split
      nMatrix(1,1) = -2

      ! end IF statement for n < nodeSize or all pts are censored
    END IF

    ! initialize an array "cstat" with 0s; nmber of times a variable was used in split

    cstat(:,1) = 0

    ! start and finish locations of indices in node

    ! "stm" array used to define start and finish locations in current node
    ! first index is node number; second index has 2 elements: start and finish locations of indices for cases in that node
    !! (1, 1) == 1 set start location of first node to 1
    !! (1, 2) == n set finish location of first node to n (number of cases in the node)
    stm(1,1) = 1
    stm(1,2) = n

    ! location of most recent storage location in matrices/vectors
    ! ncur is incremented when a node successfully splits indicating the
    ! location in the nodes list where base information for the each daughter
    ! is stored

    ncur = 1

    ! (still in loop through all the trees)
    ! loop through each of the nodes in the decision tree
    ! nrNodes = maximum number of nodes in a tree

    DO k = 1, nrNodes
      ! if k is beyond current node count, means loop is iterating over nodes that havent been created yet
      ! current node count at limit, break from loop

      ! if the iteration is < current node count OR
      ! if the current node count > total nodes - 2 (means no more nodes left to split)
      ! EXIT the loop & proceed to the next part of the code after the loop

      IF (k .GT. ncur .OR. ncur .GT. (nrNodes - 2)) EXIT

      ! if the current node can't be split anymore (-1) value, cycle to the next node count

      ! if node is not to be split, cycle to next node
      IF (nint(nMatrix(k,1)) .EQ. -1) CYCLE



      ! indices for cases contained in current node
      ! stm(k, 1) and stm(k, 2) are the start and finish locations of indices in a node
      ! then take the array of indices representing all cases
      ind = jdex(stm(k,1):stm(k,2))

    !!!!!!!
    !!!!!!! Adding print statement for indices in current node
    !!!!!!!

    ! Print values of variables
    !print *, "Indices for cases in current node:"
    !print *, "indices =", ind

    !!!!!!!
    !!!!!!!
    !!!!!!!

      ! compares count of times a variable ws used in a split cstat(:,k) with a threshold floor(srs * sum(cStat(:,k))
      ! cstat(:, k) = column of array for current node "k"
      ! represents count of how many times each varaible has been used in a split for the entire tree we are in up to the current time point (k-th node)
      ! holds a logical to identify deficient variables
      ! filters out variables that have been used more than a certain threshold

      ! if there are deficient variables, use only these variables
      cand = cstat(:,k) .LT. floor(srs * sum(cStat(:,k)))

      ! get indices of covariates used on the split based on deficient variables

      pind = pack(pindices, cand)

      ! if there are no deficient variables, then use all of the variables
      IF (size(pind) .EQ. 0) pind = pindices

      ! split cases
      indOut = ind

      ! call subroutine "tfindSplit" to find the best split for the current node based on the input variables
      ! might be all of them if there are no deficient variables, otherwise only use deficient variables
      ! use of deficient variables is to introduce diversity among trees in the forest
      ! you don't want all nodes to be constructed using the same variables at each split

      !! input parameters:
!   size(ind) = nCases : integer, the number of elements in input casesIn
!   ind = casesIn : integer(:), the indices of the cases in this node
!   size(pind) nv : integer, the number of covariates to include in search
!   pind: varsIn : integer(:), covariates to include in search

!! outputs
!   splitVar : integer, the index of the selected variable for splitting
!   cutoffBest : real(:), the cutoff (<= go to 'left' node)
!   splitFound : integer, 0 = no split; 1 = found a split
!   casesOut : integer(:), elements of casesIn that go left; ind if yes, 0
!     otherwise
!   nCuts : integer, the number of cutoff values returned
!   lft : integer, the number of cases in the left node

! we note, that propensity scores are incorporated WITHIN tfindSplit based on the probabilities


      CALL tfindSplit(size(ind), ind, size(pind), pind, splitVar, cutoffBest, &
                    & splitFound, indOut, nc, lft)

      ! if there is no split found, set the current node as a terminal node (-1)

      IF (splitFound .EQ. 0 ) THEN
        ! if no split available, set node k as terminal node
        nMatrix(k,1) = -1

        ! CYCLE to the next node
        CYCLE


      END IF

      ! set the status of node k to indicate it's an interior node with a split
      ! this means that the node is an "internal" node within a tree and has been split into two values
      ! -1: terminal node, cannot be split further
      ! -2: node can be split further but hasn't yet
      ! -3: interior node: has already been split successfuly into left and right child nodes

      ! set node k to be interior (i.e. has split)

      ! nMatrix(:,1) is indicator about node status
      nMatrix(k,1) = -3

      ! records the variable index used for the split

      ! add split information to node
      nMatrix(k,4) = pindices(splitVar)

      ! number of categories for categorical variables
      nMatrix(k,5) = nc

      ! cutoff values for the split
      ! number of columns here depends on variable type and number of categories
      nMatrix(k,6:(6+nc-1)) = cutoffBest(1:nc)

      !! set indices of the left daughter node
      nMatrix(k,2) = ncur + 1

      ! right daughter node

      nMatrix(k,3) = ncur + 2

      ! update the count of how many times a variable was used in a split

      ! increment the times the variable was used in a split by getting current count for kth node
      newstat = cStat(:,k)

      ! then increase this k-th count value
      ! nMatrix(k,4) records the variable index used for the split, so we go to that varaible index and increase the count

      newstat(nint(nMatrix(k,4))) = newstat(nint(nMatrix(k,4))) + 1

      ! update the order of cases in the "jdex" array based on the split
      ! stm(k, 1) and stm(k, 2) are the start and finish locations of indices in a node
      ! then take the array of indices representing all cases
      ! indOut: output of subroutine "tfindSplit": elements of casesIn that go left; ind if yes, 0 otherwise

      ! store new case order in jdex
      jdex(stm(k,1):stm(k,2)) = indOut

      !! left node

      !! increment the current node count (AKA a new daughter node is being added)

      ncur = ncur + 1


      ! index boundaries for cases in left node

      ! set the start and end indices for cases in the left daughter node

      stm(ncur,1) = stm(k,1)

      ! lft = the number of cases in the left node
      ! output from subroutine "tfindSplit"
      stm(ncur,2) = stm(k,1) + lft - 1

      ! extract indixes cases corresponding to the left daughter node from "jdex" array

      leftCases = jdex(stm(ncur,1):stm(ncur,2))

    !!!!!!!
    !!!!!!! Adding print statement for indices in current left node
    !!!!!!!

    ! Print values of variables
    !print *, "Indices for cases in current left node:"
    !print *, "indices =", leftCases

    !!!!!!!
    !!!!!!!
    !!!!!!!


      ! get basic node information for left daughter

      !! call "calcValueSingle" subroutine to calculate survival function & mean survival time for left daughter node


    ! size(leftCases) = nCases: as number of elements in indices of subset (number of cases)
    ! leftCases = casesIn: indices of subset for which value is calculated (number of events)

    ! --------- NOTE: we need to input the selected propensity scores from the correct indices (same number as nCases) ------------ !

      CALL calcValueSingle(size(leftCases), leftCases, survFunc(:,ncur), &
                         & mean(ncur))


    !!!!!!!
    !!!!!!! Adding print statement for indices in surival function for left daughter node
    !!!!!!!

    ! Print values of variables
    !print *, "survival function and mean survival time for left daughter node:"
    !print *, "values =", survFunc(:,ncur)

    !!!!!!!
    !!!!!!!
    !!!!!!!

      ! if isSurvival is true,

      IF (isSurvival) THEN

      ! use linear interpolation to calculate survival probability of the left node between two adjacent time points

        ! estimate survival probability at SurvivalTime
        ! we get a KM estimate of the survival curve, but this might be to estimate survival probabilities that falls between two time events
        ! recall KM estimator is a step function jumping when events occur
                ! sIndex is an input from setUpBasics
            ! this is input into dtrSurv.R for mean survival and mean survival probability
            ! the index of the last time point <= survival time
        ! sFraction (setUpBasics): the fractional location of survival time between the sIndex and the next index time point


        survProb(ncur) = survFunc(sIndex,ncur) * (1.d0 - sFraction) + &

                       & survFunc(sIndex+1,ncur) * sFraction

        ! if the survival probability is very close to 0 with a certain tolerance, set it equal to 0
        IF (survProb(ncur) .LT. 1d-8) survProb(1) = 0.d0

        ! end of using "isSurvival" flag
      END IF


      ! if the number of cases in left node <= node size (minimum number of cases in a node (global) OR
      ! if there are no uncensored cases

      IF (size(leftCases) .LE. nodeSize .OR. sum(delta(leftCases)) .LE. 1) THEN
        ! if the number of cases in the node is at or below the minimum required
        ! or the number of uncensored event is only 1
        ! status is terminal

        ! mark the current daughter node as terminal
        nMatrix(ncur,1) = -1
      ELSE

        ! otherwise, mark it as -2 to indicate it has the potential to be split (but hasn't been split yet)
        nMatrix(ncur,1) = -2

        ! end statement about minimum node size
      END IF

      ! copies update count of variable usage "newstat" to current node's ("ncur") column in variable usage count ("cstat")
      ! recall: this is a count across all nodes in a forest

      cstat(:,ncur) = newstat

      !!!!!! ------- right node ------- !!!!!!

      ! increment the node counter for the right node

      ncur = ncur + 1

      ! index boundaries for cases in right node

      ! set starting index for cases in right node to be starting index in parent node + lft
      ! lft = the number of cases in the left node

      stm(ncur,1) = stm(k,1) + lft

      ! ending index for cases in the right node same as ending index of the parent node
      stm(ncur,2) = stm(k,2)


      ! retrieve left and right cases
      ! extract indixes cases corresponding to the right daughter node from "jdex" array

      rightCases = jdex(stm(ncur,1):stm(ncur,2))


    !!!!!!!
    !!!!!!! Adding print statement for indices in right daughter node
    !!!!!!!

    ! Print values of variables
    !print *, "indices of cases in right daughter node:"
    !print *, "indices =", rightCases

    !!!!!!!
    !!!!!!!
    !!!!!!!


      ! delete input propensityright since calcvaluesingle already adjusts for propensity scores
      ! calculate survival function and mean survival time for right node


          ! --------- NOTE: we need to input the selected propensity scores from the correct indices (same number as nCases) ------------ !


      CALL calcValueSingle(size(rightCases), rightCases, survFunc(:,ncur), &
                         & mean(ncur))


    !!!!!!!
    !!!!!!! Adding print statement for survival function in right node
    !!!!!!!

    ! Print values of variables
    !print *, "survival function and mean survival in right daughter node:"
    !print *, "values =", survFunc(:,ncur)

    !!!!!!!
    !!!!!!!
    !!!!!!!

      ! if "isSurvival" is TRUE

      IF (isSurvival) THEN

        ! use linear interpolation to calculate survival probability of the left node between two adjacent time points
        ! estimate survival probability at SurvivalTime
        ! we get a KM estimate of the survival curve, but this might be to estimate survival probabilities that falls between two time events
        ! recall KM estimator is a step function jumping when events occur
        ! sIndex is an input from setUpBasics
            ! this is input into dtrSurv.R for mean survival and mean survival probability
            ! the index of the last time point <= survival time
        ! sFraction (setUpBasics): the fractional location of survival time between the sIndex and the next index time point

        survProb(ncur) = survFunc(sIndex,ncur) * (1.d0 - sFraction) + &
                       & survFunc(sIndex+1,ncur) * sFraction

        ! if the value is very close to 0, set it equal to 0

        IF (survProb(ncur) .LT. 1d-8) survProb(1) = 0.d0

        ! end the IF about using "isSurvival"
      END IF


      ! if the number of cases in right node <= node size (minimum number of cases in a node (global) OR
      ! if there are no uncensored cases

      IF (size(rightCases) .LE. nodeSize .OR. &
        & sum(delta(rightCases)) .LE. 1) THEN
        ! if the number of cases in the node is at or below the minimum required
        ! or the number of uncensored event is only 1
        ! status is terminal

        ! mark the current daughter node as terminal
        nMatrix(ncur,1) = -1

      ELSE

        ! otherwise, mark the node with potential to be split further

        nMatrix(ncur,1) = -2

        ! end IF statement about minimum node size
      END IF

      ! copies update count of variable usage "newstat" to current node's ("ncur") column in variable usage count ("cstat")
      ! recall: this is a count across all nodes in a forest, so the count will be updated for daughter nodes

      cstat(:,ncur) = newstat


      ! retrieve the variable on which data is split
      ! retrieve the index of the variable on which the data is split for the current node k
      m = nint(nMatrix(k,4))

      ! retrieve the covariate
      ! we just retrieved the index, and now we retrieve the values of the covariate itself that was used in the split
      xm = xAll(:,m)

      ! if there are fewer than 1 categories for the variable (then it's numeric),

      IF (nCat(m) .LE. 1) THEN
        ! if a numeric variable, use the cutoff value to
        ! identify if individual i goes left

        ! use a cutoff value to determine if each individual goes left...
        ! xm holds values of a covariate that was used in a split for all individuals
        ! we see if it's less than the corresponding cutoff values for the split calculated earlier and get logical
        ! tst is logical where TRUE = left node, FALSE = right node

        tst = xm .LE. nMatrix(k,6)

      ELSE

      ! otherwise, if the variable is categorical,

        ! if an unordered factor, use category to identify if individual
        ! i goes left

        ! we initialize the value of tst to FALSE

        tst = .FALSE.

        ! loop through all of the category cutoffs of that covariate (held in nMatrix(k, 5))

        DO j = 1, nint(nMatrix(k,5))

        ! if the integer representation for the covariate value nint(xm) = the specific category value from that split, then they should go to the left branch so their tst values will be TRUE
        ! the cutoff is nMatrix(k, category level + 5)
        ! adding 5 accounts for the fact that cutoff values are sorted starting from the 6th column in nMatrix

          tst = tst .OR. nint(xm) .EQ. nint(nMatrix(k,j+5))

          ! end this loop from j = 1 to nint
        END DO

        ! end if statement checking if variable is numeric or categorical
      END IF

      ! now, loop through all ofcases

      DO j = 1, nAll

      ! if the individual's status doesn't match the value of the current node "k", then CYCLE to the next case

        IF (allStatus(j) .NE. k) CYCLE

        ! if the j-th case falls into the left node (tst(j) == TRUE)

        IF (tst(j)) THEN

        ! update the status of case "j" to the indices of the left daughter node

          allStatus(j) = nint(nMatrix(k,2))
        ELSE

        ! otherwise, if the j-th case falls into the right node,
        ! update the status of case "j" to the indices of the right daughter node
          allStatus(j) = nint(nMatrix(k,3))

        ! end IF statement about whether "tst" is TRUE or FALSE
        END IF

        ! end the loop over j to nAll
      END DO

    ! end the loop from k = 1, nrNodes

    END DO

    ! if there are any values where the first column of nMatrix (node status) == -2, this means further splitting can happen
    ! where this is true, since these are not interior nodes, we just turn them into terminal nodes to ensure no further splitting is done

    ! ensure that all nodes that are not "interior" are "terminal"
    WHERE (nMatrix(:,1) .EQ. -2) nMatrix(:,1) = -1


    !! ------- Now, we want to loop through each of the terminal nodes to create duplicated observations ------ !!

    ! 1. loop through each of the terminal nodes for processing
    ! 2. Retrieve observations for the terminal node
    ! 3. Retrieve propensity scores for the observations
    ! 4. Create duplicated observations
    ! 5. call calcvaluesingle on all the observations within the node & output this into ncurr

    ! ---------------------------------------------------------------------------------------------------------- !!

    ! store survival function values for each node of current tree

    trees(iTree)%survFunc = survFunc(:,1:ncur)

    !!!!!!!
    !!!!!!! Adding print statement for each node's survival function
    !!!!!!!

    ! Print values of variables
    !print *, "survival function values for each node of current tree:"
    !print *, "survival func for each node =", survFunc(:,1:ncur)

    !!!!!!!
    !!!!!!!
    !!!!!!!

    ! store mean survival time for each node of the current tree
    trees(iTree)%mean = mean(1:ncur)

    ! store survival probability for each node of the current tree
    trees(iTree)%survProb = survProb(1:ncur)

    ! stores splie information for each node of current tree
    trees(iTree)%matrix = nMatrix(1:ncur,:)

    ! stores total number of nodes in current tree
    trees(iTree)%nNode = ncur

    ! these lines compute aggregated metrics across all trees as they are looped through

    ! accumulates survival function values across all trees
    ! allStatus: contains indixes of terminal nodes to which each case in data belongs to
    ! therefore, for each case, the survival function from that terminal node is added as we iterature through each tree

    tforestSurvFunc = tforestSurvFunc + survFunc(:,allStatus)

    ! accumulates mean survival time across all trees

    forestMean = forestMean + mean(allStatus)

    ! accumulates survival probability across all trees
    forestSurvProb = forestSurvProb + survProb(allStatus)

    ! end the loop over each tree

  END DO

  ! This change is to eliminate a strange lto warning from R
!  forestSurvFunc = reshape(tforestSurvFunc, (/nt*nAll/)) / nTree
  j = 0
  DO i = 1, SIZE(tforestSurvFunc,2)
    forestSurvFunc(j+1:j+SIZE(tforestSurvFunc,1)) = tforestSurvFunc(:,i)
    j = j + SIZE(tforestSurvFunc,1)

    ! end this loop from 1 to SIZE
  END DO

  ! now, get the output arrays, which are the average values across all trees

  forestSurvFunc = forestSurvFunc / nTree
  forestMean = forestMean / nTree
  forestSurvProb = forestSurvProb / nTree

END SUBROUTINE tSurvTree


  ! **************************************************************** !




END MODULE IH_INNERS


  ! ******************************* subroutine: predictSurvTree ********************************* !

! n, integer, number of cases in data
! np, integer, number of covariates in data
! xt, integer(:), covariates
! nCat, integer(:), number of levels in each covariate (0 = continuous,
!   1 = ordered factors, 2+ = unordered factor)
! nt, integer, number of time points
! nNodes, integer, number of nodes
! tsurvFunc, real(:), survival functions for each node
! mean, real(:), mean survival of each node
! survProb, real(:), survival probability of each node
! nCols, integer, number of columns in node matrix
! tnodes, real(:), node information
! predSurvFunc, real(:), predicted survival functions
! predMean, real(:), predicted mean survival
! predSurvProb, real(:), predicted survival probabilities



SUBROUTINE predictSurvTree(n, np, xt, nCat, nt, nNodes, tsurvFunc, mean, &
  & survProb, nCols, tnodes, predSurvFunc, predMean, predSurvProb)

  ! all input parameters must be defined
  IMPLICIT NONE


  ! initialize a constant, "dp", returning a number with 15 decimal places of precision
  ! and a range representing up to 10^307
  INTEGER, PARAMETER :: dp = selected_real_kind(15,307)

  ! take input parameters: n, np, xt, nCat, nt, nNodes, tsurvFunc, mean, survProb, nCols, tnodes


  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: np
  REAL(dp), DIMENSION(1:n*np), INTENT(IN) :: xt
  INTEGER, DIMENSION(1:np), INTENT(IN) :: nCat
  INTEGER, INTENT(IN) :: nt
  INTEGER, INTENT(IN) :: nNodes
  REAL(dp), DIMENSION(1:nt*nNodes), INTENT(IN) :: tsurvFunc
  REAL(dp), DIMENSION(1:nNodes), INTENT(IN) :: mean
  REAL(dp), DIMENSION(1:nNodes), INTENT(IN) :: survProb
  INTEGER, INTENT(IN) :: nCols
  REAL(dp), DIMENSION(1:nNodes*nCols), INTENT(IN) :: tnodes

  ! output parameters: predSurvFunc, predMean, predSurvProb

  REAL(dp), DIMENSION(1:n*nt), INTENT(OUT) :: predSurvFunc
  REAL(dp), DIMENSION(1:n), INTENT(OUT) :: predMean
  REAL(dp), DIMENSION(1:n), INTENT(OUT) :: predSurvProb

  ! declaring local variables to be used within the subroutine

  ! index counters
  INTEGER :: i, j, m
  ! array to store terminal node to which each case belongs
  INTEGER, DIMENSION(1:n) :: stat

  ! array to store covariate values
  REAL(dp), DIMENSION(1:n) :: xm
  ! array to store survival functions
  REAL(dp), DIMENSION(1:nt, 1:n) :: tsurv

  ! reshaped array to store covariates
  REAL(dp), DIMENSION(1:n, 1:np) :: x
  ! reshaped array to store node info
  REAL(dp), DIMENSION(1:nNodes, 1:nCols) :: nodes
  ! reshaped array to store survival functions
  REAL(dp), DIMENSION(1:nt, 1:nNodes) :: survFunc

  LOGICAL, DIMENSION(1:n) :: sti, tst

  ! reshape the input arrays for covariates ("xt = x")
  ! node information (tnodes = nodes)
  ! survival functions for each node (tsurvFunc == survFunc): pt ends up in terminal node,
  !!!!    then fit a KM with all the pts from that terminal node to estimate node survival function

  x = reshape(xt, (/n, np/))
  nodes = reshape(tnodes, (/nNodes, nCols/))
  survFunc = reshape(tsurvFunc, (/nt, nNodes/))

  ! initialize array to store survival functions with 0

  tsurv = 0.d0

  ! column 1 is indicator of interior/terminal
  ! column 2 is the index of left
  ! column 3 is the index of right
  ! column 4 is the split variable
  ! column 5 is the number of cutoffs
  ! column 6:nCol are the cutoff values

  ! begin with every case in node 1
  ! stat will contain the terminal node to which each case belongs

  ! initialize to the root node (node 1) where every case starts
  stat = 1

  ! loop through each node in the treee

  DO i = 1, nNodes
    ! if terminal cycle

    ! if the node is terminal (nodes(i, 1) == -1.) cycle to next node iteration

    IF (nint(nodes(i,1)) .EQ. -1) CYCLE

    ! identify individuals in this node
    ! identify individuals in this node... stat == 1 initializes all individuals to be in the first node
    ! sti is true only if the individual is currently in the node being processed

    sti = stat .EQ. i

    ! retrieve the variable index on which data is split
    m = nint(nodes(i,4))

    ! retrieve the covariate using the index to subset the list of covariates
    xm = x(:,m)

    ! if the number of categories < 1,the variable is numeric, then

    IF (nCat(m) .LE. 1) THEN
      ! if a numeric variable, use the cutoff value to
      ! identify if individual i goes left

        ! xm holds values of a covariate that was used in a split for all individuals
        ! we see if it's <= the corresponding cutoff values for the split cutoff and get logical
        ! tst is logical where TRUE = left node, FALSE = right node


      tst = xm .LE. nodes(i,6)
    ELSE

    ! otherwise, if a variable is a factor,

      ! if an unordered factor, use category to identify if individual
      ! i goes left

      ! initialize all patients to right node

      tst = .FALSE.

      ! loop through number of cutoff values for the categorical variables at the current node
      DO j = 1, nint(nodes(i,5))

      ! if the value of the covariate == the j-th category value for categorical var, tst = TRUE
      ! meaning, move left
      ! checks if patient's categorical covariate value matches the specified categories for going left

        tst = tst .OR. nint(xm) .EQ. nint(nodes(i,j+5))

        ! end do loop across all cutoff values
      END DO

      ! end if statement for numeric/categorical vars
    END IF

    ! update the "stat" variable (which terminal node) based on if patient goes left or right
    ! tst = TRUE if patient is in left node. sti = TRUE if the person is currently in the node being considered
    ! if both tst and sti are true, update "stat" array with the value from the left daughter node index
    ! this will update whether a patient is present in a given node

    WHERE (tst .AND. sti) stat = nint(nodes(i,2))


    ! where both tst == FALSE and sti is true (in right node & present),
    ! update "stat" array with value from right daughter node index


    WHERE ( (.NOT. tst) .AND. sti) stat = nint(nodes(i,3))

    ! end loop across each node

  END DO

  ! retrieve appropriate values based on terminal node
  ! survFunc is an array storing survival functions for each node
  ! extract survival function for each person's terminal node
  tsurv = survFunc(:,stat)


  ! This change is to eliminate a strange lto warning from R
!  predSurvFunc = reshape(tsurv,(/n*nt/))
  j = 0
  DO i = 1, SIZE(tsurv,2)
    predSurvFunc(j+1:j+SIZE(tsurv,1)) = tsurv(:,i)
    j = j + SIZE(tsurv,1)
  END DO

  ! based on the terminal node, retrieve the predicted mean and survival probability
  ! mean = mean survival of each node
  ! survProb, survival probability of each node
  ! we retrieve these values for each patient based on the terminal node they belong to

  predMean = mean(stat)
  predSurvProb = survProb(stat)

END SUBROUTINE

  ! ***************************** subroutine: setUpBasics *********************************** !


  ! ** not in module IH_INNERS

! set up basic information for the module
! t_nt, integer, the number of time points
! t_dt, real(:), the time differences between time points
! t_rs, real, the probability for a random split
! t_ERT, integer, the indicator of extremely randomized trees
! t_uniformSplit, integer, the indicator of method for determining cut-off
!   when using ERT
! t_nodeSize, integer, the minimum number of cases in each node
! t_minEvent, integer, the minimum number of events in each node
! t_rule, integer, 0 = mean, 1 = logrank
! t_sIndex, integer, the indices of time points that is closest to the
!   requested survival time
! t_sFraction, real, the fractional distance between time points the the
!   requested survival time
! t_stratifiedSplit, real, the coefficient for determining stratification
! t_replace, integer, indicator of sampling with replacement


SUBROUTINE setUpBasics(t_nt, t_dt, t_rs, t_ERT, t_uniformSplit, t_nodeSize, &
                     & t_minEvent, t_rule, t_sIndex, t_sFraction, &
                     & t_stratifiedSplit, t_replace)

  ! using the inners module
  ! this allows you to have access to the "INNERS" subroutines, variables,types, etc

  USE IH_INNERS

  ! input variables must be explicitly defined

  IMPLICIT NONE

  ! declare input and output parameters

  INTEGER, INTENT(IN) :: t_nt
  REAL(dp), DIMENSION(1:t_nt), INTENT(IN) :: t_dt
  REAL(dp), INTENT(IN) :: t_rs
  INTEGER, INTENT(IN) :: t_ERT
  INTEGER, INTENT(IN) :: t_uniformSplit
  INTEGER, INTENT(IN) :: t_nodeSize
  INTEGER, INTENT(IN) :: t_minEvent
  INTEGER, INTENT(IN) :: t_rule
  INTEGER, INTENT(IN) :: t_sIndex
  REAL(dp), INTENT(IN) :: t_sFraction
  REAL(dp), INTENT(IN) :: t_stratifiedSplit
  INTEGER, INTENT(IN) :: t_replace

  ! assign value of input nt = t_nt (the number of time points)
  nt = t_nt

  ! assign value of input sIndex = t_sIndex (indices of time points closest to requested survival time)
  sIndex = t_sIndex

  ! assign value of input sFraction = t_sFraction (fractional distance between time pts & requested survival time)
  sFraction = t_sFraction

  ! sets "isSurvival" to TRUE it sIndex > 0 (AKA there are indices for itme points closest to survival time)
  ! this would be false if there are no time points, invalid survival time?
  isSurvival = sIndex > 0

  ! deallocates the "dt" array if it has been previously allocated to free up memory before reallocation
  ! this will hold the time differences between the time points

  IF (dtAllocated) DEALLOCATE(dt)

  ! allocate memory for the "dt" array with size "nt" (number of time points)

  ALLOCATE(dt(1:nt))

  ! update a boolean flag to indicate that "dt" array has been allocated

  dtAllocated = .TRUE.


  ! assign value of input dt = t_dt (time differences between time points)
  dt = t_dt

  ! assign other values ot input parameters to corresponding variables in "INNERS" module

  rs = t_rs
  ERT = t_ERT
  uniformSplit = t_uniformSplit
  nodeSize = t_nodeSize
  minEvent = t_minEvent
  rule = t_rule
  stratifiedSplit = t_stratifiedSplit
  replace = t_replace

END SUBROUTINE setUpBasics



  ! *********************************** subroutine: setUpInners ***************************** !

  ! ** not in module IH_INNERS

! set up basic information for the module that is step dependent
! t_n, integer, the number of cases under consideration
! t_np, integer, the number of covariates
! t_x, real(:), the covariates
! t_pr, real(:), the probability mass vector of survival function
! t_propensity real(:) the propensity scores
! t_delta, integer(:), the indicator of censoring
! t_mTry, integer, the maximum number of covariates to try for splitting
! t_nCat, integer(:), the number of categories in each covariate
! t_sampleSize, integer, the number of cases to sample for each tree
! t_ntree, integer, the number of trees in the forest
! t_nrNodes, integer, the maximum number of nodes


SUBROUTINE setUpInners(t_n, t_np, t_x, t_pr, t_delta, t_propensity, t_mTry, t_nCat,  t_sampleSize, t_nTree, t_nrNodes)

  ! use the "IH_INNERS" module to allow access to subroutines defined there

  USE IH_INNERS

  ! all input variables must be explicitly defined

  IMPLICIT NONE

  ! initialize input parameters

  INTEGER, INTENT(IN) :: t_n
  INTEGER, INTENT(IN) :: t_np
  REAL(dp), DIMENSION(1:t_n*t_np), INTENT(IN) :: t_x
  REAL(dp), DIMENSION(1:nt*t_n), INTENT(IN) :: t_pr

  ! the propensity score should be the same dimension as delta: one at each time point
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_propensity
  INTEGER, DIMENSION(1:t_n), INTENT(IN) :: t_delta
  INTEGER, INTENT(IN) :: t_mTry
  INTEGER, DIMENSION(1:t_np), INTENT(IN) :: t_nCat
  INTEGER, INTENT(IN) :: t_sampleSize
  INTEGER, INTENT(IN) :: t_nTree
  INTEGER, INTENT(IN) :: t_nrNodes

  ! assign input nAll = t_n (the number of cases under consideration)

  nAll = t_n

  ! assign input np = t_np (the number of covariates)
  np = t_np

  ! check if memory has already been allocated for these arrays (global)
  ! if so, then deallocate them them to free up memory

! note this should have propensityAll
  IF (isAllocated) THEN
    DEALLOCATE(xAll, prAll, deltaAll, propensityALL, nCat, forest%survFunc, forest%mean,  &
             & forest%survProb, trees)

             ! end if statement about allocation
  END IF

  ! next, allocate memory for a few of the arrays we just deallocated

  ! xAll will hold covariates
  ALLOCATE(xAll(1:nAll,1:np))

  ! prAll will hold probability mass vector of survival function
  ALLOCATE(prAll(1:nAll, 1:nt))

  ! propensityAll will hold propensity scores for all patients & stages
  ALLOCATE(propensityAll(1:nAll))

  ! deltaAll will hold indicator of censoring
  ALLOCATE(deltaAll(1:nAll))

  ! nCat will hold number of categories of each covariate
  ALLOCATE(nCat(1:np))

  ! update a boolean flag to indicate that arrays have been alocated for:
  ! xAll, prAll, deltaAll, nCat

  isAllocated = .TRUE.

  ! now, assign vaulues of input into these arrays we just allocated memory for
  ! reshape them first
  ! covariates from input t_x
  xAll = reshape(t_x, (/nAll,np/))

  ! probability mass vector from input t_pr
  prAll = reshape(t_pr, (/nAll,nt/))


    ! propensity score from input t_propensity
  propensityAll = t_propensity


  ! censoring indicator from t_delta
  deltaAll = t_delta

  ! number of categories of each covar from input t_nCat
  nCat = t_nCat

  ! maximum number of levels from input (nCat) for each covariate
  nLevs = max(maxval(nCat),1)

  ! assign input maximum number of covariates to try for splitting (t_mTry)
  mTry = t_mTry

  ! assign input number of cases to sample for each tree (t_sampleSize)
  sampleSize = t_sampleSize

  ! then, allocate memory for the rest of the arrays we previously deallocated for
  ! "forest" itself is a structure of type "ForestValues", and we have arrays to hold components within this structure defined at the beginning of the module


  ALLOCATE(forest%survFunc(1:nt, 1:nAll))
  ALLOCATE(forest%mean(1:nAll))
  ALLOCATE(forest%survProb(1:nAll))

  ! initialize these allocated arrays for survival function, mean survival, and survival probability

  forest%survFunc = 0.d0
  forest%mean = 0.d0
  forest%survProb = 0.d0

  ! assign input number of trees in forest (t_nTree)

  nTree = t_nTree

  ! allocate memory for the trees array with dimensions nTree

  ALLOCATE(trees(1:nTree))

  ! assign input maximum number of nodes (t_nrNodes)

  nrNodes = t_nrNodes

END SUBROUTINE setUpInners

  ! ************************************* subroutine: survTree *************************** !

    ! ** not in module IH_INNERS

! access function for calculating forest

SUBROUTINE survTree(tSurvFunc, mean, survProb)

! allows access to subroutines defined in the "IH_INNERS" module
  USE IH_INNERS

  ! all input variables must be explicitly declared
  IMPLICIT NONE

  ! declare nodes to store OUTPUT values

  REAL(dp), DIMENSION(1:nrNodes*nt), INTENT(OUT) :: tSurvFunc
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nrNodes), INTENT(OUT) :: survProb

  ! call the "tsurvTree" subroutine, passing declared output parameters as arguments
  ! in the tsurvTree, these are the inputs, which are declared as outputs
  ! (basically the input names here will be the outputs)
  ! tSurvFunc = forestSurvFunc, real(:), survival function averaged over forest
  ! mean = forestMean, real, mean survival averaged over forest
  ! survProb = forestSurvProb, real, survival probability averaged over forest

  CALL tsurvTree(tSurvFunc, mean, survProb)

END SUBROUTINE survTree

  ! *************************************** suroutine: treeDim ************************* !

! retrieve dimensions of node matrix for the iTree-th tree
!! ** note: "trees" is of type "NODE" defined in the global statements

SUBROUTINE treeDim(iTree, nr, nc)

  ! calls on the "INNERS" module to access the other subroutines in there
  USE IH_INNERS

  ! all variables must be explicitly declared
  IMPLICIT NONE

  ! input: index of the tree that the dimensions will be retrieved for

  INTEGER, INTENT(IN) :: iTree

  ! outputs: the number of rows (nr), number of columns (nc)
  INTEGER, INTENT(OUT) :: nr
  INTEGER, INTENT(OUT) :: nc

  !retrieve the number of rows from the "trees" object, access "matrix component"
  ! indexes to the ith- tree
  ! retrieve number of rows (number of nodes) since each node is represented by a row in the node matrix

  nr = size(trees(iTree)%matrix,1)

  ! retrieve number of columns in node matrix
  nc = size(trees(iTree)%matrix,2)

END SUBROUTINE treeDim

 ! ************************************* subroutine: getTree *************************** !

! retrieve the nodes, survival function, mean survival, and survival probability
! for the iTree-th tree

SUBROUTINE getTree(iTree, nr, nc, nodes, survFunc, mean, survProb)

  ! calls on the "INNERS" module to access the other subroutines in there
  USE IH_INNERS

    ! all variables must be explicitly declared
  IMPLICIT NONE

  ! inputs: index of the tree to be retrieved, number of rows, number of columns

  INTEGER, INTENT(IN) :: iTree
  INTEGER, INTENT(IN) :: nr
  INTEGER, INTENT(IN) :: nc

  ! output: the nodes, survival function, mean survival, survival prob

  REAL(dp), DIMENSION(1:nr*nc), INTENT(OUT) :: nodes
  REAL(dp), DIMENSION(1:nt*nr), INTENT(OUT) :: survFunc
  REAL(dp), DIMENSION(1:nr), INTENT(OUT) :: mean
  REAL(dp), DIMENSION(1:nr), INTENT(OUT) :: survProb

  ! create local loop counters

  INTEGER :: i, j

  ! This change is to eliminate a strange lto warning from R
!  nodes = reshape(trees(iTree)%matrix, (/nr*nc/))

  ! initiate the loop counter as 0
  j = 0

  ! loop through columns of the "matrix" part
  DO i = 1, SIZE(trees(iTree)%matrix,2)

    ! copies value from the i-th column of matrix into part of the "nodes" array
    ! the part of the nodes array starts at index j + 1, ends at j+SIZE(trees(iTree)%matrix,1)


    nodes(j+1:j+SIZE(trees(iTree)%matrix,1)) = trees(iTree)%matrix(:,i)

    ! then ncrement j by the nuber of elements copied in the current iteration
    j = j + SIZE(trees(iTree)%matrix,1)

    ! end loop through the columns of the tree matrix
  END DO

  ! This change is to eliminate a strange lto warning from R
!  survFunc = reshape(trees(iTree)%survFunc, (/nt*nr/))

! initiate loop counter as 0
  j = 0

  ! iterate over columns of "survFunc" column of the specified tree
  DO i = 1, SIZE(trees(iTree)%survFunc,2)

  ! gets value from ith column of survFunc into survFunc array
  ! start at j+1 then ends at j+SIZE(trees(iTree)%survFunc,1)

    survFunc(j+1:j+SIZE(trees(iTree)%survFunc,1)) = trees(iTree)%survFunc(:,i)

    ! increment j counter by number of elements copied in the current iteration
    j = j + SIZE(trees(iTree)%survFunc,1)

  END DO

  ! then, directly assign the values of mean and survProb from the specified trees
  mean = trees(iTree)%mean
  survProb = trees(iTree)%survProb

END SUBROUTINE getTree

  ! **************************************************************** !

