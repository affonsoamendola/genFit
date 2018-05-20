program genFit

    implicit none

    DOUBLE PRECISION, DIMENSION(2) :: r_firstGuess
    DOUBLE PRECISION, DIMENSION(3,2) :: r_Values 

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: r_Parameters
    !First dimension is the index of the life
    !Second dimension is the number of parameters + 1, last index is reserved for life fitness
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: r_Winners
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: r_Temp

    DOUBLE PRECISION :: r_mutChance
    !Chance of a mutation happening when getting gene from parent
    DOUBLE PRECISION :: r_mutConstant
    !Constant for mutation strength, might make this variable
    DOUBLE PRECISION :: r_totalMutChance
    !Chance of completely changing a gene, instead of taking one from parent, creates a random one
    DOUBLE PRECISION :: r_totalMutConstant
    !Constant of total mutation strength, basically, the limits for the new gene
    DOUBLE PRECISION :: r_parentBias
    !Chance of taking gene from parent 1 or 2, 0.5 is 50%, 0.2 is 20% for parent 1, 80"for parent 2

    DOUBLE PRECISION :: r_criteria
    !Criteria for stopping process

    INTEGER :: i_liveNumber
    !Number of generations

    INTEGER :: i_parameterNumber
    !Number of adjustable parameters

    INTEGER :: i_winGen
    !Number of Winner Generations, that get to reproduce

    INTEGER :: i_Time

    INTEGER :: i

    INTEGER :: i_generation
    INTEGER :: i_maxGenerations

    LOGICAL :: b_foundAnswer

    INTERFACE
        SUBROUTINE calculateFitness(r_Parameters, r_Values)
            DOUBLE PRECISION, DIMENSION(:,:) :: r_Values
            DOUBLE PRECISION, DIMENSION(:) :: r_Parameters
        END SUBROUTINE
       DOUBLE PRECISION function y(r_x,r_Parameters) 
            DOUBLE PRECISION :: r_x
            DOUBLE PRECISION, DIMENSION(:) :: r_Parameters
        END function
        subroutine generateNewParameters(r_Parent1, r_Parent2, r_Offspring, r_parentBias, r_mutChance, r_mutConstant, &
                                 r_totalMutChance, r_totalMutConstant)
            DOUBLE PRECISION, DIMENSION(:) :: r_Parent1
            DOUBLE PRECISION, DIMENSION(:) :: r_Parent2
            DOUBLE PRECISION, DIMENSION(:) :: r_Offspring

            DOUBLE PRECISION :: r_mutChance
            DOUBLE PRECISION :: r_mutConstant
            DOUBLE PRECISION :: r_totalMutChance
            DOUBLE PRECISION :: r_totalMutConstant
            DOUBLE PRECISION :: r_parentBias
        END SUBROUTINE
        SUBROUTINE SORT_DOUBLE_ARRAY(r_Array)
            DOUBLE PRECISION, INTENT(INOUT) :: r_Array(:,:)
        END SUBROUTINE
    END INTERFACE

    r_Values = reshape((/ 0., 1., 2., 0., 1., 2. /),shape(r_Values))
    r_firstGuess = reshape((/1.1,0./), shape(r_firstGuess))

    b_foundAnswer = .false.
    r_mutChance = 1d0
    r_mutConstant = 1d-20
    r_totalMutChance = 0.1d0
    r_totalMutConstant = 1d-5

    r_parentBias = 0.5d0

    i_parameterNumber = SIZE(r_firstGuess)
    i_liveNumber = 50
    i_winGen = 15

    i_generation = 0
    i_maxGenerations = 1000000

    ALLOCATE(r_Parameters(i_liveNumber, i_parameterNumber + 1))
    !Last number of the r_Parameters list is fitness of that instance
    ALLOCATE(r_Winners(i_winGen, i_parameterNumber + 1))
    ALLOCATE(r_Temp(i_parameterNumber))

    CALL SYSTEM_CLOCK(i_Time)
    CALL SRAND(i_Time)

    r_Parameters(1,:SIZE(r_firstGuess)) = r_firstGuess(:)
    r_Parameters(1,SIZE(r_firstGuess)+1) = 50.
    r_Parameters(3:,:) = 0.

    CALL SORT_DOUBLE_ARRAY(r_Parameters)

    DO i=i_liveNumber-1, 1,-1
        CALL generateNewParameters(r_Parameters(i_liveNumber,:),r_Parameters(i_liveNumber,:),&
                                   r_Parameters(i,:), 1.0d0, 1.0d0, r_mutConstant, -1.0d0, 0.0d0)
        
    END DO
    
    DO i_generation=1, i_maxGenerations
        !PRINT *,"EVERYONE"
        DO i=1, i_liveNumber
            CALL calculateFitness(r_Parameters(i,:),r_Values)
            !PRINT *,r_Parameters(i,:)
        END DO
        !PRINT *,"EVERYONE---"
        CALL SORT_DOUBLE_ARRAY(r_Parameters)
        !PRINT *,"WINNERS"
        DO i=1, i_winGen
            r_Winners(i,:)=r_Parameters(i,:)
            
            !PRINT *,r_Winners(i,:)
            
        END DO
        !PRINT *,"WINNERS---"

        DO i=1, i_liveNumber

            CALL generateNewParameters(r_Winners(INT(RAND()*(DBLE(i_winGen)-0.001)+1.),:), &
                                       r_Winners(INT(RAND()*(DBLE(i_winGen)-0.001)+1.),:), &
                                       r_Parameters(i,:),r_parentBias,r_mutChance,r_mutConstant,r_totalMutChance,r_totalMutConstant)

        END DO

        PRINT *, i_generation
    END DO

    PRINT *,"ULTIMATE WINNERS"
        DO i=1, i_winGen
            r_Winners(i,:)=r_Parameters(i,:)
            
            PRINT *,r_Winners(i,:)
            
        END DO
    PRINT *,"ULTIMATE WINNERS---"

    DEALLOCATE(r_Parameters)
    DEALLOCATE(r_Temp)
    DEALLOCATE(r_Winners)
end program

subroutine generateNewParameters(r_Parent1, r_Parent2, r_Offspring, r_parentBias, r_mutChance, r_mutConstant, &
                                 r_totalMutChance, r_totalMutConstant)
    implicit none
    DOUBLE PRECISION, DIMENSION(:) :: r_Parent1
    DOUBLE PRECISION, DIMENSION(:) :: r_Parent2
    DOUBLE PRECISION, DIMENSION(:) :: r_Offspring

    DOUBLE PRECISION :: r_mutChance
    DOUBLE PRECISION :: r_mutConstant
    DOUBLE PRECISION :: r_totalMutChance
    DOUBLE PRECISION :: r_totalMutConstant
    DOUBLE PRECISION :: r_parentBias

    DOUBLE PRECISION :: r_mutation

    INTEGER :: i

    DO i=1, size(r_Parent1)-1
        IF(RAND() <= r_totalMutChance) THEN
            r_Offspring(i) = r_Offspring(i) + (RAND()-0.5)*r_totalMutConstant    
        ELSE
            r_mutation = 0.
            IF(RAND() <= r_mutChance) THEN
                r_mutation = (RAND()-0.5)*r_mutConstant
            END IF
            IF(RAND() <= r_parentBias) THEN
                r_Offspring(i) = r_Parent1(i) + r_mutation 
                !PRINT *, "PARENT1"
            ELSE
                r_Offspring(i) = r_Parent2(i) + r_mutation
                !PRINT *, "PARENT2"
            END IF
        END IF
    END DO
END subroutine

SUBROUTINE calculateFitness(r_Parameters, r_Values)
    implicit none
    DOUBLE PRECISION, DIMENSION(:,:) :: r_Values
    DOUBLE PRECISION, DIMENSION(:) :: r_Parameters

    DOUBLE PRECISION :: r_chi2

    INTEGER :: i

    INTERFACE
        DOUBLE PRECISION function y(r_x, r_Parameters)
            DOUBLE PRECISION :: r_x
            DOUBLE PRECISION, DIMENSION(:) :: r_Parameters
        END function
    END INTERFACE 

    r_chi2 = 0.

    DO i=1, SIZE(r_Values,1)
        r_chi2 = r_chi2 + (r_Values(i,2)-y(r_Values(i,1),r_Parameters))**2 
    END DO
    r_Parameters(SIZE(r_Parameters)) = r_chi2
    return
end SUBROUTINE

DOUBLE PRECISION function y(r_x,r_Parameters)
    implicit none
    DOUBLE PRECISION :: r_x
    DOUBLE PRECISION, DIMENSION(:) :: r_Parameters

    y = r_Parameters(1)*r_x + r_Parameters(2)

    return
end function

SUBROUTINE SORT_DOUBLE_ARRAY(r_Array)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT) :: r_Array(:,:)

    DOUBLE PRECISION :: r_currentSmaller
    DOUBLE PRECISION, ALLOCATABLE :: r_tempHolder(:)
    INTEGER :: i_smallerIndex
    
    INTEGER :: i
    INTEGER :: j
    !Alocação de variaveis
    !---------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ALLOCATE(r_tempHolder(SIZE(r_Array, 2)))

    DO i=1, SIZE(r_Array,1)
        r_currentSmaller = r_Array(i,SIZE(r_Array,2))
        i_smallerIndex = i
        !Inicializa os menores valores encontrados como os primeiros valores, se eles não
        !forem de fato os menores, eles seram substituidos pelo menor valor no loop a seguir

        DO j = i, SIZE(r_Array,1)
        !Começa o loop a partir da posição de i, para não substituir valores ja ordenados
            !Considerando uma comparação como um passo, incrementa o valor de passos levados

            IF (r_currentSmaller > r_Array(j,SIZE(r_Array,2))) THEN

                !Considerando uma troca de valores como um passo, incrementa o valor de passos levados
                i_smallerIndex = j
                r_currentSmaller = r_Array(j,SIZE(r_Array,2))
            END IF
            !Se o valor sendo checado atualmente for menor que o suposto menor valor, troca o suposto 
            !menor valor pelo valor sendo checado atualmente, e o indice do suposto menor
            !valor pelo indice do checado atualmente

        END DO

        !Considerando uma troca de valores como um passo, incrementa o valor de passos levados 

        r_tempHolder(:) = r_Array(i,:)
        r_Array(i,:) = r_Array(i_smallerIndex,:)
        r_Array(i_smallerIndex,:) = r_tempHolder(:)
        !Troca o valor do indice i atual pelo menor valor encontrado no array inteiro.

    END DO
    !Para cada indice, a subrotina checa todos os valores a partir dele por um valor menor que ele, se encontrar
    !ela troca o valor contido no indice pelo menor valor encontrado
    !basicamente uma implementação de Selection Sort

    DEALLOCATE(r_tempHolder)
END SUBROUTINE