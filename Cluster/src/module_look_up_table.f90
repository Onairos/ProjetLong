MODULE module_look_up_table
IMPLICIT NONE

 
  PUBLIC :: look_up_init
  PUBLIC :: look_up_push
  PUBLIC :: look_up_get
  PUBLIC :: look_up_set
  PUBLIC :: look_up_print

 TYPE look_up_table
  
  CHARACTER(128), ALLOCATABLE :: keys(:)
  DOUBLE PRECISION,                  ALLOCATABLE :: values(:)
  INTEGER                            :: look_up_index
 END TYPE look_up_table


CONTAINS
 
  TYPE(look_up_table) FUNCTION look_up_init(no_items)
    INTEGER, INTENT(in)  :: no_items
    INTEGER              :: status
    ALLOCATE(look_up_init%keys(no_items), stat=status)
    ALLOCATE(look_up_init%values(no_items), stat=status)
    look_up_init%look_up_index = 0
    RETURN
  END 
 
  SUBROUTINE look_up_push(lu_table, key, value)
    TYPE(look_up_table) :: lu_table
    CHARACTER(*), INTENT(IN)     :: key
    DOUBLE PRECISION        , INTENT(IN)     :: value
    lu_table%look_up_index = lu_table%look_up_index + 1
    IF(lu_table%look_up_index > Size(lu_table%keys, 1)) CALL print_error("Error: look_up table is already full")
    lu_table%keys(lu_table%look_up_index) = key
    lu_table%values(lu_table%look_up_index) = value
  END SUBROUTINE look_up_push
 
  SUBROUTINE look_up_set(lu_table,key, value)
    TYPE(look_up_table) :: lu_table
    CHARACTER(*), INTENT(IN)     :: key
    DOUBLE PRECISION        , INTENT(IN)     :: value
    INTEGER                      :: local_index
    LOGICAL                      :: found
    found = .FALSE.
    IF(lu_table%look_up_index.LT.Size(lu_table%keys, 1)) CALL print_error("Error: the look_uptable is not yet full")
    DO local_index = 1,Size(lu_table%keys,1)
      IF(TRIM(lu_table%keys(local_index)) == TRIM(key)) THEN
        lu_table%values(local_index) = value
        found = .TRUE.
      ENDIF
    ENDDO
    IF(.NOT.found) CALL print_error("Unknown key")
  END SUBROUTINE look_up_set
 
  DOUBLE PRECISION FUNCTION look_up_get(lu_table, key)
    TYPE(look_up_table) :: lu_table
    CHARACTER(*), INTENT(IN)     :: key
    INTEGER                      :: local_index
    LOGICAL                      :: found
    found = .FALSE.
    IF(lu_table%look_up_index.LT.Size(lu_table%keys, 1)) CALL print_error("Error: the look_uptable is not yet full")
    DO local_index = 1,Size(lu_table%keys,1)
      IF(TRIM(lu_table%keys(local_index)) == TRIM(key)) THEN
        look_up_get = lu_table%values(local_index)
        RETURN
      ENDIF
    ENDDO
  END FUNCTION look_up_get

  LOGICAL FUNCTION look_up_exists(lu_table, key)
    TYPE(look_up_table) :: lu_table
  CHARACTER(*), INTENT(IN)     :: key
  INTEGER                      :: local_index
  IF(lu_table%look_up_index.LT.Size(lu_table%keys, 1)) CALL print_error("Error: the look_uptable is not yet full")  
  look_up_exists = .FALSE.
  DO local_index = 1,Size(lu_table%keys,1)
    IF(TRIM(lu_table%keys(local_index)) == TRIM(key)) THEN
      look_up_exists = .TRUE.
    ENDIF
  ENDDO
  RETURN
  END FUNCTION look_up_exists

  SUBROUTINE look_up_print(lu_table)
    TYPE(look_up_table) :: lu_table
    INTEGER  :: local_index
    PRINT*, "Contents of the look_uptable:"
    DO local_index = 1,Size(lu_table%keys,1)
      PRINT*, TRIM(lu_table%keys(local_index)), " = ", lu_table%values(local_index)
    ENDDO
  END SUBROUTINE look_up_print
 
  SUBROUTINE print_error(text)
    CHARACTER(*) :: text
    PRINT*, text
    STOP
  END SUBROUTINE print_error
END MODULE module_look_up_table
