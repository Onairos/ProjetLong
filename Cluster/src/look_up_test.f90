PROGRAM look_up_test
 USE module_look_up_table

  IMPLICIT NONE 

  REAL :: res = 0.0
  TYPE(look_up_table)::lt
  

   lt=look_up_init(3)
 
  CALL look_up_push(lt,"one", 1.0)
  CALL look_up_push(lt,"two", 2.0)
  CALL look_up_push(lt,"three", 3.0)
 
  CALL look_up_print(lt)
  PRINT*, ""
  res=look_up_get(lt,"one")
  PRINT*, "one =>",res
  res=look_up_get(lt,"two")
  PRINT*, "one =>", res
 
  CALL look_up_set(lt,"one", 60.0)
  PRINT*, "now one must have changed to 60.2"
  PRINT*,""
  CALL look_up_print(lt)
  res=look_up_get(lt,"one")
  PRINT*,""
  PRINT*,"one =>", res
  PRINT*,""
  PRINT*,"Exist one ?", look_up_exists(lt, "one")
  PRINT*,""
  PRINT*,"Exist two ?", look_up_exists(lt, "two")
  PRINT*,""
  PRINT*,"Exist three ?", look_up_exists(lt, "three")
  PRINT*,""
  PRINT*,"Exist four ?", look_up_exists(lt, "four")
END PROGRAM look_up_test
