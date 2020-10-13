MODULE ERR

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: errore, infomsg

CONTAINS

    SUBROUTINE errore(routine, msg, ierr)
        !
        USE io_global, ONLY : stdout
        IMPLICIT NONE
        CHARACTER(*),    INTENT(in) :: routine, msg
        INTEGER,         INTENT(in) :: ierr
        !
        WRITE( stdout, FMT = '(/,1X,78("*"))')
        WRITE( stdout, FMT = '(5X,"from ",A," : error #",I10)' ) routine, ierr
        WRITE( stdout, FMT = '(5X,A)' ) msg
        WRITE( stdout, FMT = '(1X,78("*"),/)' )
        !
        STOP
        !
        RETURN
        !
    END SUBROUTINE errore
    !
    !
    !----------------------------------------------------------------------
    !
    !
    SUBROUTINE infomsg( routine, message )
        !
        USE io_global, ONLY : stdout
        IMPLICIT NONE
        CHARACTER (LEN=*) :: routine, message
        !
        ! the name of the calling routine
        ! the output message
        WRITE( stdout, '(5X,"Message from routine ",A,":")' ) routine
        WRITE( stdout, '(5X,A)' ) message
        RETURN
        !
    END SUBROUTINE infomsg

END MODULE ERR