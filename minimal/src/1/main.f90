PROGRAM QE_MINIMAL

    USE error_handler
    USE kinds
    USE constants
    USE io_global
    USE control_flags
    USE parameters
    USE matrix_inversion
    USE wannier_new
    USE symm_base
    USE cell_base
    USE ions_base
    USE random_numbers

END PROGRAM QE_MINIMAL


SUBROUTINE get_symm_base_s( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, cosbc_, s_out )
    !
    USE cell_base
    USE symm_base
    USE kinds,     ONLY: DP
    !
    IMPLICIT NONE
    !
    EXTERNAL :: dgetrf, dgetri, zgetrf, zgetri, DGEMV, DGER
    !
    INTEGER, INTENT(IN) :: ibrav_
    !
    REAL(DP), INTENT(IN) :: celldm_ (6)
    !  traditional crystallographic cell parameters (alpha=cosbc and so on)
    REAL(DP) :: rd_ht (3,3)
    !
    !    CELL_PARAMETERS (cell_option)
    !      HT(1,1) HT(1,2) HT(1,3)
    !      HT(2,1) HT(2,2) HT(2,3)
    !      HT(3,1) HT(3,2) HT(3,3)
    !
    !      card_cell_parameters() in read_cards.f90
    !     
    !      DO i = 1, 3
    !           CALL read_line( input_line )
    !           READ(input_line,*) ( rd_ht( i, j ), j = 1, 3 )
    !      ENDDO
    !      !
    !      trd_ht = .true.
    !      tcell  = .true.
    !
    REAL(DP), INTENT(IN) :: a_ , b_ , c_ , cosab_, cosac_, cosbc_
    !
    INTEGER, INTENT(OUT) :: s_out(3,3,48)
    !INTEGER :: i, j, k

    rd_ht  = 0.0_DP

    CALL cell_base_init( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
                         cosbc_, .false., rd_ht, 'alat' )
    CALL set_sym_bl()

    s_out = s + 0.0

    RETURN
END
!
!----------------------------------------------------------------------
!
!SUBROUTINE errore(routine, msg, ierr)
!    !
!    USE io_global, ONLY : stdout
!    IMPLICIT NONE
!    CHARACTER(*),    INTENT(in) :: routine, msg
!    INTEGER,         INTENT(in) :: ierr
!    !
!    IF( ierr <= 0 ) THEN
!        RETURN
!    ELSE
!        WRITE( stdout, FMT = '(/,1X,78("*"))')
!        WRITE( stdout, FMT = '(5X,"from ",A," : error #",I10)' ) routine, ierr
!        WRITE( stdout, FMT = '(5X,A)' ) msg
!        WRITE( stdout, FMT = '(1X,78("*"),/)' )
!        !
!        STOP
!        !
!        RETURN
!    ENDIF
!    !
!END SUBROUTINE errore
!
!
!----------------------------------------------------------------------
!
!
!SUBROUTINE infomsg( routine, message )
!    !
!    USE io_global, ONLY : stdout
!    IMPLICIT NONE
!    CHARACTER (LEN=*) :: routine, message
!    !
!    ! the name of the calling routine
!    ! the output message
!    WRITE( stdout, '(5X,"Message from routine ",A,":")' ) routine
!    WRITE( stdout, '(5X,A)' ) message
!    RETURN
    !
!END SUBROUTINE infomsg