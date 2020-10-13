PROGRAM QE_MINIMAL

    USE kinds
    USE constants
    USE io_global
    USE ERR
    USE control_flags
    USE parameters
    USE matrix_inversion
    USE wannier_new
    USE symm_base
    USE cell_base
    USE ions_base
    USE random_numbers

END PROGRAM QE_MINIMAL


SUBROUTINE get_symm_base_s( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
                            cosbc_, trd_ht, rd_ht, cell_units_, s_out )
    !
    USE cell_base, ONLY: cell_base_init
    USE symm_base, ONLY: s, set_sym_bl
    USE kinds,     ONLY: DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ibrav_
    REAL(DP), INTENT(IN) :: celldm_ (6)
    LOGICAL, INTENT(IN) :: trd_ht
    REAL(DP), INTENT(IN) :: rd_ht (3,3)
    CHARACTER(LEN=*), INTENT(IN) :: cell_units_
    REAL(DP), INTENT(IN) :: a_ , b_ , c_ , cosab_, cosac_, cosbc_

    INTEGER, INTENT(OUT) :: s_out(3,3,48)
    INTEGER :: i, j, k

    CALL cell_base_init( ibrav_, celldm_, a_, b_, c_, cosab_, cosac_, &
                         cosbc_, trd_ht, rd_ht, cell_units_ )
    CALL set_sym_bl()

    DO i = 1, 3
        DO j = 1, 3
            DO k = 1, 48
                s_out(i,j,k) = s(i,j,k)
            ENDDO
        ENDDO
    ENDDO
    
    RETURN
END