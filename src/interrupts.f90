MODULE interrupts
    IMPLICIT NONE
    LOGICAL,SAVE :: is_interrupted, is_inited = .FALSE.
    PRIVATE      :: is_interrupted, is_inited
  CONTAINS
    SUBROUTINE init_signal_handler()
        USE IFPORT
        INTEGER(4) :: catch_signal, sigret
        EXTERNAL   :: catch_signal
        IF (.NOT. is_inited) THEN
            is_inited = .TRUE.
            sigret = SIGNAL(SIGTERM, catch_signal, -1)
        END IF
    END SUBROUTINE
    SUBROUTINE set_is_interrupted(state)
        LOGICAL :: state
        is_interrupted = state
    END SUBROUTINE
    FUNCTION get_is_interrupted()
        LOGICAL :: get_is_interrupted
        get_is_interrupted = is_interrupted
    END FUNCTION
END MODULE interrupts

INTEGER(4) FUNCTION catch_signal(sig_num)
    USE interrupts
    INTEGER(4)         :: sig_num

    PRINT *, 'Detected SIGTERM signal.'
    CALL set_is_interrupted(.TRUE.)
    catch_signal = 1
END FUNCTION catch_signal
