    !  Stagr_Console.f90
    !
    !  FUNCTIONS:
    !  Stagr_Console      - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Stagr_Console
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************

    program Stagr_Console
    USE RivStagr4Mod_jmn
    use RivStagrMod_bmi
    implicit none
    INTEGER(4) count, num, i, status, cptArg
    logical :: lookForBMI=.FALSE.
    logical :: fileExist
    CHARACTER(LEN=250) buf
    CHARACTER(LEN=250) file
    
    !Check if arguments are found
    count = COMMAND_ARGUMENT_COUNT()

    if(count > 0) then
        !      CALL GETARG(1, buf, status)
        do cptArg=1,count
        CALL get_command_argument(cptArg, buf) !gnu fortran only take 2 args
        select case(buf)
        case('--BMI')
            lookForBMI = .TRUE.
            write(*,*)'BMI = ', lookForBMI
        case default
            inquire(file=adjustl(buf), exist=fileExist)
            if(.not.fileExist)then
                write(*,*)'file ', buf, ' not found'
                stop
            else
                file = adjustl(buf)
            endif
            !call stagr4(buf)
            write(*,*)'file name: ', file
        end select
        enddo
        if(lookForBMI) then
            inquire(file=file, exist=fileExist)
            if(.not.fileExist)then
                write(*,*)'file ', buf, ' not found'
                stop
            endif
            call STAGRBMI(file)
        else
            Call stagr4(file)
        endif
    endif

    ! Variables


    ! Body of Stagr_Console

    end program Stagr_Console

