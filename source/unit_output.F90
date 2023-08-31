!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module to:
! - set output unit
! - print array of messages to OUTPUT   
! - abort the code by 1) priting an error message and 2) closing files
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author - i.scivetti Feb 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module unit_output

  Use, Intrinsic :: iso_fortran_env, Only: error_unit,&
                                           input_unit,&
                                           output_unit
  Use numprec,                       Only: wi,&
                                           wp

  Implicit None
  Private

  Integer(Kind=wi), Save :: ounit

  Public :: close_unit, &
            error_stop, &
            info,       &
            set_output_unit

Contains

  Subroutine set_output_unit(n)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to set the unit for the output file 
    !
    ! author    - i.scivetti June 2020 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: n
    ounit=n

  End Subroutine set_output_unit

  Subroutine info(message, n)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print messages to OUT_EQCM file
    !
    ! author    - i.scivetti June 2020 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: n
    Character(Len=*),  Intent(In   ) :: message(n)

    Integer(Kind=wi) :: i

    Do i = 1, n
      Write (ounit, '(a)') Trim(message(i))
    End Do

  End Subroutine info

  Subroutine error_stop(message)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to:
    ! - print an error message
    ! - close all opened files previously 
    ! - stop the execution 
    !
    ! author    - i.scivetti June 2020 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Optional, Intent(In   ) :: message

    If (message /= ' ') Then
      Write (ounit,*)      ' '
    End If

    Write (ounit, '(a)') Trim(message)
    ! Print error to screen
    Write (output_unit,*) '*******************************************************'
    Write (output_unit,*) '** ERROR !!! ERROR !!! ERROR !!! ERROR !!! ERROR !!! **'
    Write (output_unit,*) '*******************************************************'
    Write (output_unit,*) '**  Please see details in the OUTPUT. If the         **'
    Write (output_unit,*) '**  error message makes no sense, check that the     **'
    Write (output_unit,*) '**  input files are free of non-ASCII characters     **'
    Write (output_unit,*) '*******************************************************'
    ! close all I/O channels
    Call close_unit(ounit)
    Error Stop

  End Subroutine error_stop

  Subroutine close_unit(i)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Close all opened files
    !
    ! author    - Ivan Scivetti June 2020 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(InOut) :: i

    Integer(Kind=wi) :: ierr
    Logical :: has_name, is_open

    Inquire (i, opened=is_open, named=has_name, iostat=ierr)

    If (is_open .and. has_name .and. ierr == 0 .and. All(i /= [-1, ERROR_UNIT, INPUT_UNIT, OUTPUT_UNIT])) Then
      Close (i)
    End If

  End Subroutine close_unit

End Module unit_output
