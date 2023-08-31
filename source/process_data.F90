!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutines to perform various different operations with strings
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! author       - i.scivetti  Feb 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module process_data 

  Use numprec,     Only : wi,&
                          wp

  Use unit_output, Only : info, &
                          error_stop

  Implicit None
  Private

  Public :: capital_to_lower_case, &
            check_for_symbols,     &  
            check_for_rubbish,     &
             detect_rubbish,       & 
            get_word_length,       &
            remove_symbols,        &
            remove_front_tabs
Contains

  Subroutine get_word_length(word,length)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Obtain the numer of characters in a string
    ! 
    ! author  - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: word
    Integer(Kind=wi), Intent(  Out) :: length

    Logical                         :: flag 

    length = 0
    flag = .true.

    ! Start transferring
    Do While (flag)
      If (word(length+1:length+1) == ' ') Then
        flag = .false.
      Else
        length=length+1 
      End If
    End Do

  End Subroutine get_word_length


  Subroutine capital_to_lower_case(string)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Transform string from capital to lower case 
    ! 
    ! author  - i.scivetti June 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(InOut) :: string

    Integer(Kind=wi) :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    ! Capitalize each letter if it is lowecase
    Do i = 1, Len_Trim(string)
      ic = Index(cap, string(i:i))
      If (ic > 0) Then
        string(i:i) = low(ic:ic)
      End If
    End Do

  End Subroutine capital_to_lower_case

  Subroutine remove_symbols(string, list)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Eliminates from string any of the symbols of array list
    ! 
    ! author  - i.scivetti June 2020
    ! refact  - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(InOut) :: string
    Character(Len=*), Intent(In   ) :: list

    Integer(Kind=wi) :: ic, i, iadd

    iadd=0

    ! Find if there is any of the symbols defined in "list"
    Do i = 1, Len_Trim(string)
        ic = Index(list, string(i:i))
        If (ic > 0) Then 
          iadd=iadd+1
          string(i:i) = ' '
        End If 
    End Do

    string=Adjustl(Trim(string))

  End Subroutine remove_symbols

  Subroutine check_for_symbols(string, list, error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finds if any of the characters of array list is present in string
    !
    ! author  - i.scivetti June 2020
    ! refact  - i.scivetti July 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(InOut) :: string
    Character(Len=*), Intent(In   ) :: list   
    Logical,          Intent(InOut) :: error                                 

    Integer(Kind=wi) :: ic, i

    i=1 
    ! Find if there is a is any of the symbols defined in "list"
    Do While (i <= Len_Trim(string) .And. (.Not. error)) 
      ic = Index(list, string(i:i))
      If (ic > 0) Then
        error=.True.
      End If
      i=i+1
    End Do

  End Subroutine check_for_symbols


  Subroutine remove_front_tabs(string)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finds if there are tabs in front of the string and remove them
    !
    ! author  - i.scivetti  Oct 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(InOut) :: string
    
    Integer(Kind=wi) :: i
    Logical :: change

    change=.False.
    i=1

    Do While (i <= Len_Trim(string) .And. (.Not. change))
      If(string(i:i) == achar(9)) Then
        string(i:i)=' '
      Else
        change=.True.
      End If
      i=i+1
    End Do

  End Subroutine remove_front_tabs

  Subroutine check_for_rubbish(iunit, error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finds if there are "wrong" characters in the definition of sentences
    !
    ! author  - i.scivetti December 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer,          Intent(In   ) :: iunit
    Character(Len=*), Intent(In   ) :: error

    Character(Len=256) :: string

    Backspace iunit
    Read (iunit, Fmt='(a)') string
    !string=Trim(Adjustl(string))

    Call detect_rubbish(string, error)
    Backspace iunit

  End Subroutine check_for_rubbish

  Subroutine detect_rubbish(string, error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Detect wrong characters 
    !
    ! author  - i.scivetti December 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: string
    Character(Len=*), Intent(In   ) :: error

    Character(Len=5)   :: list
    Character(Len=256) :: messages(7)
    Integer(Kind=wi) :: ic_hash, ic_rubbish, i
    Logical :: hash_found, rubbish_found, trigger, fread

    list='/\,;'//"'"
    hash_found=.False.
    rubbish_found=.False.
    trigger=.False.
    fread=.True.
            
    ! Find if # is defined in "string"
    i=1
    Do While (i <= Len_Trim(string) .And. (.Not. hash_found))
      ic_hash = Index('#', string(i:i))
      If (ic_hash > 0) Then
        ic_hash=i
        hash_found=.True.
      End If
      i=i+1
    End Do

    ! Find if any of the characters in "list" is found in "string"
    i=1
    Do While (i <= Len_Trim(string) .And. (.Not. rubbish_found))
      ic_rubbish = Index(list, string(i:i))
      If (ic_rubbish > 0) Then
        ic_rubbish=i 
        rubbish_found=.True.
      End If
      i=i+1
    End Do

    If (rubbish_found) Then
      trigger=.True.
      If (hash_found) Then
        If (ic_rubbish>ic_hash) Then
          trigger=.False.
        End If
      End If
    End If

    If (trigger) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,2a)') '*** ERROR in ', Trim(error) 
      Write (messages(2),'(1x,3a)') '*** At least one of these characters (', Trim(list), ') is found in the following line:'     
      Write (messages(3),'(a)')     ' '
      Write (messages(4),'(1x,a)')   Trim(string)
      Write (messages(5),'(a)')     ' '
      Write (messages(6),'(1x,a)')  'Please remove the wrong character(s) and rerun. If defined, real numbers&
                                    & MUST use dots (no commas).'
      Write (messages(7),'(1x,a)')  'IMPORTANT: There is no need to remove "wrong" characters when they are&
                                    & part of a comment. Comments MUST start with "#".'
      Call info(messages, 7) 
      Call error_stop(' ')
    End If

  End Subroutine detect_rubbish

End module process_data



