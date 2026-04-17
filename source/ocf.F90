!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the Orientational Correlation Function (OCF)
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module ocf

  Use atomic_model, Only: model_type, &
                          check_PBC
                          
  Use input_types,  Only: in_integer, &
                          in_logic,   &
                          in_string
                    
  Use numprec,      Only: wi,&
                          wp 
                    
  Use fileset,      Only: file_type,  &
                          refresh_out, &
                          FILE_SET, & 
                          FILE_OCF_ALL, &
                          FILE_OCF_AVG, &
                          FILE_CHEM_OCF_ALL, &                          
                          FILE_CHEM_OCF_AVG
                           
  Use process_data, Only: set_read_status, &
                          get_word_length, &
                          capital_to_lower_case, &
                          check_for_rubbish

  Use trajectory,   Only: traj_type,&
                          within_region, &
                          average_segments

  Use unit_output,  Only: error_stop, &
                          info

  Implicit None
  Private
  
  Type, Public :: chemocf_type
    Private
    Type(in_string), Public  :: invoke
    Type(in_string), Public  :: variable
    Type(in_logic),  Public  :: print_all_intervals
  End Type
    
  Type, Public :: ocf_type
    Private
    Type(in_string),  Public :: invoke
    Type(in_integer), Public :: legendre_order
    Type(in_string),  Public :: u_definition
    Type(in_logic),   Public :: print_all_intervals
  End Type  ocf_type

  Public :: read_ocf, check_ocf, orientational_correlation_function
  Public :: read_orientational_chemistry, check_orientational_chemistry, compute_orientational_chemistry
  
Contains 

  Subroutine read_ocf(iunit, ocf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information for orientational correlation function (OCF)
    ! analysis from the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(ocf_type),    Intent(InOut) :: ocf_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &OCF block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_ocf" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_ocf') Exit
      Call check_for_rubbish(iunit, '&ocf')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='u_definition') Then
        Read (iunit, Fmt=*, iostat=io) word, ocf_data%u_definition%type
        Call set_read_status(word, io, ocf_data%u_definition%fread, ocf_data%u_definition%fail,&
                           & ocf_data%u_definition%type)

      Else If (Trim(word)=='legendre_order') Then
         Read (iunit, Fmt=*, iostat=io) word, ocf_data%legendre_order%value
         Call set_read_status(word, io, ocf_data%legendre_order%fread,&
                            & ocf_data%legendre_order%fail)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, ocf_data%print_all_intervals%stat
       Call set_read_status(word, io, ocf_data%print_all_intervals%fread, ocf_data%print_all_intervals%fail)
                            
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_ocf"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_ocf


  Subroutine check_ocf(files, ocf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(ocf_type),     Intent(InOut) :: ocf_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the &OCF block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (ocf_data%legendre_order%fread) Then
      If (ocf_data%legendre_order%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "legendre_order" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (ocf_data%legendre_order%value < 1 .Or. ocf_data%legendre_order%value>4) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                &'Input value for "legendre_order" must be a value between 1 and 4 (polynomial order).'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "legendre_order" directive'
       Call info(messages, 1)
       Call error_stop(' ')
    End If

    If (ocf_data%u_definition%fread) Then
      If (ocf_data%u_definition%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "u_definition" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(ocf_data%u_definition%type)/='bond_12'  .And. &
            Trim(ocf_data%u_definition%type)/='bond_13'  .And. &
            Trim(ocf_data%u_definition%type)/='bond_123' .And. &
            Trim(ocf_data%u_definition%type)/='bond_12-13'  .And. &
            Trim(ocf_data%u_definition%type)/='plane') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "u_definition". Valid options: "bond_12", "bond_13",&
                                    & "bond_12-13", "bond_123" or "plane"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "u_definition" directive'
       Call info(messages, 1)
       Call error_stop(' ')
    End If
    
    If (ocf_data%print_all_intervals%fread) Then
      If (ocf_data%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      ocf_data%print_all_intervals%stat=.False.
    End If
    
  End Subroutine check_ocf


  Subroutine orientational_correlation_function(files, traj_data, ocf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the orientational correlation function (OCF)
    ! Different possible flavours are available depending on the settings
    ! of the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),    Intent(InOut) :: ocf_data

    Integer(Kind=wi)   :: i, j, k, l
    Integer(Kind=wi)   :: Nini_species, iunit
    Real(Kind=wp)      :: suma_i 
    Real(Kind=wp)      :: time
    Real(Kind=wp)      :: base_time
    Logical            :: set_u0
    Character(Len=256) :: message
    
    Logical           :: terminated(traj_data%Nmax_species)

    If (ocf_data%print_all_intervals%stat) Then
      ! Print header
      Open(Newunit=files(FILE_OCF_ALL)%unit_no, File=files(FILE_OCF_ALL)%filename, Status='Replace')
      iunit=files(FILE_OCF_ALL)%unit_no
      Write (iunit,'(a)') '#  OCF Analysis for all the intervals' 
      Write (iunit,'(a)') '#  Time (ps)         OCF' 
    End If
    
    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1

    Do k= 1, traj_data%analysis%N_seg
      set_u0=.True.
      l=0
      ! Initialise terminated tag
      Do j = 1, traj_data%Nmax_species
        terminated(j)=.False.
      End Do
      base_time=(traj_data%analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        l=l+1
        time=(i-1)*traj_data%timestep%value
        If (.Not. set_u0) Then
          Do j=1,3
            traj_data%species(i,:)%u0(j,1)=traj_data%species(i-1,:)%u0(j,1)
            If (Trim(ocf_data%u_definition%type) == 'bond_12-13') Then
              traj_data%species(i,:)%u0(j,2)=traj_data%species(i-1,:)%u0(j,2)
            End If
          End Do
        Else
          Nini_species=0
          Do j = 1, traj_data%Nmax_species
            If (traj_data%species(i,j)%alive) Then
              Call rotation_vector_monitored_species(traj_data, ocf_data, i, j)
              traj_data%species(i,j)%u0(:,1)=traj_data%species(i,j)%u(:,1)
              If (Trim(ocf_data%u_definition%type) == 'bond_12-13') Then
                traj_data%species(i,j)%u0(:,2)=traj_data%species(i,j)%u(:,2)
              End If
              Nini_species=Nini_species+1
            Else
              terminated(j)=.True.
            End If
          End Do
          set_u0=.False.
          If (Nini_species==0) Then
            Write (message,'(1x,a,2x,i6,a)') '***PROBLEMS: the code could not identify a single monitored species for frame ', i,&
                                            & '. Please review the settings for the &monitored_species block'
            Call info(message, 1)
            Call error_stop(' ')
          End If
        End If
      
        suma_i=0.0_wp
        traj_data%N_species=0
        Do j = 1, traj_data%Nmax_species
          If (.Not. terminated(j)) Then
            If (traj_data%species(i,j)%alive) Then
              Call rotation_vector_monitored_species(traj_data, ocf_data, i, j)
              Call orientational_correlation_term_monitored_species(traj_data, ocf_data, i, j, suma_i)  
            Else
              !terminated(j)=.True.
            End If
          End If  
        End Do
    
        If (i==traj_data%analysis%seg_indx(2,k)) Then
          If (traj_data%N_species /= 0) Then
            suma_i=suma_i/traj_data%N_species
            traj_data%analysis%variable(l,k)=suma_i
            If (ocf_data%print_all_intervals%stat) Then
              Write (iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
            End If
          End If  
          terminated=.False.
          set_u0=.True.
          If (ocf_data%print_all_intervals%stat) Then
            If ((traj_data%analysis%N_seg /=1) .And. (k /= traj_data%analysis%N_seg)) Then
              If (k /= traj_data%analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)         OCF' 
              End If
            End If  
          End If                      
        Else
          If ((traj_data%N_species) /= 0) Then
            suma_i=suma_i/traj_data%N_species
          Else  
            suma_i=0.0_wp
          End If
          traj_data%analysis%variable(l,k)=suma_i
          If (ocf_data%print_all_intervals%stat) Then
            Write (iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If
        End If
      End Do
    End Do
    
    If (ocf_data%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The OCF analysis for the multiple time intervals was printed to the "'&
                                 &//Trim(files(FILE_OCF_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The OCF analysis was printed to the "'//Trim(files(FILE_OCF_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time interval.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    
    Call average_segments(files, traj_data, FILE_OCF_AVG, 'OCF')
    If (traj_data%analysis%N_seg ==1 ) Then
      Write (message,'(1x,a)') 'WARNING: A single time interval was used to compute the average OCF! The computed STD&
                                & is zero. Use/Check the &data_analysis block to improve the statistics.'
      Call info(message, 1)
    End If
    
    If (.Not. ocf_data%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the OCF analysis for all time intervals,&
                                & the "print_all_intervals" directive (within the &ocf block) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%analysis%N_seg ==1 .And. (.Not. traj_data%analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_OCF_ALL)%filename)//'" and "'//Trim(files(FILE_OCF_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If
    
    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine orientational_correlation_function

  Subroutine rotation_vector_monitored_species(traj_data, ocf_data, i, j)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the rotation vector of the monitored species
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),    Intent(In   ) :: ocf_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j

    Integer(Kind=wi) :: indx1, indx2, indx3, k
    Logical          :: modified
    Real(Kind=wp), Dimension(3)  :: u12, u13


    indx1=traj_data%species(i,j)%list(1)
    indx2=traj_data%species(i,j)%list(2)
    indx3=traj_data%species(i,j)%list(3)
    
    If (Trim(ocf_data%u_definition%type) == 'bond_12') Then
      traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      Call check_PBC(traj_data%species(i,j)%u(:,1), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(ocf_data%u_definition%type) == 'bond_13') Then
      traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(traj_data%species(i,j)%u(:,1), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(ocf_data%u_definition%type) == 'bond_12-13') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Do k=1,3
        traj_data%species(i,j)%u(k,1)=u12(k)
        traj_data%species(i,j)%u(k,2)=u13(k)
      End Do
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
      traj_data%species(i,j)%u(:,2)=traj_data%species(i,j)%u(:,2)/norm2(traj_data%species(i,j)%u(:,2))
    Else If (Trim(ocf_data%u_definition%type) == 'bond_123') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Do k=1, 3
        traj_data%species(i,j)%u(k,1)=u12(k)+u13(k)
      End Do
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(ocf_data%u_definition%type) == 'plane') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call cross_product(u12, u13, traj_data%species(i,j)%u(:,1))
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    End If    
    
  End Subroutine rotation_vector_monitored_species  

  
  Subroutine orientational_correlation_term_monitored_species(traj_data, ocf_data, i, j, suma_i)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation for the relevant 
    ! species j at the MD frame i (cij)
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(ocf_type),    Intent(In   ) :: ocf_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j
    Real(Kind=wp),     Intent(InOut) :: suma_i
    
    Real(Kind=wp)     :: x, cij, c2ij, c0ij
    Logical           :: flag
    Integer(Kind=wi)  :: m
    
    If (traj_data%region%define%fread) Then
      m=traj_data%species(i,j)%list(1)
      Call within_region(traj_data, i, m, flag)
    Else
      flag=.True.
    End If
    
    If (flag) Then
      x=Dot_product(traj_data%species(i,j)%u(:,1),traj_data%species(i,j)%u0(:,1))
      traj_data%N_species=traj_data%N_species+1
      Select Case (ocf_data%legendre_order%value)  
        Case (1)
          cij=x
        Case (2)
          cij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
        Case (3)
          cij=(5.0_wp*(x)**3-3.0_wp*x)/2.0_wp
        Case (4)
          cij=(35.0_wp*(x)**4-30.0_wp*x**2+3.0_wp)/8.0_wp
      End Select  
      
      If (Trim(ocf_data%u_definition%type) == 'bond_12-13') Then
        x=Dot_product(traj_data%species(i,j)%u(:,2),traj_data%species(i,j)%u0(:,2))
        c0ij=cij
        Select Case (ocf_data%legendre_order%value)  
          Case (1)
            c2ij=x
          Case (2)
            c2ij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
          Case (3)
            c2ij=(5.0_wp*(x)**3-3.0_wp*x)/2.0_wp
          Case (4)
            c2ij=(35.0_wp*(x)**4-30.0_wp*x**2+3.0_wp)/8.0_wp
        End Select  
        cij=(c0ij+c2ij)/2.0_wp
      End If
      
      suma_i=suma_i+cij

    End If

  End Subroutine orientational_correlation_term_monitored_species    
  
  Subroutine cross_product(a, b, cross)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the cross_product 
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3)
    Real(Kind=wp), Intent(In   ) :: b(3)
    Real(Kind=wp), Intent(  Out) :: cross(3) 

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)

  End Subroutine cross_product   


  Subroutine read_orientational_chemistry(iunit, chemocf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &orientational_chemistry block
    !
    ! author    - i.scivetti February 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(chemocf_type), Intent(InOut)  :: chemocf_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &orientational_chemistry block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_orientational_chemistry" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_orientational_chemistry') Exit
      Call check_for_rubbish(iunit, '&orientational_chemistry')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='variable') Then
        Read (iunit, Fmt=*, iostat=io) word, chemocf_data%variable%type
        Call set_read_status(word, io, chemocf_data%variable%fread,&
                                     & chemocf_data%variable%fail,&
                                     & chemocf_data%variable%type)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, chemocf_data%print_all_intervals%stat
       Call set_read_status(word, io, chemocf_data%print_all_intervals%fread,&
                         & chemocf_data%print_all_intervals%fail)
                            
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings. See the "use_code.md" file.&
                                & Have you properly closed the block with "&end_orientational_chemistry"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_orientational_chemistry
  
  Subroutine compute_orientational_chemistry(files, traj_data, model_data, chemocf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the orientational chemistry along the trajectory
    !
    ! author    - i.scivetti Febraury 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data
    Type(chemocf_type), Intent(InOut) :: chemocf_data

    Integer(Kind=wi)   :: i, j, k, l, m
    Integer(Kind=wi)   :: iunit, ini_indx, ncouples
    Real(Kind=wp)      :: suma_i 
    Real(Kind=wp)      :: time
    Real(Kind=wp)      :: base_time
    Real(Kind=wp)      :: u(3,model_data%chem%N0%value), u0(3,model_data%chem%N0%value)

    Logical            :: set_s0, modified
    Logical            :: first_change(model_data%chem%N0%value)
    Integer(Kind=wp)   :: first_index(model_data%chem%N0%value)

    Character(Len=256) :: message
    
    Integer(Kind=wi)   :: indexes(2,model_data%chem%N0%value)
    Integer(Kind=wi)   :: ref_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: sites(4,traj_data%analysis%Ninterval,model_data%chem%N0%value)
    Integer(Kind=wi)   :: s1, s2
 
    If (chemocf_data%print_all_intervals%stat) Then
      ! Print header
      Open(Newunit=files(FILE_CHEM_OCF_ALL)%unit_no, File=files(FILE_CHEM_OCF_ALL)%filename,  Status='Replace')
      iunit=files(FILE_CHEM_OCF_ALL)%unit_no
      Write (iunit,'(a)') '#  Orientational chemistry for all the intervals' 
      Write (iunit,'(a)') '#  Time (ps)         OCF' 
    End If

    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1

    Do k= 1, traj_data%analysis%N_seg
      set_s0=.True.
      first_change=.True.
      l=0
      ini_indx=traj_data%analysis%seg_indx(1,k)
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        l=l+1
        If (Trim(chemocf_data%variable%type)=='special_pair') Then
        ! Obtain the index of the closest acceptor
          Do m = 1, model_data%chem%N0%value
            sites(1,l,m)=traj_data%track_chem%config(i,m)%indx
            sites(2,l,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
            sites(3,l,m)=traj_data%track_chem%config(i,m)%nn_indx(2)
            sites(4,l,m)=traj_data%track_chem%config(i,m)%nn_indx(3)
          End Do
        Else If (Trim(chemocf_data%variable%type)=='acceptor_donor_transfer_couple') Then
          If (set_s0) Then
            Do m = 1, model_data%chem%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%indx
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
            End Do  
            set_s0=.False.
          Else
            Do m = 1, model_data%chem%N0%value
              If (traj_data%track_chem%config(i,m)%indx/=indexes(2,m)) Then
                indexes(1,m)=indexes(2,m)
                indexes(2,m)=traj_data%track_chem%config(i,m)%indx
                If (first_change(m)) Then
                  Do j=1, l
                    sites(1,j,m)=ref_indx(m)
                    sites(2,j,m)=indexes(2,m)
                  End Do
                  first_change(m)=.False.
                  first_index(m)=l
                Else
                  sites(1,l,m)=indexes(1,m)
                  sites(2,l,m)=indexes(2,m)
                End If
              Else
                If (.Not. first_change(m)) Then
                  sites(1,l,m)=indexes(1,m)
                  sites(2,l,m)=indexes(2,m)
                End If
              End If
            End Do
          End If
        End If  
      End Do
      
      ! Set initial vector for the transfer couple at the start of the time interval
      Do m=1, model_data%chem%N0%value
        s1=sites(1,1,m)
        s2=sites(2,1,m)
        
        u0(:,m)=traj_data%config(ini_indx,s2)%r-traj_data%config(ini_indx,s1)%r
        Call check_PBC(u0(:,m), traj_data%box(ini_indx)%cell, traj_data%box(ini_indx)%invcell, 0.5_wp, modified)
        u0(:,m)=u0(:,m)/norm2(u0(:,m))
      End Do

      l=0
      ! Compute the orientational correlation function from the chaning chemistry
      base_time=(traj_data%analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        l=l+1
        time=(i-1)*traj_data%timestep%value
        suma_i=0.0_wp
        ncouples=0
        Do m=1, model_data%chem%N0%value
          s1=sites(1,l,m)
          s2=sites(2,l,m)

          u(:,m)=traj_data%config(i,s2)%r-traj_data%config(i,s1)%r
          Call check_PBC(u(:,m), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
          u(:,m)=u(:,m)/norm2(u(:,m))

          If (Trim(chemocf_data%variable%type)=='special_pair') Then
            Call orientational_correlation_term_transfer_couple(traj_data, i, s1, u(:,m), u0(:,m), suma_i, ncouples)
          Else If (Trim(chemocf_data%variable%type)=='acceptor_donor_transfer_couple') Then
            If (l<first_index(m)) Then
              Call orientational_correlation_term_transfer_couple(traj_data, i, s1, u(:,m), u0(:,m), suma_i, ncouples)
            Else
              Call orientational_correlation_term_transfer_couple(traj_data, i, s2, u(:,m), u0(:,m), suma_i, ncouples) 
            End If
          End If
        End Do 
        
        If (i==traj_data%analysis%seg_indx(2,k)) Then
          If (ncouples /= 0) Then
            suma_i=suma_i/ncouples
            traj_data%analysis%variable(l,k)=suma_i
            If (chemocf_data%print_all_intervals%stat) Then
              Write (iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
            End If
          End If
          If (chemocf_data%print_all_intervals%stat) Then
            If ((traj_data%analysis%N_seg /=1) .And. (k /= traj_data%analysis%N_seg)) Then
              If (k /= traj_data%analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)         OCF' 
              End If
            End If
           End If 
        Else
          If (ncouples /= 0) Then
            suma_i=suma_i/ncouples
          Else
            suma_i=0.0_wp
          End If
          traj_data%analysis%variable(l,k)=suma_i
          If (chemocf_data%print_all_intervals%stat) Then
            Write (iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If  
        End If
        
      End Do
      
    End Do  

    If (chemocf_data%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The CHEM_OCF analysis for the multiple time intervals was printed to the "'&
                                 &//Trim(files(FILE_CHEM_OCF_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The CHEM_OCF analysis was printed to the "'//Trim(files(FILE_CHEM_OCF_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time interval.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    ! Compute average
    Call average_segments(files, traj_data, FILE_CHEM_OCF_AVG, 'CHEM_OCF')
    If (traj_data%analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time interval was used to compute the average CHEM_OCF! The computed STD&
                              & is zero. Use/Check the &data_analysis block to improve the statistics.'
      Call info(message, 1)                        
    End If
    
    If (.Not. chemocf_data%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the CHEM_OCF analysis for all time intervals,&
                                & the "print_all_intervals" directive (within the &orientational_chemistry) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%analysis%N_seg ==1 .And. (.Not. traj_data%analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_CHEM_OCF_ALL)%filename)//'" and "'//Trim(files(FILE_CHEM_OCF_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If

    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine compute_orientational_chemistry 

  Subroutine orientational_correlation_term_transfer_couple(traj_data, i, s2, u, u0, suma_i, ncouples)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation for the relevant 
    ! transfer couple at the MD frame i (cij)
    !
    ! author    - i.scivetti Feb 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: s2
    Real(Kind=wp),     Intent(In   ) :: u(3)
    Real(Kind=wp),     Intent(In   ) :: u0(3)
    Real(Kind=wp),     Intent(InOut) :: suma_i
    Integer(Kind=wi),  Intent(InOut) :: ncouples
    
    Real(Kind=wp)     :: x, cij
    Logical           :: flag
    
    If (traj_data%region%define%fread) Then
      Call within_region(traj_data, i, s2, flag)
    Else
      flag=.True.
    End If

    If (flag) Then
      x=Dot_product(u,u0)
      ncouples=ncouples+1
      cij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
      suma_i=suma_i+cij
    End If

  End Subroutine orientational_correlation_term_transfer_couple  

  Subroutine check_orientational_chemistry(files, chemocf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &orientational_chemistry block
    !
    ! author    - i.scivetti Febraury 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(chemocf_type), Intent(InOut) :: chemocf_data

    Character(Len=256)  :: error_set
    Character(Len=256)  :: messages(2)

    error_set = '***ERROR in the &orientational_chemistry block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (chemocf_data%variable%fread) Then
      If (chemocf_data%variable%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "variable" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(chemocf_data%variable%type)/='special_pair'     .And. &
            Trim(chemocf_data%variable%type)/='acceptor_donor_transfer_couple')  Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    & 'Wrong input for "variable". Valid options:&
                                    & "special_pair" or "acceptor_donor_transfer_couple"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "variable" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "special_pair" or "acceptor_donor_transfer_couple"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If

    If (chemocf_data%print_all_intervals%fread) Then
      If (chemocf_data%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      chemocf_data%print_all_intervals%stat=.False.
    End If
    
  End Subroutine check_orientational_chemistry
  
  
End Module ocf
