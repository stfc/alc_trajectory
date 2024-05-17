!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that:
! - reads the SETTINGS file and defines the settings for analysis
! - checks correctness of defined directives
!
! Copyright   2023-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author      - i.scivetti   March 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module settings 

  Use atomic_model,      Only : model_type, &
                                geo_param_type, & 
                                geo_spec_type, &
                                short_dist_type, &
                                check_length_directive,&
                                check_model_settings,&
                                print_model_settings
                                
  Use constants,         Only : Bohr_to_A,  &
                                chemsymbol, & 
                                NPTE, &
                                max_at_species,&
                                max_components,&
                                max_unchanged_atoms
                                
  Use fileset,           Only : file_type, &
                                FILE_SET, &  
                                FILE_OUT, & 
                                refresh_out
  Use numprec,           Only : wi, &
                                wp
  Use process_data,      Only : capital_to_lower_case, &
                                check_for_rubbish, &
                                get_word_length
  Use unit_output,       Only : error_stop,&
                                info
  Use trajectory,        Only : traj_type, &
                                check_trajectory_settings

  Implicit None
  
  Private
  Public :: read_settings, check_settings

Contains

  Subroutine duplication_error(directive)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Aborts execution when duplication for
   ! a directive is found
   !
   ! author - i. scivetti  Feb 2023
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: directive

    Character(Len=256)  :: message

    Write (message,'(4a)') '***ERROR - Directive "', Trim(directive), '" is duplicated!'
    Call error_stop(message)

  End Subroutine duplication_error  

  Subroutine set_read_status(word, io, fread, fail, string)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to:
    !  - prevent duplication
    !  - define input directive is read by setting fread=.True. 
    !  - test if there was a problem with reading a directive, indicated by io/=0. This sets fail=.True.
    !
    ! author    - i.scivetti Febe 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=*), Intent(In   ) :: word
    Integer(Kind=wi), Intent(In   ) :: io
    Logical,          Intent(  Out) :: fread 
    Logical,          Intent(InOut) :: fail
    Character(Len=*), Optional, Intent(InOut) :: string

    If (fread)then
      Call duplication_error(word)
    Else
      fread=.True.
      If (io /= 0) Then
        fail=.True.
      End If
    End If

    If (Present(string)) then
      Call capital_to_lower_case(string)
    End If

  End Subroutine set_read_status 

  Subroutine read_settings(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read settings from SETTINGS file.
    ! Lines starting with # are ignored and assumed as comments. 
    ! If a directive is identified during the reading of the file, subroutine "set_read_fail" 
    ! assigns fread=.True. On the contrary, the subroutine assigns fail=.True. (fail=.False.) 
    ! if the format/syntax for the directive is correct (incorrect)
    ! If the directive is repeated the execution is aborted via subroutine duplication 
    ! 
    ! author        - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(InOut) :: files(:)
    Type(model_type),   Intent(InOut) :: model_data 
    Type(traj_type),    Intent(InOut) :: traj_data 
 
    Logical            :: safe
    Character(Len=256) :: word
    Integer(Kind=wi)   :: length, io, iunit
  
    Character(Len=256)  :: message

    Character(Len=32 )  :: set_file
    Character(Len=32 )  :: set_error

    set_file = Trim(files(FILE_SET)%filename)
    set_error = '***ERROR in the '//Trim(set_file)//' file.'

    ! Open the SETTINGS file with settings
    Inquire(File=files(FILE_SET)%filename, Exist=safe)
    
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(set_file), '(settings for analysis) not found'
      Call error_stop(message)
    Else
      Open(Newunit=files(FILE_SET)%unit_no, File=Trim(set_file), Status='old')
      iunit=files(FILE_SET)%unit_no 
    End If

     Read (iunit, Fmt=*, iostat=io) word
     ! If nothing is found, complain and abort
     If (is_iostat_end(io)) Then
       Write (message,'(3(1x,a))') Trim(set_error), Trim(set_file), 'file seems to be empty?. Please check'
       Call error_stop(message)
     End If
     ! Check header has "#" as the first character 
     If (word(1:1)/='#') Then
       Write (message,'(4(1x,a))') Trim(set_error), 'Heading comment in file', Trim(set_file), & 
                                  'is required and MUST be preceded with the symbol "#"'
       Call error_stop(message)
     End If

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Exit
      end If
      Call check_for_rubbish(iunit, Trim(set_file)) 
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
        ! Do nothing if line is a comment of we have an empty line
        Read (iunit, Fmt=*, iostat=io) word

      ! Model related variable 
      Else If (word(1:length) == 'change_chemistry') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%change_chemistry%stat
        Call set_read_status(word, io, model_data%change_chemistry%fread, model_data%change_chemistry%fail)

      Else If (word(1:length) == 'input_geometry_format') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%input_geometry_format%type
        Call set_read_status(word, io, model_data%input_geometry_format%fread, model_data%input_geometry_format%fail, &
                           & model_data%input_geometry_format%type)

      Else If (word(1:length) == 'cell_units') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%config%cell_units%type
        Call set_read_status(word, io, model_data%config%cell_units%fread, model_data%config%cell_units%fail)
        
      Else If (word(1:length) == 'position_units') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%config%position_units%type
        Call set_read_status(word, io, model_data%config%position_units%fread, model_data%config%position_units%fail)
        
      Else If (word(1:length) == '&input_composition') Then
        Read (iunit, Fmt=*, iostat=io) model_data%input_composition%invoke%type
        Call set_read_status(word, io, model_data%input_composition%invoke%fread, model_data%input_composition%invoke%fail)
        ! Read information inside the block
        Call read_input_composition(iunit, model_data)

      Else If (word(1:length) == '&simulation_cell') Then
        Read (iunit, Fmt=*, iostat=io) model_data%config%simulation_cell%type
        Call set_read_status(word, io, model_data%config%simulation_cell%fread, model_data%config%simulation_cell%fail)
        ! Read information inside the block
        Call read_input_cell(iunit, model_data)  

      Else If (word(1:length) == '&search_chemistry') Then
        Read (iunit, Fmt=*, iostat=io) model_data%search_chemistry%type
        Call set_read_status(word, io, model_data%search_chemistry%fread, model_data%search_chemistry%fail)
        ! Read information inside the block
        Call read_chemistry_criteria(iunit, model_data)

      ! Trajectory related variables
      Else If (word(1:length) == 'ensemble') Then 
        Read (iunit, Fmt=*, iostat=io) word, traj_data%ensemble%type
        Call set_read_status(word, io, traj_data%ensemble%fread, traj_data%ensemble%fail, &
                           & traj_data%ensemble%type)

      Else If (word(1:length) == 'timestep') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%timestep%tag, traj_data%timestep%value,&
                                      &traj_data%timestep%units
        Call set_read_status(word, io, traj_data%timestep%fread, traj_data%timestep%fail)

      Else If (word(1:length) == '&data_analysis') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%analysis%invoke%type
        Call set_read_status(word, io, traj_data%analysis%invoke%fread, traj_data%analysis%invoke%fail)
        Call read_settings_for_analysis(iunit, traj_data)
        
      Else If (word(1:length) == 'print_retagged_trajectory') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%print_retagged_trajectory%stat
       Call set_read_status(word, io, traj_data%print_retagged_trajectory%fread, traj_data%print_retagged_trajectory%fail)
        
      Else If (word(1:length) == '&monitored_species') Then
        Read (iunit, Fmt=*, iostat=io) model_data%config%monitored_species%type
        Call set_read_status(word, io, model_data%config%monitored_species%fread, model_data%config%monitored_species%fail)
        !Read information inside the block
        Call read_monitored_species(iunit, model_data)

      Else If (word(1:length) == '&ocf') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%ocf%invoke%type
        Call set_read_status(word, io, traj_data%ocf%invoke%fread, traj_data%ocf%invoke%fail)
        !Read information inside the block
        Call read_ocf(iunit, traj_data)

      Else If (word(1:length) == '&region') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%region%define%type
        Call set_read_status(word, io, traj_data%region%define%fread, traj_data%region%define%fail)
        !Read information inside the block
        Call read_region(iunit, traj_data)

      Else If (word(1:length) == '&msd') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%msd%invoke%type
        Call set_read_status(word, io, traj_data%msd%invoke%fread, traj_data%msd%invoke%fail)
        !Read information inside the block
        Call read_msd(iunit, traj_data)

      Else If (word(1:length) == '&selected_nn_distances') Then
        Read (iunit, Fmt=*, iostat=io) model_data%nndist%invoke%type
        Call set_read_status(word, io, model_data%nndist%invoke%fread, model_data%nndist%invoke%fail)
        !Read information inside the block
        Call read_selected_nn_distances(iunit, model_data%nndist)

      Else If (word(1:length) == '&rdf') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%rdf%invoke%type
        Call set_read_status(word, io, traj_data%rdf%invoke%fread, traj_data%rdf%invoke%fail)
        !Read information inside the block
        Call read_rdf(iunit, traj_data)

      Else If (word(1:length) == '&coord_distrib') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%coord_distrib%invoke%type
        Call set_read_status(word, io, traj_data%coord_distrib%invoke%fread, traj_data%coord_distrib%invoke%fail)
        !Read information inside the block
        Call read_coord_distrib(iunit, traj_data)

      Else If (word(1:length) == '&track_unchanged_chemistry') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%unchanged%invoke%type
        Call set_read_status(word, io, traj_data%unchanged%invoke%fread, traj_data%unchanged%invoke%fail)
        !Read information inside the block
        Call read_track_unchanged_chemistry(iunit, traj_data)
        
      Else If (word(1:length) == '&lifetime') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%lifetime%invoke%type
        Call set_read_status(word, io, traj_data%lifetime%invoke%fread, traj_data%lifetime%invoke%fail)
        !Read information inside the block
        Call read_lifetime(iunit, traj_data)

      Else If (word(1:length) == '&orientational_chemistry') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%chem_ocf%invoke%type
        Call set_read_status(word, io, traj_data%chem_ocf%invoke%fread, traj_data%chem_ocf%invoke%fail)
        !Read information inside the block
        Call read_orientational_chemistry(iunit, traj_data)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Directive not recognised. Inform and kill 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      Else
        If (word(1:1)=='&') Then
          Write (message,'(1x,4a)') Trim(set_error), ' Unknown directive found: "', Trim(word),&
                                  &'. Do you use "&" to define a block? If so,&
                                  & make sure the block is valid and has right syntax.'
        Else
          Write (message,'(1x,a)') Trim(set_error)//' Unknown directive found: "'//Trim(word)//'".&
                                  & Have you correctly defined the previous directives? Have you forgotten something maybe?'
        End If 
        Call error_stop(message)
      End If

    End Do
    ! Close file
    Close(files(FILE_SET)%unit_no)

  End Subroutine read_settings

   Subroutine check_end(io, string)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Subroutine to check if there is missing data and the end of the file
     ! has been reached
     !
     ! author    - i.scivetti Dec 2021
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Integer,          Intent(In   ) :: io
     Character(Len=*), Intent(In   ) :: string
 
     Character(Len=256) :: messages(2)
 
     If (is_iostat_end(io))Then
       Call info(' ', 1)
       Write (messages(1),'(1x,2a)') '*** ERROR in ', Trim(string)
       Write (messages(2),'(1x,2a)') 'End of file is detected. It seems there is missing data or the block is not&
                                   & closed properly. Please check'
       Call info(messages, 2)
       Call error_stop(' ')
     End If
 
   End Subroutine check_end

  Subroutine read_input_cell(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the simulation cell vectors of the input structure
    ! defined in &simulation_cell
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),   Intent(In   ) :: iunit
    Type(model_type),   Intent(InOut) :: model_data
   
    Integer(Kind=wi)   :: io, i, j
    Character(Len=64 ) :: error_simulation_cell
    Character(Len=256) :: messages(2), word
    Logical            :: endblock

    error_simulation_cell = '***ERROR in &simulation_cell of SETTINGS file'
    Write (messages(1),'(a)') error_simulation_cell

    i=1
    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&simulation_cell')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call check_for_rubbish(iunit, '&simulation_cell') 
          Read (iunit, Fmt=*, iostat=io) (model_data%config%cell(i,j), j=1,3)
          If (io/=0) Then
            Write (messages(2),'(a,i2)') 'Problems with the definition of cell vector', i
            Call info(messages, 2)
            Call error_stop(' ')
          End If
          i=i+1
        Else
          Write (messages(2),'(1x,a)') 'End of block found! Not all the cell vectors for the&  
                                     & input structure have been defined. Please check.'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End If
    End Do
 
    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&simulation_cell')
      Call capital_to_lower_case(word)
      If (word /= '&end_simulation_cell') Then
        If (word(1:1) /= '#') Then
          Write (messages(2),'(a)') 'Block for cell vectors must be closed with&
                                  & sentence &end_simulation_cell.'
          Call info(messages,2)
          Call error_stop(' ')
        End If
      Else
          endblock=.True.
      End If
    End Do

  End Subroutine read_input_cell

  Subroutine read_chemistry_criteria(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the definitions to search and identify changes chemisty
    ! Information is read from the &search_chemistry block
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

   Integer(Kind=wi)   :: io, length
   Character(Len=256) :: message, word
   Character(Len=256) :: set_error
    

    set_error = '***ERROR in &search_chemistry of SETTINGS file.'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_search_chemistry" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_search_chemistry') Exit
      Call check_for_rubbish(iunit, '&search_chemistry')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'total_number') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%chem%N0%value
        Call set_read_status(word, io, model_data%chem%N0%fread, model_data%chem%N0%fail)

       Else If (word(1:length) == '&acceptor_criteria') Then
         Read (iunit, Fmt=*, iostat=io) word
         Call set_read_status(word, io, model_data%chem%acceptor%criteria%fread, model_data%chem%acceptor%criteria%fail)
         Call read_environment_settings(iunit, model_data)

       Else If (Trim(word)=='&bonding_criteria') Then
            Read (iunit, Fmt=*, iostat=io) word
            Call set_read_status(word, io, model_data%chem%bonds%criteria%fread,&
                  & model_data%chem%bonds%criteria%fail)
            model_data%chem%bonds%criteria%stat = .True.
         Call read_bonding_criteria(iunit, model_data)
 
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&End_chemistry"?'
        Call error_stop(message)
      End If

    End Do
   
  End Subroutine read_chemistry_criteria

  Subroutine read_bonding_criteria(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settings from the &bonding_criteria block and 
    ! identify the species to be tracked
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: messages(3), word
    Character(Len=256) :: set_error
     

    set_error = '***ERROR in the "&bonding_criteria" sub-block (within &search_chemistry).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (messages(1),'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_bonding_criteria" to close the block.&
                                  & Check if directives are set correctly.'         
        Call info(messages, 1)                          
        Call error_stop(' ') 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_bonding_criteria') Exit
      Call check_for_rubbish(iunit, '&bonding_criteria')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'only_element') Then 
        Read (iunit, Fmt=*, iostat=io) word, model_data%chem%bonds%species%type
        Call set_read_status(word, io, model_data%chem%bonds%species%fread, model_data%chem%bonds%species%fail)

      Else If (Trim(word)=='cutoff') Then
          Read (iunit, Fmt=*, iostat=io) model_data%chem%bonds%cutoff%tag, model_data%chem%bonds%cutoff%value,&
                                         model_data%chem%bonds%cutoff%units
         Call set_read_status(word, io, model_data%chem%bonds%cutoff%fread,  model_data%chem%bonds%cutoff%fail)

      Else If (Trim(word)=='number_of_bonds') Then
          Read (iunit, Fmt=*, iostat=io) word, model_data%chem%bonds%N0%value
          Call set_read_status(word, io, model_data%chem%bonds%N0%fread,&
               & model_data%chem%bonds%N0%fail)
        
       Else If (word(1:length) == '&possible_extra_bonds') Then
          Read (iunit, Fmt=*, iostat=io) model_data%extra_bonds%invoke%type
          Call set_read_status(word, io, model_data%extra_bonds%invoke%fread, model_data%extra_bonds%invoke%fail)
          ! Read information inside the block
          Call read_extra_bonding(iunit, model_data)
        
      Else
        Write (messages(1),'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//'" is not recognised as a valid settings.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the sub-block with "&end_bonding_criteria"?'
        Write (messages(3),'(1x,a)') 'Have you included the units for "cutoff"?'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

    End Do
   
  End Subroutine read_bonding_criteria
  
  Subroutine read_environment_settings(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the defined criteria to explore the  environment of 
    ! selected atoms
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)   :: io, j, length
    Character(Len=256) :: messages(3), word
    Character(Len=256) :: set_error
    

    set_error = '***ERROR in the sub-block &acceptor_criteria within &search_chemistry (SETTINGS file).'
    Write (messages(2),'(1x,a)') 'Have you properly closed the sub-block with "&end_acceptor_criteria"?'
    Write (messages(3),'(1x,a)') 'Have you set the units for directive "cutoff"?'                         

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (messages(1),'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_acceptor_criteria" to close the block.&
                                  & Check if directives are set correctly.' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_acceptor_criteria') Exit
      Call check_for_rubbish(iunit, '&acceptor_criteria')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='include_tags') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%chem%acceptor%N0_incl
        Call prevent_segmentation(iunit, io, word, model_data%chem%acceptor%N0_incl,'max_components', max_components, set_error)
        Read (iunit, Fmt=*, iostat=io) word, model_data%chem%acceptor%N0_incl,&
                                          & (model_data%chem%acceptor%tg_incl(j), j = 1, model_data%chem%acceptor%N0_incl)
        Call set_read_status(word, io, model_data%chem%acceptor%info_include%fread, model_data%chem%acceptor%info_include%fail)
        model_data%chem%acceptor%info_include%stat = .True.

      Else If (Trim(word)=='exclude_pairs') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%chem%acceptor%N0_excl 
        Call prevent_segmentation(iunit, io, word, model_data%chem%acceptor%N0_excl,'max_components', max_components, set_error)
        Read (iunit, Fmt=*, iostat=io) word, model_data%chem%acceptor%N0_excl,&
                                          & (model_data%chem%acceptor%tg_excl(j), j = 1, model_data%chem%acceptor%N0_excl)
        Call set_read_status(word, io, model_data%chem%acceptor%info_exclude%fread, model_data%chem%acceptor%info_exclude%fail)
        model_data%chem%acceptor%info_exclude%stat = .True.

      Else If (Trim(word)=='cutoff') Then
         Read (iunit, Fmt=*, iostat=io) model_data%chem%acceptor%cutoff%tag, model_data%chem%acceptor%cutoff%value,&
                                        model_data%chem%acceptor%cutoff%units 
         Call set_read_status(word, io, model_data%chem%acceptor%cutoff%fread, model_data%chem%acceptor%cutoff%fail)

      Else
        Write (messages(1),'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//&
                                    &'" is not recognised as a valid settings.'
        Call info(messages, 3)
        Call error_stop(' ') 
      End If
    End Do
   
  End Subroutine read_environment_settings

  Subroutine read_extra_bonding(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to extra bond settings
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)  ::  io, i

    Character(Len=256)  :: word, messages(3)
    Character(Len=256)  :: error_block
    Logical  :: error, endblock, fread 

    error= .False.
    error_block = '***ERROR in &possible_extra_bonds (inside &bonding_criteria).' 
    Write (messages(1),'(a)') error_block 

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&possible_extra_bonds')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&possible_extra_bonds')
      End If
    End Do

    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, model_data%extra_bonds%N0

    If (Trim(word) /= 'types_of_bonds') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "types_of_bonds" is expected.'
      error=.True.
    End If 

    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "type_of_bonds"'
      error=.True.
    Else
      If (model_data%extra_bonds%N0<1) Then
        Write (messages(2),'(a)') 'The "type_of_bonds" directive MUST BE >= 1'
        error=.True.
      ElseIf (model_data%extra_bonds%N0>max_components) Then
        Write (messages(2),'(a,i3,a)') 'Are you sure you want to consider more than ', max_components,&
                                    & ' for "type_of_bonds"? Please check'
        error=.True.
      End If
    End If

    If (error) Then
      Call info(messages,2) 
      Call error_stop(' ')
    End If

    i=1
    Do While (i <= model_data%extra_bonds%N0)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&possible_extra_bonds')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&possible_extra_bonds')
  
        Read (iunit, Fmt=*, iostat=io) model_data%extra_bonds%tg1(i), model_data%extra_bonds%tg2(i),   &
                                       model_data%extra_bonds%bond(i)%value, model_data%extra_bonds%bond(i)%units

        If (io/=0) Then
          If (Trim(model_data%extra_bonds%tg1(i)) == '&end') Then
            Write (messages(2),'(2(a,i2),a)') 'Missing specification for extra bonds. Only ',  i-1,&
                                    &' species set out of ', model_data%extra_bonds%N0, ' (types_of_bonds)'
            Write (messages(3),'(a)') 'Please check. What is the value set for "type_of_bonds"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          Else  
            Write (messages(2),'(a,i3)') 'Problems to read bonding criteria ',  i
            Write (messages(3),'(a)') 'Please check. What is the value set for "type_of_bonds"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          End If
        Else
          Write (messages(1),'(a,i3)') Trim(error_block)//' Problems to read bonding criteria ',  i
          model_data%extra_bonds%bond(i)%fread= .True.
          Call check_length_directive(model_data%extra_bonds%bond(i), messages(1), .True., 'inblock')
        End If

        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&possible_extra_bonds')
      Call capital_to_lower_case(word)
      If (word /= '&end_possible_extra_bonds') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=model_data%extra_bonds%N0) Then 
            Write (messages(2),'(a)') 'Number of extra bonds specified is larger than&
                                     & the value given by directive "type_of_bonds"'
          Else
            Write (messages(2),'(a)') 'Block must be closed with sentence &end_possible_extra_bonds. Please check. Is the&
                                     & number of defined bond criteria the same as set for directive "type_of_bonds"?'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If  
    End Do
    
  End Subroutine read_extra_bonding

  Subroutine read_monitored_species(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the defined criteria to explore the selected species
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type), Intent(InOut) :: model_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: messages(3), word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &monitored_species block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (messages(1),'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_monitored_species" to close the block.&
                                  & Check if directives are set correctly.'         
        Call info(messages, 1)                          
        Call error_stop(' ') 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_monitored_species') Exit
      Call check_for_rubbish(iunit, '&monitored_species')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (word(1:length) == 'compute_amount') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%species_definition%compute_amount%stat
        Call set_read_status(word, io, model_data%species_definition%compute_amount%fread,&
                           & model_data%species_definition%compute_amount%fail)
      
      Else If (Trim(word)=='name') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%species_definition%name%type
        Call set_read_status(word, io, model_data%species_definition%name%fread, model_data%species_definition%name%fail)

      Else If (Trim(word)=='reference_tag') Then
        Read (iunit, Fmt=*, iostat=io) word, model_data%species_definition%reference_tag%type
        Call set_read_status(word, io, model_data%species_definition%reference_tag%fread,&
                           & model_data%species_definition%reference_tag%fail)

      Else If (Trim(word)=='bond_cutoff') Then
        Read (iunit, Fmt=*, iostat=io) model_data%species_definition%bond_cutoff%tag,  &
                                       model_data%species_definition%bond_cutoff%value, &
                                       model_data%species_definition%bond_cutoff%units
                                       
        Call set_read_status(word, io, model_data%species_definition%bond_cutoff%fread,&
                           & model_data%species_definition%bond_cutoff%fail)

      Else If (word(1:length) == '&atomic_components') Then
        Read (iunit, Fmt=*, iostat=io) word
        Call set_read_status(word, io, model_data%species_definition%atomic_components%fread,&
                           & model_data%species_definition%atomic_components%fail)
        Call read_components_monitored_species(iunit, model_data)

      Else If (word(1:length) == '&intramol_stat_settings') Then
        Read (iunit, Fmt=*, iostat=io) model_data%species_definition%intra_geom%invoke%type
        Call set_read_status(word, io, model_data%species_definition%intra_geom%invoke%fread, &
                            & model_data%species_definition%intra_geom%invoke%fail, &
                            & model_data%species_definition%intra_geom%invoke%type)
        model_data%species_definition%intra_geom%tag='intramol_stat_settings'
        Call read_geom_param_monitored_species(iunit, model_data%species_definition%intra_geom)

      Else If (word(1:length) == '&intermol_stat_settings') Then
        Read (iunit, Fmt=*, iostat=io) model_data%species_definition%inter_geom%invoke%type
        Call set_read_status(word, io, model_data%species_definition%inter_geom%invoke%fread, &
                            & model_data%species_definition%inter_geom%invoke%fail, &
                            & model_data%species_definition%inter_geom%invoke%type)
        model_data%species_definition%inter_geom%tag='intermol_stat_settings'
        Call read_geom_param_monitored_species(iunit, model_data%species_definition%inter_geom)

      Else
        Write (messages(1),'(1x,a)') Trim(set_error)//' Directive "'//Trim(word)//&
                                  &'" is not recognised as a valid settings. See the "use_code.md" file.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the block with "&end_monitored_species"?'
        Write (messages(3),'(1x,a)') 'Have you included the units for "bond_cutoff"?'
        Call info(messages, 3)
        Call error_stop(' ')
      End If

    End Do
    
  End Subroutine read_monitored_species

  Subroutine read_geom_param_monitored_species(iunit, T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read blocks with parameters
    ! for the statistical analysis of geometry quantitites 
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),      Intent(In   ) :: iunit
    Type(geo_spec_type),   Intent(InOut) :: T 
    
    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the "&'//Trim(T%tag)//'" block (SETTINGS file).'
    
    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_'//Trim(T%tag)//'" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_'//Trim(T%tag)) Exit
      Call check_for_rubbish(iunit, '&'//Trim(T%tag))

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If ((Trim(word)=='only_ref_tag_as_nn') .And. (Trim(T%tag)=='intermol_stat_settings')) Then
        Read (iunit, Fmt=*, iostat=io) word, T%only_ref_tags_as_nn%stat
        Call set_read_status(word, io, T%only_ref_tags_as_nn%fread, T%only_ref_tags_as_nn%fail)

      Else If (Trim(word)=='&distance_parameters') Then
        Read (iunit, Fmt=*, iostat=io) T%dist%invoke%type
        Call set_read_status(word, io, T%dist%invoke%fread, T%dist%invoke%fail, T%dist%invoke%type)
        T%dist%name='distance_parameters'
        Call read_geom_param(iunit, T%tag, T%dist)

      Else If (Trim(word)=='&angle_parameters') Then
        Read (iunit, Fmt=*, iostat=io) T%angle%invoke%type
        Call set_read_status(word, io, T%angle%invoke%fread, T%angle%invoke%fail, T%angle%invoke%type)
        T%angle%name='angle_parameters'
        Call read_geom_param(iunit, T%tag, T%angle)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_'//Trim(T%tag)//'"?'
        Call error_stop(message)
      End If

    End Do
  
  End Subroutine read_geom_param_monitored_species
  
  
  Subroutine read_geom_param(iunit, inblock, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read distance and angle settings
    ! for statistics of the monitored species
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),     Intent(In   ) :: iunit
    Character(*),         Intent(In   ) :: inblock
    Type(geo_param_type), Intent(InOut) :: M 
    
    Integer(Kind=wi)   :: io, length, i
    Character(Len=256) :: message, word
    Character(Len=256) :: messages(2)
    Character(Len=256) :: set_error
    
    M%delta%tag='delta'
    
    set_error = '***ERROR in "&'//Trim(M%name)//'" within the "&'//Trim(inblock)//'" block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_'//Trim(M%name)//'" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_'//Trim(M%name)) Exit
      Call check_for_rubbish(iunit, '&'//Trim(M%name))

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word
 
      Else If (Trim(word)=='species') Then
        If (Trim(M%name)=='distance_parameters') Then
          M%nspecies=2
        Else If (Trim(M%name)=='angle_parameters') Then
          M%nspecies=3
        End If
        Read (iunit, Fmt=*, iostat=io) M%tag_species%type, (M%species(i), i=1, M%nspecies) 
        Call set_read_status(word, io, M%tag_species%fread, M%tag_species%fail, M%tag_species%type)

      Else If (Trim(word)=='lower_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%lower_bound%tag, M%lower_bound%value, M%lower_bound%units 
         Call set_read_status(word, io, M%lower_bound%fread, M%lower_bound%fail)

      Else If (Trim(word)=='upper_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%upper_bound%tag, M%upper_bound%value, M%upper_bound%units 
         Call set_read_status(word, io, M%upper_bound%fread, M%upper_bound%fail)

      Else If (Trim(word)=='delta') Then
         Read (iunit, Fmt=*, iostat=io) M%delta%tag, M%delta%value, M%delta%units 
         Call set_read_status(word, io, M%delta%fread, M%delta%fail)

      Else
        Write (messages(1),'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the block with "&end_'//Trim(M%name)//'"? &
                                & Have you defined the directives correctly? See the "use_code.md" file'
        Call info (messages, 2)
        Call error_stop(' ')
      End If
    End Do
  
  End Subroutine read_geom_param

  Subroutine read_selected_nn_distances(iunit, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read parameters from the
    ! &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),      Intent(In   ) :: iunit
    Type(short_dist_type), Intent(InOut) :: M 
    
    Integer(Kind=wi)   :: io, length, i
    Character(Len=256) :: message, word
    Character(Len=256) :: messages(2)
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in "&selected_nn_distances" block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_selected_nn_distances" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_selected_nn_distances') Exit
      Call check_for_rubbish(iunit, '&end_selected_nn_distances')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word
 
      Else If (Trim(word)=='reference_species') Then
        Read (iunit, Fmt=*, iostat=io) M%tag_reference_species%type, M%reference_species
        Call set_read_status(word, io, M%tag_reference_species%fread, M%tag_reference_species%fail, M%tag_reference_species%type)

      Else If (Trim(word)=='nn_species') Then
        Read (iunit, Fmt=*, iostat=io) M%tag_nn_species%type, M%num_nn_species
        Call prevent_segmentation(iunit, io, M%tag_nn_species%type, M%num_nn_species,&
                                & 'max_components', max_components, set_error)
        M%nn_species=' '
        Read (iunit, Fmt=*, iostat=io) M%tag_nn_species%type, M%num_nn_species, (M%nn_species(i), i=1, M%num_nn_species) 
        Call set_read_status(word, io, M%tag_nn_species%fread, M%tag_nn_species%fail, M%tag_nn_species%type)
        
      Else If (Trim(word)=='lower_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%lower_bound%tag, M%lower_bound%value, M%lower_bound%units 
         Call set_read_status(word, io, M%lower_bound%fread, M%lower_bound%fail)

      Else If (Trim(word)=='upper_bound') Then
         Read (iunit, Fmt=*, iostat=io) M%upper_bound%tag, M%upper_bound%value, M%upper_bound%units 
         Call set_read_status(word, io, M%upper_bound%fread, M%upper_bound%fail)

      Else If (Trim(word)=='dr') Then
         Read (iunit, Fmt=*, iostat=io) M%dr%tag, M%dr%value, M%dr%units 
         Call set_read_status(word, io, M%dr%fread, M%dr%fail)

      Else
        Write (messages(1),'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.'
        Write (messages(2),'(1x,a)') 'Have you properly closed the block with "&end_selected_nn_distances"? &
                                & Have you defined the directives correctly? See the "use_code.md" file'
        Call info (messages, 2)
        Call error_stop(' ')
      End If
    End Do
  
  End Subroutine read_selected_nn_distances
  
  
  Subroutine read_components_monitored_species(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to extra bond settings
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data 

    Integer(Kind=wi)  ::  io, i

    Character(Len=256)  :: word, messages(3)
    Character(Len=256)  :: error_block
    Logical             :: endblock, fread 

    error_block = '***ERROR in the &atomic_components sub-block (wihtin &monitored_species)'
    Write (messages(1),'(a)') error_block 

    fread= .True.
    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&atomic_components (within &monitored_species)')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&atomic_components')
      End If
    End Do

    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, model_data%species_definition%num_components

    If (Trim(word) /= 'number_components') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "number_components" is expected.'
      Call info(messages, 2) 
      Call error_stop(' ')
    End If 

    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number_components"'
      Call info(messages, 2) 
      Call error_stop(' ')
    Else
      If (model_data%species_definition%num_components<1) Then
        Write (messages(2),'(a)') 'The "number_components" directive MUST BE >= 1'
        Call info(messages, 2) 
        Call error_stop(' ')
      ElseIf (model_data%species_definition%num_components>max_at_species) Then
        Write (messages(2),'(a,i3,a)') 'Are you sure you want to consider more than ', max_at_species,&
                                    & ' for "number_components"? Please check'
        Write (messages(3),'(a)') 'If you are sure of what you are doing, look for the parameter "max_at_species" in the code,&
                                  & increase its value as needed and recompile.'
        Call info(messages, 3)
        Call error_stop(' ')
      End If
    End If

    i=1
    Do While (i <= model_data%species_definition%num_components)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&atomic_components (within &monitored_species)')
      If (word(1:1)/='#') Then
        Call check_for_rubbish(iunit, '&atomic_components')
  
        Read (iunit, Fmt=*, iostat=io) model_data%species_definition%element(i), model_data%species_definition%N0_element(i)

        If (io/=0) Then
          If (Trim(model_data%species_definition%element(i)) == '&end') Then
            Write (messages(2),'(2(a,i2),a)') 'Missing specification for extra bonds. Only ',  i-1,&
                                    &' species set out of ', model_data%species_definition%num_components,&
                                    & ' (number_components)'
            Write (messages(3),'(a)') 'Please check. What is the value set for "number_components"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          Else  
            Write (messages(2),'(a,i3)') 'Problems to read the species component ',  i      
            Write (messages(3),'(a)') 'Please check. What is the value set for "number_components"?'                       
            Call info(messages, 3) 
            Call error_stop(' ')
          End If 
        End If

        i=i+1
      End If
    End Do 

    endblock=.False.

    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&atomic_components (within &monitored_species)')
      Call capital_to_lower_case(word)
      If (word /= '&end_atomic_components') Then
        If (word(1:1) /= '#') Then
          If ((i-1)/=model_data%species_definition%num_components) Then 
            Write (messages(2),'(a)') 'Number of extra bonds specified is larger than&
                                     & the value given by directive "number_components"'
          Else
            Write (messages(2),'(a)') 'Block must be closed with sentence &end_atomic_components. Please check. Is the&
                                     & number of defined bond criteria the same as set for directive "number_components"?'
          End If   
          Call info(messages,2) 
          Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If  
    End Do
    
  End Subroutine read_components_monitored_species

  Subroutine read_region(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the &region block. This block defines the portion of the
    ! system to be analysed.
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(traj_type), Intent(InOut)  :: traj_data 

    Integer(Kind=wi)   :: io, length, k
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &region block (SETTINGS file).'
 
    traj_data%region%number=0
    
    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_region" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_region') Exit
      Call check_for_rubbish(iunit, '&region')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='delta_x') Then
        traj_data%region%number(1)=traj_data%region%number(1)+1
        k=traj_data%region%number(1)
        Read (iunit, Fmt=*, iostat=io) traj_data%region%invoke(1,k)%type, &
                                    & traj_data%region%domain(1,1,k),     &
                                    & traj_data%region%domain(1,2,k),     &
                                    & traj_data%region%inout(1,k)
        Call set_read_status(word, io, traj_data%region%invoke(1,k)%fread, &
                            & traj_data%region%invoke(1,k)%fail, traj_data%region%invoke(1,k)%type)
         
      Else If (Trim(word)=='delta_y') Then
        traj_data%region%number(2)=traj_data%region%number(2)+1
        k=traj_data%region%number(2)
        Read (iunit, Fmt=*, iostat=io) traj_data%region%invoke(2,k)%type, &
                                    & traj_data%region%domain(2,1,k),     &
                                    & traj_data%region%domain(2,2,k),     &
                                    & traj_data%region%inout(2,k)
        Call set_read_status(word, io, traj_data%region%invoke(2,k)%fread, &
                            & traj_data%region%invoke(2,k)%fail, traj_data%region%invoke(2,k)%type)

      Else If (Trim(word)=='delta_z') Then
        traj_data%region%number(3)=traj_data%region%number(3)+1
        k=traj_data%region%number(3)
        Read (iunit, Fmt=*, iostat=io) traj_data%region%invoke(3,k)%type, &
                                    & traj_data%region%domain(3,1,k),     &
                                    & traj_data%region%domain(3,2,k),     &
                                    & traj_data%region%inout(3,k)
        Call set_read_status(word, io, traj_data%region%invoke(3,k)%fread, &
                            & traj_data%region%invoke(3,k)%fail, traj_data%region%invoke(3,k)%type)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_region"?'
        Call error_stop(message)
      End If

    End Do
    
    ! Assing to 1 if not read
    Do k = 1, 3
      If (traj_data%region%number(k)==0) Then
        traj_data%region%number(k)=1
      End If
    End Do
    
  End Subroutine read_region

  Subroutine read_track_unchanged_chemistry(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settings to track chemically unchanged species 
    ! along the trajectory
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(traj_type),  Intent(InOut) :: traj_data 

    Integer(Kind=wi)   :: io, length, j
    Character(Len=256) :: message, messages(2)
    Character(Len=256) :: word
    Character(Len=256) :: set_error
    Logical :: error, fread
    
    set_error = '***ERROR in the &track_unchanged_chemistry block (SETTINGS file).'
    error=.False.
    fread= .True.

    Do While (fread)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&track_unchanged_chemistry')
      If (word(1:1)/='#') Then
        fread=.False.
        Call check_for_rubbish(iunit, '&track_unchanged_chemistry')
      End If
    End Do

    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, traj_data%unchanged%N0
    If (Trim(word) /= 'number') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "number" is expected to be defined first'
      error=.True.
    End If 

    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "number"'
      error=.True.
    Else
      If (traj_data%unchanged%N0<1) Then
        Write (messages(2),'(a)') 'The "number" directive MUST BE >= 1'
        error=.True.
      End If  
      If (traj_data%unchanged%N0>max_components) Then
        Write (messages(2),'(a,i3,a)') 'Directive number: are you sure you want to consider more than ', max_components,&
                                    & '? Please check'
        error=.True.
      End If
      If (traj_data%unchanged%N0>max_unchanged_atoms) Then
        Write (messages(2),'(a,i3,a)') 'Directive "number": the user cannot track more than ', max_unchanged_atoms,&
                                       &' per simulation. In case a larger number is needed, run the code several times'
        error=.True.
      End If
    End If
    ! print erro if any
    If (error) Then
      Call info(messages,2) 
      Call error_stop(' ')
    End If
    
    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_track_unchanged_chemistry" to close the block.&
                                  & Check if directives "tag" and "list_indexes" are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_track_unchanged_chemistry') Exit
      Call check_for_rubbish(iunit, '&track_unchanged_chemistry')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='list_indexes') Then
        traj_data%unchanged%indexes=-1
        Read (iunit, Fmt=*, iostat=io) traj_data%unchanged%list_indexes%type,&
                                       (traj_data%unchanged%indexes(j), j= 1, traj_data%unchanged%N0)
        Call set_read_status(word, io, traj_data%unchanged%list_indexes%fread, traj_data%unchanged%list_indexes%fail,&
                                     & traj_data%unchanged%list_indexes%type)

      Else If (Trim(word)=='tag') Then
         Read (iunit, Fmt=*, iostat=io) word, traj_data%unchanged%tag%type 
         Call set_read_status(word, io, traj_data%unchanged%tag%fread, traj_data%unchanged%tag%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with&
                                & "&end_track_unchanged_chemistry"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_track_unchanged_chemistry
  
  Subroutine read_settings_for_analysis(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the time settings from the &data_analysis block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(traj_type),  Intent(InOut) :: traj_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &data_analysis block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_data_analysis" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_data_analysis') Exit
      Call check_for_rubbish(iunit, '&data_analysis')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='time_interval') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%analysis%time_interval%tag, &
                                      & traj_data%analysis%time_interval%value,& 
                                      & traj_data%analysis%time_interval%units
        Call set_read_status(word, io, traj_data%analysis%time_interval%fread,&
                                      & traj_data%analysis%time_interval%fail)

      Else If (Trim(word)=='ignore_initial') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%analysis%ignore_initial%tag, &
                                      & traj_data%analysis%ignore_initial%value,& 
                                      & traj_data%analysis%ignore_initial%units
        Call set_read_status(word, io, traj_data%analysis%ignore_initial%fread,&
                                      & traj_data%analysis%ignore_initial%fail)

      Else If (Trim(word)=='overlap_time') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%analysis%overlap_time%tag, &
                                      & traj_data%analysis%overlap_time%value,& 
                                      & traj_data%analysis%overlap_time%units
        Call set_read_status(word, io, traj_data%analysis%overlap_time%fread,&
                                      & traj_data%analysis%overlap_time%fail)

      Else If (word(1:length) == 'normalise_at_t0') Then
        Read (iunit, Fmt=*, iostat=io) word, traj_data%analysis%normalise_at_t0%stat
       Call set_read_status(word, io, traj_data%analysis%normalise_at_t0%fread, traj_data%analysis%normalise_at_t0%fail)
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_data_analysis"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_settings_for_analysis
  
  Subroutine read_rdf(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for Radial Distribution Function (RDF)
    ! analysis from the &RDF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(traj_type),  Intent(InOut) :: traj_data 

    Integer(Kind=wi)   :: io, length, j
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &RDF block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_rdf" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_rdf') Exit
      Call check_for_rubbish(iunit, '&rdf')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='tags_species_a') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%rdf%tags_species_a%type, traj_data%rdf%num_type_a
        Call prevent_segmentation(iunit, io, traj_data%rdf%tags_species_a%type, traj_data%rdf%num_type_a,&
                                & 'max_components', max_components, set_error)
        traj_data%rdf%type_a=' '
        Read (iunit, Fmt=*, iostat=io) traj_data%rdf%tags_species_a%type, traj_data%rdf%num_type_a,&
                                       (traj_data%rdf%type_a(j), j= 1, traj_data%rdf%num_type_a)
        Call set_read_status(word, io, traj_data%rdf%tags_species_a%fread, traj_data%rdf%tags_species_a%fail,&
                           & traj_data%rdf%tags_species_a%type)

      Else If (Trim(word)=='tags_species_b') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%rdf%tags_species_b%type, traj_data%rdf%num_type_b
        Call prevent_segmentation(iunit, io, traj_data%rdf%tags_species_b%type, traj_data%rdf%num_type_b,&
                                & 'max_components', max_components, set_error)
        traj_data%rdf%type_b=' '
        Read (iunit, Fmt=*, iostat=io) traj_data%rdf%tags_species_b%type, traj_data%rdf%num_type_b,&
                                       (traj_data%rdf%type_b(j), j= 1, traj_data%rdf%num_type_b)
        Call set_read_status(word, io, traj_data%rdf%tags_species_b%fread, traj_data%rdf%tags_species_b%fail,&
                           & traj_data%rdf%tags_species_b%type)

      Else If (Trim(word)=='dr') Then
         Read (iunit, Fmt=*, iostat=io) traj_data%rdf%dr%tag, traj_data%rdf%dr%value, traj_data%rdf%dr%units 
         Call set_read_status(word, io, traj_data%rdf%dr%fread, traj_data%rdf%dr%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_rdf"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_rdf

  Subroutine read_coord_distrib(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for computing the coordinate distribution
    ! of selective species. Information must be provided in the 
    ! &coord_distrib block 
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(traj_type),  Intent(InOut) :: traj_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &coord_distrib block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly.&
                                  & Use "&end_coord_distrib" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_coord_distrib') Exit
      Call check_for_rubbish(iunit, '&coord_distrib')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='species') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%coord_distrib%species_dir%type, traj_data%coord_distrib%species
        Call set_read_status(word, io, traj_data%coord_distrib%species_dir%fread,&
                           & traj_data%coord_distrib%species_dir%fail,traj_data%coord_distrib%species_dir%type)

      Else If (Trim(word)=='delta') Then
         Read (iunit, Fmt=*, iostat=io) traj_data%coord_distrib%delta%tag, &
                                      & traj_data%coord_distrib%delta%value,&
                                      & traj_data%coord_distrib%delta%units 
         Call set_read_status(word, io, traj_data%coord_distrib%delta%fread, traj_data%coord_distrib%delta%fail)

      Else If (Trim(word)=='coordinate') Then
        Read (iunit, Fmt=*, iostat=io) word, traj_data%coord_distrib%coordinate%type
        Call set_read_status(word, io, traj_data%coord_distrib%coordinate%fread,& 
                           & traj_data%coord_distrib%coordinate%fail,&
                           & traj_data%coord_distrib%coordinate%type)
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_coord_distrib"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_coord_distrib
  
  Subroutine read_msd(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for mean square displacement (MSD)
    ! analysis from the &MSD block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(traj_type), Intent(InOut)  :: traj_data 

    Integer(Kind=wi)   :: io, length, j
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &MSD block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_msd" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_msd') Exit
      Call check_for_rubbish(iunit, '&msd')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='select') Then
        Read (iunit, Fmt=*, iostat=io) word, traj_data%msd%select%type
        Call set_read_status(word, io, traj_data%msd%select%fread, traj_data%msd%select%fail,&
                           & traj_data%msd%select%type)

      Else If (Trim(word)=='pbc_xyz') Then
         Read (iunit, Fmt=*, iostat=io) word, (traj_data%msd%pbc(j), j= 1, 3)
         Call set_read_status(word, io, traj_data%msd%pbc_xyz%fread, traj_data%msd%pbc_xyz%fail)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%msd%print_all_intervals%stat
       Call set_read_status(word, io, traj_data%msd%print_all_intervals%fread, traj_data%msd%print_all_intervals%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_msd"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_msd

  Subroutine read_orientational_chemistry(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &orientational_chemistry block
    !
    ! author    - i.scivetti February 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(traj_type), Intent(InOut)  :: traj_data 

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
        Read (iunit, Fmt=*, iostat=io) word, traj_data%chem_ocf%variable%type
        Call set_read_status(word, io, traj_data%chem_ocf%variable%fread,&
                                     & traj_data%chem_ocf%variable%fail,&
                                     & traj_data%chem_ocf%variable%type)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%chem_ocf%print_all_intervals%stat
       Call set_read_status(word, io, traj_data%chem_ocf%print_all_intervals%fread,&
                         & traj_data%chem_ocf%print_all_intervals%fail)
                            
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings. See the "use_code.md" file.&
                                & Have you properly closed the block with "&end_orientational_chemistry"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_orientational_chemistry
  
  
  Subroutine read_lifetime(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &lifetime block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(traj_type), Intent(InOut)  :: traj_data 

    Integer(Kind=wi)   :: io, length
    Character(Len=256) :: message, word
    Character(Len=256) :: set_error
    
    set_error = '***ERROR in the &lifetime block (SETTINGS file).'

    Do
      Read (iunit, Fmt=*, iostat=io) word
      If (io /= 0) Then
        Write (message,'(2(1x,a))') Trim(set_error), 'It appears the block has not been closed correctly. Use&
                                  & "&end_lifetime" to close the block.&
                                  & Check if directives are set correctly.'         
        Call error_stop(message) 
      End If  
      
      Call get_word_length(word,length)
      Call capital_to_lower_case(word)
      If (Trim(word)=='&end_lifetime') Exit
      Call check_for_rubbish(iunit, '&lifetime')

      If (word(1:1) == '#' .Or. word(1:3) == '   ') Then
      ! Do nothing if line is a comment of we have an empty line
      Read (iunit, Fmt=*, iostat=io) word

      Else If (Trim(word)=='method') Then
        Read (iunit, Fmt=*, iostat=io) word, traj_data%lifetime%method%type
        Call set_read_status(word, io, traj_data%lifetime%method%fread,&
                                     & traj_data%lifetime%method%fail,&
                                     & traj_data%lifetime%method%type)

      Else If (Trim(word)=='rattling_wait') Then
         Read (iunit, Fmt=*, iostat=io) traj_data%lifetime%rattling_wait%tag, &
                                      & traj_data%lifetime%rattling_wait%value,& 
                                      & traj_data%lifetime%rattling_wait%units
         Call set_read_status(word, io, traj_data%lifetime%rattling_wait%fread,&
                                      & traj_data%lifetime%rattling_wait%fail)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%lifetime%print_all_intervals%stat
       Call set_read_status(word, io, traj_data%lifetime%print_all_intervals%fread, traj_data%lifetime%print_all_intervals%fail)
                            
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_lifetime"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_lifetime
!   
  Subroutine read_ocf(iunit, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information for orientational correlation function (OCF)
    ! analysis from the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(traj_type), Intent(InOut)  :: traj_data 

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
        Read (iunit, Fmt=*, iostat=io) word, traj_data%ocf%u_definition%type
        Call set_read_status(word, io, traj_data%ocf%u_definition%fread, traj_data%ocf%u_definition%fail,&
                           & traj_data%ocf%u_definition%type)

      Else If (Trim(word)=='legendre_order') Then
         Read (iunit, Fmt=*, iostat=io) word, traj_data%ocf%legendre_order%value
         Call set_read_status(word, io, traj_data%ocf%legendre_order%fread,&
                            & traj_data%ocf%legendre_order%fail)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, traj_data%ocf%print_all_intervals%stat
       Call set_read_status(word, io, traj_data%ocf%print_all_intervals%fread, traj_data%ocf%print_all_intervals%fail)
                            
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_ocf"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_ocf
  
  Subroutine prevent_segmentation(iunit, io, in_name, input, ref_name, reference, error)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to prevent segmentation fault in case the user wants to define
    ! settings beyonf the reference number. 
    !
    ! author    - i.scivetti April 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Integer(Kind=wi), Intent(In   ) :: io
    Character(Len=*), Intent(In   ) :: in_name
    Integer(Kind=wi), Intent(In   ) :: input
    Character(Len=*), Intent(In   ) :: ref_name
    Integer(Kind=wi), Intent(In   ) :: reference
    Character(Len=*), Intent(In   ) :: error
    
    Character(Len=256) :: messages(3)
    Character(Len=256) :: word, default

    Write (messages(1),'(a,i3,a)')  Trim(error)
    If (io == 0) Then
      If (input>reference) Then
        Write(word,*)    input
        Write(default,*) reference
        Write (messages(2),'(a)') 'Are you sure you want to consider '//Trim(Adjustl(word))//' components for&
                                 & "'//Trim(in_name)//'"? The maximum default value is '//Trim(Adjustl(default))
        Write (messages(3),'(a)') 'If you are sure of what you are doing, look for the parameter "'//Trim(ref_name)//&
                                  &'" in the code, increase its value as needed and recompile.'
        Call info(messages, 3)
        Call error_stop(' ')
      Else If (input < 1) Then
        Write (messages(2),'(a)') ' The number associated with "'//Trim(in_name)//'" must be positive!'
        Call info(messages, 2)
        Call error_stop(' ')
      Else
        Backspace iunit
      End If
    Else
      Write (messages(2),'(a)') 'Problems in the settings of "'//Trim(in_name)//'". Please check.'
      Call info(messages, 2)
      Call error_stop(' ')
    End If
  
  End Subroutine prevent_segmentation
  
  Subroutine read_input_composition(iunit, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the amount of atoms and chemical elements for each atomic
    ! species. Information must be defined in the &input_composition block
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(model_type),  Intent(InOut) :: model_data
 
    Logical  :: endblock, loop, error
    Logical  :: header(3), error_duplication
 
    Integer(Kind=wi)   :: io, i, j, k, ilist, ic
    Character(Len=256) :: messages(10), word
    Character(Len=64 ) :: error_input_composition
 
    error_input_composition = '***ERROR in &input_composition of SETTINGS file'
    Write (messages(1),'(a)') error_input_composition
 
    header=.False.
    error_duplication=.False.
    error=.False.
    ilist=10
 
    Write (messages(3),'(1x,a)')    'The correct structure for the block must be:'
    Write (messages(4),'(1x,a)')    '&input_composition'
    Write (messages(5),'(1x,a)')    '  atomic_species    Nsp'
    Write (messages(6),'(1x,a)')    '  tags      tg1    tg2    tg3   .... tgNsp'
    Write (messages(7),'(1x,a)')    '  elements  E_tg1  E_tg2  E_tg3 .... E_tgNsp'
    Write (messages(8),'(1x,a)')    '  amounts   N_tg1  N_tg2  N_tg3 .... N_tgNsp'
    Write (messages(9),'(1x,a)')    '&end_input_composition'
    Write (messages(10),'(1x,a)')    'See the "use_code.md" file for details'
    
    ! Read number of extra bonds
    Read (iunit, Fmt=*, iostat=io) word, model_data%input_composition%atomic_species
    Call capital_to_lower_case(word) 
    If (Trim(word) /= 'atomic_species') Then
      Write (messages(2),'(3a)') 'Directive "', Trim(word), &
                         & '" has been found, but directive "atomic_species" is expected.'
      error=.True.
    End If 
    If (io /= 0) Then
      Write (messages(2),'(a)') 'Wrong (or missing) specification for directive "atomic_species"'
      error=.True.
    Else
      If (model_data%input_composition%atomic_species<1) Then
        Write (messages(2),'(a)') 'The "atomic_species" directive MUST BE >= 1'
        error=.True.
      End If  
    End If
   
    If (error) Then
      Call info(messages, ilist) 
      Call error_stop(' ')
    End If

    Call model_data%init_input_composition()
   
    i=1
    Do While (i <= 3)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(iunit, '&input_composition')
      If (word(1:1)/='#') Then
        If (word(1:1)/='&') Then
          Call capital_to_lower_case(word) 
          If (Trim(word)=='tags') Then
            If (.Not. header(1)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&input_composition') 
              Read (iunit, Fmt=*, iostat=io) word,&
                            & (model_data%input_composition%tag(j), j = 1, model_data%input_composition%atomic_species)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read tags for atoms'
                Call info(messages, 2)
                Call error_stop(' ')
              End If
              header(1)=.True.
            Else
              error_duplication=.True.
            End If
          Else If (Trim(word)=='amounts') Then
            If (.Not. header(2)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&input_composition') 
              Read (iunit, Fmt=*, iostat=io) word,&
                               & (model_data%input_composition%N0(j), j = 1, model_data%input_composition%atomic_species)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the "amount" directive for each atomic tag'
                Call info(messages, ilist)
                Call error_stop(' ')
              End If  
              header(2)=.True.
            Else
              error_duplication=.True.
            End If  
          Else If (Trim(word)=='elements') Then
            If (.Not. header(3)) Then 
              i=i+1
              Call check_for_rubbish(iunit, '&input_composition') 
              Read (iunit, Fmt=*, iostat=io) word, &
                         & (model_data%input_composition%element(j), j = 1, model_data%input_composition%atomic_species)
              If (io /= 0) Then
                Write (messages(2),'(1x,a)') 'Problems to read the chemical element for each atomic tag (directive "element")'
                Call info(messages, ilist)
                Call error_stop(' ')
              End If  
              header(3)=.True.        
            Else
              error_duplication=.True.
            End If
          Else
            Write (messages(2),'(1x,3a)') 'Wrong descriptor "', Trim(word), '". Please chek the amount of atomic species&
                                          & defined and the input species.'
            Call info(messages, ilist)
            Call error_stop(' ')
          End If
        Else
          Write (messages(2),'(1x,a)')    ' '
          Call info(messages, ilist)
          Call error_stop(' ')
        End If
      End If
 
      If (error_duplication) Then
        Write (messages(2),'(1x,3a)') 'Descriptor "', Trim(word), '" is duplicated within the block.&
                                     & If there is no duplication, then there is an inconsistency&
                                     & in the info provided for value of "atomic_species"'
        Call info(messages, ilist)
        Call error_stop(' ')
      End If
 
    End Do 
 
    endblock=.False.
    Do While (.Not. endblock)
      Read (iunit, Fmt=*, iostat=io) word
      Call check_end(io, '&input_composition') 
      Call capital_to_lower_case(word)
      If (word /= '&end_input_composition') Then
        If (word(1:1) /= '#') Then
            Write (messages(2),'(3a)') 'All info have already been defined. Directive "',&
                                    & Trim(word), '" is not valid. Block must be&
                                    & closed with sentence &end_input_composition.' 
            Call info(messages, ilist)
            Call error_stop(' ')
        End If
      Else
        endblock=.True.
      End If
    End Do

    ! Check if species tags contain asterix
    Do i=1, model_data%input_composition%atomic_species
      ic= Index(Trim(model_data%input_composition%tag(i)), '*') 
      If (ic > 0) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%input_composition%tag(i)), &
                                 '" contains an asterisk. Defined species MUST NOT contain asterisks. Please correct.'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do
    
    ! Check if the number of atoms are correct
    Do i=1, model_data%input_composition%atomic_species
      If (model_data%input_composition%N0(i)< 0) Then
        Write (messages(2),'(3a)') 'Tag "', Trim(model_data%input_composition%tag(i)), '" CANNOT be associated with&
                                 & negative number of atoms within the input structure! Please correct'
        Call info(messages,2)
        Call error_stop(' ')
      End If
    End Do
 
    ! Calculate the number of total atoms set in the block
    model_data%input_composition%numtot=0
    Do i=1, model_data%input_composition%atomic_species
      model_data%input_composition%numtot=model_data%input_composition%numtot+model_data%input_composition%N0(i)
    End Do
 
    ! Assing atomic numbers
    Do i=1, model_data%input_composition%atomic_species
      loop=.True.
      k=1
      Do While (k <= NPTE .And. loop)
        If (Trim(chemsymbol(k))==Trim(model_data%input_composition%element(i))) Then
          loop=.False.
        End If
        k=k+1
      End Do
      If (loop) Then
         Write (messages(2),'(1x,5a)') 'Wrong chemical element "', Trim(model_data%input_composition%element(i)),&
                                      & '" defined for species tag "', Trim(model_data%input_composition%tag(i)),&
                                      & '". Please check.' 
         Call info(messages, 2)
         Call error_stop(' ')
      End If
    End Do
 
  End Subroutine read_input_composition

  Subroutine check_settings(files, model_data, traj_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of the defined directive
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(model_type),     Intent(InOut) :: model_data 
    Type(traj_type),      Intent(InOut) :: traj_data
 
    Call check_model_settings(files, model_data)
    Call check_trajectory_settings(files, traj_data, model_data)    
    ! Print model related settings
    Call print_model_settings(files, model_data)
    Call refresh_out(files)

  End Subroutine check_settings

End module settings

