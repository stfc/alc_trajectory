!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to analyse the trajectory. If the directive change_chemistry
! is set to .True., the algorithm searches and tracks changes of 
! chemical species based on the information of the &search_chemistry
! block. This modelue executes all the implemented options for 
! RDF, OCF, MSD and TCF analyses.
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:        i.scivetti  Feb 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module trajectory 
  Use atomic_model, Only : model_type, &
                           about_cell, &
                           atomistic_model, & 
                           check_definition_bonds, &
                           check_PBC,&
                           check_cell_consistency,&
                           check_length_directive, &
                           check_orthorhombic_cell,&
                           read_model,&
                           obtain_maximum_number_species,&
                           identify_monitored_indexes

  Use constants,    Only : max_components, &
                           max_at_species, &
                           max_unchanged_atoms,&
                           pi

  Use fileset,      Only : file_type, &
                           FILE_TCF, &
                           FILE_TCF_AVG, &
                           FILE_MSD, &
                           FILE_MSD_AVG, &
                           FILE_OCF, &
                           FILE_OCF_AVG, &
                           FILE_RDF, &
                           FILE_RES_TIMES, &
                           FILE_SET, & 
                           FILE_TRAJECTORY,&
                           FILE_TRACK_CHEMISTRY,&
                           FILE_TAGGED_TRAJ, &
                           FILE_UNCHANGED_CHEM, & 
                           refresh_out

  Use input_types,  Only : in_integer, &
                           in_logic,   &
                           in_scalar,  &
                           in_param,   & 
                           in_param_array, & 
                           in_string

  Use numprec,      Only : wi,&
                           wp 
  
  Use process_data, Only : capital_to_lower_case, &
                           detect_rubbish,        &
                           remove_symbols,        &
                           remove_front_tabs
  Use unit_output,  Only : error_stop, &
                           info

  Implicit None
  Private

  Type :: atom_type
     Real(Kind=wp)    :: r(3)
     Integer(Kind=wi) :: indx
     Character(Len=8) :: tag
     Character(Len=2) :: element
  End Type 

  Type :: box_type
     Real(Kind=wp)    :: cell(3,3)
     Real(Kind=wp)    :: invcell(3,3)
     Real(Kind=wp)    :: volume
     Real(Kind=wp)    :: cell_length(3)
  End Type 

  Type :: ocf_type
     Type(in_string)   :: invoke
     Type(in_integer)  :: legendre_order
     Type(in_string)   :: u_definition
  End Type 
  
  Type :: analysis_type
    Type(in_string)   :: invoke
    Type(in_param)    :: time_interval
    Type(in_param)    :: ignore_initial
    Type(in_param)    :: overlap_time
    Integer(Kind=wi)  :: N_seg
    Integer(Kind=wi)  :: Ninterval
    Integer(Kind=wi)  :: frame_ini
    Integer(Kind=wi), Allocatable :: seg_indx(:,:)
    Real(Kind=wp),    Allocatable :: variable(:,:)
    Integer(Kind=wi), Allocatable :: max_points(:) 
  End Type
  
  ! Type to describe species
  Type :: species_type
    Integer(Kind=wi) :: list(max_at_species)
    Logical          :: alive
    Real(Kind=wp)    :: u(3,2)
    Real(Kind=wp)    :: u0(3,2)
  End Type
  
  !Type to describe the region where to constrain the analysis
  Type :: region_type
    Type(in_string)  :: define
    Logical          :: belong(3,max_components)
    Type(in_string)  :: invoke(3,max_components)
    Character(Len=8) :: inout(3,max_components)
    Logical          :: inside(3,max_components)
    Real(Kind=wp)    :: domain(3,2,max_components)
    Integer(Kind=wi) :: number(3)
  End Type
  
  !Type to describe the region where to constrain the analysis
  Type :: msd_type
    Type(in_string)  :: invoke
    Type(in_string)  :: pbc_xyz
    Type(in_string)  :: select
    Logical          :: pbc(3)
    Real(Kind=wp)    :: r2
  End Type

  !Type to describe the region where to constrain the analysis
  Type :: rdf_type
    Type(in_string)  :: invoke
    Type(in_string)  :: tags_species_a
    Type(in_string)  :: tags_species_b
    Type(in_param)   :: dr
    Real(Kind=wp)    :: rmax
    Character(Len=8) :: type_a(max_components)
    Character(Len=8) :: type_b(max_components)
    Integer(Kind=wi) :: num_type_a
    Integer(Kind=wi) :: num_type_b
  End Type

  !Type to print the position of selected atoms, whose content remain unchanged
  Type :: unchanged_type
    Type(in_string)  :: invoke
    Type(in_string)  :: tag
    Integer(Kind=wi) :: N0
    Integer(Kind=wi) :: indexes(max_unchanged_atoms)
    Type(in_string)  :: list_indexes
  End Type

  !Type to describe the region where to constrain the analysis
  Type :: track_type
    Type(atom_type),    Allocatable :: config(:,:)
  End Type

  !Type to describe the region where to constrain the analysis
  Type :: life_type
    Type(in_string)   :: invoke
    Type(in_string)   :: method
    Type(in_param)    :: rattling_wait
  End Type
  
  ! Type for eqcm data and analysis
  Type, Public :: traj_type
    Private
    Type(species_type),   Allocatable :: species(:,:)
    Type(atom_type),      Allocatable :: config(:,:)
    Type(box_type),       Allocatable :: box(:)
    Type(region_type),    Public      :: region
    Type(rdf_type),       Public      :: rdf
    Type(ocf_type),       Public      :: ocf
    Type(msd_type),       Public      :: msd
    Type(life_type),      Public      :: lifetime
    Type(unchanged_type), Public      :: unchanged
    Type(analysis_type),  Public      :: analysis
    Type(track_type)                  :: track_chem              
    Type(in_param),       Public      :: timestep
    Type(in_logic),       Public      :: print_retagged_trajectory
    Integer(Kind=wi)                  :: frames
    Integer(Kind=wi)                  :: Nmax_species
    Integer(Kind=wi)                  :: N_species
    Logical                           :: reload_trajectory
    Type(in_string),      Public      :: ensemble
  Contains
    Private
      Procedure         :: alloc_trajectory  => allocate_trajectory_arrays
      Procedure         :: alloc_analysis    => allocate_analysis_arrays
      Final             :: cleanup
  End Type traj_type

  Public :: trajectory_analysis, research_trajectory  
  Public :: check_trajectory_settings
  

Contains

  Subroutine allocate_trajectory_arrays(T, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate trajectory arrays
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(traj_type),   Intent(InOut)  :: T
    Type(model_type),   Intent(In   )  :: model_data
    
    Integer(Kind=wi)    :: fail(2)
    Character(Len=256)  :: message
    Logical             :: error_alloc

    error_alloc=.False.
    fail=0
    
    Write (message,'(1x,1a)') '***ERROR: Allocation problems for the trajectory&
                                & (subroutine allocate_trajectory_arrays). It is likely that the&
                                & trajectory and/or the system is exceedingly large.'

    Allocate(T%config(T%frames, model_data%config%num_atoms), Stat=fail(1))
    Allocate(T%box(T%frames),                                 Stat=fail(2))
    If (Any(fail > 0)) Then
       error_alloc=.True.
    End If

    If(model_data%change_chemistry%stat) Then 
      Allocate(T%track_chem%config(T%frames, max_components), Stat=fail(1))
      If (Any(fail > 0)) Then
        error_alloc=.True.
      End If
    End If

    If (model_data%config%monitored_species%fread) Then
      Allocate(T%species(T%frames, model_data%config%Nmax_species), Stat=fail(1))
      If (Any(fail > 0)) Then
        error_alloc=.True.
      End If
      T%Nmax_species=model_data%config%Nmax_species
    End If

    If (error_alloc) Then
      Call error_stop(message)
    End If
      
  End Subroutine allocate_trajectory_arrays

  Subroutine allocate_analysis_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate arrays for date analysis
    !
    ! author    - i.scivetti Mrch 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(traj_type),   Intent(InOut)  :: T
    
    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message

    fail=0
    Allocate(T%analysis%seg_indx(2,T%analysis%N_seg),                     Stat=fail(1))
    Allocate(T%analysis%variable(T%analysis%Ninterval, T%analysis%N_seg), Stat=fail(2))
    Allocate(T%analysis%max_points(T%analysis%N_seg), Stat=fail(3))

    If (Any(fail > 0)) Then
       Write (message,'(1x,1a)') '***ERROR: Allocation problems for trajectory analysis &
                               & (subroutine allocate_analysis_arrays). Please review&
                               & settings of the &data_analysis block'
       Call info(message, 1)
       Call error_stop(' ')
    End If

  End Subroutine allocate_analysis_arrays
  
  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type) :: T

    If (Allocated(T%config)) Then
      Deallocate(T%config)
    End If 

    If (Allocated(T%track_chem%config)) Then
      Deallocate(T%track_chem%config)
    End If 

    If (Allocated(T%analysis%seg_indx)) Then
      Deallocate(T%analysis%seg_indx)
    End If 

    If (Allocated(T%analysis%variable)) Then
      Deallocate(T%analysis%variable)
    End If 

    If (Allocated(T%analysis%max_points)) Then
      Deallocate(T%analysis%max_points)
    End If 
    
  End Subroutine cleanup
    
  Subroutine research_trajectory(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to obtain the trajectory. Atomic positions and elements must
    ! be provided in the TRAJECTORY file. If there are changes in the chemistry,
    ! the initial atomic tags asociated to each atoms are redefined according to
    ! the settings of &seach_chemistry (and &possible_extra_bonds, if defined) 
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    
    Logical            :: safe, loop_traj, fortho
    Character(Len=256) :: message
    Character(Len=32 ) :: input_file, set_error
    Integer(Kind=wi)   :: i, j
    
    input_file=Trim(files(FILE_TRAJECTORY)%filename)
    set_error = '***ERROR -'

    Inquire(File=input_file, Exist=safe)
    If (.not.safe) Then
      Call info(' ', 1)
      Write (message,'(4(1x,a))') Trim(set_error), 'File', Trim(input_file), 'not found'
      Call error_stop(message)
    Else
      traj_data%reload_trajectory=.False. 
    End If

    ! Further checking to make sure it is all set
    Call cross_checking(files, model_data, traj_data)     

    ! Initialise flags 
    model_data%config%allocated_model_geo=.False.
    
    ! Open the TRAJECTORY file
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=Trim(input_file),Status='old')
    Call read_model(files, model_data, 1, traj_data%ensemble%type)
    Close(files(FILE_TRAJECTORY)%unit_no) 
    ! Scale simulation cell 
    model_data%config%cell=model_data%config%cell_scaling * model_data%config%cell
    If (Trim(model_data%config%position_units%type) == 'bohr') Then  
      Do j=1,3 
        model_data%config%atom(:)%r(j)=model_data%config%position_scaling* model_data%config%atom(:)%r(j) 
      End Do
    End If
    ! Compute cell related quantities for first checks 
    Call about_cell(model_data%config%cell,model_data%config%invcell,&
                  & model_data%config%cell_length, model_data%config%volume)
                    
    ! Check consistency between the system and the input cell
    Call check_cell_consistency(model_data)
    Call check_orthorhombic_cell(model_data%config%cell, fortho) 
      If (.Not.fortho) Then
        Call info(' ', 1)
        Write (message,'(1x,1a)') '***WARNING: the atomic model is not orthorhombic. This code has only been tested&
                                  & for orthorhombic cells. We do not guarantee a correct functioning.'
        Call info(message, 1)
     End If        
    
    If (model_data%config%monitored_species%fread) Then
      ! Compute the maximum amount of species to be monitored and allocate them      
      Call obtain_maximum_number_species(model_data)
      Call model_data%init_species() 
      Call identify_monitored_indexes(model_data)
      model_data%config%species(:)%alive=.False.
    End If

    !Check the labelling against info of the &track_unchanged_chemistry block
    If (traj_data%unchanged%invoke%fread) Then
      Call check_initial_unchanged_labels(files, traj_data, model_data)
    End If
    
    ! Check bonds againts the size of the simulation cell
    If(model_data%change_chemistry%stat) Then 
      Call check_definition_bonds(model_data, 1)
    End If
    
    ! Check if the defined region for analysis is within the simulation cell
    If (traj_data%region%define%fread) Then
      Call check_region_domain(model_data, traj_data, 1)
    End If 
    
    ! Search for the number of frames
    Call obtain_number_frames(files, model_data, traj_data)
    ! Allocate trajectory arrays
    Call traj_data%alloc_trajectory(model_data)
    ! Check time settings (from &data_analysis) against the trajectory 
    Call check_time_settings(files, traj_data)
    ! Print trajectory related settings
    Call print_trajectory_settings(traj_data, model_data)
    ! Refresh 
    Call refresh_out(files)
    
    ! Open the TRAJECTORY file
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=Trim(input_file),Status='old')
    Call info(' ', 1)
    Call info('Start of the analysis', 1)
    Call info('=====================', 1)
    Write (message,'(1x,a,i6,a)') 'The code has identified a total of ', traj_data%frames, ' frames'
    Call info(message, 1)
    Call info(' Reading trajectory from the "'//Trim(input_file)//'" file...', 1)
    i=1
    loop_traj=.True.
    Do While (i <= traj_data%frames)
      Call read_model(files, model_data, i, traj_data%ensemble%type)
      If (Trim(model_data%config%position_units%type) == 'bohr') Then
        Do j=1,3 
          model_data%config%atom(:)%r(j)=model_data%config%position_scaling* model_data%config%atom(:)%r(j) 
        End Do
      End If
      If (Trim(traj_data%ensemble%type) == 'npt') Then
        model_data%config%cell=model_data%config%cell_scaling * model_data%config%cell
        Call about_cell(model_data%config%cell,model_data%config%invcell,&
                        model_data%config%cell_length, model_data%config%volume)
        Call check_definition_bonds(model_data, i)
      End If
      ! Identify the components of the model
      Call atomistic_model(model_data, i)
      ! Copy to trajectory arrays for later analysis
      Call copy_to_trajectory(traj_data, model_data, i)
      i=i+1
    End Do

    Close(files(FILE_TRAJECTORY)%unit_no) 
    Call info(' INFO: The trajectory has been defined successfully', 1)
    Call info(' ', 1)
    Call refresh_out(files)

    If (traj_data%frames==1) Then
      Call info(' **********************************************************', 1)
      Call info(' ** WARNING: ONLY ONE FRAME WAS DETECTED IN THE TRAJECTORY!', 1)
      Call info(' **********************************************************', 1)
      Call info(' ', 1)
    End If

    If(model_data%change_chemistry%stat) Then
     Call residence_percentage(traj_data, model_data)
     Call refresh_out(files) 
    End If
    
    If(traj_data%print_retagged_trajectory%stat) Then 
      input_file=Trim(files(FILE_TAGGED_TRAJ)%filename)
      Call print_tagged_trajectory(files, model_data, traj_data)
      Write (message,'(1x,a)') 'A copy of the trajectory with modified tags for the atomic species was printed&
                              & to the "'//Trim(input_file)//'" file'
      Call info(message, 1)
      Call refresh_out(files)
    End If 
  
  End Subroutine research_trajectory
  
  Subroutine print_tagged_trajectory(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to analyse the trajectory
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
  
    Integer(Kind=wi)   :: iunit, i, l, k 
  
    ! Print tracked species
      If(model_data%change_chemistry%stat) Then
        Open(Newunit=files(FILE_TAGGED_TRAJ)%unit_no, File=files(FILE_TAGGED_TRAJ)%filename, Status='Replace')
        iunit=files(FILE_TAGGED_TRAJ)%unit_no
        Do l = 1, traj_data%frames
          Write(iunit,*) model_data%config%num_atoms, traj_data%frames
          Write(iunit,*) ' ' 
          Do i = 1, model_data%config%num_atoms 
            Write(iunit,'(a, 4x, 3(f11.3), 4x, a)') Trim(traj_data%config(l,i)%element),&
                                                 & (traj_data%config(l,i)%r(k), k=1, 3),&
                                                 &  Trim(traj_data%config(l,i)%tag) 
          End Do               
        End Do
        Close(iunit)
      End If
  
  End Subroutine
  
  Subroutine trajectory_analysis(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to analyse the trajectory depending on the options of the
    ! SETTINGS file
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
  
    Character(Len=256)  :: message

    If (model_data%config%monitored_species%fread .And. model_data%species_definition%compute_amount%stat) Then
      Call compute_number_monitored_species(traj_data, model_data)
    End If
    
    If(model_data%change_chemistry%stat) Then 
      Call print_tracking_species(files, traj_data, model_data)
      If (traj_data%lifetime%invoke%fread) Then 
        Call transfer_correlation_function(files, traj_data, model_data)
        Call residence_times(files, traj_data, model_data)
      End If
    End If

    If (traj_data%ocf%invoke%fread) Then
      Call orientational_correlation_function(files, traj_data)
    End If

    If (traj_data%msd%invoke%fread) Then
      Call mean_squared_displacement(files, traj_data, model_data)
    End If

    If (traj_data%rdf%invoke%fread) Then
      Call radial_distribution_function(files, traj_data, model_data)
    End If

    If (traj_data%unchanged%invoke%fread) Then
      Call print_unchanged_chemistry(files, traj_data)
    End If

    If (Trim(traj_data%ensemble%type)/='nve') Then
      If (traj_data%lifetime%invoke%fread .Or.&
          traj_data%ocf%invoke%fread      .Or.&
          traj_data%msd%invoke%fread) Then   
        Call info(' ', 1)
        Call info(' ****************************************************************************************', 1)
        Write (message,'(1x,a)') 'IMPORTANT: The user should bear in mind that the computed properties&
                                 & might be influenced'
        Call info(message, 1)                         
        If (Trim(traj_data%ensemble%type)=='nvt') Then
          Write (message,'(12x,a)') 'by the "thermostat" used to generate the trajectory'  
        Else If (Trim(traj_data%ensemble%type)=='npt') Then
          Write (message,'(12x,a)') 'by the "thermostat" and "barostat" used to generate the trajectory'  
        End If
        Call info(message, 1)
        Call info(' ****************************************************************************************', 1)
      End If
    End If
    
  End Subroutine trajectory_analysis

  Subroutine residence_percentage(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the residence percentage of changing chemistry 
    ! species along the trajectory
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, j, k, l, m
    Integer(Kind=wi)   :: counts(model_data%chem%acceptor%N0_incl)
    Real(Kind=wp)      :: amount
    Character(Len=8)   :: word
    Character(Len=256) :: messages(5)
    
    l=0
    counts=0
    Do i = traj_data%analysis%frame_ini, traj_data%frames
      l=l+1
      k=0
      Do j = 1, model_data%chem%N0%value
        Do m = 1, model_data%chem%acceptor%N0_incl
          word=Trim(model_data%chem%acceptor%tg_incl(m))//'*'
          If (Trim(traj_data%track_chem%config(i,j)%tag)==Trim(word)) Then
            counts(m)=counts(m)+1
            k=k+1 
          End If
        End Do
      End Do
    End Do

    Write (messages(1),'(1x,a)') 'Population probabilities of formed chemical species'
    Write (messages(2),'(1x,a)') '(from the species defined in the "include_tags" directive in "&enviroment_criteria")'
    Write (messages(3),'(1x,a)') '-----------------------'
    Write (messages(4),'(1x,a)') 'Fraction (%)    Species'
    Write (messages(5),'(1x,a)') '-----------------------'
    Call info(messages, 5)                            
    Do m = 1, model_data%chem%acceptor%N0_incl
      word=Trim(model_data%chem%acceptor%tg_incl(m))//'*' 
      amount= 100.0_wp * Real(counts(m),Kind=wp)/(l*model_data%chem%N0%value)
      Write (messages(1),'(6x,f7.3,4x,a)') amount, Trim(word)
      Call info(messages, 1)                                                                    
    End Do                    
    Write (messages(1),'(1x,a)') '-----------------------'
    Call info(messages, 1)
    If (model_data%chem%acceptor%N0_incl==1) Then
      Write (messages(1),'(1x,a)') 'NOTE: 100% population is consistent with having defined a single species in "include_tags"'
      Call info(messages, 1)
    End If
    Call info(' ', 1)
    
  End Subroutine residence_percentage
  
  Subroutine print_tracking_species(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the positions of those species that change their 
    ! chemsitry along the trajectory
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(InOut) :: model_data
  
    Integer(Kind=wi)   :: iunit, i, l 
    Character(Len=256) :: num_species
    Character(Len=256) :: message
  
    ! Print tracked species
    Open(Newunit=files(FILE_TRACK_CHEMISTRY)%unit_no, File=files(FILE_TRACK_CHEMISTRY)%filename, Status='Replace')
    iunit=files(FILE_TRACK_CHEMISTRY)%unit_no
    Write(num_species,*) model_data%chem%N0%value 
    Write (iunit,'(a,9x,2a)') '# Time (ps)', 'XYZ_Species_1 .... XYZ_Species_', Adjustl(Trim(num_species)) 
    Do i = traj_data%analysis%frame_ini, traj_data%frames
       Write(iunit,'(f10.4, 4x, *(f11.3))') (i-traj_data%analysis%frame_ini)*traj_data%timestep%value/1000_wp,&
                                       & (traj_data%track_chem%config(i,l)%r(:), l=1, model_data%chem%N0%value)
    End Do
    Write (message,'(1x,a)') 'The tracking of the changing chemical species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_TRACK_CHEMISTRY)%filename)//'" file'
    Call info(message, 1)
    Call refresh_out(files)
    Close(iunit)
  
  End Subroutine print_tracking_species
  
  Subroutine print_unchanged_chemistry(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the positions of those atomic indexes defined in the
    ! &track_unchanged_chemistry block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
  
    Integer(Kind=wi)   :: iunit, i, k, l 
    Character(Len=256) :: num_species, frame
    Character(Len=256) :: message, messages(3), spec
    Logical            :: flag
  
    flag=.True. 
  
    Open(Newunit=files(FILE_UNCHANGED_CHEM)%unit_no, File=files(FILE_UNCHANGED_CHEM)%filename, Status='Replace')
    iunit=files(FILE_UNCHANGED_CHEM)%unit_no
    spec=Trim(traj_data%unchanged%tag%type)
    If(traj_data%unchanged%N0==1) Then
      Write (iunit,'(a,7x,a)') '#Time (ps)', 'XYZ_'//Trim(spec)//'_1'
    Else
      Write(num_species,*) traj_data%unchanged%N0 
      Write (iunit,'(a,7x,a)') '#Time (ps)', 'XYZ_'//Trim(spec)//'_1 .... XYZ_'//Trim(spec)//&
                              &'_'// Adjustl(Trim(num_species))
    End If 
    
    i=1
    Do While (i <= traj_data%frames .And. flag)
      l =1
      Do While (l<= traj_data%unchanged%N0 .And. flag)
        k=traj_data%unchanged%indexes(l)
        If (Trim(traj_data%config(i,k)%tag) /= Trim(traj_data%unchanged%tag%type)) Then
          flag=.False.
        End If
        l=l+1
      End Do
      If (flag) Then
        Write(iunit,'(f10.4, 1x, *(f11.3))') (i-1)*traj_data%timestep%value/1000.0_wp,&
                & (traj_data%config(i,traj_data%unchanged%indexes(l))%r(:), l=1, traj_data%unchanged%N0)
      Else
        Write (messages(1),'(1x,a)') '**********************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: Problems with tracking species defined in the&
                                        & &track_unchanged_chemistry block'
        Write(num_species,*) k
        Write(frame,*)       i 
        Write (messages(2),'(1x,a)') '   Requested index "'//Trim(Adjustl(num_species))//'" has changed&
                                      & its chemistry at frame: '//Trim(Adjustl(frame))
        Write (messages(3),'(1x,a)') '   Please review the settings. The tracking was printed to the "'&
                                    &//Trim(files(FILE_UNCHANGED_CHEM)%filename)//'" file just up to this frame'
        Call info(messages, 3)
        Write (messages(1),'(1x,a)') '**********************************************'
        Call info(messages, 1)
      End If
      i=i+1
    End Do
    
    If (flag) Then
      Write (message,'(1x,a)') 'The tracking of unchanged chemical species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_UNCHANGED_CHEM)%filename)//'" file'
      Call info(message, 1)
    End If
    Call refresh_out(files)
    Close(iunit)

  End Subroutine print_unchanged_chemistry

  Subroutine compute_number_monitored_species(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the average number (and STD) of the monitored species 
    !
    ! author    - i.scivetti June 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi)  :: i, j
    Integer(Kind=wi)  :: num_at_a, net_frames
    Integer(Kind=wi)  :: accum_a
    
    Character(Len=256) :: messages(3)
    Logical            :: flag

    Real(Kind=wp) :: average, std, sum_i
    
    ! counting
    Real(Kind=wp), Allocatable  :: nat(:)
   
    ! In case &region is defined
    flag=.True.
    net_frames=0
    accum_a=0
    
    Allocate(nat(traj_data%frames-traj_data%analysis%frame_ini+1))
    
    ! Compute the histogram for atoms of type a and b
    Do i = traj_data%analysis%frame_ini, traj_data%frames
      ! Define the number and list of indexes for type of species "a"
      num_at_a=0
      net_frames=net_frames+1 
      Do j = 1, model_data%config%num_atoms
        If (Trim(model_data%species_definition%reference_tag%type)==Trim(traj_data%config(i,j)%tag)) Then
          If (traj_data%region%define%fread) Then
             Call within_region(traj_data, i, j, flag)
          End If
          If (flag) Then
            num_at_a=num_at_a+1
          End If
        End If
      End Do
      ! Accummulators
      nat(net_frames)=num_at_a 
      accum_a=accum_a+num_at_a
    End Do
      
    average= Real(accum_a,Kind=wp)/net_frames
    
    ! Compute average
    If (net_frames > 1) Then
      sum_i=0.0_wp
      j=0 
      Do i = traj_data%analysis%frame_ini, traj_data%frames
        j=j+1 
        sum_i=sum_i+(Real(nat(j),Kind=wp)-average)**2
      End Do
      std=sqrt(sum_i/(net_frames-1))
    Else
      std=0.0_wp
    End If

    If (traj_data%region%define%fread) Then
       Write (messages(1),'(1x,a)') 'Amount of monitored species "'//Trim(model_data%species_definition%name%type)//&
                                   &'" within the selected region as specified in the &region block'
    Else
       Write (messages(1),'(1x,a)') 'Amount of monitored species "'//Trim(model_data%species_definition%name%type)//&
                                   &'" within the simulation cell'
    End If
    Write (messages(2),'(1x,f8.2,5x,a,f8.2)')  average, '+/-', STD
    Call info(messages, 2)
    Call info(' ', 1)
    
  End Subroutine compute_number_monitored_species
  
  Subroutine radial_distribution_function(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the Mean Squared Displacement (MSD) based on the
    ! settings of the &MSD block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(InOut) :: model_data

    Integer(Kind=wi)  :: i, j, k, m, iunit, indx_a, indx_b
    Integer(Kind=wi)  :: num_at_a, num_at_b, nbins, net_frames
    Integer(Kind=wi)  :: accum_a, accum_b 
    
    Real(Kind=wp)     :: rmax, r_bin, rho_b, dV, r
    Real(Kind=wp)     :: rj(3), rk(3), rjk(3) 
    
    Integer(Kind=wi)  :: list_indx_a(model_data%config%num_atoms)
    Integer(Kind=wi)  :: list_indx_b(model_data%config%num_atoms)
    
    Character(Len=256) :: messages(3), message
    Character(Len=256) :: type_error
    Logical            :: modified, falloc, flag
    Logical            :: counted(model_data%config%num_atoms)
    Integer(Kind=wi)   :: fail(4)  
   
    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: gr(:)
    Real(Kind=wp),    Allocatable  :: nn(:)
    Real(Kind=wp),    Allocatable  :: cn(:)
    
    ! Search for the value of rmax 
    Do i = traj_data%analysis%frame_ini, traj_data%frames
      rmax=-Huge(1.0_wp)
      Do j = 1, 3
        If (traj_data%box(i)%cell_length(j) > rmax) Then
           rmax = traj_data%box(i)%cell_length(j)
        End If
      End Do
    End Do
    rmax=rmax/2.0_wp  
    
    ! Define number of bins
    nbins=Floor(rmax/traj_data%rdf%dr%value)

    ! In case &region is defined
    flag=.True.
    
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(gr(nbins), Stat=fail(2))
    Allocate(nn(nbins), Stat=fail(3))
    Allocate(cn(nbins), Stat=fail(4))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for RDF arrays. RDF analysis will not be executed.'
      falloc=.False.
    Else
      falloc=.True.
    End If
    
    If (falloc) Then
      gr=0.0_wp
      nn=0.0_wp
      ! Initiate Accumulators
      accum_a=0
      accum_b=0
      net_frames=0
      
      ! Compute the histogram for atoms of type a and b
      Do i = traj_data%analysis%frame_ini, traj_data%frames
        ! Define the number and list of indexes for type of species "a"
        num_at_a=0
        list_indx_a=0
        Do j = 1, model_data%config%num_atoms
          If (ANY(traj_data%rdf%type_a==Trim(traj_data%config(i,j)%tag))) Then
            If (traj_data%region%define%fread) Then
               Call within_region(traj_data, i, j, flag)
            End If
            If (flag) Then
              num_at_a=num_at_a+1
              list_indx_a(num_at_a)=j
            End If
          End If
        End Do
      
        ! Define the number and list of indexes for type of species "b"
        num_at_b=0
        list_indx_b=0
        Do j = 1, model_data%config%num_atoms
          If (ANY(traj_data%rdf%type_b==Trim(traj_data%config(i,j)%tag))) Then
            num_at_b=num_at_b+1
            list_indx_b(num_at_b)=j
          End If
        End Do
        
        ! Accummulators
        accum_b=accum_b+num_at_b
        accum_a=accum_a+num_at_a
      
        !Define rho_b 
        rho_b= num_at_b/traj_data%box(i)%volume   
        
        ! Calculate the histogram for this particular frame of the trajectory
        If (num_at_a /=0 .And. num_at_b/=0) Then
          h=0
          counted=.False.
          Do j=1, num_at_a 
            indx_a=list_indx_a(j)
            rj=traj_data%config(i,indx_a)%r
            If (traj_data%region%define%fread) Then
              Call within_region(traj_data, i, indx_a, flag)
            Else
              flag=.True.
            End If
            If (flag) Then
              Do k=1, num_at_b
                indx_b=list_indx_b(k)
                If (indx_a /= indx_b .And. (.Not. counted(indx_a)) .And. (.Not. counted(indx_b))) Then
                  rk=traj_data%config(i,indx_b)%r
                  rjk=rj-rk
                  Call check_PBC(rjk, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                  m=Floor(norm2(rjk)/traj_data%rdf%dr%value)+1
                  If (m <= nbins) Then
                    h(m)=h(m)+2
                  End If
                End If
              End Do
              counted(indx_a)=.True.
            End If
          End Do 
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do m=1, nbins 
            gr(m)= gr(m)+Real(h(m),Kind=wp)/(num_at_a*rho_b)
            nn(m)= nn(m)+Real(h(m),Kind=wp)/(num_at_b)
          End Do
        End If
        
      End Do
      
      ! Compute the radial distribution function gr
      If (accum_a /= 0 .And. accum_b /= 0) Then
        ! Print File
        Open(Newunit=files(FILE_RDF)%unit_no, File=files(FILE_RDF)%filename, Status='Replace')
        iunit=files(FILE_RDF)%unit_no
        Write (iunit,'(*(a,2x))') '#  Tags for type species "a":',&
                                & (Trim(traj_data%rdf%type_a(j)), j= 1, traj_data%rdf%num_type_a) 
        Write (iunit,'(*(a,2x))') '#  Tags for type species "b":',&
                                & (Trim(traj_data%rdf%type_b(j)), j= 1, traj_data%rdf%num_type_b) 
        Write (iunit,'(a)') '#  Radius [A]      RDF [1/A^3]       Coordination number' 
        
        cn=0.0_wp
        Do m=1, nbins
          r_bin=Real(m, Kind=wp)*traj_data%rdf%dr%value
          dV=4.0_wp*pi*r_bin**2*traj_data%rdf%dr%value 
          gr(m)=gr(m)/dV/Real(net_frames,Kind=wp)
          Do k= 1, m
           r=(Real(k,Kind=wp))*traj_data%rdf%dr%value
           cn(m)=cn(m)+nn(k)
          End Do
          cn(m)=cn(m)/Real(net_frames,Kind=wp)
          Write(iunit,'(2x,(3(f11.3, 6x)))') (Real(m,Kind=wp)-0.5)*traj_data%rdf%dr%value, gr(m), cn(m)
        End Do
        Write (message,'(1x,a)') 'The RDF analysis was printed to the "'//Trim(files(FILE_RDF)%filename)//'" file.'
        Call info(message, 1)
      Else
        If (accum_a == 0) Then
          type_error='"tags_species_a"'
        End If
        If (accum_b == 0) Then
          type_error='"tags_species_b"'
        End If
        If (accum_a == 0 .And. accum_b == 0) Then
          type_error='"tags_species_a" and "tags_species_b"'
        End If
        
        Write (messages(1),'(1x,a)') '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: RDF analysis could not be executed'
        If (traj_data%region%define%fread) Then
          Write (messages(2),'(1x,a)') '   Requested species as specified for '//Trim(type_error)//' in the &RDF&
                                  & block could not be identified along the trajectory for the selected region of&
                                  & the space (&region block).'
         Write (messages(3),'(1x,a)') '   Please verify the settings for the &RDF and &region blocks. The user should also&
                                & look at the file '//Trim(files(FILE_TAGGED_TRAJ)%filename)                        
        Else
          Write (messages(2),'(1x,a)') '   Requested species as specified for '//Trim(type_error)//' in the &RDF&
                                  & block could not be identified along the trajectory.'
          Write (messages(3),'(1x,a)') '   Please verify the settings for the &RDF block. The user should also&
                                & look at the file '//Trim(files(FILE_TAGGED_TRAJ)%filename)                        
        End If
        Call info(messages, 3)
        Write (messages(1),'(1x,a)') '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Close file
      Close(iunit)
      ! Deallocate arrays   
      Deallocate(cn, nn, gr, h)
    End If
    
  End Subroutine radial_distribution_function

  Subroutine mean_squared_displacement(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the Mean Squared Displacement (MSD) based on the
    ! settings of the &MSD block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)  :: i, j, k, l
    Integer(Kind=wi)  :: Nini_species, iunit, indx
    Real(Kind=wp)     :: time
    Real(Kind=wp)     :: base_time
    Logical           :: set_u0, fzero

    Character(Len=256) :: message, num_species
    Logical           :: terminated(traj_data%Nmax_species)

    ! Print tracked species
    Open(Newunit=files(FILE_MSD)%unit_no, File=files(FILE_MSD)%filename, Status='Replace')
    iunit=files(FILE_MSD)%unit_no
    Write(num_species,*) model_data%chem%N0%value 
    Write (iunit,'(a,8x,a)') '#  Time (ps)', '"'//Trim(traj_data%msd%select%type)//'"-MSD for species "'&
                              &//Trim(model_data%species_definition%name%type)//'" [A^2]' 

    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1
                          
    Do k= 1, traj_data%analysis%N_seg
      set_u0=.True.
      fzero=.False.
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
          End Do
        Else
          Nini_species=0
          Do j = 1, traj_data%Nmax_species
            If (traj_data%species(i,j)%alive) Then
              indx=traj_data%species(i,j)%list(1)
              traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx)%r
              traj_data%species(i,j)%u0(:,1)=traj_data%species(i,j)%u(:,1)
              Nini_species=Nini_species+1
            End If
          End Do
          set_u0=.False.
          If (Nini_species==0) Then
            Write (message,'(1x,a,2x,i6,a)') '***PROBLEMS: the code could not identify a single monitored species for frame ', i,&
                                            & '. Plase review the settings for the &monitored_species block'
            Call info(message, 1)
            Call error_stop(' ')
          End If
        End If
      
        traj_data%msd%r2=0.0_wp
        traj_data%N_species=0
        Do j = 1, traj_data%Nmax_species
          If(.Not. terminated(j)) Then
            If (traj_data%species(i,j)%alive) Then
              indx=traj_data%species(i,j)%list(1)
              traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx)%r
              Call msd_vector_difference(traj_data, i, j)
            Else
              terminated(j)=.True.
            End If
          End If  
        End Do
        
        If (i==traj_data%analysis%seg_indx(2,k)) Then
          If (traj_data%N_species /= 0) Then
            traj_data%msd%r2=traj_data%msd%r2/traj_data%N_species
            traj_data%analysis%variable(l,k)=traj_data%msd%r2
            Write(iunit,'(f11.3, 4x, f11.4)') (time-base_time)/1000.0_wp, traj_data%msd%r2
          End If  
          terminated(:)=.False.
          set_u0=.True.
          If ((traj_data%analysis%N_seg /=1) .And. (k /= traj_data%analysis%N_seg)) Then
            If (k /= traj_data%analysis%N_seg) Then
              Write (iunit,*) ' '
              Write (iunit,'(a)') '# Reseting the MSD analysis....' 
              Write (iunit,'(a,8x,a)') '#  Time (ps)', '"'//Trim(traj_data%msd%select%type)//'"-MSD for species "'&
                                &//Trim(model_data%species_definition%name%type)//'" [A^2]' 
            End If                    
          End If                      
        Else
          If ((traj_data%N_species) /= 0 .And. (.Not. fzero)) Then
            traj_data%msd%r2=traj_data%msd%r2/traj_data%N_species
            traj_data%analysis%variable(l,k)=traj_data%msd%r2
            Write(iunit,'(f11.3, 4x, f11.4)') (time-base_time)/1000.0_wp, traj_data%msd%r2
          Else  
            fzero=.True.
            traj_data%analysis%max_points(k)=i-1
          End If
        End If  
      End Do
    End Do

    Write (message,'(1x,a)') 'The MSD analysis was printed to the "'//Trim(files(FILE_MSD)%filename)//'" file.'
    Call info(message, 1)

    ! Close file
    Close(iunit)    

    If (traj_data%analysis%N_seg /=1 ) Then
      Call average_segments(files, traj_data, FILE_MSD_AVG, 'MSD')    
    End If

  End Subroutine mean_squared_displacement

  Subroutine msd_vector_difference(traj_data, i, j)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation from species j
    ! for the frame i (cij)
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j
    
    Logical           :: modified
    Real(Kind=wp)     :: du(3)
    Logical           :: flag
    Integer(Kind=wi)  :: m
    
    If (traj_data%region%define%fread) Then
      m=traj_data%species(i,j)%list(1)
      Call within_region(traj_data, i, m, flag)
    Else
      flag=.True.
    End If
    
    If (flag) Then
      traj_data%N_species=traj_data%N_species+1
      du=traj_data%species(i,j)%u(:,1)-traj_data%species(i,j)%u0(:,1)
      Call check_PBC(du, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      ! Recover if requested PBC
      If (traj_data%msd%pbc_xyz%fread) Then
        Do m= 1, 3
          If (.Not. traj_data%msd%pbc(m)) Then
            du(m)=traj_data%species(i,j)%u(m,1)-traj_data%species(i,j)%u0(m,1)   
          End If
        End Do
      End If

      Select Case (Trim(traj_data%msd%select%type))  
        Case ('x')
          traj_data%msd%r2 = traj_data%msd%r2 + du(1)**2
        Case ('y')
          traj_data%msd%r2 = traj_data%msd%r2 + du(2)**2
        Case ('z')
          traj_data%msd%r2 = traj_data%msd%r2 + du(3)**2
        Case ('xy')
          traj_data%msd%r2 = traj_data%msd%r2 + du(1)**2 + du(2)**2
        Case ('xz')
          traj_data%msd%r2 = traj_data%msd%r2 + du(1)**2 + du(3)**2
        Case ('yz')
          traj_data%msd%r2 = traj_data%msd%r2 + du(2)**2 + du(3)**2
        Case ('xyz')
          Do m= 1, 3
           traj_data%msd%r2 = traj_data%msd%r2 + du(m)**2
          End Do
      End Select  

    End If
    
  End Subroutine msd_vector_difference

  Subroutine residence_times(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the residence times for the changing species
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, m
    Integer(Kind=wi)   :: iunit
    Real(Kind=wp)      :: time
    Logical            :: set_u0
    Character(Len=256) :: message
    
!    Integer(Kind=wi)   :: indexes(2,model_data%chem%N0%value)
    Integer(Kind=wi)   :: ref_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: icount(model_data%chem%N0%value)    
    Real(Kind=wp)      :: values(traj_data%frames, model_data%chem%N0%value)
    Character(Len=8)   :: tag(traj_data%frames, 2, model_data%chem%N0%value)
    Character(Len=8)   :: ref_tag(model_data%chem%N0%value)
    Logical            :: hold(model_data%chem%N0%value)
    
    Real(Kind=wp)      :: rattling
    Real(Kind=wp)      :: tchange, base_time
    
    rattling=traj_data%lifetime%rattling_wait%value

    Do m = 1, model_data%chem%N0%value
      base_time=(traj_data%analysis%frame_ini-1)*traj_data%timestep%value
      set_u0=.True.
      icount(m)=0
      i = traj_data%analysis%frame_ini
      hold(m)=.False.
      Do While (i <= traj_data%frames)
        time=(i-1)*traj_data%timestep%value
        If (set_u0) Then
          ref_indx(m)=traj_data%track_chem%config(i,m)%indx
          ref_tag(m)=traj_data%track_chem%config(i,m)%tag
          set_u0=.False.
        Else
          If (traj_data%track_chem%config(i,m)%indx /= ref_indx(m)) Then
            If (.Not. hold(m)) Then
              tchange=time
              hold(m)=.True.
            End If
          Else 
            hold(m)=.False.
          End If  
          If (hold(m)) Then
            If (((time-tchange) > rattling .Or. i==traj_data%frames)) Then
              hold(m)=.False.
              icount(m)=icount(m)+1
              values(icount(m),m)=(tchange-base_time)/1000.0_wp
              tag(icount(m),1,m)=ref_tag(m)
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
              ref_tag(m)=traj_data%track_chem%config(i,m)%tag
              base_time=tchange
            End If
          End If
        End If
        i=i+1
      End Do
    End Do
    
    Open(Newunit=files(FILE_RES_TIMES)%unit_no, File=files(FILE_RES_TIMES)%filename, Status='Replace')
    iunit=files(FILE_RES_TIMES)%unit_no
    
    Do m = 1, model_data%chem%N0%value
    Write (iunit,'(a,i3)') '#  Species', m 
    Write (iunit,'(a)') '#  Residence Time (ps)     Tag for site' 
      Do i =1, icount(m)
        Write(iunit,'(f11.3,15x,a)') values(i,m), Trim(tag(i,1,m))
      End Do
    Write (iunit,'(a)') ' ' 
    End Do 
    
    Write (message,'(1x,a)') 'The Residence Times for each changing chemical species were&
                             & printed to the "'//Trim(files(FILE_RES_TIMES)%filename)//'" file.'
    Call info(message, 1)
    Close(iunit)
    
  End Subroutine residence_times

  Subroutine transfer_correlation_function(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the transfer correlation function (TCF) involving
    ! the changing chemistry species
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, j, k, l, m, n, ntop
    Integer(Kind=wi)   :: iunit, ini_indx
    Real(Kind=wp)      :: suma_i 
    Real(Kind=wp)      :: time
    Real(Kind=wp)      :: base_time
    Logical            :: set_u0
    Character(Len=256) :: message
    
    Logical            :: terminated(traj_data%analysis%Ninterval, model_data%chem%N0%value)
    Integer(Kind=wi)   :: indexes(2,model_data%chem%N0%value)
    Integer(Kind=wi)   :: ref_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: ref_indx2(model_data%chem%N0%value)
    Integer(Kind=wi)   :: time_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: Nnet
    
    Logical            :: follow(model_data%chem%N0%value)
    Logical            :: flag, first(model_data%chem%N0%value)
    Logical            :: hold(model_data%chem%N0%value)

    Real(Kind=wp)      :: tchange(model_data%chem%N0%value)
    Character(Len=256) :: method
    Real(Kind=wp)      :: rattling
    
    method=Trim(traj_data%lifetime%method%type)
    rattling=traj_data%lifetime%rattling_wait%value
    
    ! Print header
    Open(Newunit=files(FILE_TCF)%unit_no, File=files(FILE_TCF)%filename, Status='Replace')
    iunit=files(FILE_TCF)%unit_no
    Write (iunit,'(a)') '#  Time (ps)          TCF' 

    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1
    
    Do k= 1, traj_data%analysis%N_seg
      set_u0=.True.
      follow=.True.
      first=.True.
      tchange=0.0_wp
      hold=.False.
      ! Initialise terminated tag
      Do m = 1, model_data%chem%N0%value
        terminated(:,m)=.False.
      End Do
      
      ini_indx=traj_data%analysis%seg_indx(1,k)
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        time=(i-1)*traj_data%timestep%value
        If (Trim(method)=='hicf' .Or. Trim(method)=='hdcf') Then
          If (set_u0) Then
            Do m = 1, model_data%chem%N0%value
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
            End Do  
            set_u0=.False.
          Else
            Do m = 1, model_data%chem%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%indx/=ref_indx(m)) Then
                  If (.Not. hold(m)) Then
                    tchange(m)=time
                    hold(m)=.True.
                    time_indx(m)=i
                  End If
                Else 
                  hold(m)=.False.
                End If  
                If (((time-tchange(m)) > rattling) .And. hold(m)) Then
                  If (Trim(method)=='hicf') Then
                    ntop=i
                  Else If (Trim(method)=='hdcf') Then
                    ntop=traj_data%analysis%seg_indx(2,k)
                    follow(m)=.False.
                  End If 
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                  hold(m)=.False.
                End If
                If ((i==traj_data%analysis%seg_indx(2,k)) .And. hold(m)) Then
                  ntop=i
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                End If
              End If
            End Do
          End If
        End If
        
        If (Trim(method)=='hdcf*' .Or. Trim(method)=='hdcf-2s') Then
          If (set_u0) Then
            Do m = 1, model_data%chem%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%indx
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
            End Do  
            set_u0=.False.
          Else
            Do m = 1, model_data%chem%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%indx/=indexes(2,m)) Then
                  indexes(1,m)=indexes(2,m)
                  indexes(2,m)=traj_data%track_chem%config(i,m)%indx
                  If (first(m)) Then
                    ref_indx2(m)=traj_data%track_chem%config(i,m)%indx
                    first(m)=.False.
                  End If
                  If (Trim(method)=='hdcf*') Then
                    If (indexes(1,m) /= ref_indx(m) .And. indexes(2,m) /= ref_indx(m)) Then
                      flag=.True.
                    Else
                      flag=.False.
                    End If
                  Else If (Trim(method)=='hdcf-2s') Then
                   If (indexes(2,m) /= ref_indx2(m)) Then
                     flag=.True.
                   Else
                     flag=.False.
                   End If
                  End If
                  
                  If (flag) Then
                    If (.Not. hold(m)) Then
                     tchange(m)=time
                     hold(m)=.True.
                     time_indx(m)=i
                    End If
                  Else
                    hold(m)=.False.
                  End If
                End If
                If (((time-tchange(m)) > rattling) .And. hold(m)) Then
                  Do n=time_indx(m), traj_data%analysis%seg_indx(2,k)
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                  follow(m)=.False.
                End If
                If ((i==traj_data%analysis%seg_indx(2,k)) .And. hold(m)) Then
                  ntop=i
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                End If
              End If
            End Do
          End If
        End If
   
      End Do
      
      flag=.True.
      l=0
      base_time=(traj_data%analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        time=(i-1)*traj_data%timestep%value
        l=l+1
        If (flag) Then
          If (All(terminated(l,:))) Then
            suma_i=0.0_wp
!             flag=.False.
!             traj_data%analysis%max_points(k)=l-1
          Else
            suma_i=0.0_wp
            Nnet=0
            Do j = 1, model_data%chem%N0%value
              If(.Not. terminated(l,j)) Then
                suma_i=suma_i+1.0_wp    
                Nnet=Nnet+1
              End If  
            End Do
            suma_i=suma_i/model_data%chem%N0%value
          End If
         
          If (flag) Then
            traj_data%analysis%variable(l,k)=suma_i
            Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If
      
        End If
        If (i==traj_data%analysis%seg_indx(2,k) .And. (traj_data%analysis%N_seg /=1)) Then
           If (k /= traj_data%analysis%N_seg) Then
             Write (iunit,*) ' '
             Write (iunit,'(a)') '# Reseting the TCF analysis....' 
             Write (iunit,'(a)') '#  Time (ps)          TCF' 
           End If
        End If
      End Do
    End Do
    
    Write (message,'(1x,a)') 'The TFC analysis was&
                             & printed to the "'//Trim(files(FILE_TCF)%filename)//'" file.'
    Call info(message, 1)
    Close(iunit)
    
    If (traj_data%analysis%N_seg /=1 ) Then
      Call average_segments(files, traj_data, FILE_TCF_AVG, 'TCF')    
    End If

  End Subroutine transfer_correlation_function
  
  Subroutine average_segments(files, traj_data, file_number, what)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to average physical quantities computed for each time segment
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: file_number
    Character(Len=*),  Intent(In   ) :: what
  
    Integer(Kind=wi) :: i, k, iunit, Nnet 
    Real(Kind=wp)    :: sum_i, average, std, maximum, minimum
    Logical          :: flag
    Character(Len=256) :: message
    ! Print header
    Open(Newunit=files(file_number)%unit_no, File=files(file_number)%filename, Status='Replace')
    iunit=files(file_number)%unit_no
    Write (iunit,'(a)') '#  Time (ps)      '//Trim(what)//'         '//Trim(what)//' (max)     '&
                       &//Trim(what)//' (min)      STD' 
    i=1
    flag=.True.    
    Do While ((i<=traj_data%analysis%Ninterval) .And. flag)
      sum_i=0.0_wp
      Nnet=0
      Do k= 1, traj_data%analysis%N_seg
        If(i<=traj_data%analysis%max_points(k)) Then
          Nnet=Nnet+1
          sum_i=sum_i+traj_data%analysis%variable(i,k)
        End If
      End Do

      ! Compute average
      If (Nnet > 0) Then
        average=sum_i/Nnet
        If (Nnet>1) Then
          sum_i=0.0_wp
          Do k= 1, traj_data%analysis%N_seg
            If(i<=traj_data%analysis%max_points(k)) Then
              sum_i=sum_i+(traj_data%analysis%variable(i,k)-average)**2
            End If
          End Do
          std=sqrt(sum_i/(Nnet-1))
        Else
          std=0.0_wp
        End If
      Else
        flag=.False.
      End If
      
      If (flag) Then
        If(average+std > 1.0_wp) Then
          maximum=1.0_wp
        Else
          maximum=average+std
        End If
        minimum=average-std
        Write(iunit,'(5(f10.3, 3x))') (i-1)*traj_data%timestep%value/1000.0_wp, average, maximum, minimum, std
      End If
      i=i+1
    End Do
    Write (message,'(1x,a)') 'The average '//Trim(what)//' was printed to the "'//&
                             &Trim(files(file_number)%filename)//'" file.'
    Call info(message, 1)
    Close(iunit)
  
  End Subroutine average_segments
  
  Subroutine orientational_correlation_function(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the orientational correlation function (OCF)
    ! Different possible flavours are available depending on the settings
    ! of the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data

    Integer(Kind=wi)  :: i, j, k, l
    Integer(Kind=wi)  :: Nini_species, iunit
    Real(Kind=wp)     :: suma_i 
    Real(Kind=wp)     :: time
    Real(Kind=wp)     :: base_time
    Logical           :: set_u0, fzero
    Character(Len=256) :: message
    
    Logical           :: terminated(traj_data%Nmax_species)

    ! Print header
    Open(Newunit=files(FILE_OCF)%unit_no, File=files(FILE_OCF)%filename, Status='Replace')
    iunit=files(FILE_OCF)%unit_no
    Write (iunit,'(a)') '#  Time (ps)           OCF' 

    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1
    
    Do k= 1, traj_data%analysis%N_seg
      set_u0=.True.
      l=0
      fzero=.False.
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
            If (Trim(traj_data%ocf%u_definition%type) == 'bond_12-13') Then
              traj_data%species(i,:)%u0(j,2)=traj_data%species(i-1,:)%u0(j,2)
            End If
          End Do
        Else
          Nini_species=0
          Do j = 1, traj_data%Nmax_species
            If (traj_data%species(i,j)%alive) Then
              Call compute_rotation_vector(traj_data, i, j)
              traj_data%species(i,j)%u0(:,1)=traj_data%species(i,j)%u(:,1)
              If (Trim(traj_data%ocf%u_definition%type) == 'bond_12-13') Then
                traj_data%species(i,j)%u0(:,2)=traj_data%species(i,j)%u(:,2)
              End If
              Nini_species=Nini_species+1
            End If
          End Do
          set_u0=.False.
          If (Nini_species==0) Then
            Write (message,'(1x,a,2x,i6,a)') '***PROBLEMS: the code could not identify a single monitored species for frame ', i,&
                                            & '. Plase review the settings for the &monitored_species block'
            Call info(message, 1)
            Call error_stop(' ')
          End If
        End If
      
        suma_i=0.0_wp
        traj_data%N_species=0
        Do j = 1, traj_data%Nmax_species
          If(.Not. terminated(j)) Then
            If (traj_data%species(i,j)%alive) Then
              Call compute_rotation_vector(traj_data, i, j)
              Call evaluate_correlation_term(traj_data, i, j, suma_i)  
            Else
              terminated(j)=.True.
            End If
          End If  
        End Do

        If (i==traj_data%analysis%seg_indx(2,k)) Then
          If (traj_data%N_species /= 0) Then
           suma_i=suma_i/traj_data%N_species
           traj_data%analysis%variable(l,k)=suma_i
           Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If  
          terminated(:)=.False.
          set_u0=.True.
          If ((traj_data%analysis%N_seg /=1) .And. (k /= traj_data%analysis%N_seg)) Then
            If (k /= traj_data%analysis%N_seg) Then
              Write (iunit,*) ' '
              Write (iunit,'(a)') '# Reseting the OCF analysis....' 
              Write (iunit,'(a)') '#  Time (ps)           OCF' 
            End If
          End If                      
        Else
          If ((traj_data%N_species) /= 0 .And. (.Not. fzero)) Then
            suma_i=suma_i/traj_data%N_species
            traj_data%analysis%variable(l,k)=suma_i
            Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          Else  
            fzero=.True.
            traj_data%analysis%max_points(k)=i-1
          End If
        End If
      End Do
    End Do
    
    Write (message,'(1x,a)') 'The orientational correlation function analysis was printed to the "'//&
                             &Trim(files(FILE_OCF)%filename)//'" file.'
    Call info(message, 1)

    Close(iunit)
    
    If (traj_data%analysis%N_seg /=1 ) Then
      Call average_segments(files, traj_data, FILE_OCF_AVG, 'OCF')    
    End If

  End Subroutine orientational_correlation_function

  Subroutine compute_rotation_vector(traj_data, i, j)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the rotation vector
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j

    Integer(Kind=wi) :: indx1, indx2, indx3, k
    Logical          :: modified
    Real(Kind=wp), Dimension(3)  :: u12, u13


    indx1=traj_data%species(i,j)%list(1)
    indx2=traj_data%species(i,j)%list(2)
    indx3=traj_data%species(i,j)%list(3)
    
    If (Trim(traj_data%ocf%u_definition%type) == 'bond_12') Then
      traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      Call check_PBC(traj_data%species(i,j)%u(:,1), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(traj_data%ocf%u_definition%type) == 'bond_13') Then
      traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(traj_data%species(i,j)%u(:,1), traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(traj_data%ocf%u_definition%type) == 'bond_12-13') Then
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
    Else If (Trim(traj_data%ocf%u_definition%type) == 'bond_123') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Do k=1, 3
        traj_data%species(i,j)%u(k,1)=u12(k)+u13(k)
      End Do
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    Else If (Trim(traj_data%ocf%u_definition%type) == 'plane') Then
      u12=traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
      u13=traj_data%config(i,indx3)%r-traj_data%config(i,indx1)%r
      Call check_PBC(u12, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call check_PBC(u13, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
      Call cross_product(u12, u13, traj_data%species(i,j)%u(:,1))
      traj_data%species(i,j)%u(:,1)=traj_data%species(i,j)%u(:,1)/norm2(traj_data%species(i,j)%u(:,1))
    End If    
    
  End Subroutine compute_rotation_vector  

  Subroutine within_region(traj_data, i, m, flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if the position of the 
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer,           Intent(In   ) :: m
    Logical,           Intent(  Out) :: flag        

    Integer(Kind=wi)  :: k, j

    Logical :: fpass(3)
   
    Do k = 1, 3
      Do j = 1, traj_data%region%number(k)
        If (traj_data%region%inside(k,j)) Then
          If (traj_data%region%domain(k,1,j) <= traj_data%config(i,m)%r(k) .And. &
              traj_data%region%domain(k,2,j) >= traj_data%config(i,m)%r(k)) Then
            traj_data%region%belong(k,j) = .True.
          Else
            traj_data%region%belong(k,j) = .False.
          End If       
        Else
          If (traj_data%region%domain(k,1,j) >  traj_data%config(i,m)%r(k) .Or. &
              traj_data%region%domain(k,2,j) <  traj_data%config(i,m)%r(k)) Then
            traj_data%region%belong(k,j) = .True.
          Else
            traj_data%region%belong(k,j) = .False.
          End If       
        End If
        If (j==1) Then
          fpass(k)=traj_data%region%belong(k,j)
        Else
          fpass(k)=fpass(k) .Or. traj_data%region%belong(k,j)
        End If
      End Do
    End Do 

    flag=fpass(1) .And. fpass(2) .And. fpass(3)
    
  End Subroutine within_region
  
  Subroutine evaluate_correlation_term(traj_data, i, j, suma_i)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation from species j
    ! for the frame i (cij)
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
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
      Select Case (traj_data%ocf%legendre_order%value)  
        Case (1)
          cij=x
        Case (2)
          cij=(3.0_wp*(x)**2-1.0_wp)/2.0_wp
        Case (3)
          cij=(5.0_wp*(x)**3-3.0_wp*x)/2.0_wp
        Case (4)
          cij=(35.0_wp*(x)**4-30.0_wp*x**2+3.0_wp)/8.0_wp
      End Select  
      
      If (Trim(traj_data%ocf%u_definition%type) == 'bond_12-13') Then
        x=Dot_product(traj_data%species(i,j)%u(:,2),traj_data%species(i,j)%u0(:,2))
        c0ij=cij
        Select Case (traj_data%ocf%legendre_order%value)  
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

  End Subroutine evaluate_correlation_term  
  
  Subroutine obtain_number_frames(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to obtain the number of frames recorded in the TRAJECTORY file 
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data

    Logical            :: loop_traj
    Character(Len=256) :: check
    Integer(Kind=wi)   :: i, stat

    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=Trim(files(FILE_TRAJECTORY)%filename),Status='old')

    i=1
    loop_traj=.True.
    Do While (loop_traj)
      Call read_model(files, model_data, i, traj_data%ensemble%type)
      Read(files(FILE_TRAJECTORY)%unit_no, Fmt= *, iostat=stat) check
      If (is_iostat_end(stat)) Then
        loop_traj=.False.
      Else 
        backspace files(FILE_TRAJECTORY)%unit_no
      End If
      i=i+1
    End Do
    traj_data%frames=i-1
    
    Close(files(FILE_TRAJECTORY)%unit_no) 
    
  End Subroutine obtain_number_frames

  Subroutine check_trajectory_settings(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of trajectory-related directives
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%print_retagged_trajectory%fread) Then
      If (traj_data%print_retagged_trajectory%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_retagged_trajectory" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
      If((.Not. model_data%change_chemistry%stat) .And. traj_data%print_retagged_trajectory%stat) Then 
        Write (messages(1),'(2(1x,a))') Trim(error_set), ' The user has set "print_retagged_trajectory" to .True. but&
                                      & "change_chemistry" is set to .False. Why do you want to retag the trajectory?&
                                      & Please change'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      traj_data%print_retagged_trajectory%stat=.False.
    End If

    ! Check timestep
    Call check_time_directive(traj_data%timestep, error_set, .True.)
    
    ! Check ensemble
    If (traj_data%ensemble%fread) Then
      If (traj_data%ensemble%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "ensemble" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%ensemble%type)/='nve'  .And. &
            Trim(traj_data%ensemble%type)/='nvt'  .And. &
            Trim(traj_data%ensemble%type)/='npt') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "ensemble". Valid options: "NVE", "NVT" and "NPT"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "ensemble" directive'
       Call info(messages, 1)
       Call error_stop(' ')
    End If

   ! Check directives for data analysis
    Call check_data_analysis(files, traj_data)
    
    ! Check info &OCF block 
    If (traj_data%ocf%invoke%fread) Then
      Call check_ocf(files, traj_data)
    End If 

    ! Check settings &region block 
     Call check_region(files, traj_data)
    
    ! Check info &MSD block 
    If (traj_data%msd%invoke%fread) Then
      Call check_msd(files, traj_data)
    End If 

    ! Check info &RDF block 
    If (traj_data%rdf%invoke%fread) Then
      Call check_rdf(files, traj_data, model_data)
    End If 

    ! Check info &track_unchanged_chemistry block 
    If (traj_data%unchanged%invoke%fread) Then
      Call check_unchanged_chemistry(files, traj_data, model_data)
    End If 

    ! Check info &lifetime block 
    If (traj_data%lifetime%invoke%fread) Then
      Call check_lifetime(files, traj_data)
    End If 
    
  End Subroutine check_trajectory_settings

  Subroutine check_region(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &egion block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer             :: k, m, j
    
    error_set = '***ERROR in the &region block of file '//Trim(files(FILE_SET)%filename)//' -'

    m=0
    If (traj_data%region%define%fread) Then
      Do k=1,3
        Do j = 1, traj_data%region%number(k)
          If (traj_data%region%invoke(k,j)%fread) Then
            If (traj_data%region%invoke(k,j)%fail) Then
              Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the '//&
                                             &Trim(traj_data%region%invoke(k,j)%type)//' directive.'
              Call info(messages, 1)
              Call error_stop(' ')
            Else
              If (traj_data%region%domain(k,1,j) > traj_data%region%domain(k,2,j)) Then
                If (traj_data%region%number(k)==1) Then
                  Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                            &'The lower value of the domain for defined "'&
                                            &//Trim(traj_data%region%invoke(k,j)%type)//&
                                            &'" is larger than the upper value!!! Please change.'
                Else
                  Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                            &'The lower value of the domain for one of the defined "'&
                                            &//Trim(traj_data%region%invoke(k,j)%type)//&
                                            &'" is larger than the upper value!!! Please change.'
                End If
                Call info(messages, 1)
                Call error_stop(' ')
              End If
              If (Abs(traj_data%region%domain(k,1,j) - traj_data%region%domain(k,2,j))<epsilon(1.0_wp)) Then
                If (traj_data%region%number(k)==1) Then
                  Write (messages(1),'(2(1x,a))') Trim(error_set),& 
                                          &'The lower and upper values of the domain for the defined "'&
                                          &//Trim(traj_data%region%invoke(k,j)%type)//&
                                          &'" are exaclty the same! Please change.'
                Else
                  Write (messages(1),'(2(1x,a))') Trim(error_set),&
                                          &'The lower and upper values of the domain for one of the defined "'&
                                          &//Trim(traj_data%region%invoke(k,j)%type)//&
                                          &'" are exaclty the same! Please change.'
                End If
                Call info(messages, 1)
                Call error_stop(' ')
              End If
              Call capital_to_lower_case(traj_data%region%inout(k,j))    
              If (Trim(traj_data%region%inout(k,j)) /= 'inside' .And. Trim(traj_data%region%inout(k,j)) /= 'outside') Then
                 If (traj_data%region%number(k)==1) Then
                   Write (messages(1),'(2(1x,a))')  Trim(error_set),'The last argument of the defined directive "'&
                                                    &//Trim(traj_data%region%invoke(k,j)%type)//&
                                                    &'" must be either "inside" or "outside", referring&
                                                    & to the region defined by the limits. Please change.'
                 Else
                   Write (messages(1),'(2(1x,a))')  Trim(error_set),'The last argument for one of the defined directives "'&
                                                    &//Trim(traj_data%region%invoke(k,j)%type)//&
                                                    &'" must be either "inside" or "outside", referring&
                                                    & to the region defined by the limits. Please change.'
                 End If
                 Call info(messages, 1)
                 Call error_stop(' ')
              Else
                If (Trim(traj_data%region%inout(k,j)) == 'inside') Then
                  traj_data%region%inside(k,j)=.True.
                Else If (Trim(traj_data%region%inout(k,j)) == 'outside') Then
                  traj_data%region%inside(k,j)=.False.
                End If
              End If
            End If
          Else
            m=m+1
            traj_data%region%inside(k,j)=.True.
            traj_data%region%domain(k,1,j)=-Huge(1.0_wp)
            traj_data%region%domain(k,2,j)= Huge(1.0_wp)
          End If
        End Do
      End Do
    End If
    
    If (m==3) Then
       Write (messages(1),'(1x,a)') 'WARNING: the &region block contains no data!' 
       Call info(messages, 1)
       traj_data%region%define%fread=.False.
    End If
    
  End Subroutine check_region

  Subroutine check_unchanged_chemistry(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &track_unchanged_chemistry block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data

    Character(Len=256)  :: messages(2), word
    Character(Len=256)  :: error_set
    Integer(Kind=wi)    :: j, k
    Logical             :: flag
    
    error_set = '***ERROR in the &track_unchanged_chemistry block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%unchanged%tag%fread) Then
      If (traj_data%unchanged%tag%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "tag" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      ! Check if all tags correspond to the same element (type a)
      j=1
      flag=.True.
      Do While (j <= model_data%input_composition%atomic_species .And. flag)
        If (Trim(model_data%input_composition%tag(j))==Trim(traj_data%unchanged%tag%type)) Then
          flag=.False.
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'The atomic tag "'//Trim(traj_data%unchanged%tag%type)//&
                                       &'" (defined for the "tag" directive) has not been defined&
                                       & in the &input_composition block! Please review the settings' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If 
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must the "tag" (atomic tag)&
                                    & to track along the trajectory'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    If (traj_data%unchanged%list_indexes%fread) Then
      If (traj_data%unchanged%list_indexes%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "list_indexes" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If

      Do j=1, traj_data%unchanged%N0-1
        Do k=j+1, traj_data%unchanged%N0
          If (traj_data%unchanged%indexes(j)==traj_data%unchanged%indexes(k)) Then
            Write(word,*) traj_data%unchanged%indexes(j)
            Write (messages(1),'(2(1x,a))') Trim(error_set), 'Index "'//Trim(Adjustl(word))//' is repeated in the list!'
            Write (messages(2),'((1x,a))') 'Values in the "list_indexes" must be  different'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End Do
      End Do 
      
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define "list_indexes" for&
                                    & all those atoms that the user wants to print'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
  End Subroutine check_unchanged_chemistry
  
  Subroutine check_initial_unchanged_labels(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the atomic tages for each component of the 
    ! list_indexes (&track_unchanged_chemistry block) is the same as the "tag"
    ! directive defined
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data

    Character(Len=256)  :: messages(2), word
    Character(Len=256)  :: error_set
    Integer(Kind=wi)    :: j, k
    
    error_set = '***ERROR in the &track_unchanged_chemistry block of file '//Trim(files(FILE_SET)%filename)//' -'

    Do j=1, traj_data%unchanged%N0
      k=traj_data%unchanged%indexes(j) 
      If (Trim(model_data%config%atom(k)%tag)/=Trim(traj_data%unchanged%tag%type)) Then
        Call info(' ', 1)
        Write(word,*) k
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Index "'//Trim(Adjustl(word))//'" (defined in list_indexes)&
                                       & does not correspond to the atomic tag "'//Trim(traj_data%unchanged%tag%type)//'".'  
        Write (messages(2),'((1x,a))') 'According to the &input_composition block, this index corresponds to atomic&
                                       & tag "'//Trim(model_data%config%atom(k)%tag)//'". Please review the labels for the model'
        Call info(messages, 2)
        Call error_stop(' ')
      End If
    End Do 
      
  End Subroutine check_initial_unchanged_labels
  
  Subroutine check_rdf(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &RDF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer(Kind=wi)    :: j, k
    Logical             :: flag

    Character(Len=8)  :: tg(max_components)
    Character(Len=8)  :: el(max_components)

    error_set = '***ERROR in the &RDF block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (.Not. traj_data%rdf%dr%fread) Then
      traj_data%rdf%dr%tag='dr'
    End If
    Call check_length_directive(traj_data%rdf%dr, error_set, .True., 'directive')
    
    If (traj_data%rdf%tags_species_a%fread) Then
      If (traj_data%rdf%tags_species_a%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "tags_species_a" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      Do j=1, traj_data%rdf%num_type_a-1
        Do k=j+1, traj_data%rdf%num_type_a
          If (Trim(traj_data%rdf%type_a(j))==Trim(traj_data%rdf%type_a(k))) Then
            Write (messages(1),'(4(1x,a))') Trim(error_set), 'Tag', Trim(traj_data%rdf%type_a(j)), 'is repeated in the list!'
            Write (messages(2),'((1x,a))') 'The tags defined in "tags_species_a" must be  different'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End Do
      End Do 
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "tags_species_a" directive for RDF analysis.&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    If (traj_data%rdf%tags_species_b%fread) Then
      If (traj_data%rdf%tags_species_b%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "tags_species_b" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      Do j=1, traj_data%rdf%num_type_b-1
        Do k=j+1, traj_data%rdf%num_type_b
          If (Trim(traj_data%rdf%type_b(j))==Trim(traj_data%rdf%type_b(k))) Then
            Write (messages(1),'(4(1x,a))') Trim(error_set), 'Tag', Trim(traj_data%rdf%type_b(j)), 'is repeated in the list!'
            Write (messages(2),'((1x,a))') 'The tags defined in "tags_species_b" must be  different'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End Do
      End Do 
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "tags_species_b" directive for RDF analysis.&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If

    ! Check if all tags correspond to the same element (type a)
    Do k=1, traj_data%rdf%num_type_a
      tg(k)=Trim(traj_data%rdf%type_a(k))
      Call remove_symbols(tg(k),'*')
      flag=.True.
      j=1
      Do While (j <= model_data%input_composition%atomic_species .And. flag)
        If (Trim(model_data%input_composition%tag(j))==Trim(tg(k))) Then
          flag=.False.
          el(k)=Trim(model_data%input_composition%element(j))
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'The tag " '//Trim(traj_data%rdf%type_a(k))//' " of the "tags_species_a"&
                                       & is not a valid option. Please review the definition of the &input_composition block' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If 
    End Do
    
    Do k= 1, traj_data%rdf%num_type_a-1
      Do j= k, traj_data%rdf%num_type_a
        If (Trim(el(k)) /= Trim(el(j))) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Tags " '//Trim(traj_data%rdf%type_a(k))//' " and " '&
                                      &//Trim(traj_data%rdf%type_a(j))//' " defined in "tags_species_b" correspond to two different&
                                      & chemical elements!' 
          Call info(messages, 1)
          Call error_stop(' ') 
        End If
      End Do
    End Do   
    
    ! Check if all tags correspond to the same element (tybe b)
    Do k=1, traj_data%rdf%num_type_b
      tg(k)=Trim(traj_data%rdf%type_b(k))
      Call remove_symbols(tg(k),'*')
      flag=.True.
      j=1
      Do While (j <= model_data%input_composition%atomic_species .And. flag)
        If (Trim(model_data%input_composition%tag(j))==Trim(tg(k))) Then
          flag=.False.
          el(k)=Trim(model_data%input_composition%element(j))
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'The tag " '//Trim(traj_data%rdf%type_b(k))//' " of the "tags_species_b"&
                                       & is not a valid option. Please review the definition of the &input_composition block' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If 
    End Do
    
    Do k= 1, traj_data%rdf%num_type_b-1
      Do j= k, traj_data%rdf%num_type_b
        If (Trim(el(k)) /= Trim(el(j))) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Tags " '//Trim(traj_data%rdf%type_b(k))//' " and " '&
                                      &//Trim(traj_data%rdf%type_b(j))//' " defined in "tags_species_b" correspond to two different&
                                      & chemical elements!' 
          Call info(messages, 1)
          Call error_stop(' ') 
        End If
      End Do
    End Do   
    
  End Subroutine check_rdf

  Subroutine check_data_analysis(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &data_analysis block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set

    error_set = '***ERROR in the &data_analysis block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%analysis%invoke%fread) Then
      Call check_time_directive(traj_data%analysis%time_interval,  error_set, .False.)
      Call check_time_directive(traj_data%analysis%ignore_initial, error_set, .False.)
      Call check_time_directive(traj_data%analysis%overlap_time,  error_set, .False.)
    End If

  End Subroutine check_data_analysis

  Subroutine check_lifetime(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &lifetime block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set
    Character(Len=256)  :: messages(2)

    error_set = '***ERROR in the &lifetime block of file '//Trim(files(FILE_SET)%filename)//' -'

    Call check_time_directive(traj_data%lifetime%rattling_wait,  error_set, .False.)

    If (.Not. traj_data%lifetime%rattling_wait%fread) Then
      traj_data%lifetime%rattling_wait%value= 0.0_wp
      traj_data%lifetime%rattling_wait%units= 'fs'
    End If
    
    If (traj_data%lifetime%method%fread) Then
      If (traj_data%lifetime%method%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "method" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%lifetime%method%type)/='hicf'  .And. &
            Trim(traj_data%lifetime%method%type)/='hdcf'  .And. &
            Trim(traj_data%lifetime%method%type)/='hdcf*'  .And. &
            Trim(traj_data%lifetime%method%type)/='hdcf-2s') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "method". Valid options: "HICF", "HDCF", "HDCF*" and "HDCF-2s"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "method" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "HICF", "HDCF", "HDCF*" and "HDCF-2s"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If
    
  End Subroutine check_lifetime

  Subroutine check_time_directive(T, error_set, kill)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check time related directivesd
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: error_set
    Logical,                 Intent(In   )  :: kill

    Character(Len=256)  :: messages(2)
    
    If (T%fread) Then
      If (T%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "'&
                                      &//Trim(T%tag)//'" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (T%value < epsilon(1.0_wp)) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Input value for "'//Trim(T%tag)//&
                                    &'" MUST be larger than zero'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        Call capital_to_lower_case(T%units)
        If (Trim(T%units) /= 'fs' .And. &
           Trim(T%units) /= 'ps') Then
           Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                                    & 'Units for directive "'//Trim(T%tag)//&
                                    &'" must be "fs" or "ps". Have you included the units?'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
        ! Transform to fs
        If (Trim(T%units) == 'ps') Then
           T%value=1000_wp* T%value
        End If
      End If
    Else 
      If (kill)then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(T%tag)//'" directive'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_time_directive  
  
  Subroutine check_time_settings(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &data_analysis block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set, message
    Integer(Kind=wi)    :: i, k, j, l, kref, kini
    Logical             :: flag, one_seg
    Real(Kind=wp)       :: time, tend, teff, tini, tref
    
    error_set = '***ERROR in the &data_analysis block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (.Not. traj_data%analysis%ignore_initial%fread) Then
       traj_data%analysis%ignore_initial%value=-traj_data%timestep%value
    End If
    
    tend=(traj_data%frames-1)*traj_data%timestep%value 
    If (tend <= traj_data%analysis%ignore_initial%value) Then
      Call info(' ', 1)
      Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%ignore_initial%tag)//&
                               &'" is larger than the total time for the trajectory. Please check&
                               & the settings and the value for the "timestep" directive.'
      Call info(message, 1)
      Call error_stop(' ')
    End If   

    i=1
    flag=.True.
    Do While (i <= traj_data%frames .And. flag)
      time=(i-1)*traj_data%timestep%value
      If (time >= traj_data%analysis%ignore_initial%value) Then
        traj_data%analysis%frame_ini = i
        tini=(i-1)*traj_data%timestep%value
        flag=.False.
      End If
      i=i+1
    End do 
    
    If (.Not. traj_data%analysis%time_interval%fread) Then
       traj_data%analysis%time_interval%value=tend
    End If

    teff=tend-tini
    If (.Not. traj_data%analysis%overlap_time%fread) Then
       traj_data%analysis%overlap_time%value=teff+traj_data%timestep%value 
    End If

    ! Compare timestep with other time settings of &data_analysis
    If (traj_data%timestep%value>=traj_data%analysis%time_interval%value) Then
      Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%time_interval%tag)//&
                               &'" must be larger that the timestep for the trajectory. Please check&
                               & the value (and units) for the "timestep" directive.'
      Call info(message, 1) 
      Call error_stop(' ')
    End If

    If (traj_data%analysis%overlap_time%fread) Then
      If (traj_data%timestep%value > traj_data%analysis%overlap_time%value) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%overlap_time%tag)//&
                                 &'" must be larger that the timestep for the trajectory. Please check values (and units).'
        Call info(message, 1) 
        Call error_stop(' ')
      End If
    End If
    
    ! Calculate the number of segments
    i=0; j=0; l=0
    one_seg=.True.
    k=traj_data%analysis%frame_ini
    tref=tini; kini=k; flag=.True.
    Do While (k <= traj_data%frames)
      time=(k-1)*traj_data%timestep%value
      If (time>=(tref+traj_data%analysis%time_interval%value)) Then
        i=i+1  
        l=k-kini+1; j=0
        If (time>=(tref+traj_data%analysis%overlap_time%value)) Then
          k=kref+1; 
          kini=k
        Else
          kref=k
          kini=k+1
        End If
        tref=(kref)*traj_data%timestep%value
        j=0
        one_seg=.False.
        flag=.True.
      Else
        If (time>=(tref+traj_data%analysis%overlap_time%value) .And. flag) Then
          kref=k-1
          flag=.False.
        End If
        j=j+1
      End If
      k=k+1
    End Do

    If(one_seg) Then
      traj_data%analysis%N_seg=1
      traj_data%analysis%Ninterval=j
    Else
      traj_data%analysis%N_seg=i
      traj_data%analysis%Ninterval=l
    End If
    
    ! Allocate arrays
    Call traj_data%alloc_analysis()
    
    ! Calculate the number of segments
    If (traj_data%analysis%N_seg /= 1) Then
      i=0; j=0; l=0
      one_seg=.True.
      k=traj_data%analysis%frame_ini
      tref=tini; kini=k; flag=.True.
      Do While (k <= traj_data%frames)
        time=(k-1)*traj_data%timestep%value
        If (time>=(tref+traj_data%analysis%time_interval%value)) Then
          i=i+1
          l=k-kini+1; j=0
          traj_data%analysis%seg_indx(1,i)=kini
          traj_data%analysis%seg_indx(2,i)=k
          If (time>=(tref+traj_data%analysis%overlap_time%value)) Then
            k=kref+1; 
            kini=k
          Else
            kref=k
            kini=k+1
          End If
          tref=(kref)*traj_data%timestep%value
          k=kref+1; kini=k
          one_seg=.False.
          flag=.True.
        Else
          If (time>=(tref+traj_data%analysis%overlap_time%value) .And. flag) Then
            kref=k-1
            flag=.False.
          End If
          j=j+1
        End If
        k=k+1
      End Do
    Else
      traj_data%analysis%seg_indx(1,1)=traj_data%analysis%frame_ini
      traj_data%analysis%seg_indx(2,1)=traj_data%analysis%Ninterval+traj_data%analysis%frame_ini-1
    End If
    
  End Subroutine check_time_settings   
  
  Subroutine check_msd(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the &MSD block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%msd%select%fread) Then
      If (traj_data%msd%select%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "select" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%msd%select%type)/='x'  .And. &
            Trim(traj_data%msd%select%type)/='y'  .And. &
            Trim(traj_data%msd%select%type)/='z'  .And. &
            Trim(traj_data%msd%select%type)/='xy' .And. &
            Trim(traj_data%msd%select%type)/='xz' .And. &
            Trim(traj_data%msd%select%type)/='yz' .And. &
            Trim(traj_data%msd%select%type)/='xyz') Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "select". Valid options: "x", "y", "z", "xy",&
                                    & "xz", "yz" or "xyz"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "select" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "x", "y", "z", "xy", "xz", "yz" or "xyz"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If

    If (traj_data%msd%pbc_xyz%fread) Then
      If (traj_data%msd%pbc_xyz%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "pbc_xyz" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      traj_data%msd%pbc=.True.
    End If
    
  End Subroutine check_msd
  
  Subroutine check_ocf(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &OCF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the &OCF block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%ocf%legendre_order%fread) Then
      If (traj_data%ocf%legendre_order%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "legendre_order" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (traj_data%ocf%legendre_order%value < 1 .Or. traj_data%ocf%legendre_order%value>4) Then
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

    If (traj_data%ocf%u_definition%fread) Then
      If (traj_data%ocf%u_definition%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "u_definition" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%ocf%u_definition%type)/='bond_12'  .And. &
            Trim(traj_data%ocf%u_definition%type)/='bond_13'  .And. &
            Trim(traj_data%ocf%u_definition%type)/='bond_123' .And. &
            Trim(traj_data%ocf%u_definition%type)/='bond_12-13'  .And. &
            Trim(traj_data%ocf%u_definition%type)/='plane') Then
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
    
  End Subroutine check_ocf
  
  Subroutine print_trajectory_settings(traj_data, model_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print a summary of the trajectory settings
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data
    
    Character(Len=256) :: messages(4), word
    Integer(Kind=wi)   :: k, j

    Call info(' ', 1) 
    Call info('Trajectory settings', 1) 
    Call info('===================', 1) 

    Write (messages(1),'(1x,a)') '- by specification, the trajectory corresponds to the "'&
                                 &//Trim(traj_data%ensemble%type)//'" ensemble and was recorded in "'&
                                 &//Trim(model_data%input_geometry_format%type)//'" format'
    Call info(messages, 1)
    Write (messages(1),'(1x,a,f5.2,a)') '- the time step between recorded configurations is ', traj_data%timestep%value, ' fs'
    Call info(messages, 1)

    If (traj_data%ocf%invoke%fread .Or. traj_data%msd%invoke%fread & 
                                 & .Or. traj_data%lifetime%invoke%fread) Then
      If (traj_data%analysis%time_interval%fread) Then
        Write(word,'(f8.3)') traj_data%analysis%time_interval%value/1000.0_wp
        If (traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') '- the data analysis will be executed using time intervals of '&
                                     &//Trim(Adjustl(word))//' ps (this does not apply to RDF)' 
        Else
          Write (messages(1),'(1x,a)') '- the data analysis will be executed using time intervals of '&
                                     &//Trim(Adjustl(word))//' ps' 
        End If                             
        Call info(messages, 1)
      End If
      If (traj_data%analysis%ignore_initial%fread) Then
        Write(word,'(f8.3)') traj_data%analysis%ignore_initial%value/1000.0_wp
        Write (messages(1),'(1x,a)') '- the initial '//Trim(Adjustl(word))//' ps of the trajectory&
                                     & will be discarded' 
        Call info(messages, 1)
      End If
      If (traj_data%analysis%overlap_time%fread  .And. (traj_data%analysis%N_seg /=1)) Then
        Write(word,'(f8.3)') traj_data%analysis%overlap_time%value/1000.0_wp
        If (traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') '- the starting points of segments for analysis are separated&
                                       & by '//Trim(Adjustl(word))//' ps (this does not apply to RDF)'
        Else
          Write (messages(1),'(1x,a)') '- the starting points of segments for analysis are separated&
                                       & by '//Trim(Adjustl(word))//' ps'
        End If
        Call info(messages, 1)
      End If
    End If
    
    If (traj_data%ocf%invoke%fread) Then
      Call info(' ', 1)
      If (model_data%species_definition%atoms_per_species /= 1) Then
        If (model_data%species_definition%atoms_per_species == 2) Then
          If (Trim(traj_data%ocf%u_definition%type)/='bond_12') Then
            Write (messages(1),'(1x,a)')  '**WARNING: since the monitored species is diatomic, the&
                                          & method to compute the rotating unit vector (u_definition) &
                                          & has been reset to "bond_12". Methods "bond_13", "bond_12-13"&
                                          & "bond_123" and "plane" are meaningless for this case!'
            traj_data%ocf%u_definition%type='bond_12'                             
            Call info(messages, 1)
          End If
        Else
          If (Trim(traj_data%ocf%u_definition%type)/='bond_12-13' .And. & 
            Trim(traj_data%ocf%u_definition%type)/='bond_123') Then 
            Write (messages(1),'(1x,a)')  '**WARNING: the "'//Trim(traj_data%ocf%u_definition%type)//'"& 
                                           & option for the rotating unit vector (u_definition)&
                                           & is not recommended for most studies.'
            Write (messages(2),'(1x,a)')  '           Unless the user is fully certain, either the&
                                           & "bond_12-13" or "bond_123" option should be used instead.'
            Call info(messages, 2)
          End If
        End If
      End If
      Write (messages(1),'(1x,a)') 'The definition of the &OCF block will execute a Rotational&
                                  & Correlation funtion (OCF) analysis, using the species "'&
                                  &//Trim(model_data%species_definition%name%type)//&
                                  & '" defined in the &monitored_species block and the following settings:'
      Write (messages(2),'(1x,a)') '- the method to compute the attached rotating unit vector is: '//&
                                  & Trim(traj_data%ocf%u_definition%type)
      Write (messages(3),'(1x,a, i2)') '- the correlation terms are computed using a legendre polynomial of order ',&
                                  & traj_data%ocf%legendre_order%value
      Call info(messages, 3)
    End If

    If (traj_data%msd%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the &MSD block will execute a Mean Square&
                                  & Displacement analysis, using the species "'&
                                  &//Trim(model_data%species_definition%name%type)//&
                                  & '" defined in the &monitored_species block and the following settings:'
      Write (messages(2),'(1x,a)') '- the values will be computed for the coordinates(s): '//&
                                  & Trim(traj_data%msd%select%type)
      Call info(messages, 2)

      If (traj_data%msd%pbc_xyz%fread) Then
        Write (messages(1),'(1x,a)') '- the settings for the "pbc_xyz" directive uses (or not) periodic&
                                     & boundary condition for each coordinate' 
        Call info(messages, 1)                             
      End If
    End If

    If (traj_data%lifetime%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the &lifetime block will compute:'
      Write (messages(2),'(1x,a)') '- the transfer correlation function (TCF) for the changing species&
                                  & using the method (see manual): '//Trim(traj_data%lifetime%method%type)
      Write (messages(3),'(1x,a)') '- the residence times for each species separately (file RES_TIMES)'
      If (traj_data%lifetime%rattling_wait%fread) Then
        Write(word,'(f8.3)') traj_data%lifetime%rattling_wait%value/1000.0_wp
        Write (messages(4),'(1x,a)') 'Rattling lifetimes lower than '//Trim(Adjustl(word))//' ps will be discarded.'  
      Else
        Write (messages(4),'(1x,a)') 'Rattling effects are included in the calculation.'
      End If
      Call info(messages, 4)
    End If
    
    If (traj_data%rdf%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the &RDF block will compute the Radial&
                                  & Distribution Funtion (RDF) and the Coordination Numbers (CNs) using:'
      Write (messages(2),'(1x,a)') '- the tags defined in "tags_species_a"  and "tags_species_b"'
      Write (messages(3),'(1x,a,f6.3,a)') '- a discretization of ', traj_data%rdf%dr%value, ' Angstrom for the radius'
      Call info(messages, 3)
    End If

    If (traj_data%region%define%fread) Then
      Call info(' ', 1)
      If (traj_data%ocf%invoke%fread .Or. traj_data%msd%invoke%fread .Or. traj_data%rdf%invoke%fread) Then
        If (traj_data%ocf%invoke%fread .And. traj_data%msd%invoke%fread .And. traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'OCF, RDF and MSD analyses will be carried out for the region:'
        Else If ((.Not. traj_data%ocf%invoke%fread) .And. traj_data%msd%invoke%fread .And. traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'MSD and RDF analyses will be carried out for the region:'
        Else If (traj_data%ocf%invoke%fread .And. (.Not. traj_data%msd%invoke%fread) .And. traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'OCF and RDF analyses will be carried out for the region:'
        Else If (traj_data%ocf%invoke%fread .And. traj_data%msd%invoke%fread .And. (.Not. traj_data%rdf%invoke%fread)) Then
          Write (messages(1),'(1x,a)') 'OCF and MSD analyses will be carried out for the region:'
        Else If (traj_data%ocf%invoke%fread .And. (.Not. traj_data%msd%invoke%fread) .And. (.Not. traj_data%rdf%invoke%fread)) Then
          Write (messages(1),'(1x,a)') 'OCF analysis will be carried out for the region:'
        Else If ((.Not. traj_data%ocf%invoke%fread) .And. traj_data%msd%invoke%fread .And. (.Not. traj_data%rdf%invoke%fread)) Then
          Write (messages(1),'(1x,a)') 'MSD analysis will be carried out for the region:'
        Else If ((.Not. traj_data%ocf%invoke%fread) .And. (.Not. traj_data%msd%invoke%fread) .And. traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'RDF analysis will be carried out for the region:'
        End If
        Call info(messages, 1)
        Do k = 1, 3
          Do j = 1, traj_data%region%number(k)
            If (traj_data%region%invoke(k,j)%fread) Then
              Write (messages(1),'(3x,a,2f9.2)') '* '//Trim(traj_data%region%inout(k,j))//' the "'&
                                        //Trim(traj_data%region%invoke(k,j)%type)//'" region with&
                                        & lower and upper value: ', traj_data%region%domain(k,1,j), &
                                        & traj_data%region%domain(k,2,j)
              Call info(messages, 1)
            End If
          End Do 
        End Do
      End If 
    End If

  End Subroutine print_trajectory_settings

  Subroutine cross_checking(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the compatibility of model and trajectory data 
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(model_type),   Intent(In   ) :: model_data
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR: Problems in file '//Trim(files(FILE_SET)%filename)//' -'

    If (Trim(model_data%input_geometry_format%type)=='xyz') Then
      If (Trim(traj_data%ensemble%type)/='nve' .And. Trim(traj_data%ensemble%type)/='nvt') Then
         Call info(' ', 1)
         Write (messages(1),'(1x,a)') Trim(error_set)//' To date, trajectories in "xyz" format cannot be&
                                    & processed in the "'//Trim(traj_data%ensemble%type)//'" ensemble&
                                    & (this is part of future implementation).'  
         Call info(messages,1)           
         Call error_stop(' ')
      End If
    Else If (Trim(model_data%input_geometry_format%type)=='vasp') Then      
      If(model_data%config%simulation_cell%fread) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' Trajectories in "vasp" format contain the definition&
                                    & of the simulation cell within the file.'
         Write (messages(2),'(4x,a)') 'Thus, definition of the "&simulation_cell" block is not needed and can cause&
                                    & problems. Please remove/comment "&simulation_cell".'  
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If
   
    If (traj_data%ocf%invoke%fread) Then
      If (.Not. model_data%config%monitored_species%fread) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//'The computation of the orientational correlation&
                                     & function (OCF) requires the definition of the "&monitored_species" block'
         Call info(messages,1)           
         Call error_stop(' ')
      Else
        If (model_data%species_definition%atoms_per_species == 1) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The computation of the orientational correlation&
                                     & function (OCF) requires that the species defined in the&
                                     & "&monitored_species" block is a molecule. Please review the settings'
         Call info(messages,1)           
         Call error_stop(' ')
        End If
      End If
    End If
    
    If (traj_data%lifetime%invoke%fread) Then
      If (.Not. model_data%change_chemistry%stat) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The user has defined the &lifetime block but&
                                     & the &search_chemistry block is missing.'
         Write (messages(2),'(4x,a)') 'The computation of residence times and transfer correlation functions is&
                                     & only possible for systems with changing chemical species'
         Call info(messages, 2)           
         Call error_stop(' ')
      End If
    End If
   
  End Subroutine cross_checking

  Subroutine copy_to_trajectory(traj_data, model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to copy model arrays to each 
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data
    Integer(Kind=wi),   Intent(In   ) :: frame

    Integer(Kind=wi) :: l
    
    traj_data%config(frame,:)%tag=model_data%config%atom(:)%tag
    traj_data%config(frame,:)%element=model_data%config%atom(:)%element
    traj_data%box(frame)%cell=model_data%config%cell
    traj_data%box(frame)%invcell=model_data%config%invcell
    traj_data%box(frame)%volume=model_data%config%volume
    traj_data%box(frame)%cell_length=model_data%config%cell_length
    Do l = 1,3
      traj_data%config(frame,:)%r(l)=model_data%config%atom(:)%r(l)
    End Do
    
    ! Copy tracked species only if change_chemistry is set to True
    If(model_data%change_chemistry%stat) Then 
      Do l = 1, model_data%chem%N0%value
        traj_data%track_chem%config(frame,l)%r(:)=model_data%track_chem(l)%r(:)
        traj_data%track_chem%config(frame,l)%indx=model_data%track_chem(l)%indx
        traj_data%track_chem%config(frame,l)%tag=model_data%track_chem(l)%tag
      End Do
    End If

    ! Copy to species arrays
    If (model_data%config%monitored_species%fread) Then
      Do l = 1, model_data%config%Nmax_species
         traj_data%species(frame,l)%alive=model_data%config%species(l)%alive
         traj_data%species(frame,l)%list=model_data%config%species(l)%list
      End Do
    End If
      
  End Subroutine copy_to_trajectory

  Subroutine check_region_domain(model_data, traj_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks the definition of the &region block against the size of the 
    ! simulation cell
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),  Intent(InOut) :: traj_data
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: frame
  
    Character(Len=256) :: messages(2)  
    Integer(Kind=wi)   :: k, j
    Real(Kind=wp)      :: min_cell, max_cell, vector(3)
    Logical            :: flag1, flag2

    Write (messages(1),'(1x,a,i6)') '***ERROR: inconsistency between the size of the simulation cell and&
                                & the sepecifications of the &region block for frame: ', frame     

    Do k = 1, 3
      Do j = 1, traj_data%region%number(k)
        If (traj_data%region%invoke(k,j)%fread) Then
          vector(:)=model_data%config%cell(:,k)
          min_cell=Minval(vector)
          max_cell=Maxval(vector)
          If (traj_data%region%inside(k,j)) Then
            flag1 = (traj_data%region%domain(k,1,j) <= min_cell) .And.&
                    (traj_data%region%domain(k,2,j) <= min_cell)
            flag2 = (traj_data%region%domain(k,1,j) >= max_cell) .And.&
                    (traj_data%region%domain(k,2,j) >= max_cell)
            If (flag1 .Or. flag2) Then
               Write (messages(2),'(1x,a)') 'There are NO atoms inside the domain range defined for "'//&
                                          &Trim(traj_data%region%invoke(k,j)%type)//'". Please change' 
               Call info(messages,2)
               Call error_stop(' ')
            End If
          Else  
            flag1 = (traj_data%region%domain(k,1,j) <= min_cell) .And.&
                    (traj_data%region%domain(k,2,j) >= max_cell)
            If (flag1) Then
               Write (messages(2),'(1x,a)') 'There are NO atoms outside the domain range defined for "'//&
                                          &Trim(traj_data%region%invoke(k,j)%type)//'". Please change' 
               Call info(messages,2)
               Call error_stop(' ')
            End If
          End If
        End If
      End Do
    End Do
  
  End Subroutine check_region_domain
  
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
  
End Module trajectory
