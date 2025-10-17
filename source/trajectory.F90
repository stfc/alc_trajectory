!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to analyse the trajectory. If the directive change_chemistry
! is set to .True., the algorithm searches and tracks changes of 
! chemical species based on the information of the &search_chemistry
! block. 
!
! Copyright   2023-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Febđ 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module trajectory 
  Use atomic_model, Only : model_type, &
                           geo_param_type, &
                           about_cell, &
                           min_intra, &
                           atomistic_model, & 
                           check_definition_bonds, &
                           check_PBC,&
                           check_cell_consistency,&
                           check_length_directive, &
                           check_orthorhombic_cell,&
                           compute_distance_pbc,&
                           read_model,&
                           obtain_maximum_number_species,&
                           identify_monitored_indexes

  Use constants,    Only : max_components, &
                           max_at_species, &
                           max_unchanged_atoms, &
                           initial_tolerance, &
                           Rads_to_degrees, &
                           pi

  Use fileset,      Only : file_type, &
                           FILE_COORD_DISTRIB, &
                           FILE_INTERMOL_DISTANCES_NN1, &
                           FILE_INTERMOL_DISTANCES_NN2, &
                           FILE_INTERMOL_ANGLES, &
                           FILE_INTRAMOL_DISTANCES, &
                           FILE_INTRAMOL_ANGLES, &
                           FILE_SELECTED_NN_DISTANCES, &
                           FILE_MSD_ALL, &
                           FILE_MSD_AVG, &
                           FILE_OCF_ALL, &
                           FILE_OCF_AVG, &
                           FILE_CHEM_OCF_ALL, &
                           FILE_CHEM_OCF_AVG, &
                           FILE_SPCF_ALL, &
                           FILE_SPCF_AVG, &
                           FILE_RDF, &
                           FILE_RES_TIMES, &
                           FILE_SET, & 
                           FILE_TCF_ALL, &
                           FILE_TCF_AVG, &
                           FILE_TRAJECTORY,&
                           FILE_TRACK_CHEMISTRY,&
                           FILE_TAGGED_TRAJ, &
                           FILE_UNCHANGED_CHEM, & 
                           refresh_out

  Use input_types,  Only : in_integer, &
                           in_logic,   &
                           in_scalar,  &
                           in_param,   & 
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
    Integer(Kind=wi)  :: nn_indx(3)
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
     Type(in_logic)    :: print_all_intervals
  End Type 
  
  Type :: analysis_type
    Type(in_string)   :: invoke
    Type(in_logic)    :: normalise_at_t0
    Type(in_param)    :: time_interval
    Type(in_param)    :: end_time
    Type(in_param)    :: ignore_initial
    Type(in_param)    :: overlap_time
    Integer(Kind=wi)  :: N_seg
    Integer(Kind=wi)  :: Ninterval
    Integer(Kind=wi)  :: frame_ini
    Integer(Kind=wi)  :: frame_last
    Logical           :: normalised
    Integer(Kind=wi), Allocatable :: seg_indx(:,:)
    Real(Kind=wp),    Allocatable :: variable(:,:)
    Integer(Kind=wi), Allocatable :: max_points(:) 
  End Type
  
  ! Type to describe species
  Type :: species_type
    Integer(Kind=wi) :: list(max_at_species)
    Integer(Kind=wi) :: nn(2)
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
  
  !Type to describe the msd
  Type :: msd_type
    Type(in_string)  :: invoke
    Type(in_string)  :: pbc_xyz
    Type(in_string)  :: select
    Type(in_logic)   :: print_all_intervals
    Logical          :: pbc(3)
    Real(Kind=wp)    :: r2
  End Type

  !Type to describe the rdf
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
  
  !Type to describe the coordinate distribution
  Type :: coord_distrib_type
    Type(in_string)  :: invoke
    Type(in_string)  :: species_dir
    Character(Len=8) :: species
    Type(in_string)  :: coordinate
    Type(in_param)   :: delta
    Integer(Kind=wi) :: indx
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
    Type(in_string)  :: invoke
    Type(in_string)  :: method
    Type(in_param)   :: rattling_wait
    Type(in_logic)   :: print_all_intervals
  End Type
  
  !Type to describe the region where to constrain the analysis
  Type :: chemocf_type
    Type(in_string)  :: invoke
    Type(in_string)  :: variable
    Type(in_logic)   :: print_all_intervals
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
    Type(chemocf_type),   Public      :: chem_ocf
    Type(unchanged_type), Public      :: unchanged
    Type(analysis_type),  Public      :: analysis
    Type(track_type)                  :: track_chem              
    Type(in_param),       Public      :: timestep
    Type(in_logic),       Public      :: print_retagged_trajectory
    Integer(Kind=wi)                  :: frames
    Integer(Kind=wi)                  :: Nmax_species
    Integer(Kind=wi)                  :: N_species
    Logical                           :: reload_trajectory
    Logical                           :: active_bonds_computed
    Type(in_string),      Public      :: ensemble
    Type(coord_distrib_type), Public  :: coord_distrib
  Contains
    Private
      Procedure         :: alloc_trajectory  => allocate_trajectory_arrays
      Procedure         :: alloc_analysis    => allocate_analysis_arrays
      Final             :: cleanup
  End Type traj_type

  Public :: trajectory_analysis, extract_trajectory  
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
    
  Subroutine extract_trajectory(files, model_data, traj_data)
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
    Character(Len=256) :: message, nframes, md_length, net_md_length
    Character(Len=256) :: input_file, set_error
    Integer(Kind=wi)   :: i, j
    
    input_file=(files(FILE_TRAJECTORY)%filename)
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
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=input_file,Status='old')
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

    Call info(' ', 1) 
    Call info('Trajectory settings', 1) 
    Call info('===================', 1) 
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
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=(input_file),Status='old')
    Write(md_length,'(f12.3)') (traj_data%frames-1)*traj_data%timestep%value/1000.0_wp
    Write(nframes,'(i8)') traj_data%frames
    Call info(' ', 1)
    Call info('Start of the analysis', 1)
    Call info('=====================', 1)
    Write (message,'(1x,a)') 'The code has identified a total of '//Trim(Adjustl(nframes))//' frames. From&
                                 & the setting of the "timestep" directive, the recorded MD trajectory is '&
                                 &//Trim(Adjustl(md_length))//' ps long.'
    Call info(message, 1)
    
    If (traj_data%analysis%end_time%fread) Then
      If ((traj_data%analysis%end_time%value-(traj_data%frames-1)*traj_data%timestep%value)<0.0_wp) Then
        Write(net_md_length,'(f8.3)') traj_data%analysis%end_time%value/1000.0_wp
        Write (message,'(1x,a)') 'Nevertheless, from the set value of the "'//Trim(traj_data%analysis%end_time%tag)//&
                                    &'" directive (&data_analysis block), the analysis will consider up to '&
                                    &//Trim(Adjustl(net_md_length))//' ps of the MD trajectory.' 
        Call info(message, 1)
      End If
    End If    
         
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
    Call info(' The trajectory has been defined successfully!', 1)
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
      input_file=(files(FILE_TAGGED_TRAJ)%filename)
      Call print_tagged_trajectory(files, model_data, traj_data)
      Write (message,'(1x,a)') 'A copy of the trajectory with modified tags for the atomic species was printed&
                              & to the "'//Trim(input_file)//'" file'
      Call info(message, 1)
      Call refresh_out(files)
    End If 
  
  End Subroutine extract_trajectory
  
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
        Write(iunit,*) model_data%config%num_atoms, traj_data%frames, ' # number of total atoms and trajectory frames'
        Do l = 1, traj_data%frames
          Write(iunit,*) 'Frame=', l 
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
    Type(model_type),  Intent(In   ) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
  
    Character(Len=256)  :: message
    Logical             :: flag_inter_geo_stat

    traj_data%active_bonds_computed=.False.
    
    ! Compute the total average of monitored species within the system
    If (model_data%config%monitored_species%fread .And. model_data%species_definition%compute_amount%stat) Then
      Call compute_number_monitored_species(traj_data, model_data)
    End If

    If(model_data%change_chemistry%stat) Then 
      Call print_tracking_species(files, traj_data, model_data)
      ! Compute the lifetime and residence time for changing chemical species
      If (traj_data%lifetime%invoke%fread) Then
        Call compute_lifetime_related_quantities(files, traj_data, model_data)
      End If
      
      ! Compute OCF for the changing chemical species
      If (traj_data%chem_ocf%invoke%fread) Then
         If (.Not. traj_data%active_bonds_computed) Then
           Call find_active_bonds(traj_data, model_data)
           traj_data%active_bonds_computed=.True.
         End If
         Call compute_orientational_chemistry(files, traj_data, model_data)
      End If
    End If
    
    ! Compute OCF for monitored species
    If (traj_data%ocf%invoke%fread) Then
      Call orientational_correlation_function(files, traj_data)
    End If

    ! Compute MSD for monitores species
    If (traj_data%msd%invoke%fread) Then
      Call mean_squared_displacement(files, traj_data, model_data)
    End If

    ! Print coordinate distribution for the species under consideration
    If (traj_data%coord_distrib%invoke%fread) Then
      Call coordinate_distribution(files, traj_data, model_data)
    End If

    ! Compute RDF 
    If (traj_data%rdf%invoke%fread) Then
      Call radial_distribution_function(files, traj_data, model_data)
    End If

    ! Print coordinates for selected unchanged species along the MD trajectory
    If (traj_data%unchanged%invoke%fread) Then
      Call print_unchanged_chemistry(files, traj_data)
    End If

    ! Compute the distribution between the shortest distance of a selected pair
    If (model_data%nndist%invoke%fread) Then
      Call compute_shortest_distance_distribution(files, traj_data, model_data)
    End If
    
    ! Compute intramolecular properties for monitored species
    If (model_data%species_definition%intra_geom%invoke%fread) Then
      If (model_data%species_definition%intra_geom%dist%invoke%fread) Then
        Call obtain_intramol_geom_stat(files, traj_data, model_data%species_definition%atoms_per_species,&
                                    & model_data%species_definition%intra_geom%dist)
      End If
      If (model_data%species_definition%intra_geom%angle%invoke%fread) Then
        Call obtain_intramol_geom_stat(files, traj_data, model_data%species_definition%atoms_per_species,&
                                    & model_data%species_definition%intra_geom%angle)
      End If
    End If

    ! Compute intermolecular properties for monitored species
    If (model_data%species_definition%inter_geom%invoke%fread) Then
      Call find_neighbours_monitored_species(traj_data, model_data, flag_inter_geo_stat)
      If (model_data%species_definition%inter_geom%dist%invoke%fread .And. flag_inter_geo_stat)  Then
        Call obtain_intermol_geom_stat(files, traj_data, model_data%species_definition%inter_geom%dist, 1)
        Call obtain_intermol_geom_stat(files, traj_data, model_data%species_definition%inter_geom%dist, 2)
        Write (message,'(1x,a)') 'The probability distribution of the intermolecular distances were printed to files "'&
                                &//Trim(files(FILE_INTERMOL_DISTANCES_NN1)%filename)//'" and "'&
                                &//Trim(files(FILE_INTERMOL_DISTANCES_NN2)%filename)//'"'
        Call info(message, 1)
        Write (message,'(1x,a)') 'which separetely consider the first and the second nearest monitored species, respectively.'
        Call info(message, 1)
        Call info(' ', 1)
      End If
      If (model_data%species_definition%inter_geom%angle%invoke%fread .And. flag_inter_geo_stat) Then
        Call obtain_intermol_geom_stat(files, traj_data, model_data%species_definition%inter_geom%angle)
      End If
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

  Subroutine find_neighbours_monitored_species(traj_data, model_data, flag_exec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the two closest monitored species to a
    ! monitored species
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
    Logical,           Intent(  Out) :: flag_exec    

    Integer(Kind=wi)  :: accum, net_frames
    Integer(Kind=wi)  :: i, j, k, mk, indx1, indx2, mindx1, mindx2
    
    Character(Len=256) :: messages(3)
    Logical            :: flag, flag1, flag2

    Real(Kind=wp)  :: min_dist1, min_dist2, u(3)
    Logical        :: modified, finclude
    
    net_frames=0
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
      accum=0
      Do  j= 1, traj_data%Nmax_species
        min_dist1 = Huge(1.0_wp)
        min_dist2 = Huge(1.0_wp)
        If (traj_data%region%define%fread) Then
          mk=traj_data%species(i,j)%list(1)
          Call within_region(traj_data, i, mk, flag)
        Else
          flag=.True.
        End If

        If (traj_data%species(i,j)%alive .And. flag) Then
          accum=accum+1
          indx1=traj_data%species(i,j)%list(1)
          mindx1=indx1
          mindx2=indx1
          Do  k= 1, traj_data%Nmax_species
            If (model_data%species_definition%inter_geom%only_ref_tags_as_nn%stat) Then
              If (traj_data%species(i,k)%alive) Then
                finclude=.True.
              Else
                finclude=.False.
              End If
            Else
              finclude=.True.
            End If
            If (k /= j .And. finclude) Then  
              indx2= traj_data%species(i,k)%list(1)
              u= traj_data%config(i,indx2)%r-traj_data%config(i,indx1)%r
              Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              flag1= norm2(u) < min_dist1
              flag2= norm2(u) < min_dist2
              If (flag1 .And. flag2) Then
                min_dist2=min_dist1
                mindx2=mindx1
                min_dist1=norm2(u)
                mindx1=indx2
              Else If ((.Not. flag1) .And. flag2) Then
                min_dist2=norm2(u)
                mindx2=indx2
              End If
            End If  
          End Do
          traj_data%species(i,j)%nn(1)=mindx1
          traj_data%species(i,j)%nn(2)=mindx2
        End If
      End Do
      If (accum > 2 ) Then
        net_frames=net_frames+1
      End If
    End Do
    
    If (net_frames==0) Then
      Write (messages(1),'(1x,a)') '*************************************************************************************'
      Call info(messages, 1)
      Write (messages(1),'(1x,a)') '   WARNING: it looks the system has two or less monitored species along the trajectory (!?)'
      Write (messages(2),'(1x,a)') '   The intermolecular analysis of geometry parameters will not be executed.' 
      Write (messages(3),'(1x,a)') '   Please review the systems and the settings.'                        
      Call info(messages, 3)
      Write (messages(1),'(1x,a)') '************************************************************************************'
      Call info(messages, 1)
      flag_exec=.False.
    Else
      flag_exec=.True.
    End If
    
  End Subroutine find_neighbours_monitored_species 

  Subroutine compute_shortest_distance_distribution(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of the shortest distance between
    ! the "reference_species" and those defined by the "nn_species" directive
    ! in the distance domain defined by the lower_bound and upper_bound
    ! Analysis is performed using the definitions of the &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),        Intent(InOut) :: files(:)
    Type(traj_type),        Intent(InOut) :: traj_data
    Type(model_type),       Intent(In   ) :: model_data

    Integer(Kind=wi)  :: nbins, num_var, net_frames, accum
    Integer(Kind=wi)  :: fail(2) 

    Integer(Kind=wi)  :: i, j, k, mk
    Integer(Kind=wi)  :: num_at(2)
    Integer(Kind=wi)  :: list_indx(model_data%config%num_atoms,2)
    
    Integer(Kind=wi)  :: iunit
    
    Character(Len=256) :: messages(5), message
    Logical            :: falloc, flag, flag1, found

    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
 
    Real(Kind=wp)  :: rj(3), rk(3), rjk(3)
    Real(Kind=wp)  :: rmin
    Logical        :: modified
    
    ! Define number of bins
    nbins=Nint(Abs(model_data%nndist%upper_bound%value-model_data%nndist%lower_bound%value)/model_data%nndist%dr%value)
   
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for obtaining the statistics&
                                & of the shortest distances. Analysis will not be executed.'
      Call info(message, 1)                          
      falloc=.False.
    Else
      falloc=.True.
    End If
 
    If (falloc) Then
      d=0.0_wp
      ! Initiate Accumulators
      accum=0
      net_frames=0

      ! Compute the histogram for the selected coordinate of the selected species
      Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
        h=0
        num_var=0
        num_at=0
        list_indx=0
        Do j = 1, model_data%config%num_atoms
          If (model_data%nndist%reference_species==traj_data%config(i,j)%tag) Then
              num_at(1)=num_at(1)+1
              list_indx(num_at(1),1)=j
          End If 
          k=1
          flag=.False.
          Do While (k <= model_data%nndist%num_nn_species .And. (.Not. flag))
            If (model_data%nndist%nn_species(k)==traj_data%config(i,j)%tag) Then
              num_at(2)=num_at(2)+1
              list_indx(num_at(2),2)=j
              flag=.True.  
            End If
            k=k+1
          End Do
        End Do
        
        If (num_at(1) /= 0 .And. num_at(2) /=0) Then
          Do j = 1, num_at(1)
            rj=traj_data%config(i,list_indx(j,1))%r
            rmin=Huge(1.0_wp)
            found=.False.
            Do k= 1, num_at(2)
              If (list_indx(k,2) /= list_indx(j,1)) Then
                rk=traj_data%config(i,list_indx(k,2))%r
                If (traj_data%region%define%fread) Then
                  Call within_region(traj_data, i, list_indx(j,1), flag1)
                  If (flag1) Then
                    flag=.True.
                  Else
                    flag=.False. 
                  End If
                Else
                   flag=.True.
                End If   
                If (flag) Then
                  rjk=rj-rk
                  Call check_PBC(rjk, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                  If (norm2(rjk) < rmin) Then
                    found=.True.
                    rmin=norm2(rjk)
                  End If
                End If
              End If
            End Do
            If (found) Then
              If (rmin >= model_data%nndist%lower_bound%value .And. rmin <= model_data%nndist%upper_bound%value) Then
                mk=Floor((rmin-model_data%nndist%lower_bound%value)/model_data%nndist%dr%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            End If
          End Do
        End If
        
        If(num_var /= 0) Then
          accum=accum+num_var
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do mk= 1, nbins 
            d(mk)= d(mk)+Real(h(mk),Kind=wp)/num_var
          End Do
        End If
        
      End Do
      
      ! Print results
      If (accum /= 0) Then
         Do mk=1, nbins 
           d(mk)=d(mk)/net_frames/model_data%nndist%dr%value
         End Do
       
        ! Print File
        Open(Newunit=files(FILE_SELECTED_NN_DISTANCES)%unit_no, File=files(FILE_SELECTED_NN_DISTANCES)%filename,&
                          &Status='Replace')
        iunit=files(FILE_SELECTED_NN_DISTANCES)%unit_no
        Write (iunit,'(a)') '#  Probability distribution of the shortest distances between'
        Write (iunit,'(a)') '#  the reference species "'//Trim(model_data%nndist%reference_species)//&
                             &'" and the species defined in the "nn_species" directive'    
        Write (iunit,'(a)') '#  Distance [Angstrom]     Probability [1/Angstrom]' 
        Do mk=1, nbins
          Write(iunit,'(2x,f12.4,6x,f14.5)') (Real(mk,Kind=wp)-0.5_wp)*model_data%nndist%dr%value+&
                                            & model_data%nndist%lower_bound%value, d(mk)
        End Do
        Write (message,'(1x,a)') 'The probability distribution of shortest distances between the "reference_species"&
                                  & and the species defined in the "nn_species" directive was printed to the "'&
                                  &//Trim(files(FILE_SELECTED_NN_DISTANCES)%filename)//'" file.'
        Call info(message, 1)

      Else
        Write (messages(1),'(1x,a)')   '*************************************************************************************'
        Write (messages(2),'(1x,a)')   '   WARNING: the statistics for the shortest distances between "reference_species" and'
        Write (messages(3),'(1x,a)')   '   "nn_species" could not be executed. File "'&
                                         &//Trim(files(FILE_SELECTED_NN_DISTANCES)%filename)//'" was not generated.'
        Write (messages(4),'(1x,a)')   '   Please check the settings of the "&selected_nn_distances" block.       '
        Write (messages(5),'(1x,a)')   '************************************************************************************'
        Call info(messages, 5)
       End If

      Deallocate(d,h)
    End If
    
    Call refresh_out(files)
    
  End Subroutine compute_shortest_distance_distribution

  Subroutine obtain_intermol_geom_stat(files, traj_data, M, num_nn)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of geometrical
    ! parameters (distance and angles) between three closest monitored
    ! species, with the criteria defined in the &intermol_stat_settings
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),            Intent(InOut)  :: files(:)
    Type(traj_type),            Intent(InOut)  :: traj_data
    Type(geo_param_type),       Intent(In   )  :: M
    Integer(Kind=wi), Optional, Intent(In   )  :: num_nn

    Integer(Kind=wi)  :: nbins, num_var, net_frames, accum
    Integer(Kind=wi)  :: fail(2) 

    Integer(Kind=wi)  :: i, j, k, k1, k2, mk
    Integer(Kind=wi)  :: iunit
    
    Character(Len=256) :: messages(2), message
    Logical            :: falloc, flag

    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
 
    Real(Kind=wp)  :: u(3), u2(3), angle 
    Logical        :: modified
    
    ! Define number of bins
    nbins=Nint(Abs(M%upper_bound%value-M%lower_bound%value)/M%delta%value)
    
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for obtaining geometry statistics&
                                & of monitored species. Analysis will not be executed.'
      Call info(message, 1)                          
      falloc=.False.
    Else
      falloc=.True.
    End If
 
    If (falloc) Then
      d=0.0_wp
      ! Initiate Accumulators
      accum=0
      net_frames=0
      
      ! Compute the histogram for the selected coordinate of the selected species
      Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
        h=0
        num_var=0
        Do  j= 1, traj_data%Nmax_species
          If (traj_data%region%define%fread) Then
            mk=traj_data%species(i,j)%list(1)
            Call within_region(traj_data, i, mk, flag)
          Else
            flag=.True.
          End If

          If (traj_data%species(i,j)%alive .And. flag) Then
            If (Trim(M%invoke%type) == '&distance_parameters') Then
              k=traj_data%species(i,j)%list(1)
              k2=traj_data%species(i,j)%nn(num_nn)
              u=traj_data%config(i,k2)%r-traj_data%config(i,k)%r          
              Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              If (norm2(u) >= M%lower_bound%value .And. norm2(u) <= M%upper_bound%value) Then
                mk=Floor((norm2(u)-M%lower_bound%value)/M%delta%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            Else If (Trim(M%invoke%type) == '&angle_parameters') Then
              k=traj_data%species(i,j)%list(1)
              k1=traj_data%species(i,j)%nn(1)
              k2=traj_data%species(i,j)%nn(2)
              u =traj_data%config(i,k1)%r-traj_data%config(i,k)%r
              u2=traj_data%config(i,k2)%r-traj_data%config(i,k)%r
              Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              Call check_PBC(u2, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
              angle=Acos((u(1)*u2(1)+u(2)*u2(2)+u(3)*u2(3))/norm2(u)/norm2(u2))*Rads_to_degrees
              If (angle >= M%lower_bound%value .And. angle <= M%upper_bound%value) Then
                mk=Floor((angle-M%lower_bound%value)/M%delta%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            End If
          End If
        End Do
        
        If(num_var /= 0) Then
          accum=accum+num_var
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do mk= 1, nbins 
            d(mk)= d(mk)+Real(h(mk),Kind=wp)/num_var
          End Do
        End If
      End Do

      ! Print results
      If (accum /= 0) Then
         Do mk=1, nbins 
           d(mk)=d(mk)/net_frames/M%delta%value
         End Do
       
        ! Print File
        If (Trim(M%invoke%type) == '&distance_parameters') Then
          If (num_nn == 1) Then
            Open(Newunit=files(FILE_INTERMOL_DISTANCES_NN1)%unit_no, File=files(FILE_INTERMOL_DISTANCES_NN1)%filename,&
                              &Status='Replace')
            iunit=files(FILE_INTERMOL_DISTANCES_NN1)%unit_no
            Write (iunit,'(a)') '#  Probability distribution of the intermolecular distances&
                               & using only the first nearest monitored species and the settings of '//Trim(M%invoke%type)  
            Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
          Else If (num_nn == 2) Then
            Open(Newunit=files(FILE_INTERMOL_DISTANCES_NN2)%unit_no, File=files(FILE_INTERMOL_DISTANCES_NN2)%filename,&
                              &Status='Replace')
            iunit=files(FILE_INTERMOL_DISTANCES_NN2)%unit_no
            Write (iunit,'(a)') '#  Probability distribution of the intermolecular distances&
                               & using only the second nearest monitored species and the settings of '//Trim(M%invoke%type)  
            Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
          End If
        Else If (Trim(M%invoke%type) == '&angle_parameters') Then
          Open(Newunit=files(FILE_INTERMOL_ANGLES)%unit_no, File=files(FILE_INTERMOL_ANGLES)%filename, Status='Replace')
          iunit=files(FILE_INTERMOL_ANGLES)%unit_no
          Write (iunit,'(a)') '#  Probability distribution of the intermolecular angles using the settings of '//Trim(M%invoke%type)
          Write (iunit,'(a)') '#  Value [Degrees]      Probability [1/Degrees]' 
          Write (message,'(1x,a)') 'The probability distribution of the intermolecular angles was printed to the "'&
                                  &//Trim(files(FILE_INTERMOL_ANGLES)%filename)//'" file.'
          Call info(message, 1)
          Call info(' ', 1)
        End If
          Do mk=1, nbins
            Write(iunit,'(2x,f12.4,6x,f14.5)') (Real(mk,Kind=wp)-0.5)*M%delta%value+M%lower_bound%value, d(mk)
          End Do
      Else
        Write (messages(1),'(1x,a)')   '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)')   '   WARNING: the statistics for the requested intermolecular geometry could not be executed'
        Write (messages(2),'(1x,a)') '   Please verify the settings for the '//Trim(M%invoke%type)//' in &intermol_stat_settings'
        Call info(messages, 2)
        Write (messages(1),'(1x,a)')   '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Deallocate arrays   
      Deallocate(d,h)
    End If
    
    Call refresh_out(files)
    
  End Subroutine obtain_intermol_geom_stat

  
  Subroutine obtain_intramol_geom_stat(files, traj_data, atoms_per_species, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of internal geometrical
    ! parameters (distance and angles) for monitored species 
    ! as defined in the &intramol_stat_settings
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut)  :: files(:)
    Type(traj_type),      Intent(InOut)  :: traj_data
    Integer(Kind=wi),     Intent(In   )  :: atoms_per_species
    Type(geo_param_type), Intent(In   )  :: M

    Integer(Kind=wi)  :: nbins, num_var, net_frames, accum
    Integer(Kind=wi)  :: fail(2) 

    Integer(Kind=wi)  :: i, j, k1, k2, k3, mk, l, l1, l2
    Integer(Kind=wi)  :: ni(3), nj(2)
    Integer(Kind=wi)  :: iunit
    
    Character(Len=256) :: messages(2), message
    Logical           :: falloc, flag, flag1, flag2

    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
 
    Real(Kind=wp)  :: u(3), u2(3), angle 
    Logical        :: modified
    
    ! Define number of bins
    nbins=Nint(Abs(M%upper_bound%value-M%lower_bound%value)/M%delta%value)
    
    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for obtaining geometry statistics&
                                & of monitored species. Analysis will not be executed.'
      Call info(message, 1)                          
      falloc=.False.
    Else
      falloc=.True.
    End If
 
    If (falloc) Then
      d=0.0_wp
      ! Initiate Accumulators
      accum=0
      net_frames=0
      
      ! Compute the histogram for the selected coordinate of the selected species
      Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
        h=0
        num_var=0
        Do  j= 1, traj_data%Nmax_species
          If (traj_data%region%define%fread) Then
            mk=traj_data%species(i,j)%list(1)
            Call within_region(traj_data, i, mk, flag)
          Else
            flag=.True.
          End If

          If (traj_data%species(i,j)%alive .And. flag) Then
            If (Trim(M%invoke%type) == '&distance_parameters') Then 
              Do k1= 1, atoms_per_species
                ni(1)=traj_data%species(i,j)%list(k1)
                Do k2= k1+1, atoms_per_species
                  ni(2)=traj_data%species(i,j)%list(k2)
                  flag1=(traj_data%config(i,ni(1))%element==M%species(1)) .And.&
                        (traj_data%config(i,ni(2))%element==M%species(2))
                  flag2=(traj_data%config(i,ni(1))%element==M%species(2)) .And.&
                        (traj_data%config(i,ni(2))%element==M%species(1))      
                  If (flag1 .Or. flag2) Then
                    u=traj_data%config(i,ni(1))%r-traj_data%config(i,ni(2))%r          
                    Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                    If (norm2(u) >= M%lower_bound%value .And. norm2(u) <= M%upper_bound%value) Then
                      mk=Floor((norm2(u)-M%lower_bound%value)/M%delta%value)+1
                      If (mk <= nbins) Then
                        h(mk)=h(mk)+1
                        num_var=num_var+1
                      End If
                    End If
                  End If
                End Do
              End Do
              
            Else If (Trim(M%invoke%type) == '&angle_parameters') Then
              Do k1= 1, atoms_per_species
                ni(1)=traj_data%species(i,j)%list(k1)
                Do k2= k1+1, atoms_per_species
                  ni(2)=traj_data%species(i,j)%list(k2)
                  Do k3= k2+1, atoms_per_species
                    ni(3)=traj_data%species(i,j)%list(k3)
                    Do l= 1, 3
                      If (traj_data%config(i,ni(l))%element==M%species(2)) Then
                        Do l1= 1, atoms_per_species
                          nj(1)=traj_data%species(i,j)%list(l1)
                          Do l2= l1+1, atoms_per_species
                            nj(2)=traj_data%species(i,j)%list(l2)
                            If (l1 /= l .And. l2 /= l) Then
                              flag1=(traj_data%config(i,nj(1))%element==M%species(1)) .And.&
                                    (traj_data%config(i,nj(2))%element==M%species(3))
                              flag2=(traj_data%config(i,nj(1))%element==M%species(3)) .And.&
                                    (traj_data%config(i,nj(2))%element==M%species(1))
                              If (flag1 .Or. flag2) Then
                                u =traj_data%config(i,nj(1))%r-traj_data%config(i,ni(l))%r
                                u2=traj_data%config(i,nj(2))%r-traj_data%config(i,ni(l))%r
                                Call check_PBC(u, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                                Call check_PBC(u2, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                                angle=Acos((u(1)*u2(1)+u(2)*u2(2)+u(3)*u2(3))/norm2(u)/norm2(u2))*Rads_to_degrees
                                If (angle >= M%lower_bound%value .And. angle <= M%upper_bound%value) Then
                                  mk=Floor((angle-M%lower_bound%value)/M%delta%value)+1
                                  If (mk <= nbins) Then
                                    h(mk)=h(mk)+1
                                    num_var=num_var+1
                                  End If
                                End If
                              End If
                            End If
                          End Do
                        End Do
                      End If
                    End Do  
                  End Do
                End Do
              End Do
            End If
          End If
        End Do

        If(num_var /= 0) Then
          accum=accum+num_var
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do mk= 1, nbins 
            d(mk)= d(mk)+Real(h(mk),Kind=wp)/num_var
          End Do
        End If
      End Do

      ! Print results
      If (accum /= 0) Then
        Do mk=1, nbins 
          d(mk)=d(mk)/net_frames/M%delta%value
        End Do
      
        ! Print File
        If (Trim(M%invoke%type) == '&distance_parameters') Then
          Open(Newunit=files(FILE_INTRAMOL_DISTANCES)%unit_no, File=files(FILE_INTRAMOL_DISTANCES)%filename, Status='Replace')
          iunit=files(FILE_INTRAMOL_DISTANCES)%unit_no
          Write (iunit,'(a)') '#  Probability distribution of the intramolecular distances&
                             & using the settings of '//Trim(M%invoke%type)  
          Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
          Write (message,'(1x,a)') 'The probability distribution of the intramolecular distances was printed to the "'&
                                  &//Trim(files(FILE_INTRAMOL_DISTANCES)%filename)//'" file.'
          Call info(message, 1)
          Call info(' ', 1)
        Else If (Trim(M%invoke%type) == '&angle_parameters') Then
          Open(Newunit=files(FILE_INTRAMOL_ANGLES)%unit_no, File=files(FILE_INTRAMOL_ANGLES)%filename, Status='Replace')
          iunit=files(FILE_INTRAMOL_ANGLES)%unit_no
          Write (iunit,'(a)') '#  Probability distribution of the intramolecular angles using the settings of '//Trim(M%invoke%type)
          Write (iunit,'(a)') '#  Value [Degrees]      Probability [1/Degrees]' 
          Write (message,'(1x,a)') 'The probability distribution of the intramolecular angles was printed to the "'&
                                  &//Trim(files(FILE_INTRAMOL_ANGLES)%filename)//'" file.'
          Call info(message, 1)
          Call info(' ', 1)
        End If
          Do mk=1, nbins
            Write(iunit,'(2x,f12.4,6x,f13.5)') (Real(mk,Kind=wp)-0.5)*M%delta%value+M%lower_bound%value, d(mk)
          End Do
      Else
        Write (messages(1),'(1x,a)')   '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)')   '   WARNING: the statistics for the requested intramolecular geometry could not be executed'
        Write (messages(2),'(1x,a)') '   Please verify the settings for the '//Trim(M%invoke%type)//' in &intramol_stat_settings'
        Call info(messages, 2)
        Write (messages(1),'(1x,a)')   '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Deallocate arrays   
      Deallocate(d,h)
    End If

    Call refresh_out(files)
    
  End Subroutine obtain_intramol_geom_stat
  
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
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
      l=l+1
      k=0
      Do j = 1, model_data%chem%N0%value
        Do m = 1, model_data%chem%acceptor%N0_incl
          word=Trim(model_data%chem%acceptor%tg_incl(m))//'*'
          If (traj_data%track_chem%config(i,j)%tag==word) Then
            counts(m)=counts(m)+1
            k=k+1 
          End If
        End Do
      End Do
    End Do

    Write (messages(1),'(1x,a)') 'Population probabilities of donor species along MD trajectory'
    Write (messages(2),'(1x,a)') '(species defined in the "include_tags" directive of the "&acceptor_criteria" block)'
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
    Type(model_type),  Intent(In   ) :: model_data
  
    Integer(Kind=wi)   :: iunit, i, l 
    Character(Len=256) :: num_species
    Character(Len=256) :: message
  
    ! Print tracked species
    Open(Newunit=files(FILE_TRACK_CHEMISTRY)%unit_no, File=files(FILE_TRACK_CHEMISTRY)%filename, Status='Replace')
    iunit=files(FILE_TRACK_CHEMISTRY)%unit_no
    If (traj_data%analysis%frame_ini==1) Then
      Write(iunit,'(a)') '# Tracking the change of chemical species over the whole trajectory'    
    Else
      Write(iunit,'(a,1x,f10.4,1x,a)') '# Tracking the change of chemical species ignoring the first',& 
                                   &  traj_data%analysis%frame_ini*traj_data%timestep%value/1000_wp,&
                                   & 'ps of the whole trajectory. This value is set to time zero below.'
    End If

    If (model_data%chem%N0%value==1) Then
      Write (iunit,'(a,9x,a)') '# Time (ps)', 'XYZ_Species_1' 
    Else
      Write(num_species,*) model_data%chem%N0%value
      Write (iunit,'(a,9x,2a)') '# Time (ps)', 'XYZ_Species_1 .... XYZ_Species_', Adjustl(Trim(num_species)) 
    End If
    
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
       Write(iunit,'(f10.4, 4x, *(f11.3))') (i-traj_data%analysis%frame_ini)*traj_data%timestep%value/1000_wp,&
                                       & (traj_data%track_chem%config(i,l)%r(:), l=1, model_data%chem%N0%value)
    End Do
    Write (message,'(1x,a)') 'The tracking of the changing chemical species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_TRACK_CHEMISTRY)%filename)//'" file'
    Call info(message, 1)
    Call refresh_out(files)
    Close(iunit)
    Call info(' ', 1)
  
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
    Character(Len=256) :: frame
    Character(Len=256) :: num_species, message, messages(3), spec, num1, num2
    Logical            :: flag
  
    flag=.True. 
  
    Open(Newunit=files(FILE_UNCHANGED_CHEM)%unit_no, File=files(FILE_UNCHANGED_CHEM)%filename, Status='Replace')
    iunit=files(FILE_UNCHANGED_CHEM)%unit_no
    
    If (traj_data%analysis%frame_ini==1) Then
      Write(iunit,'(a)') '# Tracking unchanged chemical species over the whole trajectory'    
    Else  
      Write(iunit,'(a,1x,f10.4,1x,a)') '# Tracking the unchanged chemical species ignoring the first',& 
                                   &  traj_data%analysis%frame_ini*traj_data%timestep%value/1000_wp,&
                                   & 'ps of the whole trajectory. This value is set to time zero below.' 
    End If
    
    If(traj_data%unchanged%N0==1) Then
      Write(iunit,'(a)') '# The label and number for the species is consistent with the settings&
                                   & of the "&track_unchanged_chemistry" block.'
    Else                               
      Write(iunit,'(a)') '# The species labelling, ordering and numbering is consistent with the settings&
                                   & of the "&track_unchanged_chemistry" block.'    
    End If
    
    spec=Trim(traj_data%unchanged%tag%type)
    Write (num1,*) traj_data%unchanged%indexes(1)
    If(traj_data%unchanged%N0==1) Then
      Write (iunit,'(a,5x,a)') '#  Time (ps)', 'XYZ_'//Trim(spec)//'_'//Adjustl(Trim(num1))
    Else
      Write (num2,*) traj_data%unchanged%indexes(traj_data%unchanged%N0)
      Write (iunit,'(a,5x,a)') '#  Time (ps)', 'XYZ_'//Trim(spec)//'_'//Adjustl(Trim(num1))//&
                              &'.... XYZ_'//Trim(spec)//'_'//Adjustl(Trim(num2))
    End If 
    
    i=traj_data%analysis%frame_ini
    Do While (i <= traj_data%analysis%frame_last .And. flag)
      l =1
      Do While (l<= traj_data%unchanged%N0 .And. flag)
        k=traj_data%unchanged%indexes(l)
        If (Trim(traj_data%config(i,k)%tag) /= Trim(traj_data%unchanged%tag%type)) Then
          flag=.False.
        End If
        l=l+1
      End Do
      If (flag) Then
        Write(iunit,'(f10.4, 1x, *(f11.3))') (i-traj_data%analysis%frame_ini)*traj_data%timestep%value/1000.0_wp,&
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
      Write (message,'(1x,a)') 'The tracking of the selected, unchanged chemical species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_UNCHANGED_CHEM)%filename)//'" file'
      Call info(message, 1)
    End If
    
    Close(iunit)
    Call refresh_out(files)

  End Subroutine print_unchanged_chemistry

  Subroutine compute_number_monitored_species(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the average number (and STD) of the monitored species 
    !
    ! author    - i.scivetti June 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

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
    
    Allocate(nat(traj_data%analysis%frame_last-traj_data%analysis%frame_ini+1))
    
    ! Compute the histogram for atoms of type a and b
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
      ! Define the number and list of indexes for type of species "a"
      num_at_a=0
      net_frames=net_frames+1 
      Do j = 1, model_data%config%num_atoms
        If (model_data%species_definition%reference_tag%type==traj_data%config(i,j)%tag) Then
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
      Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
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

  Subroutine coordinate_distribution(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the distribution of the coordinates (x, y or z) of
    ! the species selected in the &coord_distrib block
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)  :: i, j, m, iunit, indx_ini
    Integer(Kind=wi)  :: num_at, nbins, net_frames
    Integer(Kind=wi)  :: accum
    
    Real(Kind=wp)     :: clim(2), coord_value
    Real(Kind=wp)     :: vector(3) 

    Integer(Kind=wi)  :: list_indx(model_data%config%num_atoms)
    
    Character(Len=256) :: ctap
    
    Character(Len=256) :: messages(3), message
    Logical            :: falloc

    Character(Len=256) :: type_error
    Integer(Kind=wi)   :: fail(2)  
    
    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: d(:)
    
    ! Search for the value of cmax and cmin 
    indx_ini=traj_data%analysis%frame_ini
    vector=0.0_wp
    Do i = 1, 3
      vector(:)=vector(:)+traj_data%box(indx_ini)%cell(i,:)
    End Do
    If (vector(traj_data%coord_distrib%indx) > 0.0_wp) Then
      clim(2)=vector(traj_data%coord_distrib%indx)
      clim(1)=0.0_wp
      ctap='top'
    Else
      clim(1)=vector(traj_data%coord_distrib%indx) 
      clim(2)=0.0_wp 
      ctap='bottom'
    End If
    
    Do i = traj_data%analysis%frame_ini+1, traj_data%analysis%frame_last
      vector=0.0_wp
      Do j = 1, 3
        vector(:)=vector(:)+traj_data%box(i)%cell(j,:)
      End Do
      If (Trim(ctap)=='bottom') Then
        If (vector(traj_data%coord_distrib%indx) < clim(1)) Then
          clim(1)=vector(traj_data%coord_distrib%indx)
        End If
      Else If (Trim(ctap)=='top') Then
        If (vector(traj_data%coord_distrib%indx) > clim(2)) Then
          clim(2)=vector(traj_data%coord_distrib%indx)
        End If
      End If
    End Do
     
    ! Define number of bins
    nbins=Floor(Abs(clim(1)-clim(2))/traj_data%coord_distrib%delta%value)

    !Allocate arrays
    Allocate(h(nbins),  Stat=fail(1))
    Allocate(d(nbins),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for coordinate distribution.&
                                & Analysis will not be executed.'
      falloc=.False.
    Else
      falloc=.True.
    End If
     
    If (falloc) Then
      d=0.0_wp
      ! Initiate Accumulators
      accum=0
      net_frames=0
      
      ! Compute the histogram for the selected coordinate of the selected species
      Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
        ! Define the number and list of indexes
        num_at=0
        list_indx=0
        Do j = 1, model_data%config%num_atoms
          If (traj_data%coord_distrib%species==traj_data%config(i,j)%tag) Then
            num_at=num_at+1
            list_indx(num_at)=j
          End If
        End Do
      
        ! Calculate the histogram for this particular frame of the trajectory
        If (num_at/=0) Then
          h=0
          Do j=1, num_at
            coord_value=Abs(traj_data%config(i,list_indx(j))%r(traj_data%coord_distrib%indx))
            If (Trim(ctap)=='top') Then
              m=Floor(coord_value/traj_data%coord_distrib%delta%value)+1
            Else
              m=nbins+1-(Floor(coord_value/traj_data%coord_distrib%delta%value)+1)
            End If
            If (m <= nbins) Then
              h(m)=h(m)+1
            End If
          End Do 
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do m=1, nbins 
            d(m)= d(m)+Real(h(m),Kind=wp)/num_at
          End Do
        End If
        accum=accum+num_at
      End Do

      Do m=1, nbins 
        d(m)=d(m)/net_frames/traj_data%coord_distrib%delta%value      
      End Do
      
      ! Print results
      If (accum /= 0) Then
        ! Print File
        Open(Newunit=files(FILE_COORD_DISTRIB)%unit_no, File=files(FILE_COORD_DISTRIB)%filename, Status='Replace')
        iunit=files(FILE_COORD_DISTRIB)%unit_no
        Write (iunit,'(a)') '#  Distribution of the '//Trim(traj_data%coord_distrib%coordinate%type)//'-coordinate&
                           & for the "'//Trim(traj_data%coord_distrib%species)//'" species'
        Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
        
        Do m=1, nbins
          Write(iunit,'(2x,f12.4,6x,f14.5)') (Real(m,Kind=wp)-0.5)*traj_data%coord_distrib%delta%value, d(m)
        End Do
        Write (message,'(1x,a)') 'The distribution of the '//Trim(traj_data%coord_distrib%coordinate%type)//&
                                &'-coordinate for the "'//Trim(traj_data%coord_distrib%species)//'" species was&
                                & printed to the "'//Trim(files(FILE_COORD_DISTRIB)%filename)//'" file.'
        Call info(message, 1)
      Else
        type_error=Trim(traj_data%coord_distrib%species)
        Write (messages(1),'(1x,a)') '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: coordinate distribution analysis could not be executed'
          Write (messages(2),'(1x,a)') '   Requested species '//Trim(type_error)//' as specified in the &coord_distrib&
                                  & block could not be identified along the trajectory.'
          Write (messages(3),'(1x,a)') '   Please verify the settings for the &coord_distrib block. The user should also&
                                & look at the file '//Trim(files(FILE_TAGGED_TRAJ)%filename)                        
        Call info(messages, 3)
        Write (messages(1),'(1x,a)') '************************************************************************************'
        Call info(messages, 1)
      End If
      
      ! Close file
      Close(iunit)
      ! Deallocate arrays   
      Deallocate(d,h)
   End If
    
   Call refresh_out(files) 
    
  End Subroutine coordinate_distribution


  Subroutine radial_distribution_function(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the Radial Distribution Function (RDF) based on the
    ! settings of the &RDF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)  :: i, j, k, m, iunit, indx_a, indx_b, itype
    Integer(Kind=wi)  :: num_at_a, num_at_b, nbins, net_frames
    Integer(Kind=wi)  :: accum_a, accum_b 
    
    Real(Kind=wp)     :: rmax, r_bin, rho_b, dV
    Real(Kind=wp)     :: rj(3), rk(3), rjk(3) 
    
    Integer(Kind=wi)  :: list_indx_a(model_data%config%num_atoms)
    Integer(Kind=wi)  :: list_indx_b(model_data%config%num_atoms)
    
    Character(Len=256) :: messages(3), message
    Character(Len=256) :: type_error
    Logical            :: modified, falloc, flag, ftype
    Logical            :: counted(model_data%config%num_atoms)
    Integer(Kind=wi)   :: fail(4)  
   
    Integer(Kind=wi), Allocatable  :: h(:)
    Real(Kind=wp),    Allocatable  :: gr(:)
    Real(Kind=wp),    Allocatable  :: nn(:)
    Real(Kind=wp),    Allocatable  :: cn(:)
    
    ! Search for the value of rmax 
    rmax=-Huge(1.0_wp)
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
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
      Call info(message, 1)  
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
      Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
        ! Define the number and list of indexes for type of species "a"
        num_at_a=0
        list_indx_a=0
        Do j = 1, model_data%config%num_atoms
          itype=1
          ftype=.True.
          Do While (itype <= traj_data%rdf%num_type_a .And. ftype)
            If (traj_data%rdf%type_a(itype)==traj_data%config(i,j)%tag) Then
              ftype=.False.      
              If (traj_data%region%define%fread) Then
                 Call within_region(traj_data, i, j, flag)
              End If
              If (flag) Then
                num_at_a=num_at_a+1
                list_indx_a(num_at_a)=j
              End If
            End If
            itype=itype+1
          End Do
        End Do
      
        ! Define the number and list of indexes for type of species "b"
        num_at_b=0
        list_indx_b=0
        Do j = 1, model_data%config%num_atoms
          itype=1
          ftype=.True.
          Do While (itype <= traj_data%rdf%num_type_b .And. ftype)
            If (traj_data%rdf%type_b(itype)==traj_data%config(i,j)%tag) Then
              ftype=.False.      
              num_at_b=num_at_b+1
              list_indx_b(num_at_b)=j
            End If
            itype=itype+1
          End Do
        End Do
        
        ! Accummulators
        accum_b=accum_b+num_at_b
        accum_a=accum_a+num_at_a
      
        !Define rho_b 
        rho_b= num_at_b/(traj_data%box(i)%volume)  
        
        ! Calculate the histogram for this particular frame of the trajectory
        If (num_at_a /=0 .And. num_at_b/=0) Then
          h=0
          counted=.False.
          Do j=1, num_at_a 
            indx_a=list_indx_a(j)
            rj=traj_data%config(i,indx_a)%r
            Do k=1, num_at_b
              indx_b=list_indx_b(k)
              If (indx_a /= indx_b) Then 
                rk=traj_data%config(i,indx_b)%r
                rjk=rj-rk
                Call check_PBC(rjk, traj_data%box(i)%cell, traj_data%box(i)%invcell, 0.5_wp, modified)
                m=Floor(norm2(rjk)/traj_data%rdf%dr%value)+1
                If (m <= nbins) Then
                  h(m)=h(m)+1
                End If
              End If
            End Do
            counted(indx_a)=.True.
          End Do 
          ! Count net frame
          net_frames=net_frames+1
          ! Normalise
          Do m=1, nbins 
            gr(m)= gr(m)+Real(h(m),Kind=wp)/(num_at_a*rho_b)
            nn(m)= nn(m)+Real(h(m),Kind=wp)/(num_at_a)
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
          gr(m)=gr(m)/(dV*net_frames)
          Do k= 1, m
           cn(m)=cn(m)+nn(k)
          End Do
          cn(m)=cn(m)/net_frames 
          !Write(iunit,'(2x, f11.3,(2(6x,f11.6)))') (Real(m,Kind=wp)-0.5)*traj_data%rdf%dr%value, gr(m), cn(m)
          Write(iunit,'(2x, f11.3,(2(6x,f11.6)))') (1.0_wp*m-0.5)*traj_data%rdf%dr%value, gr(m), cn(m)
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
    
    Call refresh_out(files)
    
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
    Logical           :: set_u0

    Character(Len=256) :: message, num_species
    Logical           :: terminated(traj_data%Nmax_species)

    If (traj_data%msd%print_all_intervals%stat) Then
      ! Print tracked species
      Open(Newunit=files(FILE_MSD_ALL)%unit_no, File=files(FILE_MSD_ALL)%filename, Status='Replace')
      iunit=files(FILE_MSD_ALL)%unit_no
      Write(num_species,*) model_data%chem%N0%value
      Write (iunit,'(a)') '#  MSD analysis for all the intervals' 
      Write (iunit,'(a,8x,a)') '#  Time (ps)', '"'//Trim(traj_data%msd%select%type)//'"-MSD for species "'&
                                &//Trim(model_data%species_definition%name%type)//'" [Angstrom^2]' 
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
          End Do
        Else
          Nini_species=0
          Do j = 1, traj_data%Nmax_species
            If (traj_data%species(i,j)%alive) Then
              indx=traj_data%species(i,j)%list(1)
              traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx)%r
              traj_data%species(i,j)%u0(:,1)=traj_data%species(i,j)%u(:,1)
              Nini_species=Nini_species+1
            Else
              terminated(j)=.True. 
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
            If (traj_data%msd%print_all_intervals%stat) Then
              Write(iunit,'(f11.3, 4x, f11.4)') (time-base_time)/1000.0_wp, traj_data%msd%r2
            End If  
          End If  
          terminated=.False.
          set_u0=.True.
          If (traj_data%msd%print_all_intervals%stat) Then
            If ((traj_data%analysis%N_seg /=1) .And. (k /= traj_data%analysis%N_seg)) Then
              If (k /= traj_data%analysis%N_seg) Then
                Write (iunit,'(a,8x,a)') '#  Time (ps)', '"'//Trim(traj_data%msd%select%type)//'"-MSD for species "'&
                                  &//Trim(model_data%species_definition%name%type)//'" [Angstrom^2]' 
              End If                    
            End If
          End If
        Else
          If ((traj_data%N_species) /= 0) Then
            traj_data%msd%r2=traj_data%msd%r2/traj_data%N_species
          Else  
            traj_data%msd%r2=0.0_wp
          End If
            traj_data%analysis%variable(l,k)=traj_data%msd%r2
            If (traj_data%msd%print_all_intervals%stat) Then
              Write(iunit,'(f11.3, 4x, f11.4)') (time-base_time)/1000.0_wp, traj_data%msd%r2
            End If
        End If  
      End Do
    End Do

    If (traj_data%msd%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The MSD analysis for the multiple time intervals was printed to the "'&
                                 &//Trim(files(FILE_MSD_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The MSD analysis was printed to the "'//Trim(files(FILE_MSD_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time interval.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If

    Call average_segments(files, traj_data, FILE_MSD_AVG, 'MSD')
    If (traj_data%analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time interval was used to compute the average MSD! The computed STD&
                              & is zero. Use/Check the &data_analysis block to improve the statistics.'
      Call info(message, 1)
    End If

    If (.Not. traj_data%msd%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the MSD analysis for all time intervals,&
                                & the "print_all_intervals" directive (within the &msd block) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%analysis%N_seg ==1) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_MSD_ALL)%filename)//'" and "'//Trim(files(FILE_MSD_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If
    
    Call info(' ', 1)
    Call refresh_out(files)
    
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
  
  Subroutine compute_lifetime_related_quantities(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the transfer correlation functions and residence
    ! times
    !
    ! author    - i.scivetti April 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
    
    Call transfer_correlation_function_sites(files, traj_data, model_data)
    Call residence_times_sites(files, traj_data, model_data)

    If (.Not. traj_data%active_bonds_computed) Then
      Call find_active_bonds(traj_data, model_data)
      traj_data%active_bonds_computed=.True.
    End If
    
    Call special_pair_correlation_function(files, traj_data, model_data)
  
  
  End Subroutine compute_lifetime_related_quantities

  Subroutine residence_times_sites(files, traj_data, model_data)
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
      Do While (i <= traj_data%analysis%frame_last)
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
            If (((time-tchange) > rattling .Or. i==traj_data%analysis%frame_last)) Then
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
    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine residence_times_sites
  
  Subroutine compute_orientational_chemistry(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the orientational chemistry along the trajectory
    !
    ! author    - i.scivetti Febraury 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

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
 
    If (traj_data%chem_ocf%print_all_intervals%stat) Then
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
        If (Trim(traj_data%chem_ocf%variable%type)=='special_pair') Then
        ! Obtain the index of the closest acceptor
          Do m = 1, model_data%chem%N0%value
            sites(1,l,m)=traj_data%track_chem%config(i,m)%indx
            sites(2,l,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
            sites(3,l,m)=traj_data%track_chem%config(i,m)%nn_indx(2)
            sites(4,l,m)=traj_data%track_chem%config(i,m)%nn_indx(3)
          End Do
        Else If (Trim(traj_data%chem_ocf%variable%type)=='acceptor_donor_transfer_couple') Then
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

          If (Trim(traj_data%chem_ocf%variable%type)=='special_pair') Then
            Call orientational_correlation_term_transfer_couple(traj_data, i, s1, u(:,m), u0(:,m), suma_i, ncouples)
          Else If (Trim(traj_data%chem_ocf%variable%type)=='acceptor_donor_transfer_couple') Then
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
            If (traj_data%chem_ocf%print_all_intervals%stat) Then
              Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
            End If
          End If
          If (traj_data%chem_ocf%print_all_intervals%stat) Then
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
          If (traj_data%chem_ocf%print_all_intervals%stat) Then
            Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If  
        End If
        
      End Do
      
    End Do  

    If (traj_data%chem_ocf%print_all_intervals%stat) Then
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
    
    If (.Not. traj_data%chem_ocf%print_all_intervals%stat) Then
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

  Subroutine compute_closest_pairs(traj_data, model_data, frame, nchem, s1, s2, s3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the closest possible acceptor
    !
    ! author    - i.scivetti Feb 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(In   ) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: frame
    Integer(Kind=wi),  Intent(In   ) :: nchem
    Integer(Kind=wi),  Intent(Out  ) :: s1
    Integer(Kind=wi),  Intent(Out  ) :: s2
    Integer(Kind=wi),  Intent(Out  ) :: s3
  
    Integer(Kind=wi)   :: i, j, k, l, mindx(3)
    Real(Kind=wp)      :: dist, min_dist(3)
    Logical            :: match_j, fexcl, flag(3)
    Character(Len=8)   :: tgexcl 

    i=traj_data%track_chem%config(frame,nchem)%indx
    min_dist(1)=Huge(1.0_wp) 
    min_dist(2)=Huge(1.0_wp)
    min_dist(3)=Huge(1.0_wp)
    
    If (model_data%chem%acceptor%info_exclude%fread) Then
      tgexcl=traj_data%config(frame,i)%tag  
      Call remove_symbols(tgexcl, '*')
      fexcl=.False.
      l=1
      Do While (l <= model_data%chem%acceptor%N0_excl .And. (.Not. fexcl))
        If (tgexcl==model_data%chem%acceptor%tg_excl(l)) Then
           fexcl=.True.
        End If
        l=l+1
      End Do  
    Else 
      fexcl=.False. 
    End If

    mindx(1)=i
    mindx(2)=i
    mindx(3)=i
    
    j=1
    Do While (j <= model_data%config%num_atoms)
      If (i/=j) Then
        match_j=.False.
        k=1
        Do While (k <= model_data%chem%acceptor%N0_incl .And. (.Not. match_j))
          If (traj_data%config(frame,j)%tag==model_data%chem%acceptor%tg_incl(k)) Then
            match_j=.True.
            If (fexcl) Then
              l=1
              Do While (l <= model_data%chem%acceptor%N0_excl .And. match_j)
                If (traj_data%config(frame,j)%tag==model_data%chem%acceptor%tg_excl(l)) Then
                   match_j=.False.
                End If
                l=l+1
              End Do
            End If
          End If
          k=k+1
        End Do

        If(match_j) Then
          Call compute_distance_PBC(traj_data%config(frame,i)%r, traj_data%config(frame,j)%r,&
                                  & traj_data%box(frame)%cell, traj_data%box(frame)%invcell, dist)
          flag(1)= dist < min_dist(1)
          flag(2)= dist < min_dist(2)
          flag(3)= dist < min_dist(3)
          If (flag(1) .And. flag(2) .And. flag(3)) Then
            min_dist(3)=min_dist(2)
            min_dist(2)=min_dist(1)
            min_dist(1)=dist
            mindx(3)=mindx(2)
            mindx(2)=mindx(1)
            mindx(1)=j
          Else If ((.Not. flag(1)) .And. flag(2) .And. flag(3)) Then
            min_dist(3)=min_dist(2)
            min_dist(2)=dist
            mindx(3)=mindx(2)
            mindx(2)=j
          Else If ((.Not. flag(1)) .And. (.Not. flag(2)) .And. flag(3)) Then
            min_dist(3)=dist
            mindx(3)=j
          End If
        End If
      End If
      j=j+1
    End Do
    
    s1=mindx(1)
    s2=mindx(2)
    s3=mindx(3)
    
    If (s1==s2) Then
      call error_stop('ERRROOR')
    End If
    
    If (s1==s3) Then
      call error_stop('ERRROOR')
    End If
    
    If (s2==s3) Then
      call error_stop('ERRROOR')
    End If
  
  End Subroutine compute_closest_pairs

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
 
  Subroutine find_active_bonds(traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to identify the active bond for the changing sites along the
    ! trajectory
    !
    ! author    - i.scivetti April 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, m
    Integer(Kind=wi)   :: s1, s2, s3
    
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
      Do  m= 1, model_data%chem%N0%value
        Call compute_closest_pairs(traj_data, model_data, i, m, s1, s2, s3)
        traj_data%track_chem%config(i,m)%nn_indx(1)=s1
        traj_data%track_chem%config(i,m)%nn_indx(2)=s2
        traj_data%track_chem%config(i,m)%nn_indx(3)=s3
      End Do
    End Do
  
  End Subroutine find_active_bonds
 
  Subroutine transfer_correlation_function_sites(files, traj_data, model_data)
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
    Logical            :: set_t0
    Character(Len=256) :: message
    
    Logical            :: terminated(traj_data%analysis%Ninterval, model_data%chem%N0%value)
    Integer(Kind=wi)   :: indexes(2,model_data%chem%N0%value)
    Integer(Kind=wi)   :: ref_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: time_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: Nnet
    
    Logical            :: follow(model_data%chem%N0%value)
    Logical            :: flag, first(model_data%chem%N0%value)
    Logical            :: hold(model_data%chem%N0%value)

    Real(Kind=wp)      :: tchange(model_data%chem%N0%value)
    Character(Len=256) :: method
    
    method=Trim(traj_data%lifetime%method%type)
    
    If (traj_data%lifetime%print_all_intervals%stat) Then
     ! Print header
     Open(Newunit=files(FILE_TCF_ALL)%unit_no, File=files(FILE_TCF_ALL)%filename, Status='Replace')
     iunit=files(FILE_TCF_ALL)%unit_no
     Write (iunit,'(a)') '#  TCF Analysis for all the intervals' 
     Write (iunit,'(a)') '#  Time (ps)         TCF' 
    End If

    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1

    Do k= 1, traj_data%analysis%N_seg
      set_t0=.True.
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
          If (set_t0) Then
            Do m = 1, model_data%chem%N0%value
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
            End Do  
            set_t0=.False.
          Else
            Do m = 1, model_data%chem%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%indx /= ref_indx(m)) Then
                  If (.Not. hold(m)) Then
                    tchange(m)=time
                    hold(m)=.True.
                    time_indx(m)=i
                  End If
                Else 
                  hold(m)=.False.
                End If  
                If (hold(m)) Then
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
              End If
            End Do
          End If
        End If
        
        If (Trim(method)=='hicf*' .Or. Trim(method)=='hdcf*') Then
          If (set_t0) Then
            Do m = 1, model_data%chem%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%indx
              ref_indx(m)=traj_data%track_chem%config(i,m)%indx
            End Do  
            set_t0=.False.
          Else
    
            Do m = 1, model_data%chem%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%indx/=indexes(2,m)) Then
                  indexes(1,m)=indexes(2,m)
                  indexes(2,m)=traj_data%track_chem%config(i,m)%indx
                End If
                
                If (indexes(1,m) /= ref_indx(m) .And. indexes(2,m) /= ref_indx(m)) Then
                  If (.Not. hold(m)) Then
                   tchange(m)=time
                   hold(m)=.True.
                   time_indx(m)=i
                  End If
                Else
                  hold(m)=.False.
                End If
    
                If (hold(m)) Then
                  If (Trim(method)=='hicf*') Then
                    ntop=i
                  Else If (Trim(method)=='hdcf*') Then
                    ntop=traj_data%analysis%seg_indx(2,k)
                    follow(m)=.False.
                  End If
                  
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                  hold(m)=.False.
                End If
              End If
            End Do
          End If
        End If
    
      End Do
      
      l=0
      base_time=(traj_data%analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        time=(i-1)*traj_data%timestep%value
        l=l+1
        Nnet=0
        Do j = 1, model_data%chem%N0%value
          If (traj_data%region%define%fread) Then
            Call within_region(traj_data, i, traj_data%track_chem%config(i,j)%indx, flag)
          Else
            flag=.True.
          End If
          If (flag) Then
            Nnet=Nnet+1
          End If
        End Do
    
        suma_i=0.0_wp
        If (.Not. All(terminated(l,:))) Then
          suma_i=0.0_wp
          Do j = 1, model_data%chem%N0%value
            If(.Not. terminated(l,j)) Then
              If (traj_data%region%define%fread) Then
                Call within_region(traj_data, i, traj_data%track_chem%config(i,j)%indx, flag)
              Else
                flag=.True.
              End If
              If (flag) Then
                suma_i=suma_i+1.0_wp    
              End If
            End If  
          End Do
          If (Nnet > 0) Then
            suma_i=suma_i/Nnet
          End If
        End If
       
        traj_data%analysis%variable(l,k)=suma_i
        If (traj_data%lifetime%print_all_intervals%stat) Then
          Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          If (i==traj_data%analysis%seg_indx(2,k) .And. (traj_data%analysis%N_seg /=1)) Then
             If (k /= traj_data%analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)        TCF' 
             End If
          End If   
        End If
      End Do
    End Do
    
    If (traj_data%lifetime%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The TCF analysis for the multiple time intervals was printed to the "'&
                                 &//Trim(files(FILE_TCF_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The TCF analysis was printed to the "'//Trim(files(FILE_TCF_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time interval.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    
    ! Compute average
    Call average_segments(files, traj_data, FILE_TCF_AVG, 'TCF')
    If (traj_data%analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time interval was used to compute the average TCF! The computed STD&
                              & is zero. Use/Check the &data_analysis block to improve the statistics.'
      Call info(message, 1)
    End If
    
    If (.Not. traj_data%lifetime%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the TCF analysis for all time intervals,&
                                & the "print_all_intervals" directive (within the &lifetime block) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%analysis%N_seg ==1 .And. (.Not. traj_data%analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_TCF_ALL)%filename)//'" and "'//Trim(files(FILE_TCF_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If

    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine transfer_correlation_function_sites
  
  Subroutine special_pair_correlation_function(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the correlation function of the special pair (SPCF),
    ! which is defined by the atomic pair of the chemical site and the closest NN, 
    ! either donor or acceptor
    !
    ! author    - i.scivetti April 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data

    Integer(Kind=wi)   :: i, j, k, l, m, n, ntop
    Integer(Kind=wi)   :: iunit, ini_indx
    Real(Kind=wp)      :: suma_i 
    Real(Kind=wp)      :: time
    Real(Kind=wp)      :: base_time
    Logical            :: set_t0
    Character(Len=256) :: message
    
    Logical            :: terminated(traj_data%analysis%Ninterval, model_data%chem%N0%value)
    Integer(Kind=wi)   :: indexes(2,model_data%chem%N0%value)
    Integer(Kind=wi)   :: ref_indx(2,model_data%chem%N0%value)
    
    Integer(Kind=wi)   :: time_indx(model_data%chem%N0%value)
    Integer(Kind=wi)   :: Nnet
    
    Logical            :: follow(model_data%chem%N0%value)
    Logical            :: first(model_data%chem%N0%value)
    Logical            :: flag
    Logical            :: hold(model_data%chem%N0%value)

    Real(Kind=wp)      :: tchange(model_data%chem%N0%value)
    Character(Len=256) :: method
    
    method=traj_data%lifetime%method%type
    
    If (traj_data%lifetime%print_all_intervals%stat) Then
     ! Print header
     Open(Newunit=files(FILE_SPCF_ALL)%unit_no, File=files(FILE_SPCF_ALL)%filename, Status='Replace')
     iunit=files(FILE_SPCF_ALL)%unit_no
     Write (iunit,'(a)') '#  SPCF Analysis for all the intervals' 
     Write (iunit,'(a)') '#  Time (ps)         SPCF' 
    End If

    !Set max_points to beyond the interval
    traj_data%analysis%max_points=traj_data%analysis%Ninterval+1

    Do k= 1, traj_data%analysis%N_seg
      set_t0=.True.
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
          If (set_t0) Then
            Do m = 1, model_data%chem%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
              ref_indx(1,m)=indexes(1,m)
              ref_indx(2,m)=indexes(2,m)
            End Do  
            set_t0=.False.
          Else
    
            Do m = 1, model_data%chem%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%nn_indx(1)/=indexes(2,m)) Then
                   indexes(1,m)=traj_data%track_chem%config(i,m)%indx
                   indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
                End If
                
                If ((indexes(1,m) /= ref_indx(1,m) .Or. indexes(2,m) /= ref_indx(2,m)) .And. &
                    (indexes(1,m) /= ref_indx(2,m) .Or. indexes(2,m) /= ref_indx(1,m))) Then
                  If (.Not. hold(m)) Then
                   tchange(m)=time
                   hold(m)=.True.
                   time_indx(m)=i
                  End If
                Else
                  hold(m)=.False.
                End If
    
                If (hold(m)) Then
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
              End If
            End Do
          End If
        End If

        If (Trim(method)=='hicf*' .Or. Trim(method)=='hdcf*') Then
          If (set_t0) Then
            Do m = 1, model_data%chem%N0%value
              indexes(1,m)=traj_data%track_chem%config(i,m)%indx
              indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
              ref_indx(1,m)=indexes(1,m)
            End Do  
            set_t0=.False.
          Else
    
            Do m = 1, model_data%chem%N0%value
              If (follow(m)) Then
                If (traj_data%track_chem%config(i,m)%nn_indx(1)/=indexes(2,m)) Then
                   indexes(1,m)=traj_data%track_chem%config(i,m)%indx
                   indexes(2,m)=traj_data%track_chem%config(i,m)%nn_indx(1)
                End If
                
                If (indexes(1,m) /= ref_indx(1,m) .And. indexes(2,m) /= ref_indx(1,m)) Then
                  If (.Not. hold(m)) Then
                   tchange(m)=time
                   hold(m)=.True.
                   time_indx(m)=i
                  End If
                Else
                  hold(m)=.False.
                End If
    
                If (hold(m)) Then
                  If (Trim(method)=='hicf*') Then
                    ntop=i
                  Else If (Trim(method)=='hdcf*') Then
                    ntop=traj_data%analysis%seg_indx(2,k)
                    follow(m)=.False.
                  End If
                  
                  Do n=time_indx(m), ntop
                    terminated(n-ini_indx+1,m)=.True.
                  End Do
                  hold(m)=.False.
                End If
              End If
            End Do
          End If
        End If
        
      End Do
      
      l=0
      base_time=(traj_data%analysis%seg_indx(1,k)-1)*traj_data%timestep%value
      Do i = traj_data%analysis%seg_indx(1,k), traj_data%analysis%seg_indx(2,k)
        time=(i-1)*traj_data%timestep%value
        l=l+1
        Nnet=0
        Do j = 1, model_data%chem%N0%value
          If (traj_data%region%define%fread) Then
            Call within_region(traj_data, i, traj_data%track_chem%config(i,j)%indx, flag)
          Else
            flag=.True.
          End If
          If (flag) Then
            Nnet=Nnet+1
          End If
        End Do
    
        suma_i=0.0_wp
        If (.Not. All(terminated(l,:))) Then
          suma_i=0.0_wp
          Do j = 1, model_data%chem%N0%value
            If(.Not. terminated(l,j)) Then
              If (traj_data%region%define%fread) Then
                Call within_region(traj_data, i, traj_data%track_chem%config(i,j)%indx, flag)
              Else
                flag=.True.
              End If
              If (flag) Then
                suma_i=suma_i+1.0_wp    
              End If
            End If  
          End Do
          If (Nnet > 0) Then
            suma_i=suma_i/Nnet
          End If
        End If
       
        traj_data%analysis%variable(l,k)=suma_i
        If (traj_data%lifetime%print_all_intervals%stat) Then
          Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          If (i==traj_data%analysis%seg_indx(2,k) .And. (traj_data%analysis%N_seg /=1)) Then
             If (k /= traj_data%analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)        SPCF' 
             End If
          End If   
        End If
      End Do
    End Do
    
    If (traj_data%lifetime%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'The SPCF analysis for the multiple time intervals was printed to the "'&
                                 &//Trim(files(FILE_SPCF_ALL)%filename)//'" file.'
      Else
        Write (message,'(1x,a)') 'The SPCF analysis was printed to the "'//Trim(files(FILE_SPCF_ALL)%filename)//'" file&
                                 & and corresponds to a single (only one) time interval.'
      End If
      Call info(message, 1)
      Close(iunit)
    End If
    
    ! Compute average
    Call average_segments(files, traj_data, FILE_SPCF_AVG, 'SPCF')
    If (traj_data%analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time interval was used to compute the average SPCF! The computed STD&
                              & is zero. Use/Check the &data_analysis block to improve the statistics.'
      Call info(message, 1)
    End If
    
    If (.Not. traj_data%lifetime%print_all_intervals%stat) Then
      If (traj_data%analysis%N_seg /=1 ) Then 
        Write (message,'(1x,a)') 'In case the user wants to print the SPCF analysis for all time intervals,&
                                & the "print_all_intervals" directive (within the &lifetime block) must be set to .True.'
        Call info(message, 1)
      End If
    Else
      If (traj_data%analysis%N_seg ==1 .And. (.Not. traj_data%analysis%normalised)) Then
        Write (message,'(1x,a)') 'WARNING: Files "'&
                               &//Trim(files(FILE_SPCF_ALL)%filename)//'" and "'//Trim(files(FILE_SPCF_AVG)%filename)//&
                               &'" contain redundant results.'
        Call info(message, 1)
      End If
    End If

    Call info(' ', 1)
    Call refresh_out(files)
    
  End Subroutine special_pair_correlation_function
  
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
    Real(Kind=wp)    :: sum_i, average, average_net, std, norma
    Logical          :: flag, fcompute
    Character(Len=256) :: message
    Character(Len=256) :: quantity

    If (Trim(what) == 'CHEM_OCF') Then
      quantity='OCF'
    Else
      quantity=what
    End If
    
    fcompute=.True.
    traj_data%analysis%normalised=.False.
    
    ! Check the value at t=0
    If (Trim(what) /= 'MSD') Then
      sum_i=0.0_wp
      Nnet=0
      Do k= 1, traj_data%analysis%N_seg
        Nnet=Nnet+1
        sum_i=sum_i+traj_data%analysis%variable(1,k)
      End Do
      
      If (Nnet > 0) Then
        norma=sum_i/Nnet
        If (Abs(norma-1.0_wp)>initial_tolerance)then
          If (Abs(norma)< initial_tolerance) Then
            Write (message,'(1x,a)') '*** WARNING: Problems with the computation of the average '//Trim(what)
            Call info(message, 1)
            Write (message,'(1x,a)') '             This is likely due to poor statistics'
            Call info(message, 1)
            If (traj_data%region%define%fread) Then
              Write (message,'(1x,a)') '             Please check the settings: the &region block might be too small.'
            Else
              Write (message,'(1x,a)') '             Please check the settings'
            End If
            Call info(message, 1)
            Call info('***', 1)
          Else
            If (traj_data%analysis%normalise_at_t0%stat) Then
              Write (message,'(1x,a)') '*** INFO: The average '//Trim(what)//' has been normalised at t=0.'
              Call info(message, 1)
              traj_data%analysis%variable=traj_data%analysis%variable/norma
              traj_data%analysis%normalised=.True.
            Else
              Write (message,'(1x,a)') '*** WARNING: The average '//Trim(what)//' is NOT normalised at t=0.'
              Call info(message, 1)
              If (traj_data%region%define%fread) Then
                Write (message,'(1x,a)') '             Please check the settings: the region defined in the &region&
                                         & block for analysis might be too small.'
              End If
              
              If (traj_data%analysis%invoke%fread) Then
                Write (message,'(1x,a)') '    To normalise, set the "normalise_at_t0" directive to .True. in the&
                                        & &data_analysis block.'
              Else 
                Write (message,'(1x,a)') '    To normalise, use the &data_analysis block and set the "normalise_at_t0"&
                                        & directive to .True.'
              End If
              Call info(message, 1)                       
              Call info(' ***', 1) 
            End If
          End If
        End If
      Else
        Write (message,'(1x,a)') '**** PROBLEMS: The average '//Trim(what)//' could not be computed.'
        Call info(message, 1)
        Write (message,'(1x,a)') '             This is likely due to poor statistics'
        Call info(message, 1)
         If (traj_data%region%define%fread) Then
           Write (message,'(1x,a)') '             Please check the settings: the region defined in the &region&
                                   & block for analysis might be too small.'
         End If
        Call info(message, 1)
        Call info(' ***', 1)
      End If
    End If
    
    If (fcompute) Then             
      ! Print header
      Open(Newunit=files(file_number)%unit_no, File=files(file_number)%filename, Status='Replace')
      iunit=files(file_number)%unit_no
      If (Trim(what) == 'MSD') Then
        Write (iunit,'(a)') '#  Average MSD and STD (in Angstrom^2) for the coordinate(s) "'//&
                        &Trim(traj_data%msd%select%type)//'" of the monitored species'
      End If
      Write (iunit,'(a)') '#  Time (ps)      '//Trim(quantity)//'          STD' 
      
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
          If (Trim(what)=='MSD') Then
            average_net=average
          Else
            If(average > 1.0_wp) Then
              average_net=1.0_wp
            Else
              average_net=average
            End If
          End If  
          Write(iunit,'(3(f10.3, 3x))') (i-1)*traj_data%timestep%value/1000.0_wp, average_net, std
        End If
        i=i+1
        
      End Do
      Write (message,'(1x,a)') 'The average '//Trim(what)//' was printed to the "'//&
                               &Trim(files(file_number)%filename)//'" file.'
      Call info(message, 1)
      Close(iunit)
    End If 
    
    
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
    Logical           :: set_u0
    Character(Len=256) :: message
    
    Logical           :: terminated(traj_data%Nmax_species)

    If (traj_data%ocf%print_all_intervals%stat) Then
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
            If (Trim(traj_data%ocf%u_definition%type) == 'bond_12-13') Then
              traj_data%species(i,:)%u0(j,2)=traj_data%species(i-1,:)%u0(j,2)
            End If
          End Do
        Else
          Nini_species=0
          Do j = 1, traj_data%Nmax_species
            If (traj_data%species(i,j)%alive) Then
              Call rotation_vector_monitored_species(traj_data, i, j)
              traj_data%species(i,j)%u0(:,1)=traj_data%species(i,j)%u(:,1)
              If (Trim(traj_data%ocf%u_definition%type) == 'bond_12-13') Then
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
          If(.Not. terminated(j)) Then
            If (traj_data%species(i,j)%alive) Then
              Call rotation_vector_monitored_species(traj_data, i, j)
              Call orientational_correlation_term_monitored_species(traj_data, i, j, suma_i)  
            Else
              !terminated(j)=.True.
            End If
          End If  
        End Do
    
        If (i==traj_data%analysis%seg_indx(2,k)) Then
          If (traj_data%N_species /= 0) Then
            suma_i=suma_i/traj_data%N_species
            traj_data%analysis%variable(l,k)=suma_i
            If (traj_data%ocf%print_all_intervals%stat) Then
              Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
            End If
          End If  
          terminated=.False.
          set_u0=.True.
          If (traj_data%ocf%print_all_intervals%stat) Then
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
          If (traj_data%ocf%print_all_intervals%stat) Then
            Write(iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          End If
        End If
      End Do
    End Do
    
    If (traj_data%ocf%print_all_intervals%stat) Then
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
    
    If (.Not. traj_data%ocf%print_all_intervals%stat) Then
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

  Subroutine rotation_vector_monitored_species(traj_data, i, j)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the rotation vector of the monitored species
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
    
  End Subroutine rotation_vector_monitored_species  

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
  
  Subroutine orientational_correlation_term_monitored_species(traj_data, i, j, suma_i)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation for the relevant 
    ! species j at the MD frame i (cij)
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

  End Subroutine orientational_correlation_term_monitored_species  
  
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

    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=files(FILE_TRAJECTORY)%filename,Status='old')

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
    Call check_time_directive(traj_data%timestep, 'timestep', error_set, .True.)
    
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

    ! Check info &coord_distrib block 
    If (traj_data%coord_distrib%invoke%fread) Then
      Call check_coord_distrib(files, traj_data, model_data)
    End If 
    
    ! Check info &track_unchanged_chemistry block 
    If (traj_data%unchanged%invoke%fread) Then
      Call check_unchanged_chemistry(files, traj_data, model_data)
    End If 

    ! Check info &lifetime block 
    If (traj_data%lifetime%invoke%fread) Then
      Call check_lifetime(files, traj_data)
    End If 

    ! Check info &orientational_chemistry 
    If (traj_data%chem_ocf%invoke%fread) Then
      Call check_orientational_chemistry(files, traj_data)
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
       Write (messages(1),'(1x,a)') 'ERROR: the &region block contains no settings!' 
       Call info(messages, 1)
       Call error_stop(' ')
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
        If (model_data%input_composition%tag(j)==traj_data%unchanged%tag%type) Then
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
      If ((model_data%config%atom(k)%tag)/=(traj_data%unchanged%tag%type)) Then
        Call info(' ', 1)
        Write(word,*) k
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Index "'//Trim(Adjustl(word))//'" (defined in list_indexes)&
                                       & does not correspond to the atomic tag "'//Trim(traj_data%unchanged%tag%type)//'".'  
        Write (messages(2),'((1x,a))') 'According to the &input_composition block, this index corresponds to atomic&
                                       & tag "'//Trim(model_data%config%atom(k)%tag)//'".&
                                       & Please review the labels and indexes of the atomic model'
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
          If (traj_data%rdf%type_a(j)==traj_data%rdf%type_a(k)) Then
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
          If (traj_data%rdf%type_b(j)==traj_data%rdf%type_b(k)) Then
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
    
  End Subroutine check_rdf

  Subroutine check_coord_distrib(files, traj_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &coord_distrib block
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data
    Type(model_type),   Intent(In   ) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer(Kind=wi)    :: j
    Logical             :: flag

    Character(Len=8)  :: tg
    Character(Len=8)  :: coord(3)
    
    ! Define coordinates to check directive "coordinate"
    coord(1)='x'
    coord(2)='y'
    coord(3)='z'
    
    error_set = '***ERROR in the &coord_distrib block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (.Not. traj_data%coord_distrib%delta%fread) Then
      traj_data%coord_distrib%delta%tag='delta'
    End If
    Call check_length_directive(traj_data%coord_distrib%delta, error_set, .True., 'directive')

    ! Check definition of "species" directive
    If (traj_data%coord_distrib%species_dir%fread) Then
      If (traj_data%coord_distrib%species_dir%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "species" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "species" directive to&
                                    & compute the coordinate distribution.&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
   ! Check if the definition of "species" is valid
    tg=Trim(traj_data%coord_distrib%species)
    Call remove_symbols(tg,'*')
    flag=.True.
    j=1
    Do While (j <= model_data%input_composition%atomic_species .And. flag)
      If (Trim(model_data%input_composition%tag(j))==Trim(tg)) Then
        flag=.False.
      End If  
      j=j+1
    End Do
    If (flag) Then
      Write (messages(1),'(2(1x,a))') Trim(error_set), '"'//Trim(traj_data%coord_distrib%species)//'"&
                                     & defined for the "species" directive is not a valid species.&
                                     & Please review the definition of the &input_composition block' 
      Call info(messages, 1)
      Call error_stop(' ') 
    End If 

    ! Check definition of "species" directive
    If (traj_data%coord_distrib%species_dir%fread) Then
      If (traj_data%coord_distrib%species_dir%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "species" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "species" directive to&
                                    & compute the coordinate distribution.&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
    ! Check definition of "coordinate" directive
    If (traj_data%coord_distrib%coordinate%fread) Then
      If (traj_data%coord_distrib%coordinate%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "coordinate" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        flag=.True. 
        Do j=1, 3
          If (Trim(traj_data%coord_distrib%coordinate%type)==Trim(coord(j))) Then
            flag=.False.
            traj_data%coord_distrib%indx=j
          End If
        End Do
        If (flag) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Definition for the "coordinate" directive&
                                    & is not valid. Valid options: "x", "y" or "z".&
                                    & Check correctness of the directives within the block.'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
      
      
    Else
      Write (messages(1),'(2(1x,a))') Trim(error_set), 'The user must define the "coordinate" directive to&
                                    & compute the coordinate distribution. Valid options: "x", "y" or "z".&
                                    & Check if the other directives have been defined correctly'
      Call info(messages, 1)
      Call error_stop(' ')
    End If
    
  End Subroutine check_coord_distrib
  
  
  Subroutine check_data_analysis(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &data_analysis block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set, message

    error_set = '***ERROR in the &data_analysis block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%analysis%invoke%fread) Then
      Call check_time_directive(traj_data%analysis%time_interval, 'time_interval',  error_set, .False.)
      Call check_time_directive(traj_data%analysis%end_time, 'end_time',  error_set, .False.)      
      Call check_time_directive(traj_data%analysis%ignore_initial, 'ignore_initial', error_set, .False.)
      Call check_time_directive(traj_data%analysis%overlap_time, 'overlap_time' ,error_set, .False.)

      If (traj_data%analysis%normalise_at_t0%fread) Then
        If (traj_data%analysis%normalise_at_t0%fail) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                    & "normalise_at_t0" (choose either .True. or .False.)'
          Call info(message,1)
          Call error_stop(' ')
        End If
      Else
        traj_data%analysis%normalise_at_t0%stat=.False.
      End If

    Else
      traj_data%analysis%normalise_at_t0%stat=.False.
    End If
    
  End Subroutine check_data_analysis

  Subroutine check_orientational_chemistry(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &orientational_chemistry block
    !
    ! author    - i.scivetti Febraury 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(traj_type),    Intent(InOut) :: traj_data

    Character(Len=256)  :: error_set
    Character(Len=256)  :: messages(2)

    error_set = '***ERROR in the &orientational_chemistry block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (traj_data%chem_ocf%variable%fread) Then
      If (traj_data%chem_ocf%variable%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "variable" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%chem_ocf%variable%type)/='special_pair'     .And. &
            Trim(traj_data%chem_ocf%variable%type)/='acceptor_donor_transfer_couple')  Then
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

    If (traj_data%chem_ocf%print_all_intervals%fread) Then
      If (traj_data%chem_ocf%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      traj_data%chem_ocf%print_all_intervals%stat=.False.
    End If
    
  End Subroutine check_orientational_chemistry
  
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

    Call check_time_directive(traj_data%lifetime%rattling_wait, 'rattling_wait' ,error_set, .False.)

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
        If (Trim(traj_data%lifetime%method%type)/='hicf'     .And. &
            Trim(traj_data%lifetime%method%type)/='hdcf')  Then
             Write (messages(1),'(2(1x,a))') Trim(error_set), &
                                    & 'Wrong input for "method". Valid options:&
                                    & "HICF" and "HDCF"'
          Call info(messages, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "method" directive'
       Write (messages(2),'( (1x,a))') 'Valid options: "HICF" and "HDCF"'
       Call info(messages, 2)
       Call error_stop(' ')
    End If

    If (traj_data%lifetime%print_all_intervals%fread) Then
      If (traj_data%lifetime%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      traj_data%lifetime%print_all_intervals%stat=.False.
    End If
    
  End Subroutine check_lifetime

  Subroutine check_time_directive(T, tag, error_set, kill)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check time related directivesd
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: tag 
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
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(tag)//'" directive'
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

    Character(Len=256)  :: error_set, warn_set, message
    Integer(Kind=wi)    :: i, k, j, l, kref, kini, frames
    Logical             :: flag, one_seg
    Real(Kind=wp)       :: time, tend, teff, tini, tref
    
    error_set = '***ERROR in the &data_analysis block of file '//Trim(files(FILE_SET)%filename)//' -'
    warn_set = '***WARNING from the &data_analysis block of file '//Trim(files(FILE_SET)%filename)//' -'

   
    If (traj_data%analysis%end_time%fread) Then
      If (traj_data%analysis%end_time%value > (traj_data%frames-1)*traj_data%timestep%value) Then
        Write (message,'(2(1x,a))') Trim(warn_set), 'The value assigned to "'//Trim(traj_data%analysis%end_time%tag)//&
                                 &'" is larger than the total time for the trajectory. The analysis will be performed&
                                 & up to largest recorded time.'
        Call info(message, 1)          
        tend=(traj_data%frames-1)*traj_data%timestep%value
        frames=traj_data%frames
      Else

        If (traj_data%timestep%value>=traj_data%analysis%end_time%value) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%end_time%tag)//&
                               &'" must be larger that the timestep for the trajectory. Please check the directives.'
          Call info(message, 1) 
          Call error_stop(' ')
        End If      
        
        If (traj_data%analysis%end_time%value > (traj_data%frames-2)*traj_data%timestep%value) Then
          Write (message,'(2(1x,a))') Trim(warn_set), 'The value assigned to "'//Trim(traj_data%analysis%end_time%tag)//&
                                   &'" is in between the last two recorded times. The analysis will be performed&
                                   & up to largest recorded time.'
          Call info(message, 1)          
          tend=(traj_data%frames-1)*traj_data%timestep%value
          frames=traj_data%frames
        Else
          tend=traj_data%analysis%end_time%value
          i=1
          flag=.True.
          Do While (i <= traj_data%frames .And. flag)
            time=(i-1)*traj_data%timestep%value
            If (time >= tend) Then
              frames=i
              flag=.False.
            End If
            i=i+1
          End do           
        End If
      End If
    Else
      tend=(traj_data%frames-1)*traj_data%timestep%value
      frames=traj_data%frames
    End If
    
    ! Set the net number of frames
    traj_data%analysis%frame_last=frames
    
    If (.Not. traj_data%analysis%ignore_initial%fread) Then
       traj_data%analysis%ignore_initial%value=-traj_data%timestep%value
       traj_data%analysis%frame_ini = 1
       tini=0.0_wp
    Else
      If (tend <= traj_data%analysis%ignore_initial%value) Then
        Call info(' ', 1)
        If (.Not. traj_data%analysis%end_time%fread) Then
          Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%ignore_initial%tag)//&
                                 &'" is larger than (or equal) the total time for the trajectory. Please check&
                                 & the settings and the value for the "timestep" directive.'
        Else
          Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%ignore_initial%tag)//&
                                 &'" is larger than (or equal) the value set for "'//Trim(traj_data%analysis%end_time%tag)//&
                                 &'". Please check settings'        
        End If
        Call info(message, 1)
        Call error_stop(' ')
      Else  
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
      End If   
    End If
  
  
    teff=tend-tini 

    ! Compare timestep with other time settings of &data_analysis
    If (traj_data%analysis%time_interval%fread) Then
      If (traj_data%timestep%value>=traj_data%analysis%time_interval%value) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%time_interval%tag)//&
                               &'" must be larger that the timestep for the trajectory. Please check the "timestep" directive.'
        Call info(message, 1) 
        Call error_stop(' ')
      End If
      If (teff<traj_data%analysis%time_interval%value) Then
           Write (message,'(2(1x,a))') Trim(warn_set), 'The input value for the "'//Trim(traj_data%analysis%time_interval%tag)//&
                                   &'" directive was too large and has been redefined to comply with the rest&
                                   & of the settings and the length of the trajectory.'
          Call info(message, 1)
          traj_data%analysis%time_interval%value=teff
      End If
    Else
      traj_data%analysis%time_interval%value=teff
    End If
 

    If (traj_data%analysis%overlap_time%fread) Then
      If (traj_data%timestep%value > traj_data%analysis%overlap_time%value) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'The value assigned to "'//Trim(traj_data%analysis%overlap_time%tag)//&
                                 &'" must be larger that the timestep for the trajectory. Please check values (and units).'
        Call info(message, 1) 
        Call error_stop(' ')
      End If
      If (teff<traj_data%analysis%overlap_time%value) Then
           Write (message,'(2(1x,a))') Trim(warn_set), 'The input value for the "'//Trim(traj_data%analysis%overlap_time%tag)//&
                                   &'" directive was too large and has been redefined to comply with the rest of the&
                                   & directive and the length of the trajectory.'
          Call info(message, 1)
          traj_data%analysis%time_interval%value=teff
      End If      
      
    Else
      traj_data%analysis%overlap_time%value=teff+traj_data%timestep%value      
    End If
    
    ! Calculate the number of segments
    i=0; j=0; l=0
    one_seg=.True.
    k=traj_data%analysis%frame_ini
    tref=tini; kini=k; flag=.True.
    Do While (k <= frames)
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
      Do While (k <= frames)
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

    If (traj_data%msd%print_all_intervals%fread) Then
      If (traj_data%msd%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      traj_data%msd%print_all_intervals%stat=.False.
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
    
    If (traj_data%ocf%print_all_intervals%fread) Then
      If (traj_data%ocf%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      traj_data%ocf%print_all_intervals%stat=.False.
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

    Write (messages(1),'(1x,a)') '- by specification, the trajectory corresponds to the "'&
                                 &//Trim(traj_data%ensemble%type)//'" ensemble and was recorded in "'&
                                 &//Trim(model_data%input_geometry_format%type)//'" format'
    Call info(messages, 1)
    Write(word,'(f10.2)') traj_data%timestep%value
    Write (messages(1),'(1x,a,f5.2,a)') '- the time step between recorded configurations is '//Trim(Adjustl(word))//' fs'
    Call info(messages, 1)

    If (traj_data%analysis%end_time%fread) Then
      If ((traj_data%analysis%end_time%value-(traj_data%frames-1)*traj_data%timestep%value)<0.0_wp) Then
        Write(word,'(f8.3)') traj_data%analysis%end_time%value/1000.0_wp
        Write (messages,'(1x,a)') '- the analysis will consider up to '//Trim(Adjustl(word))//' ps of the trajectory' 
        Call info(messages, 1)
      End If
    End If        
    
    If (traj_data%analysis%ignore_initial%fread) Then
      Write(word,'(f8.3)') traj_data%analysis%ignore_initial%value/1000.0_wp
      Write (messages(1),'(1x,a)') '- the initial '//Trim(Adjustl(word))//' ps of the trajectory&
                                   & will be discarded' 
      Call info(messages, 1)
    End If
    
    
    If (traj_data%ocf%invoke%fread .Or. traj_data%chem_ocf%invoke%fread .Or. &
        traj_data%msd%invoke%fread .Or. traj_data%lifetime%invoke%fread) Then
      If (traj_data%analysis%time_interval%fread .Or. &
          (traj_data%analysis%overlap_time%fread  .And. (traj_data%analysis%N_seg /=1))) Then
        Write (messages(1),'(1x,a)') 'Instructions for analysis in time segments (&data_analysis) applied to:' 
        Call info(messages, 1)
        If (traj_data%ocf%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* OCF'
          Call info(messages, 1)
        End If  
        If (traj_data%lifetime%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* TCF'
          Write (messages(2),'(5x,a)') '* SPCF'
          Call info(messages, 2)
        End If  
        If (traj_data%msd%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* MSD'
          Call info(messages, 1)
        End If
        If (traj_data%chem_ocf%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* CHEM_OCF'
          Call info(messages, 1)
        End If
        If (traj_data%analysis%time_interval%fread) Then
          Write(word,'(f8.3)') traj_data%analysis%time_interval%value/1000.0_wp
          If(traj_data%analysis%N_seg /= 1) Then
            Write (messages(1),'(2x,a)') '- the data analysis will be executed using time intervals of '&
                                       &//Trim(Adjustl(word))//' ps' 
          Else
            Write (messages(1),'(2x,a)') '- the data analysis will be executed without the use of time intervals'           
          End If
          Call info(messages, 1)
        End If
        If (traj_data%analysis%overlap_time%fread  .And. (traj_data%analysis%N_seg /=1)) Then
          Write(word,'(f8.3)') traj_data%analysis%overlap_time%value/1000.0_wp
          Write (messages(1),'(2x,a)') '- the starting points of segments for analysis are separated&
                                         & by '//Trim(Adjustl(word))//' ps'
          Call info(messages, 1)
        End If
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
              Trim(traj_data%ocf%u_definition%type)/='plane' .And. & 
              Trim(traj_data%ocf%u_definition%type)/='bond_123') Then 
            Write (messages(1),'(1x,a)')  '**WARNING: the "'//Trim(traj_data%ocf%u_definition%type)//'"& 
                                           & option for the rotating unit vector (u_definition)&
                                           & is not recommended for most studies.'
            Write (messages(2),'(1x,a)')  '           Unless the user is fully certain, either the&
                                           & "bond_12-13", "plane" or "bond_123" option should be used instead.'
            Call info(messages, 2)
          End If
        End If
      End If
      Write (messages(1),'(1x,a)') 'The definition of the "&OCF" block will compute the Orientational&
                                  & Correlation Funtion (OCF) of the species "'&
                                  &//Trim(model_data%species_definition%name%type)//&
                                  & '" (defined in the &monitored_species block) as follows:'
      Write (messages(2),'(1x,a)') '- the attached rotating unit vector is defined with the method: '//&
                                  & Trim(traj_data%ocf%u_definition%type)
      Write (messages(3),'(1x,a, i2)') '- the correlation terms are obtained using the Legendre polynomial of order ',&
                                  & traj_data%ocf%legendre_order%value
      Call info(messages, 3)
    End If

    If (traj_data%msd%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&MSD" block will execute a Mean Square&
                                  & Displacement analysis of the species "'&
                                  &//Trim(model_data%species_definition%name%type)//&
                                  & '" (defined in the &monitored_species block) as follows:'
      Write (messages(2),'(1x,a)') '- the values will be computed for the coordinates(s): '//&
                                  & Trim(traj_data%msd%select%type)
      Call info(messages, 2)

      If (traj_data%msd%pbc_xyz%fread) Then
        Write (messages(1),'(1x,a)') '- the "pbc_xyz" directive specifies which coordinate uses (or not) periodic&
                                     & boundary conditions' 
        Call info(messages, 1)                             
      End If
    End If

    If (traj_data%lifetime%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&lifetime" block will compute:'
      Write (messages(2),'(1x,a)') '- the Transfer Correlation Function (TCF) and the Special Pair Correlation&
                                  & Function (SPCF) for the changing chemical species&
                                  & using the method: '//Trim(traj_data%lifetime%method%type)
      Write (messages(3),'(1x,a)') '- the residence times for each changing species (file RES_TIMES):'
      If (traj_data%lifetime%rattling_wait%fread) Then
        Write(word,'(f8.3)') traj_data%lifetime%rattling_wait%value/1000.0_wp
        Write (messages(4),'(3x,a)') 'Rattling times lower than '//Trim(Adjustl(word))//' ps will&
                                    & be discarded for the computation of residence times'  
      Else
        Write (messages(4),'(3x,a)') 'Rattling effects are included in the calculation&
                                     & for the computation of residence times'
      End If
      Call info(messages, 4)
    End If

    If (traj_data%chem_ocf%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&orientational_chemistry" block will compute&
                                  & OCF for the changing chemical species (CHEM_OCF) using the "'&
                                  &//Trim(traj_data%chem_ocf%variable%type)//'" as the orientational vector'
      Call info(messages, 1)
    End If
    
    If (traj_data%unchanged%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&track_unchanged_chemistry" block will print the positions&
                                  & of the selected atomic indexes with unchanged chemistry along the trajectory.'
      Call info(messages, 1)
    End If
    
    If (traj_data%rdf%invoke%fread) Then
      Call info(' ', 1)
      Write(word,'(f10.2)') traj_data%rdf%dr%value
      Write (messages(1),'(1x,a)') 'The definition of the "&RDF" block will compute the Radial&
                                  & Distribution Funtion (RDF) and the Coordination Numbers (CN) using:'
      Write (messages(2),'(1x,a)') '- the tags defined in "tags_species_a"  and "tags_species_b"'
      Write (messages(3),'(1x,a)') '- a discretization of '//Trim(Adjustl(word))//' Angstrom'
      Call info(messages, 3)
    End If

    If (traj_data%coord_distrib%invoke%fread) Then
      Call info(' ', 1)
      Write(word,'(f10.2)') traj_data%coord_distrib%delta%value
      Write (messages(1),'(1x,a)') 'The definition of the "&coord_distrib block" will compute the distribution&
                                  & of the '//Trim(traj_data%coord_distrib%coordinate%type)//'-values for all&
                                  & the "'//Trim(traj_data%coord_distrib%species)//'" species in the whole system'
      Write (messages(2),'(1x,a)') 'with a selected discretization of '//Trim(Adjustl(word))//' Angstrom for the coordinate.'
      Call info(messages, 2)
    End If
    
    If (model_data%species_definition%intra_geom%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&intramol_stat_settings" will compute the probability&
                                  & distribution for the intramolecular:'
      Call info(messages, 1)                            
      If (model_data%species_definition%intra_geom%dist%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- distances, using the settings of "&distance_parameters"'
        Call info(messages, 1)                            
      End If
      If (model_data%species_definition%intra_geom%angle%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- angles, using the settings of "&angle_parameters"'
        Call info(messages, 1)                            
      End If
      Write (messages(1),'(1x,a)') 'corresponding to the species "'//Trim(model_data%species_definition%name%type)//&
                                  & '" (defined in the &monitored_species block).'
      Call info(messages, 1)                          
    End If

    If (model_data%nndist%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&selected_nn_distances" block will compute probability distribution&
                               & of the shortest distance of the selected pair of species (this is not RDF)'
      Call info(messages, 1)                             
    End If
    
    If (model_data%species_definition%inter_geom%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&intermol_stat_settings" will compute the probability&
                                  & distribution for the intermolecular:'
      Call info(messages, 1)                            
      If (model_data%species_definition%intra_geom%dist%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- distances, using the settings of "&distance_parameters"'
        Call info(messages, 1)                            
      End If
      If (model_data%species_definition%intra_geom%angle%invoke%fread) Then
        Write (messages(1),'(1x,a)') '- angles, using the settings of "&angle_parameters"'
        Call info(messages, 1)                            
      End If
      Write (messages(1),'(1x,a)') 'by considering the two closest "'//Trim(model_data%species_definition%name%type)//& 
                                  &'" species to each "'//Trim(model_data%species_definition%name%type)//'" species&
                                  & (see the &monitored_species block).'
      Call info(messages, 1)
    End If
    
    If (traj_data%region%define%fread) Then
      Call info(' ', 1)
      If (traj_data%ocf%invoke%fread .Or. &
          traj_data%msd%invoke%fread .Or. &
          traj_data%rdf%invoke%fread .Or. &
          traj_data%chem_ocf%invoke%fread .Or. &
          model_data%nndist%invoke%fread .Or. &
          model_data%species_definition%intra_geom%invoke%fread .Or. &
          model_data%species_definition%inter_geom%invoke%fread) Then
        Write (messages(1),'(1x,a)') 'From the definition of the "&region" block, the computation of'
        Call info(messages, 1)
        If (traj_data%ocf%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- OCF'
          Call info(messages, 1)
        End If  
        If (traj_data%lifetime%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- TCF'
          Write (messages(2),'(3x,a)') '- SPCF'
          Call info(messages, 2)
        End If  
        If (traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- RDF'
          Call info(messages, 1)
        End If  
        If (traj_data%msd%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- MSD'
          Call info(messages, 1)
        End If
        If (traj_data%chem_ocf%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- CHEM_OCF'
          Call info(messages, 1)
        End If
        If (model_data%species_definition%intra_geom%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- Intramolecular parameters (monitored species)'
          Call info(messages, 1)
        End If
        If (model_data%species_definition%inter_geom%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- Intermolecular parameters (monitored species)'
          Call info(messages, 1)
        End If
        If (model_data%nndist%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- shortest distance distribution for the selected pair'
          Call info(messages, 1)
        End If
        Write (messages(1),'(1x,a)') 'will be only carried out:'
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
        If (traj_data%rdf%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'IMPORTANT: For the RDF analysis, the definition of the &region block&
                                      & only applies to the species listed in "tags_species_a" (&rdf block)'
          Call info(messages, 1)
        End If  
        If (model_data%nndist%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'IMPORTANT: For the analysis of the shortest distance distribution,&
                                      & the definition of the &region block only applies to the species& 
                                      & listed in "reference_species" (&selected_nn_distances block)'
          Call info(messages, 1)
        End If  
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

    If (traj_data%chem_ocf%invoke%fread) Then
      If (.Not. model_data%change_chemistry%stat) Then
         Write (messages(1),'(1x,a)') Trim(error_set)//' The user has defined the &orientational_chemistry block but&
                                     & the &search_chemistry block is missing.'
         Write (messages(2),'(4x,a)') 'The computation of the orientational chemistry is&
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
        traj_data%track_chem%config(frame,l)%r=model_data%track_chem(l)%r
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
