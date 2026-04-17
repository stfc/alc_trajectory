!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the trajectory settings and identification. 
! If the directive change_chemistry
! is set to .True., the algorithm searches and tracks changes of 
! chemical species based on the information of the &search_chemistry
! block. 
!
! Copyright   2023-2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2023
! Refact      -  i.scivetti  Feb 2026           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module trajectory 

  Use atomic_model, Only : model_type,&
                           about_cell,&
                           atomistic_model, & 
                           check_definition_bonds, &
                           check_PBC,&
                           check_cell_consistency,&
                           check_orthorhombic_cell,&
                           compute_distance_pbc,&
                           identify_monitored_indexes,&
                           obtain_maximum_number_species,&
                           read_model

  Use constants,    Only : max_components, &
                           max_at_species, &
                           max_unchanged_atoms, &
                           initial_tolerance

  Use fileset,      Only : file_type, &
                           FILE_SET, & 
                           FILE_TRAJECTORY,&
                           FILE_TRACK_CHEMISTRY,&
                           FILE_TAGGED_TRAJ, &
                           FILE_UNCHANGED_CHEM, & 
                           refresh_out

  Use input_types,  Only : in_logic,   &
                           in_param,   & 
                           in_string

  Use numprec,      Only : wi,&
                           wp 
  
  Use process_data, Only : capital_to_lower_case, &
                           check_end, &
                           check_for_rubbish,&
                           get_word_length, &
                           remove_symbols, &
                           set_read_status
                           
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

  ! Type for eqcm data and analysis
  Type, Public :: traj_type
    Private
    Logical                           :: reload_trajectory
    Logical,                  Public  :: active_bonds_computed
    Type(in_logic),           Public  :: print_retagged_trajectory
    Type(in_param),           Public  :: timestep
    Integer(Kind=wi),         Public  :: frames
    Integer(Kind=wi),         Public  :: Nmax_species
    Integer(Kind=wi),         Public  :: N_species
    Type(in_string),          Public  :: ensemble
    Type(region_type),        Public  :: region
    Type(unchanged_type),     Public  :: unchanged
    Type(analysis_type),      Public  :: analysis
    Type(track_type),         Public  :: track_chem                      
    Type(species_type),   Public, Allocatable :: species(:,:)
    Type(atom_type),      Public, Allocatable :: config(:,:)
    Type(box_type),       Public, Allocatable :: box(:)
  Contains
    Private
      Procedure         :: alloc_trajectory  => allocate_trajectory_arrays
      Procedure         :: alloc_analysis    => allocate_analysis_arrays
      Final             :: cleanup
  End Type traj_type

  Public :: extract_trajectory, trajectory_setup  
  Public :: check_trajectory_settings, check_time_directive
  Public :: read_region, read_track_unchanged_chemistry, read_settings_segment_analysis
  Public :: within_region, compute_closest_pairs, define_trajectory_segments
  Public :: find_neighbours_monitored_species, find_active_bonds, average_segments
  Public :: print_unchanged_chemistry, print_tracking_species
    
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

    If (model_data%change_chemistry%stat) Then 
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

    If (Allocated(T%species)) Then
      Deallocate(T%species)
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
    ! Subroutine to extract the trajectory from file
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    
    Logical            :: loop_traj
    Character(Len=256) :: message, nframes, md_length, net_md_length
    Character(Len=256) :: input_file
    Integer(Kind=wi)   :: i, j
    
    input_file=(files(FILE_TRAJECTORY)%filename)
    
    ! Open the TRAJECTORY file
    Open(Newunit=files(FILE_TRAJECTORY)%unit_no, File=(input_file),Status='old')
    Write (md_length,'(f12.3)') (traj_data%frames-1)*traj_data%timestep%value/1000.0_wp
    Write (nframes,'(i8)') traj_data%frames
    Call info(' ', 1)
    Call info('Start of the analysis', 1)
    Call info('=====================', 1)
    Write (message,'(1x,a)') 'The code has identified a total of '//Trim(Adjustl(nframes))//' frames. From&
                                 & the setting of the "timestep" directive, the recorded MD trajectory is '&
                                 &//Trim(Adjustl(md_length))//' ps long.'
    Call info(message, 1)
    
    If (traj_data%analysis%end_time%fread) Then
      If ((traj_data%analysis%end_time%value-(traj_data%frames-1)*traj_data%timestep%value)<0.0_wp) Then
        Write (net_md_length,'(f8.3)') traj_data%analysis%end_time%value/1000.0_wp
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

    If (model_data%change_chemistry%stat) Then
     Call residence_percentage(model_data, traj_data)
     Call refresh_out(files) 
    End If
    
    If (traj_data%print_retagged_trajectory%stat) Then 
      input_file=(files(FILE_TAGGED_TRAJ)%filename)
      Call print_tagged_trajectory(files, model_data, traj_data)
      Write (message,'(1x,a)') 'A copy of the trajectory with modified tags for the atomic species was printed&
                              & to the "'//Trim(input_file)//'" file'
      Call info(message, 1)
      Call refresh_out(files)
    End If 
  
  End Subroutine extract_trajectory
  
  Subroutine trajectory_setup(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to setup variables and check variables against the trajectory
    !
    ! refact    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
    
    Logical            :: safe, fortho
    Character(Len=256) :: message
    Character(Len=256) :: input_file, set_error
    Integer(Kind=wi)   :: j
    
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
                                  & for orthorhombic cells. We do not guarantee a correct analysis.'
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
    If (model_data%change_chemistry%stat) Then 
      Call check_definition_bonds(model_data, 1)
    End If
    
    ! Check if the defined region for analysis is within the simulation cell
    If (traj_data%region%define%fread) Then
      Call check_consistency_spatial_domain(model_data, traj_data, 1)
    End If 

    ! Search for the number of frames
    Call obtain_number_frames(files, model_data, traj_data)
    ! Allocate trajectory arrays
    Call traj_data%alloc_trajectory(model_data)
  
  End Subroutine trajectory_setup  
  
  Subroutine find_neighbours_monitored_species(model_data, traj_data, flag_exec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the two closest monitored species to a
    ! monitored species
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(In   ) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data
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
  
  Subroutine residence_percentage(model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the residence percentage of changing chemistry 
    ! species along the trajectory
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(In   ) :: model_data
    Type(traj_type),   Intent(InOut) :: traj_data

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

        If (match_j) Then
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
      call error_stop('ERROR')
    End If
    
    If (s1==s3) Then
      call error_stop('ERROR')
    End If
    
    If (s2==s3) Then
      call error_stop('ERROR')
    End If
  
  End Subroutine compute_closest_pairs
 
  Subroutine average_segments(files, traj_data, file_number, what, msd_xyz)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to average physical quantities computed for each time segment
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),              Intent(InOut) :: files(:)
    Type(traj_type),              Intent(InOut) :: traj_data
    Integer(Kind=wi),             Intent(In   ) :: file_number
    Character(Len=*),             Intent(In   ) :: what
    Character(Len=256), Optional, Intent(In   ) :: msd_xyz
  
    Integer(Kind=wi) :: i, k, iunit, Nnet 
    Real(Kind=wp)    :: sum_i, average, average_net, std, norma
    Logical          :: flag, fcompute
    Character(Len=256) :: message
    Character(Len=256) :: quantity
    !Character(Len=256),  :: quantity

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
        If (Abs(norma-1.0_wp)>initial_tolerance) Then
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
                        &Trim(msd_xyz)//'" of the monitored species'
      End If
      Write (iunit,'(a)') '#  Time (ps)      '//Trim(quantity)//'          STD' 
      
      i=1
      flag=.True.    
      Do While ((i<=traj_data%analysis%Ninterval) .And. flag)
        sum_i=0.0_wp
        Nnet=0
        Do k= 1, traj_data%analysis%N_seg
          If (i<=traj_data%analysis%max_points(k)) Then
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
              If (i<=traj_data%analysis%max_points(k)) Then
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
            If (average > 1.0_wp) Then
              average_net=1.0_wp
            Else
              average_net=average
            End If
          End If  
          Write (iunit,'(3(f10.3, 3x))') (i-1)*traj_data%timestep%value/1000.0_wp, average_net, std
        End If
        i=i+1
        
      End Do
      Write (message,'(1x,a)') 'The average '//Trim(what)//' was printed to the "'//&
                               &Trim(files(file_number)%filename)//'" file.'
      Call info(message, 1)
      Close(iunit)
    End If 
    
    
  End Subroutine average_segments
  
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
      Read (files(FILE_TRAJECTORY)%unit_no, Fmt= *, iostat=stat) check
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
        Write (word,*) k
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
      If (kill) Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(tag)//'" directive'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_time_directive  
  
  Subroutine define_trajectory_segments(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &data_analysis block and define
    ! segments for analysis
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
          traj_data%analysis%overlap_time%value=teff+traj_data%timestep%value
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
        l=k-kini+1 
        j=0
        If (time>=(tref+traj_data%analysis%overlap_time%value)) Then
          If (flag) Then
            kref=k
          End If
        Else
           If (traj_data%analysis%overlap_time%fread) Then
             kref=Nint((tref+traj_data%analysis%overlap_time%value)/traj_data%timestep%value)+1
           Else
             kref=k
           End If
        End If
        tref=(kref-1)*traj_data%timestep%value
        k=kref
        kini=k
        one_seg=.False.
        flag=.True.
      Else
        If (time>=(tref+traj_data%analysis%overlap_time%value) .And. flag) Then
          kref=k
          flag=.False.
        End If
        j=j+1
      End If
      k=k+1
    End Do

    If (one_seg) Then
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
      k=traj_data%analysis%frame_ini
      tref=tini; kini=k; flag=.True.
      Do While (k <= frames)
        time=(k-1)*traj_data%timestep%value
        If (time>=(tref+traj_data%analysis%time_interval%value)) Then
          i=i+1
          l=k-kini+1
          j=0
          traj_data%analysis%seg_indx(1,i)=kini
          traj_data%analysis%seg_indx(2,i)=k
          If (time>=(tref+traj_data%analysis%overlap_time%value)) Then
            If (flag) Then
              kref=k
            End If
          Else
            If (traj_data%analysis%overlap_time%fread) Then
              kref=Nint((tref+traj_data%analysis%overlap_time%value)/traj_data%timestep%value)+1
            Else
              kref=k
            End If
          End If
          tref=(kref-1)*traj_data%timestep%value
          k=kref
          kini=k
          flag=.True.
        Else
          If (time>=(tref+traj_data%analysis%overlap_time%value) .And. flag) Then
            kref=k
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
    
  End Subroutine define_trajectory_segments   

  Subroutine check_consistency_spatial_domain(model_data, traj_data, frame)
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
  
  End Subroutine check_consistency_spatial_domain
  
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
      Write (iunit,'(a)') '# Tracking the change of chemical species over the whole trajectory'    
    Else
      Write (iunit,'(a,1x,f10.4,1x,a)') '# Tracking the change of chemical species ignoring the first',& 
                                   &  traj_data%analysis%frame_ini*traj_data%timestep%value/1000_wp,&
                                   & 'ps of the whole trajectory. This value is set to time zero below.'
    End If

    If (model_data%chem%N0%value==1) Then
      Write (iunit,'(a,9x,a)') '# Time (ps)', 'XYZ_Species_1' 
    Else
      Write (num_species,*) model_data%chem%N0%value
      Write (iunit,'(a,9x,2a)') '# Time (ps)', 'XYZ_Species_1 .... XYZ_Species_', Trim(Adjustl(num_species)) 
    End If
    
    Do i = traj_data%analysis%frame_ini, traj_data%analysis%frame_last
       Write (iunit,'(f10.4, 4x, *(f11.3))') (i-traj_data%analysis%frame_ini)*traj_data%timestep%value/1000_wp,&
                                       & (traj_data%track_chem%config(i,l)%r(:), l=1, model_data%chem%N0%value)
    End Do
    Write (message,'(1x,a)') 'The tracking of the changing chemical species in xyz format was printed& 
                              & to the "'//Trim(files(FILE_TRACK_CHEMISTRY)%filename)//'" file'
    Call info(message, 1)
    Call refresh_out(files)
    Close(iunit)
    Call info(' ', 1)
  
  End Subroutine print_tracking_species

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
      If (model_data%change_chemistry%stat) Then
        Open(Newunit=files(FILE_TAGGED_TRAJ)%unit_no, File=files(FILE_TAGGED_TRAJ)%filename, Status='Replace')
        iunit=files(FILE_TAGGED_TRAJ)%unit_no
        Write(iunit,'(2i8,1x,a)') model_data%config%num_atoms, traj_data%frames, ' # number of total atoms and trajectory frames'
        Do l = 1, traj_data%frames
          Write(iunit,'(a,1x,i8)') 'Frame=', l 
          Do i = 1, model_data%config%num_atoms 
            Write (iunit,'(a, 4x, 3(f11.3), 4x, a)') Trim(traj_data%config(l,i)%element),&
                                                 & (traj_data%config(l,i)%r(k), k=1, 3),&
                                                 &  Trim(traj_data%config(l,i)%tag) 
          End Do               
        End Do
        Close(iunit)
      End If
  
  End Subroutine
  
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
      Write (iunit,'(a)') '# Tracking unchanged chemical species over the whole trajectory'    
    Else  
      Write (iunit,'(a,1x,f10.4,1x,a)') '# Tracking the unchanged chemical species ignoring the first',& 
                                   &  traj_data%analysis%frame_ini*traj_data%timestep%value/1000_wp,&
                                   & 'ps of the whole trajectory. This value is set to time zero below.' 
    End If
    
    If (traj_data%unchanged%N0==1) Then
      Write (iunit,'(a)') '# The label and number for the species is consistent with the settings&
                                   & of the "&track_unchanged_chemistry" block.'
    Else                               
      Write (iunit,'(a)') '# The species labelling, ordering and numbering is consistent with the settings&
                                   & of the "&track_unchanged_chemistry" block.'    
    End If
    
    spec=Trim(traj_data%unchanged%tag%type)
    Write (num1,*) traj_data%unchanged%indexes(1)
    If (traj_data%unchanged%N0==1) Then
      Write (iunit,'(a,5x,a)') '#  Time (ps)', 'XYZ_'//Trim(spec)//'_'//Trim(Adjustl(num1))
    Else
      Write (num2,*) traj_data%unchanged%indexes(traj_data%unchanged%N0)
      Write (iunit,'(a,5x,a)') '#  Time (ps)', 'XYZ_'//Trim(spec)//'_'//Trim(Adjustl(num1))//&
                              &'.... XYZ_'//Trim(spec)//'_'//Trim(Adjustl(num2))
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
        Write (iunit,'(f10.4, 1x, *(f11.3))') (i-traj_data%analysis%frame_ini)*traj_data%timestep%value/1000.0_wp,&
                & (traj_data%config(i,traj_data%unchanged%indexes(l))%r(:), l=1, traj_data%unchanged%N0)
      Else
        Write (messages(1),'(1x,a)') '**********************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: Problems with tracking species defined in the&
                                        & &track_unchanged_chemistry block'
        Write (num_species,*) k
        Write (frame,*)       i 
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
    If (model_data%change_chemistry%stat) Then 
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

  Subroutine read_settings_segment_analysis(iunit, traj_data)
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

      Else If (Trim(word)=='end_time') Then
        Read (iunit, Fmt=*, iostat=io) traj_data%analysis%end_time%tag, &
                                      & traj_data%analysis%end_time%value,& 
                                      & traj_data%analysis%end_time%units
        Call set_read_status(word, io, traj_data%analysis%end_time%fread,&
                                      & traj_data%analysis%end_time%fail)                                      
                                      
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
    
  End Subroutine read_settings_segment_analysis 
  
  Subroutine check_trajectory_settings(files, model_data, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the correctness of trajectory-related directives
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(In   ) :: files(:)
    Type(model_type),     Intent(In   ) :: model_data
    Type(traj_type),      Intent(InOut) :: traj_data

    Character(Len=256)  :: message
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    
    If (traj_data%print_retagged_trajectory%fread) Then
      If (traj_data%print_retagged_trajectory%fail) Then
        Write (message,'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_retagged_trajectory" (choose either .True. or .False.)'
        Call info(message, 1)
        Call error_stop(' ')
      End If
      If ((.Not. model_data%change_chemistry%stat) .And. traj_data%print_retagged_trajectory%stat) Then 
        Write (message,'(2(1x,a))') Trim(error_set), ' The user has set "print_retagged_trajectory" to .True. but&
                                      & "change_chemistry" is set to .False. Why do you want to retag the trajectory?&
                                      & Please change'
        Call info(message, 1)
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
        Write (message,'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "ensemble" directive.'
        Call info(message, 1)
        Call error_stop(' ')
      Else
        If (Trim(traj_data%ensemble%type)/='nve'  .And. &
            Trim(traj_data%ensemble%type)/='nvt'  .And. &
            Trim(traj_data%ensemble%type)/='npt') Then
             Write (message,'(2(1x,a))') Trim(error_set), &
                                    &'Wrong input for "ensemble". Valid options: "NVE", "NVT" and "NPT"'
          Call info(message, 1)
          Call error_stop(' ')
        End If
      End If
    Else
       Write (message,'(2(1x,a))')  Trim(error_set), 'The user must define the "ensemble" directive'
       Call info(message, 1)
       Call error_stop(' ')
    End If

   ! Check directives for data analysis
    Call check_data_analysis_block(files, traj_data)

    ! Check info &track_unchanged_chemistry block 
    If (traj_data%unchanged%invoke%fread) Then
      Call check_unchanged_chemistry(files, traj_data, model_data)
    End If 
    
   ! Check settings &region block 
    Call check_region(files, traj_data)
     
  End Subroutine check_trajectory_settings
 
 
  Subroutine check_region(files, traj_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &region block
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
            Write (word,*) traj_data%unchanged%indexes(j)
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
  
  Subroutine check_data_analysis_block(files, traj_data)
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
    
  End Subroutine check_data_analysis_block  
 
  
End Module trajectory
