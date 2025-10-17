!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module that reads and evaluates each atomic configuration of the trajectory
!
! Copyright   2023-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti  Feb  2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module atomic_model 

  Use constants,          Only : max_components, &
                                 max_at_species, &
                                 Bohr_to_A, &
                                 Rads_to_degrees, &
                                 chemsymbol, &
                                 NPTE
                                 
  Use fileset,            Only : file_type,              &
                                 FILE_TRAJECTORY,        &
                                 FILE_SET
  Use input_types,        Only : in_integer, &
                                 in_integer_array, & 
                                 in_logic,   &
                                 in_string,  &
                                 in_param,   &
                                 in_scalar
  Use numprec,            Only : li, &
                                 wi, &
                                 wp
  Use process_data,       Only : capital_to_lower_case,&
                                 remove_symbols
  Use unit_output,        Only : error_stop,&
                                 info 

  Implicit None
  Private

  ! Maximum number of bonded neighbours 
  Integer(Kind=wi), Parameter :: max_neighbours = 12  
  ! Tolerance for length
  Real(Kind=wp), Parameter, Public :: length_tol = 1.0E-8
  ! Limit for the minimum bond distance
  Real(Kind=wp), Parameter, Public :: min_intra = 0.70_wp  

  ! Type to describe the atoms
  Type :: tracking_type
    Real(Kind=wp)    :: r(3)
    Real(Kind=wp)    :: r0(3)
    Integer(Kind=wi) :: indx
    Character(Len=8) :: tag
  End Type 
  
  ! Type to describe the atoms
  Type :: atom_type
    Real(Kind=wp)    :: r(3)
    Real(Kind=wp)    :: r0(3)
    Character(Len=8) :: tag
    Character(Len=8) :: tag_0
    Character(Len=2) :: element
    Integer(Kind=wi) :: Nbonds
    Integer(Kind=wi) :: Nbonds0
    Integer(Kind=wi) :: bonds(max_neighbours)
    Integer(Kind=wi) :: bonds0(max_neighbours)
    Logical          :: dynamics(3)
    Logical          :: identified
  End Type 

  ! Type for the  &input_composition
  Type :: type_input_composition
    Type(in_string)  :: invoke
    Integer(Kind=wi) :: numtot
    Integer(Kind=wi) :: atomic_species
    Character(Len=8), Allocatable :: tag(:)
    Character(Len=2), Allocatable :: element(:)
    Integer(Kind=wi), Allocatable :: N0(:)
  End Type 

  ! Type for the definition of extra bonds
  Type :: extra_bonding 
    Type(in_string)    :: invoke
    Character(Len=8)   :: tg1(max_components)
    Character(Len=8)   :: tg2(max_components)
    Type(in_param)     :: bond(max_components)
    Logical            :: flag(max_components)
    Integer(Kind=wi)   :: N0
  End Type 

  Type :: bonding
    Type(in_logic)     :: criteria
    Type(in_string)    :: species
    Type(in_integer)   :: N0
    Type(in_param)     :: cutoff 
    Integer(Kind=wi)   :: list(max_components)
  End Type         

  Type :: environment
    Logical           :: check
    Type(in_logic)    :: criteria
    Type(in_param)    :: cutoff 
    Type(in_logic)    :: info_include
    Type(in_logic)    :: info_exclude
    Character(Len=8)  :: tg_incl(max_components)
    Integer(Kind=wi)  :: N0_incl
    Character(Len=8)  :: tg_excl(max_components)
    Integer(Kind=wi)  :: N0_excl
    Integer(Kind=wi)  :: list(max_components)
    Integer(Kind=wi)  :: num_search
    Integer(Kind=wi)  :: num_nn
    Integer(Kind=wi)  :: accum_indxs(max_components)
    Integer(Kind=wi)  :: order_indxs(max_components)
    Integer(Kind=wi)  :: missed(max_components)
  End Type        

  ! Type for the definition of new chemistry
  Type :: chemistry_type 
    Type(in_integer)   :: N0
    Integer(Kind=wi)   :: indx_new(max_components)
    Integer(Kind=wi)   :: indx_prev(max_components)
    Type(environment)  :: acceptor 
    Type(bonding)      :: bonds
  End Type 

  ! Type for geometrical paremeter
  Type, Public :: short_dist_type
    Type(in_string)    :: invoke
    Type(in_string)    :: tag_reference_species
    Type(in_string)    :: tag_nn_species
    Character(Len=8)   :: reference_species
    Integer(Kind=wi)   :: num_nn_species
    Character(Len=8)   :: nn_species(max_at_species)
    Type(in_param)     :: lower_bound
    Type(in_param)     :: upper_bound
    Type(in_param)     :: dr 
  End Type
  
  ! Type for geometrical paremeter
  Type, Public :: geo_param_type
    Type(in_string)      :: invoke
    Character(Len=256)   :: name
    Type(in_string)      :: tag_species
    Integer(Kind=wi)     :: nspecies
    Character(Len=8)     :: species(max_at_species)
    Integer(Kind=wi)     :: num_spec(max_at_species)
    Type(in_param)       :: lower_bound
    Type(in_param)       :: upper_bound
    Type(in_param)       :: delta 
  End Type

  ! Type for computation of the internal geometry
  Type, Public :: geo_spec_type
    Type(in_string)       :: invoke
    Character(Len=256)    :: tag
    Type(geo_param_type)  :: dist
    Type(geo_param_type)  :: angle
    Type(in_logic)        :: only_ref_tags_as_nn
  End Type 
  
  ! Type to describe the definition of the monitored species
  Type :: spec_def_type
    Type(in_string)       :: name
    Integer(Kind=wi)      :: num_components
    Type(in_string)       :: reference_tag
    Type(in_param)        :: bond_cutoff
    Integer(Kind=wi)      :: atoms_per_species 
    Type(in_string)       :: atomic_components
    Character(Len=8)      :: element(max_at_species)
    Integer(Kind=wi)      :: N0_element(max_at_species)
    Type(in_logic)        :: compute_amount
    Type(geo_spec_type)   :: intra_geom
    Type(geo_spec_type)   :: inter_geom
  End Type 

  ! Types related to format for VASP files 
  Type :: list_type
    Character(Len=8)  :: tag(max_components)
    Character(Len=2)  :: element(max_components) 
    Integer(Kind=wi)  :: N0(max_components)   
    Integer(Kind=wi)  :: num_elements
    Integer(Kind=wi)  :: net_elements
    Character(Len=32) :: coord_type 
  End Type 
 
  ! Type to describe species
  Type :: type_species
    Integer(Kind=wi) :: list(max_at_species)
    Logical          :: alive
  End Type
  
  ! Type for the configuration 
  Type :: config_type
    !!! Size of the model 
    !!!!!!!!!!!!!!!!!!!!
    ! Indicator for the simulation cell 
    Type(in_string),  Public :: simulation_cell 
    ! Cell vectors 
    Real(Kind=wp)      :: cell(3,3)
    ! Inverse cell vectors
    Real(Kind=wp)      :: invcell(3,3)
    ! Length of cell vectors
    Real(Kind=wp)      :: cell_length(3)
    ! Cell volume
    Real(Kind=wp)      :: volume
    ! Cell units
    Type(in_string)     :: cell_units
    !!! Atomic configuration
    !!!!!!!!!!!!!!!!!!!!!!!!
    ! Total number of atoms  
    Integer(Kind=wi) :: num_atoms
    ! Allocation flags
    Logical, Public     :: allocated_model_geo
    ! Position units
    Type(in_string)     :: position_units
    ! size scaling
    Real(Kind=wp)       :: cell_scaling 
    ! position scaling
    Real(Kind=wp)       :: position_scaling
    ! Scale factor vasp
    Real(Kind=wp)       :: scale_factor_vasp
    ! Atoms in the model
    Type(atom_type), Allocatable :: atom(:)
    ! Arrays for reading VASP trajectories
    Character(Len= 2), Allocatable :: element_file(:)
    Integer(Kind=wi),  Allocatable :: amount_file(:)
    ! list elements
    Type(list_type)  :: list
    ! Constrained dynamics
    Logical          :: selective_dyn
    !!! Variables related to species
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Indicator for the monitored species
    Type(in_string),  Public    :: monitored_species 
    ! Maximum species number
    Integer(Kind=wi)   :: Nmax_species
    ! Species
    Type(type_species), Allocatable :: species(:)
  End Type  

  ! Type for the modelling related variables 
  Type, Public :: model_type
    Private
    ! Flag to specify if chemistry is changed wrt initial species
    Type(in_logic), Public                :: change_chemistry
    ! Extra bonds    
    Type(extra_bonding), Public           :: extra_bonds
    ! Indicator to search for chemistry 
    Type(in_string), Public               :: search_chemistry 
    ! new chemistry
    Type(chemistry_type), Public          :: chem
    ! Number of extra species found
    Integer(Kind=wi), Public              :: number_species_found
    ! Type for those components defined in input_composition
    Type(type_input_composition), Public  ::  input_composition 
    ! Format of input trajectory 
    Type(in_string), Public               :: input_geometry_format
    ! Input model 
    Type(config_type), Public             :: config 
    ! definition of species to be monitored 
    Type(spec_def_type), Public           :: species_definition
    ! Tracked species
    Type(tracking_type), Public           :: track_chem(max_components)
    ! Shortest pair
    Type(short_dist_type), Public         :: nndist

   Contains
     Private
       Procedure          :: init_atomic_arrays       => allocate_input_atomic_arrays
       Procedure          :: init_vasp_arrays         => allocate_arrays_vasp_trajectory
       Procedure, Public  :: init_species             => allocate_species_arrays
       Procedure, Public  :: init_input_composition   => allocate_input_composition_arrays
       Final              :: cleanup
  End Type model_type

  Public :: read_model, check_model_settings, print_model_settings
  Public :: atomistic_model, check_definition_bonds, about_cell
  Public :: obtain_maximum_number_species, identify_monitored_indexes
  Public :: check_PBC, check_cell_consistency, check_orthorhombic_cell, check_length_directive 
  Public :: compute_distance_pbc
  
Contains

  Subroutine allocate_input_composition_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for the input model
    !
    ! author    - i.scivetti April 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),  Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(3)
    Character(Len=256)  :: message

    fail=0

    Allocate(T%input_composition%tag(T%input_composition%atomic_species),       Stat=fail(1))
    Allocate(T%input_composition%N0(T%input_composition%atomic_species),        Stat=fail(2))
    Allocate(T%input_composition%element(T%input_composition%atomic_species),   Stat=fail(3))

    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the input composition&
                                & (subroutine allocate_input_composition_arrays).'
      Call error_stop(message)
    End If

  End Subroutine allocate_input_composition_arrays

  Subroutine allocate_input_atomic_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for the input model
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),  Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(1)
    Character(Len=256)  :: message

    fail=0

    Allocate(T%config%atom(T%config%num_atoms),   Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the atomic arrays of input model&
                                & (subroutine allocate_input_atomic_arrays). Model is too large'
      Call error_stop(message)
    End If

  End Subroutine allocate_input_atomic_arrays

  Subroutine allocate_arrays_vasp_trajectory(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate atomic arrays for trajectory in vasp format
    !
    ! author    - i.scivetti Apr 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),  Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(2)
    Character(Len=256)  :: message

    fail=0

    Allocate(T%config%element_file(T%config%list%num_elements), Stat=fail(1))
    Allocate(T%config%amount_file(T%config%list%num_elements),  Stat=fail(2))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for arrays in subroutine&
                                   &"allocate_arrays_vasp_trajectory". THIS SHOULD NOT HAPPEN! Please&
                                   & contact the developers.'
      Call error_stop(message)
    End If

  End Subroutine allocate_arrays_vasp_trajectory
  
  Subroutine allocate_species_arrays(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate species arrays as part of the model
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Class(model_type),  Intent(InOut)  :: T

    Integer(Kind=wi)    :: fail(1)
    Character(Len=256)  :: message

    fail=0

    Allocate(T%config%species(T%config%Nmax_species),        Stat=fail(1))
    If (Any(fail > 0)) Then
      Write (message,'(1x,1a)') '***ERROR: Allocation problems for the species of input model&
                                & (subroutine allocate_species_arrays). Maybe the model is too large.'
      Call error_stop(message)
    End If

  End Subroutine allocate_species_arrays

  Subroutine cleanup(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate variables
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Type(model_type) :: T

    If (Allocated(T%config%atom)) Then
      Deallocate(T%config%atom)
    End If 

    If (Allocated(T%config%species)) Then
      Deallocate(T%config%species)
    End If 

    If (Allocated(T%config%element_file)) Then
      Deallocate(T%config%element_file)
    End If 
    
    If (Allocated(T%config%amount_file)) Then
      Deallocate(T%config%amount_file)
    End If 

    If (Allocated(T%input_composition%tag)) Then
      Deallocate(T%input_composition%tag)
    End If 

    If (Allocated(T%input_composition%N0)) Then
      Deallocate(T%input_composition%N0)
    End If 

    If (Allocated(T%input_composition%element)) Then
      Deallocate(T%input_composition%element)
    End If 
    
  End Subroutine cleanup 

  Subroutine atomistic_model(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subrotine to:
    ! - evaluate if there are changing species and redefine tags
    ! - identify the set of atoms that belong to the species to be monitored
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),    Intent(InOut) :: model_data
    Integer(Kind=wi),    Intent(In   ) :: frame

    Integer(Kind=wi) :: j
    
    If(model_data%change_chemistry%stat) Then 
      ! Reset the condition to search search and identify sites 
      model_data%config%atom(:)%identified=.False.
      model_data%config%atom(:)%Nbonds=0
      model_data%number_species_found=0
      If (frame == 1) Then 
         Call identify_initial_chemistry(model_data, frame)
      Else
         Call search_chemistry_changes(model_data, frame)   
      End If
      ! Track the targeted species 
      Call tracking_chemistry(model_data, frame)
    End If
    
    ! Reset identified according to what has changed or not
    Do j=1, model_data%config%num_atoms 
      If (model_data%config%atom(j)%tag==model_data%config%atom(j)%tag_0) Then
        model_data%config%atom(j)%identified=.False. 
      Else
        model_data%config%atom(j)%identified=.True.
      End If
    End Do

    If (model_data%config%monitored_species%fread) Then
      Call check_monitored_species(model_data, frame)
    End If

  End Subroutine atomistic_model

  Subroutine check_monitored_species(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify molecular species to be monitored along the 
    ! trajectory
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: frame
  
    Integer(Kind=wi) :: j, k
    
    Do j = 1, model_data%config%Nmax_species
      k = model_data%config%species(j)%list(1)
      If(model_data%config%atom(k)%tag==model_data%species_definition%reference_tag%type) Then
        If (model_data%config%species(j)%alive) Then
          model_data%config%species(j)%alive=.True.
        Else
          If (model_data%species_definition%atoms_per_species /= 1) Then 
            Call identify_species(model_data, k, j, frame)
          Else
            model_data%config%atom(k)%identified=.True.
          End If
          model_data%config%species(j)%alive=.True.
        End If
      Else
        model_data%config%species(j)%alive=.False.
      End If
    End Do

  End Subroutine check_monitored_species
  
  Subroutine identify_species(model_data, kin, jin, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to find the all the neighbours-bonded atoms that form a species 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: kin
    Integer(Kind=wi),  Intent(In   ) :: jin
    Integer(Kind=wi),  Intent(In   ) :: frame

    Integer(Kind=wi) :: m, k, ka, il
    Integer(Kind=wi) :: knn, nacc, num_nn, num_tot, nlist, knni
    Logical             :: loop, loop_sp, is_nn, is_list
    Real(Kind=wp)       :: a(3), b(3)
    Real(Kind=wp)       :: dist
    Character(Len=256)  :: messages(6)

    Real(Kind=wp)       :: bond_cutoff
    Integer(Kind=wi), Dimension(max_neighbours) :: list_nn, list_nn_next

    bond_cutoff=model_data%species_definition%bond_cutoff%value    
    
    num_nn=1
    list_nn(1)=kin

    loop=.True. 
    num_tot=0
    list_nn_next=0 
    nlist=1

    Do While (loop)
      nacc=0
      Do knn=1, num_nn
        ka=list_nn(knn)
        Do k=1, model_data%config%num_atoms
          loop_sp=.True.
          m=1
          Do While (m <= model_data%species_definition%num_components .And. loop_sp)
            If ((model_data%config%atom(k)%element==model_data%species_definition%element(m)) .And. k/=ka .And.&
               (.Not. model_data%config%atom(k)%identified) .And. (.Not. model_data%config%atom(ka)%identified) ) Then
               ! Compute distances
               a=model_data%config%atom(k)%r
               b=model_data%config%atom(ka)%r
               Call compute_distance_PBC(a, b, model_data%config%cell, model_data%config%invcell, dist)
               ! Check if distances are within a given cutoff
               If (dist<bond_cutoff) Then 
                 If (dist<min_intra) Then
                   Write (messages(1),'(a,i6)') '***ERROR in the configuration structure of frame: ', frame
                   Write (messages(2),'(3(a,i4),3a)') 'Intermolecular distance between atoms ', k, ' and ', ka, &
                                                    ' in unit ', jin, ' of species "',&
                                                    & Trim(model_data%species_definition%name%type),'"'
                   Write (messages(3),'(a,f4.2,a)') 'is shorter than the input minimum distance criteria of ', &
                                                   &  min_intra,&
                                                   & ' Angstrom for bonding. Please review the input geometry.'
                   Call info(messages,3)
                   Call error_stop(' ')
                 Else
                   loop_sp=.False.
                   is_nn=.True.
                   il=1
                   Do While (il <= max_neighbours .And. is_nn) 
                     If (k==list_nn_next(il)) Then
                       is_nn=.False.
                     End If
                     il=il+1
                   End Do 
                   If (is_nn) Then    
                     nacc=nacc+1
                     If (nacc==max_neighbours+1) Then
                        Write (messages(1),'(a,i4,3a,i6)') '***ERROR: Trouble to identify unit ', jin, ' of species "',&
                                                 & Trim(model_data%species_definition%name%type),&
                                                 & '" for the configuration of frame: ', frame    
                        Write (messages(2),'(a)') ' It is likely the species represents does not represent a molecule!&
                                                 & Please review the definition of "&monitored_species".'
                        Call info(messages, 2)
                        Call error_stop(' ')
             
                     End If
                     list_nn_next(nacc)=k
                     il=1
                     is_list=.True.
                     Do While (il <= model_data%species_definition%atoms_per_species .And. is_list) 
                       If (k==model_data%config%species(jin)%list(il)) Then
                         is_list=.False.
                       End If 
                       il=il+1
                     End Do
                     If (is_list) Then
                       nlist=nlist+1
                       model_data%config%species(jin)%list(nlist)=k
                     End If
                   End If
                 End If  
               End If
            End If
            m=m+1
          End Do
        End Do
        If (.Not. model_data%config%atom(ka)%identified) Then
          model_data%config%atom(ka)%identified=.True.
          num_tot=num_tot+1
          If (num_tot+nacc == model_data%species_definition%atoms_per_species) Then
            Do knni= 1, nacc
              model_data%config%atom(list_nn_next(knni))%identified=.True.
            End Do
            loop=.False.
            num_tot=model_data%species_definition%atoms_per_species
          End If
        End If
      End Do
      num_nn=nacc  
      list_nn=list_nn_next
      list_nn_next=0 
      If (nacc==0) loop=.False.
    End Do

    If (model_data%species_definition%atoms_per_species/=num_tot) Then
      Write (messages(1),'(a,i4,3a,i6)') '***ERROR: Trouble to identify unit ', jin, ' of the monitored species "',&
                               & Trim(model_data%species_definition%name%type), '" for frame: ', frame
      Write (messages(2),'(a)')  'Possible reasons:' 
      Write (messages(3),'(a)')  ' 1) highly unstable configuration'

      Write (messages(4),'(a,f4.2,a)')  ' 2) the value for the bond cutoff (', bond_cutoff, &
                                    & ' Angstrom) as specified in &monitored_species&
                                    & is not adequate. The user must adjust this value'

      Write (messages(5),'(a)') ' 3) the cell vectors are not strictly consistent with the geometry of the input structure.'
      Write (messages(6),'(a)') ' 4) wrong definition of the monitored species.&
                                & Check the "reference_tag" and the &atomic_components block.'
      Call info(messages, 6)
      If (model_data%chem%acceptor%info_exclude%fread)Then
        Write (messages(1),'(a)') ' 5) Inadequate use of the "exclude_pairs" directive. Try removing this directive.'  
        Call info(messages, 1)
      End If 
      If (model_data%change_chemistry%stat) Then
         Write (messages(1),'(a)') '************************************************************************************'
         Write (messages(2),'(a)') '*** IMPORTANT: also check the settings of the &search_chemistry block, particularly '
         Write (messages(3),'(a)') '***            the "cutoff" directive for the &bonding_criteria (might be too large)'
         Write (messages(4),'(a)') '************************************************************************************'
         Call info(messages, 4)
      End If
      
      Call error_stop(' ')
    End If

  End Subroutine identify_species
  

  Subroutine obtain_maximum_number_species(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Obtain the maximum amount of specis
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
  
    Integer(Kind=wi) :: j
    
    ! For relevent indexes, find the initial index list
    model_data%config%Nmax_species=0
    Do j = 1, model_data%config%num_atoms
      If(model_data%config%atom(j)%tag_0==model_data%species_definition%reference_tag%type) Then
       model_data%config%Nmax_species=model_data%config%Nmax_species+1 
      End If
    End Do 

  End Subroutine obtain_maximum_number_species
  
  Subroutine identify_monitored_indexes(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Find the initial index list
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
  
    Integer(Kind=wi) :: j, icount
  
    icount=0
    Do j = 1, model_data%config%num_atoms
      If(model_data%config%atom(j)%tag_0==model_data%species_definition%reference_tag%type) Then
       icount=icount+1 
       model_data%config%species(icount)%list(1)=j
      End If
    End Do 

  End Subroutine identify_monitored_indexes
  
  Subroutine tracking_chemistry(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Track and reorder indexes for a better visualization 
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi),    Intent(In   ) :: frame

    Integer(Kind=wi) :: i, j, indx_new
    Real(Kind=wp)    :: dist, dist_min
   Character(Len=256) :: messages(2)

    If (frame == 1) Then
      Do i = 1, model_data%chem%N0%value
        model_data%track_chem(i)%indx=model_data%chem%indx_new(i)
      End Do
    Else        
       Do i = 1, model_data%chem%N0%value
        dist_min=Huge(1.0_wp)
        Do  j = 1, model_data%chem%N0%value
          indx_new  = model_data%chem%indx_new(j)
          Call compute_distance_PBC(model_data%config%atom(indx_new)%r, model_data%track_chem(i)%r0,&
                                & model_data%config%cell, model_data%config%invcell, dist)  
          If (dist < dist_min) Then
             dist_min=dist
             model_data%track_chem(i)%indx=indx_new
          End If                
        End Do
      End Do
    End If

    Do i=1, model_data%chem%N0%value-1
      Do j=i+1, model_data%chem%N0%value
         If (model_data%track_chem(i)%indx == model_data%track_chem(j)%indx) Then
           Write (messages(1), '(1x,a)') '*** ERROR: problems to identify chemical species!'
           Write (messages(2), '(1x,a)') '    Please reduce the value of the cutoff directive in &bonding_criteria'
           Call info(messages, 2) 
           Call error_stop(' ')
         End If
      End Do
    End Do 

    ! Update indexes
    Do i = 1, model_data%chem%N0%value
      model_data%track_chem(i)%r=model_data%config%atom(model_data%track_chem(i)%indx)%r
      model_data%track_chem(i)%r0=model_data%config%atom(model_data%track_chem(i)%indx)%r
      model_data%track_chem(i)%tag=model_data%config%atom(model_data%track_chem(i)%indx)%tag
    End Do
    model_data%chem%indx_prev = model_data%chem%indx_new

  End Subroutine tracking_chemistry

  Subroutine identify_initial_chemistry(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify the set of atomic species from the information of input settings 
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi),    Intent(In   ) :: frame

    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256) :: messages(2)

    ! Initialise  
    Do k=1, max_neighbours 
      model_data%config%atom(:)%bonds(k)=0
    End Do

    Write (messages(1), '(1x,a,i5,a)') '*** ERROR: problems to identify the chemistry for the starting&
                                      & configuration (frame ', frame, ') of the recorded trajectory.'

    ! Check for bonds according to the criteria set in the &possible_extra_bonds block 
    If (model_data%extra_bonds%invoke%fread) Then
      Call initial_extra_bonding(model_data, 1)
    End If

    ! Find those sites with the type and amount of bonds as specified in the &search_chemistry block
    If (model_data%number_species_found /= model_data%chem%N0%value) Then
      i=1
      Do While (i <= model_data%config%num_atoms-1)
        j=i+1
        Do While (j <= model_data%config%num_atoms)
          Call find_species_bonds(model_data, i, j)
          j=j+1
        End Do
        i=i+1
      End Do
      
      Do i=1, model_data%config%num_atoms
        If (model_data%config%atom(i)%Nbonds == model_data%chem%bonds%N0%value) Then
           model_data%number_species_found = model_data%number_species_found+1 
           model_data%chem%indx_new(model_data%number_species_found) = i
           model_data%config%atom(i)%tag=Trim(model_data%config%atom(i)%tag)//'*'
           Do k=1, model_data%chem%bonds%N0%value
             m=model_data%config%atom(i)%bonds(k)
             model_data%config%atom(m)%tag=Trim(model_data%config%atom(m)%tag)//'*'       
           End Do
         End If
      End Do
    End If

    ! If there are less or more target species that those specified, complain and abort
    If (model_data%number_species_found /= model_data%chem%N0%value) Then
      Call info(messages, 1)                     
      Call error_chemistry(model_data, 'failed_initial_counting')                     
    End If

    ! Copy initial configuration for environment checking
    Do i=1, model_data%config%num_atoms
      model_data%config%atom(i)%r0=model_data%config%atom(i)%r
    End Do

    ! Reset bonds indexes for retag  
    Do i=1, model_data%config%num_atoms
      Do k=1, max_neighbours 
        model_data%config%atom(i)%bonds0(k)=model_data%config%atom(i)%bonds(k)
      End Do
      model_data%config%atom(i)%Nbonds0=model_data%config%atom(i)%Nbonds
    End Do

  End Subroutine identify_initial_chemistry

  Subroutine error_chemistry(model_data, condition) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Error to identify the chemistry
    ! 
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Character(Len=*), Intent(In   ) :: condition

    Character(Len=256) :: messages(2)

    Write (messages(1), '(1x,a)') 'Please review the settings. Possible reasons:'
    Call info (messages, 1)

    If (condition=='failed_initial_counting') Then
      
      If (model_data%number_species_found > model_data%chem%N0%value) Then
        Write (messages(1), '(1x,a)') '- value for the "cutoff" directive of the "&bonding_criteria" sub-block is too large'
        Call info (messages, 1)
      Else
        Write (messages(1), '(1x,a)') '- value for the cutoff directive of the "&bonding_criteria" sub-block&
                                       & is not large enough'
        Call info (messages, 1)
        If (model_data%extra_bonds%invoke%fread) Then
          Write (messages(1), '(1x,a)') '- the information provided in the "&possible_extra_bonds" is not correct (wrong& 
                                       & species involved and/or incorrect cutoff values'
          Call info (messages, 1)
        Else
          Write (messages(1), '(1x,a)') '- the system could also form bonds with other species. The user&
                                       & should consider using the "&possible_extra_bonds" block'
          Call info (messages, 1)
        End If        
      End If
      Write (messages(1), '(1x,a)') '- the value assigned to "total_number" in the  "&search_chemistry" block is incorrect'
      Write (messages(2), '(1x,a)') '- the specification for "number_of_bonds" and/or "only_element" in the&
                                    & "&bonding_criteria" block is wrong'
      Call info (messages, 2)

    Else If (condition=='failed_counting') Then
      If (model_data%number_species_found > model_data%chem%N0%value) Then
        Write (messages(1), '(1x,a)') '- the "cutoff" directive of the "&acceptor_criteria" sub-block is too large'
      Else  
        Write (messages(1), '(1x,a)') '- the "cutoff" directive of the "&acceptor_criteria" sub-block is too low'
      End If
      If (model_data%chem%acceptor%info_exclude%fread) Then
        Write (messages(1), '(1x,a)') '- inadequate use of the "exclude_pairs" directive. Try by removing this directive.'
        Call info (messages, 1)
      End If
    End If

    If (model_data%input_geometry_format%type=='xyz') Then
       Write (messages(1),'(1x,a)') '- vectors defined in the &simulation_cell block do not correspond to the generated trajectory'
       Call info (messages, 1)
    End If

    Write (messages(1), '(1x,a)') '- the atomic configuration is highly unstable (non-equilibrated trajectories)'
    Call info (messages, 1)
    Call error_stop(' ')

  End Subroutine error_chemistry

  Subroutine search_chemistry_changes(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Identify changes of chemistry
    ! 
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),   Intent(InOut) :: model_data
    Integer(Kind=wi),   Intent(In   ) :: frame
  
    Integer(Kind=wi)  :: i, j, l, m
    Integer(Kind=wi)  :: k1, k2
    Integer(Kind=wi)  :: n1, n2, ni
    Integer(Kind=wi)  :: icount, ihit
    Integer(Kind=wi)  :: icc, indx
    Integer(Kind=wi)  :: N0_counted, N0_missed

    Logical            :: match, fenv, fprev
    Character(Len=256) :: messages(2)
  
    If (model_data%chem%acceptor%check) Then
      Call check_environment_size(model_data, frame) 
    End If

    ! Initialize arrays
    model_data%chem%acceptor%accum_indxs=0
    model_data%extra_bonds%flag=.True.
    Do j=1, max_neighbours 
      model_data%config%atom(:)%bonds(j)=0
    End Do

    ! First check if any of the previous indexes keeps or form extra bonds
    If (model_data%extra_bonds%invoke%fread) Then
      icount=1
      Do l=1, model_data%chem%N0%value
        match=.False.
        i=model_data%chem%indx_prev(l)
          j=1
          Do While (j <= model_data%config%num_atoms .And. (.Not. match))
            If (i/=j) Then
              Call search_extra_bonding(model_data, i, j, icount, match)
            End If
            j=j+1
          End Do

          If (match) Then
            model_data%extra_bonds%flag(l)=.False.
            icount=icount+1
            model_data%number_species_found=model_data%number_species_found+1
          End If
      End Do
    End If
    
    ! Now check the environment around each of the previously found species
    N0_counted=0
    N0_missed=0
    Do l=1, model_data%chem%N0%value
      If (model_data%extra_bonds%flag(l)) Then
        Call define_environment(model_data, l)
        k2=1
        ihit=0
        fprev=.True.
        Do While (k2 <= model_data%chem%acceptor%num_nn .And. fprev)
          m=model_data%chem%acceptor%list(k2)
          fenv=.True.
          ! Check if a species of each environment has formend a new extra bond
          If (model_data%extra_bonds%invoke%fread) Then
            match=.False.
            j=1
            Do While (j <= model_data%config%num_atoms .And. (.Not. match))
              If (m/=j) Then
                If ((.Not. model_data%config%atom(j)%identified) .And. (.Not. model_data%config%atom(m)%identified))  Then
                  k1=model_data%number_species_found+1      
                  Call search_extra_bonding(model_data, m, j, k1, match)
                End If
              End If
              j=j+1
            End Do
            If (match) Then
              fenv=.False.      
              model_data%number_species_found=model_data%number_species_found+1
              ihit=ihit+1  
            End If
          End If
  
          ! Now check the chemical bonds as spcified in the &bonding_criteria block  
          If (fenv) Then
            j=1
            Do While (j <= model_data%config%num_atoms)
              If (m /= j) Then
                Call find_species_bonds(model_data, m, j)
              End If
              j=j+1
            End Do
            ! Count those sites of the environment that become a target species with N0 bonds
            If (model_data%config%atom(m)%Nbonds == model_data%chem%bonds%N0%value) Then
              Do n1=1, model_data%chem%bonds%N0%value
                n2=model_data%config%atom(m)%bonds(n1)
                !model_data%config%atom(n2)%tag=Trim(model_data%config%atom(n2)%tag)//'*'       
              End Do
              If (k2/=1) Then
                If ( ( .Not. (Any(model_data%chem%acceptor%accum_indxs==m)) ) .And. &
                     ( .Not. (Any(model_data%chem%indx_prev==m)) ) ) then
                  ihit=ihit+1
                  N0_counted=N0_counted+1
                  model_data%chem%acceptor%accum_indxs(N0_counted)=m
                End If
              Else
                If (.Not. (Any(model_data%chem%acceptor%accum_indxs==m)) ) then
                  ihit=ihit+1
                  N0_counted=N0_counted+1
                  model_data%chem%acceptor%accum_indxs(N0_counted)=m
                  fprev=.False.
                End If
              End If
            End If    
          End If    
          k2=k2+1 
        End Do 
        ! If none of the sites within the environment has N0 bonds, we allocate the centre of the environment in a list,
        ! in case there are missing species at the end
        If (ihit == 0) Then
           icc=0     
           Do indx = 1, N0_counted     
             If ( .Not. (model_data%chem%acceptor%accum_indxs(indx)==model_data%chem%acceptor%list(1))) then     
               icc=icc+1
             End If
           End Do
           If (icc == N0_counted) Then
             N0_missed=N0_missed+1
             model_data%chem%acceptor%missed(N0_missed)=model_data%chem%acceptor%list(1)
           End If    
        End If  
      End If
    End Do

    ! Add to the list those sites that become/preserve a target species
    If (model_data%number_species_found < model_data%chem%N0%value) Then
      If (N0_counted/=0) Then
        model_data%chem%acceptor%order_indxs(1)=model_data%chem%acceptor%accum_indxs(1)
        icount=1
        Do j = icount+1, N0_counted
          match=.True. 
          i=1
          Do While (i <= icount .And. match)
             If (model_data%chem%acceptor%accum_indxs(j) == model_data%chem%acceptor%order_indxs(i)) Then 
                match=.False.
             End If  
             i=i+1   
          End Do
          If (match) Then
            icount=icount+1 
            model_data%chem%acceptor%order_indxs(icount)=model_data%chem%acceptor%accum_indxs(j)     
          End If        
        End Do
        Do i=1, icount
          model_data%chem%indx_new(i+model_data%number_species_found) = model_data%chem%acceptor%order_indxs(i)
          ! Retag
          ni=model_data%chem%acceptor%order_indxs(i)
          model_data%config%atom(ni)%tag=Trim(model_data%config%atom(ni)%tag)//'*'
          Do n1=1, model_data%config%atom(ni)%Nbonds 
            n2=model_data%config%atom(ni)%bonds(n1)
            model_data%config%atom(n2)%tag=Trim(model_data%config%atom(n2)%tag)//'*'       
          End Do
        End Do 
        model_data%number_species_found=model_data%number_species_found+icount
      End If
    End If        
    
    ! If it is not enough, one shall consider including the centre atom of those
    ! environments that had NONE species with N0 bonds 
    If (model_data%number_species_found < model_data%chem%N0%value) Then
      If (N0_missed /=0 ) Then
        Do i=1, N0_missed
          model_data%chem%indx_new(i+model_data%number_species_found) = model_data%chem%acceptor%missed(i)
          ! Retag
          ni = model_data%chem%acceptor%missed(i)
          model_data%config%atom(ni)%tag=Trim(model_data%config%atom(ni)%tag)//'*'
          model_data%config%atom(ni)%Nbonds=model_data%config%atom(ni)%Nbonds0
          Do n1=1, model_data%config%atom(ni)%Nbonds 
            model_data%config%atom(ni)%bonds(n1)=model_data%config%atom(ni)%bonds0(n1)
            n2=model_data%config%atom(ni)%bonds(n1)
            model_data%config%atom(n2)%tag=Trim(model_data%config%atom(n2)%tag)//'*'       
          End Do
        End Do
        model_data%number_species_found=model_data%number_species_found+N0_missed 
      End If  
    End If

    ! If there was no success, complain and abort
    If (model_data%number_species_found /=  model_data%chem%N0%value) Then
      Write (messages(1), '(1x,a,i6, a)') '*** ERROR: problems to identify the chemistry for the configuration ',&
                                         & frame,' of the trajectory.'
      Call info(messages, 1)
      Call error_chemistry(model_data, 'failed_counting')                     
    End If

    ! Further check to make sure there was no double counting, only for extreme cases  
    icount=0
    Do k1=1, model_data%chem%N0%value
      i=model_data%chem%indx_new(k1)
      Do k2= 1,  model_data%config%atom(i)%Nbonds 
        icount=icount+1
        model_data%chem%bonds%list(icount)=model_data%config%atom(i)%bonds(k2)
      End Do
    End Do  

    match=.False.
    Do k1= 1, icount
      Do k2= 1, icount
        If(k1/=k2 .And. (model_data%chem%bonds%list(k1)==model_data%chem%bonds%list(k2))) Then
          match=.True.
        End If        
      End Do
    End Do

    If(match)Then
      Write (messages(1), '(1x,a,i6, a)') '*** ERROR: problems to identify the chemistry for the configuration ',&
                                         & frame,' of the trajectory.'
      Call info(messages, 1)
      Call error_chemistry(model_data, 'failed_counting')                     
    End If         

    ! Reset bonds indexes for retag  
    Do j=1, max_neighbours 
      model_data%config%atom(:)%bonds0(j)=model_data%config%atom(:)%bonds(j)
    End Do
    model_data%config%atom(:)%Nbonds0=model_data%config%atom(:)%Nbonds

  End Subroutine search_chemistry_changes
  
  Subroutine find_species_bonds(model_data, i, j)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Determine is there is a bond between species i and j
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: i
    Integer(Kind=wi),  Intent(In   ) :: j
    
    Integer(Kind=wi) :: m, indx, N0
    Real(Kind=wp)    :: dist
    Logical          :: match1, match2 
  
    match1=.False.
    match2=.False.
    m=1

    Do While (m <= model_data%chem%acceptor%N0_incl .And. (.Not. match1))
      If(model_data%config%atom(i)%tag==model_data%chem%acceptor%tg_incl(m) .Or.&
         model_data%config%atom(j)%tag==model_data%chem%acceptor%tg_incl(m)) Then
        match1=.True. 
        If((model_data%config%atom(i)%tag==model_data%chem%acceptor%tg_incl(m)) .And. &
           (.Not. model_data%config%atom(i)%identified))Then
          indx=i
          If (((model_data%config%atom(j)%element==model_data%chem%bonds%species%type) .And. &
              (.Not. model_data%config%atom(j)%identified))) Then
            match2=.True.  
          End If
        
        Else If((model_data%config%atom(j)%tag==model_data%chem%acceptor%tg_incl(m)) .And.&
                (.Not. model_data%config%atom(j)%identified)) Then
          indx=j
          If ((model_data%config%atom(i)%element==model_data%chem%bonds%species%type) .And. &
              (.Not. model_data%config%atom(i)%identified)) Then
            match2=.True.  
          End If
          
        End If  
      End If
      m=m+1
    End Do

    If (match1 .And. match2) Then
      Call compute_distance_PBC(model_data%config%atom(i)%r, model_data%config%atom(j)%r,&
                              & model_data%config%cell, model_data%config%invcell, dist)
      If (dist < model_data%chem%bonds%cutoff%value ) Then
        model_data%config%atom(indx)%Nbonds=model_data%config%atom(indx)%Nbonds+1 
        N0=model_data%config%atom(indx)%Nbonds
        If(model_data%config%atom(i)%element==model_data%chem%bonds%species%type) Then
          model_data%config%atom(i)%identified=.True.
          model_data%config%atom(j)%bonds(N0)=i
        Else If(model_data%config%atom(j)%element==model_data%chem%bonds%species%type) Then
          model_data%config%atom(j)%identified=.True.
          model_data%config%atom(i)%bonds(N0)=j
        End If
      End If
    End If

  End Subroutine find_species_bonds
 
  Subroutine check_environment_size(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if the defined environment criteria is suitable to the system 
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: frame

    Integer(Kind=wi)   :: l, icount
    Character(Len=256) :: messages(3)

    icount=0 
    Do l=1, model_data%chem%N0%value
      Call define_environment(model_data, l)
      If (model_data%chem%acceptor%num_nn/=1) Then
        icount=icount+1      
      End If        
    End Do

    If (icount==0) Then
      Write (messages(1), '(1x,a,i6,a)') '*** ERROR: problems to identify atomic enviroments for the requested initial frame&
                                        & (configuration ' , frame-1, ' of the trajectory).' 
      Write (messages(2), '(1x,a)') 'Either the enviroment cutoff set is too small or the atomic configuration&
                                    & is unrealistic to identify the chemical enviroment with the defined criteria.'
      Write (messages(3), '(1x,a)') 'Please review the settings of the &search_chemistry block and this particular&
                                    & configuration in the TRAJECTORY file' 
      Call info (messages, 3)
      Call error_stop(' ')
    End If
   
    ! Turn off checking from now on
    model_data%chem%acceptor%check=.False. 

  End Subroutine check_environment_size

  Subroutine define_environment(model_data, kindx)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Define the set of non-bonded atoms (environment) around a particular atom
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: kindx

    Integer(Kind=wi)   :: i, j, k, m, icount
    Real(Kind=wp)      :: dist
    Logical            :: match_self, match_i, match_j
    
    Character(Len=256) :: message
    
    icount=0
    model_data%chem%acceptor%num_search=1
    model_data%chem%acceptor%num_nn=1
    i=model_data%chem%indx_prev(kindx)
    model_data%chem%acceptor%list(model_data%chem%acceptor%num_nn)=i
        
    j=1
    Do While (j <= model_data%config%num_atoms)
      If (i/=j) Then
        match_self=.False.
        If (model_data%chem%acceptor%info_exclude%fread)Then
          m=1
          Do While (m <= model_data%chem%acceptor%N0_excl .And. (.Not. match_self))
            If(model_data%config%atom(i)%tag==model_data%chem%acceptor%tg_excl(m) .And.&
               model_data%config%atom(j)%tag==model_data%chem%acceptor%tg_excl(m)) Then
               match_self=.True.
            End If
            m=m+1
          End Do
        End If
        
        If (.Not. match_self) Then
          match_i=.False.
          match_j=.False.
          Do k = 1, model_data%chem%acceptor%N0_incl
            If (model_data%config%atom(i)%tag==model_data%chem%acceptor%tg_incl(k)) Then  
              match_i=.True.
            End If
            If (model_data%config%atom(j)%tag==model_data%chem%acceptor%tg_incl(k)) Then  
              match_j=.True.
            End If
          End Do
          If(match_i .And. match_j) Then
            If (model_data%chem%acceptor%check) Then
              Call compute_distance_PBC(model_data%config%atom(i)%r0, model_data%config%atom(j)%r0,&
                                    & model_data%config%cell, model_data%config%invcell, dist)
            Else                
              Call compute_distance_PBC(model_data%config%atom(i)%r, model_data%config%atom(j)%r,&
                                    & model_data%config%cell, model_data%config%invcell, dist)
            End If                
            If (dist < model_data%chem%acceptor%cutoff%value) Then
              icount=icount+1
              model_data%chem%acceptor%num_nn=model_data%chem%acceptor%num_nn+1
              If (model_data%chem%acceptor%num_nn > max_components) Then
                Write (message,'(1x,a)') '***ERROR - The number of enviroment neighbours has reached a default maximum.' 
                Call info (message, 1)
                Write (message,'(4x,a)') 'Either review the settings of "&acceptor_criteria" or increase the&
                                      & "max_components" variable in the code and recompile'
                Call info (message, 1)                      
                Call error_stop(' ')
                End If
              model_data%chem%acceptor%list(model_data%chem%acceptor%num_nn)=j
            End If
          End If
        End If
      End If
      j=j+1
    End Do
    
  End Subroutine define_environment

  Subroutine initial_extra_bonding(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the formation of bonding for those species defined
    ! in &possible_extra_bonds.
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: frame
    
    Integer(Kind=wi)   :: i, j, l
    Logical            :: match 
    Character(Len=256) :: messages(3)
   
    l=1
    i=1 
    Do While (i <= model_data%config%num_atoms-1)
      match=.False.
      j=i+1
      Do While (j <= model_data%config%num_atoms .And. (.Not. match))
        Call search_extra_bonding(model_data, i, j, l, match)
        If (match) Then
          l=l+1      
          model_data%number_species_found=model_data%number_species_found+1  
        End If
        j=j+1
      End Do
      i=i+1
    End Do

    If (model_data%number_species_found > model_data%chem%N0%value) Then
      Write (messages(1), '(1x,a,i6,a)') '*** ERROR: problems to identify extra bonds (&possible_extra_bonds sub-block)&
                                    & for configuration ', frame, ' of the trajectory.' 
      Write (messages(2), '(1x,a)') 'Either the distance cutoffs are exeedingly large or the atomic configuration&
                                    & is highly unstable and unrealistic.'
      Write (messages(3), '(1x,a)') 'Please review the settings of the &possible_extra_bonds block and this particular&
                                    & configuration in the TRAJECTORY file' 
      Call info (messages, 3)
      Call error_stop(' ')
    End If
    
  End Subroutine initial_extra_bonding

  Subroutine search_extra_bonding(model_data, i, j, l, match)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the formation of bonding for those species defined
    ! in &extra_bonding.
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: i 
    Integer(Kind=wi), Intent(In   ) :: j
    Integer(Kind=wi), Intent(In   ) :: l
    Logical,          Intent(InOut) :: match

    Integer(Kind=wi)   :: k
    Real(Kind=wp)      :: dist
    Logical            :: f1, f2

    Do k = 1, model_data%extra_bonds%N0
      If(.Not. match) Then
        f1=(model_data%config%atom(i)%tag==model_data%extra_bonds%tg1(k) .And.&
            model_data%config%atom(j)%tag==model_data%extra_bonds%tg2(k)) 
        f2=(model_data%config%atom(j)%tag==model_data%extra_bonds%tg1(k) .And.&
            model_data%config%atom(i)%tag==model_data%extra_bonds%tg2(k)) 

        If(f1 .Or. f2) Then
          Call compute_distance_PBC(model_data%config%atom(i)%r, model_data%config%atom(j)%r,&
                                   & model_data%config%cell, model_data%config%invcell, dist)
          If (dist < model_data%extra_bonds%bond(k)%value) Then
             match=.True.
             model_data%config%atom(i)%tag=Trim(model_data%config%atom(i)%tag)//'*'
             model_data%config%atom(j)%tag=Trim(model_data%config%atom(j)%tag)//'*'
          End If
        End If   

        If (match) Then
          If((model_data%config%atom(i)%element==model_data%chem%bonds%species%type) .And. &
             (.Not. model_data%config%atom(i)%identified)) Then
             model_data%chem%indx_new(l) = j
             model_data%config%atom(i)%identified=.True.
             model_data%config%atom(j)%bonds(1)=i
             model_data%config%atom(j)%Nbonds=1
          Else If((model_data%config%atom(j)%element==model_data%chem%bonds%species%type) .And. &
             (.Not. model_data%config%atom(j)%identified)) Then
             model_data%chem%indx_new(l) = i
             model_data%config%atom(j)%identified=.True.
             model_data%config%atom(i)%bonds(1)=j
             model_data%config%atom(i)%Nbonds=1
          End If
        End If
      End If
    End Do
          
  End Subroutine search_extra_bonding
   
  Subroutine read_model(files, model_data, frame, ensemble)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read input models from TRAJECTORY file 
    ! Format Options: 
    ! - xyz 
    ! - vasp
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: frame
    Character(Len=*),  Intent(In   ) :: ensemble

    If (model_data%input_geometry_format%type == 'xyz') Then
      Call read_input_xyz_format(files, model_data)
    Else If (model_data%input_geometry_format%type == 'vasp') Then
      Call read_input_vasp_format(files, model_data, frame, ensemble)
    End If

    ! Copy tags to tag_0
    model_data%config%atom(:)%tag_0=model_data%config%atom(:)%tag
    
  End Subroutine read_model
  
  Subroutine read_input_xyz_format(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read file INPUT_STRUCTURE according to xyz format 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Integer(Kind=wi) :: io, iunit
    Integer(Kind=wi) :: i, j, k, m
    Character(Len=256)    :: title
    Character(Len=256)    :: messages(5)
    
    Character(Len=256)    :: set_error, error_format(10), input_file

    input_file = Trim(files(FILE_TRAJECTORY)%filename) 
    iunit=files(FILE_TRAJECTORY)%unit_no

    set_error = '***ERROR in file '//Trim(input_file)//' (inconsistency with xyz format):'
    error_format(1) ='In ALC_TRAJECTORY, the structure of the '//Trim(input_file)//' file in "xyz" format must be:'
    error_format(2) = ' '
    error_format(3) = 'Nat (number of atoms)'
    error_format(4) = 'Description from MD run such as time and/or energy'
    error_format(5) = 'Element_1         X_1      Y_1      Z_1'
    error_format(6) = 'Element_2         X_2      Y_2      Z_2'
    error_format(7) = '...........       .....    .....    .....'
    error_format(8) = 'Element_Nat       X_Nat    Y_Nat    Z_Nat'
    error_format(9) = ' ' 
    error_format(10) = 'Please check consistency between the structure of the file and&
                       & directive "input_geometry_format".'

    ! Start reading the file
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! Read number of atoms 
    Read (iunit, Fmt=*, iostat=io) model_data%config%num_atoms
    If (io/=0) Then
      Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the number of atoms, which must be&
                               & specified in the first line.'
      Call info(messages,1)
      Call info(error_format,10)
      Call error_stop(' ')
    End If    

    If (model_data%config%num_atoms /= model_data%input_composition%numtot) Then
     Write (messages(1),'(3a)') 'ERROR***: Inconsistency between the number of atoms in file ', Trim(input_file),&
                             &' and the amount of atoms defined in &input_composition. Please check.'   
     Call error_stop(messages(1))
    End If

    ! Read title
    Read (iunit, Fmt=*, iostat=io) title

    ! Allocate atomic arrays for the input model
    If (.Not. model_data%config%allocated_model_geo) Then
      Call model_data%init_atomic_arrays()
      model_data%config%allocated_model_geo=.True.
    End If
    
    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%config%num_atoms)
      If (model_data%input_composition%N0(j)/=0) Then
        i=i+1
        Read (iunit, Fmt=*, iostat=io) model_data%config%atom(i)%element, (model_data%config%atom(i)%r(m), m=1,3)
        If (io/=0) Then
          If (i== model_data%config%num_atoms) Then
            Write (messages(1),'(2a)') Trim(set_error), ' Missing input coordinates somewhere in the list.'
          Else
            Write (messages(1),'(2a,i5)') Trim(set_error), ' Wrong specification for atom', i
          End If
          Call info(messages, 1)
          Call info(error_format,10)
          Call error_stop(' ')
        Else
          If (Trim(model_data%config%atom(i)%element)/=Trim(model_data%input_composition%element(j))) Then
            Write (messages(1),'(a,i5,5a)') '***ERROR: Chemical element of atom ', i, &
                                          & ' is set to "', Trim(model_data%config%atom(i)%element), '", but& 
                                          & according to the definition of &input_composition it must be "',&
                                          & Trim(model_data%input_composition%element(j)), '".'
            Write (messages(2),'(3a)')  'Check the consistency between the data in ', Trim(input_file),&
                                          & ' and &input_composition.'
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
        k=k+1
        model_data%config%atom(i)%tag=model_data%input_composition%tag(j)
        If (k==model_data%input_composition%N0(j)) Then
          k=0
          j=j+1
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

  End Subroutine read_input_xyz_format
  
  Subroutine read_input_vasp_format(files, model_data, frame, ensemble)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read atomic models according to VASP/POSCAR format 
    ! Adapted from ALC_EQCM
    !
    ! author    - i.scivetti April 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(model_type),  Intent(InOut) :: model_data
    Integer(Kind=wi),  Intent(In   ) :: frame
    Character(Len=*),  Intent(In   ) :: ensemble

    Character(Len=256)  :: title, word, word_frame
    Character(Len=256)  :: messages(5)
    Character(Len=  2)  :: buffer
    Character(Len=256)  :: input_file, set_error, webpage, error_ensemble
    Integer(Kind=wi)    :: io, iunit
    Integer(Kind=wi)    :: i, j, k, m
    Logical             :: loop, error_vasp, error_elements
    Real(Kind=wp)       :: v_cart(3)

    error_vasp = .False.
    input_file = Trim(files(FILE_TRAJECTORY)%filename) 
    iunit=files(FILE_TRAJECTORY)%unit_no
    Write(word_frame,*) frame

    set_error = '***ERROR in file '//Trim(input_file)//' (inconsistency with POSCAR format) at frame '&
               &//Trim(Adjustl(word_frame))//'. '
    webpage='For correct format and details see the following link: https://www.vasp.at/wiki/index.php/POSCAR'
    error_ensemble='IMPORTANT: are you sure that the '//Trim(input_file)//' file corresponds to the "'&
                  &//Trim(ensemble)//'" ensemble?'
    
    If (.Not. model_data%config%allocated_model_geo) Then
      ! Determine vasp elements and amounts of atoms from the data in &input_composition and &block_species_components
      model_data%config%list%num_elements=1 
      j=1
      model_data%config%list%N0(j)=model_data%input_composition%N0(j)
      model_data%config%list%element(j)=model_data%input_composition%element(j)   
      buffer=model_data%input_composition%element(j)
      Do i = 2,  model_data%input_composition%atomic_species
        If (Trim(model_data%input_composition%element(i)) /= Trim(buffer)) Then
          buffer=model_data%input_composition%element(i)
          If (model_data%input_composition%N0(i)/=0) Then
            model_data%config%list%num_elements=model_data%config%list%num_elements+1
            j=j+1
            model_data%config%list%N0(j)=model_data%input_composition%N0(i)
            model_data%config%list%element(j)=model_data%input_composition%element(i)
          End If
        Else
          model_data%config%list%N0(j)=model_data%config%list%N0(j)+model_data%input_composition%N0(i)
        End If
      End Do 
      
      model_data%config%num_atoms=0
      Do i=1, model_data%config%list%num_elements
         model_data%config%num_atoms=model_data%config%num_atoms + model_data%config%list%N0(i)
      End Do
     ! Allocate atomic arrays for the input model
      Call model_data%init_atomic_arrays() 
      ! Allocate vasp related lists
      Call model_data%init_vasp_arrays()
      model_data%config%allocated_model_geo=.True.
    End If

    If (Trim(ensemble)=='npt' .Or. frame==1) Then
        error_elements=.True.
        Read (iunit, Fmt=*, iostat=io) title

      ! Read scale factor
      Read (iunit, Fmt=*, iostat=io) model_data%config%scale_factor_vasp
      If (io/=0) Then
        Write (messages(1),'(2a)') Trim(set_error), ' Scale factor for cell vectors is invalid.'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call info(error_ensemble,1)
        Call error_stop(' ')
      End If 
  
      Do i= 1, 3
        Read (iunit, Fmt=*, iostat=io) (model_data%config%cell(i,j), j=1,3)
        If (io/=0) Then
          Write (messages(1),'(2a,i1,a)') Trim(set_error), ' Definition for cell vector ', i, ' is incorrect.'
          Write (messages(2),'(a)') webpage
          Call info(messages,2)
          Call info(error_ensemble,1)
          Call error_stop(' ')
        End If 
      End Do
  
      If (error_elements) Then
        Read (iunit, Fmt=*, iostat=io) (model_data%config%element_file(i) , i=1, model_data%config%list%num_elements)
        If (io/=0) Then
          Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the list of atomic species.'
          Write (messages(2),'(a)') webpage
          Call info(messages,2)
          Call info(error_ensemble,1)
          Call error_stop(' ') 
        End If
      End If
  
      Do i=1, model_data%config%list%num_elements
        loop=.True.
        k=1
        Do While (k <= NPTE .And. loop)
          If (chemsymbol(k)==model_data%config%element_file(i)) Then
            loop=.False.
          End If
          k=k+1
        End Do
        If (loop) Then
          Write (messages(1),'(4a,i2,a)') Trim(set_error), ' Chemical element "' , Trim(model_data%config%element_file(i)), &
                                         & '" defined for component ',  i, ' of the list does not correspond to an element&
                                         & of the Periodic Table. Please use a valid element.'
          Write (messages(2),'(a)')   webpage
          Write (messages(3),'(3a)') 'IMPORTANT: The user should also check for inconsistencies between the list/number of atoms&
                                  & in file ', Trim(input_file),' and the settings in &Block_input_composition.&
                                  & Have you missed info?'
          Call info(messages,3)
          Call info(error_ensemble,1)
          Call error_stop(' ')
        End If
  
        If (Trim(model_data%config%element_file(i)) /= Trim(model_data%config%list%element(i))) Then
          error_vasp=.True.
        End If 
      End Do
  
      Read (iunit, Fmt=*, iostat=io) (model_data%config%amount_file(i) , i=1, model_data%config%list%num_elements)
      If (io/=0) Then
        Read (iunit, Fmt=*, iostat=io) (model_data%config%amount_file(i) , i=1, model_data%config%list%num_elements)
        If (io/=0) Then
          Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the list with the number of atoms&
                                                 & for each atomic species'
          Write (messages(2),'(a)') webpage
          Call info(messages,2)
          Call info(error_ensemble,1)
          Call error_stop(' ')
        End If
      End If
  
      Write (messages(1),'(3a)') '***ERROR: Inconsistent between the list of atomic species in ',  Trim(input_file),&
                               & ' (VASP format) and the data specified in &input_composition.'
      Write (messages(2),'(a)') 'The list of elements and amount of atoms per element resulting from the definition&
                               & in &input_composition should be:'
            
      Write (messages(3),'(*(4x,a2))') (model_data%config%list%element(j), j=1, model_data%config%list%num_elements)
      Write (messages(4),'(*(1x,i5))') (model_data%config%list%N0(j), j=1, model_data%config%list%num_elements)
      Write (messages(5),'(3a)') 'Which DO NOT MATCH the order/amount of atoms in file ', Trim(input_file), & 
                              & '. Please review the settings of &input_composition (or the geometry of the input structure)'
  
      Do i=1, model_data%config%list%num_elements
        If (model_data%config%amount_file(i) /= model_data%config%list%N0(i)) Then
          error_vasp=.True.
        End If
      End Do
  
      If (error_vasp) Then
        Call info(messages, 5)
        Call info(error_ensemble,1)
        Call error_stop(' ')
      End If
      
    End If

      Read (iunit, Fmt=*, iostat=io) word
      Call capital_to_lower_case(word)
      If (word(1:4)=='sele') Then
        model_data%config%selective_dyn=.True. 
      Else
        model_data%config%selective_dyn=.False. 
        Backspace iunit
      End If
  
      Read (iunit, Fmt=*, iostat=io) model_data%config%list%coord_type
      If (io/=0) Then
        Write (messages(1),'(2a)') Trim(set_error), ' Invalid specification for the type of coordinates. Valid options are either& 
                                               & "Cartesian" or "Direct"'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call info(error_ensemble,1)
        Call error_stop(' ')
      End If
      Call capital_to_lower_case(model_data%config%list%coord_type)   
    
      If (Trim(model_data%config%list%coord_type) /= 'cartesian' .And. &
        Trim(model_data%config%list%coord_type) /= 'direct') Then
        Write (messages(1),'(2a)') Trim(set_error), ' Wrong option for the type of coordinates. Valid options are either& 
                                               & "Cartesian" or "Direct"'
        Write (messages(2),'(a)') webpage
        Call info(messages,2)
        Call info(error_ensemble,1)
        Call error_stop(' ')
      End If
    
    ! Read atomic coordinates
    j=1; k=0; i=0
    Do While (i < model_data%config%num_atoms)
      If (model_data%input_composition%N0(j)/=0) Then
        i=i+1
        If (model_data%config%selective_dyn) Then
          Read (iunit, Fmt=*, iostat=io) (model_data%config%atom(i)%r(m), m=1,3), &
                                              (model_data%config%atom(i)%dynamics(m), m=1,3)
        Else 
          Read (iunit, Fmt=*, iostat=io) (model_data%config%atom(i)%r(m), m=1,3)
        End If
        If (io/=0) Then
          If (i== model_data%config%num_atoms) Then
            Write (messages(1),'(a)') Trim(set_error)//' Missing input coordinates somewhere in the list'
          Else
            Write (messages(1),'(2a,i5,a)') Trim(set_error), ' Wrong specification for the input coordinates of atom ', i
          End If
          Call info(messages, 1)
          Call info(error_ensemble,1)
          Call error_stop(' ')
        End If
        k=k+1
        model_data%config%atom(i)%tag=model_data%input_composition%tag(j)
        model_data%config%atom(i)%element=model_data%input_composition%element(j)
        If (k==model_data%input_composition%N0(j)) Then
          k=0
          j=j+1 
        End If
      Else
        k=0
        j=j+1
      End If
    End Do

    ! Transform using the scaling factor
    If (Trim(model_data%config%list%coord_type) == 'cartesian') Then
      Do i = 1 , model_data%config%num_atoms
        model_data%config%atom(i)%r=model_data%config%scale_factor_vasp*model_data%config%atom(i)%r
      End Do
    Else If (Trim(model_data%config%list%coord_type) == 'direct') Then
      Do i = 1, model_data%config%num_atoms 
        v_cart(1:3)=MatMul(model_data%config%atom(i)%r(1:3), model_data%config%cell(1:3,1:3))
        model_data%config%atom(i)%r=model_data%config%scale_factor_vasp*v_cart
      End Do
    End If
    
    model_data%config%cell= model_data%config%scale_factor_vasp * model_data%config%cell
    
  End Subroutine read_input_vasp_format
  

  Subroutine check_definition_bonds(model_data, frame)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Checks the definition of bond related quantities against the dimension of
    ! the simulation cell
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data
    Integer(Kind=wi), Intent(In   ) :: frame
  
    Character(Len=256) :: messages(2)  
    Integer(Kind=wi)   :: i, k
    Real(Kind=wp)      :: min_vector_length
    Real(Kind=wp)      :: length_from_volume
    Real(Kind=wp)      :: bond_criterion
    
    min_vector_length=Huge(1.0_wp)
    Do i = 1, 3
      If (model_data%config%cell_length(i)< min_vector_length) Then
         min_vector_length = model_data%config%cell_length(i)
      End If
    End Do
    
    length_from_volume=(model_data%config%volume)**(1.0_wp/3.0_wp)
    
    If (length_from_volume < min_vector_length) Then
      bond_criterion = length_from_volume/2.0_wp
    Else
      bond_criterion = min_vector_length/2.0_wp
    End If
    
    If (model_data%extra_bonds%invoke%fread) Then
      Do k = 1, model_data%extra_bonds%N0
        If (model_data%extra_bonds%bond(k)%value > bond_criterion) Then
          If (frame ==1) Then
            Write (messages(1),'(1x,a,i6)') '***ERROR: problems with the settings of "&possible_extra_bonds" for the intial&
                                           & configuration.' 
          Else 
            Write (messages(1),'(1x,a,i6,a)') '***ERROR: problems with the settings of "&possible_extra_bonds" for frame',&
                                          & frame, ' of the trajectory.'
          End If
          Write (messages(2),'(1x,a,i3,a)') 'Bonding cutoff for the type of bond ', k, &
                                     &' is larger than the size of the simulation cell. Please review the values.'
          Call info(messages, 2)
          Call error_stop(' ')
        End If
      End Do
    End If
    
    If (model_data%chem%acceptor%cutoff%value > bond_criterion) Then
      If (frame ==1) Then
        Write (messages(1),'(1x,a,i6)') '***ERROR: problems with the settings of "&acceptor_criteria" for the intial&
                                       & configuration.' 
      Else 
        Write (messages(1),'(1x,a,i6,a)') '***ERROR: problems with the settings of "&acceptor_criteria" for frame',&
                                      & frame, ' of the trajectory.'
      End If
      Write (messages(2),'(1x,a,i3,a)') 'Cutoff value to define the environment is larger than the size of the&
                                       & simulation cell. Please review the value.'
      Call info(messages, 2)
      Call error_stop(' ')
    End If        

    If (model_data%chem%bonds%cutoff%value > bond_criterion) Then
      If (frame ==1) Then
        Write (messages(1),'(1x,a,i6)') '***ERROR: problems with the settings of "&bonding_criteria" for the intial&
                                       & configuration.' 
      Else 
        Write (messages(1),'(1x,a,i6,a)') '***ERROR: problems with the settings of "&bonding_criteria" for frame',&
                                      & frame, ' of the trajectory.'
      End If
      Write (messages(2),'(1x,a,i3,a)') 'Cutoff value to define the bonds is larger than the size of the&
                                       & simulation cell. Please review the value.'
      Call info(messages, 2)
      Call error_stop(' ')
    End If        

  End Subroutine check_definition_bonds

  Subroutine check_cell_consistency(model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the consistency between the simulation cell and the
    ! atomic coordinates. Adapted from ALC_EQCM  
    !
    ! author    - i.scivetti April 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(model_type), Intent(InOut) :: model_data

    Integer(Kind=wi) :: i, j
    Real(Kind=wp)    :: dist, dist0
    Logical          :: error 

    Real(Kind=wp)       :: a(3), b(3), r(3)
    Character(Len=256)  :: message, messages(2)
    Real(Kind=wp)       :: r0(model_data%config%num_atoms,3)
    Logical             :: changed_geo(model_data%config%num_atoms)

    Do i=1, model_data%config%num_atoms
      r0(i,:)=model_data%config%atom(i)%r(:)
    End Do

    error=.False.

    Write (messages(1),'(a)') '***ERROR: Inconsistency found between the atomic coordinates&
                               & and the simulation cell vectors.' 
                               
    If (model_data%config%cell_units%fread .Or. model_data%config%position_units%fread) Then
      If (model_data%config%cell_units%fread .And. .Not. model_data%config%position_units%fread) Then
        Write (messages(2),'(3x, a)') 'Please also verify input coordinates,&
                                     & as well as the setting for "cell_units".'
      Else If (.Not. model_data%config%cell_units%fread .And. model_data%config%position_units%fread) Then
        Write (messages(2),'(3x, a)') 'Please also verify input coordinates,&
                                     & as well as the setting for "position_units".'
      Else If (model_data%config%cell_units%fread .And. model_data%config%position_units%fread) Then
        Write (messages(2),'(3x, a)') 'Please also verify input coordinates,&
                                     & as well as the settings for "cell_units" and "position_units".'
      End If
    Else
      Write (messages(2),'(3x, a)') 'Please also verify input coordinates.'
    End If
    ! Move all inside the ce
    Do i = 1, model_data%config%num_atoms
      r=model_data%config%atom(i)%r
      Call check_PBC(r, model_data%config%cell, model_data%config%invcell, 1.0_wp, changed_geo(i))
      If (changed_geo(i)) Then
        model_data%config%atom(i)%r=r
      End If
    End Do

    ! Try to move atoms again....if any atom is now moved, there is an inconsistency
    ! between the coordinates and the cell 
    Do i =  1, model_data%config%num_atoms
      If (changed_geo(i)) Then
        r=model_data%config%atom(i)%r
        Call check_PBC(r, model_data%config%cell, model_data%config%invcell, 1.0_wp, changed_geo(i))
        If (changed_geo(i)) Then
          ! Once again. If any atom is now moved, there is an inconsistency
          Call check_PBC(r, model_data%config%cell, model_data%config%invcell, 1.0_wp, changed_geo(i))
        End If  
        If (changed_geo(i)) Then
          Write (message,'(a,i4)') '***PROBLEMS with the position of atom',  i
          Call info(message, 1)
          error=.True.
        End If
      End If
    End Do

    If (error) Then
       Call info(messages, 2)
       Call error_stop(' ')
    End If

    Do i = 1, model_data%config%num_atoms-1
      Do j = i+1, model_data%config%num_atoms
        a=model_data%config%atom(i)%r
        b=model_data%config%atom(j)%r
        Call compute_distance_PBC(a, b, model_data%config%cell, model_data%config%invcell, dist)
        a(:)=r0(i,:)
        b(:)=r0(j,:)
        Call compute_distance_PBC(a, b, model_data%config%cell, model_data%config%invcell, dist0)

        If (Abs(dist-dist0)>length_tol) Then
          Write (message,'(a,2(i7,a))') '***PROBLEMS: Distance between atom ', i, ' and ', j, &
                                   &' does not comply with the crystal symmetry imposed by the cell vectors.'
          Call info(message,1)
          Call info(messages, 2)
          Call error_stop(' ')
        Else  
          If (dist < min_intra) Then
            Write (message,'(a,2(i7,a))') '***PROBLEMS: Distance between atom ', i, ' and ', j, &
                                     &' is too short. It is likely that the cell dimensions are&
                                     & shorter than the system size.'
            Call info(message,1)
            If (model_data%config%simulation_cell%fread) Then
              Write (message,'(13x, a)')  'The user should double check the settings for the &simulation_cell block.'
            Else
              Write (message,'(13x, a)')  'Either the input model is UNPHYSICAL or a HUGE error was&
                                        & introduced in the generation of the trajectory.' 
            End If
            Call info(message,1)
            Call info(messages, 2)
            Call error_stop(' ')
          End If
        End If
      End Do
    End Do
 
  End Subroutine check_cell_consistency
  
  
  Subroutine compute_distance_PBC(a, b, cell, invcell, dist)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the distance between atoms with PCB 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   ) :: a(3), b(3)
    Real(Kind=wp), Intent(In   ) :: cell(3,3)
    Real(Kind=wp), Intent(In   ) :: invcell(3,3)
    Real(Kind=wp), Intent(  Out) :: dist

    Real(Kind=wp) :: Dr_cart(3)
    Logical :: modified
    Integer :: ir

    ! Vector difference
    Do ir=1,3
      Dr_cart(ir)=a(ir)-b(ir) 
    End Do

    ! Find the vector difference for the nearest neighbours (NN)
    Call check_PBC(Dr_cart, cell, invcell, 0.5_wp, modified)
    ! Calculate norm
    dist=norm2(Dr_cart)

  End Subroutine compute_distance_PBC

  Subroutine check_PBC(v_cart, basis, inv_basis, ratio, changed_geo)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the components of given vector "v_cart" 
    ! (given in cartesian coordinates of the Euclidean space) in terms of a general
    ! 3D vector basis, specified by the 3x3 matrix "basis".
    ! The inverse of the basis matrix is named "inv_basis".
    ! The factor "ratio" is used to evaluate how large these components are, 
    ! and modify the "vect" accordingly to accound for PBC.
    !
    ! The input value of "ratio" depends on the quantity to evalue. Thus,
    ! - ratio=0.5 is used for nearest-neighbour distances 
    ! - ratio=1.0 is used to evaluate if atomic positions lie within the volume 
    !   defined by the basis.
    !
    ! Logical variable "changed_geo" is used to check is the vector has been modified 
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(InOut) :: v_cart(3)
    Real(Kind=wp), Intent(In   ) :: basis(3,3)
    Real(Kind=wp), Intent(In   ) :: inv_basis(3,3)
    Real(Kind=wp), Intent(In   ) :: ratio 
    Logical,       Intent(  Out) :: changed_geo                     

    Real(Kind=wp) :: v_direct(3), limit1, limit2
    Integer(Kind=wi) :: ir, i
    Logical          :: flag, flag2

    If (Abs(ratio-0.5_wp) < epsilon(ratio)) Then
      limit1=ratio
      limit2=-ratio
    Else If (Abs(ratio-1.0_wp) < epsilon(ratio)) Then
      limit1=ratio !+length_tol
      limit2=0.0_wp !-length_tol
    End If

    ! Express vector difference in terms of the cell vectors
    v_direct(1:3)= MatMul(v_cart(1:3), inv_basis(1:3,1:3))

    ! PCB effect
    i=1
    flag=.True.
    changed_geo=.False.
    Do While (i< 4 .And. flag)
      flag2=.False.
      Do ir = 1, 3
        If (v_direct(ir) > limit1) Then
           v_cart(:)= v_cart(:) - basis(ir,:)
           changed_geo=.True.
           flag2=.True.
        Else If (v_direct(ir) < limit2) Then
           v_cart(:)= v_cart(:) + basis(ir,:)
           changed_geo=.True.
           flag2=.True.
        End If
      End Do
      If (flag2) Then
        v_direct(1:3)= MatMul(v_cart(1:3), inv_basis(1:3,1:3))
      Else
        flag=.False.
      End If
      i=i+1
    End Do

  End Subroutine check_PBC 
 
  Subroutine check_orthorhombic_cell(A, flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check if the simulation cell 'A' is orthorhombic or not. 
    ! vectors. If A is not orthorhombic, flag will be set to False
    !
    ! author    - i.scivetti Aug 2022
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Real(Kind=wp), Intent(In   )  :: A(3,3)
    Logical,       Intent(  Out)  :: flag
   
    Integer(Kind=wi) :: i, j

    i=1
    flag=.True.
    Do i = 1, 2
     Do j = i+1, 3
       If (Abs(Dot_product(A(i,:), A(j,:)))>epsilon(1.0_wp)) Then
         flag=.False.      
       End If        
     End Do
    End Do 

  End Subroutine check_orthorhombic_cell
 
  Subroutine about_cell(A, invA, length, volume)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to invert a general 3x3 matrix and compute the length of the cell
    ! vectors.
    ! Matrix A is the input matrix.
    ! Matrix invA is the output matrix (inverse of A)
    ! Volume
    !
    ! author    - i.scivetti Jan 2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Real(Kind=wp), Intent(In   )  :: A(3,3)
    Real(Kind=wp), Intent(  Out)  :: invA(3,3)
    Real(Kind=wp), Intent(  Out)  :: length(3)
    Real(Kind=wp), Intent(  Out)  :: volume 
    
    Real(Kind=wp) :: Det
    Real(Kind=wp) :: Cofactor(3,3)
   
    Integer(Kind=wi) :: i, j

    length = 0.0_wp     

    Do i = 1, 3 
      Do j= 1, 3
        length(i) = length(i)+ A(i,j)**2   
      End Do
      length(i)=sqrt(length(i))
    End Do 

    Det =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)

    If (Abs(Det) <= epsilon(det)) Then
      Call error_stop('***ERROR: The determinant of the simulation cell is zero. Please check the definition of the cell vectors')
    End If

    Volume=Det
    
    Cofactor(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    Cofactor(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    Cofactor(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    Cofactor(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    Cofactor(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    Cofactor(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    Cofactor(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    Cofactor(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    Cofactor(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    invA = Transpose(Cofactor) / Det
    
  End Subroutine about_cell

  Subroutine check_model_settings(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives related to the atomistic models
    ! as well as chemistry changes
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(3)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'

    ! Check "input_geometry_format"
    If (model_data%input_geometry_format%fread) Then
      If (model_data%input_geometry_format%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set), ' Wrong specification for directive "input_geometry_format"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        If (Trim(model_data%input_geometry_format%type) /= 'xyz' .And. &
            Trim(model_data%input_geometry_format%type) /= 'vasp') Then
          Write (messages(1),'(2a)') Trim(error_set), ' Valid options for directive "input_geometry_format":'
          Write (messages(2),'(1x,a)') '- xyz'
          Write (messages(2),'(1x,a)') '- vasp'
          Call info(messages, 2)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(1),'(a,1x,a)')  Trim(error_set), 'Specification of the "input_geometry_format" directive&
                                     & is required.'
      Call info(messages,1)
      Call error_stop(' ') 
    End If
    
    If (.Not. model_data%input_composition%invoke%fread) Then
      Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                              & 'Specification of the "&input_composition" block is missing!'
      Call info(messages,1)
      Call error_stop(' ')
    End If

    If (.Not. model_data%config%simulation_cell%fread) Then
      If ((Trim(model_data%input_geometry_format%type) == 'xyz')) Then
        Write (messages(1),'(2(1x,a))')  Trim(error_set), 'Specification of the "&simulation_cell" block&
                                      & is required for trajectories in xyz format'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    End If

    ! Chek directive for cell units
    If (model_data%config%cell_units%fread) Then
      If (model_data%config%cell_units%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set), ' Wrong specification for directive "cell_units"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        If (Trim(model_data%input_geometry_format%type) == 'vasp') Then
           Write (messages(1),'(2a)') Trim(error_set), ' Definition of "cell_units" is incompatible with the&
                                       & option "vasp" for "input_geometry_format". Please review'
           Call info(messages, 1)
           Call error_stop(' ') 
        End If
        Call capital_to_lower_case(model_data%config%cell_units%type)
        If (Trim(model_data%config%cell_units%type) /= 'angstrom' .And.&
            Trim(model_data%config%cell_units%type) /= 'bohr') Then
          Write (messages(1),'(2a)') Trim(error_set), ' Specification for directive "cell_units"&
                                      & should either be "Angstrom" or "Bohr"'
          Call info(messages, 1)
          Call error_stop(' ') 
        Else
          If (Trim(model_data%config%cell_units%type) == 'angstrom') Then
            model_data%config%cell_scaling=1.0_wp
          Else If (Trim(model_data%config%cell_units%type) == 'bohr') Then  
            If (Trim(model_data%input_geometry_format%type) /= 'xyz') Then
              Write (messages(1),'(2a)') Trim(error_set), ' Definition of "cell_units" is incompatible with the&
                                       & choice of "input_geometry_format". Please review'
             Call info(messages, 1)
             Call error_stop(' ') 
            End If
            model_data%config%cell_scaling=Bohr_to_A
          End If
        End If
      End If
    Else
      If (Trim(model_data%input_geometry_format%type) == 'xyz') Then
        Write (messages(1),'(a,1x,a)')  Trim(error_set), 'Specification of the "cell_units" directive&
                                     & is required.'
        Call info(messages,1)
        Call error_stop(' ')
      Else
        model_data%config%cell_scaling=1.0_wp
      End If
    End If

    ! Chek directive for position units
    If (model_data%config%position_units%fread) Then
      If (model_data%config%position_units%fail) Then
        Write (messages(1),'(2a)')  Trim(error_set), ' Wrong specification for directive "position_units"'
        Call info(messages,1)
        Call error_stop(' ') 
      Else
        Call capital_to_lower_case(model_data%config%position_units%type)
        If (Trim(model_data%config%position_units%type) /= 'angstrom' .And.&
            Trim(model_data%config%position_units%type) /= 'bohr') Then
          Write (messages(1),'(2a)') Trim(error_set), ' Specification for directive "position_units"&
                                      & should either be "Angstrom" or "Bohr"'
          Call info(messages, 1)
          Call error_stop(' ') 
        Else
          If (Trim(model_data%config%position_units%type) == 'bohr') Then  
            If (Trim(model_data%input_geometry_format%type) == 'xyz') Then
              Write (messages(1),'(2a)') Trim(error_set), ' Definition of "position_units" is incompatible with the&
                                       & choice of "input_geometry_format". Please review'
             Call info(messages, 1)
             Call error_stop(' ') 
            End If
            model_data%config%position_scaling=Bohr_to_A
          End If
        End If
      End If
    Else
      If (Trim(model_data%input_geometry_format%type) /= 'xyz'   .And. &
          Trim(model_data%input_geometry_format%type) /= 'vasp')  Then
        Write (messages(1),'(a,1x,a)')  Trim(error_set), 'Specification of the "position_units" directive&
                                     & is required.'
        Call info(messages,1)
        Call error_stop(' ')
      Else
        model_data%config%position_scaling=1.0_wp
      End If
    End If

    ! Change chemistry
    If (model_data%change_chemistry%fread) Then
      If (model_data%change_chemistry%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "change_chemistry" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    End If
    
    If (model_data%search_chemistry%fread) Then
      If (model_data%change_chemistry%stat) Then
        Call check_chemistry_criteria(files, model_data)
      Else
        Write (messages(1),'(1x,a)') 'WARNING: the "&search_chemistry" block is defined but it is not relevant& 
                                  & to for the analysis of the trajectory, as the user has set the "change_chemistry"&
                                  & directive to "False"'
        Call info(messages, 1)         
      End If
    End If

    If (model_data%extra_bonds%invoke%fread) Then
      If (model_data%change_chemistry%stat) Then
        Call check_extra_bonding(files, model_data)
      Else
        Write (messages(1),'(1x,a)') 'WARNING: the "&possible_extra_bonds" block is defined but it is not relevant&
                                  & to for the analysis of the trajectory, as the user has set the "change_chemistry"&
                                  & directive to "False"'
        Call info(messages, 1)
      End If
    End If
    
    ! This needs to be re-evaluated under the condition of analysis when defined 
    If (model_data%config%monitored_species%fread) Then
      Call check_definition_monitored_species(files, model_data)
    End If

    ! Check &selected_nn_distances
    If (model_data%nndist%invoke%fread) Then
      Call check_selected_nn_distances(files, model_data%nndist, model_data)
    End If
    
  End Subroutine check_model_settings

  Subroutine check_chemistry_criteria(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives of the &search_chemistry block
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&search_chemistry" block'

    If (model_data%chem%N0%fread) Then
      If (model_data%chem%N0%fail) Then
        Write (messages(2),'(1x,a)')  'Wrong specification for directive "total_number"'
        Call info(messages,2)
        Call error_stop(' ') 
      Else
        If (model_data%chem%N0%value < 1) Then
          Write (messages(2),'(1x,a)')  'Directive "total_number" must be positive!'
          Call info(messages,2)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(2),'(1x,a)')  'The "total_number" directive has not been defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    End If

    If (model_data%chem%bonds%criteria%fread) Then
      Call check_bonding_criteria(files, model_data) 
    Else
      Write (messages(2),'(1x,a)')  'The "&bonding_criteria" sub-block MUST BE defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    End If

    If (model_data%chem%acceptor%criteria%fread) Then
      Call check_acceptor_criteria(files, model_data) 
    Else
      Write (messages(2),'(1x,a)')  'The "&enviroment_criteria" sub-block MUST BE defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    End If
    
  End Subroutine check_chemistry_criteria

  Subroutine check_acceptor_criteria(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives of the &enviroment_criteria 
    ! sub-block
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(3)
    Character(Len=64 )  :: error_set
    Integer             :: j, k
    Logical             :: found
    
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&enviroment_criteria" sub-block.'

    ! environment cutoff
    If (.Not. model_data%chem%acceptor%cutoff%fread) Then
      model_data%chem%acceptor%cutoff%tag='cutoff'
    End If
    Call check_length_directive(model_data%chem%acceptor%cutoff, messages(1), .True., 'directive')

    ! Exclude pairs
    If (model_data%chem%acceptor%info_exclude%fread) Then
      If (model_data%chem%acceptor%info_exclude%fail) Then    
        Write (messages(2),'(1x,a)')  'Wrong specification for directive "exclude_pairs". See the "use_code.md" file.'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If

      If (model_data%chem%acceptor%N0_excl < 1) Then
        Write (messages(2),'(1x,a)')  'The number of "exclude_pairs" must be positive'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
      
      Do k=1, model_data%chem%acceptor%N0_excl
        found=.False.
        Do j=1, model_data%input_composition%atomic_species
          If (Trim(model_data%chem%acceptor%tg_excl(k))== Trim(model_data%input_composition%tag(j))) Then
            found=.True.
          End If
        End Do
        If (.Not. found) Then
          Write (messages(2),'(3(1x,a))') 'The tag "', Trim(model_data%chem%acceptor%tg_excl(k)), &
                                      & '" in the "exclude_pairs" directive is not defined as a tag in the&
                                      & "&input_composition" block. Please review'
          Call info(messages, 2)
          Call error_stop(' ') 
        End If
      End Do

      Do j=1, model_data%chem%acceptor%N0_excl-1
        Do k=j+1, model_data%chem%acceptor%N0_excl
          If (Trim(model_data%chem%acceptor%tg_excl(j))==Trim(model_data%chem%acceptor%tg_excl(k))) Then
            Write (messages(2),'(3(1x,a))') 'Tag', Trim(model_data%chem%acceptor%tg_excl(j)), 'is repeated in the list!'
            Write (messages(3),'((1x,a))') 'The tags defined in "exclude_pairs" must be  different'
            Call info(messages, 3)
            Call error_stop(' ')
          End If
        End Do
      End Do 

    End If

    If (model_data%chem%acceptor%info_include%fread) Then
      If (model_data%chem%acceptor%info_include%fail) Then    
        Write (messages(2),'(1x,a)')  'Wrong specification for directive "include_tags". See the "use_code.md" file.'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If

      Do k=1, model_data%chem%acceptor%N0_incl
        found=.False.
        Do j=1, model_data%input_composition%atomic_species
          If (Trim(model_data%chem%acceptor%tg_incl(k))== Trim(model_data%input_composition%tag(j))) Then
            found=.True.
          End If
        End Do
        If (.Not. found) Then
          Write (messages(2),'(3(1x,a))') 'The tag "', Trim(model_data%chem%acceptor%tg_incl(k)), &
                                      & '" in the "include_tags" directive is not defined as a tag in the&
                                      & "&input_composition" block. Please review'
          Call info(messages, 2)
          Call error_stop(' ') 
        End If
      End Do

      Do j=1, model_data%chem%acceptor%N0_incl-1
       Do k=j+1, model_data%chem%acceptor%N0_incl
         If (Trim(model_data%chem%acceptor%tg_incl(j))==Trim(model_data%chem%acceptor%tg_incl(k))) Then
           Write (messages(2),'(3(1x,a))') 'Tag', Trim(model_data%chem%acceptor%tg_incl(j)), 'is repeated in the list!'
           Write (messages(3),'((1x,a))') 'The tags defined in "include_tags" must be  different'
           Call info(messages, 3)
           Call error_stop(' ')
         End If
       End Do
      End Do  
      
    Else
      Write (messages(2),'(1x,a)')  'The "include_tags" directive has not been defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    End If
    
    If (model_data%chem%acceptor%info_exclude%fread .And. model_data%chem%acceptor%info_include%fread) Then
      Do k=1, model_data%chem%acceptor%N0_excl
        found=.False.
        Do j=1, model_data%chem%acceptor%N0_incl
          If (Trim(model_data%chem%acceptor%tg_excl(k))== Trim(model_data%chem%acceptor%tg_incl(j))) Then
            found=.True.
          End If
        End Do
        If (.Not. found) Then
          Write (messages(2),'(3(1x,a))') 'The tag "', Trim(model_data%chem%acceptor%tg_excl(k)), &
                                      & '" of the "exclude_pairs" directive is not defined as a tag in the&
                                      & "include_tags" list. Please review.'
          Call info(messages, 2)
          Call error_stop(' ') 
        End If
      End Do
    End If    
    
    If (model_data%chem%acceptor%info_exclude%fread .And. model_data%chem%acceptor%info_exclude%fread) Then
      If (model_data%chem%acceptor%N0_excl==model_data%chem%acceptor%N0_incl .And. model_data%chem%acceptor%N0_excl==1) Then
          Write (messages(2),'(a)') 'The definition of "include_tags" and "exclude_pairs" are incosistent: only one&
                                    & tag is defined and the same tag is requested to be excluded. Please review.'
          Call info(messages, 2)
          Call error_stop(' ') 
      End If
    End If
   
    ! Intialise initial checking of environment
    model_data%chem%acceptor%check=.True.

  End Subroutine check_acceptor_criteria
  
  Subroutine check_bonding_criteria(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives of the &bonding_criteria sub-block
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer(Kind=wi)    :: j
    Logical             :: found
    
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&bonding_criteria" sub-block.'

    If (model_data%chem%bonds%species%fread) Then
      If (model_data%chem%bonds%species%fail) Then
        Write (messages(2),'(1x,a)')  'Wrong specification for directive "only_element"'
        Call info(messages,2)
        Call error_stop(' ') 
      End If
      found=.False.
      Do j=1, model_data%input_composition%atomic_species
        If (Trim(model_data%chem%bonds%species%type)== Trim(model_data%input_composition%element(j))) Then
          found=.True.
        End If
      End Do
      If (.Not. found) Then
        Write (messages(2),'(3(1x,a))') 'The definition of the "only_element" directive is "',&
                                      & Trim(model_data%chem%bonds%species%type), &
                                      & '", which is not defined as a valid CHEMICAL ELEMENT&
                                      & in the "&input_composition" block. Please review'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    Else
      Write (messages(2),'(1x,a)')  'The "only_element" directive has not been defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    End If

    ! Number of bonds 
    If (model_data%chem%bonds%N0%fread) Then
      If (model_data%chem%bonds%N0%fail) Then
        Write (messages(2),'(1x,a)')  'Wrong specification for directive "number_of_bonds"'
        Call info(messages,2)
        Call error_stop(' ') 
      Else
        If (model_data%chem%bonds%N0%value < 1) Then
          Write (messages(2),'(1x,a)')  'Directive "number_of_bonds" must be positive!'
          Call info(messages,2)
          Call error_stop(' ') 
        End If
      End If
    Else
      Write (messages(2),'(1x,a)')  'The "number_of_bonds" directive has not been defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    End If
    
    ! Bond cutoff
    If (.Not. model_data%chem%bonds%cutoff%fread) Then
      model_data%chem%bonds%cutoff%tag='cutoff'
    End If
    Call check_length_directive(model_data%chem%bonds%cutoff, messages(1), .True., 'directive')
    
  End Subroutine check_bonding_criteria
  
  Subroutine check_extra_bonding(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives the &possible_extra_bonds block
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256) :: messages(3)
    Character(Len=256) :: error_set
    Character(Len= 4)  :: tg1_element, tg2_element
    Integer            :: j, k
    Logical            :: a1, a2, b1, b2, found
    
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' - ' 
    Write (messages(1),'(a)')  Trim(error_set)//'Wrong specification in the &possible_extra_bonds" block.'
    
    Do k = 1, model_data%extra_bonds%N0
      ! Loop over tg1
      found = .False.
      Do j = 1, model_data%input_composition%atomic_species
        If (Trim(model_data%extra_bonds%tg1(k)) == Trim(model_data%input_composition%tag(j))) Then
          tg1_element=Trim(model_data%input_composition%element(j))
          found = .True.   
        End If
      End Do
      
      If (.Not. found) Then      
        Write (messages(2),'(1x,3a,i3,a)') 'Species "', Trim(model_data%extra_bonds%tg1(k)), '" defined for the first tag& 
                                          & of the type of bond ', k, ' has NOT been defined in the &input_composition&
                                          & block. Please review the settings'
        Call info(messages, 2)         
        Call error_stop(' ')
      End If
      
      ! Loop over tg2
      found = .False.
      Do j = 1, model_data%input_composition%atomic_species
        If (Trim(model_data%extra_bonds%tg2(k)) == Trim(model_data%input_composition%tag(j))) Then
          tg2_element=Trim(model_data%input_composition%element(j))
          found = .True.   
        End If
      End Do
      
      If (.Not. found) Then      
        Write (messages(2),'(1x,3a,i3,a)') 'Species "', Trim(model_data%extra_bonds%tg2(k)), '" defined for the second tag& 
                                          & of the type of bond ', k, ' has NOT been defined in the &input_composition&
                                          & block. Please review the settings'
        Call info(messages, 2)         
        Call error_stop(' ')
      End If

    ! Check intrinsic tags wihtin the bond
      If(Trim(model_data%extra_bonds%tg1(k)) == Trim(model_data%extra_bonds%tg2(k))) Then
        Write (messages(2),'(1x,a,i3,a)') 'The tags for the specification of the type of bond ', k, &
                                          & ' must be different. Please review the settings'
        Call info(messages, 2)         
        Call error_stop(' ')
      End If
      
      ! Check consistent with the "&bonding_criteria" block
      If (Trim(tg1_element) /= Trim(model_data%chem%bonds%species%type) .And. &
          Trim(tg2_element) /= Trim(model_data%chem%bonds%species%type)) Then
          Write (messages(2),'(1x,a,i3,a)') 'None of the chemical element associated with the defined&
                                          & species of the type of bond ', k, &
                                          &' correspond to the setting for "only_element" of the&
                                          & "&bonding_criteria" block. Please correct'
        Call info(messages, 2)         
        Call error_stop(' ')
      End If
      If (Trim(tg1_element) == Trim(model_data%chem%bonds%species%type) .And. &
          Trim(tg2_element) == Trim(model_data%chem%bonds%species%type)) Then
          Write (messages(2),'(1x,a,i3,a)') 'The two chemical element associated with the defined&
                                          & species of the type of bond ', k, &
                                          &' correspond to the setting for "only_element" of the&
                                          & "&bonding_criteria" block. Such a bond is not allowed. Please correct.'
        Call info(messages, 2)         
        Call error_stop(' ')
      End If

      found=.False. 
      Do j = 1, model_data%chem%acceptor%N0_incl 
        If ((Trim(model_data%chem%acceptor%tg_incl(j)) == Trim(model_data%extra_bonds%tg1(k))) .Or. &
            (Trim(model_data%chem%acceptor%tg_incl(j)) == Trim(model_data%extra_bonds%tg2(k)))) Then
            found=.True.
        End If
      End Do

      If (.Not. found) Then      
          Write (messages(2),'(1x,a,i3,5a)') 'None of the two species tags of the type of bond ', k, &
                                          &' (', Trim(model_data%extra_bonds%tg1(k))//' and '//Trim(model_data%extra_bonds%tg2(k)),&
                                          & ') is defined in the "include_tags" directive of the "&acceptor_criteria"&
                                          & block. One of the tags must be part of "include_tags". Please correct.'
        Call info(messages, 2)         
        Call error_stop(' ')
      End If

    End Do

    ! Check for repetition
    Do k = 1, model_data%extra_bonds%N0-1
      Do j = k+1, model_data%extra_bonds%N0
        a1 = Trim(model_data%extra_bonds%tg1(k)) == Trim(model_data%extra_bonds%tg1(j))
        a2 = Trim(model_data%extra_bonds%tg2(k)) == Trim(model_data%extra_bonds%tg2(j))
        b1 = Trim(model_data%extra_bonds%tg1(k)) == Trim(model_data%extra_bonds%tg2(j))
        b2 = Trim(model_data%extra_bonds%tg2(k)) == Trim(model_data%extra_bonds%tg1(j))
        If ((a1 .And. a2) .Or. (b1 .And. b2)) Then 
          Write (messages(2),'(1x,2(a,i2),a)') 'Information for the type of bonds ', k, &
                                   & ' and ', j, ' is redundant. Please check'
          Call info(messages, 2)
          Call error_stop(' ')
        End If    
      End Do
    End Do

  End Subroutine check_extra_bonding

  Subroutine check_definition_monitored_species(files, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the format and directives of the &monitored block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set, ref_element

    Integer(Kind=wi)    :: j, k
    Logical             :: flag, found
    
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&monitored_species" block.'

    If (model_data%species_definition%atomic_components%fread) Then
      ! Check information within the &atomic_components sub-block
      Do k = 1, model_data%species_definition%num_components
        ! Loop over elements
        found = .False.
        Do j = 1, model_data%input_composition%atomic_species
          If (Trim(model_data%species_definition%element(k)) == Trim(model_data%input_composition%element(j))) Then
            found = .True.   
          End If
        End Do
        
        If (.Not. found) Then      
          Write (messages(2),'(1x,3a,i3,a)') '(Sub-block &atomic_components): Atomic element "', &
                                       &Trim(model_data%species_definition%element(k)), &
                                       &'" defined for the ',& 
                                       & k, ' component has NOT been defined in the &input_composition&
                                       & block. Please review the settings'
          Call info(messages, 2)         
          Call error_stop(' ')
        End If
  
        If (model_data%species_definition%N0_element(k)<1) Then
          Write (messages(2),'(1x,3a,i3,a)') '(Sub-block &atomic_components): Atomic element "',&  
                                      & Trim(model_data%species_definition%element(k)),&
                                      & '" defined for the ',&
                                      & k, ' component cannot be lower than 1.'
          Call info(messages, 2)         
          Call error_stop(' ')
        
        ElseIf (model_data%species_definition%num_components>max_components) Then
          Write (messages(2),'(1x,3a,2(i3,a))') '(Sub-block &atomic_components): Atomic element "',& 
                                     & Trim(model_data%species_definition%element(k)),&
                                     &'" defined for the ',& 
                                     & k, ' has more than ', max_components,&
                                     & ' components". Are you sure? Please check' 
          Call info(messages, 2)         
          Call error_stop(' ')
        End If
      End Do
  
     ! Check there is no repetition of elements in the &atomic_components
     Do k = 1, model_data%species_definition%num_components
       Do j=1, model_data%species_definition%num_components
         If (k /= j) Then
           If (Trim(model_data%species_definition%element(k))==Trim(model_data%species_definition%element(j))) Then
             Write (messages(2),'(1x,a)') 'Problems in "&atomic_components": definition for the atomic element "'& 
                                         &//Trim(model_data%species_definition%element(k))//'" is repeated.'
             Call info(messages, 2)
             Call error_stop(' ') 
           End If 
         End If
       End Do
     End Do
  
     ! Calculate total number of atoms per species
      model_data%species_definition%atoms_per_species=0
      Do j=1, model_data%species_definition%num_components
        model_data%species_definition%atoms_per_species = model_data%species_definition%atoms_per_species + &
                                                          model_data%species_definition%N0_element(j)
      End Do
      
      If (model_data%species_definition%atoms_per_species>max_at_species) Then
          Write (messages(2),'(1x,a,i2,a)')  'The &atomic_components block indicates that each species has ',&
                                   & model_data%species_definition%atoms_per_species, ' atoms. Are you sure&
                                   & you want to monitor species with such a large number of atoms?'
         Call info(messages, 2)
         Call error_stop(' ') 
      End If
  
      If (model_data%species_definition%atoms_per_species == 1) Then
          Write (messages(2),'(1x,a,i2,a)')  'The &atomic_components block indicates that each species has ',&
                                   & model_data%species_definition%atoms_per_species, ' atoms. This block is&
                                   & only valid for molecular species.'
         Call info(messages, 2)
         Call error_stop(' ') 
      End If
      
      ! Check that the element of the reference_tag has only one component in the definition of species 
      Do j=1, model_data%species_definition%num_components
        If (Trim(ref_element)==Trim(model_data%species_definition%element(j)))  Then
           If (model_data%species_definition%N0_element(j) /= 1) Then
             Write (messages(2),'(1x,a)')  'The chemical element "'//Trim(ref_element)//'" of "reference_tag",&
                                      & must be have only one component as part of the species composition'
            Call info(messages, 2)
            Call error_stop(' ') 
           End If
        End If
      End Do
    Else
      model_data%species_definition%atoms_per_species=1
    End If
   
    If (.Not. model_data%species_definition%name%fread) Then
      If (model_data%species_definition%atoms_per_species /= 1 ) Then
        Write (messages(2),'(1x,a)')  'The "name" directive has not been defined.'
        Call info(messages, 2)
        Call error_stop(' ') 
      Else
        model_data%species_definition%name%type=Trim(model_data%species_definition%reference_tag%type)
      End If
    Else
      If (model_data%species_definition%atoms_per_species == 1 ) Then
        Write (messages(2),'(1x,a)')  'The "name" directive must only be defined when the number of atoms&
                                    & for each species is larger than 2'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If       
    End If

    If (.Not. model_data%species_definition%reference_tag%fread) Then
      Write (messages(2),'(1x,a)')  'The "reference_tag" directive has not been defined.'
      Call info(messages, 2)
      Call error_stop(' ') 
    Else
      flag=.True.
      j=1
      Do While (j <= model_data%input_composition%atomic_species .And. flag)
        If (Trim(model_data%input_composition%tag(j))==Trim(model_data%species_definition%reference_tag%type)) Then
          flag=.False.
          ref_element=Trim(model_data%input_composition%element(j))
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(2),'(1x,a)')  'The setting for "reference_tag" does not match any of&
                                    & tags defined in the "&input_composition" block'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If 
    End If

    ! Check cutoff criterion
    If (.Not. model_data%species_definition%bond_cutoff%fread) Then
      model_data%species_definition%bond_cutoff%tag='bond_cutoff'
    End If
    Call check_length_directive(model_data%species_definition%bond_cutoff, messages(1), .True., 'directive')
    
    If (model_data%species_definition%compute_amount%fread) Then    
      If (model_data%species_definition%compute_amount%fail) Then   
        Write (messages(2),'(1x,a)')  'Wrong or missing definition for "compute_amount" directive.'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    Else
      model_data%species_definition%compute_amount%stat=.False.
    End If

    If (.Not. model_data%species_definition%name%fread) Then
      If (model_data%species_definition%atoms_per_species /= 1 ) Then
        Write (messages(2),'(1x,a)')  'The "name" directive has not been defined.'
        Call info(messages, 2)
        Call error_stop(' ') 
      Else
        model_data%species_definition%name%type=Trim(model_data%species_definition%reference_tag%type)
      End If
    Else
      If (model_data%species_definition%atoms_per_species == 1 ) Then
        Write (messages(2),'(1x,a)')  'The "name" directive must only be defined when the number of atoms&
                                    & for each species is larger than 2'
        Call info(messages, 2)
        Call error_stop(' ') 
      End If       
    End If

    ! Check intramol_stat_settings 
    If (model_data%species_definition%intra_geom%invoke%fread) Then
      If (model_data%species_definition%intra_geom%angle%invoke%fread) Then
        If (model_data%species_definition%atoms_per_species < 3) Then
        Write (messages(2),'(1x,a)')  'Problems in "'//Trim(model_data%species_definition%intra_geom%invoke%type)//&
                                     &'": it is not possible to define an internal angle when monitored species are diatomic.& 
                                     & Remove "'//Trim(model_data%species_definition%intra_geom%angle%invoke%type)//'"'
        Call info(messages, 2)
        Call error_stop(' ') 
        End If
      End If
      If (model_data%species_definition%intra_geom%dist%invoke%fread) Then
        Call check_settings_geom_param(messages(1), model_data%species_definition%intra_geom%invoke%type,&
                                    & model_data%species_definition%intra_geom%dist)
        Call check_intramol_stat_species(messages(1),model_data%species_definition,model_data%species_definition%intra_geom%dist)
      End If
      If (model_data%species_definition%intra_geom%angle%invoke%fread) Then
        Call check_settings_geom_param(messages(1), model_data%species_definition%intra_geom%invoke%type,&
                                   & model_data%species_definition%intra_geom%angle)
        Call check_intramol_stat_species(messages(1),model_data%species_definition,model_data%species_definition%intra_geom%angle)  
      End If
      If ((.Not. model_data%species_definition%intra_geom%dist%invoke%fread) .And. &
          (.Not. model_data%species_definition%intra_geom%angle%invoke%fread)) Then
           Write (messages(2),'(1x,a)')  'Empty "'//Trim(model_data%species_definition%intra_geom%invoke%type)//&
                                     &'" block! Please define "&distance_parameters" and/or "&angle_parameters",&
                                     & or remove the block.'
           Call info(messages, 2)
           Call error_stop(' ') 
      End If    
    End If

    ! Check intermol_stat_settings 
    If (model_data%species_definition%inter_geom%invoke%fread) Then
      ! Check if only reference tags will be condired as NNs
      If (model_data%species_definition%inter_geom%only_ref_tags_as_nn%fread) Then
        If (model_data%species_definition%inter_geom%only_ref_tags_as_nn%fail) Then
          Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                    & "only_ref_tags_as_nn" (choose either .True. or .False.)'
          Call info(messages,1)
          Call error_stop(' ')
        End If
      Else
        model_data%species_definition%inter_geom%only_ref_tags_as_nn%stat=.False.
      End If
    
      If (model_data%species_definition%inter_geom%dist%invoke%fread) Then
        Call check_settings_geom_param(messages(1), model_data%species_definition%inter_geom%invoke%type,&
                                    & model_data%species_definition%inter_geom%dist)
      End If
      If (model_data%species_definition%inter_geom%angle%invoke%fread) Then
        Call check_settings_geom_param(messages(1), model_data%species_definition%inter_geom%invoke%type,&
                                   & model_data%species_definition%inter_geom%angle)
      End If
      If ((.Not. model_data%species_definition%inter_geom%dist%invoke%fread) .And. &
          (.Not. model_data%species_definition%inter_geom%angle%invoke%fread)) Then
           Write (messages(2),'(1x,a)')  'Empty "'//Trim(model_data%species_definition%inter_geom%invoke%type)//&
                                     &'" block! Please define "&distance_parameters" and/or "&angle_parameters",&
                                     & or remove the block.'
           Call info(messages, 2)
           Call error_stop(' ') 
      End If    
    End If
    
    
  End Subroutine check_definition_monitored_species

  Subroutine check_intramol_stat_species(error_set, T, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check species defined to 
    ! compute the statistices of intramolecular
    ! geometry
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256),   Intent(In   ) :: error_set
    Type(spec_def_type),  Intent(In   ) :: T
    Type(geo_param_type), Intent(InOut) :: M   
  
    Integer(Kind=wi)   ::  k, j, n
    Character(Len=1)   :: num
    Character(Len=256) :: messages(3)
    

    messages(1)=error_set
    Write (messages(2),'(1x,a)') 'Problems in "'//Trim(M%invoke%type)//'" of "&intramol_stat_settings"'

    Do j = 1, M%nspecies
      M%num_spec(j)=0
      Do k = 1, T%num_components
        If (Trim(M%species(j))==Trim(T%element(k)))Then
          M%num_spec(j)=M%num_spec(j)+1
        End If
      End Do      
      If (M%num_spec(j)==0) Then
        Write(num,'(i1)') j 
        Write (messages(3),'(1x,a)') 'Argument '//Trim(num)//' of the "species" directive does not&
                                     & correspond to the elements defined in "&atomic_components"' 
        Call info(messages, 3)
        Call error_stop(' ') 
      End If
    End Do
    
    
    Do k = 1, T%num_components
      n=0
      Do j = 1, M%nspecies
        If (Trim(M%species(j))==Trim(T%element(k)))Then
          n=n+1
        End If
      End Do
      If (n>T%N0_element(k))Then
        Write(num,'(i1)') n 
        Write (messages(3),'(1x,a)') 'The number of times the element "'//Trim(T%element(k))//'" is listed in the&
                                   & "species" directive ('//Trim(num)//' times) exceeds the value set in "&atomic_components"' 
        Call info(messages, 3)
        Call error_stop(' ') 
      End If
    End Do
    
  End Subroutine check_intramol_stat_species
  
  Subroutine check_settings_geom_param(error_set, inblock, M)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check distance and angle settings
    ! for statistics of the monitored species
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Character(Len=256),   Intent(In   ) :: error_set
    Character(Len=256),   Intent(In   ) :: inblock 
    Type(geo_param_type), Intent(InOut) :: M   

    Character(Len=256)  :: messages(3), error
    
    messages(1)=error_set
    
    If (Trim(inblock)=='&intramol_stat_settings') Then
      If (M%tag_species%fread) Then
         If (M%tag_species%fail) Then
            Write (messages(2),'(1x,a)') 'Check "'//Trim(inblock)//'" block: Problems to read the "'//Trim(M%tag_species%type)//&
                                        &'" directive within "'//Trim(M%invoke%type)//'".' 
            Call info(messages, 2)
            Call error_stop(' ') 
         End If
      Else
        Write (messages(2),'(1x,a)')  'Problems in "'//Trim(inblock)//'": The user must define the species involved&
                                    & inside "'//Trim(M%invoke%type)//'" using the "species" directive, which is missing.' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    Else If (Trim(inblock)=='&intermol_stat_settings') Then 
      If (M%tag_species%fread) Then
         Write (messages(2),'(1x,a)') 'Check "'//Trim(inblock)//'" block: the definition of the "'//Trim(M%tag_species%type)//&
                                     &'" directive within "'//Trim(M%invoke%type)//'" is not necessary.' 
         Write (messages(3),'(1x,a)') 'The statistical analysis is carried out using the "reference_tag" defined in&
                                     & "&monitored_species". Please remove "'//Trim(M%tag_species%type)//'" from this block.' 

         Call info(messages, 3)
         Call error_stop(' ') 
      End If
      If (M%tag_species%fread) Then
         Write (messages(2),'(1x,a)') 'Check "'//Trim(inblock)//'" block: the definition of the "'//Trim(M%tag_species%type)//&
                                     &'" directive within "'//Trim(M%invoke%type)//'" is not required.' 
         Write (messages(3),'(1x,a)') 'Please remove "'//Trim(M%tag_species%type)//'" from this block.' 

         Call info(messages, 3)
         Call error_stop(' ') 
      End If
    End If 
 
    ! Error message just in case....
    error=Trim(messages(1))//' Check "'//Trim(M%invoke%type)//'" inside "'//Trim(inblock)//'".'
    
    !Check lower_bound, upper_bound and delta
    If (.Not. M%lower_bound%fread) Then
      M%lower_bound%tag='lower_bound'
    End If
    If (.Not. M%upper_bound%fread) Then
      M%upper_bound%tag='upper_bound'
    End If
    If (.Not. M%delta%fread) Then
      M%delta%tag='delta'
    End If
    
    If (Trim(M%invoke%type) == '&distance_parameters') Then
      Call check_length_directive(M%lower_bound, error, .True., 'directive')
      Call check_length_directive(M%upper_bound, error, .True., 'directive')
      Call check_length_directive(M%delta, error, .True., 'directive')
      If (M%lower_bound%value >= M%upper_bound%value) Then
        Write (messages(2),'(1x,a)')  'Problems with "'//Trim(M%invoke%type)//'" in "'//Trim(inblock)//'": The value of&
                                    & "upper_bound" must be larger than "lower_bound" (make sure this is the case if&
                                    & you use different units)' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    Else If (Trim(M%invoke%type) == '&angle_parameters') Then
      Call check_angle_directive(M%lower_bound, error, .True., 'directive')
      Call check_angle_directive(M%upper_bound, error, .True., 'directive')
      Call check_angle_directive(M%delta, error, .True., 'directive')
      If (M%lower_bound%value >= M%upper_bound%value) Then
        Write (messages(2),'(1x,a)')  'Problems with "'//Trim(M%invoke%type)//'" in "'//Trim(inblock)//'": The value of&
                                    & "upper_bound" must be larger than "lower_bound" (make sure this is the case if&
                                    & you use different units)' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If
    End If  
     
  End Subroutine check_settings_geom_param
    
  Subroutine check_selected_nn_distances(files, M, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the definition of the
    ! parameters defined in the &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),       Intent(InOut) :: files(:)
    Type(short_dist_type), Intent(InOut) :: M 
    Type(model_type),      Intent(In   ) :: model_data

    Character(Len=256)  :: messages(3), error_set
    Character(Len=8)    :: tagj, tagk
    Logical             :: flag
    Integer(Kind=wi)    :: j, k
    
    ! Error message just in case....
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&selected_nn_distances" block.'

    If(M%tag_nn_species%fread) Then
      If (M%tag_nn_species%fail) Then
        Write (messages(2),'(1x,a)')  'Problems to define the "nn_species" directive'  
        call info(messages, 2)
        call error_stop(' ')
      End If
    Else
      Write (messages(2),'(1x,a)')  'The user must define the "nn_species" directive'  
      call info(messages, 2)
      call error_stop(' ')
    End If

    If(M%tag_reference_species%fread) Then
      If (M%tag_reference_species%fail) Then
        Write (messages(2),'(1x,a)')  'Problems to define the "reference_species" directive'  
        call info(messages, 2)
        call error_stop(' ')
      End If
    Else
      Write (messages(2),'(1x,a)')  'The user must define the "reference_species" directive'  
      call info(messages, 2)
      call error_stop(' ')
    End If

    !Check if the reference_species is defined in the &input_composition block  
    tagk=Trim(M%reference_species)
    Call remove_symbols(tagk,'*')
    flag=.True.
    j=1
    Do While (j <= model_data%input_composition%atomic_species .And. flag)
      If (Trim(model_data%input_composition%tag(j))==Trim(tagk)) Then
        flag=.False.
      End If  
      j=j+1
    End Do
    If (flag) Then
      Write (messages(2),'(1x,a)')   'The tag "'//Trim(M%reference_species)//&
                                     &'" defined in the "reference_species" directive is not a valid option.&
                                     & Please check the definition of the &input_composition block' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If 
    
    !Check if tags in nn_species are defined in the &input_composition block  
    Do k=1, M%num_nn_species
      tagk=Trim(M%nn_species(k))
      Call remove_symbols(tagk,'*')
      flag=.True.
      j=1
      Do While (j <= model_data%input_composition%atomic_species .And. flag)
        If (Trim(model_data%input_composition%tag(j))==Trim(tagk)) Then
          flag=.False.
        End If  
        j=j+1
      End Do
      If (flag) Then
        Write (messages(2),'(1x,a)')   'The tag "'//Trim(M%nn_species(k))//&
                                       &'" defined in the "nn_species" directive is not a valid option.&
                                       & Please check the definition of the &input_composition block' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If 
    End Do

    !Check if tags defined in nn_species are repeated
    Do j=1, M%num_nn_species-1
      tagj=Trim(M%nn_species(j))
      Do k=j+1, M%num_nn_species 
        tagk=Trim(M%nn_species(k))
        If (Trim(tagj)==Trim(tagk)) Then
          Write (messages(2),'(1x,a)')   'The tag "'//Trim(tagj)//&
                                         &'" is repeated in the specification of the "nn_species" directive.&
                                         & Please remove this duplication.' 
          Call info(messages, 2)
          Call error_stop(' ') 
        End If 
      End Do
    End Do

    !Check lower_bound, upper_bound and delta
    If (.Not. M%lower_bound%fread) Then
      M%lower_bound%tag='lower_bound'
    End If
    If (.Not. M%upper_bound%fread) Then
      M%upper_bound%tag='upper_bound'
    End If
    If (.Not. M%dr%fread) Then
      M%dr%tag='delta'
    End If
    
    Call check_length_directive(M%lower_bound, messages(1), .True., 'directive')
    Call check_length_directive(M%upper_bound, messages(1), .True., 'directive')
    Call check_length_directive(M%dr,          messages(1), .True., 'directive')
    If (M%lower_bound%value >= M%upper_bound%value) Then
      Write (messages(2),'(1x,a)')  'The value of "upper_bound" must be larger than "lower_bound"&
                                  & (make sure this is the case if you use different units)' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If
     
  End Subroutine check_selected_nn_distances
  
  
  Subroutine check_angle_directive(T, error_set, kill, type_directive)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check angle related directives
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: error_set
    Logical,                 Intent(In   )  :: kill
    Character(Len=*),        Intent(In   )  :: type_directive
  
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
        If (Trim(T%units) /= 'rads' .And. Trim(T%units) /= 'degrees') Then
           If (Trim(type_directive) /= 'inblock') Then
             Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                                      & 'Units for directive "'//Trim(T%tag)//'" must be "Degrees" or "Rads".&
                                      & Have you defined the units? Please review.'
           Else
             Write (messages(1),'(2(1x,a))')  Trim(error_set), '. Units for angles must be "Degrees" or "Rads".&
                                           & Have you defined the units? Please review'
           End If
           Call info(messages, 1)
           Call error_stop(' ')
        End If
        ! Transform to Angstrom
        If (Trim(T%units) == 'rads') Then
           T%value=Rads_to_degrees * T%value
        End If
      End If
    Else 
      If (kill) Then
        If (Trim(type_directive) /= 'inblock') Then
          Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(T%tag)//'" directive'
        Else
          Write (messages(1),'(1x,a)')  Trim(error_set)
        End If  
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_angle_directive  
  
  
  Subroutine check_length_directive(T, error_set, kill, type_directive)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check length related directives
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(in_param),          Intent(InOut)  :: T
    Character(Len=*),        Intent(In   )  :: error_set
    Logical,                 Intent(In   )  :: kill
    Character(Len=*),        Intent(In   )  :: type_directive
  
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
        If (Trim(T%units) /= 'angstrom' .And. Trim(T%units) /= 'bohr') Then
           If (Trim(type_directive) /= 'inblock') Then
             Write (messages(1),'(2(1x,a))')  Trim(error_set),&
                                      & 'Units for directive "'//Trim(T%tag)//'" must be "Angstrom" or "Bohr".&
                                      & Have you defined the units? Please review.'
           Else
             Write (messages(1),'(2(1x,a))')  Trim(error_set), '. Units for bonds must be "Angstrom" or "Bohr".&
                                           & Have you defined the units? Please review'
           End If
           Call info(messages, 1)
           Call error_stop(' ')
        End If
        ! Transform to Angstrom
        If (Trim(T%units) == 'bohr') Then
           T%value=Bohr_to_A * T%value
        End If
      End If
    Else 
      If (kill) Then
        If (Trim(type_directive) /= 'inblock') Then
          Write (messages(1),'(2(1x,a))')  Trim(error_set), 'The user must define the "'//Trim(T%tag)//'" directive'
        Else
          Write (messages(1),'(1x,a)')  Trim(error_set)
        End If  
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    End If
    
  End Subroutine check_length_directive  
  
  Subroutine print_model_settings(files, model_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print a summary of the atomic model settings
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(InOut) :: model_data

    Character(Len=256) :: messages(4), word, N0chem, nbonds, ncutoff
    Integer(Kind=wi)   :: j

    Call info(' ', 1) 
    Call info('Model settings', 1) 
    Call info('==============', 1) 
    Write (word,*) model_data%input_composition%numtot
    Write (messages(1),'(1x,a)') '- the atomic model contains '//Trim(Adjustl(word))//' atoms.&    
                                 & The initial tag and element for each atomic species is defined&
                                 & in the "&input_composition" block.'

                              
    If ((Trim(model_data%input_geometry_format%type) == 'xyz')) Then
      Write (messages(2),'(1x,a)') '- the size of the model is defined by the lattice vectors&
                                     & as specified in the "&simulation_cell" block.'
    Else
      Write (messages(2),'(1x,a)') '- the size of the model is specified in the '//&
                                    &Trim(files(FILE_TRAJECTORY)%filename)//' file.'
      
    End If
    Call info(messages, 2)

    If (model_data%config%monitored_species%fread) Then
      Write (messages(1),'(1x,a)') '- the "&monitored_species" block specifies the composition of the&
                                 & non-reactive part of the model that will be analysed.'
      Call info(messages, 1)
    End If

    If (model_data%change_chemistry%stat) Then
      Write (N0chem,*) model_data%chem%N0%value
      Write (messages(1),'(1x,a)') 'The procedure is set to identify (and track) '//Trim(Adjustl(N0chem))//' changing chemical& 
                            & species along the trajectory. The composition and bonding of such species is defined as follows:'
      Write (messages(2),'(1x,a)') '- species will contain bonds to "'//Trim(model_data%chem%bonds%species%type)//'" atoms' 
      Write (messages(3),'(1x,a,(*(3x,a)))') '- a bond with a "'//Trim(model_data%chem%bonds%species%type)//'" atom&
                                        & is only possible with following atomic tag(s): ', &
                                        & (Trim(model_data%chem%acceptor%tg_incl(j)), j = 1, model_data%chem%acceptor%N0_incl)
      Write (nbonds,*) model_data%chem%bonds%N0%value
      Write (ncutoff,'(f10.2)') model_data%chem%bonds%cutoff%value
      Write (messages(4),'(1x,a)') '- the number of bonds for each species is '//Trim(Adjustl(nbonds))//&
                                  & '. A bond is created when the interactomic distance (cutoff) is lower than '&
                                  &//Trim(Adjustl(ncutoff))//' Angstroms'
      Call info(messages, 4)
      
      If (model_data%extra_bonds%invoke%fread) Then
        Write (messages(1),'(1x,a)') 'From the definitions in "&possible_extra_bonds" (within the &bonding_criteria block),&
                                    & the existence of the following possible bonds will also be considered:'
        Write (messages(2),'(2x,a)') '-----------------------------'                            
        Write (messages(3),'(2x,a)') 'Tag1  Tag2  Cutoff (Angstrom)'
        Write (messages(4),'(2x,a)') '-----------------------------'
        Call info(messages, 4)                            
        Do j = 1, model_data%extra_bonds%N0
           Write (messages(1),'(2x,(2(a,3x),f6.2))') Trim(model_data%extra_bonds%tg1(j)),&
                                                   & Trim(model_data%extra_bonds%tg2(j)),&
                                                   & model_data%extra_bonds%bond(j)%value
          Call info(messages, 1)                                                                    
        End Do                    
        Write (messages(1),'(2x,a)') '-----------------------------'
        Call info(messages, 1)                            
      End If

      Write (messages(1),'(1x,a)') 'Atomic species of donors will be tagged with the "*" symbol'
      Write (messages(2),'(1x,a)') 'The algorithm tracks bond breaking and formation with "'//&
                                   &Trim(model_data%chem%bonds%species%type)//'" atoms by considering neighbouring acceptors&
                                   & within the environment around each donor species:'
      Write (ncutoff,'(f10.2)') model_data%chem%acceptor%cutoff%value                             
      Write (messages(3),'(1x,a)') '- the environment around each donor is defined as the region within a cutoff radius of '&
                                    &//Trim(Adjustl(ncutoff))//' Angstroms'
      Write (messages(4),'(1x,a,(*(3x,a)))') '- the search for new donors is restricted to the following&
                                    & atomic tag(s): ',&
                                    & (Trim(model_data%chem%acceptor%tg_incl(j)), j = 1, model_data%chem%acceptor%N0_incl)
      Call info(messages, 4)
      If (model_data%chem%acceptor%info_exclude%fread) Then
        Write (messages(1),'(1x,a,(*(3x,a)))') '- the search will also exclude the following possible acceptor--donor pairs(s): ', &
         &  (Trim(model_data%chem%acceptor%tg_excl(j))//'*--'&
         &//Trim(model_data%chem%acceptor%tg_excl(j)), j = 1, model_data%chem%acceptor%N0_excl)
        Call info(messages, 1)
      End If
    End If

  End Subroutine print_model_settings

End Module atomic_model
