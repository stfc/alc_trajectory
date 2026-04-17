!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to obtain statistics for:
!  * intermolecular angles and distances
!  * intramolecular angles and distances  
!  * distribution of the NN distances
!  * coordinate distribution for a given species
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:   -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module geo_stat

  Use atomic_model,  Only: model_type, &
                           check_PBC,&
                           check_length_directive, &
                           geo_param_type

  Use constants,     Only: Rads_to_degrees, &
                           max_components, &
                           max_at_species
                            
  Use fileset,       Only: file_type,  &
                           refresh_out, &
                           FILE_INTERMOL_DISTANCES_NN1, &
                           FILE_INTERMOL_DISTANCES_NN2, &
                           FILE_INTERMOL_ANGLES, &
                           FILE_INTRAMOL_DISTANCES, &
                           FILE_INTRAMOL_ANGLES, &
                           FILE_SELECTED_NN_DISTANCES, &
                           FILE_COORD_DISTRIB, &
                           FILE_SET

  Use input_types,  Only : in_param,   & 
                           in_string                           

  Use numprec,       Only: wi,& 
                           wp
                         
  Use trajectory,    Only: traj_type, &
                           within_region

  Use process_data,  Only: remove_symbols, &
                           set_read_status, &
                           check_for_rubbish,&
                           capital_to_lower_case,&
                           get_word_length, &
                           prevent_segmentation

  Use unit_output,   Only: info, &
                           error_stop 
                           
  ! Type for nn_distance
  Type :: nndist_type
    Private
    Type(in_string), Public  :: invoke
    Type(in_string)  :: tag_reference_species
    Type(in_string)  :: tag_nn_species
    Character(Len=8) :: reference_species
    Integer(Kind=wi) :: num_nn_species
    Character(Len=8) :: nn_species(max_at_species)
    Type(in_param)   :: lower_bound
    Type(in_param)   :: upper_bound
    Type(in_param)   :: dr 
  End Type
  
  !Type to desc`ribe the coordinate distribution
  Type :: coord_distr_type
    Type(in_string), Public  :: invoke
    Type(in_string)  :: species_dir
    Character(Len=8) :: species
    Type(in_string)  :: coordinate
    Type(in_param)   :: delta
    Integer(Kind=wi) :: indx
  End Type
    
  Type, Public :: geo_stat_type
    Private
    Type(nndist_type),      Public :: nndist
    Type(coord_distr_type), Public :: coord_distr
  End Type geo_stat_type
  
  Public :: read_coord_distrib, read_selected_nn_distances
  Public :: check_coord_distrib, check_selected_nn_distances
  Public :: obtain_intermol_geom_stat, obtain_intramol_geom_stat, compute_number_monitored_species
  Public :: compute_coordinate_distribution, compute_nn_distance_distribution
  
Contains

  Subroutine read_selected_nn_distances(iunit, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read parameters from the
    ! &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   ) :: iunit
    Type(geo_stat_type), Intent(InOut) :: geo_stat_data 
    
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
        Read (iunit, Fmt=*, iostat=io) geo_stat_data%nndist%tag_reference_species%type, geo_stat_data%nndist%reference_species
        Call set_read_status(word, io, geo_stat_data%nndist%tag_reference_species%fread, & 
                           & geo_stat_data%nndist%tag_reference_species%fail, geo_stat_data%nndist%tag_reference_species%type)

      Else If (Trim(word)=='nn_species') Then
        Read (iunit, Fmt=*, iostat=io) geo_stat_data%nndist%tag_nn_species%type, geo_stat_data%nndist%num_nn_species
        Call prevent_segmentation(iunit, io, geo_stat_data%nndist%tag_nn_species%type, geo_stat_data%nndist%num_nn_species,&
                                & 'max_components', max_components, set_error)
        geo_stat_data%nndist%nn_species=' '
        Read (iunit, Fmt=*, iostat=io) geo_stat_data%nndist%tag_nn_species%type, geo_stat_data%nndist%num_nn_species,&
                                     & (geo_stat_data%nndist%nn_species(i), i=1, geo_stat_data%nndist%num_nn_species) 
        Call set_read_status(word, io, geo_stat_data%nndist%tag_nn_species%fread, geo_stat_data%nndist%tag_nn_species%fail,&
                                     & geo_stat_data%nndist%tag_nn_species%type)
        
      Else If (Trim(word)=='lower_bound') Then
         Read (iunit, Fmt=*, iostat=io) geo_stat_data%nndist%lower_bound%tag, geo_stat_data%nndist%lower_bound%value,&
                                     & geo_stat_data%nndist%lower_bound%units 
         Call set_read_status(word, io, geo_stat_data%nndist%lower_bound%fread, geo_stat_data%nndist%lower_bound%fail)

      Else If (Trim(word)=='upper_bound') Then
         Read (iunit, Fmt=*, iostat=io) geo_stat_data%nndist%upper_bound%tag, geo_stat_data%nndist%upper_bound%value,&
                                      & geo_stat_data%nndist%upper_bound%units 
         Call set_read_status(word, io, geo_stat_data%nndist%upper_bound%fread, geo_stat_data%nndist%upper_bound%fail)

      Else If (Trim(word)=='dr') Then
         Read (iunit, Fmt=*, iostat=io) geo_stat_data%nndist%dr%tag, geo_stat_data%nndist%dr%value, geo_stat_data%nndist%dr%units 
         Call set_read_status(word, io, geo_stat_data%nndist%dr%fread, geo_stat_data%nndist%dr%fail)

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
  
  Subroutine check_selected_nn_distances(files, model_data, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the definition of the
    ! parameters defined in the &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(In   ) :: files(:)
    Type(model_type),    Intent(In   ) :: model_data
    Type(geo_stat_type), Intent(InOut) :: geo_stat_data 

    Character(Len=256)  :: messages(3), error_set
    Character(Len=8)    :: tagj, tagk
    Logical             :: flag
    Integer(Kind=wi)    :: j, k
    
    ! Error message just in case....
    error_set = '***ERROR in file '//Trim(files(FILE_SET)%filename)//' -'
    Write (messages(1),'(1x,2a)')  Trim(error_set), ' "&selected_nn_distances" block.'

    If (geo_stat_data%nndist%tag_nn_species%fread) Then
      If (geo_stat_data%nndist%tag_nn_species%fail) Then
        Write (messages(2),'(1x,a)')  'Problems to define the "nn_species" directive'  
        call info(messages, 2)
        call error_stop(' ')
      End If
    Else
      Write (messages(2),'(1x,a)')  'The user must define the "nn_species" directive'  
      call info(messages, 2)
      call error_stop(' ')
    End If

    If (geo_stat_data%nndist%tag_reference_species%fread) Then
      If (geo_stat_data%nndist%tag_reference_species%fail) Then
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
    tagk=Trim(geo_stat_data%nndist%reference_species)
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
      Write (messages(2),'(1x,a)')   'The tag "'//Trim(geo_stat_data%nndist%reference_species)//&
                                     &'" defined in the "reference_species" directive is not a valid option.&
                                     & Please check the definition of the &input_composition block' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If 
    
    !Check if tags in nn_species are defined in the &input_composition block  
    Do k=1, geo_stat_data%nndist%num_nn_species
      tagk=Trim(geo_stat_data%nndist%nn_species(k))
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
        Write (messages(2),'(1x,a)')   'The tag "'//Trim(geo_stat_data%nndist%nn_species(k))//&
                                       &'" defined in the "nn_species" directive is not a valid option.&
                                       & Please check the definition of the &input_composition block' 
        Call info(messages, 2)
        Call error_stop(' ') 
      End If 
    End Do

    !Check if tags defined in nn_species are repeated
    Do j=1, geo_stat_data%nndist%num_nn_species-1
      tagj=Trim(geo_stat_data%nndist%nn_species(j))
      Do k=j+1, geo_stat_data%nndist%num_nn_species 
        tagk=Trim(geo_stat_data%nndist%nn_species(k))
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
    If (.Not. geo_stat_data%nndist%lower_bound%fread) Then
      geo_stat_data%nndist%lower_bound%tag='lower_bound'
    End If
    If (.Not. geo_stat_data%nndist%upper_bound%fread) Then
      geo_stat_data%nndist%upper_bound%tag='upper_bound'
    End If
    If (.Not. geo_stat_data%nndist%dr%fread) Then
      geo_stat_data%nndist%dr%tag='delta'
    End If
    
    Call check_length_directive(geo_stat_data%nndist%lower_bound, messages(1), .True., 'directive')
    Call check_length_directive(geo_stat_data%nndist%upper_bound, messages(1), .True., 'directive')
    Call check_length_directive(geo_stat_data%nndist%dr,          messages(1), .True., 'directive')
    If (geo_stat_data%nndist%lower_bound%value >= geo_stat_data%nndist%upper_bound%value) Then
      Write (messages(2),'(1x,a)')  'The value of "upper_bound" must be larger than "lower_bound"&
                                  & (make sure this is the case if you use different units)' 
      Call info(messages, 2)
      Call error_stop(' ') 
    End If
     
  End Subroutine check_selected_nn_distances

  Subroutine read_coord_distrib(iunit, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for computing the coordinate distribution
    ! of selective species. Information must be provided in the 
    ! &coord_distrib block 
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   ) :: iunit
    Type(geo_stat_type), Intent(InOut) :: geo_stat_data 

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
        Read (iunit, Fmt=*, iostat=io) geo_stat_data%coord_distr%species_dir%type, geo_stat_data%coord_distr%species
        Call set_read_status(word, io, geo_stat_data%coord_distr%species_dir%fread,&
                           & geo_stat_data%coord_distr%species_dir%fail,geo_stat_data%coord_distr%species_dir%type)

      Else If (Trim(word)=='delta') Then
         Read (iunit, Fmt=*, iostat=io) geo_stat_data%coord_distr%delta%tag, &
                                      & geo_stat_data%coord_distr%delta%value,&
                                      & geo_stat_data%coord_distr%delta%units 
         Call set_read_status(word, io, geo_stat_data%coord_distr%delta%fread, geo_stat_data%coord_distr%delta%fail)

      Else If (Trim(word)=='coordinate') Then
        Read (iunit, Fmt=*, iostat=io) word, geo_stat_data%coord_distr%coordinate%type
        Call set_read_status(word, io, geo_stat_data%coord_distr%coordinate%fread,& 
                           & geo_stat_data%coord_distr%coordinate%fail,&
                           & geo_stat_data%coord_distr%coordinate%type)
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_coord_distrib"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_coord_distrib

  Subroutine check_coord_distrib(files, model_data, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &coord_distrib block
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(In   ) :: files(:)
    Type(model_type),     Intent(In   ) :: model_data
    Type(geo_stat_type),  Intent(InOut) :: geo_stat_data

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

    If (.Not. geo_stat_data%coord_distr%delta%fread) Then
      geo_stat_data%coord_distr%delta%tag='delta'
    End If
    Call check_length_directive(geo_stat_data%coord_distr%delta, error_set, .True., 'directive')

    ! Check definition of "species" directive
    If (geo_stat_data%coord_distr%species_dir%fread) Then
      If (geo_stat_data%coord_distr%species_dir%fail) Then
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
    tg=Trim(geo_stat_data%coord_distr%species)
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
      Write (messages(1),'(2(1x,a))') Trim(error_set), '"'//Trim(geo_stat_data%coord_distr%species)//'"&
                                     & defined for the "species" directive is not a valid species.&
                                     & Please review the definition of the &input_composition block' 
      Call info(messages, 1)
      Call error_stop(' ') 
    End If 

    ! Check definition of "species" directive
    If (geo_stat_data%coord_distr%species_dir%fread) Then
      If (geo_stat_data%coord_distr%species_dir%fail) Then
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
    If (geo_stat_data%coord_distr%coordinate%fread) Then
      If (geo_stat_data%coord_distr%coordinate%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "coordinate" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        flag=.True. 
        Do j=1, 3
          If (Trim(geo_stat_data%coord_distr%coordinate%type)==Trim(coord(j))) Then
            flag=.False.
            geo_stat_data%coord_distr%indx=j
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
  
  Subroutine compute_coordinate_distribution(files, model_data, traj_data, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the distribution of the coordinates (x, y or z) of
    ! the species selected in the &coord_distrib block
    !
    ! author    - i.scivetti Oct 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(model_type),    Intent(In   ) :: model_data
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(geo_stat_type), Intent(InOut) :: geo_stat_data

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
    If (vector(geo_stat_data%coord_distr%indx) > 0.0_wp) Then
      clim(2)=vector(geo_stat_data%coord_distr%indx)
      clim(1)=0.0_wp
      ctap='top'
    Else
      clim(1)=vector(geo_stat_data%coord_distr%indx) 
      clim(2)=0.0_wp 
      ctap='bottom'
    End If
    
    Do i = traj_data%analysis%frame_ini+1, traj_data%analysis%frame_last
      vector=0.0_wp
      Do j = 1, 3
        vector(:)=vector(:)+traj_data%box(i)%cell(j,:)
      End Do
      If (Trim(ctap)=='bottom') Then
        If (vector(geo_stat_data%coord_distr%indx) < clim(1)) Then
          clim(1)=vector(geo_stat_data%coord_distr%indx)
        End If
      Else If (Trim(ctap)=='top') Then
        If (vector(geo_stat_data%coord_distr%indx) > clim(2)) Then
          clim(2)=vector(geo_stat_data%coord_distr%indx)
        End If
      End If
    End Do
     
    ! Define number of bins
    nbins=Floor(Abs(clim(1)-clim(2))/geo_stat_data%coord_distr%delta%value)

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
          If (geo_stat_data%coord_distr%species==traj_data%config(i,j)%tag) Then
            num_at=num_at+1
            list_indx(num_at)=j
          End If
        End Do
      
        ! Calculate the histogram for this particular frame of the trajectory
        If (num_at/=0) Then
          h=0
          Do j=1, num_at
            coord_value=Abs(traj_data%config(i,list_indx(j))%r(geo_stat_data%coord_distr%indx))
            If (Trim(ctap)=='top') Then
              m=Floor(coord_value/geo_stat_data%coord_distr%delta%value)+1
            Else
              m=nbins+1-(Floor(coord_value/geo_stat_data%coord_distr%delta%value)+1)
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
        d(m)=d(m)/net_frames/geo_stat_data%coord_distr%delta%value      
      End Do
      
      ! Print results
      If (accum /= 0) Then
        ! Print File
        Open(Newunit=files(FILE_COORD_DISTRIB)%unit_no, File=files(FILE_COORD_DISTRIB)%filename, Status='Replace')
        iunit=files(FILE_COORD_DISTRIB)%unit_no
        Write (iunit,'(a)') '#  Distribution of the '//Trim(geo_stat_data%coord_distr%coordinate%type)//'-coordinate&
                           & for the "'//Trim(geo_stat_data%coord_distr%species)//'" species'
        Write (iunit,'(a)') '#  Value [Angstrom]      Probability [1/Angstrom]' 
        
        Do m=1, nbins
          Write (iunit,'(2x,f12.4,6x,f14.5)') (Real(m,Kind=wp)-0.5)*geo_stat_data%coord_distr%delta%value, d(m)
        End Do
        Write (message,'(1x,a)') 'The distribution of the '//Trim(geo_stat_data%coord_distr%coordinate%type)//&
                                &'-coordinate for the "'//Trim(geo_stat_data%coord_distr%species)//'" species was&
                                & printed to the "'//Trim(files(FILE_COORD_DISTRIB)%filename)//'" file.'
        Call info(message, 1)
      Else
        type_error=Trim(geo_stat_data%coord_distr%species)
        Write (messages(1),'(1x,a)') '*************************************************************************************'
        Call info(messages, 1)
        Write (messages(1),'(1x,a)') '   WARNING: coordinate distribution analysis could not be executed'
          Write (messages(2),'(1x,a)') '   Requested species '//Trim(type_error)//' as specified in the &coord_distrib&
                                  & block could not be identified along the trajectory.'
          Write (messages(3),'(1x,a)') '   Please verify the settings for the &coord_distrib block.'                 
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
    
  End Subroutine compute_coordinate_distribution

  Subroutine compute_nn_distance_distribution(files, model_data, traj_data, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the statistics of the shortest distance between
    ! the "reference_species" and those defined by the "nn_species" directive
    ! in the distance domain defined by the lower_bound and upper_bound
    ! Analysis is performed using the definitions of the &selected_nn_distances block
    !
    ! author    - i.scivetti Nov 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),       Intent(InOut) :: files(:)
    Type(model_type),      Intent(In   ) :: model_data
    Type(traj_type),       Intent(InOut) :: traj_data
    Type(geo_stat_type),   Intent(InOut) :: geo_stat_data

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
    nbins=Nint(Abs(geo_stat_data%nndist%upper_bound%value-geo_stat_data%nndist%lower_bound%value)/geo_stat_data%nndist%dr%value)
   
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
          If (geo_stat_data%nndist%reference_species==traj_data%config(i,j)%tag) Then
              num_at(1)=num_at(1)+1
              list_indx(num_at(1),1)=j
          End If 
          k=1
          flag=.False.
          Do While (k <= geo_stat_data%nndist%num_nn_species .And. (.Not. flag))
            If (geo_stat_data%nndist%nn_species(k)==traj_data%config(i,j)%tag) Then
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
              If (rmin >= geo_stat_data%nndist%lower_bound%value .And. rmin <= geo_stat_data%nndist%upper_bound%value) Then
                mk=Floor((rmin-geo_stat_data%nndist%lower_bound%value)/geo_stat_data%nndist%dr%value)+1
                If (mk <= nbins) Then
                  h(mk)=h(mk)+1
                  num_var=num_var+1
                End If
              End If
            End If
          End Do
        End If
        
        If (num_var /= 0) Then
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
           d(mk)=d(mk)/net_frames/geo_stat_data%nndist%dr%value
         End Do
       
        ! Print File
        Open(Newunit=files(FILE_SELECTED_NN_DISTANCES)%unit_no, File=files(FILE_SELECTED_NN_DISTANCES)%filename,&
                          &Status='Replace')
        iunit=files(FILE_SELECTED_NN_DISTANCES)%unit_no
        Write (iunit,'(a)') '#  Probability distribution of the shortest distances between'
        Write (iunit,'(a)') '#  the reference species "'//Trim(geo_stat_data%nndist%reference_species)//&
                             &'" and the species defined in the "nn_species" directive'    
        Write (iunit,'(a)') '#  Distance [Angstrom]     Probability [1/Angstrom]' 
        Do mk=1, nbins
          Write (iunit,'(2x,f12.4,6x,f14.5)') (Real(mk,Kind=wp)-0.5_wp)*geo_stat_data%nndist%dr%value+&
                                            & geo_stat_data%nndist%lower_bound%value, d(mk)
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
    
  End Subroutine compute_nn_distance_distribution

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
        
        If (num_var /= 0) Then
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
            Write (iunit,'(2x,f12.4,6x,f14.5)') (Real(mk,Kind=wp)-0.5)*M%delta%value+M%lower_bound%value, d(mk)
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

        If (num_var /= 0) Then
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
            Write (iunit,'(2x,f12.4,6x,f13.5)') (Real(mk,Kind=wp)-0.5)*M%delta%value+M%lower_bound%value, d(mk)
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
  
End Module geo_stat
