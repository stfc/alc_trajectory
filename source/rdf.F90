!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for Radial Distribution Function (RDF)
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module rdf

  Use atomic_model,  Only: model_type, &
                           check_PBC, &
                           check_length_directive

  Use constants,     Only: max_components,&
                           pi
                             
                           
  Use fileset,       Only: file_type,  &
                           refresh_out, &
                           FILE_RDF,&
                           FILE_SET
                    
  Use input_types,   Only: in_param,   & 
                           in_string
                     
  Use numprec,       Only: wi,&
                           wp
                     
  Use process_data,  Only: capital_to_lower_case, &
                           check_for_rubbish,&
                           get_word_length, &
                           prevent_segmentation, &
                           remove_symbols, &
                           set_read_status

  Use trajectory,    Only: traj_type, &
                           within_region

  Use unit_output,   Only: info, &
                           error_stop                            
                           
  !Type to describe the rdf
  Type, Public :: rdf_type
    Private
    Type(in_string), Public  :: invoke
    Type(in_string)  :: tags_species_a
    Type(in_string)  :: tags_species_b
    Type(in_param),  Public  :: dr
    Real(Kind=wp)    :: rmax
    Character(Len=8) :: type_a(max_components)
    Character(Len=8) :: type_b(max_components)
    Integer(Kind=wi) :: num_type_a
    Integer(Kind=wi) :: num_type_b
  End Type

  Public :: radial_distribution_function, read_rdf, check_rdf
  
Contains

  Subroutine read_rdf(iunit, rdf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for Radial Distribution Function (RDF)
    ! analysis from the &RDF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi), Intent(In   ) :: iunit
    Type(rdf_type),   Intent(InOut) :: rdf_data 

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
        Read (iunit, Fmt=*, iostat=io) rdf_data%tags_species_a%type, rdf_data%num_type_a
        Call prevent_segmentation(iunit, io, rdf_data%tags_species_a%type, rdf_data%num_type_a,&
                                & 'max_components', max_components, set_error)
        rdf_data%type_a=' '
        Read (iunit, Fmt=*, iostat=io) rdf_data%tags_species_a%type, rdf_data%num_type_a,&
                                       (rdf_data%type_a(j), j= 1, rdf_data%num_type_a)
        Call set_read_status(word, io, rdf_data%tags_species_a%fread, rdf_data%tags_species_a%fail,&
                           & rdf_data%tags_species_a%type)

      Else If (Trim(word)=='tags_species_b') Then
        Read (iunit, Fmt=*, iostat=io) rdf_data%tags_species_b%type, rdf_data%num_type_b
        Call prevent_segmentation(iunit, io, rdf_data%tags_species_b%type, rdf_data%num_type_b,&
                                & 'max_components', max_components, set_error)
        rdf_data%type_b=' '
        Read (iunit, Fmt=*, iostat=io) rdf_data%tags_species_b%type, rdf_data%num_type_b,&
                                       (rdf_data%type_b(j), j= 1, rdf_data%num_type_b)
        Call set_read_status(word, io, rdf_data%tags_species_b%fread, rdf_data%tags_species_b%fail,&
                           & rdf_data%tags_species_b%type)

      Else If (Trim(word)=='dr') Then
         Read (iunit, Fmt=*, iostat=io) rdf_data%dr%tag, rdf_data%dr%value, rdf_data%dr%units 
         Call set_read_status(word, io, rdf_data%dr%fread, rdf_data%dr%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_rdf"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_rdf

  Subroutine check_rdf(files, rdf_data, model_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &RDF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(rdf_type),     Intent(InOut) :: rdf_data
    Type(model_type),   Intent(In   ) :: model_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set
    Integer(Kind=wi)    :: j, k
    Logical             :: flag

    Character(Len=8)  :: tg(max_components)
    Character(Len=8)  :: el(max_components)

    error_set = '***ERROR in the &RDF block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (.Not. rdf_data%dr%fread) Then
      rdf_data%dr%tag='dr'
    End If
    Call check_length_directive(rdf_data%dr, error_set, .True., 'directive')
    
    If (rdf_data%tags_species_a%fread) Then
      If (rdf_data%tags_species_a%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "tags_species_a" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      Do j=1, rdf_data%num_type_a-1
        Do k=j+1, rdf_data%num_type_a
          If (rdf_data%type_a(j)==rdf_data%type_a(k)) Then
            Write (messages(1),'(4(1x,a))') Trim(error_set), 'Tag', Trim(rdf_data%type_a(j)), 'is repeated in the list!'
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
    
    If (rdf_data%tags_species_b%fread) Then
      If (rdf_data%tags_species_b%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "tags_species_b" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
      Do j=1, rdf_data%num_type_b-1
        Do k=j+1, rdf_data%num_type_b
          If (rdf_data%type_b(j)==rdf_data%type_b(k)) Then
            Write (messages(1),'(4(1x,a))') Trim(error_set), 'Tag', Trim(rdf_data%type_b(j)), 'is repeated in the list!'
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
    Do k=1, rdf_data%num_type_a
      tg(k)=Trim(rdf_data%type_a(k))
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
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'The tag " '//Trim(rdf_data%type_a(k))//' " of the "tags_species_a"&
                                       & is not a valid option. Please review the definition of the &input_composition block' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If 
    End Do
    
    ! Check if all tags correspond to the same element (tybe b)
    Do k=1, rdf_data%num_type_b
      tg(k)=Trim(rdf_data%type_b(k))
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
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'The tag " '//Trim(rdf_data%type_b(k))//' " of the "tags_species_b"&
                                       & is not a valid option. Please review the definition of the &input_composition block' 
        Call info(messages, 1)
        Call error_stop(' ') 
      End If 
    End Do
    
  End Subroutine check_rdf

  Subroutine radial_distribution_function(files, traj_data, model_data, rdf_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the Radial Distribution Function (RDF) based on the
    ! settings of the &RDF block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
    Type(rdf_type),    Intent(InOut) :: rdf_data

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
    nbins=Floor(rmax/rdf_data%dr%value)

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
          Do While (itype <= rdf_data%num_type_a .And. ftype)
            If (rdf_data%type_a(itype)==traj_data%config(i,j)%tag) Then
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
          Do While (itype <= rdf_data%num_type_b .And. ftype)
            If (rdf_data%type_b(itype)==traj_data%config(i,j)%tag) Then
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
                m=Floor(norm2(rjk)/rdf_data%dr%value)+1
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
                                & (Trim(rdf_data%type_a(j)), j= 1, rdf_data%num_type_a) 
        Write (iunit,'(*(a,2x))') '#  Tags for type species "b":',&
                                & (Trim(rdf_data%type_b(j)), j= 1, rdf_data%num_type_b) 
        Write (iunit,'(a)') '#  Radius [A]      RDF [1/A^3]       Coordination number' 
        
        cn=0.0_wp
        Do m=1, nbins
          r_bin=Real(m, Kind=wp)*rdf_data%dr%value
          dV=4.0_wp*pi*r_bin**2*rdf_data%dr%value 
          gr(m)=gr(m)/(dV*net_frames)
          Do k= 1, m
           cn(m)=cn(m)+nn(k)
          End Do
          cn(m)=cn(m)/net_frames 
          !Write (iunit,'(2x, f11.3,(2(6x,f11.6)))') (Real(m,Kind=wp)-0.5)*rdf_data%dr%value, gr(m), cn(m)
          Write (iunit,'(2x, f11.3,(2(6x,f11.6)))') (1.0_wp*m-0.5)*rdf_data%dr%value, gr(m), cn(m)
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
         Write (messages(3),'(1x,a)') '   Please verify the settings for the &RDF and &region blocks'                   
        Else
          Write (messages(2),'(1x,a)') '   Requested species as specified for '//Trim(type_error)//' in the &RDF&
                                  & block could not be identified along the trajectory.'
          Write (messages(3),'(1x,a)') '   Please verify the settings for the &RDF block.'                        
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

  
End Module rdf
