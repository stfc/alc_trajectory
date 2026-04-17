!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the Mean Square Displacement (MSD)
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module msd

  Use atomic_model,  Only: model_type, &
                           check_PBC
                           
  Use fileset,       Only: file_type,  &
                           refresh_out, &
                           FILE_MSD_ALL, &
                           FILE_MSD_AVG, &
                           FILE_SET

  
  Use input_types,   Only: in_logic, &
                           in_string
                           
  Use numprec,       Only: wi,& 
                           wp
 
  Use process_data, Only : set_read_status, &
                           get_word_length, &
                           capital_to_lower_case, &
                           check_for_rubbish  
 
  Use trajectory,    Only: traj_type, &
                           average_segments, &
                           within_region
                           
  Use unit_output,   Only: info, &
                           error_stop 

  !Type to describe the msd
  Type, Public :: msd_type
    Private
    Type(in_string), Public  :: invoke
    Type(in_string), Public  :: pbc_xyz
    Type(in_string), Public  :: select
    Type(in_logic),  Public  :: print_all_intervals
    Logical,         Public  :: pbc(3)
    Real(Kind=wp),   Public  :: r2
  End Type

  Public :: read_msd, check_msd, mean_squared_displacement 
  
Contains

  Subroutine read_msd(iunit, msd_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the settigns for mean square displacement (MSD)
    ! analysis from the &MSD block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),  Intent(In   ) :: iunit
    Type(msd_type),    Intent(InOut) :: msd_data
 

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
        Read (iunit, Fmt=*, iostat=io) word, msd_data%select%type
        Call set_read_status(word, io, msd_data%select%fread, msd_data%select%fail,&
                           & msd_data%select%type)

      Else If (Trim(word)=='pbc_xyz') Then
         Read (iunit, Fmt=*, iostat=io) word, (msd_data%pbc(j), j= 1, 3)
         Call set_read_status(word, io, msd_data%pbc_xyz%fread, msd_data%pbc_xyz%fail)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, msd_data%print_all_intervals%stat
       Call set_read_status(word, io, msd_data%print_all_intervals%fread, msd_data%print_all_intervals%fail)

      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_msd"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_msd
  
  Subroutine check_msd(files, msd_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &msd block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),    Intent(In   ) :: files(:)
    Type(msd_type),     Intent(InOut) :: msd_data

    Character(Len=256)  :: messages(2)
    Character(Len=64 )  :: error_set

    error_set = '***ERROR in the &MSD block of file '//Trim(files(FILE_SET)%filename)//' -'

    If (msd_data%select%fread) Then
      If (msd_data%select%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "select" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(msd_data%select%type)/='x'  .And. &
            Trim(msd_data%select%type)/='y'  .And. &
            Trim(msd_data%select%type)/='z'  .And. &
            Trim(msd_data%select%type)/='xy' .And. &
            Trim(msd_data%select%type)/='xz' .And. &
            Trim(msd_data%select%type)/='yz' .And. &
            Trim(msd_data%select%type)/='xyz') Then
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

    If (msd_data%pbc_xyz%fread) Then
      If (msd_data%pbc_xyz%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "pbc_xyz" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      End If
    Else
      msd_data%pbc=.True.
    End If

    If (msd_data%print_all_intervals%fread) Then
      If (msd_data%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      msd_data%print_all_intervals%stat=.False.
    End If
    
  End Subroutine check_msd

  
  Subroutine mean_squared_displacement(files, traj_data, model_data, msd_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the Mean Squared Displacement (MSD) based on the
    ! settings of the &MSD block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),   Intent(InOut) :: files(:)
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(model_type),  Intent(In   ) :: model_data
    Type(msd_type),    Intent(InOut) :: msd_data

    Integer(Kind=wi)  :: i, j, k, l
    Integer(Kind=wi)  :: Nini_species, iunit, indx
    Real(Kind=wp)     :: time
    Real(Kind=wp)     :: base_time
    Logical           :: set_u0

    Character(Len=256) :: message, num_species
    Logical           :: terminated(traj_data%Nmax_species)

    If (msd_data%print_all_intervals%stat) Then
      ! Print tracked species
      Open(Newunit=files(FILE_MSD_ALL)%unit_no, File=files(FILE_MSD_ALL)%filename, Status='Replace')
      iunit=files(FILE_MSD_ALL)%unit_no
      Write (num_species,*) model_data%chem%N0%value
      Write (iunit,'(a)') '#  MSD analysis for all the intervals' 
      Write (iunit,'(a,8x,a)') '#  Time (ps)', '"'//Trim(msd_data%select%type)//'"-MSD for species "'&
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
      
        msd_data%r2=0.0_wp
        traj_data%N_species=0
        Do j = 1, traj_data%Nmax_species
          If (.Not. terminated(j)) Then
            If (traj_data%species(i,j)%alive) Then
              indx=traj_data%species(i,j)%list(1)
              traj_data%species(i,j)%u(:,1)=traj_data%config(i,indx)%r
              Call msd_vector_difference(traj_data, msd_data, i, j)
            Else
              terminated(j)=.True.
            End If
          End If  
        End Do
        
        If (i==traj_data%analysis%seg_indx(2,k)) Then
          If (traj_data%N_species /= 0) Then
            msd_data%r2=msd_data%r2/traj_data%N_species
            traj_data%analysis%variable(l,k)=msd_data%r2
            If (msd_data%print_all_intervals%stat) Then
              Write (iunit,'(f11.3, 4x, f11.4)') (time-base_time)/1000.0_wp, msd_data%r2
            End If  
          End If  
          terminated=.False.
          set_u0=.True.
          If (msd_data%print_all_intervals%stat) Then
            If ((traj_data%analysis%N_seg /=1) .And. (k /= traj_data%analysis%N_seg)) Then
              If (k /= traj_data%analysis%N_seg) Then
                Write (iunit,'(a,8x,a)') '#  Time (ps)', '"'//Trim(msd_data%select%type)//'"-MSD for species "'&
                                  &//Trim(model_data%species_definition%name%type)//'" [Angstrom^2]' 
              End If                    
            End If
          End If
        Else
          If ((traj_data%N_species) /= 0) Then
            msd_data%r2=msd_data%r2/traj_data%N_species
          Else  
            msd_data%r2=0.0_wp
          End If
            traj_data%analysis%variable(l,k)=msd_data%r2
            If (msd_data%print_all_intervals%stat) Then
              Write (iunit,'(f11.3, 4x, f11.4)') (time-base_time)/1000.0_wp, msd_data%r2
            End If
        End If  
      End Do
    End Do

    If (msd_data%print_all_intervals%stat) Then
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

    Call average_segments(files, traj_data, FILE_MSD_AVG, 'MSD', msd_data%select%type)
    If (traj_data%analysis%N_seg ==1 ) Then 
      Write (message,'(1x,a)') 'WARNING: A single time interval was used to compute the average MSD! The computed STD&
                              & is zero. Use/Check the &data_analysis block to improve the statistics.'
      Call info(message, 1)
    End If

    If (.Not. msd_data%print_all_intervals%stat) Then
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

  Subroutine msd_vector_difference(traj_data, msd_data, i, j)  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the contribution to the correlation from species j
    ! for the frame i (cij)
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(traj_type),   Intent(InOut) :: traj_data
    Type(msd_type),    Intent(InOut) :: msd_data
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
      If (msd_data%pbc_xyz%fread) Then
        Do m= 1, 3
          If (.Not. msd_data%pbc(m)) Then
            du(m)=traj_data%species(i,j)%u(m,1)-traj_data%species(i,j)%u0(m,1)   
          End If
        End Do
      End If

      Select Case (Trim(msd_data%select%type))  
        Case ('x')
          msd_data%r2 = msd_data%r2 + du(1)**2
        Case ('y')
          msd_data%r2 = msd_data%r2 + du(2)**2
        Case ('z')
          msd_data%r2 = msd_data%r2 + du(3)**2
        Case ('xy')
          msd_data%r2 = msd_data%r2 + du(1)**2 + du(2)**2
        Case ('xz')
          msd_data%r2 = msd_data%r2 + du(1)**2 + du(3)**2
        Case ('yz')
          msd_data%r2 = msd_data%r2 + du(2)**2 + du(3)**2
        Case ('xyz')
          Do m= 1, 3
           msd_data%r2 = msd_data%r2 + du(m)**2
          End Do
      End Select  

    End If
    
  End Subroutine msd_vector_difference

End Module msd
