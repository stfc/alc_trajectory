!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for calculation related to transfer of reactive species
!
! Copyright   2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module transfer_species

  Use atomic_model, Only : model_type

  Use fileset,       Only: file_type,  &
                           refresh_out, &
                           FILE_SPCF_ALL, &
                           FILE_SPCF_AVG, & 
                           FILE_TCF_ALL, &
                           FILE_TCF_AVG,&                                                 
                           FILE_RES_TIMES, &
                           FILE_SET

  Use input_types,   Only: in_logic, &
                           in_string, &
                           in_param
                           
  Use numprec,       Only: wi,& 
                           wp                           

  Use process_data,  Only: capital_to_lower_case, &
                           check_for_rubbish, &
                           set_read_status, &
                           get_word_length
                           
  Use trajectory,    Only: traj_type, &
                           average_segments,&
                           check_time_directive,&
                           find_active_bonds, &
                           within_region

  Use unit_output,   Only: info, &
                           error_stop                            
                           
  !Type to describe the region where to constrain the analysis
  Type :: transfer_type
    Type(in_string)  :: invoke
    Type(in_string)  :: method
    Type(in_param)   :: rattling_wait
    Type(in_logic)   :: print_all_intervals
  End Type

  Public :: read_transfer_species, check_transfer_species
  Public :: compute_lifetime_related_quantities
  
Contains  

  Subroutine read_transfer_species(iunit, transfer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to read the information from the &lifetime block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Integer(Kind=wi),    Intent(In   )  :: iunit
    Type(transfer_type), Intent(InOut)  :: transfer_data 

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
        Read (iunit, Fmt=*, iostat=io) word, transfer_data%method%type
        Call set_read_status(word, io, transfer_data%method%fread,&
                                     & transfer_data%method%fail,&
                                     & transfer_data%method%type)

      Else If (Trim(word)=='rattling_wait') Then
         Read (iunit, Fmt=*, iostat=io) transfer_data%rattling_wait%tag, &
                                      & transfer_data%rattling_wait%value,& 
                                      & transfer_data%rattling_wait%units
         Call set_read_status(word, io, transfer_data%rattling_wait%fread,&
                                      & transfer_data%rattling_wait%fail)

      Else If (word(1:length) == 'print_all_intervals') Then
       Read (iunit, Fmt=*, iostat=io) word, transfer_data%print_all_intervals%stat
       Call set_read_status(word, io, transfer_data%print_all_intervals%fread, transfer_data%print_all_intervals%fail)
                            
                                      
      Else
        Write (message,'(1x,5a)') Trim(set_error), ' Directive "', Trim(word),&
                                & '" is not recognised as a valid settings.',&
                                & ' See the "use_code.md" file. Have you properly closed the block with "&end_lifetime"?'
        Call error_stop(message)
      End If

    End Do
    
  End Subroutine read_transfer_species

  Subroutine check_transfer_species(files, transfer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to check the settings of the &lifetime block
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(In   ) :: files(:)
    Type(transfer_type), Intent(InOut) :: transfer_data

    Character(Len=256)  :: error_set
    Character(Len=256)  :: messages(2)

    error_set = '***ERROR in the &lifetime block of file '//Trim(files(FILE_SET)%filename)//' -'

    Call check_time_directive(transfer_data%rattling_wait, 'rattling_wait' ,error_set, .False.)

    If (.Not. transfer_data%rattling_wait%fread) Then
      transfer_data%rattling_wait%value= 0.0_wp
      transfer_data%rattling_wait%units= 'fs'
    End If
    
    If (transfer_data%method%fread) Then
      If (transfer_data%method%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Wrong (or missing) settings for the "method" directive.'
        Call info(messages, 1)
        Call error_stop(' ')
      Else
        If (Trim(transfer_data%method%type)/='hicf'     .And. &
            Trim(transfer_data%method%type)/='hdcf')  Then
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

    If (transfer_data%print_all_intervals%fread) Then
      If (transfer_data%print_all_intervals%fail) Then
        Write (messages(1),'(2(1x,a))') Trim(error_set), 'Missing (or wrong) specification for directive&
                                  & "print_all_intervals" (choose either .True. or .False.)'
        Call info(messages,1)
        Call error_stop(' ')
      End If
    Else
      transfer_data%print_all_intervals%stat=.False.
    End If
    
  End Subroutine check_transfer_species
  
  Subroutine compute_lifetime_related_quantities(files, traj_data, model_data, transfer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the transfer correlation functions and residence
    ! times
    !
    ! author    - i.scivetti April 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(model_type),    Intent(In   ) :: model_data
    Type(transfer_type), Intent(InOut) :: transfer_data    
    
    Call transfer_correlation_function_sites(files, traj_data, model_data, transfer_data)
    Call residence_times_sites(files, traj_data, model_data, transfer_data)

    If (.Not. traj_data%active_bonds_computed) Then
      Call find_active_bonds(traj_data, model_data)
      traj_data%active_bonds_computed=.True.
    End If
    
    Call special_pair_correlation_function(files, traj_data, model_data, transfer_data)
  
  End Subroutine compute_lifetime_related_quantities

  Subroutine transfer_correlation_function_sites(files, traj_data, model_data, transfer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the transfer correlation function (TCF) involving
    ! the changing chemistry species
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(model_type),    Intent(In   ) :: model_data
    Type(transfer_type), Intent(InOut) :: transfer_data

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
    
    method=Trim(transfer_data%method%type)
    
    If (transfer_data%print_all_intervals%stat) Then
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
            If (.Not. terminated(l,j)) Then
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
        If (transfer_data%print_all_intervals%stat) Then
          Write (iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          If (i==traj_data%analysis%seg_indx(2,k) .And. (traj_data%analysis%N_seg /=1)) Then
             If (k /= traj_data%analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)        TCF' 
             End If
          End If   
        End If
      End Do
    End Do
    
    If (transfer_data%print_all_intervals%stat) Then
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
    
    If (.Not. transfer_data%print_all_intervals%stat) Then
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
  
  Subroutine special_pair_correlation_function(files, traj_data, model_data, transfer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the correlation function of the special pair (SPCF),
    ! which is defined by the atomic pair of the chemical site and the closest NN, 
    ! either donor or acceptor
    !
    ! author    - i.scivetti April 2024
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(model_type),    Intent(In   ) :: model_data
    Type(transfer_type), Intent(InOut) :: transfer_data
    
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
    
    method=transfer_data%method%type
    
    If (transfer_data%print_all_intervals%stat) Then
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
            If (.Not. terminated(l,j)) Then
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
        If (transfer_data%print_all_intervals%stat) Then
          Write (iunit,'(f11.3, 4x, 1(f11.3))') (time-base_time)/1000.0_wp, suma_i
          If (i==traj_data%analysis%seg_indx(2,k) .And. (traj_data%analysis%N_seg /=1)) Then
             If (k /= traj_data%analysis%N_seg) Then
               Write (iunit,'(a)') '#  Time (ps)        SPCF' 
             End If
          End If   
        End If
      End Do
    End Do
    
    If (transfer_data%print_all_intervals%stat) Then
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
    
    If (.Not. transfer_data%print_all_intervals%stat) Then
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
  
  Subroutine residence_times_sites(files, traj_data, model_data, transfer_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to compute the residence times for the changing species
    !
    ! author    - i.scivetti March 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),     Intent(InOut) :: files(:)
    Type(traj_type),     Intent(InOut) :: traj_data
    Type(model_type),    Intent(In   ) :: model_data
    Type(transfer_type), Intent(InOut) :: transfer_data
    
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
    
    rattling=transfer_data%rattling_wait%value

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
          Write (iunit,'(f11.3,15x,a)') values(i,m), Trim(tag(i,1,m))
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
  
  
End Module transfer_species
