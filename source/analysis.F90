!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module related to the analysis of the generated trajectory
!
! Copyright   2023-2026 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     -  i.scivetti  Feb 2026           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module analysis

  Use atomic_model,     Only: model_type
                        
  Use fileset,          Only: file_type,  &
                              refresh_out, &
                              FILE_INTERMOL_DISTANCES_NN1, &
                              FILE_INTERMOL_DISTANCES_NN2
                              
  Use geo_stat,         Only: geo_stat_type, &
                              compute_coordinate_distribution, & 
                              compute_nn_distance_distribution, &
                              obtain_intermol_geom_stat, &
                              obtain_intramol_geom_stat,& 
                              compute_number_monitored_species  
                           
  Use transfer_species, Only: transfer_type,&
                              compute_lifetime_related_quantities                         

  Use numprec,          Only: wi,& 
                              wp
                        
  Use msd,              Only: msd_type, & 
                              mean_squared_displacement
                              
  Use ocf,              Only: ocf_type, &                            
                              chemocf_type, &
                              compute_orientational_chemistry, &
                              orientational_correlation_function
                        
  Use rdf,              Only: rdf_type, &
                              radial_distribution_function   
                              
  Use trajectory,       Only: traj_type, &
                              define_trajectory_segments, &
                              find_active_bonds,&
                              find_neighbours_monitored_species,&
                              print_tracking_species, &
                              print_unchanged_chemistry
                              
  Use unit_output,      Only: info, &
                              error_stop 

                           
  Public :: trajectory_analysis, print_settings_for_trajectory_analysis
  
Contains

  Subroutine print_settings_for_trajectory_analysis(files, model_data, traj_data, ocf_data, chemocf_data,&
                                                  & msd_data, rdf_data, transfer_data, geo_stat_data) 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print a summary of the trajectory settings
    !
    ! author    - i.scivetti Feb 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(model_type),     Intent(In   ) :: model_data
    Type(traj_type),      Intent(InOut) :: traj_data
    Type(ocf_type),       Intent(InOut) :: ocf_data
    Type(chemocf_type),   Intent(InOut) :: chemocf_data
    Type(msd_type),       Intent(InOut) :: msd_data
    Type(rdf_type),       Intent(InOut) :: rdf_data
    Type(transfer_type),  Intent(InOut) :: transfer_data
    Type(geo_stat_type),  Intent(InOut) :: geo_stat_data
    
    Character(Len=256) :: messages(4), word
    Integer(Kind=wi)   :: k, j
    
    Call info(' ', 1) 
    Call info('Trajectory settings', 1) 
    Call info('===================', 1) 

    ! Check time settings (from &data_analysis) against the trajectory 
    Call define_trajectory_segments(files, traj_data)

    Write (messages(1),'(1x,a)') '- by specification, the trajectory corresponds to the "'&
                                 &//Trim(traj_data%ensemble%type)//'" ensemble and was recorded in "'&
                                 &//Trim(model_data%input_geometry_format%type)//'" format'
    Call info(messages, 1)
    Write (word,'(f10.2)') traj_data%timestep%value
    Write (messages(1),'(1x,a,f5.2,a)') '- the time step between recorded configurations is '//Trim(Adjustl(word))//' fs'
    Call info(messages, 1)

    If (traj_data%analysis%end_time%fread) Then
      If ((traj_data%analysis%end_time%value-(traj_data%frames-1)*traj_data%timestep%value)<0.0_wp) Then
        Write (word,'(f8.3)') traj_data%analysis%end_time%value/1000.0_wp
        Write (messages,'(1x,a)') '- the analysis will consider up to '//Trim(Adjustl(word))//' ps of the trajectory' 
        Call info(messages, 1)
      End If
    End If        
    
    If (traj_data%analysis%ignore_initial%fread) Then
      Write (word,'(f8.3)') traj_data%analysis%ignore_initial%value/1000.0_wp
      Write (messages(1),'(1x,a)') '- the initial '//Trim(Adjustl(word))//' ps of the trajectory&
                                   & will be discarded' 
      Call info(messages, 1)
    End If
    
    
    If (ocf_data%invoke%fread .Or. chemocf_data%invoke%fread .Or. &
        msd_data%invoke%fread .Or. transfer_data%invoke%fread) Then
      If (traj_data%analysis%time_interval%fread .Or. &
          (traj_data%analysis%overlap_time%fread  .And. (traj_data%analysis%N_seg /=1))) Then
        Write (messages(1),'(1x,a)') 'Instructions for analysis in time segments (&data_analysis) applied to:' 
        Call info(messages, 1)
        If (ocf_data%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* OCF'
          Call info(messages, 1)
        End If  
        If (transfer_data%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* TCF'
          Write (messages(2),'(5x,a)') '* SPCF'
          Call info(messages, 2)
        End If  
        If (msd_data%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* MSD'
          Call info(messages, 1)
        End If
        If (chemocf_data%invoke%fread) Then
          Write (messages(1),'(5x,a)') '* CHEM_OCF'
          Call info(messages, 1)
        End If
        If (traj_data%analysis%time_interval%fread) Then
          Write (word,'(f8.3)') traj_data%analysis%time_interval%value/1000.0_wp
          If (traj_data%analysis%N_seg /= 1) Then
            Write (messages(1),'(2x,a)') '- the data analysis will be executed using time intervals of '&
                                       &//Trim(Adjustl(word))//' ps' 
          Else
            Write (messages(1),'(2x,a)') '- the data analysis will be executed without the use of time intervals'           
          End If
          Call info(messages, 1)
        End If
        If (traj_data%analysis%overlap_time%fread  .And. (traj_data%analysis%N_seg /=1)) Then
          Write (word,'(f8.3)') traj_data%analysis%overlap_time%value/1000.0_wp
          Write (messages(1),'(2x,a)') '- the starting points of segments for analysis are separated&
                                         & by '//Trim(Adjustl(word))//' ps'
          Call info(messages, 1)
        End If
      End If
    End If
    
    If (ocf_data%invoke%fread) Then
      Call info(' ', 1)
      If (model_data%species_definition%atoms_per_species /= 1) Then
        If (model_data%species_definition%atoms_per_species == 2) Then
          If (Trim(ocf_data%u_definition%type)/='bond_12') Then
            Write (messages(1),'(1x,a)')  '**WARNING: since the monitored species is diatomic, the&
                                          & method to compute the rotating unit vector (u_definition) &
                                          & has been reset to "bond_12". Methods "bond_13", "bond_12-13"&
                                          & "bond_123" and "plane" are meaningless for this case!'
            ocf_data%u_definition%type='bond_12'                             
            Call info(messages, 1)
          End If
        Else
          If (Trim(ocf_data%u_definition%type)/='bond_12-13' .And. &
              Trim(ocf_data%u_definition%type)/='plane' .And. & 
              Trim(ocf_data%u_definition%type)/='bond_123') Then 
            Write (messages(1),'(1x,a)')  '**WARNING: the "'//Trim(ocf_data%u_definition%type)//'"& 
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
                                  & Trim(ocf_data%u_definition%type)
      Write (messages(3),'(1x,a, i2)') '- the correlation terms are obtained using the Legendre polynomial of order ',&
                                  & ocf_data%legendre_order%value
      Call info(messages, 3)
    End If

    If (msd_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&MSD" block will execute a Mean Square&
                                  & Displacement analysis of the species "'&
                                  &//Trim(model_data%species_definition%name%type)//&
                                  & '" (defined in the &monitored_species block) as follows:'
      Write (messages(2),'(1x,a)') '- the values will be computed for the coordinates(s): '//&
                                  & Trim(msd_data%select%type)
      Call info(messages, 2)

      If (msd_data%pbc_xyz%fread) Then
        Write (messages(1),'(1x,a)') '- the "pbc_xyz" directive specifies which coordinate uses (or not) periodic&
                                     & boundary conditions' 
        Call info(messages, 1)                             
      End If
    End If

    If (transfer_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&lifetime" block will compute:'
      Write (messages(2),'(1x,a)') '- the Transfer Correlation Function (TCF) and the Special Pair Correlation&
                                  & Function (SPCF) for the changing chemical species&
                                  & using the method: '//Trim(transfer_data%method%type)
      Write (messages(3),'(1x,a)') '- the residence times for each changing species (file RES_TIMES):'
      If (transfer_data%rattling_wait%fread) Then
        Write (word,'(f8.3)') transfer_data%rattling_wait%value/1000.0_wp
        Write (messages(4),'(3x,a)') 'Rattling times lower than '//Trim(Adjustl(word))//' ps will&
                                    & be discarded for the computation of residence times'  
      Else
        Write (messages(4),'(3x,a)') 'Rattling effects are included in the calculation&
                                     & for the computation of residence times'
      End If
      Call info(messages, 4)
    End If

    If (chemocf_data%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&orientational_chemistry" block will compute&
                                  & OCF for the changing chemical species (CHEM_OCF) using the "'&
                                  &//Trim(chemocf_data%variable%type)//'" as the orientational vector'
      Call info(messages, 1)
    End If
    
    If (traj_data%unchanged%invoke%fread) Then
      Call info(' ', 1)
      Write (messages(1),'(1x,a)') 'The definition of the "&track_unchanged_chemistry" block will print the positions&
                                  & of the selected atomic indexes with unchanged chemistry along the trajectory.'
      Call info(messages, 1)
    End If
    
    If (rdf_data%invoke%fread) Then
      Call info(' ', 1)
      Write (word,'(f10.2)') rdf_data%dr%value
      Write (messages(1),'(1x,a)') 'The definition of the "&RDF" block will compute the Radial&
                                  & Distribution Funtion (RDF) and the Coordination Numbers (CN) using:'
      Write (messages(2),'(1x,a)') '- the tags defined in "tags_species_a"  and "tags_species_b"'
      Write (messages(3),'(1x,a)') '- a discretization of '//Trim(Adjustl(word))//' Angstrom'
      Call info(messages, 3)
    End If

    If (geo_stat_data%coord_distr%invoke%fread) Then
      Call info(' ', 1)
      Write (word,'(f10.2)') geo_stat_data%coord_distr%delta%value
      Write (messages(1),'(1x,a)') 'The definition of the "&coord_distrib block" will compute the distribution&
                                  & of the '//Trim(geo_stat_data%coord_distr%coordinate%type)//'-values for all&
                                  & the "'//Trim(geo_stat_data%coord_distr%species)//'" species in the whole system'
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

    If (geo_stat_data%nndist%invoke%fread) Then
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
      If (ocf_data%invoke%fread .Or. &
          msd_data%invoke%fread .Or. &
          rdf_data%invoke%fread .Or. &
          chemocf_data%invoke%fread .Or. &
          geo_stat_data%nndist%invoke%fread .Or. &
          model_data%species_definition%intra_geom%invoke%fread .Or. &
          model_data%species_definition%inter_geom%invoke%fread) Then
        Write (messages(1),'(1x,a)') 'From the definition of the "&region" block, the computation of'
        Call info(messages, 1)
        If (ocf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- OCF'
          Call info(messages, 1)
        End If  
        If (transfer_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- TCF'
          Write (messages(2),'(3x,a)') '- SPCF'
          Call info(messages, 2)
        End If  
        If (rdf_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- RDF'
          Call info(messages, 1)
        End If  
        If (msd_data%invoke%fread) Then
          Write (messages(1),'(3x,a)') '- MSD'
          Call info(messages, 1)
        End If
        If (chemocf_data%invoke%fread) Then
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
        If (geo_stat_data%nndist%invoke%fread) Then
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
        If (rdf_data%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'IMPORTANT: For the RDF analysis, the definition of the &region block&
                                      & only applies to the species listed in "tags_species_a" (&rdf block)'
          Call info(messages, 1)
        End If  
        If (geo_stat_data%nndist%invoke%fread) Then
          Write (messages(1),'(1x,a)') 'IMPORTANT: For the analysis of the shortest distance distribution,&
                                      & the definition of the &region block only applies to the species& 
                                      & listed in "reference_species" (&selected_nn_distances block)'
          Call info(messages, 1)
        End If  
      End If 
    End If
    
    ! Refresh 
    Call refresh_out(files)

  End Subroutine print_settings_for_trajectory_analysis

  Subroutine trajectory_analysis(files, model_data, traj_data, ocf_data, chemocf_data, msd_data,&
                               & rdf_data, transfer_data, geo_stat_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to analyse the trajectory depending on the options of the
    ! SETTINGS file
    !
    ! author    - i.scivetti Jan 2023
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type),      Intent(InOut) :: files(:)
    Type(model_type),     Intent(In   ) :: model_data
    Type(traj_type),      Intent(InOut) :: traj_data
    Type(ocf_type),       Intent(InOut) :: ocf_data
    Type(chemocf_type),   Intent(InOut) :: chemocf_data
    Type(msd_type),       Intent(InOut) :: msd_data
    Type(rdf_type),       Intent(InOut) :: rdf_data
    Type(transfer_type),  Intent(InOut) :: transfer_data
    Type(geo_stat_type),  Intent(InOut) :: geo_stat_data
  
    Character(Len=256)  :: message
    Logical             :: flag_inter_geo_stat

    traj_data%active_bonds_computed=.False.
    
    ! Compute the total average of monitored species within the system
    If (model_data%config%monitored_species%fread .And. model_data%species_definition%compute_amount%stat) Then
      Call compute_number_monitored_species(traj_data, model_data)
    End If

    If (model_data%change_chemistry%stat) Then 
      Call print_tracking_species(files, traj_data, model_data)
      ! Compute the lifetime and residence time for changing chemical species
      If (transfer_data%invoke%fread) Then
        Call compute_lifetime_related_quantities(files, traj_data, model_data, transfer_data)
      End If
      
      ! Compute OCF for the changing chemical species
      If (chemocf_data%invoke%fread) Then
         If (.Not. traj_data%active_bonds_computed) Then
           Call find_active_bonds(traj_data, model_data)
           traj_data%active_bonds_computed=.True.
         End If
         Call compute_orientational_chemistry(files, traj_data, model_data, chemocf_data)
      End If
    End If
    
    ! Compute OCF for monitored species
    If (ocf_data%invoke%fread) Then
      Call orientational_correlation_function(files, traj_data, ocf_data)
    End If

    ! Compute MSD for monitores species
    If (msd_data%invoke%fread) Then
      Call mean_squared_displacement(files, traj_data, model_data, msd_data)
    End If

    ! Print coordinate distribution for the species under consideration
    If (geo_stat_data%coord_distr%invoke%fread) Then
      Call compute_coordinate_distribution(files, model_data, traj_data, geo_stat_data)
    End If
    
    ! Compute RDF 
    If (rdf_data%invoke%fread) Then
      Call radial_distribution_function(files, traj_data, model_data, rdf_data)
    End If

    ! Print coordinates for selected unchanged species along the MD trajectory
    If (traj_data%unchanged%invoke%fread) Then
      Call print_unchanged_chemistry(files, traj_data)
    End If

    ! Compute the distribution between the shortest distance of a selected pair
    If (geo_stat_data%nndist%invoke%fread) Then
      Call compute_nn_distance_distribution(files, model_data, traj_data, geo_stat_data)
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
      Call find_neighbours_monitored_species(model_data, traj_data, flag_inter_geo_stat)
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
      If (transfer_data%invoke%fread .Or.&
          ocf_data%invoke%fread      .Or.&
          msd_data%invoke%fread) Then   
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
   
End Module analysis
