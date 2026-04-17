!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Welcome to ALC_TRAJECTORY: a ALC software to analyse MD trajectories
! both for reactive and non-reactive systems. The main purpose of this
! code is to offer the posibility to compute simultaneously:
!
! * Radial Distribution Functions
! * Transfer Correlation Functions (only for reactive systems)
! * Mean Square Dislpacements
! * Orientational Correlation Functions
! * Special Pair Orientational and Transfer Correlation Functions
!
! Please refer to file "use_code.md" for a detailed explanation of the capabilities
! Example cases can be found in the examples folder
!
! This code is available under the BSD 3-Clause License.
!
! Copyright   2023-2026 Ada Lovelace Centre   (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)  
!               
! Author:            Ivan Scivetti (i.scivetti)
! Project support:   Gilberto Teobaldi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program alc_trajectory

  Use analysis,          Only: print_settings_for_trajectory_analysis,&
                               trajectory_analysis  
                         
  Use atomic_model,      Only: model_type, &
                               atomistic_model, &
                               read_model, &
                               print_model_settings
                         
  Use fileset,           Only: file_type, &
                               NUM_FILES, &
                               print_header_out, &
                               set_system_files, &
                               wrapping_up
                               
  Use geo_stat,          Only: geo_stat_type                              
                           
  Use msd,               Only: msd_type                         
                         
  Use numprec,           Only: wi,& 
                               wp
                               
  Use ocf,               Only: ocf_type, &
                               chemocf_type
                               
  Use rdf,               Only: rdf_type                         
                               
  Use settings,          Only: read_settings, &
                               check_settings
                         
  Use trajectory,        Only: traj_type, &
                               extract_trajectory,&
                               trajectory_setup   

  Use transfer_species,  Only: transfer_type                         
                               
  Use unit_output,       Only: info
 

Implicit None
! Definition of variables
  Type(file_type)      :: files(NUM_FILES)
  Type(model_type)     :: model_data
  Type(traj_type)      :: traj_data
  Type(ocf_type)       :: ocf_data
  Type(chemocf_type)   :: chemocf_data
  Type(msd_type)       :: msd_data
  Type(rdf_type)       :: rdf_data
  Type(transfer_type)  :: transfer_data
  Type(geo_stat_type)  :: geo_stat_data

  !Time related variables
  Integer(kind=wi)   :: start,finish,rate

  ! Array to print information
  Character(Len=256) :: message

  ! Start of the code 
  !!!!!!!!!!!!!!!!!!!
  ! Record initial time
  Call system_clock(count_rate=rate)
  Call system_clock(start)
  ! Initialise settings for input/output files
  Call set_system_files(files)
  ! Print header of OUT
  Call print_header_out(files) 
  ! Read settings from SET
  Call read_settings(files, model_data, traj_data, ocf_data, chemocf_data, msd_data,&
                   & rdf_data, transfer_data, geo_stat_data)
  ! Check the specification of settings in SET
  Call check_settings(files, model_data, traj_data, ocf_data, chemocf_data, msd_data,&
                   & rdf_data, transfer_data, geo_stat_data)
  ! Print model related settings
  Call print_model_settings(files, model_data)  
  ! Prepare trajectory
  Call trajectory_setup(files, model_data, traj_data)
  ! Print trajectory settings
  Call print_settings_for_trajectory_analysis(files, model_data, traj_data, ocf_data, chemocf_data,&
                                            & msd_data, rdf_data, transfer_data, geo_stat_data)
  ! Read and define trajectory
  Call extract_trajectory(files, model_data, traj_data)
  ! Perform the requested analysis
  Call trajectory_analysis(files, model_data, traj_data, ocf_data, chemocf_data, msd_data,&
                                           & rdf_data, transfer_data, geo_stat_data)

  ! Record final time
  Call system_clock(finish)

  ! Print execution time
  Call info(' ', 1)
  Call info(' ==========================================', 1)
  Write (message, '(1x,a,f9.3,a)') 'Total execution time = ',  Real(finish-start,Kind=wp)/rate,  ' seconds.' 
  Call info(message, 1)
  Call info(' ==========================================', 1)

  ! Print appendix to OUTPUT file
  Call wrapping_up(files)

End Program alc_trajectory
