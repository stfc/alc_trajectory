!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! Welcome to ALC_TRAJECTORY: a ALC software to analyse MD trajectories
! both for reactive and non-reactive systems. The main purpose of this
! code is to offer the posibility to compute simultaneously:
!
! * Radial Distribution Functions
! * Transfer Correlation Functions (only for reactive systems)
! * Mean Square Dislpacements
! * Orientational Correlation Functions
!
! This code is available under the BSD 3-Clause License.
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)  
!               
! Author:            Ivan Scivetti (i.scivetti)
! Project support:   Gilberto Teobaldi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program alc_trajectory

  Use atomic_model,  Only : model_type, &
                            atomistic_model, &
                            read_model

  Use fileset,       Only: file_type, &
                           NUM_FILES, &
                           print_header_out, &
                           set_system_files, &
                           wrapping_up

  Use numprec,       Only: wi,& 
                           wp
                           
  Use settings,      Only: read_settings, &
                           check_settings

  Use trajectory,    Only: traj_type, &
                           research_trajectory,&
                           trajectory_analysis

  Use unit_output,   Only: info
  
  

Implicit None
! Definition of variables
  Type(file_type)      :: files(NUM_FILES)
  Type(model_type)     :: model_data
  Type(traj_type)      :: traj_data

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
  Call read_settings(files, model_data, traj_data)
  ! Check the specification of settings in SET
  Call check_settings(files, model_data, traj_data)
  ! Read and define trajectory
  Call research_trajectory(files, model_data, traj_data)
  ! To be implemented: Perform the requested analysis
  Call trajectory_analysis(files, model_data, traj_data)

  ! Record final time
  Call system_clock(finish)

  ! Print execution time
  Call info(' ', 1)
  Call info(' ==========================================', 1)
  Write (message, '(1x,a,f9.3,a)') 'Total execution time = ',  Real(finish-start,Kind=wp)/rate,  ' seconds.' 
  Call info(message, 1)
  Call info(' ==========================================', 1)

  ! Print appendix to OUT_EQCM file
  Call wrapping_up(files)

End Program alc_trajectory
