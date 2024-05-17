!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Module for input/output files and related subroutines
!
! Copyright - 2023-2024 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author      -  i.scivetti  Feb 2023
!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module fileset

  Use constants,    Only: code_name,&
                          code_VERSION, &
                          date_RELEASE
  Use numprec,      Only: wi
  use unit_output,  Only: info, &
                          set_output_unit 

  Implicit None
  Private

  ! File data
  Type, Public :: file_type
    Private
    ! Filename
    Character(Len=256), Public :: filename
    ! Fortran unit number, set with newunit=T%unit_no
    Integer(Kind=wi), Public   :: unit_no = -2
  Contains
    Procedure, Public :: init => file_type_init
    Procedure, Public :: Close => close_file
  End Type file_type

  ! SET file
  Integer(Kind=wi), Parameter, Public :: FILE_SET=1
  ! OUT file
  Integer(Kind=wi), Parameter, Public :: FILE_OUT=2 
  ! TRAJECTORY file
  Integer(Kind=wi), Parameter, Public :: FILE_TRAJECTORY=3 
  ! TRACK_CHEMISTRY
  Integer(Kind=wi), Parameter, Public :: FILE_TRACK_CHEMISTRY=4 
  ! TAGGED_TRAJECTORY
  Integer(Kind=wi), Parameter, Public :: FILE_TAGGED_TRAJ=5
  ! Unchanged chemistry
  Integer(Kind=wi), Parameter, Public :: FILE_UNCHANGED_CHEM=6
  ! Radial Distribution Function RDF
  Integer(Kind=wi), Parameter, Public :: FILE_RDF=7
  ! Orientational correlation function OCF
  Integer(Kind=wi), Parameter, Public :: FILE_OCF_ALL=8
  ! Orientational correlation function OCF (averaged)
  Integer(Kind=wi), Parameter, Public :: FILE_OCF_AVG=9
  ! Mean Squared Displacement MSD
  Integer(Kind=wi), Parameter, Public :: FILE_MSD_ALL=10 
  ! Mean Squared Displacement MSD (averaged)
  Integer(Kind=wi), Parameter, Public :: FILE_MSD_AVG=11
  ! Residence Time Correlation Function
  Integer(Kind=wi), Parameter, Public :: FILE_TCF_ALL=12 
  ! Residence Time Correlation Function (average)
  Integer(Kind=wi), Parameter, Public :: FILE_TCF_AVG=13 
  ! Residence times
  Integer(Kind=wi), Parameter, Public :: FILE_RES_TIMES=14 
  ! Coordinate distribution
  Integer(Kind=wi), Parameter, Public :: FILE_COORD_DISTRIB=15 
  ! Intramol distances
  Integer(Kind=wi), Parameter, Public :: FILE_INTRAMOL_DISTANCES=16 
  ! Intramol angles
  Integer(Kind=wi), Parameter, Public :: FILE_INTRAMOL_ANGLES=17 
  ! Intermol distances to the first NN
  Integer(Kind=wi), Parameter, Public :: FILE_INTERMOL_DISTANCES_NN1=18 
  ! Intermol distances to the second NN
  Integer(Kind=wi), Parameter, Public :: FILE_INTERMOL_DISTANCES_NN2=19 
  ! Intermol angles
  Integer(Kind=wi), Parameter, Public :: FILE_INTERMOL_ANGLES=20 
  ! Special Pair Correlation Function
  Integer(Kind=wi), Parameter, Public :: FILE_SPCF_ALL=21 
  ! Special Pair Correlation Function (averaged)
  Integer(Kind=wi), Parameter, Public :: FILE_SPCF_AVG=22
  ! Orientational correlation function for shortest bond TBOCF
  Integer(Kind=wi), Parameter, Public :: FILE_CHEM_OCF_ALL=23
  ! Orientational correlation function OCF (averaged)
  Integer(Kind=wi), Parameter, Public :: FILE_CHEM_OCF_AVG=24
  ! Shortest distance for selected pair of species
  Integer(Kind=wi), Parameter, Public :: FILE_SELECTED_NN_DISTANCES=25
  
  ! Size of filename array
  Integer(Kind=wi), Parameter, Public :: NUM_FILES = 25

  Public :: set_system_files, print_header_out, wrapping_up, refresh_out

Contains

  Subroutine refresh_out(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to refresh the output
    !
    ! author    - i.scivetti Oct 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES) 

    Call files(FILE_OUT)%close ()
    Open (Newunit=files(FILE_OUT)%unit_no, File=files(FILE_OUT)%filename, Position='Append')
    Call set_output_unit(files(FILE_OUT)%unit_no)

  End Subroutine refresh_out

  Subroutine file_type_init(T, filename)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to initialise files
    !
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type)                :: T
    Character(Len=*), Intent(In   ) :: filename

    T%filename = Trim(filename)
  End Subroutine file_type_init


  Subroutine set_names_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set default names for files
    !
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=256), Dimension(NUM_FILES)   :: set_names
    Integer(Kind=wi)                           :: file_no

    ! Default file names array
    ! Populate default names array
    set_names(FILE_SET)             = "SETTINGS"
    set_names(FILE_OUT)             = "OUTPUT"
    set_names(FILE_TRAJECTORY)      = "TRAJECTORY"
    set_names(FILE_TRACK_CHEMISTRY) = "TRACK_CHEMISTRY"
    set_names(FILE_TAGGED_TRAJ)     = "TAGGED_TRAJECTORY"
    set_names(FILE_UNCHANGED_CHEM)  = "UNCHANGED_CHEMISTRY"
    set_names(FILE_RDF)             = "RDF"
    set_names(FILE_OCF_ALL)         = "OCF_ALL"
    set_names(FILE_OCF_AVG)         = "OCF_AVG"
    set_names(FILE_MSD_ALL)         = "MSD_ALL"
    set_names(FILE_MSD_AVG)         = "MSD_AVG"
    set_names(FILE_TCF_ALL)         = "TCF_ALL"
    set_names(FILE_TCF_AVG)         = "TCF_AVG"
    set_names(FILE_RES_TIMES)       = "RES_TIMES"
    set_names(FILE_COORD_DISTRIB)   = "COORD_DISTRIBUTION"
    set_names(FILE_SPCF_ALL)        = "SPCF_ALL"
    set_names(FILE_SPCF_AVG)        = "SPCF_AVG"
    set_names(FILE_CHEM_OCF_ALL)    = "CHEM_OCF_ALL"
    set_names(FILE_CHEM_OCF_AVG)    = "CHEM_OCF_AVG"
    set_names(FILE_INTRAMOL_DISTANCES) = "INTRAMOL_DISTANCES"
    set_names(FILE_INTRAMOL_ANGLES)    = "INTRAMOL_ANGLES"
    set_names(FILE_INTERMOL_ANGLES)    = "INTERMOL_ANGLES_NN"
    set_names(FILE_INTERMOL_DISTANCES_NN1) = "INTERMOL_DISTANCES_NN1"
    set_names(FILE_INTERMOL_DISTANCES_NN2) = "INTERMOL_DISTANCES_NN2"
    set_names(FILE_SELECTED_NN_DISTANCES)  = "SELECTED_NN_DISTANCES"
    
    ! Set default filenames
    Do file_no = 1, NUM_FILES
      Call files(file_no)%init(set_names(file_no))
    End Do

  End Subroutine set_names_files


  Subroutine close_file(T)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to close files
    !
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    Class(file_type) :: T

    Logical :: is_open

    Inquire (T%unit_no, opened=is_open)
    If (is_open) Then
      Close (T%unit_no)
      T%unit_no = -2
    End If

  End Subroutine close_file

  Subroutine set_system_files(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to open OUTPUT file 
    ! 
    ! author    - i.scivetti April 2020
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Call set_names_files(files)   
    Open (Newunit=files(FILE_OUT)%unit_no, File=files(FILE_OUT)%filename, Status='replace')
    Call set_output_unit(files(FILE_OUT)%unit_no)

  End Subroutine set_system_files   

  Subroutine print_header_out(files)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Subroutine to print the header to OUTPUT file 
    !  
    ! author        - i.scivetti July 2020
    ! contribution  - i.scivetti Oct  2021
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)

    Character(Len=*), Parameter :: fmt1 = '(a)'
    Character(Len=*), Parameter :: fmt2 = '(3a)'
    Character(Len=*), Parameter :: fmt3 = '(4a)'
    Character(Len=128)          :: header(14)

    Write (header(1), fmt1)   Repeat("#", 74)
    Write (header(2), fmt2)  "#                      WELCOME TO ", Trim(code_name),  Repeat(" ", 25)//"#"
    Write (header(3), fmt1)  "#  An ALC software for MD analysis of reactive and non-reactive systems  #"
    Write (header(4), fmt1)  "#                                                                        #"
    Write (header(5), fmt3)  "#  version:  ", Trim(code_VERSION), Repeat(' ',57),                     "#"
    Write (header(6), fmt3)  "#  release:  ", Trim(date_RELEASE), Repeat(' ',52),                     "#"
    Write (header(7), fmt1)  "#                                                                        #"
    Write (header(8), fmt1)  "#  Copyright:  2024  Ada Lovelace Centre (ALC)                           #"
    Write (header(9), fmt1)  "#              Scientific Computing Department (SCD)                     #"
    Write (header(10), fmt1) "#              Science and Technology Facilities Councils (STFC)         #"
    Write (header(11), fmt1) "#                                                                        #"
    Write (header(12), fmt1) "#  Author:            Ivan Scivetti (SCD/STFC)                           #"
    Write (header(13), fmt1) "#  Project support:   Gilberto Teobaldi (SCD/STFC)                       #"
    Write (header(14), fmt1)  Repeat("#", 74)
    Call info(header, 14)
    
    ! Refresh OUT_EQCM
    Call refresh_out(files)

  End Subroutine print_header_out

  Subroutine wrapping_up(files)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine to print final remarks to OUT_EQCM file 
  ! and close the file 
  !  
  ! author    - i.scivetti July 2020
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Type(file_type) :: files(NUM_FILES)
 
    Character(Len=*), Parameter :: fmt1 = '(1x,a)'
    Character(Len=*), Parameter :: fmt2 = '(1x,3a)'
    Character(Len=128)          :: appex(7)
     
    Write (appex(1), fmt1)   Repeat(" ", 1)
    Write (appex(2), fmt1)   Repeat("#", 38)
    Write (appex(3), fmt1)  "#                                    #" 
    Write (appex(4), fmt1)  "#  Job has finished successfully     #"
    Write (appex(5), fmt2)  "#  Thanks for using ", Trim(code_name), "!  #"
    Write (appex(6), fmt1)  "#                                    #" 
    Write (appex(7), fmt1)   Repeat("#", 38)
    Call info(appex, 7)

    Close(files(FILE_OUT)%unit_no)    

  End Subroutine wrapping_up  

End Module fileset
