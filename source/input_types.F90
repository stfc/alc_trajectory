!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Definition of common types   
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author:     i.scivetti Feb 2023 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module input_types
  Use numprec,    Only: wi,&
                        wp

  Implicit None
  Private

  Type, Public :: in_string
    Character(Len=256) :: type=repeat(' ',256) 
    Logical            :: fread= .False.
    Logical            :: fail = .False.
    Logical            :: warn = .False. 
  End Type

  Type, Public :: in_integer
    Integer(Kind=wi)   :: value
    Character(Len=256) :: tag
    Logical            :: fread= .False.
    Logical            :: fail = .False.
    Logical            :: warn = .False. 
  End Type

  Type, Public :: in_integer_array
    Integer(Kind=wi), Allocatable :: value(:)
    Character(Len=256) :: tag
    Logical            :: fread= .False.
    Logical            :: fail = .False.
    Logical            :: warn = .False. 
  End Type

  Type, Public :: in_logic
    Logical           :: stat = .False. 
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_scalar
    Real(Kind=wp)     :: value
    Logical           :: fread= .False.
    Logical           :: fail = .False.
    Logical           :: warn = .False. 
  End Type

  Type, Public :: in_param
    Character(Len=256) :: tag
    Character(Len=16)  :: units
    Real(Kind=wp)      :: value
    Real(Kind=wp)      :: convert
    Logical            :: fread= .False.
    Logical            :: fail = .False.
    Logical            :: warn = .False. 
  End Type

End Module input_types
