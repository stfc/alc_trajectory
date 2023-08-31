!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module to define real and integer precisions
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author: i.scivetti Feb 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module numprec

  Use, Intrinsic :: iso_fortran_env, Only: int32, &
                                           int64, & 
                                           real64
  Implicit None

  Private
  ! Double real
  Integer, Parameter, Public :: dp = real64
  ! Working real
  Integer, Parameter, Public :: wp = dp
  ! 32-bit integer 
  Integer, Parameter, Public :: ni = int32
  ! Working integer
  Integer, Parameter, Public :: wi = ni
  ! Long integer
  Integer, Parameter, Public :: li = int64

End Module numprec
