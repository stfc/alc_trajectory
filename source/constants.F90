!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing constants and parameters for computation
!
! Copyright - 2023 Ada Lovelace Centre (ALC)
!             Scientific Computing Department (SCD)
!             The Science and Technology Facilities Council (STFC)
!
! Author    - i.scivetti Feb 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Module constants

  Use numprec, Only: wi, &
                     wp

  Implicit None

  ! Code reference 
  Character(Len=16), Parameter, Public  :: code_name    = "ALC_TRAJECTORY" 
  Character(Len=16), Parameter, Public  :: code_VERSION = "1.2"
  Character(Len=16), Parameter, Public  :: date_RELEASE = "Dec 2023"

  ! FIXED PARAMETERS
  Real(Kind=wp), Parameter, Public  :: pi    = 3.14159265358979312e0_wp 
  Real(Kind=wp), Parameter, Public  :: twopi = 6.28318530717958623e0_wp 
  Real(Kind=wp), Parameter, Public  :: Farad = 96485332.900_wp
  Real(Kind=wp), Parameter, Public  :: Bohr_to_A = 0.529177249_wp
  Real(Kind=wp), Parameter, Public  :: Rads_to_degrees = 57.2957795130823_wp
  Real(Kind=wp), Parameter, Public  :: avogadro  = 6.02214076E+23
  Real(Kind=wp), Parameter, Public  :: K_to_eV   = 8.61732814974056E-05_wp
  Real(Kind=wp), Parameter, Public  :: Ha_to_eV  =27.211386245988_wp
  
  ! UNITS CONVERSION
  Real(Kind=wp), Parameter, Public  :: g_to_ng  = 1.0E+9 
  Real(Kind=wp), Parameter, Public  :: cm_to_Ang= 1.0E+8

  ! Number of Periodic Table Elements
  Integer(Kind=wi), Parameter, Public   :: NPTE = 118

  Character(Len=2), Dimension(NPTE), Parameter, Public :: chemsymbol=&
              (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc',  &
               'Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo',  &
               'Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',  &
               'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po',  &
               'At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',  &
               'Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/) 

  ! Maximum number of components
  Integer(Kind=wi), Parameter, Public   :: max_components= 50
  ! Maximum number of atoms per species
  Integer(Kind=wi), Parameter, Public   :: max_at_species= 50
  ! Maximum number of unchanged species to track along the trajectory
  Integer(Kind=wi), Parameter, Public   ::  max_unchanged_atoms = 10  

End Module constants
