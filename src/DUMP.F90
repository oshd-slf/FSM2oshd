!-----------------------------------------------------------------------
! Write out state variables at end of run
!-----------------------------------------------------------------------
subroutine DUMP

use MODCONF, only: CANMOD,SNFRAC

use MODE_WRITE, only: WRITE_2D

use IOUNITS, only : &
  udmp                ! Dump file unit number

use FILES, only : dump_file

use STATE_VARIABLES, only : &
  albs,              &! Snow albedo
  Ds,                &! Snow layer thicknesses (m)
  fsnow,             &! Snow cover fraction
  Nsnow,             &! Number of snow layers 
  Qcan,              &! Canopy air space humidity
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  snowdepthmin,      &! min Snow at time step of swemin (m)
  snowdepthmax,      &! max Snow at time step of swemax (m)
  snowdepthhist,     &! history of Snow depth during last 14 days (m)  
  swemin,            &! min swe during season (mm)
  swemax,            &! max swe during season (mm)  
  swehist,           &! history of Snow depth during last 14 days (kg/m^2)
  Sveg,              &! Canopy snow mass (kg/m^2)
  Tcan,              &! Canopy air space temperature (K)
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  theta,             &! Volumetric moisture content of soil layers
  Tsrf,              &! Surface skin temperature (K)
  Tveg                ! Vegetation temperature (K)

implicit none

if (dump_file /= 'none') then
  open(udmp, file = dump_file)  
  write(udmp,*) albs(:,:)
  write(udmp,*) Ds(:,:,:)
  write(udmp,*) Nsnow(:,:)
  write(udmp,*) Qcan(:,:)
  write(udmp,*) Sice(:,:,:)
  write(udmp,*) Sliq(:,:,:)
  write(udmp,*) Sveg(:,:)
  write(udmp,*) Tcan(:,:)
  write(udmp,*) theta(:,:,:)
  write(udmp,*) Tsnow(:,:,:)
  write(udmp,*) Tsoil(:,:,:)
  write(udmp,*) Tsrf(:,:)
  write(udmp,*) fsnow(:,:)
  write(udmp,*) Tveg(:,:)
  write(udmp,*) snowdepthmin(:,:)
  write(udmp,*) snowdepthmax(:,:)
  write(udmp,*) snowdepthhist(:,:,:)
  write(udmp,*) swemin(:,:)
  write(udmp,*) swemax(:,:)
  write(udmp,*) swehist(:,:,:)
  close(udmp)
end if
    
  end subroutine DUMP

