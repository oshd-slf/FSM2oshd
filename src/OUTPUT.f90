!-----------------------------------------------------------------------
! Write output
!-----------------------------------------------------------------------
subroutine OUTPUT

use DRIVING, only: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour                ! Hour of day

use IOUNITS, only : &
  uout                ! Output file unit number

use STATE_VARIABLES, only: &
  Ds,                &! Snow layer thicknesses (m)
  Nsnow,             &! Number of snow layers
  fsnow,             &! Snow cover fraction 
  rgrn,              &! Snow layer grain radius (m)
  Sice,              &! Ice content of snow layers (kg/m^2)
  Sliq,              &! Liquid content of snow layers (kg/m^2)
  Tsnow,             &! Snow layer temperatures (K)
  Tsoil,             &! Soil layer temperatures (K)
  Tsrf                ! Surface skin temperature (K)

! use DIAGNOSTICS, only : &
!   diags,             &! Cumulated diagnostics
!   Nave,              &! Number of timesteps in average outputs
!   SWint,             &! Cumulated incoming solar radiation (J/m^2)
!   SWout               ! Cumulated reflected solar radiation (J/m^2)

! implicit none

! real*8 :: &
!   alb                 ! Effective albedo

! if (SWint > 0) then
!   alb = SWout / SWint
! else
!   alb = -9
! end if

! ! Averages
! diags(:) = diags(:) / Nave

write(uout,*) year,month,day,hour,sum(Ds),fsnow,sum(Sice+Sliq),Tsrf,Nsnow

! diags(:) = 0
! SWint = 0
! SWout = 0

! 100 format(i4,2(2x,i2),7(2x,f10.3))

end subroutine OUTPUT
