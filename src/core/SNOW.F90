!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
subroutine SNOW(Esrf,G,ksnow,ksoil,Melt,unload,Gsoil,Roff,meltflux_out,Sbsrf)

use MODCONF, only: HYDROL,DENSTY,OSHDTN,HN_ON,SNFRAC

use MODTILE, only: tthresh
 
use CONSTANTS, only: &
  grav,              &! Acceleration due to gravity (m/s^2)
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  rho_ice,           &! Density of ice (kg/m^3)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use DRIVING, only: &
  dt,                &! Timestep (s)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)

use GRID, only: &
  Dzsnow,            &! Minimum snow layer thicknesses (m)
  Dzsoil,            &! Soil layer thicknesses (m)
  Nsmax,             &! Maximum number of snow layers
  Nsoil,             &! Number of soil layers
  Nx,Ny               ! Grid dimensions

use PARAMETERS, only: &
  a_eta,             &! Temperature factor for Crocus B92 compaction (K^-1)
  b_eta,             &! First density factor for Crocus B92 compaction (m^3/kg)
  c_eta,             &! Second density factor for Crocus B92 compaction (kg/m^3)
  eta0,              &! Reference snow viscosity (Pa s)
  eta1,              &! Reference snow viscosity for Crocus B92 compaction (Pa s)
  rgr0,              &! Fresh snow grain radius (m)
  rho0,              &! Fixed snow density (kg/m^3)
  rhob,              &! Temperature factor in fresh snow density (kg/m^3/K)
  rhoc,              &! Wind factor in fresh snow density (kg s^0.5/m^3.5)
  rhof,              &! Fresh snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  snda,              &! Thermal metamorphism parameter (1/s)
  trho,              &! Snow compaction timescale (s)
  Wirr                ! Irreducible liquid water content of snow

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

use LANDUSE, only: &
  dem,               &! Terrain elevation (m)
  tilefrac            ! Grid cell tile fraction

implicit none

real, intent(in) :: &
  Esrf(Nx,Ny),       &! Moisture flux from the surface (kg/m^2/s)
  G(Nx,Ny),          &! Heat flux into surface (W/m^2)
  ksnow(Nsmax,Nx,Ny),&! Thermal conductivity of snow (W/m/K)
  ksoil(Nsoil,Nx,Ny),&! Thermal conductivity of soil (W/m/K)
  Melt(Nx,Ny),       &! Surface melt rate (kg/m^2/s)
  unload(Nx,Ny)       ! Snow mass unloaded from canopy (kg/m^2)

real, intent(out) :: &
  Gsoil(Nx,Ny),      &! Heat flux into soil (W/m^2)
  Roff(Nx,Ny),       &! Total runoff (kg/m^2)
  meltflux_out(Nx,Ny),  &! Runoff from snowmelt at base of snow (kg/m^2)
  Sbsrf(Nx,Ny)        ! Sublimation from the snow surface (kg/m^2)

integer :: &
  i,j,               &! Point counters
  k,                 &! Level counter
  knew,              &! New snow layer pointer
  kold,              &! Old snow layer pointer
  Nold                ! Previous number of snow layers

real :: &
  coldcont,          &! Layer cold content (J/m^2)
  dnew,              &! New snow layer thickness (m)
  dSice,             &! Change in layer ice content (kg/m^2)
  Esnow,             &! Snow sublimation rate (kg/m^2/s)
  eta,               &! Snow viscosity (Pa s)
  f1,                &! Wet snow settlement factor
  f2,                &! Grain size snow settlement factor
  fold,              &! Previous snow cover fraction
  ggr,               &! Grain area growth rate (m^2/s)
  mass,              &! Mass of overlying snow (kg/m^2)
  phi,               &! Porosity
  rhonew,            &! Density of new snow (kg/m^3)
  rhos,              &! Density of snow layer (kg/m^3)
  SliqCap,           &! Liquid water holding capacity of snow
  SliqMax,           &! Maximum liquid content for layer (kg/m^2)
  snowdepth,         &! Snow depth (m)
  t_decompaction,    &! Decompaction time for fresh snow density (h)
  Tsnow0,            &! Falling snow temperature (K)
  wt                  ! Layer weighting

real :: &
  a(Nsmax),          &! Below-diagonal matrix elements
  b(Nsmax),          &! Diagonal matrix elements
  c(Nsmax),          &! Above-diagonal matrix elements
  csnow(Nsmax),      &! Areal heat capacity of snow (J/K/m^2)
  dTs(Nsmax),        &! Temperature increments (K)
  D(Nsmax),          &! Layer thickness before adjustment (m)
  E(Nsmax),          &! Energy contents before adjustment (J/m^2)
  Gs(Nsmax),         &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nsmax),        &! Matrix equation rhs
  R(Nsmax),          &! Snow grain radii before adjustment (kg/m^2)
  S(Nsmax),          &! Ice contents before adjustment (kg/m^2)
  U(Nsmax),          &! Layer internal energy contents (J/m^2)
  W(Nsmax),          &! Liquid contents before adjustment (kg/m^2)
  SWEtmp              ! Temporary snow water equivalent for snow covered fraction calculation (kg/m^2)
  
real :: &
  fsnow_thres(Nx,Ny),&! Thresholded snow cover fraction
  Roff_snow(Nx,Ny),  &! Runoff at base of snow (kg/m^2)
  Roff_bare(Nx,Ny)    ! Bare soil runoff (kg/m^2)

Gsoil(:,:) = G(:,:)
Roff(:,:) = 0
meltflux_out(:,:) = 0
Roff_bare(:,:) = Rf(:,:) * dt * (1 - fsnow(:,:))
Roff_snow(:,:) = Rf(:,:) * dt * fsnow(:,:)

! Points with existing snowpack
do j = 1, Ny
do i = 1, Nx

  if (tilefrac(i,j) < tthresh) goto 1 ! exclude points outside tile of interest

  Sbsrf(i,j) = 0

  if (fsnow(i,j) > epsilon(fsnow)) then ! This condition should be equivalent to Nsnow(i,j) > 0

    ! Except for point case, apply a minimum threshold of 0.1 to fsnow
    ! to avoid 'long tails' in SWE due to slowing down depletion rates
    if (SNFRAC == 3) then
      fsnow_thres(i,j) = fsnow(i,j)
    else
      fsnow_thres(i,j) = max(fsnow(i,j),0.1)
    end if

    ! Heat conduction
    do k = 1, Nsnow(i,j)
      csnow(k) = (Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat) / fsnow(i,j)
    end do
    if (Nsnow(i,j) == 1) then
      Gs(1) = 2 / (Ds(1,i,j)/ksnow(1,i,j) + Dzsoil(1)/ksoil(1,i,j))
      dTs(1) = (G(i,j) + Gs(1)*(Tsoil(1,i,j) - Tsnow(1,i,j)))*dt /  &
               (csnow(1) + Gs(1)*dt)
    else
      do k = 1, Nsnow(i,j) - 1
        Gs(k) = 2 / (Ds(k,i,j)/ksnow(k,i,j) + Ds(k+1,i,j)/ksnow(k+1,i,j))
      end do
      a(1) = 0
      b(1) = csnow(1) + Gs(1)*dt
      c(1) = - Gs(1)*dt
      rhs(1) = (G(i,j) - Gs(1)*(Tsnow(1,i,j) - Tsnow(2,i,j)))*dt
      do k = 2, Nsnow(i,j) - 1
        a(k) = c(k-1)
        b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
        c(k) = - Gs(k)*dt
        rhs(k) = Gs(k-1)*(Tsnow(k-1,i,j) - Tsnow(k,i,j))*dt  &
                 + Gs(k)*(Tsnow(k+1,i,j) - Tsnow(k,i,j))*dt 
      end do
      k = Nsnow(i,j)
      Gs(k) = 2 / (Ds(k,i,j)/ksnow(k,i,j) + Dzsoil(1)/ksoil(1,i,j))
      a(k) = c(k-1)
      b(k) = csnow(k) + (Gs(k-1) + Gs(k))*dt
      c(k) = 0
      rhs(k) = Gs(k-1)*(Tsnow(k-1,i,j) - Tsnow(k,i,j))*dt  &
               + Gs(k)*(Tsoil(1,i,j) - Tsnow(k,i,j))*dt
      call TRIDIAG(Nsnow(i,j),Nsmax,a,b,c,rhs,dTs)
    end if 
    do k = 1, Nsnow(i,j)
      Tsnow(k,i,j) = Tsnow(k,i,j) + dTs(k)
      if (HN_ON) then
        Tsnow(k,i,j) = max(Tsnow(k,i,j), (Tm-40))
      endif
    end do
    k = Nsnow(i,j)
    Gsoil(i,j) = Gs(k)*(Tsnow(k,i,j) - Tsoil(1,i,j))

    ! Convert melting ice to liquid water
    dSice = Melt(i,j)*fsnow_thres(i,j)*dt

    meltflux_out(i,j) = dSice
    
    do k = 1, Nsnow(i,j)
      coldcont = csnow(k)*(Tm - Tsnow(k,i,j))
      if (coldcont < 0) then
        dSice = dSice - fsnow(i,j)*coldcont/Lf
        Tsnow(k,i,j) = Tm
      end if
      if (dSice > epsilon(dSice)) then
        if (dSice > Sice(k,i,j)) then  ! Layer melts completely
          dSice = dSice - Sice(k,i,j)
          Ds(k,i,j) = 0
          Sliq(k,i,j) = Sliq(k,i,j) + Sice(k,i,j)
          Sice(k,i,j) = 0
        else                       ! Layer melts partially
          Ds(k,i,j) = (1 - dSice/Sice(k,i,j))*Ds(k,i,j)
          Sice(k,i,j) = Sice(k,i,j) - dSice
          Sliq(k,i,j) = Sliq(k,i,j) + dSice
          dSice = 0                ! Melt exhausted
        end if
      end if
    end do

    ! Remove snow by sublimation 
    dSice = max(Esrf(i,j)*fsnow_thres(i,j), 0.)*dt
    if (dSice > epsilon(dSice)) then
      do k = 1, Nsnow(i,j)
        if (dSice > Sice(k,i,j)) then  ! Layer sublimates completely
          dSice = dSice - Sice(k,i,j)
          Ds(k,i,j) = 0
          Sbsrf(i,j) = Sbsrf(i,j) + Sice(k,i,j)
          Sice(k,i,j) = 0
        else                       ! Layer sublimates partially
          Ds(k,i,j) = (1 - dSice/Sice(k,i,j))*Ds(k,i,j)
          Sice(k,i,j) = Sice(k,i,j) - dSice
          Sbsrf(i,j) = Sbsrf(i,j) + dSice
          dSice = 0                ! Sublimation exhausted
        end if
      end do
    end if

    ! Snow hydraulics
    ! First, unloading snow is added to liquid water if Ta above freezing point
    if (Ta(i,j) >= Tm) then
      Roff_bare(i,j) = Roff_bare(i,j) + unload(i,j) * (1 - fsnow(i,j)) ! Bare soil fraction
      Roff_snow(i,j) = Roff_snow(i,j) + unload(i,j) * fsnow(i,j) ! Snow covered ground fraction
    end if
    if (HYDROL == 0) then
      ! Free-draining snow 
      meltflux_out(i,j) = 0
      do k = 1, Nsnow(i,j)
        Roff_snow(i,j) = Roff_snow(i,j) + Sliq(k,i,j)
        meltflux_out(i,j) =  meltflux_out(i,j) + Sliq(k,i,j)
        Sliq(k,i,j) = 0
      end do
    elseif (HYDROL == 1) then
      ! Bucket storage 
      do k = 1, Nsnow(i,j)
        phi = 0
        if (Ds(k,i,j) > epsilon(Ds)) phi = 1 - Sice(k,i,j)/(rho_ice*Ds(k,i,j)*fsnow(i,j))
        SliqMax = fsnow(i,j)*rho_wat*Ds(k,i,j)*phi*Wirr
        Sliq(k,i,j) = Sliq(k,i,j) + Roff_snow(i,j)
        Roff_snow(i,j) = 0
        if (Sliq(k,i,j) > SliqMax) then       ! Liquid capacity exceeded
          Roff_snow(i,j) = Sliq(k,i,j) - SliqMax   ! so drainage to next layer
          Sliq(k,i,j) = SliqMax
        end if
        coldcont = csnow(k)*(Tm - Tsnow(k,i,j))
        if (coldcont > 0) then       ! Liquid can freeze
          dSice = min(Sliq(k,i,j), fsnow(i,j)*coldcont/Lf) 
          Sliq(k,i,j) = Sliq(k,i,j) - dSice
          Sice(k,i,j) = Sice(k,i,j) + dSice
          meltflux_out(i,j) = meltflux_out(i,j) - dSice     
          Tsnow(k,i,j) = Tsnow(k,i,j) + Lf*dSice/csnow(k)/fsnow(i,j)
        end if
      end do
      
      if (meltflux_out(i,j) < 0) then
        meltflux_out(i,j) = 0
      end if
      
    else  ! HYDROL == 2
      ! Density-dependent bucket storage
      do k = 1, Nsnow(i,j)
        if (Ds(k,i,j) > epsilon(Ds)) then
          rhos = Sice(k,i,j) / Ds(k,i,j) / fsnow(i,j)
          SliqCap = 0.03 + 0.07*(1 - rhos/200)
          SliqCap = max(SliqCap, 0.03)
        end if
        SliqMax = SliqCap*Sice(k,i,j)
        Sliq(k,i,j) = Sliq(k,i,j) + Roff_snow(i,j)
        Roff_snow(i,j) = 0
        if (Sliq(k,i,j) > SliqMax) then       ! Liquid capacity exceeded
          Roff_snow(i,j) = Sliq(k,i,j) - SliqMax   ! so drainage to next layer
          Sliq(k,i,j) = SliqMax
        end if
        coldcont = csnow(k)*(Tm - Tsnow(k,i,j))
        if (coldcont > epsilon(coldcont)) then       ! Liquid can freeze
          dSice = min(Sliq(k,i,j), fsnow(i,j)*coldcont/Lf) 
          Sliq(k,i,j) = Sliq(k,i,j) - dSice
          Sice(k,i,j) = Sice(k,i,j) + dSice
          ! to account for refreezing of melt
          meltflux_out(i,j) = meltflux_out(i,j) - dSice          
          Tsnow(k,i,j) = Tsnow(k,i,j) + Lf*dSice/csnow(k)/fsnow(i,j)
        end if
      end do
      
      if (meltflux_out(i,j) < 0) then
        meltflux_out(i,j) = 0
      end if
      
    endif

    ! Snow compaction
    if (DENSTY == 0) then
      ! Fixed snow density
      do k = 1, Nsnow(i,j)
        Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rho0 / fsnow(i,j)
      end do
    elseif (DENSTY == 1) then
      ! Snow compaction with age
      do k = 1, Nsnow(i,j)
        if (Ds(k,i,j) > epsilon(Ds)) then
          rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)
          if (Tsnow(k,i,j) >= Tm) then
              if (rhos < rmlt) rhos = rmlt + (rhos - rmlt)*exp(-dt/trho)
          else
              if (rhos < rcld) rhos = rcld + (rhos - rcld)*exp(-dt/trho)
          end if
          Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rhos / fsnow(i,j)
        end if
      end do
    else if (DENSTY == 2) then
      ! Snow compaction by overburden
      mass = 0
      do k = 1, Nsnow(i,j)
        mass = mass + 0.5*(Sice(k,i,j) + Sliq(k,i,j)) / fsnow(i,j)
        if (Ds(k,i,j) > epsilon(Ds)) then
          rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)
          rhos = rhos + (rhos*grav*mass*dt/(eta0*exp(-(Tsnow(k,i,j) - Tm)/12.4 + rhos/55.6))   &
                      + dt*rhos*snda*exp((Tsnow(k,i,j) - Tm)/23.8 - max(rhos - 150, 0.)/21.7))
          Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rhos / fsnow(i,j)
        end if
        mass = mass + 0.5*(Sice(k,i,j) + Sliq(k,i,j)) / fsnow(i,j)
      end do
    else ! DENSTY == 3
      ! Snow compaction by overburden, dependent on liquid water content (Crocus B92)
      mass = 0
      do k = 1, Nsnow(i,j)
        mass = mass + 0.5*(Sice(k,i,j) + Sliq(k,i,j)) / fsnow(i,j)
        if (Ds(k,i,j) > epsilon(Ds)) then
          rhos = (Sice(k,i,j) + Sliq(k,i,j)) / Ds(k,i,j) / fsnow(i,j)
          f1 = 1 / (1 + 600 * Sliq(k,i,j) / (rho_wat * Ds(k,i,j) * fsnow(i,j)))
          f2 = 1.0
          eta = f1 * f2 * eta1 * (rhos / c_eta) * exp(a_eta * (Tm - Tsnow(k,i,j)) + b_eta * rhos)
          rhos = rhos + rhos*grav*mass*dt/eta
          Ds(k,i,j) = (Sice(k,i,j) + Sliq(k,i,j)) / rhos / fsnow(i,j)
        end if
        mass = mass + 0.5*(Sice(k,i,j) + Sliq(k,i,j)) / fsnow(i,j)
      end do
    endif

    ! Snow grain growth --> *GM for now, this code feature is not functional because the state variable rgrn is not tracked (no bin output)
    do k = 1, Nsnow(i,j)
      ggr = 2e-13
      if (Tsnow(k,i,j) < Tm) then
        if (rgrn(k,i,j) < 1.50e-4) then
          ggr = 2e-14
        else
          ggr = 7.3e-8*exp(-4600/Tsnow(k,i,j))
        end if
      end if
      rgrn(k,i,j) = rgrn(k,i,j) + dt*ggr/rgrn(k,i,j)
    end do
  end if  ! Existing snowpack

  ! Important check: meltflux_out doesn't track liquid water retention during percolation, while Roff_snow does, so here we use Roff_snow to 
  ! constrain meltflux_out; relevant for HYDROL = [1,2]. Minor inconsistencies remain in case of rainfall and percolation through the entire 
  ! snowpack in the same timestep. Verified by LQ and GM, Nov 2021
  if (meltflux_out(i,j) > Roff_snow(i,j)) then
    meltflux_out(i,j) = Roff_snow(i,j)
  end if

  ! Add bare soil runoff to snowmelt runoff for total runoff
  Roff(i,j) = Roff_snow(i,j) + Roff_bare(i,j)
    
  1 continue 
  
end do
end do

! Add snow, recalculate layers
do j = 1, Ny
do i = 1, Nx

  if (tilefrac(i,j) < tthresh) goto 2 ! exclude points outside tile of interest

  ! Falling snow temperature
  Tsnow0 = min(Ta(i,j), Tm)
  if (HN_ON) then
    Tsnow0 = max(Tsnow0, (Tm-40))
  end if

  ! Add snowfall and frost to layer 1 with fresh snow density and grain size
  Esnow = 0
  if (Esrf(i,j) < 0 .and. Tsrf(i,j) < Tm) then
    Esnow = fsnow(i,j) * Esrf(i,j)
    Sbsrf(i,j) = Esnow*dt
  end if 
  dSice = (Sf(i,j) - Esnow)*dt  ! Think about how to scale for fsnow...

  ! Catch to round infinitesimally small new snow amounts.
  ! The small amounts were due to EnKF-assimilated daily precip being downscaled to hourly
  ! based on the hourly (COSMO) to daily (COSMO) mass ratio.  Since min daily EnKF was 1mm,
  ! when 1mm was downscaled to hourly this mass could become
  ! quite small.  When this occurred on bare ground, Tsnow calculation would blow up,
  ! place NaN in Tsnow and snow would never again melt through the season
  ! Do we still need this Catch in FSM??
  if (Nsnow(i,j) <= 1 .and. dSice < .001 .and. Sice(1,i,j) < .001) then
    dSice = FLOAT(INT(dSice * 1000 + 0.5)) / 1000
  end if
  
  if (DENSTY == 0) then
    rhonew = rho0
  else
    if (OSHDTN == 0) then 
      ! Initial formulation
      rhonew = max(rhof + rhob*(Ta(i,j) - Tm) + rhoc*Ua(i,j)**0.5, 50.)
    else ! OSHDTN == 1
      ! New formulation with decompaction
      rhonew = rhof + rhob*(Ta(i,j) - Tm) + rhoc*Ua(i,j)**0.5
      if (dem(i,j) <= 1000) then
        t_decompaction = 24
      else if (dem(i,j) > 2000) then
        t_decompaction = 0
      else
        t_decompaction = 24 + (dem(i,j) - 1000) / (2000 - 1000) * (0 - 24)
      end if
      rhonew = 300 + (rhonew - 300)*exp(t_decompaction/100)
      rhonew = max(rhonew, 50.)
    end if
  endif
  if (Sice(1,i,j) + dSice > epsilon(Sice)) then
    rgrn(1,i,j) = (Sice(1,i,j)*rgrn(1,i,j) + dSice*rgr0) / (Sice(1,i,j) + dSice)
  endif
  Sice(1,i,j) = Sice(1,i,j) + dSice
  snowdepth = sum(Ds(:,i,j)) * fsnow(i,j) + dSice / rhonew

  ! Add canopy unloading to layer 1 with bulk snow density and grain size
  rhos = rhof
  if (Ta(i,j) < Tm) then  ! only if it's cold enough
    mass = sum(Sice(:,i,j)) + sum(Sliq(:,i,j))
    if (snowdepth > epsilon(snowdepth)) then
      rhos = mass / snowdepth
    endif
    Sice(1,i,j) = Sice(1,i,j) + unload(i,j)
    snowdepth = snowdepth + unload(i,j) / rhos
  end if

  ! Store previous snow cover fraction
  fold = fsnow(i,j)
  ! Updated Fractional Snow-Covered Area 
  SWEtmp = sum(Sice(:,i,j) + Sliq(:,i,j))
  call SNOWCOVERFRACTION(snowdepth,SWEtmp,i,j)

  ! Rescale Ds with new snow cover fraction
  if (fsnow(i,j) > epsilon(fsnow)) then
    ! Update surface layer thickness based on new fsnow
    Ds(1,i,j) = Ds(1,i,j) * fold / fsnow(i,j) + dSice / rhonew / fsnow(i,j)
    if (Ta(i,j) < Tm) then
      Ds(1,i,j) = Ds(1,i,j) + unload(i,j) / rhos / fsnow(i,j)
    end if
  end if
  if (Nsnow(i,j) > 1) then
    do k = 2, Nsnow(i,j)
      Ds(k,i,j) = Ds(k,i,j) * fold / fsnow(i,j)
    end do
  end if

  ! New snowpack
  if (Nsnow(i,j) == 0 .and. Sice(1,i,j) > epsilon(Sice(1,i,j))) then
    Nsnow(i,j) = 1
    Tsnow(1,i,j) = Tsnow0
  end if

  ! Store state of old layers
  D(:) = Ds(:,i,j)
  R(:) = rgrn(:,i,j)
  S(:) = Sice(:,i,j)
  W(:) = Sliq(:,i,j)
  if (fsnow(i,j) > epsilon(fsnow)) then
    csnow(1) = (Sice(1,i,j)*hcap_ice + Sliq(1,i,j)*hcap_wat) / fsnow(i,j)
    E(1) = csnow(1)*(Tsnow(1,i,j) - Tm) + ((dSice + unload(i,j)) * hcap_ice / fsnow(i,j)) * (Tsnow0 - Tsnow(1,i,j)) ! Adjustment given that csnow(1) already includes the new snow
  else
    E(:) = 0
  end if
  if (Nsnow(i,j) > 1) then
    do k = 2, Nsnow(i,j)
      csnow(k) = (Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat) / fsnow(i,j)
      E(k) = csnow(k)*(Tsnow(k,i,j) - Tm)
    end do
  end if
  Nold = Nsnow(i,j)

  ! Initialise new layers
  Ds(:,i,j) = 0
  rgrn(:,i,j) = 0
  Sice(:,i,j) = 0
  Sliq(:,i,j) = 0
  Tsnow(:,i,j) = Tm
  U(:) = 0
  Nsnow(i,j) = 0

  if (fsnow(i,j) > epsilon(fsnow)) then  ! Existing or new snowpack

    ! Re-assign and count snow layers
    dnew = snowdepth / fsnow(i,j)
    Ds(1,i,j) = dnew
    k = 1
    if (Ds(1,i,j) > Dzsnow(1)) then 
      do k = 1, Nsmax
        Ds(k,i,j) = Dzsnow(k)
        dnew = dnew - Dzsnow(k)
        if (dnew <= Dzsnow(k) .or. k == Nsmax) then
          Ds(k,i,j) = Ds(k,i,j) + dnew
          exit
        end if
      end do
    end if
    Nsnow(i,j) = k

    ! Fill new layers from the top downwards
    knew = 1
    dnew = Ds(1,i,j)
    do kold = 1, Nold
      do
        if (D(kold) < dnew) then
          ! All snow from old layer partially fills new layer
          rgrn(knew,i,j) = rgrn(knew,i,j) + S(kold)*R(kold)
          Sice(knew,i,j) = Sice(knew,i,j) + S(kold)
          Sliq(knew,i,j) = Sliq(knew,i,j) + W(kold)
          U(knew) = U(knew) + E(kold)
          dnew = dnew - D(kold)
          exit
        else
          ! Some snow from old layer fills new layer
          wt = dnew / D(kold)
          rgrn(knew,i,j) = rgrn(knew,i,j) + wt*S(kold)*R(kold)
          Sice(knew,i,j) = Sice(knew,i,j) + wt*S(kold) 
          Sliq(knew,i,j) = Sliq(knew,i,j) + wt*W(kold)
          U(knew) = U(knew) + wt*E(kold)
          D(kold) = (1 - wt)*D(kold)
          E(kold) = (1 - wt)*E(kold)
          S(kold) = (1 - wt)*S(kold)
          W(kold) = (1 - wt)*W(kold)
          knew = knew + 1
          if (knew > Nsnow(i,j)) exit
          dnew = Ds(knew,i,j)
        end if
      end do
    end do
  
    ! Diagnose snow layer temperatures
    do k = 1, Nsnow(i,j)
      csnow(k) = (Sice(k,i,j)*hcap_ice + Sliq(k,i,j)*hcap_wat) / fsnow(i,j)
      Tsnow(k,i,j) = Tm + U(k) / csnow(k)
      if (HN_ON) then
        Tsnow(k,i,j) = max(Tsnow(k,i,j), (Tm-40))
      endif      
      rgrn(k,i,j) = rgrn(k,i,j) / Sice(k,i,j)
    end do 
  end if ! Existing or new snowpack

  2 continue
  
end do
end do

end subroutine SNOW
