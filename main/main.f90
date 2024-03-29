!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Main program
!
! - Simulate emissions and chemical reactions of gases, aerosol processes as well as 
!   transport of gases and aerosol particles within the planetary boundary layer with a
!   column model.
! - Check Fortran conventions at http://www.fortran90.org/src/best-practices.html
! - Check code conventions at
!   http://www.cesm.ucar.edu/working_groups/Software/dev_guide/dev_guide/node7.html
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program main

use chemistry_mod
use aerosol_mod

implicit none

!-----------------------------------------------------------------------------------------
! Control variables (can be moved to an input file in future)
!-----------------------------------------------------------------------------------------
logical :: use_emission   = .true.
logical :: use_chemistry  = .true.
logical :: use_deposition = .false.
logical :: use_aerosol    = .true.
integer, parameter :: model_number = 3 ! Which K-model to use.
character(len=255), parameter :: input_dir  = './input'
character(len=255), parameter :: output_dir = './output'

!-----------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------
! Double precision
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format
integer, parameter :: dp = selected_real_kind(15, 307)

! Physics constants
real(dp), parameter :: PI     = 2*asin(1.0_dp)                  ! the constant pi
real(dp), parameter :: grav   = 9.81_dp                         ! [m s-2], gravitation
real(dp), parameter :: Rgas   = 8.3144598_dp                    ! [J mol-1 K-1], universal gas constant
real(dp), parameter :: NA     = 6.022140857e23_dp               ! [molec mol-1], Avogadro's number 
real(dp), parameter :: mm_air = 28.96e-3_dp                     ! [kg mol-1], mean molar mass of air
real(dp), parameter :: kb     = 1.38064852e-23_dp               ! [m2 kg s-2 K-1], Boltzmann constant
real(dp), parameter :: Cp     = 1012.0_dp                       ! [J kg-1 K-1], air specific heat at constant pressure,
real(dp), parameter :: p00    = 1.01325e5_dp                    ! [Pa], reference pressure at surface
real(dp), parameter :: nu_air = 1.59e-5_dp                      ! [m2 s-1], kinematic viscosity of air
real(dp), parameter :: Omega  = 2*PI/(24.0_dp*60.0_dp*60.0_dp)  ! [rad s-1], Earth angular speed
real(dp), parameter :: lambda = 300.0_dp                        ! maximum mixing length, meters
real(dp), parameter :: vonk   = 0.4_dp                          ! von Karman constant, dimensionless
real(dp), parameter :: ppb = 1e-9_dp                            ! Parts per billion

real(dp), parameter :: ug = 10.0d0, vg = 0.0d0  ! [m s-1], geostrophic wind

! Latitude and longitude of Hyytiala
! real(dp), parameter :: latitude_deg  = 61.8455d0  ! [degN]
real(dp), parameter :: latitude_deg  = 56.1d0  ! [degN]
!real(dp), parameter :: longitude_deg = 24.2833d0  ! [degE]
real(dp), parameter :: longitude_deg = 13.42d0  ! [degE]
real(dp), parameter :: latitude      = latitude_deg  * PI/180.0d0  ! [rad]
real(dp), parameter :: longitude     = longitude_deg * PI/180.0d0  ! [rad]

real(dp), parameter :: fcor = 2*Omega*sin(latitude)  ! Coriolis parameter at Hyytiala

!-----------------------------------------------------------------------------------------
! Grid parameters
!-----------------------------------------------------------------------------------------
integer, parameter :: nz = 50  ! [-], number of height levels

! Model height levels, [m]
real(dp), parameter, dimension(nz) :: &
  hh = (/    0,   10,   20,   30,   40,   50,   60,   70,   80,   90, &
           100,  120,  140,  160,  180,  200,  230,  260,  300,  350, &
           400,  450,  500,  550,  600,  650,  700,  800,  900, 1000, &
          1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, &
          2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 /)

real(dp), parameter :: hc = 10.0_dp  ! [m], canopy height

!-----------------------------------------------------------------------------------------
! Time variables
!-----------------------------------------------------------------------------------------
integer, parameter :: one_hour = 60*60  ! [s], one hour in seconds

real(dp) :: time                  ! [s], current time
real(dp) :: time_start, time_end  ! [s], start and end time

real(dp) :: dt         ! [s], time step for main loop, usually is equal to meteorology time step
real(dp) :: dt_emis    ! [s], time step for emission calculation
real(dp) :: dt_chem    ! [s], time step for chemistry calculation
real(dp) :: dt_depo    ! [s], time step for deposition calculation
real(dp) :: dt_aero    ! [s], time step for aerosol calculation
real(dp) :: dt_output  ! [s], time step for output

real(dp) :: time_start_emission    ! [s], time to start calculating emission
real(dp) :: time_start_chemistry   ! [s], time to start calculating chemistry
real(dp) :: time_start_deposition  ! [s], time to start calculating deposition
real(dp) :: time_start_aerosol     ! [s], time to start calculating aerosol

integer :: daynumber_start  ! [day], start day of year
integer :: daynumber        ! [day], current day of year
integer :: days

integer :: counter  ! [-], counter of time steps

!-----------------------------------------------------------------------------------------
! Meteorology variables
!-----------------------------------------------------------------------------------------
real(dp), dimension(nz  ) :: uwind, &  ! [m s-1], u component of wind
                             vwind, &  ! [m s-1], v component of wind
                             theta     ! [K], potential temperature
real(dp), dimension(nz  ) :: temp, &   ! [K], air temperature
                             pres      ! [Pa], air pressure

! Defining K-parameters, both integer and half-values. There are only 49 intermediate values
! For each value K_m(i) the half-value above it is K_m_half(i)
real(dp), dimension(nz) :: K_m = 0.0_dp ! [m/s] K_m-parameter for the integer values
real(dp), dimension(nz-1) :: K_m_half = 0.0_dp ! [m/s] K_m parameter for the half-values
real(dp), dimension(nz) :: K_h = 0.0_dp ! [m/s] K_h-parameter for the integer values
real(dp), dimension(nz-1) :: K_h_half = 0.0_dp ! [m/s] K_h-parameter for the integer values
real(dp), dimension(nz-1) :: richard = 0.0_dp ! Richardson number for thermal stability, I only define this for half-integer values of height.
real(dp):: richards_nr10m = 0.0_dp ! Richardson number for the second lowest layer (10) where aerosols are removed.

! Final vectors to store the changes
real(dp), dimension(nz) :: du_dt ! Derivative of u with regards to time at integer values
real(dp), dimension(nz) :: dv_dt ! Derivative of v with regards to time at integer values
real(dp), dimension(nz) :: dtheta_dt ! Derivative of theta with regards to time
real(dp), dimension(nz) :: dc_dt ! Derivative of scalar concentration with regards to time  
real(dp) :: emission_isoprene = 0.0_dp ! Emission rate of isoprene for a given time step
real(dp) :: emission_monoterpene = 0.0_dp ! Emission rate of alpha-pinene for a given time step.   

! Chemistry vectors and matrixes:
real(dp), dimension(nz) :: Mair ! [molecules/cm3] Air molecules concentration at each level.
real(dp), dimension(nz) :: O2 ! [molecules/cm3] Oxygen concentration at each level
real(dp), dimension(nz) :: N2 ! [molecules/cm3] Nitrogen concentration at each level
real(dp), dimension(nz) :: H2O ! [molecules/cm3] Water concentration at each level
real(dp), dimension(nz,neq) :: concentration ! [molecules/cm3] Chemical species concentrations (row) for each level (column)

! Aerosol vectors and matrixes:
integer, parameter :: nr_bins = 100 ! Number of particle sizes used in the aerosol code
integer, parameter :: nr_cond = 2 ! Number of condensable vapours.
real(dp), dimension(nz) :: PN,& ! [# m^-3] Total number concentration for each height level
                           PM,& ! [kg m^-3] Total mass concentration for each height level
                           PV ! [um^3 cm^-3] Total volume concentration for each height level
real(dp), dimension(nz,nr_bins) :: aerosol_conc ! Aerosol concentrations
! where columns denote the height level and rows denote the size bins.
real(dp), dimension(nr_cond) :: cond_vapour ! [molec/m^3] Vector of concentrations of condensable vapours
real(dp), dimension(nz) :: daero_dt ! Derivative of aerosol with regards to time 
real(dp), dimension(nz,nr_cond) :: cond_sink ! Condensation sink for each height level for H2SO4 and ELVOC
real(dp) :: DSWF ! [W m-2], Downward Shortwave Radiation Flux
real(dp) :: wind_speed10m ! Total windspeed at 10 m height

integer :: i, j  ! used for loops


!-----------------------------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------------------------

call time_init()                 ! initialize time
call meteorology_init()          ! initialize meteorology
call chemistry_init()            ! initialize chemistry
! I don't need to do this 50 times sinces the initial aerosol concentrations are identical in each layer
! So I could have done it once and then copied over aerosol_conc to all the other layers
do i = 1, nz
  call aerosol_init(diameter, particle_mass, particle_volume, aerosol_conc(i,:), &
  particle_density, nucleation_coef, molecular_mass, molar_mass, &
  molecular_volume, molecular_dia, mass_accomm) ! initialize aerosol variables and values
end do

cond_sink = 0.001d0

! Calculate initial values of PN, PM, PV
do i = 1, nz
  PN(i) = sum(aerosol_conc(i,:)) * 1D-6                 ! [# cm-3], total particle number concentration
  PM(i) = sum(aerosol_conc(i,:)*particle_mass) * 1D9    ! [ug m-3], total particle mass concentration
  PV(i) = sum(aerosol_conc(i,:)*particle_volume) * 1D12 ! [um^3 cm^-3] Total particle volume
end do

call open_files()        ! open output files
call write_files(time)   ! write initial values

!-----------------------------------------------------------------------------------------
! Start main loop
!-----------------------------------------------------------------------------------------
do while (time <= time_end)
  !write(*,*) cond_sink
  ! Passage of time:
  !days = aint(time / (24*60*60)) ! truncates time / 24 to the integer, i.e. only increments from 0 when 24 hours passes
  !daynumber = daynumber_start + days ! increases the day number if 24 hours passes.
  !---------------------------------------------------------------------------------------
  ! Meteorology
  !---------------------------------------------------------------------------------------
  ! Set lower boundary condition
  call surface_values(theta(1), time+dt)  ! theta = temperature at the surface

  ! Update meteorology

  ! Determine the K_values at this time
  call K_values(K_m,K_m_half,K_h,K_h_half,richard)
  ! First the change of wind speeds, du_dt, dv_dt, du_dt_half, dv_dt_half
  call wind_derivatives(uwind,vwind,du_dt,dv_dt)

  ! Then the update of temperature diffusion
  ! Once again boundary conditions are defined and constant.
  dtheta_dt(1) = 0.0_dp
  dtheta_dt(nz) = 0.0_dp
  do i = 2, nz - 1
    dtheta_dt(i) = (K_h_half(i) * ((theta(i + 1) - theta(i)) / (hh(i + 1) - hh(i))) - & 
    K_h_half(i - 1) * ((theta(i) - theta(i - 1)) / (hh(i) - hh(i - 1)))) / ((hh(i + 1) - hh(i - 1)) / 2)
  end do

  ! Update the values
  ! uwind(time+1) = uwind(time)+delta(time)*du_dt
  do i = 2, nz-1
    uwind(i) = uwind(i)+dt*du_dt(i)
    vwind(i) = vwind(i)+dt*dv_dt(i)
    theta(i) = theta(i)+dt*dtheta_dt(i)
  end do
  !---------------------------------------------------------------------------------------
  ! Emission
  !---------------------------------------------------------------------------------------
  ! Start to calculate emission after time_start_emission
  ! Compute emission part every dt_emis, multiplying 1000 to convert s to ms to make mod easier
  if ( use_emission .and. time >= time_start_emission ) then
    if ( mod( nint((time - time_start_emission)*1000.0d0), nint(dt_emis*1000.0d0) ) == 0 ) then
      ! Calculate emission rates
      call calc_emission_rates(emission_monoterpene,emission_isoprene,theta,time,daynumber)
      ! Add emissions to concentrations at level 2:
      !concentration(2,13) = concentration(2,13) + emission_isoprene ! Isoprene
      !concentration(2,23) = concentration(2,23) + emission_monoterpene ! Monoterpene (alpha-pinene)
    end if
  end if

  if ( use_emission .and. (.not. use_chemistry) ) then
    ! Add emission to the number concentrations of compounds
  end if

  !---------------------------------------------------------------------------------------
  ! Deposition
  !---------------------------------------------------------------------------------------
  ! Start to calculate gas dry deposition velocity after time_start_deposition
  ! Compute deposition part every dt_depo, multiplying 1000 to convert s to ms to make mod easier
  if ( use_deposition .and. time >= time_start_deposition ) then
    if ( mod( nint((time - time_start_deposition)*1000.0d0), nint(dt_depo*1000.0d0) ) == 0 ) then
      ! First calculate Richardssons number for the second layer (10m)
      richards_nr10m = (hh(3)-hh(2)) * grav / ((theta(3) - theta(2)) / 2.0_dp) * (theta(3) - theta(2)) / ((uwind(3) - uwind(2))**2 + (vwind(3) - vwind(2))**2)
      ! Calculate temp and pressure
      temp = theta - (grav/Cp)*hh
      pres = barometric_law(p00, temp, hh)
      ! Calculate DSWF
      DSWF = 6D2 * get_exp_coszen(time,daynumber,latitude)
      ! Wind speed
      wind_speed10m = sqrt(uwind(2)**2 + vwind(2)**2)
      ! Calculate deposition velocity, remove the deposited concentration at
      ! level 2 which includes canopy and soil:
      !write(*,*) "Depositing"
      !write(*,*) aerosol_conc(2,10)
      call dry_dep_velocity(diameter,particle_density,temp(2),pres(2),DSWF, & 
      richards_nr10m,wind_speed10m,aerosol_conc(2,:),concentration(2,:),dt_depo,hh(2))
      !write(*,*) aerosol_conc(2,10)
    end if
  end if

  !---------------------------------------------------------------------------------------
  ! Chemistry
  !---------------------------------------------------------------------------------------
  ! Start to calculate chemical reactions only after some time to save the computation time
  ! Compute chemistry part every dt_chem, multiplying 1000 to convert s to ms to make mod easier
  if ( use_chemistry .and. time >= time_start_chemistry ) then
    if ( mod( nint((time - time_start_chemistry)*1000.0d0), nint(dt_chem*1000.0d0) ) == 0 ) then
      ! Solve chemical equations for each layer except boundaries
      call chemistry_1D(time,theta,Mair,O2,N2,H2O,concentration,dt_chem,emission_isoprene,emission_monoterpene,cond_sink)
    end if  ! every dt_chem
  end if

  ! Update concentrations of gas phase compounds if any of these processes are considered
  ! Deposition should not be used alone because it calculates nothing in that case
  if (use_emission .or. use_chemistry) then
    ! Trick to make bottom flux zero, from boundary conditions
    ! Setting concentrations at level 1 equal to level 2:
    concentration(1,:) = concentration(2,:)

    ! Concentrations can not be lower than 0
    concentration(:,:) = max(concentration(:,:),0.0_dp)

    ! Mixing of chemical species
    ! Flux at the bottom is guaranteed to be zero due to the above.
    ! Loop over each chemical, and inside we loop over each height-level
    ! except for the lowest and highest

    do i = 1, neq
      ! Quickly initialize the concentration change vector for the new chemical:
      dc_dt = 0.0_dp
      do j = 2, nz-1
        dc_dt(j) = (K_h_half(j) * ((concentration(j+1,i) - concentration(j,i)) / (hh(j + 1) - hh(j))) - & 
        K_h_half(j - 1) * ((concentration(j,i) - concentration(j-1,i)) / (hh(j) - hh(j - 1)))) / ((hh(j + 1) - hh(j - 1)) / 2.0_dp)
      end do
      do j = 2, nz-1
        concentration(j,i) = concentration(j,i) + dt*dc_dt(j)
      end do
    end do

    ! Set the boundary_conditions above again for output
    ! And set top concentrations to zero, to simulate chemicals
    ! leaving the box.
    concentration(1,:) = concentration(2,:)
    concentration(nz,:) = 0.0_dp
  
  end if

  !---------------------------------------------------------------------------------------
  ! Aerosol
  !---------------------------------------------------------------------------------------
  ! Start to calculate aerosol processes only after some time to save the computation time
  ! Compute aerosol part every dt_aero, multiplying 1000 to convert s to ms to make mod easier
  if ( use_aerosol .and. time >= time_start_aerosol ) then
    if ( mod( nint((time - time_start_aerosol)*1000.0d0), nint(dt_aero*1000.0d0) ) == 0 ) then
      ! Nucleation, condensation, coagulation and deposition of particles
      ! Make sure temperature and pressure is updated:
      temp = theta - (grav/Cp)*hh
      pres = barometric_law(p00, temp, hh)
      ! First reset the condensation sink
      cond_sink = 0.001d0
      do i = 2, nz-1
        ! Find the concentrations of the condensable vapours:
        cond_vapour(1) = concentration(i,21) * 1d6
        cond_vapour(2) = concentration(i,25) * 1d6
        call nucleation(nucleation_coef,cond_vapour(1),aerosol_conc(i,:),dt_aero)
        call condensation(dt_aero, temp(i), pres(i), mass_accomm, molecular_mass, &
        molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
        aerosol_conc(i,:), diameter, cond_vapour, cond_sink(i,:))
        call coagulation(dt_aero, aerosol_conc(i,:), diameter, temp(i), pres(i), particle_mass)
      end do
    end if

    ! Trick to make bottom flux zero
    aerosol_conc(1,:) = aerosol_conc(2,:)
    ! Concentrations can not be lower than 0 [molec m-3]
    aerosol_conc(:,:) = max(aerosol_conc(:,:),0.0_dp)
    ! Mixing of aerosol particles
    do i = 1, nr_bins
      ! Quickly initialize the concentration change vector for the new chemical:
      daero_dt = 0.0_dp
      do j = 2, nz-1
        daero_dt(j) = (K_h_half(j) * ((aerosol_conc(j+1,i) - aerosol_conc(j,i)) / (hh(j + 1) - hh(j))) - & 
        K_h_half(j - 1) * ((aerosol_conc(j,i) - aerosol_conc(j-1,i)) / (hh(j) - hh(j - 1)))) / ((hh(j + 1) - hh(j - 1)) / 2.0_dp)
      end do
      do j = 2, nz-1
        aerosol_conc(j,i) = aerosol_conc(j,i) + dt*daero_dt(j)
      end do
    end do
    ! Set the top to be zero to simulate aerosols leaving the box:
    !aerosol_conc(nz,:) = 0.0_dp
    ! Set the constraints above again for output
    aerosol_conc(1,:) = aerosol_conc(2,:)
    aerosol_conc(:,:) = max(aerosol_conc(:,:),0.0_dp)
    ! Update related values, e.g., total number concentration, total mass concentration

    do i = 1, nz
      PN(i) = sum(aerosol_conc(i,:)) * 1D-6                 ! [# cm-3], total particle number concentration
      PM(i) = sum(aerosol_conc(i,:)*particle_mass) * 1D9    ! [ug m-3], total particle mass concentration
      PV(i) = sum(aerosol_conc(i,:)*particle_volume) * 1D12 ! [um^3 cm^-3] Total particle volume
    end do
  end if

  !---------------------------------------------------------------------------------------
  ! Ending loop actions
  !---------------------------------------------------------------------------------------
  ! Advance to next time step
  time = time + dt

  ! Write data every dt_output [s]
  if ( mod( nint((time - time_start)*1000.0d0), nint(dt_output*1000.0d0) ) == 0 ) then
    write(*, '(a8,f8.3,a8)') 'time = ', time/one_hour, '   hours'
    call write_files(time)
  end if

  ! Count loop number
  counter = counter + 1

end do

!-----------------------------------------------------------------------------------------
! Finalization
!-----------------------------------------------------------------------------------------
! Close all the opened files
call close_files()

! Count total time steps
write(*,*) counter,'time steps'


contains


!-----------------------------------------------------------------------------------------
! subroutine open_files()
!
! Open needed files
!-----------------------------------------------------------------------------------------
subroutine open_files()
  logical :: dir_exist

  ! Create a new directory if it does not exist
  inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
  if (.not. dir_exist) then
    ! This line may change for different operating systems
    call system('mkdir ' // trim(adjustl(output_dir)))
  end if

  ! Open files to write output results
  open( 8,file=trim(adjustl(output_dir))//'/time.dat' ,status='replace',action='write')
  open( 9,file=trim(adjustl(output_dir))//'/hh.dat'   ,status='replace',action='write')
  open(10,file=trim(adjustl(output_dir))//'/uwind.dat',status='replace',action='write')
  open(11,file=trim(adjustl(output_dir))//'/vwind.dat',status='replace',action='write')
  open(12,file=trim(adjustl(output_dir))//'/theta.dat',status='replace',action='write')
  
  open(13,file=trim(adjustl(output_dir))//'/K_m.dat',status='replace',action='write')
  open(14,file=trim(adjustl(output_dir))//'/K_h.dat',status='replace',action='write')
  open(15,file=trim(adjustl(output_dir))//'/richard.dat',status='replace',action='write')

  open(16,file=trim(adjustl(output_dir))//'/emission_isoprene.dat',status='replace',action='write')
  open(17,file=trim(adjustl(output_dir))//'/emission_monoterpene.dat',status='replace',action='write')

  open(18,file=trim(adjustl(output_dir))//'/alpha_pinene.dat',status='replace',action='write')
  open(19,file=trim(adjustl(output_dir))//'/isoprene.dat',status='replace',action='write')
  open(20,file=trim(adjustl(output_dir))//'/oh_radical.dat',status='replace',action='write')
  open(21,file=trim(adjustl(output_dir))//'/ho2_radical.dat',status='replace',action='write')
  open(22,file=trim(adjustl(output_dir))//'/h2so4.dat',status='replace',action='write')
  open(23,file=trim(adjustl(output_dir))//'/elvoc.dat',status='replace',action='write')

  open(24,file=trim(adjustl(output_dir))//'/PN.dat',status='replace',action='write')
  open(25,file=trim(adjustl(output_dir))//'/PM.dat',status='replace',action='write')
  open(26,file=trim(adjustl(output_dir))//'/PV.dat',status='replace',action='write')

  open(27,file=trim(adjustl(output_dir))//'/aerosol_conc_1.dat',status='replace',action='write')
  open(28,file=trim(adjustl(output_dir))//'/diameter.dat',status='replace',action='write')
end subroutine open_files


!-----------------------------------------------------------------------------------------
! subroutine write_files(time)
!
! Write data to files at time
!-----------------------------------------------------------------------------------------
subroutine write_files(time)
  real(dp) :: time  ! current time
  character(255) :: outfmt_one_scalar, outfmt_two_scalar, outfmt_level, outfmt_mid_level, outfmt_size_bins

  ! Output real data with scientific notation with 16 decimal digits
  outfmt_one_scalar = '(es25.16e3)'                               ! for scalar
  write(outfmt_level     , '(a, i3, a)') '(', nz  , 'es25.16e3)'  ! for original levels
  write(outfmt_mid_level , '(a, i3, a)') '(', nz-1, 'es25.16e3)'  ! for middle levels
  write(outfmt_two_scalar, '(a, i3, a)') '(', 2   , 'es25.16e3)'  ! for two scalars
  write(outfmt_size_bins, '(a, i3, a)') '(', nr_bins , 'es25.16e3)' ! For size distributions

  ! Only output hh once
  if (time == time_start) then
    write(9, outfmt_level) hh
  end if

  ! Output every output time step
  write( 8, outfmt_one_scalar) time/(24*one_hour)  ! [day]
  write(10, outfmt_level     ) uwind
  write(11, outfmt_level     ) vwind
  write(12, outfmt_level     ) theta

  write(13, outfmt_mid_level     ) K_m_half
  write(14, outfmt_mid_level     ) K_h_half
  write(15, outfmt_mid_level     ) richard

  write(16, outfmt_level     ) emission_isoprene
  write(17, outfmt_level     ) emission_monoterpene

  write(18,outfmt_level) concentration(:,23) ! Alpha pinene concentrations
  write(19,outfmt_level) concentration(:,13) ! Isoprene concentrations
  write(20,outfmt_level) concentration(:,3) ! OH concentrations
  write(21,outfmt_level) concentration(:,8) ! HO2 concentrations
  write(22,outfmt_level) concentration(:,21) ! H2SO4 concentrations
  write(23,outfmt_level) concentration(:,25) ! ELVOC concentrations

  write(24, outfmt_level) PN
  write(25, outfmt_level) PM
  write(26, outfmt_level) PV

  write(27, outfmt_size_bins) aerosol_conc(1,:) ! Aerosol number concentration in the 1st model layer
  write(28, outfmt_size_bins) diameter ! Diameters corresponding to the size bins.
end subroutine write_files


!-----------------------------------------------------------------------------------------
! subroutine Close_Files()
!
! Close files
!-----------------------------------------------------------------------------------------
subroutine close_files()
  close(8)
  close(9)
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  close(26)
  close(27)
  close(27)
  close(28)
end subroutine close_files


!-----------------------------------------------------------------------------------------
! subroutine time_init()
!
! Time initiation
!-----------------------------------------------------------------------------------------
subroutine time_init()
  ! Basic time variables
  time_start = 0.0d0
  time_end   = 5.0d0 * 24.0d0 * one_hour
  !time_end   = 4.0d0 * 24.0d0 * one_hour
  time       = time_start

  ! Time steps
  dt        = 0.5d0
  dt_emis   = 0.5d0
  dt_chem   = 10.0d0
  dt_depo   = 10.0d0
  dt_aero   = 10.0d0
  dt_output = 3600.0d0

  ! Day number
  !daynumber_start = 31+28+31+30+31+30+10  ! day is July 10th, 2011
  daynumber_start = 31+28+31+6 ! 6th of april
  daynumber       = daynumber_start

  ! Start time for each process
  !time_start_emission   = 3*24*one_hour
  time_start_emission   = 4.625d0*24*one_hour
  time_start_chemistry  = 3*24*one_hour
  time_start_deposition = 3*24*one_hour
  !time_start_aerosol    = 3*24*one_hour
  time_start_aerosol    = 4*24*one_hour

  ! Loop number
  counter = 0
end subroutine time_init


!-----------------------------------------------------------------------------------------
! subroutine meteorology_init()
!
! Meteorology initiation
!-----------------------------------------------------------------------------------------
subroutine meteorology_init()
  ! Wind velocity
  uwind         = 0.0d0
  uwind(nz)     = ug
  uwind(2:nz-1) = uwind(nz) * hh(2:nz-1)/hh(nz)

  vwind = 0.0d0
  vwind(nz) = vg

  ! Potential temperature
  !theta     = 273.15d0 + 25.0d0
  !theta(nz) = 273.15d0 + 30.0d0
  theta     = 273.15d0 + 0.0d0
  theta(nz) = 273.15d0 + 5.0d0

  ! Air temperature and pressure
  temp = theta - (grav/Cp)*hh
  pres = barometric_law(p00, temp, hh)
end subroutine meteorology_init

subroutine chemistry_init()
  temp = theta - (grav/Cp)*hh
  pres = barometric_law(p00, temp, hh)
  concentration = 0.0_dp
  ! Determine the concentration of air for each height
  do i = 1, nz
    Mair(i) = pres(i)*NA / (Rgas*temp(i)) * 1d-6
  end do
  ! Set initial concentrations under the assumption that the mixing ratio
  ! is homogenous throughout the air parcel
  do i = 1, nz
    O2(i) = 0.21_dp * Mair(i) ! Concentrations of O2 at the different levels
    N2(i) = 0.78_dp * Mair(i) ! Concentrations of N2 at the different levels
    !H2O(i) = 1.0D16 * (Mair(i)/Mair(1)) ! Water vapour concentratoin at the different levels
    H2O(i) = 1.0D16 ! Water vapour concentratoin at the different levels
    !concentration(i,1) = 24.0_dp * Mair(i) * ppb ! O3 concentration
    concentration(i,1) = 40.0_dp * Mair(i) * ppb ! O3 concentration
    concentration(i,5) = 0.2_dp * Mair(i) * ppb ! NO2 concentration
    concentration(i,6) = 0.07_dp * Mair(i) * ppb ! NO concentration
    !concentration(i,9) = 100.0_dp * Mair(i) * ppb ! CO concentration
    concentration(i,9) = 200.0_dp * Mair(i) * ppb ! CO concentration
    concentration(i,11) = 1759.0_dp * Mair(i) * ppb ! CH4 concentration
    !concentration(i,20) = 0.5_dp * Mair(i) * ppb ! SO2 concentration
    concentration(i,20) = 2.0_dp * Mair(i) * ppb ! SO2 concentration
  end do

end subroutine

!-----------------------------------------------------------------------------------------
! Get the surface values from the input data file
! Now only the temperature is used.
!-----------------------------------------------------------------------------------------
subroutine surface_values(temperature, time)

  ! (Note: can also get water concentrantion, in ppt, if modify this
  ! subroutine to also use column 8)
  !
  ! Data is taken from:
  ! http://avaa.tdata.fi/web/smart

  real(dp), intent(in)            :: time ! input, in seconds
  real(dp), intent(out)           :: temperature ! output, in Kelvin
  logical, save                   :: first_time = .true.
  real(dp), dimension(8,50), save :: surface_data
  real(dp), dimension(50), save   :: temperature_data
  real(dp), parameter             :: seconds_in_day = 24*60*60
  real(dp), parameter             :: seconds_in_30min = 30*60
  integer                         :: index
  real(dp) :: time24h, time30min, time24plus15, temp1, temp2, x

  ! Only when called for the first time, read in data from file
  ! With this trick, we don't need to open the file in the main program
  IF (first_time) THEN
     !open(30, file=trim(adjustl(input_dir))//'/hyytiala_20110710-t_h2o.dat', status='old')
     open(30, file=trim(adjustl(input_dir))//'/hyltemossa_2018_4_06_t_h2o.dat', status='old')
     read(30, *) surface_data
     temperature_data(1:50) = surface_data(7,1:50) ! in Celcius
     first_time = .false.
  end IF

  time24h = modulo(time, seconds_in_day) ! time modulo 24 hours
  time24plus15 = time24h + 15*60 ! time since 23:45 previous day
  time30min = modulo(time24plus15, seconds_in_30min)
  index = 1 + floor(time24plus15/seconds_in_30min)

  temp1 = temperature_data(index)
  temp2 = temperature_data(index + 1)
  x = time30min/seconds_in_30min

  ! linear interpolation between previous and next temperature data value
  temperature = temp1 + x*(temp2 - temp1) + 273.15_dp  ! now in Kelvin

end subroutine surface_values


subroutine K_values(K_m,K_m_half,K_h,K_h_half,richard)
  real(dp), dimension(nz) :: K_m ! [m/s] K_m-parameter for the integer values
  real(dp), dimension(nz-1) :: K_m_half ! [m/s] K_m parameter for the half-values
  real(dp), dimension(nz) :: K_h ! [m/s] K_h-parameter for the integer values
  real(dp), dimension(nz-1) :: K_h_half ! [m/s] K_h-parameter for the integer values
  real(dp) :: dv_bar_dt ! Absolute value of the derivatives of the wind values at a given height
  real(dp) :: L ! The L factor at a given height, used in model 2 and 3.
  real(dp) :: du_dz, dv_dz ! Derivatives of the wind speeds with regards to height
  real(dp) :: dtheta_dz ! Derivative of the temperature with regards to height.
  real(dp) :: half_height ! Height in the in-between level
  real(dp), dimension(nz-1) :: richard ! Richardsons number for thermal stability.
  real(dp) :: f_m ! Function value of Richardsons number for wind
  real(dp) :: f_h ! Function value of Richardsons number for diffusion
  real(dp) :: theta_extrapolated ! Extrapolated temperature between two integer levels

  select case (model_number)
  case (1)
    ! Model 1 assumes that the K values are constant at 5 m/s
    K_m = 5.0_dp
    K_m_half = 5.0_dp
    K_h = 5.0_dp
    K_h_half = 5.0_dp
  case (2)
    ! First determine the absolute value of the wind speed vector for the different half-heights
    ! and the value of the L factor at the half-heights
    ! I calculate only the K values for the half_values
    do i = 1, nz-1
      du_dz = (uwind(i + 1) - uwind(i)) / (hh(i + 1) - hh(i))
      dv_dz = (vwind(i + 1) - vwind(i)) / (hh(i + 1) - hh(i))
      dv_bar_dt = sqrt(du_dz**2 + dv_dz**2)
      ! Quickly get the heigth at the half-height:
      half_height = half_z(hh(i), hh(i + 1))
      L = vonk * half_height / (1 + (vonk * half_height / lambda))
      K_m_half(i) = L**2 * dv_bar_dt
      ! They are defined to be identical for identical heights:
      K_h_half(i) = K_m_half(i)
    end do
  case (3)
    do i = 1, nz-1
      ! First find the derivatives of the parameters with respect to height:
      du_dz = (uwind(i + 1) - uwind(i)) / (hh(i + 1) - hh(i))
      dv_dz = (vwind(i + 1) - vwind(i)) / (hh(i + 1) - hh(i))
      dtheta_dz = (theta(i + 1) - theta(i)) / (hh(i + 1) - hh(i))
      dv_bar_dt = sqrt(du_dz**2 + dv_dz**2)
      ! Find the height at the half-height:
      half_height = half_z(hh(i), hh(i + 1))
      L = vonk * half_height / (1 + (vonk * half_height / lambda))
      ! Find potential temperature at the half-height by extrapolation:
      theta_extrapolated = (theta(i + 1) + theta(i)) / 2
      ! Calculate Richards Numbers
      richard(i) = (grav / theta_extrapolated) * (dtheta_dz / (du_dz**2 + dv_dz**2))
      ! Determine the case based on the magnitude and sign and calculate the correction factor:
      if (richard(i) < 0.0) then
        f_m = sqrt(1 - 16 * richard(i))
        f_h = (1 - 16 * richard(i))**(3.0/4.0)
      else if (richard(i) < 0.2 .and. richard(i) >= 0.0) then
        f_m = max((1 - 5 * richard(i))**2, 0.1_dp)
        f_h = f_m
      else if (richard(i) >= 0.2) then
        f_m = 0.1_dp
        f_h = f_m
      end if
      K_m_half(i) = L**2 * dv_bar_dt * f_m
      K_h_half(i) = L**2 * dv_bar_dt * f_h
    end do
  end select

end subroutine K_values

subroutine calc_emission_rates(emission_monoterpene,emission_isoprene,theta,time,daynumber)
  implicit none
  real(dp) :: emission_isoprene ! Emission rate of isoprene for a given time step
  real(dp) :: emission_monoterpene ! Emission rate of alpha-pinene for a given time step.
  real(dp), dimension(nz) :: theta ! [K], potential temperature
  real(dp) :: time ! [s], current time
  integer :: daynumber ! [day], current day of year
  real(dp) :: PAR ! [unit!!] Photosynthesis Active Radiation 
  real(dp) :: gamma_isoprene
  real(dp) :: gamma_monoterpene
  real(dp) :: C_L
  real(dp) :: C_T
  ! Empirical constants for the equations:
  real(dp), parameter :: alpha = 0.0027_dp ! Unitless
  real(dp), parameter :: beta = 0.09_dp ! [K^-1]
  real(dp), parameter :: C_L1 = 1.006_dp ! Unitless
  real(dp), parameter :: C_T1 = 95.0_dp*1000.0_dp ! [J/mol]
  real(dp), parameter :: C_T2 = 230.0_dp*1000.0_dp ! [J/mol]
  real(dp), parameter :: T_S = 303.15_dp ! [K]
  real(dp), parameter :: T_M = 314.0_dp ! [K]
  real(dp), parameter :: D_m = 0.0538_dp ! [g(dry mass)/cm^2] Standing leaf biomass, distributed between levels 1 and 2.
  real(dp), parameter :: emission_factor = 100_dp ! [ng(VOC) / ( g(needle dry mass) * h ) ]
  real(dp), parameter :: delta = 1.0_dp ! Emission activity factor for long term control

  real(dp), parameter :: molar_mass_isoprene = 68.12_dp ! [g/mol] molarmass of C5H8
  real(dp), parameter :: molar_mass_monoterpene = 136.24_dp ! [g/mol] molarmass of C10H16

  ! First get the PAR value:
  PAR = 1000*get_exp_coszen(time,daynumber,latitude)

  ! Find the real air temperatures, using the simplified equation and the originally initiated temp vector:
  temp = theta - (grav/Cp)*hh 
  ! We use the temperature at level 2 (between 5-15 m)

  ! First I calculate the isoprene emission:
  C_L = alpha * C_L1 * PAR / (sqrt(1 + alpha**2 * PAR**2))
  C_T = exp(C_T1 * (temp(2) - T_S) / (Rgas * temp(2) * T_S)) / (1 + exp(C_T2 * (temp(2) - T_M) / (Rgas * temp(2) * T_S)))
  gamma_isoprene = C_L*C_T
  
  ! Calculate the emission flux, convert ng to g, h to s, g to mol, mol to molecules.
  ! Finally divide it by the height interval (between level 1 and 2) to get it per volume (# cm^-3 s^-1)
  emission_isoprene = D_m * emission_factor * gamma_isoprene * delta * 10.0**(-9) / (60 * 60 * molar_mass_isoprene) &
                      * (NA / (100 * (hh(2) - hh(1))))
  
  ! Then monoterpene:
  gamma_monoterpene = exp(beta * (temp(2) - T_S))
  emission_monoterpene = D_m * emission_factor * gamma_monoterpene * delta * 10.0**(-9) / (60 * 60 * molar_mass_monoterpene) &
                        * (NA / (100 * (hh(2) - hh(1))))
end subroutine

subroutine wind_derivatives(uwind,vwind,du_dt,dv_dt)
real(dp), dimension(nz) :: uwind, &  ! [m s-1], u component of wind
                           vwind  ! [m s-1], v component of wind
real(dp), dimension(nz) :: du_dt ! Derivative of u with regards to time at integer values
real(dp), dimension(nz) :: dv_dt ! Derivative of v with regards to time at integer values

! First I calculate for the integer values:
! Loop over the second z value to the second last, the excluded values are
! constant and determined by the boundary conditions:
du_dt(1) = 0
du_dt(nz) = 0
dv_dt(1) = 0
dv_dt(nz) = 0
do i = 2, nz - 1
  ! First update du_dt
  du_dt(i) = fcor*(vwind(i)-vg) + &
  (K_m_half(i) * ((uwind(i + 1) - uwind(i)) / (hh(i + 1) - hh(i))) - & 
  K_m_half(i - 1) * ((uwind(i) - uwind(i - 1)) / (hh(i) - hh(i - 1)))) / ((hh(i + 1) - hh(i - 1)) / 2)
  ! Then the same for dv_dt
  dv_dt(i) = -fcor*(uwind(i)-ug) + &
  (K_m_half(i) * ((vwind(i + 1) - vwind(i)) / (hh(i + 1) - hh(i))) - & 
  K_m_half(i - 1) * ((vwind(i) - vwind(i - 1)) / (hh(i) - hh(i - 1)))) / ((hh(i + 1) - hh(i - 1)) / 2)
end do

end subroutine wind_derivatives

subroutine chemistry_1D(time,theta,Mair,O2,N2,H2O,concentration,dt_chem,emission_isoprene,emission_monoterpene,cond_sink)
  real(dp) :: time ! [s] Current time since start
  real(dp) :: dt_chem ! [s] Chemistry time step
  real(dp) :: emission_isoprene ! [molecules/cm3/s] Emission rate of isoprene
  real(dp) :: emission_monoterpene ! [molecules/cm3/s] Emission rate of monoterpene
  real(dp), dimension(nz) :: theta ! [K], potential temperature
  real(dp), dimension(nz) :: temp ! [K], air temperature
  real(dp), dimension(nz) :: Mair ! [molecules/cm3] Air molecules concentration at each level.
  real(dp), dimension(nz) :: O2 ! [molecules/cm3] Oxygen concentration at each level
  real(dp), dimension(nz) :: N2 ! [molecules/cm3] Nitrogen concentration at each level
  real(dp), dimension(nz) :: H2O ! [molecules/cm3] Water concentration at each level
  real(dp), dimension(nz,neq) :: concentration ! [molecules/cm3] Chemical species concentrations (row) for each level (column)
  real(dp) :: exp_coszen
  real(dp), dimension(nz) :: emi_iso_height_dep ! Height-dependent emission of isoprene
  real(dp), dimension(nz) :: emi_mono_height_dep ! Height-dependent emission of monoterpene
  real(dp), dimension(nz,nr_cond) :: cond_sink

  temp = theta - (grav/Cp)*hh
  !write(*,*) "Temperature at level 2 is:", temp(2)
  pres = barometric_law(p00, temp, hh)
  exp_coszen = get_exp_coszen(time, daynumber, latitude)
  !write(*,*) exp_coszen

  ! Determine the concentration of air for each height
  do i = 1, nz
    Mair(i) = pres(i)*NA / (Rgas*temp(i)) * 1d-6
  end do
  ! Set concentrations under the assumption that the mixing ratio
  ! is homogenous throughout the air parcel
  ! Set initial concentrations under the assumption that the mixing ratio
  ! is homogenous throughout the air parcel
  do i = 1, nz
    O2(i) = 0.21_dp * Mair(i) ! Concentrations of O2 at the different levels
    N2(i) = 0.78_dp * Mair(i) ! Concentrations of N2 at the different levels
    !H2O(i) = 1.0D16 * (Mair(i)/Mair(1)) ! Water vapour concentratoin at the different levels
    H2O(i) = 1.0D16
    concentration(i,1) = 24.0_dp * Mair(i) * ppb ! O3 concentration
    concentration(i,5) = 0.2_dp * Mair(i) * ppb ! NO2 concentration
    concentration(i,6) = 0.07_dp * Mair(i) * ppb ! NO concentration
    concentration(i,9) = 100.0_dp * Mair(i) * ppb ! CO concentration
    concentration(i,11) = 1759.0_dp * Mair(i) * ppb ! CH4 concentration
    concentration(i,20) = 0.5_dp * Mair(i) * ppb ! SO2 concentration
  end do

  ! Constructing the vectors describing the height-dependence fo the emissions:
  do i = 1, nz
    if (i == 2) then
      emi_iso_height_dep(i) = emission_isoprene
      emi_mono_height_dep(i) = emission_monoterpene
    else 
      emi_iso_height_dep(i) = 0.0_dp
      emi_mono_height_dep(i) = 0.0_dp
    end if
  end do

  ! Doing the chemistry for all the heights
  do i = 2, nz-1
    call chemistry_step(concentration(i,:),time,time+dt_chem,O2(i),N2(i),Mair(i),H2O(i),temp(i),exp_coszen,emi_iso_height_dep(i),emi_mono_height_dep(i),cond_sink(i,:))
    !write(*,*) O2(i),N2(i),H2O(i)
  end do

end subroutine
!-----------------------------------------------------------------------------------------
! Calculate the radiation related quantities
!-----------------------------------------------------------------------------------------
real(dp) function get_exp_coszen(time,daynumber,latitude)
  real(dp), intent(in) :: time,latitude
  INTEGER, intent(in) :: daynumber
  real(dp) :: hourangle,zenith,coszen
  hourangle = get_hourangle(time)
  zenith = solar_zenith_angle(hourangle,daynumber,latitude)
  coszen = cos(zenith)
  IF (coszen > 0) THEN  ! sun is above horizon
     get_exp_coszen = exp(-0.575_dp/coszen)
  ELSE
     get_exp_coszen = 0.0_dp
  endIF
end function get_exp_coszen


real(dp) function get_hourangle(time)
  real(dp), intent(in) :: time
  real(dp), parameter :: one_day = 24*one_hour
  get_hourangle = modulo(time,one_day)/one_day * 2 * pi - pi
end function get_hourangle


real(dp) function solar_zenith_angle(hourangle,daynumber,latitude)
  ! http://en.wikipedia.org/wiki/Solar_elevation_angle
  ! http://en.wikipedia.org/wiki/Position_of_the_Sun
  INTEGER, intent(in) :: daynumber
  real(dp), intent(in) :: hourangle,latitude
  real(dp) :: declination,elevation
  real(dp), parameter :: to_rad = pi/180.0_dp

  declination = -23.44_dp * to_rad * cos(2 * pi * (daynumber + 10)/365.0_dp)
  elevation = cos(hourangle)*cos(declination)*cos(latitude) &
       + sin(declination)*sin(latitude)
  solar_zenith_angle = pi/2.0_dp - elevation
  ! Notes:
  ! - Not tested near equador or on the southern hemisphere.
  ! - solar_zenith_angle can be larger than pi/2, it just means
  !   the sun is below horizon.
  ! - solar_zenith_angle assumes time is in local solar time, which
  !   is usually not exactly true
end function solar_zenith_angle


!-----------------------------------------------------------------------------------------
! Other functions
!-----------------------------------------------------------------------------------------
function barometric_law(p00, tempK, h) result(p)
  real(dp), intent(in) :: p00, tempK(nz), h(nz)
  real(dp) :: p(nz)
  real(dp) :: dh(nz)

  dh(2:nz) = h(2:nz) - h(1:nz-1)

  p(1) = p00
  do i=2, nz
    p(i) = p(i-1)*exp(-mm_air*grav/(Rgas*(tempK(i-1)+tempK(i))/2.0d0)*dh(i))
  end do
end function barometric_law

function half_z(z1,z2) result(z_half)
  real(dp) :: z1
  real(dp) :: z2
  real(dp) :: z_half
  z_half = (z1 + z2) / 2
end function

end program main
