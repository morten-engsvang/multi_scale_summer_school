program chemistry_box

  use chemistry_mod

  implicit none

  integer, parameter :: dp = selected_real_kind(15,307)

  real(dp) :: conc(neq) ! Remember to make this 2D when changing to the 1D model
  real(dp) :: O2, N2, Mair, H2O, TEMP, pres
  real(dp) :: exp_coszen

  ! time management
  real(dp) :: time, time_start, time_end
  real(dp) :: dt                          ! time step, in seconds
  integer  :: time_in_ms                  ! time in milliseconds
  integer, parameter :: one_hour = 60*60  ! in seconds
  integer :: daynumber_start, daynumber

  ! Other parameters
  real(dp), parameter :: pi = 2*asin(1.0d0)
  real(dp), parameter :: Rgas = 8.3144598d0                   ! universal gas constant [J mol-1 K-1]
  real(dp), parameter :: NA   = 6.022140857d23                ! Avogadro's number [molec mol-1]
  real(dp), parameter :: ppb = 1e-9_dp

  ! latitude and longitude of SMEAR II
  real(dp), parameter :: latitude_deg = 61.8455d0
  real(dp), parameter :: latitude = latitude_deg * pi/180.0d0
  real(dp), parameter :: longitude_deg = 24.2873d0
  real(dp), parameter :: longitude = longitude_deg * pi/180.0d0

  ! Output
  character(len=255), parameter :: output_dir = './output'
  logical :: dir_exist

  ! Start program

  ! Atmospheric oxygen, N2 and H2O are kept constant:
  pres = 1.01325e5_dp                    ! [Pa], reference pressure at surface
  TEMP = 300.0d0                         ! [K], temperature
  Mair = pres*NA / (Rgas*temp) * 1d-6    ! Air molecules concentration [molecules/cm3]
  O2   = 0.21d0*Mair                     ! Oxygen concentration [molecules/cm3]
  N2   = 0.78d0*Mair                     ! Nitrogen concentration [molecules/cm3]
  H2O  = 1.0D16                          ! Water molecules [molecules/cm3]
  ! In the 1D model we also need to define the number concentration
  ! as being height dependent / pressure dependent

  ! initial state:
  conc     = 0.0d0
  conc(1)  = 24.0d0   * Mair * ppb     ! O3 concentration
  conc(5)  = 0.2d0    * Mair * ppb     ! NO2
  conc(6)  = 0.07d0   * Mair * ppb     ! NO
  conc(9)  = 100.0d0  * Mair * ppb     ! CO
  conc(11) = 1759.0d0 * Mair * ppb     ! CH4
  conc(13) = 2.2d0    * Mair * ppb     ! C5H8 (isoprene)
  conc(20) = 0.5d0    * Mair * ppb     ! SO2
  conc(23) = 2.2d0    * Mair * ppb     ! alpha-pinene (monoterpene)

  daynumber_start = 31+28+31+30+31+30+10 ! day is July 10th, 2011
  daynumber = daynumber_start

  time_start = 0.0d0
  time_end = 5 * 24 * one_hour
  dt = 10.0d0

  time = time_start

  exp_coszen = get_exp_coszen(time, daynumber, latitude, longitude)

  ! Create a new directory if it does not exist
  inquire(file=trim(adjustl(output_dir)), exist=dir_exist)
  if (.not. dir_exist) then
    ! This line may change for different operating systems
    call system('mkdir ' // trim(adjustl(output_dir)))
  end if

  ! Open the output files and write the initial values
  open(11,file=trim(adjustl(output_dir)) // "/concentrations.dat",status='replace',action='write')
  open(12,file=trim(adjustl(output_dir)) // "/exp_coszen.dat",status='replace',action='write')
  open(13,file=trim(adjustl(output_dir)) // "/time.dat",status='replace',action='write')
  write(11,*) conc
  write(12,*) exp_coszen
  write(13,*) time

  ! Start the main loop
  do while (time < time_end)

    exp_coszen = get_exp_coszen(time, daynumber, latitude, longitude)
    
    ! Calculate chemical reactions
    call chemistry_step(conc,time,time+dt,O2,N2,Mair,H2O,TEMP,exp_coszen)

    ! Advance time
    time = time + dt
    daynumber = daynumber_start + floor(time/(24 * one_hour))
    
    ! Write output every hour
    time_in_ms = floor(1000*time)
    if (modulo(time_in_ms, 1000*one_hour) == 0) then
       write(*,'(A8,F6.3,A6)') 'time = ', time/(24*one_hour), '  days'
       write(11,*) conc
       write(12,*) exp_coszen
       write(13,*) time/(24*one_hour)
    end if

  end do

  ! Close the files
  close(11)
  close(12)


contains


  real(dp) function get_hourangle(time)
     real(dp), intent(in) :: time
     real(dp), parameter :: one_day = 24*one_hour
     get_hourangle = modulo(time,one_day)/one_day * 2 * pi - pi
  end function get_hourangle


  real(dp) function solar_zenith_angle(hourangle,daynumber,latitude,longitude)
     ! http://en.wikipedia.org/wiki/Solar_elevation_angle
     ! http://en.wikipedia.org/wiki/Position_of_the_Sun
     integer, intent(in) :: daynumber
     real(dp), intent(in) :: hourangle,latitude,longitude
     real(dp) :: declination,elevation
     real(dp), parameter :: to_rad = pi/180.0d0

     declination = -23.44d0 * to_rad * cos(2 * pi * (daynumber + 10)/365.0d0)
     elevation = cos(hourangle)*cos(declination)*cos(latitude) &
                 + sin(declination)*sin(latitude)
     solar_zenith_angle = pi/2.0d0 - elevation
     ! Notes:
     ! - Not tested near equador or on the southern hemisphere.
     ! - solar_zenith_angle can be larger than pi/2, it just means that the sun is below horizon.
     ! - solar_zenith_angle assumes time is in local solar time, which is usually not exactly true

  end function solar_zenith_angle


  real(dp) function get_exp_coszen(time,daynumber,latitude,longitude)
     real(dp), intent(in) :: time,latitude,longitude
     integer, intent(in) :: daynumber
     real(dp) :: hourangle,zenith,coszen
     hourangle = get_hourangle(time)
     zenith = solar_zenith_angle(hourangle,daynumber,latitude,longitude)
     coszen = cos(zenith)
     if (coszen > 0) then  ! sun is above horizon
        get_exp_coszen = exp(-0.575d0/coszen)
     else
        get_exp_coszen = 0
     endif
  end function get_exp_coszen

end program chemistry_box
