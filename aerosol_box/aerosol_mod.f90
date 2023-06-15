MODULE aerosol_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: pi
PUBLIC :: dp, cond_vapour, diameter, particle_mass, particle_volume, particle_conc, &
          particle_density, nucleation_coef, molecular_mass, molar_mass, &
          molecular_volume, molecular_dia, mass_accomm, &
          PN, PM, PV
PUBLIC :: aerosol_init, nucleation, condensation, coagulation, dry_dep_velocity

!====================== Definition of variables =====================================================================!
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,300)
! so that numbers will be in 64bit floating point
! http://en.wikipedia.org/wiki/Double_precision_floating-point_format

REAL(dp), PARAMETER :: pi   = 2D0*ASIN(1D0)  ! constant pi
REAL(dp), PARAMETER :: ka   = 0.4D0          ! [-], von Karman constant, dimensionless
REAL(dp), PARAMETER :: g    = 9.81D0         ! [m s-2], gravitation const
REAL(dp), PARAMETER :: Rg   = 8.3145D0       ! Universal gas constant J mol^-1 K^-1
REAL(dp), PARAMETER :: Na   = 6.022D23       ! Avogadro's number 
REAL(dp), PARAMETER :: Mair = 28.96D-3       ! Mean molecular weight of air
REAL(dp), PARAMETER :: kb   = 1.381d-23      ! Boltzmann constant [(m2*kg)/(s2*K)]

INTEGER, PARAMETER ::  nr_bins = 100           ! Number of particle size bins
INTEGER, PARAMETER ::  nr_cond = 2             ! Number of condensable vapours
  
REAL(dp), DIMENSION(nr_bins) :: diameter       , &  ! Diameter of each size bin
                                particle_mass  , &  ! mass of one particle in each size bin
                                particle_conc  , &  ! number concentration in each size bin
                                particle_volume, &  ! volume concentration in each size bin 
                                coag_loss      , &  ! coagulation loss rate of particles in each size bin
                                v_dep               ! Dry deposition velocity of particles
     
REAL(dp), DIMENSION(nr_cond) :: molecular_mass  , &  ! molecular mass of the condensing vapours [kg/#]
                                molecular_volume, &  ! Molecule volume of condensable vapours [m^3]
                                molecular_dia   , &  ! Molecule diameter of condensable vapours [m]
                                molar_mass      , &  ! Molar mass of condensable vapours [kg/m^3]
                                cond_vapour          ! Concentration of condensable vapours [molec/m^3] 
                                      
REAL(dp), DIMENSION(nr_cond) :: cond_sink = 1.0d-3  ! Assumed initial condensation sink of vapours [s^-1]

REAL(dp) :: PN, PM, PV  ! Total particle number [# m-3] and mass concentration [kg m-3], particle volume [um^3 cm^-3]
     
REAL(dp) :: vd_SO2, vd_O3, vd_HNO3  ! [m s-1], dry deposition velocity of SO2, O3 & HNO3     

REAL(dp) :: particle_density, &  ! [kg]
            nucleation_coef , &  ! Nucleation coefficient 
            mass_accomm          ! mass accomodation coefficient 

REAL(dp) :: nucleation_rate  ! [# m-3 s-1]


CONTAINS


SUBROUTINE aerosol_init(diameter, particle_mass, particle_volume, particle_conc, &
                        particle_density, nucleation_coef, molecular_mass, molar_mass, &
                        molecular_volume, molecular_dia, mass_accomm)

  !======================================!
  ! Definition of variables
  !======================================!

  REAL(dp), DIMENSION(nr_bins), INTENT(OUT) :: diameter       , &  ! [m], diamter of each size bin
                                               particle_mass  , &  ! [kg], mass of one particle
                                               particle_volume, &  ! [m3], volume of one particle
                                               particle_conc       ! [# m-3], number concentration
  
  REAL(dp), DIMENSION(nr_cond), INTENT(OUT) :: molecular_mass  , & ! [kg], molecular mass of the condensing vapours
                                               molecular_volume, & ! [m3]
                                               molecular_dia   , & ! [m]
                                               molar_mass          ! [kg mol-1], molar mass of the condensing vapours
  
  REAL(dp), INTENT(OUT) :: nucleation_coef, &  ! [m3 molec-1], nucleation coefficient
                           mass_accomm         ! [-], mass accomodation coefficient
  
  REAL(dp), DIMENSION(nr_cond) :: density  ! [kg m-3], Bulk density of condensing vapours
  
  REAL(dp) :: particle_density  ! [kg m-3], particle density
  
  INTEGER :: i

  nucleation_coef = 1D-20
  mass_accomm = 1D0
  
  !===== Particle properties =====!

  ! Particle diameters between 2D-9 and 2.5D-6 m:
  diameter(1)=2D-9 
  DO i=2,nr_bins
    diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1))
  END DO
    
  particle_conc = 1D0 ! Assume an initial particle number concentration of 1 [# m-3]
  where((abs(diameter-2D-7)-MINVAL(abs(diameter-2D-7)))<1D-20)  particle_conc=2D8  ! add 200 [# cm-3] to 200 nm sized accumulation mode particles
  
  particle_density = 1.4D3                                       ! [kg m-3], assumed fixed particle density
  particle_volume = 1D0/6D0 * pi * diameter**3                   ! [m-3], single particle volume
  particle_mass=  1D0/6D0 * pi * diameter**3 * particle_density  ! [kg], single particle mass
  
  !===== Condensable vapor properties =====!

  density = (/1.84D3, 1.4D3/)                                ! density of sulphuric acid and ELVOC
  molar_mass = (/0.098D0, 0.3D0/)                            ! H2SO4 and ELVOC
  molecular_mass = molar_mass / Na                           ! molecular mass [kg]
  molecular_volume = molecular_mass / density                ! molecular volume [m-3]
  molecular_dia = (6D0 * molecular_volume / pi )**(1D0/3D0)  ! molecular diameter [m]

END SUBROUTINE aerosol_init
  
SUBROUTINE nucleation(nucleation_coef,conc_h2so4,particle_conc,dt_aerosol) ! (Add input and output variables here)
  REAL(dp) :: nucleation_coef  ! [m3 molec-1], nucleation coefficient
  REAL(dp), DIMENSION(nr_bins) :: particle_conc ! [# m-3], number concentration
  REAL(dp) :: conc_h2so4 ! [molecules/cm^3] Concentration of sulphuric acid.
  REAL(dp) :: J_2nm ! [molecules^-3 * s^-1] formation rate of 2 nm aerosol particles
  REAL(dp) :: dt_aerosol ! [s] Time step size
  ! Consider how kinetic H2SO4 nucleation influence the number concentrations of particles 
  ! in the fist size bin particle_conc(1) within one model time step
  J_2nm = nucleation_coef * conc_h2so4**2
  where((abs(diameter-2D-9)-MINVAL(abs(diameter-2D-9)))<1D-20) particle_conc = particle_conc + J_2nm*dt_aerosol ! Update the 2 nm size bin
  ! This code above is stolen from the initialization code, not sure exactly how this works.

END SUBROUTINE nucleation

SUBROUTINE condensation(timestep, temperature, pressure, mass_accomm, molecular_mass, &
                        molecular_volume, molar_mass, molecular_dia, particle_mass, particle_volume, &
                        particle_conc, diameter, cond_vapour) ! Add more variables if you need it
  
  REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: diameter, particle_mass
  REAL(dp), DIMENSION(nr_cond), INTENT(IN) :: molecular_mass, molecular_dia, &
                                              molecular_volume, molar_mass
  REAL(dp), INTENT(IN) :: timestep, temperature, pressure, mass_accomm
  
  REAL(dp), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc
  
  REAL(dp), DIMENSION(2), INTENT(IN) :: cond_vapour  ! [molec m-3], condensing vapour concentrations, which is H2SO4 and organics (ELVOC)
  
  REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: particle_volume
  
  REAL(dp), DIMENSION(nr_bins) :: slip_correction, diffusivity, speed_p, &
                                  particle_conc_new
  
  REAL(dp), DIMENSION(nr_cond) :: diffusivity_gas, speed_gas
  
  REAL(dp) :: dyn_visc, l_gas, dens_air

  REAL(dp), DIMENSION(nr_cond) :: fuchs_beta ! Fuchs Sutugin correction for the different vapours

  REAL(dp), DIMENSION(nr_cond) :: knudsen ! Knudsen number for the different vapours

  REAL(dp), DIMENSION(nr_cond) :: CR ! Collision rates for the different vapours

  REAL(dp), DIMENSION(nr_cond) :: lambda ! Mean free path of the vapours

  REAL(dp), DIMENSION(nr_bins) :: particle_volume_new = 0.0_dp ! Vector to store the increase in particle volume

  REAL(dp) :: x_1, x_2 ! Variables to store fractions for the full-stationary structure method
  
  INTEGER :: j, i
  
  ! Add more variabels as you need it...
  
  dyn_visc = 1.8D-5*(temperature/298D0)**0.85D0  ! dynamic viscosity of air
  dens_air=Mair*pressure/(Rg*temperature)        ! Air density
  l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature))) ! Gas mean free path in air (m)
  
  slip_correction = 1D0+(2D0*l_gas/(diameter))*&
  (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter))) ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34) 
  
  diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)   ! Diffusivity for the different particle sizes m^2/s
  speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                     ! speed of particles (m/s)
  
  diffusivity_gas=5D0/(16D0*Na*molecular_dia**2D0*dens_air)*&
  SQRT(Rg*temperature*Mair/(2D0*pi)*(molar_mass+Mair)/molar_mass)            ! Diffusivity of condensable vapours (m^2 s^-1)
  
  ! Thermal velocity of vapour molecule
  speed_gas=SQRT(8D0*kb*temperature/(pi*molecular_mass)) ! speed of H2SO4 molecule
  
  ! For all the particle diameters:
  do j = 1, nr_bins
    do i = 1, nr_cond
      ! Calculate the Fuchs-Sutugin correction factor:
      ! First the mean free path
      lambda(i) = 3 * (diffusivity_gas(i) + diffusivity(j)) / (sqrt(speed_gas(i)**2 + speed_p(j)**2))
      ! Then Knudsen for the condensable vapour
      knudsen(i) = 2 * lambda(i) / (diameter(j) + molecular_dia(i))
      ! Then Fuchs correction
      fuchs_beta(i) = (0.75_dp * mass_accomm * (1 + knudsen(i)))&
      / (knudsen(i)**2 + knudsen(i) + 0.283_dp * knudsen(i) * mass_accomm + 0.75_dp * mass_accomm)
      ! Calculate the Collision rate (CR [m^3/2]) between gas molecules (H2SO4 and ELVOC) and the particles:
      CR(i) = 2 * pi * (diameter(j) + molecular_dia(i)) * (diffusivity(j) + diffusivity_gas(i)) * fuchs_beta(i)
    end do
    ! Calculate the new single particle volume after condensation (particle_volume_new):
    particle_volume_new(j) = particle_volume(j)&
    + timestep * CR(1) * cond_vapour(1) * molecular_volume(1)&
    + timestep * CR(2) * cond_vapour(2) * molecular_volume(2)
  end do
  
  
  ! Use the full-stationary method to divide the particles between the existing size bins (fixed diameter grid):
  particle_conc_new=0D0 ! Initialise a new vector with the new particle concentrations
  particle_conc_new(nr_bins)=particle_conc(nr_bins) ! To ensure that they don't grow out of the box ,
  ! Otherwise we're not going to preserve the number of particles
  
  DO j = 1,nr_bins-1
  ! Add equations that redistributes the particle number concentration 
  ! in size bin 1 to nr_bins-1 to the fixed volume (diameter) grid 
    x_1 = (particle_volume(j+1) - particle_volume_new(j)) / (particle_volume(j+1) - particle_volume(j)) ! Fraction to leave in bin j
    x_2 = 1.0_dp - x_1 ! fraction to pass onto bin j+1
    particle_conc_new(j) = particle_conc_new(j) + x_1 * particle_conc(j)
    particle_conc_new(j + 1) = particle_conc_new(j + 1) + x_2 * particle_conc(j) 
  END DO
  ! Update the particle concentration in the particle_conc vector:
  particle_conc=particle_conc_new
    
END SUBROUTINE condensation

SUBROUTINE coagulation(timestep, particle_conc, diameter, &
                       temperature,pressure,particle_mass) ! Add more variables if you need it
  
  REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: diameter
  REAL(dp), DIMENSION(nr_bins), INTENT(INOUT) :: particle_conc
  REAL(dp), INTENT(IN) :: timestep
  REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: particle_mass       ! mass of one particle                                 
  REAL(dp), INTENT(IN) :: temperature, pressure
  
  REAL(dp), DIMENSION(nr_bins,nr_bins) :: coagulation_coef        ! coagulation coefficients [m^3/s]
  
  REAL(dp), DIMENSION(nr_bins) :: slip_correction, diffusivity, dist, speed_p, &
                                  Beta_Fuchs, free_path_p
  
  REAL(dp) :: dyn_visc, &  ! dynamic viscosity, kg/(m*s)
              l_gas        ! Gas mean free path in air
  
  REAL(dp), DIMENSION(nr_bins) :: loss_self
  REAL(dp), DIMENSION(nr_bins) :: loss_other
  
  INTEGER  :: i,j
  
  ! The coagulation coefficient is calculated according to formula 13.56 in Seinfield and Pandis (2006), Page 603
  
  dyn_visc = 1.8D-5*(temperature/298.0d0)**0.85                                              ! Dynamic viscosity of air
  
  l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature)))                        ! Gas mean free path in air (m)
  
  slip_correction = 1D0+(2D0*l_gas/(diameter))*&
  (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))                                        ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)
  
  diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)                 ! Diffusivity for the different particle sizes m^2/s
  
  speed_p = SQRT(8D0*kb*temperature/(pi*particle_mass))                                   ! Speed of particles (m/s)
  
  free_path_p = 8D0*diffusivity/(pi*speed_p)                                              ! Particle mean free path (m)
  
  dist = (1D0/(3D0*diameter*free_path_p))*((diameter+free_path_p)**3D0 &
  -(diameter**2D0+free_path_p**2D0)**(3D0/2D0))-diameter                    ! mean distance from the center of a sphere reached by particles leaving the sphere's surface (m)
  
  DO i = 1,nr_bins
     Beta_Fuchs = 1D0/((diameter+diameter(i))/(diameter+diameter(i)+&
     2D0*(dist**2D0+dist(i)**2D0)**0.5D0)+8D0*(diffusivity+diffusivity(i))/&
     (((speed_p**2D0+speed_p(i)**2D0)**0.5D0)*(diameter+diameter(i))))                    ! Fuchs correction factor from Seinfeld and Pandis, 2006, p. 600
  
     coagulation_coef(i,:) = 2D0*pi*Beta_Fuchs*(diameter*diffusivity(i)+&
     diameter*diffusivity+diameter(i)*diffusivity+diameter(i)*diffusivity(i))             ! coagulation rates between two particles of all size combinations  (m^3/s)    
  END DO
  
  ! Write equations that considers how the particle number concentration in each size bin 
  !(particle_conc) is influenced by the coagulation sink (loss of smaller particles when
  ! they collide with larger ones)

  ! You can first calculate the loss (loss1) do to self-coagulation between particles in the same size bin
  ! and then calculate the loss (loss2) due to coagulation with larger particles
  ! Then add the two loss terms together loss = loss1 + loss2 
  loss_self = 0.0_dp
  loss_other = 0.0_dp
  DO i = 1, nr_bins-1
    ! First the self-coagulation for the i'th particle size
    ! Except for the largest particles size, which is handled
    ! separately
    loss_self(i) = - coagulation_coef(i,i) * particle_conc(i)**2
    DO j = i + 1, nr_bins
      ! Then the loss to all the particles larger than the i'th particle size
      loss_other(i) = loss_other(i) - coagulation_coef(i,j) * particle_conc(i) * particle_conc(j)
    END DO
  END DO
  ! Handling of the largest particle size:
  loss_self(nr_bins) = - coagulation_coef(nr_bins,nr_bins) * particle_conc(nr_bins)**2
  loss_other(nr_bins) = 0.0_dp
  ! Updating the concentrations:
  DO i = 1, nr_bins
    particle_conc(i) = particle_conc(i) + timestep * (loss_self(i) + loss_other(i))
  END DO
END SUBROUTINE coagulation

  
SUBROUTINE dry_dep_velocity(diameter,particle_density,temperature,pressure,DSWF, & 
                            Richards_nr10m,wind_speed10m,particle_conc,vapour_conc,timestep,mixing_height) ! Add more variables if you need it
   
  REAL(dp), DIMENSION(nr_bins) :: particle_conc
  REAL(dp), DIMENSION(2) :: vapour_conc
  REAL(dp), DIMENSION(nr_bins), INTENT(IN) :: diameter
  REAL(dp), INTENT(IN) :: timestep
  REAL(dp), INTENT(IN) :: mixing_height
  
  REAL(dp), INTENT(IN) :: temperature, pressure, Richards_nr10m, DSWF, &
  wind_speed10m, particle_density
        
  REAL(dp) :: z0m, r_coll, a_landuse, j_landuse, v_kinematic,dyn_visc,l_gas,Pr,beta,&
  gam,zr,u_friction,dens_air, L_Ob, raO3, raSO2, raHNO3, raisoprene, raapinene
  
  ! Specific parameters for the surface resistance of gases:
  REAL(dp) :: rj,rlu,rac, &
              DiffusivityH2O, D_ratio_SO2, D_ratio_O3, D_ratio_HNO3, D_ratio_isoprene, D_ratio_apinene, &
              DiffusivitySO2, DiffusivityO3, DiffusivityHNO3, Diffusivityisoprene, Diffusivityapinene,&
              z_roughSO2, z_roughO3, z_roughHNO3, z_roughisoprene, z_roughapinene, &
              ScSO2, ScO3, ScHNO3, Scisoprene, Scapinene, &
              rbSO2, rbO3, rbHNO3, rbisoprene, rbapinene, &
              H_effSO2, H_effO3, H_effHNO3, H_effisoprene, H_effapinene, &
              f0_SO2, f0_O3, f0_HNO3, f0_isoprene, f0_apinene, &
              rclSO2, rclO3, rgsSO2, rgsO3

  ! Bulk canopy
  REAL(dp) :: rst_h2o ! Bulk canopy stomatal resistance (for H2O)
  REAL(dp), DIMENSION(5) :: r_sm_i ! Combination of bulk canopy stomatal and mesophyll resistances
  ! Outer canopy and lower canopy
  REAL(dp), DIMENSION(5) :: r_lu_i ! Outer canopy surface resistances
  REAL(dp) :: r_dc ! Resitance in the lower canopy to transfer by bouyant convection
  REAL(dp), DIMENSION(5) :: r_cl_i ! Resistance in the lower canopy to uptake
  ! Soil uptake
  ! rac is already defined and independent of the gas, transfer resistance to the soil
  REAL(dp), DIMENSION(5) :: r_gs_i ! Resitance to uptake by the soil
  REAL(dp), DIMENSION(5) :: r_c_i ! Combined canopy resistances for gasses
  
  REAL(dp), DIMENSION(5) :: deposition_gas ! Gas deposition velocity

  REAL(dp), DIMENSION(nr_bins) :: v_sedimentation ! Sedimentation velocity
  REAL(dp), DIMENSION(nr_bins) :: slip_correction ! Cunninghams slip correction factor
  REAL(dp), DIMENSION(nr_bins) :: diffusivity ! Diffusivity of different particle sizes
  REAL(dp), DIMENSION(nr_bins) :: z0_particle ! Surface roughness length scale for molecules and particles
  REAL(dp), DIMENSION(nr_bins) :: ra_particle ! Aerodynamic resistance for the particles
  REAL(dp), DIMENSION(nr_bins) :: rb_particle ! Quasi-laminar resistance for the particles-
  REAL(dp) :: u_star ! Friction velocity
  REAL(dp), DIMENSION(nr_bins) :: Scp ! Dimension-less Schmidt number for particles
  REAL(dp), DIMENSION(nr_bins) :: St ! Stokes number for vegetation
  REAL(dp), DIMENSION(nr_bins) :: R1 ! Sticking factor
  REAL(dp), DIMENSION(nr_bins) :: deposition_particle ! Particle deposition velocity

  INTEGER :: i

  dens_air = Mair*pressure/(Rg*temperature)    ! Air density (kg/m^3)
  dyn_visc = 1.8D-5*(temperature/298.)**0.85   ! dynamic viscosity of air (kg/(m*s))
  v_kinematic = dyn_visc/dens_air              ! kinematic viscosity of air (m^2/s)
  
  ! Code stolen from coagulation code in order to calculate Cunninghams slip correction factor
  l_gas=2D0*dyn_visc/(pressure*SQRT(8D0*Mair/(pi*Rg*temperature)))                        ! Gas mean free path in air (m)
  slip_correction = 1D0+(2D0*l_gas/(diameter))*&
  (1.257D0+0.4D0*exp(-1.1D0/(2D0*l_gas/diameter)))                                        ! Cunninghams slip correction factor (Seinfeld and Pandis eq 9.34)

  ! Particle diffusivity for later :)
  diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)                 ! Diffusivity for the different particle sizes m^2/s

  zr=10D0                 ! Reference height [m]
  z0m = 0.9D0             ! Surface roughness length for momentum evergreen, needleleaf trees (m)     
  u_friction=ka*wind_speed10m/(log(zr/z0m))  ! Friction velocity (Eq. 16.67 from Seinfeld and Pandis, 2006)

  ! Land use category paramaters from Seinfeld and Pandis, 2006 Table 19.2: 
  r_coll = 2D-3 ! radius of collector evergreen, needleleaf trees

  ! coefficients based on land use categories (evergreen, needleleaf trees)
  a_landuse = 1D0
  j_landuse = 0.56D0

  Pr = 0.95D0   ! Turbulent Prandtl number (when ka = 0.4 (Hogstrom, 1988))
  beta = 7.8D0  ! When ka = 0.4 (Hogstrom, 1988)
  gam = 11.6D0  ! When ka = 0.4 (Hogstrom, 1988)

  ! Calculate the particle sedimentation velocity:
  DO i = 1, nr_bins
    v_sedimentation(i) = (diameter(i) ** 2) * (particle_density - dens_air) * g * slip_correction(i) / (18.0_dp * dyn_visc)
  END DO

  ! Particle diffusivity for later :)
  diffusivity = slip_correction*kb*temperature/(3D0*pi*dyn_visc*diameter)! Diffusivity for the different particle sizes m^2/s
  ! Calculation of aerodynamic resistance for particles for:
  ! stable boundary layer (Ri>1D-6)
  ! neutral boundary layer (abs(Ri)<1D-6
  ! unstable boundary layer Ri<-1D-6
  ! Precalculation of z0 and u_star
  u_star = ka * wind_speed10m / (log(zr / z0m))
  z0_particle = diffusivity / (ka * u_star)
  ! Check if neutral
  IF (abs(Richards_nr10m) < 1D-6) THEN
    do i = 1, nr_bins
      ra_particle(i) = Pr * log(zr / z0_particle(i)) / (ka * u_star)
    end do
  ! Check if stable:
  ELSE IF (Richards_nr10m > 0) THEN
    do i = 1, nr_bins
      ra_particle(i) = (Pr * log(zr / z0_particle(i)) + beta / ( (zr / Richards_nr10m) * (zr - z0_particle(i)))) / (ka * u_star)
    end do
  ! Else it has to be unstable:
  ELSE
    L_Ob=zr/Richards_nr10m ! Monin-Obukhov length scale
    do i = 1, nr_bins
      ra_particle(i) = Pr * log(phi(zr,L_Ob,gam) / phi(z0_particle(i),L_Ob,gam)) / (ka * u_star)
    end do
  END IF

  ! Calculate the quasi-laminar resistance (rb) for particles:
  do i = 1, nr_bins
    Scp(i) = v_kinematic / diffusivity(i)
    St(i) = v_sedimentation(i) * u_star / (g * r_coll)
    R1(i) = exp(-sqrt(St(i)))
  end do

  do i = 1, nr_bins
    rb_particle(i) = (3 * u_star * R1(i) * (Scp(i)**(-j_landuse) + (St(i)&
    / (a_landuse + St(i)))**2 + 0.5_dp * (diameter(i) / r_coll)**2))**(-1.0_dp)
  end do

  ! Calculate the dry deposition velocity for particles:
  do i = 1, nr_bins
    deposition_particle(i) = 1.0_dp / (ra_particle(i) + rb_particle(i)&
    + ra_particle(i) * rb_particle(i) * v_sedimentation(i)) + v_sedimentation(i)
  end do
  write(*,*) ra_particle(1), rb_particle(1), v_sedimentation(1), deposition_particle(1)
  ! Calculate the dry deposition velocity for O3, SO2, HNO3, isoprene and a-pinene: 

  ! Resistance components used when calculating the surface resistance for gases, 
  ! table 19.3 Seinfeld and Pandis, 2006: 

  ! The minimum, bulk canopy stomatal resistance for water vapor:
  rj = 130D0 ! (s/m) Summer, evergreen, needleleaf
  
  ! The resistance of the outer surfaces in the upper canopy
  rlu = 2000D0 ! (s/m) Summer, evergreen, needleleaf
  
  ! transfer resistance on the ground (that depends only on canopy height)
  rac = 2000D0 ! (s/m) Summer, evergreen, needleleaf
  
  ! resistance for uptake by soil, leaf litter, and so on at the ground SO2
  rgsSO2 = 500D0 ! (s/m) Summer, evergreen, needleleaf
  
  ! restistance for uptake by soil, leaf litter, and so on at the ground, O3
  rgsO3 = 200D0 ! (s/m) Summer, evergreen, needleleaf
  
  ! resistance for uptake by leaves,twigs, and other exposed surfaces, SO2
  rclSO2 = 2000D0 ! (s/m) Summer, evergreen, needleleaf  
  
  ! resistance for uptake by leaves, twigs, and other exposed surfaces, O3
  rclO3 = 1000D0 ! (s/m) Summer, evergreen, needleleaf    
    
  ! Diffusion coefficients of selected gases
  DiffusivityH2O = 0.234D-4 ! Diffusion coefficient of water vapor in air (m^2/s), table 16.2 Seinfeld and Pandis
  
  ! ratio between diffusivity of water vapor and SO2, O3 or HNO3 from table 19.4 Seinfeld and Pandis
  D_ratio_SO2 = 1.89D0
  D_ratio_O3 = 1.63D0
  D_ratio_HNO3 = 1.87D0
  D_ratio_isoprene = 2.7D0 ! Estimated 
  D_ratio_apinene = 4D0  ! Estimated 
  
  DiffusivitySO2 = DiffusivityH2O/D_ratio_SO2    ! Diffusivity of SO2 (m^2/s)
  DiffusivityO3 = DiffusivityH2O/D_ratio_O3      ! Diffusivity of O3 (m^2/s)
  DiffusivityHNO3 = DiffusivityH2O/D_ratio_HNO3  ! Diffusivity of HNO3 (m^2/s)
  Diffusivityisoprene = DiffusivityH2O/D_ratio_isoprene  ! Diffusivity of isoprene (m^2/s)
  Diffusivityapinene = DiffusivityH2O/D_ratio_apinene  ! Diffusivity of apinene (m^2/s)
  
  ! Calculate the aerodynamic resistance for O3, SO2, HNO3, isoprene & a-pinene (ra) in similar way as
  ! for particles:
  ! Precalculation of surface roughness, u_star is the same as the previous value
  z_roughSO2 = DiffusivitySO2 / (ka * u_star)
  z_roughO3 = DiffusivityO3 / (ka * u_star)
  z_roughHNO3 = DiffusivityHNO3 / (ka * u_star)
  z_roughisoprene = Diffusivityisoprene / (ka * u_star)
  z_roughapinene = Diffusivityapinene / (ka * u_star)
  ! Check if neutral
  IF (abs(Richards_nr10m) < 1D-6) THEN
    raSO2 = Pr * log(zr / z_roughSO2)
    raO3 = Pr * log(zr / z_roughO3)
    raHNO3 = Pr * log(zr / z_roughHNO3)
    raisoprene = Pr * log(zr / z_roughisoprene)
    raapinene = Pr * log(zr / z_roughapinene)
  ! Check if stable:
  ELSE IF (Richards_nr10m > 0) THEN
    raSO2 = calcZ0(zr,Richards_nr10m,Pr,ka,u_star,z_roughSO2)
    raO3 = calcZ0(zr,Richards_nr10m,Pr,ka,u_star,z_roughO3)
    raHNO3 = calcZ0(zr,Richards_nr10m,Pr,ka,u_star,z_roughHNO3)
    raisoprene = calcZ0(zr,Richards_nr10m,Pr,ka,u_star,z_roughisoprene)
    raapinene = calcZ0(zr,Richards_nr10m,Pr,ka,u_star,z_roughapinene)
  ! Else it has to be unstable:
  ELSE
    L_Ob=zr/Richards_nr10m
    raSO2 = Pr * log(phi(zr,L_Ob,gam) / phi(z_roughSO2,L_Ob,gam)) / (ka * u_star)
    raO3 = Pr * log(phi(zr,L_Ob,gam) / phi(z_roughO3,L_Ob,gam)) / (ka * u_star)
    raHNO3 = Pr * log(phi(zr,L_Ob,gam) / phi(z_roughHNO3,L_Ob,gam)) / (ka * u_star)
    raisoprene = Pr * log(phi(zr,L_Ob,gam) / phi(z_roughisoprene,L_Ob,gam)) / (ka * u_star)
    raapinene = Pr * log(phi(zr,L_Ob,gam) / phi(z_roughapinene,L_Ob,gam)) / (ka * u_star)
  END IF

  ! Calculate the quasi-laminar resistance for O3, SO2, HNO3, isoprene & a-pinene (rb):
  rbSO2 = (5.0_dp / u_star) * (v_kinematic / DiffusivitySO2)**(2.0_dp/3.0_dp) 
  rbO3 = (5.0_dp / u_star) * (v_kinematic / DiffusivityO3)**(2.0_dp/3.0_dp)
  rbHNO3 = (5.0_dp / u_star) * (v_kinematic / DiffusivityHNO3)**(2.0_dp/3.0_dp)
  rbisoprene = (5.0_dp / u_star) * (v_kinematic / Diffusivityisoprene)**(2.0_dp/3.0_dp)
  rbapinene = (5.0_dp / u_star) * (v_kinematic / Diffusivityapinene)**(2.0_dp/3.0_dp) 
  ! Calculation of surface resistance for O3, SO2, HNO3, isoprene & a-pinene (rc)
  
  ! Effective Henry's lay const:
  H_effSO2 = 1D5   ! M atm^-1
  H_effO3 = 1D-2   ! M atm^-1
  H_effHNO3 = 1D14 ! M atm^-1
  H_effisoprene = 1.2D-2 ! M atm^-1
  H_effapinene = 3D-2 ! M atm^-1
  
  ! Noramlized reactivity, table 19.4 from Seinfeld and Pandis, 2006:
  f0_SO2 = 0D0
  f0_O3 = 1D0
  f0_HNO3 = 0D0
  f0_isoprene = 0D0
  f0_apinene = 0D0
  
  ! Calculate the bulk canopy stomatal resistance (rst)
  rst_h2o = rj * (1 + (200.0_dp / (DSWF + 0.1_dp))**2 * 400.0_dp / (temperature * (40.0_dp - temperature)))

  ! Calculate the combined stomatal and mesophyll resistance (rsm):
  r_sm_i(1) = rst_h2o * D_ratio_SO2 + 1.0_dp / (3.3D-4 * H_effSO2 + 100*f0_SO2)
  r_sm_i(2) = rst_h2o * D_ratio_O3 + 1.0_dp / (3.3D-4 * H_effO3 + 100*f0_O3)
  r_sm_i(3) = rst_h2o * D_ratio_HNO3 + 1.0_dp / (3.3D-4 * H_effHNO3 + 100*f0_HNO3)
  r_sm_i(4) = rst_h2o * D_ratio_isoprene + 1.0_dp / (3.3D-4 * H_effisoprene + 100*f0_isoprene)
  r_sm_i(5) = rst_h2o * D_ratio_apinene + 1.0_dp / (3.3D-4 * H_effapinene + 100*f0_apinene)

  ! Calculate the resistance of the outer surfaces in the upper canopy (rlu):
  r_lu_i(1) = rlu / (1D-5 * H_effSO2 + f0_SO2)
  r_lu_i(2) = rlu / (1D-5 * H_effO3 + f0_O3)
  r_lu_i(3) = rlu / (1D-5 * H_effHNO3 + f0_HNO3)
  r_lu_i(4) = rlu / (1D-5 * H_effisoprene + f0_isoprene)
  r_lu_i(5) = rlu / (1D-5 * H_effapinene + f0_apinene)

  ! Calculate the resistance to transfer by buoyant convection (rdc):
  r_dc = 100.0_dp * (1.0_dp + 1000.0_dp / (DSWF + 10.0_dp))

  ! Calculate the resistance of the exposed surfaces in the lower portions of 
  ! structures of the canopy (rcl): 
  r_cl_i(1) = ((1d-5 * H_effSO2 / rclSO2) + (f0_SO2 / rclO3))**(-1.0_dp)
  r_cl_i(2) = ((1d-5 * H_effO3 / rclSO2) + (f0_O3 / rclO3))**(-1.0_dp)
  r_cl_i(3) = ((1d-5 * H_effHNO3 / rclSO2) + (f0_HNO3 / rclO3))**(-1.0_dp)
  r_cl_i(4) = ((1d-5 * H_effisoprene / rclSO2) + (f0_isoprene / rclO3))**(-1.0_dp)
  r_cl_i(5) = ((1d-5 * H_effapinene / rclSO2) + (f0_apinene / rclO3))**(-1.0_dp)
          
  ! Calculate the resistance of the exposed surfaces on the groud 
  !(soil,leaf litter, ground) (rgs):
  r_gs_i(1) = ((1d-5 * H_effSO2 / rgsSO2) + (f0_SO2 / rgsO3))**(-1.0_dp)
  r_gs_i(2) = ((1d-5 * H_effO3 / rgsSO2) + (f0_O3 / rgsO3))**(-1.0_dp)
  r_gs_i(3) = ((1d-5 * H_effHNO3 / rgsSO2) + (f0_HNO3 / rgsO3))**(-1.0_dp)
  r_gs_i(4) = ((1d-5 * H_effisoprene / rgsSO2) + (f0_isoprene / rgsO3))**(-1.0_dp)
  r_gs_i(5) = ((1d-5 * H_effapinene / rgsSO2) + (f0_apinene / rgsO3))**(-1.0_dp)
  
  ! Combine all resistances in order to get the total surface resistance 
  ! for O3, SO2, HNO3, isoprene and a-pinene (rc):
  do i = 1, 5
    r_c_i(i) = (1.0_dp / r_sm_i(i) + 1.0_dp / r_lu_i(i) + 1.0_dp / (r_dc + r_cl_i(i)) + 1.0_dp / (rac + r_gs_i(i)))**(-1.0_dp)
  end do

  ! Finally calculate the dry deposition velocity of SO2, O3, HNO3, isoprene and a-pinene:
  deposition_gas(1) = 1.0_dp / (raSO2 + rbSO2 + r_c_i(1))
  deposition_gas(2) = 1.0_dp / (raO3 + rbO3 + r_c_i(2))
  deposition_gas(3) = 1.0_dp / (raHNO3 + rbHNO3 + r_c_i(3))
  deposition_gas(4) = 1.0_dp / (raisoprene + rbisoprene + r_c_i(4))
  deposition_gas(5) = 1.0_dp / (raapinene + rbapinene + r_c_i(5))

  ! Finally do the deposition
  do i = 1, nr_bins
    particle_conc(i) = particle_conc(i)*exp(-deposition_particle(i) * timestep / mixing_height)
  end do

END SUBROUTINE dry_dep_velocity

function phi(input,L_MO,gam) result(retval)
  real(dp), intent(in) :: input
  real(dp), intent(in) :: L_MO
  real(dp), intent(in) :: gam
  real(dp) :: retval
  retval = (sqrt(1 - gam * input / L_MO) - 1) / (sqrt(1 - gam * input / L_MO) + 1)
end function phi

function calcZ0(zr,richards,pr,ka,u_star,z0) result(retval)
  real(dp), intent(in) :: zr
  real(dp), intent(in) :: richards
  real(dp), intent(in) :: pr
  real(dp), intent(in) :: ka
  real(dp), intent(in) :: u_star
  real(dp), intent(in) :: z0
  real(dp) :: retval
  real(dp), parameter :: beta = 7.8
  retval = (pr * log(zr / z0) + beta / ( (zr / richards) * (zr - z0))) / (ka * u_star)
end function


END MODULE aerosol_mod
