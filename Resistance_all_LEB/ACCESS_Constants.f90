!============================================================================================================
!    Program:    Simplified ACCESS Model
!                Simplified (Dry Deposition only) Atmospheric Chemistry and Canopy Exchange Simulation System
!
!    Developed by : Dr. Beiming Tang based on Dr. Rick Saylor's ACCESS model version 3.1.0
!============================================================================================================

module ACCESS_Constants
    use, intrinsic :: iso_fortran_env, only: dp=>real64, i4=>int32
    implicit none

!===========================================
!part 1.define constants used in calculation
!===========================================

    integer(kind=i4), parameter             :: ninteg = 31        !# of integrated species, depend on chemical mechanism used in ACCESS
    integer(kind=i4), parameter             :: ntotal = 31        !# of total species in ACCESS based on chosen chemical mechanim, including transported species
    integer(kind=i4), parameter             :: npts   = 8        !# of vertical layers
    integer(kind=i4), parameter             :: ncnpy  = 8        !# of vertical layers


!==================================================
!part 2.define global variables used in calculation
!==================================================

!    real(kind=dp), parameter                :: tsoilk = 298      !soil/litter temperature(K)  
    real(kind=dp), parameter                :: stheta = 0.01     !volumetric soil water content (m^3/m^3)
    real(kind=dp), parameter                :: sattheta = 0.1    !saturation volumetric soil water content (m^3/m^3)
    real(kind=dp), parameter                :: dsoil = 4         !depth of topsoil (cm)
    real(kind=dp), parameter                :: rtheta = 0.05     !residual volumetic soil water content (m^3/m^3)
    real(kind=dp), parameter                :: sbcoef = 0.2      !clapp and hornberger exponent
    real(kind=dp), parameter                :: pi = 3.1415926535
    real(kind=dp), parameter                :: kd = 0.6
    real(kind=dp), parameter                :: x = 1.0
    real(kind=dp), parameter                :: tsoilk = 294.0    !soil temperature,(K)

end module ACCESS_Constants
