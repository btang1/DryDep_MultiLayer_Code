!=======================================================================================================
!    Program: Simplied ACCESS model
!             Simplied (Dry deposition only) Atmospheric Chemistry and Canopy Exchange Simulation System
!
!    developed by : Dr. Beiming Tang based on Dr. Rick D. Saylor's ACCESS model version 3.1.0
!=======================================================================================================

module DryDep
    use ACCESS_Constants
    use ACCESS_Modules

    implicit none

    public  GetDryDepExCoeffs, GetSoilDryDepExCoeffs

contains

!============================================================================
!PART 1. CALCULATE LEAF-SCALE DRY DEPOSITION VELOCITIES & COMPENSATION POINTS
!============================================================================

subroutine GetDryDepExCoeffs(rb,rc,rm,rs)
    integer(kind=i4)                                       :: i                       !i is layer
    integer(kind=i4)                                       :: l                       !l is species
    real(kind = dp)                                        :: mdiffl                  !molecular diffusivity of species l in air (cm^2/s)
    real(kind = dp)                                        :: relhumi                 !relative humidity (%)
    real(kind = dp)                                        :: hstarl                  !effective Henry's Law coefficient (M/atm)
    real(kind=dp)                                          :: f01                     !reactivity parameter (0-1) !Wesely's reactivity parameter (dimensionless)
    real(kind = dp)                                        :: srad                    !Solar radiation at canopy top (W/m^2)
    character(len=10)                                      :: rs_select               !Selection of stomatal resistance algorithm
    real(kind=dp)                                          :: hc                      !hc is canopy height (cm)
    real(kind=dp),dimension(npts)                          :: ppfd                    !vertical profile of photosynthetic photon flux (umol/m^2-s) at current simulation time
    real(kind=dp),dimension(npts)                          :: lai                     !layer leaf area index of canopy (cm^2*leaf/cm^2)
    real(kind=dp),dimension(npts)                          :: clai
    real(kind=dp),dimension(npts)                          :: tk                      !vertical profile of temperature at current simulation time (K)
    real(kind=dp),dimension(npts)                          :: pmb                     !vertical profile of pressure at current simulation time (mb)
    real(kind=dp),dimension(npts)                          :: qh                      !vertical profile of specific humidity at current simulation time (g/kg)
    real(kind=dp),dimension(npts)                          :: ubar                    !vertical profile of mean wind speed at current simulation time (cm/s)
    real(kind = dp),dimension(npts,ninteg),intent(inout)                 :: rb                      !leaf boundary resistance (s/cm)
    real(kind = dp),dimension(npts,ninteg),intent(inout)                 :: rc                      !cuticular resistance (s/cm)
    real(kind = dp),dimension(npts,ninteg),intent(inout)                 :: rm                      !mesophyll resistance (s/cm)
    real(kind = dp),dimension(npts,ninteg),intent(inout)                 :: rs                      !stomatal resistance (s/cm)
    real(kind=dp),dimension(npts,ninteg)    :: vd                      !dry deposition exchange coefficient (cm/s)
    real(kind = dp),dimension(npts,ninteg)                 :: gp                      !dry deposition compensation point (mole/cm^3)
    real(kind = dp)                                        :: rnum
    real(kind = dp)                                        :: rden
    real(kind = dp)                                        :: rlx
    real(kind = dp)                                        :: vdlx
    real(kind=dp),dimension(npts)                          :: fshd                    !shaded fraction of canopy
    real(kind=dp),dimension(npts)                          :: fsun                    !sunlit fraction of canopy
    real(kind=dp)                                          :: gs_sun_i                !leaf stomatal conductance in sublit fraction (mol/m^2-s)
    real(kind=dp)                                          :: gs_shd_i 
    real(kind = dp)                                        :: tsoilk                   !soil temperature over simulation (K)
    real(kind = dp)                                        :: kd                      !diffuse radiation extinction coefficient 
    real(kind = dp)                                        :: zadeg
    real(kind = dp)                                        :: zarad
    real(kind = dp)                                        :: tin
    real(kind = dp)                                        :: slat
    real(kind = dp)                                        :: slon
    real(kind = dp)                                        :: pmbzref
    real(kind = dp)                                        :: ppfdzref
    real(kind = dp)                                        :: sradzref
    real(kind = dp)                                        :: ppfd_direct
    real(kind = dp)                                        :: ppfd_diffus
    real(kind = dp)                                        :: nir_direct
    real(kind = dp)                                        :: nir_diffus
    real(kind = dp)                                        :: tkzref
    real(kind = dp)                                        :: eatm
!    real(kind = dp),parameter                              :: x=1.0
!    real(kind = dp),parameter                              :: kd=0.6
!    real(kind = dp),parameter                              :: tsoilk=294.0
    real(kind=dp)                                          :: lwdnzref
    real(kind=dp),dimension(npts)                          :: ppfd_sun
    real(kind=dp),dimension(npts)                          :: ppfd_shd
    real(kind=dp),dimension(npts)                          :: rabs_sun
    real(kind=dp),dimension(npts)                          :: rabs_shd
    real(kind = dp)                                        :: ubarms
    real(kind=dp)                                          :: gb_i
    real(kind=dp)                                          :: anet_sun_i
    real(kind=dp)                                          :: anet_shd_i
    real(kind=dp)                                          :: tl_sun_i
    real(kind=dp)                                          :: tl_shd_i
    real(kind=dp)                                          :: gl_sun_i
    real(kind=dp)                                          :: gl_shd_i
    real(kind=dp)                                          :: rs_wgt_i
    
    !nighttime 1am, get input data from 5m to 40m, every 5m interval
!    pmb    = (/1012.66     , 1012.07     , 1011.48     , 1010.891    , 1010.292    , 1009.702    , 1009.111    , 1008.522    /)              !calc use barometric formula
!    clai    = (/0.444595E+01, 0.395672E+01, 0.331488E+01, 0.211422E+01, 0.148396E+01, 0.392690E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
!    tk     = (/0.293843E+03, 0.293680E+03, 0.293552E+03, 0.293149E+03, 0.292633E+03, 0.292133E+03, 0.292153E+03, 0.292183E+03/)              !data from CANACC
!    ppfd   = (/0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC
!    qh     = (/5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         /)              !calc use RH=30% T=20C
!    ubar   = (/0.540616E+01, 0.688076E+01, 0.898461E+01, 0.121621E+02, 0.174182E+02, 0.278164E+02, 0.861199E+02, 0.861199E+02/)              !data from CANACC

    !daytime 11am, get input data from 5m to 40m, every 5m interval
    pmb    = (/1012.665, 1012.09     , 1011.515    , 1010.945    , 1010.37     , 1009.783    , 1009.165    , 1008.567    /)              !calc use barometric formula
    lai    = (/0.49328E+00,0.64184E+00, 1.20066E+00, 0.63026E+00,   1.09027E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00/)              !data calculated from clai
    tk     = (/2.96E+02, 0.298247E+03, 0.299576E+03, 0.300502E+03, 0.300629E+03, 0.299304E+03, 0.296860E+03, 0.295961E+03/)              !data from CANACC
    ppfd   = (/1.79E+02, 0.230477E+03, 0.317002E+03, 0.558829E+03, 0.743638E+03, 0.119511E+04, 0.140834E+04, 0.225315E+04/)              !data from CANACC
    qh     = (/5.2     , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         , 5.2         /)              !calc use RH=30% T=20C
    ubar   = (/1.77E+01, 0.225845E+02, 0.294898E+02, 0.399192E+02, 0.571712E+02, 0.913008E+02, 0.282668E+03, 0.282688E+03/)              !data from CANACC
    clai   = (/4.45E+00, 0.395672E+01, 0.331488E+01, 0.211422E+01, 0.148396E+01, 0.393690E+00, 0.000000E+00, 0.000000E+00/)              !data from CANACC

    tin = (16.0*365.0+4.0)*24.0*3600.0+(31.0+29.0+31.0+30.0+31.0+30.0+27.0)*24.0*3600.0+(11.0+5.0)*3600.0          !this is 07/28/2016 11am, seconds since Jan 1,2000 NOON england
    slat = 35.0597  !this is observation site
    slon = -83.4306 !this is observation site

    rs_select = 'gs_medlyn'                                             !To Do: make this selection via input file
    srad = ppfd(int(hc)+1)/4.57_dp                                      !cannopy top PPFD converted from (umol/m^2-s) to (W/m^2)
    do l = 1,ninteg                                                     !ninteg is integrated # of species in ACCESS model
        do i = 1,npts                                                   !npts is # of vertical layers
            if (lai(i) > 0.0) then                                      !within canopy
                mdiffl  = MolecDiff(l,tk(i),pmb(i))                     !calculate molecular diffusivity (cm^2/s)
                relhumi = RelativeHumidity(tk(i),pmb(i),qh(i))          !calculate relative humidity at layer i
                rb(i,l) = rbl(mdiffl, ubar(i))                          !leaf boundary layer resistance (s/cm)
                
                hstarl  = EffHenrysLawCoeff(l)
                f01     = ReactivityParam(l)
                rc(i,l) = rcl(hstarl, f01)                              !leaf cuticular resistance (s/cm)
                rm(i,l) = rml(hstarl, f01)                              !leaf mesophyll resistance (s/cm)
                select case (rs_select)
                    case('zhang_df')                                    !Zhang et al., (2002,2003) with hardcodes deciduous forest parameters
                        rs(i,l) = rs_zhang_df(mdiffl,tk(i),pmb(i),ppfd(i),srad,relhumi)
                    case('gs_medlyn')                                   !Medlyn et al., (2011)
                        !to get ppfd_direct,ppfd_diffus,nir_direct,nir_diffus 
                        zadeg = SolarZenithAngle(tin/(24.0*3600.0),slat,slon)
                        zarad = zadeg*pi/180.0
                        pmbzref = 1000.0          !air pressure at measurement height (mb)
                        ppfdzref = 2000.0         !PPFD at measurement height (umol/m^2-s)
                        sradzref = 700.0         !shortwave solar radiation at measurement height (W/m^2)
                        Call PartitionRAD(zarad,pmbzref,ppfdzref,sradzref,ppfd_direct,ppfd_diffus,nir_direct,nir_diffus) 

                        !to get rabs_sun,rabs_shd,ppfd_sun,ppfd_shd,fsun,fshd
                        tkzref = 290.0           !air temperature at measurement height (K)
                        eatm = esat(tkzref)*relhumi/100.0   !esat(kPa),this unit is from CANACC
!                        x = 1.0
!                        kd = 0.60               !readin from input files
!                        tsoilk = 294.0           !readin temp of soil at 5cm (K)
                        lwdnzref = 4.14472E+02   !W/m^2, from CANACC
                        Call CalcRadProfiles(clai,lai,eatm,lwdnzref,nir_direct,nir_diffus,ppfd_direct,ppfd_diffus,tkzref,x,zarad,  &
                                             kd,tk,tsoilk,ppfd_sun,ppfd_shd,rabs_sun,rabs_shd,fsun,fshd)
 
                        !get gs_sun_i & !get gs_sud_i
                        ubarms = 0.01 *ubar(i)                          !change from cm/s to m/s
                        Call CalcOneLevelLEB(rabs_sun(i),ppfd_sun(i),tk(i),pmb(i),ubarms,relhumi,gs_sun_i,  &
                                             gl_sun_i,gb_i,anet_sun_i,tl_sun_i) 
                        Call CalcOneLevelLEB(rabs_shd(i),ppfd_shd(i),tk(i),pmb(i),ubarms,relhumi,gs_shd_i,  &
                                             gl_shd_i,gb_i,anet_shd_i,tl_shd_i) 

                        !get rs_wgt_i
                        Call CalcWeightedProfiles(i,rs_wgt_i,fsun(i),fshd(i),gs_sun_i,gs_shd_i,pmb(i),tk(i))   

                        !convert rs_wgt from (s/m) to (s/cm)
                        rs(i,l) = (mdiffh2o(tk(i),pmb(i))/mdiffl)*rs_wgt_i* 0.01
 
                    case default                                        !Zhang et la., (2002,2003) wiht hardcodes deciduous forest parameters
                        rs(i,l) = rs_zhang_df(mdiffl,tk(i),pmb(i),ppfd(i),srad,relhumi)
                end select
                rnum = rc(i,l) * (rs(i,l) + rm(i,l))
                rden = rc(i,l) + 2.0 * (rs(i,l) + rm(i,l))
                rlx   = rb(i,l) + (rnum/rden)
                vdlx  = 1.0/rlx
                vd(i,l) = vdlx                                           !calcualte deposition velocity (cm/s)
                gp(i,l) = gpla(i,l)                                      !compensation point concentration (mole/cm^3)
           
            else                                                         !out of canopy
                rb(i,l) = 0.0_dp
                rc(i,l) = 0.0_dp
                rm(i,l) = 0.0_dp
                rs(i,l) = 0.0_dp
                vd(i,l) = 0.0_dp
                gp(i,l) = 0.0_dp

            end if
        end do
    end do

    return
end subroutine GetDryDepExCoeffs

!==================================================================
! PART 2. CALCULATE DRY DEPOSITION VELOCITIES TO THE GROUND SURFACE
!==================================================================

subroutine GetSoilDryDepExCoeffs()
    integer(kind = i4)                :: l                       !l is species
    real(kind = dp)                   :: mdiffl                  !molecular diffusivity (cm^2/s)
    real(kind = dp)                   :: tsoilk                  !soil/litter temperature (K)
    real(kind = dp),dimension(ninteg) :: rsoill                  !resistance to diffusion thru soil pore space for chemical species (s/cm)
    real(kind = dp)                   :: rbg                     !ground boundary layer resistance (s/cm)
    real(kind = dp),dimension(ninteg) :: vs                      !soil exchange coefficients (cm/s)
    real(kind = dp)                   :: ubarg
    real(kind=dp),dimension(npts)     :: pmb
    real(kind=dp),dimension(npts)     :: ubar

    pmb    = (/1.01266     , 1.01207     , 1.01148     , 1.010891    , 1.010292    , 1.009702    , 1.009111    , 1.008522    /)              !calc use barometric formula
    ubar   = (/0.540616E+01, 0.688076E+01, 0.898461E+01, 0.121621E+02, 0.174182E+02, 0.278164E+02, 0.861199E+02, 0.861199E+02/)              !data from CANACC

    do l = 1,ninteg    !ninteg is integrated # of species in ACCESS
        mdiffl = MolecDiff(l, tsoilk, pmb(1))
        rsoill(l) = SoilResist(mdiffl)

        ubarg = ubar(1)                                          !get ground level mean wind speed
        rbg = SoilRbg(ubarg)                                     !Rbg(ground boundary layer resistance, s/cm)
                                                                 !Rbg is invariant to species not layers
        vs(l) = 1.0/(rbg+rsoill(l))                              !deposition velocity to ground surface (cm/s)
    end do

    return
end subroutine GetSoilDryDepExCoeffs

end module DryDep
