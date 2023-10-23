!===========================================================================================================
!    Program:    Simplified ACCESS Model
!                Simplified (Dry Deposition only) Atmospheric Chemistry and Canopy Exchange Simulation System
!
!    Developed by : Dr. Beiming Tang based on Dr. Rick Saylor's ACCESS model version 3.1.0
!============================================================================================================

module ACCESS_Modules
    use ACCESS_Constants
    implicit none
    
    public MolecDiff, EffHenrysLawCoeff,ReactivityParam, &
           SoilResist, rs_zhang_df, SoilRbg, &
           mdiffh2o, RelativeHumidity,gpla,  &
           rbl, rcl, rml, CalcWeightedProfiles,  &
           CalcOneLevelLEB,PartitionRad,     &
           CalcRadProfiles, SolarZenithAngle,  &
           esat 
         
      

contains

!===============================================
!PART 1. Functions directly called by DryDep.f90
!
!
!===============================================

!=========================================================================
!function RelativeHumidity - calculate Relative Humidity from the supplied 
!                            specific humidity, temperature and pressure
!=========================================================================
function RelativeHumidity(tki,pmbi,qhi)
    real(kind = dp), intent(in) :: tki
    real(kind = dp), intent(in) :: pmbi
    real(kind = dp), intent(in) :: qhi
    real(kind = dp)             :: RelativeHumidity
    real(kind = dp)             :: e
    real(kind = dp)             :: es
    real(kind = dp)             :: qhd
    real(kind = dp)             :: pkpa
    real(kind = dp)             :: rhi
    real(kind = dp),parameter   :: rhmin = 0.1
    real(kind = dp),parameter   :: rhmax = 99.0

    qhd  = 0.001*qhi
    pkpa = 0.1*pmbi
    e    = pkpa*qhd/(0.622+qhd)
    es   = esat(tki)
    rhi  = max(rhmin,min(rhmax,100.0*e/es))   !bound RH to (rhmin, rhmax)
    RelativeHumidity = rhi

    return
end function RelativeHumidity


!==========================================================================
!function MolecDiff - calculate molecular diffusivities (cm^2/s) at a given
!                     temperature and pressure
!==========================================================================
function MolecDiff(ispec, tkx, pmbx)
    integer(kind=i4), intent(in) :: ispec     !dummy id for species
    real(kind=dp), intent(in)    :: tkx
    real(kind=dp), intent(in)    :: pmbx
    real(kind=dp)                :: MolecDiff
    real(kind=dp),dimension(ntotal)  ::mdiffstp
    
    call SetMolecDiffSTP(mdiffstp)
    MolecDiff = mdiffstp(ispec)*(1013.25_dp/pmbx)*((tkx/298.15_dp)**1.81)

    return
end function MolecDiff


!=========================================================================
!function mdiffh2o - calculate molecular diffusivity of water vapor in air
!
!source: Tracy (1980)
!=========================================================================
function mdiffh2o(tki,pmbi)
    real(kind = dp), intent(in) :: tki
    real(kind = dp), intent(in) :: pmbi
    real(kind = dp)             :: mdiffh2o

    mdiffh2o = 0.226_dp*((tki/273.15_dp)**1.81_dp)*(1000.0_dp/pmbi)

    return
end function mdiffh2o


!==========================================================================
!function ReactivityParam - calculate reactivity parameters (dimensionless)
!==========================================================================
function ReactivityParam(ispec)
    integer(kind=i4), intent(in)    :: ispec   !dummy id for species
    real(kind=dp)                   :: ReactivityParam
    real(kind=dp),dimension(ntotal) :: f0

    call SetReactivityParams(f0)
    ReactivityParam = f0(ispec)

    return
end function ReactivityParam


!========================================================================
!function EffHenrysLawCoeff - calculate Henry's Law coefficient (M/atm) 
!                             have not include temperature dependence yet
!========================================================================
function EffHenrysLawCoeff(ispec)
    integer(kind=i4), intent(in)    :: ispec  !dummy id for species
    real(kind=dp)                   :: EffHenrysLawCoeff
    real(kind=dp),dimension(ntotal) :: hstar     

    call SetEffHenrysLawCoeffs(hstar)
    EffHenrysLawCoeff = hstar(ispec)

    return
end function EffHenrysLawCoeff


!===============================================================================================================
!function SoilResist - calculate the resistance to diffusion of a species from the free warer surface in the soil
!                      to the soil-atmosphere interface. Rsoil 
!Source - Sakagichi & Zeng (2009)
!================================================================================================================
function SoilResist(mdiffl)
    real(kind=dp), intent(in)  :: mdiffl       !molecular diffusivity of species in air (cm^2/s)
    real(kind=dp)              :: SoilResist   !Soil resistance (s/cm)
    real(kind=dp)              :: xe           !temporary variable
    real(kind=dp)              :: ldry         !diffusion distance through the soil (cm)
    real(kind=dp)              :: mdiffp       !effective diffusivity of species through the soil (cm^2/s)

    xe = (1.0-(stheta/sattheta))**5.0
    ldry = dsoil*(exp(xe)-1.0)/1.7183
    ldry = max(0.0,ldry)

    mdiffp = mdiffl*sattheta*sattheta*(1.0-(rtheta/sattheta))**(2.0+3.0/sbcoef)
    SoilResist = ldry/mdiffp

    return
end function SoilResist


!=====================================================================================
!function SoilRbg - calculate the boundary layer resistance at the ground surface. Rbg
!
!Source - Schuepp (1977) 
!=====================================================================================
function SoilRbg(ubarg)
    real(kind=dp), intent(in)  :: ubarg           !mean wind speed in the 1st model layer (cm/s)
    real(kind=dp)              :: SoilRbg         !Boundary layer resistance at ground surface (s/cm)
    real(kind=dp), parameter   :: rbgmax = 1.67   !maximum ground surface boundary layer resistance (s/cm)
    real(kind=dp)              :: rbg             !temporary variable for boundary layer resistance (s/cm)

    rbg = 11.534/(0.14*ubarg)                     !assume Sc=0.7,del0/zl = 0.02 and ustar = 0.14*ubar (Weber,1999)
    SoilRbg = min(rbgmax,rbg)

    return
end function SoilRbg 


!===================================================================
!function rbl - calculate leaf boundary resistance for trace species
!
!Source - rb formulation from Wu et al., (2003)
!       - ustar = 0.14*ubar from Weber (1999)
!===================================================================
function rbl(mdiffl, ubari)
    real(kind=dp)               :: rbl            !leaf boundary layer resistance (s/cm)
    real(kind=dp), intent(in)   :: mdiffl         !molecular diffusivity of species in air (cm^2/s)
    real(kind=dp), intent(in)   :: ubari          !mean wind speed at layer i (cm/s)
    
    rbl = 10.53_dp/((mdiffl**0.666667)*ubari)

    return
end function rbl


!===============================================================
!function rcl - calculate cuticular resistance for trace species
!
!Source - Wesely (1989)
!===============================================================
function rcl(hstarl,f01)
    real(kind=dp), intent(in)   :: hstarl         ! effective Henry's law coefficient (M/atm)
    real(kind=dp), intent(in)   :: f01            ! reactivity parameter (0-1)
    real(kind=dp)               :: rcl            ! cuticular resistance (s/cm)
    real(kind=dp), parameter    :: rcref=20.0     ! rc for ozone (s/cm) for deciduous forest

    rcl = rcref/((hstarl*1.0D-05)+f01)

    return
end function rcl


!===============================================================
!function rml - calculate mesophyll resistance for trace species
!
!Source - Wesely (1989)
!===============================================================
function rml(hstarl,f01)
    real(kind=dp), intent(in)   :: hstarl          ! effective Henry's law coefficient (M/atm)
    real(kind=dp), intent(in)   :: f01             ! reactivity parameter (0-1)
    real(kind=dp)               :: rml             ! mesophyll resistance (s/cm)

    rml = 1.0_dp/((hstarl/3000.0_dp)+100.0_dp*f01)

    return
end function rml


!==================================================================================
!function gpla - calculate stomatal compensation point for ACCESS-ORG trace species
!
!Source - ???
!==================================================================================
function gpla(ilev,lsp)
    integer(kind=i4), intent(in) :: ilev           !ilev is layer
    integer(kind=i4), intent(in) :: lsp            !lsp is species
    real(kind=dp)                :: gpla

    gpla = 0.0_dp                                  !dummy stub for now

    return
end function gpla


!======================================================================
!function rs_zhang_df - calcualte stomatal resistance for trace species
!
!Source - Zhang et al., (2002 & 2003)
!       - 'rsmin' by Wesely et al (1989)
!       - 'bvpd' by Wolfe and Thornton (2011)
!======================================================================
function rs_zhang_df(mdiffl, tki, pmbi, ppfdi, srad, relhumi)
    real(kind = dp), intent(in) :: mdiffl            !molecular diffusivity of trace species in air (cm^2/s)
    real(kind = dp), intent(in) :: tki               !temperature
    real(kind = dp), intent(in) :: pmbi              !air pressure (mb)
    real(kind = dp), intent(in) :: ppfdi             !photosynthetic photon flux (umol/m^2-s)
    real(kind = dp), intent(in) :: srad              !solar irradiation (W/m^2)
    real(kind = dp), intent(in) :: relhumi           !relative humidity (%)
    real(kind = dp)             :: rs_zhang_df       !stomatal resistance (s/cm)
    real(kind = dp), parameter  :: rsmin =1.0        !minimum leaf stomatal resistance (s/cm) for deciduous forest
    real(kind = dp), parameter  :: rsmax =10000.     !maximum leaf stomatal resistance (s/cm) (stoma are closed)
    real(kind = dp), parameter  :: brsp = 196.5      !empirical constant (umol/m^2-s) for deciduous forest
    real(kind = dp), parameter  :: tmin = 0.0        !temperature correction parameter-deciduous forest
    real(kind = dp), parameter  :: tmax = 45.0       !temperature correction parameter-deciduous forest
    real(kind = dp), parameter  :: topt = 27.0       !temperature correction parameter-deciduous forest
    real(kind = dp), parameter  :: bvpd = 0.10       !empirical constant for VPD correction-deciduous forest
    real(kind = dp), parameter  :: phic1 = -1.9      !empirical constant for water stress correction-deciduous forest
    real(kind = dp), parameter  :: phic2 = -2.5      !empirical constant for water stress correction-deciduous forest
    real(kind = dp)             :: cft
    real(kind = dp)             :: cfvpd
    real(kind = dp)             :: cfphi
    real(kind = dp)             :: tcel
    real(kind = dp)             :: ft1
    real(kind = dp)             :: ft2
    real(kind = dp)             :: et
    real(kind = dp)             :: vpd
    real(kind = dp)             :: phi

    !temperature correction
    tcel = tki - 273.15
    et   = (tmax-topt)/(topt-tmin)
    ft1  =(tcel-tmin)/(topt-tmin)
    ft2  = (tmax-tcel)/(tmax-topt)
    cft  = ft1*(ft2**et)

    !water vapor pressure defit correction
    vpd  = esat(tki)*(1.0_dp - (relhumi/100.0_dp))
    cfvpd= 1.0_dp - bvpd*vpd

    !water stress correction
    phi  = -0.72_dp - 0.0013_dp*srad
    cfphi= (phi-phic2)/(phic1-phic2)

    if (ppfdi >0.0) then
        rs_zhang_df = rsmin*(1.0_dp+brsp/ppfdi)*mdiffh2o(tki,pmbi)/(mdiffl*cft*cfvpd*cfphi)
    else
        rs_zhang_df = rsmax                         !nighttime, stoma are closed
    endif

    return
end function rs_zhang_df

!==================================================================
!CalWeightedProfiles - calculate sun/shade weighted canopy profiles
!
!==================================================================

subroutine CalcWeightedProfiles(i,rs_wgt_i,fsun_i,fshd_i,gs_sun_i,gs_shd_i,pmbi,tki)    !Beiming modify this code
    integer(kind=i4), intent(in) :: i                     !i is layer
    real(kind=dp),intent(out)    :: rs_wgt_i
    real(kind=dp),intent(in)     :: fsun_i
    real(kind=dp)                :: rs_sun_i
    real(kind=dp),intent(in)     :: fshd_i
    real(kind=dp)                :: rs_shd_i
    real(kind=dp),intent(in)     :: gs_sun_i
    real(kind=dp),intent(in)     :: gs_shd_i    
    real(kind=dp),intent(in)     :: pmbi 
    real(kind=dp),intent(in)     :: tki    
   
    if (i <= ncnpy) then            !ncnpy is # of vertical layers
        rs_sun_i = gtor(gs_sun_i,pmbi,tki)
        rs_shd_i = gtor(gs_shd_i,pmbi,tki)
    endif
    rs_wgt_i = fsun_i*rs_sun_i + fshd_i*rs_shd_i

    return
end subroutine CalcWeightedProfiles

 
!=============================================================
!PART 2. 2nd level functions called by functions in DryDep.f90
!
!
!=============================================================

!==============================================================================================
!function esat - calculate the saturation vaport pressure (kPa) of water at a given temperature
!
!source - Rogers et al., (1989)
!notes - valid over -30 C <= T <= 35 C
!==============================================================================================
function esat(tki)
    real(kind=dp),intent(in) :: tki    !temperature,unit = K
    real(kind=dp)            :: esat   !satuaration vapor pressure,unit = kPa
    real(kind=dp)            :: tc     !temperature, unit = C

    tc = tki - 273.15
    esat = 0.6112*exp(17.67*tc/(tc+243.5))

    return
end function esat

!=====================================================================================
!function mdiffstp - set molecular diffusivity data (cm^2/s) for all species at
!                           0 deg C and 1 atm 
!=====================================================================================
subroutine SetMolecDiffSTP(mdiffstp)
    real(kind = dp),dimension(ntotal),intent(out)  :: mdiffstp       !molecular diffusivities of species in air at 0degC and 1 atm [cm^2/s]
!    real(kind = dp), parameter         :: mdiffstp_default = 0.100   !default value of mdiffstp (cm^2/s) with no reliable data
!    integer(kind=i4)                   :: l                          !l is species
    
!    do l=1,ntotal      !ntotal is total species in ACCESS based on chosen chemical mechanim, including transported species
!        mdiffstp(l) = mdiffstp_default
!    end do

!Insert MolecDiffSTP Coefficients for RACM2_plus mechanism
!species in array are:
!         = (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
!            MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
!            PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
!           ISNP /)  
    mdiffstp =(/0.1802, 0.1361, 0.1444, 0.1349, 0.1041, 0.1041, 0.0808, 0.1807, 0.1300, 0.1952,  &
                0.1297, 0.1200, 0.1297, 0.1153, 0.2773, 0.2773, 0.2543, 0.2000, 0.1340, 0.1060,  &
                0.1040, 0.0892, 0.0845, 0.0837, 0.0837, 0.0834, 0.0750, 0.0750, 0.0745, 0.0745,  &
                0.0712 /)
!species (NO2, O3)
!    mdiffstp =(/0.1361, 0.1444/) 

    return
end subroutine SetMolecDiffSTP

!========================================================
!function f0 - set reactivity parameters for all species
!
!source - Wesely et al., (1989) and Nguyen et al., (2015)
!========================================================
subroutine SetReactivityParams(f0)
!    integer(kind=i4)                              :: l               !l is species
    real(kind = dp),dimension(ntotal),intent(out) :: f0
!    real(kind=dp),parameter                       :: f0_default = 0.0
    
!    do l=1,ntotal     !ntotal is total species in ACCESS based on chosen chemical mechanim, including transported species
!        f0(l) = f0_default
!    end do

!Insert ReactivityParams Coefficients for RACM2_plus mechanism
!species in array are:
!         = (/NO,   NO2,    O3,  HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
!            MO2,   OP1,   MOH,   NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
!            PAA, DHMOB, HPALD,  ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
!           ISNP /)  
    f0 = (/0.0, 0.1, 1.0, 0.1, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  &
           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0,  &
           0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  &
           0.0 /)

!species (NO2, O3)
!    f0 = (/0.1, 1.0 /)

    return
end subroutine SetReactivityParams 

!====================================================================
!function hstar - set henry's law coefficient for all species (M/atm)
!
!source - Nguyen et al., (2015)
!====================================================================
subroutine SetEffHenrysLawCoeffs(hstar)
!    integer(kind=i4)                            :: l                  !l is species
    real(kind=dp),dimension(ntotal),intent(out) :: hstar
!    real(kind=dp),parameter                     :: hstar_default = 1.000
    
!    do l = 1,ntotal    !ntotal is total species in ACCESS based on chosen chemical mechanim, including transported species
!        hstar(l) = hstar_default
!    End do

!Insert EffHenryLaw Coefficients for RACM2_plus mechanism
!species in array are:
!         = (/NO,   NO2,    O3, HONO,  HNO4,   HNO3,   N2O5,     CO,  H2O2,  CH4,
!            MO2,   OP1,   MOH,  NO3,   O3P,    O1D,     HO,    HO2,  ORA1,  HAC,
!            PAA, DHMOB, HPALD, ISHP, IEPOX, PROPNN, ISOPNB, ISOPND, MACRN, MVKN,
!           ISNP /)    

!    hstar = (/0.0019, 0.0120, 0.0103, 2.6D+05, 1.0D+07, 3.2D+13, 1.0D+14, 9.8D-04, 8.4D+04, 1.4D-03,  &
!             6.9D+02, 3.0D+02, 2.2D+02, 0.0380, 0.0380, 0.0380, 3.9D+01, 6.9D+02, 5.6D+03, 2.0D+03,  & 
!             5.2D+02, 2.0D+03, 4.0D+04, 7.0D+07, 7.0D+07, 1.0D+04, 5.0D+03, 5.0D+03, 6.0D+03, 6.0D+03,  &
!             5.0D+03 /)

    hstar = (/1.9D-03, 1.2D-02,1.03D-02, 2.6D+05, 1.0D+07, 3.2D+13, 1.0D+14, 9.8D-04, 8.4D+04, 1.4D-03,  &
              6.9D+02, 3.0D+02, 2.2D+02, 3.8D-02, 3.8D-02, 3.8D-02, 3.9D+01, 6.9D+02, 5.6D+03, 2.0D+03,  &
              5.2D+02, 2.0D+03, 4.0D+04, 7.0D+07, 7.0D+07, 1.0D+04, 5.0D+03, 5.0D+03, 6.0D+03, 6.0D+03,  &
              5.0D+03 /)

!species (NO2, O3)
!    hstar = (/0.0120, 0.0103/)
    
    return
end subroutine SetEffHenrysLawCoeffs

!==================================================================
!function gtor - convert conductance (mol/m^2-s) to resistance (s/m)
!
!notes - assumes non-zero value for gz!
!==================================================================
function gtor(gz, pmbi,tki)
    real(kind=dp), intent(in)     :: gz                 !conductance, mol/m^2-s
    real(kind=dp), intent(in)     :: pmbi               !air pressure, mb
    real(kind=dp), intent(in)     :: tki                !temperature, K
    real(kind=dp)                 :: gtor               !resistance, s/s
    real(kind=dp)                 :: gms                !conductance, m/s
    real(kind=dp), parameter      :: rgas=8.205D-05     !ideal gas constant, m^3-atm/K-mol

    gms = gz*(rgas*tki)/(pmbi/1013.)
    gtor = 1.0/gms

    return
end function gtor

!========================================================================
!function gs_medlyn - calculate the stomatal conductance for water vapor,
!                      mol/m^2-s
!
!Source -  Medlyn, B. E., et al., (2011) 
!========================================================================
subroutine gs_medlyn_s(tki, relhumi, anet, cci,gs_medlyn)
    real(kind=dp), intent(in)    :: tki     !air temperature, K
    real(kind=dp), intent(in)    :: relhumi  !relative humidity
    real(kind=dp), intent(in)    :: anet    !net photosynthesis assimilation rate, umol/m^2-s
    real(kind=dp), intent(out)   :: cci    !intercellular CO2 concentration,umol/mol
    real(kind=dp), intent(out)   :: gs_medlyn  !stomatal conductance for water vapor, mol/m^2-s
    real(kind=dp)                :: vpd     ! vapor pressure dificit
    real(kind=dp)                :: gs
    real(kind=dp),parameter      :: g0 = 0.02    !this is a readin value,0.02 is deciduous forest
    real(kind=dp),parameter      :: g1 = 4.00    !this is a readin value, 4.00 is deciduous forest
    real(kind=dp),parameter      :: cca = 405.0  !atmospheric CO2 mixing ratio (umol/mol)
    real(kind=dp),parameter      :: gsmax = 10.0 !maximum stomatal conductance parameter, mol/m^2-s

    vpd = esat(tki)*(1.0-relhumi/100.)
    if (vpd >0.0) then
        cci = cca*(g1/(g1 + sqrt(vpd)))
        gs  = g0 + 1.6*anet*(1.0 + g1/sqrt(vpd))/cca
    else
        cci = cca
        gs  = gsmax
    end if

    gs_medlyn = gs

    return
end subroutine gs_medlyn_s

!================================================================================================================
!subroutine CalcOneLevelLEB - calculate leaf energy balance and leaf temperature, net photosynthesis assimilation
!                             rate and stomatal conductance for one level in the canopy
!
!iterate through equations until stomatal conductance converges to specified tolerance
!================================================================================================================
subroutine CalcOneLevelLEB(rabs, ppfdi, tki, pmbi, ubarms, relhumi, gs, gv, gb, anet, tleaf)
    real(kind=dp), intent(in)    :: rabs
    real(kind=dp), intent(in)    :: ppfdi
    real(kind=dp), intent(in)    :: tki
    real(kind=dp), intent(in)    :: pmbi
    real(kind=dp), intent(in)    :: ubarms          !mean wind speed, m/s
    real(kind=dp), intent(in)    :: relhumi
    real(kind=dp), intent(out)    :: gs             !leaf stomatal conductance, mol/m^2-s
    real(kind=dp), intent(out)    :: gb             !boundary layer conductance,mol/m^2-s
    real(kind=dp), intent(out)    :: gv             !overall leaf conductance, mol/m^2-s
    real(kind=dp), intent(out)    :: anet           !net photosynthesis assimilation rate, umol/m^2-s
    real(kind=dp), intent(out)    :: tleaf          !leaf temperature, K
    real(kind=dp)                 :: cci            !intercellular CO2 mixing ratio, umol/mol
    real(kind=dp)                 :: delgs          !absolute difference in gs between iterations
    real(kind=dp)                 :: gh             !heat conductance, mol/m^2-s
    real(kind=dp), parameter      :: gr=0.15        !radiative conductance, mol/m^2-s (Rr=300 s/m)
    real(kind=dp)                 :: gs_new          !updated value of gs, mol/m^2-s
    integer(kind=i4), parameter   :: maxiter = 10
    integer(kind=i4)              :: iter
    real(kind=dp), parameter      :: dleaf= 0.06    !characteristic leaf dimension = 0.75*(leaf width) (m),read-in from input, deciduous forest is 0.06
    real(kind=dp), parameter      :: ltol=0.001     !iteration tolerance used in calculating leaf temperature
    real(kind=dp)                 :: gs_medlyn   
 
    !initial guesses
    cci = 200.0
    gs  = 0.2

    !heat conductance with radiative conductance
    gh  = gh_cn(ubarms, dleaf, gr)

    !boundary layer conductance
    gb  = gb_cn(ubarms, dleaf)

    delgs = gs
    iter  = 0
    do while (delgs/gs > ltol .and. iter <= maxiter)
    
        !overall leaf conductance
        gv = gv_cn(gb,gs)

        !calculate leaf temperature
        tleaf = LeafTemp(tki,pmbi,relhumi,rabs,gh,gv)

        !calculate photosynthesis assimilation
        anet = NetPhoto(tleaf,cci,ppfdi)

        !calculate new stomatal conductance
        Call gs_medlyn_s(tleaf, relhumi, anet, cci, gs_medlyn)
        gs_new = gs_medlyn

        delgs = abs(gs_new-gs)

        gs = gs_new

        iter = iter + 1
    end do
    return
end subroutine CalcOneLevelLEB

!==========================================================================================================
!function gh_cn - calculate the heat conductance using the forced convection formulation of
!                 Campbell & Norman, Table 7.6, p.109, plus an assumed value for the radiative conductance.
!
!source - Campbell, G. S., and Norman, J. M. (1998)
!==========================================================================================================
function gh_cn(ubar, dleaf, gr)
    real(kind=dp), intent(in) :: ubar     !mean wind speed, m/s
    real(kind=dp), intent(in) :: dleaf    !leaf dimension, m
    real(kind=dp), intent(in) :: gr       !radiative conductance, mol/m^2-s
    real(kind=dp)             :: gh_cn    !total heat conductance, mol/m^2-s

    gh_cn = 1.4*0.135*sqrt(ubar/dleaf)
    gh_cn = gh_cn*gr/(gh_cn+gr)

    return
end function gh_cn

!===============================================================================================================
!function gb_cn - calculate the boundary layer conductance using the forced convection formulation of Campbell & 
!                 Norman, Table 7.6, p. 109
!source - Campbell, G. S., and Norman, J. M. (1998)
!===============================================================================================================
function gb_cn(ubar, dleaf)
    real(kind=dp), intent(in)  :: ubar    !mean wind speed, m/s
    real(kind=dp), intent(in)  :: dleaf   !leaf dimension, m
    real(kind=dp)              :: gb_cn   !boundary layer con

    gb_cn = 1.4*0.147*sqrt(ubar/dleaf)

    return
end function gb_cn

!===========================================================================================================
!function gv_cn - calculate the overall leaf conductance using approximate formulation of Campbell & Norman,
!                 Table 7.6, p. 109
!source - Campbell, G. S., and Norman, J. M. (1998)
!===========================================================================================================
function gv_cn(gb,gs)
    real(kind=dp), intent(in)  :: gb   !boundary layer conductance, mol/m^2-s
    real(kind=dp), intent(in)  :: gs   !stomatal conductance, mol/m^2-s
    real(kind=dp)              :: gv_cn   !overall leaf conductance, mol/m^2-s

    gv_cn = gb*gs/(gb+gs)

    return
end function gv_cn

!=====================================================================================================================
!function LeafTemp -calculate the leaf temperature
!                   Uses simple Newton-Raphson method to compute tleaf based on code from Press et al., (1992)
!                   Numerical Recipes in Fortran 77, Chap. 9. Root Finding and Nonlinear Sets of Equations
!                   No approximations are made in this formulation, which should produce the most accurate estimate of 
!                   leaf temperature
!=====================================================================================================================
function LeafTemp(tki, pmbi, relhumi, rabs, gh, gv)
    real(kind=dp), intent(in)  :: tki      !air temperature, K
    real(kind=dp), intent(in)  :: pmbi     !air pressure, mb
    real(kind=dp), intent(in)  :: relhumi    !relative humidity, %
    real(kind=dp), intent(in)  :: rabs       !absorbed radiation, W/m^2
    real(kind=dp), intent(in)  :: gh         !heat conductance for leaf-atmosphere
    real(kind=dp), intent(in)  :: gv         !leaf water vapor conductance
    real(kind=dp)              :: LeafTemp   !leaf temperature, K
    integer(kind=i4)           :: i
    integer(kind=i4), parameter  :: maxiter=10   !maximum number of iterations
    real(kind=dp), parameter     :: rtol=0.1     !convergence tolerance
    real(kind=dp)                :: tair
    real(kind=dp)                :: tleaf
    real(kind=dp)                :: dt
    real(kind=dp)                :: f
    real(kind=dp)                :: df

    tleaf = tki
    tair  = tki

    do i=1,maxiter
        call LeafBalEq(tleaf, tair, pmbi, relhumi, rabs, gh, gv, f, df)
        dt = f/df
        tleaf = tleaf - dt
        if (abs(dt) < rtol) then
            LeafTemp = tleaf
            return
        end if
    end do
    write(*,*) 'LeafTemp non-convergence!'
    LeafTemp = tki

    return
end function LeafTemp

!=============================================================================================
!subroutine LeafBalEq - computes value for leaf energy balance equation and its 1st derivative
!
!used for full non-linear version of the solution empolying the Newton_Raphson method
!=============================================================================================
subroutine LeafBalEq(tleaf, tair, pmbi, relhumi, rabs, gh, gv, f, df)
    real(kind=dp), intent(in)  :: tleaf     !leaf temperature, K
    real(kind=dp), intent(in)  :: tair      !air temperature, K
    real(kind=dp), intent(in)  :: pmbi      !air pressure, mb
    real(kind=dp), intent(in)  :: relhumi   !relative humidity, %
    real(kind=dp), intent(in)  :: rabs      !absorbed radiation, W/m^2
    real(kind=dp), intent(in)  :: gh        !heat conductance for leaf-atmosphere
    real(kind=dp), intent(in)  :: gv        !leaf water vapor conductance
    real(kind=dp)              :: vpd       !vapor pressure deficit, kPa
    real(kind=dp)              :: patm      !air pressure, kPa
    real(kind=dp)              :: s         !first derivative of esat parameter
    real(kind=dp),intent(out)  :: f
    real(kind=dp),intent(out)  :: df
    real(kind=dp), parameter   :: epsleaf=0.97   
    real(kind=dp), parameter   :: sbsig = 5.67D-08  ! Stefan-Boltzmann constant (W/m2-K4)
    real(kind=dp), parameter   :: cpair=29.3        !heat capacity of air, J/mol-K



    vpd = esat(tleaf) -esat(tair)*relhumi/100.0
    patm = pmbi*0.1   !convert mb to kPa
    s = desdt(tleaf)/patm

    f = -rabs +2.0*epsleaf*sbsig*tleaf**4.0 + cpair*gh*(tleaf-tair) + lambda(tleaf)*gv*vpd/patm
    df = 8.0*epsleaf*sbsig*tleaf**3.0 +cpair*gh +lambda(tleaf)*gv*s

    return
end subroutine LeafBalEq

!================================================================================
!function desdt - calculate the 1st derivative of es wrt T at a given temperature
!
!    Derived from the esat expression in ...
!    Rogers, R. R., and Yau, M. K. (1989) A short Course in Cloud Physics,
!    3rd Ed., Butterworth-Heinemann, Burlington, MA. (p16)
!    Valid over -3C <= T <= 35C
!================================================================================
function desdt(tki)
    real(kind=dp), intent(in)  :: tki   !temperature (K)
    real(kind=dp)              :: desdt   !des/dt (kPa/C)
    real(kind=dp)              :: tc      !temperature (C)

    tc = tki -273.15
    desdt = esat(tki)*4302.6/((tc+243.5)**2.0)
    return
end function desdt

!=========================================================================================
!function lambda - calculate the latent heat of evaporation for H2O at a given temperature
!
!    Polynominal fit of data in Table 2.1 (p.16) of ...
!    Rogers, R. R., and Yau, M. K. (1989) A Short Course in Cloud Physics,
!    3rd Ed., Butterworth-Heinemann, Burlington, MA. (p16)
!    Valid over -40C <=T <=40C
!==========================================================================================
function lambda(tki)
    real(kind=dp), intent(in)  :: tki   !temperature (K)
    real(kind=dp)              :: lambda     !latent heat of evaporation (J/mol)
    real(kind=dp)              :: tc         !temperature (C)

    tc = tki - 273.15
    lambda = (2500.8 - 2.36*tc + 0.0016*tc*tc - 0.00006*tc*tc*tc)*18.0
    return
end function lambda

!====================================================================================================================
!function NetPhoto - calculates net photosynthesis assimilation rate
!      
!Use the formulation of ...
!    Collatz et al., (1991) Physiological and environmental regulation of stomatal conductance, photosynthesis, and 
!    transpiration: a model that includes a laminar boundary layer, Agric. For. Meteorol., 54, 107-136.
!as expressed in ...
!    Campbell, G. S., and Norman, J. M. (1998) An introduction to Environmental Biophysics, 2nd Ed. Springer, New York
!=====================================================================================================================
function NetPhoto(tleaf,cci,ppfd)
    implicit none
    real(kind=dp), intent(in) :: tleaf    !leaf temperature, K
    real(kind=dp), intent(in) :: cci      !intercellular CO2 concentration, umol/mol
    real(kind=dp), intent(in) :: ppfd     !photosynthetic photon flux density, umol/m^2-s
    real(kind=dp)             :: NetPhoto !net photosynthesis assimilation rate, umol/m^2-s
    real(kind=dp)             :: tlc      !leaf temperature, C
    real(kind=dp)             :: agross   !gross photosynthesis assimilation rate, umol/m^2-s
    real(kind=dp)             :: rd       !leaf respiration rate, umol/m^2-s
    real(kind=dp)             :: jp       !tmp variable, colimitation assimilation, umol/m^2-s
    real(kind=dp)             :: js       !photosynthesis assimilation rate imposed by sucrose synthesis,umol/m^2-s
    real(kind=dp)             :: jc       !Rubisco-limited photosynthesis assimilation rate,  umol/m^2-s
    real(kind=dp)             :: je       !light-limited photosynthesis assimilation rate, umol/m^2-s
    real(kind=dp)             :: gamstr   !light compensation point -CO2 conc at which assimilation is 0,umol/m^2-s
    real(kind=dp)             :: tau      !CO2/O2 specificity ratio, mmol/umol
    real(kind=dp)             :: vm       !max Rubisco capacity per unit leaf area, umol/m^2-s
    real(kind=dp)             :: kc       !Michaelis constant for CO2, umol/mol
    real(kind=dp)             :: k0       !inhibition constant for O2, umol/mol
    real(kind=dp), parameter  :: tau25=2600.0  !CO2/O2 specificity ratio @25degC, mmol/mmol
    real(kind=dp), parameter  :: vm25=54.0   !max Rubisco capacity per unit leaf area @25degC,umol/m^2-s
                                             !from Houborg et al., (2009)
    real(kind=dp), parameter  :: kc25=300.0  !Michaelis constant for CO2 @25degC, umol/mol
    real(kind=dp), parameter  :: k025=300000.0 !inhibition constant for O2 @25degC, umol/mol
    real(kind=dp), parameter  :: coa=210000.0  !oxygen mole fraction, umol/mol
    real(kind=dp), parameter  :: alfp=0.8      !leaf absorptivity for PPFD
    real(kind=dp), parameter  :: em=0.08       !maximum quantum efficiency, mol/mol
    real(kind=dp), parameter  :: beta=0.98     !transition parameter between Jp and Js
    real(kind=dp), parameter  :: theta=0.95    !transition parameter between Jp and Js
    real(kind=dp), parameter  :: rd25=1.5      !leaf respiration rate @25degC, umol/m^2-s

    !equation require leaf temperature in degC
    tlc = tleaf-273.15

    !ratio describing partitioning of RuBP to the carboxylase and oxygenase
    !reactions of Rubisco
    tau = tau25*exp(-0.056*(tlc-25.0))

    !CO2 conc where assimilation is zero
    gamstr = coa/(2.0*tau)

    !light-limited photosynthesis assimilatiion rate
    je = alfp*em*ppfd*(cci-gamstr)/(cci+2.0*gamstr)

    !max Rubisco capacity per unit leaf area
    vm = vm25*exp(0.088*(tlc-25.0))/(1.0+exp(0.29*(tlc-41.0)))

    !Michaelis constant for CO2
    kc = kc25*exp(0.074*(tlc-25.0))

    !inhibition constant for O2
    k0 = k025*exp(0.018*(tlc-25.0))

    !Rubisco-limited photosynthesis assimilation rate
    jc = vm*(cci-gamstr)/(cci+kc*(1.0+coa/k0))

    !colimitation between Je and Jc
    jp = (je+jc -sqrt((je+jc)*(je+jc) -4.0*theta*je*jc))/(2.0*theta)

    !sucrose-limited photosynthesis assimilation rate
    js = vm*0.5    

    !gross photosynthesis assimilationrate
    agross = (jp + js -sqrt((jp+js)*(jp+js) -4.0*beta*jp*js))/(2.0*beta)

    !leaf respiration rate
    rd = rd25*exp(0.069*(tlc-25.0))/(1.0+exp(1.3*(tlc-55.0)))
!    rd = 0.015*vm   !used by Collatz et al.

    !net photosynthesis assimilation rate
    NetPhoto = max(0.0, agross -rd)

    return
end function NetPhoto

!===============================================================================================================
!Subroutine CalcRadProfiles - compute radiation profiles within the canopy
!
!Uses algorithms of Bodin & Franklin (2012) Efficient modeling of sun/shade canopy radiation dynamics
!                   explicitly accounting for scattering, Geoscientific Model Development, 5, 535-541.
!  ... with elliptical kb from Campbell & Norman (1998) pp. 247-259
!  ... clear-sky effective emissivity of Prata (1996) Q.J.R. Meteorol. Soc., 122, 1127-1151.
!  ... cloudy-sky correction algorithm of Crawford and Duchon (1999) J. Appl. Meteorol., 48, 474-480.
!  ... based on evaluation of Flerchinger et al., (2009) Water Resour. Res., 45, W03423, doi:10.1029/2008WR007394
!================================================================================================================
subroutine CalcRadProfiles(clai,lai,eatm,lwdnzref,nir_direct,nir_diffus,ppfd_direct, ppfd_diffus,tkzref,x,zarad,  &
                           kd,tk,tsoilk,ppfd_sun,ppfd_shd,rabs_sun,rabs_shd,fsun,fshd)
    integer(kind=i4)                 :: i
    real(kind=dp)                    :: kbexza   !canopy extinction coefficient for given x and za values
    real(kind=dp)                    :: sky_ir   !effective sky thermal radiation (W/m^2)
    real(kind=dp)                    :: grnd_ir  !ground surface thermal radiation (W/m^2)
    real(kind=dp)                    :: epssky   !effective emissivity for sky thermal radiation
    real(kind=dp)                    :: w        !precipitable water (cm)
    real(kind=dp), parameter         :: rfl_par=0.11   !leaf reflectance coefficient for PAR
    real(kind=dp), parameter         :: trns_par=0.16  !leaf transmission coefficient for PAR
    real(kind=dp), parameter         :: sigma_par = rfl_par + trns_par  !leaf scattering coefficient for PAR
    real(kind=dp), parameter         :: rfl_nir=0.43   !leaf reflectance coefficient for NIR
    real(kind=dp), parameter         :: trns_nir=0.26  !leaf transmission coefficient for NIR
    real(kind=dp), parameter         :: sigma_nir=rfl_nir+trns_nir  !leaf scattering coefficient for NIR
    real(kind=dp), parameter         :: alfpar=0.85    !leaf absorptance for PAR
    real(kind=dp), parameter         :: alfnir=0.25    !leaf absorptance for NIR
    real(kind=dp), parameter         :: alfir=1.0      !leaf absorptance for IR
    real(kind=dp)                    :: taudt_ir
    real(kind=dp)                    :: s
    real(kind=dp)                    :: kdsig_par
    real(kind=dp)                    :: kdrfl_par
    real(kind=dp)                    :: kdtrs_par
    real(kind=dp)                    :: rho_par
    real(kind=dp)                    :: kdsig_nir
    real(kind=dp)                    :: kdrfl_nir
    real(kind=dp)                    :: kdtrs_nir
    real(kind=dp)                    :: rho_nir
    real(kind=dp)                    :: rdpar
    real(kind=dp)                    :: rdnir
    real(kind=dp)                    :: rsuppar
    real(kind=dp)                    :: rsupnir
    real(kind=dp)                    :: rsdnpar
    real(kind=dp)                    :: rsdnnir
    real(kind=dp), dimension(npts), intent(in)        :: clai      !top-down cumulative leaf area index of canopy (cm^2 leaf/cm^2)
    real(kind=dp), dimension(npts), intent(in)        :: lai       !layer leaf area index of canopy (cm^2 leaf/cm^2)
    real(kind=dp), dimension(npts), intent(out)       :: ppfd_sun
    real(kind=dp), dimension(npts)                    :: rt_sun    !total absorbed PPFD and NIR radiation in sunlit fractions (W/m^2)
    real(kind=dp), dimension(npts)                    :: nir_sun
    real(kind=dp), dimension(npts), intent(out)       :: rabs_sun
    real(kind=dp), dimension(npts)                    :: lw_dn
    real(kind=dp), dimension(npts)                    :: lw_up
    real(kind=dp), dimension(npts), intent(out)       :: ppfd_shd
    real(kind=dp), dimension(npts)                    :: rt_shd
    real(kind=dp), dimension(npts)                    :: nir_shd
    real(kind=dp), dimension(npts), intent(out)       :: rabs_shd
    real(kind=dp), intent(in)                         :: eatm    !water vapor pressure (kPa) 
    real(kind=dp), parameter                          :: epsgrnd = 0.95    ! soil/litter surface emissivity
    real(kind=dp)                                     :: ESTIMATED
    real(kind=dp)                                     :: MEASURED
    real(kind=dp)                                     :: SKYLW
    real(kind=dp), intent(in)                         :: lwdnzref
    real(kind=dp), intent(in)                         :: nir_direct
    real(kind=dp), intent(in)                         :: nir_diffus
    real(kind=dp), intent(in)                         :: ppfd_direct
    real(kind=dp), intent(in)                         :: ppfd_diffus
    real(kind=dp), parameter                          :: sbsig =5.67D-08
    real(kind=dp), intent(in)                         :: tkzref     ! air temperature (K)
    real(kind=dp), intent(in)                         :: x       ! leaf angle distribution parameter in ellipsoidal kb function
    real(kind=dp), intent(in)                         :: zarad   ! zenith angle at current meteorological data time step (radians)
    real(kind=dp), dimension(npts), intent(out)       :: fsun
    real(kind=dp), dimension(npts), intent(out)       :: fshd
    real(kind=dp)                                     :: epsleaf=0.97   ! leaf emissivity
    real(kind=dp), intent(in)                         :: kd      ! diffuse radiation extenction coefficient
    real(kind=dp),dimension(npts), intent(in)         :: tk      ! vertical profile of temperature at current simulation time (K)
    real(kind=dp),intent(in)                          :: tsoilk


    !canopy extinction coefficient for ellipitical leaf angle distribution
    kbexza = kbe(x,zarad)

    !determine downwelling sky longwave radiation
    SKYLW = MEASURED
    if (SKYLW == ESTIMATED) then
        !effective downwelling sky thermal radiation
        w = 4.65*eatm/tk(npts)
        epssky = 1.0 -(1+w)*dexp(-(1.2+3.0*w)**0.5)
        s = (((ppfd_direct +ppfd_diffus)/4.6) + nir_direct+nir_diffus)/1000.0
        epssky = (1.0-s) + s*epssky
    else
        !directly from measured downwelling LW radiation
        epssky = lwdnzref/(sbsig*tkzref**4.0)
    end if

    sky_ir = epssky*sbsig*tkzref**4.0

    !longwave upwelling radiation
    grnd_ir = epsgrnd*sbsig*tsoilk**4.0
    lw_up(1) = grnd_ir
    do i=2,npts
        if (lai(i) >0.0) then
            taudt_ir = dexp(-kd*lai(i))
            lw_up(i) = lw_up(i-1)*taudt_ir + (1.0-taudt_ir)*epsleaf*sbsig*tk(i-1)**4.0
        else
            lw_up(i) = lw_up(i-1)
        end if
    end do

    !diffuse extinction coefficients and canopy reflectances
    kdsig_par = kd/dsqrt(1.0-sigma_par)
    kdrfl_par = kd/dsqrt(1.0-rfl_par)
    kdtrs_par = kd/dsqrt(1.0-trns_par)
    rho_par   = ((1.0-dsqrt(1.0-sigma_par))/(1.0+dsqrt(1.0-sigma_par)))*(2.0/(1.0+1.6*dcos(zarad)))

    kdsig_nir = kd/dsqrt(1.0-sigma_nir)
    kdrfl_nir = kd/dsqrt(1.0-rfl_nir)
    kdtrs_nir = kd/dsqrt(1.0-trns_nir)
    rho_nir   = ((1.0-dsqrt(1.0-sigma_nir))/(1.0+dsqrt(1.0-sigma_nir)))*(2.0/(1.0+1.6*dcos(zarad)))

    !loop over the domain from the top
    do i=npts,1,-1

   
      ! within or above canopy?
      if (i > ncnpy) then         !above canopy is all sunny
          ppfd_sun(i) = ppfd_direct +ppfd_diffus
          ppfd_shd(i) = 0.0
          nir_sun(i)  = nir_direct + nir_diffus
          nir_shd(i)  = 0.0
          rt_sun(i)   = ppfd_sun(i)/4.6+nir_sun(i)
          rt_shd(i)   = 0.0

          ! day or night?
          if (rt_sun(i) >0.0) then
              fsun(i) = 1.0
              fshd(i) = 0.0
          else
              fsun(i) = 0.0
              fshd(i) = 1.0
          end if

          ! no canopy, so no absorption
          rabs_sun(i) = 0.0
          rabs_shd(i) = 0.0

          lw_dn(i) = sky_ir
      else      !within canopy
      
          ! day or night?
          if (rt_sun(ncnpy+1) >0.0) then
              fsun(i) = dexp(-kbexza*clai(i))
              fshd(i) = 1.0 - fsun(i)
          else
              fsun(i) = 0.0
              fshd(i) = 1.0
          end if

          ! downwelling radiation
          if (lai(i) >0.0) then
              taudt_ir = dexp(-kd*lai(i))
              lw_dn(i) = lw_dn(i+1)*taudt_ir + (1.0-taudt_ir)*epsleaf*sbsig*tk(i+1)**4.0
          else
              lw_dn(i) = lw_dn(i+1)
          end if

          ! diffuse radiation
          rdpar = ppfd_diffus*(1.0 - rho_par)*dexp(-kd*clai(i))
          rdnir = nir_diffus*(1.0 - rho_nir)*dexp(-kd*clai(i))

          ! upwelling scattered radiation
          rsuppar = ppfd_direct*rfl_par*(dexp(-kbexza*clai(i)) - dexp(kd*clai(i)-(kbexza+kd)*clai(1)))/(kd+kbexza)
          rsupnir = nir_direct*rfl_nir*(dexp(-kbexza*clai(i))-dexp(kd*clai(i)-(kbexza+kd)*clai(1)))/(kd+kbexza)

          ! downwelling scattered radiation
          rsdnpar = ppfd_direct*trns_par*(dexp(-kbexza*clai(i))-dexp(-kd*clai(i)))/(kd-kbexza)
          rsdnnir = nir_direct*trns_nir*(dexp(-kbexza*clai(i))-dexp(-kd*clai(i)))/(kd-kbexza)

          ! sunlit leaves
          ppfd_sun(i) = kdsig_par*rdpar + kdrfl_par*rsuppar +kdtrs_par*rsdnpar + kbexza*ppfd_direct
          nir_sun(i)  = kdsig_nir*rdnir + kdrfl_nir*rsupnir +kdtrs_nir*rsdnnir + kbexza*nir_direct
          rt_sun(i)   = alfpar*ppfd_sun(i)/4.6 + alfnir*nir_sun(i)
          rabs_sun(i) = rt_sun(i) + alfir*(lw_up(i) + lw_dn(i))
      
          ! shaded leaves
          ppfd_shd(i) = kdsig_par*rdpar +kdrfl_par * rsuppar + kdtrs_par*rsdnpar
          nir_shd(i)  = kdsig_nir*rdnir +kdrfl_nir * rsupnir + kdtrs_nir*rsdnnir
          rt_shd(i)   = alfpar*ppfd_shd(i)/4.6 + alfnir*nir_shd(i)
          rabs_shd(i) = rt_shd(i) + alfir*(lw_up(i) +lw_dn(i))
      end if
    end do
  
    return
end subroutine CalcRadProfiles

!===========================================================================
!function kbe - canopy extinction coefficient
!
!Use the ellipsoidal leaf distribution function of ...
!    Campbell and Norman (1998) An introduction to Environmental Biophysics,
!    2nd Ed., Springer, New York. (p. 251)
!===========================================================================
function kbe(x, za)
    implicit none
    integer, parameter :: dp=kind(0.d0)
    integer, parameter :: i4=kind(1)
    real(kind=dp), intent(in)   :: x   !leaf angle distribution parameter  (canopy specific)
    real(kind=dp), intent(in)   :: za  !zenith angle (radians)
    real(kind=dp)               :: kbe   !canopy extinction coefficient
    real(kind=dp)               :: numer   !tmps
    real(kind=dp)               :: denom   !tmps

    numer = (x*x +dtan(za)*dtan(za))**0.5
    denom = x + 1.774*((x + 1.182)**(-0.733))

    kbe = numer/denom
  
    return
end function kbe

!=================================================================================================
!Subroutine PartitionRAD - partitions measured global and income PAR into diffuse and direct beams
!
!Uses the algorithm of ...
!    Weiss & Norman (1985) Ag. & Forest Meteor., 34, 205-213
!with corrections included in CANOAK by Baldocchi 
!    Eq (3) of Weiss & Norman corrected to
!        Rdv = 0.4*(600*cos*theta) - RDV)
!    Eq (5) of Weiss & Norman corrected to
!        Rdv = 0.6*(720-RDN/cos(theta) -w)*cos(theta)
!
!    and adjusting all of the equations as neccessary for a solar constant value of
!    1373 W/m^2 instead of the 1320 W/m^2 used by Weiss & Norman
!==================================================================================================
subroutine PartitionRAD(zarad,pmbzref,ppfdzref,sradzref,ppfd_direct,ppfd_diffus,nir_direct,nir_diffus)
    implicit none
    integer, parameter  :: dp=kind(0.d0)
    integer, parameter  :: i4=kind(1)
    real(kind=dp)       :: cosza     ! cosine of zenith angle
    real(kind=dp)       :: rdpar     ! potential direct beam PAR (W/m^2)
    real(kind=dp)       :: rspar     ! potential diffuse PAR (W/m^2)
    real(kind=dp)       :: rdnir     ! potential direct beam NIR (W/m^2)
    real(kind=dp)       :: rsnir     ! potential diffuse NIR (W/m^2)
    real(kind=dp)       :: rtpar     ! potential total PAR (W/m^2)
    real(kind=dp)       :: rtnir     ! potential total NIR (W/m^2)
    real(kind=dp)       :: wa        ! water absoption in NIR
    real(kind=dp)       :: radratio  ! ratio of actual to potential total radiation
    real(kind=dp)       :: rr        ! temp variable for radratio
    real(kind=dp)       :: nirx      ! IR inferred from sradzref and ppfdzref (W/m^2) 
    real(kind=dp)       :: rfac      ! computation variable 
    real(kind=dp)       :: fdpar     ! calculated direct beam fraction of PAR
    real(kind=dp)       :: fspar     ! calculated diffuse fraction of PAR
    real(kind=dp)       :: fdnir     ! calculated direct beam fraction of NIR
    real(kind=dp)       :: fsnir     ! calculated diffuse fraction of NIR
    real(kind=dp)       :: pc        ! pressure and optical factor
    real(kind=dp), parameter  :: coszamin=0.017    ! cos(89deg)
    real(kind=dp),intent(in)       :: pmbzref
    real(kind=dp),intent(out)       :: ppfd_direct
    real(kind=dp),intent(out)       :: ppfd_diffus
    real(kind=dp),intent(out)       :: nir_direct
    real(kind=dp),intent(out)       :: nir_diffus
    real(kind=dp),intent(in)       :: ppfdzref
    real(kind=dp),intent(in)       :: sradzref
    real(kind=dp),intent(in)       :: zarad     ! zenith angle at current meteorological data time step (radians)

    cosza=cos(zarad)
    if (cosza < coszamin) then
        ppfd_direct = 0.0
        ppfd_diffus = 0.0
        nir_direct = 0.0
        nir_diffus = 0.0
        return
    end if

    pc = (pmbzref/1013.25)/cosza

    ! potential direct beam PAR
    rdpar = 624.0*exp(-0.185*pc)*cosza

    ! potential diffuse PAR
    rspar = 0.4*(624.0*cosza-rdpar)

    ! water absorption in NIR for 10 mm precipitable water
    ! (expression used here is that used in CANOAK, but is different than that given in Weis & Norman)
    wa = 1773.0* 0.077*(2.0*pc) ** 0.30

    ! potential direct beam NIR
    rdnir = (748.0*exp(-0.06*pc) - wa)* cosza

    ! potential diffuse NIR
    rsnir = 0.6*(748.0 - (rdnir/cosza)- wa)*cosza

    ! total potential PAR &NIR
    rtpar = rdpar + rspar
    rtnir = rdnir + rsnir

    ! ratio of observed and potential total radiation
    radratio = sradzref/(rtpar + rtnir)

    !measured IR
    nirx = sradzref - (ppfdzref/4.6)

    ! PAR direct and diffuse (umol/m^2-s)
    if (radratio >0.9) then
        rr = 0.9
    else
        rr = radratio
    end if
    rfac = 1.0 - ((0.9 - rr)/0.7)**0.667
    fdpar = (rdpar/rtpar)*rfac

    if (fdpar <0.0) then
        fdpar = 0.0
    end if

    if (fdpar >1.0) then
        fdpar = 1.0
    end if

    fspar = 1.0 - fdpar

    ppfd_direct = fdpar*ppfdzref
    ppfd_diffus = fspar*ppfdzref

    if (ppfd_direct < 0.0) then
        ppfd_direct = 0.0
        ppfd_diffus = ppfdzref
    end if

    if (ppfdzref == 0.0) then
        ppfd_direct = 0.001
        ppfd_diffus = 0.001
    end if

    ! NIR direct and diffuse (W/m^2)
    if (radratio > 0.88) then
        rr = 0.88
    else
        rr = radratio
    end if
    rfac = 1.0 - ((0.88 - rr)/0.68)**0.667
    fdnir = (rdnir/rtnir)*rfac

    if (fdnir <0.0) then
        fdnir =0.0
    end if

    if (fdnir >1.0) then
        fdnir = 1.0
    end if

    fsnir = 1.0 -fdnir

    nir_direct = fdnir * nirx
    nir_diffus = fsnir * nirx

    if (nir_direct < 0.0) then
        nir_direct = 0.0
        nir_diffus = nirx
    end if

    if (nirx <= 0.0) then
        nir_direct = 0.1
        nir_diffus = 0.1
    end if

    return
end subroutine PartitionRAD

!===============================================================================================================
!SolarZenithAngle .. calculate the topocentric (i.e., local) zenith angle (in degrees) based on the algorithm of 
!                    Michalsky,J.J. (1988) 
!                    for elevation refraction correction, algorithm assumes surface T=288K and P=1013.25mb
!                    
!                    input parameters:
!                        time21 = modified Julian Day from Noon on Jan 1, 2000
!                        slat=latitude in degree (north is positive)
!                        slon=longitude in degree (east is positive)
!      
!                    internal parameters:
!                        el= sun elevation angle(degs)
!                        ha = solar hour angle
!                        dec = declination
!
!                    return
!                        za= zenith angle (degs)
!================================================================================================================
function SolarZenithAngle(time21,slat,slon)
    real(kind=dp), intent(in)  :: time21 !time in days since Jan 1, 2000
    real(kind=dp), intent(in)  :: slat   !latitude (deg)
    real(kind=dp), intent(in)  :: slon   !longitude (deg)
    real(kind=dp)  :: SolarZenithAngle
    real(kind=dp)  :: twopi
    real(kind=dp)  :: pi
    real(kind=dp)  :: deg2rad
    real(kind=dp)  :: mnlong
    real(kind=dp)  :: mnanom
    real(kind=dp)  :: eclong
    real(kind=dp)  :: oblqec
    real(kind=dp)  :: num
    real(kind=dp)  :: den
    real(kind=dp)  :: ra
    real(kind=dp)  :: gmst
    real(kind=dp)  :: lmst
    real(kind=dp)  :: latrad
!    real(kind=dp)  :: elc
    real(kind=dp)  :: refrac
!    real(kind=dp)  :: soldia
    real(kind=dp)  :: el
    real(kind=dp)  :: dec
    real(kind=dp)  :: ha
    

    pi = 2.0*acos(0.0)
    twopi = 2.0*pi
    deg2rad = pi/180.0

    !force mean longitude between 0 and 360 degs
    mnlong=280.460+0.9856474*time21
    mnlong=mod(mnlong,360.0)
    if (mnlong.lt.0.0) then
        mnlong = mnlong+360.0
    end if

    !mean anomaly in radians between 0 and 2*pi
    mnanom=357.528+0.9856003*time21
    mnanom=mod(mnanom,360.0)
    if (mnanom.lt.0.0) then
        mnanom = mnanom+360.0
    end if
    mnanom=mnanom*deg2rad
    
    !compute the ecliptic longitude and obliquity of ecliptic in radians
    eclong=mnlong+1.915*sin(mnanom)+0.020*sin(2.0*mnanom)
    eclong = mod(eclong,360.0)
    if (eclong.lt.0.0) then
        eclong=eclong+360.0
    end if 
    oblqec = 23.439 - 0.0000004*time21
    eclong = eclong*deg2rad
    oblqec = oblqec*deg2rad
  
    !calcualte right ascension and declination
    num = cos(oblqec)*sin(eclong)
    den = cos(eclong)
    ra = atan(num/den)

    !force ra between 0 and 2*pi
    if (den.lt.0.0) then
        ra=ra+pi
    elseif (num.lt.0.0) then
        ra=ra+twopi
    end if

    !dec in radians
    dec = asin(sin(oblqec)*sin(eclong))

    !calculate Greenwich mean sidereal time in hours
    !subtitute the appoximate gmst formula from the U.S. Naval Observatory web site
    !http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time
    gmst = 18.697374558 + 24.06570982441908 * time21
    gmst = mod(gmst,24.0)
    if (gmst.lt.0.0) then
        gmst=gmst+24.0
    end if

    !calcualte local mean sidereal time in radians
    lmst = gmst+slon/15.0
    lmst = mod(lmst,24.0)
    if (lmst.lt.0.0) then
        lmst=lmst+24.0
    end if
    lmst = lmst*15.0*deg2rad

    !calculate hour angle in radians between -pi and pi
    ha = lmst-ra
    if (ha.lt.-pi) then
        ha=ha+twopi
    end if
    if (ha.gt.pi) then
        ha=ha-twopi
    end if

    !change latitude to radians
    latrad=slat*deg2rad

    !calculate elevation
    el=asin(sin(dec)*sin(latrad)+cos(dec)*cos(latrad)*cos(ha))

    !calcualte refraction correction for US stan. atmosphere
    !need to have el in degs before calculating correction
    el=el/deg2rad

    !note that 3.51823=1013.25 mb/288K
    if (el.ge.19.225) then
        refrac=0.00452*3.51823/tan(el*deg2rad)
    else if (el.gt.-0.766.and.el.lt.19.225) then
        refrac=3.51823*(0.1594+0.0196*el+0.00002*el**2)/(1.0+0.505*el+0.0845*el**2)
    else if (el.le.-0.766) then
        refrac = 0.0
    end if

    !elevation in degs
    el = el+refrac
 
    !zenith angle in degs
    SolarZenithAngle = 90.0-el

end function SolarZenithAngle 


end module ACCESS_Modules


