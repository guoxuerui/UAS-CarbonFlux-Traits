Module Plants

  Use Datatypes
  Use TimeData, only: TimeFactor
  Implicit None

  Integer, parameter :: WinterWheat = 1   ! Winterweizen
  Integer, parameter :: SummerWheat = 2   ! Sommerweizen
  Integer, parameter :: Maize       = 3   ! Mais
  Integer, parameter :: Potatoe     = 4   ! Kartoffel
  Integer, parameter :: SugarBeet   = 5   ! Zuckerruebe
  Integer, parameter :: Gras        = 6   ! C3 Gras
  Integer, parameter :: Beech       = 7   ! Buche
  Integer, parameter :: Spruce      = 8   ! Fichte
  Integer, parameter :: Oak         = 9   ! Eiche
  Integer, parameter :: MixedForest =10   ! Mischwald
  Integer, parameter :: Alder       =11   ! Erle
  Integer, parameter :: Barley      =12   ! Winter Gerste
  Integer, parameter :: c4Gras      =13   ! C4 Gras
  Integer, parameter :: MaxPlant    =13
  Character(16), Parameter :: PlantName(MaxPlant) = (/ &
       'Winter Wheat    ',&
       'Summer Wheat    ',&
       'Maize           ',&
       'Potatoes        ',&
       'Sugar Beet      ',&
       'C3 Gras         ',&
       'Beech           ',&
       'Spruce          ',&
       'Oak             ',&
       'Mixed Forest    ',&
       'Alder           ',&
       'Barley          ',&
       'C4 Gras         '/)
  ! factor for abovegroundResidues
  Real(dp) :: agr_fact(MaxPlant)=(/ 0.25,0.25,0.1,1.0,1.0, &
       1.0,1.0,1.0,1.0,1.0,1.0,0.25,1.0 /)
  Integer :: PlantTabs=12
  !raus Integer, Parameter :: SimplePlantTabs=6
  Integer, Parameter :: RootTab=12        ! distribution of roots
  !raus Integer, Parameter :: SimpleRootTab=5   ! distribution of roots
  !raus Real(dp) :: dtPlants
  Logical :: harvestresidues
  Real(dp) :: senescence_start_temp ! only used for potatoe (Grad C)

  ! climate data
  Type ClimateDataType
     Real(dp) :: time
     Real(dp) :: et0          ! potential evapotranspiration of reference crop
     Real(dp) :: prec         ! precipitation
     !Real(dp) :: intc
     Real(dp) :: temp
     Real(dp) :: temp_avg
     Real(dp) :: temp_min     ! daily min temperature
     Real(dp) :: temp_max     ! daily max temperature
     Real(dp) :: global_rad   ! "timestep" total global radiation
  End Type ClimateDataType
  Type(ClimateDataType) :: ClimateData

  Type FarquharType
     Logical :: used
     Real(dp) :: w_c, w_j_sl, w_j_sh, w_e   ! ( μ mol CO_2 ) /  m**2 / s 
     Real(dp) :: vcmax        ! maximum rate of carboxylation
     Real(dp) :: vcmax25      ! maximum rate of carboxylation at 25° C
     Real(dp) :: vcmax25N     ! vcmax25 * f(N)   nitrogen availability factor f(N)
     Real(dp) :: alfa         ! quantum efficiency ( µmol CO2 per µmol photons)
     Real(dp) :: ftv          ! function that mimics thermal breakdown of metabolic processes
     Real(dp) :: m            ! Ball/Berry coefficient m (dimensionless)
     Real(dp) :: c_i          ! internal leaf CO_2 partial pressure (Pa)
     Real(dp) :: c_s          ! CO_2 partial pressure at leaf surface (Pa)
     Real(dp) :: o_i          ! O_2 partial pressure (Pa)
     Real(dp) :: gamma        ! CO_2 compensation point gamma (Pa)
     Real(dp) :: K_c, K_o     ! K_c and K_o are the Michaelis-Menten constants (Pa) for CO_2 and O_2
     Real(dp) :: K_c25, K_o25 ! values (Pa) at 25Grad
     Real(dp) :: a_kc, a_ko   ! relative changes in K_c25 and K_o25 , respectively, for a 10Grad change in temperature
     Real(dp) :: temp_fw      ! freezing temperature of water
     Real(dp) :: Boltzmann    ! J / K / molecule
     Real(dp) :: Avogadro     ! molecule / kmol
     Real(dp) :: R_gas        ! J / K / kmol
     Real(dp) :: CN_L         ! CN_L is the leaf carbon-to-nitrogen ratio (g C / g N)
     Real(dp) :: F_LNR        ! fraction of leaf nitrogen in Rubisco (g N in Rubisco / g N)
     Real(dp) :: F_NR         ! mass ratio of total Rubisco molecular mass to nitrogen in Rubisco (g Rubisco / g N in Rubisco)
     Real(dp) :: a_R25        ! specific activity of Rubisco (µmol CO2 / g Rubisco / s)
     Real(dp) :: N_a          ! area-based leaf nitrogen concentration (g N / m**2 leaf area)
     Real(dp) :: daylength_max! max daylength for latitude
     Real(dp) :: SLA_0        ! value for SLA at the top of the canopy (m**2 leaf area / g C)
     Real(dp) :: PSI_O        ! soil water potential (mm) when stomata are fully open
     Real(dp) :: PSI_C        ! soil water potential (mm) when stomata are fully closed
  End Type FarquharType
  Logical :: Farquhar_used

  ! plant data for all plants
  Type PlantDataType
     Integer :: NoTypes          ! Number of plants
     Real(dp) :: dayOfYear       ! day of the year
     Real(dp) :: latitude        ! latitude of the site
     Integer :: InterceptModel   ! 1=Bormann  2=Hoyningen-Huene
     Real(dp) :: UnitFactor      ! factor to convert from trace (mm) to plants: 1 for mm, 10 for cm, ...
     Real(dp) :: Area_ha_to_L2   ! factor to convert from ha to L2
     Real(dp) :: Area_L2_to_m2   ! factor to convert from L2 to m2
     Real(dp) :: Unit_L_to_m     ! factor to convert from L to m
     Real(dp) :: root_max        ! max rooting depth for all plants
     Real(dp) :: wcwp_phead      ! pressure head for water content at willting point
  End Type PlantDataType
  Type(PlantDataType) :: PlantData

  ! plant data for each plant type for all models
  Type PlantType
     Integer :: type      ! plant type
     Integer :: season    ! ==0 if plant is not between planting and harvest
     Logical :: c3 ! true for c3 plants
     ! pressure values p0 < p1 <p2h < p2l < p3
     Real(dp) :: p0, p1, p2h, p2l, p3
     ! threshold collar water potential, Equivalent conductance, Compensatory RWU conductance
     Real(dp) :: hx_min, Krs, Kcomp
     Real(dp) :: rna       ! actual depth above there is no root water uptake
     Real(dp) :: so_length ! depth of storage organ for underground so
     Type(FarquharType), Pointer :: Far
     Type(TabType) :: root_distribution   ! relative root depth against root density
  End Type PlantType
  Type(PlantType), Allocatable :: Plant(:)

  ! plant data for the Sucros plant model
  Type SucrosPlantType
     Integer :: akctype      ! actual Kc against 1=dvs or 2=time (day of the year) 3=Kc calculation
     Real(dp) :: akcMin      ! minimum of Kc for akctype=3
     Real(dp) :: Tmin_v1     ! min temperature for development, dependent on developmental phase, EM-TS
     Real(dp) :: Topt_v1     ! optimum temperature for development, dependent on developmental phase, EM-TS
     Real(dp) :: Tmax_v1     ! max temperature for development, dependent on developmental phase, EM-TS
     Real(dp) :: Tmin_v2     ! min temperature for development, dependent on developmental phase, TS-AN
     Real(dp) :: Topt_v2     ! optimum temperature for development, dependent on developmental phase, TS-AN
     Real(dp) :: Tmax_v2     ! max temperature for development, dependent on developmental phase, TS-AN
     Real(dp) :: Tmin_r      ! min temperature for development, dependent on developmental phase, AN-PM
     Real(dp) :: Topt_r      ! optimum temperature for development, dependent on developmental phase, AN-PM
     Real(dp) :: Tmax_r      ! max temperature for development, dependent on developmental phase, AN-PM
     Real(dp) :: Tmin_vn     ! min temperature for development, dependent on developmental phase, vernalization
     Real(dp) :: Topt_vn     ! optimum temperature for development, dependent on developmental phase, vernalization
     Real(dp) :: Tmax_vn     ! max temperature for development, dependent on developmental phase, vernalization
     Real(dp) :: Popt        ! ==0 or for short- or long-day plants ==optimal photoperiod -> different calculation of omega
     Real(dp) :: Pcrit       ! critical photoperiod below which no development occurs, cultivar dependent
     Real(dp) :: omega       ! photoperiod sensitivity coefficient, cultivar dependent
     Real(dp) :: rmax_v1     ! max daily development rate in EM-TS phase, cultivar dependent
     Real(dp) :: rmax_v2     ! max daily development rate in TS-AN phase, cultivar dependent
     Real(dp) :: rmax_r      ! max daily development rate in AN-PM phase, cultivar dependent
     Real(dp) :: rna_max     ! max depth above there is no root water uptake
     Real(dp) :: root_max    ! max rooting depth
     Real(dp) :: root_init   ! initial rooting depth
     ! vaiables used for plant model Sucros
     Integer :: no_dates     ! no of dates for planting and harvest
     ! start/end time for the reduction of extraction because of senescence
     Real(dp) :: senescence_start, senescence_end
     Real(dp) :: nsl         ! number of seedlings per m**2
     Real(dp) :: rgr         ! relative growth rate during exponential leaf area growth
     Real(dp) :: tempbase    ! base temperature for juvenile leaf area growth
     Real(dp) :: sla         ! specific area of new leaves
     Real(dp) :: rsla        ! change of specific leaf area per unit thermal time
     Real(dp) :: amx         ! potential CO2-assimilation rate of a unit leaf area for light saturation
     Real(dp) :: eff         ! initial light use efficiency
     Real(dp) :: rkdf        ! extinction coefficient for diffuse PAR
     Real(dp) :: scp         ! scattering coefficient of leaves for PAR
     Real(dp) :: rmainso     ! maintenance demand rate for storage organs per unit dry matter
     Real(dp) :: asrqso      ! conversion efficiency coefficient (assimilate requirement for DM for storage organs)
     !Real(dp) :: ear        ! ear area ratio (HA EAR/SHOOT) (not used)
     Real(dp) :: tempstart   ! start temperature for the plant growth
     Real(dp) :: debr_fac    ! dead LAI debris factor (added by D.Farber)
     Real(dp) :: ls          ! LAI as switch from temperature to radiation-limited LAI expansion
     Real(dp) :: rlaicr      ! critical LAI for leaf death due to self shading
     Real(dp) :: eai         ! initial value of the ear area index (2sided) (crop 1-3,5)
     Real(dp) :: rmatr       ! initial value of the maturity class (crop 4)
     Real(dp) :: ssl         ! leaf area of one seedling
     Real(dp) :: srw         ! specific root weight
     Real(dp) :: slaid_off   ! slaid for outside the season
     Integer, Pointer :: emergence(:)       ! day of crop emergence
     Integer, Pointer :: harvest(:)         ! harvest date
     Integer :: lCeres                      ! 0=original or 1=new approach for simulation of dvs (added by D.Farber)
     Real(dp) :: rlv, rst, rso, rcrn, rrt   ! respiration for each plant organ (growth + maintenance)
     Real(dp) :: rgrowthlv, rgrowthst, rgrowthso, rgrowthcrn, rgrowthrt   ! growth respiration for each plant organ
     Real(dp) :: deathfacMax                ! max factor used for factor of rootdeath
     Real(dp) :: dvsReset, dvsLow, dvsHigh, storage, storagerat, storageFacMax, storageMax, Tstorage, cutLai, temp_sum_crit   ! LINGRA (added by D.Farber)
     Real(dp) :: exu_fac                    !added by D.Farber
     Real(dp) :: bp(3)                      ! phosphorus uptake parameter
     Integer :: harvestType                 ! type of last harvest in year: 1=continuous grassland  2=remove plant
     Logical, Allocatable :: last_Date_in_year(:) ! only used for grass
     Type(TabType), Allocatable :: tab(:)
  End Type SucrosPlantType
  Type(SucrosPlantType), Allocatable :: SucrosPlant(:)

  Type PlantGeometryType
     Integer :: index        ! index for plant, 0 for no plant
     Integer :: output       ! index for output of plant data
     Real(dp) :: rootdepth   ! actual length of roots
     Real(dp) :: akc         ! actual crop factor Kc
     Real(dp) :: p2          ! pressure value that depends on transpiration/reduction factor of the potential root water uptake
     !raus   Real(dp) :: area        ! area of the surface element side
     ! accumulative variables
     Real(dp) :: temp_sum    ! effective temperature sum
     Real(dp) :: temp_sum_till   ! effective temperature sum
     !Real(dp) :: temp_eff_sen
     Real(dp) :: dvs         ! plant development stage
     Real(dp) :: VD          ! effective vernalization days
     Real(dp) :: cgphot      ! cumulative glucose assimilation
     Real(dp) :: pond_plant  ! interception
     Real(dp) :: pond_soil   ! ponding on soil
     Real(dp) :: pond_plant_new   ! actual changing of interception
     Real(dp) :: rain        ! precipitation, later rain on soil
     ! total leaf area lai = slaid + slaig
     Real(dp) :: slaig       ! total green leaf area
     Real(dp) :: slaid       ! total dead leaf area
     ! weight in the different plant organs (DM)
     ! weight of the leaves wlv = wlvg + wlvd
     Real(dp) :: wlvg        ! dry weight of green leaves per unit area
     Real(dp) :: wlvd        ! dry weight of dead leaves per unit area
     Real(dp) :: wso         ! dry weight of storage organs per unit area
     Real(dp) :: wst         ! dry weight of stems per unit area
     Real(dp) :: wrt         ! dry weight of roots per unit area
     Real(dp) :: wrtd        ! dry weight of dead roots per unit area
     Real(dp) :: ratwrtd        ! dry weight of dead roots per unit area
     Real(dp) :: wcrn        ! dry weight of the crown (crop 5: sugar beets) per unit area
     Integer :: count        !(added by A.Klosterhalfen)
     ! end accumulative variables
     ! surface flux related stuff
     Real(dp) :: eva_soil        ! potential soil evaporation
     Real(dp) :: eva_plants      ! potential evapotranspiration of plants (transpiration + evaporation of interception)
     Real(dp) :: eva_soil_intercept   ! evaporation of soil minus evaporation of interception
     Real(dp) :: transpiration   ! potential trancpiration (eva_plants - evaporation of interception)
     Real(dp) :: transpiration_act   ! actual transpiration
     Real(dp) :: flux            ! flux=precipitation-evaporation
     Real(dp) :: fluxsum         ! accumulated flux over 1 day
     Real(dp) :: z0              ! z value of the surface
     Real(dp) :: plai            ! 
  End Type PlantGeometryType
  Type(PlantGeometryType) :: PlantGeometry

  Type PlantGeometryArrayType
     ! arrays with dimension NumNP
     ! numbering starts from top, not from bottom!
     Real(dp), Pointer :: z(:)         ! z value, z from 0 to -drz
     Real(dp), Pointer :: dz(:)        ! z value, z from 0 to -drz
     Real(dp), Pointer :: rrd(:)       ! relative root distribution (added by D.Farber)
     Real(dp), Pointer :: rsod(:)      ! relative distribution for storage organs in soil
     Real(dp), Pointer :: rex(:)       ! potential root extraction rate
     Real(dp), Pointer :: rex_act(:)   ! actual root extraction rate
     Real(dp), Pointer :: wcwp(:)      ! water content at willting point
  End Type PlantGeometryArrayType
  Type(PlantGeometryArrayType) :: PlantGeometryArray

  ! soilco2 coupling.
  Real(dp), Parameter :: MolCO2=0.0440098 ! kg/mol
  ! nitrogen coupling
  Real(dp) :: NitReduct=1, PhoReduct=1
  ! local subroutines
  Private :: between, CalculateLAI, SetPlantTime, ResetPlant, NormalizeRootDistribution

Contains

  Subroutine s2bool(s, b)
    Character(len=1), Intent(in) :: s
    Logical, Intent(out) :: b
    If(s=='t' .Or. s=='T') Then
       b=.True.
    Else
       b=.False.
    Endif
  End Subroutine s2bool

  Subroutine s2int(s, i)
    Character(len=1), Intent(in) :: s
    Integer, Intent(out) :: i
    If(s == 'f' .Or. s == 'F') Then
       i=0
    Else If(s == 't' .Or. s == 'T') Then
       i=1
    Else
       i=Index('0123456789', s)-1
    End If
  End Subroutine s2int
  
  Subroutine ApplyRootSourceSink(wc, sink, hRoot,vRoot,cRoot)
    Use Geometry, Only: NumNP,hNew,coord,coord_eps, zSurf, hTot
    Use Timedata, Only: SimTime
    Use Variables, Only: lCO2, waterstress, alphaAvg, PlantsExist
    Use Carbon, Only: CO2
    Implicit None
    Real(dp), Intent(in)  :: wc(:)
    Real(dp), Intent(out) :: sink(:),hRoot,vRoot,cRoot
    Integer :: e,ep, ip
    Real(dp) :: rex, rexpot, dx, ARoot,act_depth, rexarea, rexareaPot
    Real(dp) :: alphaRel, phwc, minus
    Real(dp) :: hSeq, Tpot, PHcollar, Tact, Kcomp

    ! reduce the root extraction term for the plants for pressure
    ! near to saturation (anaerobic) or for dry conditions (water stress).
    If(.not.PlantsExist) Return
    ip=PlantGeometry%index
    If(ip==0) Return
    sink=0
    vRoot=0.
    hRoot=0.
    cRoot=0.
    act_depth=0.
    alphaRel = 0.
    rexarea = 0.
    rexareaPot = 0.
    rex=0
    If(waterstress==2) Then
       minus=-1.0
    Else
       minus=1.0
    End If
    PlantGeometryArray%rex_act=null_dp
    If(waterstress<3) Then
       Do ep=1,NumNP
          e=NumNP-ep+1   ! for plants the elements are numbered from top to bottom
          If(e==1) Then
             dx=0.5*(coord(2)-coord(1))
          Elseif(e==NumNP) Then
             dx=0.5*(coord(NumNP)-coord(NumNP-1))
          Else
             dx=0.5*(coord(e+1)-coord(e-1))
          Endif
          rexpot=PlantGeometryArray%rex(ep)
          If(rexpot/=null_dp) Then
             rexAreaPot = rexareaPot + rexpot*dx
             ! average pressure head (ph) in element
             If(waterstress==0) Then
                rex = rexpot
             Else
                If(waterstress==2) Then
                   phwc=wc(e)
                Else
                   phwc=PlantData%UnitFactor*hNew(e)
                End If
                If(between(phwc, Plant(ip)%p0, Plant(ip)%p3)) Then
                   ! rex=rex * a
                   ! with a: 0-1 between p0 and p1
                   !           1 between p1 and p2
                   !         1-0 between p2 and p3
                   If(between(phwc, Plant(ip)%p0, Plant(ip)%p1)) Then
                      ! linear interpolation between p0 and p1
                      rex=rexpot * (Plant(ip)%p0-phwc)/(Plant(ip)%p0-Plant(ip)%p1)
                   Else If(between(phwc, PlantGeometry%p2, Plant(ip)%p3)) Then
                      ! interpolation between p2 and p3
                      rex=rexpot * 10**(minus*(PlantGeometry%p2-phwc)/Plant(ip)%p3)
                   Else
                      rex = rexpot
                   End If
                End If ! between p0 and p3
             End If ! waterstress==0
             ! Root extraction rate is limited to avoid
             ! extraction below wilting point (WCWP)
             rex=Max(null_dp, &
                  Min( rex, (wc(e)-PlantGeometryArray%wcwp(ep))/SimTime%delta_t ) )
             PlantGeometryArray%rex_act(ep)=rex
             ! modify source/sink
             sink(e)=rex
             rexarea = rexarea + rex*dx
             vRoot=vRoot+sink(e)*dx
             hRoot=hRoot+hNew(e)*dx
             If(lCO2) cRoot=cRoot+CO2(e)*dx
          Else ! rexpot==0
             !modification N.Prolingheuer (hRoot and cRoot when root extraction is zero but roots exist)
             If(PlantGeometryArray%rrd(ep) .Ne. null_dp) Then
                hRoot=hRoot+hNew(e)*dx
                If(lCO2) cRoot=cRoot+CO2(e)*dx
             End If
          End If
       End Do
    Else ! waterstress==3
       hSeq=0
       Do ep=1,NumNP
          e=NumNP-ep+1   ! for plants the elements are numbered from top to bottom
          ! hTot(n)=hNew(n)-(zSurf-coord(n))
          ! The water potential is the sum of
          ! the matric potential and the gravitational potential.
          ! The "DOT_PRODUCT" of hTot by SSF:
          hTot(e)=hNew(e)-(zSurf-coord(e))
          hSeq=hSeq+PlantGeometryArray%rrd(ep)*hTot(e)
       End Do
       Tpot=PlantGeometry%transpiration/PlantData%UnitFactor
       PHcollar=-Abs(Tpot)/Plant(ip)%Krs+hSeq
       If(PHcollar<Plant(ip)%hx_min) Then
          !We check if we are under the threshold collar water potential.
          PHcollar=Plant(ip)%hx_min            ! In case of stress, the collar water potential stays blocked at the treshold value, and will be used to calculate the actual transpiration rate.
          Tact=Plant(ip)%Krs*Max(0.0_dp, hSeq-PHcollar)   ! [L/T] Convert the actual collar water potential into the actual transpiration (Couvreur jul 2011)
          !Write(*,'(A,1P,5e13.5)') 'hx_min,hseq,PHcollar=',Plant(ip)%hx_min,hSeq,PHcollar,hSeq-PHcollar
       Else
          Tact=Abs(Tpot)             ! [L/T]
       End If
       !Write(99,*) simtime%time,Plant(ip)%hx_min,PHcollar,Tact,Tpot,hSeq
       Do ep=1,NumNP
          e=NumNP-ep+1   ! for plants the elements are numbered from top to bottom
          If(e==1) Then
             dx=0.5*(coord(2)-coord(1))
          Elseif(e==NumNP) Then
             dx=0.5*(coord(NumNP)-coord(NumNP-1))
          Else
             dx=0.5*(coord(e+1)-coord(e-1))
          Endif
          Kcomp=Plant(ip)%Kcomp*(hTot(e)-hSeq)
          Sink(e)=PlantGeometryArray%rrd(ep) * &
               (Tact + Kcomp) / (PlantGeometryArray%dz(ep)/PlantData%UnitFactor)
          rexpot=PlantGeometryArray%rex(ep)
          PlantGeometryArray%rex_act(ep)=Sink(e)
          If(rexpot/=null_dp) Then
             rexAreaPot = rexareaPot + rexpot*dx
             rexarea = rexarea + Sink(e)*dx
             vRoot=vRoot+sink(e)*dx
             hRoot=hRoot+hNew(e)*dx
             If(lCO2) cRoot=cRoot+CO2(e)*dx
          Else ! rexpot==0
             !modification N.Prolingheuer (hRoot and cRoot when root extraction is zero but roots exist)
             If(PlantGeometryArray%rrd(ep) .Ne. null_dp) Then
                hRoot=hRoot+hNew(e)*dx
                If(lCO2) cRoot=cRoot+CO2(e)*dx
             End If
          End If
       End Do
    End If ! waterstress

    ! modified by M.Herbst, sets alphaAvg to zero for pressure heads lower than p3
    If(waterstress==0 .Or. rexareaPot==0) Then
       alphaAvg = 1.0
    Else
       alphaAvg = rexarea/rexareaPot
    Endif

    ARoot=(Plant(ip)%rna-PlantGeometry%rootdepth)/PlantData%UnitFactor
    if(ARoot>coord_eps) then
       hRoot=hRoot/ARoot
       cRoot=cRoot/ARoot
    end if
  End Subroutine ApplyRootSourceSink


  Function ConvertToTraceUnit(x, unit)
    Implicit None
    Real(dp) :: ConvertToTraceUnit
    Real(dp), Intent(in) :: x
    Integer, Intent(in) :: unit
    ConvertToTraceUnit=Units_L(unit)/PlantData%UnitFactor*x
  End Function ConvertToTraceUnit


  Subroutine InitPlantData()
    Use Geometry, Only: NumNP,coord,MatNum
    Use Material, Only: FQ,Par
    Use Variables, Only: PlantsExist
    Implicit None
    Integer :: i, err, ip
    Integer :: nsucros !raus nsimple
    Real(dp) :: rootdepth, z0
    Real(dp) :: pressure

    If(.not.PlantsExist) Return
    nsucros=PlantData%NoTypes
    !raus nsimple=PlantData%NoSimpleTypes
    Do i=0,PlantData%NoTypes
       Plant(i)%season=0 ! init
    Enddo
    !raus  If(nsimple>0) Plant(nsucros+1:nsucros+nsimple)%season=1
    ! convert unit for length, for plants mm are used
    !                 mm     cm      dm       m         km
    ! from plants (mm) to trace: /UnitFactor
    ! from trace to plants (mm): *UnitFactor
    PlantData%UnitFactor=Units_L(Unit_L_Input)
    PlantData%Area_ha_to_L2 = 0.0001/(1000.0/PlantData%UnitFactor)**2   ! ha -> L2
    PlantData%Area_L2_to_m2=(1000.0_dp/PlantData%UnitFactor)**2
    PlantData%Unit_L_to_m=1000.0_dp/PlantData%UnitFactor
    ! get the max root depth for all plant types and
    rootdepth=null_dp
    If(nsucros>0) rootdepth=Minval(SucrosPlant%root_max)
    PlantData%root_max=rootdepth
    Write(*,*) 'Maximum depth of roots: ',rootdepth,' mm'
    If(Abs(rootdepth)/PlantData%UnitFactor>coord(NumNP)-coord(1)) &
         Call WriteError('Geometry smaller than maximum root depth.')
    ! normalise the root distribution curves
    Do i=1, PlantData%NoTypes
       Call NormalizeRootDistribution(Plant(i)%root_distribution)
    End Do

    ! compute the max z-value
    PlantGeometry%z0=PlantData%UnitFactor * coord(NumNP)

    ! set pressure head for water content at wilting point: -10**4.2 cm
    PlantData%wcwp_phead=ConvertToTraceUnit(-10.0_dp**4.2_dp,2)
    pressure=PlantData%wcwp_phead
    ! compute the number of elements in z-direction for the roots
    !raus   ! calculate the area top element side
    !raus   PlantGeometry%area=1 !HH
    z0=PlantGeometry%z0
    PlantGeometry%transpiration_act=null_dp
    Allocate (PlantGeometryArray%z(NumNP), stat=err)
    If(err/=0) Call AllocateError('PlantGeometryArray%z')
    Allocate (PlantGeometryArray%dz(NumNP), stat=err)
    If(err/=0) Call AllocateError('PlantGeometryArray%dz')
    Allocate (PlantGeometryArray%rex(NumNP), stat=err)
    If(err/=0) Call AllocateError('PlantGeometryArray%rex')
    Allocate (PlantGeometryArray%rrd(NumNP), stat=err)   !(added by D.Farber)
    If(err/=0) Call AllocateError('PlantGeometryArray%rrd')
    Allocate (PlantGeometryArray%rsod(NumNP), stat=err)
    If(err/=0) Call AllocateError('PlantGeometryArray%rsod')
    Allocate (PlantGeometryArray%rex_act(NumNP), stat=err)
    If(err/=0) Call AllocateError('PlantGeometryArray%rex_act')
    Allocate (PlantGeometryArray%wcwp(NumNP), stat=err)
    If(err/=0) Call AllocateError('PlantGeometryArray%wcwp')
    ! reverse the z-values and subtract height of the surface
    Do i=1,NumNP
       PlantGeometryArray%rex(i)=null_dp
       PlantGeometryArray%rrd(i)=null_dp
       PlantGeometryArray%rsod(i)=null_dp
       PlantGeometryArray%z(i)=PlantData%UnitFactor * coord(NumNP-i+1) - z0
    End Do

    ! calculate dz:
    ! first and last node
    PlantGeometryArray%dz(1)=0.5_dp*(PlantGeometryArray%z(1)-PlantGeometryArray%z(2))
    PlantGeometryArray%dz(NumNP)=0.5_dp*(PlantGeometryArray%z(NumNP-1)-PlantGeometryArray%z(NumNP))
    ! inner nodes
    Do i=2,NumNP-1
       PlantGeometryArray%dz(i)=0.5_dp*(PlantGeometryArray%z(i-1)-PlantGeometryArray%z(i+1))
    End Do

    ! calculate the water content at wilting point
    Do i=1,NumNP
       PlantGeometryArray%wcwp(i)=FQ(pressure,Par(:,MatNum(NumNP-i+1)))
    End Do

    ! output=global node number if output is wanted, 0 else
    PlantGeometry%output=NumNP
    PlantGeometry%flux=null_dp
    PlantGeometry%fluxsum=null_dp
    PlantGeometry%eva_plants=null_dp
    PlantGeometry%pond_soil=null_dp
    PlantGeometry%pond_plant=null_dp
    PlantGeometry%pond_plant_new=null_dp
    PlantGeometry%rain=null_dp
    PlantGeometry%rootdepth=null_dp



    !PlantGeometry%slaig=null_dp
    PlantGeometry%slaid=null_dp
    PlantGeometry%temp_sum=null_dp
    PlantGeometry%temp_sum_till=null_dp
    !PlantGeometry%temp_eff_sen=null_dp
    PlantGeometry%dvs=null_dp
    PlantGeometry%cgphot=null_dp
    PlantGeometry%wrt=null_dp
    PlantGeometry%wrtd=null_dp
    PlantGeometry%wso=null_dp
    PlantGeometry%wst=null_dp
    PlantGeometry%wcrn=null_dp
    PlantGeometry%slaid=null_dp
    PlantGeometry%wlvg=null_dp
    PlantGeometry%wlvd=null_dp
    PlantGeometry%index=0
    PlantGeometry%vd=null_dp
    PlantGeometry%count=0   !(added by A.Klosterhalfen)

    SucrosPlant%storage = 0.      ! LINGRA
    SucrosPlant%storagerat = 0.   ! LINGRA

    ! init akc
    Do ip=1,PlantData%NoTypes
       If(SucrosPlant(ip)%akctype==3) Then
          PlantGeometry%akc=calculate_akc(ip)
       Else
          ! set akc to first value in table
          PlantGeometry%akc=SucrosPlant(ip)%tab(11)%Data(1,2)
       Endif
    End Do
    ! convert parameters of Couvreur model from  plant to soilco2 unit
    Plant%hx_min=Plant%hx_min/PlantData%UnitFactor
    Plant%Krs=Plant%Krs*TimeFactor
    Plant%Kcomp=Plant%Kcomp*TimeFactor
  End Subroutine InitPlantData

    !modification D.Farber - harvestresidues
  Subroutine ResetPlant(ip, simTime, oldseason)
    ! Called after each season
    Use Timedata, only: simtime_to_iday, iday_to_date, resetPlantTime
    Use Variables, Only: anlv, anst, anrt, anso, ancrn, aplv, apst, aprt, apso, apcrn, aprt, ap_tot
    Use Geometry, Only: NumNP, coord, rNodeResC, rNodeResN, rNodeResP, delta_z
    Use Carbon, only: co2pdepth
    Implicit None
    Integer, Intent(in) :: ip, oldseason
    Real(dp), Intent(in) :: simTime
    Integer :: ptype, ios
    Real(dp) :: laiFactor, laiBoundary
    Real(dp) :: diffWso, diffWst, diffWcrn, diffWlvg
    character(len=*), parameter :: harvestFile = "harvest.out"
    Logical, save :: firstCall = .true.
    Real(dp) :: t
    Integer :: yy, mm, dd,i,j, init
    Real(dp) :: abovegroundResiduesC, abovegroundResiduesN, abovegroundResiduesP, depth
    Real(dp), Parameter :: wst_fact=0.494_dp
    Real(dp), Parameter :: wrt_fact=0.467_dp
    Real(dp), Parameter :: wlv_fact=0.459_dp
    Real(dp), Parameter :: wso_fact=0.446_dp

    ptype=Plant(ip)%type
    Call WriteOutput('Resetting season of plant: '//PlantName(ptype))
    
    ! reset variables used outside the season of the plant
    resetPlantTime=SimTime
    If(ptype/=Gras .And. ptype/=c4Gras) Then
       ! no lingra
       If(harvestresidues) Then
          Select Case(ptype)
          Case(WinterWheat,Barley,Maize,SummerWheat)
             abovegroundResiduesC = agr_fact(ptype)*PlantGeometry%wst*wst_fact   !kg C/L2/T
             abovegroundResiduesN = agr_fact(ptype)*anst   !kg C/L2/T
             abovegroundResiduesP = agr_fact(ptype)*apst   !kg C/L2/T
          Case(Potatoe,Sugarbeet)
             abovegroundResiduesC=(PlantGeometry%wst*wst_fact)+&
                  (PlantGeometry%wlvg+PlantGeometry%wlvd)*wlv_fact
             abovegroundResiduesN=anst+anlv
             abovegroundResiduesP=apst+aplv
          Case default
             abovegroundResiduesC=0
             abovegroundResiduesN=0
             abovegroundResiduesP=0
          End Select
          depth=coord(NumNP)-co2pdepth
          Do i=1,NumNP
             j=NumNP-i+1
             rNodeResC(j) = PlantGeometryArray%rrd(i)*PlantGeometry%wrt*wrt_fact/delta_z(j)
             rNodeResN(j) =  PlantGeometryArray%rrd(i)*anrt/delta_z(j)
             rNodeResP(j) =  PlantGeometryArray%rrd(i)*aprt/delta_z(j)
             If(coord(j).Ge.depth) Then
                rNodeResC(j) = rNodeResC(j) + abovegroundResiduesC/co2pdepth
                rNodeResN(j) = rNodeResN(j) + abovegroundResiduesN/co2pdepth
                rNodeResP(j) = rNodeResP(j) + abovegroundResiduesP/co2pdepth
             Endif
          End Do
       End If
       Call ResetPlantData(ip, 0, 2)
    Else   ! LINGRA 
       If(SucrosPlant(ip)%last_date_in_year(oldseason)) Then
          If(harvestresidues) Then
             abovegroundResiduesC=(PlantGeometry%wst*wst_fact)+&
                  (PlantGeometry%wlvg+PlantGeometry%wlvd)*wlv_fact+(PlantGeometry%wso*wso_fact)
             abovegroundResiduesN=anst+anlv+anso
             abovegroundResiduesP=apst+aplv+apso
             depth=coord(NumNP)-co2pdepth
             Do i=1,NumNP
                j=NumNP-i+1
                rNodeResC(j) = PlantGeometryArray%rrd(i)*PlantGeometry%wrt*wrt_fact/delta_z(j)
                rNodeResN(j) =  PlantGeometryArray%rrd(i)*anrt/delta_z(j)
                rNodeResP(j) =  PlantGeometryArray%rrd(i)*aprt/delta_z(j)
                If(coord(j).Ge.depth) Then
                   rNodeResC(j) = rNodeResC(j) + abovegroundResiduesC/co2pdepth
                   rNodeResN(j) = rNodeResN(j) + abovegroundResiduesN/co2pdepth
                   rNodeResP(j) = rNodeResP(j) + abovegroundResiduesP/co2pdepth
                Endif
             End Do
          Endif
          If(SucrosPlant(ip)%harvestType==1) Then
             init=1
          Else
             init=0
          Endif
          Call ResetPlantData(ip, init, SucrosPlant(ip)%harvestType)
       Else ! not last date in year
          !       rootdepth=init tempsum=0 dvs=0 cgphot=0 slaig=0 slaid=0 wlvg=wlvd=wso=wst=wrtd=0 wrt=0.5*wrt
          PlantGeometry%count=0   !(added by A.Klosterhalfen)
          laiBoundary = SucrosPlant(ip)%cutLai
          If(PlantGeometry%slaig>laiBoundary) Then

             laiFactor = PlantGeometry%slaig/laiBoundary
             PlantGeometry%slaig = laiBoundary
             PlantGeometry%slaid = 0

             diffWso = PlantGeometry%wso - PlantGeometry%wso/laiFactor
             diffWst = PlantGeometry%wst - PlantGeometry%wst/laiFactor
             diffWcrn = PlantGeometry%wcrn - PlantGeometry%wcrn/laiFactor
             diffWlvg = PlantGeometry%wlvg - PlantGeometry%wlvg/laiFactor
             if(firstCall) then
                firstCall = .false.
                Open(62,file=Trim(harvestFile),Status='REPLACE',ACTION='WRITE',&
                     FORM='FORMATTED',iostat=ios)
                If(ios/=0) Call WriteError( 'Error while opening plants output file: '&
                     // Trim(harvestFile))
                write(62,*) "wso/wst/wcrn/wlvg"
             Else
                Open(62,file=harvestFile,Status='OLD',ACTION='WRITE',&
                     FORM='FORMATTED',POSITION='APPEND')
             End If

             t=simtime_to_iday(simTime)
             Call iday_to_date(simtime_to_iday(SimTime),yy,mm,dd)
             write(62,'(I5, 2I3, F8.1, 4E13.5)') yy, mm, dd, simTime, diffWso, diffWst, diffWcrn, diffWlvg
             close(62)

             PlantGeometry%wso = PlantGeometry%wso/laiFactor
             PlantGeometry%wst = PlantGeometry%wst/laiFactor
             PlantGeometry%wcrn = PlantGeometry%wcrn/laiFactor
             PlantGeometry%wlvg = PlantGeometry%wlvg/laiFactor
             PlantGeometry%wlvd=null_dp

             If(PlantGeometry%dvs > SucrosPlant(ip)%dvsReset) Then
                PlantGeometry%dvs = SucrosPlant(ip)%dvsReset
             End If
             ! reduce N masses
             anlv=anlv/laiFactor
             anst=anst/laiFactor
             anso=anso/laiFactor
             ancrn=anso/laiFactor
             ! reduce P masses
             ap_tot=ap_tot*(aplv+apst+apso+apcrn)/(aplv+apst+apso+apcrn+aprt)/laiFactor
             aplv=aplv/laiFactor
             apst=apst/laiFactor
             apso=apso/laiFactor
             apcrn=apso/laiFactor
          Else
             Print *,'slaig: ',PlantGeometry%slaig,' less equal laiBoundary:',laiBoundary
          End If
       End If ! harvestresidues
    End If ! grass
    !end modification
  End Subroutine ResetPlant

  
  Subroutine ResetPlantData(ip, init, mode)
    Use Variables, Only: GPP, NPP, aboveground_respiration
    Use Geometry, Only: rnodert, rnodexu
    Implicit None
    Integer, Intent(in) :: ip, init, mode
    Real(dp) :: x
    Print *,"ResetGras"
    !modification N.Prolingheuer - variable akc
    If(SucrosPlant(ip)%akctype==3) Then
       PlantGeometry%akc=calculate_akc(ip)
    Else
       If(SucrosPlant(ip)%akctype==1) Then
          x=PlantGeometry%dvs
       Else
          x=PlantData%dayOfYear
       End If
       PlantGeometry%akc=InterpolateTab(SucrosPlant(ip)%tab(11),x)
    Endif
    !PlantGeometry%akc=one_dp
    !end modification
    If(init==0) Then
       PlantGeometry%rootdepth=null_dp
       PlantGeometry%slaig=null_dp
       PlantGeometry%slaid=null_dp
       PlantGeometry%wlvg=null_dp
    Else
       PlantGeometry%rootdepth=SucrosPlant(ip)%root_init
       PlantGeometry%slaig=SucrosPlant(ip)%nsl*SucrosPlant(ip)%ssl
       PlantGeometry%slaid=Max(SucrosPlant(ip)%slaid_off-PlantGeometry%slaig,null_dp)
       PlantGeometry%wlvg=PlantGeometry%slaig/SucrosPlant(ip)%sla*PlantData%Area_ha_to_L2
    Endif
    PlantGeometry%pond_plant=null_dp
    PlantGeometry%pond_plant_new=null_dp
    PlantGeometry%eva_plants=null_dp
    PlantGeometry%transpiration = null_dp
    PlantGeometryArray%rex=null_dp
    ! reset variables not used outside the season of the plant
    ! (for cosmetic reasons)
    PlantGeometry%temp_sum=null_dp
    PlantGeometry%temp_sum_till=null_dp
    !PlantGeometry%temp_eff_sen=null_dp
    PlantGeometry%dvs=null_dp
    PlantGeometry%cgphot=null_dp
    PlantGeometry%wrtd=null_dp
    PlantGeometry%wso=null_dp
    PlantGeometry%wst=null_dp
    PlantGeometry%wcrn=null_dp
    PlantGeometry%wlvd=null_dp
    PlantGeometry%vd=null_dp
    !modification D.Farber
    ! respiration
    SucrosPlant(ip)%rlv = 0.
    SucrosPlant(ip)%rst = 0.
    SucrosPlant(ip)%rso = 0.
    SucrosPlant(ip)%rcrn = 0.
    SucrosPlant(ip)%rrt = 0.
    SucrosPlant(ip)%rgrowthrt = 0.
    SucrosPlant(ip)%rgrowthlv = 0.
    SucrosPlant(ip)%rgrowthst = 0.
    SucrosPlant(ip)%rgrowthso = 0.
    SucrosPlant(ip)%rgrowthcrn = 0.
    SucrosPlant(ip)%storage=0
    rnodert = 0.
    rnodexu = 0.
    GPP = 0
    NPP = 0
    aboveground_respiration = 0
    If(mode==1) Then
       PlantGeometry%wrt=0.05_dp*PlantGeometry%wrt
    Else
       PlantGeometry%wrt=0
    Endif
  End Subroutine ResetPlantData

  Subroutine ResetGrass()
    Implicit None
    PlantGeometry%dvs=null_dp
    PlantGeometry%temp_sum=null_dp
    PlantGeometry%temp_sum_till=null_dp
    PlantGeometry%cgphot=null_dp
  End Subroutine ResetGrass

  Subroutine InitPlant(ip)
    ! called when season starts
    Implicit None
    Integer, Intent(in) :: ip
    Integer :: ptype
    Real(dp) :: rna=null_dp
    ptype=Plant(ip)%type
    Call WriteOutput('Starting season of plant: '//PlantName(ptype))
    PlantGeometry%index=ip
    !raus  If(Plant(ip)%model==1) Then
    PlantGeometry%temp_sum=null_dp
    PlantGeometry%temp_sum_till=null_dp
    PlantGeometry%dvs=null_dp
    PlantGeometry%cgphot=null_dp
    PlantGeometry%wrt=null_dp
    PlantGeometry%wso=null_dp
    PlantGeometry%wst=null_dp
    PlantGeometry%wcrn=null_dp
    PlantGeometry%slaig=SucrosPlant(ip)%nsl*SucrosPlant(ip)%ssl
    PlantGeometry%slaid=Max(SucrosPlant(ip)%slaid_off-PlantGeometry%slaig,null_dp)
    PlantGeometry%wlvg=PlantGeometry%slaig/SucrosPlant(ip)%sla*PlantData%Area_ha_to_L2
    PlantGeometry%wlvd=null_dp
    PlantGeometry%pond_plant=null_dp
    PlantGeometry%rootdepth=SucrosPlant(ip)%root_init
    PlantGeometry%vd=null_dp
    PlantGeometry%count=null_dp   !(added by A.Klosterhalfen)
    If(ptype==SugarBeet) Then
       Plant(ip)%so_length=-350.0 ! mm
    Else If(ptype==Potatoe) Then
       Plant(ip)%so_length=-200.0 ! mm
    Else
       Plant(ip)%so_length=0.0
    End If
    If(ptype==SugarBeet .Or. ptype==Potatoe) Then
       Call DistributeRoots(Plant(ip)%root_distribution, PlantGeometryArray%rsod, rna,&
            Plant(ip)%so_length, PlantGeometryArray%z)
!       Print *,'roots=',PlantGeometryArray%rsod
!       Print *,'root_dist=',Plant(ip)%root_distribution%data
!       Print *,'z=',PlantGeometryArray%z
!       Print *,'length=',Plant(ip)%so_length
    End If
    !raus End If
  End Subroutine InitPlant


  Subroutine CalculatePlantGrowth(SimTime)
    Use Variables, Only: dailyCalculation
    Use TimeData, Only: EpsTime
    Implicit None
    Real(dp), Intent(in) :: SimTime
    Real(dp) :: t

    Call SetPlantTime(SimTime)
    if(dailyCalculation) then
       t = 0   ! Daily Calculation
    Else
       !calculate daytime (from 0 to 23)
       t = Nint(((SimTime-1)*TimeFactor-Int((SimTime-1)*TimeFactor+epsTime))/TimeFactor)
    end if
    Call CalculateLAI(t)

  End Subroutine CalculatePlantGrowth


  Subroutine CalculateLAI(t)
    ! calculate the leaf area index (LAI)
    Use Variables, Only: aboveground_respiration, GPP,NPP, waterstress, dailyCalculation, rootexu, &
         rootdeath, dailyCalculation, alphaAvg, RelHum, fluorescence_755nm, Nit_ancrt, Pho_ancrt
    Use Geometry, Only: rnodert, rnodexu, NumNP, rnodedeadw, rnodedeadwN, rnodedeadwP
    Use Carbon, Only: Patm, CO2Top
    Use Timedata, Only: SimTime
    Implicit None
    Real(dp), Intent(in) :: t
    Real(dp) :: temp, temp_avg, temp_avg2, temp_eff, temp_min, temp_max
    Real(dp) :: dec, rad         ! solar declination, rad=PI/180
    Real(dp) :: sinld, cosld     ! intermediate variables in calculating solar declination
    Real(dp) :: dl               ! daylength
    Real(dp) :: detr             ! daily extra terrestrial radiation (W/m2)
    Real(dp) :: dpar             ! daily photosynthetically active radiation
    Real(dp) :: dsinb, dsinbe, sinb=null_dp
    Real(dp) :: atmtr            ! atmospheric transmission coefficient
    Real(dp) :: frdf             ! diffuse radiation as a fraction of total solar radiation
    Real(dp) :: amdvs, amtmp, eff, dtga
    Real(dp) :: apar, pardf, pardr, parlsh, parlpp
    Real(dp) :: fslla, asssh, asssl, amax
    Real(dp) :: fgros, rlaic, refl, scp, sqrt_scp, rkdf, clustf, rkbl
    Real(dp) :: ratdvs=null_dp, wlv, rmaint, gphot
    Real(dp) :: dvs, plai, x
    Real(dp) :: fsh, frt=null_dp, flv, fso=null_dp, fst, fcrn=null_dp, asrq
    Real(dp) :: ratwtot, ratwrt, ratwrtd, ratwsh, ratwlvg, ratwst, ratwso, ratwcrn, ratgrtot, avCcont, ratexu
    Real(dp) :: rdrdv=null_dp, rdrsh, rdrwlvd
    Real(dp) :: ratlaid, ratwlvd, ratlaig, sla, slaig
    Real(dp) :: rlnew
    Real(dp) :: fT=null_dp, fP, fV
    Real(dp) :: alpha, fvn=null_dp, omega, rmax=null_dp, Tmin=null_dp, Topt=null_dp, Tmax=null_dp,VD
    Real(dp) :: rmaintlv, rmaintst, rmaintso, rmaintrt, rmaintcrn
    Integer :: i, j,k, ri, ip, ptype
    Real(dp) :: sunrise, sunset, time, hour
    Real(dp), Parameter :: secPerH=3600.0_dp
    Real(dp), Save :: sinbold=0.0_dp
    Real(dp) :: deathfac, deathfacMax
    Logical :: gaussIntegration   ! for daily calculation an additional gaussintegration is needed
    Real(dp), Dimension(3), Parameter :: &
         gsdst=(/ 0.112702_dp, 0.5_dp , 0.887298_dp /), &      ! distance in Gaussian integration
         gswt =(/ 0.277778_dp, 0.444444_dp, 0.277778_dp /)     ! weighting factor in Gaussian integration
    Real(dp) :: storagerat, storageFac, storageDecrease, glucosestor   !added by D.Farber
    Real(dp) :: Patm_pa, fact, ci_gamma=0
    Real(dp) :: photosynthesis, photosynthesis_sh, photosynthesis_sl
    Type(FarquharType), Pointer :: fq=>null()
    ! SIF
    Real(dp) :: A,Jo,Je,Phi_p,Phi_po,Phi_f,kf,kd,kn,fluorescence,apar_pho

    apar=0
    fluorescence_755nm=0
    atmtr=0
    !modification D.Farber
    ! average temperature and average daytime temperature for the current day
    temp_avg = ClimateData%temp_avg
    if(dailyCalculation) then
       temp_min = ClimateData%temp_min
       temp_max = ClimateData%temp_max
       temp_avg = (temp_min+temp_max)/2.0_dp
       temp = temp_avg
       gaussIntegration = .true.
    else
       temp = ClimateData%temp
       temp_min = ClimateData%temp_min
       temp_max = ClimateData%temp_max
    end if

    temp_avg2=0.25_dp*temp_min+0.75_dp*temp_max
    If(temp_min<null_dp) temp_min=2.0_dp+0.805_dp*temp_min   ! assessment of soil temperature for subsequent calculations
    If(temp_max<null_dp) temp_max=2.0_dp+0.805_dp*temp_max

    ! calculation of the daylength and daily extral terrestrial
    ! radiation from daynumber and latitude
    dec = -23.4_dp*cos(2.0_dp*PI*(PlantData%dayOfYear+10.0_dp)/365.0_dp)
    rad=PI/180.0_dp
    sinld= Sin(dec*rad)*Sin(PlantData%latitude*rad)
    cosld= Cos(dec*rad)*Cos(PlantData%latitude*rad)
    ! daylength (dl=h)
    sunrise = 12-Acos(-tan(PlantData%latitude*rad)*tan(dec*rad))/15/rad
    sunset = 12+Acos(-tan(PlantData%latitude*rad)*tan(dec*rad))/15/rad
    dl = sunset-sunrise
    if(.not. between(t, sunrise, sunset+1)) then
       sinbOld = null_dp
    end if
    
    ip=PlantGeometry%index
    If(ip==0) Return

    if(dailyCalculation) then
       ! daily integral of sine of solar inclination above the horizon (dsinb)
       dsinb=secPerH*(dl*sinld+24.0_dp*cosld*Sqrt(one_dp-(sinld/cosld)**2)/PI)
       ! daily integral of sinb  with a correction of lower atmospheric transmission
       ! at lower solar elevations (dsinbe)
       dsinbe=secPerH*(dl*(sinld+0.4_dp*(sinld*sinld+0.5_dp*cosld*cosld))+ &
            12.0_dp*cosld*(2.0_dp+3.0_dp*0.4_dp*sinld)*Sqrt(one_dp-(sinld/cosld)**2)/PI)
    else
       time = t
       If(t>sunset .And. t<sunset+one_dp) Then
          time = sunset
       End If
!       averaging deactivated
!       sinb = 0.5_dp*(sinbOld+(sin(dec*rad)*sin(PlantData%latitude*rad)+ &
!            Cos(dec*rad)*Cos(15._dp*(time-12._dp)*rad)*Cos(PlantData%latitude*rad)))
       sinb = sin(dec*rad)*sin(PlantData%latitude*rad)+ &
            Cos(dec*rad)*Cos(15._dp*(time-12._dp)*rad)*Cos(PlantData%latitude*rad)
       sinb=Max(sinb, 0.0_dp)
       sinbOld = sinb
       dsinb = sinb*secPerH
       ! ("timestep") integral of sinb  with a correction of lower  atmospheric transmittance
       ! at lower solar elevations
       dsinbe = sin(asin(sinb)+0.4_dp*sinb)*secPerH
    end if

    ! daily extra terrestrial radiation (detr) from corrected solar constant sc (1370)
    detr=1370.0_dp*(one_dp+0.033_dp*Cos(2.0_dp*PI*PlantData%dayOfYear/365.0_dp)) ! (W/m2)
    ! calculation of the daily radiation above the canopy
    ! calculation of the daily photosynthetically active radiation (dpar) and
    ! diffusive fraction of incoming radiation (frdf) from measured global
    ! radiation (global_rad) and atmospheric transmittance
    ! daily photosynthetically active radiation (J/m2/d)
    dpar=0.5_dp*ClimateData%global_rad*PlantData%Area_L2_to_m2  ! W/L2 -> W/m2
    ! diffuse radiation as a fraction of total solar radiation (frdf) from atmospheric transmittance
    If(dsinb>null_dp) Then
       atmtr=ClimateData%global_rad*PlantData%Area_L2_to_m2/(detr*dsinb) ! dimensionless
       If(atmtr.Le.0.07_dp) Then
          frdf=one_dp
       Else If(atmtr.Le.0.35_dp) Then
          frdf=one_dp-2.3_dp*(atmtr-0.07_dp)**2
       Else If(atmtr.Le.0.75_dp) Then
          frdf=1.333_dp -1.46_dp*atmtr
       Else  
          frdf=0.23_dp
       Endif
    Else
       frdf=null_dp
    Endif
    !raus If(Plant(ip)%model==2) Then
    !raus    PlantGeometry%akc=SimplePlant(ip)%akc
    !raus    PlantGeometry%slaig=SimplePlant(ip)%slaig
    !raus    PlantGeometry%slaid=SimplePlant(ip)%slaid
    !raus    PlantGeometry%rootdepth=SimplePlant(ip)%rootdepth
    !raus    Return
    !raus End If
    If(Plant(ip)%season==0) Return
    ptype=Plant(ip)%type
    dvs=PlantGeometry%dvs
    If(dvs>2.0_dp) Then
       If(ptype==WinterWheat .or. ptype==SummerWheat &
            .or. ptype==Maize .or. ptype==Gras .or. ptype==c4Gras .or. ptype==Barley) Return
    End If
    slaig=PlantGeometry%slaig
    plai=slaig+0.5_dp*SucrosPlant(ip)%eai
    PlantGeometry%plai=plai
    ! daily effective temperature for growth
    temp_eff=Max(temp-SucrosPlant(ip)%tempbase, null_dp)

    if( t .eq. null_dp) then
       PlantGeometry%temp_sum=PlantGeometry%temp_sum+Max(temp_avg-SucrosPlant(ip)%tempbase, null_dp)
    end if
    If(PlantGeometry%temp_sum<SucrosPlant(ip)%tempstart) Return
    ! calculation of actual CO2-assimilation rate of a unit leaf area
    ! for light saturation (amax) as a function of dvs/senescence (amdvs)
    ! and temperature (amtmp)
    If(ptype==Potatoe) Then
       amdvs=one_dp
    Else if(ptype==SugarBeet) then
       amdvs=InterpolateTab(SucrosPlant(ip)%tab(1),PlantGeometry%temp_sum)
    Else
       amdvs=InterpolateTab(SucrosPlant(ip)%tab(1),dvs)
    End If
    if(dailyCalculation) then
       amtmp=InterpolateTab(SucrosPlant(ip)%tab(2),temp_avg2)
    else
       amtmp=InterpolateTab(SucrosPlant(ip)%tab(2),temp)
    end if
    amax=SucrosPlant(ip)%amx*amdvs*amtmp*PlantData%Area_ha_to_L2   ! kg CO2/ha leaf/h -> kg CO2/L2 leaf/h
    ! computes potential daily total gross assimilation of the canopy (dtga, kg CO2/ha/d)
    eff=SucrosPlant(ip)%eff*PlantData%Area_ha_to_L2 ! (kg CO2/ha leaf/h)/(J/m2/s)) -> (kg CO2/L2 leaf/h)/(J/m2/s))
    rkdf=SucrosPlant(ip)%rkdf
    scp=SucrosPlant(ip)%scp
    sqrt_scp=Sqrt(one_dp-scp)
    dtga=null_dp
    If(Farquhar_used) fq => Plant(ip)%far ! abbreviation
    fq => Plant(ip)%far ! testing
    If(sinb > 1E-10 .Or. dailyCalculation) Then

       Do i = 1,3   ! if no gaussintegration is needed (more than one calculation per day)
                    ! this loop will stop after the first interation
          if(dailyCalculation) then
             hour=12.0_dp+dl*0.5_dp*gsdst(i)
             ! calculation of the instantaneous radiation above the canopy on this hour
             ! sin of solar inclination above the horizon
             sinb=Max(null_dp,sinld+cosld*Cos(2.0_dp*PI*(hour+12.0_dp)/24.0_dp))
          end if
          ! diffuse PAR (pardf) and direct PAR (pardr) (W/m2)
          apar=dpar*sinb*(1.0_dp+0.4_dp*sinb)/dsinbe   !(D.Farber) Photosynthetically active radiation flux
          pardf=Min(apar,frdf*dpar*sinb/dsinb)
          pardr=apar-pardf

          ! calculation of hourly CO2 assimilation rate (frgos) at different canopy depths
          ! with partial cumulated LAI at various canopy depth (rlaic)
          fgros=null_dp
          ! selection of the canopy depths
          Do j=1,3
             rlaic=plai*gsdst(j)
             ! calculates the radiation profile in the canopy and gives instantaneous
             ! values of absorbed radiation for succesive leaf layers
             ! canopy reflection coefficient for PAR (refl)
             ! as a function of the leaf scattering coefficient (scp)
             refl=(one_dp-sqrt_scp)/(one_dp+sqrt_scp)         
             ! extinction coefficient for direct component of PAR (rkbl) and total direct flux 
             ! (akdrt) / cluster factor and ratio between emperical and theoretical value of (rkdf)
             clustf=rkdf/(0.8_dp*sqrt_scp)
             rkbl=(0.5_dp/sinb)*clustf
             ! absorbed radiation fluxes per unit leaf area (W/m2 leaf)
             ! diffuse flux, total direct flux and direct component of direct flux
             ! absorbed fluxes (W/m2 leaf) for shaded (parlsh) and sunlit leaves per unit leaf area (parlsl)
             parlsh=(one_dp-refl)*pardf*rkdf*Exp(-rkdf*rlaic)+ &
                  (one_dp-refl)*pardr*rkbl*sqrt_scp*Exp(-rkbl*sqrt_scp*rlaic)- &
                  (one_dp-scp)*pardr*rkbl*Exp(-rkbl*rlaic)
             ! direct PAR absorbed by leaves perpendicular to direct beam (parlpp) (W/m2 leaf)
             parlpp=pardr*(one_dp-scp)/sinb
             ! fraction of sunlit leaf area (fslla)
             fslla=Exp(-rkbl*rlaic)*clustf
             If(Farquhar_used) Then
                ! Patm from kg/L/T**2 to kg/m/s**2
                Patm_pa=Patm*PlantData%Unit_L_to_m/(86400_dp*TimeFactor)**2
                fq%c_s=CO2Top*Patm_pa
                fq%o_i=0.209*Patm_pa ! O2 partial pressure (Pascal)
                fq%K_c = fq%K_c25 * fq%a_kc ** ((temp - 25_dp)/10_dp) ! Pa
                fq%K_o = fq%K_o25 * fq%a_ko ** ((temp - 25_dp)/10_dp) ! Pa
                fq%gamma=0.5 * fq%K_c/fq%K_o * 0.21*fq%o_i
                fq%ftv=one_dp/(one_dp+Exp((-220000_dp+710_dp*(temp+fq%temp_fw))&
                     /(0.001*fq%R_gas*(temp+fq%temp_fw))))
                fq%vcmax=fq%vcmax25N * 2.4_dp**((temp - 25_dp)/10_dp) * fq%ftv &
                     * (dl/fq%daylength_max)**2 * alphaAvg
                If(Plant(ip)%c3) Then
                   fq%c_i = Max(fq%gamma,fq%c_s*(1.0-1.6/(fq%m*RelHum*Max(alphaAvg,0.0001_dp))))
                Else
                   fq%c_i = Max(9.9e-1*Patm_pa,fq%c_s*(1-1.6/(fq%m*RelHum*Max(alphaAvg,0.0001_dp))))
                End If
                ! export limited rate of carboxylation for C 3 plant asnd the PEP carboxylase limited
                ! rate of carboxylation for C 4 plants w_e ( µ mol CO2 / m**2 / s )
                If(Plant(ip)%c3) Then
                   fq%w_e=0.5*fq%vcmax
                Else
                   fq%w_e=4000_dp*fq%vcmax*fq%c_i/Patm_pa ! or 18000 or 20000 instead of 4000
                End If
                photosynthesis=fq%w_e
                ci_gamma=fq%c_i - fq%gamma
                !The RuBP carboxylase (Rubisco) limited rate of carboxylation:
                If(Plant(ip)%c3) Then
                   If(ci_gamma>null_dp) Then
                      fq%w_c=fq%vcmax * ci_gamma / (  fq%c_i + fq%K_c*(1_dp+fq%o_i/fq%K_o) )
                      photosynthesis=Min(photosynthesis, fq%w_c)
                   Else
                      fq%w_c=null_dp
                   End If
                Else
                   fq%w_c=fq%vcmax
                   photosynthesis=Min(photosynthesis, fq%w_c)
                End If
                ! The maximum rate of carboxylation allowed by the capacity to regenerate RuBP
                ! (i.e., the light-limited rate) w_j ( µ mol CO2 / m**2 / s )
                !  shaded leaves
                If(Plant(ip)%c3) Then
                   If(ci_gamma>null_dp) Then
                      fq%w_j_sh=4.6*parlsh*fq%alfa * ci_gamma / (fq%c_i+2.0*fq%gamma)
                      photosynthesis_sh=Min(photosynthesis, fq%w_j_sh)
                   Else
                      fq%w_j_sh=null_dp
                      photosynthesis_sh=photosynthesis
                   End If
                Else
                   fq%w_j_sh=4.6*parlsh*fq%alfa
                   photosynthesis_sh=Min(photosynthesis, fq%w_j_sh)
                End If
                fq%w_j_sl=null_dp  !  sunlit leaves 
                Do k=1,3
                   If(Plant(ip)%c3) Then
                      If(ci_gamma>null_dp) Then
                         fq%w_j_sl=fq%w_j_sl+4.6*(parlsh+parlpp*gsdst(k))*fq%alfa &
                              * ci_gamma / (fq%c_i+2.0*fq%gamma) * gswt(k)
                      End If
                   Else
                      fq%w_j_sl=fq%w_j_sl+4.6*(parlsh+parlpp*gsdst(k))*fq%alfa*gswt(k)
                   End If
                End Do
                If(Plant(ip)%c3) Then
                   If(ci_gamma>null_dp) Then
                      photosynthesis_sl=Min(photosynthesis, fq%w_j_sl)
                   Else
                      photosynthesis_sl=photosynthesis
                   End If
                Else
                   photosynthesis_sl=Min(photosynthesis, fq%w_j_sl)
                End If
                ! µmol CO2 / m**2 / s = 0.0001584 kg CO2 / cm**2 / h
                ! 1 mol co2 = 44 g
                ! 1 umol co2 = 44e-6 g
                ! 1 umol co2 = 44e-9 kg
                ! m -> L : Units(4) = m   Unit(4)/unit(UnitinTrace)
                fact=44e-9*(Units_L(Unit_L_Input)/Units_l(4))**2 * secPerH
                asssh=photosynthesis_sh*fact
                asssl=photosynthesis_sl*fact
             Else ! end Farquhar
                ! CO2 assimilation rate of shaded leaf area (asssh, kg CO2/L2 leaf/h)
                asssh=amax*(one_dp-Exp(-eff*parlsh/amax))
                ! CO2 assimilation rate of sunlit leaf area (asssl, kg CO2/L2 leaf/h)
                asssl=null_dp 
                Do k=1,3
                   asssl=asssl+amax*(one_dp-Exp(-eff*(parlsh+parlpp*gsdst(k))/amax))*gswt(k)
                End Do
                If(debug) Then
                   Write(34,*) simtime%time,parlsh,parlpp,plai,clustf,rkbl,sinb,asssh,asssl,eff,amax,pardr,pardf,apar
                End If
             End If
             ! hourly CO2 assimilation rate of the crop (kg CO2/L2 soil/h) 
             fgros=fgros+((one_dp-fslla)*asssh+fslla*asssl)*plai*gswt(j)
          End Do
          If(dailyCalculation) Then
             ! integration of instantaneous assimilation to a daily total (dtga)
             dtga=dtga+fgros*dl*gswt(i)
          else
             ! "timestep" instantaneous assimilation
             dtga=fgros   ! kg CO2/L2/h
             Exit
          End If
       End Do

       ! solar induced chlorophyll fluorescence (SIF)
       fluorescence_755nm=null_dp
       If(Farquhar_used .And. apar>1e-8) Then
!          A=fslla*photosynthesis_sl+(one_dp-fslla)*photosynthesis_sh ! A=photosysnthesis (umol co2/m2/s)
          If(plai>null_dp) Then
             A=dtga/plai/fact ! A=photosysnthesis (umol co2/m2/s)
          Else
             A=null_dp
          Endif
          apar_pho=4.6*apar ! J/m2/s -> umol photons/m2/s
          Jo=apar_pho*fq%alfa ! maximum possible electron transport rate Jo (umol co2/m2/s)
          ! actual electron transport rate Je
          If(Plant(ip)%c3) Then
             If(ci_gamma>1e-8) Then
                Je=A*(fq%c_i+2.0*fq%gamma)/(fq%c_i-fq%gamma) ! (umol co2/m2/s)
             Else
                Je=null_dp
             End If
          Else
             Je=A ! (umol co2/m2/s)
          End If
          Phi_Po = 0.8 ! (according to Björkman and Demmig, 1987) is the efficiency of photochemical trapping in the dark adapted state (dimless)
          If(jo<1e-10) Then
             Phi_p=null_dp
          Else
             Phi_p=Phi_po*Min(Je/Jo,one_dp) ! photochemical quantum yield (dimless)
          End If
          kf=0.05 ! rate coefficient fluorescence (1/s)
          kd=Max(0.03*temp+0.0773, 0.087) ! rate coefficient dark-adapted (1/s)
          x=1.0-Phi_p/Phi_Po
          kn=(6.2473*x-0.5944)*x ! rate coefficient light-adapted (1/s)
          Phi_f=kf/(kf+kd+kn)*(1.0-Phi_p) ! fluorescence yield (dimless)
          ! fluorescence=APAR*Phi_f
          fluorescence=apar_pho*Phi_f ! leaf level SIF (umol photons/m2/s)
          If(fq%vcmax25N<=70.0) Then
             ! k= conversion factor from leaf level fluorescence to spectrometer oserved fluorescence, see Fig. 1 Lee et al.
             x=0.047716*fq%vcmax25N+7.70092
          Else
             x=0.032686*fq%vcmax25N+8.75302
          Endif
          fluorescence_755nm=fluorescence/x ! SIF at 755nm (W/m2/sr/um)
          If(debug) Then
             Write(35,*) simtime%time,A,apar_pho,Jo,Je,Phi_p,kf,kd,kn,fluorescence,fluorescence_755nm, fact
          End If
       End If
    End If
    If(debug) Then
       Write(33,*) simtime%time,fq%w_c,fq%w_j_sh,fq%w_j_sl,fq%w_e,&
            fq%K_c,fq%K_o,fq%c_i,fq%gamma,parlsh,parlpp,plai,clustf,rkbl,sinb,asssh,asssl,fslla,&
            temp,patm_pa,fq%o_i,Units_L(4)/Units_L(Unit_L_Input),fact,&
            alphaAvg,frdf,dsinb,apar,pardf,pardr,atmtr,detr,ClimateData%global_rad,detr*dsinb,sinb,dsinbe
    End If
    If(Farquhar_used) Then
       Write(64,12) Nint(simtime%time),fq%w_c,fq%w_j_sh,fq%w_j_sl,fq%w_e,&
            fq%K_c,fq%K_o,fq%c_i,fq%gamma,parlpp,parlpp,plai,clustf,rkbl,sinb,asssh,asssl,fslla,&
            temp,patm_pa,fq%o_i,Units_L(4)/Units_L(Unit_L_Input),fact,&
            alphaAvg,frdf,dsinb,apar,pardf,pardr,atmtr,detr,ClimateData%global_rad,detr*dsinb,sinb,dsinbe 
       If(sinb > 1E-10 .And. apar>1e-8) Then
          Write(63,1) Nint(simtime%time),A,apar_pho,Jo,Je,Phi_p,kf,kd,kn,fluorescence,fluorescence_755nm, fact
       Else
          Write(63,1) Nint(simtime%time),null_dp,null_dp,null_dp,null_dp,null_dp,null_dp,null_dp,null_dp,null_dp,null_dp, null_dp
       End If
1      Format(I8,1P,1X,32E21.6)
12     Format(I8,1P,1X,34E21.6)
    End If
    ! conversion from CO2 assimilation to glucose production
    ! and reduction for water stress
    ! daily total gross assimilation (gphot, kg CH2O/ha/d)

    If (waterstress > 0 .AND. Farquhar_used .EQV. .FALSE.) Then
       gphot=dtga*(30.0_dp/44.0_dp)*Min(alphaAvg,NitReduct,PhoReduct)   ! kg CH2O/L2/"timestep"
    else
       gphot=dtga*(30.0_dp/44.0_dp)*Min(NitReduct,PhoReduct)
    Endif

    if( t .eq. null_dp) then
       ! calculation of the crop development rate
       Select Case(ptype)
       Case(SugarBeet)
          ratdvs=Min(19.0_dp,Max(temp_avg-2.0_dp,null_dp))
       Case(Potatoe)
          If(temp_avg<13.0_dp) Then
             x=temp_avg-2.0_dp
          Else
             x=29.0_dp-temp_avg
          Endif
          ratdvs=Min(11.0_dp,Max(x,null_dp))
       Case(WinterWheat,SummerWheat,Maize,Gras,c4Gras,Barley)
          !if-else added by D.Farber
          if( SucrosPlant(ip)%lCeres == 1 ) then
             !modification N.Prolingheuer
             ! Calculating the development stage according to CERES-Wheat (Wang & Engel, 1998, Streck et al., 2003)
             ! Depending on developmental stages (EM-TS-AN-PM) and temperature (fT), photoperiod (fP) and vernalization (fV) 
             ! EM-TS: fT, fP, FV
             ! TS-AN: fT, fP
             ! AN-PM: fT
             ! Temperature reduction function (fT)
             if(dvs<0.4) Then                         ! EM-TS
                Tmin=SucrosPlant(ip)%Tmin_v1
                Topt=SucrosPlant(ip)%Topt_v1
                Tmax=SucrosPlant(ip)%Tmax_v1
             else if(dvs.ge.0.4 .and. dvs<1.0) Then   ! TS-AN
                Tmin=SucrosPlant(ip)%Tmin_v2
                Topt=SucrosPlant(ip)%Topt_v2
                Tmax=SucrosPlant(ip)%Tmax_v2
             else if(dvs.ge.1.0) Then                 ! AN-PM
                Tmin=SucrosPlant(ip)%Tmin_r
                Topt=SucrosPlant(ip)%Topt_r
                Tmax=SucrosPlant(ip)%Tmax_r
             end if

             if(Tmin.le.temp_avg .and. temp_avg.le.Tmax) Then
                alpha=log(2.0)/(log((Tmax-Tmin)/(Topt-Tmin)))
                fT=(2.0*(temp_avg-Tmin)**alpha*(Topt-Tmin)**alpha-(temp_avg-Tmin)**(2.0*alpha)) &
                     /((Topt-Tmin)**(2.0*alpha))
             else if(temp_avg<Tmin .or. temp_avg>Tmax) Then
                fT=null_dp
             end if

             ! Vernalization reduction function (fV) (Streck et al., 2003)

             ! effective vernalization days (VD)
             ! daily vernalization rate function (fvn)
             if(dvs<0.4) Then   ! EM-TS
                if(SucrosPlant(ip)%Tmin_vn.le.temp_avg .and. temp_avg.le.SucrosPlant(ip)%Tmax_vn) Then
                   alpha=log(2.0)/(log((SucrosPlant(ip)%Tmax_vn-SucrosPlant(ip)%Tmin_vn) &
                        /(SucrosPlant(ip)%Topt_vn-SucrosPlant(ip)%Tmin_vn)))
                   fvn=(2.0*(temp_avg-SucrosPlant(ip)%Tmin_vn)**alpha*(SucrosPlant(ip)%Topt_vn &
                        -SucrosPlant(ip)%Tmin_vn)**alpha-(temp_avg-SucrosPlant(ip)%Tmin_vn)**(2.0*alpha)) &
                        /((SucrosPlant(ip)%Topt_vn-SucrosPlant(ip)%Tmin_vn)**(2*alpha))
                else if(temp_avg<SucrosPlant(ip)%Tmin_vn .or. temp_avg>SucrosPlant(ip)%Tmax_vn) Then
                   fvn=null_dp
                end if
                VD=PlantGeometry%VD+fvn
                PlantGeometry%VD=VD
                fV=(VD**5.0)/(22.5**5.0+VD**5.0)
             else   ! TS-AN and AN-PM
                fV=one_dp
             endif

             ! Photoperiod reduction function (fP)
             if(dvs<1.0) Then   ! EM-TS and TS-AN
                if(SucrosPlant(ip)%Popt/=null_dp) then
                   omega=4.0/(SucrosPlant(ip)%Popt-SucrosPlant(ip)%Pcrit)
                else
                   omega=SucrosPlant(ip)%omega
                end if
                fP=1-exp(-omega*(dl-SucrosPlant(ip)%Pcrit))
             else   ! AN-PM
                fP=one_dp
             endif

             ! rate of development (ratdvs)
             if(dvs<0.4) Then   ! EM-TS
                rmax=SucrosPlant(ip)%rmax_v1
             else if(dvs.ge.0.4 .and. dvs<1.0) Then   ! TS-AN
                rmax=SucrosPlant(ip)%rmax_v2
             else if(dvs.ge.1.0) Then   ! AN-PM
                rmax=SucrosPlant(ip)%rmax_r
             end if

             ratdvs=rmax*fT*fP*fV
          else
             !end modification (original model:)
             ! Development rate (ratdvs)
             If(dvs<one_dp) Then
                k=3
             Else
                k=4
             Endif
             ratdvs=InterpolateTab(SucrosPlant(ip)%tab(k),temp_avg)
          end if
       case default
          write(*,*) "error in ptype"
       End Select
    else
       ratdvs = 0._dp
    end if

    ! maintenance glucose requirements
    ! maintenance respiration of crop (rmaint)
    wlv=PlantGeometry%wlvg+PlantGeometry%wlvd   ! dry weight of leaves per unit area
    If(wlv<=null_dp) Then 
       x=null_dp
    Else
       x=PlantGeometry%wlvg/wlv
    Endif
    rmaint=0.03_dp*TimeFactor*PlantGeometry%wlvg &   !modified by D.Farber (original: 0.03_dp*wlv)
         +0.015_dp*TimeFactor*PlantGeometry%wst &
         +0.015_dp*TimeFactor*PlantGeometry%wrt &
         +SucrosPlant(ip)%rmainso*TimeFactor*PlantGeometry%wso &
         +0.015_dp*TimeFactor*PlantGeometry%wcrn   ! kg CH2O/L2/"timestep"

    !modification D.Farber
    ! relative fractions of maintenance respiration for each plant organ
    rmaintlv = 0.03_dp*TimeFactor*PlantGeometry%wlvg/rmaint
    rmaintst = 0.015_dp*TimeFactor*PlantGeometry%wst/rmaint
    rmaintrt = 0.015_dp*TimeFactor*PlantGeometry%wrt/rmaint
    rmaintso = SucrosPlant(ip)%rmainso*TimeFactor*PlantGeometry%wso/rmaint
    rmaintcrn = 0.015_dp*TimeFactor*PlantGeometry%wcrn/rmaint
    !end modification

    rmaint = rmaint*x*2.0_dp**((temp-25.0_dp)/10.0_dp)
    ! Q10=2.0, base coefficient (accounting for increase in rmaint with a 10C rise in temperature)

    If (dailyCalculation) Then
       rmaint=Min(rmaint, gphot)
    end if

    ! calculation of the dry matter fraction for each plant organ
    ! calculate cumulative values for the gross photosynthesis
    PlantGeometry%cgphot=PlantGeometry%cgphot+gphot
    If(ptype==Potatoe) Then
       x=(dvs-one_dp/(0.0015_dp+0.00079_dp*SucrosPlant(ip)%rmatr))/430.0_dp
       fsh=Min(one_dp,Max(0.8_dp+0.2_dp*x, 0.8_dp))
       frt=one_dp-fsh
       flv=Min(0.75_dp,Max(0.75_dp-x,null_dp))
       fso=Min(1.0D0,Max(x,null_dp))
       fst=one_dp-flv-fso
       fcrn=null_dp
    ! Modified by M.Herbst 7.7.2016 for sugarbeet not related to DVS but to temp_sum
    Else if (ptype==SugarBeet) then
       fsh=InterpolateTab(SucrosPlant(ip)%tab(5),PlantGeometry%temp_sum)
       flv=InterpolateTab(SucrosPlant(ip)%tab(6),PlantGeometry%temp_sum)
       fst=InterpolateTab(SucrosPlant(ip)%tab(7),PlantGeometry%temp_sum)
       fcrn=one_dp-flv-fst
       frt=InterpolateTab(SucrosPlant(ip)%tab(8),PlantGeometry%temp_sum)*(one_dp-fsh)
       fso=one_dp-fsh-frt
    Else    ! no potatoes or SugarBeet
       fsh=InterpolateTab(SucrosPlant(ip)%tab(5),dvs)
       flv=InterpolateTab(SucrosPlant(ip)%tab(6),dvs)
       fst=InterpolateTab(SucrosPlant(ip)%tab(7),dvs)
       Select Case(ptype)
       Case(WinterWheat,SummerWheat,Barley,Gras,c4Gras)
          fso=one_dp-flv-fst
          frt=one_dp-fsh
          fcrn=null_dp
       Case(Maize)
          fcrn=InterpolateTab(SucrosPlant(ip)%tab(8),dvs)
          fso=one_dp-flv-fst-fcrn
          frt=one_dp-fsh
       End Select
    End If

    ! global glucose production to dry matter production
    ! conversion efficiency coefficient (asrq) (assimilate requirement for dry matter production)
    asrq=fsh*(1.46*flv+1.51*fst+1.51*fcrn)+1.44*frt

    !modification D.Farber
    ! according to SUCROS97: Simulation of crop growth for potential and water-limited production situations
    ! as applied to spring wheat, Editors: H.H. van Laar, J. Goudriaan, H. van Keulen
    avCcont = fsh*(0.459*flv+0.494*fst+0.491*fcrn)+0.467*frt   ! organ-specific carbon content: g C/g DM * organ fraction
    If(ptype==SugarBeet) Then
       asrq=asrq+fso*SucrosPlant(ip)%asrqso
       avCcont = avCcont+0.446*fso
    Else
       asrq=asrq+fso*fsh*SucrosPlant(ip)%asrqso
       avCcont = avCcont+0.471*fso*fsh
    End If
    ! net dry matter growth rate of the crop per unit area (rawtot)
    ! and for each plant organ per unit area


    glucosestor = gphot-rmaint

    If(ptype == Gras .Or. ptype == c4Gras) Then   ! LINGRA
       storageFac = 0.
       if(dvs >  SucrosPlant(ip)%dvsLow) then
          storageFac = &
               SucrosPlant(ip)%storageFacMax/(SucrosPlant(ip)%dvsHigh-SucrosPlant(ip)%dvsLow)*(dvs-SucrosPlant(ip)%dvsLow)
       else if(dvs > SucrosPlant(ip)%dvsHigh) then
          storageFac = SucrosPlant(ip)%storageFacMax
       end if
       storagerat = glucosestor*storageFac
       if(SucrosPlant(ip)%storage+storagerat > SucrosPlant(ip)%storageMax) then
          storagerat = SucrosPlant(ip)%storageMax - SucrosPlant(ip)%storage
          SucrosPlant(ip)%storage = SucrosPlant(ip)%storageMax
       else
          SucrosPlant(ip)%storage = SucrosPlant(ip)%storage + storagerat
       end if
       glucosestor = glucosestor - storagerat
    end if

    ! dry matter growth rates
    ratwtot=glucosestor/asrq

    ratgrtot = (gphot-rmaint)-ratwtot*avCcont*30.0_dp/12.0_dp   ! kg CH2O/L2/timestep -> total growth respiration
    ratgrtot = max(ratgrtot, null_dp)   !set ratgrtot to 0 during night if ratwtot is negative due to maintenance respiration

    If(rootexu) Then
       ratwrt=ratwtot*frt*(1._dp-SucrosPlant(ip)%exu_fac)
       ratexu=ratwtot*frt*SucrosPlant(ip)%exu_fac *0.467   ! kg C/L2/"timestep"
       Do ri=1,NumNP
          rnodexu(ri) = PlantGeometryArray%rrd(NumNP-ri+1)/PlantGeometryArray%dz(NumNP-ri+1)*PlantData%UnitFactor*ratexu   ! kg C/L3/"timestep"
       End Do
    Else
       ratwrt=ratwtot*frt
    Endif

    If(rootdeath) Then
       deathfacMax = SucrosPlant(ip)%deathfacMax
       !       if(ptype==WinterWheat .or. ptype==SummerWheat .or. ptype==Barley .or. ptype==Gras .or. ptype==c4Gras) then
       if(dvs < 0.2) then
          deathfac = 0.   ! no rootdeath
       else if(dvs > 0.5) then
          deathfac = deathfacMax
       else
          deathfac = deathfacMax/(0.5-0.2)*(dvs-0.2)   ! lower rootdeath during juvenile growth
       end if
       !       else
!          deathfac = 0.   ! needs to be modified for other crops
!       end if
       ratwrtd = ratwrt*deathfac
       ratwrt = ratwrt*(1-deathfac)
       Do ri=1,NumNP
          rnodedeadw(ri) = PlantGeometryArray%rrd(NumNP-ri+1)/PlantGeometryArray%dz(NumNP-ri+1)&
               *PlantData%UnitFactor*ratwrtd*0.467   ! kg C/L3/timestep
          rnodedeadwN(ri) = PlantGeometryArray%rrd(NumNP-ri+1)/PlantGeometryArray%dz(NumNP-ri+1)&
               *PlantData%UnitFactor*ratwrtd*Nit_ancrt   ! kg C/L3/timestep
          rnodedeadwP(ri) = PlantGeometryArray%rrd(NumNP-ri+1)/PlantGeometryArray%dz(NumNP-ri+1)&
               *PlantData%UnitFactor*ratwrtd*Pho_ancrt   ! kg C/L3/timestep
       End Do
    Else
       ratwrtd = 0.
    End If
    PlantGeometry%ratwrtd=ratwrtd ! for use in nit_upt()

    ratwsh=ratwtot*fsh
    if((ptype == Gras .or. ptype==c4Gras) .and. dvs < SucrosPlant(ip)%dvsLow) then   ! LINGRA
       storageDecrease = log(100.)/SucrosPlant(ip)%Tstorage*TimeFactor
       ratwsh = ratwsh + SucrosPlant(ip)%storage*storageDecrease/(fsh*(1.46*flv+1.51*fst+1.51*fcrn))
       SucrosPlant(ip)%storage = SucrosPlant(ip)%storage*(1-storageDecrease)
    end if
    ratwlvg=ratwsh*flv
    ratwst=ratwsh*fst
    ratwso=ratwsh*fso
    ratwcrn=ratwsh*fcrn

    If(ptype==SugarBeet) ratwso=ratwtot*fso
    ! growth respiration of leaves, stems, root and storage organs
    SucrosPlant(ip)%rgrowthlv = ratgrtot*fsh*flv/(30.0_dp/44.0_dp)/molCO2   ! mol CO2/L2/"timestep"
    SucrosPlant(ip)%rgrowthst = ratgrtot*fsh*fst/(30.0_dp/44.0_dp)/molCO2
    SucrosPlant(ip)%rgrowthso = ratgrtot*fsh*fso/(30.0_dp/44.0_dp)/molCO2
    SucrosPlant(ip)%rgrowthrt = ratgrtot*frt/(30.0_dp/44.0_dp)/molCO2
    SucrosPlant(ip)%rgrowthcrn = ratgrtot*fsh*fcrn/(30.0_dp/44.0_dp)/molCO2

    ! maintenance respiration + growth respiration
    SucrosPlant(ip)%rlv = rmaintlv*rmaint/(30.0_dp/44.0_dp)/molCO2 +SucrosPlant(ip)%rgrowthlv   ! mol CO2/L2/d
    SucrosPlant(ip)%rst = rmaintst*rmaint/(30.0_dp/44.0_dp)/molCO2 +SucrosPlant(ip)%rgrowthst
    SucrosPlant(ip)%rso = rmaintso*rmaint/(30.0_dp/44.0_dp)/molCO2 +SucrosPlant(ip)%rgrowthso
    SucrosPlant(ip)%rrt = rmaintrt*rmaint/(30.0_dp/44.0_dp)/molCO2 +SucrosPlant(ip)%rgrowthrt
    SucrosPlant(ip)%rcrn = rmaintcrn*rmaint/(30.0_dp/44.0_dp)/molCO2 +SucrosPlant(ip)%rgrowthcrn

    ! root respiration at a specific node
    Do ri=1,NumNP
       rnodert(ri) = PlantGeometryArray%rrd(NumNP-ri+1)/PlantGeometryArray%dz(NumNP-ri+1)*PlantData%UnitFactor*SucrosPlant(ip)%rrt   ! mol CO2/L3/d
    End Do
    If(ptype==SugarBeet .Or. ptype==potatoe) Then
       Do ri=1,NumNP
          rnodert(ri) = rnodert(ri)+PlantGeometryArray%rsod(NumNP-ri+1) / &
                        PlantGeometryArray%dz(NumNP-ri+1) * &
                        PlantData%UnitFactor*SucrosPlant(ip)%rso   ! mol CO2/L3/d
       End Do
    End If

    GPP= -1.0_dp*gphot/(30.0_dp/44.0_dp)/MolCO2   ! mol CO2/L2/timestep
    If(ptype==SugarBeet .Or. ptype==potatoe) Then
       aboveground_respiration = SucrosPlant(ip)%rlv+SucrosPlant(ip)%rst+SucrosPlant(ip)%rcrn   ! mol CO2/L2/timestep
       NPP=GPP+aboveground_respiration+SucrosPlant(ip)%rrt+SucrosPlant(ip)%rso   ! mol CO2/L2/timestep
    Else
       aboveground_respiration = SucrosPlant(ip)%rlv+SucrosPlant(ip)%rst+SucrosPlant(ip)%rso+SucrosPlant(ip)%rcrn   ! mol CO2/L2/timestep
       NPP=GPP+aboveground_respiration+SucrosPlant(ip)%rrt   ! mol CO2/L2/timestep
    End If
    !end modification

    ! calculation of the death rates for leaf weight and LAI
    ! senescence component of the growth rate coefficient of dead leaf area (rdrdv)
    Select Case(ptype)
    Case(WinterWheat,SummerWheat,Barley)
       If(dvs<one_dp) Then
          rdrdv=null_dp
       Else
          rdrdv=InterpolateTab(SucrosPlant(ip)%tab(10),temp_eff)
       End If
       If(ptype==WinterWheat .or. ptype==Barley) Then
            rdrdv=rdrdv*InterpolateTab(SucrosPlant(ip)%tab(9),dvs)
       End If
    Case(Gras,c4Gras)   ! LINGRA
       If(dvs<0.3) Then
          rdrdv=null_dp
       Else
          rdrdv=InterpolateTab(SucrosPlant(ip)%tab(10),temp_eff)
       End If
    Case(SugarBeet)
       rdrdv=ratdvs*InterpolateTab(SucrosPlant(ip)%tab(10),PlantGeometry%temp_sum)
    Case(Maize)
       If(dvs<1.35_dp) Then
          rdrdv=0.0005_dp
       Else
          rdrdv=0.003_dp
       Endif
       rdrdv=rdrdv*(temp_avg-8.0_dp)
    Case(Potatoe)
       x=PlantGeometry%temp_sum-senescence_start_temp
       If(x<null_dp) Then
          rdrdv=null_dp
       Else 
          rdrdv=Max(temp_avg-2.0_dp,8.0_dp) &
               *Exp(-11.7_dp+0.68_dp*SucrosPlant(ip)%rmatr) &
               *Exp(x*(0.0068_dp-0.0006_dp*SucrosPlant(ip)%rmatr))
       Endif
    End Select
    ! shade component of the growth rate coefficient of the dead leaf area (rdrsh)
    x=0.03_dp*(plai-SucrosPlant(ip)%rlaicr)/SucrosPlant(ip)%rlaicr
    rdrsh=Max(null_dp,Min(x,0.03_dp))
    ! relative death rates of leaves (rdrwlvd)
    rdrwlvd=Max(rdrsh,rdrdv)
    ! coefficient for chilling temperatures (maize only)
    If(ptype==Maize) Then
       If(dvs<one_dp) Then
          rdrwlvd=null_dp
       Else
          If(dvs<1.25_dp) Then
             x=null_dp
          Else
             x=Min(one_dp,Max(null_dp,(6.0_dp-temp_avg)/6.0_dp))
          Endif
          rdrwlvd=Max(0.001_dp,x,rdrwlvd)
       Endif
    else if (ptype==Gras .or. ptype==c4Gras) then   ! LINGRA
       If(dvs<0.3) Then
          rdrwlvd=null_dp
       Else
          If(dvs<2.0_dp) Then
             x=null_dp
          Else
             x=Min(one_dp,Max(null_dp,(6.0_dp-temp_avg)/6.0_dp))
          Endif
          rdrwlvd=Max(0.001_dp,x,rdrwlvd)
       Endif
    Endif
    ! dead leaf area growth rate 
    ratlaid=slaig*(Exp(rdrwlvd*TimeFactor)-one_dp)
    ! net dry matter growth rate for dead leaves per unit area (ratwlvd) 
    If(plai<=null_dp .Or. slaig<=null_dp  ) Then
       ratwlvd=null_dp
    Else
       ratwlvd=PlantGeometry%wlvg*ratlaid/slaig
    Endif
    ! specific leaf area for maize
    If(ptype==Maize .or. ptype==Gras .or. ptype==c4Gras) Then   !(ptype Gras added by A.Klosterhalfen)
       sla=InterpolateTab(SucrosPlant(ip)%tab(9),dvs)
    Else
       sla=SucrosPlant(ip)%sla
    End If
    ! green leaf area growth rate (ratlaig)
    Select Case(ptype)
    Case(Potatoe,SugarBeet)
       If(PlantGeometry%temp_sum<=450.0_dp .And. plai<=0.75_dp) Then
          ratlaig=slaig*(Exp(SucrosPlant(ip)%rgr*TimeFactor*temp_eff)-one_dp)      
       Else
          ratlaig=sla/PlantData%Area_ha_to_L2*max(0.0,ratwlvg)
       Endif
    Case Default
       If(SucrosPlant(ip)%lCeres==1) then        !if-else added by D.Farber
          ! Temperature-limited leaf expansion rate
          !modification N.Prolingheuer
          !If(dvs<0.3_dp .And. plai<0.75_dp) Then
          If(plai<SucrosPlant(ip)%ls) Then
             if( t .eq. null_dp .and. (dvs>0.2_dp.or.PlantGeometry%temp_sum_till/=null_dp))then
                PlantGeometry%temp_sum_till=PlantGeometry%temp_sum_till+temp_avg
             endif
             if(PlantGeometry%temp_sum_till>null_dp) then
                ratlaig=slaig*(Exp(SucrosPlant(ip)%rgr*temp_eff+(SucrosPlant(ip)%rsla*temp_eff) &
                     *(PlantGeometry%temp_sum_till+0.5*temp_eff))-one_dp)
             else
                ratlaig=slaig*(Exp(SucrosPlant(ip)%rgr*temp_eff)-one_dp)*TimeFactor
             end if
          !end modification
          ! Radiation-limited leaf expansion rate 
          else
             ratlaig=sla/PlantData%Area_ha_to_L2*max(0.0, ratwlvg)
          end if
       else
          If(dvs<0.3_dp .And. plai<0.75_dp) Then
             ratlaig=slaig*(Exp(SucrosPlant(ip)%rgr*temp_eff*TimeFactor)-one_dp)
          Else
             ratlaig=sla/PlantData%Area_ha_to_L2*max(0.0,ratwlvg)
          Endif
       end if
    End Select
    ! root penetration rates and root density growth rates
    rlnew=1000.0_dp*ratwrt*SucrosPlant(ip)%srw
    If(rlnew/=null_dp) Then
       ! calculate effective temperature for root growth
       !x=0.5_dp*(temp_min+temp_max)
       x = temp   !modification D.Farber (temp equal "0.5_dp*(temp_min+temp_max)" for daily calculation)
       If(temp_max<null_dp) Then
          x=null_dp
       Else If(temp_min<null_dp) Then
          x=Max(null_dp,&
               x+Max(null_dp,0.1583_dp*(temp_max-temp_min)-Abs(x)*0.4043_dp))
       Endif
       ! calculate the rooting depth 
       ! rootdepth = cm beneath  the soil surface
       !           ==> now in mm beneath  the soil surface
       PlantGeometry%rootdepth=Max(PlantGeometry%rootdepth &
            -Min(x*2.2_dp,18.0_dp)*TimeFactor, SucrosPlant(ip)%root_max)
    Endif
    ! calculation of the integrals
    ! plant development stage
    PlantGeometry%dvs=dvs+ratdvs
    !start modification A.Klosterhalfen   ! LINGRA
    If(ptype==Gras .or. ptype==c4Gras) Then
       !write (*,*) PlantGeometry%count
       if(PlantGeometry%temp_sum < SucrosPlant(ip)%temp_sum_crit) then
          PlantGeometry%count=0
       end if
       if(PlantGeometry%temp_sum > SucrosPlant(ip)%temp_sum_crit) then
          PlantGeometry%count=PlantGeometry%count+1
       end if
       !write (*,*) PlantGeometry%count
       if(PlantGeometry%temp_sum > SucrosPlant(ip)%temp_sum_crit .and. PlantGeometry%count==1) then
          PlantGeometry%dvs = PlantGeometry%dvs + 0.6
       end if
    End If
    !end modification
    If(ptype/=potatoe) PlantGeometry%dvs = Min(PlantGeometry%dvs, 2.0)
    ! weight in the different plant organs
    PlantGeometry%wlvg=PlantGeometry%wlvg+(ratwlvg-ratwlvd)
    if(PlantGeometry%wlvg < 0.) PlantGeometry%wlvg = 0.
    PlantGeometry%wlvd=PlantGeometry%wlvd+ratwlvd
    if(PlantGeometry%wlvd < 0.) PlantGeometry%wlvd = 0.
    PlantGeometry%wst=PlantGeometry%wst+ratwst
    if(PlantGeometry%wst < 0.) PlantGeometry%wst = 0.
    PlantGeometry%wso=PlantGeometry%wso+ratwso
    if(PlantGeometry%wso < 0.) PlantGeometry%wso = 0.
    PlantGeometry%wrt=PlantGeometry%wrt+ratwrt
    if(PlantGeometry%wrt < 0.) PlantGeometry%wrt = 0.
    PlantGeometry%wrtd=PlantGeometry%wrtd+ratwrtd
    if(PlantGeometry%wrtd < 0.) PlantGeometry%wrtd = 0.
    PlantGeometry%wcrn=PlantGeometry%wcrn+ratwcrn
    if(PlantGeometry%wcrn < 0.) PlantGeometry%wcrn = 0.
    ! leaf area characteristics 
    !PlantGeometry%slaig=Max(slaig+(ratlaig-ratlaid),null_dp)
    PlantGeometry%slaig=calculate_glai(ip)
    PlantGeometry%slaid=PlantGeometry%slaid+ratlaid*(1-SucrosPlant(ip)%debr_fac)   !added by D.Farber
    !PlantGeometry%slaid=PlantGeometry%slaid+InterpolateTab(SucrosPlant(ip)%tab(18),x)
    ! calculation of the actual Kc factor
    If(SucrosPlant(ip)%akctype==3) Then
       PlantGeometry%akc=calculate_akc(ip)
    Else
       If(SucrosPlant(ip)%akctype==1) Then
          x=PlantGeometry%dvs
       Else
          x=PlantData%dayOfYear
       End If
       PlantGeometry%akc=InterpolateTab(SucrosPlant(ip)%tab(11),x)
    Endif
  End Subroutine CalculateLAI


  Subroutine CalculateSurfaceFlux()
    ! calculate the flux through the surface
    Implicit None
    Real(dp) :: x, rain,rain2, intercap, lai, intercept, et0
    Integer :: ip

    rain=ClimateData%prec

    et0=ClimateData%et0
    rain2=rain   ! precipitation

    If(et0<null_dp) Then   ! dew
       rain2=rain2-et0
       et0=null_dp
    End If
    ip=PlantGeometry%index
    If(ip==0) Then   ! no plants
       PlantGeometry%eva_soil=et0
       PlantGeometry%eva_plants=null_dp
    Else
       ! calculation of the potential soil evaporation and the potential evapotranspiration of plants
       x=Max(null_dp, PlantGeometry%akc*et0)
       PlantGeometry%eva_soil=x * Exp(-0.6_dp*PlantGeometry%slaig)
       PlantGeometry%eva_plants=x - PlantGeometry%eva_soil
    End If
    ! calculation of interception (intercept) depending on capacity and total leaf area (lai)
    lai=PlantGeometry%slaig+PlantGeometry%slaid
    intercap=InterceptCapacity(rain2,lai)
    intercept=Min(intercap,PlantGeometry%pond_plant+rain2)
    PlantGeometry%pond_plant_new=intercept-PlantGeometry%pond_plant
    ! rain on soil (rain)
    PlantGeometry%rain=Max(null_dp,rain2-PlantGeometry%pond_plant_new)
    ! differentiation between transpiration and evaporation from interception
    x=Min(intercept, PlantGeometry%eva_plants)
    ! evaporation plants = transpiration + ponding plant
    PlantGeometry%transpiration=PlantGeometry%eva_plants - x
    intercept=intercept - x
    If(intercept>null_dp) Then
       ! try to evaporate the remaining interception "with energy leftover" from eva_soil
       x=Min(intercept, PlantGeometry%eva_soil)
       PlantGeometry%eva_soil_intercept=PlantGeometry%eva_soil-x
       intercept=intercept - x
    Else
       PlantGeometry%eva_soil_intercept=PlantGeometry%eva_soil
    End If
    PlantGeometry%pond_plant=intercept
    ! flux=precipitation-evaporation
    PlantGeometry%flux=PlantGeometry%eva_soil_intercept - PlantGeometry%rain
    If(ip>0) Then   ! plants
       ! calculation of p2 for the reduction
       ! of the potential root water uptake (smx)
       x=PlantGeometry%transpiration
       If(x < 0.1*TimeFactor) Then
          PlantGeometry%p2 = Plant(ip)%p2l
       Else If(x < 0.5*TimeFactor) Then
          PlantGeometry%p2 = Plant(ip)%p2h + &
               ((0.5*TimeFactor-x)/(0.4*TimeFactor))*(Plant(ip)%p2l-Plant(ip)%p2h)
       Else
          PlantGeometry%p2 = Plant(ip)%p2h
       End If
    End If
  End Subroutine CalculateSurfaceFlux


  Subroutine SetPlantTime(SimTime)
    Use TimeData, Only: simtime_to_iday, day_of_year
    Implicit None
    Real(dp), Intent(in) :: SimTime
    Integer :: i, j, t, ptype
    Integer :: oldseason

    ! calculate day of the year
    PlantData%dayOfYear=day_of_year(simtime_to_iday(simtime))
    ! end of year events
    If(equal(PlantData%dayOfYear,365.0_dp)) Then
       i=PlantGeometry%index
       If(Plant(i)%season>0) Then
          ptype=Plant(i)%type
          If(ptype==Gras .Or. ptype==c4Gras) Then
             call ResetGrass()
          Endif
       Endif
    Endif
    ! calculate season number
    Do i=1, PlantData%NoTypes
       oldseason=Plant(i)%season
       Plant(i)%season=0
       Do j=1, SucrosPlant(i)%no_dates
          t=simtime_to_iday(simTime)
          ! for hourly-calc harvest is at 12am. Due to FP-aritm the upper bound is 12.1
          If(t>=SucrosPlant(i)%emergence(j) .and. t<=SucrosPlant(i)%harvest(j)+int(12.1/24./Timefactor)) Then
             Plant(i)%season=j
             Exit
          End If
       End Do
       If(oldseason>0 .and. Plant(i)%season==0) Then
          ! after harvest
          Call ResetPlant(i, SimTime, oldseason)
       Else If(oldseason==0 .and. Plant(i)%season>0) Then
          ! plant or emergence
          Call InitPlant(i)
       else if(oldseason>0 .and. Plant(i)%season>oldseason) then
          call ResetPlant(i, SimTime, oldseason)
       End If
    End Do
  End Subroutine SetPlantTime


  Subroutine CalculateRootExtraction()
    Implicit None
    PlantGeometryArray%rex = PlantGeometry%transpiration*PlantGeometryArray%rrd/PlantGeometryArray%dz
  End Subroutine CalculateRootExtraction

  Subroutine CalculateRelativeRootdist()
    Use Geometry, Only: beta
    Implicit None
    Real(dp) :: SimTime, rna, rootdepth
    Integer :: i, ip

    PlantGeometryArray%rrd=0 ! init
    ! determine zone from 0 to rna where is no extraction because of senescence
    SimTime=PlantData%dayOfYear
    Do i=1,PlantData%NoTypes
       If( between(SimTime,SucrosPlant(i)%senescence_start,SucrosPlant(i)%senescence_end) ) Then
          ! linear interpolation between start and end time of senescence
          ! the max value for the inactive zone is a plant parameter
          Plant(i)%rna=SucrosPlant(i)%rna_max * &
               (SimTime-SucrosPlant(i)%senescence_start) / &
               (SucrosPlant(i)%senescence_end-SucrosPlant(i)%senescence_start)
       Else
          Plant(i)%rna=null_dp
       End If
    End Do
    !raus  ip=PlantData%NoTypes
    !raus  Do i=1,PlantData%NoSimpleTypes
    !raus     ip=ip+1
    !raus     Plant(ip)%rna=InterpolateTab(SimplePlant(i)%tab(6),SimTime)
    !raus  End Do
    ! calculate the root extraction rates
    ip=PlantGeometry%index
    If(ip==0) Return
    If(Plant(ip)%season==0) Return
    rootdepth=PlantGeometry%rootdepth
    rna=Plant(ip)%rna
    If(rna<=rootdepth) Return   ! rna is below depth of the roots
    Call DistributeRoots(Plant(ip)%root_distribution, PlantGeometryArray%rrd, rna, rootdepth, &
         PlantGeometryArray%z)
    beta=PlantGeometryArray%rrd
  End Subroutine CalculateRelativeRootdist


  Subroutine DistributeRoots(tab, roots, rna, rootdepth, z)
    Use Geometry, only: NumNP
    Implicit None
    Type(TabType), intent(in) :: tab
    Real(dp), intent(inout) :: roots(:)
    Real(dp), intent(in) :: rna, rootdepth
    Real(dp), intent(in) :: z(:)
    Real(dp) :: rl, z1,z2, dz,dz2
    Integer :: i,i1

    ! from surface to rna root extraction is 0,
    ! from rna to rootdepth root extraction is the transpiration
    ! distributed by the root density distribution (w) divided by dz
    ! rex = transpiration * w / dz

    ! determine compartment for which the lower boundary is below rna
    !
    !  | . .  |  .  .   .    |      .       .      .        |
    ! top    rna          rootdepth                       bottom
    !         <    rl        >
    !  0     -50           -300                           -1000
    roots=0 ! init
    Do i=1,NumNP
       If(rna>z(i)) Exit
    End Do
    If(i>NumNP) Return
    i1=Min(Max(i,2),NumNP-1)
    rl=rna-rootdepth
    z1=null_dp
    dz=0.5*(z(i1-1)-z(i1))
    z2=(rna-z(i1-1)+dz)/rl
    ! first node
    If(rna>0.5*(z(i1)+z(i1-1))) Then
       roots(i1-1)=IntegrateTab(tab,z1,z2)
    Endif
    Do i=i1,NumNP-1
       z1=z2
       dz=0.5*(z(i-1)-z(i+1))
       dz2=0.5*(z(i)-z(i+1))
       z2=( rna-z(i)+dz2 )/rl
       roots(i)=IntegrateTab(tab,z1,z2)
       If(rootdepth>=z(i+1)) Exit
    End Do
    ! last node
    If(i<NumNP) Then
       If(rootdepth<0.5*(z(i)+z(i+1))) Then
          z1=z2
          dz=0.5*(z(i)-z(i+1))
          z2=(rna-z(i+1))/rl
          roots(i+1)=IntegrateTab(tab,z1,z2)
       Endif
    Endif
  End Subroutine DistributeRoots


  Subroutine CalculateActualTranspiration(delta_t,vTop)
    Use Geometry, Only: NumNP
    Implicit None
    Real(dp), Intent(in) :: delta_t,vTop
    Integer :: i,ip

    If(debug) Print *,'CalculateActualTranspiration'
    ! calculate the actual transpiration = sum over the root sink terms
    ip=PlantGeometry%index
    If(ip/=0) then
       If(Plant(ip)%season/=0) then
          Do i=1,NumNP
             !mass=rex*volume*delta_t
             If(PlantGeometryArray%rex(i)>null_dp) Then
                PlantGeometry%transpiration_act= &
                     PlantGeometry%transpiration_act+ &
                     PlantGeometryArray%rex_act(i)*delta_t*PlantGeometryArray%dz(i)
             End If
          End Do
       end if
     end if
    PlantGeometry%fluxsum=PlantGeometry%fluxsum+delta_t*vTop*PlantData%UnitFactor

  End Subroutine CalculateActualTranspiration


  Subroutine ResetActualTranspiration()
    Implicit None
    ! reset accumulated variables
    PlantGeometry%transpiration_act=null_dp
    PlantGeometry%fluxsum=null_dp
  End Subroutine ResetActualTranspiration


  Logical Function between(x, a, b)
    ! test if x is between a and b
    Real(dp), Intent(in) :: x, a, b
    between=( (a<=x .and. x<=b) .or. (b<=x .and. x<=a) )
  End Function between

  
  Real(dp) Function daylength(day, lat)
    Implicit None
    Real(dp), Intent(in) :: day, lat
    Real(dp), Parameter :: pi1=4.0_dp*Atan(1.0_dp)
    Real(dp), Parameter :: rad=pi1/180.0_dp
    Real(dp) dec,sinld,cosld,sunrise,sunset
    
    dec = -23.4_dp*Cos(2.0_dp*PI*(day+10.0_dp)/365.0_dp)
    sinld= Sin(dec*rad)*Sin(lat*rad)
    cosld= Cos(dec*rad)*Cos(lat*rad)
    sunrise = 12-Acos(-Tan(lat*rad)*Tan(dec*rad))/15/rad
    sunset = 12+Acos(-Tan(lat*rad)*Tan(dec*rad))/15/rad
    daylength = sunset-sunrise
  End Function daylength


  Subroutine ReadPlantData(fn)
    Use TimeData, Only: SetStartDate, date_to_iday, iday_to_date
    Use Variables, Only: co2_fluxes, respiration, maint_growth, waterstress, dailyCalculation,&
         rootexu, rootdeath, PlantsExist, lNitrogen, lPhosphorus, P_uptake_method
    Implicit None
    Character(*), Intent(in) :: fn
    Integer :: ios, err, i, j, k, yy, mm, dd
    !raus   Integer :: ik
    Character :: line*256
    Character(len=1), Allocatable :: switches(:)
    Integer, Parameter :: lbuffer=50
    Real(dp) :: buffer(lbuffer)
    Integer :: lCeresTmp   !added by D.Farber
    Integer :: libuffer
    Integer, Allocatable :: ibuffer(:)
    Integer, allocatable :: datebuffer(:)
    Integer :: version, ptype, year
    Integer, Allocatable :: years(:)
    Type(FarquharType), Pointer :: fq=>Null()

    If(.Not.PlantsExist) Then
       PlantData%NoTypes=0
       !raus     PlantData%NoSimpleTypes=0
       Call WriteOutput('Plants disabled')
       Return
    End If
    If(lPhosphorus) Then
       PlantTabs=18
    Elseif(lNitrogen) Then
       PlantTabs=17
    Else
       PlantTabs=12
    Endif
    libuffer=PlantTabs+4
    Allocate(ibuffer(libuffer))
    Call WriteOutput('Plants enabled')
    Open(32,file=Trim(fn),status='OLD',iostat=ios)
    If(ios/=0) Call WriteError( 'Error while opening plant input file: ' // Trim(fn))
    Read(32,2,end=1) line
    If(Trim(line)/='soilco2 plant input') &
         Call WriteError( 'Wrong file format for the plant input file: ' // Trim(fn))
    Call WriteOutput('Reading plant input file: '//Trim(fn))
    Read(32,*,End=1) version
    ! version:  2=alt  3=Farquhar
    If(version<2) Call WriteError('Ungueltige Versionsnummer in: '//Trim(fn))

    Read(32,*,end=1)
    ! co2_fluxes, respiration, maint_growth, waterstress, rootexu, rootdeath, harvestresidues
    Read(32,2,end=1) line
    i=wordsf(line)
    If(i<7 .Or. i>9) &
         Call WriteError('Invalid number of switches in line 4 of plant.in')
    Allocate(switches(i))
    Read(line,*,End=1) switches
    Call s2bool(switches(1), co2_fluxes)
    Call s2bool(switches(2), respiration)
    Call s2bool(switches(3), maint_growth)
    Call s2int(switches(4), waterstress)
    Call s2bool(switches(5), rootexu)
    Call s2bool(switches(6), rootdeath)
    Call s2bool(switches(7), harvestresidues)
    If(i>7) Then 
       Call s2bool(switches(8), Farquhar_used)
       If(i>8) Then
          Call s2int(switches(9), P_uptake_method)
          If(P_uptake_method/=1 .And. P_uptake_method/=2) &
               Call WriteError('Invalid P_uptake_method in plant.in')
       Else
          P_uptake_method=1
       End If
    Else
       Farquhar_used=.False.
    Endif
    read(32,*, end=1) dailyCalculation
    if(dailyCalculation) then
       TimeFactor = 1.0_dp
    else
       TimeFactor = 1.0_dp/24.0_dp
    end if
    Read(32,*,end=1) yy,mm,dd
    Call SetStartDate(yy,mm,dd)
    Read(32,*,end=1) PlantData%NoTypes
    If(PlantData%NoTypes<0 .Or. PlantData%NoTypes>MaxPlant) &
         Call WriteError( 'Wrong number of plant types in file: ' // Trim(fn))
    If(PlantData%NoTypes>0) Then
       Allocate (SucrosPlant(PlantData%NoTypes), stat=err)
       If(err/=0) Call AllocateError('SucrosPlant(PlantData%NoTypes)')
    End If
    Do i=1,PlantData%NoTypes
       Allocate (SucrosPlant(i)%tab(PlantTabs), stat=err)
       If(err/=0) Call AllocateError('SucrosPlant(i)%tab()')
    Enddo
    Read(32,*,end=1) 
    Read(32,*,end=1) PlantData%InterceptModel
    If(PlantData%InterceptModel<1 .or. PlantData%InterceptModel>2) &
         Call WriteError( 'Wrong input for the intercept model in file: ' // Trim(fn))
    !raus  Read(32,*,end=1) PlantData%NoSimpleTypes
    !raus  If(PlantData%NoSimpleTypes<0 .or. PlantData%NoSimpleTypes>MaxPlant) &
    !raus     Call WriteError( 'Wrong number of plants for the simple model in file: ' // Trim(fn))
    Read(32,*,end=1) PlantData%latitude
    !raus  PlantData%NoTypes=PlantData%NoTypes+PlantData%NoSimpleTypes
    PlantData%NoTypes=PlantData%NoTypes
    Allocate (Plant(0:PlantData%NoTypes), stat=err)
    If(err/=0) Call AllocateError('Plant(PlantData%NoTypes)')
    If(PlantData%NoTypes>0) Write(*,*) &
         'Number of plants for the SUCROS model: ', PlantData%NoTypes

    Do i=1,PlantData%NoTypes
       Read(32,2,end=1)
       Read(32,*,end=1) j
       If(j<1 .or. j>MaxPlant) &
            Call WriteError( 'Wrong plant type in file: ' // Trim(fn))
       ibuffer(PlantTabs+3)=j
       Read(32,*,end=1) ibuffer(1:PlantTabs)
       Read(32,*,end=1) ibuffer(PlantTabs+1)
       Read(32,*,end=1) ibuffer(PlantTabs+2)
       Read(32,*,end=1) ibuffer(PlantTabs+4)
       Do j=1,PlantTabs
          Allocate (SucrosPlant(i)%tab(j)%data(ibuffer(j),2), stat=err)
          If(err/=0) Call AllocateError('SucrosPlant(i)%tab(j)%data')
          SucrosPlant(i)%tab(j)%nrows=ibuffer(j)
          SucrosPlant(i)%tab(j)%ncols=2
       End Do
       If(ibuffer(PlantTabs+2)/=lbuffer) &
            Call WriteError( 'Wrong number of parameters in file: ' // Trim(fn))
       ptype=ibuffer(PlantTabs+3)
       Plant(i)%type=ptype
       If(ptype==Maize .Or. ptype==c4Gras) Then
          Plant(i)%c3=.False.
       Else
          Plant(i)%c3=.True.
       End If
       !raus     Plant(i)%model=1
       SucrosPlant(i)%no_dates=ibuffer(PlantTabs+1)
       SucrosPlant(i)%akctype=ibuffer(PlantTabs+4)
       If(SucrosPlant(i)%akctype<1 .or. SucrosPlant(i)%akctype>3) &
            Call WriteError( 'Wrong Kc type in file: ' // Trim(fn))
       Allocate (SucrosPlant(i)%emergence(SucrosPlant(i)%no_dates), stat=err)
       If(err/=0) Call AllocateError('SucrosPlant(i)%emergence')
       Allocate (SucrosPlant(i)%harvest(SucrosPlant(i)%no_dates), stat=err)
       If(err/=0) Call AllocateError('SucrosPlant(i)%harvest')
       Allocate (datebuffer(1:3*2*(SucrosPlant(i)%no_dates)), stat=err)
       If(err/=0) Call AllocateError('datebuffer')
       Allocate(SucrosPlant(i)%last_Date_in_year(SucrosPlant(i)%no_dates))
       SucrosPlant(i)%last_Date_in_year=.False.
       If(ptype==potatoe) Then
          Read(32,*,End=1) buffer(1:2),senescence_start_temp
       Else
          Read(32,*,End=1) buffer(1:2)
       Endif
       If(waterstress==0) Then
          Read(32,*,End=1)
          buffer(3:7)=0
       Else If(waterstress==3) Then
          Read(32,*,End=1) buffer(3:5)
          buffer(6:7)=0
       Else
          Read(32,*,End=1) buffer(3:7)
       Endif
       !modification N.Prolingheuer/D.Farber
       Read(32,*,end=1) lCeresTmp, buffer(8:19)   ! Ceres: temperatures
       SucrosPlant(i)%lCeres = lCeresTmp          ! original or new approach? (added by D.Farber)
       Read(32,*,end=1) buffer(20:22)             ! Ceres: photoperiod

       If(ptype==Gras .Or. ptype==c4Gras) Then   ! LINGRA
          Read(32,*,End=1) SucrosPlant(i)%dvsReset, SucrosPlant(i)%dvsLow, SucrosPlant(i)%dvsHigh, &
                  SucrosPlant(i)%storageFacMax, SucrosPlant(i)%storageMax, SucrosPlant(i)%Tstorage, &
                  SucrosPlant(i)%cutLai, SucrosPlant(i)%temp_sum_crit
          Read(32,*,End=1) SucrosPlant(i)%harvestType
          If(SucrosPlant(i)%harvestType<1 .Or. SucrosPlant(i)%harvestType>2) &
             Call WriteError('Error in plant.in: Invalid grass input.')
       Endif

       Read(32,*,end=1) buffer(23:25)             ! Ceres: development rates
       !end modification

       Do j=26,lbuffer
          Read(32,*,End=1) buffer(j)
       End Do
       SucrosPlant(i)%senescence_start=buffer(1)
       SucrosPlant(i)%senescence_end=buffer(2)
       If(waterstress==3) Then
          ! hx_min [L] is the threshold collar water potential that triggers
          ! a reduction of Tact, as compared to Tpot.
          Plant(i)%hx_min=buffer(3)
          Plant(i)%Krs=buffer(4)   ! Equivalent conductance of the whole root system [cm³/hPa/day]
          Plant(i)%Kcomp=buffer(5) ! Compensatory RWU conductance of the root system [cm³/hPa/day]
       Else If(waterstress==2) Then
          Plant(i)%p0=buffer(3)
          Plant(i)%p1=buffer(4)
          Plant(i)%p2h=buffer(5)
          Plant(i)%p2l=buffer(6)
          Plant(i)%p3=buffer(7)
       Else
          Plant(i)%p0=-Abs(buffer(3))
          Plant(i)%p1=-Abs(buffer(4))
          Plant(i)%p2h=-Abs(buffer(5))
          Plant(i)%p2l=-Abs(buffer(6))
          Plant(i)%p3=-Abs(buffer(7))
       End If
       SucrosPlant(i)%Tmin_v1=buffer(8)
       SucrosPlant(i)%Topt_v1=buffer(9)
       SucrosPlant(i)%Tmax_v1=buffer(10)
       SucrosPlant(i)%Tmin_v2=buffer(11)
       SucrosPlant(i)%Topt_v2=buffer(12)
       SucrosPlant(i)%Tmax_v2=buffer(13)
       SucrosPlant(i)%Tmin_r=buffer(14)
       SucrosPlant(i)%Topt_r=buffer(15)
       SucrosPlant(i)%Tmax_r=buffer(16)
       SucrosPlant(i)%Tmin_vn=buffer(17)
       SucrosPlant(i)%Topt_vn=buffer(18)
       SucrosPlant(i)%Tmax_vn=buffer(19)
       SucrosPlant(i)%Popt=buffer(20)
       SucrosPlant(i)%Pcrit=buffer(21)
       SucrosPlant(i)%omega=buffer(22)
       SucrosPlant(i)%rmax_v1=buffer(23)
       SucrosPlant(i)%rmax_v2=buffer(24)
       SucrosPlant(i)%rmax_r=buffer(25)
       SucrosPlant(i)%rna_max=-Abs(buffer(26))
       If(SucrosPlant(i)%senescence_start>=SucrosPlant(i)%senescence_end) Then
          SucrosPlant(i)%senescence_start=null_dp
          SucrosPlant(i)%senescence_end=one_dp
          SucrosPlant(i)%rna_max=null_dp
       End If
       SucrosPlant(i)%root_max=-Abs(buffer(27))
       SucrosPlant(i)%root_init=-Abs(buffer(28))
       SucrosPlant(i)%exu_fac = buffer(29)
       SucrosPlant(i)%deathfacMax = buffer(30)
       SucrosPlant(i)%nsl=buffer(31)
       SucrosPlant(i)%rgr=buffer(32)
       SucrosPlant(i)%tempbase=buffer(33)
       SucrosPlant(i)%sla=buffer(34)
       SucrosPlant(i)%rsla=-abs(buffer(35))   !added by N.Prolingheuer
       SucrosPlant(i)%amx=buffer(36)
       SucrosPlant(i)%eff=buffer(37)
       SucrosPlant(i)%rkdf=buffer(38)
       SucrosPlant(i)%scp=buffer(39)
       SucrosPlant(i)%rmainso=buffer(40)
       SucrosPlant(i)%asrqso=buffer(41)
       SucrosPlant(i)%tempstart=buffer(42)
       SucrosPlant(i)%debr_fac=buffer(43)
       SucrosPlant(i)%ls=buffer(44)
       SucrosPlant(i)%rlaicr=buffer(45)
       SucrosPlant(i)%eai=buffer(46)
       SucrosPlant(i)%rmatr=buffer(47)
       SucrosPlant(i)%ssl=buffer(48)
       SucrosPlant(i)%srw=buffer(49)
       SucrosPlant(i)%slaid_off=buffer(50)
       Allocate(fq)
       Plant(i)%Far => fq
       ! read input that can also be missing
       Read(32,2, End=1) line ! comment line
       line=to_lower(line)
       If(line(1:9)=='#farquhar' .Or. line(1:10)=='# farquhar') Then
          If(Farquhar_used) Then
             Read(32,*, End=1) fq%vcmax25N
             Read(32,*, End=1) fq%m
             ! Farquhar constants
             fq%K_c25=25.0
             fq%K_o25=30000.0
             fq%a_kc=2.1
             fq%a_ko=1.2
             fq%temp_fw=273.15_dp
             fq%Boltzmann=1.38065E-23_dp
             fq%Avogadro=6.02214E26_dp
             fq%R_gas=fq%Boltzmann*fq%Avogadro
             fq%CN_L=25_dp
             fq%a_R25=60_dp
             fq%F_NR=7.16_dp
             fq%SLA_0=0.03_dp
             fq%PSI_O= -74000_dp
             fq%PSI_C=-275000_dp
             If(ptype==C4Gras .Or. ptype==Gras) Then
                fq%F_LNR=0.09_dp
             Else
                fq%F_LNR=0.1_dp
             End If
             ! the quantum efficiency ( µ mol CO2 per µ mol photons) (Table 8.1)
             If(ptype==C4Gras) Then
                fq%alfa=0.04_dp
             Else
                fq%alfa=0.06_dp
             End If
             ! 21.6. is day 172 of year
             fq%daylength_max=daylength(172.5_dp,PlantData%latitude)
          Else
             Read(32,2, End=1)
             Read(32,2, End=1)
          End If
       Else
          If(Farquhar_used) Then
             Call WriteError('No input for Farquhar in plant.in')
          Else
             Backspace(32)
          Endif
       Endif
       Read(32,2, End=1) line ! comment line
       line=to_lower(line)
       If(line(1:11)=='#phosphorus' .or. line(1:12)=='# phosphorus') Then
          If(P_uptake_method==2) Then
             Read(32,*,End=1)  SucrosPlant(i)%bp
          Else
             Read(32,*,End=1)
          Endif
       Else
          If(P_uptake_method==2) Then
             Call WriteError('Missing phosphorus input in plant.in')
          Else
             Backspace(32)
          Endif
       Endif
       ! emergence und harvest date(s)
       Read(32,2, End=1)
       If(ptype==Gras .Or. ptype==c4Gras)Then   ! LINGRA
          ! for lingra there is only one emergence date, the 1st of January of the first year
          ! all of the dates are hrvest dates.
          Read(32,*,end=1) datebuffer(1:(3*SucrosPlant(i)%no_dates))
          Allocate(years(SucrosPlant(i)%no_dates))
          Do j=1,SucrosPlant(i)%no_dates
             year=datebuffer(1+3*(j-1))
             years(j)=year
             SucrosPlant(i)%emergence(j)=date_to_iday(year, 1, 1) ! only 1 emergence date on 1th of January
             SucrosPlant(i)%harvest(j)=date_to_iday(year,datebuffer(2+3*(j-1)),datebuffer(3+3*(j-1)))
          End Do
          SucrosPlant(i)%last_date_in_year(SucrosPlant(i)%no_dates)=.True. ! last year
          Do j=SucrosPlant(i)%no_dates-1,1,-1
             !write(*,*) "em", datebuffer(1), SucrosPlant(i)%emergence(j), SucrosPlant(i)%harvest(j)
             If(years(j)<year) Then
                SucrosPlant(i)%last_date_in_year(j)=.True.
                If(year-1 == years(j) .And. SucrosPlant(i)%harvestType==1) &
                     SucrosPlant(i)%emergence(j+1)=SucrosPlant(i)%emergence(j)
                year=years(j)
             Endif
          End Do
          Write(*,*) "emergence", SucrosPlant(i)%emergence
          Write(*,*) "harvest  ", SucrosPlant(i)%harvest
          Write(*,*) "lastdate ", SucrosPlant(i)%last_date_in_year
          Deallocate(years)
       Else
          Read(32,*,End=1) datebuffer(1:(3*2*SucrosPlant(i)%no_dates))
          Do j=1,SucrosPlant(i)%no_dates
             SucrosPlant(i)%emergence(j)=date_to_iday(datebuffer(1+6*(j-1)),datebuffer(2+6*(j-1)),datebuffer(3+6*(j-1)))
             SucrosPlant(i)%harvest(j)=date_to_iday(datebuffer(4+6*(j-1)),datebuffer(5+6*(j-1)),datebuffer(6+6*(j-1)))
          End Do
       End If
       ! Tables
       Do j=1,PlantTabs
          Read(32,2,end=1)
          Do k=1,SucrosPlant(i)%tab(j)%nrows
             Read(32,*,end=1) SucrosPlant(i)%tab(j)%data(k,1:2)
          End Do
       End Do
       ! copy table with root distribution into Plant
       Plant(i)%root_distribution%nrows=SucrosPlant(i)%tab(RootTab)%nrows
       Plant(i)%root_distribution%ncols=2
       Plant(i)%root_distribution%data=>SucrosPlant(i)%tab(RootTab)%data
       Deallocate(datebuffer)
       If(SucrosPlant(i)%akctype==3) Then
          SucrosPlant(i)%akcMin=Minval(SucrosPlant(i)%tab(11)%Data(:,2))
          Print *,'akcmin=',SucrosPlant(i)%akcMin
       Endif
    End Do

    Do i=1,PlantData%NoTypes
       If( Plant(i)%p0<Plant(i)%p1 .or. Plant(i)%p1<Plant(i)%p2h .or. &
            Plant(i)%p2h<Plant(i)%p2l .or. Plant(i)%p2l<Plant(i)%p3 ) &
            Call WriteError( 'Wrong input pressure head values in plant file: ' // Trim(fn))
    End Do

    Close(32)

    Return
1   Call WriteError( 'Unexpected end of file in: '// Trim(fn))
2   Format(A)
  End Subroutine ReadPlantData


  Real(dp) Function InterpolateTab(tab, x)
    Implicit None
    Type(TabType), Intent(in) :: tab
    Real(dp), Intent(in) :: x
    Integer :: i
    If(x<=tab%data(1,1)) Then
       InterpolateTab=tab%data(1,2)
    Else If(x>=tab%data(tab%nrows,1)) Then
       InterpolateTab=tab%data(tab%nrows,2)
    Else
       Do i=2,tab%nrows
          If(x<=tab%data(i,1)) Exit
       End Do
       InterpolateTab=(x-tab%data(i-1,1))/(tab%data(i,1)-tab%data(i-1,1))* &
            (tab%data(i,2)-tab%data(i-1,2))+tab%data(i-1,2)
    End If
  End Function InterpolateTab


  Real(dp) Function IntegrateTab(tab, a,b)
    Implicit None
    Type(TabType), Intent(in) :: tab
    Real(dp), Intent(in) :: a,b

    Integer :: i
    Do i=1,tab%nrows
       If(a<=tab%data(i,1)) Exit
    End Do
    If(i>tab%nrows .or. b==a) Then
       IntegrateTab=null_dp
    Else If(b<=tab%data(i,1)) Then
       If(i<=1) Then
          IntegrateTab=null_dp
       Else
          IntegrateTab=0.5_dp*(InterpolateTab(tab,a)+InterpolateTab(tab,b))*(b-a)
       End If
    Else
       If(i>1) Then
          IntegrateTab=0.5_dp*(tab%data(i,2)+InterpolateTab(tab,a))*(tab%data(i,1)-a)
       Else
          IntegrateTab=null_dp
       End If
       Do i=i+1,tab%nrows
          If(b<=tab%data(i,1)) Exit
          IntegrateTab=IntegrateTab+0.5_dp*(tab%data(i,2)+tab%data(i-1,2)) *&
               (tab%data(i,1)-tab%data(i-1,1))
       End Do
       If(i<=tab%nrows) IntegrateTab=IntegrateTab+0.5_dp*&
            (InterpolateTab(tab,b)+tab%data(i-1,2)) * (b-tab%data(i-1,1))
    End If
  End Function IntegrateTab


  Subroutine NormalizeRootDistribution(tab)
    Implicit None
    Type(TabType), Intent(inout) :: tab
    Integer :: j, n1
    Real(dp) :: z, xmin

    ! scale the x-values of the root distribution to 0-1
    n1=tab%nrows
    If(n1<=1) Call WriteError('Invalid distribution of roots')
    xmin=tab%data(1,1)
    z=tab%data(n1,1)-xmin
    If(z==null_dp) Call WriteError('Invalid distribution of roots')
    Do j=1,n1
       tab%data(j,1)=(tab%data(j,1)-xmin)/z
    End Do
    ! Integrate the root distribution curve
    ! integral = 0.5 * sum( (y(i)+y(i-1)) * (x(i)-x(i-1)) )
    z=null_dp
    Do j=2,n1
       z=z+(tab%data(j,2)+tab%data(j-1,2)) * (tab%data(j,1)-tab%data(j-1,1))  
    End Do
    If(z==null_dp) Call WriteError('Invalid integral of the roots distribution')
    z=0.5_dp*z
    ! set integral to 1
    Do j=1,n1
       tab%data(j,2)=tab%data(j,2)/z
    End Do
  End Subroutine NormalizeRootDistribution


  Function InterceptCapacity(rain, lai)
    Implicit None
    Real(dp) :: InterceptCapacity
    Real(dp), Intent(in) :: rain, lai

    If(PlantData%InterceptModel==1) Then
       InterceptCapacity=0.2_dp*lai
    Else
       InterceptCapacity=Max(null_dp, -0.42_dp+0.245_dp*rain+0.2_dp*lai+&
            0.0271_dp*rain*lai-0.0111_dp*rain*rain-0.0109_dp*lai*lai)
    End If
  End Function InterceptCapacity


  Subroutine WritePlantOutput(SimTime)

    Use TimeData, Only: iday_to_date, simtime_to_iday
    Use Datatypes, Only: dp, null_dp
    Use Variables, Only: alphaAvg
    Implicit none
    Real(dp), Intent(in) :: SimTime
    Real(dp), Save :: lasttime=null_dp
    Real(dp) :: flux, evapOut
    Logical, Save :: firstCall=.true.
    Integer :: ios,yy,mm,dd, ptype


    character(len=*), parameter :: oPlants = "plants.out"
    character(len=*), parameter :: oPlants2 = "plants2.out"

    If(debug) Print *,'WritePlantOutput'
    if(firstCall) Then
       firstCall=.false.
       Open(60,file=Trim(oPlants),Status='REPLACE',ACTION='WRITE',&
            FORM='FORMATTED',iostat=ios)
       If(ios/=0) Call WriteError( 'Error while opening plants output file: '&
            // Trim(oPlants))
       Write(60,2) '#plants'
       Write(60,2) '# 1=year','# 2=month','# 3=day',&
            '# 4=simulation time',&
            '# 5=depth of the roots [mm]','# 6=actual Kc',&
            '# 7=effective temperature sum [dC]','# 8=development stage (dvs)',&
            '# 9=total gross assimilation [kg CH2O/L2]','#10=ponding on plants/interception [L/T]',&
            '#11=ponding on soil [L]',&
            '#12=total green leaf area (slaig)',&
            '#13=total dead leaf area (slaid)',&
            '#14=dry weight of green leaves (wlvg) [kg DM/L2]',&
            '#15=dry weight of dead leaves (wlvd) [kg DM/L2]',&
            '#16=dry weight of storage organs (wso) [kg DM/L2]',&
            '#17=dry weight of stems (wst) [kg DM/L2]',&
            '#18=dry weight of roots (wrt) [kg DM/L2]',&
            '#19=dry weight of dead roots (wrtd) [kg DM/L2]',&
            '#20=dry weight of crowns (wcrn) [kg DM/L2]',&
            '#21=C storage of grass [kg C/L2]',&
            '#22=plant_type'
       Write(60,3) 'year month','day time',&
            'rootdepth','akc','temp_sum',&
            'dvs','cgphot','pond_plant','pond_soil','slaig','slaid','wlvg',&
            'wlvd','wso','wst','wrt','wrtd','wcrn','storage','plant_type'
       Open(61,file=Trim(oPlants2),Status='REPLACE',ACTION='WRITE',&
            FORM='FORMATTED',iostat=ios)
       If(ios/=0) Call WriteError( 'Error while opening plants output file: '&
            // Trim(oPlants2))
       Write(61,2) '#plants2'
       Write(61,2) '# 1=year','# 2=month','# 3=day',&
            '# 4=simulation time',&
            '# 5=potential soil evaporation [mm/T]',&
            '# 6=potential transpiration of plants (including evaporation of interception) [mm/T]',&
            '# 7=evaporation of interception due "with energy leftover" from eva_soil [mm/T]',&
            '# 8=evaporation of interception plants [mm/T]',&
            '# 9=actual soil evaporation [mm/T]',&
            '#10=actual transpiration [mm/T]',&
            '#11=potential flux [mm/T]',&
            '#12=actual flux [mm/T]',&
            '#13=alphaAvg'
       Write(61,3) 'year month','day time',&
            'eva_pot','trans_pot','i_soil','i_plants',&
            'eva_act','transp_act','flux_pot','flux_act', 'alphaAvg'
       If(farquhar_used) Then
          Open(63,file='sif.out',Status='REPLACE',ACTION='WRITE',&
               FORM='FORMATTED',iostat=ios)
          If(ios/=0) Call WriteError( 'Error while opening plants output file: sif.out')
          Write(63,22) '','time','A','apar_pho','Jo','Je','Phi_p','kf','kd','kn','fluorescence','fluorescence_755nm','fact'
!          Write(63,22) '#','(t)','(umol co2/m2/s)','(umol photons/m2/s)','(umol co2/m2/s)','(umol co2/m2/s)',&
!               '(-)','(1/s)','(1/s)','(1/s)','(umol photons/m2/s)','(W/m2/sr/um)'
       
          Open(64,file='faquhar_Photosynthese.out',Status='REPLACE',ACTION='WRITE',&
               FORM='FORMATTED',iostat=ios)
          If(ios/=0) Call WriteError( 'Error while opening plants output file: faquhar_Photosynthese.out')

          Write(64,37) 'time','fq%w_c','fq%w_j_sh','fq%w_j_sl','fq%w_e',&
            'fq%K_c','fq%K_o','fq%c_i','fq%gamma','parlsh','parlpp','plai','clustf','rkbl','sinb','asssh','asssl','fslla',&
            'temp','patm_pa','fq%o_i','Units_L(4)/Units_L(Unit_L_Input)','fact','vcmax',&
            'alphaAvg','frdf','dsinb','apar','pardf','pardr','atmtr','detr','ClimateData%global_rad','detr*dsinb','sinb','dsinbe'
22        Format(A1,A8,20A21)
37        Format(A4,33A13)  
       Endif
    Else
       Open(60,file=oPlants,Status='OLD',ACTION='WRITE',&
            FORM='FORMATTED',POSITION='APPEND')
       Open(61,file=oPlants2,Status='OLD',ACTION='WRITE',&
            FORM='FORMATTED',POSITION='APPEND')
    End If
    Call iday_to_date(simtime_to_iday(SimTime),yy,mm,dd)

    !write(*,*) "time", dd, mm, yy

    ptype=PlantGeometry%index
    If(ptype>0) ptype=Plant(ptype)%type
    flux=PlantGeometry%fluxsum/(Simtime-lasttime)   !raus /PlantGeometry%area
    Write(60,'(I4,2I3,F9.2,1P,17E13.5,0P,I12)') yy,mm,dd, SimTime,&
         PlantGeometry%rootdepth,PlantGeometry%akc,&
         PlantGeometry%temp_sum,PlantGeometry%dvs,&
         PlantGeometry%cgphot,PlantGeometry%pond_plant,&
         PlantGeometry%pond_soil,PlantGeometry%slaig,&
         PlantGeometry%slaid,PlantGeometry%wlvg,&
         PlantGeometry%wlvd,PlantGeometry%wso,&
         PlantGeometry%wst,PlantGeometry%wrt,PlantGeometry%wrtd,&
         PlantGeometry%wcrn, sum(SucrosPlant%storage), ptype

    evapOut = flux+PlantGeometry%rain   ! evaporation from soil

    if(abs(PlantGeometry%flux-flux) > abs(PlantGeometry%flux)/100 +0.001 .and. PlantGeometry%flux<0 .and. flux<0) then
       !write(*,*) "evap",  PlantGeometry%flux-flux, PlantGeometry%flux
       evapOut = PlantGeometry%flux+PlantGeometry%rain
    end if
    Write(61,'(I4,2I3,F9.2,1P,22E13.5)') yy,mm,dd, SimTime,&
         PlantGeometry%eva_soil,&
         PlantGeometry%eva_plants,&
         PlantGeometry%eva_soil - PlantGeometry%eva_soil_intercept,&
         PlantGeometry%eva_plants - PlantGeometry%transpiration,&
         evapOut,&
         PlantGeometry%transpiration_act/(Simtime-lasttime),&
         PlantGeometry%flux, flux, alphaAvg
    Close(60)
    Close(61)
    lasttime=SimTime
2   Format(A)
3   Format(A10,A10,A11,22A13)

  End Subroutine WritePlantOutput

  Function calculate_akc(ip)
    Implicit None
    Integer, Intent(in) :: ip
    Real(dp) :: calculate_akc
    Real(dp) :: x
    x=SucrosPlant(ip)%akcMin
    calculate_akc=x+PlantGeometry%slaig/(2.88_dp*one_dp/(one_dp-x))
  End Function calculate_akc

  Function calculate_glai(ip)
    Implicit None
    Integer, Intent(in) :: ip
    Real(dp) :: calculate_glai
    Real(dp) :: x
    x=PlantData%dayOfYear
    calculate_glai=InterpolateTab(SucrosPlant(ip)%tab(8),x)
  End Function calculate_glai

End Module Plants
