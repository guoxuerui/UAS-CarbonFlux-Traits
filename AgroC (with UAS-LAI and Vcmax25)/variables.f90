Module Variables

  Use datatypes
  Implicit None
  ! Basic information
  Integer :: inputVersion   ! version number of the format of the input data
  Integer :: MaxIt,KodTop,KodBot,KTOld,KBOld
  Real(dp) :: TolTh,TolH,CosAlf,rTop,rRoot,rBot,hCritS,hCritA,GWL0L,Aqh,Bqh,TopBotH(2)
  Logical :: ExitConv,TopInF,BotInF,ShortF,lWat,lCO2,lRoot,lTemp,SinkF
  Logical :: WLayer,qGWLF,FreeD,SeepF,AtmBC,lAmpl,InitWC
  ! nodal information
  Real(dp) :: hTop,hBot
  ! sink input
  Real(dp) :: P0,H50,fR6
  ! root input
  Integer :: kRoot,kBeta
  Real(dp) :: T1,T2,T3,tRMin,tRHarv,xRMin,xRMax,RGR,RDDMax,alpha
  ! temperature input
  Integer :: kTopT,kBotT
  Real(dp) :: Ampl,tPeriod,tTop,tBot
  ! printing
  Integer :: MaxTPrint   ! number of t-level output
  Real(dp), Allocatable :: TPrint(:)   ! dimension MaxTPrint

  Real(dp) :: root_respiration
  Real(dp) :: heterotrophic_respiration
  Real(dp) :: aboveground_respiration
  Real(dp) :: GPP
  Real(dp) :: TER
  Real(dp) :: NPP
  Real(dp) :: NEE
  Real(dp) :: RelHum
  Real(dp) :: fluorescence_755nm

  Real(dp) :: rorad      ! [L] mean root radius
  Real(dp) :: UnitFactorPlants ! factor to convert from plants (kg) to agroc
  Real(dp) :: w0_dens    ! L/L3 root density at top

  Logical :: co2_fluxes, respiration, maint_growth   ! write output files?
  Integer :: waterstress                             ! calculate gphot with respect to waterstress?
  Logical :: rootexu                                 ! calculate rootexudation?
  Logical :: dailyCalculation=.True.                 ! only one caluclation per day? (=> atmosph.in has tMin and tMax)
  Logical :: rootdeath
  Logical :: Transport    ! solute transport
  Logical :: PlantsExist  ! plants
  Logical :: LNitrogen, LNitrogenReady=.False., LNPools
  Logical :: lPhosphorus=.False.
  Integer :: P_uptake_method ! 1=Wave 2=Apex_Epic
  !Real(dp) :: rlv, rst, rso, rrt, rcrn
  Real(dp) :: anlv=0, anst=0, anrt=0, anso=0, ancrn=0, aplv=0, apst=0, aprt=0, apso=0, apcrn=0, ap_tot=0
  Real(dp) :: alphaAvg=1.0
  Real(dp) :: Nit_ancrt=0, Pho_ancrt=0

Contains

  Subroutine init_variables()
    CosAlf=1.0
    fR6=1.0
  End Subroutine init_variables

  Subroutine Pause(text)
    Character(len=*) :: text
    Print *,text
    call sleep(1)
  End Subroutine Pause

  !Subroutine reset_Variables()
    !rlv = 0.
    !rst = 0.
    !rso = 0.
    !rcrn = 0.
    !rrt = 0.
  !End Subroutine reset_Variables

End Module Variables
