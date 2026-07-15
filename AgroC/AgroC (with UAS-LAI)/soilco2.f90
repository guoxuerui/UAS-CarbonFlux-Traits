! Source file SOILCO2.FOR |||||||||||||||||||||||||||||||||||||||||||||*
!                                                                      *
!     SOILCO2 - Numerical model of one-dimensional water flow, heat    *
!               transport and carbon dioxide production and transport  *
!               in a variably saturated porous medium                  *
!                                                                      *
!     J. Simunek and D. L. Suarez, May, 1994                           *
!                                                                      *
!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

program SOILCO2

  Use datatypes
  Use Geometry
  Use Material
  Use Carbon
  Use Input
  Use Output
  Use Variables
  Use timedata
  Use Watflow
  Use SourceSink
  Use Temperature
  Use Plants, only: InitPlantData, &
    CalculateActualTranspiration, ResetActualTranspiration, &
    CalculatePlantGrowth, WritePlantOutput, ReadPlantData, ApplyRootSourceSink
  Use Nitrogen, Only: WriteNitrogenOutput, NitrogenIn, nit_sink, NITinpFromOC, NitrogenTop
  Use Phosphorus, Only: WritePhosphorusOutput, PhosphorusIn, PhosphorusTop, PinpFromOC, p_sink
  Use Solute, Only: solutes, ChemIn, Conc,SorbOut, SetChemBC, cTop, cBot
  Implicit None

  Integer :: PLevel=1,ALevel=1,TLevel=1
  logical :: ConvgF
  real(dp) :: t, dt
  Real(dp) ::dtMaxC
  logical :: outp
  logical :: lOutput(14)

  real(dp) :: COCumT=0,COCumA=0
  real(dp) :: COVolI,CumT=0,wCumT=0,cvTop,cvBot
  real(dp) :: dtOld, plantinp=0, co2Sink, wCumA=0, wVolI, Yield=0
  real(dp) :: fET=0
  Real(dp) :: cRoot=0,hRoot=0,vRoot=0,vTop  ! sink
  real(dp) :: vProd,vProdh,vProdr, molProdr, molProdh, belowground_respiration
  Integer :: ItCum=0,Iter
  Integer :: i,kk,err=0
  Integer :: timestep
  Character(len=256), Dimension(:), Allocatable :: AtmosLines
  Character(len=32) :: arg
  ! solute
  Logical :: lVapor
  Real(dp) :: tConv,xConv
  Integer :: tnew
  !CHARACTER(len=32) :: atmosphFile

  ! COMMON outputdirectory
  ! character*12 :: outputdirectory = "outputfiles/"

  lOutput=.true.
  Call get_command_Argument(1, arg)
  If(Trim(arg) == 'debug') Then
     debug=.True.
  Else
     debug=.False.
  End If
  outp=.false.
  call init_variables()

  !     Read input data
  open(30,file='selector.in',  status='old')
  call BasInf(lOutput)
  if(TopInF.or.BotInF) &
       open(31,file='atmosph.in',   status='OLD')
  call NodInf(lOutput)
  call MatIn(lOutput)
  call ReadPlantData("plants.in")
  Call InitPlantData()

  ! hh   convert wc to h if required
  If(InitWC) Call wc2h()
  call SetMat(ThOld)
  Call TmIn(t,dt)
  call InitSimTime(t,dt)
  Allocate (AtmosLines(AtmPerDay))

  If(TopInF.Or.BotInF) Then
     Call initAtmosph(lAmpl)
     call readAtmosph(AtmosLines, AtmPerDay)
     !timestep = dtplants
     timestep = 1
     call SetBC(rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L, &
          TopInF,BotInF,tTop,tBot,fET,plantinp,lAmpl,Ampl, AtmosLines(1),tnew)
  endif

  Call RootIn(lOutput)
  if(lRoot) then
     If(PlantsExist) Then 
        lRoot=.false.
     else 
        Call SetRG(t,CumT,fET)
     End If
  end if

  Call SinkIn()
  if(PlantsExist) then
     Call ApplyRootSourceSink(ThOld, Sink, hRoot,vRoot,cRoot)
  Else
     If(SinkF) Then
        Call SetSnk(hRoot,vRoot,rRoot,cRoot)
     End If
  end if

  Call TempIn(lOutput)
  Call CO2In(lOutput)
  if(lCO2) then
     call OCinit()
     Call NitrogenIn(lOutput)
     Call OCinp(plantinp)
     Call NITinpFromOC()
     Call PhosphorusIn(lOutput)
     Call Produc(tOld,dt,cumT,vProd,vProdh,vProdr, molProdr, molProdh, belowground_respiration,outp,lOutput)
  end if
  Call ChemIn()
  Call SetChemBC(tnew, Prec, rSoil, dt, hNew(NumNP))
  Call NitrogenTop(tnew, Prec, rSoil, cTop, cBot)
  Call PinpFromOC()
  Call PhosphorusTop(tnew, Prec, rSoil, cTop, cBot)

  Close(30) ! selector.in
  If(lOutput(1)) Call Profil()
  If(lOutput(6)) Then
     Call NodOut(ThOld,tInit)
     If(Transport) Call ConcOut(Conc,SorbOut,tInit)
  End If
  call SubReg(ThOld,ThOld,t-dt,dt,0,wCumA,wCumT,wVolI,COCumT, &
       COVolI,Yield, lOutput, COCumA)

  if(lWat.or.lTemp.or.lCO2) &
       call Veloc(NumNP,hOld,Con,coord,CosAlf,vOld,ThOld(NumNP), &
       ThOld(NumNP),Sink(NumNP),dt)
  if(.not.lWat) then
     do i=1,NumNP
        vNew(i)=vOld(i)
        ThNew(i)=ThOld(i)
     end do
  end if
  if(lCO2 .and. .not.lTemp) then
     TempO=TempN
  end if
  if(PlantsExist) call CalculatePlantGrowth(t)

  !     Solve water movement
  TimeLoop: Do
     If(lWat) Then
        call calculateWatFlow(t,dt,HelpP,HelpR,HelpS,Iter,ItCum,ConvgF)
     else
        iter=1
        ItCum=ItCum+1
     end if

     !     To calculate the moistures and the velocities
     if(lWat.or.lTemp.or.lCO2) &
          call Veloc (NumNP,hNew,Con,coord,CosAlf,vNew,ThNew(NumNP), &
          ThOld(NumNP),Sink(NumNP),dt)
     !     Calculation of heat transport
     if(lTemp) &
          Call Temper(t,dt,HelpP,HelpR,HelpS,HelpQ,HelpT,Disp)

     !     Calculation of nitrogen
     Call nit_sink(dt)
     ! caclulation of phosphorus
     Call p_sink(dt)

     !     Calculations of the CO2 transport
     if(lCO2) &
          Call Gas(dt,cvTop,cvBot,HelpP,HelpR,HelpS,HelpQ,HelpT,co2Sink)

     !     Calculations of the solute transport -----------------------------
     If(Transport) Then
        call Solutes(coord,dt,t,MatNum,ThOld,ThNew,vOld,vNew,Disp,&
             HelpP,HelpR,HelpS,HelpQ,dtMaxC,TempO,TempN,ths,hTemp,&
             Sink,TLevel,dtMin,dtOpt,lWat,xConv,tConv,lVapor,rBot,err,&
             Beta,AtmBC,SinkF)
        If(err.Ne.0) Call WriteError('Does not converge in the solute transport module !')
     End If
     !     T-level information
     kk=1
     If(TopInF.Or.BotInF) kk=0
     call TLInf(kk,Con(1),Con(2),Con(NumNP),Con(NumNP-1),coord(1),coord(2), &
          coord(NumNP),coord(NumNP-1),CosAlf,t,dt,Iter,TLevel,ShortF, &
          TPrint(PLevel),rTop,rRoot,vTop,vRoot,hNew(NumNP), &
          hNew(NumNP-1),hRoot,hNew(1),hNew(2),CumQ,ItCum,KodTop, &
          KodBot,ConvgF,lCO2,cRoot,CO2(NumNP),CO2(1),cvTop,cvBot, &
          vProdh,vProd,CumT,TempN(NumNP),co2Sink,ThNew(NumNP), &
          ThOld(NumNP),Sink(NumNP),T1,T2,T3,wCumT,wCumA,COCumT, &
          Yield,lRoot,lOutput,COCumA, molProdr, molProdh, belowground_respiration)

     if(NObs.gt.0.and.kk.eq.1 .and. lOutput(8)) &
          call ObsNod(t, cvTop, lOutput, belowground_respiration)

     !     P-level information
     if(abs(TPrint(PLevel)-t).lt.EpsTime) then
        If(lOutput(6)) Then
           Call NodOut(ThNew, TPrint(PLevel))
        End If
        call SubReg(ThNew,ThOld,t,dt,PLevel,wCumA,wCumT, &
             wVolI,COCumT,COVolI,Yield, lOutput,COCumA)
        If(PLevel<MaxTPrint) PLevel=PLevel+1
     end if
     ! write output for the plants.
     If(PlantsExist) Then
        Call CalculateActualTranspiration(dt,vTop)
        If(SavePlantOutput(t)) Then
           Call WritePlantOutput(t)
           Call ResetActualTranspiration()
        End If
     End If
     If( (TopInF.Or.BotInF) .And. Abs(t-tAtm)<EpsTime ) Then
        If( timestep<AtmPerDay) Then
            timestep = timestep+1
        else
           call readAtmosph(AtmosLines, AtmPerDay)
           timestep = 1
        End If
     end if
     !     A-level information
     If(Abs(t-tAtm)<EpsTime.And.(TopInF.Or.BotInF)) Then
        call ALInf(t,CumQ,hNew(NumNP),hRoot,hNew(1),ALevel,cvTop, &
             cvBot,CO2(NumNP),cRoot,CO2(1),vProdh,vProd,lCO2,lOutput)
        call SetBC(rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L, &
             TopInF,BotInF,tTop,tBot,fET,plantinp,lAmpl,Ampl, AtmosLines(timestep),tnew)
        Call SetChemBC(tnew,Prec,rSoil,dt,hNew(NumNP))
        Call NitrogenTop(tnew, Prec, rSoil, cTop, cBot)
        Call PhosphorusTop(tnew, Prec, rSoil, cTop, cBot)
        If(lCO2) Then
           Call OCinp(plantinp)
           Call NITinpFromOC()
           Call PinpFromOC()
        Endif
        if(NObs.gt.0) call ObsNod(t, cvTop, lOutput, belowground_respiration)
        call PoolOut(t, lOutput)
        If(Transport) Then
           Call MassOut(Conc,SorbOut,t)
           If(lOutput(6)) Call ConcOut(Conc,SorbOut,t)
        Endif
        If(lOutput(11)) Call ProdOut(t)
        Call WriteNitrogenOutput(t)
        Call WritePhosphorusOutput(t)
        ALevel=ALevel+1
        outp=.true.
     end if
     If(PlantsExist) Then
        Call ApplyRootSourceSink(ThNew, Sink, hRoot,vRoot,cRoot)
     End If
     !     Time governing
     if(abs(t-tMax).le.EpsTime) then
        if(lOutput(3)) write(71,'(''end'')')
        if(lOutput(4)) write(72,'(''end'')')
        if(lOutput(5)) write(73,'(''end'')')
        if(lOutput(8)) write(77,'(''end'')')
        Exit TimeLoop
     else
        dtOld=dt
        KTOld=KodTop
        KBOld=KodBot
        Call TmCont(t,dt,Iter,TPrint(PLevel),NumNP,coord,ThNew,Disp,MatNum,thS,dtMaxC,Transport)
        TLevel=TLevel+1
     end if
     if(.not.PlantsExist) then
        If(lRoot) &
             Call SetRG(t,CumT,fET)
        If(SinkF) &
             Call SetSnk(hRoot,vRoot,rRoot,cRoot)
     end if

     If(PlantsExist) Then
        ! calculate plant growth only once per day
        If(SimTime%newSimDay) then
           Call CalculatePlantGrowth(t)
        end if
     End If

     if(lCO2) &
          Call Produc(tOld,dt,cumT,vProd,vProdh,vProdr, molProdr, molProdh, belowground_respiration, outp,lOutput)
     !     New pressure heads
     if(lWat) then
        If(.Not.ConvgF .And. ExitConv) Exit TimeLoop
        do i=1,NumNP
           hHelp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
           hOld(i) =hNew(i)
           hNew(i) =hHelp(i)
           ThOld(i)=ThNew(i)
           vOld(i) =vNew(i)
        end do
     end if

  End Do TimeLoop
!  If(warnings) Print *,'Warnings in agroc.warnings...'

end program SOILCO2
