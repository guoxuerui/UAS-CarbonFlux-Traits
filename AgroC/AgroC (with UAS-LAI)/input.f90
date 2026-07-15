Module Input

  Use datatypes
  Implicit None
  Logical :: lSurf
  Real(dp) :: rSoil,Prec
  Integer :: iTemp=0,jTemp=0  ! column index of temperature in atmosph.in
Contains
  
  Subroutine BasInf(output)
    Use datatypes
    Use geometry
    Use Variables
    Implicit None
    Integer, Parameter :: NUnitD=8
    Integer :: iUnit(NUnitD)=(/ 50,70,71,72,73,75,76,77 /)

    Integer :: i
    character*72 Hed
    logical :: output(:)
    Character :: line*256,unit*4

    Read(30, '(a)') line
    i = Index(line,'Pcp_File_Version=')
    If(i==1) Then
       Read(line(18:),*) inputVersion
       Read(30,*)
    Else
       inputVersion=1
    End If
    read(30,*)
    read(30,1) Hed
    read(30,1) line
    if(line(1:15).eq.'i_check run_inf') then
       read(30,*) output
       read(30,*)
    endif

    If(wordsf(line).Eq.6) Then
       read(30,*) LUnit,TUnit,MUnit,MaxIt,TolTh,TolH
    else
       read(30,*) LUnit,TUnit,MUnit,MaxIt,TolTh,TolH,ExitConv
    Endif
    unit=to_lower(LUnit)
    If(unit=='mm') Then
       Unit_L_Input=1
    Elseif(unit=='cm') Then
       Unit_L_Input=2
    Elseif(unit=='dm') Then
       Unit_L_Input=3
    Elseif(unit=='m') Then
       Unit_L_Input=4
    Elseif(unit=='km') Then
       Unit_L_Input=5
    Else
       Call WriteError('Invalid unit for length in selector.in')
    Endif
    unit=to_lower(TUnit)
    If(unit=='hour') Then
       Unit_T_Input=1
    Elseif(unit=='day') Then
       Unit_T_Input=2
    Else
       Call WriteError('Invalid unit for time in selector.in')
    Endif
    unit=to_lower(MUnit)
    If(unit=='mg') Then
       Unit_M_Input=1
    Elseif(unit=='g') Then
       Unit_M_Input=2
    Elseif(unit=='kg') Then
       Unit_M_Input=3
    Elseif(unit=='t') Then
       Unit_M_Input=4
    Else
       Call WriteError('Invalid unit for mass in selector.in')
    Endif
    read(30,*)
    read(30,1) line
    PlantsExist=.false.
    LNitrogen=.False.
    Transport=.False.
    i=wordsf(line)
    If(i.Eq.5) Then
       read(line,*) ShortF,lWat,lCO2,lRoot,lTemp
       InitWC=.false.
    Else If(i.Eq.6) Then
       read(line,*) ShortF,lWat,lCO2,lRoot,lTemp,InitWC
    Else If(i.Eq.7) Then
       read(line,*) ShortF,lWat,lCO2,lRoot,lTemp,InitWC,PlantsExist
    Else If(i.Eq.8) Then
       Read(line,*) ShortF,lWat,lCO2,lRoot,lTemp,InitWC,PlantsExist,LNitrogen
    Else If(i.Eq.9) Then
       Read(line,*) ShortF,lWat,lCO2,lRoot,lTemp,InitWC,PlantsExist,LNitrogen,Transport
    Else If(i.Eq.10) Then
       Read(line,*) ShortF,lWat,lCO2,lRoot,lTemp,InitWC,PlantsExist,LNitrogen,lPhosphorus,Transport
    Else
       Call WriteError('Invalid number of logical switches in line 9 of selector.in')
    Endif
    If(lNitrogen) Then
       lCO2=.True.
       NSfrom=NSnit
       Transport=.True.
    Endif
    If(Lphosphorus) NSfrom=NSfrom+1 ! + phosphorus
    if(output(1)) open(50,file='i_check.out',Status='REPLACE')
    if(output(2)) open(70,file='run_inf.out',Status='REPLACE')
    if(output(3)) open(71,file='t_level.out',Status='REPLACE')
    if(output(4)) open(72,file='a_level.out',Status='REPLACE')
    if(output(5)) open(73,file='co2_inf.out',Status='REPLACE')
    If(output(6)) Then
       Open(75,file='nod_inf.out',Status='REPLACE')
       if(Transport) Open(86,file='conc.out',Status='REPLACE')
    Endif
    if(output(7)) open(76,file='balance.out',Status='REPLACE')
    if(output(8)) open(77,file='point.out',Status='REPLACE')
    if(output(9)) open(78,file='nod_pool.out',Status='REPLACE')
    if(output(10)) open(79,file='reduction.out',Status='REPLACE')
    if(output(11)) open(80,file='nod_prod.out',Status='REPLACE')
    if(output(12)) open(81,file='matlab.out',Status='REPLACE')
    if(output(13)) open(82,file='invers.out',Status='REPLACE')

    read(30,*)
    read(30,1) line
    if(wordsf(line).eq.4) then
       read(line,*) TopInF,AtmBC,SinkF,WLayer
       lAmpl=.false.
    else
       read(line,*) TopInF,AtmBC,SinkF,WLayer,lAmpl
    endif
    read(30,*)
    read(30,*) BotInF,qGWLF,FreeD,SeepF
    read(30,*)
    read(30,1) line
    if(wordsf(line).eq.2) then
       read(line,*) KodTop,KodBot
       TopBotH(1)=1e35
       TopBotH(2)=1e35
    else
       read(line,*) KodTop,KodBot,TopBotH
    Endif
    read(30,*)
    read(30,*) rTop,rRoot,rBot,hCritS,hCritA
    read(30,*)
    read(30,*) zSurf,GWL0L,Aqh,Bqh

    ! Input modifications
    rRoot=abs(rRoot)
    hCritA=-abs(hCritA)
    if(TopInF) KodTop=isign(3,KodTop)
    if(BotInF) Then
       If(KodBot.eq.6) Then
          ! bottom boundary condition for a specified pressure head.
          ! The pressure head value is the hB value read in Atmosph.in.
          ! The pressure head BC switches to a zero flux BC,
          ! if there is a flux into the domain.
          ! BotInF must be set to true.
       Else
          KodBot=isign(3,KodBot)
       Endif
    Endif
    if(AtmBC) then
       hCritS=0
       KodTop=-4
    end if
    If(WLayer) KodTop=-Abs(KodTop)
    if(qGWLF)  KodBot=-7
    if(FreeD)  KodBot=-5
    if(SeepF)  KodBot=-2
    kTOld=KodTop
    kBOld=KodBot
    if(PlantsExist .and. (lRoot.or.SinkF)) &
         Print *,'Warning: PlantsExist is activated, therefore lRoot/SinkF are disactivated.'

    write(*,*) '-----------------------------------------------------'
    write(*,*) '|                                                   |'
    write(*,*) '|                     SOILCO2                       |'
    write(*,*) '|                                                   |'
    write(*,*) '|  Code for simulating water flow, heat transport   |'
    write(*,*) '|  and carbon dioxide production and transport in   |'
    write(*,*) '|  one-dimensional variably saturated porous media  |'
    write(*,*) '|                                                   |'
    write(*,*) '|                   version 1.2                     |'
    write(*,*) '|              Last modified: May, 1994             |'
    write(*,*) '|                                                   |'
    write(*,*) '-----------------------------------------------------'
    write(*,*)
    write(*,*) Hed
    Write(*,*)
    do i=1,NUnitD
       if(output(i)) then
          write(iUnit(i),*) '******* Program SOILCO2'
          write(iUnit(i),*) '******* ',Hed
          write(iUnit(i),*) 'Units: L = ',LUnit,', T = ',TUnit, &
               ', M = ',MUnit
       endif
    end do
    if(output(1)) then
       write(50,*)
       write(50,*) 'CosAlf,MaxIt,TolTh,TolH'
       write(50,110) CosAlf,MaxIt,TolTh,TolH
       write(50,*)
       write(50,*) 'TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,FreeD,SeepF&
            &,lWat,lCO2,lRoot,lTemp,InitWC'
       write(50,120) TopInF,BotInF,AtmBC,SinkF,WLayer,qGWLF,FreeD, &
            SeepF,lWat,lCO2,lRoot,lTemp,InitWC
       If(PlantsExist) Write(50,*) 'using plants'
    endif
    if(.not.(TopInF.or.BotInF) .and. output(4)) write(72,130)

1   format(A)
110 format(f6.3,i5,f8.5,f8.3)
120 format(13l6)
130 format(////)
  End Subroutine BasInf


  !***********************************************************************
  subroutine NodInf(output)
    Use Geometry
    Use Material, only: NMat, allocate_material_data
    Use Carbon, only: CO2, CO2old,allocate_co2data
    Use Variables, Only: hBot,hTop, Transport, LNitrogen, Lphosphorus
    Use Solute, Only: NS, Conc, allocate_chem_data
    Implicit None
    Logical :: output(:)
    Integer :: i,j,n,nOld,nw
    Real(dp) :: dx,SBeta,SCO2,ShOld,STemp
    Real(dp), Allocatable :: SConc(:)
    Character :: line*120
    Character(len=16), Save :: field(5)

    write(*,*)'reading nodal information'
    read(30,*)
    read(30,*)
    Read(30,'(A)') line
    nw=read_words(line,field)
    Read(field(1),*) NumNP
    Read(field(2),*) NMat
    Read(field(3),*) NLay
    Read(field(4),*) NObs
    NS=0
    If(nw>4 .And. (LNitrogen.Or.Transport)) Read(field(5),*) NS
    If(LNitrogen .And. NS<3) Call WriteError('At least 3 solutes required for nitrogen.')
    If(Lphosphorus .And. NS<4) Call WriteError('At least 4 solutes required for phosphorus.')
    call allocate_geometry_data()
    call allocate_co2data()
    Call allocate_material_data()
    Call allocate_chem_data(NumNP)
    read(30,*)
    read(30,*) (Node(i),i=1,NObs)
    read(30,*)
    read(30,*) (Coord(i),i=NumNP-1,1,-1)
    read(30,*)

    ! Transformation from increments to coordinates,
    ! highest index is top (surface), lowest bottom
    Coord(NumNP)=zSurf
    do i=NumNP-1,1,-1
       Coord(i)=Coord(i+1)-Coord(i)
    end do

    !     Read nodal point information
    Allocate(SConc(NS))
    j=NumNP
    nOld=1
    do while(j>0)
       Read(30,*) n,hOld(n),MatNum(n),LayNum(n),Beta(n),CO2(n),TempN(n),Conc(:,n)
       If(j<n) Then
          write(*,*)'ERROR in NodInf at node =', n
          stop
       Else If(j>n) Then
          dx=Coord(nOld)-Coord(n)
          ShOld=(hOld(nOld)-hOld(n))/dx
          SBeta=(Beta(nOld)-Beta(n))/dx
          SCO2=(CO2(nOld)-CO2(n))/dx
          STemp=(TempN(nOld)-TempN(n))/dx
          If(Transport) SConc=(Conc(:,nOld)-Conc(:,n))/dx
          do i=nOld-1,n+1,-1
             dx=Coord(nOld)-Coord(i)
             hOld(i)=hOld(nOld)-ShOld*dx
             Beta(i)=Beta(nOld)-SBeta*dx
             CO2(i)=CO2(nOld)-SCO2*dx
             CO2old(i)=CO2(i)   !modification N.Prolingheuer
             TempN(i)=TempN(nOld)-STemp*dx
             If(Transport) Conc(:,i)=Conc(:,nOld)-SConc*dx
             MatNum(i)=MatNum(i+1)
             LayNum(i)=LayNum(i+1)
          end do
          j=n
       End If
       nOld=n
       j=j-1
    End Do
    Deallocate(SConc)
    coord_eps=0.001*(coord(NumNP)-coord(NumNP-1))
    SBeta=Beta(NumNP)*(coord(NumNP)-coord(NumNP-1))/2.
    do i=2,NumNP-1
       SBeta=SBeta+Beta(i)*(coord(i+1)-coord(i-1))/2.
    end do
    do i=2,NumNP
       if(SBeta.gt.0.) then
          Beta(i)=Beta(i)/SBeta
       else
          Beta(i)=0.
       end if
    end do

    ! calculate volume of nodes
    delta_z(    1)=0.5*(coord(2)-coord(1))
    delta_z(NumNP)=0.5*(coord(NumNP)-coord(NumNP-1))
    Do n=2,NumNP-1
       delta_z(n)=0.5*(coord(n+1)-coord(n-1))
    Enddo
    !     Print nodal information
    if(output(1)) write(50,110)
    Do n=NumNP,1,-1
       If(Transport) Then
          If(output(1)) Write(50,120) n,Coord(n),hOld(n),MatNum(n),LayNum(n), &
               Beta(n),CO2(n),TempN(n),Conc(:,n)
       Else
          If(output(1)) Write(50,120) n,Coord(n),hOld(n),MatNum(n),LayNum(n), &
               Beta(n),CO2(n),TempN(n)
       Endif
       hNew(n) =hOld(n)
       hHelp(n)=hOld(n)
    end do
    if(output(1)) write(50,'(''end'')')
    hBot=hNew(1)
    hTop=hNew(NumNP)

110 format (/'Nodal point information'//'node      x         ', &
         'hOld    MatN LayN  Beta      CO2       Temp     Conc...'/)
120 Format (i4,2f11.3,2i5,f8.3,e12.4,f8.3,99e12.4)

  end subroutine NodInf


  !***********************************************************************
  subroutine MatIn(output)
    Use Material
    Implicit None
    Integer, parameter :: iMax=12
    Real(dp) :: Qe
    logical :: output(:)
    Integer :: i,M,NParTemp
    Character :: line*256

    write(*,*) 'reading material information'
    read(30,*)
    read(30,*)
    Read(30,'(A)') line
    i=wordsf(line)
    NTab=100
    If(i<4) Then
       Read(line,*) hTab1,hTabN,NParTemp
    Elseif(i<5) Then
       Read(line,*) hTab1,hTabN,NParTemp,iModel
    Else
       Read(line,*) hTab1,hTabN,NParTemp,iModel,NTab
    Endif
    If(iModel==1) Then
       NPar=10
    Elseif(iModel==2) Then
       NPar=13
    Else
       Write(*,*) 'Input error: Wrong model for material parameters.'
       Stop
    Endif
    If(.Not.(iModel==1 .And. NParTemp==9)) NParTemp=NPar
    hTab1=-min(abs(hTab1),abs(hTabN))
    hTabN=-max(abs(hTab1),abs(hTabN))
    read(30,*)
    If(output(1)) Then
       Write(50,*) 'MatNum, Param. array:'
       Write(50,*)
       Write(50,150) 'Mat','Qr','Qs','Qa','Qm','Alfa','n','Ks','Kk','Qk','BPar','w2','alfa2','n2'
    Endif
    do M=1,NMat
       Par(10,M)=0.5
       read(30,*)      (Par(i,M),i=1,NPar)
       if(output(1)) write(50,120) M,(Par(i,M),i=1,NPar)
    end do

    call allocate_material_tab_data()
    Call GenMat()
    If(output(1)) Then
       Write(50,*)
       Write(50,'(A6,5A16)') 'MatNum','Qe','Q','h','C','K'
       do M=1,NMat
          Write(50,*)
          Do i=1,NTab
             Write(*,'(1p,10e14.5)') ths(M),TheTab(i,M),thr(M),ths(M)-thr(M),ths(M)-TheTab(i,M)
             Qe=(TheTab(i,M)-thr(M))/(ths(M)-thr(M))
             Write(50,140) M,Qe,TheTab(i,M),hTab(i),CapTab(i,M),ConTab(i,M)
          End Do
       End Do
    Endif
120 Format(i5,1P,50e12.3)
140 Format(i6,1P,5e16.8)
150 Format(A5,50A12)
  end subroutine MatIn


  !***********************************************************************
  subroutine GenMat()
    Use Material
    Implicit None
    Integer :: i,M
    Real(dp) :: dlh,alh,Hr

    write(*,*) 'generating materials'
    dlh=(log10(-hTabN)-log10(-hTab1))/(NTab-1)
    do i=1,NTab
       alh=log10(-hTab1)+(i-1)*dlh
       hTab(i)=-10**alh
    end do
    do M=1,NMat
       Hr       =FH(0.0_dp,Par(:,M))
       hSat(M)  =FH(1.0_dp,Par(:,M))
       ConSat(M)=FK(0.0_dp,Par(:,M))
       thr(M)   =FQ(Hr   ,Par(:,M))
       ths(M)   =FQ(0.0_dp,Par(:,M))
       do i=1,NTab
          ConTab(i,M)=FK(hTab(i),Par(:,M))
          CapTab(i,M)=FC(hTab(i),Par(:,M))
          TheTab(i,M)=FQ(hTab(i),Par(:,M))
       end do
    end do

  end subroutine GenMat


  !***********************************************************************
  ! hh   convert initial watercontent to initial pressure head
  subroutine wc2h()
    Use Geometry, only: NumNP,MatNum,hNew,hOld,hHelp
    Use Material, only: Par,FH
    Use Variables, Only: hBot,hTop,TopBotH,KodTop,KodBot
    Implicit None
    Integer :: n,m
    Real(dp) :: missing

    missing=1e30
    do n=1,NumNP
       m=MatNum(n)
       hOld(n)=FH( (hOld(n)-Par(1,m))/(Par(2,m)-Par(1,m)), Par(:,m) )
       hNew(n) =hOld(n)
       hHelp(n)=hOld(n)
    end do
    if(TopBotH(2).lt.missing .and. KodBot.gt.0) Then
       hBot=TopBotH(2)
       hNew(1)=TopBotH(2)
       hOld(1)=TopBotH(2)
       hHelp(1)=TopBotH(2)
    else
       hBot=hNew(1)
    endif
    if(TopBotH(1).lt.missing .and. KodTop.gt.0) Then
       hTop=TopBotH(1)
       hNew(NumNP)=TopBotH(1)
       hOld(NumNP)=TopBotH(1)
       hHelp(NumNP)=TopBotH(1)
    else
       hTop=hNew(NumNP)
    endif
  end subroutine wc2h


  !***********************************************************************
  Subroutine TmIn(t,dt)
    Use Variables, Only: TPrint,MaxTPrint,TopInF,BotInF,PlantsExist
    Use timedata
    Implicit None
    Integer :: i,MPL
    Real(dp) :: t,dt

    write(*,*)'reading time information'
    read(30,*)
    read(30,*)
    read(30,*) dt,dtMin,dtMax,DMul,DMul2,MPL
    read(30,*)

    MaxTPrint=MPL+1
    write(*,*) "MaxTPrint", MaxTPrint
    Allocate (TPrint(MaxTPrint))
    read(30,*) (TPrint(i),i=1,MPL)
    DtOpt=dt
    if(TopInF.or.BotInF) then
       do i=1,8
          read (31,*)
       end do
       read(31,*) tInit,lSurf
       read(31,*)
       read(31,*)
       read(31,*)
       tMax=1e20
       if(PlantsExist) lSurf=.false.
    else
       tInit=0.0
       lSurf=.false.
       tMax=TPrint(MPL)
       tAtm=tMax
    end if
    TPrint(MPL+1)=tMax
    tOld=tInit
    t=tInit+dt
    dtInit=dt
    tAtmOld=0
    read(30,'(A)') headline

  end subroutine TmIn


  !***********************************************************************
  subroutine RootIn(output)
    Use Variables
    Implicit None
    Integer :: iRFak
    Real(dp) :: rtm,tRMed,xRMed
    logical :: output(:)

    If(.Not.lRoot) Then
       Print *,'Skipping block E'
       If(headline(1:12) == '*** BLOCK E:') Call SkipBlock('E')
       Return
    Endif
    write(*,*)'reading root growth information'
    if(output(1)) write(50,110)
    read(30,*)
    read(30,*)
    read(30,*) kRoot,kBeta
    if(kRoot.eq.1) then
       read(30,*)
       read(30,*) T1,T2,T3
    end if
    read(30,*)
    read(30,*) iRFak,tRMin,tRMed,tRHarv,xRMin,xRMed,xRMax,RDDMax
    if(kRoot.eq.1) then
       RGR=-log(xRMin/99./(xRMax-xRMin))
       if(output(1)) write(50,120) tRMin,tRHarv,xRMax,RDDMax,RGR
    else
       if(iRFak.eq.0) then
          tRMed=(tRHarv-tRMin)/2.
          xRMed=(xRMax-xRMin)/2.
       end if
       rtm=tRMed-tRMin
       RGR=-(1./rtm)*log((xRMin*(xRMax-xRMed))/(xRMed*(xRMax-xRMin)))
       if(output(1)) write(50,130) tRMin,tRHarv,xRMin,xRMax,RGR
    end if
    if(kBeta.eq.1) then
       read(30,*)
       read(30,*) alpha
    end if
    read(30,'(A)') headline

110 format(//' Root growth information'/1x,23('='))
120 format(/' Degree day concept in combination with Verhulst-Pearl', &
         ' logistic growth function'//'  tRMin = ',f10.3, &
         ' tRHarv = ',f10.3,' xRMax = ',f10.3,' RDDMax = ',f10.3// &
         ' Root growth rate = ',e11.3)
130 format(/' Verhulst-Pearl logistic growth function'// &
         ' tRMin = ',f10.3,' tRHarv = ',f10.3/' xRMin = ',f10.3, &
         ' xRMax = ',f10.3/' Root growth rate = ',e11.3)

  end subroutine RootIn


  !***********************************************************************
  Subroutine SinkIn()
    Use Variables, Only: P0,H50,fR6, SinkF
    Implicit None
    If(.Not.SinkF) Then
       Print *,'Skipping block F'
       If(headline(1:12) == '*** BLOCK F:') Call SkipBlock('F')
       Return
    Endif
    write(*,*)'reading sink information'
    read(30,*)
    read(30,*) P0,H50,fR6
    read(30,'(A)') headline
  End Subroutine SinkIn

  Subroutine initAtmosph(lAmpl)
    Use Variables, Only: dailyCalculation,PlantsExist
    Logical, Intent(in) :: lAmpl
    If(PlantsExist) Then
       If(lAmpl) Then
          iTemp=14
       Else
          iTemp=13
       Endif
       If(dailyCalculation) Then
          jTemp=iTemp+1
       Else
          jTemp=iTemp
       Endif
    Endif
  End Subroutine initAtmosph
  
  !***********************************************************************
  subroutine readAtmosph(lines, timesteps)
    Use Plants, Only: ClimateData
    Use Variables, Only: dailyCalculation,PlantsExist
    Implicit None

    Integer, Intent(in) :: timesteps
    Character(len=*), Intent(out) :: lines(:)
    real(dp), dimension(15,timesteps) :: input
    Integer :: i, iend

    If(debug) Print *,'readAtmosph'
    iend = 0
    do i = 1, timesteps
       Read(31,'(A)') lines(i)
       if(lines(i)(1:3).eq.'end') then
          iend = 1
          exit
       end if
    end do
    If(PlantsExist .And. iend==0) Then
       Do i = 1, timesteps
          Read(lines(i), *) input(1:jTemp,i)
       End Do
       
       If(dailyCalculation) Then
          ClimateData%temp_avg = (input(iTemp,1)+input(iTemp+1,1))/2.0
          ClimateData%temp_min = input(iTemp,1)
          ClimateData%temp_max = input(iTemp+1,1)
       Else
          ClimateData%temp_avg = Sum(input(itemp,:))/timesteps
          ClimateData%temp_min = Minval(input(itemp,:))
          ClimateData%temp_max = Maxval(input(itemp,:))
       End If
    End If
  end subroutine readAtmosph

  !***********************************************************************
  subroutine SetBC(rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L, &
       TopInF,BotInF,tTop,tBot,fET,plantin,lAmpl,Ampl,line,tnew)
    Use Plants, Only: ClimateData, PlantData, PlantGeometry, &
         CalculateSurfaceFlux, CalculateRootExtraction, CalculateRelativeRootdist, Farquhar_used
    Use timedata, Only: tAtm,tAtmOld,tMax,lMinStep
    Use Variables, Only: dailyCalculation, RelHum, PlantsExist, Transport, KodTop
    Use Solute, Only: cTop, cTopOld, cBot
    Implicit None

    Integer, Intent(out) :: tnew
    Logical :: TopInF,BotInF,lAmpl
    Real(dp) :: rTop,rRoot,rBot,hCritA,hBot,hTop,GWL0L,tTop,tBot,fET,plantin,Ampl
    ! local variables
    character(256) :: line
    Character(len=32), Save :: field(100)
    Real(dp) :: hb,hCA,hT,rB,rR
    Integer :: i,nw
    Real(dp) :: rTopOld

    If(line(1:3).Eq.'end') Then
       Goto 1
    End If
    nw=read_words(line,field)
    Read(field(1),*) tAtm
    If(tAtm<tAtmOld .Or. equal(tAtm,tAtmOld)) Then
       Write(*,'(A,F10.1,A,F10.1)') 'Error in atmosph.in. tAtm=',tAtm,' <=',tAtmOld
       Call WriteError('Execution stopped because of input error for time data.')
    Elseif(Abs(tAtm-tAtmOld)>one_dp+eps) Then
       Write(*,'(A,F10.1)') 'Warning in atmosph.in. tAtm-tAtmOld>1 at time:',tAtm
    Endif
    tAtmOld=tAtm
    tnew=Nint(tatm)
    Read(field(2),*) Prec
    Read(field(3),*) rSoil
    Read(field(4),*) rR
    Read(field(5),*) hCA
    Read(field(6),*) rB
    Read(field(7),*) hB
    Read(field(8),*) hT
    Read(field(9),*) tTop
    Read(field(10),*) tBot
    Read(field(11),*) plantin
    If(lAmpl) Then
       Read(field(12),*) Ampl
       i=13
    Else
       i=12
    Endif
    If(PlantsExist) Then
       if(i>nw) call WriteError('Missing column et0 in atmosph.in.')
       Read(field(i),*) ClimateData%et0
       i=i+1
       If(dailyCalculation) Then
          If(i>nw) Call WriteError('Missing column temp_min in atmosph.in.')
          Read(field(i),*) ClimateData%temp_min
          i=i+1
          If(i>nw) Call WriteError('Missing column temp_max in atmosph.in.')
          Read(field(i),*) ClimateData%temp_max
          i=i+1
       Else
          If(i>nw) Call WriteError('Missing column temp in atmosph.in.')
         Read(field(i),*) ClimateData%temp
          i=i+1
       Endif
       if(i>nw) call WriteError('Missing column rad in atmosph.in.')
       Read(field(i),*) ClimateData%global_rad
       i=i+1
       ClimateData%time=tAtm
       ClimateData%prec=Prec*PlantData%UnitFactor
    Endif
    If(Farquhar_used) Then
       if(i>nw) call WriteError('Missing column RelHum in atmosph.in.')
       Read(field(i),*) RelHum
       i=i+1
       RelHum=RelHum/100_dp
    Endif
    If(Transport) Then
       cTopOld=cTop
       cTop=0
       cBot=0
    Endif

    If(lSurf) Then
       rSoil=rSoil+(1-fET)*rR
       rR=fET*rR
    end if

    !     Top of the profile
    if(PlantsExist) Then
       Call CalculateSurfaceFlux()
       rTop=PlantGeometry%flux / PlantData%UnitFactor
       hTop=hT
       rRoot=abs(PlantGeometry%transpiration) / PlantData%UnitFactor
       hCritA=-abs(hCA)
       Call CalculateRelativeRootdist()
       Call CalculateRootExtraction()
       rSoil=PlantGeometry%eva_soil_intercept / PlantData%UnitFactor
       Prec=PlantGeometry%rain / PlantData%UnitFactor
    else if(TopInF) then
       rTopOld=rTop
       hCritA=-abs(hCA)
       rTop=abs(rSoil)-abs(Prec)
       If(Abs(rTopOld-rTop).Gt.Abs(rTop)*0.2.And.rTop.Lt.0.) lMinStep=.True.
       If(KodTop.Eq.3) Then
          If(Abs(hTop-hT).Gt.Abs(hTop)*0.2) lMinStep=.True.
          hTop=hT
       End If
       rRoot=abs(rR)
    end if

    !     Bottom of the profile
    if(BotInF) then
       If(Abs(rBot-rB).Gt.Abs(rBot)*0.2) lMinStep=.True.
       rBot=rB
       If(Abs(hBot-hB-GWL0L).Gt.Abs(hBot)*0.2) lMinStep=.True.
       hBot=hB+GWL0L
    end if

    return
    ! EOF
1   tMax=tAtm
  end subroutine SetBC


  !***********************************************************************
  subroutine TempIn(output)
    Use Variables
    Use Material, Only: tempParam
    Implicit None
    logical :: output(:)
    Integer :: i
    Real(dp) :: tampl,tT,tB

    If(.Not.lTemp) Then
       Print *,'Skipping block G'
       If(headline(1:12) == '*** BLOCK G:') Call SkipBlock('G')
       Return
    Endif
    write(*,*) 'reading heat transport information'
    read(30,*)
    read(30,*) tampl,tPeriod
    if(.not.lAmpl) Ampl=tampl
    read(30,*)
    read(30,*) kTopT,tT,kBotT,tB
    if(.not.TopInf) tTop=tT
    if(.not.BotInf) tBot=tB
    read(30,*)
    if(output(1)) write(50,110)  &
         ' Heat transport information', &
         ' ==========================', &
         ' ample = ',Ampl,'  period=',tPeriod,'day/tUnits', &
         'Beta','Qn','Qo','B1','B2','B3','Cn','Co','Cw'
110 format(//,A,/,A,//,A,f10.3,A,f8.5,A,//,3A7,6A11)
    Do i=1,Ubound(tempParam,2)
       read(30,*) tempParam(:,i)
       if(output(1)) write(50,120) tempParam(:,i)
    end do
    read(30,'(A)') headline

120 format(3f7.3,1P,6e11.3)

  end subroutine TempIn


  !***********************************************************************
  subroutine CO2In(output)
    Use Geometry, only: NumNP
    Use Material, only: NMat
    Use timedata, Only: co2dt,co2dtMin
    Use Carbon
    Use Variables, Only: lCO2
    Implicit None
    logical output(:)
    character line*256
    integer i,j,nfrom,nto,nby

    If(.Not.lCO2) Then
       Print *,'Skipping block H'
       If(headline(1:12) == '*** BLOCK H:') Call SkipBlock('H')
       Return
    Endif
    write(*,*)'reading CO2 transport and production information'
    if(output(1)) write(50,110)
    read(30,*)
    Read(30,1) line
    i=wordsf(line)
    if(i.eq.1) then
       read(line,*) lStagn
       iGasDiff=1
    else if(i.eq.2) then
       read(line,*) lStagn,iGasDiff
    Else
       Read(line,*) lStagn,iGasDiff,co2dt,co2dtMin
    Endif
    if(output(1)) write(50,*) 'gas diffusion model =',iGasDiff
    if(output(1)) write(50,*)
    read(30,*)
    read(30,*) kTopCO,CO2Top,kBotCO,CO2Bot
    if(output(1)) write(50,120) kTopCO,CO2Top,kBotCO,CO2Bot
    read(30,*)
    if(output(1)) write(50,130)
    do i=1,NMat
       if(iGasDiff.eq.3 .or. iGasDiff.eq.5) then
          read(30,*) (co2parm(3+j,i),j=1,6)
          if(output(1)) write(50,140) (co2parm(3+j,i),j=1,6)
       else 
          read(30,*) (co2parm(3+j,i),j=1,3) 
          if(output(1)) write(50,140) (co2parm(3+j,i),j=1,3)
       end if
    end do
    read(30,*)
    read(30,*) gamS0,gamR0,PDDMax,kProd
    if(output(1)) write(50,150) gamS0,gamR0,PDDMax,kProd
    read(30,*)
    if(kProd.ne.0) then
       read(30,*) xR
       if(output(1)) write(50,160) xR
    else
       read(30,*) alf
       if(output(1)) write(50,170) alf
    end if
    read(30,*)
    read(30,1) line
    if(wordsf(line).eq.7) then
       read(line,*) B1,B2,cM1,cM2,HB1,HB2,fS6
       Patm=7.56387e12
       gasconst=6.20667e14
       MolCO2=0.0440098
       print *, 'Using standard units for constants'
    else
       read(line,*) B1,B2,cM1,cM2,HB1,HB2,fS6,Patm,gasconst,MolCO2
    endif
    if(output(1)) then
       write(50,180) B1,B2,cM1,cM2,HB1,HB2,fS6
       write(50,'(3A12)') 'Patm','gasconst','MolCO2'
       write(50,'(1P,3E12.5)') Patm,gasconst,MolCO2
    endif
    if(kProd.eq.3) then
       ! input for the RothC parameters
       read(30,*)
       read(30,*) co2rdpm,co2rrpm,co2rbio,co2rhum, co2pdepth,co2alf
       if(output(1)) write(50,190) co2rdpm,co2rrpm,co2rbio,co2rhum, &
            co2pdepth,co2alf

       read(30,*)
       read(30,1) line
       if(wordsf(line).eq.5) then
          read(line,*) iredW,iredT,iredCO2,iReftemp,Reftemp
          hCentury=-1000.0
          B1DPM=B1
          B1RPM=B1
          B1Bio=B1
          B1Hum=B1
       else if(wordsf(line).eq.6) then
          read(line,*) iredW,iredT,iredCO2,iReftemp,Reftemp,hCentury
          B1DPM=B1
          B1RPM=B1
          B1Bio=B1
          B1Hum=B1
       else
          read(line,*) iredW,iredT,iredCO2,iReftemp,Reftemp,hCentury, &
               B1DPM, B1RPM, B1Bio, B1Hum   
       endif
       if(output(1)) write(50,'(A,I3,A,I3/A,I3/A,I3)') &
            'Reduction function: Water:',iredW,'  Temperature:',iredT, &
            'Reference temperature model number:',iReftemp, &
            'Reduction function: CO2:',iredCO2
       If(iReftemp.eq.0 .and. output(1)) &
            Write(50,'(A,F6.2)') 'Reference temperature:',Reftemp
       if(iredW.eq.10) then
          ! read moisture intervals
          read(30,*)
          read(30,*) NumWCInt
          if(NumWCInt.lt.1) then
             print *,'Invalid number of intervals for moisture'
             stop 
          endif
          Allocate(wcint(NumWCInt-1))
          Allocate(wcfact(NumWCInt))
          read(30,*)
          read(30,*) wcint
          read(30,*)
          read(30,*) wcfact
       else if(iredW.eq.11) then
          ! read wc parameters for stepwise linear finction
          read(30,*)
          read(30,*) WCPar
       endif
       if(iredT.eq.10) then
          ! read temperature intervals
          read(30,*)
          read(30,*) NumTempInt
          if(NumTempInt.lt.1) then
             print *,'Invalid number of intervals for temperature'
             stop 
          endif
          Allocate(tempint(NumTempInt-1))
          Allocate(tempfact(NumTempInt))
          read(30,*)
          read(30,*) tempint
          read(30,*)
          read(30,*) tempfact
       else if(iredT.eq.11 .or. iredT.eq.12 .or. iredT.eq.13 ) then
          ! read parameter for temperature function
          read(30,*)
          read(30,*) TempPar
       endif
       
       read(30,*)
       if(output(1)) write(50,'(A11)') 'Clay'
       read(30,*) (co2parm(1,i),i=1,NMat)
       if(output(1)) write(50,'(F8.5)') (co2parm(1,i),i=1,NMat)
       read(30,*)
       if(output(1)) write(50,'(5A11)') 'DPM','RPM','BIO','HUM','IOM'
       if(co2alf.gt.null_dp) then
          nfrom=1
          nto=1
          nby=1
       Else If(equal(co2alf,-1.0_dp)) Then
          nfrom=NumNP
          nto=1
          nby=-1
       else 
          nfrom=1
          nto=NMat
          nby=1
       endif
       do i=nfrom,nto,nby
          read(30,*) (co2init(j,i),j=1,5)
          if(output(1)) write(50,200) (co2init(j,i),j=1,5)
       Enddo
190    format (/' kdpm = ',e11.3,' krpm = ',e11.3/' kbio = ',e11.3, &
            ' khum = ',e11.3,' Pdepth =',e11.3,' alf =',e11.3)
200    format(6e11.3)
    Endif
    read(30,'(A)') headline

1   format(A)
110 format (//' CO2 transport information'/1X,25('=')/)
120 format ('kTopCO       CO2Top       kBotCO     CO2Bot'/ &
         i4,6x,e10.3,6x,i4,4x,e10.3)
130 format (/'      Da         Dw       l      eps     campb     Xm')
140 Format (1P,2e11.3,3f8.3,e11.3)
150 format (//' Production of CO2 information'/1x,29('=')// &
         '    gamS0       gamR0       PDDMax     kProd'/2e12.3,f10.3,i10)
160 format (/' xR  = ',f10.3)
170 format (/' alf = ',f10.3)
180 format (/' B1  = ',f15.3,'  B2  = ',f15.3/' cM1 = ',f15.3, &
         '  cM2 = ',f15.3/' HB2 = ',F15.3,'  HB2 = ',F15.1/ &
         ' fS6 = ',f15.3)

  end subroutine CO2In


  !***********************************************************************
  subroutine Profil()
    Use datatypes
    Use Geometry
    Use Material, Only: thR,thS,ConSat,hSat
    implicit None
    Integer :: i, M
    Real(dp) :: ConSatN

    write(*,*)'printing profile information'
    write(50,110)
    ConSatN=ConSat(MatNum(NumNP))
    do i=NumNP,1,-1
       M=MatNum(i)
       write(50,120) i,zSurf-coord(i),thR(M),thS(M),hSat(M),ConSat(M), &
            ConSat(M)/ConSatN,Beta(i)
    end do
    write(50,'(''end'')')
    Close(50)
110 format(///' Profile information'/1x,20('=')// &
         '    n      depth     THr       THs       hs       Ks', &
         '        Ks/KsTop    Beta'/)
120 format(i5,f10.2,2f10.3,f10.1,e12.3,2f10.3)
    return
  end subroutine Profil



End Module Input
