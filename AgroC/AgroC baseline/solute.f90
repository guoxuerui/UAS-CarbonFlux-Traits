Module Solute

  Use Datatypes
  Implicit None
  Integer, Parameter :: maxChPar=20
  Integer :: NS=0     ! number of solutes
  Integer :: iDualPor, iMoistDep, IterC, iTort, kBotCh, kTopCh, iBact, iNonEqul
  Real(dp) :: cAtm, CMin, Courant, cTolA, cTolR, dSurf, epsi
  Real(dp), Allocatable :: cBot(:)     ! dim NS
  Real(dp), Allocatable :: cBotOld(:)     ! dim NS
  Real(dp), Allocatable :: cTop(:)     ! dim NS
  Real(dp), Allocatable :: cTopOld(:)     ! dim NS
  Real(dp), Allocatable :: cvch0(:),cvch1(:),cvchR(:),cvchIm(:) ! dim NS
  Real(dp), Allocatable :: ChPar(:,:,:)  ! dim maxChPar x NMat x NS
  Real(dp), Allocatable :: Conc(:,:)   ! dim NS x NumNP
  Real(dp), Allocatable :: DMoist(:,:,:,:)   ! dim NMat x NS x 13 x 6
  Real(dp), Allocatable :: CNew(:)     ! dim NumNP
  Real(dp), Allocatable :: cPrevO(:)   ! dim NumNP
  Real(dp), Allocatable :: cTemp(:)    ! dim NumNP
  Real(dp), Allocatable, Dimension(:) :: g0,g1 ! dim NumNP
  Logical, Allocatable, Dimension(:) :: lLinear ! dim NS
  Logical, Allocatable, Dimension(:) :: freundlich ! dim NS
  Logical, Allocatable, Dimension(:) :: langmuir ! dim NS
  Logical, Allocatable, Dimension(:) :: lMobIm  ! dim NMAT
  Logical :: lTort,lUpW,lArtD,lBact,lEqInit,lActRSU,lFlux,lEquil,lDensity,lMassIni,lfiltr,ldualneq
  Integer :: MaxItC,iConcType
  Integer, Parameter :: maxCumCh=10  ! first dimension of CumCh
  Real(dp) :: OmegaS,OmegaW,Peclet,PeCr,rKM,TPulse
  Real(dp) :: SPot=0
  Real(dp), Allocatable, Dimension(:) :: q0,q1,Retard,SinkIm,SorbN,SorbN2,STrans
  Real(dp), Allocatable, Dimension(:,:) :: Sorb,Sorb2,SorbOut,CumCh
  Real(dp), Allocatable, Dimension(:,:,:) :: WDep
  Real(dp), Allocatable, Dimension(:,:) :: TDep
  Real(dp), Allocatable, Dimension(:) :: ThOldIm,ThNewIm,sSink
  Real(dp), Allocatable, Dimension(:) :: chem_vTop,chem_vBot,wc
  ! SubReg
  Real(dp), Allocatable, Dimension(:) :: cTot,ConVol,cVolI ! dim NS
  Real(dp), Allocatable, Dimension(:) :: cTotIm,ConVolIm,ConVolIm2,cCumT,cCumA, cGWL ! dim NS
  Real(dp), Allocatable, Dimension(:,:) :: ConSub,cMean ! dim NS x NLay
  Character(len=8), Allocatable, Dimension(:) :: solute_name
  ! top and bottom atmosheric boundary
  Integer :: ATM_dim, ATM_ind
  Integer(dp), Allocatable :: ATM_Time(:)
  Real(dp), Allocatable :: ATM_Top(:,:), ATM_Bot(:,:)

Contains

  ! allocate data that needs to be allocated before ChemIn
  Subroutine allocate_chem_data(NumNP)
    Integer, Intent(in) :: NumNP
    Allocate(solute_name(NS))
    Allocate(Conc(NS,NumNP))
    Allocate(Sorb(NS,NumNP))
    Allocate(cTop(NS))
    Allocate(cBot(NS))
    Allocate(cTopOld(NS))
    Allocate(cBotOld(NS))
    Conc=0
    cTop=0
    cBot=0
    cTopOld=0
    cBotOld=0
  End Subroutine allocate_chem_data

  
  ! Read information about solute transport
  Subroutine ChemIn()
    Use Variables, Only: inputVersion, Transport, AtmBC
    Use Geometry, Only: NumNP, NLay, NSfrom, cRootMax
    Use Material, Only: NMat, Par, FQ, iModel
    Implicit None
    Logical :: lTDep, lMoistDep=.False., lEqInit, lVar=.False.
    Integer :: j,jj,M,nPar2,nw,ns2
    Real(dp), Allocatable :: cTopTemp(:),cBotTemp(:)
    Character :: line*4096
    Character(len=32) :: field(maxChPar)
    
!     New Options - Supported by GUI
!     lBact   - Virus transport, ka,kd concept
!     lFiltr  - Filtration theory
!     lVapor  - Vapor flow
!     lEnter  - End the run with pushing Enter key
!     nPrStep - Print to the screen and T_Level files at each nPrStep
!     lPrintD - Print at a daily (given) interval
!     tPrintInt- Print interval
!     nTabMod - Kode for the input of the soil hydraulic properties tables
!     iDualPor- Dual porosity model; 
!               = 1: transfer proportional to difference in water contents
!               = 2: transfer proportional to difference in pressure heads
!     lDualNEq- Both physical and chemical nonequilibrium are considered simultaneously, two-site sorption in the mobile zone
!               (fraction of equilibrium sites in mobile zone - ChPar(13), Rate of kinetic sorption - ChPar(16)
!     lMeteo  - Meterorological input to calculate ET
!     lLAI    - Distribution of pET based on LAI
!     lDayVar - Daily variations in root water uptake and evaporation
!     lSinPrec- Sinusoidal distribution of precipitation
!     lEnBal  - Evaporation and heat flux is calculated from energy balance
!     iSunSh  - =0 Sunshine hours; =1 Cloudeness; =2 Transmission coeff.
!     iRelHum - =0 Relative Humidity; =1 Vapor Pressure
!     lMetDaily - Daily variations of meteorological variables are generated from daily average, max, and min data
!     lSnow   - Snow accumulation at the soil surface
!     SnowMF  - Amount of snow melted per 1 degree [cm3/cm2/K/d]
!     SnowLayer- Thickness of the snow layer
!     iTort   - tortuosity in solute transport, =0 for Millington and =1 for Moldrup
!     lSeep   - Seepage face initiated by different pressure head
!     hSeep   - seepage face with a different bottom pressure
!     lEqInit - initial noequilibrium phase is in equilibrium with liquid phase
!     lMassIni- initial condition is given in the total concentration [M_solute/M_soil]
!     iMoistDep- reaction rates are dependent on the water content (iMoistDep=2)
!     OmegaC  - Compensated root water uptake
!     lFlux   - Print fluxes in observation nodes instead of temperatures
!     lActRSU - Active root solute uptake
!     OmegaS  - Compensated root solute active uptake
!     SPot    - Potential root solute uptake
!     lOmegaW - Reduction of the potential root solute uptake due to reduction of the root water uptake

    If(.Not.Transport) Then
       Print *,'Skipping block K'
       If(headline(1:12) == '*** BLOCK K:') Call SkipBlock('K')
       Return
    Endif
    iModel=1 ! see material
    iConcType=0
    lBact=.False.
    lEqInit=.False.
    lActRSU=.False.
    lFlux=.False.
    lEquil=.False.
    lDensity=.False.
    lTort=.False.
    lBact=.False.
    lFiltr=.False.
    lMassIni=.False.
    ldualneq=.False.
    iDualPor=0
    iMoistDep=0
    IterC=0
    iTort=0
    iBact=0
    iNonEqul=0
    Write(*,*) 'reading solute transport information'
    Read(30,*,err=901)
    If(inputVersion<3) Then
       Read(30,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,PeCr,lTort
    Else If(inputVersion<5) Then
       Read(30,*,err=901) epsi,lUpW,lArtD,lTDep,cTolA,cTolR,MaxItC,PeCr,lTort,iBact,lFiltr
       If(iBact.Eq.1) lBact=.True.
    Else
       Call WriteError('Invalid input version in file: selector.in')
    End If
    If(inputVersion==4) Then
       Read(30,*,err=901)
       Read(30,*,err=901) iNonEqul,lMoistDep,lDualNEq,lMassIni,lEqInit,lVar
       If(lMoistDep) iMoistDep=1
       If(lVar) iTort=1
    End If
    ! allocate arrays
    Allocate(cRootMax(NS))
    Allocate(cvch0(NS))
    Allocate(cvch1(NS))
    Allocate(cvchR(NS))
    Allocate(cvchIm(NS))
    Allocate(cTot(NS),ConVol(NS),cVolI(NS))
    Allocate(cTotIm(NS),ConVolIm(NS),ConVolIm2(NS),cCumT(NS),cCumA(NS),cGWL(NS))
    Allocate(ConSub(NS,NLay))
    Allocate(cMean(NS,NLay))
    Allocate(ChPar(maxChPar,NMat,NS))
    Allocate(DMoist(NMat,NS,13,6))
    Allocate(cNew(NumNP))
    Allocate(cTemp(NumNP))
    Allocate(cPrevO(NumNP))
    Allocate(g0(NumNP))
    Allocate(g1(NumNP))
    Allocate(q0(NumNP))
    Allocate(q1(NumNP))
    Allocate(wc(NumNP))
    wc=0
    Sorb=0
    Allocate(Sorb2(NS,NumNP))
    Allocate(SorbN(NumNP))
    Allocate(SorbN2(NumNP))
    Allocate(SorbOut(NS,NumNP))
    SorbOut=0
    Allocate(STrans(NumNP))
    Allocate(Retard(NumNP))
    Retard=0
    Allocate(SinkIm(NumNP))
    Allocate(ThOldIm(NumNP))
    Allocate(ThNewIm(NumNP))
    Allocate(sSink(NumNP))
    sSink=0
    Allocate(lLinear(NS))
    lLinear=.True.
    Allocate(freundlich(NS))
    freundlich=.False.
    Allocate(langmuir(NS))
    langmuir=.False.
    Allocate(chem_vBot(NS))
    Allocate(chem_vTop(NS))
    Allocate(lMobIm(NMAT))
    Allocate(TDep(maxChPar,NS))
    Allocate(WDep(9,NMAT+2,NS))
    Allocate(CumCh(maxCumCh,NS))
    cTemp=0
    cCumT=0.0
    cCumA=0.0
    PeCr=Max(PeCr,0.1)
    Write(50,*)
    Write(50,1) 'Solute transport input'
    Write(50,1) '======================'
    if(lUpW) then
       write(50,120,err=902)
    else
       write(50,130,err=902)
       if(lArtD) write(50,140,err=902) PeCr
    end if
    write(50,150,err=902) lTDep,iMoistDep,cTolA,cTolR,MaxItC
    Read(30,*,err=901)
    ! init data
    lEquil=.True.
    ChPar=0
    Do M=1,NMat
       Read(30,*,err=901) (ChPar(j,M,1),j=1,4)
       Write(50,160,err=901) M,(ChPar(j,M,1),j=1,4)
       If(ChPar(3,M,1).Lt.1.Or.ChPar(4,M,1).Gt.0..Or.lBact) lEquil=.False.
       lMobIm(M)=.false.
       If(.Not.lBact.And.ChPar(4,M,1).Gt.0.) lMobIm(M)=.True.
       If(.Not.lEquil.And.ChPar(1,M,1).Eq.0.) Goto 903
    End Do
    Do jj=2,NS
       ChPar(1:4,:,jj)=ChPar(1:4,:,1)
    Enddo
!    If(lNitrogen) Then
!       ! set parameters for nitrogen
!       ! same values for all materials?
!       ChPar( 5,:,1)=0.18  ! DifW
!       ChPar( 7,:,1)=0.001 ! Ks
!       ChPar( 9,:,1)=1     ! Beta
!       ChPar(14,:,1)=0.005 ! SnkL1'
!       ChPar(15,:,1)=0.005 ! SnkS1'
!       j1=4
!    Else
!       j1=1
!    Endif
    Do jj=1,NS
       write(50,170,err=902) jj
       read(30,*,err=901) solute_name(jj)
       read(30,*,err=901)
       Read(30,*,err=904) (ChPar(j,1,jj),j=5,6)
       Write(50,180,err=902) (ChPar(j,1,jj),j=5,6)
       read(30,*,err=901)
       Do M=1,NMat
          If(m>1) ChPar(5:6,M,jj)=ChPar(5:6,1,jj)
          Read(30,'(A)',err=904) line
          nw=read_words(line,field)
          nw=Min(maxChPar, nw)+6
          Read(line,*,err=904) (ChPar(j,M,jj),j=7,nw)
          Write(50,190,err=902) M,(ChPar(j,M,jj),j=7,maxChPar)
          If(ChPar(8,M,jj).Gt.eps) Then
             write(50,200,err=902) M
          Else If(ChPar(9,M,jj).Gt.1+eps) Then
             write(50,210,err=902) M
          else
             write(50,220,err=902) M
          end if
          If(.Not.lEquil) Then
             if(lMobIm(M)) then
                write(50,222,err=902)
             else
                write(50,224,err=902)
             end if
          End If
          If(ChPar(9,M,jj)>1+eps) freundlich(jj)=.True.
          If(ChPar(8,M,jj)>eps) langmuir(jj)=.True.
          If(langmuir(jj) .Or. freundlich(jj)) lLinear(jj)=.False.
          If(lBact.And.(ChPar(18,M,jj).Gt.eps.Or.ChPar(15,M,jj).Gt.eps)) lLinear(jj)=.False.
       End Do
    End Do
    TDep=0.
    WDep(:,1,:)=1.
    WDep(:,2,:)=0.
    CumCh=0.0
    if(lTDep) then
       Read(30,*,err=901)
       Do jj=1,NS
          read(30,*,err=901)
          Read(30,*,err=901) (TDep(j,jj),j=5,6)
          read(30,*,err=901)
          Read(30,*,err=901) (TDep(j,jj),j=7,20)
       end do
    end if
    if(iMoistDep.eq.1) then
       Read(30,*,err=901)
       do jj=1,NS
          read(30,*,err=901)
          read(30,*,err=901) nPar2
          read(30,*,err=901)
          Read(30,*,err=901) WDep(:,1,jj)
          Read(30,*,err=901) WDep(:,2,jj)
          Do M=1,NMat
             do j=1,9
                ! WDep(j,2+M,jj)=FQ(iModel,WDep(j,2,jj),Par(1,M))
                WDep(j,2+M,jj)=FQ(WDep(j,2,jj),Par(:,M))
             end do
          end do
       end do
    end if

    Allocate(cTopTemp(NS))
    Allocate(cBotTemp(NS))
    read(30,*,err=901)
    Read(30,*,err=901) kTopCh,(cTopTemp(jj),jj=1,NS)
    read(30,*,err=901)
    Read(30,*,err=901) kBotCh,(cBotTemp(jj),jj=1,NS)
    cTopOld=cTopTemp(1:NS)
    cBotOld=cBotTemp(1:NS)
    Deallocate(cTopTemp)
    Deallocate(cBotTemp)
    If(.Not.AtmBC) Then
       ! don't overwrite data from atmosph.in
       cTop=cTopOld
       cBot=cBotOld
    Endif
    if(kTopCh.eq.-2) then 
       read(30,*,err=901)
       read(30,*,err=901) dSurf,cAtm
    end if
    write(50,230,err=902) kTopCh,(cTopOld(jj),jj=1,NS)
    Write(50,240,err=902) kBotCh,(cBotOld(jj),jj=1,NS)
    read(30,*,err=901)
    read(30,*,err=901) cRootMax
    Write(50,'(A,1P,99E12.5)',err=902) ' cRootMax = ',cRootMax
    Read(30,*,err=901)
    Read(30,*,err=901) tPulse
    Write(50,250,err=902) tPulse
    ! top and bottom boundaries
    ns2=NS-NSfrom
    Read(30,1) line
    ATM_ind=1 ! index of next time for boundary
    If(line(1:1)=='*') Then
       ATM_dim=0
       Allocate(ATM_Time(ATM_dim))
       Allocate(ATM_Top(ATM_dim,ns2))
       Allocate(ATM_Bot(ATM_dim,ns2))
       headline=Trim(line)
    Else
       Read(30,*,err=901)
       Read(30,1) line
       ATM_dim=wordsf(line)
       If(ATM_dim<1)  Call WriteError('Invalid number of solute boundaries.')
       Allocate(ATM_Time(ATM_dim))
       Allocate(ATM_Top(ATM_dim,ns2))
       Allocate(ATM_Bot(ATM_dim,ns2))
       Read(line,*,err=901) ATM_Time
       Do j=1,ns2
          Read(30,*,err=901)
          Read(30,*,err=901) ATM_Top(:,j)
          Read(30,*,err=901)
          Read(30,*,err=901) ATM_Bot(:,j)
       Enddo
       Read(30,1) headline
    Endif
    ! end of input
    Return ! normal end of subroutine
901 Call WriteError('Error when reading from input file selector.in')
902 Call WriteError('Error when writing into an output file')
903 Call WriteError('Bulk Density is equal to zero')
904 Call WriteError('Missing solute input in selector.in')

1   Format(A)
120 format(/' Upstream weighting finite-element method')
130 format(/' Galerkin finite-element method')
140 format (/' Artificial dispersion is added when Peclet number is higher than',f10.3)
150 format(//' lTDep     lWDep     cTolA     cTolR   MaxItC'&
         &        /l3,6x,i3,e13.3,f10.4,i7/&
         &        //' Mat.     Bulk.D.    DispL    Fraction  Immobile WC')
160 format(i3,f13.4,3f10.4)
170 format(/'    Dif.w.      Dif.g.   ',50('-'),' (',i2,'.solute)')
180 format(2e12.4/' Mat.     KS         Nu         Beta      Henry&
         &  SinkL1     SinkS1     SinkG1     SinkL1`    SinkS1`    SinkG1`&
         &  SinkL0     SinkS0     SinkG0      Alfa')
190 Format(i4,1P,14e11.4)
200 format(/' Langmuir nonlinear adsorption isotherm for material ',i2)
210 format(/' Freundlich nonlinear adsorption isotherm for material ',i2)
220 format(/' No adsorption or linear adsorp. isotherm for material ',i2)
222 format(/' Physical non-equilibrium solute transport with mobile and imobile water.')
224 format(/' Chemical non-equilibrium solute transport with kinetic and equilibrium sorption sites.')
230 format(/' kTopCh      cTop(1...NS)'/i4,7x,20e10.3)
240 format(/' kBotCh      cBot(1...NS)'/i4,7x,20e10.3)
250 format(/' tPulse =   ',f15.3)
  End Subroutine ChemIn


  ! Source file SOLUTE.FOR |||||||||||||||||||||||||||||||||||||||||||||||

  !     To assemble and solve the solute transport equation
  !     Mass-lumping finite elements

  subroutine Solutes(x,dt,t,MatNum,thO,thN,vO,vN,Disp,B,D,E,F,&
       dtMaxC,TempO,TempN,thSat,vCorr,Sink,TLevel,dtMin,dtOpt,lWat,&
        xConv,tConv,lVapor,rBot,ierr,Beta,AtmBC,SinkF)

    Use Geometry, Only: NumNP, cRootMax
    Use Variables, Only: lNitrogen
    Logical :: lConv,lWat,lVapor,SinkF,lNEquil,AtmBC
    Real(dp) :: t,dt,dtMin,dtOpt,xConv,tConv,rBot
    Integer :: TLevel, MatNum(:),ierr
    Real(dp) :: dtMaxC
    Real(dp), Dimension(:) :: B,D,E,F,x,thO,thN,vO,vN,Disp,&
         TempO,TempN,thSat,vCorr,Sink,Beta
    Real(dp) :: alf,BN,DN,FN,dtMxC,dtOld,E1,D1,F1,Henry,Cour,Dg,dSurfT,Pecl,R,rMin,&
         Tr,TT
    Integer :: i,Iter,jS,Level,M
    Integer :: N,NLevel

    If(debug) Print *,'solutes'
    N=NumNP
    alf=1.-epsi
    IterC=1.
    NLevel=2
    Peclet=0.
    Courant=0.
    dtMaxC=1.e+30
    rMin=1.e-30
    Tr=293.15
    R=8.314
    !     Sequential first order decay goes into equilibrium phase (lNEquil=.false.) or nonequilbrium phase (lNEquil=.true.) 
    lNEquil=.false.

10  continue

    !     Loop on species in the chain

    Do jS=1,NS
       Iter=0
       chem_vTop(jS)=0.
       chem_vBot(jS)=0.
       cvCh0(jS)=0.
       cvCh1(jS)=0.
       cvChR(jS)=0.
       cvChIm(jS)=0.
       If(t-tPulse.Gt.dtMin.And..Not.AtmBC) Then
          cTop(jS)=0.
          cBot(jS)=0.
       end if
       if(kBotCh.lt.0) then
          if(vO(1).ge.0.) chem_vBot(jS)=alf*cBot(jS)  *vO(1)
          if(vO(1).lt.0.) chem_vBot(jS)=alf*Conc(jS,1)*vO(1)
          if(lVapor.and.rBot.eq.0.) chem_vBot(jS)=0.
       else if(kBotCh.eq.0) then
          chem_vBot(jS)=alf*Conc(jS,1)*vO(1)
       end if
       if(kTopCh.lt.0..and.TLevel.ne.1) then 
          if(vO(N).lt.0.) chem_vTop(jS)=alf*cTop(jS)*vO(N)
       end if
       if(kTopCh.eq.-2) then
          M=MatNum(N)
          Tr=293.15
          R=8.314
          TT=(TempO(N)+273.15-Tr)/R/(TempO(N)+273.15)/Tr
          Dg   =ChPar(6,M,jS)*Exp(TDep(6,jS)*TT)
          Henry=ChPar(10,M,jS)*Exp(TDep(10,jS)*TT)
          dSurfT=dSurf*Exp(TDep(9,jS)*TT)
          chem_vTop(jS)=chem_vTop(jS)+alf*Dg/dSurfT*Henry*Conc(jS,N)-Dg/dSurfT*cAtm
       end if
       if(.not.lLinear(jS)) then
          Do i=1,N
             cNew(i)=Conc(jS,i)
             if(.not.lEquil)       SorbN(i) =Sorb(jS,i)
             if(lBact.or.lDualNEq) SorbN2(i)=Sorb2(jS,i)
          End Do
       end if

       !       Root Solute Uptake
       If(SinkF .or. lNitrogen) Then
            call SetSSnk(jS,x,Beta,Sink,sSink,OmegaW,&
                         cRootMax(jS),lActRSU,OmegaS,SPot,rKM,cMin) !,ThN)
         Endif
       !       Iterative loop for a nonlinear adsorption isotherm
12     Iter=Iter+1
       if(.not.lLinear(jS)) then
          Do i=1,N
             cTemp(i)=cNew(i)
          End Do
       end if

       !       To construct the matrix equation
       do Level=1,NLevel

          !         Calculate the dispersion coefficients, retardation factors, source/
          !         decay coefficients, Peclet and Courant numbers, upstream weighting
          !         factors
          call Coeff(jS,Level,NLevel,x,Disp,vO,vN,thO,thN,&
               &               thSat,MatNum,TempN,TempO,&
               &               dt,Pecl,Cour,dtMxC,lLinear,lEquil,&
               &               lArtD,Iter,vCorr,&
               &               sSink,lBact,&
               &               iDualPor,xConv,tConv,&
               &               iMoistDep,lNEquil)
          Peclet=Max(Peclet,Pecl)
          Courant=Max(Courant,Cour)
          dtMaxC=Min(dtMaxC,dtMxC)

          !         Set up the matrix equation
          call MatSet(jS,Level,alf,dt,cBot,&
               x,thO,thN,vO,vN,Disp,B,&
               D,E,F,E1,D1,F1,BN,DN,FN,TempO,TempN,&
               dSurfT,MatNum,iDualPor,lVapor,&
               rBot,lBact)
          do i=1,N
             if(Level.eq.1) vO(i)=vO(i)+vCorr(i)
             if(Level.eq.2) vN(i)=vN(i)+vCorr(i)
          end do

          !         Calculate mass-transfer fluxes at the beginning of the time interval
          If(Level.Eq.1.And.Iter.Eq.1) &
               Call MassTran(jS,MatNum,TempO,lEquil,x,alf,&
               sSink,lBact,ThO,vO,iDualPor,xConv,tConv,lLinear)
       end do

       !       Solve matrix equation
       call BanSol(N,B,D,E,F)
       
       !       Test for convergence for nonlinear problem
       lConv=.true.
       do i=1,N
          if((NS.gt.1.and.Iter.eq.1).or.lDensity) cPrevO(i)=Conc(jS,i)
          if(lLinear(jS)) then
             Conc(jS,i)=max(F(i),0.0)
             if(Conc(jS,i).lt.1.e-30.and.Conc(jS,i).gt.0.) Conc(jS,i)=0.
          else
             cNew(i)=F(i)
             if(cNew(i).lt.1.0e-30) cNew(i)=0.
             if(abs(cNew(i)-cTemp(i)).gt.cTolA+cTolR*Conc(jS,i)) lConv=.false.
          end if
       end do
       If(.Not.lLinear(jS)) Then
          if(.not.lConv) then
             if(iter.lt.MaxItC) then
                goto 12
             else if(dt.gt.dtMin.and..not.lWat) then
                !              ierr=1
                dtOld=dt
                dt=Max(dt/3.,dtMin)
                dtOpt=dt
                t=t-dtOld+dt
                goto 10
             else
                ierr=1
             end if
          end if
          Do i=1,N
             Conc(jS,i)=cNew(i)
             if(.not.lEquil)       Sorb(jS,i) =SorbN(i)
             if(lBact.or.lDualNEq) Sorb2(jS,i)=SorbN2(i)
          End Do
       end if

       !       Calculate sorbed concentration for linear noneq. adsorption or 
       !       concentration in the imobile water.
       if(.not.lEquil.and.lLinear(jS)) &
            &    call SorbConc(jS,MatNum,TempN,dt,lBact,thN,vN,iDualPor,&
            &                  iMoistDep,xConv,tConv)

       !       Calculate mass-transfer fluxes at the end of the time interval
       call MassTran(jS,MatNum,TempN,lEquil,x,epsi,sSink,lBact,ThN,vN,iDualPor,&
            &                xConv,tConv,lLinear)

       !       Set up mass fluxes
       if(kTopCh.lt.0) then
          if(TLevel.ne.1) then
             if(vN(N).lt.0.) chem_vTop(jS)=chem_vTop(jS)+epsi*vN(N)*cTop(jS)
          else
             if(vN(N).lt.0.) chem_vTop(jS)=chem_vTop(jS)+     vN(N)*cTop(jS)
          end if
       else 
          chem_vTop(jS)=FN-BN*Conc(jS,N-1)-DN*Conc(jS,N)
       end if
       if(kTopCh.eq.-2) then
          M=MatNum(N)
          TT=(TempN(N)+273.15-Tr)/R/(TempN(N)+273.15)/Tr
          Dg   =ChPar( 6,M,jS)*Exp(TDep( 6,jS)*TT)
          Henry=ChPar(10,M,jS)*Exp(TDep(10,jS)*TT)
          chem_vTop(jS)=chem_vTop(jS)+epsi*Dg/dSurfT*Henry*Conc(jS,N)-Dg/dSurfT*cAtm
       end if
       if(kBotCh.lt.0) then
          if(vN(1).ge.0.) chem_vBot(jS)=chem_vBot(jS)+epsi*cBot(jS)  *vN(1)
          if(vN(1).lt.0.) chem_vBot(jS)=chem_vBot(jS)+epsi*Conc(jS,1)*vN(1)
          if(lVapor.and.rBot.eq.0.) chem_vBot(jS)=0.
       else if(kBotCh.eq.0) then 
          chem_vBot(jS)=chem_vBot(jS)+epsi*Conc(jS,1)*vN(1)
       else
          chem_vBot(jS)=D1*Conc(jS,1)+E1*Conc(jS,2)-F1
       end if
       IterC=max0(IterC,Iter)
       if(abs(chem_vTop(jS)).lt.rMin) chem_vTop(jS)=0.
       if(abs(chem_vBot(jS)).lt.rMin) chem_vBot(jS)=0.
       If(lLinear(jS)) Then
          ! Calculate linear Sorption
          Do i=1,n
             m=MatNum(i)
             SorbOut(jS,i)=ChPar(7,m,jS)*Conc(jS,i)
          Enddo
       Else If(freundlich(jS)) Then
          ! Calculate freundlich Sorption
          Do i=1,n
             m=MatNum(i)
             SorbOut(jS,i)=ChPar(7,m,jS)*Conc(jS,i)**ChPar(9,m,jS)
          Enddo
       Else If(langmuir(jS)) Then
          ! Calculate langmuir Sorption
          Do i=1,n
             m=MatNum(i)
             SorbOut(jS,i)=ChPar(7,m,jS)*Conc(jS,i) / (one_dp+ChPar(8,m,jS)*Conc(jS,i))
          Enddo
       Endif
    end do
    !     Calculate flux concentrations
    if(iConcType.eq.2)&
         Call FluxConc(x,vN,thN,thSat,MatNum,TempN,iDualPor,ThNewIm,1)    
  End subroutine Solutes

  !***********************************************************************    

  !     Calculate the dispersion coefficients, retardation factors, source/
  !     decay coefficients, Peclet and Courant numbers, upstream weighting
  !     factors

  Subroutine Coeff(jS,Level,NLevel,x,Disp,vO,vN,thO,&
       &                 thN,thSat,MatNum,TempN,TempO,&
       &                 dt,Peclet,Courant,dtMaxC,&
       &                 lLinear,lEquil,lArtD,Iter,vCorr,&
       &                 sSink,lBact,&
       &                 iDualPor,&
       &                 xConv,tConv,iMoistDep,&
       &                 lNEquil)
    Use Variables, Only: lNitrogen
    Implicit None
    Real(dp) :: dtMaxC,dt,Peclet,Courant,xConv,tConv
    Logical :: lLinear(:),lEquil,lArtD,lBact,lNEquil
    Real(dp), Dimension(:) :: x,Disp,vO,vN,thO,thN,thSat,TempO,TempN,sSink,vCorr
    Integer :: MatNum(:),jS,Level,NLevel,Iter,iDualPor,iMoistDep
    Real(dp) :: aa,Alfa1,Alfa2,cc,cG,cG1,cMid,CourMax,cprev,dConc,dConcS,&
         ddExp,derK,Dg,dHenry,dKs,Dc,DMobI,dNu,dRetard,dRetardS,dSConc,dSConcS,&
         Dw,dx,f1,f_em,fExp,fExpN,fExpO,fExpP,FlMacro,Frac,GamG,GamG1,GamG1i,&
         GamG1P,GamGi,GamL,GamL1,GamL1i,GamL1O,GamL1Oi,GamL1P,GamL1Pi,GamLi,&
         GamLO,GamLOi,GamS,GamS1,GamS1i,GamS1O,GamS1Oi,GamS1P,GamS1Pi,GamSi,&
         GamSO,GamSOi,Henry,Henryi,henryj,HenryN,HenryO,HenryP,&
         Omega,OmegaO,psi1,psi1O,psi2,psi2O,R,rKa1,rKa1O,rKa2,rKa2O,rKd1,rKd1O,&
         rKd2,rKd2O,ro,SConc,SConcO,SConcOS,SConcP,SConcPS,SConcS,SMax1,SMax1O,&
         SMax2,SMax2O,sMid,ss,ss1,ss2,SSorb,SSorb2,TauG,ThG,ThImob,ThImobO,&
         Thj,ThW,ThWO,Tr,TT,TTi,TTj,TTN,TTO,v,vj,xKs,xKsN,xKsO,xKsP,xMuG,xMuGO,&
         xMuL,xMuLO,xMuS,xMuSO,xNu,xNuN,xNuO,xNuP,Dep
    Integer :: i,iPsi1,iPsi2,j,jS1,k,M,nmat

    If(debug) Print *,'coeff'
    !     Inicialization
    XNuO=0.0
    HenryP=0.0
    HenryI=0.0
    HenryJ=0.0
    GamG1P=0.0
    GamL1P=0.0
    GamS1P=0.0
    fExpO=0.0
    fExpP=0.0
    TTi=0.0
    TTj=0.0
    XNuP=0.0
    cprev=0.0
    jS1=jS-1
    Peclet=0.
    Courant=0.
    CourMax=1.
    dtMaxC=1.e+30
    Tr=293.15
    R=8.314
    nmat=Ubound(MatNum,1)
    do i=nmat,1,-1
       j=i+1
       k=i-1
       M=MatNum(i)
       if(Level.eq.NLevel) then
          ThW=ThN(i)
          ThWO=ThO(i)
          ThG=Max(0.,thSat(M)-ThW)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then 
             ThImob=ChPar(4,M,jS)
             ThImobO=ThImob
             ThW=Max(ThW-ThImob,0.001)
          end if
          if(iDualPor.gt.0) then
             ThImob=ThNewIm(i)
             ThImobO=ThOldIm(i)
          end if
          v=vN(i)
          if(i.ne.nmat) then
             vj=vN(j)
             Thj=ThN(j)
             if(lMobIm(M).and.iDualPor.eq.0.or.lBact) Thj=Max(Thj-ThImob,0.001)
          end if
          TT=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          if(jS.gt.1) cPrev=Conc(jS-1,i)
       else
          ThW=ThO(i)
          ThG=Max(0.,thSat(M)-ThW)
          if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
             ThImob=ChPar(4,M,jS)
             ThW=Max(ThW-ThImob,0.001)
          end if
          if(iDualPor.gt.0) then
             ThImob=ThOldIm(i)
          end if
          v=vO(i)
          if(i.ne.nmat) then
             vj=vO(j)
             Thj=ThO(j)
             if(lMobIm(M).and.iDualPor.eq.0.or.lBact) &
                  &                                       Thj=Max(Thj-ThImob,0.001)
          end if
          TT=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          if(jS.gt.1) cPrev=cPrevO(i)
       end if
       If(lNitrogen) cPrev=0
       !       Temperature dependence
       f1=1.
       ro   =ChPar(1,M,jS)*Exp(TDep(1,jS)*TT)
       Frac =ChPar(3,M,jS)*Exp(TDep(3,jS)*TT)
       Dw   =ChPar(5,M,jS)*Exp(TDep(5,jS)*TT)
       Dg   =ChPar(6,M,jS)*Exp(TDep(6,jS)*TT)
       xKs  =ChPar(7,M,jS)*Exp(TDep(7,jS)*TT)
       xNu  =ChPar(8,M,jS)*Exp(TDep(8,jS)*TT)
       fExp =ChPar(9,M,jS)
       Henry=ChPar(10,M,jS)*Exp(TDep(10,jS)*TT)
       if(iMoistDep.gt.0)  &
            &    f1=rMD(M,jS,1,1,ThW,iMoistDep)
       GamL  =ChPar(11,M,jS)*Exp(TDep(11,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,10,1,ThImob,iMoistDep)
       GamLi =ChPar(11,M,jS)*Exp(TDep(11,jS)*TT)*f1    ! reaction in the immobile phase
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,2,2,ThW,iMoistDep)
       GamS  =ChPar(12,M,jS)*Exp(TDep(12,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,11,2,ThImob,iMoistDep)
       GamSi =ChPar(12,M,jS)*Exp(TDep(12,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,3,3,ThW,iMoistDep)
       GamG  =ChPar(13,M,jS)*Exp(TDep(13,jS)*TT)*f1
       GamGi=GamG
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,4,4,ThW,iMoistDep)
       GamL1 =ChPar(14,M,jS)*Exp(TDep(14,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,12,4,ThImob,iMoistDep)
       GamL1i=ChPar(14,M,jS)*Exp(TDep(14,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,5,5,ThW,iMoistDep)
       GamS1 =ChPar(15,M,jS)*Exp(TDep(15,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,13,5,ThImob,iMoistDep)
       GamS1i=ChPar(15,M,jS)*Exp(TDep(15,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,6,6,ThW,iMoistDep)
       GamG1 =ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,7,7,ThW,iMoistDep)
       xMuL  =ChPar(17,M,jS)*Exp(TDep(17,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,8,8,ThW,iMoistDep)
       xMuS  =ChPar(18,M,jS)*Exp(TDep(18,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,9,9,ThW,iMoistDep)
       xMuG  =ChPar(19,M,jS)*Exp(TDep(19,jS)*TT)*f1
       Omega =ChPar(20,M,jS)*Exp(TDep(20,jS)*TT)
       f_em=1.
       If(lDualNEq) f_em  =ChPar(13,M,jS)*Exp(TDep(13,jS)*TT)
       If(lDualNEq) OmegaS=ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)
       if(lBact) then
          Dg=0. 
          GamG =0.
          GamGi=0.
          GamL1=0.
          GamL1i=0.
          GamS1=0.
          GamS1i=0.
          GamG1=0.
          GamG1i=0.
          xMuL =0.
          xMuS =0.
          xMuG =0.
          Omega=0.
          SMax2 =ChPar(15,M,jS)*Exp(TDep(15,jS)*TT)
          rKa2  =ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)
          rKd2  =ChPar(17,M,jS)*Exp(TDep(17,jS)*TT)
          SMax1 =ChPar(18,M,jS)*Exp(TDep(18,jS)*TT)
          rKa1  =ChPar(19,M,jS)*Exp(TDep(19,jS)*TT)
          rKd1  =ChPar(20,M,jS)*Exp(TDep(20,jS)*TT)
          iPsi1=0
          iPsi2=0
          if(.not.lFiltr) iPsi2=int(ChPar(13,M,jS))
          if(.not.lFiltr) iPsi1=int(ChPar(14,M,jS))
          if(iPsi1.eq.0.and.SMax1.gt.0.) iPsi1=1
          if(iPsi2.eq.0.and.SMax2.gt.0.) iPsi2=1
          if(iPsi1.ge.3.or.iPsi2.ge.3) &
               &      Dc=ChPar(6,M,jS)*Exp(TDep(6,jS)*TT)
          if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(15,M,jS)
          if(Level.eq.NLevel) then
             ss1=SorbN(i)
             ss2=SorbN2(i)
          else
             ss1=Sorb(jS,i)
             ss2=Sorb2(jS,i)
          end if
          psi1=1.
          psi2=1.
          if(iPsi1.gt.0) call Blocking(iPsi1,SMax1,psi1,x(i),ss1,dc,aa)
          if(iPsi2.gt.0) call Blocking(iPsi2,SMax2,psi2,x(i),ss2,dc,aa)

          !         recalculate ka1 and ka2 based on filtration theory
          if(lFiltr) then
             GamG =0.
             GamL1=0.
             Dc=ChPar(13,M,jS)*Exp(TDep(13,jS)*TT)
             Dep=ChPar(14,M,jS)*Exp(TDep(14,jS)*TT)
             Alfa1=rKa1
             Alfa2=rKa2
             call Deposit(rKa1,rKa2,Dc,Dep,Alfa1,Alfa2,ThW,v,TempN(i),&
                  &                   xConv,tConv)
          end if
       end if
       if(jS.gt.1) then
          xKsP  =ChPar( 7,M,jS1)*Exp(TDep(7,jS-1)*TT)
          xNuP  =ChPar( 8,M,jS1)*Exp(TDep(8,jS-1)*TT)
          fExpP =ChPar( 9,M,jS1) !*exp(TDep(9,jS-1)*TT)
          HenryP=ChPar(10,M,jS1)*Exp(TDep(10,jS-1)*TT)
          f1=1.
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS-1,4,4,ThW,iMoistDep)
          GamL1P =ChPar(14,M,jS1)*Exp(TDep(14,jS-1)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS-1,12,4,ThImob,iMoistDep)
          GamL1Pi=ChPar(14,M,jS1)*Exp(TDep(14,jS-1)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS-1,5,5,ThW,iMoistDep)
          GamS1P =ChPar(15,M,jS1)*Exp(TDep(15,jS-1)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS-1,13,5,ThImob,iMoistDep)
          GamS1Pi=ChPar(15,M,jS1)*Exp(TDep(15,jS-1)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS-1,6,6,ThW,iMoistDep)
          GamG1P=ChPar(16,M,jS1)*Exp(TDep(16,jS-1)*TT)*f1
          if(lBact) then
             GamL1P=0.
             GamL1Pi=0.
             GamS1P=0.
             GamS1Pi=0.
             GamG1P=0.
          end if
       end if
       if(Level.eq.NLevel) then
          TTO=(TempO(i)+273.15-Tr)/R/(TempO(i)+273.15)/Tr
          xKsO  =ChPar( 7,M,jS)*Exp(TDep( 7,jS)*TTO)
          xNuO  =ChPar( 8,M,jS)*Exp(TDep( 8,jS)*TTO)
          fExpO =ChPar( 9,M,jS)
          HenryO=ChPar(10,M,jS)*Exp(TDep(10,jS)*TTO)
          f1=1.
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,1,1,ThO(i),iMoistDep)
          GamLO =ChPar(11,M,jS)*Exp(TDep(11,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,10,1,ThImobO,iMoistDep)
          GamLOi=ChPar(11,M,jS)*Exp(TDep(11,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,2,2,ThO(i),iMoistDep)
          GamSO =ChPar(12,M,jS)*Exp(TDep(12,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,11,2,ThImobO,iMoistDep)
          GamSOi=ChPar(12,M,jS)*Exp(TDep(12,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,4,4,ThO(i),iMoistDep)
          GamL1O =ChPar(14,M,jS)*Exp(TDep(14,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,12,4,ThImobO,iMoistDep)
          GamL1Oi=ChPar(14,M,jS)*Exp(TDep(14,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,5,5,ThO(i),iMoistDep)
          GamS1O=ChPar(15,M,jS)*Exp(TDep(15,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,13,5,ThImobO,iMoistDep)
          GamS1Oi=ChPar(15,M,jS)*Exp(TDep(15,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,7,7,ThO(i),iMoistDep)
          xMuLO =ChPar(17,M,jS)*Exp(TDep(17,jS)*TTO)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,8,8,ThO(i),iMoistDep)
          xMuSO =ChPar(18,M,jS)*Exp(TDep(18,jS)*TTO)*f1
          OmegaO=ChPar(20,M,jS)*Exp(TDep(20,jS)*TTO)
          dKs   =(xKs  -  xKsO)/dt
          dNu   =(xNu  -  xNuO)/dt
          ddExp =(fExp - fExpO)/dt
          dHenry=(Henry-HenryO)/dt
          if(i.ne.1)     TTi=(TempN(k)+273.15-Tr)/R/(TempN(k)+273.15)/Tr
          if(i.ne.nmat) TTj=(TempN(j)+273.15-Tr)/R/(TempN(j)+273.15)/Tr
          if(lBact) then
             GamS1O=0.
             GamS1Oi=0.
             xMuLO =0.
             xMuSO =0.
             xMuGO =0.
             OmegaO=0.
             SMax2O =ChPar(15,M,jS)*Exp(TDep(15,jS)*TTO)
             rKa2O  =ChPar(16,M,jS)*Exp(TDep(16,jS)*TTO)
             rKd2O  =ChPar(17,M,jS)*Exp(TDep(17,jS)*TTO)
             SMax1O =ChPar(18,M,jS)*Exp(TDep(18,jS)*TTO)
             rKa1O  =ChPar(19,M,jS)*Exp(TDep(19,jS)*TTO)
             rKd1O  =ChPar(20,M,jS)*Exp(TDep(20,jS)*TTO)
             iPsi1=0
             iPsi2=0
             if(.not.lFiltr) iPsi2=int(ChPar(13,M,jS))
             if(.not.lFiltr) iPsi1=int(ChPar(14,M,jS))
             if(iPsi1.eq.0.and.SMax1O.gt.0.) iPsi1=1
             if(iPsi2.eq.0.and.SMax2O.gt.0.) iPsi2=1
             if(iPsi1.ge.3.or.iPsi2.ge.3) &
                  &        Dc=ChPar(6,M,jS)*Exp(TDep(6,jS)*TT)
             if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(15,M,jS)
             psi1O=1.
             psi2O=1.
             if(iPsi1.gt.0) &
                  &        call Blocking(iPsi1,SMax1O,psi1O,x(i),Sorb(jS,i),dc,aa)
             if(iPsi2.gt.0) &
                  &        call Blocking(iPsi2,SMax2O,psi2O,x(i),Sorb2(jS,i),dc,aa)
             if(lFiltr) then
                GamL1O=0.
                Dc=ChPar(13,M,jS)*Exp(TDep(13,jS)*TTO)
                Dep=ChPar(14,M,jS)*Exp(TDep(14,jS)*TTO)
                Alfa1=rKa1O
                Alfa2=rKa2O
                call Deposit(rKa1O,rKa2O,Dc,Dep,Alfa1,Alfa2,ThWO,vO(i),&
                     &                     TempO(i),xConv,tConv)
             end if
          end if
       else
          TTN=(TempN(i)+273.15-Tr)/R/(TempN(i)+273.15)/Tr
          xKsN  =ChPar( 7,M,jS)*Exp(TDep( 7,jS)*TTN)
          xNuN  =ChPar( 8,M,jS)*Exp(TDep( 8,jS)*TTN)
          fExpN =ChPar( 9,M,jS) !*exp(TDep( 9,jS)*TTN)
          HenryN=ChPar(10,M,jS)*Exp(TDep(10,jS)*TTN)
          dKs   =(xKsN  -  xKs)/dt
          dNu   =(xNuN  -  xNu)/dt
          ddExp =(fExpN - fExp)/dt
          dHenry=(HenryN-Henry)/dt
          if(i.ne.1)     TTi=(TempO(k)+273.15-Tr)/R/(TempO(k)+273.15)/Tr
          if(i.ne.nmat) TTj=(TempO(j)+273.15-Tr)/R/(TempO(j)+273.15)/Tr
       end if
       If(i.Ne.1) Henryi=ChPar(10,MatNum(k),jS)*Exp(TDep(10,jS)*TTi)
       If(i.Ne.nmat) Henryj=ChPar(10,MatNum(j),jS)*Exp(TDep(10,jS)*TTj)

       dSConc=1.
       dConc=1.
       SConcP=1.
       SConc=1.
       SConcO=1.
       dRetard=0.

       SConcS=1.
       SConcOS=1.
       dSConcS=1.
       dConcS=1.
       SConcPS=1.
       dRetardS=0.
       cMid=0
       !       Effects of nonlinear adsorption
       if(.not.lLinear(jS)) then
          cc=Conc(jS,i)
          cMid=(Conc(jS,i)+cNew(i))/2.
          if(Level.eq.NLevel) cc=cNew(i)
          if(cc.gt.0.) then
             dSConc=fExp*cc**(fExp-1.)/(1.+xNu*cc**fExp)**2
             SConc =     cc**(fExp-1.)/(1.+xNu*cc**fExp)
          end if
          if(cMid.gt.0.) then
             dConc=fExp*cMid**(fExp-1.)/(1.+xNu*cMid**fExp)**2
             dRetard=cMid**fExp/(1.+xNu*cMid**fExp)*dKs-&
                  &           xKs*cMid**(2.*fExp)/(1.+xNu*cMid**fExp)**2*dNu+&
                  &           xKs*log(cMid)*cMid**fExp/(1.+xNu*cMid**fExp)**2*ddExp
          end if
          if(Level.eq.NLevel.and..not.lEquil.and.Conc(jS,i).gt.0.)&
               &      SConcO=Conc(jS,i)**(fExpO-1.)/(1.+xNuO*Conc(jS,i)**fExpO)
          if(lMobIm(M).or.iDualPor.gt.0) then     ! mobile-immobile model
             ss=Sorb(jS,i)
             sMid=(Sorb(jS,i)+SorbN(i))/2.
             if(Level.eq.NLevel) ss=SorbN(i)
             if(ss.gt.0.) then
                dSConcS=fExp*ss**(fExp-1.)/(1.+xNu*ss**fExp)**2
                SConcS =     ss**(fExp-1.)/(1.+xNu*ss**fExp)
             end if
             if(sMid.gt.0.) then
                dConcS=fExp*sMid**(fExp-1.)/(1.+xNu*sMid**fExp)**2
                dRetardS=sMid**fExp/(1.+xNu*sMid**fExp)*dKs-&
                     &            xKs*sMid**(2.*fExp)/(1.+xNu*sMid**fExp)**2*dNu+&
                     &            xKs*log(sMid)*sMid**fExp/(1.+xNu*sMid**fExp)**2*ddExp
             end if
             if(Level.eq.NLevel.and..not.lEquil.and.Sorb(jS,i).gt.0.)&
                  &        SConcOS=Sorb(jS,i)**(fExpO-1.)/(1.+xNuO*Sorb(jS,i)**fExpO)
          end if
       else
          if(Conc(jS,i).gt.0.) dRetard=Conc(jS,i)*dKs
          if(lMobIm(M).or.iDualPor.gt.0) then     ! mobile-immobile model
             if(Sorb(jS,i).gt.0) dRetardS=Sorb(jS,i)*dKs
          end if
       end if
       if(jS.gt.1) then
          if(.not.lLinear(jS-1)) then
             if(cPrev.gt.0.)&
                  &        SConcP=cPrev**(fExpP-1.)/(1.+xNuP*cPrev**fExpP)
             if(Sorb(jS-1,i).gt.0.)&
                  &        SConcPS=Sorb(jS-1,i)**(fExpP-1.)/&
                  &                                     (1.+xNuP*Sorb(jS-1,i)**fExpP)
          end if
       end if

       !       Calculate the retardation factors
       Retard(i)=(ro*Frac*f_em*xKs*dConc+ThG*Henry)/ThW+1.
       !       Calculate the dispersion coefficients
       call Disper(i,M,dt,lArtD,&
            &              iDualPor,Level,NLevel,Disp,thSat,&
            &              ThImob,ThW,ThG,v,Dw,Dg,Henry,ro,Frac,xKs,fExp,xNu,&
            &              cMid,dSConc,TauG,lBact)

       !       Calculate the adsorbed concentration on kinetic sites or
       !       the concentration in an imobile zone, before solving matrix equation
       if(.not.lEquil)&
            &    call NEquil(i,jS,M,&
            &                SSorb,SSorb2,lLinear,lBact,Level,&
            &                NLevel,dt,ro,xKs,xKsO,cc,SConc,SConcO,SConcS,&
            &                SConcOS,dSConcS,xMuL,xMuLO,xMuS,xMuSO,dRetardS,&
            &                GamLi,GamL1i,GamLOi,GamL1Oi,GamSi,GamS1i,GamSOi,&
            &                GamS1Oi,Omega,OmegaO,rKa1,rKa1O,rKa2,rKa2O,rKd1,&
            &                rKd1O,rKd2,rKd2O,ThW,ThWO,psi1,psi1O,psi2,psi2O,&
            &                DMobI,Frac,ThImob,ThImobO,iDualPor,FlMacro,&
            &                lNEquil,GamL1Pi,GamS1Pi,xKsP,SConcPS,&
            &                f_em,OmegaS)

       !       Calculate zero-order coefficient g0
       g0(i)=xMuL*ThW+Frac*f_em*ro*xMuS+ThG*xMuG-sSink(i)
       q0(i)=xMuL*ThW+          ro*xMuS+ThG*xMuG
       if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
             g0(i)=g0(i)+Omega*SSorb
             if(iDualPor.gt.0.and.SinkIm(i).le.0) g0(i)=g0(i)-FlMacro
             if(lDualNEq) g0(i)=g0(i)+OmegaS*ro*SSorb2
          else if(.not.lBact) then
             g0(i)=g0(i)+Omega*ro*SSorb
          else if(lBact) then
             g0(i)=g0(i)+rKd1*ro*SSorb+rKd2*ro*SSorb2
          end if
       end if
       if(jS.gt.1) then
          cG=cPrev*(GamL1P*ThW+ro*Frac*f_em*xKsP*GamS1P*SConcP+ThG*HenryP*GamG1P)
          cG1=cG
          if(.not.lEquil) then
             if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
                aa=Sorb(jS-1,i)*(ThImob*GamL1Pi+&
                     &                               (1.-Frac)*ro*GamS1Pi*xKsP*SConcPS)
                if(.not.lNEquil) cG=cG+aa
                cG1=cG1+aa    
                if(lDualNEq) then
                   aa=GamS1Pi*ro*Sorb2(jS-1,i)
                   if(.not.lNEquil) cG=cG+aa
                   cG1=cG1+aa
                end if
             else if(.not.lBact) then
                aa=GamS1Pi*ro*Sorb(jS-1,i)
                if(.not.lNEquil) cG=cG+aa
                cG1=cG1+aa
             else if(lBact) then
                write(*,*) 'Attachment/dettachment model is implemented &
                     &only for one solute'
                write(*,*)'Press Enter to continue'
                read(*,*)
                stop
             end if
          end if
          g0(i)=g0(i)+cG
          q0(i)=q0(i)+cG1
       end if
       if(cMid.gt.0.) g0(i)=g0(i)-ro*Frac*f_em*dRetard
       !       Calculate first-order coefficient g1
       g1(i)=-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*f_em*xKs*SConc-&
            &         (GamG+GamG1)*ThG*Henry
       !        if(Level.eq.NLevel) g1(i)=g1(i)-ThG*dHenry-Henry*(ThWO-ThW)/dt
       if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then ! mobile-immobile model
             g1(i)=g1(i)-Omega
             if(iDualPor.gt.0.and.SinkIm(i).gt.0) g1(i)=g1(i)-SinkIm(i)
             if(Level.eq.NLevel.and.lLinear(jS))&
                  &         g1(i)=g1(i)+Omega*dt*Omega/DMobI
             if(lDualNEq) then
                g1(i)=g1(i)-OmegaS*ro*Frac*(1.-f_em)*SConc*xKs
                if(Level.eq.NLevel.and.lLinear(jS)) g1(i)=g1(i)+OmegaS*ro*&
                     &                  (dt*OmegaS*Frac*(1.-f_em)*xKs/&
                     &                  (2.+dt*(OmegaS+GamSi+GamS1i)))
             end if
          else if(.not.lBact) then                             ! two-site sorption model
             g1(i)=g1(i)-Omega*ro*(1.-Frac)*SConc*xKs
             if(Level.eq.NLevel.and.lLinear(jS)) g1(i)=g1(i)+Omega*ro*&
                  &             (dt*Omega*(1.-Frac)*xKs/(2.+dt*(Omega+GamSi+GamS1i)))
          else if(lBact) then                                  ! filtration model
             g1(i)=g1(i)-ThW*(rKa1*psi1+rKa2*psi2)
             if(Level.eq.NLevel.and.lLinear(jS)) g1(i)=g1(i)+dt*ThW*&
                  &               (rKd1*rKa1/(2.+dt*(rKd1+GamSi))+&
                  &                rKd2*rKa2/(2.+dt*(rKd2+GamSi)))
          end if
       end if
       q1(i)=(-(GamL+GamL1)*ThW-(GamS+GamS1)*ro*Frac*f_em*xKs*SConc-&
            &          (GamG+GamG1)*ThG*Henry)*Conc(jS,i)
       if(.not.lEquil) then
          if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
             q1(i)=q1(i)-Sorb(jS,i)*(ThImob*(GamLi+GamL1i)+&
                  &                  (1.-Frac)*ro*xKs*SConcS*(GamSi+GamS1i))
             if(lDualNEq) q1(i)=q1(i)-(GamSi+GamS1i)*ro*Sorb2(jS,i)
          else if(.not.lBact) then
             q1(i)=q1(i)-(GamSi+GamS1i)*ro*Sorb(jS,i)
          else if(lBact) then
             q1(i)=q1(i)-ro*(GamSi+GamS1i)*(Sorb(jS,i)+Sorb2(jS,i))
          end if
       end if

       !       Velocity corrections
       if(i.eq.1) then
          dx=x(2)-x(1)
          derK=(Henryj-Henry)/dx
       else if(i.eq.nmat) then
          dx=x(nmat)-x(nmat-1)
          derK=(Henry-Henryi)/dx
       else
          dx=(x(j)-x(k))/2.
          derK=(Henryj-Henryi)/dx
       end if
       vCorr(i)=ThG*Dg*TauG*derK
       if(Level.eq.1)      vO(i)=vO(i)-vCorr(i)
       if(Level.eq.NLevel) vN(i)=vN(i)-vCorr(i)

       !       Calculate the maximum local Peclet and Courant numbers
       call PeCour(i,j,nmat,Level,NLevel,lArtD,dt,x,v,ThW,vj,&
            &              Thj,Disp,Peclet,Courant,CourMax,dtMaxC,&
            &              Iter)
    end do
  End Subroutine Coeff

  !***********************************************************************

  Subroutine MatSet(jS,Level,alf,dt,&
       &                  cBot,x,thO,thN,vO,vN,Disp,&
       &                  B,D,E,F,E1,D1,F1,BN,DN,FN,&
       &                  TempO,TempN,dSurf,MatNum,&
       &                  iDualPor,lVapor,rBot,lBact)
    Use geometry, Only: NumNP
    Logical :: lVapor,lBact
    Real(dp), Dimension(:) :: B,D,E,F,cBot,x,thO,thN,vO,vN,Disp,TempO,TempN
    Integer :: MatNum(:),jS,Level,iDualPor
    Real(dp) :: alf,dt,E1,D1,F1,BN,DN,FN,rBot,dSurf
    Real(dp) :: a1,b1,Dg,dx,F2,FE,Henry,R,ThImob,Tr,TT
    Integer :: i,M,N
    Real(dp),save :: x1,x2

    If(debug) Print *,'matset'
    N=NumNP
    do i=1,N
       M=MatNum(i)
       if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
          ThImob=ChPar(4,M,jS)
          if(ThImob.gt.thO(i)) write(*,*) "Warning !!! ThImob > Theta"
          thN(i)=max(thN(i)-ThImob,0.001)
          thO(i)=max(thO(i)-ThImob,0.001)
       end if
    end do

    !     Lower boundary condition
    b1=x(2)-x(1)
    If(Level.Eq.1) Then
       F1=           Conc(jS,1)*&
            &     (b1/2./dt*thO(1)*Retard(1)+&
            &      alf*(-(thO(1)*Disp(1)+thO(2)*Disp(2))/b1/2.-&
            &           ((2.+3.*wc(1))*vO(1)+vO(2))/6.+&
            &           b1/12.*(3.*g1(1)+g1(2))))+&
            &                Conc(jS,2)*&
            &     alf*((thO(1)*Disp(1)+thO(2)*Disp(2))/b1/2.-&
            &          (vO(1)+(2.-3.*wc(1))*vO(2))/6.+b1/12.*(g1(1)+g1(2)))+&
            &                alf*b1/6.*(2.*g0(1)+g0(2))

       !       3. type  BC
       if(kBotCh.eq.-1) F(1)=F1+alf*cBot(jS)*vO(1)
    else
       E1=epsi*(-(thN(1)*Disp(1)+thN(2)*Disp(2))/b1/2.+&
            &           (vN(1)+(2.-3.*wc(1))*vN(2))/6.-b1/12.*(g1(1)+g1(2)))
       D1=b1/2./dt*thN(1)*Retard(1)+&
            &     epsi*((thN(1)*Disp(1)+thN(2)*Disp(2))/b1/2.+&
            &           ((2.+3.*wc(1))*vN(1)+vN(2))/6.-b1/12.*(3.*g1(1)+g1(2)))
       F2=epsi*b1/6.*(2.*g0(1)+g0(2))
       F1=F1+F2

       !       1.type BC
       if(kBotCh.eq.1) then
          D(1)=1.
          E(1)=0.
          F(1)=cBot(jS)
       end if

       !       3. type  BC
       if(kBotCh.eq.-1) then
          if(vN(1).gt.0..or.(lVapor.and.rBot.eq.0.)) then
             E(1)=E1
             D(1)=D1
             F(1)=F(1)+F2+epsi*cBot(jS)*vN(1)
          else
             D(1)=-1.
             E(1)=1.
             F(1)=0.
          end if
       end if

       !       Free drainage
       if(kBotCh.eq.0) then
          D(1)=-1.
          E(1)=1.
          F(1)=0.
       end if
    end if

    do i=2,N-1
       a1=b1
       b1=x(i+1)-x(i)
       dx=(x(i+1)-x(i-1))/2.
       if(Level.eq.1) then
          F(i)=       Conc(jS,i-1)*&
               &       alf*((thO(i-1)*Disp(i-1)+thO(i)*Disp(i))/a1/2.+&
               &            ((2.+3.*wc(i-1))*vO(i-1)+vO(i))/6.+&
               &            a1/12.*(g1(i-1)+g1(i)))+&
               &                Conc(jS,i)*&
               &      (dx/dt*thO(i)*Retard(i)+&
               &      alf*(-(thO(i-1)*Disp(i-1)+thO(i)*Disp(i))/a1/2.-&
               &           (thO(i+1)*Disp(i+1)+thO(i)*Disp(i))/b1/2.-&
               &           (vO(i+1)+3.*(wc(i-1)+wc(i))*vO(i)-vO(i-1))/6.+&
               &           (a1*(g1(i-1)+3.*g1(i))+b1*(3.*g1(i)+g1(i+1)))/12.))+&
               &                Conc(jS,i+1)*&
               &      alf*((thO(i+1)*Disp(i+1)+thO(i)*Disp(i))/b1/2.-&
               &           (vO(i)+(2.-3.*wc(i))*vO(i+1))/6.+&
               &           b1/12.*(g1(i)+g1(i+1)))+&
               &              alf*(a1*(g0(i-1)+2.*g0(i))+b1*(2.*g0(i)+g0(i+1)))/6.
       else
          B(i)=epsi*(-(thN(i-1)*Disp(i-1)+thN(i)*Disp(i))/a1/2.-&
               &               ((2.+3.*wc(i-1))*vN(i-1)+vN(i))/6.-&
               &               a1/12.*(g1(i-1)+g1(i)))
          D(i)=dx/dt*thN(i)*Retard(i)+&
               &         epsi*((thN(i-1)*Disp(i-1)+thN(i)*Disp(i))/a1/2.+&
               &               (thN(i+1)*Disp(i+1)+thN(i)*Disp(i))/b1/2.+&
               &               (vN(i+1)+3.*(wc(i-1)+wc(i))*vN(i)-vN(i-1))/6.-&
               &               (a1*(g1(i-1)+3.*g1(i))+b1*(3.*g1(i)+g1(i+1)))/12.)
          E(i)=epsi*(-(thN(i+1)*Disp(i+1)+thN(i)*Disp(i))/b1/2.+&
               &               (vN(i)+(2.-3.*wc(i))*vN(i+1))/6.-&
               &               b1/12.*(g1(i)+g1(i+1)))
          F(i)=F(i)+epsi*(a1*(g0(i-1)+2.*g0(i))+b1*(2.*g0(i)+g0(i+1)))/6.
       end if
    end do

    !     Upper boundary condition
    If(Level.Eq.1) Then
       x1=alf*((thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+&
            ((2.+3.*wc(N-1))*vO(N-1)+vO(N))/6.+b1/12.*(g1(N-1)+g1(N)))
       x2=(b1/2./dt*thO(N)*Retard(N)+&
            &      alf*(-(thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+&
            &           (vO(N-1)+(2.-3.*wc(N-1))*vO(N))/6.+&
            &           b1/12.*(g1(N-1)+3*g1(N))))+&
            &                alf*b1/6.*(g0(N-1)+2.*g0(N))

       FN=           Conc(jS,N-1)*&
            &    alf*((thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+&
            &       ((2.+3.*wc(N-1))*vO(N-1)+vO(N))/6.+b1/12.*(g1(N-1)+g1(N)))+&
            &                Conc(jS,N)*&
            &     (b1/2./dt*thO(N)*Retard(N)+&
            &      alf*(-(thO(N-1)*Disp(N-1)+thO(N)*Disp(N))/b1/2.+&
            &           (vO(N-1)+(2.-3.*wc(N-1))*vO(N))/6.+&
            &           b1/12.*(g1(N-1)+3*g1(N))))+&
            &                alf*b1/6.*(g0(N-1)+2.*g0(N))

       !       3. type BC
       if(kTopCh.le.0) then
          F(N)=FN
          if(vO(N).lt.0.) F(N)=F(N)-alf*vO(N)*cTop(jS)
          if(kTopCh.eq.-2) then
             M=MatNum(N)
             Tr=293.15
             R=8.314
             TT=(TempO(N)+273.15-Tr)/R/(TempO(N)+273.15)/Tr
             Dg=ChPar(6,M,jS)*Exp(TDep(6,jS)*TT)
             Henry=ChPar(10,M,jS)*Exp(TDep(10,jS)*TT)
             F(N)=F(N)-alf*Dg/dSurf*Henry*Conc(jS,N)+Dg/dSurf*cAtm
          end if
       end if
    else
       BN=epsi*(-(thN(N-1)*Disp(N-1)+thN(N)*Disp(N))/b1/2.-&
            &        ((2.+3.*wc(N-1))*vN(N-1)+vN(N))/6.-b1/12.*(g1(N-1)+g1(N)))
       DN=b1/2./dt*thN(N)*Retard(N)+&
            &     epsi*((thN(N-1)*Disp(N-1)+thN(N)*Disp(N))/b1/2.-&
            &           (vN(N-1)+(2.-3.*wc(N-1))*vN(N))/6.-&
            &           b1/12.*(g1(N-1)+3.*g1(N)))
       FE=epsi*b1/6.*(g0(N-1)+2.*g0(N))
       FN=FN+FE

       !       1. type BC
       if(kTopCh.gt.0) then
          B(N)=0.
          D(N)=1.
          F(N)=cTop(jS)

          !       3. type BC
       else
          B(N)=BN
          D(N)=DN
          F(N)=F(N)+FE
          if(vN(N).lt.0.) F(N)=F(N)-epsi*vN(N)*cTop(jS)
          if(kTopCh.eq.-2) then
             M=MatNum(N)
             Tr=293.15
             R=8.314
             TT=(TempN(N)+273.15-Tr)/R/(TempN(N)+273.15)/Tr
             Dg=ChPar(6,M,jS)*Exp(TDep(6,jS)*TT)
             Henry=ChPar(10,M,jS)*Exp(TDep(10,jS)*TT)
             D(N)=D(N)+epsi*Dg/dSurf*Henry
          end if
       end if
    end if

    Do i=1,N
       M=MatNum(i)
       if(lMobIm(M).and.iDualPor.eq.0.or.lBact) then
          ThImob=ChPar(4,M,jS)
          thN(i)=thN(i)+ThImob
          thO(i)=thO(i)+ThImob
       end if
    End Do
  End Subroutine MatSet

  !************************************************************************

  !     Solve matrix equation

  Subroutine BanSol(N,A,B,C,F)
    integer :: N
    Real(dp) :: A(:),B(:),C(:),F(:)
    Integer :: i,j

    do i=2,N
       B(i)=B(i)-A(i)*C(i-1)/B(i-1)
       F(i)=F(i)-A(i)*F(i-1)/B(i-1)
    end do
    F(N)=F(N)/B(N)
    do i=2,N
       j=N-i+1
       ! checking for underflow
       !If(abs(F(j+1))<1e-300_dp) F(j+1)=0
       !Write(*,'(a,1p,99e13.5)') 'f,c,f1,b,cf=',F(j),C(j),F(j+1),B(j),C(j)*F(j+1)
       F(j)=(F(j)-C(j)*F(j+1))/B(j)
    end do
    return
  End Subroutine BanSol

  !************************************************************************

  !     Calculate the maximum local Peclet and Courant numbers

  Subroutine PeCour(i,j,NumNP,Level,NLevel,lArtD,dt,x,v,ThW,&
       &                  vj,Thj,Disp,Peclet,Courant,CourMax,&
       &                  dtMaxC,Iter)

    Integer :: i,j,NumNP,Level,NLevel,Iter
    Real(dp) :: dtMaxC,dt,v,ThW,vj,Thj,Peclet,Courant,CourMax
    logical :: lArtD
    Real(dp) :: x(:),Disp(:)
    Real(dp) :: TanH,z,Cour,Cour1,DD,dtMax,dx,Pe2,Pec,RMin,vMax,vv,vv1

    TanH(z)=(exp(z)-exp(-z))/(exp(z)+exp(-z))

    If(debug) Print *,'pecour'
    if(i.ne.NumNP) then
       dx=x(j)-x(i)
       vv=0.
       if(ThW.gt.1.e-6.and.Thj.gt.1.e-6) vv=(abs(v)/ThW+abs(vj)/Thj)/2.
       vv1=0.
       if(ThW.gt.1.e-6.and.Thj.gt.1.e-6) vv1=(v/ThW+vj/Thj)/2.
       DD=(Disp(i)+Disp(j))/2.
       if(Level.eq.NLevel) then
          Pec=99999.
          dtMax=1.e+30
          !          vMax=amax1(abs(v)/ThW,abs(vj)/Thj)
          vMax=(abs(v)+abs(vj))/(ThW+Thj)
          RMin=Min(Retard(i),Retard(j))
          if(DD.gt.0.) Pec=abs(vv)*dx/DD
          Cour=vMax*dt/dx/RMin
          Peclet=Max(Peclet,Pec)
          Courant=Max(Courant,Cour)
          Cour1=CourMax
          if(.not.lUpW.and..not.lArtD) then
             If(Pec.Ne.99999.) Cour1=Min(1.,PeCr/Max(0.5,Pec))
          end if
          if(epsi.lt.1..and.vMax.gt.1.e-20) dtMax=Cour1*dx*RMin/vMax
          !         the von Neumann time step limit
          !          RThE=(ThW+thj)/2.*RMin
          !          if(abs(DD).gt.1.e-20)dtMax=Min(dtMax,10.*RThE*dx*dx/2./DD)
          dtMaxC=Min(dtMaxC,dtMax)

          !       Calculate upstream weighting factors
       else if(lUpW.and.Iter.eq.1) then
          Pe2=11.
          if(DD.gt.0.) Pe2=dx*vv1/DD/2.
          if(abs(vv).lt.1.e-30) then
             wc(i)=0.
          else if(abs(Pe2).gt.10.) then
             if(vv1.gt.0.) wc(i)=1.
             if(vv1.lt.0.) wc(i)=-1
          else
             wc(i)=1./TanH(Pe2)-1./Pe2
             wc(i)=Min( 1.,wc(i))
             wc(i)=Max(-1.,wc(i))
          end if
       end if
    end if

  End Subroutine PeCour

  !************************************************************************

  !     Calculate the dispersion coefficients

  Subroutine Disper(i,M,dt,lArtD,&
       &                  iDualPor,Level,NLevel,Disp,thSat,&
       &                  ThImob,ThW,ThG,v,Dw,Dg,Henry,ro,Frac,xKs,fExp,&
       &                  xNu,cMid,dSConc,TauG,lBact)
    Integer :: i,M,iDualPor,Level,NLevel
    logical :: lArtD,lBact
    Real(dp) :: thSat(:),Disp(:),dt,ThImob,ThW,ThG,v,Dw,Dg,Henry,ro,Frac,xKs,fExp,&
         xNu,cMid,dSConc,TauG,DD,DPom,fi,TauW,ThS

    If(debug) Print *,'disper'
    if(lTort) then
       ThS=thSat(M)
       if(lMobIm(M).and.iDualPor.eq.0.or.lBact) ThS=max(thSat(M)-ThImob,0.001)
       if(              iDualPor.gt.0)          ThS=    thSat(M)+ThImob
       If(iTort.Eq.0) Then
          TauW=ThW**(7./3.)/ThS**2
          TauG=ThG**(7./3.)/ThS**2
       else
          TauW=0.66*(ThW/ThS)**(8./3.)
          TauG=ThG**1.5/ThS
       End If
    else
       TauW=1.
       TauG=1.
    end if
    Disp(i)=ChPar(2,M,1)*Abs(v)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
    if(.not.lArtD.and..not.lUpW) then
       fi=0.
       if(cMid.gt.0.)&
            &    fi=6.*ThW*ro*xKs*cMid**(fExp-1.)*&
            &         (fExp/(1.+xNu*cMid**fExp)**2-1./(1.+xNu*cMid**fExp))
       DPom=Max(dt/(6.*ThW*(ThW+ro*Frac*xKs*dSConc+ThG*Henry)+fi),0.)
       if(Level.ne.NLevel) then
          Disp(i)=Disp(i)+v*v*DPom
       else
          Disp(i)=Max(Disp(i)-v*v*DPom,Disp(i)/2.)
       end if
    end if
    if(lArtD) then
       DD=0.
       if(PeCr.ne.0.and.abs(v).gt.1.e-15) DD=v*v*dt/thW/thW/&
            &       Retard(i)/PeCr
       if(DD.gt.Disp(i)) Disp(i)=DD
    end if

  End Subroutine Disper

  !************************************************************************

  !     Calculate the adsorbed concentration on kinetic sites or
  !     the concentration in an imobile zone, before solving matrix equation

  Subroutine NEquil(i,jS,M,&
       &                  SSorb,SSorb2,lLinear,lBact,&
       &                  Level,NLevel,dt,ro,xKs,xKsO,cc,SConc,&
       &                  SConcO,SConcS,SConcOS,dSConcS,xMuL,xMuLO,xMuS,&
       &                  xMuSO,dRetardS,GamL,GamL1,GamLO,GamL1O,GamS,&
       &                  GamS1,GamSO,GamS1O,Omega,OmegaO,rKa1,rKa1O,rKa2,&
       &                  rKa2O,rKd1,rKd1O,rKd2,rKd2O,ThW,ThWO,psi1,psi1O,&
       &                  psi2,psi2O,DMobI,Frac,ThImob,ThImobO,iDualPor,&
       &                  FlMacro,lNEquil,GamL1Pi,GamS1Pi,xKsP,&
       &                  SConcPS,f_em,OmegaS)
    Integer :: i,jS,M,Level,NLevel,iDualPor
    Logical :: lLinear(:),lBact,lNEquil
    Real(dp) :: SSorb,SSorb2,dt,ro,xKs,xKsO,cc,SConc,&
         SConcO,SConcS,SConcOS,dSConcS,xMuL,xMuLO,xMuS,&
         xMuSO,dRetardS,GamL,GamL1,GamLO,GamL1O,GamS,&
         GamS1,GamSO,GamS1O,Omega,OmegaO,rKa1,rKa1O,rKa2,&
         rKa2O,rKd1,rKd1O,rKd2,rKd2O,ThW,ThWO,psi1,psi1O,&
         psi2,psi2O,DMobI,Frac,ThImob,ThImobO,&
         FlMacro,GamL1Pi,GamS1Pi,xKsP,&
         SConcPS,f_em,OmegaS,AMobI,BMobI,BMobIO,dTheta,EMobI,EMobIO,GMobI0,cS
    If(debug) Print *,'nequil'
    FlMacro=0.
    if(iDualPor.gt.0) then            ! mobile-immobile model
       if(SinkIm(i).gt.0) then
          FlMacro=SinkIm(i)*Conc(jS,i)
       else
          FlMacro=SinkIm(i)*Sorb(jS,i)
       end if
    end if
    SSorb =Sorb (jS,i)
    SSorb2=Sorb2(jS,i)
    if(Level.eq.NLevel) then

       !       mobile-immobile model
       if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then
          AMobI=(ThImob+ThImobO)/2.+(1.-Frac)*ro*xKs*dSConcS
          dTheta=ThImob-ThImobO
          EMobI =ThImob*xMuL+(1.-Frac)*ro*xMuS-(1.-Frac)*ro*dRetardS
          EMobIO=ThImobO*xMuLO+(1.-Frac)*ro*xMuSO-(1.-Frac)*ro*dRetardS
          EMobI =EMobI +FlMacro
          EMobIO=EMobIO+FlMacro
          if(lNEquil) then
             cS=Sorb(jS-1,i)*(ThImob*GamL1Pi+&
                  &                               (1.-Frac)*ro*GamS1Pi*xKsP*SConcPS)
             EMobI =EMobI +cS
             EMobIO=EMobIO+cS
          end if
          BMobI =ThImob *(GamL +GamL1 )+&
               &                       (1.-Frac)*ro*(GamS +GamS1 )*xKs *SConcS
          BMobIO=ThImobO*(GamLO+GamL1O)+&
               &                       (1.-Frac)*ro*(GamSO+GamS1O)*xKsO*SConcOS
          if(lLinear(jS)) then
             DMobI=  2.*AMobI+dt*(Omega +BMobI )+dTheta
             GMobI0=(2.*AMobI-dt*(OmegaO+BMobIO)-dTheta)/DMobI
             Sorb(jS,i)=Sorb(jS,i)*GMobI0+&
                  &                      dt*(OmegaO*Conc(jS,i)+EMobI+EMobIO)/DMobI
             SSorb=Sorb(jS,i)
          else
             SorbN(i)=Sorb(jS,i)+dt/AMobI*  (   epsi *(Omega *&
                  &              (cNew(i)   -SorbN(i))  -BMobI *SorbN(i)  +EMobI)+&
                  &                                     (1.-epsi)*(OmegaO*&
                  &              (Conc(jS,i)-Sorb(jS,i))-BMobIO*Sorb(jS,i)+EMobIO))
             SSorb=SorbN(i)
          end if
          if(lDualNEq) then
             cS=0.
             if(lNEquil) cS=GamS1Pi*Sorb2(jS-1,i)*dt*2.
             if(lLinear(jS)) then
                Sorb2(jS,i)=((2.-(OmegaS+GamSO+GamS1O)*dt)*Sorb2(jS,i)+&
                     &                  dt*Frac*(1.-f_em)*OmegaS*xKsO*Conc(jS,i)+&
                     &                  dt*(1.-f_em)*(xMuSO+xMuS)+cS)/&
                     &                  (2.+dt*(OmegaS+GamS+GamS1))
                SSorb2=Sorb2(jS,i)
             else
                SorbN2(i)=Sorb2(jS,i)+dt*&
                     &    (epsi* (OmegaS*(Frac*(1.-f_em)*SConc *xKs *cc     -SorbN2(i))-&
                     &                         (GamS+GamS1)*SorbN2(i)+(1.-f_em)*xMuS)+&
                     & (1.-epsi)*(OmegaS*(Frac*(1.-f_em)*SConcO*xKsO*Conc(jS,i)-SSorb2)-&
                     &                         (GamSO+GamS1O)*SSorb2+(1.-f_em)*xMuSO))
                SSorb2=SorbN2(i)
             end if
          end if

          !       two-site sorption model
       else if(.not.lBact) then
          cS=0.
          if(lNEquil) cS=GamS1Pi*Sorb(jS-1,i)*dt*2.
          if(lLinear(jS)) then
             Sorb(jS,i)=((2.-(OmegaO+GamSO+GamS1O)*dt)*Sorb(jS,i)+&
                  &                  dt*(1.-Frac)*OmegaO*xKsO*Conc(jS,i)+&
                  &                  dt*(1.-Frac)*(xMuSO+xMuS)+cS)/&
                  &                  (2.+dt*(Omega+GamS+GamS1))
             SSorb=Sorb(jS,i)
          else
             SorbN(i)=Sorb(jS,i)+dt*&
                  &          (epsi* (Omega* ((1.-Frac)*SConc *xKs *cc     -SorbN(i))-&
                  &                         (GamS+GamS1)*SorbN(i)+(1.-Frac)*xMuS)+&
                  &       (1.-epsi)*(OmegaO*((1.-Frac)*SConcO*xKsO*Conc(jS,i)-SSorb)-&
                  &                         (GamSO+GamS1O)*SSorb+(1.-Frac)*xMuSO))
             SSorb=SorbN(i)
          end if

          !       filtration model
       else if(lBact) then
          if(lLinear(jS)) then
             Sorb(jS,i)=((2.-dt*(rKd1O+GamSO+GamS1O))*Sorb(jS,i)+&
                  &                      dt*rKa1O*ThW*Conc(jS,i)/ro)/&
                  &                      (2.+dt*(rKd1+GamS+GamS1))
             SSorb=Sorb(jS,i)
             Sorb2(jS,i)=((2.-dt*(rKd2O+GamSO+GamS1O))*Sorb2(jS,i)+&
                  &                      dt*rKa2O*ThW*Conc(jS,i)/ro)/&
                  &                      (2.+dt*(rKd2+GamS+GamS1))
             SSorb2=Sorb2(jS,i)
          else
             SorbN(i)=Sorb(jS,i)+dt*&
                  &               (epsi*    (rKa1*ThW/ro*psi1*cc-&
                  &                         (rKd1+GamS+GamS1)*SorbN(i))+&
                  &               (1.-epsi)*(rKa1O*ThWO/ro*psi1O*Conc(jS,i)-&
                  &                         (rKd1O+GamSO+GamS1O)*Sorb(jS,i)))
             SSorb=SorbN(i)
             SorbN2(i)=Sorb2(jS,i)+dt*&
                  &               (epsi*    (rKa2*ThW/ro*psi2*cc-&
                  &                         (rKd2+GamS+GamS1)*SorbN2(i))+&
                  &               (1.-epsi)*(rKa2O*ThWO/ro*psi2O*Conc(jS,i)-&
                  &                         (rKd2O+GamSO+GamS1O)*Sorb2(jS,i)))
             SSorb2=SorbN2(i)
          end if
       end if
    end if

  End Subroutine NEquil

  !************************************************************************

  !     Calculate sorbed concentration for linear noneq. adsorption or 
  !     concentration in the imobile water. 
  !     At the end of the time step. After solving matrix equation

  Subroutine SorbConc(jS,MatNum,Temp,dt,lBact,ThW,Veloc,&
       iDualPor,iMoistDep,xConv,tConv)

    Logical lBact
    Integer :: MatNum(:),jS,iDualPor,iMoistDep
    Real(dp), Dimension(:) :: Temp,ThW,Veloc
    Real(dp) :: dt,xConv,tConv
    Integer :: i,N,M
    Real(dp) :: Tr,R,TT,Frac,xKs,f1,Alfa1,Alfa2,AMobI,BMobI,Dc,DMobI,dTheta,f_em,FlMacro,&
         GamL,GamL1,GamS,GamS1,Omega,rKa1,rKa2,rKd1,rKd2,ro,Theta,ThImob,ThImobO,Dep

    If(debug) Print *,'sorbconc'
    ThImobO=0.0
    Tr=293.15
    R=8.314
    ThImob=0.0
    N=Ubound(MatNum,1)
    do i=1,N
       M=MatNum(i)
       TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
       Frac  =ChPar( 3,M,jS)*Exp(TDep(3,jS)*TT)
       xKs   =ChPar( 7,M,jS)*Exp(TDep(7,jS)*TT)
       f1=1.
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,11,2,ThW(i),iMoistDep)
       GamS  =ChPar(12,M,jS)*Exp(TDep(12,jS)*TT)*f1
       if(iMoistDep.gt.0) &
            &    f1=rMD(M,jS,13,5,ThW(i),iMoistDep)
       GamS1 =ChPar(15,M,jS)*Exp(TDep(15,jS)*TT)*f1
       Omega =ChPar(20,M,jS)*Exp(TDep(20,jS)*TT)

       if((lMobIm(M).or.iDualPor.gt.0).and..not.lBact) then ! mobile-immobile model
          ro    =ChPar(1,M,jS)*Exp(TDep(1,jS)*TT)
          if(iDualPor.eq.0) then
             ThImob=ChPar(4,M,jS)*Exp(TDep(4,jS)*TT)
             ThImobO=ThImob
          else if(iDualPor.gt.0) then
             ThImob =ThNewIm(i)
             ThImobO=ThOldIm(i)
          end if
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,10,1,ThImob,iMoistDep)
          GamL  =ChPar(11,M,jS)*Exp(TDep(11,jS)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,12,4,ThImob,iMoistDep)
          GamL1 =ChPar(14,M,jS)*Exp(TDep(14,jS)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,11,2,ThImob,iMoistDep)
          GamS  =ChPar(12,M,jS)*Exp(TDep(12,jS)*TT)*f1
          if(iMoistDep.gt.0) &
               &      f1=rMD(M,jS,13,5,ThImob,iMoistDep)
          GamS1 =ChPar(15,M,jS)*Exp(TDep(15,jS)*TT)*f1
          dTheta=ThImob-ThImobO
          AMobI=(ThImob+ThImobO)/2.+(1.-Frac)*ro*xKs
          BMobI=ThImob*(GamL+GamL1)+(GamS+GamS1)*ro*(1.-Frac)*xKs
          DMobI=2.*AMobI+dt*(Omega+BMobI)+dTheta
          Sorb(jS,i)=Sorb(jS,i)+dt*Omega*Conc(jS,i)/DMobI
          FlMacro=0.
          if(iDualPor.gt.0) then
             if(SinkIm(i).gt.0) then
                FlMacro=SinkIm(i)*Conc(jS,i)
             else
                FlMacro=SinkIm(i)*Sorb(jS,i)
             end if
          end if
          if(jS.eq.1) STrans(i)=Omega*(Conc(jS,i)-Sorb(jS,i))+FlMacro
          if(lDualNEq) then
             f_em  =ChPar(13,M,jS)*Exp(TDep(13,jS)*TT)
             OmegaS=ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)
             Sorb2(jS,i)=&
                  &              Sorb2(jS,i)+dt*OmegaS*Frac*(1.-f_em)*xKs*Conc(jS,i)/&
                  &                         (2.+dt*(OmegaS+GamS+GamS1))
          end if

       else if(.not.lBact) then                 ! two-site sorption model
          Sorb(jS,i)=Sorb(jS,i)+dt*Omega*(1.-Frac)*xKs*Conc(jS,i)/&
               &                       (2.+dt*(Omega+GamS+GamS1))

       else if(lBact) then                      ! filtration model
          ro    =ChPar(1,M,jS)*Exp(TDep(1,jS)     *TT)
          ThImob=ChPar(4,M,jS)
          Theta=ThW(i)-ThImob
          GamS  =ChPar(12,M,jS)*Exp(TDep(12,jS)*TT)
          GamS1 =0.
          rKa1  =ChPar(19,M,jS)*Exp(TDep(19,jS)*TT)
          rKd1  =ChPar(20,M,jS)*Exp(TDep(20,jS)*TT)
          rKa2  =ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)
          rKd2  =ChPar(17,M,jS)*Exp(TDep(17,jS)*TT)
          if(lFiltr) then
             Dc=ChPar(13,M,jS)*Exp(TDep(13,jS)*TT)
             Dep=ChPar(14,M,jS)*Exp(TDep(14,jS)*TT)
             Alfa1=rKa1
             Alfa2=rKa2
             call Deposit(rKa1,rKa2,Dc,Dep,Alfa1,Alfa2,Theta,Veloc(i),&
                  &                   Temp(i),xConv,tConv)
          end if
          Sorb(jS,i) =Sorb(jS,i) +dt*rKa1*Theta*Conc(jS,i)/ro/&
               &                          (2.+dt*(rKd1+GamS+GamS1))
          Sorb2(jS,i)=Sorb2(jS,i)+dt*rKa2*Theta*Conc(jS,i)/ro/&
               &                          (2.+dt*(rKd2+GamS+GamS1))
       end if
    end do

  End Subroutine SorbConc

  !************************************************************************

  !     Calculate mass-transfer fluxes at the end of the time interval

  Subroutine MassTran(jS,MatNum,Temp,lEquil,&
       &                    x,epsi,&
       &                    sSink,lBact,theta,&
       &                    Veloc,iDualPor,xConv,tConv,&
       &                    lLinear)
    Use Geometry, Only: NumNP
    Logical :: lEquil,lBact,lLinear(:)
    Integer :: MatNum(:),jS,iDualPor
    Real(dp), dimension(:) :: x,Temp,sSink,theta,Veloc
    Real(dp) :: epsi
    Real(dp) :: xConv,tConv,Tr,R,TTi,TTj,dx,Omegai,Omegaj,FlMacroi,FlMacroj,&
         aa,Alfa1,Alfa2,cci,ccj,Dc,f_emi,f_emj,fExpi,fExpj,fraci,fracj,&
         OmegaSi,OmegaSj,psi1i,psi1j,psi2i,psi2j,rKa1,rKa2,rka1i,rka1j,rKa2i,rKa2j,&
         rKd1i,rKd1j,rkd2i,rkd2j,roi,roj,SMax1i,SMax1j,SMax2i,SMax2j,SorbEi,SorbEj,&
         Thetai,Thetaj,thimobi,thimobj,xksi,xksj,xnui,xnuj,Dep
    Integer :: i,Mi,j,Mj,iPsi1,iPsi2

    If(debug) Print *,'masstran'
    Tr=293.15
    R=8.314
    Do i=1,NumNP
       Mi=MatNum(i)
       TTi=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
       If(i.Eq.NumNP) Cycle
       j=i+1
       dx=x(j)-x(i)
       cvCh0(jS)=cvCh0(jS)+epsi*dx*(q0(i)+q0(j))/2.
       cvCh1(jS)=cvCh1(jS)+epsi*dx*(q1(i)+q1(j))/2.
       cvChR(jS)=cvChR(jS)+epsi*dx*(sSink(i)+sSink(j))/2.
       if(.not.lEquil) then
          Mj=MatNum(j)
          TTj=(Temp(j)+273.15-Tr)/R/(Temp(j)+273.15)/Tr
          Omegai=ChPar(20,Mi,jS)*Exp(TDep(20,jS)*TTi)
          Omegaj=ChPar(20,Mj,jS)*Exp(TDep(20,jS)*TTj)

          !         mobile-immobile model
          if((lMobIm(Mi).or.iDualPor.gt.0).and..not.lBact) then
             FlMacroi=0.
             FlMacroj=0.
             if(iDualPor.gt.0) then
                if(SinkIm(i).gt.0) then
                   FlMacroi=SinkIm(i)*Conc(jS,i)
                else
                   FlMacroi=SinkIm(i)*Sorb(jS,i)
                end if
                if(SinkIm(j).gt.0) then
                   FlMacroj=SinkIm(j)*Conc(jS,j)
                else
                   FlMacroj=SinkIm(j)*Sorb(jS,j)
                end if
             end if
             cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*&
                  &                  (Omegai*(Conc(jS,i)-Sorb(jS,i))+&
                  &                   Omegaj*(Conc(jS,j)-Sorb(jS,j))+&
                  &                   FlMacroi+FlMacroj)
             if(jS.eq.1.and..not.lLinear(jS)) &
                  &          STrans(i)=epsi*(Omegai*(Conc(jS,i)-Sorb(jS,i))+&
                  &                          FlMacroi)
             if(lDualNEq) then
                roi   =ChPar( 1,Mi,jS)*Exp(TDep( 1,jS)*TTi)
                Fraci =ChPar( 3,Mi,jS)*Exp(TDep( 3,jS)*TTi)
                xKsi  =ChPar( 7,Mi,jS)*Exp(TDep( 7,jS)*TTi)
                xNui  =ChPar( 8,Mi,jS)*Exp(TDep( 8,jS)*TTi)
                fExpi =ChPar( 9,Mi,jS) !*exp(TDep( 9,jS)*TTi)
                roj   =ChPar( 1,Mj,jS)*Exp(TDep( 1,jS)*TTj)
                Fracj =ChPar( 3,Mj,jS)*Exp(TDep( 3,jS)*TTj)
                xKsj  =ChPar( 7,Mj,jS)*Exp(TDep( 7,jS)*TTj)
                xNuj  =ChPar( 8,Mj,jS)*Exp(TDep( 8,jS)*TTj)
                fExpj =ChPar( 9,Mj,jS) !*exp(TDep( 9,jS)*TTj)
                f_emi =ChPar(13,Mi,jS)*Exp(TDep(13,jS)*TTi)
                f_emj =ChPar(13,Mj,jS)*Exp(TDep(13,jS)*TTj)
                OmegaSi=ChPar(16,Mi,jS)*Exp(TDep(16,jS)*TTi)
                OmegaSj=ChPar(16,Mj,jS)*Exp(TDep(16,jS)*TTj)
                cci   =Conc(jS,i)
                ccj   =Conc(jS,j)
                SorbEi=0.
                SorbEj=0.
                if(cci.gt.0.) SorbEi=&
                     &             Fraci*(1.-f_emi)*xKsi*cci**fExpi/(1.+xNui*cci**fExpi)
                if(ccj.gt.0.) SorbEj=&
                     &             Fracj*(1.-f_emj)*xKsj*ccj**fExpj/(1.+xNuj*ccj**fExpj)
                cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*&
                     &                      (roi*OmegaSi*(SorbEi-Sorb2(jS,i))+&
                     &                       roj*OmegaSj*(SorbEj-Sorb2(jS,j)))
                if(jS.eq.1.and..not.lLinear(jS)) &
                     &         STrans(i)=STrans(i)+epsi*roi*OmegaSi*(SorbEi-Sorb2(jS,i))
             end if

             !         two-site sorption model
          else if(.not.lBact) then
             roi   =ChPar(1,Mi,jS)*Exp(TDep(1,jS)*TTi)
             Fraci =ChPar(3,Mi,jS)*Exp(TDep(3,jS)*TTi)
             xKsi  =ChPar(7,Mi,jS)*Exp(TDep(7,jS)*TTi)
             xNui  =ChPar(8,Mi,jS)*Exp(TDep(8,jS)*TTi)
             fExpi =ChPar(9,Mi,jS) !*exp(TDep(9,jS)*TTi)
             roj   =ChPar(1,Mj,jS)*Exp(TDep(1,jS)*TTj)
             Fracj =ChPar(3,Mj,jS)*Exp(TDep(3,jS)*TTj)
             xKsj  =ChPar(7,Mj,jS)*Exp(TDep(7,jS)*TTj)
             xNuj  =ChPar(8,Mj,jS)*Exp(TDep(8,jS)*TTj)
             fExpj =ChPar(9,Mj,jS) !*exp(TDep(9,jS)*TTj)
             cci   =Conc(jS,i)
             ccj   =Conc(jS,j)
             SorbEi=0.
             SorbEj=0.
             if(cci.gt.0.)&
                  &          SorbEi=(1.-Fraci)*xKsi*cci**fExpi/(1.+xNui*cci**fExpi)
             if(ccj.gt.0.)&
                  &          SorbEj=(1.-Fracj)*xKsj*ccj**fExpj/(1.+xNuj*ccj**fExpj)
             cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*&
                  &                      (roi*Omegai*(SorbEi-Sorb(jS,i))+&
                  &                       roj*Omegaj*(SorbEj-Sorb(jS,j)))
             if(jS.eq.1) &
                  &        STrans(i)=epsi*roi*Omegai*(SorbEi-Sorb(jS,i))

             !         filtration model
          else if(lBact) then
             roi   =ChPar(1, Mi,jS)*Exp(TDep( 1,jS)*TTi)
             roj   =ChPar(1, Mj,jS)*Exp(TDep( 1,jS)*TTj)
             SMax1i=ChPar(18,Mi,jS)*Exp(TDep(18,jS)*TTi)
             SMax1j=ChPar(18,Mj,jS)*Exp(TDep(18,jS)*TTj)
             rKa1i =ChPar(19,Mi,jS)*Exp(TDep(19,jS)*TTi)
             rKa1j =ChPar(19,Mj,jS)*Exp(TDep(19,jS)*TTj)
             rKd1i =ChPar(20,Mi,jS)*Exp(TDep(20,jS)*TTi)
             rKd1j =ChPar(20,Mj,jS)*Exp(TDep(20,jS)*TTj)
             SMax2i=ChPar(15,Mi,jS)*Exp(TDep(15,jS)*TTi)
             SMax2j=ChPar(15,Mj,jS)*Exp(TDep(15,jS)*TTj)
             rKa2i =ChPar(16,Mi,jS)*Exp(TDep(16,jS)*TTi)
             rKa2j =ChPar(16,Mj,jS)*Exp(TDep(16,jS)*TTj)
             rKd2i =ChPar(17,Mi,jS)*Exp(TDep(17,jS)*TTi)
             rKd2j =ChPar(17,Mj,jS)*Exp(TDep(17,jS)*TTj)
             ThImobi=ChPar(4,Mi,jS)
             ThImobj=ChPar(4,Mj,jS)
             Thetai=theta(i)-ThImobi
             Thetaj=theta(j)-ThImobj
             iPsi1=0
             iPsi2=0
             if(.not.lFiltr) iPsi2=int(ChPar(13,Mi,jS))
             if(.not.lFiltr) iPsi1=int(ChPar(14,Mj,jS))
             if(iPsi1.eq.0.and.SMax1i.gt.0.) iPsi1=1
             if(iPsi2.eq.0.and.SMax2i.gt.0.) iPsi2=1
             if(iPsi1.ge.3.or.iPsi2.ge.3) &
                  Dc=ChPar(6,Mi,jS)*Exp(TDep(6,jS)*TTi)
             if(iPsi1.eq.5.or.iPsi2.eq.5) aa=ChPar(15,Mi,jS)
             psi1i=1.
             psi1j=1.
             psi2i=1.
             psi2j=1.
             if(iPsi1.gt.0) then
                call Blocking(iPsi1,SMax1i,psi1i,x(i),Sorb(jS,i),Dc,aa)
                call Blocking(iPsi1,SMax1j,psi1j,x(j),Sorb(jS,j),Dc,aa)
             end if
             if(iPsi2.gt.0) then
                call Blocking(iPsi2,SMax2i,psi2i,x(i),Sorb2(jS,i),Dc,aa)
                call Blocking(iPsi2,SMax2j,psi2j,x(j),Sorb2(jS,j),Dc,aa)
             end if
             if(lFiltr) then
                Dc=ChPar(13,Mi,jS)*Exp(TDep(13,jS)*TTi)
                Dep=ChPar(14,Mi,jS)*Exp(TDep(14,jS)*TTi)
                Alfa1=rKa1i
                Alfa2=rKa2i
                call Deposit(rKa1,rKa2,Dc,Dep,Alfa1,Alfa2,Thetai,&
                     &                     Veloc(i),Temp(i),xConv,tConv)
                rKa1i=rKa1
                rKa1j=rKa1
                rKa2i=rKa2
                rKa2j=rKa2
             end if
             cvChIm(jS)=cvChIm(jS)+epsi*dx/2.*&
                  &            (Conc(jS,i)*Thetai*(psi1i*rKa1i+psi2i*rKa2i)+&
                  &             Conc(jS,j)*Thetaj*(psi1j*rKa1j+psi2j*rKa2j)-&
                  &             roi*(Sorb(jS,i)*rKd1i+Sorb2(jS,i)*rKd2i)-&
                  &             roj*(Sorb(jS,j)*rKd1j+Sorb2(jS,j)*rKd2j))
             if(jS.eq.1) STrans(i)=epsi*&
                  &                     (Conc(jS,i)*Thetai*(psi1i*rKa1i+psi2i*rKa2i)-&
                  &                      roi*(Sorb(jS,i)*rKd1i+Sorb2(jS,i)*rKd2i))
          end if
       end if
    End Do
  End Subroutine MassTran

  !************************************************************************

  !     Calculate flux concentration for the first solute

  Subroutine FluxConc(x,v,theta,thSat,MatNum,Temp,iDualPor,ThIm,jS)
    Integer :: MatNum(:),iDualPor,jS
    Real(dp) :: x(:),v(:),theta(:),thSat(:),Temp(:),ThIm(:)
    Integer :: i,M,NumNP
    Real(dp) :: Tr,R,ThW,ThG,cGrad,Dg,Disp,Dw,Henry,qW,TauG,TauW,ThImob,ThS,TT
    If(debug) Print *,'fluxconc'
    Tr=293.15
    R=8.314
    ThImob=0.0
    NumNP=Ubound(MatNum,1)
    Do i=1,NumNP
       M=MatNum(i)
       ThW=Theta(i)
       ThG=Max(0.,thSat(M)-ThW)
       if(lMobIm(M)) then
          if(iDualPor.eq.0) then
             ThImob=ChPar(4,M,jS)
             ThW=max(ThW-ThImob,0.001)
          else if(iDualPor.gt.0) then
             ThImob=ThIm(i)
          end if
       end if
       TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
       Dw   =ChPar( 5,M,jS)*Exp(TDep( 5,jS)*TT)
       Dg   =ChPar( 6,M,jS)*Exp(TDep( 6,jS)*TT)
       Henry=ChPar(10,M,jS)*Exp(TDep(10,jS)*TT)
       if(lTort) then
          ThS=thSat(M)
          if(lMobIm(M).and.iDualPor.eq.0) ThS=max(thSat(M)-ThImob,0.001)
          if(              iDualPor.gt.0) ThS=    thSat(M)+ThImob
          TauW=ThW**(7./3.)/ThS**2
          TauG=ThG**(7./3.)/ThS**2
       else
          TauW=1.
          TauG=1.
       end if
       qW=v(i)
       Disp=ChPar(2,M,jS)*abs(qW)/ThW+Dw*TauW+ThG/ThW*Dg*Henry*TauG
       cGrad=0.
       if(i.eq.1) then
          cGrad=(Conc(jS,i+1)-Conc(jS,i))/(x(i+1)-x(i))
       else if(i.eq.NumNP) then
          cGrad=(Conc(jS,i)-Conc(jS,i-1))/(x(i)-x(i-1))
       else
          cGrad=(Conc(jS,i+1)-Conc(jS,i-1))/(x(i+1)-x(i-1))
       end if
       cNew(i)=Conc(jS,i)
       if(qW.ne.0.) cNew(i)=Conc(jS,i)-Disp*ThW/qW*cGrad
    end do

  End Subroutine FluxConc

  !************************************************************************

  !     Calculate blocking coefficient for the attachment process

  Subroutine Blocking(iPsi,SMax,psi,x,ss,Dc,SMax2)

    integer :: iPsi
    Real(dp) :: SMax,psi,x,ss,Dc,SMax2
    Real(dp) :: Minf,Binf,Sinf,const

    If(debug) Print *,'blocking'
    psi=1.
    if(iPsi.eq.1) then
       if(SMax.gt.0.) psi=1.-ss/SMax
    else if(iPsi.eq.2) then
       if(SMax.gt.0.) psi=max(ss**SMax,psi)
    else if(iPsi.eq.3) then
       Binf=1./SMax
       Sinf=.546
       Minf=Dc
       const=Sinf*Binf*ss  
       if(ss.le.(0.8*SMax))&
            &    psi=1.-(4.*const)+(3.08*const**2.)+(1.4069*const**3.)
       if(ss.gt.(0.8*SMax))&
            &    psi=((1.-Binf*ss)**3.)/(2.*(Minf**2.)*(Binf**3.))
    else if(iPsi.eq.4) then
       if(SMax.gt.0..and.Dc.gt.0.) psi=((abs(x)+Dc)/Dc)**(-SMax)
    else if(iPsi.eq.5) then
       if(SMax.gt.0..and.Dc.gt.0.) psi=((abs(x)+Dc)/Dc)**(-SMax)
       if(SMax2.gt.0.) psi=psi*(1.-ss/SMax2)
    end if

  End Subroutine Blocking

  !************************************************************************

  !     Calculate the deposition coefficient for the bacteria transport,
  !     All calculations within this subroutines are in meters and seconds
  !     Conversions are needed

  Subroutine Deposit(Ka1,Ka2,Dc1,Dp1,Alfa1,Alfa2,Theta,q,Temp,xConv,tConv)
    Real(dp) :: Dc1,Dp1,Alfa1,Alfa2,Theta,q,Temp,xConv,tConv
    Real(dp) :: Ka1,Ka2,mu,N_Pe,N_Lo,N_R,N_G,As,Bk,Dc,e_diff,e_grav,e_inter
    Real(dp) :: eta,g,gamma,H,PVeloc,rof,rop,Veloc,PI,Dep

    If(debug) Print *,'deposit'
    !     Ka       - deposition coefficient (output) [1/T]
    !     Dc       - diameter of the sand grains (m)
    !     Dp       - diameter of the bacteria (0.95 microm) (m)
    !     Alfa     - sticking efficiency (-)
    !     Theta    - porosity (-)
    !     q        - Darcy flux [L/T]
    !     Temp     - Temperature in Celcius
    PI=3.1415                ! Ludolf's number
    mu=0.00093               ! fluid viscosity (Pa s)
    Bk=1.38048e-23           ! Boltzman constatnt (J/K)
    H=1.e-20                 ! Hamaker constant (J)
    g=9.81                   ! gravitational acceleration (m/s2)
    rop=1080.                ! bacterial density (kg/m3)
    rof=998.                 ! fluid density (kg/m3)
    Veloc=abs(q/xConv*tConv) ! absolute value of Darcy flux (converted to m/s)
    PVeloc=Veloc/Theta       ! pore velocity (converted to m/s)
    Dc=Dc1/xConv             ! conversion to m
    Dep=Dp1/xConv             ! conversion to m
    If(Dep.Le.0.And.Dc.Le.0.) Call WriteError('Both Dp and Dc are equal to zero !!!')

    if(Veloc.gt.0.) then
       gamma=(1.-Theta)**(1./3.)
       As=2.*(1.-gamma**5)/(2.-3.*gamma+3.*gamma**5-2.*gamma**6) ! Correct.factor
       N_Pe=3.*PI*mu*Dep*Dc*Veloc/(Bk*(Temp+273.15)) ! Peclet number
       e_diff=4.*As**(1./3.)*N_Pe**(-2./3.)         ! removal by diffusion

       N_Lo=4.*H/(9.*PI*mu*Dep**2*Veloc)             ! London number
       N_R=Dep/Dc                                    ! Interception number
       e_inter=As*N_Lo**(1./8.)*N_R**(15./8.)       ! removal interception

       N_G=g*(rop-rof)*Dep**2/(18.*mu*Veloc)         ! Gravitation number
       e_grav=0.00338*As*N_G**1.2*N_R**(-0.4)       ! removal by gravitational
       !                                                      sedimentation
    else
       e_diff =0.
       e_inter=0.
       e_grav =0.
    end if
    eta=e_diff+e_inter+e_grav                      ! single-collector efficiency

    !     Original Filtration Theory
    Ka1=3.*(1.-Theta)/2./dc*eta*Alfa1*PVeloc
    Ka1=Ka1/tConv
    Ka2=3.*(1.-Theta)/2./dc*eta*Alfa2*PVeloc
    Ka2=Ka2/tConv

  End Subroutine Deposit

  !************************************************************************

  !     Nonequilibrium phase is initially in equilibrium with liquid phase

  Subroutine NonEqInit(MatNum,Temp,lLinear,iDualPor,lBact,Theta)

    Logical :: lLinear(:),lBact
    Integer :: MatNum(:),iDualPor
    Real(dp) :: Temp(:),Theta(:)
    Integer :: jS,i,M
    Real(dp) :: Tr,R,TT,ro,Frac,xKs,xNu,SConc,cc,rKa1,rKd1,xKs1,rKa2,rKd2,xKs2,fExp
    If(debug) Print *,'noneqinit'
    Tr=293.15
    R=8.314
    Do jS=1,Ubound(Sorb,1)
       Do i=1,Ubound(Sorb,2)
          M=MatNum(i)
          TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
          if(lMobIm(M).or.iDualPor.gt.0) then
             Sorb(jS,i)=Conc(jS,i)
          else
             ro   =ChPar(1,M,jS)*Exp(TDep(1,jS)*TT)
             Frac =ChPar(3,M,jS)*Exp(TDep(3,jS)*TT)
             xKs  =ChPar(7,M,jS)*Exp(TDep(7,jS)*TT)
             xNu  =ChPar(8,M,jS)*Exp(TDep(8,jS)*TT)
             fExp =ChPar(9,M,jS) !*exp(TDep(9,jS)*TT)
             SConc=1.
             cc=Conc(jS,i)
             if(.not.lLinear(jS).and.cc.gt.0.) &
                  &        SConc=cc**(fExp-1.)/(1.+xNu*cc**fExp)
             Sorb(jS,i)=(1.-Frac)*SConc*xKs*cc
             if(lBact) then
                Frac=0.
                rKa1=ChPar(19,M,jS)*Exp(TDep(19,jS)*TT)
                rKd1=ChPar(20,M,jS)*Exp(TDep(20,jS)*TT)
                xKs1=Theta(i)*rKa1/ro/rKd1 
                Sorb(jS,i)=xKs1*Conc(jS,i)
                rKa2=ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)
                rKd2=ChPar(17,M,jS)*Exp(TDep(17,jS)*TT)
                xKs2=Theta(i)*rKa1/ro/rKd1 
                Sorb2(jS,i)=xKs2*Conc(jS,i)
             end if
          end if
       End Do
    End Do
  End Subroutine NonEqInit

  !***********************************************************************

  !     Subroutine calculating accesible water content for colloids and colloid velocity
  !     Only for steady-state water flow and homogeneous soil profile

  Subroutine Exclusion(Par,ThNew,vNew,ThOld,vOld)

    Real(dp) :: thr,ths,swr,alpha,vgn1,vgn2,vgm1,vgm2,Pc_c,&
         &      sw_c_eff,sw_c,sw,sw_eff,kr_c,krw,r_c
    Real(dp) :: Par(:,:),ThNew(:),vNew(:),ThOld(:),vOld(:)
    Integer :: i,M
    Real(dp) :: PorVelSolute,PorVelColloid
    Real(dp) :: th_c,xL,ThC,VelC

    If(debug) Print *,'exclusion'
    M=1
    th_c=ChPar(4,M,1)   ! Water Content from which colloids are excluded
    if(abs(th_c).lt.1.e-20) return
    ChPar(4,M,:)=0.
    thr=Par(1,M)
    ths=Par(2,M)
    swr=thr/ths
    sw_c=th_c/ths
    alpha=Par(3,M)
    vgn1=Par(4,M)
    xL=Par(6,M)
    vgm1=1.0-1.0/vgn1
    vgn2=vgn1+1
    vgm2=1.0-2.0/vgn2
    if(sw_c.lt.0.or.sw_c.gt.1.) then
       write(*,*) 'Problem with size exclusion!'
       write(*,*) 'Press Enter to continue'
       read(*,*)
    end if

    !     Accessible Water Content to Colloid
    ThC=ThNew(1)-ths*sw_c

    sw=ThNew(1)/ths
    sw_eff=(sw-swr)/(1.0-swr)
    sw_c_eff=(sw_c-swr)/(1.0-swr)

    !     Colloid Permeability according to Burdine Model
    if(sw_eff.gt.sw_c_eff) then
       kr_c=(sw_eff**2.0)*(((1.0-sw_c_eff**(1.0/vgm2))**vgm2)-&
            &                      ((1.0-sw_eff  **(1.0/vgm2))**vgm2)) 
    else
       kr_c=0.0
    end if

    !     Mualem Water Relative Permeability
    if(sw_eff.gt.0.0) then
       krw=(sw_eff**xL)*(1-(1-sw_eff**(1/vgm1))**vgm1)**2.0
       VelC=vNew(1)*kr_c/krw
    else
       krw=0.0
       VelC=vNew(1)
    end if

    !     Convert sw_c to r_c
    Pc_c=(1.0/alpha)*(sw_c_eff**(-1./vgm1)-1)**(1.0/vgn1) 
    r_c=(2.0*72.0)/(Pc_c*981.0)
    PorVelSolute=vNew(1)/ThNew(1)
    PorVelColloid=VelC/ThC
    !     write(*,*) "r_c microns", 10000.*r_c

    Do i=1,Size(ThNew,1)
       ThNew(i)=ThC
       ThOld(i)=ThC
       vNew(i)=VelC
       vOld(i)=VelC
    end do
  End Subroutine Exclusion

  !***********************************************************************

  !     Reads parameters for the function expressing reaction rate dependence 
  !     on the water content

  Subroutine MoistDepIn(cDataPath,cFileName,iMoistDep)
    Integer :: iMoistDep
    Character :: cFileName*260,cDataPath*260
    Integer :: i,iLengthPath,M,jReact,jS
    If(debug) Print *,'moistdepin'
    iLengthPath = Len_Trim(cDataPath)
    cFileName = cDataPath(1:iLengthPath)//'\MoistDep.in'
    open(15,file=cFileName, status='unknown',err=901)

    read(15,*,err=902)
    Do M=1,Ubound(DMoist,1)
       read(15,*,err=902)
       do jS=1,Ubound(DMoist,2)
          read(15,*,err=902)
          Do jReact=1,Ubound(DMoist,3)
             read(15,*,err=902) (DMoist(M,jS,jReact,i),i=1,6)
          end do
       end do
    end do
    close(15)
    return

    !     Error opening an input file 
901 if(iMoistDep.eq.2) iMoistDep=0
    return
    !     Error reading from an input file 
902 if(iMoistDep.eq.2) iMoistDep=0
    write(*,*) 'Error reading from an input file MoistDep.in !!!!'
    close(15)
  End Subroutine MoistDepIn

  !***********************************************************************

  Real(dp) Function rMD(M,jS,jReact,iReact,Theta,iMoistDep)
    Integer :: M,jS,jReact,iReact,iMoistDep
    Real(dp) :: Theta
    Real(dp) :: ReacMin0,Theta0,Theta1,Theta2,Theta3,ReacMin1
    !     Function expressing reaction rate dependence on the water content
    !     ReacMin0  - relative minimum rate of reaction at low water contents
    !     Theta0    - water content at which reaction rate start increasing
    !     Theta1    - water content at which reaction rate stops increasing
    !     Theta2    - water content at which reaction rate start decreasing  
    !     Theta3    - water content at which reaction rate stops decreasing 
    !     ReacMin1  - relative minimum rate of reaction at high water contents
    !     If theta2=theta3=thetaS -> Anaerobic process
    !     If theta0=theta1=0      -> Aerobic process
    !     If theta2=0 -> no reduction

    rMD=1.
    if(iMoistDep.eq.2) then
       if(jReact.eq.0) return
       ReacMin0=DMoist(M,jS,jReact,1)
       Theta0  =DMoist(M,jS,jReact,2)
       Theta1  =DMoist(M,jS,jReact,3)
       Theta2  =DMoist(M,jS,jReact,4)
       Theta3  =DMoist(M,jS,jReact,5)
       ReacMin1=DMoist(M,jS,jReact,6)
       if(abs(Theta2).lt.0.001) return
       if     (Theta.le.Theta0) then
          rMD=ReacMin0
       else if(Theta.le.Theta1) then
          rMD=ReacMin0+(Theta-Theta0)/(Theta1-Theta0)*(1.-ReacMin0)
       else if(Theta.le.Theta2) then
          rMD=1.
       else if(Theta.le.Theta3) then
          rMD=ReacMin1+(Theta-Theta3)/(Theta2-Theta3)*(1.-ReacMin1)
       else
          rMD=ReacMin1
       end if
    else if(iMoistDep.eq.1) then ! Walker's formula
       If(WDep(iReact,2+M,jS).Gt.Theta.And.WDep(iReact,2+M,jS).Gt.0.) &
            rMD=(Theta/WDep(iReact,2+M,jS))**WDep(iReact,1,jS)
    end if

  End Function rMD

  !***********************************************************************

  Real(dp) Function rMD1(M,jS,jReact,Theta)
    Integer :: M,jS,jReact
    Real(dp) :: Theta
    Real(dp) :: rType,Theta0,Theta1,ReacMin
    !     Function expressing reaction rate dependence on the water content
    !     Type  =1 or -1: Rate increases or decreases with water content, respectively)
    !     Theta0 - water content at which reaction rate start increasing or decreasing  
    !     Theta1 - water content at which reaction rate stops increasing or decreasing 
    !     ReacMin  - relative minimum rate of reaction


    rType  =DMoist(M,jS,jReact,1)
    Theta0 =DMoist(M,jS,jReact,2)
    Theta1 =DMoist(M,jS,jReact,3)
    ReacMin=DMoist(M,jS,jReact,4)
    rMD1=1.
    if(abs(rType).lt.0.1) return
    if(rType.gt.0) then                 ! increasing rate
       if(Theta.ge.Theta1) then
          rMD1=1.
       else if(Theta.le.Theta0) then
          rMD1=ReacMin
       else
          if(abs(Theta1-Theta0).gt.0.)&
               &      rMD1=ReacMin+(Theta-Theta0)/(Theta1-Theta0)*(1.-ReacMin)
       end if
    else                                ! decreasing rate
       if(Theta.ge.Theta1) then
          rMD1=ReacMin
       else if(Theta.le.Theta0) then
          rMD1=1.
       else
          if(abs(Theta1-Theta0).gt.0.)&
               &      rMD1=ReacMin+(Theta-Theta1)/(Theta0-Theta1)*(1.-ReacMin)
       end if
    end if

  End Function rMD1

  !************************************************************************

  !     Distribute mass into different phases

  Subroutine MassInit(NumNP,MatNum,Temp,Theta,ThetaIm,ThSat,lLinear,lBact)

    logical :: lLinear(:),lBact
    Integer :: NumNP,MatNum(:)
    Real(dp) :: Temp(:),Theta(:),ThSat(:),ThetaIm(:)
    Integer :: jS,i,M
    Real(dp) :: Tr,R,Par(10),rKa1,rKa2,rKd1,rKd2,TT,xKs1,xKs2
    If(debug) Print *,'massinit'
    Tr=293.15
    R=8.314
    Do jS=1,NS
       Do i=1,NumNP
          if(Conc(jS,i).gt.0.) then
             M=MatNum(i)
             TT=(Temp(i)+273.15-Tr)/R/(Temp(i)+273.15)/Tr
             Par(1)=ChPar(1,M,jS)*Exp(TDep(1,jS)*TT) !ro
             Par(2)=ChPar(3,M,jS)*Exp(TDep(3,jS)*TT) !frac
             Par(3)=ChPar(7,M,jS)*Exp(TDep(7,jS)*TT) !xKs
             Par(4)=ChPar(8,M,jS)*Exp(TDep(8,jS)*TT) !xNu
             Par(5)=ChPar(9,M,jS) !*exp(TDep(9,jS)*TT) !fExp
             Par(6)=ChPar(10,M,jS)*Exp(TDep(7,jS)*TT) !xKH
             Par(7)=Theta(i)+ThetaIm(i)
             Par(8)=ThSat(M)-Par(7)
             if(lBact) then
                rKa1=ChPar(19,M,jS)*Exp(TDep(19,jS)*TT)
                rKd1=ChPar(20,M,jS)*Exp(TDep(20,jS)*TT)
                xKs1=Theta(i)*rKa1/Par(1)/rKd1 
                rKa2=ChPar(16,M,jS)*Exp(TDep(16,jS)*TT)
                rKd2=ChPar(17,M,jS)*Exp(TDep(17,jS)*TT)
                xKs2=Theta(i)*rKa1/Par(1)/rKd1
                Par(3)=xKs1+xKs2 
             end if
             if(lLinear(jS)) then
                Conc(jS,i)=Conc(jS,i)/(Par(7)+Par(1)*Par(3)+Par(8)*Par(6))
             else
                Conc(jS,i)=cInit(Conc(jS,i),Par)
             end if
          end if
       End Do
    End Do

  End Subroutine MassInit


  !***********************************************************************

  !     Evaluate Liquid concentration from the total solute mass

  real(dp) function cInit(xMass,Par)
    Real(dp) :: Par(:),xMass
    Real(dp) :: x1,x2,XB1,XB2

    x1=1.e-3
    x2=1.e+3
    call ZBRAK1(X1,X2,XB1,XB2,xMass,Par)
    cInit=ZBRENT1(XB1,XB2,xMass,Par)

  end function cInit

  !***********************************************************************

  Real(dp) Function SolMass(Conc,xMass,Par)

    !     Calculate total solute mass for concentration Conc
    Real(dp) :: Par(:),xMass,Conc
    Real(dp) :: ro,frac,xKs,xNu,fExp,xKh,Theta,ThetaA,yMass
    ro=Par(1)
    frac=Par(2)
    xKs=Par(3)
    xNu=Par(4)
    fExp=Par(5)
    xKH=Par(6)
    Theta=Par(7)
    ThetaA=Par(8)

    yMass=Theta*Conc+ro*xKs*Conc**fExp/(1.+xNu*Conc**fExp)+&
         &      ThetaA*xKH*Conc
    SolMass=xMass-yMass

  End Function SolMass

  !***********************************************************************

  !     Bracketing of the root, Numerical recepies (345)

  Subroutine ZBRAK1(X1,X2,XB1,XB2,xMass,Par)
    Integer :: i,NB,NBB
    Real(dp) :: Par(:),X1,X2,XB1,XB2,xMass
    Real(dp) :: dlh,FC,FP,dx2

    If(debug) Print *,'zbrak1'
    NBB=1
    NB=1000
    dlh=(log10(X2)-log10(X1))/(NB-1)
    FP=SolMass(X1,xMass,Par)
    Do i=1,NB
       dx2=log10(X1)+(i)*dlh
       X2=10**dx2
       FC=SolMass(X2,xMass,Par)
       if(FC*FP.lt.0.) then
          NBB=NBB+1
          XB1=X1
          XB2=X2
          return
       end if
       FP=FC
       X1=X2
       if(NBB.eq.NB) return
    end do

  End Subroutine ZBRAK1

  !***********************************************************************

  !     Brent method of finding root that lies between x1 and x2, 
  !     Numerical recepies (354)

  Real(dp) Function ZBRENT1(X1,X2,xMass,Par)

    Integer, Parameter ::ITMAX=100
    Real(dp), Parameter :: EPS=3.E-8,TOL=1.e-6
    Real(dp) :: X1,X2,xMass,Par(:)
    Integer :: Iter
    Real(dp) :: A,B,C=0,D=0,E=0,FA,FB,FC,P,Q,R,S,Tol1,xm

    A=X1
    B=X2
    FA=SolMass(A,xMass,Par)
    FB=SolMass(B,xMass,Par)
    If(FB*FA.Gt.0.) Call WriteOutput('Root must be bracketed for ZBRENT1.')
    FC=FB
    Do ITER=1,ITMAX
       IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
       ENDIF
       IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
       ENDIF
       TOL1=2.*EPS*ABS(B)+0.5*TOL
       XM=.5*(C-B)
       IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT1=B
          RETURN
       ENDIF
       If(Abs(E).Ge.TOL1 .And. Abs(FA).Gt.Abs(FB)) Then
          S=FB/FA
          IF(A.EQ.C) THEN
             P=2.*XM*S
             Q=1.-S
          ELSE
             Q=FA/FC
             R=FB/FC
             P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
             Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.0*P .LT. MIN(3.0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
             E=D
             D=P/Q
          ELSE
             D=XM
             E=D
          ENDIF
       ELSE
          D=XM
          E=D
       ENDIF
       A=B
       FA=FB
       IF(ABS(D) .GT. TOL1) THEN
          B=B+D
       ELSE
          B=B+SIGN(TOL1,XM)
       ENDIF
       FB=SolMass(B,xMass,Par)
    End Do
    Call WriteOutput('ZBRENT1 exceeding maximum iterations.')
    ZBRENT1=B

  End Function ZBRENT1


  !     Subroutine calculating root solute uptake with and without compensation
  Subroutine SetSSnk(jS,x,Beta,Sink,SinkS,OmegaW,&
       cRootMax,lActRSU,OmegaS,SPot,rKM,cMin) !,ThN)
    Use Variables, Only: lNitrogen, lPhosphorus
    Use Geometry, Only: NumNP, NSnit, NitSink, NitSinkRoot, Pho,PhoSink,PhoSinkRoot
    Integer :: jS
    Real(dp) :: OmegaW,cRootMax,OmegaS,SPot,rKM,cMin
    Real(dp), Dimension(:) :: x,Beta,Sink,SinkS !,ThN
    Logical :: lActRSU,lLast
    Real(dp) :: AUptakeA,cc,Compen,dxM,Omega,Omega1,SAUptakeA,SAUptakeAN,SAUptakeP,&
         SPUptake
    Integer :: i,iStep,nStep
!    Real(dp) :: sinksum1,sinksum2
    If(debug) Print *,'setssink'
    !     Inputs:
    !     SPot      - potential root solute uptake
    !     OmegaS    - solute stress index
    !     rKM       - Michaelis-Menten constant
    !     lActRSU   - consider active root solute uptake
    !     cRootMax  - maximum concentration for the passive solute uptake

    !     From Water Flow
    !     Sink(i)   - Root water uptake
    !     OmegaW    - ratio of actual and potential transpiration

    !     SPUptake  - passive root solute uptake (step 1)
    !     SAUptakeP - potential active solute uptake (step 1)
    !     SAUptakeA - uncompensated actual active solute uptake (step 2)
    !     SAUptakeA - compensated actual active solute uptake (step 3)
    !     SinkS(i)  - local active solute uptake

    !     Initialization
    Compen=1.
    nStep=1
    if(lActRSU)                  nStep=2
    if(lActRSU.and.OmegaS.lt.1.) nStep=3
    !     step 1: Passive uptake
    !     step 2: Active uptake without compensation
    !     step 3: Active uptake with compensation
    lLast=.false.                 ! Active uptake only for the last solute
    if(lLast.and.jS.lt.NS) nStep=1
    Omega=0.
    Omega1=0.
    SPUptake=0.
    SinkS=0.
    AUptakeA=0.

    do iStep=1,nStep
       SAUptakeA=0.
       SAUptakeP=0.
       !sinksum1=0
       !sinksum2=0
       Do i=1,NumNP
          if(i.eq.NumNP) then
             dxM=(x(i)-x(i-1))/2.
          else if(i.eq.1) then
             dxM=(x(i)-x(i+1))/2.
          else
             dxM=(x(i+1)-x(i-1))/2.
          end if
          If( (lNitrogen .And. jS<=NSnit ) .Or. (lPhosphorus .And. jS==Pho)) Then
             If( jS<=NSnit ) Then
                ! NitSink from soil volume to liquid phase concentration
                SinkS(i)=-(NitSink(i,jS)+NitSinkRoot(i,jS)) ! /ThN(i)
                !sinksum1=sinksum1+Sink(i)
                !sinksum2=sinksum2+NitSinkRoot(i,jS)
             Else
                SinkS(i)=-(PhoSink(i)+PhoSinkRoot(i))
             End If
             SPUptake=SPUptake+SinkS(i)*dxM
             !  This is needed only for the last node, but that node may not have beta
             SAUptakeP=Max(SPot*OmegaW-SPUptake,0.)
          Elseif(Beta(i).Gt.0.0) Then
             cc=Max(Conc(jS,i)-cMin,0.)
             if(iStep.eq.1) then
                SinkS(i)=Sink(i)*Max(Min(Conc(jS,i),cRootMax),0.)
                SPUptake=SPUptake+SinkS(i)*dxM
                !             This is needed only for the last node, but that node may not have beta
                SAUptakeP=Max(SPot*OmegaW-SPUptake,0.)
             else if(iStep.eq.2) then
                AUptakeA=cc/(rKM+cc)*Beta(i)*SAUptakeP
                Omega=Omega+AUptakeA*dxM
                if(nStep.eq.2) SinkS(i)=SinkS(i)+AUptakeA
                !             This is needed only for the last node, but that node may not have beta
                SAUptakeA =Omega
                SAUptakeAN=Omega
                if(SAUptakeP.gt.0.) Omega1=Omega/SAUptakeP
             else if(iStep.eq.3) then
                !             This is needed only for the first node, but that node may not have beta
                if(Omega1.lt.OmegaS.and.Omega1.gt.0.) Compen=OmegaS
                if(Omega1.ge.OmegaS)                  Compen=Omega1
                if(Compen.gt.0.) AUptakeA=cc/(rKM+cc)*Beta(i)*SAUptakeP/Compen
                SinkS(i)=SinkS(i)+AUptakeA
                SAUptakeA=SAUptakeA+AUptakeA*dxM
             end if
          else
             SinkS(i)=0.
          End If
       End Do
       !Print *,'sinksum=',sinksum1,sinksum2
       !if(iStep.eq.nStep.and.jS.eq.NS) &
       !     &    write(78,'(3x,e14.7,1x,4e12.4)') t,SPUptake,SAUptakeP,SAUptakeA,SAUptakeAN ! the last is uncompensated
    End Do

  End Subroutine SetSSnk


  Subroutine SetChemBC(tnew,Prec,rSoil,dt,hNewTop)
    Use Variables, Only: Wlayer,KodTop,Transport,AtmBC
    Use timedata, Only: lMinStep
    Use Geometry, Only: NSfrom
    Implicit None
    Real(dp), Intent(in) :: Prec,rSoil,dt,hNewTop
    Integer, Intent(in) :: tnew
    Integer :: j
    Logical :: plus

    If(.Not.Transport) Return
    ! other solutes than nitrogen
    If(AtmBC) Then
       plus=.False.
       Do j=NSfrom+1,NS
          cTop(j)=0
          cBot(j)=0
          If(ATM_ind<=ATM_dim) Then
             If(tnew==ATM_Time(ATM_ind)) Then
                cTop(j)=ATM_Top(ATM_ind,NS-j+1)
                cBot(j)=ATM_Bot(ATM_ind,NS-j+1)
                plus=.True.
                Write(*,'(a,i8,a,1p,e13.5,a,e13.5)')&
                     'solute boundary, time=',tnew,' ctop=',cTop(j),' cbot=',cBot(j)
             Endif
          Endif
       End Do
    Endif
    If(plus) ATM_ind=ATM_ind+1

    If(WLayer.And.hNewTop.Gt.0.0) Then ! mass balance in the surface layer
       If((hNewTop+dt*(Prec-rSoil)).Gt.0.) &
            cTop=(hNewTop*cTopOld+dt*Prec*cTop)/(hNewTop+dt*(Prec-rSoil))
    Else
       If(Abs(KodTop).Eq.4.And.kTopCh.Le.0) Then
          If(Prec-rSoil.Gt.0.) Then
             cTop=cTop*Prec/(Prec-rSoil)
          Else If(rSoil.Gt.0.) Then
             cTop=0.
          End If
       End If
    endif

    Do j=1,NS
       If(.Not.lLinear(J).And.cTop(j).Gt.0.) lMinStep=.True.
       If(lActRSU.And.NS.Eq.1) Then
          cBot(j)=0.
          SPot=cBot(j)
       End If
    Enddo
  End Subroutine SetChemBC

End Module Solute
