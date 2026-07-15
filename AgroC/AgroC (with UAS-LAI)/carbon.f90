Module Carbon

  Use datatypes
  Implicit None

  Integer, Parameter :: MaxPool=5   ! number of CO2 pools
  ! DPM = decomposable plant material
  ! RPM = resistant plant material
  ! HUM = humus
  ! IOM = inert organic matter (inert = chem. not reacting)
  ! BIO = microbial biomass

  Integer :: NumTempInt
  Integer :: NumWCInt
  Integer :: kTopCO,kBotCO
  Integer :: iGasDiff
  Integer :: nPool,iredW,iredT,iredCO2,iReftemp,kProd
  Logical :: lStagn
  Real(dp) :: CO2Top,CO2Bot
  Real(dp) :: TempPar(5),WCPar(4)
  Real(dp) :: gamS0,gamR0,alf,PDDMax,xR,B1,B2, &
       cM1,cM2,HB1,HB2,fS6,Patm,gasconst,MolCO2,B1DPM, B1RPM, B1Bio, &
       B1Hum,co2rdpm,co2rrpm,co2rbio,co2rhum,co2pdepth,co2alf, &
       Reftemp,hCentury
  Real(dp), Allocatable :: eff_rate_man_from_nit(:), eff_rate_lit_from_nit(:) ! rates from nitrogen
  Real(dp) :: OCinpMan, OCinpLit ! for N and P organic fertilizer
  ! Arrays for nodes (dimension NumNP).
  Real(dp), Allocatable :: co2init(:,:)   ! dimension MaxPool x NumNP
  Real(dp), Allocatable :: co2len(:),PoolDPM(:),PoolRPM(:)
  Real(dp), Allocatable :: PoolBio(:),PoolHum(:),PoolIOM(:),PoolCO2(:)
  Real(dp), Allocatable :: diffDPM(:),diffRPM(:),diffBio(:),diffHum(:)
  Real(dp), Allocatable :: PoolDPMOld(:),PoolRPMOld(:),PoolBioOld(:),PoolHumOld(:)
  Real(dp), Allocatable :: gamSDPM(:),gamSRPM(:), gamSBIO(:),gamSHUM(:)
  Real(dp), Allocatable :: rnodeDPM(:),rnodeRPM(:)
  Real(dp), Allocatable :: CO2(:),CO2old(:),COIn(:),g0(:),g0h(:),g0r(:)
  Real(dp), Allocatable :: sfact(:,:)   ! dimension 3 x NumNP


  ! Arrays for materials (dimension NMat).
  Real(dp), Allocatable :: co2matlen(:)
  Real(dp), Allocatable :: co2parm(:,:)   ! dimension 3 x NMat

  Real(dp), Allocatable :: tempint(:), tempfact(:)   ! NumTempInt
  Real(dp), Allocatable :: wcint(:), wcfact(:)   ! NumWCInt


Contains


  Subroutine allocate_co2data()

    Use Geometry, Only: NumNP
    Use Material, Only: NMat
    Implicit None
    Allocate(co2init(MaxPool,NumNP))
    Allocate(co2len(NumNP))
    Allocate(PoolDPM(NumNP))
    Allocate(PoolRPM(NumNP))
    Allocate(PoolBio(NumNP))
    Allocate(PoolHum(NumNP))
    Allocate(diffDPM(NumNP))
    Allocate(diffRPM(NumNP))
    Allocate(diffBio(NumNP))
    Allocate(diffHum(NumNP))
    Allocate(PoolIOM(NumNP))
    Allocate(PoolCO2(NumNP))
    Allocate(PoolDPMOld(NumNP))
    Allocate(PoolRPMOld(NumNP))
    Allocate(PoolBioOld(NumNP))
    Allocate(PoolHumOld(NumNP))
    Allocate(gamSDPM(NumNP))
    Allocate(gamSRPM(NumNP))
    Allocate(gamSBio(NumNP))
    Allocate(gamSHum(NumNP))
    Allocate(CO2(NumNP))
    Allocate(CO2old(NumNP))
    Allocate(COIn(NumNP))
    Allocate(g0(NumNP))
    Allocate(g0h(NumNP))
    Allocate(g0r(NumNP))
    Allocate(sfact(3,NumNP))
    Allocate(eff_rate_man_from_nit(NumNP))
    Allocate(eff_rate_lit_from_nit(NumNP))
    Allocate(co2matlen(NMat))
    Allocate(co2parm(9,NMat))
    Allocate(rnodeDPM(NumNP),rnodeRPM(NumNP))
    
    ! initialize data
    g0=0.0
    PoolDPM=0
    PoolRPM=0
    PoolBio=0
    PoolHUM=0
    PoolIOM=0
    PoolCO2=0
    co2matlen=0
    rnodeDPM=0
    rnodeRPM=0

  End Subroutine allocate_co2data



  !     To assemble and solve the CO2 transport equation
  Subroutine Gas(dt,cvTop,cvBot,B,D,E,F,Retard,co2Sink)
    Use Geometry, Only: NumNP,MatNum,coord,vNew,vOld,ThNew,ThOld,vAir,Disp,Sink,TempO,TempN
    Use Material, Only: thS
    Implicit None
    Real(dp) :: dt
    Real(dp) :: cvTop,cvBot,B(:),D(:),E(:),F(:),Retard(:),co2Sink
    Integer :: i,j,M,Level,N
    Real(dp) :: x1,x2,BN,FN,DN,FE,D1,E1,F1,F2,DisA,DisW,gasdiff,dx,tha,th,Sum,v
    Real(dp) :: Temper,Henry
    Data d1,e1,f1,f2,fn,disa,temper,henry,x2/9*0.0/, M/0/

    ! To calculate the velocities of the gaseous phase
    N=NumNP
    Sum=0.0
    If(lStagn) Then
       vAir=0.0
    Else
       vAir(1)=0.0
       Do i=2,NumNP
          dx=coord(i)-coord(i-1)
          Sum=Sum+dx*(Sink(i)+Sink(i-1))/2.
          vAir(i)=-vNew(i)-Sum+vNew(1)
       End Do
    End If

    !     To construct the matrix equation
    Do Level=1,2
       Do i=1,NumNP
          M=MatNum(i)
          If(Level.Eq.1) Then
             th=ThOld(i)
             v=vOld(i)
             Temper=TempO(i)+273.15
             Henry=10.**(-13.417+2299.6/Temper+0.01422*Temper)/101.3*8.314*Temper
             vOld(i)=vOld(i)*Henry+vAir(i)
             Sink(i)=-Sink(i)*Henry
          Else
             th=ThNew(i)
             v=vNew(i)
             Temper=TempN(i)+273.15
             Henry=10.**(-13.417+2299.6/Temper+0.01422*Temper)/101.3*8.314*Temper
             vNew(i)=vNew(i)*Henry+vAir(i)
          End If
          tha=thS(M)-th+0.0001
          Retard(i)=tha+Henry*th
          DisW=co2parm(5,M)*(1.+0.026942*(Temper-293.15))
          DisA=co2parm(4,M)*(1.+0.007576*(Temper-293.15))
          If(iGasDiff.Eq.1) Then
             !     gaseous diffusion after Millington and Quirk 1961
             gasdiff=tha**3.333/thS(M)**2

          Else If(iGasDiff.Eq.2) Then
             !     gaseous diffusion after Moldrup et al 2000
             gasdiff=tha**2.5/thS(M)

          !start modification A.Klosterhalfen
          Else if(iGasDiff .eq. 3) then
             !     gaseous diffusion after Moldrup, Olesen et al 2000 (Par(26)=eps100, Par(27)=campb)
             gasdiff = (2*co2parm(7,M)**3+0.04*co2parm(7,M))* &
                  (tha/co2parm(7,M))**(2.+3./co2parm(8,M))
          !end modification

          Else if(iGasDiff .eq. 4) then
             !     double linear
             If(th.Lt.0.213) Then 
                gasdiff=tha*0.93461691-0.18377
             Else 
                gasdiff=tha*0.00911683+0.00845
             Endif

          !start modification A.Klosterhalfen
          Else if(iGasDiff .eq. 5) then
             !     gaseous diffusion after Kristensen et al. 2010 (Par(26)=fracture prosity (eps*), Par(27)=fracture tortuosity factor (H), Par(28)=matrix tortuosity factor (Xm))
             If(tha.Le.co2parm(7,M)) then
                gasdiff=co2parm(8,M)*tha
             Else
                gasdiff=(co2parm(8,M)*co2parm(7,M))+(tha-co2parm(7,M))**co2parm(9,M)*((tha-co2parm(7,M))/(thS(M)-co2parm(7,M)))
             End If
         !end modification

          Else
             !     gaseous diffusion after Millington and Quirk 1961
             gasdiff=tha**3.333/thS(M)**2
          Endif
          Disp(i)=DisA*gasdiff+ &
               (DisW*th**3.333/thS(M)**2+co2parm(6,M)*Abs(v))*Henry
       End Do

       !     Lower boundary condition
       dx=coord(2)-coord(1)
       If(Level.Eq.1) Then
          F1=CO2(1)* &
               (-(Disp(1)+Disp(2))/4./dx-(2.*vOld(1)+vOld(2))/12.+ &
               dx/12./dt*(3.*Retard(1)+Retard(2))+dx/24.*(3.*Sink(1)+Sink(2))) &
               +CO2(2)* &
               ((Disp(1)+Disp(2))/4./dx-(vOld(1)+2.*vOld(2))/12.+ &
               dx/12./dt*(Retard(1)+Retard(2))+dx/24.*(Sink(1)+Sink(2)))+ &
               dx/12.*(2.*g0(1)+g0(2))
       Else
          D1=(Disp(1)+Disp(2))/4./dx+(2.*vNew(1)+vNew(2))/12.+ &
               dx/12./dt*(3.*Retard(1)+Retard(2))-dx/24.*(3.*Sink(1)+Sink(2))
          E1=-(Disp(1)+Disp(2))/4./dx+(vNew(1)+2.*vNew(2))/12.+ &
               dx/12./dt*(Retard(1)+Retard(2))-dx/24.*(Sink(1)+Sink(2))
          F2=dx/12.*(2*g0(1)+g0(2))
          F1=F1+F2
       End If
       If(kBotCO.Eq.0) Then
          D(1)=-1.
          E(1)=1.
          F(1)=0.
       Else If(kBotCO.Eq.1) Then
          D(1)=1.
          E(1)=0.
          F(1)=CO2Bot
       Else
          If(Level.Eq.1) Then
             Temper=TempO(1)+273.15
             Henry=10.**(-13.417+2299.6/Temper+0.01422*Temper)/101.3* &
                  8.314*Temper
             F(1)=F1+0.5*vOld(1)*Henry*CO2(1)
             cvBot=0.5*vOld(1)*Henry*CO2(1)
          Else
             E(1)=E1
             Temper=TempN(1)+273.15
             Henry=10.**(-13.417+2299.6/Temper+0.01422*Temper)/101.3* &
                  8.314*Temper
             D(1)=D1-0.5*vNew(1)*Henry
             F(1)=F(1)+F2
          End If
       End If

       Do i=2,NumNP-1
          x1=coord(i)-coord(i-1)
          x2=coord(i+1)-coord(i)
          dx=0.5*(coord(i+1)-coord(i-1))
          If(Level.Eq.1) Then
             F(i)=CO2(i-1)* &
                  ((Disp(i)+Disp(i-1))/4./x1+ &
                  (vOld(i)+2.*vOld(i-1))/12.+ &
                  x1/12./dt*(Retard(i-1)+Retard(i))+ &
                  x1/24.*(Sink(i-1)+Sink(i)))+ &
                  CO2(i)* &
                  ((x1*(Retard(i-1)+3.*Retard(i))+ &
                  x2*(3.*Retard(i)+Retard(i+1)))/12./dt- &
                  (vOld(i+1)-vOld(i-1))/12.- &
                  (Disp(i+1)+Disp(i))/4./x2- &
                  (Disp(i)+Disp(i-1))/4./x1+ &
                  (x1*(Sink(i-1)+3.*Sink(i))+ &
                  x2*(3.*Sink(i)+Sink(i+1)))/24.)+ &
                  CO2(i+1)* &
                  ((Disp(i+1)+Disp(i))/4./x2- &
                  (2.*vOld(i+1)+vOld(i))/12.+ &
                  x2/12./dt*(Retard(i)+Retard(i+1))+ &
                  x2/24.*(Sink(i)+Sink(i+1)))+ &
                  x1/12.*(g0(i-1)+2.*g0(i))+x2/12.*(2.*g0(i)+g0(i+1))
          Else
             B(i)=-(Disp(i)+Disp(i-1))/4./x1- &
                  (vNew(i)+2.*vNew(i-1))/12.+ &
                  x1/12./dt*(Retard(i-1)+Retard(i))- &
                  x1/24.*(Sink(i-1)+Sink(i))
             D(i)=(x1*(Retard(i-1)+3*Retard(i))+ &
                  x2*(Retard(i+1)+3*Retard(i)))/12./dt+ &
                  (Disp(i)+Disp(i-1))/4./x1+ &
                  (Disp(i+1)+Disp(i))/4./x2+ &
                  (vNew(i+1)-vNew(i-1))/12.- &
                  (x1*(Sink(i-1)+3.*Sink(i))+ &
                  x2*(3.*Sink(i)+Sink(i+1)))/24.
             E(i)=-(Disp(i+1)+Disp(i))/4./x2+ &
                  (2.*vNew(i+1)+vNew(i))/12.+ &
                  x2/12./dt*(Retard(i)+Retard(i+1))- &
                  x2/24.*(Sink(i)+Sink(i+1))
             F(i)=F(i)+x1/12.*(g0(i-1)+2.*g0(i))+ &
                  x2/12.*(2*g0(i)+g0(i+1))
          End If
       End Do

       !     Upper boundary condition
       If(Level.Eq.1) Then
          FN=CO2(N-1)* &
               ((Disp(N-1)+Disp(N))/4./x2+(2.*vOld(N-1)+vOld(N))/12.+ &
               x2/dt/12.*(Retard(N-1)+Retard(N))+x2/24.*(Sink(N-1)+Sink(N)))+ &
               CO2(N)* &
               (-(Disp(N-1)+Disp(N))/4./x2+(vOld(N-1)+2.*vOld(N))/12.+ &
               x2/12./dt*(Retard(N-1)+3.*Retard(N))+ &
               x2/24.*(Sink(N-1)+3*Sink(N)))+ &
               x2/12.*(g0(N-1)+2.*g0(N))
       Else
          BN=-(Disp(N-1)+Disp(N))/4./x2-(2.*vNew(N-1)+vNew(N))/12.+ &
               x2/dt/12.*(Retard(N-1)+Retard(N))-x2/24.*(Sink(N-1)+Sink(N))
          DN=(Disp(N-1)+Disp(N))/4./x2- &
               (vNew(N-1)+2.*vNew(N))/12.+ &
               x2/12./dt*(Retard(N-1)+3.*Retard(N))- &
               x2/24.*(Sink(N-1)+3.*Sink(N))
          FE=x2/12.*(g0(N-1)+2.*g0(N))
          FN=FN+FE
          If(kTopCO.Lt.0) Then
             DisA=co2parm(4,M)*(1.+0.007576*(Temper-293.15))
             B(N)=BN
             D(N)=DN+0.5*DisA/CO2Top
             F(N)=FN-0.5*DisA/CO2Top*CO2(N)+DisA/CO2Top*0.00033
             cvTop=0.5*DisA/CO2Top*(CO2(N)-0.00033)
          Else
             B(N)=0.
             D(N)=1.
             F(N)=CO2Top
          End If
       End If
    End Do

    !     Solve matrix equation
    Do i=2,NumNP
       D(i)=D(i)-B(i)*E(i-1)/D(i-1)
       F(i)=F(i)-B(i)*F(i-1)/D(i-1)
    End Do
    CO2(N)=F(N)/D(N)
    Do i=2,NumNP
       j=N-i+1
       CO2(j)=(F(j)-E(j)*CO2(j+1))/D(j)
    End Do
    
    !modification N.Prolingheuer - separating R into Ra and Rh
    Do i=2,NumNP
       j=N-i+1
       !print *, 'CO2',j,(CO2(j)-CO2old(j))/CO2(j)*100.0
       CO2old(j)=CO2(j)
    End Do
    !end N.Prolingheuer
    
    !     Set up mass flux
    If(kTopCO.Eq.-1) Then
       cvTop=cvTop+0.5*DisA/CO2Top*(CO2(N)-0.00033)
    Else
       cvTop=FN-BN*CO2(N-1)-DN*CO2(N)
    End If
    If(kBotCO.Eq.-1) Then
       cvBot=cvBot+0.5*vNew(1)*Henry*CO2(1)
    Else
       cvBot=D1*CO2(1)+E1*CO2(2)-F1
    End If

    co2Sink=0.
    Do i=1,NumNP
       Temper=TempO(i)+273.15
       Henry=10.**(-13.417+2299.6/Temper+0.01422*Temper)/101.3* &
            8.314*Temper
       vOld(i)=(vOld(i)-vAir(i))/Henry
       If(i.Ne.NumNP) co2Sink=co2Sink-(Sink(i)*CO2(i)+Sink(i+1)*CO2(i+1))/2.* &
            (coord(i+1)-coord(i))
       Sink(i)=-Sink(i)/Henry
       Temper=TempN(i)+273.15
       Henry=10.**(-13.417+2299.6/Temper+0.01422*Temper)/101.3* &
            8.314*Temper
       vNew(i)=(vNew(i)-vAir(i))/Henry
    End Do
    Return
  End Subroutine Gas


  !***********************************************************************
  !     Production of CO2
  !***********************************************************************
  Subroutine Produc(t,dt,cumT,vProd,vProdh,vProdr, molProdr, molProdh, belowground_respiration, outp,output)

    Use Geometry
    Use Material
    Use Variables, Only: rRoot, fR6, root_respiration, heterotrophic_respiration,&
         aboveground_respiration, TER, NPP, NEE, PlantsExist, lNitrogen, LNitrogenReady, LNPools
    Use TimeData, Only: resetPlantTime, EpsTime
    Implicit None
    Real(dp) :: t,dt,cumT,vProd,vProdh,vProdr, molProdr, molProdh, belowground_respiration
    Logical :: outp, output(:)
    ! local variables
    Integer :: i,j,m
    Real(dp) :: expo,dpmnew,rpmnew,bionew,humnew,HB0,help,scal,Tref
    Real(dp) :: TSMDmax,TSMDacc,Volm,xbh,xco
    Real(dp) :: fS1,fS2,fS3,fS4,ax,dx
    Real(dp) :: fR2,fR3,fR4,fR5
    Real(dp) :: co2mol,co2molRPM,co2molDPM,co2molBIO,co2molHUM
    Real(dp) :: gamS,gamRPM,gamDPM,gamBIO,gamHUM
    Real(dp) :: diff
    Real(dp) :: fs3DPM,fs3RPM,fs3Bio,fs3Hum


    If(debug) Print *,'Produc'
    fs3hum=0.0
    fS3=-99999.9
    heterotrophic_respiration=null_dp
    !     Time reduction for plant production
    If(PDDMax.Gt.0.) Then
       fR5=CumT/PDDMax
       If(fR5.Gt.1.) fR5=1.
    Else
       fR5=1.
    End If

    If(outp .And. output(10)) &
         Write(79,*) 'simulation time:',t


    Do i=1,NumNP
       M=MatNum(i)

       !       Space and tension reduction for plant production
       fR2=0.
       If(rRoot.GT.0.) then 
          fR2=Sink(i)/rRoot
       End If

       !       Space reduction for soil production
       If(KProd.Eq.1) Then
          If(coord(i).Lt.coord(NumNP)-xR) Then
             fS1=0.
          Else
             fS1=1/xR
          End If
       Else If(KProd.Eq.2) Then
          If(coord(i).Lt.coord(NumNP)-xR) Then
             fS1=0.
          Else If(coord(i).Lt.coord(NumNP)-0.2*xR) Then
             fS1=2.08333/xR*(1-(coord(NumNP)-coord(i))/xR)
          Else
             fS1=1.6667/xR
          End If
       Else
          ax=-alf*(coord(NumNP)-coord(i))
          If(ax.Lt.-85.) Then
             fS1=0.
          Else
             fS1=Exp(ax)*alf
          End If
       End If

       !       Tension reduction for soil production
       If(KProd.Eq.3) Then
          If (iredW.Eq.0) Then
             !          SOILCO2 (Simunek & Suarez, 1993)
             HB0=-1.0/Par(5,m)
             If(HB1.Le.hNew(i).And.hNew(i).Le.HB0) Then
                fS2=(Log10(Abs(hNew(i)))-Log10(Abs(HB0)))/ &
                     (Log10(Abs(HB1))    -Log10(Abs(HB0)))
             Else If(HB2.Le.hNew(i).And.hNew(i).Lt.HB1) Then
                fS2=(Log10(Abs(hNew(i)))-Log10(Abs(HB2)))/ &
                     (Log10(Abs(HB1))    -Log10(Abs(HB2)))
             Else
                fS2=0.0
             End If
          Else If (iredW.Eq.1) Then
             !          SOILCO2 (original version)
             If(hNew(i).Le.HB2) Then
                fS2=0.
             Else If(hNew(i).Gt.HB1) Then
                fS2=1.
             Else
                fS2=(Log10(Abs(hNew(i)))-Log10(Abs(HB2)))/ &
                     (Log10(Abs(HB1))    -Log10(Abs(HB2)))
             End If
          Else If(iredW.Eq.2) Then
             !          DAISY
             If(hNew(i).Ge.-1.0) Then
                fS2=0.6
             Else If(hNew(i).Lt.-1.0.And.hNew(i).Ge.-10**(1.5)) Then
                fS2=0.6+(0.4*Log10(-1.0*hNew(i))/1.5)
             Else If(hNew(i).Lt.-10**(1.5).And.hNew(i).Ge.-10**(2.5)) &
                  Then
                fS2=1.0
             Else If(hNew(i).Lt.-10**(2.5).And.hNew(i).Ge.-10**(6.5)) &
                  Then
                fS2=1.625-(Log10(-1.0*hNew(i))/4.0)
             Else
                fS2=0.0
             End If
          Else If(iredW.Eq.8) Then
             !          new DAISY
             If(hNew(i).Ge.-10**(2.5)) &
                  Then
                fS2=1.0
             Else If(hNew(i).Lt.-10**(2.5).And.hNew(i).Ge.-10**(6.5)) &
                  Then
                fS2=1.625-(Log10(-1.0*hNew(i))/4.0)
             Else
                fS2=0.0
             End If
          Else If(iredW.Eq.3) Then
             !          CANDY
             If(ThNew(i)/ths(m).Le.0.5) Then
                fS2=4.0*ThNew(i)/ths(m)*(1-ThNew(i)/ths(m))
             Else
                fS2=1.0 
             End If
          Else If(iredW.Eq.4) Then
             !          PATCIS - mineral soil
             help=1.0-Exp(-22.6*ThNew(i)+0.11)
             If(help.Lt.0.0) Then
                fS2=0.0
             Else
                fS2=help
             End If
          Else If(iredW.Eq.5) Then
             !          CENTURY
             !fS2=1.0/(1.0+30.*exp(-8.5*(ThNew(i)-thr(m))/(ths(m)-thr(m))))
             fS2=1.0/(1.0+30.*Exp(-8.5 * ThNew(i)/FQ(hCentury,Par(:,m))))
          Else If(iredW.Eq.6) Then
             !          ROTHC
             TSMDmax=ths(m)-thr(m)
             TSMDacc=ths(m)-ThNew(i)
             If(TSMDacc.Lt.0.444*TSMDmax) Then
                fS2=1.0
             Else
                fS2=0.2+0.8*(TSMDmax-TSMDacc)/(0.556*TSMDmax)
             End If
          Else If (iredW.Eq.7) Then
             !          no reduction
             fS2=1.0
          Else If(iredW.Eq.10) Then
             !          reduction factors for moisture intervals from input
             Do j=1,NumWCInt-1
                If(ThNew(i).Le.wcint(j)) Exit
             Enddo
             fS2=wcfact(j)
          Else If(iredW.Eq.11) Then
             !          stepwise linear function
             If(ThNew(i).Lt.WCPar(1)) Then
                fs2=0.0
             Else If(ThNew(i).Lt.WCPar(2)) Then
                fs2=(ThNew(i)-WCPar(1))/(WCPar(2)-WCPar(1))
             Else If(ThNew(i).Lt.WCPar(3)) Then
                fs2=1.0
             Else If(ThNew(i).Lt.WCPar(4)) Then
                fs2=(ThNew(i)-WCPar(4))/(WCPar(3)-WCPar(4))
             Else
                fs2=0.0
             Endif
         
             !modification A.Klosterhalfen - May 2014
          Else If(iredW.Eq.12) Then
             !          after Bauer et al. 2012, Skopp et al. 1990
             fS2=(Exp(HB1*ThNew(i)+HB2*ThNew(i)**2))/(Exp((-HB1**2)/(4*HB2)))

             !          after M.Herbst - Jan 2015
          Else If(iredW.Eq.13) Then
             fS2=HB1+HB2*ThNew(i)**2
             If(fS2.Gt.1.0) Then
                fS2=1.0
             Else If(fS2.Lt.0.0) Then
                fS2=0.0
             EndIf
             !end modification
         
          Else
             Write(*,*) 'Ungueltiger Wert fuer iredW'
             Stop            
          End If
       Else
          !          SOILCO2 (original version)
          If(hNew(i).Le.HB2) Then
             fS2=0.
          Else If(hNew(i).Gt.HB1) Then
             fS2=1.
          Else
             fS2=(Log10(Abs(hNew(i)))-Log10(Abs(HB2)))/ &
                  (Log10(Abs(HB1))    -Log10(Abs(HB2)))
          End If
       End If

       !       Temperature reduction
       If(KProd.Eq.3) Then
          !          Reference temperature for scaling of reduction functions
          If(iReftemp.Eq.0) Then
             !             reference temperature from input
             Tref=Reftemp
          Else If(iReftemp.Eq.1) Then
             !             reference temperature of SOILCO2 
             Tref=20.0
          Else If(iReftemp.Eq.2) Then
             !             reference temperature of DAISY
             Tref=10.0
          Else If(iReftemp.Eq.3) Then
             !             reference temperature of CANDY
             Tref=35.0
          Else If(iReftemp.Eq.4) Then
             !             reference temperature of PATCIS
             Tref=10.0
          Else If(iReftemp.Eq.5) Then
             !             reference temperature of CENTURY
             Tref=Tan((1-0.56)/0.465)/0.097+15.7
          Else If(iReftemp.Eq.6) Then
             !             reference temperature of ROTHC
             Tref=106.0/Log(46.9)-18.3
          Else
             Write(*,*) 'Ungueltiger Wert fuer iReftemp'
             Stop
          End If
          If(iredT.Eq.1) Then
             !             SOILCO2 with variable reference temperature
             fS3=Exp(B1*(tempN(i)-Tref)/(tempN(i)+273.15)/(273.15+Tref)) 
          Else If(iredT.Eq.8) Then
             !             SOILCO2 with variable reference temperature,
             !             pool dependant activation energy
             fS3DPM=Exp(B1DPM*(tempN(i)-Tref) / &
                  (tempN(i)+273.15)/(273.15+Tref)) 
             fS3RPM=Exp(B1RPM*(tempN(i)-Tref) / &
                  (tempN(i)+273.15)/(273.15+Tref)) 
             fS3Bio=Exp(B1Bio*(tempN(i)-Tref) / &
                  (tempN(i)+273.15)/(273.15+Tref)) 
             fS3Hum=Exp(B1Hum*(tempN(i)-Tref) / &
                  (tempN(i)+273.15)/(273.15+Tref)) 
          Else If(iredT.Eq.9) Then
             !             SOILCO2 with variable reference temperature,
             !             temperature dependant activation energy
             If(tempN(i).Lt.10.0) Then
                fS3=Exp(B1DPM*(tempN(i)-Tref) / &
                     (tempN(i)+273.15)/(273.15+Tref)) 
             Else If(tempN(i).Lt.20.0) Then
                fS3=Exp(B1RPM*(tempN(i)-Tref) / &
                     (tempN(i)+273.15)/(273.15+Tref)) 
             Else If(tempN(i).Lt.30.0) Then
                fS3=Exp(B1Bio*(tempN(i)-Tref) / &
                     (tempN(i)+273.15)/(273.15+Tref)) 
             Else 
                fS3=Exp(B1Hum*(tempN(i)-Tref) / &
                     (tempN(i)+273.15)/(273.15+Tref)) 
             Endif
          Else If(iredT.Eq.2) Then
             !             DAISY with variable reference temperature
             If(0.0.Lt.Tref.And.Tref.Le.20.0) Then
                scal=1.0/(0.1*Tref)
             Else If(Tref.Gt.20.0) Then
                scal=1.0/(Exp(0.47-0.027*Tref+0.00193*Tref**2))
             Else
                scal=1.0
             End If
             If(tempN(i).Le.0.0) Then
                fS3=scal*0.0
             Else If(tempN(i).Gt.20.0) Then
                fS3=scal*(Exp(0.47-0.027*tempN(i)+0.00193*tempN(i)**2))
             Else
                fS3=scal*(0.1*tempN(i))
             End If
          Else If(iredT.Eq.3) Then
             !             CANDY with variable reference temperature
             If(tempN(i).Le.35.0) Then
                fS3=2.1**((tempN(i)-Tref)/10.0)  
             Else
                fS3=2.1**((35.0-Tref)/10.0)
             End If
          Else If(iredT.Eq.4) Then
             !          PATCIS with variable reference temperature
             If(tempN(i).Gt.20.0) Then
                fS3=Exp(78.2*1000/8.314*(tempN(i)-Tref)/ &
                     (tempN(i)+273.15)/(Tref+273.15))
             Else If(tempN(i).Lt.10) Then
                fS3=Exp(94.9*1000/8.314*(tempN(i)-Tref)/ &
                     (tempN(i)+273.15)/(Tref+273.15))
             Else
                fS3=Exp(79.3*1000.0/8.314*(tempN(i)-Tref)/ &
                     (tempN(i)+273.15)/(Tref+273.15))
             End If
          Else If(iredT.Eq.5) Then
             !          CENTURY with variable reference temperature
             scal=1.0/(0.56+0.465*Atan(0.097*(Tref-15.7)))
             fS3=scal*(0.56+0.465*Atan(0.097*(tempN(i)-15.7)))
          Else If(iredT.Eq.6) Then
             !          ROTHC  with variable reference temperature 
             !          (assumption: air temperature = soil temperature)
             scal=1.0/(47.9/(1+Exp(106.0/(Tref+18.3))))
             fS3=scal*(47.9/(1+Exp(106.0/(tempN(i)+18.3))))
          Else If(iredT.Eq.7) Then
             !          no reduction
             fS3=1.0
          Else If(iredT.Eq.10) Then
             !          reduction factors for temperature intervals from input
             Do j=1,NumTempInt-1
                If(tempN(i).Le.tempint(j)) Exit
             Enddo
             fS3=tempfact(j)
          Else If(iredT.Eq.11) Then
             !          Kirschbaum et al 2000
             fS3=Exp(TempPar(1)+TempPar(2)*tempN(i) &
                  +TempPar(3)*tempN(i)**2)
          Else If(iredT.Eq.12) Then
             !          Parton et al 1987
             If ((tempN(i).Le.0.0)) Then
                fS3=0.0
             Else 
                fS3=TempPar(1)+TempPar(2)*tempN(i)**TempPar(4) &
                     -(tempN(i)/TempPar(3))**TempPar(5)
             End If
             If (fS3.Lt.0.0) Then
                fS3=0.0
             End If
          Else If(iredT.Eq.13) Then
             !          O'Connel, 1990
             fS3=Exp(TempPar(1)+TempPar(2)*tempN(i) &
                  *(1-0.5*tempN(i)/TempPar(3)))
          Else
             Write(*,*) 'Ungueltiger Wert fuer iredT'
             Stop
          End If

       Else
          !          SOILCO2 referred to 20C (original expression)
          fS3=Exp(B1*(tempN(i)-20.0)/(tempN(i)+273.15)/293.15)
       End If

       fR3=Exp(B2*(tempN(i)-20.0)/(tempN(i)+273.15)/293.15)

       !       Carbon dioxide concentration reduction - microorganisms
       If(KProd.Eq.3) Then
          If(iredCO2.Eq.0) Then
             !       SOILCO2  (original version)
             If(CO2(i).Lt.0.21) Then
                fS4=(0.21-CO2(i))/(0.42-CO2(i)-cM1)
             Else
                fS4=0.
             End If
          Else If (iredCO2.Eq.1) Then
             !       SOILCO2; modified by M.Herbst 21/02/2006
             If(CO2(i).Lt.0.21) Then
                fS4=(0.21-CO2(i))/(0.42-CO2(i)-cM1)+1.0-0.21/(0.42-cM1)
             Else
                fS4=0.0
             End If
          Else If (iredCO2.Eq.3) Then
             !       no reduction
             fS4=1.0
          Else
             Write(*,*) 'Ungueltiger Wert fuer iredCO2'
             Stop
          End If
       Else  
          !       SOILCO2; modified by M.Herbst 21/02/2006
          !       original: fS4=(0.21-CO2(i))/(0.42-CO2(i)-cM1)
          If(CO2(i).Lt.0.21) Then
             fS4=(0.21-CO2(i))/(0.42-CO2(i)-cM1)+1.0-0.21/(0.42-cM1)
          Else
             fS4=0.
          End If
       End If

       !       Carbon dioxide concentration reduction - plant roots
       If(CO2(i).Lt.0.21) Then
          fR4=(0.21-CO2(i))/(0.42-CO2(i)-cM2)
       Else
          fR4=0.
       End If

       !          ideal gas: molvolume=RT/P
       Volm=gasconst*(273.15+tempN(i))/Patm   ! cm^3 /mol
       If(KProd.Eq.3) Then
          !          calculate decomposition of the pools
          If(iredT.Eq.8) Then 
             expo=-fS2*fS4*dt
             dpmnew=PoolDPM(i)*Exp(fs3DPM*expo*co2rdpm)
             rpmnew=PoolRPM(i)*Exp(fs3RPM*expo*co2rrpm)
             bionew=PoolBIO(i)*Exp(fs3Bio*expo*co2rbio)
             humnew=PoolHUM(i)*Exp(fs3Hum*expo*co2rhum)
          Else
             expo=-fS2*fS3*fS4*dt
             If(lNitrogen .And. LNitrogenReady .And. LNPools) Then
                dpmnew=PoolDPM(i)*Exp(-eff_rate_lit_from_nit(i)*dt)
                rpmnew=PoolRPM(i)*Exp(-eff_rate_man_from_nit(i)*dt)
             Else
                dpmnew=PoolDPM(i)*Exp(expo*co2rdpm)
                rpmnew=PoolRPM(i)*Exp(expo*co2rrpm)
             Endif
             bionew=PoolBIO(i)*Exp(expo*co2rbio)
             humnew=PoolHUM(i)*Exp(expo*co2rhum)
          Endif
          diff=PoolDPM(i)-dpmnew + PoolRPM(i)-rpmnew +  &
               PoolBIO(i)-bionew + PoolHUM(i)-humnew
          diffDPM(i)=PoolDPM(i)-dpmnew
          diffRPM(i)=PoolRPM(i)-rpmnew
          diffBIO(i)=PoolBIO(i)-bionew
          diffHUM(i)=PoolHUM(i)-humnew
          !write(*,*) PoolHUM(50)
          xbh=co2parm(2,m)
          xco=co2parm(3,m)
          !          update pools
          !rnodexu(i) = rnodexu(i)
          rnodeDPM(i)=(rnodexu(i) + (rnodedeadw(i)+rNodeResC(i))*0.59)*dt
          rnodeRPM(i)=              (rnodedeadw(i)+rNodeResC(i))*0.41*dt
          PoolDPM(i)=dpmnew + rnodeDPM(i)
          PoolRPM(i)=rpmnew + rnodeRPM(i)
          PoolBIO(i)=bionew
          PoolHUM(i)=humnew
          PoolCO2(i)=PoolCO2(i)+xco*diff
          !          0.044 = molar mass co2, 0.012 = molar mass c
          gamS=xco*diff/dt*0.0440098/0.0120107   ! kg CO2 /L^3 /T
          gamDPM=xco*diffDPM(i)/dt*0.0440098/0.0120107
          gamRPM=xco*diffRPM(i)/dt*0.0440098/0.0120107
          gamBIO=xco*diffBIO(i)/dt*0.0440098/0.0120107
          gamHUM=xco*diffHUM(i)/dt*0.0440098/0.0120107
          !write(*,*) gasconst/Patm/MolCO2
          !          mol=mass/molar mass C02
          co2mol=gamS/MolCO2   ! mol CO2 /L^3 /T
          rnodeh(i) = co2mol

          co2molDPM=gamDPM/MolCO2   ! mol CO2 /L^3 /T
          co2molRPM=gamRPM/MolCO2
          co2molBIO=gamBIO/MolCO2
          co2molHUM=gamHUM/MolCO2
          !          volume=molvolume*mol
          gamS=Volm*co2mol   ! cm^3 /L^3 /T
          !write(*,*) gamS, co2mol
          gamSDPM(i)=Volm*co2molDPM
          gamSRPM(i)=Volm*co2molRPM
          gamSBIO(i)=Volm*co2molBIO
          gamSHUM(i)=Volm*co2molHUM
          g0h(i)=gamS
          heterotrophic_respiration=heterotrophic_respiration+gamS
          If(PlantsExist) Then
             !g0r(i)=root_respiration*Volm*fR2
             g0r(i)=rnodert(i)*Volm   ! cm^3 /L^3 /day
          else
             g0r(i)=gamR0*fR2*fR3*fR4*fR5*fR6
          Endif
          g0(i)=g0r(i)+g0h(i)
          If(outp .And. output(10)) Then
             If(iredT.Eq.8) Then 
                Write(79,'(I6,6F10.6)') &
                     i,fS2,fS3DPM,fs3RPM,fs3Bio,fs3Bio,fS4
             Else
                Write(79,'(I6,3F10.6)') i,fS2,fS3,fS4
             Endif
          Endif
          sfact(1,i)=fs2
          sfact(2,i)=fs3
          sfact(3,i)=fs4
       Else
          g0h(i)=gamS0*fS1*fS2*fS3*fS4*fS6
          If(PlantsExist) Then
             g0r(i)=root_respiration*Volm*fR2
          else
             g0r(i)=gamR0*fR2*fR3*fR4*fR5*fR6
          Endif
          g0(i)=g0r(i)+g0h(i)
       End If
    End Do
    LNitrogenReady=.True.
    If(t > resetPlantTime+1d0-EpsTime) Then
       rNodeResC = 0
       rNodeResN = 0
       rNodeResP = 0
    Endif
    !modification N.Prolingheuer
    !GPP,TER,NEE,NPP
    !end N.Prolingheuer

    !modification N.Prolingheuer -  CO2 production roots (vProdr)
    vProdr=0.
    vProdh=0.
    vProd=0.
    molProdr=0.
    molProdh=0.
    Do i=2,NumNP
       dx=coord(i)-coord(i-1)
       vProdr=vProdr+dx*(g0r(i)+g0r(i-1))/2.
       molProdr = molProdr + dx*(rnodert(i)+rnodert(i-1))/2.
       vProdh=vProdh+dx*(g0h(i)+g0h(i-1))/2.
       molProdh = molProdh + dx*(rnodeh(i)+rnodeh(i-1))/2.
       vProd=vProd+dx*(g0(i)+g0(i-1))/2.
    End Do
    TER=aboveground_respiration+molProdr+molProdh
    belowground_respiration = molProdr+molProdh
    NEE=NPP+molProdh

    outp=.False.

    !print *, (g0r(NumNP-1)+2.*g0r(NumNP))/(g0(NumNP-1)+2.*g0(NumNP)), vProdr, vProd
    !end N.Prolingheuer
  End Subroutine Produc


  !***********************************************************************
  Subroutine OCinp(plantinp)
    Use Geometry, Only: NumNP,MatNum, OCinpLit, OCinpMan
    Implicit None
    Real(dp) :: plantinp
    Integer :: i
    Real(dp) :: diff,xbh,dpmnew,rpmnew

    If(KProd.Ne.3) Return
    Do i=1,NumNP
       xbh=co2parm(2,MatNum(i))
       diff=diffDPM(i) + diffRPM(i) +  diffBIO(i) + diffHUM(i)
       PoolBIO(i)=PoolBIO(i)+xbh*0.46*diff
       PoolHUM(i)=PoolHUM(i)+xbh*0.54*diff
    End Do
    !     distribute new plant material on nodes
    OCinpLit=0
    OCinpMan=0
    If(plantinp.Gt.0.0) Then 
       plantinp=plantinp/co2pdepth
       dpmnew=plantinp*0.59
       rpmnew=plantinp*0.41
       Do i=nPool,NumNP
          PoolDPM(i)=PoolDPM(i)+dpmnew
          PoolRPM(i)=PoolRPM(i)+rpmnew
          ! calc N+P input
          OCinpLit(i)=dpmnew
          OCinpMan(i)=rpmnew
       End Do
    End If
    Do i=1,NumNP
       PoolDPMOld(i)=PoolDPM(i)
       PoolRPMOld(i)=PoolRPM(i)
       PoolBIOOld(i)=PoolBIO(i)
       PoolHUMOld(i)=PoolHUM(i)
    End Do
  End Subroutine OCinp

  !***********************************************************************
  Subroutine OCinit()
    Use Geometry, Only: NumNP,MatNum,coord
    Use Material, Only: NMat
    Implicit None
    Integer :: i,m
    Real(dp) :: depth,exp_alt,exp_neu,w,xcbh

    If(KProd.Ne.3) Return

    !     determine nodes for distributing organic carbon
    nPool=NumNP
    depth=coord(NumNP)-0.99999*co2pdepth
    Do i=NumNP-1,1,-1
       If(coord(i).Lt.depth) Exit
       nPool=nPool-1
    End Do
    !     init co2len
    co2len(1)=0.5*(coord(2)-coord(1))
    co2len(NumNP)=0.5*(coord(NumNP)-coord(NumNP-1))
    Do i=2,NumNP-1
       co2len(i)=0.5*( coord(i+1)-coord(i-1) )
       m=MatNum(i)
       co2matlen(m)=co2matlen(m)+co2len(i)
    End Do
    If(co2alf.Le.null_dp) Then
       ! equally distribute init pool material on the nodes
       Do i=1,NumNP
          If(equal(co2alf,-1.0_dp)) Then
             m=i
          Else 
             m=MatNum(i)
          Endif
          PoolDPM(i)=co2init(1,m)
          PoolRPM(i)=co2init(2,m)
          PoolBio(i)=co2init(3,m)
          PoolHum(i)=co2init(4,m)
          PoolIOM(i)=co2init(5,m)
          PoolDPMOld(i)=PoolDPM(i)
          PoolRPMOld(i)=PoolRPM(i)
          PoolBioOld(i)=PoolBio(i)
          PoolHumOld(i)=PoolHum(i)
       End Do
    Else
       ! exponentially distribute init pool material on the nodes
       exp_alt=1.0
       exp_neu=1.0
       Do i=NumNP,2,-1
          exp_neu=Exp((coord(i-1)-coord(NumNP))*co2alf)
          w=exp_alt-exp_neu
          !write(*,*) w,coord(NumNP)-coord(i-1)
          exp_alt=exp_neu
          PoolDPM(i)=co2init(1,1)*w
          PoolRPM(i)=co2init(2,1)*w
          PoolBio(i)=co2init(3,1)*w
          PoolHum(i)=co2init(4,1)*w
          PoolIOM(i)=co2init(5,1)*w
          PoolDPMOld(i)=PoolDPM(i)
          PoolRPMOld(i)=PoolRPM(i)
          PoolBioOld(i)=PoolBio(i)
          PoolHumOld(i)=PoolHum(i)
       End Do
       w=exp_neu
       !write(*,*) w,coord(NumNP),coord(1)
       PoolDPM(1)=co2init(1,1)*w
       PoolRPM(1)=co2init(2,1)*w
       PoolBio(1)=co2init(3,1)*w
       PoolHum(1)=co2init(4,1)*w
       PoolIOM(1)=co2init(5,1)*w
       PoolDPMOld(1)=PoolDPM(1)
       PoolRPMOld(1)=PoolRPM(1)
       PoolBioOld(1)=PoolBio(1)
       PoolHumOld(1)=PoolHum(1)
    Endif
    !     get partition factor for co2 - (hum+bio) depending from clay
    Do i=1,NMat
       xcbh=1.67*(1.85+1.6*Exp(-0.0786*co2parm(1,i)*100.0))
       co2parm(3,i)=xcbh/(xcbh+1.0)
       co2parm(2,i)=1.0-co2parm(3,i) ! Fe
    End Do
    !     correct co2pdepth
    co2pdepth=0.0
    Do i=nPool,NumNP
       co2pdepth=co2pdepth+co2len(i)
    End Do
  End Subroutine OCinit

End Module Carbon
