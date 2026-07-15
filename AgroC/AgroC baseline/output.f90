Module Output

Contains

  subroutine TLInf(kk,Con1,Con2,ConN,ConB,x1,x2,xN,xB,CosAlf,t,dt, &
       Iter,TLevel,ShortF,TPrint,rTop,rRoot,vTop,vRoot,hNewN, &
       hNewB,hRoot,hNew1,hNew2,CumQ,ItCum,KodTop,KodBot, &
       ConvgF,lCO2,cRoot,CO2Top,cBot,cvTop,cvBot,vProdm, &
       vProd,CumT,tTop,co2Sink,ThNew,ThOld,SinkTop,T1,T2, &
       T3,wCumT,wCumA,COCumT,Yield,lRoot,output,COCumA, molProdr, molProdh, belowground_respiration)

    Use datatypes
    Use Geometry, Only: NumNP, coord, hNew
    Use TimeData, Only: iday_to_date, simtime_to_iday, dtMin
    Use variables, Only: GPP, NPP, TER, NEE, aboveground_respiration, co2_fluxes, respiration, &
         maint_growth, fluorescence_755nm,PlantsExist, Transport
    Use Solute, Only: NS,Conc,SorbOut,solute_name, cGWL, CumCh, cCumT, cCumA, chem_vTop, chem_vBot,&
         cvCh0, cvCh1, cvChR, cvChIm
    Use Plants, Only: SucrosPlant, PlantGeometry
    Implicit None 
    Integer :: i,j,js
    ! arguments
    Integer :: kk,Iter,TLevel,ItCum,KodTop,KodBot
    logical :: ShortF,ConvgF,lCO2,lRoot,output(:)
    Real(dp) :: Con1,Con2,ConN,ConB,x1,x2,xN,xB,CosAlf,t,dt
    Real(dp) :: TPrint,rTop,rRoot,vTop,vRoot,hNewN
    Real(dp) :: hNewB,hRoot,hNew1,hNew2,CumQ(:)
    Real(dp) :: cRoot,CO2Top,cBot,cvTop,cvBot,vProdm, molProdr, molProdh, belowground_respiration
    Real(dp) :: vProd,CumT,tTop,co2Sink,ThNew,ThOld,SinkTop,T1,T2
    Real(dp) :: T3,wCumT,wCumA,COCumT,Yield,COCumA
    ! local variables
    Real(dp) :: delta,vBot,dGWL,dx
    Integer :: yy,mm,dd

    If(debug) Print *,'TLInf'
    vTop=-(ConN+ConB)/2*((hNewN-hNewB)/(xN-xB)+CosAlf)- &
         (ThNew-ThOld)*(xN-xB)/2./dt-SinkTop*(xN-xB)/2.
    vBot=-(Con1+Con2)/2*((hNew2-hNew1)/(x2-x1)+CosAlf)
    CumQ(1)=CumQ(1)+rTop *dt
    CumQ(2)=CumQ(2)+rRoot*dt
    CumQ(3)=CumQ(3)+vTop *dt
    CumQ(4)=CumQ(4)+vRoot*dt
    CumQ(5)=CumQ(5)+vBot *dt
    wCumA=wCumA+(abs(vTop)+abs(vBot)+abs(vRoot))*dt
    wCumT=wCumT+(vBot-vTop-vRoot)*dt
    if(CumQ(2).gt.0.) Yield=CumQ(4)/CumQ(2)*100
    if(lCO2) then
       CumQ(6)=CumQ(6)+cvTop*dt
       CumQ(7)=CumQ(7)+cvBot*dt
       CumQ(8)=CumQ(8)+vProd*dt
       CumQ(9)=CumQ(9)+co2Sink*dt
       COCumA=COCumA+(abs(cvBot)+abs(cvTop)+abs(vProd)+abs(co2Sink))*dt
       COCumT=COCumT+(cvBot-cvTop+vProd-co2Sink)*dt
    end if
    delta=0.
    if(rRoot.gt.0.) delta=vRoot/rRoot
    If(lRoot)&
         CumT=CumT+dt*delta*(Max(tTop-T1,0.0d0)-Max(tTop-T2,0.0d0)-Max(tTop-T3,0.0d0))
    If(output(14)) Then
       if(float((TLevel+19)/20).eq.(TLevel+19)/20.) write(*,110) &
            'Time','Iter','ItCum','vTop','SvTop','SvRoot','SvBot','hTop','hRoot','hBot'
       write(*,120) t,Iter,ItCum,vTop,CumQ(3),CumQ(4),CumQ(5),hNewN,hRoot,hNew1
    Endif
    
    If(Transport) Then
       Do jS=1,NS
          CumCh(1,jS)=CumCh(1,jS)-chem_vTop(jS)*dt
          CumCh(2,jS)=CumCh(2,jS)+chem_vBot(jS)*dt
          CumCh(3,jS)=CumCh(3,jS)+cvCh0(jS)*dt
          CumCh(4,jS)=CumCh(4,jS)+cvCh1(jS)*dt
          CumCh(5,jS)=CumCh(5,jS)+cvChR(jS)*dt
          CumCh(6,jS)=CumCh(6,jS)+cvChIm(jS)*dt
          cCumT(jS)=cCumT(jS)+(chem_vTop(jS)-chem_vBot(jS)-cvCh0(jS)-cvCh1(jS)+cvChR(jS))*dt
          cCumA(jS)=cCumA(jS)+(Abs(chem_vBot(jS))+Abs(chem_vTop(jS))+&
               Abs(cvCh0(jS))+Abs(cvCh1(jS))+Abs(cvChR(jS)))*dt
          !         Average GWL concentration
          cGWL(jS)=0.
          dGWL=0
          Do i=1,NumNP-1
             j=i+1
             dx=coord(j)-coord(i)
             if(hNew(j)>0.0) then
                cGWL(jS)=cGWL(jS)+(Conc(jS,i)+Conc(jS,j))/2.*dx
                dGWL=dGWL+dx          
             Else
                Exit
             end if
          enddo
          if(dGWL.gt.0.) cGWL(jS)=cGWL(jS)/dGWL
       Enddo
    Endif

110 Format(A8,A5,A8,7A10)
120 Format(f10.3,i3,i8,1P,7e10.2)
    if(TLevel.eq.1) then
       If(output(3)) Then
          Write(71,210) &
            'Time','t-level','rTop','rRoot','vTop','vRoot','vBot', &
            'sum(vTop)','sum(vBot)','hTop','hRoot','hBot', &
            (Trim(solute_name(i))//'_w(bot)',Trim(solute_name(i))//'_w(top)',&
             Trim(solute_name(i))//'_s(bot)',Trim(solute_name(i))//'_s(top)',i=1,NS)
          Write(71,210) &
            '[T]','[-]', '[L/T]','[L/T]','[L/T]','[L/T]', &
            '[L/T]','[L]','[L]','[L]','[L]','[L]',&
            ('[M/L**3]','[M/L**3]','[M/(M soil)]','[M/(M soil)]',i=1,NS)
210       Format(a8,a10,200A15)
       Endif
       if(output(2)) write(70,140) &
            'TLevel','Time','dt','Iter','ItCum','KodT','KodB','Convergency'
       if(output(5) .and. lCO2 .and. kk.eq.1) write(73,150) &
            'Time','cvTop','cvBot','sum(cvTop)','sum(cvBot)', &
            'cTop','cRoot','cBot','vProd','vProd', &
            'sum(vProd)','sum(Sink)','TLevel', &
            '[T]','[L3/L2/T]','[L3/L2/T]','[L3/L2]','[L3/L2]', &
            '[L3/L3]','[L3/L3]','[L3/L3]','[L3/L2/T]','[L3/L2/T]', &
            '[L3/L2]','[L3/L2]'

       !modification D.Farber - output
       If(co2_fluxes) Then
          open(83,file='co2_fluxes.out',Status='REPLACE')
          Write(83,200) 'year month','day Time','rTop', 'sum(rTop)',&
               'Rh','Rr', 'RaboveG', 'RbelowG', 'GPP', 'NPP', 'NEE', 'TER', 'F_755nm'
!          write(83,201) '# [T]', '[L]', '[L]', '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]', &
!               '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]'
       End If
       If(PlantsExist) Then
          If(respiration) Then
             Open(84, file='respiration.out', Status='REPLACE')
             Write(84, fmt='(A10, A9,4A17)') 'year month','day Time', 'rlv', 'rst', 'rso', 'rrt'
!             Write(84, fmt='(A19,4A17)') '# [T]', '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]','[molCO2/L2/T]'
          End If
          If(maint_growth) Then
             Open(85, file='maint_growth.out', Status='REPLACE')
             Write(85, fmt='(A10, A9,10A17)') 'year month','day Time', 'rgrowthlv', 'rmaintlv', 'rgrowthst', &
                  'rmaintst', 'rgrowthso', 'rmaintso', 'rgrowthcrn', 'rmaintcrn', 'rgrowthrt', 'rmaintrt'
!             Write(85, fmt='(A19,10A17)') '# [T]', '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]','[molCO2/L2/T]', &
!                  '[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]','[molCO2/L2/T]', '[molCO2/L2/T]', '[molCO2/L2/T]'
          End If
          !end D.Farber
       Endif
    end if
    if(.not.ShortF.or.abs(TPrint-t).lt. 0.5*dtMin) then   !changed by D.Farber
       If(output(3)) Then
          Write(71,160) Nint(t),TLevel,rTop,rRoot,vTop,vRoot,vBot,CumQ(3),CumQ(5),hNewN,hRoot,hNew1,&
               (Conc(i,1),Conc(i,NumNP),SorbOut(i,1),SorbOut(i,NumNP),i=1,NS)
       Endif
       if(output(2)) write(70,170) TLevel,t,dt,Iter,ItCum,KodTop,KodBot,ConvgF
       if(lCO2.and.kk.eq.1 .and. output(5)) &
            write(73,180) t,cvTop,cvBot,CumQ(6),CumQ(7),CO2Top,cRoot,cBot, &
            vProdm,vProd,CumQ(8),CumQ(9),TLevel
    end if
    !modification N.Prolingheuer - output
    if(abs(t-nint(t)) .le. 0.5*dtMin) then   !changed by D.Farber
       Call iday_to_date(simtime_to_iday(t),yy,mm,dd)
       if(co2_fluxes) then
          write(83,190) yy,mm,dd, t, rTop,CumQ(1),molProdh, molProdr, aboveground_respiration, &
               belowground_respiration, GPP, NPP, NEE, TER, fluorescence_755nm
       end if
       If(PlantsExist) Then
          i=Max(1, PlantGeometry%index)
          if(respiration) then
             Write(84, fmt='(I4, 2I3, f9.2,1p,4e17.3)') yy,mm,dd, t, &
                  SucrosPlant(i)%rlv, SucrosPlant(i)%rst, SucrosPlant(i)%rso, SucrosPlant(i)%rrt
          end if
          if(maint_growth) then
             Write(85,fmt='(I4, 2I3, f9.2, 1p, 10e17.3)') yy,mm,dd, t, &
                  SucrosPlant(i)%rgrowthlv, SucrosPlant(i)%rlv-SucrosPlant(i)%rgrowthlv, &
                  SucrosPlant(i)%rgrowthst, SucrosPlant(i)%rst-SucrosPlant(i)%rgrowthst, &
                  SucrosPlant(i)%rgrowthso, SucrosPlant(i)%rso-SucrosPlant(i)%rgrowthso, &
                  SucrosPlant(i)%rgrowthcrn, SucrosPlant(i)%rcrn-SucrosPlant(i)%rgrowthcrn, &
                  SucrosPlant(i)%rgrowthrt, SucrosPlant(i)%rrt-SucrosPlant(i)%rgrowthrt
          end if
       end if
    end if
    !end N.Prolingheuer
140 Format(//A5,2A12,A5,3A6,A6)
170 format(i5,2e12.3,i5,3i6,L6)
160 Format(i8,i10,1P,200E15.7)
150 Format(/A12, 11A11,  A8) 
180 format(f12.5,11e11.3,i8)
!modification N.Prolingheuer - output
200 format(A10, A9,11A17)
201 format(A19,10A17)
190 format(I4, 2I3, f9.2,1p,15e17.3)
!end N.Prolingheuer
    return
  end subroutine TLInf

  !***********************************************************************

  subroutine ALInf(t,CumQ,hNewN,hRoot,hNew1,ALevel,cvTop,cvBot,CO2Top, &
       cRoot,cBot,vProdm,vProd,lCO2,output)

    Use datatypes
    Use Solute, Only : NS, solute_name, CumCh
    Use Variables, Only: Transport
    Implicit None 
    ! arguments
    integer ALevel
    logical lCO2, output(*)
    Real(dp) :: t,CumQ(20),hNewN,hRoot,hNew1,cvTop,cvBot,CO2Top
    Real(dp) :: cRoot,cBot,vProdm,vProd
    ! local variables
    Real(dp), Allocatable, Save :: vtopold(:),vbotold(:)
    Integer :: i
    Logical, Save :: first_call=.True.

    If(first_call) Then
       Allocate(vtopold(NS),vbotold(NS))
       vtopold=0
       vbotold=0
       first_call=.False.
    Endif
    If(ALevel.Eq.1) Then
       If(output(4)) Then
          Write(72,*)
          Write(72,110) 'Time','sum(rTop)','sum(rRoot)','sum(vTop)','sum(vRoot)', &
               'sum(vBot)','hTop','hRoot','hBot',&
               ('v'//Trim(solute_name(i))//'_top','v'//Trim(solute_name(i))//'_bot',i=1,NS)
          Write(72,110)  '[T]','[L]','[L]','[L]','[L]','[L]','[L]','[L]','[L]', &
               ('[M/L2/T]',i=1,2*NS)
       Endif
       if(output(5)) write(73,130) &
            'Time','cvTop','cvBot','sum(cvTop)','sum(cvBot)', &
            'cTop','cRoot','cBot','vProdm','vProd', &
            'sum(vProd)','sum(Sink)','ALevel', &
            '[T]','[L3/L2/T]','[L3/L2/T]','[L3/L2]','[L3/L2]', &
            '[L3/L3]','[L3/L3]','[L3/L3]','[L3/L2/T]','[L3/L2/T]', &
            '[L3/L2]','[L3/L2]'
    end if
    If(output(4)) Then
       Write(72,120) Nint(t),(CumQ(i),i=1,5),hNewN,hRoot,hNew1,&
            (CumCh(1,i)-vtopold(i),CumCh(2,i)-vbotold(i),i=1,NS)
       If(Transport) Then
          vtopold=CumCh(1,:)
          vbotold=CumCh(2,:)
       Endif
    Endif
    if(lCO2 .and. output(5))  &
         write(73,140) t,cvTop,cvBot,CumQ(6),CumQ(7),CO2Top,cRoot, &
         cBot,vProdm,vProd,CumQ(8),CumQ(9),ALevel

110 Format(A10,99A14)
120 Format(i10,1p,99e14.6)
130 Format(A12,11A14,A8)
140 format(f12.5,11e14.6,i8)
  end subroutine ALInf

  !***********************************************************************

  subroutine SubReg(ThN,ThO,t,dt,PLevel,wCumA,wCumT,wVolI,COCumT,COVolI, &
       Yield, output,COCumA)

    Use datatypes
    Use Geometry
    Use Carbon
    Use Variables
    Use Material, Only: Par,FQ,thS,tempParam
    Use Solute, Only: NS,ConSub,cMean,cTot,ConVol,cVolI,ChPar,lDualNEq,Conc,lLinear,iDualPor,lBact,&
         TDep,lMobIm,lEquil,ConVolIm,ConVolIm2,cTotIm,cCumT,cCumA,solute_name
    Implicit None 
    ! arguments
    logical :: output(:)
    integer :: PLevel
    Real(dp) ::ThN(:), ThO(:)
    Real(dp) ::t,dt,wCumA,wCumT,wVolI,COCumT,COVolI,Yield,COCumA
    ! local variables
    Integer :: N,i,j,Mi,Mj,L,jS
    Real(dp) :: thS1,thS2,Con1,Con2,ConN,ConM,aTot,hTotal,Vol,Change,DeltW
    Real(dp) :: dx,COBalR,COBalT,CONewi,COTot,COVol,DeltCO,v1,vN,vOldi,vNewi
    Real(dp) :: wBalR,wBalT,ww,Henry,deltC,cNewi,C1,C2,cE,f_em,R,Tr,TT,cEl
    Real(dp) :: TTi,xKsi,xNui,fExpi,Henryi
    Real(dp) :: TTj,xKsj,xNuj,fExpj,Henryj
    Real(dp) :: ThWi,ThWj,ThImobi,ThImobj,ThGi,ThGj
    ! temperature
    Real(dp) :: TTot,TVol,TE,TNewE
    Real(dp) :: cBalR,cBalT,cc

    If(debug) Print *,'SubReg'
    Tr=293.15
    R=8.314
    Con1=Con(1)
    Con2=Con(2)
    ConN=Con(NumNP)
    ConM=Con(NumNP-1)
    N=NumNP
    thS2=FQ(0.0_dp,Par(:,MatNum(N)))
    do L=1,NLay
       If(lWat.Or.PLevel.Eq.0) Then
          SubCha(L)=0.
          SubVol(L)=0.
          hMean(L) =0.
       End If
       SubCO(L)=0.
       COMean(L) =0.
       Ar(L)=0.
    end do
    ATot=0.
    hTotal=0.
    Vol=0.
    Change=0.
    DeltW=0.
    ! temperature
    TTot=0
    TVol=0
    SubT=0
    TMean=0
    If(lCO2) Then
       COTot=0.
       COVol=0.
       DeltCO=0.
    Endif
    If(Transport) Then
       cTot=0.0
       ConVol=0.0
       ConVolIm=0.0
       ConVolIm2=0.0
       cTotIm=0.0
       cVolI=0.0
       ConSub=0.0
       cMean=0.0
       deltC=0.0
    End If
    Do i=N-1,1,-1
       j=i+1
       cEl=0.0
       Mi=MatNum(i)
       Mj=MatNum(j)
       L=LayNum(i)
       dx=coord(j)-coord(i)
       VNewi=dx*(ThN(i)+ThN(j))/2.
       VOldi=dx*(ThO(i)+ThO(j))/2.
       Vol=Vol+VNewi
       Change=Change+(VNewi-VOldi)/dt
       SubCha(L)=SubCha(L)+(VNewi-VOldi)/dt
       SubVol(L)=SubVol(L)+VNewi
       TT=(TempN(i)+TempN(j))/2.+273.15
       if(PLevel.eq.0) then
          WatIn(i)=VNewi
       else
          DeltW=DeltW+abs(WatIn(i)-VNewi)
       end if
       If(lTemp) Then
          TE=(TempN(i)+TempN(j))/2.
          TNewE=0.5*dx*((TempN(i)+273.15)*(TempParam(1,Mi)*TempParam(7,Mi)+&
               TempParam(2,Mi)*TempParam(8,Mi)+TempParam(9,Mi)*ThN(i))+&
               (TempN(j)+273.15)*(TempParam(1,Mj)*TempParam(7,Mj)+&
               TempParam(2,Mj)*TempParam(8,Mj)+TempParam(9,Mj)*ThN(j)))
          TVol=TVol+TNewE
          SubT(L)=SubT(L)+TNewE
          TTot=TTot+TE*dx
          TMean(L)=TMean(L)+TE*dx
       End If
       If(lCO2) Then
          thS1=thS2
          thS2=FQ(0.0_dp,Par(:,Mi))
          Henry=10.**(-13.417+2299.6/(tempN(i)+273.15)+0.01422* &
               (tempN(i)+273.15))/101.3*8.314*(tempN(i)+273.15)
          CONewi=dx*(CO2(j)*(thS1-ThN(j)+Henry*ThN(j))+ &
               CO2(i)  *(thS2-ThN(i)+Henry*ThN(i)))/2.
          COVol=COVol+CONewi
          SubCO(L)=SubCO(L)+CONewi
          COMean(L)=COMean(L)+CO2(i)*dx
          If(PLevel.Eq.0) Then
             COIn(i)=CONewi
          else
             DeltCO=DeltCO+abs(COIn(i)-CONewi)
          End If
       Endif
       If(Transport) Then
          ! solute transport
          Do jS=1,NS
             cE=(Conc(jS,i)+Conc(jS,j))/2.
             TTi=(TempN(i)+273.15-Tr)/R/TT/Tr
             xKsi  =ChPar( 7,Mi,jS)*Exp(TDep( 7,jS)*TTi)
             xNui  =ChPar( 8,Mi,jS)*Exp(TDep( 8,jS)*TTi)
             fExpi =ChPar( 9,Mi,jS)*Exp(TDep( 9,jS)*TTi)
             Henryi=ChPar(10,Mi,jS)*Exp(TDep(10,jS)*TTi)
             TTj=(TempN(j)+273.15-Tr)/R/TT/Tr
             xKsj  =ChPar( 7,Mj,jS)*Exp(TDep( 7,jS)*TTj)
             xNuj  =ChPar( 8,Mj,jS)*Exp(TDep( 8,jS)*TTj)
             fExpj =ChPar( 9,Mj,jS)*Exp(TDep( 9,jS)*TTj)
             Henryj=ChPar(10,Mj,jS)*Exp(TDep(10,jS)*TTj)
             C1=1.0
             C2=1.0
             If(.Not.lLinear(jS)) Then
                If(Conc(jS,i).Gt.0.) C1=Conc(jS,i)**(fExpi-1.)/(1.+xNui*Conc(jS,i)**fExpi)
                If(Conc(jS,j).Gt.0.) C2=Conc(jS,j)**(fExpj-1.)/(1.+xNuj*Conc(jS,j)**fExpj)
             End If
             ThWi=ThN(i)
             ThWj=ThN(j)
             ThImobi=ChPar(4,Mi,jS)
             ThImobj=ChPar(4,Mj,jS)
             ThGi=Max(0.,ths(Mi)-ThWi)
             ThGj=Max(0.,ths(Mj)-ThWj)
             If(iDualPor.Gt.0) Then
                ThImobi=ThNewIm(i)
                ThImobj=ThNewIm(j)
                !              ThGi=amax1(0.,ths(Mi)-ThWi+thSIm(Mi)-ThImobi)
                !              ThGj=amax1(0.,ths(Mj)-ThWj+thSIm(Mj)-ThImobj)
             End If
             If(lMobIm(Mi).And.iDualPor.Eq.0.Or.lBact) ThWi=Max(ThWi-ThImobi,0.001)
             If(lMobIm(Mj).And.iDualPor.Eq.0.Or.lBact) ThWj=max(ThWj-ThImobj,0.001)
             f_em=1.
             If(lDualNEq) f_em=ChPar(13,Mi,jS)
             cNewi=(Conc(jS,i)*(ThN(i)+f_em*ChPar(3,Mi,jS)*ChPar(1,Mi,jS)*xKsi*C1+thGi*Henryi)+&
                    Conc(jS,j)*(ThN(j)+f_em*ChPar(3,Mj,jS)*ChPar(1,Mj,jS)*xKsj*C2+thGj*Henryj))*0.5*dx
             ConVol(jS)=ConVol(jS)+cNewi
             ConSub(jS,L)=ConSub(jS,L)+cNewi
             cTot(jS)=cTot(jS)+cE*dx
             cMean(jS,L)=cMean(jS,L)+cE*dx
             If(jS.Eq.1) cEl=cNewi
          End Do
       End If
       hMean(L)=hMean(L)+hNew(i)*dx
       Ar(L)=Ar(L)+dx
       ATot=ATot+dx
    End Do
    Do L=1,NLay
       hTotal=hTotal+hMean(L)/ATot
       hMean(L)=hMean(L)/Ar(L)
       If(lCO2) Then
          COTot=COTot+COMean(L)/ATot
          COMean(L)=COMean(L)/Ar(L)
          If(Transport) cMean(:,L)=cMean(:,L)/Ar(L)
       End If
    End Do
    If(lCO2) COTot=COTot/ATot
    v1=-(Con1+Con2)/2.*((hNew(2)-hNew(1))/(coord(2)-coord(1))+CosAlf)
    vN=-(ConN+ConM)/2.*((hNew(N)-hNew(N-1))/(coord(N)-coord(N-1))+CosAlf)

    If(output(7)) Then
       write(76,110) t
       write(76,120) (i,i=1,NLay)
       write(76,130)
       write(76,140) ATot,  (    Ar(i),i=1,NLay)
       write(76,150) Vol,(SubVol(i),i=1,NLay)
       write(76,160) Change,(SubCha(i),i=1,NLay)
       write(76,170) hTotal,  ( hMean(i),i=1,NLay)
       If(lCO2) Then
          write(76,180) COVol,( SubCO(i),i=1,NLay)
          Write(76,190) COTot,(COMean(i),i=1,NLay)
       Endif
       If(Transport) Then
          Do jS=1,NS
             Write(76,251) Trim(solute_name(jS))//' volume  [M/L2] ',ConVol(jS),ConSub(jS,:)
             Write(76,251) Trim(solute_name(jS))//' cMean   [M/L3] ',cTot(jS),  cMean(jS,:)
251          Format(1x,a,1P,11e12.4)
          End Do
       End If
       write(76,200) vN,v1
    Endif

    !     Mass balance calculation
    If(PLevel.Eq.0) Then
       wVolI=Vol
       If(lCO2) COVolI=COVol
       If(Transport) cVolI=ConVol
    Else
       wBalT=Vol-wVolI-wCumT
       if(output(7)) write(76,210) wBalT
       ww=max(DeltW,wCumA)
       If(ww.Ge.1.e-25) Then
          wBalR=abs(wBalT)/ww*100.
          if(output(7)) write(76,220) wBalR
       End If
       If(lCO2) Then
          COBalT=COVol-COVolI-COCumT
          if(output(7)) write(76,230) COBalT
          ww=max(DeltCO,COCumA)
          If(ww.Ge.1.e-25) Then
             COBalR=abs(COBalT)/ww*100.
             if(output(7)) write(76,231) COBalR
          End If

          If(Transport) Then
             Do jS=1,NS
                cBalT=ConVol(jS)-cVolI(jS)+cCumT(jS)
                If(.Not.lEquil)       cBalT=cBalT+ConVolIm(jS)
                If(lBact.Or.lDualNEq) cBalT=cBalT+ConVolIm2(jS)
                If(output(7)) Write(76,250) jS,cBalT
                cc=amax1(DeltC,cCumA(jS))
                If(cc.Gt.1.e-25) Then
                   cBalR=Abs(cBalT)/cc*100.
                   If(output(7)) Write(76,260) jS,cBalR
                End If
             End Do
          Endif
       End If
       If(output(7)) Write(76,240) Yield
    End If
    if(output(7)) write(76,130)

110 format(//'-----------------------------------------------------'/ &
         ' Time       [T]',f12.4/ &
         '-----------------------------------------------------')
120 format( ' Sub-region num.             ',9(i7,4x))
130 format( '-----------------------------------------------------')
140 format( ' Area          [L]',1P,9e12.4)
150 format( ' W-volume      [L]',1P,9e12.4)
160 format( ' In-flow     [L/T]',1P,9e12.4)
170 format( ' h Mean        [L]',1P,9e12.4)
180 format( ' CO2-volume    [L]',1P,9e12.4)
190 Format( ' CO2-mean  [L3/L3]',1P,9e12.4)
250 Format( ' ConcVol    [M/L2]',I2,1P,11e12.4)
260 Format( ' cMean      [M/L3]',I2,1P,11e12.4)
200 format( ' Top Flux    [L/T]',1P,e12.4/ &
            ' Bot Flux    [L/T]',1P,e12.4)
210 format( ' WatBalT       [L]',1P,e12.4)
220 format( ' WatBalR       [%]',   f12.4)
230 format( ' CO2BalT       [L]',1P,e12.4)
231 format( ' CO2BalR       [%]',1P,e12.4)
240 format( ' Crop Yield    [%]',   f12.4)

  end subroutine SubReg

  !*********************************************************************

  subroutine NodOut(thN, TPrint)

    Use datatypes
    Use Geometry
    Use carbon, only: CO2,g0
    Use Variables, Only: CosAlf
    Use Solute, Only: solute_name, Conc, SorbOut
    Implicit None 
   
    ! arguments
    Real(dp) :: thN(:), TPrint
    ! local variables
    Integer :: i,j,N,ns
    Real(dp) :: v,dxA,dxB,vA,vB
    Logical, Save :: first_call=.True.

    If(debug) Print *,'NodOut'
    N=NumNP
    ns=Ubound(Conc,1)
    v=-(Con(N)+Con(N-1))/2*((hNew(N)-hNew(N-1))/(coord(N)-coord(N-1))+CosAlf)
    If(first_call) Then
       first_call=.False.
       Write(75,110) '#Node','Depth','Head','Moisture','CO2','Temp','K','C','Flux','Sink','Product',&
            (Trim(solute_name(j))//'_w',Trim(solute_name(j))//'_s',j=1,ns)
       Write(75,110) '#    ','[L]','[L]','[-]','[-]','[C]','[L/T]','[1/L]','[L/T]','[1/T]','[1/T]',&
            ('[M/L**3]','[M/(M soil)]',j=1,ns)
    Endif
    Write(75,'(a,i8)') '#Time:',Nint(TPrint)
    write(75,120) N,coord(N)-zSurf,hNew(N),thN(N),CO2(N),TempN(N),Con(N), &
         Cap(1),v,Sink(1),g0(1),(Conc(j,N),SorbOut(j,N),j=1,ns)
    do i=N-1,2,-1
       dxA=coord(i+1)-coord(i)
       dxB=coord(i)-coord(i-1)
       vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+CosAlf)
       vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+CosAlf)
       v= (vA*dxA+vB*dxB)/(dxA+dxB)
       write(75,120) i,coord(i)-zSurf,hNew(i),thN(i),CO2(i),TempN(i),Con(i), &
            Cap(1),v,Sink(1),g0(1),(Conc(j,i),SorbOut(j,i),j=1,ns)
    end do
    v=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/(coord(2)-coord(1))+CosAlf)
    write(75,120) 1,coord(1)-zSurf,hNew(1),thN(1),CO2(1),TempN(1),Con(1), &
         Cap(1),v,Sink(1),g0(1),(Conc(j,1),SorbOut(j,1),j=1,ns)

110 Format(A5,A9,A10,A9,A14,A8,5A12,99A12)
120 Format(i5,f9.2,f10.1,f9.3,e14.4,f8.3,1P,99e12.3)

  end subroutine NodOut

  
  Subroutine ConcOut(Conc, SorbOut, TPrint)

    Use datatypes
    Use Geometry, Only: coord,zSurf
    Use Solute, Only: solute_name
    Implicit None 
   
    ! arguments
    Real(dp), Intent(in) :: Conc(:,:), SorbOut(:,:), TPrint
    ! local variables
    Integer :: i,j,n,ns
    Logical, Save :: first_call=.True.

    If(debug) Print *,'ConcOut'
    n=Ubound(Conc,2)
    ns=Ubound(Conc,1)
    If(first_call) Then
       first_call=.False.
       Write(86,110) '#Node','Depth',(Trim(solute_name(j))//'_w',Trim(solute_name(j))//'_s',j=1,ns)
       Write(86,110) '#    ','[L]',('[M/L**3]','[M/(M soil)]',j=1,ns)
    Endif
    Write(86,'(A,I8)') '#Time:',Nint(TPrint)
    Do i=n,1,-1
       Write(86,120) i,coord(i)-zSurf,(Conc(j,i),SorbOut(j,i),j=1,ns)
    End Do
110 Format(A5,A10,99A15)
120 Format(0P,i5,1x,f9.2,1P,99E15.7)

  end subroutine ConcOut

  
  Subroutine MassOut(Conc, SorbOut, TPrint)

    Use datatypes
    Use Geometry, Only: delta_z, Matnum, ThNew, NSnit, Pho
    Use Solute, Only: ChPar, solute_name
    Use Nitrogen, Only: rvol, PoolNitMan,PoolNitHum,PoolNitLit, NFact
    Use Phosphorus, Only: PoolPhoActive,PoolPhoStable,PoolPhoLit,PoolPhoMan,PoolPhoHum
    Use Variables, Only: lNitrogen, lPhosphorus
    Implicit None 
   
    ! arguments
    Real(dp), Intent(in) :: Conc(:,:), SorbOut(:,:), TPrint
    ! local variables
    Integer :: i,j,n,ns
    Logical, Save :: first_call=.True.
    Real(dp), Allocatable, Save :: mass_c(:), mass_s(:)
    Real(dp) :: rvolm, mass_N, mass_P, mass_Pools, active, stable, lit, man, hum

    If(debug) Print *,'MassOut'
    n=Ubound(Conc,2)
    ns=Ubound(Conc,1)
    If(first_call) Then
       first_call=.False.
       Allocate(mass_c(ns),mass_s(ns))
       Open(51,file='mass.out',Status='REPLACE')
       Write(51,110) '#Time   ',&
            (Trim(solute_name(j))//'_w',Trim(solute_name(j))//'_s',j=1,ns),&
            'NH3_Prod','Mass_N-conc','Mass_N-Pools','Mass_N','Mass_P',&
            'PoolPhoActive','PoolPhoStable','PoolPhoLit','PoolPhoMan','PoolPhoHum'
       Write(51,110) '#       ',('[M/L**2]',j=1,2*ns),&
            '[M/L**2/T]','[M/L**2]','[M/L**2]','[M/L**2]','[M/L**2]',&
            '[M/L**2]','[M/L**2]','[M/L**2]','[M/L**2]','[M/L**2]'
    Endif
    
    ! nitrogen
    If(lNitrogen) Then
       Do j=1,NSnit
          mass_c(j)=0
          mass_s(j)=0
          Do i=1,n
             ! conc = (mass solute) / (volume water)
             mass_c(j)=mass_c(j)+Conc(j,i)*delta_z(i)*ThNew(i)
             ! sorb = (mass solute) / (mass soil)
             mass_s(j)=mass_s(j)+SorbOut(j,i)*delta_z(i)*ChPar(1,MatNum(i),j)
          End Do
          mass_N=mass_N+NFact(j)*(mass_c(j)+mass_s(j))
       Enddo
       rvolm=dot_Product(rvol,delta_z)
       mass_Pools=dot_Product(delta_z,PoolNitMan)+dot_Product(delta_z,PoolNitHum)&
            +dot_Product(delta_z,PoolNitLit)
    Else
       mass_N=null_dp
       rvolm=null_dp
       mass_Pools=null_dp
    Endif

    ! phosphorus
    If(lPhosphorus) Then
       mass_P=0
       mass_c(Pho)=0
       mass_s(Pho)=0
       If(lPhosphorus) Then
          Do i=1,n
             ! conc = (mass solute) / (volume water)
             mass_c(Pho)=mass_c(Pho)+Conc(Pho,i)*delta_z(i)*ThNew(i)
             ! sorb = (mass solute) / (mass soil)
             mass_s(Pho)=mass_s(Pho)+SorbOut(Pho,i)*delta_z(i)*ChPar(1,MatNum(i),Pho)
          End Do
          active=dot_Product(delta_z,PoolPhoActive)
          stable=dot_Product(delta_z,PoolPhoStable)
          lit=dot_Product(delta_z,PoolPhoLit)
          man=dot_Product(delta_z,PoolPhoMan)
          hum=dot_Product(delta_z,PoolPhoHum)
       Endif
       mass_P=mass_c(Pho)+mass_s(Pho)
    Else
       mass_P=null_dp
       active=null_dp
       stable=null_dp
       lit=null_dp
       man=null_dp
       hum=null_dp
    Endif

    Write(51,120) Nint(Tprint),(mass_c(j),mass_s(j),j=1,ns),&
         rvolm,mass_N,mass_Pools,mass_N+mass_Pools,&
         mass_P,active, stable, lit, man, hum
110 Format(A8,99A14)
120 Format(0P,i8,1P,99ES14.5E3)

  end subroutine MassOut

  !**********************************************************************

  subroutine ObsNod(t, cvTop, output, belowground_respiration)
    Use datatypes
    Use Geometry, only: NObs,Node,ThNew,hNew,TempN,vNew
    Use Carbon, Only: CO2,PoolDPM,PoolRPM,PoolHUM
    Use Variables, Only: aboveground_respiration,NEE,GPP,alphaAvg,lNitrogen
    Use Plants, Only: PlantGeometry
    Use Nitrogen, Only: PoolNitLit,PoolNitMan,PoolNitHUM
    Use Solute, Only: NS, Conc, SorbOut, solute_name
    Implicit None 
    logical :: output(:)
    Real(dp) :: t, cvTop
    Real(dp) :: belowground_respiration ! (added by R.Peters)
    ! local variables
    Integer :: i,j,n
    
    If(debug) Print *,'ObsNode'
    If(output(8)) Then
       Write(77,100) 'Time:',t,'node','theta','head','CO2','temp','v',&
            (Trim(solute_name(j))//'_w',Trim(solute_name(j))//'_s',j=1,NS)
       Do i=1,NObs
          n=Node(i)
          Write(77,111) n,ThNew(n),hNew(n),CO2(n),tempN(n),vNew(n),(Conc(j,n),SorbOut(j,n),J=1,NS)
       End Do
    Endif

100 Format(A,f12.3,/,A5,A8,A11,A9,A8,A11,100A15)
111 Format(0P,I5,f8.4,f11.2,f9.5,f8.3,1P,e11.3,100E15.7)

    ! (added by R.Peters 16.06.2016)
    if(output(12)) then
       If(lNitrogen) Then
          Write(81,'(1P,99999E14.5)') t, &
               (ThNew(Node(i)),i=1,NObs), &
               (hNew(Node(i)),i=1,NObs), &
               (CO2(Node(i)),i=1,NObs), &
               (tempN(Node(i)),i=1,NObs), &
               aboveground_respiration, &
               belowground_respiration, &
               NEE, &
               GPP, &
               alphaAvg, &
               PlantGeometry%rootdepth, &
               PlantGeometry%slaig, &
               PlantGeometry%wlvg, &
               PlantGeometry%wso, &
               PlantGeometry%wst, &
               PlantGeometry%wrt, &
               (PoolDPM(Node(i)),i=1,NObs), &
               (PoolRPM(Node(i)),i=1,NObs), &
               (PoolHUM(Node(i)),i=1,NObs), &
               (PoolNitLit(Node(i)),i=1,NObs), &
               (PoolNitMan(Node(i)),i=1,NObs), &
               (PoolNitHUM(Node(i)),i=1,NObs), &
               ((Conc(j,Node(i)),j=1,NS),i=1,NObs)
       Else
          Write(81,'(1P,99999E14.5)') t, &
               (ThNew(Node(i)),i=1,NObs), &
               (hNew(Node(i)),i=1,NObs), &
               (CO2(Node(i)),i=1,NObs), &
               (tempN(Node(i)),i=1,NObs), &
               aboveground_respiration, &
               belowground_respiration, &
               NEE, &
               GPP, &
               alphaAvg, &
               PlantGeometry%rootdepth, &
               PlantGeometry%slaig, &
               PlantGeometry%wlvg, &
               PlantGeometry%wso, &
               PlantGeometry%wst, &
               PlantGeometry%wrt, &
               (PoolDPM(Node(i)),i=1,NObs), &
               (PoolRPM(Node(i)),i=1,NObs), &
               (PoolHUM(Node(i)),i=1,NObs), &
               ((Conc(j,Node(i)),j=1,NS),i=1,NObs)
       Endif
    endif

    if(output(13)) write(82,'(1P,99999E14.5)') &
         t, (tempN(Node(i)),i=1,NObs), cvTop

  end subroutine ObsNod


  !**********************************************************************

  Subroutine PoolOut(TPrint, output)

    Use datatypes
    Use Carbon
    Use Geometry, Only: NumNP,hNew,thNew,coord,tempN,zSurf
    Use Nitrogen, Only: PoolNitMan,PoolNitHum,PoolNitLit
    Use Phosphorus, Only: PoolPhoActive,PoolPhoStable,PoolPhoLit,PoolPhoMan,PoolPhoHum
    Use Variables, Only: lNitrogen, lPhosphorus
    Implicit None 
    ! arguments
    logical :: output(:)
    Real(dp) :: TPrint
    ! local variables
    Integer :: i
    Logical, Save :: first_call=.True.
    
    if(KProd.ne.3) return
    
    If(output(9)) Then
       If(first_call) Then
          first_call=.False.
          Write(78,110) '#Node','Depth','Head','Moisture','temp', &
               'CO2','Product','DPM','RPM','BIO','HUM','IOM','CO2-C',&
               'NitMan','NitHum','NitLit',&
                'PoolPhoActive','PoolPhoStable','PoolPhoLit','PoolPhoMan','PoolPhoHum'
          Write(78,110) '#    ',   '[L]', '[L]',     '[-]', &
               '[-]',  '[1/T]','M/L3','M/L3','M/L3','M/L3', &
               'M/L3','M/L3',' ','M/L3','M/L3','M/L3', &
               'M/L3','M/L3','M/L3','M/L3','M/L3'
       Endif
       Write(78,'(a,i8)') '#Time:',Nint(TPrint)
       do i=1,NumNP
          If(lNitrogen) Then
             If(lPhosphorus) Then
                Write(78,120) i,coord(i)-zSurf,hNew(i),thNew(i),tempN(i), &
                     CO2(i),g0(i),PoolDPM(i),PoolRPM(i),PoolBio(i),PoolHum(i), &
                     PoolIOM(i),PoolCO2(i),PoolNitMan(i),PoolNitHum(i),PoolNitLit(i), &
                     PoolPhoActive(i),PoolPhoStable(i),&
                     PoolPhoLit(i),PoolPhoMan(i),PoolPhoHum(i)
             Else
                Write(78,120) i,coord(i)-zSurf,hNew(i),thNew(i),tempN(i), &
                     CO2(i),g0(i),PoolDPM(i),PoolRPM(i),PoolBio(i),PoolHum(i), &
                     PoolIOM(i),PoolCO2(i),PoolNitMan(i),PoolNitHum(i),PoolNitLit(i)
             Endif
          Else
             Write(78,120) i,coord(i)-zSurf,hNew(i),thNew(i),tempN(i), &
                  CO2(i),g0(i),PoolDPM(i),PoolRPM(i),PoolBio(i),PoolHum(i), &
                  PoolIOM(i),PoolCO2(i)
          End If
       end do
    Endif
    PoolCO2=0.0
       
110 format(A5,  A10,  A16,  2A9,  99A16)
120 Format(i5,f10.1,1P,e16.8,0P,f9.5,f9.3,1P,99e16.8)
  end subroutine PoolOut

  !**********************************************************************

  subroutine ProdOut(t)

    Use datatypes
    Use Geometry, only: NumNP,coord,zSurf
    Use Carbon
    Implicit None 
    ! arguments
    Real(dp) :: t
    ! local variables
    Integer :: i

    if(KProd.ne.3) return

    write(80,111) ' Time:',t,' '
111 format(//A,f13.6//A)
    write(80,110) 'Node','Depth', &
         'CO2total','CO2DPM','CO2RPM','CO2BIO','CO2HUM', 'CO2root'
    write(80,110)   ' ',   '[L]',      &
         '[1/T]','[1/T]','[1/T]','[1/T]','[1/T]'
    do i=1,NumNP
       write(80,120) i,coord(i)-zSurf,g0(i), &
            gamSDPM(i),gamSRPM(i),gamSBio(i),gamSHum(i), g0r(i)
    end do
    write(80,'(''end'')')

110 format(A5,  A10,  6A16)
120 format(i5,f10.3,6e16.8)
  end subroutine ProdOut


end Module Output
