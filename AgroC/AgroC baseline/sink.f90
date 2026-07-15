Module SourceSink

  Use datatypes

Contains

  Subroutine SetSnk(hRoot,vRoot,TPot,cRoot)

    Use Geometry, only: NumNP,MatNum,coord,Sink,Beta,hNew
    Use Variables, only: P0,H50,fR6,lCO2
    Use Carbon, only: CO2
    Implicit None
    Real(dp) :: hRoot,vRoot,TPot,cRoot
    Integer :: i,M
    Real(dp) :: aRoot,dxM,alfa

    vRoot=0.
    hRoot=0.
    ARoot=0.
    cRoot=0.
    do i=2,NumNP
       if(Beta(i).gt.0.) then
          if(i.eq.NumNP) then
             dxM=(coord(i)-coord(i-1))/2.
          else
             dxM=(coord(i+1)-coord(i-1))/2.
          end if
          M=MatNum(i)
          if(hNew(i).ge.0.) then
             Alfa=0.
          else
             Alfa=1./(1.+(hNew(i)/H50)**P0)
          end if
          Sink(i)=Alfa*Beta(i)*TPot*fR6
          vRoot=vRoot+Sink(i)*dxM
          hRoot=hRoot+hNew(i)*dxM
          if(lCO2) cRoot=cRoot+CO2(i)*dxM
          ARoot=ARoot+dxM
       else
          Sink(i)=0.
       end if
    end do
    if(ARoot.gt.0.001) then
       hRoot=hRoot/ARoot
       cRoot=cRoot/ARoot
    end if

  end subroutine SetSnk

  !***********************************************************************

  subroutine SetRG(t,CumT,fET)

    Use Geometry, only: NumNP,coord,Beta
    Use Variables, only: kRoot,kBeta,tRMin,tRHarv,xRMin,xRMax,RGR,RDDMax,alpha
    Implicit None
    Real(dp) :: t,CumT,fET
    Integer :: i
    Real(dp) :: tt,xR,coordN,SBeta

    if(t.lt.tRMin.or.t.gt.tRHarv) then
       fET=0.
       CumT=0.
       Beta=0.
       return
    end if
    if(kRoot.eq.1) then
       tt=CumT/RDDMax
       if(tt.gt.1.) tt=1.
    else
       xR=xRMax
       if(xRMin.le.0.00001) xRMin=0.00001
       tt=t-tRMin
    end if
    xR=(xRMax*xRMin)/(xRMin+(xRMax-xRMin)*exp(-RGR*tt))
    fET=xR/xRMax
    SBeta=0.
    coordN=coord(NumNP)
    do i=2,NumNP
       if(kBeta.eq.0) then
          if(coord(i).lt.coordN-xR) then
             Beta(i)=0.
          else if(coord(i).lt.coordN-0.2*xR) then
             Beta(i)=2.08333/xR*(1-(coordN-coord(i))/xR)
          else
             Beta(i)=1.66667/xR
          end if
       else
          Beta(i)=0.
          if(coord(i).gt.coordN-xR) Beta(i)=alpha*exp(-alpha*(coordN-coord(i)))
       end if
       if(i.ne.NumNP) then
          SBeta=SBeta+Beta(i)*(coord(i+1)-coord(i-1))/2.
       else
          SBeta=SBeta+Beta(i)*(coord(i)-coord(i-1))/2.
       end if
    end do
    if(SBeta.le.0.) then
       Beta(NumNP-1)=1
    else
       do i=2,NumNP
          Beta(i)=Beta(i)/SBeta
       end do
    end if

  end subroutine SetRG

End Module SourceSink
