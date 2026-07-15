Module Watflow

  Use datatypes
  Implicit None

Contains

  subroutine calculateWatFlow(t,dt,P,R,S,Iter,ItCum,ConvgF)

    Use datatypes
    Use Geometry
    Use Variables
    Use Material, only: hSat,thR,thS
    Use Timedata, only: decreaseTimeStep
    Implicit None 
    logical :: ConvgF
    Integer :: Iter,ItCum
    Real(dp) :: t,dt
    Real(dp), dimension(:) :: P,R,S
    ! local variables
    Integer :: i,m
    Logical :: ItCrit
    Real(dp) :: EpsH,EpsTh,PB,RB,SB,PT,RT,ST,Th
    Character(len=128) :: text
    Iter=0
    ConvgF=.true.


    PicardIteration: do
       !     Generate terms of matrix equation and solv by Gauss elimination
       call SetMat(ThNew)
       call Reset (NumNP,rTop,rBot,CosAlf,dt,coord,hOld,Con,Cap,WLayer,hNew, &
            Sink,P,R,S,PB,RB,SB,PT,RT,ST,FreeD,qGWLF,Aqh,Bqh, &
            GWL0L,ThNew,ThOld)
       call Shift (KodTop,rTop,rBot,hTop,hBot,hCritA,hCritS,CosAlf, &
            WLayer,Con(NumNP),Con(NumNP-1),hNew(NumNP), &
            hNew(NumNP-1),coord(NumNP),coord(NumNP-1),AtmBC,KodBot,Con(1), &
            Con(2),hNew(1),hNew(2),coord(1),coord(2),SeepF,ThNew(NumNP), &
            ThOld(NumNP),Sink(NumNP),dt)
       hHelp=hNew
       call Gauss (NumNP,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT, &
            RT,ST)
       Iter =Iter+1
       ItCum=ItCum+1

       ! Test for convergence
       ItCrit=.true.
       do i=1,NumNP
          m=MatNum(i)
          EpsTh=0.
          EpsH=0.
          if(hHelp(i).lt.hSat(m).and.hNew(i).lt.hSat(m)) then
             Th=ThNew(i)+Cap(i)*(hNew(i)-hHelp(i))/(ths(m)-thr(m))
             EpsTh=abs(ThNew(i)-Th)
             !write(*,*) i, ": ", Cap(i), " * ", (hNew(i)-hHelp(i)), " * ", (1/(ths(m)-thr(m))), " = ", EpsTh
          else
             EpsH=abs(hNew(i)-hHelp(i))
          end if
          !write(*,*) "Eps(", i, "): ", NumNP, " : ", EpsTh, " - ", EpsH
          if(EpsTh.gt.TolTh.or.EpsH.gt.TolH) then
             ItCrit=.false.
             Exit
          end if
       end do

       if(ItCrit) Exit   ! convergence

       If(Iter>=MaxIt) Then
          Write(text,'(a,1p,e12.5,a,e12.5,a)') 't=',t,' dt=',dt,&
               " maximum number of iterations reached, decrease time step."
          call WriteWarning(text)
          ConvgF=decreaseTimeStep(t,dt)
          if(ConvgF) then
             do i=1,NumNP
                hNew(i) =hOld(i)
                hHelp(i)=hOld(i)
             end do
             KodTop=KTOld
             KodBot=KBOld
             Iter=0
          else
             write(*,*) "can not decrease time step"
             return
          end if
       end if
    end do PicardIteration
  end subroutine calculateWatFlow

  !***********************************************************************

  subroutine Reset(N,rTop,rBot,CosAlf,dt,x,hOld,Con,Cap,WLayer,hNew, &
       Sink,P,R,S,PB,RB,SB,PT,RT,ST,FreeD,qGWLF,Aqh,Bqh, &
       GWL0L,ThNew,ThOld)

    Implicit None
    Integer :: N
    logical WLayer,FreeD,qGWLF
    Real(dp) :: rTop,rBot,CosAlf,dt,PB,RB,SB,PT,RT,ST,Aqh,Bqh,GWL0L
    Real(dp), Dimension(:) :: x,hOld,hNew,P,R,S,Con,Cap,Sink,ThNew,ThOld
    ! local variables
    Integer :: i
    Real(dp) :: dx,dxA,dxB,A2,A3,B,F2,ConA,ConB

    !     Finite differences

    dxB=x(2)-x(1)
    ConB=(Con(1)+Con(2))/2.
    S(1)=-ConB/dxB
    RB=ConB/dxB
    SB=-RB
    if(qGWLF) rBot=Fqh(hNew(1)-GWL0L,Aqh,Bqh)
    PB=ConB*CosAlf+rBot
    if(FreeD) PB=0
    do i=2,N-1
       dxA=x(i)-x(i-1)
       dxB=x(i+1)-x(i)
       dx=(dxA+dxB)/2.
       ConA=(Con(i)+Con(i-1))/2.
       ConB=(Con(i)+Con(i+1))/2.
       A2=ConA/dxA+ConB/dxB
       A3=-ConB/dxB
       B =(ConA-ConB)*CosAlf
       F2=Cap(i)*dx
       R(i)=A2+F2/dt
       P(i)=F2*hNew(i)/dt-(ThNew(i)-ThOld(i))*dx/dt-B-Sink(i)*dx
       S(i)=A3
    end do
    dxA=x(N)-x(N-1)
    dx=dxA/2.
    ConA=(Con(N)+Con(N-1))/2.
    RT=ConA/dxA+Cap(N)*dx/dt
    ST=-ConA/dxA
    PT=Cap(N)*dx*hNew(N)/dt-(ThNew(N)-ThOld(N))*dx/dt-Sink(N)*dx- &
         ConA*CosAlf-rTop
    if(WLayer) then
       if(hNew(N).gt.0.) then
          RT=RT+1./dt
          PT=PT+max(hOld(N),0.0d0)/dt
       else
          PT=PT+max(hOld(N),0.0d0)/dt
       end if
    end if

  end subroutine Reset

  !**********************************************************************

  subroutine Gauss(N,KodTop,KodBot,hTop,hBot,hNew,P,R,S,PB,RB,SB,PT,RT,ST)

    Implicit None 
    Integer :: N,KodTop,KodBot
    Real(dp) ::hTop,hBot,PB,RB,SB,PT,RT,ST
    Real(dp), Parameter :: eps=1e-20
    Real(dp), Dimension(:) :: hNew,P,R,S
    ! local variables
    Integer :: i

    !     Forward
    if(KodBot.ge.0) then
       P(2)=P(2)-S(1)*hBot
    else
       P(2)=P(2)-PB*S(1)/RB
       R(2)=R(2)-SB*S(1)/RB
    end if
    do i=3,N-1
       If(Abs(R(i-1)) > eps) Then
          P(i)=P(i)-P(i-1)*S(i-1)/R(i-1)
          R(i)=R(i)-S(i-1)*S(i-1)/R(i-1)
       Else
          Print *,'watflow gauss exception 1'
       Endif
    end do
    if(KodTop.ge.0) then
       P(N-1)=P(N-1)-S(N-1)*hTop
    else
       P(N)=PT-P(N-1)*ST/R(N-1)
       R(N)=RT-S(N-1)*ST/R(N-1)
    end if

    !     Back
    if(KodTop.ge.0) then
       hNew(N)=hTop
       hNew(N-1)=P(N-1)/R(N-1)
    else
       hNew(N)=P(N)/R(N)
       hNew(N-1)=(P(N-1)-S(N-1)*hNew(N))/R(N-1)
    end if
    do i=N-2,2,-1
       If(Abs(R(i)) > eps) Then
          hNew(i)=(P(i)-S(i)*hNew(i+1))/R(i)
       Else
          Print *,'watflow gauss exception 2'
       Endif
    end do
    if(KodBot.ge.0) then
       hNew(1)=hBot
    else
       hNew(1)=(PB-SB*hNew(2))/RB
    end if

  end subroutine Gauss

  !**********************************************************************

  subroutine Shift(KodTop,rTop,rBot,hTop,hBot,hCritA,hCritS,CosAlf, &
       WLayer,ConTop,ConBlw,hNT,hNB,xTop,xBlw,AtmBC, &
       KodBot,ConBot,ConAbv,hN1,hN2,xBot,xAbv,SeepF, &
       ThNew,ThOld,SinkTop,dt)

    Implicit None
    Integer :: KodTop,KodBot
    logical WLayer,AtmBC,SeepF
    Real(dp) :: rTop,rBot,hTop,hBot,hCritA,hCritS,CosAlf, &
         ConTop,ConBlw,hNT,hNB,xTop,xBlw, &
         ConBot,ConAbv,hN1,hN2,xBot,xAbv, &
         ThNew,ThOld,SinkTop,dt
    Real(dp) :: vBot,vTop

    if(SeepF) then
       ! Seepage face at the bottom
       if(KodBot.ge.0) then
          vBot=-(ConBot+ConAbv)/2.*((hN2-hN1)/(xAbv-xBot)+CosAlf)
          if(vBot.gt.0.) then
             KodBot=-2
             rBot=0.
          end if
       else
          if(hN1.ge.0.) then
             KodBot=2
             hBot=0.
          end if
       end if
    Else if(KodBot.eq.6) then
       ! bottom boundary condition: pressure head / zero flux into domain
       vBot=-0.5*(ConBot+ConAbv)*((hN2-hN1)/(xAbv-xBot)+CosAlf)
       if(vBot.gt.0.) then
          KodBot=-6
          rBot=0.
       end if
    else if(KodBot.eq.-6) then
       if(hN1.ge.hBot) then
          KodBot=6
       end if
    end if

    !     Atmospheric boundary condition
    if(AtmBC) then
       if(KodTop.ge.0) then
          vTop=-(ConTop+ConBlw)/2.*((hNT-hNB)/(xTop-xBlw)+CosAlf)- &
               (ThNew-ThOld)*(xTop-xBlw)/2./dt-SinkTop*(xTop-xBlw)/2.
          if(abs(vTop).gt.abs(rTop).or.vTop*rTop.le.0) KodTop=-4
       else
          if(.not.WLayer) then
             if(hNT.gt.hCritS) then
                KodTop=4
                hTop=hCritS
             end if
          end if
          if(hNT.le.hCritA) then
             KodTop=4
             hTop=hCritA
          end if
       end if
    end if

  end subroutine Shift

  !***********************************************************************

  subroutine SetMat(theta)

    Use Geometry, only: NumNP,MatNum,hNew,hHelp,Con,Cap
    Use Material, only: NTab,hTab,ConTab,CapTab,Par,ConSat,hSat,TheTab,thS,FK,FC,FQ
    Implicit None
    Real(dp), Intent(out) :: theta(:)
    Integer :: i,M,iT
    Real(dp), Parameter :: small=1.e-15
    Real(dp) :: alh1,dlh,hi1,hi2,hiM,Ci,S1,S2,Ti

    alh1=log10(-hTab(1))
    dlh =(log10(-hTab(NTab))-alh1)/(NTab-1)
    do i=1,NumNP
       M=MatNum(i)
       hi1=min(hSat(M),hHelp(i))
       hi2=min(hSat(M), hNew(i))
       hiM=0.1*hi1+0.9*hi2
       if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
          Ci=ConSat(M)
       else if(hiM.gt.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=min(int((log10(-hiM)-alh1)/dlh)+1, NTab-1)
          S1=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=ConTab(iT,M)+S1*(hiM-hTab(iT))
       else
          Ci=FK(hiM,Par(:,M))
       end if
       Con(i)=Ci
       Con(i)=max(Con(i),small)
       hiM=hNew(i)
       if(hiM.ge.hSat(M)) then
          Ci=0.
          Ti=ThS(M)
       else if(hiM.ge.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=min(int((log10(-hiM)-alh1)/dlh)+1, NTab-1)
          S1=(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
          S2=(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=CapTab(iT,M)+S1*(hiM-hTab(iT))
          Ti=TheTab(iT,M)+S2*(hiM-hTab(iT))
       else
          Ci=FC(hiM,Par(:,M))
          Ti=FQ(hiM,Par(:,M))
       end if
       Cap(i)=max(Ci,small)
       theta(i)=Ti
    end do

  end subroutine SetMat


  !**********************************************************************
  function Fqh(GWL,Aqh,Bqh)
    Real(dp) Fqh,GWL,Aqh,Bqh
    Fqh=Aqh*exp(Bqh*abs(GWL))
  end function Fqh


  !***********************************************************************
  !     To calculate the velocities
  subroutine Veloc(N,hNew,Con,x,CosAlf,v,ThNew,ThOld,SinkTop,dt)

    Implicit None
    Integer :: N
    Real(dp) :: hNew(:),x(:),Con(:),v(:),CosAlf,ThNew,ThOld,SinkTop,dt
    Integer :: i
    Real(dp) :: dxA,dxB,vA,vB

    v(N)=-(Con(N)+Con(N-1))/2.*((hNew(N)-hNew(N-1))/(x(N)-x(N-1))+ &
         CosAlf)-(x(N)-x(N-1))/2.*((ThNew-ThOld)/dt+SinkTop)
    do i=2,N-1
       dxA=x(i+1)-x(i)
       dxB=x(i)-x(i-1)
       vA=-(Con(i)+Con(i+1))/2.*((hNew(i+1)-hNew(i))/dxA+CosAlf)
       vB=-(Con(i)+Con(i-1))/2.*((hNew(i)-hNew(i-1))/dxB+CosAlf)
       v(i)=(vA*dxB+vB*dxA)/(dxA+dxB)
    end do
    v(1)=-(Con(1)+Con(2))/2.*((hNew(2)-hNew(1))/(x(2)-x(1))+CosAlf)

  end subroutine Veloc

end Module Watflow
