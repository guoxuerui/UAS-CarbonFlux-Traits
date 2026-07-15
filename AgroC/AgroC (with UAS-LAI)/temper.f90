Module Temperature

  ! Calculation of heat transport
  Use datatypes
  Implicit None

Contains
  
  Subroutine Temper(t,dt,B,D,E,F,Cap,Cond)
    Use Geometry, Only: NumNP,coord,MatNum,TempO,TempN,vOld,vNew,ThOld,ThNew,Sink
    Use Variables, Only: Ampl,tPeriod,kTopT,kBotT,tTop,tBot
    Use Material, Only: tempParam
    Implicit None
    Real(dp), Intent(in) :: t,dt
    Real(dp), Dimension(:), Intent(inout) :: B,D,E,F,Cap,Cond
    Integer :: N,i,j,Level,M
    Real(dp) :: Cw,dx,dxA,dxB,v,th,tTopA

    N=NumNP
    Cw=tempParam(9,1)
    tTopA=tTop
    If(tPeriod.Gt.0.) tTopA=tTop+Ampl*Sin(2.*PI*t/tPeriod-7.*PI/12.)
    TempO=TempN

    Do Level=1,2
       Do i=1,N
          M=MatNum(i)
          If(Level.Eq.1) Then
             th=ThOld(i)
             v=vOld(i)
          Else
             th=ThNew(i)
             v=vNew(i)
          End If
          Cap(i)=tempParam(7,M)*tempParam(2,M)+tempParam(8,M)*tempParam(3,M)+tempParam(9,M)*th
          Cond(i)=tempParam(4,M)+tempParam(5,M)*th+tempParam(6,M)*Sqrt(th)+ &
               tempParam(9,M)*tempParam(1,M)*Abs(v)
       End Do
       dx=coord(2)-coord(1)
       If(kBotT.Gt.0) Then
          D(1)=1.
          E(1)=0.
          F(1)=tBot
       Else If(kBotT.Lt.0) Then
          If(Level.Eq.2) Then
             D(1)=dx/12./dt*(3.*Cap(1)+Cap(2))+(Cond(1)+Cond(2))/dx/4.+ &
                  Cw*(2.*vNew(1)+vNew(2))/12.+ &
                  dx/24.*Cw*(3.*Sink(1)+Sink(2))
             E(1)=-(Cond(1)+Cond(2))/4./dx+dx/12./dt*(Cap(1)+Cap(2))+ &
                  Cw*(2.*vNew(2)+vNew(1))/12.+dx/24.*Cw*(Sink(1)+Sink(2))
          Else
             F(1)=TempO(1)*(dx/12./dt*(3.*Cap(1)+Cap(2))- &
                  (Cond(1)+Cond(2))/dx/4.- &
                  Cw*(2.*vOld(1)+vOld(2))/12.- &
                  dx/24.*Cw*(3.*Sink(1)+Sink(2)))+ &
                  TempO(2)*((Cond(1)+Cond(2))/4./dx+ &
                  dx/12./dt*(Cap(1)+Cap(2))- &
                  Cw*(2.*vOld(2)+vOld(1))/12.- &
                  dx/24.*Cw*(Sink(1)+Sink(2)))+ &
                  tBot*Cw*(vNew(1)+vOld(1))/2.
          End If
       Else
          D(1)=-1.
          E(1)=1.
          F(1)=0.
       End If
       Do i=2,N-1
          dxA=coord(i)-coord(i-1)
          dxB=coord(i+1)-coord(i)
          dx=(coord(i+1)-coord(i-1))/2.
          If(Level.Eq.2) Then
             B(i)=-(Cond(i)+Cond(i-1))/4./dxA-Cw*(vNew(i)+ &
                  2.*vNew(i-1))/12.+dxA/12./dt*(Cap(i-1)+Cap(i))+ &
                  dxA/24.*Cw*(Sink(i-1)+Sink(i))
             D(i)=(Cond(i-1)+Cond(i))/4./dxA+(Cond(i)+Cond(i+1))/4./dxB+ &
                  (dxA*(Cap(i-1)+3.*Cap(i))+ &
                  dxB*(Cap(i+1)+3.*Cap(i)))/12./dt+ &
                  Cw*(vNew(i+1)-vNew(i-1))/12.+ &
                  dxA/24.*Cw*(Sink(i-1)+3.*Sink(i))+ &
                  dxB/24.*Cw*(3.*Sink(i)+Sink(i+1))
             E(i)=-(Cond(i)+Cond(i+1))/4./dxB+ &
                  dxB/12./dt*(Cap(i+1)+Cap(i))+ &
                  Cw*(2.*vNew(i+1)+vNew(i))/12.+ &
                  dxB/24*Cw*(Sink(i+1)+Sink(i))
          Else
             F(i)=TempO(i-1)*(dxA/12./dt*(Cap(i-1)+Cap(i))+ &
                  (Cond(i)+Cond(i-1))/4./dxA+ &
                  Cw*(vOld(i)+2.*vOld(i-1))/12.- &
                  dxA/24.*Cw*(Sink(i-1)+Sink(i)))+ &
                  TempO(i)*(-Cw*(vOld(i+1)-vOld(i-1))/12.+ &
                  (dxA*(Cap(i-1)+3.*Cap(i))+ &
                  dxB*(3.*Cap(i)+Cap(i+1)))/12./dt- &
                  (Cond(i+1)+Cond(i))/4./dxB- &
                  (Cond(i)+Cond(i-1))/4./dxA- &
                  dxA/24.*Cw*(Sink(i-1)+3.*Sink(i))- &
                  dxB/24.*Cw*(3.*Sink(i)+Sink(i+1)))+ &
                  TempO(i+1)*(dxB/12./dt*(Cap(i)+Cap(i+1))+ &
                  (Cond(i+1)+Cond(i))/4./dxB- &
                  Cw*(2.*vOld(i+1)+vOld(i))/12.- &
                  dxB/24.*Cw*(Sink(i+1)+Sink(i)))
          End If
       End Do
       If(kTopT.Gt.0) Then
          B(N)=0.
          D(N)=1.
          F(N)=tTopA
       Else If(kTopT.Lt.0) Then
          dx=coord(N)-coord(N-1)
          If(Level.Eq.2) Then
             B(N)=-(Cond(N)+Cond(N-1))/4./dx- &
                  Cw*(vNew(N)+2.*vNew(N-1))/12.+ &
                  dx/12./dt*(Cap(N-1)+Cap(N))+ &
                  dx/24.*Cw*(Sink(N-1)+Sink(N))
             D(N)=dx/12./dt*(Cap(N-1)+3.*Cap(N))+ &
                  (Cond(N-1)+Cond(N))/4./dx- &
                  Cw*(2.*vNew(N)+vNew(N-1))/12.+ &
                  dx/24.*Cw*(Sink(N-1)+3.*Sink(N))
          Else
             F(N)=TempO(N-1)*((Cond(N)+Cond(N-1))/4./dx+ &
                  Cw*(vOld(N)+2.*vOld(N-1))/12.+ &
                  dx/12./dt*(Cap(N-1)+Cap(N))- &
                  dx/24.*Cw*(Sink(N-1)+Sink(N)))+ &
                  TempO(N)*(dx/12./dt*(Cap(N-1)+3.*Cap(N))- &
                  (Cond(N-1)+Cond(N))/4./dx+ &
                  Cw*(2.*vOld(N)+vOld(N-1))/12.- &
                  dx/24.*Cw*(Sink(N-1)+3.*Sink(N)))- &
                  tTopA*Cw*(vNew(N)+vOld(N))/2.
          End If
       End If
    End Do

    !     Solve matrix equation
    Do i=2,N
       D(i)=D(i)-B(i)*E(i-1)/D(i-1)
       F(i)=F(i)-B(i)*F(i-1)/D(i-1)
    End Do
    TempN(N)=F(N)/D(N)
    Do i=2,N
       j=N-i+1
       TempN(j)=(F(j)-E(j)*TempN(j+1))/D(j)
    End Do

  End Subroutine Temper

End Module Temperature
