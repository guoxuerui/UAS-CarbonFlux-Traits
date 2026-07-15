Module Timedata

  Use datatypes
  Implicit None
  Save

  ! simulation time data
  Type SimTimeType
     Real(dp) :: time
     Real(dp) :: delta_t
     Integer :: LastSimDay   ! last simulation day
     Logical :: NewSimDay    ! new simulation day
  End Type SimTimeType
  Type(SimTimeType) SimTime

  ! time information
  Real(dp) :: tInit,tMax,tAtm,tAtmOld,dtMax,DMul,DMul2,dtMin,dtOpt,tOld
  Real(dp) :: co2dt, co2dtMin
  Real(dp) :: resetPlantTime=0

  Integer :: StartDate   ! in iday format
  Real(dp), Parameter :: yearLength=365.2425_dp   !365, 242190517
  Integer, Parameter :: daysOfMonth(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  Real(dp) :: dtInit
  Logical :: lMinStep

  Type ProfileType
     ! 1 Profile = data(2, NoPoints)
     Integer :: NoPoints
     Real(dp), Pointer :: data(:,:)
     Real(dp) :: InterpolData
  End Type ProfileType

  ! NoProfiles time profiles
  Type(ProfileType), Allocatable :: Profile(:)

  ! private data
  Real(dp), Private :: OutTimestep
  Real(dp) :: EpsTime
  Real(dp) :: TimeFactor=1
  Integer :: AtmPerDay
Contains

  !-------------------------------------------
  ! Simulation time routines
  !-------------------------------------------
  Subroutine InitSimTime(t,dt)
    Implicit None
    Real(dp) :: t,dt
    EpsTime=Max(eps,0.001*dtMin)
    SimTime%time=t
    SimTime%delta_t=dt
    SimTime%LastSimDay=-1      ! last simulation day
    SimTime%NewSimDay=.true.   ! new simulation day
    OutTimestep=1
    ! CO2 time control
    co2dt=dt
    co2dtMin=0
    AtmPerDay=Nint(1.0/TimeFactor)
    lMinstep=.True.
  End Subroutine InitSimTime


  Function decreaseTimeStep(t,dt)
    Implicit None
    Logical :: decreaseTimeStep
    Real(dp) :: t,dt

    If(dt<=dtMin) Then
       decreaseTimeStep=.false.
    else
       dt=max(dt/3,dtMin)
       dtOpt=dt
       t=tOld+dt
       SimTime%time=t
       SimTime%delta_t=dt
       decreaseTimeStep=.true.
    end if
  End Function decreaseTimeStep


  Subroutine TmCont(t,dt,Iter,tPrint,NumNP,x,ThNew,Disp,MatNum,thS,&
       dtMaxC,Transport)
    implicit None
    Integer :: Iter,NumNP, MatNum(:)
    Real(dp) :: t,dt,tPrint,dtMaxC
    Real(dp), Dimension(:) :: x, ThNew, Disp, thS
    Logical :: Transport
    ! local variables
    Integer :: i, M
    Real(dp), Parameter :: eps=1d-5
    Real(dp) :: tFix,th_air1,th_air2,maxdt
    tOld=t
    if(lMinStep) then
       maxdt=Min(dtMax,dtInit,dtOpt)
       If(Transport) maxdt=Min(maxdt,dtMaxC)
       dtOpt=maxdt
       lMinStep=.false.
    else
       If(Transport) Then
          maxdt=Min(dtMax,dtMaxC)
       Else
          maxdt=dtMax
       Endif
    end if
    tFix=min(tPrint,tAtm,tMax)
    if(Iter.le.3.and.(tFix-t).ge.DMul*dtOpt) &
         dtOpt=min(maxdt,DMul*dtOpt)
    if(Iter.ge.7) &
         dtOpt=max(dtMin,DMul2*dtOpt)
    dt=Min(dtOpt,tFix-t)
    dt=min((tFix-t)/anint((tFix-t)/dt),maxdt)
    If(Transport) dt=Min(dt,dtMaxC) ! apply dtMaxC from solute transport
    if( tFix-(t+dt) .lt. dtMin/100 ) dt = tFix -t   !added by D.Farber

    ! calculate dt for CO2 transport
    if(co2dtMin.gt.0.0d0) then
       co2dt=dt
       Do i=2,NumNP
          M=MatNum(i)
          th_air1=thS(M)-ThNew(i)+eps
          M=MatNum(i-1)
          th_air2=thS(M)-ThNew(i-1)+eps
          co2dt=Min( co2dt, 0.5d0 * (x(i-1)-x(i))**2 / &
               Min(Disp(i-1),Disp(i)) * Max(th_air1, th_air2) )
       enddo
       co2dt=Max(co2dt,co2dtMin)
    endif
    ! calculating SimTime
    SimTime%LastSimDay=Int(SimTime%time+EpsTime)
    t=t+dt
    SimTime%time=t
    SimTime%delta_t=dt
    ! needed to call the plant routines once per day (... per defined timeperiod; D.Farber)
    If(Int(SimTime%time+EpsTime)>SimTime%LastSimDay) Then
       SimTime%NewSimDay=.true.
    Else
       SimTime%NewSimDay=.false.
    End If
  end subroutine TmCont


  Function SavePlantOutput(time)
    Implicit None
    Logical SavePlantOutput
    Real(dp), Intent(in) :: time
    If(OutTimestep<=null_dp) Then
       SavePlantOutput=.false.
    Else If(Mod(time+EpsTime,OutTimestep)<0.5_dp*dtMin) Then
       SavePlantOutput=.true.
    Else
       SavePlantOutput=.false.
    End If
  End Function SavePlantOutput


  !-------------------------------------------
  ! date routines
  !-------------------------------------------
  Integer Function date_to_iday(y, m, d)
  ! date format: number of timesteps from year 1.1.0001
    Implicit None
    Integer, Intent(in) :: y, m ,d
    Integer :: y1, i, iday

    iday=Int(y*yearLength+eps) + d
    y1=y-1
    Do i=1,m-1
       iday=iday+daysOfMonth(i)
    End Do
    If(m>2) iday=iday+leap_year(y)
    date_to_iday=int(iday/TimeFactor+eps)
  End Function date_to_iday


  Integer Function simtime_to_iday(simtime)
    Implicit None
    Real(dp), Intent(in) :: simtime
    simtime_to_iday=Int(simtime+EpsTime)+StartDate
  End Function simtime_to_iday


  Subroutine iday_to_date(iTime, y, m, d)
    Implicit None
    Integer, Intent(in) :: iTime
    Integer, Intent(out) :: y, m, d
    Integer :: i, iday
    !changed by N.Prolingheuer (correct results but wrong way) maybe correct now (D.Farber)
    iday = int(iTime*TimeFactor+eps)
    y=Int(Real(iday,dp)/yearLength+eps)

    d=iday-Int(y*yearLength)
    Do m=1,12
       i=d-daysOfMonth(m)

       If(m==2) i=i-leap_year(y)
       If(i<1) Exit
       d=i
    End Do
  End Subroutine iday_to_date


  Integer Function day_of_year(iTime)
    Implicit None
    Integer, Intent(in) :: iTime
    Integer :: y, iday
    iday = int(iTime*TimeFactor+eps)
    y=Int(Real(iday,dp)/yearLength+eps)
    day_of_year=iday-Int(y*yearLength+eps)
  End Function day_of_year


  Integer Function leap_year(y)
    Implicit None
    Integer, Intent(in) :: y
    If( (Mod(y,4)==0 .and. Mod(y,100)/=0) .or. Mod(y,400)==0 ) Then
       leap_year=1
    Else
       leap_year=0
    End If
  End Function leap_year


  Subroutine SetStartDate(y, m, d)
    Implicit None
    Integer, Intent(in) :: y, m ,d
    StartDate=date_to_iday(y, m, d)-1   ! -1 because simtime starts with 1 =>  simtime = 1: "startDate"+simTime = correct startDate
  End Subroutine SetStartDate


  Real(dp) Function date_to_simtime(y, m, d)
    Implicit None
    Integer, Intent(in) :: y, m ,d
    date_to_simtime=date_to_iday(y, m, d) - StartDate
  End Function date_to_simtime

End Module Timedata
