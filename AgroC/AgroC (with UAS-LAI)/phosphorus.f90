Module Phosphorus
  Use datatypes
  Use Variables, Only: lPhosphorus
  Use Geometry, Only: Pho
  Implicit None
  Logical :: lPpools, lPRootUpt
  ! 0=calcareous,1=slighlty weathered,2=moderately weathered,3=highly weathered
  Integer :: i,m,SoilGen
  Integer :: nfrom, nto, nby, ndim
  Real(dp), Allocatable :: CaC(:)   ! CACO3 conc., for SoilGen=0
  Real(dp), Allocatable :: BSA(:)   ! base saturation NH4 acetate in %, for SoilGen=2
  Real(dp), Allocatable :: pH(:)    ! soil pH, for SoilGen=2
  Real(dp) :: PhosRo ! [-] C/P ratio of the biomass(nicht sicher ob wir das brauchen)
  ! Phosphorus humification coefficient
  ! default EPIC = 0.2
  ! negative Value indicates the Use of the RothC Value= 0.54
  Real(dp) :: PhosFh ! [-]
  ! humus P mineralizaton rate constant [d-1]
  ! negative Value indicates the Use of the RothC kHUM
  Real(dp) :: RkCMPi  ! [1/T]
  ! D0 travel distance resistance between bulk soil and root for Phosphorus
  Real(dp) :: Pho_rdo ! [1/L]
  Real(dp) :: P_conv_diff_factor ! reduce convective uptake to factor*(diffusive uptake)
  Real(dp) :: P_Fert_lab ! fraction of labile (soluble) P in org. fertilizer
  Real(dp) :: P_Fert_C   ! fraction of P to C (used to calculate P in org. fertilizer)
  Integer :: P_Fert_num ! Number of mineral P fertilizer applications on top
  ! Mineral P fertilizer application on top
  ! 1 line per time: time, mass_P
  ! (dpm,rpm,hum)
  Real(dp), Allocatable :: P_init_Lit(:),P_init_Man(:),P_init_Hum(:)
  Real(dp), Allocatable :: PoolPhoLit(:),PoolPhoMan(:),PoolPhoHum(:)
  Real(dp), Allocatable :: PoolPhoActive(:),PoolPhoStable(:)
  Real(dp), Allocatable :: PhoSinkFert(:),unc(:),und(:)
  ! boundary for solute transport
  Real(dp) :: P_on_top
  Integer :: PTopDim, PTopInd=1
  Real(dp), Allocatable :: PTopTime(:), PTopMass(:)
  Type BalanceType
     Real(dp) :: PhoSinkSum
     Real(dp) :: r_active_labile_sum
     Real(dp) :: r_active_stable_sum
     Real(dp) :: sorp
  End Type BalanceType
  Type(BalanceType) :: Balance

  Type PhoSinkBalanceType
     Real(dp) :: dvs,temp_sum,xnc_tot,PlantWeight,apc_tot,&
          tunc,totdem,tund,totvr,ap_tot,PhoReduct,&
          sinksum,sinkrootsum
  End Type PhoSinkBalanceType
  Type(PhoSinkBalanceType) :: PhoSinkBalance
Contains

  Subroutine PhosphorusIn(output)
    Use Carbon, Only: co2alf, co2rhum
    Use Geometry, Only: NumNP, MatNum, coord, PhoSink, PhoSinkRoot
    Use Material, Only: NMat

    Implicit None
    Logical, Intent(in) :: output(:)
    Integer :: i, m, nfrom,nto,nby,ndim
    Real(dp) :: exp_alt, exp_neu, w
    Real(dp), Allocatable :: p_init(:,:)

    If(debug) Print *,'PhosphorusIn'
    If(.Not.lPhosphorus) Then
       Print *,'Skipping block J'
       If(headline(1:12) == '*** BLOCK J:') Call SkipBlock('J')
       Return
    Endif
    Write(*,*) 'reading phosphorus input'
    read(30,*)
    Read(30,*) lPpools, lPRootUpt, SoilGen
    If(SoilGen<0 .Or. SoilGen>3) Call WriteError('Invalid input for SoilGen')
    Read(30,*)
    Allocate(CaC(NMat), BSA(nMat), pH(nMat))
    Do m=1,NMat
       Read(30,*) CaC(m), BSA(m), pH(m)
    Enddo
    Read(30,*) 
    Read(30,*) PhosRo
    Read(30,*) PhosFh
    If(PhosFh<0.0) PhosFh=0.54_dp
    Read(30,*)
    Read(30,*) RkCMPi
    If(RkCMPi<0.0) Then
       RkCMPi=co2rhum
    Endif
   
    Read(30,*) Pho_rdo
    Read(30,*) P_conv_diff_factor
    Read(30,*)
    Read(30,*) P_Fert_lab
    Read(30,*) P_Fert_C
    Read(30,*) 
    ! read fertilizer input
    Read(30,*) PTopDim
    Allocate(PTopTime(PTopDim),PTopMass(PTopDim))
    Read(30,*)
    Do m=1,PTopDim
       Read(30,*) PTopTime(m),PTopMass(m)
    Enddo
    ! read initial pools
    if(co2alf.gt.null_dp) then
       nfrom=1
       nto=1
       nby=1
       ndim=1
    Else If(equal(co2alf, -1.0_dp)) Then
       nfrom=NumNP
       nto=1
       nby=-1
       ndim=NumNP
    else 
       nfrom=1
       nto=NMat
       nby=1
       ndim=NMat
    endif
    Allocate(p_init(5,ndim))
    Allocate(PoolPhoLit(NumNP),PoolPhoMan(NumNP),PoolPhoHum(NumNP))
    Allocate(PoolPhoActive(NumNP),PoolPhoStable(NumNP))
    Allocate(PhoSinkFert(NumNP))
    Allocate(unc(NumNP))
    Allocate(und(NumNP))
    PhoSinkFert=0
    PoolPhoActive=0
    PoolPhoStable=0
    Read(30,*)
    Do m=nfrom,nto,nby
       Read(30,*) p_init(:,m)
    Enddo
    Read(30,1) headline
    ! init pools
    If(co2alf.Le.null_dp) Then
       ! equally distribute init pool material on the nodes
       Do i=1,NumNP
          If(equal(co2alf, -1.0_dp)) Then
             m=i
          Else 
             m=MatNum(i)
          Endif
          PoolPhoLit(i)=p_init(1,m)
          PoolPhoMan(i)=p_init(2,m)
          PoolPhoHum(i)=p_init(3,m)
          PoolPhoActive(i)=p_init(4,m)
          PoolPhoStable(i)=p_init(5,m)
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
          PoolPhoLit(i)=p_init(1,1)*w
          PoolPhoMan(i)=p_init(2,1)*w
          PoolPhoHum(i)=p_init(3,1)*w
          PoolPhoActive(i)=p_init(4,1)*w
          PoolPhoStable(i)=p_init(5,1)*w
       End Do
       w=exp_neu
       !write(*,*) w,coord(NumNP),coord(1)
       PoolPhoLit(1)=p_init(1,1)*w
       PoolPhoMan(1)=p_init(2,1)*w
       PoolPhoHum(1)=p_init(3,1)*w
       PoolPhoActive(1)=p_init(4,1)*w
       PoolPhoStable(1)=p_init(5,1)*w
    Endif
    ! output
    If(output(1)) Then
       Write(50,'(/A)') 'Phosphorus input'
       Write(50,'(A/)') '=============='
       Write(50,*) 'lPpools: ',lPpools
       Write(50,*) 'lPRootUpt: ',lPRootUpt
       Write(50,*) 'Soilgen: ',SoilGen
       Write(50,3) 'CAC BSA NGH4:',(CaC(m),BSA(m),pH(m),m=1,NMat)
       Write(50,1) 'PhosRo: ',PhosRo
       Write(50,1) 'PhosFh: ',PhosFh
       Write(50,1) 'P_Fert_lab P_Fert_C: ',P_Fert_lab,P_Fert_C
       Write(50,2) 'P_Fertilizer: ',(PTopTime(m),PTopMass(m),m=1,PTopDim)
       Write(50,3) 'n_init:',p_init
    End If
    Deallocate(p_init)
    Allocate(PhoSink(NumNP))
    Allocate(PhoSinkRoot(NumNP))
    PhoSink=null_dp
    PhoSinkRoot=null_dp
    Call reset_PhoSinkBalance()
1   Format(A,1P,3E13.5)
2   Format(A,/,1P,(2E13.5))
3   Format(A,/,1P,(3E13.5))
  End Subroutine PhosphorusIn

  Subroutine WritePhosphorusOutput(SimTime)
    Use TimeData, Only: iday_to_date, simtime_to_iday
    Implicit None
    Real(dp), Intent(in) :: SimTime
    Logical, Save :: firstCall=.True.
    Integer :: ios,yy,mm,dd
    Character(len=*), Parameter :: fn_phos = "phosphorus.out"

    If(debug) Print *,'WritePhosphorusOutput'
    If(.Not.lPhosphorus) Return
    If(firstCall) Then
       firstCall=.false.
       Open(87,file=Trim(fn_phos),Status='REPLACE',ACTION='WRITE',&
            FORM='FORMATTED',iostat=ios)
       If(ios/=0) Call WriteError( 'Error while opening nitrogen output file: '&
            // Trim(fn_phos))
       Write(87,2) '#phophorus'
!            '# 5=PhoSink  (M/(L2 soil)/T)',&
       Write(87,2) '# 1=year','# 2=month','# 3=day','# 4=simulation time',&
            '# 5=rActLab (M/L2/T)',&
            '# 6=rActStab (M/L2/T)',&
            '# 7=sorp (-)',&
            '# 8=dvs (-)',&
            '# 9=temp_sum ( degree*T )',&
            '#10=xnc_tot ( (mass Phosphorus)/(mass dry matter) )',&
            '#11=Plantweight ( (mass dry matter )/L2) ) ',&
            '#12=apc_tot ( (mass Phosphorus)/(mass dry matter) )',&
            '#13=tunc ( (mass Phosphorus)/L2/T )',&
            '#14=totdem ( (mass Phosphorus)/L2/T )',&
            '#15=tund ( (mass Phosphorus)/L2/T )',&
            '#16=totvr ( (mass Phosphorus)/L2/T )',&
            '#17=ap_tot ( (mass Phosphorus)/L2 )',&
            '#18=PhoReduct (-)',&
            '#19=PhoSink  (M/(L2 soil)/T)',&
            '#20=PhoSinkRoot  (M/(L2 soil)/T)'
       Write(87,'(A,A,99A13)') '#     date','    time',&
       'rActLab','rActStab','sorp',&
       'dvs','temp_sum','xnc_tot','PlantWeight','apc_tot','tunc','totdem',&
       'tund','totvr','ap_tot','PhoReduct','PhoSink','PhoSinkRoot'
    Else
       Open(87,file=fn_phos,Status='OLD',ACTION='WRITE',&
            FORM='FORMATTED',POSITION='APPEND')
    Endif
    Call iday_to_date(simtime_to_iday(SimTime),yy,mm,dd)
    Write(87,'(I4,2I3,I8,1P,99E13.5)') yy,mm,dd, Nint(SimTime),&
         Balance%r_active_labile_sum,&
         Balance%r_active_stable_sum,Balance%sorp,&
         PhoSinkBalance%dvs, PhoSinkBalance%temp_sum,&
         PhoSinkBalance%xnc_tot,PhoSinkBalance%PlantWeight,&
         PhoSinkBalance%apc_tot,PhoSinkBalance%tunc,&
         PhoSinkBalance%totdem,PhoSinkBalance%tund,&
         PhoSinkBalance%totvr,PhoSinkBalance%ap_tot,&
         PhoSinkBalance%PhoReduct,&
         PhoSinkBalance%sinksum,PhoSinkBalance%sinkrootsum
    PhoSinkBalance%sinksum=null_dp
    PhoSinkBalance%sinkrootsum=null_dp
2   Format(A)
    Close(87)
  End Subroutine WritePhosphorusOutput


  Subroutine p_sink(dt)
    Implicit None
    Real(dp), Intent(in) :: dt
    If(debug) Print *,'p_sink'
    If(lPhosphorus) Then
       Call P_uptake(dt)
       Call Mineral_P_Cycling(dt)
       Call P_Mineralization(dt)
       Call p_correct_sink(dt)
    End If
  End Subroutine p_sink

  Subroutine Mineral_P_Cycling(dt)
    ! inorganic pools: labile (conc), active and stable
    Use Geometry, Only: NumNP,MatNum,delta_z,ThNew,PhoSink
    Use Solute, Only: Conc
    Use Carbon, Only: co2parm
    Use TimeData, Only: TimeFactor
    Implicit None
    Real(dp), Intent(in) :: dt
    Integer :: m,n
    Real(dp) :: bo,delta_t,vol_w
    Real(dp) :: sorp  ! sorption coefficient 
    Real(dp) :: r_active_labile  ! mineral P flow rate between labile and active in kg/ha/d
    Real(dp) :: r_active_stable ! mineral P flow rate between active and stable in kg/ha/d
    Real(dp) :: active_pool ! active mineral in kg/ha
    Real(dp) :: labile_pool ! labile mineral pool
    Real(dp) :: stabile_pool ! stable mineral pool
    Real(dp) :: PhoSinkSum, r_active_labile_sum, r_active_stable_sum,sorp_sum
    Real(dp) :: labile_pool_sum

    If(debug) Print *,'Mineral_P_Cycling'
    delta_t=dt*TimeFactor
    PhoSinkSum=0
    r_active_labile_sum=0
    r_active_stable_sum=0
    sorp_sum=0
    labile_pool_sum=0
    Do n=1,NumNP
       vol_w=delta_z(n)*ThNew(n)
       labile_pool=Conc(Pho,n)*vol_w ! conc->mass
       labile_pool_sum=labile_pool_sum+labile_pool*delta_z(n)
    Enddo
    Do n=1,NumNP
       m=MatNum(n)
       vol_w=delta_z(n)*ThNew(n)
       labile_pool=Conc(Pho,n)*vol_w ! conc->mass
       active_pool=PoolPhoActive(n)
       stabile_pool=PoolPhoStable(n)
       Select Case(SoilGen)
       Case(0) ! PSP=0.58-0.0061*CAC
          sorp=0.58-0.0061*CAC(m)
       Case(1) ! PSP=0.02+0.014*WPML
          ! to_ha(Unit_L_Input): convert from agroc length unit to ha
          ! Units_M(3)/Units_M(Unit_M_Input): convert from agroc mass unit to kg
          sorp=0.02+0.014*labile_pool_sum*to_ha(Unit_L_Input)*Units_M(3)/Units_M(Unit_M_Input)
       Case(2) ! PSP=0.0054*BSA+0.116*PH-0.73
          sorp=0.0054*BSA(m)+0.116*Ph(m)-0.73
       Case(3) ! PSP=0.46-0.0916*ln(CLA)
          sorp=0.46-0.0916*Log(co2parm(1,m)*100_dp) ! clay
       Case default
          Call WriteError('Invalid value of SoilGen')
       End Select

       ! labile to active
       !Print *,'sorp=',sorp
       sorp=Max(0.05_dp,Min(0.75_dp,sorp))
       sorp_sum=sorp_sum+sorp
       ! MPR=WPML-WPMA*PSP/(1.0-PSP)  Eq238
       r_active_labile=labile_pool-active_pool*sorp/(1.0_dp-sorp)*delta_t
       !Print *,'r_active_labile',r_active_labile,labile_pool,active_pool,sorp
       If(r_active_labile<null_dp) &
            r_active_labile=0.1*r_active_labile ! reverse flow is much lower
       r_active_labile_sum=r_active_labile_sum+r_active_labile*delta_z(n)
       active_pool=active_pool+r_active_labile
       ! set sink term for solute transport
       PhoSink(n)=PhoSinkFert(n)-r_active_labile/delta_t
       PhoSinkSum=PhoSinkSum+PhoSink(n)*delta_z(n)
       PhoSinkFert(n)=0
       ! stable to active
       If(SoilGen==0) Then
          bo=0.0076
       Else
          bo=exp(-1.77*sorp-7.05)
       End If
       ! ASPR=bo*(4.0*WPMA-WPMS) Eq243
       r_active_stable=bo*(4.0*active_pool-stabile_pool)*delta_t
       If(r_active_stable<null_dp) &
            r_active_stable=0.1*r_active_stable ! reverse flow is much lower
       r_active_stable_sum=r_active_stable_sum+r_active_stable*delta_z(n)
       stabile_pool=stabile_pool+r_active_stable
       active_pool=active_pool-r_active_stable
       ! copy back to pool
       PoolPhoActive(n)=active_pool
       PoolPhoStable(n)=stabile_pool
    End Do
    sorp_sum=sorp_sum/NumNP
    Call set_balance(PhoSinkSum,r_active_labile_sum,r_active_stable_sum,sorp_sum)
  End Subroutine Mineral_P_Cycling

  Subroutine P_Mineralization(dt)
    ! organic pools: Lit (DPM), Man (RPM) and Hum
    Use datatypes
    Use Geometry, Only: NumNP, MatNum, ThNew, delta_z, PhoSink
    Use Carbon, Only: PoolDPM, PoolRPM, sfact, co2parm, co2rdpm,co2rrpm,&
                      eff_rate_lit_from_nit,eff_rate_man_from_nit
    Use Solute, Only: Conc

    Implicit None
    Real(dp), Intent(in) :: dt
    Integer :: i
    Real(dp) :: eff_rate_hum, eff_rate_lit, eff_rate_man
    Real(dp) :: check_mim, fe
    Real(dp) :: pxl, pxm, red_fact, shortage, vimma, vimmb
    Real(dp) :: change_p, dx !, plus

    If(debug) Print *,'P_Mineralization'
    If(.Not.LPPools) Return
    !     pxl= ratio defining mineralisation or immobilisation for litter
    !     pxm= ratio defining mineralisation or immobilisation for manure
    !          > 0 = mineralisation
    !          < 0 = immobilisation
    Do i=1,NumNP
       dx=delta_z(i)
       red_fact = sfact(1,i)*sfact(2,i)*sfact(3,i)
       fe=co2parm(2,MatNum(i)) ! dimensionless
       !       If(PoolDPM(i)<1e-10) Then
       !          Print *,i,PoolNitLit(i),PoolDPM(i),PoolNitLit(i)/PoolDPM(i),PoolNitMan(i)/PoolRPM(i),fe/ro
       !       Endif
       pxl=PoolPhoLit(i)/PoolDPM(i)-(fe/PhosRo) ! dimensionless
       pxm=PoolPhoMan(i)/PoolRPM(i)-(fe/PhosRo) ! dimensionless
       !Print *,i,pxl,pxm
       eff_rate_lit = co2rdpm*red_fact ! 1/T
       eff_rate_man = co2rrpm*red_fact
       eff_rate_hum = RkCMPi*red_fact
       ! unit = unit of Pool = M/(L3 soil)
       check_mim = ( eff_rate_hum*PoolPhoHum(i)+ &
            pxm*eff_rate_man*PoolRPM(i)+ &
            pxl*eff_rate_lit*PoolDPM(i) )*dt ! M/(L3 soil)
       !If(i==1) Write(*,'(a,1p,99e13.5)') 'pxl,pxm,chek_mim=',pxl,pxm,check_mim,eff_rate_lit, eff_rate_man,eff_rate_hum,&
       !     PoolNitHum(i),PoolRPM(i),PoolDPM(i),red_fact, sfact(1,i),sfact(2,i),sfact(3,i)
       If(check_mim .Ge. 0.d0) Then ! P Mineralization, no immobilization
          !Print *,'miner_immob then'
          change_p = check_mim ! unit of Pool = M/(L3 soil)
       Else ! Immobilization occurs
          If(Conc(Pho,i)*ThNew(i) .Ge. (-check_mim)) Then
             change_p = check_mim
          Else ! The present P concentration is not enough to fulfill mineralization P demand
             change_p = -Conc(Pho,i)*ThNew(i)
             shortage=check_mim+Conc(Pho,i)*ThNew(i) ! negativ
             If(pxl>null_dp  .And. pxm<null_dp) Then   ! Lit pool has enough P, Man pool is missing P for mineralization
                vimmb    = pxm*eff_rate_man*PoolRPM(i)*dt  ! this is the amount of P potentially decomposed from Man(RPM) pool ! unit of Pool = M/(L3 soil)
                vimma    = vimmb-shortage            ! this is the amount of P that could ACTUALLY be decomposed from Man(RPM) pool accounting for the P gap
                eff_rate_lit = co2rdpm*red_fact      ! Lit pool decomposition rate not affected
                eff_rate_man = vimma/(pxm*PoolRPM(i)*dt) ! Man(RPM) decomp. rate is computed as ratio between P that could actually be decomposed and amount of P potentially decomposed    
             Elseif(pxl<null_dp  .And. pxm>null_dp) Then   ! Man pool has enough P, Lit pool is missing P for mineralization
                vimmb    = pxl*eff_rate_lit*PoolDPM(i)*dt  ! this is the amount of P potentially decomposed from Lit(DPM) pool
                vimma    = vimmb-shortage             ! this is the amount of P that could ACTUALLY be decomposed from Lit(DPM) pool accounting for the P gap
                eff_rate_lit = vimma/(pxl*PoolDPM(i)*dt) ! Lit(DPM) decomp. rate is computed as ratio between P that could actually be decomposed and amount of P potentially decomposed
                eff_rate_man = co2rrpm*red_fact          ! Man pool decomposition rate not affected  ! 1/T  
             Elseif(pxl<null_dp  .And. pxm<null_dp) Then   ! Lit pool and Man pool is missing P for mineralization
                vimmb    = pxl*eff_rate_lit*PoolDPM(i)*dt  ! this is the amount of P potentially decomposed from Lit(DPM) pool
                vimma    = vimmb-shortage   ! this is the amount of P that could ACTUALLY be decomposed from Lit(DPM) pool accounting for the P gap
                If(vimma .Gt. 0.d0) Then   ! Not?? Enough P to fulfill Lit(DPM) decomposition P demand, may be there is enough NO3 left for Man decomposition?
                   vimmb    = pxm*eff_rate_man*PoolRPM(i)*dt ! this is the amount of P potentially decomposed from Man(RPM) pool ! unit of Pool
                   vimma    = vimmb+vimma ! ??
                   eff_rate_lit = 0.d0  ! No Lit(DPM) decomposition
                   eff_rate_man = vimma/(pxm*PoolRPM(i)*dt) ! Man decomp. rate is reduced ??
                Else ! NOT enough NO3 to fulfill Lit(DPM) decomposition N demand
                   Print *,'troubles in P_Mineralization1'
                   eff_rate_lit = 0.0_dp
                   eff_rate_man = 0.0_dp
                Endif
             Else
                Print *,'troubles in P_Mineralization2'
                eff_rate_lit = 0.0_dp
                eff_rate_man = 0.0_dp
             Endif
          Endif
       Endif
       ! transfer to rates to carbon
       eff_rate_lit_from_nit(i)=Min(eff_rate_lit_from_nit(i),eff_rate_lit)
       eff_rate_man_from_nit(i)=Min(eff_rate_man_from_nit(i),eff_rate_man)
       !        the gain/losses of the phosphorus litter pool
       !        through mineralisation/immobilisation and phosphorus humification
       PoolPhoLit(i)=PoolPhoLit(i)+(-pxl-fe*PhosFh/PhosRo)*eff_rate_lit*PoolDPM(i)*dt  ! unit of Pool = M/(L3 soil)
       !        the gain/losses of the phosphorus litter pool
       !        through mineralisation/immobilisation and phosphorus humification
       PoolPhoMan(i)=PoolPhoMan(i)+(-pxm-fe*PhosFh/PhosRo)*eff_rate_man*PoolRPM(i)*dt
       !plus=plus+(-pxm-fe*fh/ro)*eff_rate_man*PoolRPM(i)*dt
       !        the changes in the n-humus pool through mineralisation and n-litter humification
       PoolPhoHum(i)=PoolPhoHum(i)+(fe*PhosFh/PhosRo*(eff_rate_lit*PoolDPM(i)+ &
            eff_rate_man*PoolRPM(i))-eff_rate_hum*PoolPhoHum(i))*dt
       !     eff_rate_man*PoolRPM(i))-eff_rate_hum*PoolNitHum(i))*dt
       !        retain the rate of change of the ammonia and nitrate species
       !        to be incorporated in the sequential transformation process  (mg/(day*liter))
       PhoSink(i)=PhoSink(i)+change_p/dt ! M/(L3 soil)/T add mineralization to sink term
    Enddo
  End Subroutine P_Mineralization
  
  Subroutine P_uptake(dt)
    Use Geometry, Only: NumNP, ThNew, MatNum, cRootMax, delta_z, Sink, Pho, PhoSink, PhoSinkRoot
    Use Solute, Only: ChPar, Conc
    Use Plants, Only: SucrosPlant, InterpolateTab, Plant,&
         PlantGeometry, PlantGeometryArray, Potatoe, SugarBeet, PhoReduct
    Use Material, Only: ThS
    Use Variables, Only: rorad, UnitFactorPlants, w0_dens, P_uptake_method, &
         aplv, apst, aprt, apso, apcrn, ap_tot, Pho_ancrt

    
    Implicit None
    Real(dp), Intent(in) :: dt
    Real(dp) :: xnc_tot
    Real(dp) :: totdem
    Real(dp) :: apc_tot=0
    Real(dp) :: dvs, temp_sum, tunc, hcsolo, tund, tpdnup, diffus_rm, Dw, TauW, rdens
    Real(dp) :: totvr, PlantWeight
    Integer :: ip,ptype, m, i
    Integer, Save :: oldseason=0
    Real(dp), Save :: totup=0
    Real(dp), Parameter :: rlncl=0.005

    If(debug) Print *,'nit_upt'
  ! uptake from plants
    If(.Not.lPRootUpt) Return
    ip=PlantGeometry%index
    ptype=Plant(ip)%type ! plant type
    If(Plant(ip)%season/=oldseason) Then
       oldseason=Plant(ip)%season
       ! s_ancl = 0.05 in wave
       ap_tot=0
       totup=0
       Call reset_PhoSinkBalance()
       PhoReduct=1
       PhoSink=0
       PhoSinkRoot=0
    Endif
    If(ip==0 .Or. Plant(ip)%season==0) Return
    dvs=PlantGeometry%dvs
    temp_sum=PlantGeometry%temp_sum
    If(P_uptake_method==1) Then
       ! Wave
       If(ptype==Potatoe .Or. ptype==SugarBeet) Then
          xnc_tot=InterpolateTab(SucrosPlant(ip)%tab(18),temp_sum)
       Else
          xnc_tot=InterpolateTab(SucrosPlant(ip)%tab(18),dvs)
       End If
    Else
       ! Apex_Epic
       If(ptype==Potatoe .Or. ptype==SugarBeet) &
            Call WriteError('P_uptake_method 2 not usable for Potatoe or SugarBeet.')
       xnc_tot=SucrosPlant(ip)%bp(1)+SucrosPlant(ip)%bp(2) &
            *Exp(-SucrosPlant(ip)%bp(3)*Min(2.0_dp, dvs)/2.0_dp)
    Endif
    PlantWeight=PlantGeometry%wlvg+PlantGeometry%wst+PlantGeometry%wrt+PlantGeometry%wso+PlantGeometry%wcrn
    If(PlantWeight>null_dp) Then
       apc_tot=ap_tot/PlantWeight
    Else
       apc_tot=null_dp
    Endif
    aprt=PlantGeometry%wrt*apc_tot
    apst=PlantGeometry%wst*apc_tot
    aplv=PlantGeometry%wlvg*apc_tot
    apso=PlantGeometry%wso*apc_tot
    apcrn=PlantGeometry%wcrn*apc_tot
    totdem=Max(0.0,PlantWeight*xnc_tot-ap_tot) ! kg P / L2
    tunc=0.0
    tund=0.0
    Do i=1,NumNP
       ! calculate convection
       hcsolo = Max(0.0,Min(Conc(Pho,i),cRootMax(Pho))) ! M / (L3 water) 
       unc(i) = Max(0.0,Sink(i))*hcsolo*dt*delta_z(i) ! dt*dx*(unit of NitSinkRoot)
       tunc=tunc+unc(i)
       ! calculate diffusion
       m=MatNum(i)
       Dw=ChPar(5,m,Pho) !TODO: noch nicht Temperatur abhaengig
       TauW=ThNew(i)**(7./3.)/ThS(m)**2
       hcsolo = Max(0.0,Conc(Pho,i))*ThNew(i)
       diffus_rm = Dw*TauW
       rdens=PlantGeometryArray%rrd(NumNP-i+1)*w0_dens
       und(i)=rdens*rorad*2*PI*diffus_rm*hcsolo*dt*delta_z(i)/Pho_rdo
       tund=tund+und(i)
    Enddo
    If (tund>totdem) Then
       ! reduce diffusion
       und=und*(totdem/tund)
       tund=totdem
       ! set convection to 0
       unc =0.0 
       tunc=0.0
    Else
       tpdnup=totdem-tund
       If(tunc>tpdnup) Then
          ! reduce convection
          unc=unc*(tpdnup/tunc)
          tunc=tpdnup
       Endif
    Endif
    ! reduce convective uptake
    If(tunc>null_dp) Then
       tpdnup=Min(tunc, P_conv_diff_factor*(tunc+tund))
       unc=unc*(tpdnup/tunc)
       tunc=tpdnup
    Endif
    !     calculate the total uptake (mg/m**2) (wave)
    totvr=tunc+tund
    Pho_ancrt=apc_tot
    ap_tot=Max(null_dp,  (ap_tot+totvr) - PlantGeometry%ratwrtd*Pho_ancrt) ! accumulate
!    PhoReduct = pho_reduct(totvr,totdem)
    PhoReduct = pho_reduct(apc_tot,xnc_tot)
    Call set_PhoSinkBalance(dvs,temp_sum,xnc_tot,PlantWeight,apc_tot,&
         tunc,totdem,tund,totvr,ap_tot,PhoReduct)
    ! calculate the uptake rates
    Do i=1,NumNP
       PhoSinkRoot(i)=(unc(i)+und(i))*UnitFactorPlants/(delta_z(i)*dt)
    Enddo
  End Subroutine P_uptake

  Subroutine set_PhoSinkBalance(dvs,temp_sum,xnc_tot,PlantWeight,apc_tot,&
       tunc,totdem,tund,totvr,ap_tot,PhoReduct)
    Implicit None
    Real(dp), Intent(in) :: dvs,temp_sum,xnc_tot,PlantWeight,apc_tot,&
         tunc,totdem,tund,totvr,ap_tot,PhoReduct

    PhoSinkBalance%dvs=dvs
    PhoSinkBalance%temp_sum=temp_sum
    PhoSinkBalance%xnc_tot=xnc_tot
    PhoSinkBalance%PlantWeight=PlantWeight
    PhoSinkBalance%apc_tot=apc_tot
    PhoSinkBalance%tunc=tunc
    PhoSinkBalance%totdem=totdem
    PhoSinkBalance%tund=tund
    PhoSinkBalance%totvr=totvr
    PhoSinkBalance%ap_tot=ap_tot
    PhoSinkBalance%PhoReduct=PhoReduct
  End Subroutine set_PhoSinkBalance

  Subroutine reset_PhoSinkBalance()
    Implicit None
    PhoSinkBalance%dvs=null_dp
    PhoSinkBalance%temp_sum=null_dp
    PhoSinkBalance%xnc_tot=null_dp
    PhoSinkBalance%PlantWeight=null_dp
    PhoSinkBalance%apc_tot=null_dp
    PhoSinkBalance%tunc=null_dp
    PhoSinkBalance%totdem=null_dp
    PhoSinkBalance%tund=null_dp
    PhoSinkBalance%totvr=null_dp
    PhoSinkBalance%ap_tot=null_dp
    PhoSinkBalance%PhoReduct=one_dp
    PhoSinkBalance%sinksum=null_dp
    PhoSinkBalance%sinkrootsum=null_dp
  End Subroutine reset_PhoSinkBalance

  Real(dp) Function pho_reduct(totvr,totdem)
    Use Plants, Only: PlantGeometry, Plant, Potatoe, SugarBeet
    Implicit None
    Real(dp), Intent(in) :: totvr,totdem
    Real(dp) :: PUptRat
    Integer :: ip, ptype
    
    If(debug) Print *,'pho_reduct'
    ip=PlantGeometry%index
    If(Plant(ip)%season==0 .Or. ip==0 .Or. .Not. lPRootUpt) Then
       pho_reduct=1.0
    Elseif (lPhosphorus) Then
       ptype=Plant(ip)%type   ! plant type
       If(ptype==Potatoe .Or. ptype==SugarBeet) Then
          If(PlantGeometry%temp_sum.Le.450.0 .And. PlantGeometry%plai .Lt. 0.75) Then
             pho_reduct=1.0
          else
             If(totdem<eps) Then
                pho_reduct=one_dp
             Else
                PUptRat=200.0_dp*totvr/totdem
                pho_reduct=PUptRat/(PUptRat+exp(4.065_dp-0.0535_dp*PUptRat))
                pho_reduct=Min(1.0, Max(0.0,totvr/totdem))
             Endif
          endif
       !        winter wheat, spring wheat and maize
       else
          If(PlantGeometry%dvs.Le.0.3 .And. PlantGeometry%plai.Lt.0.75) Then
             pho_reduct=1.0
          Else
             If(totdem<eps) Then
                pho_reduct=one_dp
             Else
                PUptRat=200.0_dp*totvr/totdem
                pho_reduct=PUptRat/(PUptRat+exp(4.065_dp-0.0535_dp*PUptRat))
                pho_reduct=Min(1.0, Max(0.0,totvr/totdem))
             Endif
          endif
       endif
    Else
       pho_reduct=1.0
    Endif
  End Function pho_reduct

  Subroutine p_correct_sink(dt)
    ! reduce P sink term to available P mass in conc
    Use Geometry, Only: PhoSink, PhoSinkRoot, NumNP, MatNum, ThNew, delta_z
    Use Solute, Only: ChPar, Conc
    Implicit None
    Real(dp), Intent(in) :: dt
    Real(dp) :: deficit,supply,demand
    ! correct sink term
    Do i=1,NumNP
       supply=Conc(Pho,i)*(ThNew(i)+ChPar(7,MatNum(i),Pho)*ChPar(1,MatNum(i),Pho))/dt
       If(PhoSink(i)<0.0d0) Then
          demand=-PhoSink(i)
          If(demand>supply) Then
             deficit=supply-demand
             supply=supply-deficit
             PhoSink(i)= Min(PhoSink(i)-deficit,0.0d0)
          Endif
       Endif
       If(lPRootUpt .And. PhoSinkRoot(i)>0.0d0) Then
          demand=PhoSinkRoot(i)
          If(demand>supply) demand=supply
          PhoSinkRoot(i)=-Max(demand,null_dp) ! negative means removal
       Endif
       ! accumulated over 1 timestep (day or hour), set to 0 in output
       PhoSinkBalance%sinksum=PhoSinkBalance%sinksum+Sum(PhoSink(:)*delta_z)
       PhoSinkBalance%sinkrootsum=PhoSinkBalance%sinkrootsum+Sum(PhoSinkRoot(:)*delta_z)
    Enddo
  End Subroutine p_correct_sink

  Subroutine PhosphorusTop(tnew, Prec, rSoil, cTop, cBot)
    ! inorganic (mineral) P fertilizer: goes as sink term into conc P
    Implicit None
    Integer, Intent(in) :: tnew
    Real(dp), Intent(in) :: Prec, rSoil
    Real(dp), Intent(out) :: cTop(:), cBot(:)
    Real(dp) :: nConc, vol
    
    If(debug) Print *,'PhosphorusTop'
    ! boundaries for phosphorus
    If(.Not.lPhosphorus) Return
    cTop(Pho)=0
    cBot(Pho)=0
    !Print *,'ntop',tnew,PTop_ind,PTop_dim,PTopTime
    If(PTopInd<=PTopDim) Then
       If(tnew==PTopTime(PTopInd)) Then
          P_on_top=P_on_top+PTopMass(PTopInd)
          PTopInd=PTopInd+1
       Endif
    Endif
    ! waiting for a flow into the soil
    vol=Prec-rSoil
    If(P_on_top>0 .And. vol>0.05/Units_L(Unit_L_Input)*Units_T(Unit_T_Input)) Then
       ! more than 0.05 mm/h
       nConc=P_on_top/vol
       cTop(Pho)=nconc
       Write(*,'(a,i8,a,1p,e13.5,a,e13.5,a,e13.5)')&
            'phosphorus on top, time=',tnew,' mass=',P_on_top,' flow=',vol,' conc=',nConc
       P_on_top=0
    Endif
  End Subroutine PhosphorusTop

  Subroutine PinpFromOC()
    ! organic P fertilizer: goes into Lit, Man and as sink term into conc P (labile inorganic)
    Use Geometry, Only: ThNew, OCinpLit, OCinpMan, rnodeResP, rnodedeadwP
    Implicit None
    Real(dp) :: fact

    If(debug) Print *,'PinpFromOC'
    If(.Not.lPhosphorus) Return
    fact=P_Fert_C*(1.0_dp-P_Fert_lab)
    PoolPhoLit=PoolPhoLit+OCinpLit*fact+(rnodeResP+rnodedeadwP)*0.59
    PoolPhoMan=PoolPhoMan+OCinpMan*fact+(rnodeResP+rnodedeadwP)*0.41
    ! set sink term for organic, add to mineral sink term
    fact=P_Fert_C*P_Fert_lab
    If(All(ThNew/=0)) Then
       PhoSinkFert=(OCinpLit+OCinpMan)*fact/ThNew
    Endif
    !Print *,'ocinp',fact,Sum(OCinpLit),Sum(PhoSinkFert)
  End Subroutine PinpFromOC

  Subroutine set_Balance(PhoSinkSum,r_active_labile_sum,r_active_stable_sum,sorp)
    Implicit None
    Real(dp), Intent(in) :: PhoSinkSum, r_active_labile_sum, r_active_stable_sum,sorp
    Balance%PhoSinkSum=PhoSinkSum
    Balance%r_active_labile_sum=r_active_labile_sum
    Balance%r_active_stable_sum=r_active_stable_sum
    Balance%sorp=sorp
  End Subroutine set_Balance

End Module Phosphorus
