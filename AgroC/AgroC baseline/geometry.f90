Module Geometry

  Use Datatypes
  Implicit None

  ! Constants
  Integer, Parameter :: dimCumQ=20
  Integer, Parameter :: Urea=1  ! urea index in Conc(solute,numnp)
  Integer, Parameter :: NH4=2   ! ammonia index in Conc(solute,numnp)
  Integer, Parameter :: NO3=3   ! nitrate index in Conc(solute,numnp)
  Integer, Parameter :: Pho=4   ! phosphorus index in Conc(solute,numnp)
  Integer, Parameter :: NSnit=3 ! number of solutes for nitrogen/phosphorus
  ! Variables
  Integer :: NSfrom=0           ! start index for solutes except nitrogen solutes
  Integer :: NumNP
  Integer :: NLay
  Integer :: NObs
  Integer :: NELz   ! number of elements
  Integer :: NSfirst ! first index number of solute other than nitrogen (1 or 4)

  ! real data.
  Real(dp) :: zSurf, coord_eps

  ! Arrays for nodes (dimension NumNP).
  Integer, Dimension(:), Allocatable :: MatNum,LayNum
  Real(dp), Dimension(:), Allocatable :: coord, delta_z
  Real(dp), Dimension(:), Allocatable :: hNew,hOld,hHelp,hTot
  Real(dp), Dimension(:), Allocatable :: Sink,Beta,Con,Cap, rnodert, rnodeh, rnodexu
  Real(dp), Dimension(:), Allocatable :: rnodedeadw, rnodedeadwN, rnodedeadwP
  Real(dp), Dimension(:), Allocatable :: rNodeResC, rNodeResN, rNodeResP
  Real(dp), Dimension(:), Allocatable :: HelpP,HelpR,HelpS,HelpQ,HelpT,HelpU,HelpV,hTemp
  Real(dp), Dimension(:), Allocatable :: ThNew,ThOld,vNew,vOld,Disp,ThOldIm,ThNewIm
  Real(dp), Dimension(:), Allocatable :: TempO,TempN,vAir,WatIn
  Real(dp), Allocatable :: NitSink(:,:), NitSinkRoot(:,:), PhoSink(:),PhoSinkRoot(:)
  Real(dp), Allocatable :: cRootMax(:) ! dim NS
  Real(dp), Allocatable :: OCinpMan(:), OCinpLit(:)
  ! arrays for layers
  Real(dp), Dimension(:), Allocatable :: SubVol,hMean, Ar,SubCha,SubCO,COMean,SubT,TMean
  ! other arrays.
  Integer, Dimension(:), Allocatable :: Node   ! dimension NObs
  Real(dp) :: CumQ(dimCumQ)

Contains

  Subroutine allocate_geometry_data()
    Allocate(coord(NumNP))
    Allocate(delta_z(NumNP))
    Allocate(hNew(NumNP))
    Allocate(hOld(NumNP))
    Allocate(hHelp(NumNP))
    Allocate(hTot(NumNP))
    Allocate(MatNum(NumNP))
    Allocate(Sink(NumNP))
    Allocate(rnodert(NumNP))
    Allocate(rNodeResC(NumNP))
    Allocate(rNodeResN(NumNP))
    Allocate(rNodeResP(NumNP))
    Allocate(rnodexu(NumNP))
    Allocate(rnodeh(NumNP))
    Allocate(rnodedeadw(numNP))
    Allocate(rnodedeadwN(numNP))
    Allocate(rnodedeadwP(numNP))
    Allocate(Beta(NumNP))
    Allocate(LayNum(NumNP))
    Allocate(Con(NumNP))
    Allocate(Cap(NumNP))
    Allocate(HelpP(NumNP))
    Allocate(HelpR(NumNP))
    Allocate(HelpS(NumNP))
    Allocate(HelpQ(NumNP))
    Allocate(HelpT(NumNP))
    Allocate(HelpU(NumNP))
    Allocate(HelpV(NumNP))
    Allocate(ThNew(NumNP))
    Allocate(ThOld(NumNP))
    Allocate(ThNewIm(NumNP))
    Allocate(ThOldIm(NumNP))
    Allocate(vNew(NumNP))
    Allocate(vOld(NumNP))
    Allocate(Disp(NumNP))
    Allocate(TempO(NumNP))
    Allocate(TempN(NumNP))
    Allocate(vAir(NumNP))
    Allocate(WatIn(NumNP))
    Allocate(hTemp(NumNP))
    Allocate(NitSink(NumNP,NSnit))     ! 
    Allocate(NitSinkRoot(NumNP,NSnit)) ! plant uptake
    Allocate(OCinpMan(NumNP))
    allocate(OCinpLit(Numnp))
    Allocate(SubVol(NLay))
    Allocate(hMean(NLay))
    Allocate(Ar(NLay))
    Allocate(SubCha(NLay))
    Allocate(SubCO(NLay))
    Allocate(COMean(NLay))
    Allocate(SubT(NLay))
    Allocate(TMean(NLay))
    ! observation nodes
    Allocate(Node(NObs))
    CumQ=0.0
    Sink=0.0
    NELz=NumNP-1
    rnodert = 0.0_dp
    rnodexu = 0.0_dp
    rnodeh = 0._dp
    rnodedeadw = 0_dp
    rnodedeadwN = 0_dp
    rnodedeadwN = 0_dp
    rnodedeadwP = 0_dp
    rNodeResC = 0.0_dp
    rNodeResN = 0.0_dp
    rNodeResP = 0.0_dp

  End Subroutine allocate_geometry_data

End Module Geometry
