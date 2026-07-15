Module Material

  Use datatypes
  Implicit None
  Integer, Parameter :: MaxPar=28
  Integer :: NTab=100 ! number of rows in table
  Integer :: NMat   ! number of materials.
  Integer :: NPar   ! number of parameters for each material.
  Integer :: iModel=1
  ! 1 = Mualem / van Genuchten
  ! 2 = dual-porosity function (Durner)
  
!     iModel = 0: van Genuchten
!              1: modified van Genuchten (Vogel and Cislerova)
!              2: Brooks and Corey
!              3: van Genuchte with air entry value of 2 cm
!              4: log-normal (Kosugi)
!              5: dual-porosity function (Durner)
!              6: dual-porosity system with transfer proportional to water content
!              7: dual-porosity system with transfer proportional to pressure head
!              8: dual-permeability system, not handled by this program
!              9: nTabMod: general tables (nTabMod)
  ! real data
  Real(dp) :: hTab1,hTabN

  ! Arrays for materials (dimension NMat).
  Real(dp), Allocatable :: thR(:),hSat(:),ConSat(:),thS(:)
  Real(dp), Allocatable :: Par(:,:)   ! dimension MaxPar x NMat
  ! Arrays for materials (dimension NTab x NMat).
  Real(dp), Allocatable :: ConTab(:,:),CapTab(:,:),TheTab(:,:)
  ! Arrays for materials (dimension NTab).
  Real(dp), Allocatable :: hTab(:)
  ! parameters for temperature
  Real(dp), Allocatable :: tempParam(:,:)

Contains

  Subroutine allocate_material_data()
    allocate(thR(NMat))
    allocate(hSat(NMat))
    allocate(ConSat(NMat))
    allocate(thS(NMat))
    allocate(Par(MaxPar,NMat))
    Allocate(tempParam(9,NMat))
  end Subroutine allocate_material_data

  Subroutine allocate_material_tab_data()
    allocate(ConTab(NTab,NMat))
    allocate(CapTab(NTab,NMat))
    allocate(TheTab(NTab,NMat))
    allocate(hTab(NTab))
  end Subroutine allocate_material_tab_data


  
  Function FK(h,Par)
    Use datatypes
    Implicit None
    Real(dp) :: FK,h,Par(:)
    ! local variables
    Real(dp) :: BPar,Qr,Qs,Qa,Qm,Qk,Alfa,n,m,Ks,Kr,Kk,HMin,HH,Qees,Qeek,Hs,Hk,Qee,Qe,Qek,FFQ,FFQk
    Real(dp) :: w1,w2,n2,Alfa2,m2,Sw1,Sw2,Sv1,Sv2,Sk1,Sk2,rNumer,rDenom
    Integer :: PPar

    Qr=Par(1)
    Qs=Par(2)
    Alfa=Par(5)
    n=Par(6)
    Ks=Max(Par(7),1.e-37)
    BPar=Par(10)
    FK=0
    If(iModel==1) Then
       PPar=2
       Qa=Par(3)
       Qm=Par(4)
       Kk=Par(8)
       Qk=Par(9)
       m=1.-1./n
       HMin=-1.d300**(1./n)/max(Alfa,1.d0)
       HH=max(dble(h),HMin)
       Qees=min((Qs-Qa)/(Qm-Qa),.999999999999999d0)
       Qeek=Min((Qk-Qa)/(Qm-Qa),Qees)
       !    Print *,'Qees=',Qees,' Qeek=',Qeek,' m=',m,' n=',n
       Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
       Hk=-1./Alfa*(Qeek**(-1./m)-1.)**(1./n)
       !    Print *,'Hs=',Hs,' Hk=',Hk
       If(h.Lt.Hk) Then
          Qee=(1.+(-Alfa*HH)**n)**(-m)
          Qe =(Qm-Qa)/(Qs-Qa)*Qee
          Qek=(Qm-Qa)/(Qs-Qa)*Qeek
          FFQ =1.-(1.-Qee**(1./m))**m
          FFQk=1.-(1.-Qeek**(1./m))**m
          if(FFQ.le.0.d0) FFQ=m*Qee**(1./m)
          Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
          FK=max(Ks*Kr,1.d-37)
       else if(h.lt.Hs) then
          Kr=(1.-Kk/Ks)/(Hs-Hk)*(h-Hs)+1.
          FK=Ks*Kr
       else 
          FK=Ks
       Endif
    Elseif(iModel==2) Then
       w2=Par(11)
       Alfa2=Par(12)
       n2=Par(13)
       m =1.d0-1.d0/n
       m2=1.d0-1.d0/n2
       w1=1.d0-w2
       Sw1=w1*(1.d0+(-Alfa *h)**n )**(-m )
       Sw2=w2*(1.d0+(-Alfa2*h)**n2)**(-m2)
       Qe=Sw1+Sw2
       Sv1=(-Alfa *h)**(n -1)
       Sv2=(-Alfa2*h)**(n2-1)
       Sk1=w1*Alfa *(1.d0-Sv1*(1.d0+(-Alfa *h)**n )**(-m ))
       Sk2=w2*Alfa2*(1.d0-Sv2*(1.d0+(-Alfa2*h)**n2)**(-m2))
       rNumer=Sk1+Sk2
       rDenom=w1*Alfa+w2*Alfa2
       If(rDenom.Ne.0.) Then
          Kr=Qe**BPar*(rNumer/rDenom)**2
       Else
          Kr=1
       Endif
       FK=max(Ks*Kr,1.0d-37)
    Endif
  End Function FK

  !***********************************************************************


  function FC(h,Par)
    Use datatypes
    Implicit None
    Real(dp) :: FC,h,Par(:)
    ! local variables
    Real(dp) :: Qr,Qs,Qa,Qm,Alfa,n,m,HMin,HH,Hs,C1,C2,Qees
    Real(dp) :: w1,w2,Alfa2,n2,m2,C1a,C2a,C1b,C2b

    Qr=Par(1)
    Qs=Par(2)
    Alfa=Par(5)
    n=Par(6)
    FC=0
    If(iModel==1) Then
       Qa=Par(3)
       Qm=Par(4)
       m=1.-1./n
       HMin=-1.d300**(1./n)/max(Alfa,1.d0)
       HH=max(dble(h),HMin)
       Qees=min((Qs-Qa)/(Qm-Qa),.999999999999999d0)
       Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
       if(h.lt.Hs) then
          C1=(1.+(-Alfa*HH)**n)**(-m-1.)
          C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.)*C1
          FC=max(C2,1.d-37)
          return
       else
          FC=0.0
       end if
    Else If(iModel==2) Then
       w2=Par(11)
       Alfa2=Par(12)
       n2=Par(13)
       m =1.d0-1.d0/n
       m2=1.d0-1.d0/n2
       w1=1.d0-w2
       C1a=(1.d0+(-Alfa *h)**n )**(-m -1.d0)
       C1b=(1.d0+(-Alfa2*h)**n2)**(-m2-1.d0)
       C2a=(Qs-Qr)*m *n *(Alfa **n )*(-h)**(n -1.d0)*C1a*w1
       C2b=(Qs-Qr)*m2*n2*(Alfa2**n2)*(-h)**(n2-1.d0)*C1b*w2
       FC=C2a+C2b
    Endif
  end function FC

  !***********************************************************************

  function FQ(h,Par)

    Use datatypes
    Implicit None
    Real(dp) :: FQ,h,Par(:)
    ! local variables
    Real(dp) :: Qr,Qs,Qa,Qm,Alfa,n,m,HMin,HH,Qees,Hs,Qee
    Real(dp) :: w1,w2,Alfa2,n2,m2,Sw1,Sw2,Qe

    Qr=Par(1)
    Qs=Par(2)
    Alfa=Par(5)
    n=Par(6)
    FQ=0
    If(iModel==1) Then
       Qa=Par(3)
       Qm=Par(4)
       m=1.-1./n
       HMin=-1.d300**(1./n)/max(Alfa,1.d0)
       HH=max(dble(h),HMin)
       Qees=min((Qs-Qa)/(Qm-Qa),.999999999999999d0)
       Hs=-1./Alfa*(Qees**(-1./m)-1.)**(1./n)
       if(h.lt.Hs) then
          Qee=(1.+(-Alfa*HH)**n)**(-m)
          FQ=max(Qa+(Qm-Qa)*Qee,1.d-37)
          return
       else
          FQ=Qs
       end if
    Else If(iModel==2) Then
       w2=Par(11)
       Alfa2=Par(12)
       n2=Par(13)
       m =1.d0-1.d0/n
       m2=1.d0-1.d0/n2
       w1=1.d0-w2
       Sw1=w1*(1.d0+(-Alfa *h)**n )**(-m )
       Sw2=w2*(1.d0+(-Alfa2*h)**n2)**(-m2)
       Qe=Sw1+Sw2
       FQ=Max(Qr+(Qs-Qr)*Qe,1.d-37)
    End If
  end function FQ

  !***********************************************************************

  function FH(Qe,Par)
    Use datatypes
    Implicit None
    Real(dp) :: FH,Qe,Par(:)
    ! local variables
    Real(dp) :: Qr,Qs,Qa,Qm,Alfa,n,m,HMin,QeeM,Qee
    Real(dp) :: w1,w2,Alfa2,n2,m2,h

    Qr=Par(1)
    Qs=Par(2)
    Alfa=Par(5)
    n=Par(6)
    FH=0
    If(iModel==1) Then
       Qa=Par(3)
       Qm=Par(4)
       m=1.-1./n
       HMin=-1.d300**(1./n)/max(Alfa,1.d0)
       QeeM=(1.+(-Alfa*HMin)**n)**(-m)
       Qee=min(max(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
       FH=max(-1./Alfa*(Qee**(-1./m)-1.)**(1./n),-1.d37)
    Else If(iModel==2) Then
       w2=Par(11)
       Alfa2=Par(12)
       n2=Par(13)
       m =1.d0-1.d0/n
       m2=1.d0-1.d0/n2
       w1=1.d0-w2
       Qee=Qe
       if(Qee.gt.0.9999d0) then
          FH=0.0
       else if(Qee.lt.0.00001d0) then
          FH=-1.e+8
       else
          h=xMualem(Qee,Par)
          FH=Max(h,-1.d37)
       end if
    End If
  end function FH


  !     Evaluate h for given theta_e for dual-porosity function
  Function xMualem(Se,Par)
    Use Variables, Only: Pause
    Implicit None
    Real(dp) :: Se,Par(:),xMualem
    Real(dp) :: x1,x2,xb1,xb2,hhh

    x1=-1.e-6
    x2=-1.e+6
    call ZBRAK(X1,X2,XB1,XB2,SE,Par)
    hhh=ZBRENT(XB1,XB2,SE,Par)
    xMualem=hhh    ! for calculation hh  
    if(hhh.ne.0.) then
       !        xMualem=1./hhh ! for integration
    else
       Call Pause('xMualem: h is equal to zero!')
    end if
  End Function xMualem

  
  Function DoublePor(hh,SE,Par)
    ! Double porosity function - for evaluation of h for given theta_e
    Implicit None
    Real(dp) :: hh,Se,Par(:),DoublePor
    Real(dp) :: wcr,wcs,Alpha,Alpha2,rn,rn2,rm,rm2,w1,w2,Sw1,Sw2,rwc

    wcr=Par(1)
    wcs=Par(2)
    Alpha=Par(5)
    rn=Par(6)
    rm=1.-1./rn
    w2=Par(11)
    w1=1.-W2
    Alpha2=Par(12)
    rn2=Par(13)
    rm2=1.-1./rn2

    Sw1=w1*(1.+(-Alpha *hh)**rn )**(-rm )
    Sw2=w2*(1.+(-Alpha2*hh)**rn2)**(-rm2)
    rwc=Sw1+Sw2
    DoublePor=SE-rwc

  End Function DoublePor


  !     Bracketing of the root, Numerical recepies (345)
  Subroutine ZBRAK(X1,X2,XB1,XB2,SE,Par)
    Implicit None
    Real(dp) :: x1,x2,xb1,xb2,se,Par(:)
    Integer :: nb,nbb,i,n
    Real(dp) :: dlh,fp,dx2,fc

    NB=1
    NBB=NB
    NB=0
    n=1000

    dlh=(Log10(-X2)-Log10(-X1))/(N-1)
    FP=DoublePor(X1,SE,Par)
    Do i=1,n
       dx2=Log10(-X1)+i*dlh
       X2=-10**dx2
       FC=DoublePor(X2,SE,Par)
       if(FC*FP.lt.0.) then
          XB1=X1
          XB2=X2
          return
       end if
       FP=FC
       X1=X2
       if(NBB==NB) return
    Enddo
  End Subroutine ZBRAK
  
  !     Brent method of finding root that lies between x1 and x2, 
  !     Numerical recepies (354)
  Function ZBRENT(X1,X2,SE,Par)
    Use Variables, Only: Pause
    implicit None
    Integer, Parameter :: ITMAX=100
    Real(dp), Parameter :: EPS=3.0E-8,TOL=1.0e-6
    Real(dp) :: x1,x2,se,Par(:),ZBRENT
    Real(dp) :: a,b,fa,fb,fc,c,d=0,e=0,xm,s,p,q,r,TOL1
    Integer :: iter
    
    A=X1
    B=X2
    c=0
    FA=DoublePor(A,SE,Par)
    FB=DoublePor(B,SE,Par)
    If(FB*FA.Gt.0.) Call Pause('Root must be bracketed for ZBRENT.')
    FC=FB
    DO ITER=1,ITMAX
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
       TOL1=2.*EPS*Abs(B)+0.5*TOL
       XM=.5*(C-B)
       IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          ZBRENT=B
          RETURN
       ENDIF
       IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
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
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
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
       FB=DoublePor(B,SE,Par)
    end do
    call PAUSE('ZBRENT exceeding maximum iterations.')
    ZBRENT=B
  END function ZBRENT

end Module Material
