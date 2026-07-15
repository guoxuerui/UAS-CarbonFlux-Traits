Module Datatypes

  Use ISO_FORTRAN_ENV
  Integer, Parameter :: dp = real64 ! Kind(1.0D0)
  Real(dp), Parameter :: null_dp=0_dp
  Real(dp), Parameter :: one_dp=1_dp
  Real(dp), Parameter :: PI=4.0_dp*Atan(one_dp) ! pi/4 = atan(1)
  Real(dp), Parameter :: eps=1e-12
  Logical :: debug, warnings=.False.
  Character(80) :: headline
  Character(4) ::  LUnit,TUnit,MUnit
  Integer :: Unit_L_Input      ! 1=mm 2=cm 3=dm 4=m 5=km
  Real(dp), Dimension(5) ::Units_L=(/ 1.0_dp,10.0_dp,100.0_dp,1000.0_dp,1.0e6_dp /)
  Real(dp), Dimension(5) ::to_ha=(/ 1.0e10_dp,1.0e8_dp,1.0e6_dp,1.0e4_dp,1.0e-2_dp /)
  Integer :: Unit_T_Input      ! 1=h 2=d
  Real(dp), Dimension(2) ::Units_T=(/ 1.0_dp,24.0_dp /)
  Integer :: Unit_M_Input      ! 1=mg 2=g 3=kg 4=t
  Real(dp), Dimension(4) ::Units_M=(/ 1.0_dp,1e3_dp,1e6_dp,1e9_dp /)
  Type TabType
     Integer :: nrows, ncols
     Real(dp), Pointer :: data(:,:)
  End Type TabType

Contains
  ! CCC  FUNKTION ZUR BERECHNUNG DER ANZAHL DER WOERTER AUF DEM STRING S
  Integer Function wordsf(s)
    Implicit None
    character(len=*) :: s
    integer i,j,verify
    wordsf=0
    i=1
    do
       j=verify(s(i:),' ')
       if(j.eq.0) exit
       wordsf=wordsf+1
       i=i+j-1
       j=index(s(i:),' ')
       if(j.eq.0) exit
       i=i+j-1
    end do
  End Function wordsf
  
! read words from line
  Integer Function read_words(s,fields)
    Implicit None
    Character(len=*), Intent(in) :: s
    Character(len=*), Intent(out) :: fields(:)
    Integer :: i,j,verify
    read_words=0
    i=1
    do
       j=verify(s(i:),' ')
       if(j.eq.0) exit
       read_words=read_words+1
       If(read_words>Ubound(fields,1)) Exit
       i=i+j-1
       j=index(s(i:),' ')
       If(j.Eq.0) Exit
       fields(read_words)=s(i:i+j-1)
       i=i+j-1
    end do
  End Function read_words

  Subroutine WriteError(text)
    Implicit None
    Character(*), Intent(in) :: text
    Write(*,*) Trim(text)
    Stop 1
  End Subroutine WriteError

  Subroutine WriteWarning(text)
    Implicit None
    Character(*), Intent(in) :: text
    Logical, Save :: first_call=.True.
    If(first_call) Then
       first_call=.False.
       warnings=.True.
!       Open(52,file='agroc.warnings',Status='REPLACE')
    Endif
!    Write(52,*) Trim(text)
    Print *, Trim(text)
  End Subroutine WriteWarning

  Subroutine WriteOutput(text)
    Implicit None
    Character(*), Intent(in) :: text
    Write(*,*) Trim(text)
  End Subroutine WriteOutput


  Subroutine AllocateError(object)
     Character*(*) :: object
     Write(*,*) 'Error while allocating object: ',Trim(object)
     Stop
  End Subroutine AllocateError

  ! used in input 
  Subroutine SkipBlock(Block)
    Character(*), Intent(in) :: block
    Do While(headline(1:7) /= '*** END')
       Read(30,'(A)',End=1) headline
       If(headline(1:9) == '*** BLOCK') Exit
    End Do
    Return
1   Call WriteError('Unexpected end of file (selector.in) while searching for Block ' // block)
  End Subroutine SkipBlock
  
  Function to_lower(str) Result (string)
    !   ==============================
    !   Changes a string to lower case
    !   ==============================
    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string
    Integer :: ic, i
    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    !   lowercase each letter
    string = str
    do i = 1, LEN_TRIM(str)
       ic = INDEX(cap, str(i:i))
       if (ic > 0) string(i:i) = low(ic:ic)
    end do

  End Function to_lower
  
  Function equal(x,y) Result (eq)
    ! compare floating point variables for equality
    Implicit None
    Logical :: eq
    Real(dp), Intent(in) :: x, y
    
    If(Abs(x-y)<eps) Then
       eq=.True.
    Else
       eq=.False.
    Endif
  End Function equal

End Module Datatypes
