MODULE UHF
  implicit NONE
  private
  public init,P_init,G_matrixa,G_matrixb,G_matrixt,Energy,UHF_control,F_matrix,P_matrixa,P_matrixb,P_matrixt
  public twomatrix,Pa,Pb,Pt,Ga,Gb,Gt,Fa,Fb,H,vecsa,vecsb
  public NORB,NELECa,NELECb,E,E_nucc,E_old,iters,valsa,valsb
  DOUBLE PRECISION, DImension(:,:,:,:),allocatable :: twomatrix
  DOUBLE PRECISION, Dimension(:,:),allocatable :: H,Pa,Pb,Pt,Ga,Gb,Gt,Fa,Fb,vecsa,vecsb
  DOUBLE PRECISION, Dimension(:),allocatable:: valsa,valsb
      INTEGER:: NORB,NELEC,iters,NELECa,NELECb
      DOUBLE PRECISION:: E,E_old,E_nucc

  CONTAINS
    SUBROUTINE init
      IMPLICIT NONE
      Integer, PARAMETER:: MAX_size=100
      INTEGER, PARAMETER:: LU=99
      DOUBLE PRECISION:: integral
      INTEGER:: my,ny,la,si,n,index1,index2,ierror,elecadd
      INTEGER:: MS2,ISYM
      DOUBLE PRECISION,DIMENSION(MAX_size)::ORBSYM
      Character(len=24):: filename
      Character(len=128):: ierrormsg
      NAMELIST/FCI/NORB,NELEC,MS2,ORBSYM,ISYM
      
      !Write (*,1000)
      !1000 Format (1x,'Enter Filename')
      !Read (*,1010) filename
      !1010 Format (A)
      filename='FCIDUMP'
      open (Unit=LU, File=filename,Status='old', iostat=ierror)
      errorcheck: IF (ierror>0) THen
        write(*,*) ierror
        write(*,*) 'fileerror'
      END if errorcheck

      READ(UNIT=LU,NML=FCI)
      Write (*,1001)
      1001 Format (1x,'HOW MUCH ELECTRONS TO ADD?')
      Read (*,1011) elecadd 
      1011 Format (I3)
      NELEC=NELEC+elecadd
      IF (MODULO(NELEC,2)==0) then
        NELECa=int(NELEC/2)
        NELECb=int(NELEC/2)
      ELSE
        NELECa=int((NELEC+1)/2)
        NELECb=int((NELEC-1)/2)
      END IF
      write (*,*) 'NELECa:',NELECa,'NELECb:',NELECb
      allocate(H(norb,norb),Pa(NORB,NORB),Pb(NORB,NORB),Pt(NORB,NORB),Ga(NORB,NORB),Gb(NORB,NORB),Gt(NORB,NORB),&
          Fa(NORB,NORB),Fb(NORB,NORB),vecsa(NORB,NORB),vecsb(NORB,NORB))
      allocate(valsa(NORB),valsb(NORB))
      allocate(twomatrix(norb,norb,norb,norb))
      
      DO
       READ(LU,300,IOSTAT=ierror,IOMSG=ierrormsg) integral, my, ny, la, si
      300 Format (1x,E23.16,4I4)
        error: If (ierror /= 0) then
          !write(*,*) ierrormsg
          EXit
        end if error
       integral_type: if (my==0 .and. ny==0 .and. la==0 .and. si==0) then
          E_nucc=integral
        else if (la==0 .and. si==0) Then
          H(my,ny)=integral
          H(ny,my)=integral
        else
          twomatrix(my,ny,la,si)=integral
          twomatrix(my,ny,si,la)=integral
          twomatrix(ny,my,si,la)=integral
          twomatrix(ny,my,la,si)=integral
          twomatrix(la,si,my,ny)=integral
          twomatrix(la,si,ny,my)=integral
          twomatrix(si,la,my,ny)=integral
          twomatrix(si,la,ny,my)=integral
        end if integral_type
      END DO
    END Subroutine init

    !subroutine terminate()
!deallocate(twomatrix,H,P,G,F)
    !  end subroutine terminate

    Subroutine P_init()
      IMPLICIT NONE
      integer::i=1
        Do i=1,NELECa
          Pa(i,i)=1
        END DO
        Do i=1,NELECb
          Pb(i,i)=1
        END DO
      Pt=Pa+Pb
    END Subroutine P_init

  Subroutine P_matrixa
    Implicit NONE
    integer::my,ny,a
    DO my=1,NORB
      DO ny=1,NORB
        Pa(my,ny)=0
        DO a=1,NELECa
         Pa(my,ny)=Pa(my,ny)+vecsa(my,a)*vecsa(ny,a)
        End Do
      End Do
    End Do
  End Subroutine P_matrixa

  Subroutine P_matrixb
    Implicit NONE
    integer::my,ny,a
    DO my=1,NORB
      DO ny=1,NORB
        Pb(my,ny)=0
        DO a=1,NELECb
         Pb(my,ny)=Pb(my,ny)+vecsb(my,a)*vecsb(ny,a)
        End Do
      End Do
    End Do
  End Subroutine P_matrixb

  Subroutine P_matrixt
    Pt=Pa+Pb
  End Subroutine P_matrixt

  Subroutine G_matrixa
    IMPLICIT NONE
    integer::my,ny,la,si

    Do my=1,NORB
      DO ny=1,NORB
        Ga(my,ny)=0
        DO la=1,NORB
          Do si=1,NORB
            Ga(my,ny)=Ga(my,ny)+Pa(la,si)*(-1)*twomatrix(my,la,si,ny)
          END Do
        End Do
      End Do
    END DO
  END Subroutine G_matrixa

  Subroutine G_matrixb
    IMPLICIT NONE
    integer::my,ny,la,si

    Do my=1,NORB
      DO ny=1,NORB
        Gb(my,ny)=0
        DO la=1,NORB
          Do si=1,NORB
            Gb(my,ny)=Gb(my,ny)+Pb(la,si)*(-1)*twomatrix(my,la,si,ny)
          END Do
        End Do
      End Do
    END DO
  END Subroutine G_matrixb

  Subroutine G_matrixt
    IMPLICIT NONE
    integer::my,ny,la,si

    Do my=1,NORB
      DO ny=1,NORB
        Gt(my,ny)=0
        DO la=1,NORB
          Do si=1,NORB
            Gt(my,ny)=Gt(my,ny)+Pt(la,si)*twomatrix(my,ny,la,si)
          END Do
        End Do
      End Do
    END DO
  END Subroutine G_matrixt

  Subroutine F_matrix
    Implicit NONE
    Fa=H+Ga+Gt
    Fb=H+Gb+Gt
    !write(*,*) F
  End Subroutine F_matrix

  Subroutine Energy
    IMPLICIT NONE
    INTEGER:: my,ny,la,si

    E=E_nucc
    Do my=1,NORB
      DO ny=1,NORB
        E=E+0.5*(Pt(ny,my)*H(my,ny)+Pa(ny,my)*Fa(my,ny)+Pb(ny,my)*Fb(my,ny))
      END Do
    END Do
  END Subroutine ENergy

  Subroutine matprint(title,A)
    DOUBLE PRECISION,DIMENSION(NORB,NORB)::A
    INTEGER::i
    CHARacter::title
    write(*,*) title
    Do i=1,NORB
      write (*,1500) A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),A(i,8),A(i,9),A(i,10),A(i,11),A(i,12),A(i,13)
      1500 Format (' ',13ES16.8)
    End Do
  End Subroutine matprint

  Subroutine fcimaker
    IMPLICIT NONE
    INTEGER::I,J,K,L,my,ny,la,si
    DOUBLE PRECISION:: integral
    INTEGER,DIMENSION(NORB)::ORBSYM
    INTEGER::ierror,fcidumpnew,MS2=0,ISYM=0
    CHARACTER(len=10)::filename
    CHARACTER(len=124)::imsg
    NAMELIST/FCI/NORB,NELEC,MS2,ORBSYM,ISYM
    CALL get_command_argument(2,filename)
    ORBSYM=1
    fcidumpnew=25
    OPEN(UNIT=fcidumpnew, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
    WRITE(UNIT=fcidumpnew,NML=FCI)
    !alpha-alpha
    DO I=1,NORB
      Do J=1,I
        Do K=1,I
         Do L=1,K
          if (K==I .and. L>J) then
            CYCLE
          end if
          integral=0
          DO my=1,NORB
            DO ny=1,NORB
              Do la=1,NORB
                Do si=1,NORB
                  integral=integral+vecsa(my,I)*vecsa(ny,J)*vecsa(la,K)*vecsa(si,L)*twomatrix(my,ny,la,si)
                END DO
              END DO
            END DO
          END DO
          WRITE(fcidumpnew,300,IOSTAT=ierror) integral, I, J, K, L
          300 Format (1x,E23.16,4I4)
        END DO
      END DO
    END DO
  END DO
  !beta-beta
  DO I=1,NORB
    DO J=1,I
      DO K=1,I
        DO L=1,K
          if (K==I .and. L>J) then
            CYCLE
          end if
          integral=0
          DO my=1,NORB
            DO ny=1,NORB
              DO la=1,NORB
                DO si=1,NORB
                  integral=integral+vecsb(my,I)*vecsb(ny,J)*vecsb(la,K)*vecsb(si,L)*twomatrix(my,ny,la,si)
                END DO
              END DO
            END DO
          END DO
          WRITE(fcidumpnew,300,IOSTAT=ierror) integral, I+NELECa, J+NELECa, K+NELECa, L+NELECa
        END DO
      END DO
    END DO
  END DO
  !alpha-beta
  DO I=1,NORB
    DO J=1,I
      DO K=1,I
        DO L=1,K
          if (K==I .and. L>J) then
            CYCLE
          end if
          integral=0
          DO my=1,NORB
            DO ny=1,NORB
              DO la=1,NORB
                DO si=1,NORB
                  integral=integral+vecsa(my,I)*vecsa(ny,J)*vecsb(la,K)*vecsb(si,L)*twomatrix(my,ny,la,si)
                END DO
              END DO
            END DO
          END DO
          WRITE(fcidumpnew,300,IOSTAT=ierror) integral, I+NELECa+NELECb, J+NELECa+NELECb, K+NELECa+NELECb, L+NELECa+NELECb
        END DO
      END DO
    END DO
  END DO
  K=0
  L=0
  !alpa-alpha
  DO I=1,NORB
    DO J=1,I
      integral=0
      DO my=1,NORB
        DO ny=1,NORB
          integral=integral+vecsa(my,I)*vecsa(ny,J)*H(my,ny)
        END DO
      END DO
    WRITE(fcidumpnew,300,IOSTAT=ierror) integral, I, J, K, L
    END DO
  END DO
  !beta-beta
  DO I=1,NORB
    DO J=1,I
      integral=0
      DO my=1,NORB
        DO ny=1,NORB
          integral=integral+vecsb(my,I)*vecsb(ny,J)*H(my,ny)
        END DO
      END DO
    WRITE(fcidumpnew,300,IOSTAT=ierror) integral, I+NELECa, J+NELECa, K+NELECa, L+NELECa
    END DO
  END DO
  I=0
  J=0
  integral=E_nucc
  WRITE(fcidumpnew,300,IOSTAT=ierror,IOMSG=imsg) integral, I, J, K, L
  !300 Format (1x,E23.16,4I4)
END SUbroutine fcimaker

Subroutine terminate
  Deallocate(Ga,Gb,Gt,Pa,Pb,Pt,H,twomatrix,Fa,Fb,vecsa,vecsb)
END SUBroutine terminate

Subroutine UHF_control

IMPLICIT NONE
integer::iters=1,i
DOUBLE PRECISION::E_old=0
INTEGER::INFO
Double Precision, Parameter:: tresh=1d-10
DOUBLE PRECISION, Dimension(:),allocatable::work
integer :: lwork,lwork_max

call init()
lwork_max=norb**2
allocate(work(lwork_max))
lwork=NORB**2


DO
  if (iters==1) then
    call P_init
  else
    call P_matrixa
    call P_matrixb
    call P_matrixt
  end if
  call G_matrixa
  call G_matrixb
  call G_matrixt
  call F_matrix
  call Energy
  call DSYEV('V','L',NORB,Fa,NORB,valsa,work,lwork,INFO)
  vecsa=Fa
  call DSYEV('V','L',NORB,Fb,NORB,valsb,work,lwork,INFO)
  vecsb=Fb
  write(*,1234)E,iters
  1234 Format (' ','ENERGY:',F29.10,'Iteration:',I3)
  if (abs(E_old-E)<tresh) then
    write(*,*) 'converged'
    exit
  end if
  E_old=E
  iters=iters+1
  if (iters>50) then
    exit
  end if
END DO 
write(*,*) 'DONE'
Call fcimaker
deallocate(work)
Call terminate()
end Subroutine UHF_control

END MODULE UHF

PROGRAM main
  use UHF
  call UHF_control
END Program
