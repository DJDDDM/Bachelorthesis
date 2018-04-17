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
      REAL:: E,E_old,E_nucc

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
        DO a=1,int(NELEC/2)
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
        DO a=1,int(NELEC/2)
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

Subroutine UHF_control

IMPLICIT NONE
integer::iters=1,i
DOUBLE PRECISION::E_old=0
INTEGER::INFO
DOUBLE PRECISION, Dimension(NORB,NORB)::trash,Fa_store,Fb_store

call init()


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
  !Fa_store=Fa
  !Fb_store=Fb
  call Energy
  call matprint('Pt',Pt)
  call DSYEV('V','L',NORB,Fa,NORB,valsa,trash,NORB*NORB,INFO)
  vecsa=Fa
  call DSYEV('V','L',NORB,Fb,NORB,valsb,trash,NORB*NORB,INFO)
  vecsb=Fb
  !Fa=Fa_store
  !Fb=Fb_store
  write(*,1234)E,iters
  1234 Format (' ','ENERGY:',ES29.20,'Iteration:',I3)
  if (abs(E_old-E)<0.00000001) then
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
!Call terminate()
end Subroutine UHF_control

END MODULE UHF

PROGRAM main
  use UHF
  call UHF_control
END Program
