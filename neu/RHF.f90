MODULE RHF
  implicit NONE
  private
  public dic,terminate,P_init,G_matrix,Energy,RHF_control,F_matrix
  public twomatrix,P,G,F,H
  public NORB,NELEC,E,E_nucc,E_old,iters
      Real, DImension(:,:,:,:),allocatable :: twomatrix
      Real, Dimension(:,:),allocatable :: H,P,G,F
      INTEGER:: NORB,NELEC,iters
      REAL:: E,E_old,E_nucc

  CONTAINS
    SUBROUTINE dic()
      IMPLICIT NONE
      Integer, PARAMETER:: MAX_size=100
      INTEGER, PARAMETER:: LU=99
      REAL:: integral,E_nucc
      INTEGER:: my,ny,la,si,n,index1,index2,ierror
      INTEGER:: MS2,ISYM
      REAL,DIMENSION(MAX_size)::ORBSYM
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
        Write (*,1020) filename
        1020 Format (1X, 'Error')
      END if errorcheck
      READ(UNIT=LU,NML=FCI)
      allocate(H(norb,norb),P(NORB,NORB),G(NORB,NORB),F(NORB,NORB))
      allocate(twomatrix(norb,norb,norb,norb))
      
      DO
       READ(LU,300,IOSTAT=ierror,IOMSG=ierrormsg) integral, my, ny, la, si
      300 Format (1x,E23.16,4I4)
        If (ierror /= 0) then
          write(*,*) ierrormsg
          EXit
        end if
        if (my==0 .and. ny==0 .and. la==0 .and. si==0) then
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
        end if
      END DO

    END Subroutine dic
    subroutine terminate()
deallocate(twomatrix,H,P,G,F)
      end subroutine terminate

    Subroutine P_init()
      IMPLICIT NONE
      integer::i=1
      DO i = 1,int(NELEC/2)
        P(i,i)=2
      END DO
    END Subroutine P_init

  Subroutine P_matrix
    Implicit NONE
    integer::my,ny,a

    DO my=1,NORB
      DO ny=1,NORB
        P(my,ny)=0
        DO a=1,int(NELEC/2)
         P(my,ny)=P(my,ny)+2*C(my,a)*C(ny,a)
        End Do
      End Do
    End Do
  End Subroutine P_matrix

  Subroutine G_matrix
    IMPLICIT NONE
    integer::my,ny,la,si

    Do my=1,NORB
      DO ny=1,NORB
        G(my,ny)=0
        DO la=1,NORB
          Do si=1,NORB
            G(my,ny)=G(my,ny)+P(la,si)*(twomatrix(my,ny,la,si)-0.5*twomatrix(my,la,si,ny))
          END Do
        End Do
      End Do
    END DO
  END Subroutine G_matrix

  Subroutine F_matrix
    Implicit NONE

    F=H+G
  End Subroutine F_matrix

  Subroutine Energy
    IMPLICIT NONE
    INTEGER:: my,ny,la,si

    E=E_nucc
    Do my=1,NORB
      DO ny=1,NORB
        E=E+P(my,ny)*(H(my,ny)+0.5*G(my,ny))
      END Do
    END Do
  END Subroutine ENergy

Subroutine RHF_control
IMPLICIT NONE
integer::iters=1
real::E_old=0

call dic()

DO
  if (iters==1) then
    call P_init
  else
    call P_matrix
  end if
  call G_matrix
  call F_matrix
  call diagonal(F)
  call Energy
  write(*,*)'ENERGY:',E,'Iteration:',iters
  if (E_old==E) then
    exit
  end if
  E_old=E
  iters=iters+1
END DO 
Call terminate()
end Subroutine RHF_control

END MODULE RHF

PROGRAM main
  use RHF
  call RHF_control
END Program
