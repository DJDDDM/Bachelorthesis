MODULE RHF
  implicit NONE
  private
  public dic,terminate,P_init,G_matrix,Energy,RHF_control,F_matrix
  public twomatrix,P,G,F,H,vecs
  public NORB,NELEC,E,E_nucc,E_old,iters,vals
  DOUBLE PRECISION, DImension(:,:,:,:),allocatable :: twomatrix
  DOUBLE PRECISION, Dimension(:,:),allocatable :: H,P,G,F,vecs
  DOUBLE PRECISION, Dimension(:),allocatable:: vals,ttt
      INTEGER:: NORB,NELEC,iters
      double precision:: E,E_old,E_nucc

  CONTAINS
    SUBROUTINE dic()
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
      CALL get_command_argument(1,filename)
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
      write (*,*) 'new NELEC:',NELEC,norb
      allocate(H(norb,norb))
      allocate(P(NORB,NORB))
      allocate(G(NORB,NORB))
      allocate(F(NORB,NORB))
      allocate(vecs(NORB,NORB))
      allocate(vals(NORB))
      allocate(twomatrix(norb,norb,norb,norb))
      write(*,*)'haha'
      allocate(ttt(10))
      write(*,*)'huhu'
      
      DO
       READ(LU,300,IOSTAT=ierror,IOMSG=ierrormsg) integral, my, ny, la, si
      300 Format (1x,E23.16,4I4)
        If (ierror /= 0) then
          !write(*,*) ierrormsg
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
      write(*,*)'hehe'

    END Subroutine dic
    subroutine terminate()
deallocate(twomatrix,H,P,G,F,vecs,vals)
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
    !write(*,*) vecs
    DO my=1,NORB
      DO ny=1,NORB
        P(my,ny)=0
        DO a=1,int(NELEC/2)
         P(my,ny)=P(my,ny)+2*vecs(my,a)*vecs(ny,a)
        End Do
      End Do
    End Do
    !write(*,*) P
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
    !write(*,*) F
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

  Subroutine matprint(title,A)
    DOUBLE PRECISION,DIMENSION(NORB,NORB)::A
    INTEGER::i
    CHARACTER::title
    write(*,*)title
    Do i=1,NORB
      write (*,1500) A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),A(i,8),A(i,9),A(i,10),A(i,11),A(i,12),A(i,13)
      1500 Format (' ',13F12.8)
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
    DO I=1,NORB
      Do J=1,I
        Do K=1,I
         Do L=1,K
          integral=0
          DO my=1,NORB
            DO ny=1,NORB
              Do la=1,NORB
                Do si=1,NORB
                  integral=integral+vecs(my,I)*vecs(ny,J)*vecs(la,K)*vecs(si,L)*twomatrix(my,ny,la,si)
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
  K=0
  L=0
  DO I=1,NORB
    DO J=1,I
      integral=0
      DO my=1,NORB
        DO ny=1,NORB
          integral=integral+vecs(my,I)*vecs(ny,J)*H(my,ny)
        END DO
      END DO
    WRITE(fcidumpnew,300,IOSTAT=ierror) integral, I, J, K, L
    if (ierror/=0) then
      write (*,*) imsg
    end if
    !300 Format (1x,E23.16,4I4)
    END DO
  END DO
  I=0
  J=0
  integral=E_nucc
  WRITE(fcidumpnew,300,IOSTAT=ierror,IOMSG=imsg) integral, I, J, K, L
  !300 Format (1x,E23.16,4I4)
END SUbroutine fcimaker

Subroutine RHF_control

IMPLICIT NONE
integer::iters=1,i,j
DOUBLE PRECISION::E_old=0
INTEGER::INFO
DOUBLE PRECISION, Dimension(:),allocatable::trash
integer :: lwork
integer :: lwork_max
double precision, parameter :: thresh=1d-10

allocate(trash(NORB**4))
call dic()
lwork_max=3*norb

DO
  if (iters==1) then
    call P_init
  else
    call P_matrix
  end if
  call G_matrix
  call F_matrix
  call Energy
!  if (iters==1) then 
    write(*,*)'norb',norb
    call DSYEV('V','L',norb,F,NORB,vals,trash,-1,INFO)
    lwork=min(lwork_max,int(trash(1)))
    write(*,*)'lwork',lwork
!  end if
!  call DSYEV('V','L',NORB,F,NORB,vals,trash,lwork,INFO)
  !deallocate(F)
  !write(*,*)1
  !deallocate(vecs)
  !write(*,*)2
do i=1,norb
 do j=1,norb
 !write(*,*)i,j
  vecs(i,j)=F(i,j)
!write(*,*)'jjj',i,j
 end do
end do
write(*,*)'after'
write(*,*)E
write(*,*)iters
  
  write(*,12)E,iters
  12 Format (' ','ENERGY:',F10.4,'Iteration:',I3)
  write(*,*)'test',E,iters
  if (abs(E_old-E)<(thresh)) then
    exit
  end if
  E_old=E
  iters=iters+1
  if (iters>50) then
    exit
  end if
END DO
write(*,*)-1
deallocate(ttt)

write(*,*)0
deallocate(trash)
write(*,*) 1
Call terminate()
write(*,*) 'DONE Calculating'
Call fcimaker
end Subroutine RHF_control

END MODULE RHF

PROGRAM main
  use RHF
  call RHF_control
END Program
