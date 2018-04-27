MODULE manz_RHF
  implicit NONE
  private
  public init,P_init,G_matrix,Energy,RHF_control,manz_RHF_controlinit
      DOUBLE PRECISION,allocatable,DImension(:,:,:,:):: twomatrix
      DOUBLE PRECISION,allocatable,DImension(:,:):: H
      Double Precision::E_Nucc
      INTEGER::NORB,NELEC,iunit_dump,MS2,ISYM
      Double Precision,Dimension(100)::ORBSYM
      CHARACTER(len=50)::filename_out
      INTEGER::addelec

  CONTAINS
    SUBROUTINE init
      IMPLICIT NONE
      Integer, PARAMETER:: MAX_size=100
      INTEGER:: my,ny,la,si,n,index1,index2,ierror,elecadd
      Character(len=24):: filename
      Character(len=128):: ierrormsg
      NAMELIST/FCI/NORB,NELEC,MS2,ORBSYM,ISYM
      
      rewind(Unit=iunit_dump)
      READ(UNIT=iunit_dump,NML=FCI)
      NELEC=NELEC+addelec
      write (*,*) 'new NELEC:',NELEC,norb
    END SUBROUTINE init 


  Subroutine read_integral
      INTEGER:: NORB_store, NELEC_store, MS2, ISYM, ierror
      INTEGER, PARAMETER :: MAX_size=100,LU=99
      INTEGER:: my,ny,la,si
      DOUBLE PRECISION:: integral
      Character(len=128):: ierrormsg
      DOUBLE PRECISION,DImension(MAX_Size)::ORBSYM
      Character(len=24):: filename
      NAMELIST/FCI/NORB,NELEC,MS2,ORBSYM,ISYM
      CALL get_command_argument(1,filename)
      DO
       READ(iunit_dump,300,IOSTAT=ierror,IOMSG=ierrormsg) integral, my, ny, la, si
      300 Format (1x,E23.16,4I4)
        If (ierror /= 0) then
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
      close(unit=iunit_dump)
    END Subroutine read_integral

    Subroutine P_init(P)
      IMPLICIT NONE
      integer::i=1
      DOUBLE Precision,Dimension(NORB,NORB),intent(out)::P
      P=0
      DO i = 1,int(NELEC/2)
        P(i,i)=2
      END DO
    END Subroutine P_init

  Subroutine P_matrix(P,vecs)
    Implicit NONE
    integer::my,ny,a
    Double Precision,Dimension(NORB,NORB),intent(out)::P
    Double Precision,Dimension(NORB,NORB),intent(in)::vecs
    DO my=1,NORB
      DO ny=1,NORB
        P(my,ny)=0
        DO a=1,int(NELEC/2)
         P(my,ny)=P(my,ny)+2*vecs(my,a)*vecs(ny,a)
        End Do
      End Do
    End Do
  End Subroutine P_matrix

  Subroutine G_matrix(G,P)
    IMPLICIT NONE
    integer::my,ny,la,si
    DOUBLE Precision,Dimension(NORB,NORB),Intent(out)::G
    Double Precision,Dimension(NORB,NORB), INtent(in)::P
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

  Subroutine Energy(E,P,G)
    IMPLICIT NONE
    INTEGER:: my,ny,la,si
    DOUBLE PRECISion, intent(out)::E
    Double precision, dimension(NORB,NORB),intent(in)::P,G

    E=E_nucc
    Do my=1,NORB
      DO ny=1,NORB
        E=E+P(my,ny)*(H(my,ny)+0.5*G(my,ny))
      END Do
    END Do
  END Subroutine ENergy

  Subroutine matprint(title,A)
    DOUBLE PRECISION,DIMENSION(NORB,NORB),Intent(in)::A
    INTEGER::i
    CHARACTER::title
    write(*,*)title
    Do i=1,NORB
      write (*,1500) A(i,1),A(i,2),A(i,3),A(i,4),A(i,5),A(i,6),A(i,7),A(i,8),A(i,9),A(i,10),A(i,11),A(i,12),A(i,13)
      1500 Format (' ',13F16.8)
    End Do
  End Subroutine matprint

  Subroutine fcimaker(vecs)
    IMPLICIT NONE
    INTEGER::I,J,K,L,my,ny,la,si
    DOUBLE PRECISION:: integral
    INTEGER::ierror,fcidumpnew
    CHARACTER(len=124)::imsg
    DOUBLE Precision,Dimension(NORB,NORB),intent(in)::vecs
    NAMELIST/FCI/NORB,NELEC,MS2,ORBSYM,ISYM
    fcidumpnew=25
    OPEN(UNIT=fcidumpnew, FILE=filename_out, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror)
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
    END DO
  END DO
  I=0
  J=0
  integral=E_nucc
  WRITE(fcidumpnew,300,IOSTAT=ierror,IOMSG=imsg) integral, I, J, K, L
END SUbroutine fcimaker

Subroutine RHF_control(P,G,vecs,vals)

IMPLICIT NONE
integer::iters=1,i,j
DOUBLE PRECISION::E_old=0
DOUBLE Precision::E=0
INTEGER::INFO
DOUBLE PRECISION, Dimension(:),allocatable::work
integer :: lwork,lwork_max,deallo_STATUS
double precision, parameter :: thresh=1d-10
CHARACTER(len=124)::ERRalloc
Double precision,Dimension(:,:),allocatable,intent(inout)::P,G,vecs
DOUBLE PRECISION,Dimension(:),allocatable,intent(inout)::vals
Double Precision,Dimension(:,:),allocatable::F
call init
lwork_max=norb**2
allocate(work(lwork_max))
lwork=NORB**2
allocate(P(NORB,NORB),G(NORB,NORB),F(NORB,NORB),vecs(NORB,NORB),vals(NORB),H(NORB,NORB),twomatrix(NORB,NORB,NORB,NORB))
P=0
G=0
F=0
vecs=0
vals=0
H=0
twomatrix=0
work=0
call read_integral

DO
  if (iters==1) then
    call P_init(P)
  else
    call P_matrix(P,vecs)
  end if
  call G_matrix(G,P)
  F=H+G
  call Energy(E,P,G)
  call DSYEV('V','L',NORB,F,NORB,vals,work,lwork,INFO)
do i=1,norb
 do j=1,norb
  vecs(i,j)=F(i,j)
 end do
end do
  write(*,1100)E,iters
  1100 Format (' ','ENERGY:',F20.14,'Iteration:',I3)
  if (abs(E_old-E)<(thresh)) then
    exit
  end if
  E_old=E
  iters=iters+1
  if (iters>50) then
    exit
  end if
END DO
write(*,*) 'DONE Calculating'
Call fcimaker(vecs)
write(*,*) "DONE FCIMAKER"
deallocate(work,stat=deallo_Status,ERRMSG=ERRalloc)
deallocate(twomatrix)
deallocate(H)
deallocate(P)
deallocate(G)
deallocate(F)
deallocate(vecs)
deallocate(vals)
write(*,*)'deallocated'
end Subroutine RHF_control

Subroutine manz_RHF_controlinit(iunit_dump_stor,filename_out_stor,addelec_stor)
  Double Precision, Dimension(:,:),allocatable::P,G,vecs
  Double Precision, Dimension(:),allocatable::vals
  Integer::iunit_dump_stor,addelec_stor
  CHARACTER(len=50)::filename_out_stor
  iunit_dump=iunit_dump_stor
  filename_out=filename_out_stor
  addelec=addelec_stor
  call RHF_control(P,G,vecs,vals)
END Subroutine manz_RHF_controlinit

END MODULE manz_RHF

PROGRAM RHF_main
  use manz_RHF
  INTEGER::iunit_dump=99
  Character (len=50)::filename_in,filename_out
  INTEGER::addelec
  addelec=0
  CALL get_command_argument(1,filename_in)
  Call get_command_argument(2,filename_out)
  open (Unit=iunit_dump, File=filename_in,Status='old',action='read')
  call manz_RHF_controlinit(iunit_dump,filename_out,addelec)
END PROGRAM
