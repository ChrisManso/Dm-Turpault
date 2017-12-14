module CSRconvert
implicit none

contains

subroutine Lecture(name,A,t)
  character(len=19),intent(in)::name
  integer,intent(out)::t
  real*8,dimension(:,:),allocatable,intent(out)::A
  integer::Ncols,Nlins,Nelmt,n,j,i
  real*8::a1


  open(1,file=trim(name))
  read(1,*)
  read(1,*) Ncols, Nlins,Nelmt
  Allocate(A(Ncols,Ncols))
  t=Ncols
A=0.d0
n=1
  do while (n<Nelmt)!Nelmt)
    read(1,*) i,j,a1
    A(i,j)=a1


    n=n+1
  end do
  close(1)

end subroutine Lecture



  !!! Recuperation du nombre de lignes,colonnes et du nombre d'élément du fichier Matrix
  ! subroutine NbrMat(name,nlen,Ncols,Nlins,Nelmt)
  !   integer,intent(in)::nlen
  !   character(len=nlen),intent(in)::name
  !   integer,intent(out)::Ncols,Nlins,Nelmt
  !   open(1,file=trim(name))
  !   read(1,*)
  !   read(1,*) Ncols, Nlins,Nelmt
  !   close(1)
  !
  ! end subroutine NbrMat


  !!! Lecture du fichier Matrix
  ! subroutine readMat(name,ncols,Nelmt,AA,IA,JA,nlen)
  !   integer,intent(in)::nlen
  !   character(len=nlen),intent(in)::name
  !   integer,intent(in)::ncols,Nelmt
  !   real*8,dimension(1:Nelmt),intent(inout)::AA
  !   integer,dimension(1:Nelmt),intent(inout)::JA
  !   integer,dimension(1:ncols+1),intent(inout)::IA
  !   integer::n,j,i,line,Elmt_line
  !   real*8::a
  !
  !   open(1,file=trim(name))
  !   read(1,*)
  !   read(1,*)
  !   n=1
  !   IA(1)=1
  !   line=1
  !   Elmt_line=0
  !
  !   do while (n<Nelmt)!Nelmt)
  !     read(1,*) j,i,a
  !     AA(n)=a
  !     JA(n)=j
  !     if (i==line) then
  !       Elmt_line=Elmt_line+1
  !     else
  !       line = i
  !       Elmt_line=Elmt_line+1
  !       IA(j)=Elmt_line
  !       !print*, i, Elmt_line
  !     end if
  !
  !     n=n+1
  !   end do
  !   close(1)
  ! end Subroutine readMat
  !
  !
  !
  ! !!!Donne la taille des matrices CSR (nombre d'éléménets non nul et nonbre de ligne)
  ! subroutine NbrElemt(A,nbr_elmt,n)
  !   real*8,dimension(:,:),intent(in)::A
  !   integer,intent(out)::n,nbr_elmt
  !   integer :: i,j
  !   n=size(A,1)
  !   nbr_elmt=0
  !
  !   do i=1,n
  !     do j=1,n
  !       if (abs(A(i,j))>10**(-7)) then
  !         nbr_elmt=nbr_elmt+1
  !       end if
  !     end do
  !   end do
  ! end subroutine NbrElemt
  !
  !
  ! !!! Convertion d'une matrice dense en format CSR
  ! Subroutine DenseToCSR(A,AA,JA,IA,nbr_elmt)
  !   real*8,dimension(:,:),intent(in)::A
  !   integer,intent(in)::nbr_elmt
  !   real*8,dimension(1:nbr_elmt),intent(inout)::AA
  !   integer,dimension(1:nbr_elmt),intent(inout)::JA
  !   integer,dimension(1:size(A,1)+1),intent(inout)::IA
  !   integer::i,j,n,Elmt_line
  !   n=size(A,1)
  !   Elmt_line=1 !Compte le nombre d'éléments non nul par ligne --> IA
  !   IA(1)=1
  !   do i=1,n
  !     do j = 1,n
  !       if (abs(A(i,j))>10**(-7)) then
  !         Elmt_line=Elmt_line+1
  !         AA(Elmt_line-1)=A(i,j)
  !         JA(Elmt_line-1)=j
  !       end if
  !     end do
  !     IA(i+1)=Elmt_line
  !   end do
  !
  !   !resize(AA,(/nbr_elmt/))
  !   !reshape(JA, shape=(/0:nbr_elmt/))
  !
  ! End subroutine DenseToCSR

End module CSRconvert
