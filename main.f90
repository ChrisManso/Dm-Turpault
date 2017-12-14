program main
  use fonctions
  use CSRconvert
  use fonctionCSR
  implicit none

  real*8,dimension(:,:),allocatable::A,G,Id
  real*8,dimension(:),allocatable::AA,x,x0,b
  integer,dimension(:),allocatable::JA,IA

  real*8,dimension(3,3)::A1,A2,Q,R
  real*8,dimension(4,4)::Test

  real*8::nombre,t1,t2,alpha
  integer::nlen,ncols,nlines,nelmt,t,i,j

  !!! Exemple d'uitilisation des fonctions NbrMat et readMat
  ! nlen=len('matrix/bcsstk18.mtx')
  ! call NbrMat('matrix/bcsstk18.mtx',nlen,ncols,nlines,nelmt)
  ! !print*, ncols,nlines,nelmt
  ! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1))
  ! call readMat('matrix/bcsstk18.mtx',ncols,nelmt,AA,IA,JA,nlen)


  !!! Exemple d'utilisation de la fonctions DenseToCSR
  ! Test(1,:)=(/12.,4.,0.,0./)
  ! Test(2,:)=(/0.,7.,9.,-3./)
  ! Test(3,:)=(/1.,0.,5.,3.4/)
  ! Test(4,:)=(/0.,0.,-3.9,1./)
  ! call NbrElemt(Test,nelmt,nlines)
  ! !print*,nelmt,nlines
  ! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:nlines+1))
  ! call DenseToCSR(Test,AA,JA,IA,nelmt)
  ! !print*, "AA vaut",AA
  ! !print*, "JA vaut",JA
  ! !print*, "IA vaut",IA

  !!! PARAMETRE
  t=50 !/!\ choisie la dimension
  Allocate(x(1:t),x0(1:t),b(1:t),A(1:t,1:t),G(1:t,1:t),Id(1:t,1:t))
  x0=1.


  !!! Création de la matrice An
    t=5

     !/!\ choisie la dimension

    ! Allocate(x(1:t),b(1:t),A(1:t,1:t),G(1:t,1:t),Id(1:t,1:t))
    alpha=1
    do i=1,t
       do j=1,t
         call random_number(nombre)
          G(i,j)=nombre/real(t)
    end do
  end do

  Id=reshape((/1,((0,i=1,t),1,j=1,t-1)/),(/t,t/))
  A=alpha*Id+matmul(transpose(G),G)

  !!! Initialisation de notre membre de droite, la matrice B
  b=1.
  b=b/sqrt(sum(b*b))


  !!! ALGORYTHME
  x=x0  !Réinitialisation du vecteur d'entré
  t1= wtime ( )
  call GPO(A,b,x,t)
  t2= wtime ( )
  print*,"temps de GPO =",t2-t1
  Print*," "


  x=x0  !Réinitialisation du vecteur d'entré
  t1= wtime ( )
  call residu(A,b,x,t) !! celle ci non plus
  t2=wtime ( )
  print*,"temps residu =",t2-t1
  Print*," "


  x=x0  !Réinitialisation du vecteur d'entré
  t1= wtime ( )
  call Jacobi(A,b,x,t)
  t2=wtime ( )
  print*,"temps de jacobi=",t2-t1
  Print*," "

  x=x0  !Réinitialisation du vecteur d'entré
  call GMRes(A,b,x,t)
  ! print*,"temps de jacobi=",t2-t1
  Print*, "GMres : ", x


  ! x=x0  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  ! call GMRes(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"temps de GMRes=",t2-t1
  ! Print*," "



  !!! PRECONDITIONNEUR
    ! x=x0  !Réinitialisation du vecteur d'entré
    ! t1= wtime ( )
    ! call precon_residu_SSOR(A,b,x,t)
    ! t2=wtime ( )
    ! print*,"temps de precon_residu_SSOR =",t2-t1
    ! Print*," "
    !
    !
    ! x=x0  !Réinitialisation du vecteur d'entré
    ! t1= wtime ( )
    ! call precon_residu_Jacobi(A,b,x,t)
    ! t2=wtime ( )
    ! print*,"temps de precon_residu_Jacobi = ",t2-t1
    ! Print*," "
    !
    !
    ! x=x0  !Réinitialisation du vecteur d'entré
    ! t1= wtime ( )
    ! call precon_residu_droite_Jacobi(A,b,x,t)
    ! t2=wtime ( )
    ! print*,"temps de precon_residu_droite_Jacobi =",t2-t1
    ! Print*," "
    !

  !!! Liberation de la mémoire
    deallocate(x,x0,b,A,G,Id)



!!! ALGORYTHME CSR
nlen=len('matrix/bcsstk18.mtx')
call NbrMat('matrix/bcsstk18.mtx',nlen,ncols,nlines,nelmt)
allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
call readMat('matrix/bcsstk18.mtx',ncols,nelmt,AA,IA,JA,nlen)
b=1.
b=b/sqrt(sum(b*b))
x0=1.

x=x0  !Réinitialisation du vecteur d'entré
t1= wtime ( )
call GPOCSR(AA,IA,JA,b,x,ncols)
t2=wtime ( )
print*,"temps de GPO au format CSR=",t2-t1
Print*," "


x=x0  !Réinitialisation du vecteur d'entré
t1= wtime ( )
call residuCSR(AA,JA,IA,b,x,ncols)
t2=wtime ( )
print*,"temps du ResiduMin au format CSR=",t2-t1
Print*," "


x=x0  !Réinitialisation du vecteur d'entré
t1= wtime ( )
call JacobiCSR(AA,JA,IA,b,x,ncols)
t2=wtime ( )
print*,"temps de Jacobi au format CSR=",t2-t1
Print*," "

deallocate(AA,JA,IA,b,x,x0)







end program main
