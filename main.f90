program main
  use fonctions
  use CSRconvert
  use fonctionCSR
  implicit none

  integer::t,i,j
  real*8,dimension(2,2)::A,G,Id

  real*8,dimension(3,3)::A1,A2,Q,R

  real*8,dimension(2)::b
  real*8,dimension(2)::x
  real*8::alpha

  integer::nlen,ncols,nlines,nelmt
  real*8,dimension(:),allocatable::AA
  integer,dimension(:),allocatable::JA,IA
  real*8,dimension(4,4)::Test

  t=2
  alpha=0.1
  do i=1,t
     do j=1,t
        G(i,j)=rand(1)
     end do
  end do


!!! Exemple d'uitilisation des fonctions NbrMat et readMat
  nlen=len('matrix/bcsstk18.mtx')
  call NbrMat('matrix/bcsstk18.mtx',nlen,ncols,nlines,nelmt)
  print*, ncols,nlines,nelmt
  allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1))
  call readMat('matrix/bcsstk18.mtx',ncols,nelmt,AA,IA,JA,nlen)


!!! Exemple d'utilisation de la fonctions DenseToCSR
!!!Test(1,:)=(/12.,4.,0.,0./)
!!!Test(2,:)=(/0.,7.,9.,-3./)
! Test(3,:)=(/1.,0.,5.,3.4/)
! Test(4,:)=(/0.,0.,-3.9,1./)
! call NbrElemt(Test,nelmt,nlines)
! print*,nelmt,nlines
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:nlines+1))
! call DenseToCSR(Test,AA,JA,IA,nelmt)
! print*, "AA vaut",AA
! print*, "JA vaut",JA
! print*, "IA vaut",IA!

  Id=reshape((/1,((0,i=1,2),1,j=1,1)/),(/t,t/))
  A=alpha*Id+ matmul(transpose(G),G)



  b(1)=-26.
  b(2)=0.125
  b=b/sqrt(sum(b*b))
  x=1.
  

  call GPO(A,b,x,t)

  print*,"GPO : ",x
  x=1.

  call residu(A,b,x,t) !! celle ci non plus

  print*,"residu : ",x

  x=1.

  call precon_residu_Jacobi(A,b,x,t)
  print*,"precon_residu_Jacobi : ",x
  x=1.

  call precon_residu_SSOR(A,b,x,t)
  print*,"precon_residu_SSOR : ",x
  x=1.

  call precon_residu_droite_Jacobi(A,b,x,t)
  print*,"precon_residu_droite_Jacobi: ",x
  x=1.

  call Jacobi(A,b,x,t) !! pas sur qu'elle marche cette méthode..

  print*,"Jacobi : ",x

  x=1

  call JacobiCSR(A,b,x,t)

  print*,"jacobi en CRS: ",x

  call GMRes(A,b,x,t)

  print*, "GMres : ", x

  call GPOCSR(A,b,x,t)
  print*, "GPO avec CRS: ", x

  A1(1,1)=2
  A1(1,2)=12
  A1(1,3)=-3
  A1(2,1)=1
  A1(2,2)=-6
  A1(2,3)=-3
  A1(3,1)=2
  A1(3,2)=0
  A1(3,3)=18


  call givens(A1,3,Q,R) !! essai de givens
  !!print*,"je suis Q ",Q
  !!print*,"je suis R:",R
  A2=matmul(Q,R)


  deallocate(AA,JA,IA)

end program main
