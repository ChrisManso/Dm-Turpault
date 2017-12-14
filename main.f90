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



  !!! PARAMETRE

  t=20 !/!\ choisie la dimension
  Allocate(x(1:t),x0(1:t),b(1:t),A(1:t,1:t),G(1:t,1:t),Id(1:t,1:t))
  x0=1.

  !!! Création de la matrice An

     !/!\ choisie la dimension

    ! Allocate(x(1:t),b(1:t),A(1:t,1:t),G(1:t,1:t),Id(1:t,1:t))
    alpha=1
    do i=1,t
       do j=1,t
         call random_number(nombre)
          G(i,j)=nombre/15
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
  call residu(A,b,x,t)
t2=wtime ( )
  print*,"temps residu =",t2-t1
  !!print*,"residu : ",x
  Print*," "
  x=1.





  x=x0  !Réinitialisation du vecteur d'entré
  t1= wtime ( )
  call Jacobi(A,b,x,t)
  t2=wtime ( )
  !!print*,"temps de jacobi=",t2-t1
  Print*," "


   !!Réinitialisation du vecteur d'entré
  ! call GMRes(A,b,x,t)
  !
  ! Print*, "GMres : ", x


  ! x=x0  Réinitialisation du vecteur d'entré
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

  !  deallocate(x,x0,b,A,G,Id)





 



print*,"----------------------------------------------------------------------------------------"
!!! ALGORYTHME CSR
! nlen=len('matrix/bcsstk18.mtx')
! call NbrMat('matrix/bcsstk18.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! call readMat('matrix/bcsstk18.mtx',ncols,nelmt,AA,IA,JA,nlen)

! nlen=len('matrix/fidapm37.mtx')
! call NbrMat('matrix/fidapm37.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! IA=0.
! call readMat('matrix/fidapm37.mtx',ncols,nelmt,AA,IA,JA,nlen)

! nlen=len('matrix/fs_541_4.mtx')
! call NbrMat('matrix/fs_541_4.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! call readMat('matrix/fs_541_4.mtx',ncols,nelmt,AA,IA,JA,nlen)

nlen=len('matrix/fs_760_3.mtx')
call NbrMat('matrix/fs_760_3.mtx',nlen,ncols,nlines,nelmt)
allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
call readMat('matrix/fs_760_3.mtx',ncols,nelmt,AA,IA,JA,nlen)

b=1.
b=b/sqrt(sum(b*b))
x0=1.

x=x0  !Réinitialisation du vecteur d'entré
t1= wtime ( )
call GPOCSR(AA,JA,IA,b,x,ncols)
t2=wtime ( )
print*,"temps de GPO en CSR=",t2-t1
Print*," "


x=x0  !Réinitialisation du vecteur d'entré
t1= wtime ( )
call residuCSR(AA,JA,IA,b,x,ncols)
t2=wtime ( )
print*,"temps de ResMin en CSR=",t2-t1
Print*," "

x=x0  !Réinitialisation du vecteur d'entré
t1= wtime ( )
call JacobiCSR(AA,JA,IA,b,x,ncols)
t2=wtime ( )
print*,"temps de Jacobi en CSR=",t2-t1
Print*," "

!!! PRECONDITIONNEUR


  ! x=x0  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  !
  ! call precon_residu_SSOR(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"le x de precon_residu_SSOR :",x
  ! !print*,"temps de precon_residu_SSOR =",t2-t1
  ! Print*," "
  !
  !
  !
  !
  ! x=x0  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  !
  ! call precon_residu_Jacobi(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"le x de precon_residu_Jacobi :",x
  ! !print*,"temps de precon_residu_Jacobi =",t2-t1
  ! Print*," "



  ! x=1.  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  ! call precon_residu_droite_Jacobi(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"le x de precon_residu_droite_Jacobi",x
  ! !!print*,"temps de precon_residu_droite_Jacobi =",t2-t1
  ! Print*," "
  !
  !
  ! x=x0  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  ! call precon_residu_droite_SSOR(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"le x de precon_residu_droite_SSOR",x
  ! print*,"temps de precon_residu_droite_SSOR =",t2-t1
  ! Print*," "
  !
  !
  ! x=x0  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  ! call precon_residu_droite_SSOR_FlexibleB(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"le x de precon_residu_droite_SSOR_FlexibleB",x
  ! print*,"temps de precon_residu_droite_SSOR_FlexibleB =",t2-t1
  ! Print*," "


  ! x=x0  !Réinitialisation du vecteur d'entré
  ! t1= wtime ( )
  ! call precon_residu_droite_SSOR_FlexibleC(A,b,x,t)
  ! t2=wtime ( )
  ! print*,"le x de precon_residu_droite_SSOR_FlexibleC",x
  ! print*,"temps de precon_residu_droite_SSOR_FlexibleC =",t2-t1
  ! Print*," "



end program main
