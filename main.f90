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

  real*8::nombre,t1,t2,alpha,br
  integer::nlen,ncols,nlines,nelmt,t,i,j, userChoice

  !!! PARAMETRE
  print*, "Sur quelle matrice voulez-vous travailler ?"
  print*, "1) = An (matrice aléatoire)"
  print*, "2) =  matrice de la vie réelle en dense"
  read*, userChoice



  select case(userChoice)
  case (1)
    !/!\ choisie la dimension
    print*, "Quelle taille de matrice ?"
    read*, t

    !Allocation de la mémoire
    Allocate(x(1:t),x0(1:t),b(1:t),A(1:t,1:t),Id(1:t,1:t),G(1:t,1:t))!,A(1:t,1:t))
    x0=1.

    ! Création de la matrice An
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
    x=1.


  case(2)
    call Lecture("matrix/fs_760_3.mtx",A,t)
    Allocate(x(1:t),x0(1:t),b(1:t),A(1:t,1:t))!,A(1:t,1:t))
    !!! Initialisation de notre membre de droite, la matrice B
    b=1.
    x0=1.
  case DEFAULT
    print*, "Ce choix n'est pas possible"
  end select


!-----------------------  GPO  -------------------------!

  ! x=x0  !Réinitialisation du vecteur d'entré
  !  t1= wtime ( )
  !  call GPO(A,b,x,t)
  !  t2= wtime ( )
  !
  ! print*,"temps de GPO =",t2-t1
  ! Print*," "


!---------------------- RESIDU --------------------------!

!Réinitialisation du vecteur d'entré
! x=1.
!
! t1= wtime ( )
! call residu(A,b,x,t)
! t2=wtime ( )
! print*,"temps residu =",t2-t1
! Print*," "

!-------------------- JACOBI -----------------------!
! Réinitialisation du vecteur d'entré
!   x=1.
!   t1= wtime ( )
!   call Jacobi(A,b,x,t)
!   t2=wtime ( )
!   !!print*,"temps de jacobi=",t2-t1
!   Print*," "


!-------------- GMRES -----------------!
!(Ne fonctionne pas)
!  Réinitialisation du vecteur d'entré
! x=1.
! t1= wtime ( )
! call GMRes(A,b,x,t)
! t2=wtime ( )
! print*,"temps de GMRes=",t2-t1
! Print*," "



!!! -------------------- PRECONDITIONNEUR ------------------------!!!


!--------- PRECONDIONNEUR DROITE SSORS FLEXIBLE ------- !
! Réinitialisation du vecteur d'entré
!  t1= wtime ( )
!  call precon_residu_droite_SSOR_FlexibleC(A,b,x,t)
!  t2=wtime ( )
!  print*,"le x de precon_residu_droite_SSOR_FlexibleC",x
!  print*,"temps de precon_residu_droite_SSOR_FlexibleC =",t2-t1
!  Print*," "


! Réinitialisation du vecteur d'entré
! x=1.
! t1= wtime ( )
! call precon_residu_SSOR(A,b,x,t)
! t2=wtime ( )
! print*,"temps de precon_residu_SSOR =",t2-t1
! print*,"le x de precon_residu_SSOR ",x
! Print*," "

!-----------------  Préconditionneu droite SSOR -----------------!
! x=x0  !Réinitialisation du vecteur d'entré
! t1= wtime ( )
! call precon_residu_droite_SSOR(A,b,x,t)
! t2=wtime ( )
! !print*,"le x de precon_residu_droite_SSOR",x
! print*,"temps de precon_residu_droite_SSOR =",t2-t1
! Print*," "


!---------- Préconditionneur gauche Jacobi -----------!
! Réinitialisation du vecteur d'entré
! x=1.
! t1= wtime ( )
! call precon_residu_Jacobi(A,b,x,t)
! t2=wtime ( )
! print*,"temps de precon_residu_Jacobi = ",t2-t1
! print*,"le x de precon_residu_Jacobi ",x
! Print*," "


!------------ Préconditionneur droite ----------------!
! Réinitialisation du vecteur d'entré
! t1= wtime ( )
! call precon_residu_droite_Jacobi(A,b,x,t)
! t2=wtime ( )
! print*,"temps de precon_residu_droite_Jacobi =",t2-t1
! Print*," "

!!! Liberation de la mémoire
!deallocate(x,x0,b,A,G,Id)

print*,"----------------------------------------------------------------------------------------"


!!! -----------------------ALGORYTHME CSR ------------------------!!!

!! ecuperation des mtrices et mis en forme CSR
!! 1ere matrice
! nlen=len('matrix/bcsstk18.mtx')
! call NbrMat('matrix/bcsstk18.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! call readMat('matrix/bcsstk18.mtx',ncols,nelmt,AA,IA,JA,nlen)

!!2eme matrice
! nlen=len('matrix/fidapm37.mtx')
! call NbrMat('matrix/fidapm37.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! IA=0.
! call readMat('matrix/fidapm37.mtx',ncols,nelmt,AA,IA,JA,nlen)

!!3eme matrice
! nlen=len('matrix/fs_541_4.mtx')
! call NbrMat('matrix/fs_541_4.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! call readMat('matrix/fs_541_4.mtx',ncols,nelmt,AA,IA,JA,nlen)

!! 4eme matrice
! nlen=len('matrix/fs_760_3.mtx')
! call NbrMat('matrix/fs_760_3.mtx',nlen,ncols,nlines,nelmt)
! allocate (AA(1:nelmt),JA(1:nelmt),IA(1:ncols+1),b(1:ncols),x(1:ncols),x0(1:ncols))
! call readMat('matrix/fs_760_3.mtx',ncols,nelmt,AA,IA,JA,nlen)
! print*,"hello"
! b=1.
! b=b/sqrt(sum(b*b))
! x0=1.


! ---------------GPOCSR------------------ !
! x=x0  !Réinitialisation du vecteur d'entré
! t1= wtime ( )
! print*,"hello"
! call GPOCSR(AA,JA,IA,b,x,ncols)
! t2=wtime ( )
! print*,"temps de GPO en CSR=",t2-t1
! Print*," "
!

!! ------------ RESIDU EN CSR -----------!!
! x=x0  !Réinitialisation du vecteur d'entré
! t1= wtime ( )
! call residuCSR(AA,JA,IA,b,x,ncols)
! t2=wtime ( )
! print*,"temps de ResMin en CSR=",t2-t1
! Print*," "

!! ------------- JACOBICSR -------------- !!
! x=x0  !Réinitialisation du vecteur d'entré
! t1= wtime ( )
! call JacobiCSR(AA,JA,IA,b,x,ncols)
! t2=wtime ( )
! print*,"temps de Jacobi en CSR=",t2-t1
! Print*," "

!! désallocation de la mémoire
!deallocate(AA,JA,IA,b,x,x0)


!deallocate(x,x0,b,A,G,Id)


end program main
