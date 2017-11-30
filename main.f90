program main
  use fonctions
  implicit none

  integer::t,i,j
  real*8,dimension(2,2)::A,G,Id

  real*8,dimension(3,3)::A1,A2,Q,R
  
  real*8,dimension(2)::b
  real*8,dimension(2)::x
  real*8::alpha

t=2
  alpha=0.1
  do i=1,t
     do j=1,t
        G(i,j)=rand(1)
     end do
  end do
  
Id=reshape((/1,((0,i=1,2),1,j=1,1)/),(/t,t/))        
 A=alpha*Id+ matmul(transpose(G),G)
     


  b(1)=-26.
  b(2)=0.125
  b=b/sqrt(sum(b*b))
  x=1.


  call GPO(A,b,x,t)

  print*,"GPO : ",x
  x=1.
  call residu(A,b,x,t)
  print*,"residu : ",x

  x=1.
  call Jacobi(A,b,x,t) !! pas sur qu'elle marche cette méthode.. 
  print*,"Jacobi : ",x

!!$  x=1
!!$  call GMRes(A,b,x,t)
!!$  print*, "GMres : ", x


A1(1,1)=2
A1(1,2)=12
A1(1,3)=-3
A1(2,1)=1
A1(2,2)=-6
A1(2,3)=-3
A1(3,1)=2
A1(3,2)=0
A1(3,3)=18


  call givens(A1,3,Q,R)
print*,"je suis Q ",Q
print*,"je suis R:",R
A2=matmul(Q,R)
print*,A2

end program main
    
