module fonctions

contains

!! JACOBI
  Subroutine Jacobi(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::d,Xnext
    integer ::i,j,k
    real*8::sum

    do k=0,10
       do i=1,t
          sum=0.
          do j=1,t
             if (i/=j) then
                sum=sum+A(i,j)*X(j)
             end if
          end do
          Xnext(i)=(b(i)-sum)/A(i,i)
       end do
       X=Xnext
    end do
  end subroutine Jacobi


!! GRADIENT A PAS OPTIMAL
  subroutine GPO(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::r,z
    real*8:: alpha,eps,nume,denom,max
    integer :: k, kmax,i

    call multi_mat(r,A,x,t)
    r=b-r

    k=0
    kmax=5
    eps=0.1
    nume=0.
    denom=0.
    alpha=0.
    z=0.
    do i=0,t
       if (abs(r(i))>max) then
          max=abs(r(i))
       end if
    end do
    do while (k<kmax .and. max>eps) 
       call multi_mat(z,A,r,t)

       do i=1,t
          nume=nume+r(i)**2
          denom=denom+z(i)*r(i)
       end do
       alpha = nume/denom

       x=x+alpha*r
       r=r-alpha*z

       do i=0,t
          if (abs(r(i))>max) then
             max=abs(r(i))
          end if
       end do
       k=k+1
    end do
  end subroutine GPO



!! RESIDU OPTIMAL
  subroutine residu(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::r,z
    real*8:: alpha,eps,nume,denom,max
    integer :: k, kmax,i

    kmax=5
    eps=0.1
    nume=0.
    denom=0.
    alpha=0.
    z=0.
    Call multi_mat(r,A,x,t)
    r=b-r 
    max=0
    k=0
    do i=0,t
       if (abs(r(i))>max) then
          max=r(i)
       end if
    end do

    do while (k<kmax .and.  max>eps) 
       call multi_mat(z,A,r,t)
       do i=1,t
          nume=nume+r(i)*z(i)
          denom=denom+z(i)*z(i)
       end do
       alpha=nume/denom

       x=x+alpha*r
       r=r-alpha*z
       do i=0,t
          if (abs(r(i))>max) then
             max=r(i)
          end if
       end do
       k=k+1
    end do
  end subroutine residu



!! ARNOLDI
!!$  subroutine Arnoldi(A,v,H,t)
!!$    integer,intent(in)::t !!taille des matrices
!!$    real*8,dimension(t,t),intent(in)::A
!!$    real*8,dimension(t,t),intent(inout)::v
!!$    real*8,dimension(t,t),intent(out)::H
!!$    real*8,dimension(t)::z,q
!!$    integer :: i,j 
!!$    real*8::sum
!!$
!!$    v(1,:)=v/sqrt(sum(v(1,:)*v(1,:)))   !! sqrt(sum(v*v)) revient à faire la norme :) 
!!$
!!$    do j=1,t
!!$       q=0.
!!$       do i=1,j
!!$          
!!$          h(i,j)=A*v(j,:)*v(i,:)
!!$       end do
!!$    q=0
!!$    do i=1,j
!!$       q=q+h(i,j)*v(i,:)
!!$    end do
!!$    call multi_mat(z,A,v(j,:),t)
!!$   
!!$    z=Z-q
!!$    H(j+1,j)=sqrt(sum(z*z))
!!$    if (H(j+1,j)==0) then
!!$       stop
!!$    end if 
!!$    v(j+1)=z/H(j+1,j)    
!!$    end do 
!!$
!!$end subroutine
!!$
!!$
!!$!!GMRes
!!$  subroutine GMRes(A,b,x,t)
!!$    integer,intent(in)::t !!taille des matrices
!!$    real*8,dimension(t,t),intent(in)::A
!!$    real*8,dimension(t),intent(in)::b
!!$    real*8,dimension(t),intent(inout)::x
!!$    real*8,dimension(t)::e,r,z,v
!!$    real*8,dimension(t,t)::H
!!$    real*8:: alpha,eps,nume,denom,max,beta,y
!!$    integer :: k, kmax,i
!!$    
!!$    
!!$    e=0.
!!$    e(1)=1.
!!$    call Multi_mat(z,A,x,t)
!!$    r=b-z
!!$   
!!$    beta=sqrt(sum(r*r))
!!$    k=0
!!$    eps=0.01
!!$    do while(beta>eps .and. k<kmax)
!!$       call arnoldi(A,r,H,t)
!!$       call Multi_mat(z,H,y,t) !! gros problème là, on ne sait pas ce que c'est y 
!!$       y=sqrt(sum((beta*e-H*y)**2)) !! c'est quoi le premier y ? :/ 
!!$       x=x+r*y
!!$       r=y
!!$       beta=sqrt(sum(r*r))
!!$       k=k+1
!!$    end do
!!$    
!!$    if (k>kmax) then 
!!$       print*, 'tolérence non atteinte', beta 
!!$    end if
!!$    
!!$  end subroutine GMRes
       



  Subroutine givens(A,t,Q,R) !! givens marche :D 
    integer,intent(in)::t
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t,t),intent(out)::Q,R
    real*8,dimension(t,t)::G

    integer :: i,j,L1,L2,n,k
    real*8::Norme
    R=A
    G=0.
    Q=0.
    do k=1,t
       Q(k,k)=1.
    end do

    do j=1,t-1
       do i=t,j+1,-1
          l1=i !! ligne qu'on veut annuler
          l2=i-1
          print*,R
          Norme=sqrt(R(l1,j)**2+R(l2,j)**2)
          print*,Norme,"  la norme "
          call mat_rot(t,l1,l2,R(L2,j)/Norme,R(l1,j)/Norme,G)
          !! on applique la rotation
          R=matmul(G,R)
          Q=matmul(Q,transpose(G))
          print*,j,i
       end do
    end do
    n=t-1
    print*,Q
  if (R(n,n)<0)then
     R(n,n)=-R(n,n)
     print*, "je suis là"
     Q(:,2)=-Q(:,2)
  end if
    print*,"hello"
  end subroutine givens

subroutine mat_rot(t,i,j,c,s,M)
  integer,intent(in)::i,j,t
  real*8,intent(in)::c,s
  real*8,dimension(t,t)::M
  integer::k
  m=0.
  do k=1,t
     M(k,k)=1.
  end do
  M(i,i)=c
  M(j,j)=c
  M(i,j)=-s
  M(j,i)=s
  print*,"voilà M",M

end subroutine mat_rot







!! MODULE COMPLEMENTAIRE

    subroutine multi_mat(FF,B,F,N)
      integer,intent(in)::N
      real*8,dimension(N),intent(out)::FF
      real*8,dimension(N),intent(in)::F
      real*8,dimension(N,N),intent(in)::B
      integer :: i,j
      real*8::res

      do i=1,N
         res=0.d0
         do j=1,N
            FF(i)=res+B(j,i)*F(j)
         end do
      end do
    end subroutine multi_mat


    subroutine cholesky(n,A,L,L2)
      integer,intent(in)::n
      real*8,dimension(n,n)::A,L,L2
      integer::i,j,k
      real*8 :: res,res2

      L=0d0
      L2=0d0
      do i=1,n
         res=0.
         res2=0.

         if (i>1) then
            do k=1,i-1
               res=res+L(i,k)**2
            end do
         end if

         L(i,i)=sqrt(A(i,i)-res)
         L2(i,i)=sqrt(A(i,i)-res)

         do j=i+1,n
            res2=0.
            do k=1,i-1
               res2=res2+L(i,k)*L(j,k)
            end do
            L(j,i)=(A(i,j)-res2)/L(i,i)
            L2(i,j)=(A(i,j)-res2)/L(i,i)
         end do
      end do
    end subroutine cholesky

    subroutine reso(n,L,L2,F,X)
      implicit none
      integer,intent(in)::n
      real*8,dimension(n,n),intent(in)::L,L2
      real*8,dimension(n),intent(in)::F
      real*8,dimension(n)::y
      real*8,dimension(n),intent(out)::x
      integer::i,j,k

      y(1)=f(1)/L(1,1)
      do i=2,n   
         y(i)=(F(i)-L(i,i-1)*y(i-1))/L(i,i)
      end do

      x(n)=y(n)/L2(n,n)
      do i=n-1,1,-1
         x(i)=(y(i)-L2(i,i+1)*y(i+1))/L2(i,i)
      end do
    end subroutine reso

  end module fonctions
