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


          if (abs(SUM(r*r))>max) then
             max=abs(r(i))
          end if

       k=k+1
       call write(k,sqrt(sum(r*r)),"GPO.txt")
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

    kmax=50
    eps=0.1
    nume=0.
    denom=0.
    alpha=0.
    z=0.
    Call multi_mat(r,A,x,t)

    r=b-r
    max=0
    k=0
    max=abs(sum(r*r))

    do while (k<kmax .and.  max>eps)
       call multi_mat(z,A,r,t)
       do i=1,t
          nume=nume+r(i)*z(i)
          denom=denom+z(i)*z(i)
       end do
       alpha=nume/denom

       x=x+alpha*r
       r=r-alpha*z

          if (abs(sum(r*r))>max) then
             max=r(i)
          end if

       k=k+1
    end do
  end subroutine residu




  subroutine precon_residu_Jacobi(A,b,x,t)  !!M-1Ax=M-1B aevc M diagonale
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t,t)::M
    real*8,dimension(t)::r,z,q,w
    real*8:: alpha,eps,nume,denom,max
    integer :: k, kmax,i

    call multi_mat(r,A,x,t)
    r=b-r
    m=0.
    do i=1,t
      M(i,i)=A(i,i)
      q(i)=r(i)*M(i,i)
    end do
    nume=0.
    denom=0.
    k=0
    kmax=10
    eps=0.1
    max=abs(sum(r*r))
    do while((k<kmax .and.  max>eps))

      call multi_mat(w,A,q,t)


      do i=1,t
         z(i)=w(i)*M(i,i)

         nume=nume+q(i)*z(i)
         denom=denom+z(i)*z(i)

      end do

      alpha=nume/denom

      x=x+alpha*q
      r=r-alpha*w

      q=q-alpha*z
      k=k+1


         if (abs(sum(r*r))>max) then
            max=r(i)
         end if

    end do
  end subroutine precon_residu_Jacobi



  subroutine precon_residu_SSOR(A,b,x,t)  !!M-1Ax=M-1B aevc M=(D-wR)D-1(D-wF))
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t,t)::M,L,L2,D,E,F
    real*8,dimension(t)::r,z,q,w
    real*8:: alpha,eps,nume,denom,max
    integer :: k, kmax,i

    call multi_mat(r,A,x,t)
    r=b-r
    m=0.

    do i=1,t
      do j=1,t
        if (i==j) then
          D(i,j)=A(i,i)
        else if (j<i) then
          E(i,j)=-A(i,j)
        else
          F(i,j)=-A(i,j)
        end if

      end do
    end do

    M=matmul(matmul((D-0.5*E),transpose(D)),(D-0.5*F))  !! construction du préconditionneur


    call cholesky(t,M,L,L2)
    call reso(t,L,L2,r,q)

    nume=0.
    denom=0.
    k=0
    kmax=10
    eps=0.1
    max=abs(sum(r*r))
    do while((k<kmax .and.  max>eps))
      call multi_mat(w,A,q,t)

      call cholesky(t,M,L,L2)
      call reso(t,L,L2,w,z)

      do i=1,t
         nume=nume+q(i)*z(i)
         denom=denom+z(i)*z(i)

      end do

      alpha=nume/denom
      x=x+alpha*q
      r=r-alpha*w

      q=q-alpha*z
      k=k+1

         if (abs(sum(r*r))>max) then
            max=r(i)
         end if

    end do
  end subroutine precon_residu_SSOR


 subroutine precon_residu_droite_Jacobi(A,b,x,t)
   integer,intent(in)::t !!taille des matrices
   real*8,dimension(t,t),intent(in)::A
   real*8,dimension(t),intent(in)::b
   real*8,dimension(t),intent(inout)::x
   real*8,dimension(t,t)::M
   real*8,dimension(t)::r,z,q,w
   real*8:: alpha,eps,nume,denom,max
   integer :: k, kmax,i
   M=0.
   do i=1,t
     M(i,i)=A(i,i)

     z(i)=x(i)*M(i,i)
   end do

   call multi_mat(r,matmul(A,transpose(M)),z,t)

   r=b-r

   nume=0.
   denom=0.
   k=0
   kmax=100
   eps=0.1
   max=abs(sum(r*r))


   do while((k<kmax .and.  max>eps))

     !!resoudre Mz=r
     do i=1,t
        z(i)=r(i)*M(i,i)
      end do

      w=matmul(A,z)

     do i=1,t
        nume=nume+r(i)*w(i)
        denom=denom+w(i)*w(i)
     end do

     alpha=nume/denom

     x=x+alpha*z
     r=r-alpha*w
     k=k+1


    if (abs(sum(r*r))>max) then
        max=r(i)
    end if

   end do


 end subroutine precon_residu_droite_Jacobi




!!$!! ARNOLDI
  subroutine Arnoldi(A,r,H,vm,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t,t)::v
    real*8,dimension(t,t),intent(out)::H
    real*8,dimension(t),intent(in)::r
    real*8,dimension(t), intent(out)::vm
    real*8,dimension(t)::z,q,z1
    integer :: i,j
    real*8::sum

    v(1,:)=r/sqrt(sum(r*r))   !! sqrt(sum(v*v)) revient à faire la norme :)

    do j=1,t

       do i=1,j
          call multi_mat(z1,A,v(j,:),t)
          h(i,j)=dot_product(z1,v(i,:))  !!fait le produit scalaire (à mettre de partout peut etre)
       end do
    q=0.
    do i=1,j
       q=q+h(i,j)*v(i,:)
    end do
    call multi_mat(z,A,v(j,:),t)

    z=Z-q
    H(j+1,j)=sqrt(sum(z*z))
    if (H(j+1,j)==0) then
       stop
    end if
    v(j+1,:)=z/H(j+1,j)
    vm=v(j+1,:)
    end do

end subroutine


!!GMRes
  subroutine GMRes(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::e,r,z,v,vm,sol,y
    real*8,dimension(t,t)::H,Q,Rm,G,L,L2
    real*8:: alpha,eps,nume,denom,max,beta
    integer :: k, kmax,i


    e=0.
    e(1)=1.
    call Multi_mat(z,A,x,t)
    r=b-z

    beta=sqrt(sum(r*r))
    k=0
    eps=0.01
    do while(beta>eps .and. k<kmax)
       call arnoldi(A,r,H,vm,t)
       call givens(H,t,Q,Rm) !! decompostion QR de H
!!$       !! on va calculer argmin
!!$       G=transpose(Q)*beta*e1
!!$
!!$       call cholesky(t,R,L,L2) !! resolution du systeme Rny=Gn
!!$       call solv(t,L,L2,G,y)
!!$       call Multi_mat(matmul(

       call Multi_mat(sol,H,y,t)
       y=sqrt(sum((beta*e-sol)**2))
       x=x+r*y
       r=y
       beta=sqrt(sum(r*r))
       k=k+1
    end do

    if (k>kmax) then
       print*, 'tolérence non atteinte', beta
    end if

  end subroutine GMRes





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

          Norme=sqrt(R(l1,j)**2+R(l2,j)**2)

          call mat_rot(t,l1,l2,R(L2,j)/Norme,R(l1,j)/Norme,G)
          !! on applique la rotation
          R=matmul(G,R)
          Q=matmul(Q,transpose(G))

       end do
    end do
    n=t-1

  if (R(n,n)<0)then
     R(n,n)=-R(n,n)

     Q(:,2)=-Q(:,2)
  end if

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
            res=res+B(i,j)*F(j)
         end do
         FF(i)=res
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


    subroutine write(n,x,name)
      integer,intent(in)::n
      real*8,intent(in)::x
      character*7 :: name


      open(1,file=name, form="formatted",position="append")
      print*,n,x
      write(1,*)n,x
      close(1)
    end subroutine write





  end module fonctions
