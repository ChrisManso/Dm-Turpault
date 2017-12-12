module fonctions

contains

!! JACOBI
  Subroutine Jacobi(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t,t)::D,N
    real*8,dimension(t)::r,Xnext
    integer ::i,j,k,kmax
    real*8::max,sigma, eps

    r=b-matmul(A,x)
    max=abs(sum(r*r))
    kmax=1000
    eps=0.0000001
    k=0
    xnext=0.
      do while (k<kmax .and. max>eps)
       do i=1,t
          sigma=0.
          do j=1,t
             if (i/=j) then
                sigma=sigma+A(i,j)*X(j)
             end if
          end do

          Xnext(i)=1./A(i,i)*(b(i)-sigma)
       end do
       X=Xnext

       r=b-matmul(A,X)

       if (abs(SUM(r*r))<max) then
          max=abs(sum(r*r))
       end if
       k=k+1
    end do
    call write(k,sqrt(sum(r*r)),"Jacobi.txt")
    print*,"jacobi",max,k
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
    kmax=100

    nume=0.
    denom=0.
    alpha=0.
    z=0.
    eps=0.0000001

    max= abs(sum(r*r))
    do while (k<kmax .and. max>eps)
       call multi_mat(z,A,r,t)

       do i=1,t
          nume=nume+r(i)**2
          denom=denom+z(i)*r(i)
       end do

       alpha = nume/denom

       x=x+alpha*r
       r=r-alpha*z


        if (abs(SUM(r*r))<max) then
           max=abs(sum(r*r))
        end if

       k=k+1
       call write(k,sqrt(sum(r*r)),"GPOpti.txt")
    end do
print*,"GPO = ", max,k
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

    kmax=10000
    eps=0.0000001
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

          if (abs(sum(r*r))<max) then
             max=abs(sum(r*r))
          end if

       k=k+1
       call write(k,sqrt(sum(r*r)),"ResMin.txt")
    end do
print*,"residu = ",max,k
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
    kmax=1000000
    eps=0.00001
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


         if (abs(sum(r*r))<max) then
            max=abs(sum(r*r))
         end if
         call write(k,sqrt(sum(r*r)),"ResJac.txt")
    end do
print*,"precon_residu_Jacobi = ",max,k
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

    E=0.
    D=0.
    F=0.
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


    M=matmul(matmul((D-0.8*E),transpose(D)),(D-0.8*F))  !! construction du préconditionneur

    call cholesky(t,M,L,L2)
    call reso(t,L,L2,r,q)

    nume=0.
    denom=0.
    k=0
    kmax=1000000
    eps=0.00001
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

         if (abs(sum(r*r))<max) then
            max=abs(sum(r*r))
         end if
         call write(k,sqrt(sum(r*r)),"ResSSO.txt")
    end do

print*,"precon_residu_SSOR = ",max,k
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
   kmax=1000000
   eps=0.00001
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


    if (abs(sum(r*r))<max) then
        max=abs(sum(r*r))
    end if

   end do
   print*,"precon_residu_droite_Jacobi = ",max,k

 end subroutine precon_residu_droite_Jacobi




!!$!! ARNOLDI
  subroutine Arnoldi(A,r,H,vm,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t+1,t)::v
    real*8,dimension(t+1,t),intent(out)::H
    real*8,dimension(t),intent(in)::r
    real*8,dimension(t), intent(out)::vm
    real*8,dimension(t)::z,q,z1
    integer :: i,j,k

    v(1,:)=r/sqrt(sum(r*r))   !! sqrt(sum(v*v)) revient à faire la norme :)

    do j=1,t
    z1 = matmul(A,v(j,:))
    do i=1,j
          h(i,j)=dot_product(z1,v(i,:))   !!fait le produit scalaire (à mettre de partout peut etre)
          z1=z1-h(i,j)*v(i,:)
       end do
!     q=0.
!     do k=1,j
!       q = q + h(i,j)*v(k,:)
!     end do
!     z=0.
    h(j+1,j)=sqrt(sum(z1*z1))
    if (h(j+1,j)==0) then
       stop
    end if
    v(j+1,:)=z1/h(j+1,j)
    vm=v(j+1,:)
    end do

end subroutine


!!GMRes
  subroutine GMRes(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::e,r,z,Vm,sol,y,G
    real*8,dimension(t,t)::H,Q,L,L2,Rm,jdr
    real*8,dimension(t+1,t) :: Hm
    real*8, dimension(t) :: v
    real*8:: alpha,eps,nume,denom,beta,som
    integer :: k, kmax,i,ui,uj


    e=0.
    e(1)=1.
    r=b-matmul(A,x)
    beta=sqrt(sum(r*r))
    k=0
    kmax=10
    eps=0.001
    
    do while(beta>eps .and. k<kmax)
       call Arnoldi(A,r,Hm,v,t)
       H=Hm(1:t,:)
       Vm=v(1:t)
       Rm=0.
       Q=0.
       call givens(H,t,Q,Rm) !! decompostion QR de H
       !! on va calculer argmin
       G=beta * matmul(transpose(Q),e)
!        print*,
!        L=0.
!        L2=0.
!        call cholesky(t,Rm,L,L2) !! resolution du systeme Rny=Gn
!        do ui=1,2
!         print*, jdr(i,:)
!        end do
!        print*,"oui oui"
        !print*, size(h)
       print*, "oui oui"
        y(t) = G(t)/Rm(t,t)
        do ui=t-1,1,-1
          som=0.
          do uj=ui+1,t
            som = som+G(uj)*Rm(uj,ui)
          end do
          y(ui) = (G(ui)-som)/Rm(ui,ui)
       end do
       r = beta*e-matmul(H,y)
       y = sqrt(sum(y*y))
       x = x+y*Vm
       beta = sqrt(sum(r*r))
       k = k+1
    end do

    if (k>kmax) then
       !print*, 'tolérence non atteinte', beta
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
      character*10 :: name
      if (n==1) then
        open(1,file=name,form="formatted")
      else

        open(1,file=name, form="formatted",position="append")
      end if
      !print*,n,x
      write(1,*)n,x
      close(1)
    end subroutine write





  end module fonctions
