module fonctionCSR

contains

Subroutine JacobiCSR(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(:),allocatable::AA
    integer,dimension(:),allocatable::JA,IA
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::d,Xnext
    integer ::i,j,k,l,m
    real*8::sum,diag

    call CSR(A,AA,JA,IA)
    
    do k=0,10
       do i=1,t
          sum=0.
          diag=0.
          l=IA(i)
          m=IA(i+1)-1
         do j=l,m
            if (i/=JA(j)) then
               sum=sum+AA(j)*X(JA(j))
            else
               diag=AA(j)
            end if
         end do
         Xnext(i)=(b(i)-sum)/diag
       end do
       X=Xnext
    end do
    
  end subroutine JacobiCSR


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



subroutine multi_matCSR(FF,B,F,N)
      integer,intent(in)::N
      real*8,dimension(N),intent(out)::FF
      real*8,dimension(N),intent(in)::F
      real*8,dimension(N,N),intent(in)::B
      real*8,dimension(:),allocatable::AA
      integer,dimension(:),allocatable::JA,IA
      integer :: i,j,k,l
      real*8::res

      call CSR(B,AA,JA,IA)


      do i=1,N
         res=0.d0
         k=IA(i)
         l=IA(i+1)-1
         do j=k,l
            res=res+AA(j)*F(JA(j))
         end do
         FF(i)=res
      end do
 end subroutine multi_mat


subroutine GPOCSR(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::r,z
    real*8:: alpha,eps,nume,denom,max
    integer :: k, kmax,i

    call multi_matCSR(r,A,x,t)
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
       call multi_matCSR(z,A,r,t)

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
       call write(k,sqrt(sum(r*r)),"GPO.txt")
    end do
  end subroutine GPOCSR

subroutine residuCSR(A,b,x,t)
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
    Call multi_matCSR(r,A,x,t)
    r=b-r 
    max=0
    k=0
    do i=0,t
       if (abs(r(i))>max) then
          max=r(i)
       end if
    end do

    do while (k<kmax .and.  max>eps) 
       call multi_matCSR(z,A,r,t)
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
  end subroutine residuCSR


!!$!! ARNOLDI
  subroutine ArnoldiCSR(A,r,H,vm,t)
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
          call multi_matCSR(z1,A,v(j,:),t)
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
    
end subroutine ArnoldiCSR




end module fonctionCSR