module fonctionCSR
use fonctions
use CSRconvert
contains

Subroutine JacobiCSR(AA,JA,IA,b,x,t)
    integer,intent(in)::t !!taille des matrices

    real*8,dimension(t),intent(in)::b
    real*8,dimension(:),allocatable,intent(in)::AA
    integer,dimension(:),allocatable,intent(in)::JA,IA
    real*8,dimension(t),intent(inout),intent(in)::x
    real*8,dimension(t)::d,Xnext,r
    integer ::i,j,k,l,m, Nb_elem,Nb_li,kmax,max

    real*8::sigma,diag,epsn,norme


    call NbrElemt(A,Nb_elem,Nb_li)
    Allocate(AA(1:Nb_elem), JA(1:Nb_elem),IA(1:Nb_li+1))
    call DenseToCSR(A,AA,JA,IA,Nb_elem)


    call multi_matCSR(r,AA,IA,JA,x,t)
    r=b-r
    max=abs(sum(r*r))
    kmax=1000
    eps=0.00001
    k=0
    xnext=0.
    do while (k<kmax .and. max>eps)
       do i=1,t
          sigma=0.
          diag=0.
          l=IA(i)
          m=IA(i+1)-1
         do j=l,m
            if (i/=JA(j)) then
               sigma=sigma+AA(j)*X(JA(j))
            else
               diag=AA(j)
            end if
         end do
         Xnext(i)=(b(i)-sigma)/diag
       end do
       X=Xnext
       call multi_mat(r,AA,IA,JA,x,t)
       r=b-r
       norme=abs(sum(r*r))
       if (norme)<max) then
          max=norme
       end if
       k=k+1
    end do

    deallocate(AA,IA,JA)

  end subroutine JacobiCSR



subroutine multi_matCSR(FF,B,F,N)
      integer,intent(in)::N
      real*8,dimension(N),intent(out)::FF
      real*8,dimension(N),intent(in)::F
      real*8,dimension(N,N),intent(in)::B
      real*8,dimension(:),allocatable::AA
      integer,dimension(:),allocatable::JA,IA
      integer :: i,j,k,l, Nb_elem, Nb_li
      real*8::res


      call NbrElemt(B,Nb_elem,Nb_li)
      Allocate(AA(1:Nb_elem), JA(1:Nb_elem),IA(1:Nb_li+1))
      call DenseToCSR(B,AA,JA,IA, Nb_elem)


      do i=1,N
         res=0.d0
         k=IA(i)
         l=IA(i+1)-1
         do j=k,l
            res=res+AA(j)*F(JA(j))
         end do
         FF(i)=res
      end do

      deallocate(AA,JA,IA)
 end subroutine multi_matCSR


subroutine GPOCSR(A,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t,t),intent(in)::A
    real*8,dimension(t),intent(in)::b
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::r,z
    real*8:: alpha,eps,nume,denom,max,norme
    integer :: k, kmax,i


    call multi_matCSR(r,AA,IA,JA,x,t)
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

       norme=abs(sum(r*r))
        if (norme<max) then
           max=norme
        end if

       k=k+1
       call write(k,sqrt(sum(r*r)),"GPOpti.txt")
    end do
print*,"GPO = ", max,k
  end subroutine GPOCSR

subroutine residuCSR(A,b,x,t)
  integer,intent(in)::t !!taille des matrices
  real*8,dimension(t,t),intent(in)::A
  real*8,dimension(t),intent(in)::b
  real*8,dimension(t),intent(inout)::x
  real*8,dimension(t)::r,z
  real*8:: alpha,eps,nume,denom,max,norme
  integer :: k, kmax,i

  kmax=10000
  eps=0.00001
  nume=0.
  denom=0.
  alpha=0.
  z=0.
  Call multi_matCSR(r,AA,JA,IA,x,t)

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
     norme=abs(sum(r*r))
        if (norme<max) then
           max=norme
        end if

     k=k+1
     call write(k,sqrt(sum(r*r)),"ResMin.txt")
  end do
print*,"residu = ",max,k
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
