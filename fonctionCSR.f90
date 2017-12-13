module fonctionCSR
use fonctions
use CSRconvert
contains

Subroutine JacobiCSR(AA,JA,IA,b,x,t)
    integer,intent(in)::t !!taille des matrices


    real*8,dimension(t),intent(in)::b
    real*8,dimension(:),intent(in)::AA
    integer,dimension(:),intent(in)::JA
    integer,dimension(:),intent(in)::IA
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::d,Xnext,r
    integer ::i,j,k,l,m, Nb_elem,Nb_li,kmax,max


    real*8::sigma,diag,epsn,norme




    call multi_matCSR(AA,IA,JA,x,r,t)
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
       call multi_matCSR(AA,IA,JA,x,r,t)
       r=b-r
       norme=abs(sum(r*r))
       if (norme<max) then
          max=norme
       end if
       k=k+1
    end do

  end subroutine JacobiCSR


subroutine GPOCSR(AA,JA,IA,b,x,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(t),intent(in)::b
    real*8,dimension(:),intent(in)::AA
    integer,dimension(:),intent(in)::JA
    integer,dimension(:),intent(in)::IA
    real*8,dimension(t),intent(inout)::x
    real*8,dimension(t)::r,z
    real*8:: alpha,eps,nume,denom,max,norme
    integer :: k, kmax,i



    call multi_matCSR(AA,JA,IA,x,r,t)

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

       call multi_matCSR(AA,JA,IA,r,z,t)


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


subroutine residuCSR(AA,JA,IA,b,x,t)
      integer,intent(in)::t !!taille des matrices
      real*8,dimension(:),intent(in)::AA
      integer,dimension(:),intent(in)::JA
  integer,dimension(:),intent(in)::IA
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
  Call multi_matCSR(AA,JA,IA,x,r,t)

  r=b-r
  max=0
  k=0
  max=abs(sum(r*r))

  do while (k<kmax .and.  max>eps)
     call multi_matCSR(AA,JA,IA,r,z,t)
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
  subroutine ArnoldiCSR(AA,JA,IA,r,H,vm,t)
    integer,intent(in)::t !!taille des matrices
    real*8,dimension(:),intent(in)::AA
    integer,dimension(:),intent(in)::JA
    integer,dimension(:),intent(in)::IA
    real*8,dimension(t,t)::v
    real*8,dimension(t,t),intent(out)::H
    real*8,dimension(t),intent(in)::r
    real*8,dimension(t), intent(out)::vm
    real*8,dimension(t)::z,q,z1
    integer :: i,j

    v(1,:)=r/sqrt(sum(r*r))   !! sqrt(sum(v*v)) revient à faire la norme :)

    do j=1,t

       do i=1,j
          call multi_matCSR(AA,JA,IA,v(j,:),Z1,t)
          h(i,j)=dot_product(z1,v(i,:))  !!fait le produit scalaire (à mettre de partout peut etre)
       end do
    q=0.
    do i=1,j
       q=q+h(i,j)*v(i,:)
    end do
    call multi_matCSR(AA,JA,IA,v(j,:),z,t)

    z=Z-q
    H(j+1,j)=sqrt(sum(z*z))
    if (H(j+1,j)==0) then
       stop
    end if
    v(j+1,:)=z/H(j+1,j)
    vm=v(j+1,:)
    end do

end subroutine ArnoldiCSR



subroutine multi_matCSR(AA,JA,IA,F,AF,t)
      integer,intent(in)::t
      real*8,dimension(t),intent(out)::AF
      real*8,dimension(t),intent(in)::F
      real*8,dimension(:),intent(in)::AA
      integer,dimension(:),intent(in)::IA
      integer,dimension(:),intent(in)::JA
      integer :: i,j,k,l, Nb_elem, Nb_li
      real*8::res

      do i=1,N
         res=0.d0
         k=IA(i)
         l=IA(i+1)-1
         do j=k,l
            res=res+AA(j)*F(JA(j))
         end do
         AF(i)=res
      end do

 end subroutine multi_matCSR




end module fonctionCSR
