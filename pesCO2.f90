module rep_ker
implicit none
real*8, parameter :: dk26f1 = 1.0d0/14.0d0, dk26f2 = 1.0d0/18.0d0, &
dk24f1 = 2.0d0/15.0d0, dk24f2 = 2.0d0/21.0d0, dk25f1 = 2.0d0/21.0d0,&
dk25f2 = 1.0d0/14.0d0, akf1=2.0d0/3.0d0

contains

function drker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker24, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker24 = dk24f1/xl**5 - dk24f2*xs/xl**6

end function drker24

function ddrker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker24, xl, xs

if (x .lt. xi) then
  ddrker24 = -dk24f2/xi**6
else
  ddrker24 = -5.0d0*dk24f1/x**6 + 6.0d0*dk24f2*xi/x**7
end if

end function ddrker24

function drker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker25, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker25 = dk25f1/xl**6 - dk25f2*xs/xl**7

end function drker25

function ddrker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker25, xl, xs

if (x .lt. xi) then
  ddrker25 = -dk25f2/xi**7
else
  ddrker25 = -6.0d0*dk25f1/x**7 + 7.0d0*dk25f2*xi/x**8
end if

end function ddrker25

function drker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker26, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker26 = dk26f1/xl**7 - dk26f2*xs/xl**8

end function drker26

function ddrker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker26, xl, xs

if (x .lt. xi) then
  ddrker26 = -dk26f2/xi**8
else
  ddrker26 = -7.0d0*dk26f1/x**8 + 8.0d0*dk26f2*xi/x**9
end if

end function ddrker26

function atker23(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: atker23, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

atker23 = 1.0d0 + xs*xl + 2.0d0*xs**2*xl  - akf1*xs**3

end function atker23

function datker23(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: datker23, xl, xs

if (x .lt. xi) then
  datker23 = xi + 4.0d0*x*xi  - 3.0d0*akf1*x**2
else
  datker23 = xi + 2.0d0*xi**2
end if

end function datker23


end module rep_ker

module surface1
implicit none
real*8, allocatable, dimension(:,:) :: asy_array1, asy_array2, & !asy_array3, &
darray1, darray2!, darray3
integer :: na1, na2, nda1, nda2

contains
subroutine pes3d(r12,r23,r31,totener,dvdr)
implicit none
real*8,parameter ::  m2 = 15.99491462d0, m3 = m2, m1 =12.00000000d0
real*8, intent(in) :: r12, r23, r31
real*8, dimension(:), intent(out) :: dvdr(3)
real*8, intent(out) :: totener
real*8 :: capr, theta
real*8, dimension(:,:) :: derv(3,3), deri(3,3), dw(3,3)
real*8, dimension(:) :: derj(3), w(3), ener(3), r(3), m(3), bigr(3)


!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! C---------O
! 1   r1    2

!=====================================================
!
!PES1               
!        3           
!        O           
!  r3   /            
!      /  r2         
!     /THETA            
!C___/_____O         
!1   r1    2         
!                    
!                 

    r(1)=r12
    r(2)=r23
    r(3)=r31
    m(1) = m1
    m(2) = m2
    m(3) = m3
    call crdtrf(r,m,capr,theta,derv)
    call calcener(capr,r12,theta, ener(1), derj, 1)
    ener(1)=ener(1)+0.08072566764673417d0

    deri(1,1) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
    deri(1,2) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
    deri(1,3) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

     dvdr(1) = deri(1,1)
     dvdr(2) = deri(1,2)
     dvdr(3) = deri(1,3)
    totener = ener(1)

end subroutine pes3d

subroutine crdtrf(r,m,capr,theta, derv)
implicit none
real*8, dimension(:), intent(in) :: r(3), m(3)
real*8, dimension(:,:), intent(out) :: derv(3,3)
real*8, intent(out) :: capr, theta
real*8, parameter :: pc = sqrt(epsilon(1.0d0))
real*8 :: cmu1, cmu2, rcm1, rcm2

!        3
!       /|\
!      / | \
!   r3/ R|  \r2
!    /   |th \
!   /____|____\ 
!  1    r1    2

cmu1 = m(2) / (m(1)+m(2))
cmu2 = m(1) / (m(1)+m(2))

rcm1 = r(1) * cmu1
rcm2 = r(1) * cmu2

capr = sqrt (r(2)**2/r(1)*rcm1 + r(3)**2/r(1)*rcm2 - rcm1*rcm2 )
if (abs(capr) < pc) capr = pc

theta = (rcm2**2+capr**2-r(2)**2)/2.0d0/rcm2/capr
theta=min(1.0d0,max(-1.0d0,theta))

theta = acos(theta)
if (theta==acos(-1.d0)) then
        !nprint*,"shite"
        theta=theta-1.d-5
end if 
derv(1,1) = -cmu1*cmu2*r(1)/capr   !dR/dr1
  
derv(1,2) = r(2)*cmu1/capr !dR/dr2

derv(1,3) = r(3)*cmu2/capr

derv(2,1) = 1.0d0
derv(2,2) = 0.0d0
derv(2,3) = 0.0d0

derv(3,1) = (derv(1,1)/capr*cos(theta)+cos(theta)/r(1)-(capr*derv(1,1)+rcm2*cmu2)/rcm2/capr)&
            /sqrt(1.0d0-cos(theta)**2)

derv(3,2) = (r(2)/rcm2/capr-derv(1,2)/rcm2+cos(theta)/capr*derv(1,2))/sqrt(1.0d0-cos(theta)**2)

derv(3,3) = (cos(theta)/capr*derv(1,3)-derv(1,3)/rcm2)/sqrt(1.0d0-cos(theta)**2)

return

end subroutine crdtrf


subroutine calcener(capr,smlr,theta, ener, der, sno)
use rep_ker
use RKHS            ! This module needs to be used by your code
!use param
implicit none
real*8 :: lambda
real*8, intent(out) :: ener
real*8, intent(in) :: capr, smlr, theta
integer, intent(in) :: sno
real*8, dimension(:), intent(out) :: der(3)
real*8,parameter :: pi = acos(-1.0d0), piby180 = pi/180.0d0
real*8 :: asener, anener, asder, z1, z2
real*8, dimension(:) :: ander(3), x(3)
integer :: kk, ii
type(kernel), save  :: pes11 !pes12!, pes13           ! The kernel type is needed to set up and evaluate a RKHS model
logical, save :: stored = .false., kread = .false.
logical, save :: ker1 = .false., ker2 = .false.!, ker3 = .false.
character (len=80), save :: datapath="./"

if (.not. ker1) then
  inquire(file=trim(datapath)//"pes11_ht.kernel", exist=ker1)   ! file_exists will be true if the file exists and false otherwise
end if


lambda=0.1d-19

if (.not. stored ) then

open(unit=1001,file=trim(datapath)//"asymp.coeff", status = "old")

read(1001,*)na1
allocate(asy_array1(na1,2))
do ii = 1, na1
  read(1001,*)asy_array1(ii,1), asy_array1(ii,2)
end do


read(1001,*)nda1
allocate(darray1(nda1,2))
do ii = 1, nda1
  read(1001,*)darray1(ii,1), darray1(ii,2)
end do

stored = .true.

end if

if (.not. kread) then
if (ker1) then
  call pes11%load_from_file(trim(datapath)//"pes11_ht.kernel")
  kread = .true.
else
  call pes11%read_grid(trim(datapath)//"pes11_ht.csv")
!  print*,"IAMHERE"
  call pes11%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
  call pes11%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
  call pes11%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL) 

  call pes11%calculate_coefficients_fast()
  call pes11%calculate_sums()
  call pes11%save_to_file(trim(datapath)//"pes11_ht.kernel")
  kread = .true.
end if
end if

x(1)=(1.0d0-cos(theta))/2.0d0
x(2)=capr
x(3)=smlr

if (sno==1) then
  asener = 0.0d0
  asder=0.0d0
  do kk = 1,na1
    asener = asener + drker26(smlr,asy_array1(kk,1))*asy_array1(kk,2)
    asder = asder + ddrker26(smlr,asy_array1(kk,1))*asy_array1(kk,2)
  end do

  anener = 0.0d0
  ander=0.0d0
  call pes11%evaluate_fast(x,anener,ander)

    ener = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/2.0d0
end if

return

end subroutine calcener

end module

module callpot1

contains

subroutine co21appes(tmpr, totener, dvdr)
use surface1
implicit none
real*8, dimension (:), intent (in) :: tmpr(3)
real*8, dimension (:), intent (out) :: dvdr(3)
real*8, intent (out) :: totener
real*8, dimension (:) :: dvdx(4),xp(3), tmpdvdr(3), r(3)
real*8, parameter :: dx = 0.005d0
real*8 :: d0h, d02h,ener
integer :: ii

!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! C---------O
! 1   r1    2

r(1) = tmpr(1)
r(2) = tmpr(2)
r(3) = tmpr(3)

call pes3d(r(1),r(2),r(3),totener,dvdr)
if (any(isNaN(dvdr))) then

  do ii = 1, 3

!==========================================================
!   Forward difference to calculate first derivative      =
!==========================================================

    xp=r
    xp(ii)=r(ii)-2.0d0*dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdx(1)=ener

    xp=r
    xp(ii)=r(ii)-dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdx(2)=ener

    xp=r
    xp(ii)=r(ii)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdx(3)=ener

    xp=r
    xp(ii)=r(ii)+2.0d0*dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdx(4)=ener

    d0h=(dvdx(3)-dvdx(2))/2.0d0/dx
    d02h=(dvdx(4)-dvdx(1))/4.0d0/dx
    dvdr(ii)=(4.0d0*d0h-d02h)/3.0d0
end do
end if
return

end subroutine co21appes

subroutine diatco(r,ener,der)
use surface1
use rep_ker
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: ener, der
integer :: kk

ener=0.0d0
der=0.0d0
do kk = 1, nda1
    ener = ener + drker26(r,darray1(kk,1))*darray1(kk,2)
    der = der + ddrker26(r,darray1(kk,1))*darray1(kk,2)
end do

return

end subroutine diatco

end module 
