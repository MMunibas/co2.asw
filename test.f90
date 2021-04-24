program testpes
use surface1 
use callpot1

implicit none
real*8 :: roo, roc, rco, ener1
real*8, dimension(:) :: dvdr(3), r(3)

!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! C---------O
! 1   r1    2


roo=4.4d0 !in bohr
roc=2.3d0 !in bohr
rco=2.3d0 !in bohr

r(1)=rco
r(2)=roo
r(3)=roc

call co21appes(r, ener1, dvdr) 

write(*,*)"Energy = ", ener1, "Hartree"
write(*,*)"dV/dr_i (Hartree/bohr) = ", dvdr

end
