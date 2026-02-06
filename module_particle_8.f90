module module_particle

   implicit none

   type type_particle
   
      real*8   :: m        ! mass
      real*8   :: x(3)     ! position
      real*8   :: v(3)     ! velocity
      real*8   :: a(3)     ! acceleration
      
   end type type_particle

end module module_particle