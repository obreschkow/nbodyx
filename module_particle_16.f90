module module_particle

   implicit none

   type type_particle
   
      real*8   :: m        ! mass
      real*16  :: x(3)     ! position, using 16-byte representation for very high accuracy simulations
      real*8   :: v(3)     ! velocity
      real*8   :: a(3)     ! acceleration
      
   end type type_particle

end module module_particle