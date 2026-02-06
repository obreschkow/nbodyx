module module_acceleration

   use shared_module_core
   use module_global

   private
   public   :: external_acceleration
   
   contains
   
   subroutine external_acceleration(x,t,a)
   
      implicit none
      real*8,intent(in)    :: x(3)  ! particle position
      real*8               :: t     ! time
      real*8,intent(out)   :: a(3)  ! particle acceleration
      
      ! *** user part ********************************************************
      
      ! Kepler potential
      ! Central mass = para%p(1)
      real*8   :: r
      real*8   :: ar
      call nil(t) ! to avoid compiler warnings, since t is not used here
      r = sqrt(sum(x**2))
      ar = -para%G*para%p(1)/r**2 ! radial acceleration
      a = x/r*ar
      
      ! Hernquist potential
      ! Halo mass = para%p(1)
      ! Halo characteristic radius = para%(2)
      ! real*8   :: r
      ! real*8   :: ar
      !call nil(t) ! to avoid compiler warnings, since t is not used here
      !r = sqrt(sum(x**2))
      !ar = -para%G*para%p(1)/para%p(2)**2/(1+r/para%p(2))**2 ! radial acceleration
      !a = x/r*ar
      
      ! *** end user part ****************************************************
   
   end subroutine external_acceleration
   
end module module_acceleration