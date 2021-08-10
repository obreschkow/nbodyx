module module_global

   ! subroutine define_global_variables

   implicit none

   type type_parameter

      character(255) :: parameterfile
      character(255) :: parameterset
      character(255) :: inputfile
      character(255) :: outputpath
      integer*4      :: outputformat
      real*8         :: tinitial
      real*8         :: tfinal
      real*8         :: dtmax
      real*8         :: dtmin
      real*8         :: dtout
      real*8         :: G
      real*8         :: smoothing_radius
      real*8         :: eta
      character(8)   :: integrator
   
   end type
   
   type type_particle
   
      ! particle type
      logical  :: acceleratable  ! if .false., the particle is considered part of an external background, and not accelerated
            
      ! dynamic variables
      real*8   :: m        ! mass
      real*16   :: x(3)     ! position ! for increased accuracy, it normally suffices to change the type of x to real*16
      real*8   :: v(3)     ! velocity
      real*8   :: a(3)     ! acceleration
      
   end type type_particle
   
   type(type_particle),allocatable  :: p(:)
      
   ! Simulation parameters read from the file 'parameters.txt'
   type(type_parameter)    :: para
   
   ! simulation parameters calculated internally
   integer*4            :: n           ! Total number of particles
   integer*4            :: ntest       ! Number of particles that cannot accelerate others (these have zero mass in the ICs)      
   integer*4            :: nbackground ! Number of particles that cannot be accelerated (these have negative mass in the ICs)
   
   ! simulation time variables
   real*8               :: t           ! Physical simulation time
   integer*4            :: snapshot
   integer*4            :: niterations
   integer*4            :: naccelerationevaluations

end module module_global