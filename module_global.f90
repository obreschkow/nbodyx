module module_global

   use module_particle

   implicit none

   type type_parameter

      character(255) :: parameterfile
      character(255) :: parameterset
      character(255) :: inputfile
      character(255) :: outputpath
      integer*4      :: outputformat
      integer*4      :: kind
      logical*4      :: include_bg
      logical*4      :: acceleration
      real*8         :: p(5)
      real*8         :: tinitial
      real*8         :: tfinal
      real*8         :: dtmax
      real*8         :: dtmin
      real*8         :: dtout
      real*8         :: G
      real*8         :: smoothing_radius
      real*8         :: eta
      real*8         :: box_size
      character(8)   :: integrator
   
   end type
   
   type(type_particle),allocatable  :: p(:)
      
   ! Simulation parameters read from the file 'parameters.txt'
   type(type_parameter)    :: para
   
   ! simulation parameters calculated internally
   integer*4            :: n           ! Total number of particles
   integer*4            :: ntest       ! Number of particles that cannot accelerate others (these have zero mass in the ICs)      
   integer*4            :: nbackground ! Number of particles that cannot be accelerated (these have negative mass in the ICs)
   integer,allocatable  :: set(:)      ! indices of acceleratable particles (with non-negative masses in the ICs)
   
   ! simulation time variables
   real*8               :: t           ! Physical simulation time
   integer*4            :: snapshot
   integer*4            :: niterations
   integer*4            :: naccelerationevaluations

end module module_global