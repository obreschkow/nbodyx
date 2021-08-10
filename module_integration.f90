module module_integration

   use module_global
   use module_io
   
   private
   public   :: run_simulation
   
   abstract interface
      subroutine integrator_abstract(dt)
         real*8,intent(in) :: dt
      end subroutine
   end interface
   
   integer,allocatable  :: set(:)   ! indices of acceleratable particles
   real*8               :: dtvar    ! suggested variable time-step
   procedure(integrator_abstract), pointer   :: integrator => NULL()
   
   contains
   
   subroutine run_simulation
   
      implicit none
      integer  :: i,j
      real*8   :: dt,t_next_snapshot
      
      ! choose custom integrator
      select case (trim(para%integrator))
         case ('euler')
            integrator => integrator_euler
         case ('leapfrog')
            integrator => integrator_leapfrog
         case ('yoshida')
            integrator => integrator_yoshida4
         case ('yoshida4')
            integrator => integrator_yoshida4
         case ('yoshida6')
            integrator => integrator_yoshida6
         case default
            call error('unknown integrator: '//trim(para%integrator))
      end select              
      
      ! make indices of acceleratable particles (i.e. = normal particles that are accelerated, as opposed to background particles)
      if (allocated(set)) deallocate(set)
      allocate(set(n-nbackground))
      j = 0
      do i = 1,n
         if (p(i)%acceleratable) then
            j = j+1
            set(j) = i
         end if
      end do
      if (j.ne.n-nbackground) call error('indexing error')
      
      ! initialise
      t = para%tinitial
      snapshot = 0
      niterations = 0
      naccelerationevaluations = 0
      call save_snapshot(first=.true.)
      t_next_snapshot = para%tinitial+para%dtout
      if (trim(para%integrator).ne.'euler') call evaluate_accelerations
      
      ! iterate
      do while (t<para%tfinal)
      
         ! make time step  
         dt = minval((/para%dtmax, para%tfinal-t, t_next_snapshot-t, max(para%dtmin, dtvar)/))
      
         ! advance all variables by dt
         call integrator(dt)
         t = t+dt
      
         ! save snapshot, if needed
         if ((t>=t_next_snapshot).or.(t>=para%tfinal))   then
            snapshot = snapshot+1
            call save_snapshot(last=t>=para%tfinal)
            t_next_snapshot = para%tinitial+(snapshot+1)*para%dtout
            if (abs(t_next_snapshot-para%tfinal)/para%dtout<1e-10) t_next_snapshot = para%tfinal
         end if
         
         ! count number of iterations
         niterations = niterations+1
      
      end do
      
      ! save simulation statistics
      call save_statistics
   
   end subroutine run_simulation
   
   subroutine integrator_euler(dt)
   
      implicit none
      real*8,intent(in) :: dt
      integer           :: d
   
      call evaluate_accelerations
      do d = 1,3
         p(set)%x(d) = p(set)%x(d)+p(set)%v(d)*dt+0.5*p(set)%a(d)*dt**2
         p(set)%v(d) = p(set)%v(d)+p(set)%a(d)*dt
      end do
   
   end subroutine integrator_euler
   
   subroutine integrator_leapfrog(dt)
   
      implicit none
      real*8,intent(in) :: dt
      integer           :: d
      real*8            :: ht
   
      ht = 0.5*dt
      do d = 1,3
         p(set)%v(d) = p(set)%v(d)+p(set)%a(d)*ht
         p(set)%x(d) = p(set)%x(d)+p(set)%v(d)*dt
      end do
      call evaluate_accelerations
      do d = 1,3
         p(set)%v(d) = p(set)%v(d)+p(set)%a(d)*ht
      end do
   
   end subroutine integrator_leapfrog
   
   subroutine integrator_yoshida4(dt)
   
      implicit none
      real*8,intent(in)    :: dt
      integer              :: dim,i
      
      real*8,parameter     :: w1 = 1.0_8/(2.0_8-2.0_8**(1.0_8/3.0_8))
      real*8,parameter     :: w0 = -2.0_8**(1.0_8/3.0_8)*w1
      real*8,parameter     :: c(4) = (/w1/2,(w0+w1)/2,(w0+w1)/2,w1/2/)
      real*8,parameter     :: d(3) = (/w1,w0,w1/)
      real*8,allocatable   :: dx(:,:)
      
      allocate(dx(size(set),3))
      dx = 0 ! the use of dx allows to increase the numerical accuracy significantly
      
      do i = 1,3
         do dim = 1,3
            dx(:,dim) = dx(:,dim)+c(i)*p(set)%v(dim)*dt
         end do
         call evaluate_accelerations(dx)
         do dim = 1,3
            p(set)%v(dim) = p(set)%v(dim)+d(i)*p(set)%a(dim)*dt
         end do
      end do
      do dim = 1,3
         dx(:,dim) = dx(:,dim)+c(4)*p(set)%v(dim)*dt
         p(set)%x(dim) = p(set)%x(dim)+dx(:,dim)
      end do
      
   end subroutine integrator_yoshida4
   
   subroutine integrator_yoshida6(dt)
   
      implicit none
      real*8,intent(in) :: dt
      integer           :: dim,i
      
      real*8,parameter     :: x1 = 1.0_8/(2.0_8-2.0_8**(1.0_8/3.0_8))
      real*8,parameter     :: x0 = -2.0_8**(1.0_8/3.0_8)*x1
      real*8,parameter     :: y1 = 1.0_8/(2.0_8-2.0_8**(1.0_8/5.0_8))
      real*8,parameter     :: y0 = -2.0_8**(1.0_8/5.0_8)*y1
      real*8,parameter     :: d(9) = (/x1*y1,x0*y1,x1*y1,x1*y0,x0*y0,x1*y0,x1*y1,x0*y1,x1*y1/)
      real*8,parameter     :: c(10) = (/d(1)/2.0_8,(d(1)+d(2))/2.0_8,(d(2)+d(3))/2.0_8, &
                              & (d(3)+d(4))/2.0_8,(d(4)+d(5))/2.0_8,(d(5)+d(6))/2.0_8,&
                              & (d(6)+d(7))/2.0_8,(d(7)+d(8))/2.0_8,(d(8)+d(9))/2.0_8,d(9)/2.0_8/)
      real*8,allocatable   :: dx(:,:)
      
      allocate(dx(size(set),3))
      dx = 0 ! the use of dx allows to increase the numerical accuracy significantly
      
      do i = 1,9
         do dim = 1,3
            dx(:,dim) = dx(:,dim)+c(i)*p(set)%v(dim)*dt
         end do
         call evaluate_accelerations(dx)
         do dim = 1,3
            p(set)%v(dim) = p(set)%v(dim)+d(i)*p(set)%a(dim)*dt
         end do
      end do
      do dim = 1,3
         dx(:,dim) = dx(:,dim)+c(10)*p(set)%v(dim)*dt
         p(set)%x(dim) = p(set)%x(dim)+dx(:,dim)
      end do
      
   end subroutine integrator_yoshida6
      
   subroutine evaluate_accelerations(dx)
   
      ! evaluates the accelerations of all non-fixed particles (i.e. those with mass>=0 in the ICs file) due to
      ! all particles with non-zero mass
      
      implicit none
      real*8,intent(in),optional :: dx(:,:)
      integer  :: i,j,k
      real*8   :: dist(3),rsqr,z,asqr,rsmoothsqr,f(3)
      real*8,allocatable   :: minrsqr(:)
      
      ! reset accelerations (could be extended to reset to an external field)
      do k = 1,n-nbackground
         i = set(k)
         p(i)%a = 0.0
      end do
      
      ! compute accelerations between pairs
      rsmoothsqr = para%smoothing_radius**2
      allocate(minrsqr(n))
      minrsqr = huge(minrsqr)
      
      if ((ntest==0).and.(nbackground==0)) then
      
         ! deal with the case where all particles are massive and acceleratable
      
         do i = 1,n-1 ! acceleration of particle i
            do j = i+1,n ! due to particle j
               dist = real(p(j)%x-p(i)%x,8)
               if (present(dx)) dist = dist+dx(j,:)-dx(i,:)
               rsqr = max(rsmoothsqr,sum(dist**2))
               minrsqr(i) = min(minrsqr(i),rsqr)
               minrsqr(j) = min(minrsqr(j),rsqr)
               f = dist/rsqr**1.5_8
               p(i)%a = p(i)%a+p(j)%m*f
               p(j)%a = p(j)%a-p(i)%m*f
            end do
         end do
         do i = 1,n
            p(i)%a = p(i)%a*para%G
         end do
      
      else
      
         ! deal with the case where some particles are massless and/or non-acceleratable
      
         do k = 1,n-nbackground ! acceleration of particle i
            i = set(k)
            do j = 1,n ! due to particle j
               if (j.ne.i) then
                  if (p(j)%m>0.0) then
                     dist = real(p(j)%x-p(i)%x,8)
                     if (present(dx)) dist = dist+dx(j,:)-dx(i,:)
                     rsqr = max(rsmoothsqr,sum(dist**2))
                     minrsqr(i) = min(minrsqr(i),rsqr)
                     p(i)%a = p(i)%a+p(j)%m*dist/rsqr**1.5_8
                  end if
               end if
            end do
            p(i)%a = para%G*p(i)%a
         end do
         
      end if
      
      ! new variable time step
      z = huge(z)
      do i = 1,n
         asqr = sum(p(i)%a**2)+1e-50_8
         z = min(z,minrsqr(i)/asqr)
      end do
      dtvar = para%eta*z**0.25
      
      ! count number of acceleration evaluations
      naccelerationevaluations = naccelerationevaluations+1
   
   end subroutine evaluate_accelerations
   
end module module_integration