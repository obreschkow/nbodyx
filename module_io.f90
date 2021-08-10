module module_io

use shared_module_core
use shared_module_constants
use shared_module_parameters
use module_global

contains

subroutine load_parameters

   implicit none
   
   call read_parameters
   
   call get_parameter_value(para%inputfile,'inputfile')
   call get_parameter_value(para%outputpath,'outputpath','')
   call get_parameter_value(para%outputformat,'outputformat',1,min=1,max=3)
   call get_parameter_value(para%tinitial,'tinitial')
   call get_parameter_value(para%tfinal,'tfinal',min=para%tinitial)
   call get_parameter_value(para%dtmax,'dtmax',min=0.0_8)
   call get_parameter_value(para%dtmin,'dtmin',min=0.0_8,max=para%dtmax)
   call get_parameter_value(para%dtout,'dtout',min=0.0_8)
   call get_parameter_value(para%G,'G',min=0.0_8)
   call get_parameter_value(para%smoothing_radius,'smoothing_radius',min=0.0_8)
   call get_parameter_value(para%eta,'eta',min=0.0_8)
   call get_parameter_value(para%integrator,'integrator','leapfrog',min=5,max=8)
   call require_no_parameters_left
   
   ! turn integrator parameter into lower case and remove tabs
   para%integrator = remove_tabs(lowercase(para%integrator))
   
   ! handle output path
   para%outputpath = remove_tabs(para%outputpath)
   para%outputpath = dir(para%outputpath,ispath=.true.) ! fix path
   call make_path(para%outputpath)
   call system('rm -f '//trim(para%outputpath)//'snapshot_*')
   
   ! check input file
   para%inputfile = remove_tabs(para%inputfile)
   call check_file(para%inputfile,'r')
   
end subroutine load_parameters

subroutine load_particles

   implicit none
   integer              :: status,i,d
   real                 :: empty(7)
   
   ! determine number of particles
   open(1,file=trim(para%inputfile),action='read',form='formatted')
   n = 0
   status = 0
   do while (status==0)
      read(1,*,IOSTAT=status) empty
      if (status.ne.0) exit
      n = n+1
   end do
   close(1)
   call out('Number of particles: ',n)
   
   ! allocate memory
   if (allocated(p)) deallocate(p)
   allocate(p(n))
   
   ! load particles
   open(1,file=trim(para%inputfile),action='read',form='formatted')
   do i = 1,n
      read(1,*) p(i)%m,p(i)%x,p(i)%v
   end do
   close(1)
   
   ! pre-process
   do d = 1,3
      p%v(d) = p%v(d)
   end do
   p%acceleratable = p%m>=0.0
   p%m = abs(p%m)
   ntest = count(p%m==0.0)
   nbackground = n-count(p%acceleratable)
   
end subroutine load_particles

subroutine save_snapshot(first,last)

   implicit none
   logical,optional,intent(in)   :: first,last
   logical                       :: isfirst,islast
   integer                       :: i
   character(255)                :: sntxt,fn
   integer,parameter             :: id = 1 ! file id
   
   select case(para%outputformat)
      
      ! ascii format
      case (1)
      write(sntxt,'(A,I0.5,A)') 'snapshot_',snapshot,'.txt'
      fn = dir(para%outputpath,sntxt)
      open(id,file=trim(fn),action='write',form='formatted',status='replace')
      write(id,*) n
      write(id,*) t
      do i = id,n
         write(id,*) p(i)%x,p(i)%v
      end do
      close(id)
      
      ! binary format
      case (2)
      write(sntxt,'(A,I0.5,A)') 'snapshot_',snapshot,'.bin'
      fn = dir(para%outputpath,sntxt)
      open(id,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
      write(id) int(n,8),t,p%x(1),p%x(2),p%x(3),p%v(1),p%v(2),p%v(3)
      close(id)
      
      ! single binary stream
      case (3)
      if (present(first)) then
         isfirst = first
      else
         isfirst = .false.
      end if
      if (present(last)) then
         islast = last
      else
         islast = .false.
      end if
      if (isfirst) then
         fn = dir(para%outputpath,'snapshot_all.bin')
         open(id,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
         write(id) int(n,8)
      end if
      write(id) real(t,8),real(p%x(1),8),real(p%x(2),8),real(p%x(3),8),real(p%v(1),8),real(p%v(2),8),real(p%v(3),8)
      if (islast) close(id)
      
      ! default case
      case default
      call error('unknown output format')
      
   end select     
   
end subroutine save_snapshot

subroutine save_statistics

   implicit none
   
   open(1,file=dir(para%outputpath,'statistics.txt'),action='write',form='formatted',status='replace')
   write(1,('(A,I0)')) 'Number_of_particles                ',n
   write(1,('(A,I0)')) 'Number_of_test_particles           ',ntest
   write(1,('(A,I0)')) 'Number_of_background_particles     ',nbackground
   write(1,('(A,I0)')) 'Number_of_snapshots                ',snapshot+1
   write(1,('(A,I0)')) 'Number_of_interations              ',niterations
   write(1,('(A,I0)')) 'Number_of_acceleration_evaluations ',naccelerationevaluations
   close(1)
   
end subroutine save_statistics

end module module_io