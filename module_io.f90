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
   call get_parameter_value(para%kind,'kind',8,min=4,max=8)
   call get_parameter_value(para%include_bg,'include_bg',.true.)
   call get_parameter_value(para%tinitial,'tinitial')
   call get_parameter_value(para%tfinal,'tfinal',min=para%tinitial)
   call get_parameter_value(para%dtmax,'dtmax',min=0.0_8)
   call get_parameter_value(para%dtmin,'dtmin',min=0.0_8,max=para%dtmax)
   call get_parameter_value(para%dtout,'dtout',min=0.0_8)
   call get_parameter_value(para%G,'G',1.0_8,min=0.0_8)
   call get_parameter_value(para%smoothing_radius,'smoothing_radius',0.0_8,min=0.0_8)
   call get_parameter_value(para%eta,'eta',0.01_8,min=0.0_8)
   call get_parameter_value(para%box_size,'box_size',0.0_8,min=0.0_8)
   call get_parameter_value(para%integrator,'integrator','leapfrog',min=3,max=8)
   call get_parameter_value(para%acceleration,'acceleration',.true.)
   call get_parameter_value(para%p(1),'p1',0.0_8)
   call get_parameter_value(para%p(2),'p2',0.0_8)
   call get_parameter_value(para%p(3),'p3',0.0_8)
   call get_parameter_value(para%p(4),'p4',0.0_8)
   call get_parameter_value(para%p(5),'p5',0.0_8)
   call require_no_parameters_left
   
   ! additional argument checks
   if ((para%kind.ne.4).and.(para%kind.ne.8)) then
      call error('parameter kind must be 4 or 8')
   end if
   
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
   integer              :: io,i,d,j
   character(128)       :: line
   
   ! determine number of particles
   open(1,file=trim(para%inputfile),action='read',form='formatted')
   n = 0
   do
      read(1,'(A)',IOSTAT=io) line
      if (io/=0) exit
      if (.not.(isempty(line).or.(line(1:1)=='#'))) n = n+1
   end do
   close(1)
   call out('Total number of particles: ',n)
   
   ! allocate memory
   if (allocated(p)) deallocate(p)
   allocate(p(n))
   
   ! load particles
   open(1,file=trim(para%inputfile),action='read',form='formatted')
   i = 0
   do
      read(1,'(A)',IOSTAT=io) line
      if (io/=0) exit
      if (.not.(isempty(line).or.(line(1:1)=='#'))) then
         i = i+1
         read(line,*) p(i)%m,p(i)%x,p(i)%v
      end if
   end do
   close(1)
   
   ! reset accelerations (do not try to change this loop to avoid compiler optimisation issues)
   do i = 1,n
      p(i)%a = 0.0
   end do
   
   ! count test and background particles
   ntest = count(p%m==0.0)
   nbackground = count(p%m<0.0)
   call out('Number of test particles (m=0): ',ntest)
   call out('Number of background particles: ',nbackground)
   
   ! make indices of acceleratable particles (i.e. = normal particles that are accelerated, as opposed to background particles)
   if (allocated(set)) deallocate(set)
   allocate(set(n-nbackground))
   j = 0
   do i = 1,n
      if (p(i)%m>=0.0) then
         j = j+1
         set(j) = i
      end if
   end do
   if (j.ne.n-nbackground) call error('indexing error')
   p%m = abs(p%m)
   
end subroutine load_particles

subroutine save_snapshot(first,last)

   implicit none
   logical,optional,intent(in)   :: first,last
   logical                       :: isfirst,islast
   integer                       :: i,k
   character(255)                :: sntxt,fn
   integer,parameter             :: id = 1 ! file id
   
   select case(para%outputformat)
      
      ! ascii format
      case (1)
      write(sntxt,'(A,I0.6,A)') 'snapshot_',snapshot,'.txt'
      fn = dir(para%outputpath,sntxt)
      open(id,file=trim(fn),action='write',form='formatted',status='replace')
      if (para%include_bg.or.(nbackground==0)) then
         if (para%kind==4) then
            write(id,*) int(n,4)
            write(id,*) real(t,4)
            do i = 1,n
               write(id,*) real(p(i)%x,4),real(p(i)%v,4)
            end do
         else
            write(id,*) int(n,8)
            write(id,*) real(t,8)
            do i = 1,n
               write(id,*) real(p(i)%x,8),real(p(i)%v,8)
            end do
         end if
      else
         if (para%kind==4) then
            write(id,*) int(n-nbackground,4)
            write(id,*) real(t,4)
            do k = 1,n-nbackground
               i = set(k)
               write(id,*) real(p(i)%x,4),real(p(i)%v,4)
            end do
         else
            write(id,*) int(n-nbackground,8)
            write(id,*) real(t,8)
            do k = 1,n-nbackground
               i = set(k)
               write(id,*) real(p(i)%x,8),real(p(i)%v,8)
            end do
         end if
      end if
      close(id)
      
      ! binary format
      case (2)
      write(sntxt,'(A,I0.6,A)') 'snapshot_',snapshot,'.bin'
      fn = dir(para%outputpath,sntxt)
      open(id,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
      if (para%include_bg.or.(nbackground==0)) then
         if (para%kind==4) then
            write(id) int(n,4),real(t,4),real(p%x(1),4),real(p%x(2),4),real(p%x(3),4), &
            & real(p%v(1),4),real(p%v(2),4),real(p%v(3),4)
         else
            write(id) int(n,8),real(t,8),real(p%x(1),8),real(p%x(2),8),real(p%x(3),8), &
            & real(p%v(1),8),real(p%v(2),8),real(p%v(3),8)
         end if
      else
         if (para%kind==4) then
            write(id) int(n,4),real(t,4),real(p(set)%x(1),4),real(p(set)%x(2),4),real(p(set)%x(3),4), &
            & real(p(set)%v(1),4),real(p(set)%v(2),4),real(p(set)%v(3),4)
         else
            write(id) int(n,8),real(t,8),real(p(set)%x(1),8),real(p(set)%x(2),8),real(p(set)%x(3),8), &
            & real(p(set)%v(1),8),real(p(set)%v(2),8),real(p(set)%v(3),8)
         end if
      end if
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
         if (para%include_bg.or.(nbackground==0)) then
            if (para%kind==4) then
               write(id) int(n,4)
            else
               write(id) int(n,8)
            end if
         else
            if (para%kind==4) then
               write(id) int(n-nbackground,4)
            else
               write(id) int(n-nbackground,8)
            end if
         end if
      end if
      if (para%include_bg.or.(nbackground==0)) then
         if (para%kind==4) then
            write(id) real(t,4),real(p%x(1),4),real(p%x(2),4),real(p%x(3),4), &
            & real(p%v(1),4),real(p%v(2),4),real(p%v(3),4)
         else
            write(id) real(t,8),real(p%x(1),8),real(p%x(2),8),real(p%x(3),8), &
            & real(p%v(1),8),real(p%v(2),8),real(p%v(3),8)
         end if
      else
         if (para%kind==4) then
            write(id) real(t,4),real(p(set)%x(1),4),real(p(set)%x(2),4),real(p(set)%x(3),4), &
            & real(p(set)%v(1),4),real(p(set)%v(2),4),real(p(set)%v(3),4)
         else
            write(id) real(t,8),real(p(set)%x(1),8),real(p(set)%x(2),8),real(p(set)%x(3),8), &
            & real(p(set)%v(1),8),real(p(set)%v(2),8),real(p(set)%v(3),8)
         end if
      end if
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