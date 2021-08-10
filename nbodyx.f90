program surfsuite

   use shared_module_core
   use shared_module_arguments
   use shared_module_parameters
   use module_global
   use module_io
   use module_integration

   implicit none
   
   ! start user interface
   call set_version('1.0')
   call handle_arguments(require_task=.false.)
   call start_output
   
   ! handle general options
   call get_option_value(para%parameterfile,'-parameterfile','parameters.txt')
   call set_parameterfile(para%parameterfile)
   call get_option_value(para%parameterset,'-parameterset','')
   call set_parameterset(para%parameterset)
   
   ! load parameters
   call load_parameters
   para%parameterset = parameterset ! this is in case of a default parameterset
   
   ! run tasks
   call load_particles
   call run_simulation
   
   ! finalize output on screen/logfile
   call stop_output
    
end program surfsuite