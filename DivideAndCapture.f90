program DivideAndCapture

  ! Beta version of new cascade algorithm that incorporates river capture
  ! but not the computation of precise divide locations
  ! Jean BRAUN 20 july 2009

  use definitions

  implicit none

  ! for a definition of the derived types see the file module_definitions.f90

  type (geom) geometry
  type (del) delaunay
  !type (del) rddelaunay
  type (netw) network
  type (stck) stack
  type (timr) timer
  type (parm) params
  type (lands)  landscape
  character*11  title
  integer old_ncapture, i, j, k, m, l, n
  double precision li,lj



  integer date_time(8),ncount
 character*10 b(3)

  ! double precision, dimension(:), allocatable :: oldz
  ! sets a common for timing of various subroutines
  ! the results are printed to the screen at the end of the run

  common /global_timer/ timer

  call date_and_time(b(1), b(2), b(3), date_time)
  print *,'wall clock=',date_time(5),':',date_time(6)

  call start_timing
  call floader( params, geometry)
  call initialize_parameters (params) ! initializes the model parameters
  if (params%read_restart) then
     call reinit_geometry(geometry,params,network,delaunay)
  else
     call initialize_geometry (geometry,params) ! initializes the model geometry and history
  endif
  call set_params_as_function ( geometry, params )

  if (.not.params%read_restart) then
     delaunay%ntriangles=0
     call calculate_delaunay (geometry,delaunay) ! calculates Delaunay and neighbour (triangle) information
  endif

  if (params%read_restart) then
     call find_network (geometry,delaunay,network,1) ! find neighbour list
  else
     network%nnode=0
     call find_network (geometry,delaunay,network,0) ! find neighbour list and donor/receiver list
  endif

  if (.not.params%read_restart) then
     call find_precipitation (geometry,params) ! orography (not fully working yet)
     call find_polygon_surface (geometry,network,params,delaunay) ! find surface area attached to each point
  endif

  call reboundary(geometry, params)

  stack%nnode=0
  call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms

  call find_strahler(geometry,network,stack)

  call find_catchment (geometry,network,stack,delaunay,params) ! find catchments

  call output_z_tau(params,geometry,network,stack)



    if (params%read_restart) then !need to call all these to correctly calculate divides and subgrid structure
      call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
      call find_network (geometry,delaunay,network,1) ! find neighbour list only
      call captures_and_divides (geometry,network,params,delaunay) ! called before time loop in DAC-BIO to find divide elevations
      call find_polygon_surface (geometry,network,params,delaunay) ! find surface area attached to each point
      if (params%add_nodes) deallocate(stack%order)
      call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms
      ! call find_catchment (geometry,network,stack,delaunay,params) ! find catchments
      ! call find_discharge (geometry,network,stack) ! compute discharge (improved cascade algorithm)
    end if

  call VTK (geometry,delaunay,params,network,stack)


  if (params%show_vtkfine)  call finetri (geometry, params,delaunay, network, landscape)
  if (params%show_vtkfine) call VTKfine (geometry,delaunay,params,network,stack,landscape)

  if (params%ascii) then
    call write_ascii(geometry,params,network)
  endif

  call output_geometry(geometry,params,network,delaunay)

  !more frequent global data on the system is output to these files
  if (params%ascii) then
      open (80,file='ASCII/capture_data',status='unknown')
      write(80,'(F8.1,i10)') params%time, geometry%ncapture
      old_ncapture = geometry%ncapture
  endif

  open (101,file='ASCII/deltaerosion',status='unknown')
  write(101,'(e15.7, e15.7, i10)') 0.d0, 0.d0, geometry%nnode
  open (102,file='ASCII/deltachi',status='unknown')
  write(102,'(e15.7, e15.7, i10)') 0.d0, 0.d0, geometry%nnode

  do while (params%time.lt.params%tfinal) ! start of time loop

     params%deltat=min(params%deltat,params%tfinal-params%time)
     params%time=params%time+params%deltat
     params%istep=params%istep+1
     
     call set_params_as_function ( geometry, params )

     call uplift_and_advect (geometry,params) ! move and uplift the grid according to the prescribed velocity field

     call cardinal_flow_direction (geometry,params,network) !calculate what boundary each point flows to (used in bio routines and lateral jitter)

     if (params%jitter_nodes) call lateral_jitter (geometry,network,params,delaunay) !jitter nodes

     if (params%move_points) then !with horizontal advection
       call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
       call find_network (geometry,delaunay,network,1) ! find neighbour list only
     endif

     if (params%add_nodes) then
       call add_remove_nodes(geometry,network,params,delaunay,stack)
       call check_neighbor_receiver (geometry,network,ncount)
     endif

     if (params%transient_divide) call update_erosion_rate_history(geometry,params)

     call captures_and_divides (geometry,network,params,delaunay) ! find potential captures and adjusts network accordingly

     call find_polygon_surface (geometry,network,params,delaunay) ! find surface area attached to each point

     if (params%add_nodes) deallocate(stack%order)
     call find_order (geometry,network,stack,delaunay,params) ! find stack order to perform the erosion and cascade algorithms

     call find_catchment (geometry,network,stack,delaunay,params) ! find catchments
     
     call cardinal_flow_direction (geometry,params,network) !calculate again card_flow after captures,@YWANG, May 2022
     
     call find_discharge (geometry,network,stack) ! compute discharge (improved cascade algorithm)
     call erode (geometry,network,stack,delaunay,params) ! erode


     if ((params%istep/params%freq)*params%freq.eq.params%istep) then
         call find_strahler(geometry,network,stack)
         if (params%show_vtkfine) call finetri (geometry, params,delaunay,network, landscape)
         call output_z_tau(params,geometry,network,stack)
         !call VTK (geometry,delaunay,params,network,stack)
         if (params%show_vtkfine) call VTKfine (geometry,delaunay,params,network,stack,landscape)
         if (params%ascii) call write_ascii(geometry,params,network)
         call projection_area(geometry,delaunay,params,network,stack)
         call output_geometry(geometry,params,network,delaunay)
         print*, 'max number of neighbors is:', maxval(geometry%nb)
         print*, 'nnode = ', geometry%nnode
         print*, 'time (Myr) =', params%time/1.0d6
     endif


     !frequent data is output here
     if (params%ascii) then
         if (old_ncapture.ne.geometry%ncapture.or.MOD(params%istep,100).eq.0) then
            write(80,'(e15.7,i10)') params%time, geometry%ncapture
            old_ncapture = geometry%ncapture
         endif
     endif


  enddo ! end of time loop

  if (params%ascii) then
    write(80,'(e15.7,i10)') params%time, geometry%ncapture
    close (80)
  endif

  close (101)
  close (102)

  call find_strahler(geometry,network,stack)
  if (params%show_vtkfine)  call finetri (geometry, params,delaunay,network, landscape)
  call output_z_tau(params,geometry,network,stack)
  call VTK (geometry,delaunay,params,network,stack)
  if (params%show_vtkfine) call VTKfine (geometry,delaunay,params,network,stack,landscape)
  if (params%ascii) call write_ascii(geometry,params,network)
  call output_geometry(geometry,params,network,delaunay)
  call find_hack(geometry,network,stack)
  call show_timing

!   call clean_dacmem(geometry, params,delaunay,network,stack, landscape)

  call date_and_time(b(1), b(2), b(3), date_time)
  print*, 'wall clock=',date_time(5),':',date_time(6)

end program DivideAndCapture
