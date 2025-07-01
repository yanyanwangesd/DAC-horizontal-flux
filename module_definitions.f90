MODULE definitions

   type parm
      integer n, istep,freq,num_restart
      double precision deltat,time,tfinal,h,ka,xc,tanthetac,hmn,m,tstart
      double precision rainfall_available,rainfall_minimum,rainfall_height
      double precision lmax,diffusivity,amin,min_erosion_rate,uplift_scalar1
      double precision uplift_scalar2,k_scalar1,ldivmax, max_adv, min_tan_head_slope, advection_scalarx, advection_scalary
      logical plot_triangles,plot_receiver,plot_donors,plot_no_receiver
      logical plot_no_donors,write_discharge,write_stack_order,plot_catchment
      logical plot_height,plot_precipitation
      logical move_points,capture,divide,small_divide,diffusion,jitter_nodes
      logical add_nodes
      logical read_restart, ascii, transient_divide
      logical show_vtkfine
      logical do_biology, do_elevation, do_slope, do_aspect, do_discharge, only_do_habitable, do_patches, do_split_merge, do_speciation, do_evolution, do_patch_size, do_topo_var, do_tree_ascii
      integer number_init_species, max_nspecies, max_npops, org_type
      double precision, dimension(:), pointer:: divergence_threshold, dispersal_rate, discharge_lim
      double precision, dimension(:), pointer::elev_lims, slope_lims, aspect_lims
      double precision aquatics_amin, evolution_rate, fitness_rate, min_patchsz, draw_rate, dispersal_mobility
      integer num_bins, sample_per_bin
      logical f_varies_with_xyz
      integer f_num_sets        ! number of condition blocks
      integer,dimension(:),pointer:: f_depends_on, f_variable_determined,f_polyseg, f_superpose     ! dimension f_num_sets
      double precision,dimension(:,:),pointer:: f_timebound       ! limit validity of condition block in time
      double precision,dimension(:,:,:),pointer::poly       ! polynomial coefficients
      integer,dimension(:,:),pointer:: pn   ! degree/elements per line
      !@YWANG, JAN. 2022
      double precision jitter_range  ! define how far node jitters laterally, should be between 0 and 1, recommend value 0.001
      double precision jitter_diffusivity ! define the random walk diffusivity, recommend value range 0.0001-0.01
      !@YWANG, JAN. 2022 
      
   end type parm

   type geom
      double precision xl,yl,zl
      double precision,dimension(:),pointer::x,y,z,u,v,w
      double precision,dimension(:),pointer::xdiv,ydiv,zdiv
      double precision,dimension(:),pointer::surface,discharge,precipitation
      double precision,dimension(:),pointer::k,erosion_rate,sediment_flux, chi, soil, max_divide_elev
      integer,dimension(:),pointer::strahler
      double precision,dimension(:,:),pointer::surface_share
      integer,dimension(:),pointer::fix,nb,catchment,card_flow,jitter_direction
      integer,dimension(:,:),pointer::nn
      integer,dimension(:,:),pointer::nndivtri !for each divide two triangles
      integer,dimension(:,:),pointer::nndivnode !for each divide two nodes
      integer,dimension(:),pointer::nt,nst ! number of triangles for each node, number of small triangles for each node (sub-grid structure for bio)
      integer, dimension(:,:),pointer::nnodetri ! for each node its neighboring triangles
      double precision,dimension(:,:),pointer::erosion_rate_history
      ! double precision,dimension(:,:),pointer::divide_elevation, dist_to_divide
      double precision,dimension(:,:),pointer::chanheadx,chanheady,chanheadz !(ndivide, 1:2)
      double precision,dimension(:,:,:),pointer::centroid_divides,chanhead_coords,divide_coords !for subgrid structure (node, neighbor, coordinate)
      double precision, dimension(:),pointer:: aw_elevation, aw_slope, aw_aspect, cell_area !area weighted means for grid nodes
      double precision,dimension(:,:),pointer:: triangle_area, triangle_slope, triangle_aspect, triangle_mean_elevation !for subgrid structure; multiple values per node (more than one per nb) for percent of total area, slope, and aspect
      ! double precision,dimension(:,:,:),pointer::tri_elev !for subgrid structure (node, small triangle, min/max)
      ! double precision, dimension(:,:),pointer::erosion_history_omega
      integer nx,ny,nnode,nnode_max,nnmax,ndivide,ndivide_max,ncapture
      integer boundary(4)
      integer,dimension(:,:),pointer::jitter_flag
      double precision, dimension(:),pointer:: jitter_residual
      
   end type geom

   type del
      integer ntriangles
      integer,dimension(:,:),pointer::icon,neighbours,numdivides,divides
      double precision,dimension(:,:),pointer::centers
   end type del

   type netw
      integer nnode,ndonmax,nlake
      integer,dimension(:),pointer::receiver,ndon,lakes,lakes_catch
      integer,dimension(:,:),pointer::donors
   end type netw

   type stck
      integer nnode
      integer,dimension(:),allocatable::order
   end type stck

   type timr
      sequence
      integer ntimer
      character*256 name(1024)
      real time_spent(1024)
      character*256 namein
      real timein
      real time_start,time_stop
   end type timr

   type lands
    	sequence
    	integer nfinetri, nfinenode
    	double precision, dimension(:),allocatable::xx,yy,zz,edot
    	integer, dimension(:,:), pointer::fineicon
   end type lands


end MODULE definitions
