subroutine add_remove_nodes(geometry,network,params,delaunay,stack)
! This version check crossing rivers of non-neighbour receiver connections to stablize grid.
! The check crossing river algorithm is especially useful for complex horizontal velocity field
! where relative shortening exists in the domain and near boundaries. --YWANG, October2024

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  type (stck) stack

  double precision l,xx,xd,fact3,fact2,lidon,lirec,px,py
  integer i,j,k,h,sp,tri,ncount,c, ptr, missing_node, it
  integer nrem_small_area,nrem_close_to_boundary,nrem_short_channel
  integer nadd_on_channel,nadd_river_no_connection,nadd_boundary, nadd_between_channel
  integer nnold, icount,flag,don,rec
  integer, dimension(2):: trik
  double precision xmin,xmax,ymin,ymax,min_dist
  double precision xmin_internal,xmax_internal,ymin_internal,ymax_internal
  integer num_samples, num_full_bin,samle_in_last_bin
  double precision tau,ave_erosion, x1,y1,x2,y2,l1,l2
  double precision, dimension(2) :: xcentroid

  double precision vmax
  logical flagcrossi, flagcrossj
  integer ri, rj
  integer crossriver
  
  
  call time_in ('add_remove_nodes')
    
  xmin = 0
  xmax = 0
  ymin = 0
  ymax = 0
  xmin_internal = 0
  xmax_internal = 0
  ymin_internal = 0
  ymax_internal = 0
 
 if(params%move_points .or. 1.eq.1)then
    do i = 1, geometry%nnode
    if(geometry%fix(i).ge.1)then
        xmin = min(geometry%x(i), xmin)
        xmax = max(geometry%x(i), xmax)
        ymin = min(geometry%y(i), ymin)
        ymax = max(geometry%y(i), ymax)
    else
        xmin_internal = min(geometry%x(i), xmin_internal)
        xmax_internal = max(geometry%x(i), xmax_internal)
        ymin_internal = min(geometry%y(i), ymin_internal)
        ymax_internal = max(geometry%y(i), ymax_internal)        
    endif    
    enddo   
 endif

 print*,'xmin xmax ymin ymax is:', xmin, xmax, ymin, ymax
 if(ymax.lt.ymax_internal) call VTK_debug(geometry,delaunay,params,network,stack,'nodesout',params%istep) 
 if(ymax.lt.ymax_internal) stop
 min_dist =  abs(params%max_adv*params%deltat)

  nnold=geometry%nnode !original number of nodes
  nadd_between_channel = 0
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! First stage: take care of all the cases where a new node should added !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!---------------------------- (1) ---------------------------- !!!
!   This section is to add nodes on channels where the channel connect two nodes that are not neighbors.
! a precheck and treatment of crossing rivers
      crossriver=0
  do i=1,geometry%nnode ! loop over the nodes
      if(geometry%fix(i).eq.0.and.network%receiver(i).ne.0)then ! deal with an internal non-lake node
           flag=0
           do k=geometry%nb(i),1,-1 ! loop over the neighbors
              if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
           enddo
           if(flag.eq.0)then !receiver is not a neighbor
              j=network%receiver(i)
              call check_crossrivers(geometry,network,params,delaunay, stack , i, flagcrossi, flagcrossj, ri, rj)
              if(ri.ne.0) print*,'ri and rj are:', ri, rj, geometry%x(ri), geometry%y(ri)
              if(flagcrossi.or.flagcrossj) then 
                  if(ri.ne.0) print*, 'network%receiver(i) before reroute,', network%receiver(i), geometry%x(i), geometry%y(i)
                  call reroute_donor_in_cross_river(geometry,network,params,delaunay, stack , i, ri, rj)
                  if(ri.ne.0) print*, 'network%receiver(i) after reroute,', network%receiver(i), geometry%x(i), geometry%y(i)
                  crossriver=crossriver+1                  
              endif              
           endif          
      endif!end deal with an internal non-lake node      
  enddo!end loop nodes
  print*,'number of cross rivers are:', crossriver

 !need to continuously update donors and receivers 
  network%nnode=geometry%nnode
  network%donors=0
  network%ndon=0
  do j=1,network%nnode
      k=network%receiver(j)
      if (k.ne.0) then
          network%ndon(k)=network%ndon(k)+1
          network%donors(network%ndon(k),k)=j
      endif
  enddo
  !call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
  call find_network (geometry,delaunay,network,1) ! find neighbour list only
  
  nadd_river_no_connection=0
  missing_node=0
  call check_neighbor_receiver (geometry,network,ncount)
  if (ncount.gt.0) missing_node = 1
  if(ncount.gt.0)print*,'number of non-neighbour connnections is:', ncount
  it=1
  do while (missing_node.eq.1)
      nnold = geometry%nnode
      do i=1,nnold ! loop over the nodes
          if(geometry%fix(i).eq.0.and.network%receiver(i).ne.0)then ! deal with an internal non-lake node
               flag=0
               do k=geometry%nb(i),1,-1 ! loop over the neighbors
                  if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
               enddo
               if(flag.eq.0)then !receiver is not a neighbor 
                  j=network%receiver(i)
                  call check_crossrivers(geometry,network,params,delaunay, stack , i, flagcrossi, flagcrossj, ri, rj)                   
                  if(flagcrossi.or.flagcrossj) then ! crossing river
                    crossriver=crossriver+1
                    call reroute_donor_in_cross_river(geometry,network,params,delaunay, stack , i, ri, rj)  
                  endif
                  if(ri.eq.0 .and. rj.eq.0)then ! not crossing river
                      ! assign coordinates to the newly added node
                      l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length between i and j
                     if(l.le.1.d1)then
                        network%receiver(i)=0 !make i a lake
                     else
                      nadd_river_no_connection=nadd_river_no_connection+1                
                      geometry%nnode=geometry%nnode+1
                      if (geometry%nnode.gt.geometry%nnode_max) then
                        print*,'problem while adding nodes 1'
                        STOP 'too many nodes added. Increase nnode_max'
                      endif                     
                     ! put the newly added node not at the exact midpoint, with a randomness of 5% of l, this helps for prevent flipping after re-delaunay triangulation,@YWang
                      call random_number(xx)
                      geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0+xx*l*5.d-3
                      geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0+xx*l*5.d-3
                      if (geometry%x(geometry%nnode).gt.xmax)then !stop nodes from being too close to right boundary
                        geometry%x(geometry%nnode) = xmax-xx
                      endif
                      if(geometry%x(geometry%nnode).lt.xmin)then !stop nodes from being too close to left boundary
                        geometry%x(geometry%nnode) = xmin+xx
                      endif
                      if (geometry%y(geometry%nnode).gt.ymax)then !stop nodes from being too close to top boundary
                        geometry%y(geometry%nnode) = ymax-xx
                      endif
                      if(geometry%y(geometry%nnode).lt.ymin)then !stop nodes from being too close to bottom boundary
                        geometry%y(geometry%nnode) = ymin+xx
                      endif                       
                      geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                      
                      ! update receive and donor relationship between i, j and the newly added node
                      network%receiver(geometry%nnode)=j
                      network%receiver(i)=geometry%nnode
                      geometry%fix(geometry%nnode)=0
                      
                      ! Assign physical properties to new node as average of its endmembers
                      geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                      geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                      geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                      geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0

                      geometry%surface(geometry%nnode)=0.0d0
                      geometry%discharge(geometry%nnode)=0.0d0
                      geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                      geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                      geometry%sediment_flux(geometry%nnode)=0.0d0
                      geometry%jitter_direction(geometry%nnode) = 0
                      !-------------------------------------
                      geometry%jitter_flag(geometry%nnode,1) = 2
                      geometry%jitter_flag(geometry%nnode,2) = 0
                      geometry%jitter_flag(geometry%nnode,2) = 0
                      geometry%jitter_residual(geometry%nnode) = 0.d0
                      !-------------------------------------
                      if (params%transient_divide) then
                         do h=1,params%num_bins
                            geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                 geometry%erosion_rate_history(j,h))/2.d0
                         enddo
                      endif ! finish updating all nodal physical properties 
                      endif ! end too small segment
                  endif !endif non-neighbour receiber but not crossing river
                  
               endif ! end if of non-neighbour receiver
          endif ! end if of an internal non-lake node
      enddo !end loop of nodes
  
     !need to continuously update donors and receivers
      network%nnode=geometry%nnode
      network%donors=0
      network%ndon=0
      do j=1,network%nnode
         k=network%receiver(j)
         if (k.ne.0) then          
            network%ndon(k)=network%ndon(k)+1
            if (network%ndon(k).gt.geometry%nnmax) call VTK_debug(geometry,delaunay,params,network,stack,'toomanydono_',params%istep)
            if (network%ndon(k).gt.geometry%nnmax) print*,'Problem is at node j coordinates of', geometry%x(j),geometry%y(j)
            if (network%ndon(k).gt.geometry%nnmax) print*,'Problem is at node k coordinates of', geometry%x(k),geometry%y(k)
            network%donors(network%ndon(k),k)=j
         endif
      enddo

      call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
      call find_network (geometry,delaunay,network,1) ! find neighbour list only
      call check_neighbor_receiver (geometry,network,ncount)
    
      if (ncount.eq.0) missing_node=0
      if (it.gt.50) then
          missing_node=0
          call VTK_debug(geometry,delaunay,params,network,stack,'stopflipping1_',params%istep)          
          STOP 'more than 50 iterations in add_and_remove_nodes routine!!!' 
      endif
      it=it+1
  enddo !end of while missing node loop
print*, 'number of add node for non-neighbour receiver connection is : ', nadd_river_no_connection
print*,'cross river after first non-neighbour receiver check:', crossriver
    !!!---------------------------- (2) ---------------------------- !!!
  ! add node on non-channel divides that are longer than a given distance
  ! add node to lower existing node
  ! note that this loop does not take into acount boundary nodes.
  nnold=geometry%nnode !original number of nodes
  nadd_between_channel = 0
  do i=1,nnold ! loop over the original nodes
      do k=geometry%nb(i),1,-1 ! loop over the neighbors
         j=geometry%nn(k,i)
         if (geometry%fix(i).eq.0.or.geometry%fix(j).eq.0) then !when both are boundary nodes treated with nadd_boundary
            if(network%receiver(i).ne.j.and.network%receiver(j).ne.i) then ! check non-channel
               if (geometry%z(i).lt.geometry%z(j)) then !add to lowest node
                  l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment between i and j
                  if (l.gt.params%ldivmax) then ! add node on the connection
                        nadd_between_channel=nadd_between_channel+1                        
                        geometry%nnode=geometry%nnode+1
                        if (geometry%nnode.gt.geometry%nnode_max) then
                           print*,'problem while adding nodes 2'
                           STOP 'too many nodes added. Increase nnode_max'
                        endif
                        geometry%x(geometry%nnode)=geometry%x(i)+(l-2.d0*params%xc)/2*(geometry%x(j)-geometry%x(i))/l
                        geometry%y(geometry%nnode)=geometry%y(i)+(l-2.d0*params%xc)/2*(geometry%y(j)-geometry%y(i))/l
                        geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                        ! Assign physical properties to new node from neighbor node i
                        !  sdw  Not sure if these are all the node-local properties
                        geometry%erosion_rate(geometry%nnode)=geometry%erosion_rate(i)
                        geometry%surface(geometry%nnode)=params%amin
                        geometry%u(geometry%nnode)=geometry%u(i)
                        geometry%v(geometry%nnode)=geometry%v(i)
                        geometry%w(geometry%nnode)=geometry%w(i)
                        geometry%erosion_rate(geometry%nnode)=geometry%erosion_rate(i)
                        geometry%fix(geometry%nnode)=0
                        geometry%discharge(geometry%nnode)=0.0d0
                        geometry%precipitation(geometry%nnode)=geometry%precipitation(i)
                        geometry%k(geometry%nnode)=geometry%k(i)
                        geometry%sediment_flux(geometry%nnode)=0.0d0
                        geometry%jitter_direction(geometry%nnode) = 0
                        network%ndon(i)=network%ndon(i)+1
                        network%ndon(geometry%nnode)=0
                        network%receiver(geometry%nnode)=i
                        !-------------------------------------
                        geometry%jitter_flag(geometry%nnode,1) = 2
                        geometry%jitter_flag(geometry%nnode,2) = 0
                        geometry%jitter_flag(geometry%nnode,3) = 0
                        geometry%jitter_residual(geometry%nnode) = 0.d0
                        !-------------------------------------
                        if (params%transient_divide) then
                           do h=1,params%num_bins
                              geometry%erosion_rate_history(geometry%nnode,h)=geometry%erosion_rate_history(i,h)
                           enddo
                        endif ! end updating physical properties of newly added node
                     endif ! end if internal non-connected two nodes
                  endif
               endif
            endif
      enddo! loop over the neighbours
  enddo! end loop over the original nodes

 
  !!!---------------------------- (3) ---------------------------- !!!
  !  This section is to add nodes on channels. It checks whether nodes on the drainage network
  ! are too far apart or not. It adds nodes when the distance between two nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90) - SDW
  
  ! add randomness to the added node if movepoint is set true @YWANG

  nadd_on_channel=0
  nnold = geometry%nnode
  do i=1,nnold ! loop over the nodes
     if (network%receiver(i).ne.0) then ! only for connections that are part of the drainage network
        j=network%receiver(i)
        l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length of segment between i and j
        if (l.gt.params%lmax) then ! add a node on the connection
           nadd_on_channel=nadd_on_channel+1
           geometry%nnode=geometry%nnode+1
           if (geometry%nnode.gt.geometry%nnode_max) then
              print*,'problem while adding nodes 3'
              STOP 'too many nodes added. Increase nnode_max'
           endif
           !-------add randomness to the added node,if there is advection, @YWANG
           if(params%move_points)then
               call random_number(xx)            
               geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0&
                       &+(2.d0*xx-1.d0)
               call random_number(xx)         
               geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0&
                       &+(2.d0*xx-1.d0)
           else
               geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0
               geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0               
           endif
           if (geometry%x(geometry%nnode).gt.xmax)then !stop nodes from being too close to right boundary
                geometry%x(geometry%nnode) = xmax-xx
           endif
           if(geometry%x(geometry%nnode).lt.xmin)then !stop nodes from being too close to left boundary
                geometry%x(geometry%nnode) = xmin+xx
           endif
           if (geometry%y(geometry%nnode).gt.ymax)then !stop nodes from being too close to top boundary
                geometry%y(geometry%nnode) = ymax-xx
           endif
           if(geometry%y(geometry%nnode).lt.ymin)then !stop nodes from being too close to bottom boundary
                geometry%y(geometry%nnode) = ymin+xx
           endif             
           geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
           !note every physical properties will have to be updated at each new
           !node....
           ! Assign physical properties to new node as average of its endmembers
           ! added 15.6.2010 sdw  Not sure if these are all the node-local properties
           geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
           geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
           geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
           geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
           network%receiver(geometry%nnode)=j
           network%receiver(i)=geometry%nnode
           geometry%fix(geometry%nnode)=0
           geometry%surface(geometry%nnode)=0.0d0
           geometry%discharge(geometry%nnode)=0.0d0
           geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
           geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
           geometry%sediment_flux(geometry%nnode)=0.0d0
           geometry%jitter_direction(geometry%nnode) = 0
           !-------------------------------------
           geometry%jitter_flag(geometry%nnode,1) = 2
           geometry%jitter_flag(geometry%nnode,2) = 0
           geometry%jitter_flag(geometry%nnode,3) = 0
           geometry%jitter_residual(geometry%nnode) = 0.d0
           !-------------------------------------
           if (params%transient_divide) then
              do h=1,params%num_bins
                 geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                      geometry%erosion_rate_history(j,h))/2.d0
              enddo
           endif

        endif
     endif
  enddo! loop over the nodes
  
  !!!---------------------------- (4) ---------------------------- !!!
  !   This section is to add boundary nodes when the boundaries are advected.
  ! It adds nodes when the distance between two boundary nodes is
  ! greater than a critical distance "lmax" (chosen in initialize_parameters.f90)
  ! for some reason the nb array for boundary nodes is not trustable. Need to generate a new list.
  ! Each node will point to its neighbor in anticlockwise direction.

  !divided into four section according top, bottom, right, top and left boundaries

  ! try again with the neighbors list 1/4/11
  nadd_boundary=0
  nnold = geometry%nnode
  if (params%move_points) then
     do i=1,nnold 
        if (geometry%fix(i).ge.1) then
           !if (geometry%y(i).eq.ymin) then ! the front boundary
           if (geometry%fix(i).eq.3) then ! the front boundary
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 ! neighbour j is a boundary node
                 if (geometry%fix(j).ge.1.and.geometry%y(j).eq.ymin.and.geometry%x(j).gt.geometry%x(i)) then !this is the closest neighbor to the right
                    l = dsqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes 4'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       !-------add randomness to the added boundary node, YWANG
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                                &+(2.d0*xx-1.d0)*l/8.d0      
                       geometry%y(geometry%nnode)=geometry%y(i)
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=geometry%fix(i)
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       geometry%jitter_direction(geometry%nnode) = 0
                       network%receiver(geometry%nnode)=0
                       !-------------------------------------
                       geometry%jitter_flag(geometry%nnode,1) = 2
                       geometry%jitter_flag(geometry%nnode,2) = 0
                       geometry%jitter_flag(geometry%nnode,3) = 0
                       geometry%jitter_residual(geometry%nnode) = 0.d0
                       !-------------------------------------
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif
                      
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           !if (geometry%x(i).eq.xmax) then ! the right boundary
           if (geometry%fix(i).eq.2) then ! the right boundary
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%x(j).eq.xmax.and.geometry%y(j).gt.geometry%y(i)) then !this is the closesent neighbor to the top
                    l = dsqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes 5'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       
                       !-------add randomness to the added node, YWANG
                       call random_number(xx)
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                                &+(2.d0*xx-1.d0)*l/8.d0
                       geometry%x(geometry%nnode)=geometry%x(i)
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       
                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmax and yi=ymin)
                       !geometry%fix(geometry%nnode)=geometry%fix(i)
                       if (geometry%y(j) .eq. ymax) then
                           geometry%fix(geometry%nnode)=geometry%fix(i)
                       else
                           geometry%fix(geometry%nnode)=geometry%fix(j)
                       endif
                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmax and yi=ymin)
                       
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       geometry%jitter_direction(geometry%nnode) = 0
                       network%receiver(geometry%nnode)=0
                       !-------------------------------------
                       geometry%jitter_flag(geometry%nnode,1) = 2
                       geometry%jitter_flag(geometry%nnode,2) = 0
                       geometry%jitter_flag(geometry%nnode,3) = 0
                       geometry%jitter_residual(geometry%nnode) = 0.d0
                       !-------------------------------------
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif

                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           !if (geometry%y(i).eq.ymax) then !the back boundary
           if (geometry%fix(i).eq.4) then !the back boundary
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%y(j).eq.ymax.and.geometry%x(j).lt.geometry%x(i)) then !this is the closesent neighbor to the left
                    l = dsqrt((geometry%x(i)-geometry%x(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes 6'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       !-------add randomness to the added node if there is advection, YWANG
                       call random_number(xx)
                       geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0 &
                                &+(2.d0*xx-1.d0)*l/8.d0                               
                       geometry%y(geometry%nnode)=geometry%y(i)
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0 !should be zero
                       ! print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       geometry%fix(geometry%nnode)=geometry%fix(i)
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       geometry%jitter_direction(geometry%nnode) = 0
                       network%receiver(geometry%nnode)=0
                       !-------------------------------------
                       geometry%jitter_flag(geometry%nnode,1) = 2
                       geometry%jitter_flag(geometry%nnode,2) = 0
                       geometry%jitter_flag(geometry%nnode,3) = 0
                       geometry%jitter_residual(geometry%nnode) = 0.d0
                       !-------------------------------------
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif

                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
           !if (geometry%x(i).eq.xmin) then ! the left boundary
           if (geometry%fix(i).eq.1) then ! the left boundary
              do k = 1,geometry%nb(i)
                 j = geometry%nn(k,i)
                 if (geometry%fix(j).ge.1.and.geometry%x(j).eq.xmin.and.geometry%y(j).lt.geometry%y(i)) then !this is the closesent neighbor to the bottom
                    l = dsqrt((geometry%y(i)-geometry%y(j))**2)
                    if (l.gt.params%lmax) then
                       nadd_boundary=nadd_boundary+1
                       geometry%nnode=geometry%nnode+1
                       if (geometry%nnode.gt.geometry%nnode_max) then
                          print*,'problem while adding nodes 7'
                          STOP 'too many nodes added. Increase nnode_max'
                       endif
                       
                       !-------add randomness to the added node if there is advection, YWANG
                       call random_number(xx)
                       geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0 &
                                &+(2.d0*xx-1.d0)*l/8.d0                               
                       geometry%x(geometry%nnode)=geometry%x(i)
                       geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                       ! print*, 'adding node', geometry%nnode, 'because distance between boundary nodes is too large'
                       !print*, geometry%x(geometry%nnode), geometry%y(geometry%nnode), geometry%z(geometry%nnode)
                       geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                       geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                       geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                       geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                       
                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmin and yi=ymin)
                       !geometry%fix(geometry%nnode)=geometry%fix(i)
                       if (geometry%y(j) .eq. ymin) then
                           geometry%fix(geometry%nnode)=geometry%fix(i)
                       else
                           geometry%fix(geometry%nnode)=geometry%fix(j)
                       endif
                       ! @YWANG JAN. 2022, fix situation of corner node (xi = xmin and yi=ymin)
                       
                       geometry%surface(geometry%nnode)=0.0d0
                       geometry%discharge(geometry%nnode)=0.0d0
                       geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                       geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                       geometry%sediment_flux(geometry%nnode)=0.0d0
                       geometry%jitter_direction(geometry%nnode) = 0
                       network%receiver(geometry%nnode)=0
                       !-------------------------------------
                       geometry%jitter_flag(geometry%nnode,1) = 2
                       geometry%jitter_flag(geometry%nnode,2) = 0
                       geometry%jitter_flag(geometry%nnode,3) = 0
                       geometry%jitter_residual(geometry%nnode) = 0.d0
                       !-------------------------------------
                       if (params%transient_divide) then
                          do h=1,params%num_bins
                             geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                                  geometry%erosion_rate_history(j,h))/2.d0
                          enddo
                       endif
                       
                    endif
                    exit !because the neighbor to the right was found.
                 endif
              enddo
           endif
        endif
     enddo
  endif
 
 
    ! update triangle mesh and neighbour list
  if(nadd_between_channel .gt. 0 .or. nadd_on_channel .gt. 0 .or. nadd_boundary .gt. 0 &
       &.or. nadd_river_no_connection.gt. 0 .or. crossriver .gt. 0)then
      print*, 'add between channel', nadd_between_channel
      print*, 'add on channel', nadd_on_channel
      print*, 'add boundary nodes ', nadd_boundary
      !need to continuously update donors and receivers
      network%nnode=geometry%nnode
      network%donors=0
      network%ndon=0
      do j=1,network%nnode
        k=network%receiver(j)
        if (k.ne.0) then
            network%ndon(k)=network%ndon(k)+1
            network%donors(network%ndon(k),k)=j
        endif
      enddo
      call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
      call find_network (geometry,delaunay,network,1) ! find neighbour list only
      call check_neighbor_receiver (geometry,network,ncount)
      !if(ncount .gt. 0)then
      !  call VTK_debug(geometry,delaunay,params,network,stack,'afteraddingnode',params%istep) 
      !endif
  endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! Second stage: take care of all the cases where a new node should be removed !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !!!---------------------------- (5) ---------------------------- !!!
  !this section removes nodes based on minimum surface area in headwater
  nrem_small_area=0
  nnold = geometry%nnode
  do i=1, nnold
     !because we remove nodes, need to check if we are still within the bounds of old nodes.
     if(i.le.geometry%nnode-nadd_on_channel-nadd_between_channel-nadd_river_no_connection-nadd_boundary)then
         if(geometry%fix(i).eq.0.and.network%ndon(i).eq.0.and.geometry%surface(i).lt.params%amin )then
           nrem_small_area=nrem_small_area+1
           geometry%nnode=geometry%nnode-1
           do j=1,i-1
              if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
           enddo
           if(i.le.geometry%nnode)then
              do j=i,geometry%nnode
                 geometry%x(j)=geometry%x(j+1)
                 geometry%y(j)=geometry%y(j+1)
                 geometry%z(j)=geometry%z(j+1)
                 geometry%u(j)=geometry%u(j+1)
                 geometry%v(j)=geometry%v(j+1)
                 geometry%w(j)=geometry%w(j+1)
                 geometry%fix(j)=geometry%fix(j+1)
                 geometry%surface(j)=geometry%surface(j+1)
                 geometry%discharge(j)=geometry%discharge(j+1)
                 geometry%precipitation(j)=geometry%precipitation(j+1)
                 geometry%k(j)=geometry%k(j+1)
                 geometry%nb(j)=geometry%nb(j+1)
                 geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                 geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                 geometry%jitter_direction(j) = geometry%jitter_direction(j+1)
                 network%ndon(j)=network%ndon(j+1)
                 !-------------------------------------
                 geometry%jitter_flag(j,1) = geometry%jitter_flag(j+1,1)
                 geometry%jitter_flag(j,2) = geometry%jitter_flag(j+1,2)
                 geometry%jitter_flag(j,3) = geometry%jitter_flag(j+1,3)
                 geometry%jitter_residual(j) = geometry%jitter_residual(j+1)
                 !-------------------------------------
                 if (network%receiver(j+1).gt.i) then
                    network%receiver(j)=network%receiver(j+1)-1
                 else
                    network%receiver(j)=network%receiver(j+1)
                 endif
              enddo
              if (params%transient_divide) then
                 do j=i,geometry%nnode
                    geometry%erosion_rate_history(j,:)=geometry%erosion_rate_history(j+1,:)
                 enddo
              endif

           endif
        endif
     endif
  enddo
 print*, 'number of removed small area nodes are:', nrem_small_area
!!! the next remove section operates also on nodes that were just now added.
!!!For that reason the donor receiver array needs to be updated
  if (nrem_small_area.gt.0.or.nadd_on_channel.gt.0&
       &.or.nadd_between_channel.gt.0.or.nadd_river_no_connection.gt.0.or.nadd_boundary.gt.0)then
     network%nnode=geometry%nnode
     network%donors=0 
     network%ndon=0
     do i=1,network%nnode
        k=network%receiver(i)
        if (k.ne.0) then
           network%ndon(k)=network%ndon(k)+1
           network%donors(network%ndon(k),k)=i
        endif
     enddo
  endif
  
    
  !!!---------------------------- (6) ---------------------------- !!!
! This section removes nodes that are too close to the boundary due to jitter or advection
!!! routine to remove nodes close to boundary, is necessary when jitter nodes
! Only internal nodes
  nrem_close_to_boundary = 0
  if ((params%move_points .or. params%jitter_nodes) ) then
      ! the min_dist seems need to be scaled with the jitter length (especially the jitter cell ratio, use bigger dist for big ratio)
      min_dist =  min(abs(params%max_adv*params%deltat),params%lmax) !if a distance of a node to the boundary < min_dist-->remove it      
      do i=1, geometry%nnode !also new nodes that were just added      
557     if (i.le.geometry%nnode) then ! i can grow above geometry%nnode because geometry%nnode can become smaller duringt this loop
            if( geometry%fix(i).eq.0 .or. & !internal nodes
                &geometry%y(i).eq.ymin.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.& !boundary nodes but not corners
                &geometry%y(i).eq.ymax.and.geometry%x(i).ne.xmin.and.geometry%x(i).ne.xmax.or.& !boundary nodes but not corners
                &geometry%x(i).eq.xmin.and.geometry%y(i).ne.ymin.and.geometry%y(i).ne.ymax.or.& !boundary nodes but not corners
                &geometry%x(i).eq.xmax.and.geometry%y(i).ne.ymin.and.geometry%y(i).ne.ymax) then!boundary nodes but not corners
              if (abs(geometry%x(i)-xmin).lt.min_dist.and.geometry%x(i).ne.xmin.or.& !close to left
                    &abs(geometry%x(i)-xmax).lt.min_dist.and.geometry%x(i).ne.xmax.or.& !close to right
                    &abs(geometry%y(i)-ymin).lt.min_dist.and.geometry%y(i).ne.ymin.or.& !close to bottom
                    &abs(geometry%y(i)-ymax).lt.min_dist.and.geometry%y(i).ne.ymax) then !close to top
                  nrem_close_to_boundary = nrem_close_to_boundary +1 
                  print*,'node to remove is:', geometry%x(i), geometry%y(i), geometry%fix(i)
                  geometry%nnode=geometry%nnode-1
                  ! first reorganize network
                  k = network%receiver(i)
                  if (k.ne.0) then !i is not a lake
                      do j = 1,network%ndon(i)
                          h = network%donors(j,i) ! h is a donor of i
                          network%receiver(h) = k ! now h is a donor of k
                      enddo
                  else! i is a lake
                      do j = 1,network%ndon(i)
                          h = network%donors(j,i) ! h is a donor of i
                          network%receiver(h) = 0 ! now h is a lake
                      enddo
                  endif
                  do j=1,i-1
                      if (network%receiver(j).gt.i) network%receiver(j)=network%receiver(j)-1
                  enddo
                  if(i.le.geometry%nnode)then !not the last node
                      do j=i,geometry%nnode
                          geometry%x(j)=geometry%x(j+1)
                          geometry%y(j)=geometry%y(j+1)
                          geometry%z(j)=geometry%z(j+1)
                          geometry%u(j)=geometry%u(j+1)
                          geometry%v(j)=geometry%v(j+1)
                          geometry%w(j)=geometry%w(j+1)
                          geometry%fix(j)=geometry%fix(j+1)
                          geometry%surface(j)=geometry%surface(j+1)
                          geometry%discharge(j)=geometry%discharge(j+1)
                          geometry%precipitation(j)=geometry%precipitation(j+1)
                          geometry%k(j)=geometry%k(j+1)
                          geometry%nb(j)=geometry%nb(j+1)
                          geometry%erosion_rate(j)=geometry%erosion_rate(j+1)
                          geometry%sediment_flux(j)=geometry%sediment_flux(j+1)
                          network%ndon(j)=network%ndon(j+1)
                          !-------------------------------------
                          geometry%jitter_direction(j) = geometry%jitter_direction(j+1)
                          geometry%jitter_flag(j,1) = geometry%jitter_flag(j+1,1)
                          geometry%jitter_flag(j,2) = geometry%jitter_flag(j+1,2)
                          geometry%jitter_flag(j,3) = geometry%jitter_flag(j+1,3)
                          geometry%jitter_residual(j) = geometry%jitter_residual(j+1)
                          !-------------------------------------
                          if (network%receiver(j+1).gt.i) then
                              network%receiver(j)=network%receiver(j+1)-1
                          else
                              network%receiver(j)=network%receiver(j+1)
                          endif
                      enddo
                      if (params%transient_divide) then
                          do j=i,geometry%nnode
                              geometry%erosion_rate_history(j,:)=geometry%erosion_rate_history(j+1,:)
                          enddo
                      endif
                  endif
                  !need to continuously update donors and receivers
                  network%nnode=geometry%nnode
                  network%donors=0
                  network%ndon=0
                  do j=1,network%nnode
                      k=network%receiver(j)
                      if (k.ne.0) then
                          network%ndon(k)=network%ndon(k)+1
                          network%donors(network%ndon(k),k)=j
                      endif
                  enddo
                  go to 557 !without increasing i becuase the node array is smaller.
              endif
            endif
        endif
      enddo ! end do nodes
  endif

print*,'number of removed nodes close to boundary is:', nrem_close_to_boundary

!-------------------------------------------------------------!
  !need to continuously update donors and receivers for the last section of add nodes
  network%nnode=geometry%nnode
  network%donors=0
  network%ndon=0
  do j=1,network%nnode
      k=network%receiver(j)
      if (k.ne.0) then
          network%ndon(k)=network%ndon(k)+1
          network%donors(network%ndon(k),k)=j
      endif
  enddo
  call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
  call find_network (geometry,delaunay,network,1) ! find neighbour list only
  call check_neighbor_receiver (geometry,network,ncount)
  !if (ncount.gt.0) call VTK_debug(geometry,delaunay,params,network,stack,'afterboundary',params%istep)

    
  !!!---------------------------- (7) ---------------------------- !!!
  !This section is to add nodes on channels where the channel connect two nodes that are not neighbors------- !!!
  crossriver=0
  do i=1,geometry%nnode ! loop over the nodes      
      if(geometry%fix(i).eq.0.and.network%receiver(i).ne.0)then ! deal with an internal non-lake node
           flag=0
           do k=geometry%nb(i),1,-1 ! loop over the neighbors
              if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
           enddo
           if(flag.eq.0)then !receiver is not a neighbor
              j=network%receiver(i)
              call check_crossrivers(geometry,network,params,delaunay, stack , i, flagcrossi, flagcrossj, ri, rj)
              print*,'ri and rj are:', ri, rj
              if(flagcrossi.or.flagcrossj) then
                  print*, 'network%receiver(i) before reroute,', j
                  call reroute_donor_in_cross_river(geometry,network,params,delaunay, stack , i, ri, rj)
                  crossriver=crossriver+1  
                  print*, 'network%receiver(i) after reroute,', network%receiver(i)
              endif              
           endif          
      endif!end deal with an internal non-lake node      
  enddo!end loop nodes  
print*,'last check number of cross rivers are:', crossriver

  !need to continuously update donors and receivers 
  network%nnode=geometry%nnode
  network%donors=0
  network%ndon=0
  do j=1,network%nnode
      k=network%receiver(j)
      if (k.ne.0) then
          network%ndon(k)=network%ndon(k)+1
          network%donors(network%ndon(k),k)=j
      endif
  enddo
  !call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
  call find_network (geometry,delaunay,network,1) ! find neighbour list only
  call check_neighbor_receiver (geometry,network,ncount)
  !if (ncount.gt.0) call VTK_debug(geometry,delaunay,params,network,stack,'afterboundary',params%istep)

  nadd_river_no_connection=0
  missing_node=0
  if (ncount.gt.0) missing_node = 1
  if (ncount.gt.0) print*, 'last check, number of non-neighbour connections is: ', ncount
  it=1
  do while (missing_node.eq.1)
    nnold = geometry%nnode
    do i=1,nnold ! loop over the nodes
      if(geometry%fix(i).eq.0.and.network%receiver(i).ne.0)then
            flag=0
            do k=geometry%nb(i),1,-1 ! loop over the neighbors
              if(network%receiver(i).eq.geometry%nn(k,i))flag=1  ! receiver is a neighbor
            enddo
            if(flag.eq.0)then !receiver is not a neighbor
              j=network%receiver(i)             
              !!!!!!!!!!!!!!!!!~~~~~~~~~~~~~check neighbour edge is a river, crossing edge~~~~~~~~~~~~~~~
              call check_crossrivers(geometry,network,params,delaunay, stack , i, flagcrossi, flagcrossj, ri, rj) 
              if(flagcrossi.or.flagcrossj) call reroute_donor_in_cross_river(geometry,network,params,delaunay, stack , i, ri, rj)              
              if(ri.eq.0 .and. rj.eq.0)then 
                l=dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2) ! l is length between i and j
                if(l.le.1.d1)then
                    network%receiver(i)=0
                else
                    nadd_river_no_connection=nadd_river_no_connection+1
                    geometry%nnode=geometry%nnode+1
                   if (geometry%nnode.gt.geometry%nnode_max) then
                      print*,'node xi and yi is', geometry%x(i), geometry%y(i)
                      print*,'problem while adding nodes 8'
                      STOP 'too many nodes added. Increase nnode_max'
                   endif
                  ! assign coordinates to the newly added node              
                  ! put the newly added node not at the exact midpoint, with a randomness of 5% of l, this helps for prevent flipping after re-delaunay triangulation,@YWang
                  call random_number(xx)
                  geometry%x(geometry%nnode)=(geometry%x(i)+geometry%x(j))/2.d0+xx*l*5.d-2
                  geometry%y(geometry%nnode)=(geometry%y(i)+geometry%y(j))/2.d0+xx*l*5.d-2
                  
                  if (geometry%x(geometry%nnode).gt.(xmax-min_dist))then !stop nodes from being too close to right boundary
                    geometry%x(geometry%nnode) = xmax-min_dist-xx
                  endif
                  if(geometry%x(geometry%nnode).lt.(xmin+min_dist))then !stop nodes from being too close to left boundary
                    geometry%x(geometry%nnode) = xmin+min_dist+xx
                  endif
                  if (geometry%y(geometry%nnode).gt.(ymax-min_dist))then !stop nodes from being too close to top boundary
                    geometry%y(geometry%nnode) = ymax-min_dist-xx
                  endif
                  if(geometry%y(geometry%nnode).lt.(ymin+min_dist))then !stop nodes from being too close to bottom boundary
                    geometry%y(geometry%nnode) = ymin+min_dist+xx
                  endif
              
                  print*, 'wrong connection between:', geometry%x(i), geometry%y(i), geometry%x(j), geometry%y(j)
                  print*, 'adding node to stop flipping 2 at:', geometry%x(geometry%nnode), geometry%y(geometry%nnode)
                  geometry%z(geometry%nnode)=(geometry%z(i)+geometry%z(j))/2.d0
                  ! Assign physical properties to new node as average of its endmembers
                  geometry%erosion_rate(geometry%nnode)=(geometry%erosion_rate(i)+geometry%erosion_rate(j))/2.d0
                  geometry%u(geometry%nnode)=(geometry%u(i)+geometry%u(j))/2.d0
                  geometry%v(geometry%nnode)=(geometry%v(i)+geometry%v(j))/2.d0
                  geometry%w(geometry%nnode)=(geometry%w(i)+geometry%w(j))/2.d0
                  network%receiver(geometry%nnode)=j
                  network%receiver(i)=geometry%nnode
                  geometry%fix(geometry%nnode)=0
                  geometry%surface(geometry%nnode)=0.0d0
                  geometry%discharge(geometry%nnode)=0.0d0
                  geometry%precipitation(geometry%nnode)=(geometry%precipitation(i)+geometry%precipitation(j))/2.d0
                  geometry%k(geometry%nnode)=(geometry%k(i)+geometry%k(j))/2.d0
                  geometry%sediment_flux(geometry%nnode)=0.0d0
                  geometry%jitter_direction(geometry%nnode) = 0
                  !-------------------------------------
                  geometry%jitter_flag(geometry%nnode,1) = 2
                  geometry%jitter_flag(geometry%nnode,2) = 0
                  geometry%jitter_flag(geometry%nnode,3) = 0
                  geometry%jitter_residual(geometry%nnode) = 0.d0
                  !-------------------------------------
                  if (params%transient_divide) then
                      do h=1,params%num_bins
                        geometry%erosion_rate_history(geometry%nnode,h)=(geometry%erosion_rate_history(i,h)+ &
                              geometry%erosion_rate_history(j,h))/2.d0
                      enddo
                  endif             
                endif! endif too small segment
            endif! endif flagcross
        endif!endif non-neighbour 
      endif !endif an internal non-lake node
    enddo ! enddo loop over nodes

    !need to continuously update donors and receivers
    network%nnode=geometry%nnode
    network%donors=0
    network%ndon=0
    do j=1,network%nnode
      k=network%receiver(j)
      if (k.ne.0) then
          network%ndon(k)=network%ndon(k)+1
          network%donors(network%ndon(k),k)=j
      endif
    enddo
    call calculate_delaunay (geometry,delaunay)     !redo delaunay triangulation
    call find_network (geometry,delaunay,network,1) ! find neighbour list only
    call check_neighbor_receiver (geometry,network,ncount)
    if (ncount.eq.0) missing_node=0
    if (it.gt.50) then
      missing_node=0
      call VTK_debug(geometry,delaunay,params,network,stack,'lastcheck',params%istep)    
      print*, 'more than 50 iterations in add_and_remove_nodes last check!!!'
      stop        
    endif! end 50 interations
    it = it+1    
  enddo !while loop




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! Short debuging section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(network%nnode.ne.geometry%nnode)print*, 'error in addremovenodes network nnode wrong'

  icount=0
  do i=1,geometry%nnode
     if(geometry%fix(i).ge.1)then
        icount=icount+1
        go to 555
     endif
     if(network%receiver(i).eq.0)then
        icount=icount+1
        go to 555
     endif
     j=network%receiver(i)
     if(j.gt.0.and.j.le.geometry%nnode)then
        icount=icount+1
        go to 555
     endif
     print*,' node not in network: ',i, network%receiver(i)
555  continue
  enddo

  call time_out ('add_remove_nodes')

  return
end subroutine add_remove_nodes




!!~~~~~~~~~~~~~~~~subroutine to reroute the donor of crossing river non-neighbour connections
subroutine reroute_donor_in_cross_river(geometry,network,params,delaunay, stack , i, ri, rj)
! this subroutine is to 
    use definitions
    implicit none
    type (geom) geometry
    type (netw) network
    type (parm) params
    type (del) delaunay
    type (stck) stack 
    
    integer i, j, k, ri, rj    
    integer flagi, flagi1, flagi2, flagj, flagj1, flagj2
    integer p, q, w,v
    double precision znbi

    flagi=0
    flagi1=0
    flagi2=0
    flagj=0
    flagj1=0
    flagj2=0    
    p=0
    q=0
    w=0
    v=0
    
   print*, 'reroute the two receivers of crossed rivers are:', ri, rj
    ! case1, crossed river nodes are neighbours of i
    j = network%receiver(i)
    if(ri.ne.0)then 
       if(geometry%z(i).ge.geometry%z(ri))then !connect to the lower receiver of crossed river
          network%receiver(i)=ri
          flagi = 1
          print*, 'reroute proceed to here 1', network%receiver(i)
          return
       else
          p=network%receiver(ri)
          if(p.ne.0)then
            if(geometry%z(i).ge.geometry%z(p) .and. p.ne.j .and. p.ne.i)then!connect to the 1st receiver of the lower receiver of crossed river
                network%receiver(i)=p
                flagi1=1
                print*, 'reroute proceed to here 11', network%receiver(i), geometry%x(p), geometry%y(p)
                return        
            endif               
          endif              
       endif          
    endif
     
    ! case2, crossed river nodes are neighbours of j
    if(ri.eq.0 .and. rj.ne.0)then
        if(geometry%z(i).ge.geometry%z(rj) .and. rj.ne.j)then
            network%receiver(i)=rj
            flagj=1
            print*, 'reroute proceed to here 2', network%receiver(i)
            return
        endif
        if(geometry%z(i).lt.geometry%z(rj) .and. rj.ne.j)then
            w=network%receiver(rj)
            if(w.ne.0)then
                if(geometry%z(i).ge.geometry%z(w) .and. w.ne.j .and. w.ne.i)then
                    network%receiver(i)=w
                    flagj1=1
                    print*, 'reroute proceed to here 22', network%receiver(i)
                    return       
                endif
            endif
        endif
    endif
    !case 3, make i to connect to i own lowest neighbour
    znbi=1.d10
    if(flagi.eq.0 .and. flagi1.eq.0 .and. flagj.eq.0 .and. flagj1.eq.0)then 
        do k=geometry%nb(i),1,-1 ! loop over the neighbors
            znbi = min(znbi,geometry%z(geometry%nn(k,i)))
        enddo
        if(znbi.le.geometry%z(i))then
            do k=geometry%nb(i),1,-1 ! loop over the neighbors
                if(geometry%z(geometry%nn(k,i)).eq.znbi)then
                    network%receiver(i)=geometry%nn(k,i)
                    flagi2=1
                    print*, 'reroute proceed to here 3', network%receiver(i)
                    return
                endif                    
            enddo
        endif
    endif 

    !case 4, make it a physical lake
    if(flagi.eq.0 .and. flagi1.eq.0 .and. flagj.eq.0 .and. flagj1.eq.0 .and.flagi2.eq.0)then ! can't find lower point to connect, make i a lake
        network%receiver(i)=0
        geometry%erosion_rate(i)=0.0d0
        geometry%sediment_flux(i) = 0.d0
        print*, 'reroute proceed to here 4', network%receiver(i)
        call VTK_debug(geometry,delaunay,params,network,stack,'artifactlakei_',params%istep)              
     endif 
               
end subroutine reroute_donor_in_cross_river


!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~subroutine to check crossing rivers
subroutine check_crossrivers(geometry,network,params,delaunay, stack , i, flagcrossi, flagcrossj, ri, rj)
  ! this subroutine is to check whether a non-neighbour receiver connection ("problem" connection) is crossing a river
  ! connection that the two nodes are neighbours of the problem connection nodes. --@YWang, September 2024
  
    use definitions
    implicit none
    type (geom) geometry
    type (netw) network
    type (parm) params
    type (del) delaunay
    type (stck) stack 
    
    integer i, j, k, ri, rj    
    logical flagcrossi, flagcrossj
    integer intersect, intersect2
    integer, external:: checkIntersection
    logical flag2, flag3
    
    !initialize flagcross to be false
    flagcrossi = .false.
    flagcrossj = .false.
    
    ! initialize receiver of the crossed river segment
    ri = 0
    rj = 0
    
    ! 1) check neighbours of j, see if there are two of them are connecting each other as a river
    j=network%receiver(i)
    do k=geometry%nb(j),1,-1 ! loop over the neighbors of j
        flag2 = any(network%receiver(geometry%nn(k,j)) .eq. geometry%nn(1:geometry%nb(j),j))             
        if (flag2 .and. network%receiver(geometry%nn(k,j)).ne.0 .and. geometry%nn(k,j).ne.0)then
            intersect = checkIntersection(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j), &
                                          geometry%x(geometry%nn(k,j)),geometry%y(geometry%nn(k,j)), &
                                          geometry%x(network%receiver(geometry%nn(k,j))),geometry%y(network%receiver(geometry%nn(k,j))))         
            if(intersect .eq. 1)then
                !call VTK_debug(geometry,delaunay,params,network,stack,'wrongconnection1_',params%istep)                
                flagcrossj = .true.   
                rj = network%receiver(geometry%nn(k,j)) ! the neighbour of j which is the receiver of the crossed river
                !print*, 'find non-neighbour connection crossing rivers j', rj                
            endif
        endif
    enddo 
    ! 2) check neighbours of i, see if there are two of them are connecting each other as a river
    do k=geometry%nb(i),1,-1 ! loop over the neighbors of i
        flag3 = any(network%receiver(geometry%nn(k,i)) .eq. geometry%nn(1:geometry%nb(i),i))
        if (flag3 .and. network%receiver(geometry%nn(k,i)).ne.0 .and. geometry%nn(k,i).ne.0)then
            intersect2 = checkIntersection(geometry%x(i),geometry%y(i),geometry%x(j),geometry%y(j), &
                                          geometry%x(geometry%nn(k,i)),geometry%y(geometry%nn(k,i)), &
                                          geometry%x(network%receiver(geometry%nn(k,i))),geometry%y(network%receiver(geometry%nn(k,i)))                          )   
             if(intersect2 .eq. 1)then
                !call VTK_debug(geometry,delaunay,params,network,stack,'wrongconnection2_',params%istep)                                
                flagcrossi = .true. 
                ri = network%receiver(geometry%nn(k,i))! the neighbour of i which is the receiver of the crossed river
                !print*, 'find non-neighbour connection crossing rivers i', ri                
            endif                                                 
        endif
    enddo

end subroutine check_crossrivers



  
 !!!~~~~~~~~~~~~~~~~~~~~small functions to determine whether two line segment intersect~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  function orientation(xp,yp,xq,yq,xr,yr)
    implicit none
    real(8), intent(in):: xp,yp,xq,yq,xr,yr
    integer :: orientation
    real:: val
    ! compute orientation value
    val = (yq-yp)*(xr-xq)-(xq-xp)*(yr-yq)
    if(val>0.0)then
        orientation=1 !counterclockwise
    elseif(val<0.0) then
        orientation=-1 !clockwise
    else
        orientation=0 !colinear
    endif    
  end function orientation
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function onSegment(ppx,ppy,qx,qy,rx,ry)
    implicit none
    real(8), intent(in) :: ppx,ppy,qx,qy,rx,ry
    integer :: onSegment
    logical :: condition1, condition2
    ! check if point R is on segment PQ
    condition1 = (min(ppx,qx)<=rx) .and. (rx<=max(ppx,qx))
    condition2 = (min(ppy,qy)<= ry) .and. (ry<=max(ppy,qy))
    if( condition1 .and. condition2)then
        onSegment = 1
    else
        onSegment = 0
    endif
    
  end function onSegment
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
   function checkIntersection(xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4)
    implicit none
    real(8), intent(in) :: xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4
    integer :: o1, o2, o3, o4
    integer, external :: orientation
    integer :: isOn1, isOn2, isOn3, isOn4
    integer, external :: onSegment
    integer checkIntersection
    !compute the four orientation
    o1 = orientation(xx1,yy1,xx2,yy2,xx3,yy3)
    o2 = orientation(xx1,yy1,xx2,yy2,xx4,yy4)
    o3 = orientation(xx1,yy1,xx4,yy4,xx3,yy3)
    o4 = orientation(xx2,yy2,xx4,yy4,xx3,yy3)
    
    isOn1= onSegment(xx1,yy1,xx2,yy2,xx3,yy3)
    isOn2= onSegment(xx1,yy1,xx2,yy2,xx4,yy4)
    isOn3= onSegment(xx1,yy1,xx4,yy4,xx3,yy3)
    isOn4= onSegment(xx2,yy2,xx4,yy4,xx3,yy3)
    
    ! general case
    if(o1/=o2 .and. o3/=o4)then
        checkIntersection = 1
        return
    endif
    ! special case
    if(o1==0 .and. isOn1.eq.1)then
        checkIntersection = 1
        return
    endif
    ! special case
    if(o2==0 .and. isOn2.eq.1)then
        checkIntersection = 1
        return
    endif   
    ! special case
    if(o3==0 .and.  isOn3.eq.1)then
        checkIntersection = 1
        return
    endif    
    ! special case
    if(o4==0 .and. isOn4.eq.1)then
        checkIntersection = 1
        return
    endif    
    ! no intersection
    checkIntersection = 0
  end function checkIntersection
  !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
