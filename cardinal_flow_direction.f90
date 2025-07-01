subroutine cardinal_flow_direction(geometry,params,network)
  
  use definitions

  implicit none

  type (geom) geometry
  type (parm) params
  type (netw) network
  
  integer n,i, di
  integer, dimension(geometry%nnode) :: nl, pp
  integer, dimension(:), allocatable :: shadowfix
  
  
  di = size(geometry%catchment)
  
  allocate(shadowfix(di))

  ! add pieces to deal with the special buffer zone in the back boundary when advecting to the back boundary
  ! with y+ advection, rivers flow internally close to the back boundary, @YWANG, April 21, 2022
  shadowfix = 0
  if (params%move_points ) then
      do i = 1, geometry%nnode
          shadowfix(i) = geometry%fix(i)
      enddo
  endif

  !calculate cardinal direction of flow for every pixel (i.e., which border each pixel flows towards) 1=west, 2=east, 4=north, 3=south
  nl=0
  geometry%card_flow = 0 !reset cardinal direction of flow
  do n=1,4
    if(params%move_points )then
        pp = pack(geometry%catchment, geometry%fix.eq.n, nl) !find all catchment numbers for each border
    endif
!    pp = pack(geometry%catchment, geometry%fix.eq.n, nl) !find all catchment numbers for each border
    do i=1,size(pp)
      if (pp(i).gt.0) then
        where(geometry%catchment.eq.pp(i))
          geometry%card_flow = n
        endwhere
      endif
    enddo
  enddo
  
  
  deallocate(shadowfix)


end subroutine cardinal_flow_direction
