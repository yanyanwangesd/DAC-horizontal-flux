subroutine projection_area(geometry,delaunay,params,network,stack)
  use definitions
  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  type (stck) stack

  integer i, j, k, n, l, is, icatch, ki, nbdonor, ji
  double precision, dimension(:), allocatable :: x, y, z
  integer, dimension(:), allocatable :: upi, donors
  double precision, dimension(:), allocatable :: proj_area
  double precision area
  double precision, dimension(:,:), allocatable ::  horizontal_flux 
  double precision, dimension(:), allocatable:: hull_x, hull_y
  integer hull_size

  integer stacknnode 
  integer nstack
  integer,dimension(:), allocatable :: stackorder
 
  character cs*8
  integer icount
  icount=params%istep
  write(cs,'(I8)') icount
  if (icount.lt.10)      cs(1:7)='0000000'
  if (icount.lt.100)     cs(1:6)='000000'
  if (icount.lt.1000)    cs(1:5)='00000'
  if (icount.lt.10000)   cs(1:4)='0000'
  if (icount.lt.100000)  cs(1:3)='000'
  if (icount.lt.1000000)  cs(1:2)='00'
  if (icount.lt.10000000)  cs(1:1)='0'
  
  
  allocate(proj_area(geometry%nnode))
  allocate(horizontal_flux(geometry%nnode,2))

! initialize 
  horizontal_flux=0.d0-9.99999d6 
  proj_area=0.d0

  ! calculate number of donors for each node
  allocate(donors(geometry%nnode))  
    donors = network%ndon
    do is=stack%nnode,1,-1
        j=stack%order(is)
        k=network%receiver(j)               
        if(k.ne.0)then
            donors(k)=donors(k)+donors(j)
        endif               
    enddo 
    
    ! loop through nodes to find donor nodes of each node
    do i=1, geometry%nnode  ! loop through all nodes          
        if(donors(i).ge.3 .and. geometry%card_flow(i).ne.0)then ! has a minimum number of contributing nodes to continue calculation         
            nstack = donors(i)+1
            allocate(stackorder(nstack))           
            stacknnode = 0           
            if(stacknnode .le. donors(i)+1 )then
                call add_to_stack_local(i,stacknnode,stackorder,network,nstack)
            endif
            allocate(x(nstack))
            allocate(y(nstack))           
            allocate(z(nstack))  
            nbdonor = 0
            do k=1,nstack
                ki = stackorder(k)
                nbdonor = nbdonor+1
                !print*, 'ki and its receiver', ki, network%receiver(ki)
                x(k) = geometry%x(ki)
                y(k) = geometry%y(ki)
                z(k) = geometry%z(ki)               
            enddo
   
            allocate(hull_x(nstack))
            allocate(hull_y(nstack))   
            n = nstack           
            ! determine a direction using the cardinal flow
            if (geometry%card_flow(i).eq.1 .or. geometry%card_flow(i) .eq.2 ) then ! flow to left or right boundary
                call compute_convex_hull(y, z, n, hull_x, hull_y, hull_size)                
            endif
            !print*,'found x in basin are:', x
            !print*,'found z in basin are:', z            
            if (geometry%card_flow(i).eq.3 .or. geometry%card_flow(i) .eq.4 ) then ! flow to top or bottom boundary
                call compute_convex_hull(x, z, n, hull_x, hull_y, hull_size)
            endif            
            call compute_hull_area(hull_x, hull_y, hull_size, area)
            !print*,'projection area is:', area
            proj_area(i) = area 
            
           
            deallocate(x)
            deallocate(y)
            deallocate(z)
            deallocate(hull_x)
            deallocate(hull_y)           
            
            deallocate(stackorder)
        endif !endif not channel head nodes
    enddo ! end loop through nodes
    
   
    ! calculate horizontal flux
    do i=1,geometry%nnode
        if (proj_area(i) .gt. 1.d0)then
            horizontal_flux(i,1) = (geometry%sediment_flux(i)-geometry%w(i)*geometry%discharge(i))/proj_area(i)
            horizontal_flux(i,2) = geometry%sediment_flux(i)/geometry%discharge(i)
        else
            horizontal_flux(i,1) = -9.99999d6
            horizontal_flux(i,2) = -9.99999d6
        endif            
    enddo
    ! make the small basins to have the same 
     do is=1,stack%nnode
        i = stack%order(is)
        if(horizontal_flux(i,2).lt.0.d0)then
            j = network%receiver(i)
            if(j.ne.0) then
                horizontal_flux(i,1) = horizontal_flux(j,1)
                horizontal_flux(i,2) = horizontal_flux(j,2)
            endif            
        endif
     enddo

    
    ! write to vtk files
    call VTK_horizontal_flux(geometry,delaunay,params,network,stack, horizontal_flux, proj_area, 'horizontal_flux',params%istep)
    open (900,file='ASCII/horizontal_flux'//cs,status='unknown')
    do i=1,geometry%nnode
        write(900,'(f19.7, f19.7, f15.1)') horizontal_flux(i,1), horizontal_flux(i,2), proj_area(i)
    enddo
    close (900)
    
    deallocate(proj_area)
    deallocate(horizontal_flux)
    deallocate(donors)   
    
endsubroutine projection_area



! ~~~~~ subroutine to compute the area of the convex hull
subroutine compute_hull_area(x, y, n, area)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in):: x(n), y(n) ! input point coordinates
    double precision, intent(out) :: area
    integer i
    double precision sum1, sum2
    sum1 = 0.d0
    sum2 = 0.d0
    area = 0.d0
    do i = 1, n-1
        sum1 = sum1 + x(i)*y(i+1)
        sum2 = sum2 +y(i)*x(i+1)
    enddo
    sum1 = sum1 +x(n)*y(1)
    sum2 = sum2 +y(n)*x(1)
    area = 0.5*abs(sum1-sum2)    
endsubroutine compute_hull_area



!~~~~~subroutine to find the convex hull of basin points
subroutine compute_convex_hull(x, y, n, hull_x, hull_y, hull_size)
    implicit none
    integer, intent(in):: n ! number of input points
    double precision, intent(in):: x(n), y(n) ! input point coordinates
    double precision, intent(out):: hull_x(n), hull_y(n)
    integer,intent(out) :: hull_size
    integer :: orientation2
    
    integer :: i, p, q, next
    integer :: hull_indices(n)
    
    ! use Gift Wrapping Algorithm to decide hull 
    p = 1
    do i = 2, n
        if(x(i)<x(p)) p = i ! find the leftmost point
    enddo    
    hull_size = 0
    next = p
    do        
        if(hull_size.le.n-1)then
        hull_size = hull_size + 1
        hull_indices(hull_size) = next
        q = mod(next,n)+1
        do i = 1, n
            if(orientation2(x(next), y(next), x(i), y(i), x(q), y(q)) <0) q = i
        enddo       
        next = q
        if(next == p) exit       
        endif

    enddo
    
    ! store hull points
    do i = 1, hull_size
        hull_x(i) = x(hull_indices(i))
        hull_y(i) = y(hull_indices(i))
    enddo

endsubroutine compute_convex_hull


! small function to find relative position of three points, p to q, q to r, and p to r
function orientation2(xp,yp,xq,yq,xr,yr)
    implicit none
    double precision, intent(in):: xp,yp,xq,yq,xr,yr
    integer :: orientation2
    real:: val
    ! compute orientation value
    val = (yq-yp)*(xr-xq)-(xq-xp)*(yr-yq)
    if(val>0.0)then
        orientation2=1 !counterclockwise
    elseif(val<0.0) then
        orientation2=-1 !clockwise
    else
        orientation2=0 !colinear
    endif    
end function orientation2



! build stack order for nodes locally
recursive subroutine add_to_stack_local(i,stacknnode,stackorder,network,nstack)
! Recursive routine to go through the node upstream and keeping track of confluences
use definitions
implicit none
type (netw) network

integer stacknnode
integer nstack
integer i, j
integer,dimension(nstack) :: stackorder
 
stacknnode=stacknnode+1
if (stacknnode > nstack) then
    print *, "Error: stacknnode exceeds nstack"
    return
endif

stackorder(stacknnode)=i

if (network%ndon(i).eq.0) return
do j=1,network%ndon(i)
    call add_to_stack_local(network%donors(j,i),stacknnode,stackorder,network,nstack)
enddo
return

end subroutine add_to_stack_local

