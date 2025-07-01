subroutine projection_area(geometry,delaunay,params,network,stack)
  use definitions
  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  type (stck) stack

  integer i, j, k, n, l
  real(8), dimension(:), allocatable :: x, y, z
  integer, dimension(:), allocatable :: upi
  integer, dimension(:), allocatable :: proj_area
  real(8) area
  real(8), dimension(:), allocatable ::  horizontal_flux 
  real(8) , dimension(:), allocatable:: hull_x, hull_y
  integer:: hull_size
  
  allocate(proj_area(geometry%nnode))
  allocate(horizontal_flux(geometry%nnode))
  horizontal_flux = 0.d0 ! initialize with a super big number
    ! loop through all nodes
    do i = 1, geometry%nnode
        if(network%ndon(i).gt.5)then ! only compute for none-channel head node        
            ! find all its upstream nodes, and divides?
            n = network%ndon(i)+1
            print*, 'number of upstream nodes and projection area is here:', n
            
            area = 0.d0
            allocate(x(n))
            allocate(y(n))
            allocate(z(n))
            allocate(hull_x(n))
            allocate(hull_y(n))           
            allocate(upi(n))
            
            x(1) = geometry%x(i)
            y(1) = geometry%y(i)
            z(1) = geometry%z(i)            
            do k = 1, network%ndon(i)
                j = network%donors(k,i)
                x(k+1) = geometry%x(j)
                y(k+1) = geometry%y(j)
                z(k+1) = geometry%z(j) 
                upi(k+1) = j
            enddo                      
            !print*, 'number of nodes upstream is:', n, geometry%nnode, i, geometry%x(i), geometry%y(i)
            ! determine a direction using the cardinal flow
            if (geometry%card_flow(i).eq.1 .or. geometry%card_flow(i) .eq.2) then ! flow to left or right
                call compute_convex_hull(x, z, n, hull_x, hull_y)                
            endif
            if (geometry%card_flow(i).eq.3 .or. geometry%card_flow(i) .eq.4) then
                call compute_convex_hull(y, z, n, hull_x, hull_y)
            endif
            call compute_hull_area(hull_x, hull_y, n, area)
            proj_area(i) = area
            print*, 'number of upstream nodes and projection area is:', n, proj_area(i)
            deallocate(x)
            deallocate(y)
            deallocate(z)
            deallocate(hull_x)
            deallocate(hull_y)           
            deallocate(upi)         
        endif ! endif of non-channel nodes      
    enddo ! end loop through nodes
    
    ! write to vtk files
    horizontal_flux = 0.d0
    do i=1,geometry%nnode
        if (proj_area(i) .gt. 0.d0)then
            !horizontal_flux(i) = (geometry%sediment_flux(i)-geometry%w(1)*geometry%discharge(i))/proj_area(i)
            horizontal_flux(i) = 1.d0
        endif
    enddo
    call VTK_horizontal_flux(geometry,delaunay,params,network,stack, horizontal_flux, 'horizontal_flux',params%istep)

    deallocate(proj_area)
    deallocate(horizontal_flux)
    
endsubroutine projection_area



! ~~~~~ subroutine to compute the area of the convex hull
subroutine compute_hull_area(x, y, n, area)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in):: x(n), y(n) ! input point coordinates
    real(8), intent(out) :: area
    integer i
    real(8) sum1, sum2
    
    sum1 = 0.d0
    sum2 = 0.d0
    
    do i = 1, n-1
        sum1 = sum1 + x(i)*y(i+1)
        sum2 = sum2 +y(i)*x(i+1)
    enddo
    sum1 = sum1 +x(n)*y(1)
    sum2 = sum2 +y(n)*x(1)
    area = 0.5*abs(sum1-sum2)    
endsubroutine compute_hull_area



!~~~~~subroutine to find the convex hull of basin points
subroutine compute_convex_hull(x, y, n, hull_x, hull_y)
    implicit none
    integer, intent(in):: n ! number of input points
    real(8), intent(in):: x(n), y(n) ! input point coordinates
    real(8), intent(out):: hull_x(n), hull_y(n)
    integer :: hull_size
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
        q = mod(next,n) +1
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
    real(8), intent(in):: xp,yp,xq,yq,xr,yr
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




