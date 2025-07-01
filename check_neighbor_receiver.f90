subroutine check_neighbor_receiver (geometry,network,ncount)

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network

  integer, intent(inout) :: ncount
  integer i,flag,k
  double precision md

  ncount=0
  md = sum(geometry%discharge(1:geometry%nnode))/real(geometry%nnode)

  do i=1,geometry%nnode
    if(geometry%fix(i).eq.0.and.network%receiver(i).ne.0)then
      flag=0
      do k=geometry%nb(i),1,-1 ! loop over the neighbors
        if(network%receiver(i).eq.geometry%nn(k,i)) flag=1  ! receiver is a neighbor
      enddo
      if (flag.eq.0) then
        ncount=ncount+1
        !print*,'i and its receiver are at', geometry%x(i), geometry%y(i), geometry%x(network%receiver(i)),geometry%y(network%receiver(i))
      endif
      
      
      
    endif
  enddo

!    if (ncount.gt.0) print*, 'receiver is not a neighbor in ', ncount,' cases'

end subroutine check_neighbor_receiver
