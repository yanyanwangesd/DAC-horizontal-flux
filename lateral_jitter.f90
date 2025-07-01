subroutine lateral_jitter (geometry,network,params,delaunay)
  ! gives nodes a random perturbation of length that is a function of grid spacing and scaled valley width
  !
  !The valley width(vw) of each node is a normal distribution that centered at node i (mean = 0), with variance of vw/4
  !
  ! geometry%jitter_flag is an integer array. The first dimension of this array marks nodes to jitter (0), not jitter (1), jitter not specified(-1)
  !                                           The second dimension stores accumulated counts that this node has been jittered since Time 0 of model run 


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay
  type (stck) stack


  integer n,i,j,k,r,ncount,tri,ptr,jit_dir,d,jcount,keep,outside,kk
  double precision xx, t, md, maxd, ndd
  double precision m1,m2,b1,b2 !for lines
  double precision  vw, c, kc !valley width parameters
  double precision L, Lpmax,Lqmax, Lvw !lengths normal to downstream channel connection
  double precision px, py, qx, qy, pxe, pye, qxe, qye
  double precision, dimension(2) :: tx, ty !point coordinates
  double precision, external::  r8_normal_ab
  double precision :: length, Ll !length between i and j, ratio of lengths
  
  integer ival, icount,xpcount,xqcount
  integer rdcount

  integer, dimension(:), allocatable :: flag, norder, hitwall
  double precision,dimension(:),allocatable :: jitter_L
  integer pi
  integer ik
  
  
  integer sticount
  character cs*8
  
  double precision ndt, tauv, rwsigma, rwdif



  allocate(flag(geometry%nnode), norder(geometry%nnode),jitter_L(geometry%nnode),hitwall(geometry%nnode))

  ! -------------------------- part 1 --------------------------
  ! set constants
  ! -------------------------- part 1 --------------------------
  !define coefficient and exponent of valley width-discharge scalling
  kc = 1.0d-2
  c = 0.25

  keep=0

  ncount = 0
  jcount = 0
  
  xpcount = 0
  xqcount = 0

  md = sum(geometry%discharge(1:geometry%nnode))/real(geometry%nnode)
  maxd = maxval(geometry%discharge(1:geometry%nnode))
  
  
  ! -------------------------- part 2 --------------------------
  ! select nodes for jittering, nodes that are not neighbours to each other are selected
  ! -------------------------- part 2 --------------------------
  ! initiate flag array
  flag = -1
  ! random pick geometry%nnode times of node in array and update flag
  rdcount = 0
  do while(jcount<geometry%nnode .and. rdcount<50)
      999 call random_number(xx)
      xx = xx*geometry%nnode
      pi = nint(xx)
      ! not first node
      if (pi .eq. 0)then
          go to 999
      endif
      ! not a boundary node
      if (geometry%fix(pi) .ne. 0)then
          go to 999
      endif
      
      ! initialize node pi
      if (flag(pi) .eq. -1) then
          flag(pi) = 1
          do j=1,geometry%nb(pi)
              k=geometry%nn(j,pi)
              flag(k) = 0
          enddo
          rdcount = 0
      else
          rdcount = rdcount + 1
      endif
      jcount = jcount + 1
  enddo
  
! !- - - - - -this piece was to make sure all non-neighbour nodes are selected
!   ! sort nodes based on the last pick of node pi
!   do i=1,geometry%nnode-pi+1 !loop over nodes
!       norder(i) = pi+i-1
!   enddo
!   do i=geometry%nnode-pi+2,geometry%nnode !loop over nodes
!       norder(i) = i-(geometry%nnode-pi+1)
!   enddo

!   ! continue to label the other node 
!   do ik = 2,geometry%nnode !loop over nodes
!       i = norder(ik) ! node index
!       if (flag(i) .eq. -1) then
!           flag(i) = 0
!           do j=1,geometry%nb(i)
!               k=geometry%nn(j,i)
!               flag(k) = 1
!           enddo
!       endif
!   enddo
! !- - - - - -this piece was to make sure all non-neighbour nodes are selected

  ! update the number of steps since last jitter
  do i=1,geometry%nnode
      if(geometry%jitter_flag(i,1).ne.1) geometry%jitter_flag(i,3) = geometry%jitter_flag(i,3)+1 ! accumulate
      if(flag(i).ne.1 .and. geometry%jitter_flag(i,1).eq.1) geometry%jitter_flag(i,3) = 2 
      if(flag(i).eq.1 .and. geometry%jitter_flag(i,1).eq.1) geometry%jitter_flag(i,3) = 1 !reset to 1
      if(geometry%fix(i).ne.0) geometry%jitter_flag(i,3) = 0 ! boundary nodes are set to 0
  enddo
  
  ! update flag of whether it is jittered for the current step, flag 1 is jitter, flag 0 is not jitter, -1 is not selected in the nnode times of random pick
  do i=1,geometry%nnode
      geometry%jitter_flag(i,1) = flag(i)
  enddo
  
   ! -------------------------- part 3 --------------------------
   ! Jitter nodes
   ! -------------------------- part 3 --------------------------

!   ! rewrite jitter directions for nodes not to be jittered
!   do i=1,geometry%nnode
!       if (flag(i) .ne. 1) then
!           geometry%jitter_direction(i) = 0
!       endif
!   enddo
  geometry%jitter_direction = 0
  jitter_L = 0.d0
  hitwall = 0
  ! start jitter nodes
  do i=1,geometry%nnode !loop over nodes
    ! only jitter selected nodes
      if (flag(i) .eq. 1) then
          !no boundary nodes, no headwater nodes
          if (geometry%fix(i).eq.0.and.network%receiver(i).gt.0.and.network%ndon(i).ge.1) then
              !only jitter nodes flowing to the boundaries
              if (geometry%card_flow(i).eq.3.or.geometry%card_flow(i).eq.4.or.geometry%card_flow(i).eq.1.or.geometry%card_flow(i).eq.2) then 
                  
                  ! 1) calculate lateral erosion length based on timestep length and uplift rate
!                   vw = geometry%w(i)*params%deltat
                  vw = kc*geometry%discharge(i)**c
    
                  ! 2) use a very small L to find the two triangles to jitter to
                  keep = 0
                  L = 0.5 !pick a really short distance just to determine what triangle should be used for edge
                  do k=geometry%nb(i),1,-1 ! loop over the neighbors
                      ! find the receiver
                      j=geometry%nn(k,i)
                      if (network%receiver(i).eq.j) then
                          keep=1
                          length=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length betwenn i and j
                          Ll = L/length
                          ! find point p at the right jitter direction
                          px = -Ll*(geometry%y(i)-geometry%y(j)) + geometry%x(i)
                          py = Ll*(geometry%x(i)-geometry%x(j)) + geometry%y(i)
                          ! find pint q at the left jitter direction
                          qx = Ll*(geometry%y(i)-geometry%y(j)) + geometry%x(i)
                          qy = -Ll*(geometry%x(i)-geometry%x(j)) + geometry%y(i)
                      endif
                  enddo
                  if (keep.eq.0) then
                      print*, 'error in lateral jitter: receiver is not a neighbor, i=', i
                  endif
                  
                  ! 2.2) find the two triangles that point p and q are in
                  ! check point p and find the new point p which is on the edge
                  icount = 0
                  do k=1,geometry%nt(i) ! go over all triangles that share i
                      tri=geometry%nnodetri(k,i)
                      ptr=0
                      call point_in_triangle(px, py, geometry, delaunay, tri, ptr)
                      
                      if (ptr.eq.1) then
                          icount = icount+1
                          kk = 1
                          do ncount=1,3
                              if (delaunay%icon(ncount,tri).ne.i) then ! find the vertices that are not i (should occur twice)
                                  tx(kk)=geometry%x(delaunay%icon(ncount,tri))
                                  ty(kk)=geometry%y(delaunay%icon(ncount,tri))
                                  kk = kk+1
                              endif
                          enddo
    
                         !find equation for line connecting the two vertices
                         if (tx(1) .eq. tx(2))then
                              
                              !find equation for line connecting point i with old point p (that is outside of the Voronoi cell but in the correct direction)
                             if (px .eq.geometry%x(i)) then
                                  print*, 'both edge and line ip vertical, no intersection, no jitter applied'
                                  pxe = geometry%x(i)
                                  pye = geometry%y(i)
                             else
                                  pxe = tx(1)
                                  m2 = (py-geometry%y(i))/(px-geometry%x(i))
                                  b2 = py-m2*px
                                  !find point on edge that is the intersection of two above lines
                                  pye = m2*pxe + b2
                             endif
                         else
                             if (px .eq.geometry%x(i)) then
                                 print*,'only line ip is vertical, normal jitter'
                                 m1 = (ty(1)-ty(2))/(tx(1)-tx(2))
                                 b1 = ty(1)-m1*tx(1)
                                 pxe = geometry%x(i)
                                 pye = m1*pxe + b1
                             else
                                 m1 = (ty(1)-ty(2))/(tx(1)-tx(2))
                                 b1 = ty(1)-m1*tx(1)
                                 !find equation for line connecting point i with old point p (that is outside of the Voronoi cell but in the correct direction)
                                 m2 = (py-geometry%y(i))/(px-geometry%x(i))
                                 b2 = py-m2*px
                                 !find point on edge that is the intersection of two above lines
                                 pxe = (b2-b1)/(m1-m2)
                                 pye = m1*pxe + b1
                             endif
                         endif
                      endif
                  enddo
                  
                  ! check point q and find the new point q which is on the edge
                  icount = 0
                  do k=1,geometry%nt(i) ! go over all triangles that share i
                      tri=geometry%nnodetri(k,i)
                      ptr=0
                      call point_in_triangle(qx, qy, geometry, delaunay, tri, ptr)
                      if (ptr.eq.1) then
                          icount = icount+1
                          kk = 1
                          do ncount=1,3
                              if (delaunay%icon(ncount,tri).ne.i) then ! find the vertices that are not i (should occur twice)
                                  tx(kk)=geometry%x(delaunay%icon(ncount,tri))
                                  ty(kk)=geometry%y(delaunay%icon(ncount,tri))
                                  kk = kk+1
                              endif
                          enddo
    
                         !find equation for line connecting the two vertices
                         if (tx(1) .eq. tx(2))then
                              
                              !find equation for line connecting point i and p
                             if (qx .eq.geometry%x(i)) then
                                  print*, 'both edge and line iq vertical, no intersection, no jitter applied'
                                  qxe = geometry%x(i)
                                  qye = geometry%y(i)
                             else
                                  qxe = tx(1)
                                  m2 = (qy-geometry%y(i))/(qx-geometry%x(i))
                                  b2 = qy-m2*qx
                                  !find point on edge that is the intersection of two above lines
                                  qye = m2*qxe + b2
                             endif
                         else
                             if (qx .eq.geometry%x(i)) then
                                 print*,'only line iq is vertical, normal jitter'
                                 m1 = (ty(1)-ty(2))/(tx(1)-tx(2))
                                 b1 = ty(1)-m1*tx(1)
                                 qxe = geometry%x(i)
                                 qye = m1*qxe + b1
                             else
                                 m1 = (ty(1)-ty(2))/(tx(1)-tx(2))
                                 b1 = ty(1)-m1*tx(1)
                                 !find equation for line connecting point i and q
                                 m2 = (qy-geometry%y(i))/(qx-geometry%x(i))
                                 b2 = qy-m2*qx
                                 !find point on edge that is the intersection of two above lines
                                 qxe = (b2-b1)/(m1-m2)
                                 qye = m1*qxe + b1
                             endif
                         endif
                      endif
                  enddo
                  
                  ! 2.3) calculate distance ip and iq
                  Lpmax = dsqrt((pxe-geometry%x(i))**2 + (pye - geometry%y(i))**2)
                  Lqmax = dsqrt((qxe-geometry%x(i))**2 + (qye - geometry%y(i))**2)
                  ! 2.4) use a bias to decide jitter direction
!                   call assign_jitter_direction(i,jit_dir,geometry,network,delaunay,0)

                  ! 3) generate a normal distribution jitter range
!                   ! option 1: use random walk thoery
                  if (geometry%jitter_flag(i,3).eq.0)then
                      ndt = params%deltat/2
                  else
                      ndt = geometry%jitter_flag(i,3)*params%deltat
                  endif
                   rwdif = params%jitter_diffusivity*vw ! diffusivity is a function of valley width
!                   rwdif = params%jitter_diffusivity! diffusivity is a constant
                  rwsigma = sqrt(2.d0*rwdif*ndt)
                  xx = r8_normal_ab (0.d0, rwsigma)
                  ! option 2: random jitter from a gaussian distribution
!                   xx = r8_normal_ab (0.d0, vw/4)
                  ! option 3: random jitter from a uniform distribution
!                   call random_number(xx)
!                   xx = (xx*2.d0-1.d0)*vw/4

                  ! jitter according to directional residual
!                   xx = r8_normal_ab (0.d0, vw/4)
!                   xx = xx+geometry%jitter_residual(i) ! compensate the directional jitter residual 
                  if (xx < 0.d0 ) then
                      L = min(abs(xx), Lpmax*params%jitter_range)
                      jit_dir = -1
                      !jitter length to right, is negative, respecting the jitter direction value
                      jitter_L(i) = -1.d0*L
                      if(abs(xx).gt.Lpmax*params%jitter_range) then
                          hitwall(i) = hitwall(i)+1
                          !right residual, is negative, respecting the jitter direction value
                          geometry%jitter_residual(i) = Lpmax*params%jitter_range-abs(xx)
                      else
                          ! no residual
                          geometry%jitter_residual(i) = 0.d0
                      endif
                  else
                      L = min(abs(xx), Lqmax*params%jitter_range)
                      jit_dir = 1
                      !jitter length to left, is positive, respecting the jitter direction value
                      jitter_L(i) = 1.d0*L
                      if(abs(xx).gt.Lqmax*params%jitter_range) then 
                          hitwall(i) = hitwall(i)+1
                          !left residual, is positive, respecting the jitter direction value
                          geometry%jitter_residual(i) = abs(xx)-Lqmax*params%jitter_range
                      else
                          !no residual
                          geometry%jitter_residual(i) = 0.d0
                      endif
                  endif
                  
                  
                  ! 4) calculate the new coordinates according to jitter direction and jitter range L
                  do k=geometry%nb(i),1,-1 ! loop over the neighbors
                      ! find the receiver
                      j=geometry%nn(k,i)
                      if (network%receiver(i).eq.j) then
                          length=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length between i and j
                          Ll = jitter_L(i)/length
                          geometry%x(i) = Ll*(geometry%y(i)-geometry%y(j)) + geometry%x(i)
                          geometry%y(i) = -Ll*(geometry%x(i)-geometry%x(j)) + geometry%y(i)
                          geometry%jitter_direction(i)= jit_dir
                      endif
                  enddo
                  
                  !update total jittered times
                  geometry%jitter_flag(i,2) = geometry%jitter_flag(i,2)+1

              endif !end if of nodes not internal draining
          endif ! end if not boundary and not channel head nodes
      endif ! end if of selected nodes
      
      
  enddo ! end loop of all nodes
  
  ! write jitter data to ascii
  if ((params%istep/params%freq)*params%freq.eq.params%istep) then
      sticount = params%istep
      write(cs,'(I8)') sticount
      if (sticount.lt.10)      cs(1:7)='0000000'
      if (sticount.lt.100)     cs(1:6)='000000'
      if (sticount.lt.1000)    cs(1:5)='00000'
      if (sticount.lt.10000)   cs(1:4)='0000'
      if (sticount.lt.100000)  cs(1:3)='000'
      if (sticount.lt.1000000)  cs(1:2)='00'
      if (sticount.lt.10000000)  cs(1:1)='0'
      open (888,file='ASCII/jitter_data'//cs,status='unknown')
      do i=1,geometry%nnode
         write(888,'(3f16.3, i5,i5, i8, i5, f10.5, i5, f15.0)') geometry%x(i),geometry%y(i),geometry%z(i), geometry%jitter_flag(i,1), geometry%jitter_flag(i,2), geometry%jitter_flag(i,3), geometry%jitter_direction(i), jitter_L(i), hitwall(i),geometry%discharge(i)
      enddo
      close(888)
  endif
  
  deallocate(flag, norder,jitter_L,hitwall)

end subroutine lateral_jitter





!*****************************************************************************
! A lot of small functions to solve geometry problems
!*****************************************************************************
subroutine assign_jitter_direction(i,jit_dir,geometry,network,delaunay,d)
  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay

  integer, intent(in) :: i
  integer, intent(inout) :: jit_dir
  integer k,j,old_jitter,d,a
  double precision, dimension(2) :: dsv, nv !downstream vector, normal vector
  double precision :: xx, bias_old, bias_us !random number


  bias_old=0.5d0 !bias towards the old direction
  bias_us=0.5d0 !bias towards direction receiver went (for the first time)

  if (d.eq.0) then !otherwise jit dir is previously assigned
      !assign the jitter direction
      if (geometry%jitter_direction(i).eq.0) then !this is the first time this node has been jittered
          ! if (network%ndon(i).eq.1) then !if only one donor (so probably a new node)
          do a=1,network%ndon(i)
              if (geometry%jitter_direction(network%donors(a,i)).gt.0) then !if any donor has been jittered before
                  call random_number(xx)
                  if (xx.le.bias_us) then
                      jit_dir = geometry%jitter_direction(network%donors(a,i))
                  else
                      jit_dir = -geometry%jitter_direction(network%donors(a,i))
                  endif
              endif
          enddo
          !if still here, no donors had been jittered before, so just pick direction randomly
          !pick it randomly
          call random_number(xx)
          if (xx.gt.0.5) then
              jit_dir = 1
          else
              jit_dir = -1
          endif

      !this is not the first time this node has been jittered
      else !bias jitter direction based on previous direction
          call random_number(xx)
          if (xx.le.bias_old) then
              jit_dir = geometry%jitter_direction(i)
          else
              jit_dir = -geometry%jitter_direction(i)
          endif
      endif
  endif
  
end subroutine assign_jitter_direction


!*****************************************************************************
subroutine assign_new_coords(i,px,py,jit_dir,geometry,network,delaunay,L,keep,d)

  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (parm) params
  type (del) delaunay

  integer, intent(in) :: i
  integer, intent(inout) :: jit_dir, keep
  double precision, intent(in) :: L
  double precision, intent(inout) :: px,py
  integer k,j,check,d,a
  double precision, dimension(2) :: dsv, nv !downstream vector, normal vector
  double precision :: length, Ll !length between i and j, ratio of lengths
  double precision :: xx, bias_old, bias_us !random number

  check=0
  keep=0
  px=0
  py=0
  bias_old=0.9 !bias towards the old direction
  bias_us=0.8 !bias towards direction receiver went (for the first time)

  if (d.eq.0) then !otherwise jit dir is previously assigned
    !assign the jitter direction
    if (geometry%jitter_direction(i).eq.0) then !this is the first time this node has been jittered
      ! if (network%ndon(i).eq.1) then !if only one donor (so probably a new node)
      do a=1,network%ndon(i)
        if (geometry%jitter_direction(network%donors(a,i)).gt.0) then !if any donor has been jittered before
          call random_number(xx)
          if (xx.le.bias_us) then
            jit_dir = geometry%jitter_direction(network%donors(a,i))
          else
            jit_dir = -geometry%jitter_direction(network%donors(a,i))
          endif
          go to 111
        endif
      enddo
      !if still here, no donors had been jittered before, so just pick direction randomly
        !pick it randomly
        call random_number(xx)
        if (xx.gt.0.5) then
          jit_dir = 1
        else
          jit_dir = -1
        endif

    !this is not the first time this node has been jittered
    else !bias jitter direction based on previous direction
      call random_number(xx)
      if (xx.le.bias_old) then
        jit_dir = geometry%jitter_direction(i)
      else
        jit_dir = -geometry%jitter_direction(i)
      endif
    endif
  endif

  111 continue

  !find the receiver
  do k=geometry%nb(i),1,-1 ! loop over the neighbors
    j=geometry%nn(k,i)
    if (network%receiver(i).eq.j) then
      keep=1
      length=dsqrt((geometry%x(i)-geometry%x(j))**2.d0+(geometry%y(i)-geometry%y(j))**2.d0) ! l is length betwenn i and j
      Ll = L/length
      ! using like triangles on either side of right angle
      if (jit_dir.eq.1) then !jitter node i to river right
        px = -Ll*(geometry%y(i)-geometry%y(j)) + geometry%x(i)
        py = Ll*(geometry%x(i)-geometry%x(j)) + geometry%y(i)
      else !jitter node i to river left
        px = Ll*(geometry%y(i)-geometry%y(j)) + geometry%x(i)
        py = -Ll*(geometry%x(i)-geometry%x(j)) + geometry%y(i)
      endif !pick jitter direction
    endif !check j is receiver
  enddo !neighbor loop
  if (keep.eq.0) then
    print*, 'error in lateral jitter: receiver is not a neighbor, i=', i
  endif


end subroutine assign_new_coords


!*****************************************************************************
subroutine point_in_triangle(px, py, geometry, delaunay, triangle, ptr)

  use definitions

  implicit none

  type (geom) geometry
  type (del) delaunay

  double precision ax,bx,cx,ay,by,cy,w1,w2,wa
  integer, intent(inout) :: ptr
  integer, intent(in) :: triangle
  double precision, intent(in) :: px, py


  ax = geometry%x(delaunay%icon(1,triangle))
  bx = geometry%x(delaunay%icon(2,triangle))
  cx = geometry%x(delaunay%icon(3,triangle))
  ay = geometry%y(delaunay%icon(1,triangle))
  by = geometry%y(delaunay%icon(2,triangle))
  cy = geometry%y(delaunay%icon(3,triangle))


  w1 = (ax*(cy-ay)+(py-ay)*(cx-ax)-px*(cy-ay)) / ((by-ay)*(cx-ax)-(bx-ax)*(cy-ay))
  w2 = (py-ay-w1*(by-ay)) / (cy-ay)
  wa = w1+w2

  if (w1.ge.0.and.w2.ge.0.and.wa.le.1) then
    ptr=1
  endif

end subroutine point_in_triangle


!*****************************************************************************
!!! function of generating a Gaussian random real mumber
function r8_normal_ab ( a, b )

!
!! R8_NORMAL_AB returns a scaled pseudonormal R8.
!
!  Discussion:
!
!    The normal probability distribution function (PDF) is sampled,
!    with mean A and standard deviation B.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!    06 August 2013
!  Author:
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A, the mean of the PDF.
!    Input, real ( kind = rk ) B, the standard deviation of the PDF.
!    Output, real ( kind = rk ) R8_NORMAL_AB, a sample of the normal PDF.
!
  implicit none
  integer, parameter :: rk = 8
  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) r8_normal_ab
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) x

  call random_number ( harvest = r1 )
  call random_number ( harvest = r2 )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
  r8_normal_ab = a + b * x
  return
end function r8_normal_ab





