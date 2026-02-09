program main
  use omp_lib
  use basicmod
  use mpimod
  use boundarymod
  implicit none
  real(8)::time_begin,time_end
  logical::is_final
  logical,parameter::nooutput=.false.
  logical,parameter:: forceoutput=.true., usualoutput=.false.
  data is_final /.false./
  call InitializeMPI
  if(myid_w == 0) print *, "setup grids and fields"
  if(myid_w == 0) print *, "grid size for x y z",ngrid1*ntiles(1),ngrid2*ntiles(2),ngrid3*ntiles(3)
  if(myid_w == 0 .and. nooutput ) print *, "Intermediate results are not outputed"
  call GenerateGrid
  call GenerateProblem
  call ConsvVariable
  if(myid_w == 0) print *, "entering main loop"
! main loop
  if(myid_w == 0 .and. .not. nooutput )                        print *,"step ","time ","dt"
  time_begin = omp_get_wtime()
  mloop: do nhy=1,nhymax
     call TimestepControl
     if(mod(nhy,nhydis) .eq. 0  .and. .not. nooutput .and. myid_w == 0) print *,nhy,time,dt
     call BoundaryCondition
     call StateVevtor
     call GravForce
     call EvaulateCh
     call NumericalFlux1
     call NumericalFlux2
     call NumericalFlux3
     call UpdateConsv
     call DampPsi
     call PrimVariable
     time=time+dt
     if(.not. nooutput ) call Output(usualoutput)
     if(time > timemax) exit mloop
  enddo mloop

  time_end = omp_get_wtime()
      
  if(myid_w == 0) print *, "sim time [s]:", time_end-time_begin
  if(myid_w == 0) print *, "time/count/cell", (time_end-time_begin)/(ngrid1*ngrid2*ngrid3)/nhymax
  
  is_final = .true.
  call Output(forceoutput)

  call FinalizeMPI
  if(myid_w == 0) print *, "program has been finished"
  
end program main


subroutine GenerateGrid
  use basicmod
  use mpimod
  implicit none
  real(8)::dx,dy,dz
  real(8)::x1minloc,x1maxloc
  real(8)::x2minloc,x2maxloc
  real(8)::x3minloc,x3maxloc
  integer::i,j,k
  
  ! x coordinates
      
  x1minloc = x1min + (x1max-x1min)/ntiles(1)* coords(1)
  x1maxloc = x1min + (x1max-x1min)/ntiles(1)*(coords(1)+1)
  
  dx=(x1maxloc-x1minloc)/dble(ngrid1)
  do i=1,in
     x1a(i) = dx*(i-(mgn+1))+x1minloc
  enddo
  do i=1,in-1
     x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
  enddo
 
  ! y coordinates
  x2minloc = x2min + (x2max-x2min)/ntiles(2)* coords(2)
  x2maxloc = x2min + (x2max-x2min)/ntiles(2)*(coords(2)+1)
 
  dy=(x2maxloc-x2minloc)/dble(ngrid2)
  do j=1,jn
     x2a(j) = dy*(j-(mgn+1))+x2minloc
  enddo

  do j=1,jn-1
     x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
  enddo
  
  ! z coordinates
  x3minloc = x3min + (x3max-x3min)/ntiles(3)* coords(3)
  x3maxloc = x3min + (x3max-x3min)/ntiles(3)*(coords(3)+1)
 
  dz=(x3maxloc-x3minloc)/ngrid3
  do k=1,kn
     x3a(k) = dz*(k-(mgn+1))+x3minloc
  enddo
  do k=1,kn-1
     x3b(k) = 0.5d0*(x3a(k+1)+x3a(k))
  enddo

!$acc update device (x1a,x1b)
!$acc update device (x2a,x2b)
!$acc update device (x3a,x3b)
  
  return
end subroutine GenerateGrid

subroutine GenerateProblem
  use basicmod
  use eosmod
  use mpimod
  use boundarymod
  implicit none
  integer::i,j,k

  real(8) :: RU,RD,Gra,P0,Rcenter
  data RU  / 2.0d0 /
  data RD  / 1.0d0 /
  data GRA / 0.1d0 /
  data P0  / 2.5d0 /
  data Rcenter / 0.0d0 /
  
  real(8)::d1d(in),pa1d(in),pb1d(in)

  real(8)::pi
  pi=acos(-1.0d0)

  d(:,:,:) = 1.0d0

  do i=1,in-1
     if(x1b(i) .lt. Rcenter) then
        d1d(i) = RD   
     else
        d1d(i) = RU
     endif
  enddo

  pa1d(1:is) = P0 - GRA * RD * x1a(1:is) ! x1a(1) is negative

  do i=is,in-1
     pa1d(i+1) = pa1d(i) & 
          &       - d1d(i)*GRA*(x1a(i+1)-x1a(i))
  enddo
      

  do i=1,in-1
     pb1d(i) = 0.5d0*(pa1d(i+1)+pa1d(i))
  enddo
  
  do k=ks,ke
  do j=js,je
  do i=is,ie
     d(i,j,k) = d1d(i)
     p(i,j,k) = pb1d(i)
     v1(i,j,k) = 0.0d0
     v2(i,j,k) = 0.0d0
     v3(i,j,k) = 0.0d0
  enddo
  enddo
  enddo

  do k=ks,ke
  do j=js,je
  do i=1,in-1
     gp(i,j,k) = -GRA*x1b(i)
  enddo
  enddo
  enddo

! pert
  do k=ks,ke
  do j=js,je
  do i=is,ie
     v1(i,j,k)= 0.01d0/4.0d0 &
     & *(1.0d0+cos(2.0d0*pi*(x2b(j)-(x2max+x2min)/2.0d0)/(x2max-x2min))) &
     & *(1.0d0+cos(2.0d0*pi*(x1b(i)-(x1max+x1min)/2.0d0)/(x1max-x1min)))
  enddo
  enddo
  enddo

  do k=ks,ke
  do j=js,je
  do i=is,ie
! adiabatic
     ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
     cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
! isotermal
!          ei(i,j,k) = p(i,j,k)
!          cs(i,j,k) = csiso
  enddo
  enddo
  enddo
      
  if(myid_w ==0 )print *,"initial profile is set"
  call BoundaryCondition

!$acc update device (d,v1,v2,v3)
!$acc update device (p,ei,cs)
!$acc update device (b1,b2,b3,bp)
!$acc update device (gp)
      
  return
end subroutine GenerateProblem
