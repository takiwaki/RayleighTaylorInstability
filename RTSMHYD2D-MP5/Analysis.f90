module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1b,x2b
    real(8),dimension(:),allocatable:: x1a,x2a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p,gp
    real(8),dimension(:,:,:),allocatable:: kin
    real(8):: dx,dy
end module fieldmod

program data_analysis
  use fieldmod
  implicit none
  integer:: fbeg, fend
  logical:: flag
  integer,parameter:: unitcon=100

  INQUIRE(FILE ="control.dat",EXIST = flag)
  if(flag) then
     open (unitcon,file="control.dat" &
     &        ,status='old',form='formatted')
     read (unitcon,*) fbeg,fend
     close(unitcon)
  endif

  FILENUMBER: do incr  = fbeg,fend
     write(6,*) "file number",incr
     call ReadData
     call Visualize2D
     call MixingRate
  enddo FILENUMBER

  stop
end program data_analysis

subroutine ReadData
  use fieldmod
  implicit none   
  character(20),parameter::dirname="bindata/"
  character(40)::filename
  integer,parameter::unitinp=13
  integer,parameter::unitbin=14
  character(8)::dummy
  logical,save:: is_inited
  data is_inited / .false. /

  write(filename,'(a3,i5.5,a4)')"unf",incr,".dat"
  filename = trim(dirname)//filename
  open(unitinp,file=filename,status='old',form='formatted')
  read(unitinp,*) dummy,time,dt
  read(unitinp,*) dummy,izone,igs
  read(unitinp,*) dummy,jzone,jgs
  close(unitinp)
  in=izone+2*igs
  jn=jzone+2*jgs
  kn=1

  is=1+igs
  js=1+jgs
  ie=in-igs
  je=jn-jgs

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in))
     allocate( x2b(jn),x2a(jn))
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate( p(in,jn,kn))
     is_inited = .true.
  endif

  write(filename,'(a3,i5.5,a4)')"bin",incr,".dat"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='binary')
  read(unitbin)x1b(:),x1a(:)
  read(unitbin)x2b(:),x2a(:)
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin)  p(:,:,:)
  close(unitbin)
  
  dx = x1b(2)-x1b(1)
  dy = x2b(2)-x2b(1)

  return
end subroutine ReadData

subroutine Visualize2D
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit2D=123

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     allocate( kin(in,jn,kn))
     is_inited = .true.
  endif

  k=1
  do j=js,je
  do i=is,ie 
     kin(i,j,k)= 0.5d0*d(i,j,k)*( &
               & +v1(i,j,k)*v1(i,j,k) &
               & +v2(i,j,k)*v2(i,j,k) &
               & )
  enddo
  enddo

  write(filename,'(a3,i5.5,a4)')"den",incr,".dat"
  filename = trim(dirname)//filename
  open(unit2D,file=filename,status='replace',form='formatted')

  write(unit2D,'(1a,4(1x,E12.3))') "#",time
  write(unit2D,'(1a,4(1x,a8))') "#","1:x    ","2:y     ","3:den ","4:E_kin "
  do j=js,je
  do i=is,ie
     write(unit2D,'(4(1x,E12.3))') x1b(i),x2b(j),d(i,j,k),kin(i,j,k)
  enddo
     write(unit2D,*)
  enddo

  close(unit2D)


  return
end subroutine Visualize2D

subroutine MixingRate
  use fieldmod
  implicit none
  integer::i,j,k
  real(8),dimension(:),allocatable,save:: aved
  real(8)                               :: mixd

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit1D=123
  integer,parameter::unittime=1234

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     allocate(aved(in))
     is_inited = .true.
  endif
  k=1
  do i=is,ie
     aved(i) = 0.0d0
  do j=js,je 
     aved(i) = aved(i) + d(i,j,k)*(x2a(j+1)-x2a(j))/(x2a(je+1)-x2a(js))
  enddo
  enddo

  mixd=0.0d0
  do i=is,ie
     mixd = mixd + (aved(i)-1.0d0)*(2.0d0-aved(i))*(x1a(i+1)-x1a(i))/(x1a(ie+1)-x1a(is))
  enddo

  write(filename,'(a3,i5.5,a4)')"var",incr,".dat"
  filename = trim(dirname)//filename
  open(unit1D,file=filename,status='replace',form='formatted')

  write(unit1D,'(1a,4(1x,E12.3))') "#",time
  write(unit1D,'(1a,4(1x,a8))') "#","1:x    ","2:y     ","3:den "
  do i=is,ie
     write(unit1D,'(4(1x,E12.3))') x1b(i),aved(i)
  enddo
  close(unit1D)

  write(filename,'(a3,i5.5,a4)')"tot",incr,".dat"
  filename = trim(dirname)//filename
  open(unittime,file=filename,status='replace',form='formatted')
  write(unittime,'(4(1x,E12.3))') time,mixd
  close(unittime)

  return
end subroutine MixingRate
