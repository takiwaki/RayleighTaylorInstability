module unitsmod
  implicit none
  real(8),parameter::    pc  = 3.085677581d18   ! parsec in [cm]
  real(8),parameter::    mu  = 1.660539066d-24  ! g
  real(8),parameter:: Msolar = 1.989e33         ! g
  real(8),parameter::   kbol = 1.380649d-23     ! J/K
  real(8),parameter::   year = 365.0d0*24*60*60 ! sec
  
  real(8),parameter:: mp=1.67262192369d-24  !! proton mass [g]
  real(8),parameter:: kB=1.380649d-16  !! Boltzmann constant [erg/K]
  real(8),parameter:: erg_to_keV= 6.242d8  !! erg => keV
end module unitsmod

module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1a,x1b,dvl1a
    real(8),dimension(:),allocatable:: x2a,x2b,dvl2a
    real(8),dimension(:),allocatable:: x3a,x3b,dvl3a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p,ei,gp
    real(8),dimension(:,:,:),allocatable:: b1,b2,b3,bp
    real(8),dimension(:,:,:),allocatable:: Tem,edot
    real(8):: dx
    real(8):: gam,rho0,Eexp
    real(8),dimension(:,:),allocatable:: rshock_ray
    real(8):: rshock,Msw,kTshock,Vshock,Lbol
end module fieldmod

program data_analysis
  use fieldmod
  implicit none
  integer:: fbeg, fend
  logical:: flag
  integer,parameter:: unitcon=100

  INQUIRE(FILE ="control.dat",EXIST = flag)
  if(flag) then
!     open (unitcon,file="control.dat" &
!     &        ,status='old',form='formatted')
!     read (unitcon,*) fbeg,fend
!     close(unitcon)
  endif
  fbeg=1
  fend=30
  FILENUMBER: do incr  = fbeg,fend
     write(6,*) "file index",incr
     call ReadData
     call WriteXDMF_for_Visit(incr, time, in, jn, kn)
!     call Visualize2D
!     call TimeProfle
  enddo FILENUMBER

  stop
end program data_analysis

subroutine ReadData
  use fieldmod
  implicit none   
  character(20),parameter::dirname="../bindata/"
  character(40)::filename
  integer,parameter::unitinp=13
  integer,parameter::unitbin=14
  character(8)::dummy
  logical flag
  logical,save:: is_inited
  data is_inited / .false. /

  write(filename,'(a3,i5.5,a4)')"unf",incr,".dat"
  filename = trim(dirname)//filename

  INQUIRE(FILE =filename,EXIST = flag)
  if(.not. flag) then
     write(6,*) "FILE:",filename
     stop 'Cannot Open  data'
  endif
  open(unitinp,file=filename,status='old',form='formatted')
  read(unitinp,*) dummy,time,dt
  read(unitinp,*) dummy,izone
  read(unitinp,*) dummy,jzone
  read(unitinp,*) dummy,kzone
  close(unitinp)
  in=izone
  jn=jzone
  kn=kzone
!  write(6,*)igs,jgs
  is=1
  js=1
  ks=1
  ie=in
  je=jn
  ke=kn

  if(.not. is_inited)then
     allocate( x1b(in+1),x1a(in+1),dvl1a(in+1))
     allocate( x2b(jn+1),x2a(jn+1),dvl2a(jn+1))
     allocate( x3b(kn+1),x3a(kn+1),dvl3a(kn+1))
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate(b1(in,jn,kn))
     allocate(b2(in,jn,kn))
     allocate(b3(in,jn,kn))
     allocate( p(in,jn,kn))
     allocate(gp(in,jn,kn))
     is_inited = .true.
  endif

  write(filename,'(A)')"g1dDT"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='unformatted',access="stream")
  read(unitbin)x1b(:),x1a(:)
  
  write(filename,'(A)')"g2dDT"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='unformatted',access="stream")
  read(unitbin)x2b(:),x2a(:)
  close(unitbin)
  
  write(filename,'(A)')"g3dDT"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='unformatted',access="stream")
  read(unitbin) x3b(:),x3a(:)
  close(unitbin)
  
  
  write(filename,'((A),i5.5)')"d3dDT.",incr
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='unformatted',access="stream")
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin) b1(:,:,:)
  read(unitbin) b2(:,:,:)
  read(unitbin) b3(:,:,:)
  read(unitbin) bp(:,:,:)
  read(unitbin)  p(:,:,:)
  read(unitbin) gp(:,:,:)
  close(unitbin)
  
  return
end subroutine ReadData

subroutine EstimateEmissivity
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k
  logical:: is_inited
  data is_inited / .false. /
  if(.not. is_inited) then
     allocate(Tem(in,jn,kn))
     allocate(edot(in,jn,kn))
     is_inited = .true.
  endif
  k=ks
  do j=js,je
  do i =is,ie
     ! p = n k T => T = p/(n)/k since kbol[J/K] kbol*1.0d5 [erg/K] 
     Tem(i,j,k) = p(i,j,k)/(d(i,j,k)/mp) / (kbol*1.0d5) ! [K]
     edot(i,j,k) = 1.4d-27 * (d(i,j,k)/mp)**2 *sqrt(Tem(i,j,k)) ! erg/s/cm^3
  enddo
  enddo
end subroutine EstimateEmissivity

subroutine FindShockRadius
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k
  real(8),dimension(:,:),allocatable,save:: pmax_ray,kTshock_ray,Vshock_ray
  logical::is_inited
  data is_inited /.false./

  if(.not. is_inited)then
     allocate(pmax_ray(jn,kn))
     allocate(rshock_ray ,mold=pmax_ray)
     allocate(kTshock_ray,mold=pmax_ray)
     allocate(Vshock_ray ,mold=pmax_ray)
  endif
  
  rshock_ray(:,:) = 0.0d0
  pmax_ray(:,:) = 0.0d0
  kTshock_ray(:,:) = 0.0d0
  Vshock_ray(:,:)  = 0.0d0
  k = ks
  do j=js,je 
     do i=is,ie
        if(pmax_ray(j,k) < p(i,j,k)) then
           pmax_ray(j,k) = p(i,j,k)
           rshock_ray = x1b(i)
           ! p = n T 
           kTshock_ray(j,k) = (kbol*1.0d5)*Tem(i,j,k)*erg_to_keV ! T [keV]
           Vshock_ray(j,k)  = v1(i,j,k)/1.0e5 ! cm/s => km/s
        endif
     enddo
  enddo
  rshock = 0.0d0
  do j=js,je
     if(rshock < rshock_ray(j,k) )then
        rshock  =  rshock_ray(j,k)
        kTshock = kTshock_ray(j,k)
        Vshock  =  Vshock_ray(j,k)
     endif
  enddo
  !print *, "rshock=",rshock/pc,"[pc]"
  !print *, "kTshock=",kTshock,"[keV]"
  !print *, "Vshock=",Vshock,"[km/s]"
  is_inited = .true.
end subroutine FindShockRadius

subroutine Visualize2D
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit2D=432

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
     is_inited = .true.
  endif

  k = int((ks+ke)/2)
! boundary 
  do i=is,ie
      d(i,js-1,k) =   d(i,js,k)
      p(i,js-1,k) =   p(i,js,k)
     v1(i,js-1,k) =  v1(i,js,k)

      d(i,je+1,k) =   d(i,je,k)
      p(i,je+1,k) =   p(i,je,k)
     v1(i,je+1,k) =  v1(i,je,k)
  enddo

  write(filename,'(a6,i5.5,a4)')"xy-pro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit2D,file=filename,status='replace',form='formatted')

  write(unit2D,'(1a,a6,1(1x,E12.3))') "#"," time=",time
!                                    12345678    1234567890123   1234567890123   123456789012
  write(unit2D,'(1a,2(1x,a7,i0))') "#"," Nrad= ",ie-is+1," Nthe= ",je-js+2

  write(unit2D,'(1a,6(1x,a13))') "#","1:r[cm] ","2:theta[rad] ","3:den[1/cm^3] ","4:p[erg/cm3] ","5:vel[cm/s] ","6:edot[cgs]"

  do j=js,je
  do i=is,ie
     write(unit2D,'(1x,SP,7(E13.3,1x))') x1b(i),x2b(j),d(i,j,k),p(i,j,k),v1(i,j,k),v2(i,j,k),v3(i,j,k)
  enddo
     write(unit2D,*)
  enddo

  close(unit2D)


  j = int((js+je)/2)

  write(filename,'(a6,i5.5,a4)')"zx-pro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit2D,file=filename,status='replace',form='formatted')

  write(unit2D,'(1a,a6,1(1x,E12.3))') "#"," time=",time
!                                    12345678    1234567890123   1234567890123   123456789012
  write(unit2D,'(1a,2(1x,a7,i0))') "#"," Nrad= ",ie-is+1," Nthe= ",je-js+2

  write(unit2D,'(1a,6(1x,a13))') "#","1:r[cm] ","2:theta[rad] ","3:den[1/cm^3] ","4:p[erg/cm3] ","5:vel[cm/s] ","6:edot[cgs]"

  do k=ks,ke
  do i=is,ie
     write(unit2D,'(1x,SP,7(E13.3,1x))') x1b(i),x3b(k),d(i,j,k),p(i,j,k),v1(i,j,k),v2(i,j,k),v3(i,j,k)
  enddo
     write(unit2D,*)
  enddo

  close(unit2D)

  
  return
end subroutine Visualize2D

subroutine TimeProfle
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="./"
  character(40)::filename
  integer,save::unittpr
  real(8)::Etot,pi

  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
  endif

  pi = acos(-1.0d0)
  Msw  = 0.0d0
  Etot = 0.0d0
  Lbol = 0.0d0

  Etot=0.0d0
  k=ks
  do j=js,je
  do i=is,ie
     if(x1b(i) <= rshock_ray(j,k) ) Msw  = Msw  + d(i,j,k)*dvl1a(i)*4.0d0*pi
     Lbol = Lbol + edot(i,j,k) * dvl1a(i)*4.0d0*pi
     Etot = Etot + (0.5d0*d(i,j,k)*(v1(i,j,k)**2+v2(i,j,k)**2)+ei(i,j,k))*dvl1a(i)*dvl2a(j)*2.0d0*pi
  enddo
  enddo

  if(.not. is_inited) then
     write(filename,'(A)')"t-prof.dat"
     filename = trim(dirname)//filename
     print *,"time evolution is written in", filename
     open(newunit=unittpr,file=filename,status='replace',form='formatted')
     write(unittpr,'(1a,1x,A)') "#"," time[year] rshock[pc] Msw[Ms] Tshock[keV] Vshock[km/s] Lbol[erg/s] Etot[erg]"
  endif
!  write(unittot,'(1a,4(1x,E12.3))') "#",time/year
!                                    12345678   1234567890123     1234567890123   123456789012
!  write(unittot,'(1a,4(1x,a13))') "#","1:r[pc] ","2:den[1/cm^3] ","3:p[erg/cm3] ","4:vel[km/s] "

  write(unittpr,'(1x,7(1x,E13.3))') time/year,rshock/pc,Msw/Msolar,kTshock,Vshock,Lbol,Etot
  ! close(unittot)

  is_inited = .true.
  return
end subroutine  TimeProfle

subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  write(*, *) trim(command)
  call system(command)
end subroutine makedirs


subroutine WriteXDMF_for_Visit(incr, time, in, jn, kn)
  use, intrinsic :: iso_fortran_env, only: int64,real64
  implicit none
  integer, intent(in) :: incr, in, jn, kn
  real(real64), intent(in) :: time

  character(*), parameter :: dirname = "../bindata/"
  character(256) :: xmfname, fgridx, fgridy, fgridz, fdata
  integer :: u, ncell
  integer(int64) :: bytes_per_real, bytes_per_field
  integer(int64) :: off_d, off_v1, off_v2, off_v3, off_b1, off_b2, off_b3, off_bp, off_p, off_gp

  ! ---- file names (match your ReadData) ----
  write(xmfname,'(a,i5.5,a)') "dt_", incr, ".xmf"
  xmfname = trim(dirname)//trim(xmfname)

  fgridx = trim(dirname)//"g1dDT"
  fgridy = trim(dirname)//"g2dDT"
  fgridz = trim(dirname)//"g3dDT"
  write(fdata,'(a,i5.5)') trim(dirname)//"d3dDT.", incr

  ! ---- sizes & offsets ----
  ! stream/unformatted wrote raw reals; assume real64 (8 bytes) because iso_fortran_env real64
  bytes_per_real = int(storage_size(0.0_real64)/8, int64)

  ncell = in*jn*kn
  bytes_per_field = int(ncell, int64) * bytes_per_real

  off_d  = 0_int64 * bytes_per_field
  off_v1 = 1_int64 * bytes_per_field
  off_v2 = 2_int64 * bytes_per_field
  off_v3 = 3_int64 * bytes_per_field
  off_b1 = 4_int64 * bytes_per_field
  off_b2 = 5_int64 * bytes_per_field
  off_b3 = 6_int64 * bytes_per_field
  off_bp = 7_int64 * bytes_per_field
  off_p  = 8_int64 * bytes_per_field
  off_gp = 9_int64 * bytes_per_field

  ! ---- write XDMF (XML) ----
  open(newunit=u, file=xmfname, status="replace", action="write", form="formatted")

  write(u,'(a)') '<?xml version="1.0" ?>'
  write(u,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(u,'(a)') '<Xdmf Version="2.0">'
  write(u,'(a)') '  <Domain>'
  write(u,'(a)') '    <Grid Name="Grid" GridType="Uniform">'
  write(u,'(a,es24.16,a)') '      <Time Value="', time, '"/>'

  ! Topology: Rectilinear mesh uses node dimensions = (in+1, jn+1, kn+1)
  ! Use Order="Fortran" consistently so Dimensions list matches (i,j,k).
  write(u,'(a,i0,1x,i0,1x,i0,a)') '      <Topology TopologyType="3DRectMesh" Dimensions="', &
                                   kn+1, jn+1, in+1, '"/>'

  write(u,'(a)') '      <Geometry GeometryType="VXVYVZ">'
  call write_axis(u, fgridx, in+1, int(in+1,int64)*bytes_per_real, bytes_per_real)  ! x1a(:)
  call write_axis(u, fgridy, jn+1, int(jn+1,int64)*bytes_per_real, bytes_per_real)  ! x2a(:)
  call write_axis(u, fgridz, kn+1, int(kn+1,int64)*bytes_per_real, bytes_per_real)  ! x3a(:)

  write(u,'(a)') '      </Geometry>'

  ! ---- Cell-centered attributes (in,jn,kn) ----
  call write_attr(u, "d" , fdata, in, jn, kn, off_d , bytes_per_real)
  call write_attr(u, "v1", fdata, in, jn, kn, off_v1, bytes_per_real)
  call write_attr(u, "v2", fdata, in, jn, kn, off_v2, bytes_per_real)
  call write_attr(u, "v3", fdata, in, jn, kn, off_v3, bytes_per_real)
  call write_attr(u, "b1", fdata, in, jn, kn, off_b1, bytes_per_real)
  call write_attr(u, "b2", fdata, in, jn, kn, off_b2, bytes_per_real)
  call write_attr(u, "b3", fdata, in, jn, kn, off_b3, bytes_per_real)
  call write_attr(u, "bp", fdata, in, jn, kn, off_bp, bytes_per_real)
  call write_attr(u, "p" , fdata, in, jn, kn, off_p , bytes_per_real)
  call write_attr(u, "gp", fdata, in, jn, kn, off_gp, bytes_per_real)

  write(u,'(a)') '    </Grid>'
  write(u,'(a)') '  </Domain>'
  write(u,'(a)') '</Xdmf>'

  close(u)

contains

  subroutine write_axis(u, fname, n, seek_bytes, bpr)
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer, intent(in) :: u, n
    character(*), intent(in) :: fname
    integer(int64), intent(in) :: seek_bytes, bpr

    ! XDMF Binary DataItem:
    ! - Format="Binary"
    ! - Seek="bytes" to skip within file
    ! - Precision="8" etc.
    write(u,'(a)') '        <DataItem Dimensions="'//trim(itoa(n))// &
                   '" NumberType="Float" Precision="'//trim(itoa(int(bpr)))// &
                   '" Format="Binary" Endian="Little" Seek="'//trim(i64toa(seek_bytes))// &
                   '"  >'//trim(fname)//'</DataItem>'
    
!    write(u,'(a)') '        <DataItem Dimensions="'//trim(itoa(n))// &
!                   '" NumberType="Float" Precision="'//trim(itoa(int(bpr)))// &
!                   '" Format="Binary" Endian="Little" Seek="'//trim(i64toa(seek_bytes))// &
!                   '" Order="Fortran">'//trim(fname)//'</DataItem>'
  end subroutine write_axis

  subroutine write_attr(u, name, fname, nx, ny, nz, seek_bytes, bpr)
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer, intent(in) :: u, nx, ny, nz
    character(*), intent(in) :: name, fname
    integer(int64), intent(in) :: seek_bytes, bpr

    write(u,'(a)') '      <Attribute Name="'//trim(name)//'" AttributeType="Scalar" Center="Cell">'
    write(u,'(a)') '        <DataItem Dimensions="'//trim(itoa(nz))//' '//trim(itoa(ny))//' '//trim(itoa(nx))// &
                   '" NumberType="Float" Precision="'//trim(itoa(int(bpr)))// &
                   '" Format="Binary" Endian="Little" Seek="'//trim(i64toa(seek_bytes))// &
                   '" >'//trim(fname)//'</DataItem>'
!                   '" Order="Fortran">'//trim(fname)//'</DataItem>'
    write(u,'(a)') '      </Attribute>'
  end subroutine write_attr

  function itoa(i) result(s)
    implicit none
    integer, intent(in) :: i
    character(32) :: s
    write(s,'(i0)') i
  end function itoa

  function i64toa(i) result(s)
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    integer(int64), intent(in) :: i
    character(64) :: s
    write(s,'(i0)') i
  end function i64toa

end subroutine WriteXDMF_for_Visit
