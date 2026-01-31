module unitsmod
  implicit none
  real(8),parameter::    pc  = 3.085677581d18   ! parsec in [cm]
  real(8),parameter::    mu  = 1.660539066d-24  ! g
  real(8),parameter:: Msolar = 1.989e33         ! g
  real(8),parameter:: Rsolar = 6.96340d10       ! cm
  real(8),parameter::   kbol = 1.380649d-23     ! J/K
  real(8),parameter::   year = 365.0d0*24*60*60 ! sec
  real(8),parameter::    day =         24*60*60 ! sec
  
end module unitsmod

module fieldmod
    implicit none
    integer:: incr
    real(8):: time,dt
    integer:: in,jn,kn
    integer:: izone,jzone,kzone
    integer:: igs,jgs,kgs
    integer:: is,js,ks,ie,je,ke
    real(8),dimension(:),allocatable:: x1b,x2b,dvl1a
    real(8),dimension(:),allocatable:: x1a,x2a,dvl2a
    real(8),dimension(:,:,:),allocatable:: d,v1,v2,v3,p,ei,gp
    integer,parameter::ncomp=4
    real(8),dimension(:,:,:,:),allocatable:: Xcomp

    real(8):: dx
    real(8):: gam,rho0,Eexp

    real(8),dimension(:),allocatable:: Mass
    real(8),dimension(:,:),allocatable:: r_shock
    integer,dimension(:,:),allocatable:: i_shock
    real(8):: r_shock_ave, r_shock_min, r_shock_max

    real(8),dimension(:),allocatable::d1d,p1d,v11d
    real(8),dimension(:,:),allocatable::x1d

    real(8),dimension(:,:),allocatable::d2d,p2d,v12d
    real(8),dimension(:,:,:),allocatable::x2d
    real(8),dimension(:,:),allocatable::deld2d

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
     write(6,*) "file index",incr
     call ReadData
     call Visualize1D
     call DetectShock
     call Visualize2D
     call VelocityDistibution
     call TimeProfile
  enddo FILENUMBER

  stop
end program data_analysis

subroutine ReadData
  use fieldmod
  implicit none
  integer:: n
  character(20),parameter::dirname="../bindata/"
  character(40)::filename
  integer,parameter::unitinp=13
  integer,parameter::unitbin=14
  character(8)::dummy
  logical flag
  real(8),dimension(:,:,:),allocatable,save:: tmp
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
  read(unitinp,*) dummy,izone,igs
  read(unitinp,*) dummy,jzone,jgs
  close(unitinp)
  in=izone+2*igs
  jn=jzone+2*jgs
  kn=1
!  write(6,*)igs,jgs
  is=1+igs
  js=1+jgs
  ks=1
  ie=in-igs
  je=jn-jgs
  ke=1

  if(.not. is_inited)then
     allocate( x1b(in),x1a(in),dvl1a(in))
     allocate( x2b(jn),x2a(jn),dvl2a(jn))
     allocate( d(in,jn,kn))
     allocate(v1(in,jn,kn))
     allocate(v2(in,jn,kn))
     allocate(v3(in,jn,kn))
     allocate( p(in,jn,kn))
     allocate(ei(in,jn,kn))
     allocate(Xcomp(ncomp,in,jn,kn))
     allocate(tmp(in,jn,kn))
     is_inited = .true.
  endif

  write(filename,'(a3,i5.5,a4)')"bin",incr,".dat"
  filename = trim(dirname)//filename
  open(unitbin,file=filename,status='old',form='unformatted',access="stream")
  read(unitbin)x1b(:),x1a(:),dvl1a(:)
  read(unitbin)x2b(:),x2a(:),dvl2a(:)
  read(unitbin)  d(:,:,:)
  read(unitbin) v1(:,:,:)
  read(unitbin) v2(:,:,:)
  read(unitbin) v3(:,:,:)
  read(unitbin)  p(:,:,:)
  read(unitbin) ei(:,:,:)
  do n= 1,ncomp
     read(unitbin) tmp(:,:,:)
     Xcomp(n,:,:,:) = tmp(:,:,:)
  enddo
  close(unitbin)
  
  return
end subroutine ReadData


!=======================================================================
!    \\\\\\\\\\       DATA_ANALYSIS                       //////////
!    //////////       DetectShock                         \\\\\\\\\\
!----------------------------------------------------------------------
subroutine  DetectShock
  use fieldmod
  implicit none
  integer:: i,j,k
  real*8 :: ra
  real*8 :: surface_shock
  real*8 :: dcos,dphi
  real*8,parameter:: tiny=1.0d-10

  logical, save :: is_inited
  data is_inited / .false. /
!
! Output
!
! r_shock(j,k) :          Shock Radius [cm]
! r_shock_ave  : Averaged Shock radius [cm]
! r_shock_min  : Minimum  Shock radius [cm]
! r_shock_max  : Maximum  Shock radius [cm]
!
!ccccccccccccccccccccccccccccccccccccccccccc
! begin
!ccccccccccccccccccccccccccccccccccccccccccc
  if(.not. is_inited)then
     allocate (r_shock(jn,kn))
     allocate (i_shock(jn,kn))
     is_inited = .true.
  endif

!ccccccccccccccccccccccccccccccccccccccccccc
! shock gain
!ccccccccccccccccccccccccccccccccccccccccccc
      do k=ks,ke
      do j=js,je
         i_shock(j,k) = is
         r_shock(j,k) = 0.0d0
         iloop_s: do i=is,ie+1,1
!..   kon150729 shock defined as eq. (72) in Marti & Mueller '96
            MM96: if( min(p(i+1,j,k),     p(i-1,j,k)) &
     &   .le.   1.0d0*abs(p(i+1,j,k) -    p(i-1,j,k)) &! parameter to judge the shocked is employed
     &   .and.           v1(i+1,j,k).lt. v1(i-1,j,k)) then

     ! At the very begining we assume the shock
            i_shock(j,k) = i
            r_shock(j,k) = x1b(i)
            exit iloop_s

            endif MM96
         enddo iloop_s
      enddo
      enddo

      surface_shock = 0.0d0
      r_shock_ave = 0.0d0
      k=ks
      do j=js,je
         dcos = dvl2a(j)
!         write(6,*) "test",dvl2a(j),dvl3a(k)
         if ( r_shock(j,k) .ne. 0.0d0 ) then
            surface_shock   = surface_shock +                dcos
            r_shock_ave     = r_shock_ave   + r_shock(j,k) * dcos
         endif
      enddo
      r_shock_ave = r_shock_ave /(surface_shock + tiny)
      r_shock_min = minval(r_shock(:,:))
      r_shock_max = maxval(r_shock(:,:))

  return
end subroutine DetectShock

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
     allocate( d2d(in,jn))
     allocate( p2d(in,jn))
     allocate(v12d(in,jn))
     allocate( x2d(ncomp,in,jn))

     allocate( deld2d(in,jn))

     is_inited = .true.
  endif

  k = ks
! boundary 
  do i=is,ie
      d(i,js-1,k) =  d(i,js,k)
      p(i,js-1,k) =  p(i,js,k)
     v1(i,js-1,k) = v1(i,js,k)
     Xcomp(1:ncomp,i,js-1,k) = Xcomp(1:ncomp,i,js,k)

      d(i,je+1,k) =  d(i,je,k)
      p(i,je+1,k) =  p(i,je,k)
     v1(i,je+1,k) = v1(i,je,k)
     Xcomp(1:ncomp,i,je+1,k) = Xcomp(1:ncomp,i,js,k)
  enddo

  do j=js,je+1
  do i=is,ie
       d2d(i,j) =  0.5d0*( d(i,j,k)+ d(i,j-1,k))
       p2d(i,j) =  0.5d0*( p(i,j,k)+ p(i,j-1,k))
      v12d(i,j) =  0.5d0*(v1(i,j,k)+v1(i,j-1,k))
       x2d(:,i,j) =  0.5d0*(XComp(:,i,j,k)+Xcomp(:,i,j-1,k))
  enddo
  enddo

  do j=js,je+1
  do i=is,ie
       deld2d(i,j) = (d2d(i,j)-d1d(i))/d1d(i)  
  enddo
  enddo



  write(filename,'(a6,i5.5,a4)')"twopro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit2D,file=filename,status='replace',form='formatted')
!                                     12345678901
  write(unit2D,'(1a11,1(1x,E12.3))') "#  time_sc=",time
  write(unit2D,'(1a11,1(1x,E12.3))') "#  rad_max=",x1a(ie+1)
  write(unit2D,'(1a11,1(1x,E12.3))') "#  r_shock=",r_shock_max
  write(unit2D,'(1a,2(1x,a7,i0))') "#"," Nrad= ",ie-is+1," Nthe= ",je-js+2

!                                     1234567890123   1234567890123   1234567890123    1234567890123   1234567890123   1234567890123
  write(unit2D,'(1a,10(1x,a13))') "#","1:r[Rs]      ","2:theta[rad] ","3:den[g/cm^3] ","4:p[erg/cm3] " &
                               &    ,"5:vel[cm/s]  ","6:dden       " &
                               &    ,"7:X_Ni       ","8:X_CO       ","9:X_He        ","10:X_H        "

  do j=js,je+1
  do i=is,ie
     write(unit2D,'(1x,10(1x,E13.3))') x1b(i),x2a(j),d2d(i,j),p2d(i,j) & 
                                  &   ,v12d(i,j), deld2d(i,j) &
                                  &   ,x2d(1,i,j),x2d(2,i,j),x2d(3,i,j),x2d(4,i,j)
  enddo
     write(unit2D,*)
  enddo

  close(unit2D)

  return
end subroutine Visualize2D

subroutine Visualize1D
  use unitsmod
  use fieldmod
  implicit none
  integer::i,j,k

  character(20),parameter::dirname="output/"
  character(40)::filename
  integer,parameter::unit1D=123

  real(8),save:: pi
  logical,save:: is_inited
  data is_inited / .false. /

  if(.not. is_inited)then
     call makedirs(dirname)
     is_inited = .true.
     allocate( d1d(in))
     allocate( p1d(in))
     allocate(v11d(in))
     allocate(x1d(ncomp,in))

     allocate(Mass(in))

     pi = acos(-1.0d0)
  endif


  d1d(:) = 0.0d0
  p1d(:) = 0.0d0
  v11d(:) = 0.0d0
  x1d(:,:) = 0.0d0
  Mass(:) = 0.0d0

  k=ks
  do i=is,ie
  do j=js,je
      d1d(i) =  d1d(i) +  d(i,j,k)*dvl2a(j)
      p1d(i) =  p1d(i) +  p(i,j,k)*dvl2a(j)
     v11d(i) = v11d(i) + v1(i,j,k)*dvl2a(j)
     x1d(:,i) = x1d(:,i) + d(i,j,k)*Xcomp(:,i,j,k)*dvl2a(j)
     Mass(i+1) = Mass(i+1) + d(i,j,k)*dvl1a(i)*dvl2a(j)*2.d0*pi
  enddo
      d1d(i) =  d1d(i)/sum(dvl2a(:))
      p1d(i) =  p1d(i)/sum(dvl2a(:))
     v11d(i) = v11d(i)/sum(dvl2a(:))
      x1d(:,i) =  x1d(:,i)/sum(dvl2a(:))/d1d(i)
     Mass(i+1) = Mass(i) + Mass(i+1)
  enddo

  write(filename,'(a6,i5.5,a4)')"onepro",incr,".dat"
  filename = trim(dirname)//filename
  open(unit1D,file=filename,status='replace',form='formatted')
!                                    12345678901
  write(unit1D,'(a11,1(1x,E12.3))') "#  time_sc=",time
!                                     1234567890123   1234567890123    1234567890123   1234567890123   1234567890123   1234567890123   1234567890123
  write(unit1D,'(1a,9(1x,a13))') "#","1:r[cm]      ","2:M[Ms]      ","3:den[g/cm^3] ","4:p[erg/cm3] ","5:vel[cm/s]  " &
 &                                  ,"6:X_Ni       ","7:X_CO        ","8:X_He       ","9:X_H        "

  do i=is,ie
     write(unit1D,'(1x,9(1x,E13.3))') x1b(i),Mass(i)/Msolar,d1d(i),p1d(i),v11d(i) &
                                     &,x1d(1,i),x1d(2,i),x1d(3,i),x1d(4,i)
  enddo
  close(unit1D)

  return
end subroutine Visualize1D

subroutine TimeProfile
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

  Etot=0.0d0
  k=ks
  do j=js,je
  do i=is,ie
     Etot = Etot + (0.5d0*d(i,j,k)*(v1(i,j,k)**2+v2(i,j,k)**2)+ei(i,j,k))*dvl1a(i)*dvl2a(j)*2.0d0*pi
  enddo
  enddo

  if(.not. is_inited) then
     write(filename,'(A)')"t-prof.dat"
     filename = trim(dirname)//filename
     print *,"time evolution is written in", filename
     open(newunit=unittpr,file=filename,status='replace',form='formatted')
  endif
  

!  write(unittpr,'(1a,4(1x,E12.3))') "#",time/year
!                                    12345678   1234567890123     1234567890123   123456789012
!  write(unittpr,'(1a,4(1x,a13))') "#","1:r[pc] ","2:den[1/cm^3] ","3:p[erg/cm3] ","4:vel[km/s] "

  write(unittpr,'(1x,4(1x,E13.3))') time/day,Etot
  
  is_inited = .true.

  return
end subroutine TimeProfile

subroutine VelocityDistibution
  implicit none
end subroutine VelocityDistibution

subroutine makedirs(outdir)
  implicit none
  character(len=*), intent(in) :: outdir
  character(len=256) command
  write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
  write(*, *) trim(command)
  call system(command)
end subroutine makedirs


