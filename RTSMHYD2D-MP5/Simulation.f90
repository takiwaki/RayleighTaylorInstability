      module commons
      implicit none
      integer::nhy
      integer,parameter::nhymax=6000000
      real(8)::time,dt
      real(8),parameter:: Coul=0.03d0
      data time / 0.0d0 /
      real(8),parameter:: timemax=15.0d0
      real(8),parameter:: dtout=timemax/500

      integer,parameter::izones=150
      integer,parameter::jzones=50
      integer,parameter::mgn=3
      integer,parameter::in=izones+2*mgn+1 &
     &                  ,jn=jzones+2*mgn+1 &
     &                  ,kn=1
      integer,parameter::is=mgn+1 &
     &                  ,js=mgn+1 &
     &                  ,ks=1
      integer,parameter::ie=izones+mgn &
     &                  ,je=jzones+mgn &
     &                  ,ke=1

      real(8),parameter:: x1min=-0.75d0,x1max=0.75d0
      real(8),parameter:: x2min=-0.25d0,x2max=0.25d0
      real(8),dimension(in)::x1a,x1b
      real(8),dimension(jn)::x2a,x2b
      real(8),dimension(kn)::x3a,x3b

      real(8),dimension(in,jn,kn)::d,et,mv1,mv2,mv3
      real(8),dimension(in,jn,kn)::p,ei,v1,v2,v3,cs
      real(8),dimension(in,jn,kn)::gp,gp1a,gp2a
      end module commons
     
      module eosmod
      implicit none
! adiabatic
      real(8),parameter::gam=1.4d0 !! adiabatic index
! isothermal
!      real(8)::csiso  !! isothemal sound speed
end module eosmod

      module fluxmod
      use commons, only : in,jn,kn
      implicit none
      integer,parameter::nden=1,nve1=2,nve2=3,nve3=4,nene=5,npre=6,ncsp=7
      integer,parameter::nhyd=7
      real(8),dimension(nhyd,in,jn,kn):: svc

      integer,parameter::mudn=1,muvu=2,muvv=3,muvw=4,muet=5  &
     &                  ,mfdn=6,mfvu=7,mfvv=8,mfvw=9,mfet=10 &
     &                  ,mcsp=11,mvel=12,mpre=13
      integer,parameter:: mflx=5,madd=3

      integer,parameter:: mden=1,mrv1=2,mrv2=3,mrv3=4,meto=5 &
     &                          ,mrvu=muvu,mrvv=muvv,mrvw=muvw
      real(8),dimension(mflx,in,jn,kn):: nflux1,nflux2,nflux3
      real(8),dimension(in,jn,kn):: grvsrc1,grvsrc2,grvsrc3

      integer,parameter:: nRK=3 ! third order RK
      real(8),parameter::twothirds =2.0d0/3.0d0
      real(8),parameter::RKfac(nRK) =(/ 1.0d0, 0.25d0, twothirds /)
      real(8),dimension(in,jn,kn,mflx,nRk):: RKstates

      end module fluxmod

      program main
      use commons
      use fluxmod
      implicit none
      integer:: RKstep
      write(6,*) "setup grids and fields"
      call GenerateGrid
      call GenerateProblem
      call ConsvVariable
      write(6,*) "entering main loop"
! main loop
                                  write(6,*)"step","time","dt"
      mloop: do nhy=1,nhymax
         call TimestepControl
         if(mod(nhy,100) .eq. 0 ) write(6,*)nhy,time,dt

         do RKstep=1,nRk
            call BoundaryCondition
            call StateVevtor
            call GravForce
            call RKsave(RKstep,RKstates)
            call NumericalFlux1
            call NumericalFlux2
            call UpdateConsv(RKstep)
            call PrimVariable
         enddo

         time=time+dt
         call Output
         if(time > timemax) exit mloop
      enddo mloop

      write(6,*) "program has been finished"
      end program main

      subroutine GenerateGrid
      use commons
      implicit none
      real(8)::dx,dy
      integer::i,j,k
      dx=(x1max-x1min)/izones
      do i=1,in
         x1a(i) = dx*(i-(mgn+1))+x1min
      enddo
      do i=1,in-1
         x1b(i) = 0.5d0*(x1a(i+1)+x1a(i))
      enddo

      dy=(x2max-x2min)/jzones
      do j=1,jn
         x2a(j) = dy*(j-(mgn+1))+x2min
      enddo
      do j=1,jn-1
         x2b(j) = 0.5d0*(x2a(j+1)+x2a(j))
      enddo

      return
      end subroutine GenerateGrid

      subroutine GenerateProblem
      use commons
      use eosmod
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
          ei(i,j,k) = p(i,j,k)/(gam-1.0d0)
          cs(i,j,k) = sqrt(gam*p(i,j,k)/d(i,j,k))
      enddo
      enddo
      enddo
      

      call BoundaryCondition

      return
      end subroutine GenerateProblem


      subroutine BoundaryCondition
      use commons
      implicit none
      integer::i,j,k

!reflection
      k=ks
      do j=1,jn-1
      do i=1,mgn
           d(is-i,j,k) =   d(is+i-1,j,k)
          ei(is-i,j,k) =  ei(is+i-1,j,k)
          v1(is-i,j,k) = -v1(is+i-1,j,k) ! sign change
          v2(is-i,j,k) =  v2(is+i-1,j,k)
          v3(is-i,j,k) =  v3(is+i-1,j,k)
          gp(is-i,j,k) =  gp(is+i-1,j,k)
      enddo
      enddo

      k=ks
      do j=1,jn-1
      do i=1,mgn
           d(ie+i,j,k) =   d(ie-i+1,j,k)
          ei(ie+i,j,k) =  ei(ie-i+1,j,k)
          v1(ie+i,j,k) = -v1(ie-i+1,j,k) ! sign change
          v2(ie+i,j,k) =  v2(ie-i+1,j,k)
          v3(ie+i,j,k) =  v3(ie-i+1,j,k)
          gp(ie+i,j,k) =  gp(ie-i+1,j,k)
      enddo
      enddo

! periodic
      k=ks
      do i=1,in-1
      do j=1,mgn
           d(i,js-j,k) =  d(i,je-j+1,k)
          ei(i,js-j,k) = ei(i,je-j+1,k)
          v1(i,js-j,k) = v1(i,je-j+1,k)
          v2(i,js-j,k) = v2(i,je-j+1,k)
          v3(i,js-j,k) = v3(i,je-j+1,k)
          gp(i,js-j,k) = gp(i,je-j+1,k)
      enddo
      enddo

! periodic
      k=ks
      do i=1,in-1
      do j=1,mgn
           d(i,je+j,k) =  d(i,js+j-1,k)
          ei(i,je+j,k) = ei(i,js+j-1,k)
          v1(i,je+j,k) = v1(i,js+j-1,k)
          v2(i,je+j,k) = v2(i,js+j-1,k)
          v3(i,je+j,k) = v3(i,js+j-1,k)
          gp(i,je+j,k) = gp(i,js+j-1,k)
      enddo
      enddo

      return
      end subroutine BoundaryCondition

      subroutine ConsvVariable
      use commons
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          et(i,j,k) = 0.5d0*d(i,j,k)*(   &
     &                    +v1(i,j,k)**2  &
     &                    +v2(i,j,k)**2  &
     &                    +v3(i,j,k)**2) &
     &                    +ei(i,j,k)
          mv1(i,j,k) =d(i,j,k)*v1(i,j,k)
          mv2(i,j,k) =d(i,j,k)*v2(i,j,k)
          mv3(i,j,k) =d(i,j,k)*v3(i,j,k)
      enddo
      enddo
      enddo
      
      return
      end subroutine Consvvariable

      subroutine PrimVariable
      use commons
      use eosmod
      implicit none
      integer::i,j,k
      do k=ks,ke
      do j=js,je
      do i=is,ie
          v1(i,j,k) = mv1(i,j,k)/d(i,j,k)
          v2(i,j,k) = mv2(i,j,k)/d(i,j,k)
          v3(i,j,k) = mv3(i,j,k)/d(i,j,k)

          ei(i,j,k) =  et(i,j,k)          &
     &          -0.5d0*d(i,j,k)*(         &
     &                    +v1(i,j,k)**2   &
     &                    +v2(i,j,k)**2   &
     &                    +v3(i,j,k)**2)

! adiabatic
           p(i,j,k) =  ei(i,j,k)*(gam-1.0d0)
          cs(i,j,k) =  sqrt(gam*p(i,j,k)/d(i,j,k))
! isotermal
!           p(i,j,k) =  d(i,j,k)*csiso**2
!          cs(i,j,k) =  csiso
      enddo
      enddo
      enddo

      return
      end subroutine PrimVariable

      subroutine TimestepControl
      use commons
      implicit none
      real(8)::dtl1
      real(8)::dtl2
      real(8)::dtl3
      real(8)::dtlocal
      real(8)::dtmin
      integer::i,j,k
      dtmin=1.0d90
      do k=ks,ke
      do j=js,je
      do i=is,ie
         dtl1 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
         dtl2 =(x2a(j+1)-x2a(j))/(abs(v2(i,j,k)) +cs(i,j,k))
!         dtl3 =(x1a(i+1)-x1a(i))/(abs(v1(i,j,k)) +cs(i,j,k))
         dtlocal = min (dtl1,dtl2)
         if(dtlocal .lt. dtmin) dtmin = dtlocal
      enddo
      enddo
      enddo

      dt = Coul * dtmin
!      write(6,*)"dt",dt
      return
      end subroutine TimestepControl

      subroutine StateVevtor
      use commons
      use fluxmod
      use eosmod
      implicit none
      integer::i,j,k

      k=ks
      do j=1,jn-1
      do i=1,in-1
         svc(nden,i,j,k) =  d(i,j,k)
         svc(nve1,i,j,k) = v1(i,j,k)
         svc(nve2,i,j,k) = v2(i,j,k)
         svc(nve3,i,j,k) = v3(i,j,k)
! adiabatic
         svc(nene,i,j,k) = ei(i,j,k)/d(i,j,k)
         svc(npre,i,j,k) = ei(i,j,k)*(gam-1.0d0)
         svc(ncsp,i,j,k) = sqrt(gam*(gam-1.0d0)*ei(i,j,k)/d(i,j,k))
! isotermal
!         svc(nene,i,j,k) = csiso**2
!         svc(npre,i,j,k) = d(i,j,k)*csiso**2
!         svc(ncsp,i,j,k) = csiso

 ! for output
         p(i,j,k) = svc(npre,i,j,k) 
      enddo
      enddo

      return
      end subroutine StateVevtor

      subroutine minmod(a,b,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))                &
     &                                        ,sign(1.0d0,a(n))*b(n)))
      enddo

      return
      end subroutine minmod


      subroutine vanLeer(dvp,dvm,dv)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::dvp,dvm
      real(8),dimension(nhyd),intent(out)::dv
      integer:: n

      do n=1,nhyd
         if(dvp(n)*dvm(n) .gt. 0.0d0)then
            dv(n) =2.0d0*dvp(n)*dvm(n)/(dvp(n)+dvm(n))
         else
            dv(n) = 0.0d0
         endif

      enddo

      return
      end subroutine vanLeer



      subroutine MClimiter(a,b,c,d)
      use fluxmod, only : nhyd
      implicit none
      real(8),dimension(nhyd),intent(in)::a,b,c
      real(8),dimension(nhyd),intent(out)::d
      integer:: n

      do n=1,nhyd
         d(n) = sign(1.0d0,a(n))*max(0.0d0,min(abs(a(n))         &
     &                                  ,sign(1.0d0,a(n))*b(n)   &
     &                                  ,sign(1.0d0,a(n))*c(n))) 
      enddo

      return
      end subroutine MClimiter

      subroutine WENO(vc,vmh,vph)
        implicit none
        real(8),dimension(5),intent(in) :: vc
! ph : plus half i+1/2
! mh : minus half i-1/2 
        real(8),intent(out) :: vmh,vph

        real(8),dimension(2),parameter::cs1 = (/  13.0/12.0, 1.0/4.0 /)
        real(8),dimension(2),parameter::cs2 = (/  13.0/12.0, 1.0/4.0 /)
        real(8),dimension(2),parameter::cs3 = (/  13.0/12.0, 1.0/4.0 /) 

        real(8),dimension(3),parameter::gamph = (/  1.0/10.0, 6.0/10.0, 3.0/10.0 /) 
        real(8),dimension(3),parameter::gammh = (/  3.0/10.0, 6.0/10.0, 1.0/10.0 /) 

        real(8),dimension(3),parameter::c1ph  = (/  1.0/3.0, -7.0/6.0,  11.0/6.0 /)
        real(8),dimension(3),parameter::c2ph  = (/ -1.0/6.0,  5.0/6.0,   1.0/3.0 /)
        real(8),dimension(3),parameter::c3ph  = (/  1.0/3.0,  5.0/6.0,  -1.0/6.0 /) 
        
        real(8),dimension(3),parameter::c1mh = (/  -1.0/6.0,  5.0/6.0,   1.0/3.0 /)
        real(8),dimension(3),parameter::c2mh = (/   1.0/3.0,  5.0/6.0,  -1.0/6.0 /)
        real(8),dimension(3),parameter::c3mh = (/  11.0/6.0, -7.0/6.0,   1.0/3.0 /)

        real(8),parameter:: eps = 1.0d-6

        real(8)::s1,s2,s3
        real(8)::wt1ph,wt2ph,wt3ph,w1ph,w2ph,w3ph
        real(8)::wt1mh,wt2mh,wt3mh,w1mh,w2mh,w3mh
        real(8)::v1ph,v2ph,v3ph
        real(8)::v1mh,v2mh,v3mh

! Smoothness Indices of the stencils.
!
!        SI1 -> for stensil S1 = { i-2, i-1, i   }
!        SI2 -> for stensil S2 = { i-1, i  , i+1 }
!        SI2 -> for stensil S3 = { i  , i+1, i+2 }

        s1 =    cs1(1)*( vc(1)-2.0*vc(2)+      vc(3) )**2 &
             & +cs1(2)*( vc(1)-4.0*vc(2)+3.0d0*vc(3) )**2

        s2 =    cs2(1)*( vc(2)-2.0*vc(3)+vc(4) )**2 &
             & +cs2(2)*( vc(2)         -vc(4) )**2

        s3 =    cs3(1)*(       vc(3)-2.0*vc(4)+vc(5) )**2 &
             & +cs3(2)*( 3.0d0*vc(3)-4.0*vc(4)+vc(5) )**2

        wt1ph = gamph(1)/(eps+s1)**2
        wt2ph = gamph(2)/(eps+s2)**2
        wt3ph = gamph(3)/(eps+s3)**2

        w1ph =wt1ph/(wt1ph+wt2ph+wt3ph)
        w2ph =wt2ph/(wt1ph+wt2ph+wt3ph)
        w3ph =wt3ph/(wt1ph+wt2ph+wt3ph)

        v1ph = c1ph(1)*vc(1) + c1ph(2)*vc(2) + c1ph(3)*vc(3)
        v2ph = c2ph(1)*vc(2) + c2ph(2)*vc(3) + c2ph(3)*vc(4)
        v3ph = c3ph(1)*vc(3) + c3ph(2)*vc(4) + c3ph(3)*vc(5)
        vph =  w1ph*v1ph + w2ph*v2ph + w3ph*v3ph

! 
        wt1mh = gammh(1)/(eps+s1)**2
        wt2mh = gammh(2)/(eps+s2)**2
        wt3mh = gammh(3)/(eps+s3)**2

        w1mh =wt1mh/(wt1mh+wt2mh+wt3mh)
        w2mh =wt2mh/(wt1mh+wt2mh+wt3mh)
        w3mh =wt3mh/(wt1mh+wt2mh+wt3mh)

        v1mh = c1mh(1)*vc(1) + c1mh(2)*vc(2) + c1mh(3)*vc(3)
        v2mh = c2mh(1)*vc(2) + c2mh(2)*vc(3) + c2mh(3)*vc(4)
        v3mh = c3mh(1)*vc(3) + c3mh(2)*vc(4) + c3mh(3)*vc(5)
        vmh =  w1mh*v1mh + w2mh*v2mh + w3mh*v3mh
        return
      end subroutine WENO

      subroutine MP5(vc,vmh,vph)
        implicit none
        real(8),dimension(5),intent(in) :: vc
! ph : plus half i+1/2
! mh : minus half i-1/2 
        real(8),intent(out) :: vmh,vph
        real(8),dimension(5),parameter:: ccph = (/  1.0d0/30.0d0, -13.0d0/60.0d0, 47.0d0/60.0d0,   9.0d0/20.0d0,  -1.0d0/20.0d0 /)
        real(8),dimension(5),parameter:: ccpm = (/ -1.0d0/20.0d0,   9.0d0/20.0d0, 47.0d0/60.0d0, -13.0d0/60.0d0,   1.0d0/30.0d0 /)
        real(8) :: vor,vol,vmp,vul,vmd,vlc,vmin,vmax
        real(8) :: dvp,dvm
        real(8) :: djm1,dj,djp1,dm4jph,dm4jmh
        real(8),parameter :: Alpha = 4.0d0, BC2 = 4.0d0/3d0
        real(8),parameter :: eps=1.0d-10

! interpolation for right side => 
               vor = sum(ccph(:)*vc(:))
               dvp = (vc(4)-vc(3))
               dvm = (vc(3)-vc(2))
               
               vmp = vc(3)+minmod2(Alpha*dvm,dvp)
               if( (vor-vc(3))*(vor-vmp) .le. eps)then
                  vph = vor
               else

                  vul = vc(3) + Alpha*dvm

                  djm1 = vc(1) -2d0*vc(2) + vc(3)
                  dj   = vc(2) -2d0*vc(3) + vc(4)
                  djp1 = vc(3) -2d0*vc(4) + vc(5)           

                  dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                  dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)

                  vmd = 0.5d0*(vc(3)+vc(4) - dm4jph)
                  vlc = vc(3) + 0.5d0*(vc(3)-vc(2)) + BC2*dm4jmh

                  vmin = max(min(vc(3),vc(4),vmd),min(vc(3),vul,vlc))
                  vmax = min(max(vc(3),vc(4),vmd),max(vc(3),vul,vlc))

                  vph = vor + minmod2((vmin-vor),(vmax-vor))
               endif

! interpolation for left side  <= 
               vol = sum(ccpm(:)*vc(:))

               dvm = (vc(2)-vc(3))
               dvp = (vc(3)-vc(4))

               vmp = vc(3) + minmod2(Alpha*dvp,dvm)
               if( (vol-vc(3))*(vol-vmp) .le. eps)then
!vmh:  v_i-1/2
                  vmh = vol
               else
                  vul = vc(3) + Alpha*dvp

                  djm1 = vc(1) -2d0*vc(2) + vc(3)
                  dj   = vc(2) -2d0*vc(3) + vc(4)
                  djp1 = vc(3) -2d0*vc(4) + vc(5)           

                  dm4jph = minmod4(4.0d0*dj-djp1,4.0d0*djp1-dj,dj,djp1)
                  dm4jmh = minmod4(4.0d0*dj-djm1,4.0d0*djm1-dj,dj,djm1)

                  vmd = 0.5d0*(vc(3)+vc(2) - dm4jmh)
                  vlc = vc(3) + 0.5d0*(vc(3)-vc(4)) + BC2*dm4jph
           
                  vmin = max(min(vc(3),vc(2),vmd),min(vc(3),vul,vlc))
                  vmax = min(max(vc(3),vc(2),vmd),max(vc(3),vul,vlc))
!vmh:  v_i-1/2
                  vmh = vol + minmod2((vmin-vol),(vmax-vol))
               endif

        return

      contains

      function minmod2(x,y)
      real(8), intent(in) :: x,y
      real(8) :: minmod2

      minmod2 = 0.5d0*(sign(1.0d0,x)+sign(1.0d0,y))*min(dabs(x),dabs(y))
      
      end function minmod2

      function minmod4(d1,d2,d3,d4)
      real(8), intent(in) :: d1,d2,d3,d4
      real(8) :: sign1,minmod4
      
      sign1 = sign(1d0,d1)
      minmod4 = 0.125d0*(sign1+sign(1.0d0,d2)) &
     & *dabs( (sign1+sign(1.0d0,d3)) &
     & *(sign1+sign(1.0d0,d4))) &
     & *min(dabs(d1),dabs(d2),dabs(d3),dabs(d4))
      
      end function minmod4

      end subroutine MP5

      subroutine NumericalFlux1
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn
      use fluxmod
      implicit none
      integer::i,j,k,n
      real(8),dimension(5):: vc
      real(8):: vph,vmh
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

      k=ks
      do j=js,je
      do i=is-1,ie+1
         do n=1,nhyd
            vc(1:5) = svc(n,i-2:i+2,j,k)
!            call MP5(vc,vmh,vph)
            call WENO(vc,vmh,vph)

            leftpr(n,i+1,j,k) = vph
            rigtpr(n,i  ,j,k) = vmh
!            leftpr(n,i+1,j,k) = svc(n,i,j,k)
!            rigtpr(n,i  ,j,k) = svc(n,i,j,k)
!            write(6,*)"wwc",wwc
!            write(6,*)"l,r",toleft,torigt
         enddo
      enddo
!            stop
      enddo

      do j=js,je
      do i=is,ie+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k) ! rho
         leftco(muvu,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)  ! rho v_x
         leftco(muvv,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k)  ! rho v_y
         leftco(muvw,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)  ! rho v_z
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &! e_i+ rho v^2/2
     &               +0.5d0*leftpr(nden,i,j,k)*(    &
     &                     +leftpr(nve1,i,j,k)**2   &
     &                     +leftpr(nve2,i,j,k)**2   &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve1,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve1,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve1,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k)  &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2  &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k)     &
     &                       )                                  *leftpr(nve1,i,j,k) 

         leftco(mcsp,i,j,k)= leftpr(ncsp,i,j,k)
         leftco(mvel,i,j,k)= leftpr(nve1,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(  &
     &                     +rigtpr(nve1,i,j,k)**2 &
     &                     +rigtpr(nve2,i,j,k)**2 &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve1,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve1,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve1,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2) &
     &                     +rigtpr(npre,i,j,k)     &
     &                      )                                    *rigtpr(nve1,i,j,k)

         rigtco(mcsp,i,j,k)= rigtpr(ncsp,i,j,k)
         rigtco(mvel,i,j,k)= rigtpr(nve1,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo

      do j=js,je
      do i=is,ie+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux1(mden,i,j,k)=nflux(mden)
         nflux1(mrv1,i,j,k)=nflux(mrvu)
         nflux1(mrv2,i,j,k)=nflux(mrvv)
         nflux1(mrv3,i,j,k)=nflux(mrvw)
         nflux1(meto,i,j,k)=nflux(meto)
      enddo
      enddo

      return
      end subroutine Numericalflux1

      subroutine NumericalFlux2
      use commons, only: is,ie,in,js,je,jn,ks,ke,kn
      use fluxmod
      implicit none
      integer::i,j,k,n
      real(8),dimension(5):: vc
      real(8):: vmh,vph
      real(8),dimension(nhyd,in,jn,kn):: leftpr,rigtpr
      real(8),dimension(2*mflx+madd,in,jn,kn):: leftco,rigtco
      real(8),dimension(2*mflx+madd):: leftst,rigtst
      real(8),dimension(mflx):: nflux

      k=ks
      do i=is,ie
      do j=js-1,je+1
         do n=1,nhyd
            vc(1:5) = svc(n,i,j-2:j+2,k)
            call MP5(vc,vmh,vph)
            call WENO(vc,vmh,vph)

            leftpr(n,i,j+1,k) = vph
            rigtpr(n,i,j  ,k) = vmh
!            leftpr(n,i+1,j,k) = svc(n,i,j,k)
!            rigtpr(n,i  ,j,k) = svc(n,i,j,k)
!            write(6,*)"wwc",wwc
!            write(6,*)"l,r",toleft,torigt
         enddo
       enddo
       enddo  

      do i=is,ie
      do j=js,je+1
         leftco(mudn,i,j,k)=leftpr(nden,i,j,k)
         leftco(muvw,i,j,k)=leftpr(nve1,i,j,k)*leftpr(nden,i,j,k)
         leftco(muvu,i,j,k)=leftpr(nve2,i,j,k)*leftpr(nden,i,j,k) ! rho v
         leftco(muvv,i,j,k)=leftpr(nve3,i,j,k)*leftpr(nden,i,j,k)
         leftco(muet,i,j,k)=leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(  &
     &                     +leftpr(nve1,i,j,k)**2 &
     &                     +leftpr(nve2,i,j,k)**2 &
     &                     +leftpr(nve3,i,j,k)**2)

         leftco(mfdn,i,j,k)=leftpr(nden,i,j,k)                   *leftpr(nve2,i,j,k) ! rho v
         leftco(mfvw,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve1,i,j,k)*leftpr(nve2,i,j,k)
         leftco(mfvu,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve2,i,j,k)*leftpr(nve2,i,j,k) &
     &                     +leftpr(npre,i,j,k)
         leftco(mfvv,i,j,k)=leftpr(nden,i,j,k)*leftpr(nve3,i,j,k)*leftpr(nve2,i,j,k)
         leftco(mfet,i,j,k)=(leftpr(nene,i,j,k)*leftpr(nden,i,j,k) &
     &               +0.5d0*leftpr(nden,i,j,k)*(   &
     &                     +leftpr(nve1,i,j,k)**2  &
     &                     +leftpr(nve2,i,j,k)**2  &
     &                     +leftpr(nve3,i,j,k)**2) &
     &                     +leftpr(npre,i,j,k)     &
     &                                       )*leftpr(nve2,i,j,k)

         leftco(mcsp,i,j,k)= leftpr(ncsp,i,j,k)
         leftco(mvel,i,j,k)= leftpr(nve2,i,j,k)
         leftco(mpre,i,j,k)= leftpr(npre,i,j,k)


         rigtco(mudn,i,j,k)=rigtpr(nden,i,j,k)
         rigtco(muvw,i,j,k)=rigtpr(nve1,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvu,i,j,k)=rigtpr(nve2,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muvv,i,j,k)=rigtpr(nve3,i,j,k)*rigtpr(nden,i,j,k)
         rigtco(muet,i,j,k)=rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2)

         rigtco(mfdn,i,j,k)=rigtpr(nden,i,j,k)                   *rigtpr(nve2,i,j,k)
         rigtco(mfvw,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve1,i,j,k)*rigtpr(nve2,i,j,k)
         rigtco(mfvu,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve2,i,j,k)*rigtpr(nve2,i,j,k) &
     &                     +rigtpr(npre,i,j,k)
         rigtco(mfvv,i,j,k)=rigtpr(nden,i,j,k)*rigtpr(nve3,i,j,k)*rigtpr(nve2,i,j,k)
         rigtco(mfet,i,j,k)=(rigtpr(nene,i,j,k)*rigtpr(nden,i,j,k) &
     &               +0.5d0*rigtpr(nden,i,j,k)*(   &
     &                     +rigtpr(nve1,i,j,k)**2  &
     &                     +rigtpr(nve2,i,j,k)**2  &
     &                     +rigtpr(nve3,i,j,k)**2) &
     &                     +rigtpr(npre,i,j,k)     &
     &                                       )*rigtpr(nve2,i,j,k)

         rigtco(mcsp,i,j,k)= rigtpr(ncsp,i,j,k)
         rigtco(mvel,i,j,k)= rigtpr(nve2,i,j,k)
         rigtco(mpre,i,j,k)= rigtpr(npre,i,j,k)

      enddo
      enddo

      do i=is,ie
      do j=js,je+1
         leftst(:)=leftco(:,i,j,k)
         rigtst(:)=rigtco(:,i,j,k)
!         call HLLE(leftst,rigtst,nflux)
         call HLLC(leftst,rigtst,nflux)
         nflux2(mden,i,j,k)=nflux(mden)
         nflux2(mrv1,i,j,k)=nflux(mrvw)
         nflux2(mrv2,i,j,k)=nflux(mrvu) ! mrv2=3, mrvu=2
         nflux2(mrv3,i,j,k)=nflux(mrvv)
         nflux2(meto,i,j,k)=nflux(meto)
      enddo
      enddo

      return
      end subroutine Numericalflux2

      subroutine HLLE(leftst,rigtst,nflux)
      use fluxmod
      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux
      real(8),dimension(mflx)::ul,ur,fl,fr
      real(8)::csl,csr
      real(8):: vl, vr
      real(8):: sl, sr

      ul(1:mflx) = leftst(1:mflx)
      fl(1:mflx) = leftst(mflx+1:2*mflx)
      ur(1:mflx) = rigtst(1:mflx)
      fr(1:mflx) = rigtst(mflx+1:2*mflx)
      csl=leftst(mcsp)
      csr=rigtst(mcsp)
       vl=leftst(mvel)
       vr=rigtst(mvel)

       sl = min(vl,vr) - max(csl,csr)
       sl = min(0.0d0,sl)
       sr = max(vl,vr) + max(csl,csr)
       sr = max(0.0d0,sr)

       nflux(:) = (sr*fl(:)-sl*fr(:) +sl*sr*(ur(:)-ul(:)))/(sr-sl)

      return
      end subroutine HLLE

      subroutine HLLC(leftst,rigtst,nflux)
!=====================================================================
!
! HLLC Scheme
!
! Purpose
! Calculation of Numerical Flux by HLLC method
!
! Reference
!  Toro EF, Spruce M, Speares W. (1992,1994)
!
! Input
! Output
!=====================================================================
      use fluxmod, only: mflx,madd                 &
     &                 , mudn,muvu,muvv,muvw,muet  &
     &                 , mfdn,mfvu,mfvv,mfvw,mfet  &
     &                 , mcsp,mvel,mpre            &
     &                 , mden,mrvu,mrvv,mrvw,meto

      implicit none
      real(8),dimension(2*mflx+madd),intent(in)::leftst,rigtst
      real(8),dimension(mflx),intent(out)::nflux

!----- U -----
! qql :: left state
! qqr :: right state
      real(8) :: rol,vxl,vyl,vzl,ptl,eel
      real(8) :: ror,vxr,vyr,vzr,ptr,eer
      real(8) :: rxl,ryl,rzl
      real(8) :: rxr,ryr,rzr
      real(8) :: ptst

!----- U* ----
! qqlst ::  left state
! qqrst :: right state
      real(8) :: rolst,vxlst,vylst,vzlst,eelst
      real(8) :: rorst,vxrst,vyrst,vzrst,eerst
      real(8) :: rxlst,rylst,rzlst
      real(8) :: rxrst,ryrst,rzrst

!----- flux ---
! fqql ::  left physical flux
! fqqr :: right physical flux
      real(8) :: frol,frxl,fryl,frzl,feel
      real(8) :: fror,frxr,fryr,frzr,feer

!----- wave speed ---
! sl ::  left-going fastest signal velocity
! sr :: right-going fastest signal velocity
! sm :: contact discontinuity velocity
! slst ::  left-going alfven velocity
! srst :: right-going alfven velocity
      real(8) :: sm,sl,sr

! cfl :: left-state Fast wave velocity
! cfr :: right-sate Fast wave velocity
      real(8) :: cfl,cfr

!--------------------
! temporary variables
      real(8) :: sdl,sdr,sdml,sdmr,isdml,isdmr,rosdl,rosdr
      real(8) :: temp
  
! no if
      real(8) :: sign1,maxs1,mins1
      real(8) :: msl,msr

!----- Step 0. ----------------------------------------------------------|

!---- Left state
        
        rol = leftst(mudn)
        eel = leftst(muet)
        rxl = leftst(muvu)
        ryl = leftst(muvv)
        rzl = leftst(muvw)
        vxl = leftst(muvu)/leftst(mudn)
        vyl = leftst(muvv)/leftst(mudn)
        vzl = leftst(muvw)/leftst(mudn)
        ptl = leftst(mpre)

!---- Right state
        
        ror = rigtst(mudn)
        eer = rigtst(muet)
        rxr = rigtst(muvu)
        ryr = rigtst(muvv)
        rzr = rigtst(muvw)
        vxr = rigtst(muvu)/rigtst(mudn)
        vyr = rigtst(muvv)/rigtst(mudn)
        vzr = rigtst(muvw)/rigtst(mudn)
        ptr = rigtst(mpre)
!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!
         
        cfl = leftst(mcsp)
        cfr = rigtst(mcsp)

        sl = min(vxl,vxr)-max(cfl,cfr) ! note sl is negative
        sr = max(vxl,vxr)+max(cfl,cfr)
!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxs
!
! Left value
        frol = leftst(mfdn)
        feel = leftst(mfet)
        frxl = leftst(mfvu)
        fryl = leftst(mfvv)
        frzl = leftst(mfvw)

! Right value
! Left value
        fror = rigtst(mfdn)
        feer = rigtst(mfet)
        frxr = rigtst(mfvu)
        fryr = rigtst(mfvv) 
        frzr = rigtst(mfvw)

!----- Step 4. ----------------------------------------------------------|
! compute middle and alfven wave
!
        sdl = sl - vxl
        sdr = sr - vxr
        rosdl = rol*sdl
        rosdr = ror*sdr

        temp = 1.0d0/(rosdr - rosdl)
! Eq. 45
        sm = (rosdr*vxr - rosdl*vxl - ptr + ptl)*temp
           
        sdml = sl - sm; isdml = 1.0d0/sdml
        sdmr = sr - sm; isdmr = 1.0d0/sdmr
        
!----- Step 5. ----------------------------------------------------------|
! compute intermediate states
!
! Eq. 49
        ptst = (rosdr*ptl-rosdl*ptr+rosdl*rosdr*(vxr-vxl))*temp

!----- Step 5A. ----------------------------------------------------------|
! compute Ul*
!

        rolst = rol*sdl   *isdml
        vxlst = sm
        rxlst = rolst*vxlst
           
        vylst = vyl
        rylst = rolst*vylst
        vzlst = vzl
        rzlst = rolst*vzlst

        eelst =(sdl*eel - ptl*vxl + ptst*sm  )*isdml

!----- Step 5B. ----------------------------------------------------------|
! compute Ur*
!

        rorst   = rosdr   *isdmr
        vxrst = sm
        rxrst = rorst*vxrst
        vyrst = vyr
        ryrst = rorst*vyrst
        vzrst = vzr
        rzrst = rorst*vzrst
           
        eerst = (sdr*eer - ptr*vxr  + ptst*sm  )*isdmr
              
!----- Step 6. ----------------------------------------------------------|
! compute flux
        sign1 = sign(1.0d0,sm)    ! 1 for sm>0, -1 for sm<0
        maxs1 =  max(0.0d0,sign1) ! 1 sm>0, 0 for sm<0
        mins1 = -min(0.0d0,sign1) ! 0 sm>0,-1 for sm<0

        msl   = min(sl  ,0.0d0)   ! 0 for sl > 0, sl for sl < 0
        msr   = max(sr  ,0.0d0)   ! S_R > 0

        nflux(mden) = (frol+msl*(rolst-rol))*maxs1 &
     &               +(fror+msr*(rorst-ror))*mins1
        nflux(meto) = (feel+msl*(eelst-eel))*maxs1 &
     &               +(feer+msr*(eerst-eer))*mins1
        nflux(mrvu) = (frxl+msl*(rxlst-rxl))*maxs1 &
     &               +(frxr+msr*(rxrst-rxr))*mins1
        nflux(mrvv) = (fryl+msl*(rylst-ryl))*maxs1 &
     &               +(fryr+msr*(ryrst-ryr))*mins1
        nflux(mrvw) = (frzl+msl*(rzlst-rzl))*maxs1 &
     &               +(frzr+msr*(rzrst-rzr))*mins1

      return
      end subroutine HLLC

      subroutine GravForce
      use commons
      use fluxmod
      implicit none
      integer :: i,j,k,n

      do k=ks,ke
      do j=js,je
      do i=is,ie+1
         gp1a(i  ,j,k) = gp(i,j,k) &
     & - 0.5d0*(gp(i  ,j,k)-gp(i-1,j,k))

         gp1a(i+1,j,k) = gp(i,j,k) &
     & + 0.5d0*(gp(i+1,j,k)-gp(i  ,j,k))

       grvsrc1(i,j,k) = (gp1a(i+1,j,k)-gp1a(i,j,k))/(x1a(i+1)-x1a(i))*d(i,j,k)

      enddo
      enddo
      enddo

      do k=ks,ke
      do i=is,ie
      do j=js,je+1
         gp2a(i  ,j,k) = gp(i,j,k) &
     & - 0.5d0*(gp(i  ,j,k)-gp(i,j-1,k))

         gp2a(i,j+1,k) = gp(i,j,k) &
     & + 0.5d0*(gp(i,j+1,k)-gp(i  ,j,k))

       grvsrc2(i,j,k) = (gp2a(i,j+1,k)-gp2a(i,j,k))/(x2a(j+1)-x2a(j))*d(i,j,k)

      enddo
      enddo
      enddo

       grvsrc3(:,:,:) = 0.0d0

      return
      end subroutine  GravForce

      subroutine UpdateConsv(RKstep)
      use commons
      use fluxmod
      implicit none
      integer,intent(in):: RKstep
      integer::i,j,k

      do k=ks,ke
      do j=js,je
      do i=is,ie
        
         d(i,j,k) = d(i,j,k)* RKfac(RKstep)           &
     & + RKstates(i,j,k,mden,1)*(1.0d0-RKfac(RKstep)) &
     & +dt*(                                       &
     & +(- nflux1(mden,i+1,j,k)                    &
     &   + nflux1(mden,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mden,i,j+1,k)                    &
     &   + nflux2(mden,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )*RKfac(RKstep)

         mv1(i,j,k) = mv1(i,j,k)*RKfac(RKstep)        &
     & + RKstates(i,j,k,mrv1,1)*(1.0d0-RKfac(RKstep)) &
     & +dt*(                                       &
     &      +  grvsrc1(i,j,k)                      &
     & +(- nflux1(mrv1,i+1,j,k)                    &
     &   + nflux1(mrv1,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mrv1,i,j+1,k)                    &
     &   + nflux2(mrv1,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )*RKfac(RKstep)

         mv2(i,j,k) = mv2(i,j,k)*RKfac(RKstep)        &
     & + RKstates(i,j,k,mrv2,1)*(1.0d0-RKfac(RKstep)) &
     & +dt*(                                       &
     &      +  grvsrc2(i,j,k)                      &
     & +(- nflux1(mrv2,i+1,j,k)                    &
     &   + nflux1(mrv2,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mrv2,i,j+1,k)                    &
     &   + nflux2(mrv2,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )*RKfac(RKstep)

         mv3(i,j,k) = mv3(i,j,k)*RKfac(RKstep)        & 
     & + RKstates(i,j,k,mrv3,1)*(1.0d0-RKfac(RKstep)) &
     & +dt*(                                       &
     &      +  grvsrc3(i,j,k)                      &
     & +(- nflux1(mrv3,i+1,j,k)                    &
     &   + nflux1(mrv3,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(mrv3,i,j+1,k)                    &
     &   + nflux2(mrv3,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )*RKfac(RKstep)

          et(i,j,k) = et(i,j,k)*RKfac(RKstep)         & 
     & + RKstates(i,j,k,meto,1)*(1.0d0-RKfac(RKstep)) &
     & +dt*(                                       &
     &  (- nflux1(meto,i+1,j,k)                    &
     &   + nflux1(meto,i  ,j,k))/(x1a(i+1)-x1a(i)) &
     & +(- nflux2(meto,i,j+1,k)                    &
     &   + nflux2(meto,i,j  ,k))/(x2a(j+1)-x2a(j)) &
     &      )*RKfac(RKstep)
      enddo
      enddo
      enddo

      return
      end subroutine UpdateConsv

      subroutine Output
      use commons
      implicit none
      integer::i,j,k
      character(20),parameter::dirname="bindata/"
      character(40)::filename
      real(8),save::tout
      data tout / 0.0d0 /
      integer::nout
      data nout / 1 /
      integer,parameter:: unitout=17
      integer,parameter:: unitbin=13
      integer,parameter:: gs=1
      integer,parameter:: nvar=6
      real(8)::x1out(is-gs:ie+gs,2)
      real(8)::x2out(js-gs:je+gs,2)
      real(8)::hydout(is-gs:ie+gs,js-gs:je+gs,ks,nvar)

      logical, save:: is_inited
      data is_inited /.false./

      if (.not. is_inited) then
         call makedirs("bindata")
         is_inited =.true.
      endif

      if(time .lt. tout+dtout) return

      write(filename,'(a3,i5.5,a4)')"unf",nout,".dat"
      filename = trim(dirname)//filename

      open(unitout,file=filename,status='replace',form='formatted')
      write(unitout,*) "# ",time,dt
      write(unitout,*) "# ",izones,gs
      write(unitout,*) "# ",jzones,gs
      close(unitout)

      x1out(is-gs:ie+gs,1) = x1b(is-gs:ie+gs)
      x1out(is-gs:ie+gs,2) = x1a(is-gs:ie+gs)

      x2out(js-gs:je+gs,1) = x2b(js-gs:je+gs)
      x2out(js-gs:je+gs,2) = x2a(js-gs:je+gs)

      hydout(is-gs:ie+gs,js-gs:je+gs,ks,1) =  d(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,2) = v1(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,3) = v2(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,4) = v3(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,5) =  p(is-gs:ie+gs,js-gs:je+gs,ks)
      hydout(is-gs:ie+gs,js-gs:je+gs,ks,6) = gp(is-gs:ie+gs,js-gs:je+gs,ks)

      write(filename,'(a3,i5.5,a4)')"bin",nout,".dat"
      filename = trim(dirname)//filename
      open(unitbin,file=filename,status='replace',form='binary') 
      write(unitbin) x1out(:,:)
      write(unitbin) x2out(:,:)
      write(unitbin) hydout(:,:,:,:)
      close(unitbin)

      write(6,*) "output:",nout,time

      nout=nout+1
      tout=time

      return
      end subroutine Output

      subroutine makedirs(outdir)
      implicit none
      character(len=*), intent(in) :: outdir
      character(len=256) command
      write(command, *) 'if [ ! -d ', trim(outdir), ' ]; then mkdir -p ', trim(outdir), '; fi'
      write(*, *) trim(command)
      call system(command)
      end subroutine makedirs

      subroutine RKsave(RKstep)
      use fluxmod
      use commons
      implicit none
      integer :: i,j,k,n
      integer,intent(in):: RKstep
      
      RKstates(:,:,:,mden,RKstep) =    d(:,:,:)
      RKstates(:,:,:,meto,RKstep) =   et(:,:,:)
      RKstates(:,:,:,mrv1,RKstep) =  mv1(:,:,:)
      RKstates(:,:,:,mrv2,RKstep) =  mv2(:,:,:)
      RKstates(:,:,:,mrv3,RKstep) =  mv3(:,:,:)

      return
      end subroutine RKsave
