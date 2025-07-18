! $Id: ana_grid.F 1620 2015-01-08 10:47:13Z marchesiello $
!
!======================================================================
! CROCO is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
! 
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
#include "cppdefs.h"
#ifdef ANA_GRID
      subroutine ana_grid (tile)
      implicit none
# include "param.h"
      integer tile, trd
C$    integer omp_get_thread_num
# include "compute_tile_bounds.h"
      call ana_grid_tile (Istr,Iend,Jstr,Jend)
      return
      end

      subroutine ana_grid_tile (Istr,Iend,Jstr,Jend)
!
!  Set up model grid using analytical expressions:
!-------------------------------------------------------------------
!
!  INPUT: Grid configuration parameters
!  ------------------------------------
! Length_XI      are the physical dimensions of the computational
! Length_ETA     domain [usually measured in  meters];
!
! depth          is the maximum depth [meters, positive];
! f0,beta        are coriolis parameters which set up a beta-plane
!                                                  [1/s, 1/(m*s)].
!
!  OUTPUT: stored in common blocks, see files "scalars.h" "grid.h"
!  ---------------------------------------------------------------
!  xl,el      Physical dimensions of the computational domain 
!                                    [usually measured in  meters];  
!  h          Model bathymetry [meters, positive] at RHO-points.
!  hmin,hmax  Minimum and maximum values of depth of bathymetry [m].
!  f          Coriolis parameter (1/seconds) at RHO-points.
!  pm,pn      Coordinate transformation metric "m" [1/meters]
!                    associated with the differential distances in 
!                    XI- and ETA-directions, both are at RHO-points.
!  xp,xr      XI-coordinates [m] at PSI- and RHO-points.
!  yp,yr      ETA-coordinates [m] at PSI- and RHO-points.
!-------------------------------------------------------------------
!
# ifdef AGRIF
      use Agrif_Util
# endif
      implicit none
      integer Istr,Iend,Jstr,Jend, i,j
# include "param.h"
# include "grid.h"
# include "scalars.h"
!
      real Length_XI,Length_ETA, depth, f0,beta, 
     &     cff,y, x0,y0,dx,dy, NSUB
!
# ifdef INTERNAL
      real cff1, ridge_width, ridge_height
# endif
# ifdef RIP
      real a1,a2,a3,b1,b2,b3,c1,c2,c3,h0,per,xx,yy,xs,db,alpha,lambda,eps,
     &     nodeX,nodeY,slope,xout,xtot
# endif
# ifdef FLUME
      real p1,p2,p3,p4
# endif
# if defined RIVER || defined IB09
      integer is,i0,j0
#  include "sources.h"
# endif
#ifdef MPI
# define LOCALLM Lmmpi
# define LOCALMM Mmmpi
#else
# define LOCALLM Lm
# define LOCALMM Mm
#endif
!
# include "compute_extended_bounds.h"
!
!----------------------------------------------------------------------
! Configuration parameters (h,f,beta)
!----------------------------------------------------------------------
!
# if defined BASIN
                     depth=5000.
                     f0=1.E-4 
                     beta=2.E-11
# elif defined SINGLE_COLUMN
#  if defined KATO_PHILIPS || defined WILLIS_DEARDORFF 
                     depth=50.
                     f0=0.
                     beta=0.
#  elif defined DIURNAL_CYCLE
                     depth=150.
                     f0=0.
                     beta=0.
#  elif defined FORCED_EKBBL
                     depth=1500.
                     f0=1.E-4
                     beta=0.
#  elif defined FORCED_DBLEEK
                     depth=25.0
                     f0=1.E-4
                     beta=0.
#  elif defined FORCED_NONROTBBL || defined FORCED_OSCNONROTBBL
                     depth=5.
                     f0=0.0
                     beta=0.
#  endif 
# elif defined CANYON_A || defined CANYON_B
                     depth=4000.
                     f0=1.E-4
                     beta=0.
# elif defined PLUME
                     depth=1300.
                     beta=0.
                     f0=0.5E-4 
# elif defined EQUATOR
                     depth=5000.
                     f0=0.
                     beta=2.2829e-11
# elif defined KH_INST 
                     depth=286.     
                     f0=0.
                     beta=0.
# elif defined GRAV_ADJ
#  ifdef NBQ
                     depth=0.29      ! Horn et al. (2001)
#  else
                     depth=20.
#  endif
                     f0=0.
                     beta=0.
# elif defined ACOUSTIC
                     depth=128.
                     f0=0.
                     beta=0.
# elif defined INNERSHELF
                     depth=0.
                     f0=4*pi/86400*sin(-21*pi/180)
                     beta=0.
# elif defined INTERNAL
                     depth=2000.
                     f0=1.e-4
                     beta=0.
                     ridge_width=30.e3
                     ridge_height=1600.
# elif defined TS_HADV_TEST
                     depth=10.
                     f0=0.
                     beta=0.                     
# elif defined OVERFLOW
                     depth=40.
                     f0=0.
                     beta=0.
# elif defined SEAMOUNT
                     depth=4500.
                     beta=0.
                     f0=1.E-4
# elif defined SHELFRONT
                     depth=1660.
                     f0=1.E-4
                     beta=0.
# elif defined SOLITON
                     depth=1.
                     f0=0.
                     beta=1.
# elif defined UPWELLING
                     depth=150.
                     f0=-8.26E-5
                     beta=0.
# elif defined SED_TEST2
                     depth=10.0    
                     f0=1.E-4
                     beta=0.
# elif defined RIVER
                     depth=150.
                     f0=8.26E-5
                     beta=0.
# elif defined JET
                     depth=4000.
                     f0=1.E-4
                     beta=1.6E-11
# elif defined SHOREFACE
                     depth=15.
                     f0=0.   
                     beta=0.
# elif defined FLUME
                     depth=0.6
                     f0=0.   
                     beta=0.
# elif defined SWASH
                     depth=1.
                     f0=0.   
                     beta=0.
# elif defined RIP
                     depth=7.
                     f0=1.E-4  
                     beta=0.
# elif defined THACKER
                     depth=10.
                     beta=0.
#  ifdef THACKER_2DV
                     f0=0.
#  else
                     f0=1.E-4
#  endif
# elif defined TANK
#  ifndef MOVING_BATHY 
                     depth=10.
#  else
                     depth=0.394
#  endif
                     f0=0.
                     beta=0.
# else
                     depth=???
                     f0=???
                     beta=???
# endif
!
!----------------------------------------------------------------------
! Grid dimensions (Length_XI, Length_ETA)
!----------------------------------------------------------------------
!
# ifdef AGRIF
      if (Agrif_Root()) then
# endif
# if defined BASIN
                     Length_XI =3600.0e+3
                     Length_ETA=2800.0e+3
# elif defined SINGLE_COLUMN
                     Length_XI = 10.0e+3
                     Length_ETA= 10.0e+3
# elif defined CANYON_A || defined CANYON_B
                     Length_XI =128.0e+3
                     Length_ETA=96.0e+3
# elif defined EQUATOR
                     Length_ETA=3200.0E+3
                     Length_XI =4000.0E+3
# elif defined KH_INST 
                     Length_XI =286.0
                     Length_ETA=Length_XI*float(MMm0)/float(LLm0) 
!                    Length_ETA=286.0/256.0
# elif defined ACOUSTIC
                     Length_XI =128.0
                     Length_ETA=2.0
# elif defined GRAV_ADJ
#  ifdef NBQ
#   ifdef GRAV_ADJ_SOLITON
                     Length_XI =6.0    ! Horn et al. (2001)
#   else
                     Length_XI =3.0
#   endif
                     Length_ETA=0.1
#  else
                     Length_XI =64.0E+3
                     Length_ETA=2.0E+3
#  endif
# elif defined INNERSHELF
                     Length_XI =200.0e+3
                     Length_ETA=3.0e+3
# elif defined INTERNAL
                     Length_XI =1200.0e+3
                     Length_ETA=6.0e+3
!                     Length_XI =2400.0e+3
!                     Length_ETA=100.0e+3
# elif defined TS_HADV_TEST
                     Length_XI  = 1.0E+3
                     Length_ETA = 1.0E+3
# elif defined OVERFLOW
                     Length_XI =2.0E+3
                     Length_ETA=64.0E+3
# elif defined SEAMOUNT
                     Length_XI =512.0e+3
                     Length_ETA=512.0e+3
# elif defined SHELFRONT
                     Length_XI =20.0e+3
                     Length_ETA=200.0e+3
# elif defined SOLITON
                     Length_XI =48.
                     Length_ETA=16.
# elif defined UPWELLING
                     Length_ETA=8.E+4
                     Length_XI =LLm*Length_ETA/MMm
# elif defined SED_TEST2
                     Length_XI =25000.0
                     Length_ETA=2500.0
# elif defined RIVER
                     Length_ETA=8.E+4
                     Length_XI =LLm*Length_ETA/MMm
# elif defined JET
                     Length_XI =500.e+3 
                     Length_ETA=2000.e+3 
# elif defined SHOREFACE
                    Length_XI =1180.0 
                    Length_ETA=140.0
# elif defined FLUME
                    Length_XI =32.0 
                    Length_ETA=1.0
# elif defined SWASH
                    Length_XI =110.0 
                    Length_ETA=MMm*Length_XI/LLm
# elif defined RIP
                    Length_XI =350     ! RIP XI
                    Length_ETA=500    ! RIP ETA
# elif defined THACKER
                     Length_XI =199.e+3 
#  ifdef THACKER_2DV
                     Length_ETA=3.e+3
#  else
                     Length_ETA=199.e+3
#  endif
# elif defined TANK
#  ifndef MOVING_BATHY
#   ifndef TANKY
                     Length_XI =10.0
                     Length_ETA=0.20*MMm0
#   else
                     Length_XI =0.20*LLm0
                     Length_ETA=10.0
#   endif
#  else
                     Length_XI =4.0
                     Length_ETA=0.001
#  endif
# else
                     Length_XI =???
                     Length_ETA=??? 
# endif

!----------------------------------------------------------------------
! Copy physical dimensions of the grid into globally vizible variables
!----------------------------------------------------------------------
!
        xl=Length_XI
        el=Length_ETA

# ifdef AGRIF
      else      
        xl=Agrif_Parent(xl)/Agrif_Parent(LLm)*LLm/Agrif_Rhox()
        el=Agrif_Parent(el)/Agrif_Parent(MMm)*MMm/Agrif_Rhoy()
        Length_XI=xl
        Length_ETA=el
      endif
# endif /* AGRIF */
!
!
!----------------------------------------------------------------------
! Set grid spacings for rectangular grids
!----------------------------------------------------------------------
!
      dx=Length_XI/float(LLm)    
      dy=Length_ETA/float(MMm)   
!
!----------------------------------------------------------------------
! Set reference point location
!----------------------------------------------------------------------
!
      x0=0.
      y0=0.
# ifdef EQUATOR
      y0=y0-0.5*Length_ETA  
# endif
# if defined INTERNAL || defined THACKER
      x0=x0-0.5*Length_XI
      y0=y0-0.5*Length_ETA  
# endif


CR      write(*,'(4(A3,pe15.9,1x),I3)') 'dx=',dx, 'dy=',dy,
CR     &                                'x0=',x0, 'y0=',y0, mynode
!
!----------------------------------------------------------------------
! Setup rectangulag grid:  coordinates (XI,ETA) at PSI- and RHO-points.
!----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
         xp(i,j) = x0 + dx*float(i+iminmpi-1-1)
         xr(i,j) = x0 + dx*(float(i+iminmpi-1-1)+0.5)
         yp(i,j) = y0 + dy*float(j+jminmpi-1-1)
         yr(i,j) = y0 + dy*(float(j+jminmpi-1-1)+0.5)
        enddo
      enddo
!
!----------------------------------------------------------------------
! Compute coordinate transformation metrics at RHO-points "pm" and
! "pn"  (1/m) associated with the differential distances in XI and
! ETA, respectively.
!----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          pm(i,j)=1./dx
          pn(i,j)=1./dy
        enddo
      enddo
!
!----------------------------------------------------------------------
!  Set Coriolis parameter [1/s] at RHO-points.
!----------------------------------------------------------------------
!
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef EQUATOR
          f(i,j)=f0+beta*yr(i,j)
#  else
          f(i,j)=f0+beta*(yr(i,j)-(0.5*el))
#  endif
        enddo
      enddo

# if defined UV_COR_NT || defined CROCO_QH
!
!  Horizontal (non-traditional) Coriolis parameter
!
      do j=JstrR,JendR
        do i=IstrR,IendR
          if (beta.eq.0) then     ! f-plan
            e(i,j)=2.*Erotation*cos(asin(f0/(2.*Erotation)))
          else
            e(i,j)=beta*Eradius   ! beta plan
          endif
        enddo
      enddo
# endif
!
!----------------------------------------------------------------------
!  Set bathymetry [meters; positive] at RHO-points.
!----------------------------------------------------------------------
!
# if defined BASIN || defined EQUATOR || defined GRAV_ADJ \
                   || defined SOLITON || defined JET \
                   || defined KH_INST || defined TS_HADV_TEST \
                   || defined ACOUSTIC || defined SINGLE_COLUMN
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth
        enddo
      enddo

# elif defined CANYON_A || defined CANYON_B
      do j=JstrR,JendR
        do i=IstrR,IendR
          cff=32000.-16000.*(sin(pi*xr(i,j)/Length_XI))**24
          h(i,j)=20.+0.5*(depth-20.)*(1.+tanh((yr(i,j)-cff)/10000.))
        enddo
      enddo

# elif defined INNERSHELF
      do j=JstrR,JendR
        do i=IstrR,IendR
!         --- constant slope ---
          h(i,j)=1.e-3*(xl+0.5*dx-xr(i,j))+4.
!         --- shelf/slope tanh ---
!          h(i,j)=30.+500.*(1+sinh((150.e3-xr(i,j))/40.e3)
!     &                      /cosh((150.e3-xr(i,j))/40.e3))
        enddo
      enddo

# elif defined INTERNAL
      cff1=1./(ridge_width*ridge_width)
      do j=JstrR,JendR
        do i=IstrR,IendR
          if ((xr(i,j).gt.ridge_width).or.
     &       (xr(i,j).lt.(-ridge_width))) then
            cff=0.
          else
            cff=1-(xr(i,j)*xr(i,j)*cff1);
          endif
          h(i,j)=depth-ridge_height*cff*cff;
#   ifdef INTERNALSHELF
          if (xr(i,j).ge.0.) then
            h(i,j)=depth-ridge_height
          endif
#   endif
        enddo
      enddo

# elif defined OVERFLOW
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=20.+0.5*(depth-20.)*(1.+
     &                         tanh((yr(i,j)-40000.)/5000.))
        enddo
      enddo

# elif defined SEAMOUNT
      cff=(1./50.e3)**2
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth*(1-0.6*exp(-cff*((xr(i,j)-0.5*xl)**2+
     &                                  (yr(i,j)-0.5*el)**2))) 
        enddo
      enddo

# elif defined SHELFRONT
      do j=JstrR,JendR
        do i=IstrR,IendR
          cff=yr(i,j)/1000.
          if (cff.lt.50.) then
            h(i,j)=50.+2.*cff
          elseif (cff.lt.60.) then
            h(i,j)=160.+1.5*(cff-50.)**2-0.1*(cff-60.0)**2
          elseif (cff.lt.100.) then
            h(i,j)=310.+30.*(cff-60.)
          elseif (cff.lt.110.) then
            h(i,j)=1660.-1.5*(cff-110.)**2
          else
            h(i,j)=1660.
          endif
        enddo
      enddo

# elif defined SED_TEST2
      do j=Jstr,Jend
        do i=Istr,Iend
          val1=user(3)            !ocean-side depth
          val2=user(4)            !nearshore depth
          val3=user(5)/3.0        ! ~ width of slope
          h(i,j)=0.5*(val1+
     &     (val1-val2)*TANH((7000.0-xr(i,j))/val3))  
        enddo
      enddo

# elif defined UPWELLING
      do j=JstrR,JendR
#  ifdef MPI
        y=dy*float(j+jj*Mm)
#  else
        y=dy*float(j)
#  endif
        if (y.gt.Length_ETA/2.) y=Length_ETA-y+dy
        cff=min(depth,84.5+66.526*tanh(0.00015*(y-0.125*Length_ETA)))
        do i=IstrR,IendR
          h(i,j)=cff
        enddo
      enddo

# elif defined RIVER
      do i=IstrR,IendR
#  ifdef MPI
        cff=(float(i +ii*Lm)-0.5)/float(LLm)
#  else
        cff=(float(i       )-0.5)/float(LLm)
#  endif
        if (cff.lt.0.05) then
          h(i,JstrR)=15.
        elseif (cff.lt.0.15) then
          h(i,JstrR)=15.+843.75*(cff-0.05)**2
        elseif (cff.lt.0.85) then
          h(i,JstrR)=15.+168.75*(cff-0.1)
        elseif (cff.lt.0.95) then
          h(i,JstrR)=150.-843.75*(cff-0.95)**2
        else
          h(i,JstrR)=150.
        endif
      enddo
      do j=JstrR+1,JendR
        do i=IstrR,IendR
          h(i,j)=h(i,JstrR)
        enddo
      enddo

# elif defined SHOREFACE
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
          cff=(float(i +ii*Lm) )
#  else
          cff=(float(i       ) )
#  endif
#  ifdef WET_DRY
          h(i,j)=12.0-0.0125*Length_XI/float(LLm+1)*cff ! original fonction
#  else
          h(i,j)=12.0-0.0125*Length_XI/float(LLm+1)*cff + 3.
#  endif
#  ifdef MASKING
          rmask(i,j)=1.0
          if (h(i,j).lt.0.01) then
            h(i,j)=0.01
            rmask(i,j)=0.0
            rmask(48,j)=1.0
          endif
#  endif
        enddo
      enddo

# elif defined FLUME
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
          cff=(float(i +ii*Lm) )
#  else
          cff=(float(i       ) )
#  endif
!
! Grasso et al. 2009
!          p1 = -6.953e-05
!          p2 =  0.003705 
!          p3 = -0.06435
!          p4 =  0.5716
!
! Roelvink & Stive 1989
          p1 =  -6.228e-06
          p2 =   0.0004507
          p3 =     0.01459
          p4 =     -0.6769
!
          h(i,j) = p1*xr(i,j)**3 + p2*xr(i,j)**2 + p3*xr(i,j) + p4
          h(i,j) = min(-h(i,j),depth)
        enddo
      enddo

# elif defined SWASH
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef MPI
          cff=(float(i +ii*Lm) )
#  else
          cff=(float(i       ) )
#  endif
          h(i,j)=1.0571-0.0125*Length_XI/float(LLm+1)*cff  ! GLOBEX profil
#  ifndef WET_DRY
     &                                              + 0.6
#  endif
          if (xr(i,j).lt.16.57) h(i,j)=0.85               ! GLOBEX
        enddo
      enddo

# elif defined RIP
#  ifdef GRANDPOPO
!
!  Idealized Grand Popo Beach in Benin
!
      xs=85         ! inner surf zone 
      db=50         ! distance from xs to sand bar
      alpha=0.025
      lambda=10     ! sand bar perturbation wavelength
      do j=JstrR,JendR
        do i=IstrR,IendR
#   ifdef OBC_EAST
          xx=xr(i,j)
#   else
          xx=Length_XI-xr(i,j)
#   endif
          h(i,j)=-4.6
     &       -1.5*exp(-6*(((xx-xs-db)/db)**2))
     &       +3.1*(1+tanh(0.025*(xx-xs)))
     &       +0.014*(xx+log(cosh(alpha*(xx-xs))/cosh(alpha*xs))/alpha)
          h(i,j)=min(h(i,j), 11.3-0.88e-5*(abs(xx-700))**2.17)

!          db=55.  !!! DUKE
!          h(i,j)=-4.5
!     &       -1.0*exp(-6*(((xx-xs-db)/db)**2))
!     &       +2.95*(1+tanh(0.025*(xx-xs)))
!     &       +0.01*(xx+log(cosh(alpha*(xx-xs))/cosh(alpha*xs))/alpha)
!          h(i,j)=min(h(i,j), 8.3-0.9e-5*(abs(xx-700))**2.12)

#   ifdef GRANDPOPO_PERTURB
          yy=yr(i,j)
          eps=0.01
          per=eps*cos(2*pi*yy/lambda)*exp(-5*(((xx-xs-db)/db)**2))
          h(i,j)=(1+per)*h(i,j)
#   endif
        enddo
     enddo

#  elif defined IB09
!
!  Idealized Imperial Beach in California
!
#   define HR16B2 
#   ifdef HR16B3 
! HR16 measured bathy
      a1 = 67.33
      b1 =0.01802
      c1 =0.1579
      a2 = 64.91
      b2 = 0.01844
      c2 = -2.979
      a3 = 0.1357
      b3 = 0.06453
      c3 = -1.341
      slope=2.8
      xtot=450
      xout=50
      nodeX=-149.8485
      nodeY=-4.4402

#   elif defined HR16B2 
! HR16 longshore-averaged bathy
      a1=9.509;
      b1=0.002914;
      c1=0.014;
      a2=0.5618;
      b2=0.02377;
      c2=-2.684;
      a3=0.1274;
      b3=0.05868;
      c3=1.32;
#   endif
      xs=40         ! inner surf zone 
      db=40         ! distance from xs to sand bar
      alpha=0.025
      lambda=10     ! sand bar perturbation wavelength

      do j=JstrR,JendR
        do i=IstrR,IendR
#   ifdef OBC_EAST
          xx=xr(i,j)
#   else
          xx=Length_XI-xr(i,j)
#   endif
#   ifdef HR16B3
          if (xx.gt.(xout-nodeX)) then
                xx=((xout-xx-nodeX)+(nodeX+600-xout))
     &              /(nodeX+600-xout)*(1+slope)-slope
                h(i,j)=tanh(1.0-xx)*(7+nodeY)-nodeY-0.0287
          elseif (xx.lt.(xout+0.1693)) then
                h(i,j)=-0.02*(xout-xx+0.1693)+0.0386+0.0158
          else
                h(i,j)=-(a1*sin(b1*(xout-xx)+c1) +
     &                   a2*sin(b2*(xout-xx)+c2) +
     &                   a3*sin(b3*(xout-xx)+c3))
          endif
#   elif defined HR16B2
          if (xx.ge.320) then
                h(i,j)=depth-0.05
          elseif (xx.le.20) then
                h(i,j)=-(20-xx)*0.02
          elseif ((xx.gt.20) .AND. (xx.lt.410)) then
                h(i,j)=-(a1*sin(b1*(20-xx)+c1) +
     &                   a2*sin(b2*(20-xx)+c2) +
     &                   a3*sin(b3*(20-xx)+c3))
          endif
#   endif
#   ifdef IB09_PERTURB
          yy=yr(i,j)
          eps=0.01
          per=eps*cos(2*pi*yy/lambda)*exp(-5*(((xx-xs-db)/db)**2))
          h(i,j)=(1+per)*h(i,j)
#   endif
        enddo
      enddo

#  else
!
! Weir et al. 2011 after Lippmann et al. 1999 (idealization of Duck beach)
!
      xs=150      ! inner surf zone 
      db=80       ! distance from xs to sand bar
      alpha=0.02
      lambda=256  ! sand bar perturbation wavelength
#  ifdef RIP_TOPO_2D
      eps=0.0     ! no sand bar perturbation
#  else
      eps=0.1     ! sand bar perturbation magnitude
#  endif
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef OBC_EAST
          xx=xr(i,j)
#  else
          xx=Length_XI-xr(i,j)
#  endif
          yy=yr(i,j);
          h0=-1.5*exp(-5*(((xx-xs-db)/db)**2))
     &       +1.35*(1+tanh(0.0253*(xx-xs)))
     &       +0.0032*(xx+log(cosh(alpha*(xx-xs))/cosh(alpha*xs))/alpha)
          per=eps*cos(2*pi*yy/lambda)*exp(-5*(((xx-xs-db)/db)**2))
          h(i,j)=(1+per)*h0!-0.8;
        enddo
      enddo
#  endif

# elif defined THACKER
      cff=80.e3
      do j=JstrR,JendR
        do i=IstrR,IendR
#  ifdef THACKER_2DV
          h(i,j)=depth*(1-(xr(i,j)/cff)**2)
#  else
          h(i,j)=depth*(1-(xr(i,j)/cff)**2 -
     &                    (yr(i,j)/cff)**2 )
#  endif
        enddo
      enddo

# elif defined TANK
#  ifndef MOVING_BATHY
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth
        enddo
      enddo
#  else
      cff=1./(0.03686)**2
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=depth-0.1*exp(-cff*((xr(i,j)-2.)**2+
     &                               (yr(i,j)-el*0.5)**2)) 
        enddo
      enddo
#  endif      
# else
      do j=JstrR,JendR
        do i=IstrR,IendR
          h(i,j)=???
        enddo
      enddo
# endif

!----------------------------------------------------------------------
! initialisation of z0b in case of analytical cases
!----------------------------------------------------------------------
      do j=JstrR,JendR
        do i=IstrR,IendR
          zob(i,j)=zobt
        end do
      end do
!
!----------------------------------------------------------------------
!  Set masking at RHO-points.
!----------------------------------------------------------------------
!
# ifdef MASKING
#  ifdef RIVER
      do j=JstrR,JendR             ! Set mask to all-water status,
        do i=IstrR,IendR           ! then mask out 8-point wide strip
          rmask(i,j)=1.            ! of land on the west, and finally
        enddo                      ! carve two 1-point wide channels
      enddo                        ! through that strip.
      if (WESTERN_EDGE) then
        do j=JstrR,JendR
          do i=IstrR,8
            rmask(i,j)=0.            ! <-- strip of land
          enddo
        enddo
        do is=1,Nsrc
#  ifdef MPI
          i0=Isrc_mpi(is,mynode)
          j0=Jsrc_mpi(is,mynode)
#  else
          i0=Isrc(is)
          j0=Jsrc(is)
#  endif
          if (IstrR.le.i0 .and. i0.le.IendR .and.
     &        JstrR.le.j0 .and. j0.le.JendR) then
            if (is.eq.1) then
              do j=j0,j0+6
                rmask(3,j)=1.      ! <-- upper channel with corner
              enddo
              do i=IstrR+3,8
                rmask(i,j0+6)=1.
              enddo
            else
              do i=IstrR+3,8
                rmask(i,j0)=1.     ! <-- lower channel along xi
              enddo
            endif
          endif
        enddo
      endif ! WESTERN_EDGE
#  elif defined IB09
      do j=JstrR,JendR
        do i=IstrR,IendR
          rmask(i,j)=1.
        enddo
      enddo
      do is=1,Nsrc
#   ifdef MPI
        i0=Isrc_mpi(is,mynode)
        j0=Jsrc_mpi(is,mynode)
#   else
        i0=Isrc(is)
        j0=Jsrc(is)
#   endif
        if (IstrR.le.i0 .and. i0.le.IendR .and.
     &      JstrR.le.j0 .and. j0.le.JendR
     &      .and. i0.ne.-1 .and. j0.ne.-1 ) then
          do i=i0,min(IendR,i0+6)
            rmask(i,j0)=0.
          enddo
        endif
      enddo
#  endif
# endif
!
#else
      subroutine ana_grid_empty
#endif /* ANA_GRID */
      return
      end

