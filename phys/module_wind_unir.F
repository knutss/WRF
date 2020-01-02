
!WRF:MODEL_LAYER:PHYSICS

MODULE module_wind_unir
!
!
! This module resembles the socalled "fitch" module in structure
! but is based on a different approach.
! The scientific desciption of the algorithm
! can be found in:
! "A revised mesoscale wind farm model implemented in the 
!  Weather Research and Forecasting(WRF) code" 
! by
! H Heiberg-Andersen and I Barstad (2016)
! in MWR

  USE module_driver_constants, ONLY : max_domains
  USE module_model_constants,  ONLY : piconst
  USE module_llxy
  USE module_dm,               ONLY : wrf_dm_min_real
  USE module_configure,        ONLY : grid_config_rec_type

  IMPLICIT NONE

  INTEGER, PARAMETER :: MAXVALS  = 100   
  INTEGER, PARAMETER :: MAXVALS2 = 100     
!
  INTEGER           :: nt 
  INTEGER           :: flow_induction  
  REAL              :: tiplossf,speedupstr
  INTEGER, DIMENSION(:), ALLOCATABLE   :: nkind
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival,jval
  REAL, DIMENSION(:), ALLOCATABLE :: hubheight,diameter,stc,stc2,cutin,cutout,npower
  LOGICAL           :: tke_adv
!
  REAL :: turbws(maxvals,maxvals2),turbtc(maxvals,maxvals2),turbpw(maxvals,maxvals2)
  CHARACTER*256  :: message_wind

!
CONTAINS


  SUBROUTINE  dragforceunir(                     &
       & id                                      &
       &,z_at_w,u,v                              &
       &,dx,dz,dt,qke                            &
       &,du,dv                                   &
       &,windfarm_opt,power                      &
       &,ids,ide,jds,jde,kds,kde                 &
       &,ims,ime,jms,jme,kms,kme                 &
       &,its,ite,jts,jte,kts,kte                 &
       &                                         )  
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: id,windfarm_opt 
  INTEGER, INTENT(IN) :: its,ite,jts,jte,kts,kte
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
  REAL, INTENT(IN) :: dx,dt
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz,u,v,z_at_w
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: du,dv,qke
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: power
  CHARACTER*256 message_wind                             

!
! Local
!
  REAL, ALLOCATABLE, DIMENSION(:) :: weight
  REAL     powcof,thrcof,c0,c1,c2,test_nan
  INTEGER  itf,jtf,ktf
  INTEGER  i,j,k,n
  INTEGER  kbot, ktop


  REAL :: speedhub, area, uhub, vhub, induction
  REAL :: ind_fac_sum,rad2deg,ue,ve
  INTEGER :: kt  

  REAL, PARAMETER :: nn=0.5, pi=3.14, cepsilon=0.168, beta=1.414, rho=1.225, rootfrac=0.1
  REAL :: area1,area2, smom
  REAL :: length, afrac, rw,thrust, sumweight, z1, z2


  itf=MIN0(ite,ide-1)
  jtf=MIN0(jte,jde-1)
  ktf=MIN0(kte,kde-1)

    power       = 0.
    ind_fac_sum = 0.
    rad2deg     = 45./atan(1.0)



    DO kt = 1,nt


!
! vertical layers cut by turbine blades
!
        i = ival(kt,id)
        j = jval(kt,id)

        if (i.ne.-9999.and.j.ne.-9999) then
        IF (( its .LE. i .AND. i .LE. itf ) .AND. &
            ( jts .LE. j .AND. j .LE. jtf )  ) THEN


           kbot=0 ! The lower blade tip is between kbot and kbot+1
           DO k=kts,ktf-1
             IF (hubheight(kt)-diameter(kt)/2. .GE. z_at_w(i,k,j) .AND. &
                 hubheight(kt)-diameter(kt)/2. .LE. z_at_w(i,k+1,j)) THEN
              kbot=k
              EXIT
             ENDIF
           ENDDO
           
           IF (kbot .EQ. 0) THEN
             print*, 'No affected cells found for (i,j)=(',i,j,')'
             GOTO 1000
           ENDIF
    
           ktop=kbot ! The upper blade tip is between ktop and ktop+1
           DO k=kbot+1,ktf
             IF (hubheight(kt)+diameter(kt)/2. .GE. z_at_w(i,k,j) .AND. &
                 hubheight(kt)+diameter(kt)/2. .LE. z_at_w(i,k+1,j)) THEN
              ktop=k
              EXIT
             ENDIF
           ENDDO
         
           IF (ktop .EQ. kbot) THEN
             print*, 'No affected cells found for (i,j)=(',i,j,')'
             GOTO 1000
           ENDIF

           area=piconst/4.*diameter(kt)**2.          ! area swept by turbine blades

!Compute speed and speed components at hub height at the turbine location

           uhub = 0.
           vhub = 0.
           DO k=kbot,ktop 
              uhub=uhub+u(i,k,j)*dz(i,k,j)/(z_at_w(i,ktop+1,j)-z_at_w(i,kbot,j))
              vhub=vhub+v(i,k,j)*dz(i,k,j)/(z_at_w(i,ktop+1,j)-z_at_w(i,kbot,j))
            ENDDO


            speedhub=sqrt(uhub**2+vhub**2)
            IF (speedhub < 0.01) THEN
              length=dx
            ELSE  
              IF (abs(uhub) .GE. abs(vhub)) THEN
                length=abs(uhub)*dx/speedhub
              ELSE
                length=abs(vhub)*dx/speedhub
              ENDIF
            ENDIF

            rw=beta*diameter(kt)/2.
          
            afrac=piconst*rw**2-rw**2*acos((z_at_w(i,ktop+1,j)-hubheight(kt))/rw) &
                 -(z_at_w(i,ktop+1,j)-hubheight(kt))*sqrt(rw**2-(z_at_w(i,ktop+1,j)-hubheight(kt))**2) &
                 -rw**2*acos((hubheight(kt)-z_at_w(i,kbot,j))/rw) &
                 -(hubheight(kt)-z_at_w(i,kbot,j))*sqrt(rw**2-(hubheight(kt)-z_at_w(i,kbot,j))**2)

            afrac=afrac/length/(z_at_w(i,ktop+1,j)-z_at_w(i,kbot,j))


          !Compute upstream wind velocity
            IF (flow_induction .EQ. 0) THEN
               induction=0.
            ELSE
            if (speedhub.ne.0.)then
            CALL induction_factorunir(speedhub,afrac,cutin(kt),cutout(kt),stc(kt),stc2(kt),nkind(kt),induction,tiplossf)
            else
             induction=0.
            endif
            ENDIF

            speedupstr=speedhub/(1.-afrac*induction)

            CALL dragcofunir(powcof,thrcof,               &
                           speedupstr,cutin(kt),cutout(kt),   &
                           npower(kt),diameter(kt),stc(kt),stc2(kt),nkind(kt))
          !Compute power produced by turbine


	    c2=2.*(sqrt(u(i,ktop,j)**2.+v(i,ktop,j)**2.)+sqrt(u(i,kbot,j)**2.+v(i,kbot,j)**2.)-2.*speedupstr)
            c2=c2/(diameter(kt)/2.)**2/3.

            c1=sqrt(u(i,ktop,j)**2.+v(i,ktop,j)**2.)-sqrt(u(i,kbot,j)**2.+v(i,kbot,j)**2.)
            c1=c1/diameter(kt) 
            c0=speedupstr-c2*(diameter(kt)/4.)**2.

            power(i,j)=power(i,j)+powcof*0.5*1.23*(c0**3.+(3.*c0**2.*c2+3.*c1**2.*c0)*((diameter(kt)/2.)**2.)/4.)*area
            power(i,j)=power(i,j) &
                       +powcof*0.5*1.23*((3.*c1**2.*c2+3.*c0*c2**2.)*((diameter(kt)/2.)**4.)/8.  &
	               +5.*c2**3.*((diameter(kt)/2.)**6.)/64.)*area 
	    

             thrust=thrcof**0.5*rho*(c0**2.+(c1**2.+2.*c0*c2)*((diameter(kt)/2.)**2.)/4.+c2**2.*((diameter(kt)/2.)**4.)/8.)*area 


            ALLOCATE(weight(ktop))

! -- prandtl-distribution:
	    sumweight=0.
	    DO k=kbot,ktop
	      weight(k)=0.
		z1=z_at_w(i,k,j)-hubheight(kt)
		z2=z_at_w(i,k+1,j)-hubheight(kt)
		IF (abs(z1) .LE. diameter(kt)/2. .AND. abs(z2) .LE. diameter(kt)/2.) THEN
                   weight(k)=1.
		ENDIF
	       sumweight=sumweight+weight(k)
	    ENDDO

             qke(i,kbot,j)=2.**(1./3.)*((thrust*speedhub-power(i,j))/cepsilon/nn)**(2./3.) &
                                    *(diameter(kt)/2./rho/dx/dx/dz(i,kbot,j))**(2./3.)
             qke(i,ktop,j)=2.**(1./3.)*((thrust*speedhub-power(i,j))/cepsilon/nn)**(2./3.) &
                                    *(diameter(kt)/2./rho/dx/dx/dz(i,ktop,j))**(2./3.)



            DO k=kbot,ktop
            du(i,k,j) = du(i,k,j)-thrust*weight(k)*uhub/sumweight/speedhub/rho/dx/dx/dz(i,k,j)
            dv(i,k,j) = dv(i,k,j)-thrust*weight(k)*vhub/sumweight/speedhub/rho/dx/dx/dz(i,k,j)
            test_nan=qke(i+1,k,j)-qke(i-1,k,j)-qke(i,k,j+1)-qke(i,k,j-1)
            if(test_nan.ge.0. .or. test_nan.lt. 0.)then   ! nan is not equal to itself
             du(i,k,j) = du(i,k,j)-(qke(i+1,k,j)-qke(i-1,k,j))/dx/6.
             dv(i,k,j) = dv(i,k,j)-(qke(i,k,j+1)-qke(i,k,j-1))/dx/6.
            endif 
            ENDDO                 ! k-loop

           DEALLOCATE(weight) 
       

        ENDIF                     ! its-loop
1000    CONTINUE
        endif                     ! -9999 -loop
    ENDDO                         ! kt-loop


  END SUBROUTINE dragforceunir


  SUBROUTINE dragcofunir(powcof,thrcof,speed,cutin,cutout, &
                     tpower,tdiameter,stdthrcoef,stdthrcoef2,nk)

  IMPLICIT NONE
  REAL, INTENT(IN):: speed, cutin, cutout, tpower,tdiameter,stdthrcoef,stdthrcoef2
  REAL, INTENT(OUT):: powcof,thrcof
  REAL :: power,area,mspeed,hspeed
!
! ... PAJ ...
!
   INTEGER :: nk,k,nu,nb
   LOGICAL :: vfound
   REAL :: fac1,fac2
!
  area=piconst/4.*tdiameter**2.          ! area swept by turbine blades
      nu=2
      nb=1
      vfound=.false.
      DO k=3,maxvals2
            IF(.NOT. vfound) THEN
              IF(turbws(nk,k).GT.speed) THEN
                nu=k 
                nb=k-1
                vfound=.true.
              ENDIF
            ENDIF
      ENDDO
!
  IF (speed .LE. cutin) THEN
     thrcof = stdthrcoef
  ELSE
    IF (speed .GE. cutout) THEN
     thrcof = stdthrcoef2
     ELSE
     thrcof = turbtc(nk,nb)+(turbtc(nk,nu)-turbtc(nk,nb))/(turbws(nk,nu)-turbws(nk,nb))*(speed-turbws(nk,nb))
    ENDIF
  ENDIF
!
! ... power coeficient ...
!
  IF(speed .LE. cutin .OR. speed .GE. cutout) THEN
     power=0.
     powcof=0.
  ELSE
      fac1=1000./(0.5*1.23*turbws(nk,nb)**3.*area)
      fac2=1000./(0.5*1.23*turbws(nk,nu)**3.*area)
      power = turbpw(nk,nb)+(turbpw(nk,nu)-turbpw(nk,nb))/(turbws(nk,nu)-turbws(nk,nb)) &
                               *(speed-turbws(nk,nb))
      powcof = turbpw(nk,nb)*fac1+(turbpw(nk,nu)*fac2-turbpw(nk,nb)*fac1)/(turbws(nk,nu)-turbws(nk,nb)) &
                                     *(speed-turbws(nk,nb))
  ENDIF

  END SUBROUTINE dragcofunir
!

  SUBROUTINE induction_factorunir(speedhub,afrac,cutin,cutout,stdthrcoef,stdthrcoef2,nk,a,F)
  
   IMPLICIT NONE
   REAL, INTENT(IN) :: speedhub, afrac,cutin, cutout, stdthrcoef, stdthrcoef2, F 
   REAL, INTENT(OUT) ::  a                    
   INTEGER :: nk, j, k, nu, nb
   LOGICAL :: vfound
   REAL :: speed, af, ap, b0, b1, b2 
   REAL :: thrcof
   REAL, PARAMETER :: a0 = 0.4
   CHARACTER*256 message_wind                             
   INTEGER, PARAMETER :: iterations = 20


   a=a0
  
   b1 = 4*(a0*a0*F-a0*(2*F+1.)+F)/(a0-1.)/(a0-1.)
   b2 = (2.-4*(a0-1.)*(a0-1.)*F)/(a0-1.)/(a0-1.)
   b0 = 2.-b1-b2

   DO j=1, iterations
      speed=speedhub/(1.-2.*afrac*a)
      vfound=.false.
      nb=1
      nu=2
      DO k=3,maxvals2
            IF(.NOT. vfound) THEN
              IF(turbws(nk,k).GT.speed) THEN
                nu=k
                nb=k-1
                vfound=.true.
              ENDIF
            ENDIF
      ENDDO 


      IF (speed .LE. cutin) THEN
       thrcof = stdthrcoef
      ELSE
       IF (speed .GE. cutout) THEN
        thrcof = stdthrcoef2
       ELSE
        thrcof = turbtc(nk,nb)+(turbtc(nk,nu)-turbtc(nk,nb))/(turbws(nk,nu)-turbws(nk,nb))*(speed-turbws(nk,nb))
       ENDIF
      ENDIF

 
      a = (1.-sqrt(1.-amin1(thrcof/F,1.)))/2. 



      if ( a .gt. a0 ) then
         if ( 4.*b2*(b0-thrcof) .gt. b1**2 ) then
             a = -b1/b2/2.
         else
             a = -(b1+sqrt(b1**2-4.*b2*(b0-thrcof)))/2./b2
             if ( a .lt. a0 ) then
                 a = -(b1-sqrt(b1**2-4.*b2*(b0-thrcof)))/2./b2
             endif
         endif
      endif


      if (j .eq.iterations-1) then
         af=a
      endif
      if (j .eq.iterations) then
         ap=a
      endif
     if(a.gt.0.99 .and. a.lt.1.01)goto 200
    ENDDO
      a=(af+ap)/2.

 200 continue

    if(a .ge.0.95)then
      write(message_wind,*)'induction factor is bigger than 1.0',a
      CALL wrf_debug(0,message_wind)
      a=0.95        
    elseif(a.le.0)then
      write(message_wind,*)'induction factor is smaller than 0.0',a,speedhub
      CALL wrf_debug(0,message_wind)
      a=0.
    endif


 END SUBROUTINE induction_factorunir

      REAL FUNCTION func(R1,R2,U,z,y)
      implicit none
      REAL, INTENT(IN) :: R1, R2, U, z, y
      REAL :: B, a, pi, at, r, argt, argr, w, rfunc

      B=3.
      a=0.333
      pi=3.14

      r=sqrt(z*z+y*y)


      IF (r .GE. R2) THEN
         func=0.
      ELSE
        at=(1.-3.*a)/(4.*a-1.)
        w=2.*(U/R2)*sqrt(a*(1.-a)/at/(1.+at))
        argt=-(B/2.)*((R2-r)/r)*sqrt(1+(r*w/U/(1.-a))**2)
        IF (argt .GT. 30.) THEN
           func=1.
        ELSE
           func=(2./pi)*acos(exp(argt))
        ENDIF
      ENDIF

      IF (r .LE. R1) THEN
         rfunc=0.
      ELSE
         argr=-(B/2.)*((r-R1)/r)*sqrt(1+(r*w/U/(1.-a))**2)
         IF (argr .GT. 30.) THEN
            rfunc=1.
         ELSE   
            rfunc=(2./pi)*acos(exp(argr))
         ENDIF
      ENDIF

      func=func*rfunc

      END FUNCTION func



      REAL FUNCTION ymin(R1,z)
      implicit none
      REAL, INTENT(IN) :: R1, z
      IF (z*z .GE. R1*R1) THEN
        ymin=0.
      ELSE
        ymin=sqrt(R1*R1-z*z)
      ENDIF
      END FUNCTION ymin

      REAL FUNCTION ymax(R2,z)
      implicit none
      REAL, INTENT(IN) :: R2, z
      ymax=sqrt(R2*R2-z*z)
      END FUNCTION ymax
      
      REAL FUNCTION innerint(R1,R2,U,z)
      implicit none
      REAL, INTENT(IN) :: R1, R2, U, z
      CALL simpson(R1,R2,U,ymin(R1,z),ymax(R2,z),z,innerint)
      END FUNCTION innerint



        SUBROUTINE simpson(R1,R2,U,a,b,z,s)
        implicit none
        REAL, INTENT(IN) :: R1, R2, U, a, b, z
        REAL, INTENT(OUT) :: s

        INTEGER :: j
        INTEGER, PARAMETER :: jmax=40
        REAL, PARAMETER    :: eps=1.e-3
        REAL    :: st, ost, os 

        ost=-1.e30
        os=-1.e30

        st=0.
        DO j=1,jmax
           CALL trapz(R1,R2,U,a,b,z,j,st)
           s=(4.*st-ost)/3.
           IF (j .GT. 5) THEN
              IF (abs(s-os) .LT. eps*abs(os) .OR. &
                  s .EQ. 0. .AND. os .EQ. 0.) THEN
                  RETURN
              ENDIF
           ENDIF
           os=s
           ost=st
        ENDDO 

        print *, "Too many steps in simpson"
!        STOP
        s=1.
	
        END SUBROUTINE simpson

        SUBROUTINE trapz(R1,R2,U,a,b,z,n,s)
        implicit none
        INTEGER, INTENT(IN) :: n
        REAL, INTENT(IN)  :: R1, R2, U, a, b, z
        REAL, INTENT(INOUT) :: s

        INTEGER :: i, it
        REAL    :: del, x, ssum, tmp

        !print *, 'z=', z

        IF (z .LT. -1.e10) THEN
        IF (n .EQ. 1) THEN
            s=0.5*(b-a)*(innerint(R1,R2,U,a)+innerint(R1,R2,U,b))
        ELSE
            it=2**(n-2)
            tmp=it
            del=(b-a)/tmp
            x=a+0.5*del
            ssum=0.
            DO i=1,it
               ssum=ssum+innerint(R1,R2,U,x)
               x=x+del
            ENDDO
            s=0.5*(s+(b-a)*ssum/tmp)
        ENDIF
        ENDIF

        IF (z .GE. -1.e10) THEN
        IF (n .EQ. 1) THEN
            s=0.5*(b-a)*(func(R1,R2,U,z,a)+func(R1,R2,U,z,b))
        ELSE
            it=2**(n-2)
            tmp=it
            del=(b-a)/tmp
            x=a+0.5*del
            ssum=0.
            DO i=1,it
               ssum=ssum+func(R1,R2,U,z,x)
              ! print *, 'ssum=',ssum
               x=x+del
            ENDDO
            s=0.5*(s+(b-a)*ssum/tmp)
        ENDIF
        ENDIF
        
        END SUBROUTINE trapz


  SUBROUTINE init_module_wind_unir(id,config_flags,xlong,xlat,windfarm_initialized,&
                                            ims,ime,jms,jme,its,ite,jts,jte,ids,ide,jds,jde)
!
  IMPLICIT NONE
!
   integer ims,ime,jms,jme,ids,ide,jds,jde
   integer its,ite,jts,jte
   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
   TYPE (grid_config_rec_type) :: config_flags
   TYPE (PROJ_INFO) :: ts_proj
   logical :: windfarm_initialized
! 
   CHARACTER*256 num,input
   real lat,lon,ts_rx,ts_ry
   REAL :: known_lat, known_lon
   INTEGER i,j,nval,k,id
  
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
!
 
      IF ( wrf_dm_on_monitor() ) THEN
      
      flow_induction=1
      IF (config_flags%windfarm_ind .EQ. 0) THEN
        flow_induction=0
      ENDIF
      
      
      tiplossf= config_flags%windfarm_tiplossf

      tke_adv=.false.
      if (config_flags%bl_pbl_physics .EQ. 5 )then
      tke_adv = config_flags%bl_mynn_tkeadvect
      endif
      if ( tke_adv ) then
      call wrf_error_fatal('WRF-unir_drag_algr: Do not use the tke_adv option!')
      endif

      
! ... PAJ: Opens the file with the location of the wind turbines ...
!
        if ( config_flags%windfarm_ij .eq. 1 ) then
          open(70,file='windturbines-ij.txt',form='formatted',status='old')
        else
          open(70,file='windturbines.txt',form='formatted',status='old')
        end if
!
! ... PAJ: Counts the turbines ...
!
       nt=0
 10    read(70,*,end=100) 
       nt=nt+1
       goto 10
!
 100   continue
       rewind (70)
     END IF         ! if-monitor

      CALL wrf_dm_bcast_real(tiplossf,1)
      CALL wrf_dm_bcast_logical(tke_adv,1)
      CALL wrf_dm_bcast_integer(flow_induction,1)
      CALL wrf_dm_bcast_integer(nt,1)
!
! ... PAJ: Initializes the configuration of the wind farm(s) ...
!
     if (.not. windfarm_initialized) then
       allocate (nkind(nt),ival(nt,max_domains),jval(nt,max_domains))
       allocate (hubheight(nt),stc(nt),stc2(nt),cutin(nt),cutout(nt),diameter(nt),npower(nt))
       ival=-9999
       jval=-9999
       windfarm_initialized=.true.
     endif
!
     IF ( wrf_dm_on_monitor() ) THEN
     do k=1,nt
       if ( config_flags%windfarm_ij .eq. 1 ) then
         read(70,*) ival(k,id), jval(k,id), nkind(k)
         write(message_wind,*)'WINDFARM Turbine #',k,': I, J = ',ival(k,id), jval(k,id),'; Type = ',nkind(k)
         CALL wrf_message(message_wind)

       else

         read(70,*)lat,lon,nkind(k)
         write(message_wind,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon,'; Type = ',nkind(k)
         CALL wrf_message(message_wind)

         CALL map_init(ts_proj)

         known_lat = xlat(its,jts)
         known_lon = xlong(its,jts)

      ! Mercator
      IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj,               &
                      truelat1 = config_flags%truelat1, &
                      lat1     = known_lat,             &
                      lon1     = known_lon,             &
                      knowni   = REAL(its),             &
                      knownj   = REAL(jts),             &
                      dx       = config_flags%dx)

      ! Lambert conformal
      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
         CALL map_set(PROJ_LC, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      truelat2 = config_flags%truelat2,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = REAL(its),              &
                      knownj   = REAL(jts),              &
                      dx       = config_flags%dx)
!      ! Polar stereographic
      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = REAL(its),              &
                      knownj   = REAL(jts),              &
                      dx       = config_flags%dx)
!
      END IF
!
         CALL latlon_to_ij(ts_proj, lat, lon, ts_rx, ts_ry)
!
          ival(k,id)=nint(ts_rx)
          jval(k,id)=nint(ts_ry)
          if (ival(k,id).lt.ids.and.ival(k,id).gt.ide) then
            ival(k,id)=-9999
            jval(k,id)=-9999
          endif

          write(message_wind,*)'WINDFARM Turbine #',k,': Lat, lon = ',lat,lon, &
                               ', (i,j) = (',ival(k,id),',',jval(k,id),'); Type = ',nkind(k)
          CALL wrf_debug(0,message_wind)

     end if     ! ij-test
    

     enddo
      close(70)
!
! ... PAJ: Read the tables with the turbine's characteristics ...
!
         turbws=0.
         turbtc=0.
         turbpw=0.
         DO i=1,nt
          write(num,*) nkind(i)
          num=adjustl(num)
          input="wind-turbine-"//trim(num)//".tbl"
          OPEN(file=TRIM(input),unit=19,FORM='FORMATTED',STATUS='OLD')
          READ (19,*,ERR=132)nval
          READ(19,*,ERR=132)hubheight(i),diameter(i),stc(i),npower(i)
            DO k=1,nval
              READ(19,*,ERR=132)turbws(nkind(i),k),turbtc(nkind(i),k),turbpw(nkind(i),k)
            ENDDO
          cutin(i)  = turbws(nkind(i),1)
          cutout(i) = turbws(nkind(i),nval)
          stc2(i) = turbtc(nkind(i),nval)
          close (19)
         ENDDO


 132   continue

      endif             !  if-monitor

        CALL wrf_dm_bcast_integer(ival,nt*max_domains)
        CALL wrf_dm_bcast_integer(jval,nt*max_domains)
        CALL wrf_dm_bcast_real(hubheight,nt)
        CALL wrf_dm_bcast_real(diameter,nt)
        CALL wrf_dm_bcast_real(stc,nt)
        CALL wrf_dm_bcast_real(npower,nt)
        CALL wrf_dm_bcast_real(cutin,nt)
        CALL wrf_dm_bcast_real(cutout,nt)
        CALL wrf_dm_bcast_integer(nkind,nt) 
        CALL wrf_dm_bcast_real(turbws,maxvals*maxvals2) 
        CALL wrf_dm_bcast_real(turbtc,maxvals*maxvals2) 
        CALL wrf_dm_bcast_real(turbpw,maxvals*maxvals2)
        CALL wrf_debug(100,'Leaving windfarm_init')


  END SUBROUTINE init_module_wind_unir
  
END MODULE module_wind_unir
