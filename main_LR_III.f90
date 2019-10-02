      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!																   !!!!
      !!!!                    Author = dx Andxea Boghi					   !!!!
      !!!!																   !!!!
      !!!!                            LEGEND							   !!!!
      !!!!																   !!!!
      !!!!                        C(1) = UO_2^{+2}						   !!!!
      !!!!                        C(2) = Ca^{+2}						   !!!!
      !!!!                        C(3) = Cl^-							   !!!!
	  !!!!                        C(4) = L^-							   !!!!
      !!!!                        C(5) = H_3O^+							   !!!!
	  !!!!                        C(6) = HCO_3^-						   !!!!
	  !!!!                        C(7) = LH								   !!!!
	  !!!!                        C(8) = LUO_2^+						   !!!!
      !!!!                        C(9) = Ca_2UO_2(CO_3)_3				   !!!!
      !!!!                        C(10) = CaUO_2(CO_3)_3^{-2}			   !!!!
      !!!!                        C(11) = UO_2CO_3						   !!!!
      !!!!                        C(12) = (UO_2)_2CO_3(OH)_3^-						   !!!!
      !!!!																   !!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Program LO_RISE_Third
      implicit none
      real*8,allocatable :: ct0(:),gammaA(:,:),al(:,:),cl(:,:),cs(:,:),ct(:,:),I(:),x(:)
      real*8,allocatable :: ql(:),qs(:),qt(:),Dl(:),totCl(:,:),totCs(:,:),totCt(:,:)
	  real*8,allocatable :: Phi(:),z(:),nu(:,:),Ka(:),Bd0(:),Kd(:),Kd_old(:)
	  real*8,allocatable :: totCt_old(:,:),dgammadI(:,:),dKddt(:),dKddt_old(:),gammaA0(:,:)
	  real*8,allocatable :: dphidx(:),ct_star(:)
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  real*8,allocatable :: Cul(:),Cul_old(:),SOH(:),dcsdul(:)
	  real*8,allocatable :: Scul_new(:),Scul_old(:),dUsdt_old(:)
      real*8,allocatable :: cuss(:),Scuss_new(:),Scuss_old(:),cusf(:),dUsdt(:)
	  real*8,allocatable :: cusf_old(:),cuss_old(:)
	  integer n,m,nx,nt,j,k,it,iit,nit1
	  real*8 ctU0,ctCl0,theta,rho,Lx,CFL,TC,f,res
	  real*8 rhow,nuw,V_REF,dx,dt,Td,time,To,dx2,Pe,rc,Ad,h,gammaAREF,p_atm,eps0
	  real*8 bhs,pH,kusf,kusb,ctL0
	  real*8 Du,loge,v,totHS_eq,FUT,FU1,FU2,FU3,factor
	  real*8 pCO2,SOH0,nSOH,ASOH,pH0,Up,U0,Un,Up_old,U0_old,Un_old
	  integer n_d,n_h,n_m,n_s,ind,indU
																	
      !!!character(len=80)::fname		 

	  !!!!!!!!!!! CONSTANTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      n = 12 !!!!! Number of species
      m = 5 !!!! Number of independent species
      p_atm = 1.d0 !!!! atmosp_CO2eric pressure in atm
        
      open(1,file='imp_lr_3.dat')
      read(1,*) ctU0             !!! 1d-06 !!!! ctU0 [mol/l]
      read(1,*) nSOH             !!! 1.92e-06
      read(1,*) ASOH             !!! 3e+03
      read(1,*) pH            !!! 7
      read(1,*) bHS              !!! 1.d-3  !!!! [mol/(l pH)]
      read(1,*) pCO2          !!! 3.4d-03 [atm]	
      read(1,*) theta            !!! 1.d0 !!! porosity  
      read(1,*) rho              !!! 0.d0 !!!! terrain density kg/l
      read(1,*) Lx             !!! 6.d0 !!!! Lx
      read(1,*) nx               !!! 9 !!!! nx	
      read(1,*) CFL              !!! 9 !!!! nx											
      read(1,*) To               !!! To [s] 
      read(1,*) TC               !!! Temp	
      read(1,*) f                !!! impedence fatcor
      read(1,*) v               !!! velocity
      close(1)

      ctCl0 = theta*20.d-03
      ctL0 = theta*1.d-04;        

      print*, ' ctCl0 =', ctCl0, ' [M] ' 
      print*, ' ctU0 =', ctU0, ' [M], pCO2 =', pCO2, ' [atm] '
      print*, ' theta =', theta, ' rho =', rho, ' [kg/l], T =', TC, ' [C] ' 
      print*, ' Lx =', Lx, ' [dm]'
      print*, ' v =', v, ' [dm/s] f =', f
      print*, ' nx =', nx, ' To =',To, ' [s], CFL =',CFL

      allocate(ct0(1:m),Dl(1:n),ql(1:nx),qs(1:nx),qt(1:nx),totCl(1:m,1:nx),totCs(1:m,1:nx),totCt(1:m,1:nx))
      allocate(gammaA(1:n,1:nx),al(1:n,1:nx),cl(1:n,1:nx),cs(1:n,1:nx),ct(1:n,1:nx),I(1:nx),x(1:nx))
	  allocate(Phi(1:nx),z(1:n),nu(1:m,1:n),Ka(1:n),Bd0(1:n),totCt_old(1:m,1:nx),dgammadI(1:n,1:nx))
	  allocate(dphidx(1:nx),ct_star(1:m))
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  allocate(Cul(1:nx),Cul_old(1:nx),SOH(1:nx),dcsdul(1:n))
	  allocate(Scul_new(1:nx),Scul_old(1:nx))
      allocate(cuss(1:nx),cuss_old(1:nx),Scuss_new(1:nx),Scuss_old(1:nx),cusf(1:nx))
	  allocate(cusf_old(1:nx))
	  allocate(Kd(1:nx),dKddt(1:nx),Kd_old(1:nx),dKddt_old(1:nx))
	  allocate(gammaA0(1:n,1:nx),dUsdt_old(1:nx),dUsdt(1:nx))

      indU = 10
      dcsdul(1)= 2.298514d-03
      dcsdul(indU)= 7.336760d+19
      kusf = 6.948103d-06
      kusb = 9.392974d-07

      !!! r
      dx = Lx/(nx-1)
	  do j=1,nx
	     x(j) = -Lx/2.d0 +dx*(j-1)
 	  enddo

      !!! Fixed Parameters
      call FixedParameters(n,m,pCO2,p_atm,TC,Ka,Ad,Bd0,gammaAREF,h,rhow,nuw,V_REF,z,nu,Dl)
      pH0 = 9.5d0
      
      !!!!!!!!!!!!!!!!!!!!!
      !!!!! define total initial concentration
      ct0(1) = ctU0
      ct0(2) = 0.d0
      ct0(3) = ctCl0
      ct0(4) = ctL0
      !ct0(7) = -(z(1)*ct0(1) +z(2)*ct0(3)/2.d0 +z(3)*ct0(3))

      SOH0 = nSOH*ASOH !!! [mol/kg]

      !!! Initrialization
      j = 1
      gammaA(:,j) = 1.d0
      al(:,j) = 1.d-60
      I(j) = 3.d-2 !(z(1)*z(1)*ct0(1) +z(2)*z(2)*ct0(3)/2.d0 +z(3)*z(3)*ct0(3))	
	
      call SpeciationII(pH,al(:,j),ct(:,j),cl(:,j),cs(:,j),dcsdul,SOH(j),SOH0, &
	  gammaA(:,j),dgammadI(:,j),I(j),ct0,nu,z,theta,rho,Ka,Bd0,Ad,gammaAREF,h,pH0,bHS,n,m,indU)

 	  !!! calculate Charge
	  call calculateCharge(ql(j),qs(j),qt(j),z,cl(:,j),cs(:,j),n,theta,rho)
      !!! calculate Total Concentration
	  call calculateTotalConc(totCl(:,j),totCs(:,j),totCt(:,j),cl(:,j),cs(:,j),ct(:,j),nu,n,m)

      print*, ' totCt(1,j)=',totCt(1,j),' totCl(1,j)=',totCl(1,j),' totCs(1,j)=',totCs(1,j),' Kd0=',totCs(1,j)/totCl(1,j)
      print*, ' SOH(j)=',SOH(j)
	  do j=1,nx/2
		 do k=1,n
			al(k,j) = al(k,1)
			cl(k,j) = cl(k,1)
			cs(k,j) = cs(k,1)
			ct(k,j) = ct(k,1)
			gammaA(k,j) = gammaA(k,1)
			dgammadI(k,j) = dgammadI(k,1)
		 enddo
		 do k=1,m
			totCl(k,j) = totCl(k,1)
			totCs(k,j) = totCs(k,1)
			totCt(k,j) = totCt(k,1)
		 enddo
		 I(j) = I(1)
		 ql(j) = ql(1)
		 qs(j) = qs(1)
		 qt(j) = qt(1)
		 SOH(j) = SOH(1)
	  enddo
         
      !!!!! define total initial concentration
      ct0(1) = ctU0*1.d-06

      !!! Initrialization
      j = nx
      gammaA(:,j) = 1.d0
      al(:,j) = 1.d-60
      I(j) = 3.d-2 !(z(1)*z(1)*ct0(1) +z(2)*z(2)*ct0(3)/2.d0 +z(3)*z(3)*ct0(3))	
	
      call SpeciationII(pH,al(:,j),ct(:,j),cl(:,j),cs(:,j),dcsdul,SOH(j),SOH0, &
	  gammaA(:,j),dgammadI(:,j),I(j),ct0,nu,z,theta,rho,Ka,Bd0,Ad,gammaAREF,h,pH0,bHS,n,m,indU)

 	  !!! calculate Charge
	  call calculateCharge(ql(j),qs(j),qt(j),z,cl(:,j),cs(:,j),n,theta,rho)
      !!! calculate Total Concentration
	  call calculateTotalConc(totCl(:,j),totCs(:,j),totCt(:,j),cl(:,j),cs(:,j),ct(:,j),nu,n,m)

      print*, ' totCt(1,j)=',totCt(1,j),' totCl(1,j)=',totCl(1,j),' totCs(1,j)=',totCs(1,j),' Kd0=',totCs(1,j)/totCl(1,j)
      print*, ' SOH(j)=',SOH(j)
	  do j=nx/2+1,nx
		 do k=1,n
			al(k,j) = al(k,nx)
			cl(k,j) = cl(k,nx)
			cs(k,j) = cs(k,nx)
			ct(k,j) = ct(k,nx)
			gammaA(k,j) = gammaA(k,nx)
			dgammadI(k,j) = dgammadI(k,nx)
		 enddo
		 do k=1,m
			totCl(k,j) = totCl(k,nx)
			totCs(k,j) = totCs(k,nx)
			totCt(k,j) = totCt(k,nx)
		 enddo
		 I(j) = I(nx)
		 ql(j) = ql(nx)
		 qs(j) = qs(nx)
		 qt(j) = qt(nx)
		 SOH(j) = SOH(nx)
	  enddo
         
	  totHS_eq = totCt(5,nx)

      print*, ' pCO2 =', pCO2, ' [atm], pH =', -dlog10(al(7,nx))
	  
	  !!! New parameters
	  loge = dlog10(dexp(1.d0))
	  Du = theta*f*Dl(1)
	  !print*, ' Du =', Du, ' [dm2/s]'
	  !pause
	  
	  dt = CFL*dx*dx/Du

          nt = 1 +floor(To/dt)

          call calculateTime(n_d,n_h,n_m,n_s,To)
	        print *, 'TIME'
	        print *, N_D,' DAYS',N_H,' HOURS',N_M,' MINUTES',N_S,' SECONDS'

	  print*, ' dt=', dt,' [s], To =', To, ' [s]'
          print*, ' Du=', Du,' [dm2/s]'

	  Td = Lx*Lx/Du

	  !!!! NON -DIMENSIONALIZE Dl
	  Dl = f*Dl/Du
	  x = x/Lx
	  dx = dx/Lx
	  dx2 = dx*dx
	  dt = dt/Td
	  Pe = Lx*v/Du
	  rc = dt/(2.d0*dx2)
          kusf = kusf*Td
          kusb = kusb*Td

          !!! Because of the different units
          factor = 10.d0

	  print*, ' pH=',pH
	  do j=1,nx
	     Cul(j) = totCl(1,j)
	     Cusf(j) = (kusb/(kusf +kusb))*totCs(1,j) !/(1.d0 +kusf/kusb)
	     !Cuss(j) = (kusf/(kusf +kusb))*totCs(1,j)
             Cuss(j) = (kusf/kusb)*Cusf(j)
	  enddo
	  !!!kusf = 0.d0
	  !!!kusb = 0.d0
	  !!!Cuss = 0.d0

          call calcKdCusf(Kd,Cusf,Cul,pH,al,gammaA,SOH0,dcsdul,Ka,kusf,kusb,n,nx,indU)
	 
	  Up = 0.d0
	  U0 = 0.d0
	  Un = 0.d0 
	  !call percCharge(Up,U0,Un,cl(:,1),nu(1,:),z,n)
	  !print*, ' Up =', Up,' U0 =', U0,' Un =', Un,' Cul(1) =', Cul(1)

	  !pause

	  Cul_old =	Cul
      Cuss_old = Cuss
      Cusf_old = Cusf
	  Kd_old = Kd
          dKddt_old = 0.d0
          dUsdt_old = 0.d0
      call createData(pH,Cul,Cuss,Cusf,Kd,SOH,totCt,pCO2,ct,al,x,I,qt,ql,qs,0,n,m,nx,Lx)

	  !do it=1,nt
	  !   if (mod(it,nt/10).eq.0) then
	!	    print*, ' it=',it
!		    ind = 10*it/nt
!		    print *, '10*it/nt=',ind
!		 endif
!	  enddo

	  !goto 345

	  Up_old = Up
	  U0_old = U0
	  Un_old = Un

	  open(27,file='fluxes.csv',form='formatted')
      write(27,*) '"t","FUb","FUi"'

	  !!!! TIME LOOP
	  nit1 = 4
	  do it=1,nt
	     !!!!  calculate Source old
         call sourceOld(Scuss_old,Scul_old,Cul_old,Cuss_old,Cusf_old,Kd_old,dKddt_old,dUsdt_old, &
                        theta,rho,Pe,kusf,kusb,dt,dx,dx2,Up_old,U0_old,Un_old,nx)

	  	 iit=1
         res=1.d0
         eps0 = 1.d-6
         !do while((res.ge.eps0).and.(iit.le.nit1))
         do iit=1,nit1
			call calcKdCusf(Kd,Cusf,Cul,pH,al,gammaA,SOH0,dcsdul,Ka,kusf,kusb,n,nx,indU)

		 	do j=1,nx
			   dKddt(j) = (Kd(j) -Kd_old(j))/dt
                           dUsdt(j) = 2.d0*rho*(cuss(j) -cuss_old(j) +cusf(j) -cusf_old(j))/(bHS*dt) 
                        enddo

		 !do iit=1,3
			!!!! update for convergence
		    ct_star = cul

			call sourceNew(Scuss_new,Scul_new,Cuss,Cusf,Kd,dUsdt,theta,rho,kusf,kusb,dx,dt,nx)

			!!! solve diffusion
            call solveDiffusion(Cul,Cuss,Scul_new,Scul_old,Scuss_old,Scuss_new, &
			     Kd,dKddt,Pe,theta,rho,kusf,kusb,dt,dx,dx2,Up,U0,Un,nx)

			do j=1,nx
			   totCt(1,j) = theta*Cul(j) +rho*(Cusf(j) +Cuss(j))
			   totCt(2,j) = -(z(1)*totCt(1,j) +z(3)*totCt(3,j) +z(4)*totCt(4,j) &
			   +z(5)*(totCt(5,j) +rho*bHS*(pH -pH0)))/z(2)
            enddo


			!!! Calculate Residual
			res = dabs(cul(1)/ct_star(1) -1.d0)
			!!! update iteration
			!iit=iit+1
		 enddo

		 print *, 'res=',res,' iit=',iit
		 
		 !!!!! SAVE DATA !!!!!!
	     if (mod(it,nt/10).eq.0) then

		    !!! calculate Charge
			do j=1,nx
			   call SpeciationII(pH,al(:,j),ct(:,j),cl(:,j),cs(:,j),dcsdul,SOH(j),SOH0, &
	           gammaA(:,j),dgammadI(:,j),I(j),totCt(:,j),nu,z,theta,rho,Ka,Bd0,Ad,gammaAREF,h,pH0,bHS,n,m,indU)

			   call calculateCharge(ql(j),qs(j),qt(j),z,cl(:,j),cs(:,j),n,theta,rho)
            enddo
                      !!! Unifor Aciidity Because it is influenced by uranium acid

		    ind = 10*it/nt +1

		    print *, '10*it/nt +1=',ind
      
	        call createData(pH,Cul,Cuss,Cusf,Kd,SOH,totCt,pCO2,ct,al,x,I,qt,ql,qs,ind,n,m,nx,Lx)
	        print *,'****************************'
            print *, 'lacking iterations=',nt -it

		    time = dt*Td*it
             
		    call calculateTime(n_d,n_h,n_m,n_s,time)

	        print *, 'TIME'
	        print *, N_D,' DAYS',N_H,' HOURS',N_M,' MINUTES',N_S,' SECONDS'

		 endif

		 FUT = (Cul(nx/2+1) -Cul(nx/2-1))/(2.d0*dx) -Pe*Cul(nx/2)

		 call IntegrateL(FU1,cul_old,cul,dt,dx,nx/2)
		 call IntegrateL(FU2,cusf_old,cusf,dt,dx,nx/2)
		 call IntegrateL(FU3,cuss_old,cuss,dt,dx,nx/2)

		 write(27,103) Td*dt*it,',',(Du/Lx)*FUT,',',(Du/Lx)*(theta*FU1+rho*(FU2+FU3))

		 !!!!!!!
	     Cul_old =	Cul
         Cuss_old = Cuss
         Cusf_old = Cusf
		 Kd_old = Kd
		 dKddt_old = dKddt
		 dUsdt_old = dUsdt
	     Up_old = Up
	     U0_old = U0
	     Un_old = Un

	  enddo

  	  close(27)

	  !345 continue

	  open(10,file='time.dat')
	  write(10,*) 'TIME'
      write(10,*) N_D,'DAYS'
      write(10,*) N_H,'HOURS'
	  write(10,*) N_M,'MINUTES'
	  write(10,*) N_S,'SECONDS' 					
      close(10)

  	  103 format(e18.10,A,e18.10,A,e18.10)


      end


      !!!*********************************************************************************
      !!!*						               			 * 
      !!!*                             FixedParameters 	        		     	 *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine FixedParameters(n,m,pCO2,p_atm,TC,Ka,Ad,Bd0,gammaAREF,h,rhow,nuw,V_REF,z,nu,D)
	  implicit none
	  integer, intent(in) :: n,m
      real*8, intent(in) :: TC,pCO2,p_atm
      real*8, intent(out) :: Ka(1:n),Ad,Bd0(1:n),gammaAREF,h,rhow,nuw,V_REF,D(1:n)
	  real*8, intent(out) :: nu(1:m,1:n),z(1:n)
	  real*8 logK(1:n),TK,kb,e,NA,epsilon0,epsilonx,Bd,pi,a0(1:n),nu_CO3(1:n),K12
	  integer i

	  !!! electric charge
	  z = (/+2,+2,-1,-1,+1,-1, 0,+1, 0,-2, 0,-1/)

      !!!! stoichiometric coefficients
      nu(1,:) = (/1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2/)
      nu(2,:) = (/0, 1, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0/)
      nu(3,:) = (/0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
      nu(4,:) = (/0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0/)
      nu(5,:) = (/0, 0, 0, 0, 1,-1, 1, 0,-6,-6,-2,-3/)
      nu_CO3  = (/0, 0, 0, 0, 0, 1, 0, 0, 3, 3, 1, 1/)

      pi = 4.d0*datan(1.d0)
 
      !!! Mass Action Constants (Activity)
      logK = (/0.0, 0.0, 0.0, 0.0, 0.0, 10.329, 5.00, 7.01, 30.7, 27.18, 9.94, -0.86/)

      Ka = 10.d0**logK
	  
	  !!!! K12
	  K12 = 10.d0**18.149d0			

	  !!! press
	  Ka = Ka*((pCO2/p_atm)/K12)**nu_CO3
	  !!! OH-
	  !!!!Ka(12) = Ka(12)*10^(-13.997);
      !!! Temperature (/K/)
      TK = TC +273.15
 
      !!! Boltzmann constant (/m2 kg s-2 K-1/)
      kb = 1.38064852d-23
 
      !!! Electron Charge (/C/)
      e = 1.6021766d-19
 
      !!! Avogadxo Number (/mol**-1/)
      NA = 6.022141d+23
 
      !!! Electric Permittivity of vacuum (/s**2  C**2  m**?3  kg**?1/)
      epsilon0 = 8.854187817d-12
 
      !!! Relative Electric Permittivity (/1/)
      epsilonx = 1.d0
 
      !!! Debye�H�ckel B (/mol**-1/2 m**1/2/)
      Bd = e*dsqrt((2.d0*NA)/(epsilon0*epsilonx*kb*TK))
 
      !!! Debye�H�ckel Ad (/mol**-1/2 m**3/2/) 
      Ad = e*e*Bd/(dlog(10.d0)*8*pi*epsilon0*epsilonx*kb*TK)
 
      Bd0 = Bd0*(10.d0)**1.5d0 !!! convert to liters
      Ad = Ad*(10.d0)**1.5d0 !!! convert to liters
  
      !!!! gammaA_REF (/l/mol/)
      gammaAREF = 1.d0
 
      !!!!! Davies coeff
      h = 0.24d0
 
      !!!! empjical
      Ad = 0.5d0
	  Bd0 = 1.d0
      h = 0.3d0 

      !!!! water density kg/l
      rhow = 1.d0

      !!!! water viscosity dm2/s
      nuw = 1.d-04

      !!! Ref pot
      V_REF = kb*TK/e

      !!! Debye�H�ckel a0
      a0(1) = 2.8734d-10 !!! m	  UO_2^{+2}
      a0(2) = 2.7643d-10 !!! m	  Ca^{+2}
      a0(3) = 1.0758d-10 !!! m	  Cl^-
      a0(5) = 2.3454d-11 !!! m	  H_3O^+
	  a0(6) = 1.8351d-10 !!! m	  HCO3-
      a0(9) = 4.7474d-10 !!! m	  Ca2UO2(CO3)3 (aq)
      a0(10) = 4.2820d-10 !!! m   CaUO2(CO3)3-2
      a0(11) = 3.2594d-10 !!! m	  UO2CO3 (aq) 
      a0(12) = 2.8734d-10 !!! m	  UO2OH+

      a0(4) = a0(2)
	  a0(7:8) = a0(2)

      a0 = a0*10.d0 !!! dm
      kb = kb*100.d0 !!! dm2

      !!!! Diffusivity   
      do i=1,n
         D(i) = kb*TK/(6.d0*pi*rhow*nuw*a0(i))
		 print*, ' D(i)=', D(i),' [dm2/s]'
      enddo
	!  pause
      !D(5) = 8.4d-07  !!! [dm2/s]
      !D(6) = 1.08d-07	!!! [dm2/s]

 
      end

      !!!*********************************************************************************
      !!!*						               			 * 
      !!!*                             calcConc                 		     	 *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine calcConc(cl,cs,ct,dcsdul,SOH,theta,rho,gammaA,al,SOH0,pH0,bHS,n,indU)
      implicit none
      integer, intent(in) :: n,indU
      real*8, intent(in) :: al(1:n),gammaA(1:n),dcsdul(1:n),theta,rho,SOH0,pH0,bHS
      real*8, intent(out) ::cl(1:n),cs(1:n),ct(1:n),SOH
      integer i
 
      !!!! calculate liquid concentration
      do i=1,n
         cl(i) = al(i)/gammaA(i)
      enddo

      !!!! calculate concentration in solids
      call Sorption(cs,dcsdul,SOH,al,SOH0,pH0,bHS,n,rho,indU)

      !!!! calculate total concentration
      do i=1,n
         ct(i) = theta*cl(i) +rho*cs(i)
      enddo
 
      end

      !!!*********************************************************************************
      !!!*						               			 * 
      !!!*                             Sorption                 		     	 *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine Sorption(cs,dcsdul,SOH,al,SOH0,pH0,bHS,n,rho,indU)
      implicit none
      integer, intent(in) :: n,indU
      real*8, intent(in) :: al(1:n),SOH0,pH0,bHS,rho,dcsdul(1:n)
      real*8, intent(out) ::cs(1:n),SOH

	  cs = 0.d0
   
      SOH = SOH0/(1.0 +dcsdul(1)*al(1)/al(5)**2.d0 +dcsdul(indU)*al(indU)*al(5)**2.d0);

      cs(1) = SOH*dcsdul(1)*al(1)/al(5)**2.d0;

      cs(indU) = SOH*dcsdul(indU)*al(indU)*al(5)**2.d0;

      cs(5) = bHS*(log10(al(5))+pH0)/rho
     
	  end

      !!!*********************************************************************************
      !!!*						               			 * 
      !!!*                                calculateIonic        	     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine calculateIonic(I,gammaA,dgammadI,n,z,cl,Ad,Bd0,h,gammaAREF)
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: z(1:n),cl(1:n),Ad,Bd0(1:n),h,gammaAREF
      real*8, intent(out) ::I,gammaA(1:n),dgammadI(1:n)
      real*8 sum,g,dgdI,alpha
      integer j
      
      sum = 0.d0
      do j=1,n
          sum = sum +(z(j)**2.d0)*cl(j)
      enddo

      g= I -0.5d0*sum

      sum = 0.d0
      do j=1,n
         sum = sum +(z(j)**4.d0)*cl(j)
      enddo

      dgdI = 1.d0 -0.5d0*sum*Ad*dlog(10.d0)*( 1.d0/(2.d0*dsqrt(dabs(I))*(1.d0 +Bd0(1)*dsqrt(dabs(I)))**2.d0) -h)
      I = I -g/dgdI
      I = dabs(I)
    
      !!!! l/[mol]
      alpha = 0.1d0
      gammaA = 1.d0
      do j=1,n
         if(z(j) == 0.d0)  then
            gammaA(j) = gammaAREF*10.d0**(alpha*I)
			dgammadI(j) = dlog(10.d0)*alpha*gammaA(j)
         else
            gammaA(j) = gammaAREF*10.d0**( -Ad*(z(j)**2.d0)*( dsqrt(dabs(I))/(1+Bd0(j)*dsqrt(dabs(I))) -h*I) )
			dgammadI(j) = -dlog(10.d0)*Ad*(z(j)**2.d0)*(1.d0/(2.d0*(dsqrt(dabs(I)))*(1+Bd0(j)*dsqrt(dabs(I)))**2.d0) -h)*gammaA(j)
         endif
      enddo
 
      end

      !!!*********************************************************************************
      !!!*						               											 * 
      !!!*                             ActivityBlock        			    	         *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine ActivityBlock(al,Ka,nu,n,m)
      implicit none
      integer, intent(in) :: n,m
      real*8, intent(in) :: nu(1:m,1:n),Ka(1:n)
      real*8, intent(inout) ::al(1:n)
      integer i,k

	  do i=m+1,n
         al(i) = Ka(i)
         do k=1,m
            al(i) = al(i)*(al(k)**nu(k,i))
         enddo
      enddo
 
      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                             calcCharge                 		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine calculateCharge(ql,qs,qt,z,cl,cs,n,theta,rho)
 	  implicit none
	  integer, intent(in) :: n
      real*8, intent(in) :: cl(1:n),cs(1:n),z(1:n),theta,rho
      real*8, intent(out) ::ql,qs,qt
	  integer j
 
      ql=0.d0
      do j=1,n
         ql = ql +z(j)*theta*cl(j)
      enddo
 
      qs=0
      do j=1,4
         qs = qs +(z(j)-2.d0)*rho*cs(j)
      enddo
	  do j=6,n-1
         qs = qs +(z(j)-2.d0)*rho*cs(j)
      enddo
      qs = qs +(z(12)-1.d0)*rho*cs(12)

      qt = ql + qs
 
      end
	  
	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                        calculate Total Conc             		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine calculateTotalConc(xl,xs,xt,cl,cs,ct,nu,n,m)
 	  implicit none
	  integer, intent(in) :: n,m
      real*8, intent(in) :: cl(1:n),cs(1:n),ct(1:n),nu(1:m,1:n)
      real*8, intent(out) ::xl(1:m),xs(1:m),xt(1:m)
	  integer i,j
	   
      xl = 0.d0
      do i=1,m
         do j=1,n
            xl(i) = xl(i) +nu(i,j)*cl(j)
         enddo
      enddo
 
      xs = 0.d0
      do i=1,m
         do j=1,n
            xs(i) = xs(i) +nu(i,j)*cs(j)
         enddo
      enddo
 
      xt = 0.d0
      do i=1,m
         do j=1,n
            xt(i) = xt(i) +nu(i,j)*ct(j)
         enddo
      enddo
 
      end
	  
	  !!!*********************************************************************************
      !!!*						               											 * 
      !!!*                             Calculate Time			        		     	 *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine calculateTime(n_d,n_h,n_m,n_s,time)
	  implicit none
	  real*8, intent(in) :: time
	  integer, intent(out) :: n_d,n_h,n_m,n_s
	  real*8 day,hour,minute,second

      !!! Time Parameters
      day = 24.d0*3600.d0
      hour = 3600.d0
      minute = 60.d0
      second = 1.d0

	  n_d = floor(time/day)
      n_h = floor((time -n_d*day)/hour)
      n_m = floor((time -n_d*day -n_h*hour)/minute)
      n_s = floor((time -n_d*day -n_h*hour -n_m*minute)/second)

	  end
	  
	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                          Fjst Oder Derivative        		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine fjstOrdDer(dcldx,cl,dx,nx)
 	  implicit none
	  integer, intent(in) :: nx
      real*8, intent(in) :: cl(1:nx),dx
      real*8, intent(out) ::dcldx(1:nx)
	  integer j

      dcldx(1) = (-3.d0*cl(1) +4.d0*cl(2) -cl(3))/(2.d0*dx)
      do j=2,nx-1
         dcldx(j) = (cl(j+1)-cl(j-1))/(2.d0*dx)
      enddo
      dcldx(nx) = (3.d0*cl(nx) -4.d0*cl(nx-1) +cl(nx-2))/(2.d0*dx)
 
      end

      !!!*********************************************************************************
      !!!*						          	                 * 
      !!!*                         Calculate Source Old	        		     	 *
      !!!*                                                                               *
      !!!*********************************************************************************
      subroutine sourceOld(Scuss,Scul,Cul,Cuss,Cusf,Kd,dKddt,dUsdt, &
	             theta,rho,Pe,kusf,kusb,dt,dx,dx2,Up,U0,Un,nx)
 	  implicit none
	  integer, intent(in) :: nx
      real*8, intent(in) :: Cul(1:nx),dt,dx,dx2,Up,U0,Un
	  real*8, intent(in) :: Cuss(1:nx),Cusf(1:nx),Kd(1:nx),dKddt(1:nx),dUsdt(1:nx)
	  real*8, intent(in) ::	Pe,kusf,kusb,theta,rho
	  real*8, intent(out) :: Scuss(1:nx),Scul(1:nx)
	  real*8 betaU,aule,aulw,aulp,Scul_s(1:nx)
	  integer j

	  do j=1,nx
		 Scuss(j) = Cuss(j) +(dt/2.d0)*(kusf*Cusf(j) -kusb*Cuss(j))
		 Scul_s(j) = (dt/2.d0)*rho*Cuss(j)*kusb/(theta +rho*(kusb/(kusf+kusb))*Kd(j))
	  enddo
 
	  call coeffU(betaU,aule,aulw,Kd(1),dKddt(1),Pe,theta,rho,kusf,kusb,dt,dx,dx2)
	  aulp = (1.d0 -aule -aulw -betaU)
	  Scul(1) = (aulp -aulw*(2.d0*dx*Pe))*Cul(1)+(aule +aulw)*Cul(2) +Scul_s(1)

	  do j=2,nx-1
	  	 call coeffU(betaU,aule,aulw,Kd(j),dKddt(j),Pe,theta,rho,kusf,kusb,dt,dx,dx2)
		 aulp = (1.d0 -aule -aulw -betaU)
	  	 Scul(j) = aulp*Cul(j) +aule*Cul(j+1) +aulw*Cul(j-1) +Scul_s(j) 
	  enddo

	  call coeffU(betaU,aule,aulw,Kd(nx),dKddt(nx),Pe,theta,rho,kusf,kusb,dt,dx,dx2)
	  aulp = (1.d0 -aule -aulw -betaU)
	  Scul(nx) = (aulp +aule*(2.d0*Pe*dx))*Cul(nx) +(aule +aulw)*Cul(nx-1) +Scul_s(nx) 

      end

      !!!*********************************************************************************
      !!!*						                	         * 
      !!!*                         Calculate Source New	        		         *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine sourceNew(Scuss,Scul,Cuss,Cusf,Kd,dUsdt,theta,rho,kusf,kusb,dx,dt,nx)
 	  implicit none
	  integer, intent(in) :: nx
      real*8, intent(in) :: Cuss(1:nx),Cusf(1:nx),Kd(1:nx)
	  real*8, intent(in) :: dx,dt,kusf,kusb,theta,rho,dUsdt(1:nx)
	  real*8, intent(out) :: Scuss(1:nx),Scul(1:nx)
	  integer j
 
	  do j=1,nx
		 Scul(j) = (dt/2.d0)*rho*Cuss(j)*kusb/(theta +rho*(kusb/(kusf+kusb))*Kd(j))
                 Scuss(j) = (dt/2.d0)*kusf*Cusf(j)
	  enddo

      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                         Solve Diffusion		        		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************
      subroutine solveDiffusion(Cul,Cuss,Scul_new,Scul_old, &
                 Scuss_old,Scuss_new,Kd,dKddt,Pe, &
                 theta,rho,kusf,kusb,dt,dx,dx2,Up,U0,Un,nx)
 	  implicit none
	  integer, intent(in) :: nx
	  real*8, intent(in) :: Scul_new(1:nx),Scul_old(1:nx),Up,U0,Un
	  real*8, intent(in) :: Scuss_old(1:nx),Scuss_new(1:nx),Kd(1:nx),dt,dx,dx2
	  real*8, intent(in) :: Pe,dKddt(1:nx),theta,rho,kusf,kusb
	  real*8, intent(inout) :: Cul(1:nx),Cuss(1:nx)
	  real*8 betaU,aule,aulw,aulp
	  integer j

	  !!!! Cul

	  call coeffU(betaU,aule,aulw,Kd(1),dKddt(1),Pe,theta,rho,kusf,kusb,dt,dx,dx2)
	  aulp = (1.d0 +aule +aulw +betaU)
	  Cul(1) = ((aule +aulw)*Cul(2) +Scul_new(1) +Scul_old(1))/(aulp +aulw*(2.d0*Pe*dx))

	  do j=2,nx-1
	  	 call coeffU(betaU,aule,aulw,Kd(j),dKddt(j),Pe,theta,rho,kusf,kusb,dt,dx,dx2)
         aulp = (1.d0 +aule +aulw +betaU)
	     Cul(j) = (aule*Cul(j+1) +aulw*Cul(j-1) +Scul_new(j) +Scul_old(j))/aulp
	  enddo

	  call coeffU(betaU,aule,aulw,Kd(nx),dKddt(nx),Pe,theta,rho,kusf,kusb,dt,dx,dx2)
      aulp = (1.d0 +aule +aulw +betaU)
	  Cul(nx) = ((aule +aulw)*Cul(nx-1) +Scul_new(nx) +Scul_old(nx))/(aulp -aule*(2.d0*Pe*dx))
	  
      do j=1,nx
	     Cuss(j) = (Scuss_old(j) +Scuss_new(j))/(1.d0 +(dt/2.d0)*kusb)
      enddo

      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                          Second Order Derivative        		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine secondOrdDer(d2cldx2,cl,dx2,nx)
 	  implicit none
	  integer, intent(in) :: nx
      real*8, intent(in) :: cl(1:nx),dx2
      real*8, intent(out) ::d2cldx2(1:nx)
	  integer j

      d2cldx2(1) = (2.d0*cl(1) -5.d0*cl(2) +4.d0*cl(3) -cl(4))/dx2
      do j=2,nx-1
         d2cldx2(j) = (cl(j+1) -2.d0*cl(j) +cl(j-1))/dx2
      enddo
      d2cldx2(nx) = (2.d0*cl(nx) -5.d0*cl(nx-1) +4.d0*cl(nx-2) -cl(nx-3))/dx2
 
      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                               Integrate              		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine Integrate(phi,dphi,dx,nx)
 	  implicit none
	  integer, intent(in) :: nx
      real*8, intent(in) :: dphi(1:nx),dx
      real*8, intent(out) ::phi(1:nx)
	  integer j

      phi(nx) = 0.d0
	  phi(nx-1) = phi(nx) -dx*(dphi(nx-1) +dphi(nx))/2.d0
      do j= nx-2,1,-1
         phi(j) = phi(j+2) -2.d0*dx*dphi(j+1)
      enddo
 
      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                             Create Data  	        		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine createData(pH,Cul,Cuss,Cusf,Kd,SOH,totCt,pCO2,ct,al,x,I,qt,ql,qs,ind,n,m,nx,Lx)
 	  implicit none
	  integer, intent(in) :: n,m,nx,ind
      real*8, intent(in) :: x(1:nx),ct(1:n,1:nx),al(1:n,1:nx),I(1:nx),Lx
	  real*8, intent(in) :: pH,pCO2,qt(1:nx),ql(1:nx),qs(1:nx),SOH(1:nx)
	  real*8, intent(in) :: Cul(1:nx),Cuss(1:nx),Cusf(1:nx),Kd(1:nx),totCt(1:m,1:nx)
	  character(len=80) title_I,title_II,title_R
	  character(len=80) title_cha,title_C0,title_Act,title_Act_I,title_Act_II,title_D1,title_D2,title_D3,title_D4
	  integer j

	  write(title_R,301)'rest_',ind,'.csv'
      301 format(a5,i2.2,a4)
	  write(title_I,201)'ct_first_',ind,'.csv'
      201 format(a9,i2.2,a4)
      write(title_II,202)'ct_second_',ind,'.csv'
      202 format(a10,i2.2,a4)
	  write(title_cha,204)'charge_',ind,'.csv'
      204 format(a7,i2.2,a4)
	  write(title_C0,205)'C0_',ind,'.csv'
      205 format(a3,i2.2,a4)
	  write(title_Act,201)'Activity_',ind,'.csv'
	  write(title_Act_I,201)'Ac_first_',ind,'.csv'
          write(title_Act_II,202)'Ac_second_',ind,'.csv'
	  write(title_D1,205)'D1_',ind,'.csv'
	  write(title_D2,205)'D2_',ind,'.csv'
	  write(title_D3,205)'D3_',ind,'.csv'
	  write(title_D4,205)'D4_',ind,'.csv'

	  open(61,file=title_R,form='formatted')
      write(61,*) '"x","Cul","Cuss","Cusf","Kd","SOH"'
      do j=1,nx					  
         write(61,102) x(j)*Lx,',',cul(j),',',cuss(j),',',cusf(j),',',Kd(j),',',SOH(j)
      enddo
	  close(61)

	  open(11,file=title_I,form='formatted')
      write(11,*) '"x","UO2+2","Ca+2","Cl-","L-","H+","HCO_3-"'
      do j=1,nx
         write(11,105) x(j)*Lx,',',ct(1,j),',',ct(2,j),',',ct(3,j),',',ct(4,j),',',ct(5,j),',',ct(6,j)
      enddo
	  close(11)

	  open(12,file=title_II,form='formatted')
      write(12,*) '"x","LH","LUO_2+","Ca_2UO_2(CO_3)_3","CaUO_2(CO_3)_3-2","UO_2CO_3","(UO_2)_2CO_3(OH)_3^-"'
      do j=1,nx
         write(12,105) x(j)*Lx,',',ct(7,j),',',ct(8,j),',',ct(9,j),',',ct(10,j),',',ct(11,j),',',ct(12,j)
      enddo
	  close(12)

	  open(11,file=title_cha,form='formatted')
      write(11,*) '"x","I","qt","ql","qs","pH","pCO2"'
      do j=1,nx
         write(11,105) x(j)*Lx,',',I(j),',',qt(j),',',ql(j),',',qs(j),',',pH,',',pCO2
      enddo
	  close(11)

	  open(12,file=title_Act_I,form='formatted')
      write(12,*) '"x","UO2+2","Ca+2","Cl-","L-","H+","HCO_3-"'
      do j=1,nx
         write(12,105) x(j)*Lx,',',al(1,j),',',al(2,j),',',al(3,j),',',al(4,j),',',al(5,j),',',al(6,j)
      enddo
	  close(12)

	  open(12,file=title_Act_II,form='formatted')
      write(12,*) '"x","LH","LUO_2+","Ca_2UO_2(CO_3)_3","CaUO_2(CO_3)_3-2","UO_2CO_3","(UO_2)_2CO_3(OH)_3^-"'
      do j=1,nx
         write(12,105) x(j)*Lx,',',al(7,j),',',al(8,j),',',al(9,j),',',al(10,j),',',al(11,j),',',al(12,j)
      enddo
	  close(12)

	  open(13,file=title_C0,form='formatted')
      write(13,*) '"x","UO2+2","Ca+2","Cl-","L-","DHSt"'
      do j=1,nx
         write(13,102) x(j)*Lx,',',totCt(1,j),',',totCt(2,j),',',totCt(3,j),',',totCt(4,j),',',totCt(5,j)
      enddo
	  close(13)

	  102 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)
	  105 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

	  end

      !!!*********************************************************************************
      !!!*						               			 * 
      !!!*                             Speciation	        		     	 *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine SpeciationII(pH,al,ct,cl,cs,dcsdul,SOH,SOH0,gammaA,dgammadI, &
	             I,ct0,nu,z,theta,rho,Ka,Bd0,Ad,gammaAREF,h,pH0,bHS,n,m,indU)
      implicit none
      integer, intent(in) :: n,m,indU
      real*8, intent(in) :: pH,ct0(1:m),theta,rho,nu(1:m,1:n),dcsdul(1:n)
	  real*8, intent(in) :: z(1:n),Ka(1:n),Bd0(1:n),SOH0,Ad,gammaAREF,h,pH0,bHS
      real*8, intent(out) :: ct(1:n),cl(1:n),cs(1:n),SOH
	  real*8, intent(inout) :: al(1:n),gammaA(1:n),dgammadI(1:n),I
      real*8 res,eps0,g(1:m-1)
      integer nit1,k,m1

	  m1 = m-1

	  al(5) = 10.d0**(-pH)
      nit1=240
      k=1
      res=1.d0
      eps0 = 1.d-20
      g = 1.d0
      do while((res.ge.eps0).and.(k.le.nit1))
         !!! calculate concentration
         call calcConc(cl,cs,ct,dcsdul,SOH,theta,rho,gammaA,al,SOH0,pH0,bHS,n,indU)
         !!! calculate Ionic Strength
         call calculateIonic(I,gammaA,dgammadI,n,z,cl,Ad,Bd0,h,gammaAREF)
         !!! Solve System
         call solveSystemII(al,g,ct,cl,cs,ct0,nu,z,theta,rho,m,n,m1)
         !!! calculate Activities
         call ActivityBlock(al,Ka,nu,n,m)
         !!! calculate residual
         res = maxval(dabs(g))
         k=k+1
      enddo

      !print*, ' k =',k,' m1=',m1
      !print*, ' resigual U=', g(1),' [M], residual Ca=', g(2),' [M]', ' resigual Cl=', g(3),' [M]', ' resigual L=', g(4)
	  
      !!! calculate concentration
      call calcConc(cl,cs,ct,dcsdul,SOH,theta,rho,gammaA,al,SOH0,pH0,bHS,n,indU)
      !!! calculate Ionic Strength
      call calculateIonic(I,gammaA,dgammadI,n,z,cl,Ad,Bd0,h,gammaAREF)
	  
      end

	  !!!*********************************************************************************
      !!!*						                		 * 
      !!!*                                solveSystem           	      	         *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine solveSystemII(al,g,ct,cl,cs,ct0,nu,z,theta,rho,m,n,m1)
      implicit none
      integer, intent(in) :: m,n,m1
      real*8, intent(in) :: nu(1:m,1:n),ct(1:n),cl(1:n),cs(1:n),ct0(1:m),z(1:n),theta,rho
      real*8, intent(inout) :: al(1:n)
      real*8, intent(out) :: g(1:m1)
      real*8 dg(1:m1,1:m1),dgT(1:m1,1:m1),A(1:m1,1:m1),eps
      integer i,j,k

      !!!!nu(4,4) = 0
      eps = 1e-4

      g(1) = 0.d0
      do j=1,n
         g(1) = g(1) +nu(1,j)*ct(j)
      enddo
      g(1) = g(1) -ct0(1)

      g(2) = 0.d0
      do j=1,n
         g(2) = g(2) +z(j)*theta*cl(j)
      enddo

	  do k=3,m1
         g(k) = 0.d0
         do j=1,n
            g(k) = g(k) +nu(k,j)*ct(j)
         enddo
         g(k) = g(k) -ct0(k)
	  enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !dg=0.d0
      do k=1,m1
         dg(k,1) = 0.d0
         do j=1,n
            dg(k,1) = dg(k,1) +nu(k,j)*nu(1,j)*theta*cl(j)/al(k)
         enddo
      enddo
      do k=1,m1
         dg(k,2) = 0.d0
         do j=1,n
            dg(k,2) = dg(k,2) +nu(k,j)*z(j)*theta*cl(j)/al(k)
         enddo
      enddo
	  !print*, 'dg(2,2)=',dg(2,2)
	  do i=3,m1
         do k=1,m1
         dg(k,i) = 0.d0
            do j=1,n
               dg(k,i) = dg(k,i) +nu(k,j)*nu(i,j)*theta*cl(j)/al(k)
            enddo
         enddo
	  enddo

      do j=1,n
         dg(1,1) = dg(1,1) +nu(1,j)*rho*cs(j)/al(1)
      enddo

	  do i=3,m1
	     do j=1,n
            dg(1,i) = dg(1,i) +nu(i,j)*rho*cs(j)/al(1)
         enddo
	  enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k=1,m1
         do i=1,m1
            dgT(i,k) = dg(k,i)
         enddo
      enddo

      call invert4t4M(A,dgT)
	  !call inverse(dgT,A,m1)
	  		 		 	  	  	
	 ! 	  pause
					  
      do i=1,m1
         do j=1,m1
            al(i) = al(i) -A(i,j)*g(j)
         enddo
		 !print*, ' dg(i,i)=', dg(i,i)
		 !pause
	     !al(i) = al(i) -g(i)/dg(i,i)
      enddo

      !al = dabs(al)
        
      end

      !!!*********************************************************************************
      !!!*						                	                                     * 
      !!!*                         Calculate Kd and Cusf    	        		         *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine calcKdCusf(Kd,Cusf,Cul,pH,al,gammaA,SOH0,dcsdul,Ka,kusf,kusb,n,nx,indU)
 	  implicit none
	  integer, intent(in) :: n,nx,indU
      real*8, intent(in) :: Cul(1:nx),pH,al(1:n,1:nx),SOH0
	  real*8, intent(in) :: gammaA(1:n,1:nx),dcsdul(1:n),Ka(1:n),kusf,kusb
	  real*8, intent(out) :: Kd(1:nx),Cusf(1:nx)
	  real*8 num(1:nx),den(1:nx)
	  integer j

          !SOH = SOH0/(1.0 +dcsdul(1)*al(1)/al(5)**2.d0 +dcsdul(indU)*al(indU)*al(5)**2.d0);

          !cs(1) = SOH*dcsdul(1)*al(1)/al(5)**2.d0;

          !cs(indU) = SOH*dcsdul(indU)*al(indU)*al(5)**2.d0;

	  do j=1,nx
	     num(j) = SOH0*(dcsdul(1)*10.d0**(2.d0*pH) +dcsdul(indU)*(al(indU,j)/al(1,j))*10.d0**(-2.d0*pH)) &
                      /(1.d0 +dcsdul(1)*al(1,j)*10.d0**(2.d0*pH) +dcsdul(indU)*al(indU,j)*10.d0**(-2.d0*pH))
		 den(j) = (1.d0/gammaA(1,j) +(Ka(8)/gammaA(8,j))*al(4,j)) &
				  +(Ka(12)/gammaA(12,j))*al(1,j)*10.d0**(3.d0*pH) & 
				  +(Ka(11)/gammaA(11,j))*10.d0**(2.d0*pH) &
				  +((al(2,j)**2.d0)*Ka(9)/gammaA(9,j) +al(2,j)*Ka(10)/gammaA(10,j)) &
				  *10.d0**(6.d0*pH)
		 ! al(12) = al(1)^2*al(5)^(-3)									
		 Kd(j) = num(j)/den(j)

		 Cusf(j) = (kusb/(kusf +kusb))*Kd(j)*Cul(j)
	  enddo

      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                         Coeff U		        		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************
      subroutine coeffU(betaU,aue,auw,Kd,dKddt,Pe,theta,rho,kusf,kusb,dt,dx,dx2)
 	  implicit none
	  real*8, intent(in) :: Kd,dKddt,Pe,theta,rho,kusf,kusb,dt,dx,dx2
	  real*8, intent(out) :: betaU,aue,auw
	  real*8 den

	  den = (theta +rho*(kusb/(kusf+kusb))*Kd)

      betaU = (dt/2.d0)*rho*(dKddt+kusf*Kd)*(kusb/(kusf+kusb))/den

	  aue = (dt/2.d0)*(1.d0/dx2 +Pe/(2.d0*dx))/den

	  auw = (dt/2.d0)*(1.d0/dx2 -Pe/(2.d0*dx))/den

      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                               Integrate              		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine IntegrateL(phi,phi_old,phi_new,dt,dx,nx)
 	  implicit none
	  integer, intent(in) :: nx
      real*8, intent(in) :: phi_old(1:nx),phi_new(1:nx),dt,dx
      real*8, intent(out) ::phi
	  real*8 dphi(1:nx)
	  integer j

	  dphi = 0.d0

	  phi=0.d0
	  
	  do j=1,nx
	     dphi(j) = (phi_new(j) -phi_old(j))/dt
         phi = phi +dphi(j)*dx
      enddo

      end

	  !!!*********************************************************************************
      !!!*						                										 * 
      !!!*                               Perc Charge             		     	         *
      !!!*                                                                               *
      !!!*********************************************************************************

	  subroutine percCharge(Ap,A0,An,ct,nu,z,n)
 	  implicit none
	  integer, intent(in) :: n
      real*8, intent(in) :: ct(1:n),nu(1:n),z(1:n)
      real*8, intent(out) ::Ap,A0,An
	  real*8 zp(1:n),z0(1:n),zn(1:n)
	  integer j

      do j=1,n
         if(z(j)>0.d0) then
            zp(j) = 1.d0
            z0(j) = 0.d0
            zn(j) = 0.d0
         elseif(z(j)<0.d0) then
            zp(j) = 0.d0
            z0(j) = 0.d0
            zn(j) = 1.d0
         else
            zp(j) = 0.d0
            z0(j) = 1.d0
            zn(j) = 0.d0
         endif 
     enddo

     Ap = 0.d0
     A0 = 0.d0
     An = 0.d0
     do j=1,n
        Ap = Ap +zp(j)*nu(j)*ct(j)
        A0 = A0 +z0(j)*nu(j)*ct(j)
        An = An +zn(j)*nu(j)*ct(j)
     enddo

     end

      !!!*********************************************************************************
      !!!*						               			 * 
      !!!*                                invert4t4M           	        	         *
      !!!*                                                                               *
      !!!*********************************************************************************

      subroutine invert4t4M(I,A)
      implicit none
      real*8, intent(in) :: A(1:4,1:4)
      real*8, intent(out) :: I(1:4,1:4)
      real*8 a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44
      real*8 D

      I = 0.d0
      !!!! define

      !!! Fjst row
      a11 = A(1,1)
      a12 = A(1,2)
      a13 = A(1,3)
      a14 = A(1,4)
      !!! Second row
      a21 = A(2,1)
      a22 = A(2,2)
      a23 = A(2,3)
      a24 = A(2,4)
      !!! Thjd row
      a31 = A(3,1)
      a32 = A(3,2)
      a33 = A(3,3)
      a34 = A(3,4)
      !!! Fourth row
      a41 = A(4,1)
      a42 = A(4,2)
      a43 = A(4,3)
      a44 = A(4,4)

      !!!!! Calculate determinant
      D = a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 &
        + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 &
        + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 &
        + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 &
        + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 &
        + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41

      !!!!! Fjst row
      !!!!! I(1,1)
      I(1,1) = a22*a33*a44 - a22*a34*a43 - a23*a32*a44 + a23*a34*a42 + a24*a32*a43 - a24*a33*a42
      !!!!! I(1,2)
      I(1,2) = a12*a34*a43 - a12*a33*a44 + a13*a32*a44 - a13*a34*a42 - a14*a32*a43 + a14*a33*a42
      !!!!! I(1,3)
      I(1,3) = a12*a23*a44 - a12*a24*a43 - a13*a22*a44 + a13*a24*a42 + a14*a22*a43 - a14*a23*a42
      !!!!! I(1,4)
      I(1,4) = a12*a24*a33 - a12*a23*a34 + a13*a22*a34 - a13*a24*a32 - a14*a22*a33 + a14*a23*a32
      !!!!! Second row
      !!!!! I(2,1)
      I(2,1) = a21*a34*a43 - a21*a33*a44 + a23*a31*a44 - a23*a34*a41 - a24*a31*a43 + a24*a33*a41
      !!!!! I(2,2)
      I(2,2) = a11*a33*a44 - a11*a34*a43 - a13*a31*a44 + a13*a34*a41 + a14*a31*a43 - a14*a33*a41
      !!!!! I(2,3)
      I(2,3) = a11*a24*a43 - a11*a23*a44 + a13*a21*a44 - a13*a24*a41 - a14*a21*a43 + a14*a23*a41
      !!!!! I(2,4)
      I(2,4) = a11*a23*a34 - a11*a24*a33 - a13*a21*a34 + a13*a24*a31 + a14*a21*a33 - a14*a23*a31
      !!!!! Thjd row
      !!!!! I(3,1)
      I(3,1) = a21*a32*a44 - a21*a34*a42 - a22*a31*a44 + a22*a34*a41 + a24*a31*a42 - a24*a32*a41
      !!!!! I(3,2)
      I(3,2) = a11*a34*a42 - a11*a32*a44 + a12*a31*a44 - a12*a34*a41 - a14*a31*a42 + a14*a32*a41
      !!!!! I(3,3)
      I(3,3) = a11*a22*a44 - a11*a24*a42 - a12*a21*a44 + a12*a24*a41 + a14*a21*a42 - a14*a22*a41
      !!!!! I(3,4)
      I(3,4) = a11*a24*a32 - a11*a22*a34 + a12*a21*a34 - a12*a24*a31 - a14*a21*a32 + a14*a22*a31
      !!!!! Fourth row
      !!!!! I(4,1)
      I(4,1) = a21*a33*a42 - a21*a32*a43 + a22*a31*a43 - a22*a33*a41 - a23*a31*a42 + a23*a32*a41
      !!!!! I(4,2)
      I(4,2) = a11*a32*a43 - a11*a33*a42 - a12*a31*a43 + a12*a33*a41 + a13*a31*a42 - a13*a32*a41
      !!!!! I(4,3)
      I(4,3) = a11*a23*a42 - a11*a22*a43 + a12*a21*a43 - a12*a23*a41 - a13*a21*a42 + a13*a22*a41
      !!!!! I(4,4)
      I(4,4) = a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31

      !!!! scale the matrix
      I = I/D

      end

!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
subroutine inverse(a,c,n)
implicit none 
integer, intent(in) :: n
real*8, intent(inout) :: a(n,n)
real*8, intent(out) :: c(n,n)
real*8 L(n,n), U(n,n), b(n), d(n), x(n)
real*8 coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.d0
U=0.d0
b=0.d0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.d0
end do

end !subroutine inverse
