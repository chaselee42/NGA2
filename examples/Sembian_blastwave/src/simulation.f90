!> Sembian blastwave case – high-pressure region initialization (dimensional)
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use amrgrid_class,     only: amrgrid
   use amrmpcomp_class,   only: amrmpcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use amrio_class,       only: amrio
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> Stiffened gas EOS parameters (liquid and gas)
   real(WP) :: GammaL,PinfL,CvL
   real(WP) :: GammaG,PinfG,CvG

   !> Flow parameters
   real(WP) :: Grho0,GP0          !< Pre-shock gas state
   real(WP) :: Grho1,GP1,u1       !< Post-shock gas state (from Rankine-Hugoniot)
   real(WP) :: M1                 !< Post-shock Mach number
   real(WP) :: relshockvel        !< Velocity at which shock moves
   real(WP) :: Lrho0              !< Liquid density
   real(WP) :: Xs                 !< Shock location [m]
   real(WP) :: Ms                 !< Shock Mach number
   real(WP) :: muG,muL            !< Dynamic viscosities
   real(WP) :: sigma

   !> High-pressure region parameters
   real(WP) :: HP_thickness       !< Thickness of HP region [m]
   real(WP) :: HP_center          !< Center of HP region [m]
   real(WP) :: HP_density         !< Density in HP region [kg/m^3]
   real(WP) :: HP_pressure        !< Pressure in HP region [Pa]

   !> Cylinder parameters
   real(WP) :: dcyl               !< Cylinder diameter [m]
   real(WP) :: xcyl               !< Cylinder center x location [m]

   !> Domain dimensions
   real(WP) :: Lx,Ly              !< Domain lengths [m]

   !> Tagging parameters
   real(WP) :: vorticity_tag=huge(1.0_WP)
   real(WP) :: rho_ratio_tag=huge(1.0_WP)

contains

   !> Levelset function for 2D cylinder centered at (xcyl, 0)
   function levelset_cyl(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP*dcyl-sqrt((xyz(1)-xcyl)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP*dcyl-sqrt((xyz(1)-xcyl)**2+xyz(2)**2) ! Enable quasi-2D runs
   end function levelset_cyl

   !> Liquid EOS: P=f(RHO,I) - Stiffened gas
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> Liquid EOS: T=f(RHO,P)
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> Liquid EOS: C=f(RHO,P)
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(max(0.0_WP,GammaL*(P+PinfL)/RHO))
   end function get_CL
   !> Liquid EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IL=(P+GammaL*PinfL)/(RHO*(GammaL-1.0_WP))
   end function get_IL

   !> Gas EOS: P=f(RHO,I) - Ideal gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> Gas EOS: T=f(RHO,P)
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> Gas EOS: C=f(RHO,P)
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(max(0.0_WP,GammaG*(P+PinfG)/RHO))
   end function get_CG
   !> Gas EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IG=(P+GammaG*PinfG)/(RHO*(GammaG-1.0_WP))
   end function get_IG

   !> Generalized mechanical relaxation for stiffened gas EOS pair
   !> Solves quadratic for equilibrium pressure Peq where PL+Pjump=PG=Peq,
   !> then adjusts VF and internal energies via p*dV work exchange.
   !> Conserves phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7)
   !> Enforces pressure jump provided in Pjump
   subroutine P_relax_generalized(VF,Q,Pjump)
      use amrmpcomp_class, only: VFlo,VFhi
      implicit none
      real(WP), intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP), intent(in) :: Pjump
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq
      real(WP), parameter :: RHOGmin=1.0e-2_WP
      real(WP), parameter :: phist=1.0_WP,phi0=0.0_WP   !< Temporal weighting, phist=1 should yield best results
      ! Skip if any conserved quantity is non-positive (EOS undefined)
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2))
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)**2
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG)**2
      cJ=ZL/(ZG+ZL)
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      n1=VF*phist
      n0=VF*(phi0*Pint-phist*cJ*pjump)+Q(3)
      d1=phist+1.0_WP/(GammaL-1.0_WP)
      d0=phi0*Pint-phist*cJ*pjump+GammaL/(GammaL-1.0_WP)*PinfL
      a=d1*(1.0_WP/(GammaG-1.0_WP)+phist*VF)+n1*(-1.0_WP/(GammaG-1.0_WP)-phist)
      b=d1*((GammaG*PinfG-pjump)/(GammaG-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump))+n1*(-(GammaG*PinfG-pjump)/(GammaG-1.0_WP)-phi0*Pint+phist*cJ*pjump)+d0*(1.0_WP/(GammaG-1.0_WP)+phist*VF)+n0*(-1.0_WP/(GammaG-1.0_WP)-phist)
      d=d0*((GammaG*PinfG-pjump)/(GammaG-1.0_WP)-Q(4)+VF*(phi0*Pint-phist*cJ*pjump))+n0*(-(GammaG*PinfG-pjump)/(GammaG-1.0_WP)-phi0*Pint+phist*cJ*pjump)
      ! Get equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.-PinfL.or.Peq-Pjump.le.-PinfG) return
      ! Get equilibrium volume fraction
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.VFlo.or.VFeq.gt.VFhi) return
      ! Adjust conserved quantities
      Q(3)=Q(3)-(phi0*Pint+phist*Peq)*(VFeq-VF)
      Q(4)=Q(4)+(phi0*Pint+phist*Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_generalized

   !> Implicit mechanical relaxation for stiffened gas EOS pair
   subroutine P_relax_implicit(VF,Q)
      use amrmpcomp_class, only: VFlo,VFhi
      implicit none
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP) :: a,b,d,Peq,VFeq
      real(WP) :: invG1G,invG1L,d0,d1,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-2_WP
      ! Skip if any conserved quantity is non-positive (EOS undefined)
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Precompute EOS constants
      invG1G=1.0_WP/(GammaG-1.0_WP)
      invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=GammaL*PinfL*invG1L
      d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G
      facL=invG1G+VF
      ! Quadratic coefficients: a*Peq^2 + b*Peq + d = 0
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      d=d0*(facG-Q(4))-Q(3)*facG
      ! Solve for equilibrium pressure (positive root)
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Bail if pressure is unphysical
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Equilibrium volume fraction from liquid energy constraint
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Update internal energies via p*dV work exchange
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_implicit

   !> Compute viscosity: constant gas and liquid, VF-weighted blend
   !> Contains commented-out Sutherland law for variable gas viscosity (dimensional form)
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pVisc,pBeta,pDiff
      real(WP) :: mu_g,mu_l
      real(WP), parameter :: myeps=1.0e-15_WP
      !> Sutherland's law parameters (dimensional, SI units)
      !> mu_ref = 1.716e-5 Pa·s at T_ref = 273.15 K, S = 110.4 K
      ! real(WP), parameter :: mu_ref=1.716e-5_WP  !< Reference viscosity [Pa·s]
      ! real(WP), parameter :: T_ref=273.15_WP      !< Reference temperature [K]
      ! real(WP), parameter :: S_suth=110.4_WP      !< Sutherland constant [K]
      ! real(WP) :: T_gas                           !< Local gas temperature [K]
      ! Loop over levels
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                     ! --- Gas viscosity ---
                     ! Constant gas viscosity (current default)
                     mu_g=muG
                     ! Sutherland's law for variable gas viscosity (dimensional):
                     ! T_gas=pTG(i,j,k,1)  ! Gas temperature from EOS [K]
                     ! mu_g=mu_ref*(T_gas/T_ref)**1.5_WP*(T_ref+S_suth)/(T_gas+S_suth)
                     ! --- Liquid viscosity ---
                     mu_l=muL
                     ! Mixture viscosity (harmonic averaging)
                     pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps))
                     ! Zero bulk viscosity
                     pBeta(i,j,k,1)=0.0_WP
                     ! Zero thermal diffusivity
                     pDiff(i,j,k,1)=0.0_WP
                  end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback – set Q and VF/barycenters for blastwave + cylinder
   subroutine blastwave_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG_local,pG_local,IG_local
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energy of liquid at ambient
      IEL=get_IL(Lrho0,GP0)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Compute VF and barycenters from levelset (cylinder)
                  call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
                  &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
                  &                              levelset=levelset_cyl,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
                  ! Store volume fraction
                  pVF(i,j,k,1)=myVF
                  ! Store barycenters
                  if (lvl.eq.solver%amr%maxlvl) then
                     pCL(i,j,k,:)=BL
                     pCG(i,j,k,:)=BG
                  end if
                  ! Cell-center x coordinate
                  x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
                  ! Determine gas state: high pressure region or ambient
                  if ((x_cc.gt.(HP_center-0.5_WP*HP_thickness)).and.(x_cc.lt.(HP_center+0.5_WP*HP_thickness))) then
                     ! High pressure region
                     rhoG_local=HP_density
                     pG_local  =HP_pressure
                  else
                     ! Ambient region
                     rhoG_local=Grho0
                     pG_local  =GP0
                  end if
                  IG_local=get_IG(rhoG_local,pG_local)
                  ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG, 0, 0, 0)
                  ! Zero velocity everywhere (blastwave forms from pressure release)
                  pQ(i,j,k,1)=(       myVF)*Lrho0
                  pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG_local
                  pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
                  pQ(i,j,k,4)=pQ(i,j,k,2)*IG_local
                  pQ(i,j,k,5)=0.0_WP
                  pQ(i,j,k,6)=0.0_WP
                  pQ(i,j,k,7)=0.0_WP
               end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine blastwave_init

   !> Tagger based on vorticity and density ratio (from amrcomp_drop)
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP) :: irho_cc,irho_xp,irho_xm,irho_yp,irho_ym,irho_zp,irho_zm
      real(WP) :: vort_x,vort_y,vort_z,vort_mag
      real(WP) :: rho_max,rho_min,rho_nb,rho_ratio
      integer :: i,j,k,ii,jj,kk
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      ! Recast tags
      tags=tags_ptr
      ! Compute tags
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Get local inverse densities
                  irho_cc=1.0_WP/max(sum(pQ(i  ,j,  k,  1:2)),solver%rho_floor)
                  irho_xp=1.0_WP/max(sum(pQ(i+1,j,  k,  1:2)),solver%rho_floor)
                  irho_xm=1.0_WP/max(sum(pQ(i-1,j,  k,  1:2)),solver%rho_floor)
                  irho_yp=1.0_WP/max(sum(pQ(i,  j+1,k,  1:2)),solver%rho_floor)
                  irho_ym=1.0_WP/max(sum(pQ(i,  j-1,k,  1:2)),solver%rho_floor)
                  irho_zp=1.0_WP/max(sum(pQ(i,  j,  k+1,1:2)),solver%rho_floor)
                  irho_zm=1.0_WP/max(sum(pQ(i,  j,  k-1,1:2)),solver%rho_floor)
                  ! Compute vorticity and tag based on it
                  vort_x=(pQ(i,j+1,k,7)*irho_yp-pQ(i,j-1,k,7)*irho_ym)*0.5_WP*dyi-(pQ(i,j,k+1,6)*irho_zp-pQ(i,j,k-1,6)*irho_zm)*0.5_WP*dzi
                  vort_y=(pQ(i,j,k+1,5)*irho_zp-pQ(i,j,k-1,5)*irho_zm)*0.5_WP*dzi-(pQ(i+1,j,k,7)*irho_xp-pQ(i-1,j,k,7)*irho_xm)*0.5_WP*dxi
                  vort_z=(pQ(i+1,j,k,6)*irho_xp-pQ(i-1,j,k,6)*irho_xm)*0.5_WP*dxi-(pQ(i,j+1,k,5)*irho_yp-pQ(i,j-1,k,5)*irho_ym)*0.5_WP*dyi
                  vort_mag=sqrt(vort_x**2+vort_y**2+vort_z**2)
                  if (vort_mag.gt.vorticity_tag) tagarr(i,j,k,1)=SETtag
                  ! Compute density ratio in 3x3x3 stencil and tag based on it
                  rho_max=solver%rho_floor; rho_min=huge(1.0_WP)
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                           rho_nb=sum(pQ(i+ii,j+jj,k+kk,1:2))
                           rho_max=max(rho_max,rho_nb)
                           rho_min=min(rho_min,max(rho_nb,solver%rho_floor))
                        end do; end do; end do
                  rho_ratio=rho_max/rho_min
                  if (rho_ratio.gt.rho_ratio_tag) tagarr(i,j,k,1)=SETtag
               end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         ! Set PinfG to zero (ideal gas)
         PinfG=0.0_WP
         ! Get fluid parameters
         call param_read('Liquid gamma',GammaL)
         call param_read('Gas gamma',GammaG)
         call param_read('Liquid Pinf',PinfL)
         call param_read('Liquid density',Lrho0)
         call param_read('Liquid Q',QL,default=0.0_WP)
         call param_read('Liquid B',BL,default=0.0_WP)
         call param_read('Pre-shock density',Grho0)
         call param_read('Pre-shock pressure',GP0)
         call param_read('Gas Q',QG,default=0.0_WP)
         call param_read('Gas B',BG,default=0.0_WP)
         call param_read('Mach number of shock',Ms)
         call param_read('Shock location',Xs)
         call param_read('Liquid viscosity',muL)
         call param_read('Gas viscosity',muG)
         call param_read('Surface tension coefficient',sigma)
         ! High-pressure region parameters
         call param_read('HP region thickness',HP_thickness)
         call param_read('HP center location',HP_center)
         call param_read('HP density',HP_density)
         call param_read('HP pressure',HP_pressure)
         ! Cylinder geometry
         call param_read('Cylinder diameter',dcyl)
         call param_read('Cylinder location',xcyl)
         ! Domain dimensions
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         ! Use shock relations to get post-shock numbers (for informational / tagging purposes)
         GP1 = GP0 * (2.0_WP*GammaG*Ms**2 - (GammaG-1.0_WP)) / (GammaG+1.0_WP)
         Grho1 = Grho0 * (Ms**2 * (GammaG+1.0_WP) / ((GammaG-1.0_WP)*Ms**2 + 2.0_WP))
         ! Calculate post-shock Mach number
         M1 = sqrt(((GammaG-1.0_WP)*(Ms**2)+2.0_WP)/(2.0_WP*GammaG*(Ms**2)-(GammaG-1.0_WP)))
         ! Calculate post-shock velocity
         u1 = -M1 * sqrt(GammaG*GP1/Grho1) + Ms*sqrt(GammaG*GP0/Grho0)
         ! Velocity at which shock moves
         relshockvel = -Grho1*u1/(Grho0-Grho1)
         ! Set specific heats
         CvL=1077.7_WP
         CvG=1956.45_WP
         ! Log setup
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  Grho0=",es12.5," GP0=",es12.5)') Grho0,GP0; call log(message)
         write(message,'("[Post-shock] Grho1=",es12.5," GP1=",es12.5," u1=",es12.5)') Grho1,GP1,u1; call log(message)
         write(message,'("[Shock velocity] =",es12.5)') relshockvel; call log(message)
         write(message,'("[HP region] center=",es12.5," thickness=",es12.5)') HP_center,HP_thickness; call log(message)
         write(message,'("[HP region] density=",es12.5," pressure=",es12.5)') HP_density,HP_pressure; call log(message)
         write(message,'("[Liquid] rhoL=",es12.5," GammaL=",es12.5," PinfL=",es12.5)') Lrho0,GammaL,PinfL; call log(message)
         write(message,'("[Gas]    GammaG=",es12.5," PinfG=",es12.5)') GammaG,PinfG; call log(message)
         write(message,'("[Visc]   muG=",es12.5," muL=",es12.5)') muG,muL; call log(message)
         write(message,'("[Cylinder] diameter=",es12.5," location=",es12.5)') dcyl,xcyl; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='Sembian_blastwave'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=0.0_WP;      amr%xhi=Lx
         amr%ylo=-0.5_WP*Ly;  amr%yhi=+0.5_WP*Ly
         amr%zlo=-0.5_WP*Ly;  amr%zhi=+0.5_WP*Ly
         amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zper=.true.
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves here
      handle_restart: block
         integer :: restart_step
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=1)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
         time%dt=1.0e-9_WP
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap, amrex_bc_reflect_even, amrex_bc_reflect_odd
         use amrmpcomp_class, only: BC_REFLECT
         use amrdata_class,    only: interp_face_lin
         ! Create flow solver
         call fs%initialize(amr=amr,name='Sembian_blastwave')
         ! Set surface tension coefficient
         !fs%sigma=1.0_WP/Weber
         fs%sigma=sigma
         ! Use face-linear interp if 2D (divfree requires ratio=2 in all dirs)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Provide thermodynamic model (6 EOS pointers)
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getTG=>get_TG
         ! Provide pressure relaxation model
         fs%relax=>P_relax_generalized
         ! Set initial conditions via blastwave callback
         fs%user_init=>blastwave_init
         ! Set BCs
         if (.not.amr%xper) then
            ! fs%Q%lo_bc(1,:)=amrex_bc_foextrap
            ! fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            ! Slip wall
            fs%lo_bc(1)=BC_REFLECT
            fs%Q%lo_bc(1,:)=amrex_bc_reflect_even
            fs%Q%lo_bc(1,5)=amrex_bc_reflect_odd
            fs%U%lo_bc(1,:)=amrex_bc_reflect_odd
            fs%V%lo_bc(1,:)=amrex_bc_reflect_even
            fs%W%lo_bc(1,:)=amrex_bc_reflect_even
            fs%hi_bc(1)=BC_REFLECT
            fs%Q%hi_bc(1,:)=amrex_bc_reflect_even
            fs%Q%hi_bc(1,5)=amrex_bc_reflect_odd
            fs%U%hi_bc(1,:)=amrex_bc_reflect_odd
            fs%V%hi_bc(1,:)=amrex_bc_reflect_even
            fs%W%hi_bc(1,:)=amrex_bc_reflect_even
            ! No slip wall
            ! fs%lo_bc(1)=BC_REFLECT
            ! fs%Q%lo_bc(1,1:4)=amrex_bc_reflect_even
            ! fs%Q%lo_bc(1,5:7)=amrex_bc_reflect_odd
            ! fs%U%lo_bc(1,:)=amrex_bc_reflect_odd
            ! fs%V%lo_bc(1,:)=amrex_bc_reflect_odd
            ! fs%W%lo_bc(1,:)=amrex_bc_reflect_odd
            ! fs%hi_bc(1)=BC_REFLECT
            ! fs%Q%hi_bc(1,1:4)=amrex_bc_reflect_even
            ! fs%Q%hi_bc(1,5:7)=amrex_bc_reflect_odd
            ! fs%U%hi_bc(1,:)=amrex_bc_reflect_odd
            ! fs%V%hi_bc(1,:)=amrex_bc_reflect_odd
            ! fs%W%hi_bc(1,:)=amrex_bc_reflect_odd
         end if
         if (.not.amr%yper) then
            ! fs%Q%lo_bc(2,:)=amrex_bc_foextrap
            ! fs%Q%hi_bc(2,:)=amrex_bc_foextrap
            ! Slip wall
            fs%lo_bc(2)=BC_REFLECT
            fs%Q%lo_bc(2,:)=amrex_bc_reflect_even
            fs%Q%lo_bc(2,6)=amrex_bc_reflect_odd
            fs%U%lo_bc(2,:)=amrex_bc_reflect_even
            fs%V%lo_bc(2,:)=amrex_bc_reflect_odd
            fs%W%lo_bc(2,:)=amrex_bc_reflect_even
            fs%hi_bc(2)=BC_REFLECT
            fs%Q%hi_bc(2,:)=amrex_bc_reflect_even
            fs%Q%hi_bc(2,6)=amrex_bc_reflect_odd
            fs%U%hi_bc(2,:)=amrex_bc_reflect_even
            fs%V%hi_bc(2,:)=amrex_bc_reflect_odd
            fs%W%hi_bc(2,:)=amrex_bc_reflect_even
            ! No slip wall
            ! fs%lo_bc(2)=BC_REFLECT
            ! fs%Q%lo_bc(2,1:4)=amrex_bc_reflect_even
            ! fs%Q%lo_bc(2,5:7)=amrex_bc_reflect_odd
            ! fs%U%lo_bc(2,:)=amrex_bc_reflect_odd
            ! fs%V%lo_bc(2,:)=amrex_bc_reflect_odd
            ! fs%W%lo_bc(2,:)=amrex_bc_reflect_odd
            ! fs%hi_bc(2)=BC_REFLECT
            ! fs%Q%hi_bc(2,1:4)=amrex_bc_reflect_even
            ! fs%Q%hi_bc(2,5:7)=amrex_bc_reflect_odd
            ! fs%U%hi_bc(2,:)=amrex_bc_reflect_odd
            ! fs%V%hi_bc(2,:)=amrex_bc_reflect_odd
            ! fs%W%hi_bc(2,:)=amrex_bc_reflect_odd
         end if
         if (.not.amr%zper) then
            ! fs%Q%lo_bc(3,:)=amrex_bc_foextrap
            ! fs%Q%hi_bc(3,:)=amrex_bc_foextrap
            ! Slip wall
            fs%lo_bc(3)=BC_REFLECT
            fs%Q%lo_bc(3,:)=amrex_bc_reflect_even
            fs%Q%lo_bc(3,7)=amrex_bc_reflect_odd
            fs%U%lo_bc(3,:)=amrex_bc_reflect_even
            fs%V%lo_bc(3,:)=amrex_bc_reflect_even
            fs%W%lo_bc(3,:)=amrex_bc_reflect_odd
            fs%hi_bc(3)=BC_REFLECT
            fs%Q%hi_bc(3,:)=amrex_bc_reflect_even
            fs%Q%hi_bc(3,7)=amrex_bc_reflect_odd
            fs%U%hi_bc(3,:)=amrex_bc_reflect_even
            fs%V%hi_bc(3,:)=amrex_bc_reflect_even
            fs%W%hi_bc(3,:)=amrex_bc_reflect_odd
            ! No slip wall
            ! fs%lo_bc(3)=BC_REFLECT
            ! fs%Q%lo_bc(3,1:4)=amrex_bc_reflect_even
            ! fs%Q%lo_bc(3,5:7)=amrex_bc_reflect_odd
            ! fs%U%lo_bc(3,:)=amrex_bc_reflect_odd
            ! fs%V%lo_bc(3,:)=amrex_bc_reflect_odd
            ! fs%W%lo_bc(3,:)=amrex_bc_reflect_odd
            ! fs%hi_bc(3)=BC_REFLECT
            ! fs%Q%hi_bc(3,1:4)=amrex_bc_reflect_even
            ! fs%Q%hi_bc(3,5:7)=amrex_bc_reflect_odd
            ! fs%U%hi_bc(3,:)=amrex_bc_reflect_odd
            ! fs%V%hi_bc(3,:)=amrex_bc_reflect_odd
            ! fs%W%hi_bc(3,:)=amrex_bc_reflect_odd
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=7,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging vorticity',vorticity_tag)
         call param_read('Tagging rho ratio',rho_ratio_tag)
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Fresh start
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC and reset moments
            call fs%build_plic(time%t)
            call fs%build_subVF()
            ! Initialize primitive variables
            call fs%get_primitive(Q=fs%Q)
            ! Initialize face velocities
            call fs%get_face_velocity()
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'Sembian_blastwave',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(fs%C,1,'C')
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%dPmax,'dPmax')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%add_column(fs%CFLst ,'CFLst' )
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_sl,'sl_max')
         call tfile%add_column(fs%wtmin_sl,'sl_min')
         call tfile%add_column(fs%wtmax_fv,'fv_max')
         call tfile%add_column(fs%wtmin_fv,'fv_min')
         call tfile%add_column(fs%wtmax_div,'div_max')
         call tfile%add_column(fs%wtmin_div,'div_min')
         call tfile%add_column(fs%wtmax_plicnet,'plicnet_max')
         call tfile%add_column(fs%wtmin_plicnet,'plicnet_min')
         call tfile%add_column(fs%wtmax_polygon,'polygon_max')
         call tfile%add_column(fs%wtmin_polygon,'polygon_min')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%add_column(fs%nmixed_min,'mixed_min')
         call tfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         if (.not.restarted.and.time%n.le.1000) then
            time%dt=1e-9_WP
         else
            call time%adjust_dt()
         end if
         call time%increment()

         ! Remember old state
         call fs%store_old()

         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Get most up-to-date pressure
         call fs%apply_relax(time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities
         call fs%get_face_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         ! Add surface tension term
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         ! Increment Q without pressure gradient
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Get most up-to-date pressure
         call fs%apply_relax(time=time%t)
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities
         call fs%get_face_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=time%dt)
         ! Add surface tension term
         call fs%add_surface_tension(scale=time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         ! ======================================================================================

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         if (fs%amr%amRoot) print*, 'time: ', time%t, ' dt: ', time%dt,' occurs: ',viz_evt%occurs()
         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time%t)
         end if
         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/Sembian_blastwave_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
