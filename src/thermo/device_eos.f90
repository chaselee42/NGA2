!> Flat, non-polymorphic, device-callable equation-of-state (EOS) module for NGA2.
!> Supports Ideal Gas, Stiffened Gas, Noble-Abel Stiffened Gas (NASG), and Mie-Gruneisen EOS models.
module device_eos
   use precision, only: WP
   implicit none

   ! EOS model identifiers
   integer, parameter :: EOS_IDEAL_GAS     = 1
   integer, parameter :: EOS_STIFFENED_GAS = 2
   integer, parameter :: EOS_NASG          = 3
   integer, parameter :: EOS_MIE_GRUNEISEN = 4

   ! Maximum size for tabulated reference curves (Mie-Gruneisen)
   integer, parameter :: MAX_NTAB = 2000

   !> Flat derived type containing all parameters for supported EOS models.
   type :: device_eos_params
      integer  :: model = EOS_IDEAL_GAS

      ! Ideal Gas / Stiffened Gas / NASG parameters
      real(WP) :: gamma = 0.0_WP
      real(WP) :: cv    = 0.0_WP
      real(WP) :: cp    = 0.0_WP
      real(WP) :: R     = 0.0_WP
      real(WP) :: q     = 0.0_WP
      real(WP) :: qp    = 0.0_WP

      ! Stiffened Gas extension
      real(WP) :: pinf  = 0.0_WP

      ! NASG extensions
      real(WP) :: b       = 0.0_WP
      real(WP) :: brhomax = 0.9_WP

      ! Mie-Gruneisen parameters
      real(WP) :: rho0    = 0.0_WP
      real(WP) :: c0      = 0.0_WP
      real(WP) :: s1      = 0.0_WP
      real(WP) :: s2      = 0.0_WP
      real(WP) :: s3      = 0.0_WP
      real(WP) :: gamma0  = 0.0_WP
      real(WP) :: T0      = 0.0_WP
      real(WP) :: etamin  = -1.0_WP
      real(WP) :: etamax  = 0.95_WP
      real(WP) :: Dmin    = 1.0e-3_WP
      integer  :: ntab    = 0
      real(WP) :: gr0m    = 0.0_WP
      real(WP) :: etaknee = 0.0_WP
      real(WP) :: pr_knee = 0.0_WP
      real(WP) :: er_knee = 0.0_WP
      real(WP) :: deta    = 0.0_WP
      real(WP) :: Ttab(MAX_NTAB) = 0.0_WP
   end type device_eos_params

   !> Device-resident copy of the active material's EOS parameters.
   !> Populated on the host with device_eos_init_from_material, then pushed to the
   !> device once with device_eos_push. Kernels reference this directly rather than
   !> mapping the struct per-kernel -- Ttab alone is MAX_NTAB reals, far too large to
   !> firstprivate per team.
   type(device_eos_params) :: dev_eos
   !$omp declare target to(dev_eos)

   ! OpenMP target declarations for procedures
   ! NOTE: device_eos_init_from_material is deliberately NOT listed -- it does a
   ! select type on class(material) and is host-only by construction.
   !$omp declare target to(device_mg_ode_rhs, device_mg_get_ref, device_mg_get_ref_d, device_mg_get_Tref, device_mg_get_rho_from_p_T, &
   !$omp&                  device_get_p_from_rho_e, device_get_T_from_p_rho, device_get_c_from_p_rho, &
   !$omp&                  device_get_e_from_p_rho, device_get_e_from_p_T, device_get_p_from_rho_T, &
   !$omp&                  device_get_rho_from_p_T, device_get_cv_from_rho_T, device_get_h_from_p_T, &
   !$omp&                  device_get_hk_from_p_T, device_get_s_from_p_T, device_get_g_from_p_T, &
   !$omp&                  device_get_gruneisen_from_rho_e, device_get_T_from_rho_e, device_get_c_from_rho_e, &
   !$omp&                  device_get_cv_from_rho_e, device_get_rhoe_from_p_rho, device_get_rhoe_from_p_T)

contains

   ! ============================================================================
   ! HOST INITIALIZATION HELPER
   ! ============================================================================

   !> Push the host-side contents of dev_eos to the device.
   !> Call once after device_eos_init_from_material, before any kernel that uses it.
   subroutine device_eos_push()
      implicit none
      !$omp target update to(dev_eos)
   end subroutine device_eos_push

   !> Populate a flat device_eos_params structure from a polymorphic host material object
   subroutine device_eos_init_from_material(mat, dev_eos)
      use material_class,       only: material
      use ideal_gas_class,      only: ideal_gas
      use stiffened_gas_class,  only: stiffened_gas
      use nasg_class,           only: nasg
      use mie_gruneisen_class,  only: mie_gruneisen
      implicit none
      class(material), intent(in) :: mat
      type(device_eos_params), intent(out) :: dev_eos

      select type (m => mat)
      type is (ideal_gas)
         dev_eos%model = EOS_IDEAL_GAS
         dev_eos%gamma = m%gamma
         dev_eos%cv    = m%cv
         dev_eos%cp    = m%cp
         dev_eos%R     = m%R
         dev_eos%q     = m%q
         dev_eos%qp    = m%qp
      type is (stiffened_gas)
         dev_eos%model = EOS_STIFFENED_GAS
         dev_eos%gamma = m%gamma
         dev_eos%pinf  = m%pinf
         dev_eos%cv    = m%cv
         dev_eos%cp    = m%cp
         dev_eos%R     = m%R
         dev_eos%q     = m%q
         dev_eos%qp    = m%qp
      type is (nasg)
         dev_eos%model   = EOS_NASG
         dev_eos%gamma   = m%gamma
         dev_eos%pinf    = m%pinf
         dev_eos%b       = m%b
         dev_eos%brhomax = m%brhomax
         dev_eos%cv      = m%cv
         dev_eos%cp      = m%cp
         dev_eos%R       = m%R
         dev_eos%q       = m%q
         dev_eos%qp      = m%qp
      type is (mie_gruneisen)
         dev_eos%model   = EOS_MIE_GRUNEISEN
         dev_eos%rho0    = m%rho0
         dev_eos%c0      = m%c0
         dev_eos%s1      = m%s1
         dev_eos%s2      = m%s2
         dev_eos%s3      = m%s3
         dev_eos%gamma0  = m%gamma0
         dev_eos%cv      = m%cv
         dev_eos%T0      = m%T0
         dev_eos%q       = m%q
         dev_eos%qp      = m%qp
         dev_eos%etamin  = m%etamin
         dev_eos%etamax  = m%etamax
         dev_eos%Dmin    = m%Dmin
         dev_eos%ntab    = m%ntab
         dev_eos%gr0m    = m%gr0m
         dev_eos%etaknee = m%etaknee
         dev_eos%pr_knee = m%pr_knee
         dev_eos%er_knee = m%er_knee
         dev_eos%deta    = m%deta
         if (allocated(m%Ttab)) then
            dev_eos%Ttab(1:min(m%ntab, MAX_NTAB)) = m%Ttab(1:min(m%ntab, MAX_NTAB))
         end if
      class default
         ! Unsupported or custom material model for device offload
      end select
   end subroutine device_eos_init_from_material

   ! ============================================================================
   ! MIE-GRUNEISEN INTERNAL DEVICE HELPERS
   ! ============================================================================

   pure real(WP) function device_mg_ode_rhs(eos, eta, T) result(f)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: eta, T
      real(WP) :: D, Dp, pH, dpH, B
      if (eta > eos%etaknee) then
         f = 0.0_WP
      else if (eta > 0.0_WP) then
         D   = max(1.0_WP - eos%s1*eta - eos%s2*eta**2 - eos%s3*eta**3, eos%Dmin)
         Dp  = -(eos%s1 + 2.0_WP*eos%s2*eta + 3.0_WP*eos%s3*eta**2)
         pH  = eos%rho0 * eos%c0**2 * eta / D**2
         dpH = eos%rho0 * eos%c0**2 * (D - 2.0_WP*eta*Dp) / D**3
         B   = (dpH*eta - pH) / (2.0_WP * eos%rho0)
         f   = eos%gamma0 * T + B / eos%cv
      else
         B   = -eos%c0**2 * eta / (1.0_WP - eta)
         f   = eos%gamma0 * T + B / eos%cv
      end if
   end function device_mg_ode_rhs

   pure subroutine device_mg_get_ref(eos, rho, pr, er)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in)  :: rho
      real(WP), intent(out) :: pr, er
      real(WP) :: eta, D
      eta = 1.0_WP - eos%rho0 / rho
      if (eta > eos%etaknee) then
         pr = eos%pr_knee; er = eos%er_knee
      else if (eta > 0.0_WP) then
         D  = max(1.0_WP - eos%s1*eta - eos%s2*eta**2 - eos%s3*eta**3, eos%Dmin)
         pr = eos%rho0 * eos%c0**2 * eta / D**2
         er = pr * eta / (2.0_WP * eos%rho0)
      else
         pr = eos%c0**2 * (rho - eos%rho0)
         er = 0.0_WP
      end if
   end subroutine device_mg_get_ref

   pure subroutine device_mg_get_ref_d(eos, rho, pr, er, dpr, der)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in)  :: rho
      real(WP), intent(out) :: pr, er, dpr, der
      real(WP) :: eta, D, Dp, dpH, detadrho
      eta = 1.0_WP - eos%rho0 / rho
      if (eta > eos%etaknee) then
         pr = eos%pr_knee; er = eos%er_knee; dpr = 0.0_WP; der = 0.0_WP
      else if (eta > 0.0_WP) then
         D   = max(1.0_WP - eos%s1*eta - eos%s2*eta**2 - eos%s3*eta**3, eos%Dmin)
         Dp  = -(eos%s1 + 2.0_WP*eos%s2*eta + 3.0_WP*eos%s3*eta**2)
         pr  = eos%rho0 * eos%c0**2 * eta / D**2
         er  = pr * eta / (2.0_WP * eos%rho0)
         detadrho = eos%rho0 / rho**2
         dpH = eos%rho0 * eos%c0**2 * (D - 2.0_WP*eta*Dp) / D**3
         dpr = dpH * detadrho
         der = (dpH*eta + pr) / (2.0_WP * eos%rho0) * detadrho
      else
         pr  = eos%c0**2 * (rho - eos%rho0)
         er  = 0.0_WP
         dpr = eos%c0**2
         der = 0.0_WP
      end if
   end subroutine device_mg_get_ref_d

   pure real(WP) function device_mg_get_Tref(eos, rho) result(Tr)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho
      real(WP) :: eta, w
      integer :: i
      eta = min(max(1.0_WP - eos%rho0/rho, eos%etamin), eos%etamax)
      w   = (eta - eos%etamin) / eos%deta
      i   = min(int(w) + 1, eos%ntab - 1)
      w   = w - real(i - 1, WP)
      Tr  = (1.0_WP - w) * eos%Ttab(i) + w * eos%Ttab(i+1)
   end function device_mg_get_Tref

   pure real(WP) function device_mg_get_rho_from_p_T(eos, p, T, y) result(rho)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rlo, rhi, flo, fhi, f, df, pr, er, dpr, der, dTr, eta
      integer :: it
      integer, parameter :: itmax = 100
      real(WP), parameter :: rtol = 1.0e-14_WP

      rlo = 1.0e-9_WP * eos%rho0
      rhi = eos%rho0 / (1.0_WP - eos%etaknee)
      flo = device_get_p_from_rho_T(eos, rlo, T, y) - p
      fhi = device_get_p_from_rho_T(eos, rhi, T, y) - p
      if (flo >= 0.0_WP) then; rho = rlo; return; end if
      if (fhi <= 0.0_WP) then; rho = rhi; return; end if
      rho = eos%rho0
      do it = 1, itmax
         call device_mg_get_ref_d(eos, rho, pr, er, dpr, der)
         f = pr + eos%gr0m * eos%cv * (T - device_mg_get_Tref(eos, rho)) - p
         if (f > 0.0_WP) then; rhi = rho; else; rlo = rho; end if
         if (rhi - rlo < rtol * eos%rho0) exit
         eta = 1.0_WP - eos%rho0 / rho
         dTr = 0.0_WP
         if (eta > eos%etamin .and. eta < eos%etaknee) &
            dTr = device_mg_ode_rhs(eos, eta, device_mg_get_Tref(eos, rho)) * eos%rho0 / rho**2
         df = dpr - eos%gr0m * eos%cv * dTr
         if (abs(df) > tiny(1.0_WP) .and. abs(f) < 0.5_WP * abs(df) * (rhi - rlo)) then
            rho = rho - f / df
            if (rho <= rlo .or. rho >= rhi) rho = 0.5_WP * (rlo + rhi)
         else
            rho = 0.5_WP * (rlo + rhi)
         end if
      end do
   end function device_mg_get_rho_from_p_T

   ! ============================================================================
   ! DEVICE-CALLABLE THERMO PROCEDURES (FLAT DISPATCH)
   ! ============================================================================

   pure real(WP) function device_get_p_from_rho_e(eos, rho, e, y) result(p)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         p = (eos%gamma - 1.0_WP) * rho * (e - eos%q)
      case (EOS_STIFFENED_GAS)
         p = (eos%gamma - 1.0_WP) * rho * (e - eos%q) - eos%gamma * eos%pinf
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         p = (eos%gamma - 1.0_WP) * rho * (e - eos%q) / ombm - eos%gamma * eos%pinf
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref(eos, rho, pr, er)
         p = pr + eos%gr0m * (e - eos%q - er)
      case default
         p = 0.0_WP
      end select
   end function device_get_p_from_rho_e

   pure real(WP) function device_get_T_from_p_rho(eos, p, rho, y) result(T)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         T = p / (eos%R * rho)
      case (EOS_STIFFENED_GAS)
         T = (p + eos%pinf) / (eos%R * rho)
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         T = (p + eos%pinf) * ombm / (eos%R * rho)
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref(eos, rho, pr, er)
         T = device_mg_get_Tref(eos, rho) + (p - pr) / (eos%gr0m * eos%cv)
      case default
         T = 0.0_WP
      end select
   end function device_get_T_from_p_rho

   pure real(WP) function device_get_c_from_p_rho(eos, p, rho, y) result(c)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, dpr, der, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         c = sqrt(max(0.0_WP, eos%gamma * p / rho))
      case (EOS_STIFFENED_GAS)
         c = sqrt(max(0.0_WP, eos%gamma * (p + eos%pinf) / rho))
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         c = sqrt(max(0.0_WP, eos%gamma * (p + eos%pinf) / (rho * ombm)))
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref_d(eos, rho, pr, er, dpr, der)
         c = sqrt(max(0.0_WP, dpr - eos%gr0m * der + p * eos%gr0m / rho**2))
      case default
         c = 0.0_WP
      end select
   end function device_get_c_from_p_rho

   pure real(WP) function device_get_e_from_p_rho(eos, p, rho, y) result(e)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         e = p / ((eos%gamma - 1.0_WP) * rho) + eos%q
      case (EOS_STIFFENED_GAS)
         e = (p + eos%gamma * eos%pinf) / ((eos%gamma - 1.0_WP) * rho) + eos%q
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         e = ombm * (p + eos%gamma * eos%pinf) / ((eos%gamma - 1.0_WP) * rho) + eos%q
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref(eos, rho, pr, er)
         e = eos%q + er + (p - pr) / eos%gr0m
      case default
         e = 0.0_WP
      end select
   end function device_get_e_from_p_rho

   pure real(WP) function device_get_e_from_p_T(eos, p, T, y) result(e)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         e = eos%cv * T + eos%q
      case (EOS_STIFFENED_GAS, EOS_NASG)
         e = eos%cv * T * (p + eos%gamma * eos%pinf) / (p + eos%pinf) + eos%q
      case (EOS_MIE_GRUNEISEN)
         rho = device_mg_get_rho_from_p_T(eos, p, T, y)
         e = device_get_e_from_p_rho(eos, p, rho, y)
      case default
         e = 0.0_WP
      end select
   end function device_get_e_from_p_T

   pure real(WP) function device_get_p_from_rho_T(eos, rho, T, y) result(p)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         p = eos%R * rho * T
      case (EOS_STIFFENED_GAS)
         p = eos%R * rho * T - eos%pinf
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         p = eos%R * rho * T / ombm - eos%pinf
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref(eos, rho, pr, er)
         p = pr + eos%gr0m * eos%cv * (T - device_mg_get_Tref(eos, rho))
      case default
         p = 0.0_WP
      end select
   end function device_get_p_from_rho_T

   pure real(WP) function device_get_rho_from_p_T(eos, p, T, y) result(rho)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         rho = p / (eos%R * T)
      case (EOS_STIFFENED_GAS)
         rho = (p + eos%pinf) / (eos%R * T)
      case (EOS_NASG)
         rho = (p + eos%pinf) / (eos%R * T + eos%b * (p + eos%pinf))
         if (eos%b * rho > eos%brhomax) then
            rho = (1.0_WP - eos%brhomax) * (p + eos%pinf) / (eos%R * T)
         end if
      case (EOS_MIE_GRUNEISEN)
         rho = device_mg_get_rho_from_p_T(eos, p, T, y)
      case default
         rho = 0.0_WP
      end select
   end function device_get_rho_from_p_T

   pure real(WP) function device_get_cv_from_rho_T(eos, rho, T, y) result(cv)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, T
      real(WP), dimension(:), intent(in) :: y
      cv = eos%cv
   end function device_get_cv_from_rho_T

   pure real(WP) function device_get_h_from_p_T(eos, p, T, y) result(h)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      select case (eos%model)
      case (EOS_IDEAL_GAS, EOS_STIFFENED_GAS)
         h = eos%cp * T + eos%q
      case (EOS_NASG)
         h = eos%cp * T + eos%b * p + eos%q
      case (EOS_MIE_GRUNEISEN)
         rho = device_mg_get_rho_from_p_T(eos, p, T, y)
         h = device_get_e_from_p_rho(eos, p, rho, y) + p / rho
      case default
         h = 0.0_WP
      end select
   end function device_get_h_from_p_T

   pure subroutine device_get_hk_from_p_T(eos, p, T, y, hk)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(:), intent(out) :: hk
      if (size(hk) > 0) then
         hk(1) = device_get_h_from_p_T(eos, p, T, y)
      end if
   end subroutine device_get_hk_from_p_T

   pure real(WP) function device_get_s_from_p_T(eos, p, T, y) result(s)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         s = eos%cp * log(T) - eos%R * log(p) + eos%qp
      case (EOS_STIFFENED_GAS, EOS_NASG)
         s = eos%cp * log(T) - eos%R * log(p + eos%pinf) + eos%qp
      case (EOS_MIE_GRUNEISEN)
         rho = device_mg_get_rho_from_p_T(eos, p, T, y)
         s = eos%cv * (log(T) + eos%gr0m / rho) + eos%qp
      case default
         s = 0.0_WP
      end select
   end function device_get_s_from_p_T

   pure real(WP) function device_get_g_from_p_T(eos, p, T, y) result(g)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         g = (eos%cp - eos%qp) * T - T * (eos%cp * log(T) - eos%R * log(p)) + eos%q
      case (EOS_STIFFENED_GAS)
         g = (eos%cp - eos%qp) * T - T * (eos%cp * log(T) - eos%R * log(p + eos%pinf)) + eos%q
      case (EOS_NASG)
         g = eos%cp * T + eos%b * p + eos%q - T * (eos%cp * log(T) - eos%R * log(p + eos%pinf) + eos%qp)
      case (EOS_MIE_GRUNEISEN)
         rho = device_mg_get_rho_from_p_T(eos, p, T, y)
         g = device_get_e_from_p_rho(eos, p, rho, y) + p / rho - T * (eos%cv * (log(T) + eos%gr0m / rho) + eos%qp)
      case default
         g = 0.0_WP
      end select
   end function device_get_g_from_p_T

   pure real(WP) function device_get_gruneisen_from_rho_e(eos, rho, e, y) result(gruneisen)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS, EOS_STIFFENED_GAS)
         gruneisen = eos%gamma - 1.0_WP
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         gruneisen = (eos%gamma - 1.0_WP) / ombm
      case (EOS_MIE_GRUNEISEN)
         gruneisen = eos%gr0m / rho
      case default
         gruneisen = 0.0_WP
      end select
   end function device_get_gruneisen_from_rho_e

   pure real(WP) function device_get_T_from_rho_e(eos, rho, e, y) result(T)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         T = (e - eos%q) / eos%cv
      case (EOS_STIFFENED_GAS)
         T = (e - eos%q - eos%pinf / rho) / eos%cv
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         T = (e - eos%q - eos%pinf * ombm / rho) / eos%cv
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref(eos, rho, pr, er)
         T = device_mg_get_Tref(eos, rho) + (e - eos%q - er) / eos%cv
      case default
         T = 0.0_WP
      end select
   end function device_get_T_from_rho_e

   pure real(WP) function device_get_c_from_rho_e(eos, rho, e, y) result(c)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: pr, er, dpr, der, p, ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         c = sqrt(max(0.0_WP, eos%gamma * (eos%gamma - 1.0_WP) * (e - eos%q)))
      case (EOS_STIFFENED_GAS)
         c = sqrt(max(0.0_WP, eos%gamma * (eos%gamma - 1.0_WP) * (e - eos%q - eos%pinf / rho)))
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         p = (eos%gamma - 1.0_WP) * rho * (e - eos%q) / ombm - eos%gamma * eos%pinf
         c = sqrt(max(0.0_WP, eos%gamma * (p + eos%pinf) / (rho * ombm)))
      case (EOS_MIE_GRUNEISEN)
         call device_mg_get_ref_d(eos, rho, pr, er, dpr, der)
         p = pr + eos%gr0m * (e - eos%q - er)
         c = sqrt(max(0.0_WP, dpr - eos%gr0m * der + p * eos%gr0m / rho**2))
      case default
         c = 0.0_WP
      end select
   end function device_get_c_from_rho_e

   pure real(WP) function device_get_cv_from_rho_e(eos, rho, e, y) result(cv)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: rho, e
      real(WP), dimension(:), intent(in) :: y
      cv = eos%cv
   end function device_get_cv_from_rho_e

   pure real(WP) function device_get_rhoe_from_p_rho(eos, p, rho, y) result(rhoe)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: ombm
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         rhoe = p / (eos%gamma - 1.0_WP) + rho * eos%q
      case (EOS_STIFFENED_GAS)
         rhoe = (p + eos%gamma * eos%pinf) / (eos%gamma - 1.0_WP) + rho * eos%q
      case (EOS_NASG)
         ombm = max(1.0_WP - eos%b * rho, 1.0_WP - eos%brhomax)
         rhoe = ombm * (p + eos%gamma * eos%pinf) / (eos%gamma - 1.0_WP) + rho * eos%q
      case (EOS_MIE_GRUNEISEN)
         rhoe = rho * device_get_e_from_p_rho(eos, p, rho, y)
      case default
         rhoe = 0.0_WP
      end select
   end function device_get_rhoe_from_p_rho

   pure real(WP) function device_get_rhoe_from_p_T(eos, p, T, y) result(rhoe)
      implicit none
      type(device_eos_params), intent(in) :: eos
      real(WP), intent(in) :: p, T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: rho
      select case (eos%model)
      case (EOS_IDEAL_GAS)
         rhoe = p * (eos%cv * T + eos%q) / (eos%R * T)
      case (EOS_STIFFENED_GAS)
         rhoe = ((p + eos%gamma * eos%pinf) * eos%cv * T + eos%q * (p + eos%pinf)) / (eos%R * T)
      case (EOS_NASG)
         rho = (p + eos%pinf) / (eos%R * T + eos%b * (p + eos%pinf))
         if (eos%b * rho > eos%brhomax) then
            rho = (1.0_WP - eos%brhomax) * (p + eos%pinf) / (eos%R * T)
            rhoe = (1.0_WP - eos%brhomax) * (p + eos%gamma * eos%pinf) / (eos%gamma - 1.0_WP) + rho * eos%q
         else
            rhoe = ((p + eos%gamma * eos%pinf) * eos%cv * T + eos%q * (p + eos%pinf)) / (eos%R * T + eos%b * (p + eos%pinf))
         end if
      case (EOS_MIE_GRUNEISEN)
         rho = device_mg_get_rho_from_p_T(eos, p, T, y)
         rhoe = rho * device_get_e_from_p_rho(eos, p, rho, y)
      case default
         rhoe = 0.0_WP
      end select
   end function device_get_rhoe_from_p_T

end module device_eos
