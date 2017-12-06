!-----------------------------------------------------------------------
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the FABM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Jorn Bruggeman
!-----------------------------------------------------------------------

#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_halite
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_halite
!     Variable identifiers
      type (type_state_variable_id)        :: id_Na,id_Cl,id_NaCl,id_Het
      type (type_dependency_id)            :: id_temp
      type (type_diagnostic_variable_id)   :: id_NaCl_sat

      real(rk) :: K_NaCl
      real(rk) :: K_NaCl_diss
      real(rk) :: K_NaCl_form
      real(rk) :: k_diss_bot      
      real(rk) :: w_NaCl

   contains
      procedure :: initialize
      procedure :: do
   end type

   real(rk), parameter :: u_Na = 22.9898_rk ! atomic masss of sodium in g/mol
   real(rk), parameter :: u_Cl = 35.453_rk  ! atomic masss of chloride in g/mol
   real(rk), parameter :: NaCl_density = 2.16e6_rk ! 2.16 g/cm3 = 2.16e6 g/m3

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the BROM equilibrium constant model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_halite), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
   !!!!!!!!!!!!!!!!!!!!!!!!!   37.3 mol**2/L**2 =37.3*10e6*10e6 mmol**2/m**2
   call self%get_parameter(self%K_NaCl,     'K_NaCl',     'mmol**2/m**6', 'NaCl conditional equilibrium constant',default=37.19e12_rk) ! default=37.19e12_rk)
   call self%get_parameter(self%K_NaCl_diss,'k_NaCl_diss','mmol/m**3/d',  'Specific rate of dissolution of NaCl',default=1._rk)
   call self%get_parameter(self%k_diss_bot, 'k_NaCl_diss_bot', 'mmol/m**2/d',  'maximum dissolution rate of bottom NaCl',default=0.2e6_rk)
   call self%get_parameter(self%K_NaCl_form,'k_NaCl_form','mmol/m**3/d',  'Specific rate of formation of NaCl',default=1._rk)
   call self%get_parameter(self%w_NaCl,     'w_NaCl',     'm/d',          'Sedimentation rate for particulate NaCl',default=1._rk)
!   call self%get_parameter(self%porosity,   'porosity',   '-',           'Porosity of precipitated NaCl',default=0.5_rk)

   call self%register_state_variable(self%id_Na,   'Na',   'mmol/m**3','Na',  minimum=0.0_rk)
   call self%register_state_variable(self%id_Cl,   'Cl',   'mmol/m**3','Cl',  minimum=0.0_rk)
   call self%register_state_variable(self%id_NaCl, 'NaCl', 'mmol/m**3','NaCl', minimum=0.0_rk, &
        vertical_movement=-self%w_NaCl/86400._rk)
!   call self%register_state_variable(self%id_NaCl_bot, 'NaCl_bot', 'mmol/m**2','particulate sodium chloride at bottom', minimum=0.0_rk)

    call self%register_state_dependency(self%id_Het, 'Het', 'mmol N/m**3','Het')

 !  call self%register_diagnostic_variable(self%id_psu,     'salt' ,   '-','salinity in PSU', standard_variable=standard_variables%practical_salinity)
   call self%register_diagnostic_variable(self%id_NaCl_sat,'NaCl_sat','-','NaCl saturation')
!   call self%register_diagnostic_variable(self%id_h_bot,   'h_bot',   'm','thickness of bottom salt layer')

!   call self%add_to_aggregate_variable(total_sodium, self%id_Na)
!   call self%add_to_aggregate_variable(total_sodium, self%id_NaCl)
!   call self%add_to_aggregate_variable(total_sodium, self%id_NaCl_bot)
!   call self%add_to_aggregate_variable(total_chloride, self%id_Cl)
!   call self%add_to_aggregate_variable(total_chloride, self%id_NaCl)
!   call self%add_to_aggregate_variable(total_chloride, self%id_NaCl_bot)

   ! Let Na and Cl ions contribute to salinity metric, and bottom-bound particulate NaCL to bottom thickness
!   call self%add_to_aggregate_variable(mass_concentration_of_solute, self%id_Na, scale_factor=u_Na*1e-6_rk)
!   call self%add_to_aggregate_variable(mass_concentration_of_solute, self%id_Cl, scale_factor=u_Cl*1e-6_rk)
!   call self%add_to_aggregate_variable(volume_fraction_of_particulates, self%id_NaCl_bot, scale_factor=(u_Na+u_Cl)/1000/NaCl_density/(1-self%porosity))

   call self%register_dependency(self%id_temp,standard_variables%temperature)
!Specify that rates are per day (default: per second)
    self%dt = 86400._rk
   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_halite),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) :: temp,Na,Cl,NaCl, Het
   real(rk) :: Om_NaCl,NaCl_prec,NaCl_diss
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature - not used yet

   _GET_(self%id_Na,Na)
   _GET_(self%id_Cl,Cl)
   _GET_(self%id_NaCl,NaCl)

   _GET_(self%id_Het, Het)
   
   ! NaCl saturation state (NB all units mmol/m**3)
   Om_NaCl=Na*Cl/(self%K_NaCl) & 
       *exp(2.86623-571.361/(temp+272.15) +77.128/(temp+272.15)/(temp+272.15))/2.57 !dependence on T
   
   _SET_DIAGNOSTIC_(self%id_NaCl_sat,Om_NaCl)

   ! Precipitation/dissolution
   NaCl_prec=self%K_NaCl_form*max(0._rk,(Om_NaCl-1._rk))
   if (NaCl.gt.(3.72e7))  NaCl_prec=0.0_rk ! no precipitation if cristalls occupy more than 100% of volume
   NaCl_diss=self%K_NaCl_diss*max(0._rk,(1._rk-Om_NaCl))*NaCl

   _SET_ODE_(self%id_Na,  -NaCl_prec+NaCl_diss)
   _SET_ODE_(self%id_Cl,  -NaCl_prec+NaCl_diss)
   _SET_ODE_(self%id_NaCl,+NaCl_prec-NaCl_diss)
   _SET_ODE_(self%id_Het,-(0.5_rk+0.5_rk*tanh((Cl*35.453e-6 + Na*22.990e-6)-300._rk))*0.6_rk*Het)


   ! From mmol/m3 to g/L (note: if we want g/kg, we'd need to divide by density rather than 1000)
!   _SET_DIAGNOSTIC_(self%id_psu,(Na*u_Na+Cl*u_Cl)*1e-6_rk)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !DESCRIPTION:
!
!
! !INPUT PARAMETERS:
   class (type_niva_brom_halite),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !REVISION HISTORY:
!  Original author(s):
!
! !LOCAL VARIABLES:
   real(rk) :: Na,Cl,NaCl,NaCl_bot
   real(rk) :: Om_NaCl,NaCl_diss
!   real(rk),parameter :: max_dt = 6*3600 ! Maximum time step = minimum time needed to dissolve bottom salt layer
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _HORIZONTAL_LOOP_BEGIN_

   !! Environment
   !_GET_(self%id_Na,Na)
   !_GET_(self%id_Cl,Cl)
   !_GET_(self%id_NaCl,NaCl)
   !_GET_HORIZONTAL_(self%id_NaCl_bot,NaCl_bot)
   !
   !_SET_BOTTOM_EXCHANGE_(self%id_NaCl,-self%w*NaCl)
   !_SET_BOTTOM_ODE_(self%id_NaCl_bot, +self%w*NaCl)
   !
   !Om_NaCl=Na*Cl/(self%K_NaCl)
   !NaCl_diss=min(self%k_diss_bot,NaCl_bot/max_dt)*max(0._rk,(1._rk-Om_NaCl))
   !_SET_BOTTOM_ODE_(self%id_NaCl_bot,-NaCl_diss)
   !_SET_BOTTOM_EXCHANGE_(self%id_Na,  NaCl_diss)
   !_SET_BOTTOM_EXCHANGE_(self%id_Cl,  NaCl_diss)
   !
   !! From mmol/m2 to thickness: first to g/m2 (multiply by molecuar weight), then to m by dividing by density in g/m3, accounting for porosity (void fraction)
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_h_bot,NaCl_bot*(u_Na+u_Cl)/1000/NaCl_density/(1-self%porosity))

   ! Leave spatial loops (if any)
   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!EOC

end module
