!-----------------------------------------------------------------------
! brom_carbon is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_carbon
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_carbon
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_DIC,id_Alk

    type(type_diagnostic_variable_id):: id_pCO2
    type(type_diagnostic_variable_id):: id_CO3
    type(type_diagnostic_variable_id):: id_dAlk
    !diagnostic dependencies
    type(type_dependency_id):: id_Hplus,id_Kc0,id_Kc1,id_Kc2
    !standard variable dependincies
    type(type_dependency_id):: id_dens,id_temp,id_salt
    !for do_surface
    type(type_dependency_id):: id_pCO2w
    type(type_horizontal_dependency_id):: id_windspeed,id_pCO2a
  contains
    procedure :: initialize
    procedure :: do_surface
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class(type_niva_brom_carbon),intent(inout),target :: self
    integer,                     intent(in)           :: configunit

    call self%register_state_variable(&
         self%id_DIC,'DIC','mM m^-3','DIC',minimum=0.0_rk)
    call self%register_state_variable(&
         self%id_Alk,'Alk','mM m^-3','Alk',2300._rk,minimum=1.e-4_rk)

    !register diagnostic variables
    !mol mol-1 and atm are the same
    call self%register_diagnostic_variable(&
         self%id_pCO2,'pCO2','mol mol-1','CO2, mole fraction')
    call self%register_diagnostic_variable(&
         self%id_CO3,'CO3','mol kg-1','CO3 concentration')

    !Register standard dependencies
    call self%register_dependency(&
         self%id_dens,standard_variables%density)
    call self%register_dependency(&
         self%id_temp,standard_variables%temperature)
    call self%register_dependency(&
         self%id_salt,standard_variables%practical_salinity)

    !dependencies
    !for do_surface
    call self%register_dependency(&
         self%id_pco2a,&
         standard_variables%mole_fraction_of_carbon_dioxide_in_air)
    call self%register_dependency(&
         self%id_windspeed,standard_variables%wind_speed)
    call self%register_dependency(&
         self%id_pCO2w,'pCO2','mol mol-1','CO2, mole fraction')

    !diagnostic variables
    call self%register_dependency(self%id_Hplus,'Hplus','mol kg^-1','H+')
    call self%register_dependency(&
         self%id_Kc0,'Kc0','-','Henry''s constant')
    call self%register_dependency(&
         self%id_Kc1,'Kc1','-','[H+][HCO3-]/[H2CO3]')
    call self%register_dependency(&
         self%id_Kc2,'Kc2','-','[H+][CO3--]/[HCO3-]')
  end subroutine initialize
  !
  !Sea water CO2 exchange, adapted from PML's ERSEM code
  !
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    class (type_niva_brom_carbon),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_SURFACE_
    real(rk):: dens
    real(rk):: Q_DIC
    real(rk):: temp, salt, Kc0
    real(rk):: Sc, fwind
    real(rk):: pCO2w, pCO2a
    real(rk):: windspeed

    _HORIZONTAL_LOOP_BEGIN_
      ! Environment
      _GET_(self%id_dens,dens) !density
      _GET_(self%id_temp,temp) !temperature
      _GET_(self%id_salt,salt) !salinity
      !previouss pCO2 which calculates here in the subroutine do
      _GET_(self%id_pCO2w,pCO2w) !in atm = (mol mol-1)
      !Equilibrium constants
      _GET_(self%id_Kc0,Kc0)
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)
      _GET_HORIZONTAL_(self%id_pCO2a,pCO2a) !from atmosphere, microM/M

      !PML
      !calculate the scmidt number and unit conversions(only temp<30)
      Sc = 2073.1_rk-125.62_rk*temp+3.6276_rk*temp**2._rk-0.043219_rk*&
           temp**3.0_rk
      fwind = (0.222_rk*windspeed**2_rk+0.333_rk*windspeed)*&
              (Sc/660._rk)**(-0.5_rk)
      fwind=fwind*24._rk/100._rk !convert to m/day

      !flux depends on the difference in partial pressures, wind and henry
      !we shall convert pCO2a from microM/M to M/M
      Q_DIC = fwind*Kc0*(pCO2a/1.e6_rk-max(0e0,pCO2w))
      Q_DIC = Q_DIC*1.e3_rk*dens !convert to mmol C m^-2 d^-1
      Q_DIC = Q_DIC/86400._rk !to mmol C m^-2 s^-1, self%dt is not specified

      _SET_SURFACE_EXCHANGE_(self%id_DIC,Q_DIC)
    _HORIZONTAL_LOOP_END_
  end subroutine do_surface
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_carbon),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    real(rk):: dens
    !state variables
    real(rk):: DIC
    !diagnostic variables
    real(rk):: Kc0,Kc1,Kc2
    real(rk):: H_
    real(rk):: hco3,co3,co2,pCO2

    _LOOP_BEGIN_
      ! Environment
      _GET_(self%id_dens,dens) !density
      !state variables
      _GET_(self%id_DIC,DIC)
      !turn units to mol/kg
      DIC = DIC*(1000._rk/dens)*1.e-6_rk
      !diagnostic variables
      _GET_(self%id_Hplus,H_)
      !Equilibrium constants
      _GET_(self%id_Kc0,Kc0)
      _GET_(self%id_Kc1,Kc1)
      _GET_(self%id_Kc2,Kc2)

      !Calculate all the others as a function of
      ![H+], DIC and constants, moles/kg
      hco3 = DIC/(1._rk+H_/Kc1+Kc2/H_)
      co3  = DIC/(1._rk+H_/Kc2+H_*H_/Kc1/Kc2)
      co2  = DIC/(1._rk+Kc1/H_+Kc1*Kc2/H_/H_)
      pco2 = co2/Kc0 !atm

      _SET_DIAGNOSTIC_(self%id_pCO2,pCO2)
      _SET_DIAGNOSTIC_(self%id_CO3,co3)
    _LOOP_END_
  end subroutine do
end module
