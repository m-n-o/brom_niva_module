!-----------------------------------------------------------------------
! brom_nitrogen is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_nitrogen
  use fabm_types
  use brom_functions

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_nitrogen
    !all descriptions are in the initialize subroutine
    !state variables dependencies
    type(type_state_variable_id):: id_NH4,id_NO2,id_NO3
    type(type_state_variable_id):: id_PO4,id_Si
    type(type_state_variable_id):: id_O2
    type(type_state_variable_id):: id_POM,id_DOM
    type(type_state_variable_id):: id_Alk,id_DIC

    type(type_diagnostic_variable_id):: id_anammox
    type(type_diagnostic_variable_id):: id_Nitrif1,id_Nitrif2
    type(type_diagnostic_variable_id):: id_DcPOM_NO3,id_DcDOM_NO3
    type(type_diagnostic_variable_id):: id_DcPOM_NO2,id_DcDOM_NO2
    type(type_diagnostic_variable_id):: id_dAlk
    !Model parameters
    !specific rates of biogeochemical processes
    !----N--------!
    real(rk):: K_nitrif1,K_nitrif2,O2s_nf
    real(rk):: K_annamox,O2s_dn
    real(rk):: K_omno_no3,K_omno_no2
    real(rk):: K_POM_NO3,K_POM_NO2
    real(rk):: K_DOM_NO3,K_DOM_NO2
    !---- Stoichiometric coefficients ----!
    real(rk):: c_to_n, c_to_si, c_to_p
  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_nitrogen), intent(inout), target :: self
    integer,                         intent(in)            :: configunit

    !----Model parameters----!
    !----N---------!
    call self%get_parameter(&
         self%K_nitrif1, 'K_nitrif1', '[1/day]',&
         'Spec.rate of 1st st. of nitrification',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_nitrif2, 'K_nitrif2', '[1/day]',&
         'Spec.rate of 2d st. of nitrification',&
         default=0.1_rk)
    call self%get_parameter(&
         self%K_annamox, 'K_annamox', '[1/day]',&
         'Spec. rate of Anammox',&
         default=0.8_rk)
    call self%get_parameter(&
         self%O2s_nf, 'O2s_nf', '[uM O]',&
         'half saturation for nitrification',&
         default=4.488_rk)
    call self%get_parameter(&
         self%O2s_dn, 'O2s_dn', '[uM O]',&
         'half saturation for annamox and denitrification',&
         default=10.0_rk)
    call self%get_parameter(&
         self%K_POM_NO3, 'K_POM_NO3', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of POM',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_POM_NO2, 'K_POM_NO2', '[1/day]',&
         'Spec.rate of 2 stage of denitrif of POM',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_DOM_NO3, 'K_DOM_NO3', '[1/day]',&
         'Spec.rate of 1 stage of denitrif of DOM',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_DOM_NO2, 'K_DOM_NO2', '[1/day]',&
         'Spec.rate of 2 stage of denitrif of DOM',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_omno_no3, 'K_omno_no3', '[uM N]',&
         'half sat. of no3 for OM denitr.',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_omno_no2, 'K_omno_no2', '[uM N]',&
         'half sat. of no2 for OM denitr.',&
         default=0.001_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(self%c_to_n,'c_to_n','[-]','C[uM]/N[uM]',&
                            default=106._rk/16._rk)
    call self%get_parameter(self%c_to_si,'c_to_si','[-]','C[uM]/Si[uM]',&
                            default=106._rk/15._rk)
    call self%get_parameter(self%c_to_p,'c_to_p','[-]','C[uM]/P[uM]',&
                            default=106._rk)

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3 N',&
         'ammonium')
    call self%register_state_dependency(&
         self%id_NO2,'NO2','mmol/m**3 N',&
         'nitrite')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3 N',&
         'nitrate')
    call self%register_state_dependency(&
         self%id_PO4,'PO4','mmol/m**3 P',&
         'phosphate')
    call self%register_state_dependency(&
         self%id_Si,'Si','mmol/m**3 Si',&
         'silicate')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3 C',&
         'total dissolved inorganic carbon')
    call self%register_state_dependency(&
         self%id_Alk,'Alk','mmol/m**3',&
         'alkalinity')
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3 O2',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_POM,'POM','mg C m^-3',&
         'POM')
    call self%register_state_dependency(&
         self%id_DOM,'DOM','mg C m^-3',&
         'DOM')

    !Register diagnostic variables
    call self%register_diagnostic_variable(self%id_Nitrif1,'Nitrif1','mmol/m**3 N',&
         'Nitrification 1 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_Nitrif2,'Nitrif2','mmol/m**3 N',&
         'Nitrification 2 stage',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_anammox,'Anammox','mmol/m**3 N',&
         'Anammox',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOM_NO2,'DcPOM_NO2','mg C m^-3',&
         'POM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOM_NO2,'DcDOM_NO2','mg C m^-3',&
         'DOM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOM_NO3,'DcPOM_NO3','mg C m^-3',&
         'POM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOM_NO3,'DcDOM_NO3','mg C m^-3',&
         'DOM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_dAlk,&
         'd_alk','mM m^-3',&
         'Alkalinity generation due to nitrification and denitrification',&
         output=output_time_step_integrated)
    !call self%add_to_aggregate_variable(&
    !     standard_variables%alkalinity_expressed_as_mole_equivalent,self%id_dAlk)

    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_nitrogen),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: O2,POM,DOM
    real(rk):: NO2,NO3,NH4
    !increments
    real(rk):: d_O2,d_DOM,d_POM
    real(rk):: d_NO2,d_NO3,d_NH4
    real(rk):: d_PO4,d_Si
    real(rk):: d_alk,d_DIC
    !processes
    real(rk):: Nitrif1,Nitrif2,Anammox
    real(rk):: denitrification_1, denitrification_2
    real(rk):: DcPOM_NO3,DcDOM_NO3,DcPOM_NO2,DcDOM_NO2

    real(rk):: dPOM_no3_in_m,dPOM_no2_in_m
    real(rk):: dDOM_no3_in_m,dDOM_no2_in_m

    real(rk):: kf

    _LOOP_BEGIN_
      !Retrieve current variable values
      !state
      _GET_(self%id_DOM,DOM)
      _GET_(self%id_NO2,NO2)
      _GET_(self%id_NO3,NO3)
      !solids
      _GET_(self%id_POM,POM)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_NH4,NH4)

      !
      if (o2 < 0.5_rk) then
        Nitrif1 = 0._rk
        Nitrif2 = 0._rk
      else
        !Nitrification 1st stage: NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O
        Nitrif1 = self%K_nitrif1*NH4*hyper_limiter(self%O2s_nf, o2, 1._rk)
        !Nitrification 2d stage: NO2- + 0.5 O2 -> NO3-
        Nitrif2 = self%K_nitrif2*NO2*hyper_limiter(self%O2s_nf, o2, 1._rk)
      end if
      !
      !Anammox NO2- + NH4+ -> N2 + 2H2O
      Anammox = self%K_annamox*NO2*NH4*hyper_inhibitor(self%O2s_dn, o2, 1._rk)
      !
      !OM denitrification (Richards, 1965)
      !(CH2O)106(NH3)16H3PO4 + 84.8HNO3 = 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4
      !
      !OM denitrification (1st stage) (Anderson,1982)
      !1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2
      !It should be in the units of OM, so mg C m^-3
      kf = monod_squared(self%K_omno_no3, NO3)&
          *hyper_inhibitor(self%O2s_dn, o2, 1._rk)
      !denitrification1
      DcPOM_NO3 = self%K_POM_NO3*POM*kf
      DcDOM_NO3 = self%K_DOM_NO3*DOM*kf
      !recalculate to moles
      dPOM_no3_in_m = carbon_g_to_mole(DcPOM_NO3)
      dDOM_no3_in_m = carbon_g_to_mole(DcDOM_NO3)
      !
      !OM denitrification (2d stage)
      !3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,1982)
      !It should be in the units of OM, so mg C m^-3
      kf = monod_squared(self%K_omno_no2, NO2)&
          *hyper_inhibitor(self%O2s_dn, o2, 1._rk)&
          *hyper_inhibitor(self%K_omno_no3, no3, 1._rk)
      !denitrification2
      DcPOM_NO2 = self%K_POM_NO2*POM*kf
      DcDOM_NO2 = self%K_DOM_NO2*DOM*kf
      !recalculate to moles
      dPOM_no2_in_m = carbon_g_to_mole(DcPOM_NO2)
      dDOM_no2_in_m = carbon_g_to_mole(DcDOM_NO2)

      !amount of NO3 consumed / NO2 excreted
      denitrification_1 = 2._rk*(dDOM_no3_in_m+dPOM_no3_in_m)
      !amount of NO2 consumed
      denitrification_2 = (4._rk/3._rk)*(dDOM_no2_in_m+dPOM_no2_in_m)

      !to prevent negative values of NO2 after summation outside FABM
      if (NO3 < denitrification_1/self%dt*300._rk) then
        dPOM_no3_in_m = 0._rk; DcPOM_NO3 = 0._rk
        dDOM_no3_in_m = 0._rk; DcDOM_NO3 = 0._rk
        denitrification_1 = 0._rk
      end if
      if (NO2 < denitrification_2/self%dt*300._rk) then
        dPOM_no2_in_m = 0._rk; DcPOM_NO2 = 0._rk
        dDOM_no2_in_m = 0._rk; DcDOM_NO2 = 0._rk
        denitrification_2 = 0._rk
      end if
      !
      !Set increments
      d_DOM = -DcDOM_NO3-DcDOM_NO2
      d_POM = -DcPOM_NO3-DcPOM_NO2

      d_NO2 = Nitrif1-Nitrif2-Anammox&
             +denitrification_1&
             -denitrification_2

      d_NO3 = Nitrif2&
             -denitrification_1

      d_NH4 =-Nitrif1-Anammox&
             +(dDOM_no3_in_m&
              +dPOM_no3_in_m&
              +dDOM_no2_in_m&
              +dPOM_no2_in_m)/self%c_to_n

      d_PO4 =+(dDOM_no3_in_m&
              +dPOM_no3_in_m&
              +dDOM_no2_in_m&
              +dPOM_no2_in_m)/self%c_to_p

      d_Si  =+(dDOM_no3_in_m&
              +dPOM_no3_in_m&
              +dDOM_no2_in_m&
              +dPOM_no2_in_m)/self%c_to_si

      d_O2  =-1.5_rk*Nitrif1-0.5_rk*Nitrif2

      d_DIC = dDOM_no3_in_m+dPOM_no3_in_m&
             +dDOM_no2_in_m+dPOM_no2_in_m

      d_alk = d_NH4-d_NO3-d_NO2-d_PO4

      _SET_ODE_(self%id_NO2,d_NO2)
      _SET_ODE_(self%id_NO3,d_NO3)
      _SET_ODE_(self%id_NH4,d_NH4)
      _SET_ODE_(self%id_PO4,d_PO4)
      _SET_ODE_(self%id_Si,d_Si)

      _SET_ODE_(self%id_O2,d_O2)
      _SET_ODE_(self%id_DIC,d_DIC)

      _SET_ODE_(self%id_Alk,d_Alk)

      _SET_ODE_(self%id_DOM,d_DOM)
      _SET_ODE_(self%id_POM,d_POM)

      _SET_DIAGNOSTIC_(self%id_anammox,Anammox)
      _SET_DIAGNOSTIC_(self%id_Nitrif1,Nitrif1)
      _SET_DIAGNOSTIC_(self%id_Nitrif2,Nitrif2)
      _SET_DIAGNOSTIC_(self%id_DcPOM_NO2,DcPOM_NO2)
      _SET_DIAGNOSTIC_(self%id_DcDOM_NO2,DcDOM_NO2)
      _SET_DIAGNOSTIC_(self%id_DcPOM_NO3,DcPOM_NO3)
      _SET_DIAGNOSTIC_(self%id_DcDOM_NO3,DcDOM_NO3)
      _SET_DIAGNOSTIC_(self%id_dAlk,d_alk)
    _LOOP_END_
  end subroutine do
end module fabm_niva_brom_nitrogen
