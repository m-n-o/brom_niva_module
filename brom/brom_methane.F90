!-----------------------------------------------------------------------
! BROM2 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM2 distribution.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_methane
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_methane
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_CH4
    !state dependencies
    type(type_state_variable_id):: id_DIC,id_NH4,id_PO4
    type(type_state_variable_id):: id_DON,id_PON
    type(type_state_variable_id):: id_O2,id_NO3,id_SO4

    !diagnostic variables by bacteria needed
    type(type_diagnostic_variable_id):: id_DcPM_ch4

    !Model parameters
    !specific rates of biogeochemical processes
    !Methan
    real(rk):: s_omso_o2,s_omso_no3,s_omch_so4
    real(rk):: K_DON_ch4,K_PON_ch4
    real(rk):: K_ch4_o2,K_ch4_so4
    !---- Stoichiometric coefficients ----!
    real(rk):: r_c_n, r_n_p

  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class(type_niva_brom_methane),intent(inout),target :: self
    integer,                      intent(in)           :: configunit

    !-----Model parameters------
    !Specific rates of biogeochemical processes
    !Methane
    call self%get_parameter(&
         self%s_omso_o2, 's_omso_o2', '[uM O]',&
         'threshold of o2 for OM sulfate reduction',&
         default=25.0_rk)
    call self%get_parameter(&
         self%s_omso_no3, 's_omso_no3', '[uM N]',&
         'threshold of noX for OM sulfate reduction',&
         default=5.0_rk)
    call self%get_parameter(&
         self%s_omch_so4, 's_omch_so4', '[uM S]',&
         'threshold of SO4 for methane production from OM',&
         default=15000.0_rk)
    call self%get_parameter(&
         self%K_DON_ch4, 'K_DON_ch4', '[1/day]',&
         'Specific rate of methane production from DON',&
         default=0.00014_rk)
    call self%get_parameter(&
         self%K_PON_ch4, 'K_PON_ch4', '[1/day]',&
         'Specific rate of methane production from PON',&
         default=0.00014_rk)
    call self%get_parameter(&
         self%K_ch4_o2, 'K_ch4_o2', '[1/day]',&
         'Specific rate of oxidation of CH4 with O2',&
         default=0.14_rk)
    call self%get_parameter(&
         self%K_ch4_so4, 'K_ch4_so4', '[1/day]',&
         'Specific rate of anoxic oxidation of CH4 with SO4',&
         default=0.0000274_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(&
         self%r_c_n,   'r_c_n',  '[-]',&
         'C[uM]/N[uM]',&
         default=6.625_rk)
    call self%get_parameter(&
         self%r_n_p,   'r_n_p',  '[-]',&
         'N[uM]/P[uM]',&
         default=16.0_rk)
    !register state variables
    call self%register_state_variable(&
         self%id_CH4,'CH4','mmol/m**3','CH4',&
         minimum=0.0_rk)

    !Register state dependencies
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mmol/m**3','DIC')
    call self%register_state_dependency(&
         self%id_O2,'O2','mmol/m**3',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mmol/m**3',&
         'nitrate')
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3',&
         'phosphate',required=.false.)
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mmol/m**3',&
         'ammonium')
    call self%register_state_dependency(&
         self%id_PON,'PON','mmol/m**3',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_DON,'DON','mmol/m**3',&
         'dissolved organic nitrogen')
    call self%register_state_dependency(&
         self%id_SO4,'SO4','mmol/m**3','sulphate')

    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPM_ch4,'DcPM_ch4','mmol/m**3',&
         'CH4 production from PON and DON',output=output_time_step_integrated)

    !Specify that rates are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !do we need to add OM, Alk, ammonium and PO4?
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_methane),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: DIC,DON,SO4,NO3,NH4,PO4
    real(rk):: O2,CH4
    real(rk):: PON
    !processes
    real(rk):: DcDM_CH4,DcPM_CH4
    real(rk):: ch4_o2,ch4_so4
    real(rk):: Dc_OM_total
    !increments
    real(rk):: d_SO4,d_O2,d_CH4,d_NH4,d_DIC,d_PO4

    _LOOP_BEGIN_
      !state variables
      _GET_(self%id_DON,DON)
      _GET_(self%id_SO4,SO4)
      _GET_(self%id_NO3,NO3)
      _GET_(self%id_NH4,NH4)
      _GET_(self%id_PO4,PO4)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_CH4,CH4)
      !solids
      _GET_(self%id_PON,PON)

      !C
      !
      !CH4 production from PON and DON
      !(CH2O)106(NH3)16H3PO4 -> 53 CO2 + 53 CH4 + 16 NH3 + H3PO4
      DcDM_CH4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                *(1._rk-0.5_rk*(1._rk+tanh(SO4-self%s_omch_so4))) &
                * self%K_DON_ch4*DON
      DcPM_CH4 = (1._rk-0.5_rk*(1._rk+tanh(o2-self%s_omso_o2))) &
                *(1._rk-0.5_rk*(1._rk+tanh(NO3-self%s_omso_no3))) &
                *(1._rk-0.5_rk*(1._rk+tanh(SO4-self%s_omch_so4))) &
                * self%K_PON_ch4*PON
      !
      !CH4 oxidation with O2
      !CH4 + 2O2 = CO2 + 2H2O
      ch4_o2 = self%K_ch4_o2*CH4*O2
      !
      !CH4 anoxic oxidation with SO4
      !CH4 + SO42- + 2 H+  = CO2 + H2S + 2H2O
      ch4_so4 = self%K_ch4_so4*CH4*SO4

      !Summariazed OM mineralization
      Dc_OM_total = 0.5_rk*(DcDM_ch4+DcPM_ch4)

      !Set increments
      !O2
      d_O2 = -2._rk*ch4_o2
      _SET_ODE_(self%id_O2,d_O2)
      !S
      d_SO4 = -ch4_so4
      _SET_ODE_(self%id_SO4,d_SO4)
      !DIC
      d_DIC = Dc_OM_total*self%r_c_n
      _SET_ODE_(self%id_DIC,d_DIC)
      !Methane
      d_CH4 = 0.5_rk*(DcDM_CH4+DcPM_CH4)-ch4_o2-ch4_so4
      _SET_ODE_(self%id_CH4,d_CH4)
      !NH4
      d_NH4 = Dc_OM_total
      _SET_ODE_(self%id_NH4,d_NH4)
      !P
      d_PO4 = Dc_OM_total/self%r_n_p
      _SET_ODE_(self%id_PO4,d_PO4)
      
      _SET_DIAGNOSTIC_(self%id_DcPM_ch4,DcPM_ch4)
    _LOOP_END_
  end subroutine do
end module
