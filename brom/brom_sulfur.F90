!-----------------------------------------------------------------------
! brom_sulfur is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_sulfur
  use brom_functions
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_sulfur
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_H2S
    type(type_state_variable_id):: id_S0,id_S2O3,id_SO4
    !state variables dependencies
    type(type_state_variable_id):: id_O2,id_NO3,id_NH4,id_PO4,id_Si
    type(type_state_variable_id):: id_POML,id_POMR,id_DOML,id_DOMR
    type(type_state_variable_id):: id_DIC,id_Alk
    !standard variables
    type(type_dependency_id):: id_temp
    !diagnostic variables - processes
    type(type_diagnostic_variable_id):: id_DcPOMR_so4, id_DcDOMR_so4
    type(type_diagnostic_variable_id):: id_DcPOML_so4, id_DcDOML_so4
    type(type_diagnostic_variable_id):: id_s2o3_no3,id_s0_no3
    type(type_diagnostic_variable_id):: id_s0_ox,id_s2o3_ox
    type(type_diagnostic_variable_id):: id_s0_disp,id_hs_ox
    type(type_diagnostic_variable_id):: id_hs_no3
    type(type_diagnostic_variable_id):: id_dAlk
    !Model parameters
    !specific rates of biogeochemical processes
    real(rk):: K_s0_disp,K_hs_ox,K_s0_ox,K_s0_no3
    real(rk):: K_POMR_so4,K_POML_so4
    real(rk):: K_DOMR_so4,K_DOML_so4
    real(rk):: K_s2o3_ox,K_s2o3_no3,K_hs_no3
    real(rk):: tref
    !---- Switches-------!
    real(rk):: s_omso_o2,s_omso_no3
    real(rk):: K_omso_so4,O2_sr,NO3_sr
    !---- Stoichiometric coefficients ----!
    real(rk):: c_to_n, c_to_si, c_to_p
  contains
    procedure :: initialize
    procedure :: do
  end type
  contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_sulfur), intent(inout), target :: self
    integer,                       intent(in)            :: configunit

    !-----Model parameters------
    !---- S---------!
    call self%get_parameter(&
            self%K_s0_disp,'K_s0_disp','[1/day]',&
            'Specific rate of S0 dispropotionation',&
            default=0.001_rk)
    call self%get_parameter(&
            self%K_hs_ox,'K_hs_ox','[1/day]',&
            'Specific rate of oxidation of H2S to S0 with O2',&
            default=0.5_rk)
    call self%get_parameter(&
            self%K_hs_no3, 'K_hs_no3', '[1/day]',&
            'Specific rate of thiodenitrification.',&
            default=0.8_rk)
    call self%get_parameter(&
            self%K_s0_ox,'K_s0_ox','[1/day]',&
            'Specific rate of oxidation of S0 with O2',&
            default=0.02_rk)
    call self%get_parameter(&
            self%K_s0_no3,'K_s0_no3','[1/day]',&
            'Specific rate of oxidation of S0 with NO3',&
            default=0.9_rk)
    call self%get_parameter(&
            self%K_s2o3_ox,'K_s2o3_ox','[1/day]',&
            'Specific rate of oxidation of S2O3 with O2',&
            default=0.01_rk)
    call self%get_parameter(&
            self%K_s2o3_no3,'K_s2o3_no3','[1/day]',&
            'Specific rate of oxidation of S2O3 with NO3',&
            default=0.01_rk)
    call self%get_parameter(&
            self%K_POML_so4,'K_POML_so4','[1/day]',&
            'Specific rate of OM sulfate reduction with sulfate',&
            default=0.000005_rk)
    call self%get_parameter(&
            self%K_POMR_so4,'K_POMR_so4','[1/day]',&
            'Specific rate of POMR sulfate reduction with sulfate',&
            default=0.000005_rk)
    call self%get_parameter(&
            self%K_DOML_so4,'K_DOML_so4','[1/day]',&
            'Specific rate of DOML sulfate reduction with sulfate',&
            default=0.000005_rk)
    call self%get_parameter(&
            self%K_DOMR_so4,'K_DOMR_so4','[1/day]',&
            'Specific rate of DOMR sulfate reduction with sulfate',&
            default=0.000005_rk)
    call self%get_parameter(&
         self%tref,'tref','degrees Celsius',&
         'Reference temperature at which temperature factor = 1',&
         default=0.0_rk)
    !----Inhibitors----!
    call self%get_parameter(&
            self%s_omso_o2, 's_omso_o2', '[uM O2]',&
            'Threshold of O2 for OM sulfate reduction',&
            default=25.0_rk)
    call self%get_parameter(&
            self%s_omso_no3, 's_omso_no3', '[uM N]',&
            'Threshold of NO3 for OM sulfate reduction',&
            default=5.0_rk)
    !----Limiters------!
    call self%get_parameter(&
            self%K_omso_so4, 'K_omso_so4', '[uM S]',&
            'Half-saturation of SO4 for OM sulfate reduction',&
            default=1.0_rk)
    call self%get_parameter(&
            self%O2_sr, 'O2_sr', '[uM O2]',&
            'Half-saturation of O2',&
            default=1.0_rk)
    call self%get_parameter(&
            self%NO3_sr, 'NO3_sr', '[uM N]',&
            'Half-saturation of NO3',&
            default=1.0_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(self%c_to_n,'c_to_n','[-]','C[uM]/N[uM]',&
                            default=106._rk/16._rk)
    call self%get_parameter(self%c_to_si,'c_to_si','[-]','C[uM]/Si[uM]',&
                            default=106._rk/15._rk)
    call self%get_parameter(self%c_to_p,'c_to_p','[-]','C[uM]/P[uM]',&
                            default=106._rk)

    !Register state variables
    call self%register_state_variable(&
            self%id_H2S, 'H2S', 'mmol/m**3','H2S',&
            minimum=0.0_rk)
    call self%register_state_variable(&
            self%id_S0 , 'S0',  'mmol/m**3','S0',&
            minimum=0.0_rk)
    call self%register_state_variable(&
            self%id_S2O3,'S2O3','mmol/m**3','S2O3',&
            minimum=0.0_rk)
    call self%register_state_variable(&
            self%id_SO4, 'SO4', 'mmol/m**3','SO4',&
            minimum=0.0_rk)

    !Register state dependencies
    call self%register_state_dependency(&
            self%id_DIC,'DIC','mmol/m**3',&
            'total dissolved inorganic carbon')
    call self%register_state_dependency(&
         self%id_Alk,'Alk','mmol/m**3',&
         'alkalinity')
    call self%register_state_dependency(&
         self%id_Si,'Si','mmol/m**3 Si',&
         'silicate')
    call self%register_state_dependency(&
            self%id_PO4,'PO4','mmol/m**3',&
            'phosphate')
    call self%register_state_dependency(&
            self%id_O2, 'O2', 'mmol/m**3',&
            'dissolved oxygen')
    call self%register_state_dependency(&
            self%id_NH4,'NH4','mmol/m**3',&
            'ammonium')
    call self%register_state_dependency(&
            self%id_NO3,'NO3','mmol/m**3',&
            'nitrate')
    call self%register_state_dependency(&
            self%id_POML,'POML','mmol/m**3',&
            'particulate organic nitrogen')
    call self%register_state_dependency(&
            self%id_POMR,'POMR','mmol/m**3',&
            'POM refractory')
    call self%register_state_dependency(&
            self%id_DOML,'DOML','mmol/m**3',&
            'dissolved organic nitrogen')
    call self%register_state_dependency(&
            self%id_DOMR,'DOMR','mmol/m**3',&
            'POM refractory')
    !Register standard variables
    call self%register_dependency(&
         self%id_temp,standard_variables%temperature)

    !Register diagnostic variables
    call self%register_diagnostic_variable(&
            self%id_DcPOMR_SO4,'DcPOMR_SO4','mmol/m**3',&
            'POMR sulfatereduction',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_DcDOMR_SO4,'DcDOMR_SO4','mmol/m**3',&
            'DOMR sulfatereduction',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_DcPOML_SO4,'DcPOML_SO4','mmol/m**3',&
            'POML sulfatereduction',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_DcDOML_SO4,'DcDOML_SO4','mmol/m**3',&
            'DOML sulfatereduction',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_s2o3_no3,'s2o3_no3','mmol/m**3',&
            'S2O3 with NO3 oxidation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_s0_no3,'s0_no3','mmol/m**3',&
            'S0 with NO3 oxidation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_s2o3_ox,'s2o3_ox','mmol/m**3',&
            'S2O3  with O2oxidation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_s0_ox,'s0_ox','mmol/m**3',&
            'S0 with O2 oxidation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_s0_disp,'s0_disp','mmol/m**3',&
            'S0 disproportionation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_hs_ox,'hs_ox','mmol/m**3',&
            'H2S with O2 oxidation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
            self%id_hs_no3,'hs_no3','mmol/m**3',&
            'H2S with NO3 oxidation',&
            output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_dAlk,&
         'd_alk','mM m^-3',&
         'Alkalinity generation due to sulfur transformation',&
         output=output_time_step_integrated)

    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_sulfur),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !standard variables
    real(rk):: temp
    !state variables
    real(rk):: H2S,S0,S2O3,SO4
    !state dependencies
    real(rk):: o2
    real(rk):: POML,POMR,DOML,DOMR
    real(rk):: NO3
    !increments
    real(rk):: d_S2O3,d_SO4,d_S0,d_H2S
    real(rk):: d_O2,d_DOML,d_POML,d_POMR,d_DOMR
    real(rk):: d_NH4,d_NO3,d_DIC,d_PO4,d_Si
    real(rk):: d_alk
    !processess
    real(rk):: s0_disp,hs_ox,s0_ox,s0_no3,s2o3_ox,s2o3_no3,hs_no3
    real(rk):: DcPOML_so4,DcDOML_so4,DcPOMR_so4,DcDOMR_so4
    real(rk):: dpoml_so4_in_m, ddoml_so4_in_m
    real(rk):: dpomr_so4_in_m, ddomr_so4_in_m
    !parameters
    real(rk):: thr_no3,thr_o2
    real(rk):: kf

    _LOOP_BEGIN_
      !Retrieve current state variable values
      _GET_(self%id_temp,temp) ! temperature
      !state
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_DOMR,DOMR)
      _GET_(self%id_S2O3,S2O3)
      _GET_(self%id_SO4,SO4)
      _GET_(self%id_NO3,NO3)
      !solids
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)
      _GET_(self%id_S0,S0)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_H2S,H2S)

      !S
      !S0 disportionation: 4S0 + 3H2O -> 2H2S + S2O3= + 2H+
      s0_disp = self%K_s0_disp*S0
      if (O2 < 0.1_rk) then
        hs_ox = 0._rk
        s0_ox = 0._rk
        s2o3_ox = 0._rk
      else
        !HS oxidation with O2: 2H2S + O2 -> 2S0 + 2H2O
        hs_ox = self%K_hs_ox*H2S*O2*hyper_limiter(self%O2_sr, O2, 1._rk) &
               *hyper_limiter(1e-3_rk, H2S, 1._rk)
        !S0 oxidation with O2: 2S0 + O2 + H2O -> S2O3= + 2H+
        s0_ox = self%K_s0_ox*O2*S0*hyper_limiter(self%O2_sr, O2, 1._rk)
        !S2O3 oxidation with O2: S2O3= + 2O2 + 2OH- -> 2SO4= + H2O
        s2o3_ox = self%K_s2o3_ox*O2*S2O3*hyper_limiter(self%O2_sr, O2, 1._rk)
      end if
      !
      if (NO3 < 0.1_rk) then
        s0_no3 = 0._rk
        s2o3_no3 = 0._rk
        hs_no3 = 0._rk
      else
        !S0 oxidation with NO3: 4S0 + 3NO3- + 7H2O -> 4SO4= + 3NH4+ + 2H+
        s0_no3 = self%K_s0_no3*NO3*S0*hyper_limiter(self%NO3_sr, NO3, 1._rk)
        !S2O3 oxidation with NO3: S2O3= + NO3- + 2H2O --> 2SO4= + NH4+
        s2o3_no3 = self%K_s2o3_no3*NO3*S2O3*hyper_limiter(self%NO3_sr, NO3, 1._rk)
        !Thiodenitrification: 3H2S + 4NO3- + 6OH- -> 3SO4= + 2N2 + 6H2O
        hs_no3 = self%K_hs_no3*H2S*NO3*hyper_limiter(self%NO3_sr, NO3, 1._rk) &
                *hyper_limiter(1e-3_rk, H2S, 1._rk)
      end if

      !in anoxic conditions:
      thr_O2 = hyper_inhibitor(self%s_omso_O2, O2, 1._rk)
      thr_no3 = hyper_inhibitor(self%s_omso_no3, no3, 1._rk)
      kf = monod_squared(self%K_omso_so4, SO4)*f_t(temp,2._rk,self%tref)
      !OM sulfatereduction (Boudreau, 1996)
      !(CH2O)106(NH3)16H3PO4 + 53SO42- = 106HCO3- + 16NH3 + H3PO4 + 53H2S
      !It should be in the units of OM, so mg C m^-3
      !POML sulfatereduction (1st stage):
      DcPOML_so4 = thr_O2*thr_no3 &
                  *self%K_POML_so4*POML*kf
      !DOML sulfatereduction (1st stage):
      DcDOML_so4 = thr_O2*thr_no3 &
                  *self%K_DOML_so4*DOML*kf
      !recalculate to moles
      dpoml_so4_in_m = carbon_g_to_mole(DcPOML_so4)
      ddoml_so4_in_m = carbon_g_to_mole(DcDOML_so4)

      !POMR sulfatereduction (1st stage):
      DcPOMR_so4  = thr_O2*thr_no3 &
                   *self%K_POMR_so4*POMR*kf
      !DOMR sulfatereduction (1st stage):
      DcDOMR_so4  = thr_O2*thr_no3 &
                   *self%K_DOMR_so4*DOMR*kf
      !recalculate to moles
      dpomr_so4_in_m = carbon_g_to_mole(DcPOMR_so4)
      ddomr_so4_in_m = carbon_g_to_mole(DcDOMR_so4)

      if (SO4 < 0.5_rk*(dpomr_so4_in_m+ddomr_so4_in_m)/self%dt*300._rk) then
        dpomr_so4_in_m = 0._rk; DcPOMR_so4 = 0._rk
        ddomr_so4_in_m = 0._rk; DcDOMR_so4 = 0._rk
        dpoml_so4_in_m = 0._rk; DcPOML_so4 = 0._rk
        ddoml_so4_in_m = 0._rk; DcDOML_so4 = 0._rk
      end if

      !Set increments
      d_SO4 = hs_no3+2._rk*s2o3_ox+s0_no3+2._rk*s2o3_no3 &
             -0.5_rk*(dpomr_so4_in_m+ddomr_so4_in_m)
      _SET_ODE_(self%id_SO4,d_SO4)
      d_S2O3 = 0.5_rk*s0_ox-s2o3_ox+0.25_rk*s0_disp-s2o3_no3
      _SET_ODE_(self%id_S2O3,d_S2O3)
      d_S0 = hs_ox-s0_ox-s0_disp-s0_no3
      _SET_ODE_(self%id_S0,d_S0)
      d_H2S = -hs_ox+0.5_rk*s0_disp-hs_no3 &
              +0.5_rk*(dpomr_so4_in_m+ddomr_so4_in_m)
      _SET_ODE_(self%id_H2S,d_H2S)

      d_DIC = dpomr_so4_in_m+ddomr_so4_in_m
      _SET_ODE_(self%id_DIC,d_DIC)
      d_O2 = -0.5_rk*hs_ox-0.5_rk*s0_ox-2._rk*s2o3_ox
      _SET_ODE_(self%id_O2,d_O2)

      d_DOML = -DcDOML_so4
      _SET_ODE_(self%id_DOML,d_DOML)
      d_POML = -DcPOML_so4
      _SET_ODE_(self%id_POML,d_POML)
      d_POMR = DcPOML_so4-DcPOMR_so4
      _SET_ODE_(self%id_POMR,d_POMR)
      d_DOMR = DcPOML_so4-DcDOMR_so4
      _SET_ODE_(self%id_DOMR,d_DOMR)

      d_NO3 = -(4._rk/3._rk)*hs_no3-0.75_rk*s0_no3-s2o3_no3
      _SET_ODE_(self%id_NO3,d_NO3)
      d_NH4 = 0.75_rk*s0_no3+s2o3_no3 &
             +(dpoml_so4_in_m+ddoml_so4_in_m)/self%c_to_n
      _SET_ODE_(self%id_NH4,d_NH4)
      d_PO4 = (dpoml_so4_in_m+ddoml_so4_in_m)/self%c_to_p
      _SET_ODE_(self%id_PO4,d_PO4)
      d_Si = (dpoml_so4_in_m+ddoml_so4_in_m)/self%c_to_si
      _SET_ODE_(self%id_Si,d_Si)

      d_alk = d_NH4-d_NO3-d_PO4-2._rk*d_SO4
      _SET_ODE_(self%id_Alk,d_alk)

      _SET_DIAGNOSTIC_(self%id_s2o3_no3,s2o3_no3)
      _SET_DIAGNOSTIC_(self%id_s0_no3,s0_no3)
      _SET_DIAGNOSTIC_(self%id_s2o3_ox,s2o3_ox)
      _SET_DIAGNOSTIC_(self%id_s0_ox,s0_ox)
      _SET_DIAGNOSTIC_(self%id_DcPOML_so4,DcPOML_so4)
      _SET_DIAGNOSTIC_(self%id_DcDOML_so4,DcDOML_so4)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_so4,DcPOMR_so4)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_so4,DcDOMR_so4)
      _SET_DIAGNOSTIC_(self%id_s0_disp,s0_disp)
      _SET_DIAGNOSTIC_(self%id_hs_ox,hs_ox)
      _SET_DIAGNOSTIC_(self%id_hs_no3,hs_no3)
      _SET_DIAGNOSTIC_(self%id_dAlk,d_alk)
    _LOOP_END_
  end subroutine do
end module fabm_niva_brom_sulfur
