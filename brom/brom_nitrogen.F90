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

module fabm_niva_brom_nitrogen
  use fabm_types

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_nitrogen
    !all descriptions are in the initialize subroutine
    !state variables dependencies
    type(type_state_variable_id):: id_NH4,id_NO2,id_NO3
    type(type_state_variable_id):: id_O2
    type(type_state_variable_id):: id_PON,id_DON
    type(type_state_variable_id):: id_DIC,id_Alk
    type(type_state_variable_id):: id_PO4

    type(type_diagnostic_variable_id):: id_DcPM_NOX,id_DcDM_NOX
    type(type_diagnostic_variable_id):: id_anammox
    type(type_diagnostic_variable_id):: id_Nitrif1, id_Nitrif2
    type(type_diagnostic_variable_id):: id_Denitr1_PM,id_Denitr1_DM
    type(type_diagnostic_variable_id):: id_Denitr2_PM,id_Denitr2_DM
    type(type_diagnostic_variable_id):: id_Denitr1,id_Denitr2
    !Model parameters
    !specific rates of biogeochemical processes
    !----N--------!
    real(rk):: K_nitrif1,K_nitrif2,O2s_nf
    real(rk):: K_annamox,O2s_dn
    real(rk):: K_denitr1,K_omno_no3
    real(rk):: K_denitr2,K_omno_no2
    !---- Stoichiometric coefficients ----!
    real(rk):: r_c_n,r_n_p
  contains
    procedure :: initialize
    procedure :: do
  end type
contains
  !
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
         self%K_denitr1, 'K_denitr1', '[1/day]',&
         'Spec.rate of 1 stage of denitrif',&
         default=0.20_rk)
    call self%get_parameter(&
         self%K_denitr2, 'K_denitr2', '[1/day]',&
         'Spec.rate of 2 stage of denitrif',&
         default=0.25_rk)
    call self%get_parameter(&
         self%K_omno_no3, 'K_omno_no3', '[uM N]',&
         'half sat. of no3 for OM denitr.',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_omno_no2, 'K_omno_no2', '[uM N]',&
         'half sat. of no2 for OM denitr.',&
         default=0.001_rk)
    call self%get_parameter(&
         self%K_annamox, 'K_annamox', '[1/day]',&
         'Spec.rate of Anammox',&
         default=0.8_rk)
    !----O2--------!
    call self%get_parameter(&
         self%O2s_nf, 'O2s_nf', '[uM O]',&
         'half saturation for nitrification',&
         default=4.488_rk)
    call self%get_parameter(&
         self%O2s_dn, 'O2s_dn', '[uM O]',&
         'half saturation for denitrification',&
         default=10.0_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(&
         self%r_c_n,   'r_c_n',  '[-]',&
         'C[uM]/N[uM]',&
         default=8.0_rk)
    call self%get_parameter(&
         self%r_n_p,   'r_n_p',  '[-]',&
         'N[uM]/P[uM]',&
         default=16.0_rk)

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
         self%id_DIC,'DIC','mmol/m**3 C',&
         'total dissolved inorganic carbon',required=.false.)
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    call self%register_state_dependency(&
         self%id_O2, 'O2', 'mmol/m**3 O2',&
         'dissolved oxygen')
    call self%register_state_dependency(&
         self%id_PON,'PON','mmol/m**3 N',&
         'particulate organic nitrogen')
    call self%register_state_dependency(&
         self%id_DON,'DON','mmol/m**3 N',&
         'dissolved organic nitrogen')
    call self%register_state_dependency(&
         self%id_po4,'PO4','mmol/m**3 P',&
         'phosphate',required=.false.)

    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPM_NOX,'DcPM_NOX','mmol/m**3 N',&
         'PON denitrification (1+2 stage)',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDM_NOX,'DcDM_NOX','mmol/m**3 N',&
         'DON denitrification (1+2 stage)',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Nitrif1,'Nitrif1','mmol/m**3 N',&
         'Nitrification 1 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Nitrif2,'Nitrif2','mmol/m**3 N',&
         'Nitrification 2 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_anammox,'Anammox','mmol/m**3 N',&
         'Anammox',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Denitr1_PM,'Denitr1_PM','mmol/m**3 N',&
         'PON denitrification 1 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Denitr1_DM,'Denitr1_DM','mmol/m**3 N',&
         'DON denitrification 1 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Denitr2_PM,'Denitr2_PM','mmol/m**3 N',&
         'PON denitrification 2 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Denitr2_DM,'Denitr2_DM','mmol/m**3 N',&
         'DON denitrification 2 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Denitr1,'Denitr1','mmol/m**3 N',&
         '(PON+DON) denitrification 1 stage',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Denitr2,'Denitr2','mmol/m**3 N',&
         '(PON+DON) denitrification 2 stage',&
         output=output_time_step_integrated)

    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !check alkalinity
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_nitrogen),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: O2,PON,DON
    real(rk):: NO2,NO3,NH4
    !increments
    real(rk):: d_O2,d_DON,d_PON
    real(rk):: d_NO2,d_NO3,d_NH4
    real(rk):: d_Alk,d_DIC
    real(rk):: d_PO4
    !processes
    real(rk):: DcPM_NOX,DcDM_NOX
    real(rk):: Nitrif1,Nitrif2,Anammox
    real(rk):: Denitr1_PM,Denitr1_DM
    real(rk):: Denitr2_PM,Denitr2_DM
    real(rk):: Denitr1,Denitr2
    !summiriazed OM mineralization in N units
    real(rk):: Dc_OM_total

    _LOOP_BEGIN_
      !Retrieve current variable values
      !state
      _GET_(self%id_DON,DON)
      _GET_(self%id_NO2,NO2)
      _GET_(self%id_NO3,NO3)
      !solids
      _GET_(self%id_PON,PON)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_NH4,NH4)

      !N
      !Nitrification 1st stage: NH4+ + 1.5 O2 ->
      !                         NO2- + 2H+ + H2O (Canfield,2005)
      Nitrif1 = self%K_nitrif1*NH4*o2*0.5_rk*(1._rk+tanh(o2-self%O2s_nf))
      !Nitrification 2d stage: NO2- + 0.5 O2 -> NO3- (Canfield,2005)
      Nitrif2 = self%K_nitrif2*NO2*o2*0.5_rk*(1._rk+tanh(o2-self%O2s_nf))
      !in suboxic conditions
      !Anammox NO2- + NH4+ -> N2 + 2H2O (Canfield,2005)
      Anammox = self%K_annamox*NO2*NH4 &
               *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn)))
      !OM denitrification (Richards, 1965)
      !(CH2O)106(NH3)16H3PO4 + 84.8HNO3 =
      ! 106CO2 + 42.4N2 + 148.4H2O + 16NH3 + H3PO4
      !POM and DOM denitrification (1st stage) (Anderson,1982)
      !1/2CH2O + NO3- -> NO2- + 1/2H2O + 1/2CO2
      Denitr1_PM = self%K_denitr1*PON &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO3/(NO3+self%K_omno_no3)
      Denitr1_DM = self%K_denitr1*DON &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO3/(NO3+self%K_omno_no3)
      !POM and DOM denitrification (2d stage)
      !3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 (Anderson,1982)
      Denitr2_PM = self%K_denitr2*PON &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO2/(NO2+self%K_omno_no2)
      Denitr2_DM = self%K_denitr2*DON &
                  *(1._rk-0.5_rk*(1._rk+tanh(o2-self%O2s_dn))) &
                  *NO2/(NO2+self%K_omno_no2)
      Denitr1 = Denitr1_PM + Denitr1_DM
      Denitr2 = Denitr2_PM + Denitr2_DM
      !from (Anderson, 1982) and (Richards, 1965)
      !Denitrification POM and DOM (1+2 stage)
      DcPM_NOX = 16._rk/212._rk*Denitr1_PM+16._rk/141.3_rk*Denitr2_PM
      DcDM_NOX = 16._rk/212._rk*Denitr1_DM+16._rk/141.3_rk*Denitr2_DM

      !Summariazed OM mineralization
      Dc_OM_total = DcPM_NOX+DcDM_NOX

      !Set increments
      !O2
      d_O2 = -1.5_rk*Nitrif1-0.5_rk*Nitrif2
      _SET_ODE_(self%id_O2,d_O2)
      !organic matter
      d_DON = -DcDM_NOX
      _SET_ODE_(self%id_DON,d_DON)
      !solid OM
      d_PON = -DcPM_NOX
      _SET_ODE_(self%id_PON,d_PON)
      !N
      d_NO2 = Nitrif1-Nitrif2+Denitr1-Denitr2-Anammox
      _SET_ODE_(self%id_NO2,d_NO2)
      d_NO3 = Nitrif2-Denitr1
      _SET_ODE_(self%id_NO3,d_NO3)
      !N gas
      d_NH4 = Dc_OM_total-Nitrif1-Anammox
      _SET_ODE_(self%id_NH4,d_NH4)
      !DIC
      d_DIC = Dc_OM_total*self%r_c_n
      _SET_ODE_(self%id_DIC,d_DIC)
      !P
      d_PO4 = Dc_OM_total/self%r_n_p
      _SET_ODE_(self%id_PO4,d_PO4)
      !Alkalinity changes due to redox reactions:
      d_Alk = (&
             !NH4+ + 1.5 O2 -> NO2- + 2H+ + H2O
             -2._rk*Nitrif1 &  !(Wolf-Gladrow, Zeebe, 2007)
             !3/4CH2O + H+ + NO2- -> 1/2N2 + 5/4H2O + 3/4CO2 or
             !5 CH2O + 4 H+ + 4 NO3- -> 2 N2 + 5 CO2 + 7H2O
             +1._rk*(Denitr2_PM+Denitr2_DM) &
             + Dc_OM_total &
             )
      _SET_ODE_(self%id_Alk,d_Alk)

      _SET_DIAGNOSTIC_(self%id_DcPM_NOX,DcPM_NOX)
      _SET_DIAGNOSTIC_(self%id_DcDM_NOX,DcDM_NOX)
      _SET_DIAGNOSTIC_(self%id_anammox,Anammox)
      _SET_DIAGNOSTIC_(self%id_Denitr1_PM,Denitr1_PM)
      _SET_DIAGNOSTIC_(self%id_Denitr1_DM,Denitr1_DM)
      _SET_DIAGNOSTIC_(self%id_Denitr2_PM,Denitr2_PM)
      _SET_DIAGNOSTIC_(self%id_Denitr2_DM,Denitr2_DM)
      _SET_DIAGNOSTIC_(self%id_Denitr1,Denitr1)
      _SET_DIAGNOSTIC_(self%id_Denitr2,Denitr2)
      _SET_DIAGNOSTIC_(self%id_Nitrif1,Nitrif1)
      _SET_DIAGNOSTIC_(self%id_Nitrif2,Nitrif2)
    _LOOP_END_
  end subroutine do
end module fabm_niva_brom_nitrogen
