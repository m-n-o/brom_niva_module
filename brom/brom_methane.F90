!-----------------------------------------------------------------------
! brom_methane is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_methane
  use fabm_types
  use brom_functions

  implicit none
  private
  type,extends(type_base_model),public:: type_niva_brom_methane
    !all descriptions are in the initialize subroutine
    type(type_state_variable_id):: id_CH4
    !state dependencies
    type(type_state_variable_id):: id_DIC,id_NH4,id_PO4,id_Si
    type(type_state_variable_id):: id_DOML,id_POML,id_POMR,id_DOMR
    type(type_state_variable_id):: id_O2,id_NO3,id_SO4
    !diagnostic variables by bacteria needed
    type(type_diagnostic_variable_id):: id_DcPOML_ch4,id_DcDOML_ch4
    type(type_diagnostic_variable_id):: id_DcPOMR_ch4,id_DcDOMR_ch4
    !for do_surface
    type(type_dependency_id):: id_temp,id_salt
    type(type_horizontal_dependency_id):: id_windspeed
    !Model parameters
    !specific rates of biogeochemical processes
    real(rk):: s_omso_o2,s_omso_no3,s_omch_so4,s_OM_refr
    real(rk):: K_DOML_ch4,K_POML_ch4,K_POMR_ch4,K_DOMR_ch4
    real(rk):: K_ch4_o2,K_ch4_so4
    !---- Stoichiometric coefficients ----!
    real(rk):: c_to_n, c_to_si, c_to_p

    contains
    procedure :: initialize
    procedure :: do_surface
    procedure :: do
  end type

  contains

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
        self%s_OM_refr, 's_OM_refr', '[uM N]',&
        'threshold of decay of refractory OM',&
        default=5.0_rk)
    call self%get_parameter(&
        self%K_DOML_ch4, 'K_DOML_ch4', '[1/day]',&
        'Specific rate of methane production from DOML',&
        default=0.00014_rk)
    call self%get_parameter(&
        self%K_POML_ch4, 'K_POML_ch4', '[1/day]',&
        'Specific rate of methane production from POML',&
        default=0.00014_rk)
    call self%get_parameter(&
        self%K_POMR_ch4, 'K_POMR_ch4', '[1/day]',&
        'Specific rate of methane production from POMR',&
        default=0.00014_rk)
    call self%get_parameter(&
        self%K_DOMR_ch4, 'K_DOMR_ch4', '[1/day]',&
        'Specific rate of methane production from DOMR',&
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
    call self%get_parameter(self%c_to_n,'c_to_n','[-]','C[uM]/N[uM]',&
                            default=106._rk/16._rk)
    call self%get_parameter(self%c_to_si,'c_to_si','[-]','C[uM]/Si[uM]',&
                            default=106._rk/15._rk)
    call self%get_parameter(self%c_to_p,'c_to_p','[-]','C[uM]/P[uM]',&
                            default=106._rk)
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
        self%id_PO4,'PO4','mmol/m**3',&
        'phosphate')
    call self%register_state_dependency(&
        self%id_Si,'Si','mmol/m**3',&
        'silicate')
    call self%register_state_dependency(&
        self%id_NH4,'NH4','mmol/m**3',&
        'ammonium')
    call self%register_state_dependency(&
        self%id_POML,'POML','mg m^3',&
        'particulate organic nitrogen')
    call self%register_state_dependency(&
        self%id_POMR,'POMR','mg m^3',&
        'particulate organic nitrogen')
    call self%register_state_dependency(&
        self%id_DOMR,'DOMR','mg m^3',&
        'DOMR')
    call self%register_state_dependency(&
        self%id_DOML,'DOML','mg m^3',&
        'dissolved organic nitrogen')
    call self%register_state_dependency(&
        self%id_SO4,'SO4','mmol/m**3','sulphate')

    !for do_surface
    call self%register_dependency(&
        self%id_temp,standard_variables%temperature)
    call self%register_dependency(&
        self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(&
        self%id_windspeed,standard_variables%wind_speed)
    !Register diagnostic variables
    call self%register_diagnostic_variable(&
        self%id_DcPOML_ch4,'DcPOML_ch4','mmol/m**3',&
        'CH4 production from POML',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
        self%id_DcDOML_ch4,'DcDOML_ch4','mmol/m**3',&
        'CH4 production from DOML',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
        self%id_DcPOMR_ch4,'DcPOMR_ch4','mmol/m**3',&
        'CH4 production from POMR ',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
        self%id_DcDOMR_ch4,'DcDOMR_ch4','mmol/m**3',&
        'CH4 production from DOMR ',output=output_time_step_integrated)

    !Specify that rates are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize

  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    class (type_niva_brom_methane),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_SURFACE_
    real(rk):: Q_CH4, Q_pCH4, CH4
    real(rk):: temp, salt, abs_temp
    real(rk):: Sc,k_660,k_CH4_660
    real(rk):: pCH4a, pCH4w
    real(rk):: windspeed
    real(rk):: a1,a2,a3,b1,b2,b3,bunsen,s

    _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_temp,temp) !temperature
      _GET_(self%id_salt,salt) !salinity
      _GET_(self%id_CH4, CH4) !previous CH4 which calculates here in the subroutine do
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)

      abs_temp = temp + 273._rk ![k]
      a1 = -67.1962_rk  ! #-68.8862
      a2 = 99.1624_rk  ! #101.4956
      a3 = 27.9015_rk ! #28.7314
      b1 = -0.072909_rk ! #0.076146
      b2 = 0.041674_rk ! #0.043970
      b3 = -0.0064603_rk ! #-0.0068672

      bunsen = 2.718281828459_rk  **             &
      (a1 + a2*(100._rk / abs_temp) + a3 * log(abs_temp / 100._rk) + &
      (salt * ( b1 + b2 * (abs_temp / 100._rk)  +  &
      b3 * (abs_temp / 100._rk **2._rk))))
      ! solubility
      s = bunsen / 22.4  !# mole /l
      pCH4a = 1.8E-6_rk
      pCH4w = CH4 * 0.082057 * abs_temp ![uatm] if ch4 in uM

      !calculate the scmidt number and unit conversions
      !Sc = 2073.1_rk-125.62_rk*temp+3.6276_rk*temp**2._rk-0.043219_rk*&
      !     temp**3.0_rk
      !Sc corrected for CH4
      Sc = 2101.2_rk - 131.54_rk * temp + &
          4.4931_rk * temp **2_rk - 0.08676_rk * temp **3 + &
          0.00070663_rk * temp ** 4_rk

      k_CH4_660 = (24._rk/100._rk)*(0.24_rk * windspeed**2._rk)*&
      (Sc/ 660._rk)**(-0.5_rk) !m/d
      Q_pCH4 = k_CH4_660 * (pCH4a - max(0e0,pCH4w))
      Q_CH4 = Q_pCH4 * s

      _SET_SURFACE_EXCHANGE_(self%id_CH4,Q_CH4)
    _HORIZONTAL_LOOP_END_
  end subroutine do_surface

  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_methane),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    !state variables
    real(rk):: DIC,SO4,NO3,NH4
    real(rk):: O2,CH4
    real(rk):: POML,POMR,DOML,DOMR
    !processes
    real(rk):: DcDOML_ch4,DcPOML_ch4,DcPOMR_ch4,DcDOMR_ch4
    real(rk):: DcTOM_CH4, ch4_o2,ch4_so4
    real(rk):: dpoml_ch4_in_m, ddoml_ch4_in_m
    real(rk):: dpomr_ch4_in_m, ddomr_ch4_in_m
    !increments
    real(rk):: d_SO4,d_O2,d_CH4,d_NH4,d_DIC,d_PO4,d_Si
    real(rk):: d_DOML,d_POML,d_POMR,d_DOMR
    real(rk):: thr_o2_l,thr_no3,thr_so4

    _LOOP_BEGIN_
      !state variables
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_DOMR,DOMR)
      _GET_(self%id_SO4,SO4)
      _GET_(self%id_NO3,NO3)
      _GET_(self%id_NH4,NH4)
      !gases
      _GET_(self%id_O2,O2)
      _GET_(self%id_CH4,CH4)
      !solids
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)

      !Correct for methane max solubility value
      if (CH4 > 1340._rk) then
          CH4 = 1340._rk
      end if

      !CH4 production from POML and DOML
      !(CH2O)106(NH3)16H3PO4 -> 53 CO2 + 53 CH4 + 16 NH3 + H3PO4
      thr_o2_l = hyper_inhibitor(self%s_omso_o2,o2,1._rk)
      thr_no3  = hyper_inhibitor(self%s_omso_no3,no3,1._rk)
      thr_so4  = hyper_inhibitor(self%s_omch_so4,so4,1._rk)

      !OM [mg C m^-3]
      DcDOML_ch4 = thr_o2_l*thr_no3*thr_so4*self%K_DOML_ch4*DOML
      DcPOML_ch4 = thr_o2_l*thr_no3*thr_so4*self%K_POML_ch4*POML
      DcDOMR_ch4 = thr_o2_l*thr_no3*thr_so4 &
                  *hyper_limiter(self%s_OM_refr,DOMR,0.1_rk) &
                  *self%K_DOMR_ch4*DOMR
      DcPOMR_ch4 = thr_o2_l*thr_no3*thr_so4 &
                  *hyper_limiter(self%s_OM_refr,POMR,0.1_rk) &
                  *self%K_POMR_ch4*POMR
      !recalculate to moles
      dpoml_ch4_in_m = carbon_g_to_mole(DcPOML_ch4)
      ddoml_ch4_in_m = carbon_g_to_mole(DcDOML_ch4)
      dpomr_ch4_in_m = carbon_g_to_mole(DcPOMR_ch4)
      ddomr_ch4_in_m = carbon_g_to_mole(DcDOMR_ch4)

      !CH4 oxidation with O2
      !CH4 + 2O2 = CO2 + 2H2O
      ch4_o2 = self%K_ch4_o2*CH4*hyper_limiter(self%s_omso_o2,o2,1._rk)
      !CH4 + SO42- + 2 H+  = CO2 + H2S + 2H2O
      ch4_so4 = self%K_ch4_so4*CH4*SO4*thr_o2_l

      !Calculate increments
      d_DOML = -DcDOML_ch4
      d_POML = -DcPOML_ch4
      d_DOMR = DcDOML_ch4-DcDOMR_ch4
      d_POMR = DcPOML_ch4-DcPOMR_ch4
      d_O2  = -2._rk*ch4_o2
      d_SO4 = -ch4_so4
      d_DIC = 0.5_rk*(ddomr_ch4_in_m+dpomr_ch4_in_m) &
             +ch4_o2+ch4_so4
      d_CH4 = 0.5_rk*(ddomr_ch4_in_m+dpomr_ch4_in_m) &
             -ch4_o2-ch4_so4
      d_NH4 = (dpoml_ch4_in_m+ddoml_ch4_in_m)/self%c_to_n
      _SET_ODE_(self%id_NH4,d_NH4)
      d_PO4 = (dpoml_ch4_in_m+ddoml_ch4_in_m)/self%c_to_p
      _SET_ODE_(self%id_PO4,d_PO4)
      d_Si = (dpoml_ch4_in_m+ddoml_ch4_in_m)/self%c_to_si
      _SET_ODE_(self%id_Si,d_Si)

      !Set increments
      _SET_ODE_(self%id_DOML,d_DOML)
      _SET_ODE_(self%id_POML,d_POML)
      _SET_ODE_(self%id_DOMR,d_DOMR)
      _SET_ODE_(self%id_POMR,d_POMR)
      _SET_ODE_(self%id_O2,d_O2)
      _SET_ODE_(self%id_SO4,d_SO4)
      _SET_ODE_(self%id_DIC,d_DIC)
      _SET_ODE_(self%id_CH4,d_CH4)
      _SET_ODE_(self%id_NH4,d_NH4)
      _SET_ODE_(self%id_PO4,d_PO4)
      _SET_ODE_(self%id_Si,d_Si)

      _SET_DIAGNOSTIC_(self%id_DcPOML_ch4,DcPOML_ch4)
      _SET_DIAGNOSTIC_(self%id_DcDOML_ch4,DcDOML_ch4)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_ch4,DcPOMR_ch4)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_ch4,DcDOMR_ch4)
    _LOOP_END_
  end subroutine do
end module fabm_niva_brom_methane
