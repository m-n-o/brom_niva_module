!-----------------------------------------------------------------------
! brom_bio is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

#include "fabm_driver.h"

module fabm_niva_brom_bio
  use fabm_types
  use brom_functions

  implicit none
  private
  type,extends(type_base_model),public :: type_niva_brom_bio
    !variables allocated here
    type(type_state_variable_id):: id_Phy,id_Het
    type(type_state_variable_id):: id_O2,id_POML,id_POMR,id_DOML,id_DOMR
    !state dependencies
    type(type_state_variable_id):: id_NO2,id_NO3,id_NH4,id_PO4
    !type(type_state_variable_id):: id_Baae,id_Baan,id_Bhae,id_Bhan
    type(type_state_variable_id):: id_DIC,id_Si,id_Alk
    !type(type_state_variable_id):: id_Sipart, id_H2S
    !standard variables dependencies
    type(type_dependency_id):: id_temp,id_salt,id_par
    type(type_horizontal_dependency_id):: id_windspeed
    type(type_horizontal_dependency_id):: id_lat
    type(type_global_dependency_id):: id_day
    !diagnostic variables
    !organic matter
    type(type_diagnostic_variable_id):: id_DcPOML_O2,id_DcDOML_O2
    type(type_diagnostic_variable_id):: id_DcPOMR_O2,id_DcDOMR_O2!,id_DcTOM_O2
    !type(type_diagnostic_variable_id):: id_POMTot,id_DOMTot
    !primary producers
    type(type_diagnostic_variable_id):: id_LimNH4,id_LimNO3,id_LimN
    type(type_diagnostic_variable_id):: id_LimP,id_LimSi
    type(type_diagnostic_variable_id):: id_Biorate,id_growthrate
    type(type_diagnostic_variable_id):: id_GrowthPhy,id_MortPhy,id_ExcrPhy
    type(type_diagnostic_variable_id):: id_ChlCratio,id_N_fixation
    !heterotrophs
    type(type_diagnostic_variable_id):: id_GrazPhy,id_GrazPOP!,id_GrazBact
    !type(type_diagnostic_variable_id):: id_GrazBaae,id_GrazBaan
    !type(type_diagnostic_variable_id):: id_GrazBhae,id_GrazBhan
    type(type_diagnostic_variable_id):: id_Grazing
    type(type_diagnostic_variable_id):: id_RespHet,id_MortHet
    !oxygen
    type(type_diagnostic_variable_id):: id_O2_rel_sat,id_O2_sat,id_AOU
    !Model parameters
    !specific rates of biogeochemical processes
    real(rk):: K_POML_DOML,K_POMR_DOMR !for autolysis
    real(rk):: K_DOML_ox,K_DOMR_ox,K_POML_ox,K_POMR_ox !specific decay scales
    real(rk):: tref,K_omox_o2 !for OM decay
    !---- N, P, Si--!
    real(rk):: K_nh4_lim,K_nox_lim,K_nfix,K_po4_lim,K_si_lim
    !----Phy  ----------!
    real(rk):: pbm,alpha
    real(rk):: K_phy_mrt,K_phy_exc
    !----Het -----------!
    real(rk):: K_het_phy_gro,K_het_phy_lim
    real(rk):: K_het_pom_gro,K_het_pom_lim
    !real(rk):: K_het_bac_gro,limGrazBac
    real(rk):: K_het_res,K_het_mrt,Uz,Hz
    !---- Sinking---!
    real(rk):: Wsed,Wphy,Whet
    !---- Stoichiometric coefficients ----!
    real(rk):: c_to_n, c_to_si, c_to_p
  contains
    procedure :: initialize
    procedure :: do
    procedure :: do_surface
    !procedure :: graz
  end type
contains
  !
  !
  !
  subroutine initialize(self,configunit)
    class (type_niva_brom_bio), intent(inout), target :: self
    integer,                    intent(in)            :: configunit

    !----Phy----------!
    call self%get_parameter(&
         self%alpha,'alpha','mg C (mg Chl a h)^-1 (microM quanta m^-2 s^-1)^-1',&
         'Photosynthetic efficiency at low irradiance',&
         default=.03_rk)
    call self%get_parameter(&
         self%pbm,'pbm','mg C (mg Chl a h)^-1','Maximum rate of photosynthesis',&
         default=8.0_rk)
    call self%get_parameter(&
         self%K_po4_lim,'K_po4_lim','mM P m^-3',&
         'Half-sat. constant for uptake of PO4 by Phy',&
         default=0.01_rk)
    call self%get_parameter(&
         self%K_si_lim,'K_si_lim','mM Si m^-3',&
         'Half-sat. constant for uptake of Si by Phy',&
         default=1._rk)
    call self%get_parameter(&
         self%K_nox_lim,'K_nox_lim','mM N m^-3',&
         'Half-sat. constant for uptake of NO3+NO2',&
         default=1._rk)
    call self%get_parameter(&
         self%K_nh4_lim,'K_nh4_lim','mM N m^-3',&
         'Half-sat. constant for uptake of NH4',&
         default=0.5_rk)
    call self%get_parameter(&
         self%K_nfix,'K_nfix','1 d^-1',&
         'Max. specific rate of nitrogen fixation',&
         default=0.4_rk)
    call self%get_parameter(&
         self%K_phy_mrt,'K_phy_mrt','1 d^-1','Specific rate of mortality',&
         default=0.0002_rk)
    call self%get_parameter(&
         self%K_phy_exc,'K_phy_exc','1 d^-1','Specific rate of excretion',&
         default=0.015_rk)
    !----Het----------!
    call self%get_parameter(&
         self%K_het_phy_gro,'K_het_phy_gro','1 d^-1',&
         'Maximum rate of grazing of Het on Phy',&
         default=0.4_rk)
    call self%get_parameter(&
         self%K_het_phy_lim,'K_het_phy_lim','dimensionless',&
         'Half-sat. constant for grazing of Het on Phy for Phy/Het ratio',&
         default=2._rk)
    call self%get_parameter(&
         self%K_het_pom_gro,'K_het_pom_gro','1 d^-1',&
         'Maximum rate of grazing of Het on POML (eats labile particulate organic))',&
         default=0.4_rk)
    call self%get_parameter(&
         self%K_het_pom_lim,'K_het_pom_lim','dimensionless',&
         'Half-sat. constant for grazing of Het on POM for POM/Het ratio',&
         default=2._rk)
    !call self%get_parameter(&
    !     self%K_het_bac_gro,'K_het_bac_gro','mmol/m**3',&
    !     'Max.spec.rate of grazing of Het on bacteria',&
    !     default=0.70_rk)
    !call self%get_parameter(&
    !     self%limGrazBac,'limGrazBac','mmol/m**3',&
    !     'Limiting parameter for bacteria grazing by Het',&
    !     default=2._rk)
    call self%get_parameter(&
         self%K_het_res,'K_het_res','1 d^-1',&
         'Specific respiration rate',&
         default=0.015_rk)
    call self%get_parameter(&
         self%K_het_mrt,'K_het_mrt','1 d^-1',&
         'Maximum specific rate of mortality of Het',&
         default=0.1_rk)
    call self%get_parameter(&
         self%Uz,'Uz','dimensionless',&
         'Food absorbency for Het',&
         default=0.5_rk)
    call self%get_parameter(&
         self%Hz,'Hz','dimensionless',&
         'Ratio between dissolved and particulate excretes of Het',&
         default=0.5_rk)
    !----OM-----------!
    call self%get_parameter(&
         self%K_POML_DOML, 'K_POML_DOML', '1 d^-1',&
         'Specific rate of Autolysis of POML to DOML',&
         default=0.15_rk)
    call self%get_parameter(&
         self%K_POMR_DOMR, 'K_POMR_DOMR', '1 d^-1',&
         'Specific rate of Autolysis of POMR to DOMR',&
         default=0.00001_rk)
    call self%get_parameter(&
         self%K_DOML_ox,'K_DOML_ox','1 d^-1',&
         'Specific rate of oxidation of DOML with O2',&
         default=0.1_rk)
    call self%get_parameter(&
         self%K_DOMR_ox,'K_DOMR_ox','1 d^-1',&
         'Specific rate of oxidation of DOMR with O2',&
         default=0.1_rk)
    call self%get_parameter(&
         self%K_POML_ox,'K_POML_ox','1 d^-1',&
         'Specific rate of oxidation of POML with O2',&
         default=0.05_rk)
    call self%get_parameter(&
         self%K_POMR_ox,'K_POMR_ox','1 d^-1',&
         'Specific rate of oxidation of POMR with O2',&
         default=0.05_rk)
    call self%get_parameter(&
         self%tref,'tref','degrees Celsius',&
         'Reference temperature at which temperature factor = 1',&
         default=0.0_rk)
    call self%get_parameter(&
         self%K_omox_o2,'K_omox_o2','mM O2 m^-3',&
         'Half-sat. constant of O2 for OM mineralization',&
         default=1.0_rk)
    !----Sinking--------
    call self%get_parameter(&
         self%Wsed,'Wsed','m d^-1',&
         'Rate of sinking of detritus',&
         default=5.00_rk)
    call self%get_parameter(&
         self%Wphy,'Wphy','m d^-1',&
         'Rate of sinking of Phy',&
         default=0.10_rk)
    call self%get_parameter(&
         self%Whet,'Whet','m d^-1',&
         'Rate of sinking of Het',&
         default=1.00_rk)
    !----Stoichiometric coefficients----!
    call self%get_parameter(self%c_to_n,'c_to_n','[-]','C[uM]/N[uM]',&
                            default=106._rk/16._rk)
    call self%get_parameter(self%c_to_si,'c_to_si','[-]','C[uM]/Si[uM]',&
                            default=106._rk/15._rk)
    call self%get_parameter(self%c_to_p,'c_to_p','[-]','C[uM]/P[uM]',&
                            default=106._rk)
    !Register state variables
    call self%register_state_variable(&
         self%id_Phy,'Phy','mg C m^-3','Phy',&
         minimum=0._rk,initial_value=1._rk,&
         vertical_movement=-self%Wphy/86400._rk)
    call self%register_state_variable(&
         self%id_Het,'Het','mg C m^-3','Het',&
         minimum=0._rk,initial_value=1._rk,&
         vertical_movement=-self%Whet/86400._rk)
    call self%register_state_variable(&
         self%id_POML,'POML','mg C m^-3','POML labile',&
         minimum=0._rk,initial_value=0._rk,&
         vertical_movement=-self%Wsed/86400._rk)
    call self%register_state_variable(&
         self%id_POMR,'POMR','mg C m^-3','POML refractory',&
         minimum=0.0_rk,initial_value=0._rk,&
         vertical_movement=-self%Wsed/86400._rk)
    call self%register_state_variable(&
         self%id_DOML,'DOML','mg C m^-3','DOML labile',&
         minimum=0.0_rk,initial_value=0._rk)
    call self%register_state_variable(&
         self%id_DOMR,'DOMR','mg C m^-3','DOMR refractory',&
         minimum=0.0_rk,initial_value=0._rk)
    call self%register_state_variable(&
         self%id_O2,'O2','mM O2 m^-3','O2',&
         minimum=0.0_rk,initial_value=350._rk)
    !Register state dependencies
    call self%register_state_dependency(&
         self%id_PO4,'PO4','mM P m^-3','PO4')
    call self%register_state_dependency(&
         self%id_Si,'Si','mM Si m^-3','Si')
    !call self%register_state_dependency(&
    !     self%id_Sipart,'Sipart','mmol/m**3','Si particulate')
    call self%register_state_dependency(&
         self%id_NO3,'NO3','mM N m^-3','NO3')
    call self%register_state_dependency(&
         self%id_NH4,'NH4','mM N m^-3','NH4')
    call self%register_state_dependency(&
         self%id_NO2,'NO2','mM N m^-3','NO2')
    call self%register_state_dependency(&
         self%id_DIC,'DIC','mM C m^-3','DIC')
    !call self%register_state_dependency(&
    !     self%id_H2S,'H2S','mmol/m**3','H2S')
    call self%register_state_dependency(self%id_Alk,&
         standard_variables%alkalinity_expressed_as_mole_equivalent)
    !call self%register_state_dependency(&
    !     self%id_Baae,'Baae','mmol/m**3','aerobic autotrophic bacteria')
    !call self%register_state_dependency(&
    !     self%id_Bhae,'Bhae','mmol/m**3','aerobic heterotrophic bacteria')
    !call self%register_state_dependency(&
    !     self%id_Baan,'Baan','mmol/m**3','anaerobic aurotrophic bacteria')
    !call self%register_state_dependency(&
    !     self%id_Bhan,'Bhan','mmol/m**3','anaerobic heterotrophic bacteria')
    !Register diagnostic variables
    call self%register_diagnostic_variable(&
         self%id_DcPOML_O2,'DcPOML_O2','mM N m^-3',&
         'POML with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcPOMR_O2,'DcPOMR_O2','mM N m^-3',&
         'POMR with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOMR_O2,'DcDOMR_O2','mM N m^-3',&
         'DOMR with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_DcDOML_O2,'DcDOML_O2','mM N m^-3',&
         'DOM with O2 mineralization',output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortHet,'MortHet','mM N m^-3','Mortality of Het',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_RespHet,'RespHet','mM N m^-3','Respiration rate of Het',&
         output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_GrazBhae,'GrazBhae','mmol/m**3','GrazBhae',&
    !     output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_GrazBhan,'GrazBhan','mmol/m**3','GrazBhan',&
    !     output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_GrazBaae,'GrazBaae','mmol/m**3','GrazBaae',&
    !     output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_GrazBaan,'GrazBaan','mmol/m**3','GrazBaan',&
    !     output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazPhy,'GrazPhy','mM N m^-3','Het grazing on Phy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrazPOP,'GrazPOP','mM N m^-3','Het grazing on particulate labile OM',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Grazing,'Grazing','mM N m^-3','Het total grazing',&
         output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_GrazBact,'GrazBact','mmol/m**3','GrazBact',&
    !     output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MortPhy,'MortPhy','mM N m^-3','Mortality of Phy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ExcrPhy,'ExcrPhy','mM N m^-3','Excretion of Phy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimNH4,'LimNH4','dimensionless','Limitation of Phy on NH4',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimNO3,'LimNO3','dimensionless','Limitation of Phy on NO3',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimN,'LimN','dimensionless','Limitation of Phy on total N',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimP,'LimP','dimensionless','Limitation of Phy on P',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_LimSi,'LimSi','dimensionless','Limitation of Phy on Si',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Biorate,'Biorate','mg C(mg Chl a d)^-1','Biomass specific photosynthetic rate',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_growthrate,'Growth_rate','1 d^-1','Daily Phy growth rate',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_GrowthPhy,'GrowthPhy','mM N m^-3','Daily growth of Phy',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_ChlCratio,'ChlCratio','mg Chl a (mg C)^-1','Chl to C ratio',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_N_fixation,'N_fixation','mM N m^-3','Daily N2 fixation',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(self%id_O2_rel_sat, &
         'O2_rel_sat','dimensionless','Relative oxygen saturation', &
         standard_variable=standard_variables%fractional_saturation_of_oxygen)
    call self%register_diagnostic_variable(self%id_O2_sat, &
         'O2_sat','mM O2 m^-3','Oxygen saturation concentration')
    call self%register_diagnostic_variable(self%id_AOU, &
         'AOU','mM O2 m^-3','Apparent oxigen utilization')
    !call self%register_diagnostic_variable(&
    !     self%id_DcTOM_O2,'DcTOM_O2','mmol/m**3',&
    !     'Refractory OM oxidation with O2',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_DOMTot,'DOMTot','mmol/m**3',&
    !     'DOMTot: refractory+labile',output=output_time_step_integrated)
    !call self%register_diagnostic_variable(&
    !     self%id_POMTot,'POMTot','mmol/m**3',&
    !     'POMTot: refractory+labile',output=output_time_step_integrated)
    !Register environmental dependencies
    call self%register_dependency(&
         self%id_par,&
         standard_variables%downwelling_photosynthetic_radiative_flux)
    call self%register_dependency(&
         self%id_temp,standard_variables%temperature)
    call self%register_dependency(&
         self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(&
         self%id_windspeed,standard_variables%wind_speed)
    call self%register_horizontal_dependency(self%id_lat,&
         standard_variables%latitude)
    call self%register_global_dependency(self%id_day,&
         standard_variables%number_of_days_since_start_of_the_year)
    !Specify that are rates computed in this module are per day
    !(default: per second)
    self%dt = 86400._rk
  end subroutine initialize
  !
  !
  subroutine do(self,_ARGUMENTS_DO_)
    class (type_niva_brom_bio),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_
    real(rk):: temp,salt,PAR,PAR_M_day,latitude,day
    !nutrients and limiters
    real(rk):: NH4,NO2,NO3,PO4,O2,Si!,H2S,Sipart
    real(rk):: LimT,LimP,LimNO3,LimNH4,LimN,LimSi,LimNut
    !biology
    !real(rk):: Baae,Baan,Bhae,Bhan
    real(rk):: Phy,ChlC,biorate,growthrate
    real(rk):: GrowthPhy,ExcrPhy,MortPhy,N_fixation
    real(rk):: Het,GrazPhy,GrazPOP
    !real(rk):: GrazBaae,GrazBaan,GrazBhae,GrazBhan,GrazBact
    real(rk):: Grazing,RespHet,MortHet
    !OM
    real(rk):: POML,POMR,DOML,DOMR,kf
    real(rk):: DcDOML_O2,DcPOML_O2!,DcTOM_O2
    real(rk):: DcPOMR_O2,DcDOMR_O2!,DOMTot,POMTot
    real(rk):: Autolysis_L,Autolysis_R
    !etc.
    real(rk):: O2_sat,Alk
    !increments
    real(rk):: d_NO2,d_NO3,d_PO4,d_Si,d_DIC,d_O2,d_NH4!,d_Sipart
    real(rk):: d_Phy,d_Het!,d_Baae,d_Baan,d_Bhae,d_Bhan
    real(rk):: d_alk,d_DOML,d_DOMR,d_POML,d_POMR
    !some organic variables in moles
    real(rk):: dphy_in_m,dzoo_resp_in_m
    real(rk):: ddoml_o2_in_m,dpoml_o2_in_m
    real(rk):: ddomr_o2_in_m,dpomr_o2_in_m

    ! Enter spatial loops (if any)
    _LOOP_BEGIN_
      ! Retrieve current environmental conditions.
      _GET_(self%id_par,PAR) ! local photosynthetically active radiation microM/s
      _GET_(self%id_temp,temp) ! temperature
      _GET_(self%id_salt,salt) ! salinity
      _GET_GLOBAL_(self%id_day,day)
      _GET_HORIZONTAL_(self%id_lat,latitude)
      ! Retrieve current (local) state variable values.
      !state variables
      _GET_(self%id_NO2,NO2)
      _GET_(self%id_NO3,NO3)
      _GET_(self%id_PO4,PO4)
      _GET_(self%id_Si,Si)
      _GET_(self%id_Alk,Alk)
      _GET_(self%id_NH4,NH4)
      !_GET_(self%id_H2S,H2S)
      _GET_(self%id_O2,O2)
      _GET_(self%id_Phy,Phy)
      _GET_(self%id_Het,Het)
      !_GET_(self%id_Baae,Baae)
      !_GET_(self%id_Baan,Baan)
      !_GET_(self%id_Bhae,Bhae)
      !_GET_(self%id_Bhan,Bhan)
      _GET_(self%id_POML,POML)
      _GET_(self%id_POMR,POMR)
      _GET_(self%id_DOML,DOML)
      _GET_(self%id_DOMR,DOMR)
      !_GET_(self%id_Sipart,Sipart)

      PAR_M_day = PAR*86400._rk/1000000._rk !PAR M per day
      !
      !Phy
      !dependence on NH4
      LimNH4 = monod_squared(self%K_nh4_lim, NH4)
      !dependence on NO3+NO2
      LimNO3 = monod_squared(self%K_nox_lim, NO3+NO2)&
             * monod_squared_inhibitor(self%K_nh4_lim, NH4)
      !dependence on total nitrogen
      LimN = LimNO3+LimNH4
      !dependence on PO4
      LimP = monod_squared(self%K_po4_lim, PO4)
      !dependence on Si
      LimSi = monod_squared(self%K_si_lim, Si)
      !total limiter
      LimNut = min(LimN, LimP, LimSi)
      !Chl a / Carbon ratio
      ChlC = chl_c_ratio(temp, PAR_M_day, LimNut)
      !photosynthetic rate
      biorate = photosynthetic_rate(photoperiod(latitude, day),&
                                    self%pbm, self%alpha, PAR)
      !daily growth rate
      growthrate = daily_growth(biorate, ChlC)
      GrowthPhy = growthrate*Phy
      !1 mg C m^-3, overwintering value
      if (Phy < 1.01_rk) then
        ExcrPhy = 0._rk
        MortPhy = 0._rk
      else
        !excretion of phy
        ExcrPhy = self%K_phy_exc*Phy
        !mortality of phy
        MortPhy = Phy*Phy*self%K_phy_mrt
      end if
      !
      !Het
      !Grazing of Het on phy
      if (Phy < 1.01_rk) then
        GrazPhy = 0._rk
      else
        GrazPhy = self%K_het_phy_gro*Het*&
                  monod_squared(self%K_het_phy_lim, quota(Phy, Het))
      end if
      !Grazing of Het on detritus
      GrazPOP = self%K_het_pom_gro*Het*&
                monod_squared(self%K_het_pom_lim, quota(POML, Het))
      !Grazing of Het on bacteria
      !GrazBaae = 1.0_rk*self%graz(Baae, Het)
      !GrazBaan = 0.5_rk*self%graz(Baan, Het)
      !GrazBhae = 1.0_rk*self%graz(Bhae, Het)
      !GrazBhan = 1.3_rk*self%graz(Bhan, Het)
      !GrazBact = GrazBaae+GrazBaan+GrazBhae+GrazBhan
      !Total grazing of Het
      Grazing = GrazPhy+GrazPOP!+GrazBact
      !Respiration of Het
      RespHet = self%K_het_res*Het!*hyper_limiter(20._rk, O2, 1._rk)
      !Mortality of Het
      MortHet = Het*self%K_het_mrt
                     !+ hyper_inhibitor(20._rk, O2, 1._rk)*0.75_rk)!*0.3_rk&
                     !+ hyper_limiter(10._rk, H2S, 1._rk)*0.45_rk)
      !
      !Nitrogen fixation described as appearence of NH4 available for
      !phytoplankton: N2 -> NH4 :
      N_fixation = 0._rk!self%K_nfix*LimP*&
      !1._rk/(1._rk+((NO3+NO2+NH4)/n_zero(PO4)*16._rk)**4._rk)*GrowthPhy
      !
      !POML and DOML (Savchuk, Wulff,1996)
      Autolysis_L = self%K_POML_DOML*POML
      Autolysis_R = self%K_POMR_DOMR*POMR
      !OM decay in N units for release of DIC and consumption of O2
      !(CH2O)106(NH3)16H3PO4+106O2->106CO2+106H2O+16NH3+H3PO4
      kf = monod_squared(self%K_omox_o2, O2)*f_t(temp,2._rk,self%tref)
      DcDOML_O2 = self%K_DOML_ox*DOML*kf
      DcPOML_O2 = self%K_POML_ox*POML*kf
      DcDOMR_O2 = self%K_DOMR_ox*DOMR*kf
      DcPOMR_O2 = self%K_POMR_ox*POMR*kf
      !DcTOM_O2 = DcDOMR_O2+DcPOMR_O2

      dphy_in_m = carbon_g_to_mole(GrowthPhy)
      dzoo_resp_in_m = carbon_g_to_mole(RespHet)
      ddoml_o2_in_m = carbon_g_to_mole(DcDOML_O2)
      dpoml_o2_in_m = carbon_g_to_mole(DcPOML_O2)
      ddomr_o2_in_m = carbon_g_to_mole(DcDOMR_O2)
      dpomr_o2_in_m = carbon_g_to_mole(DcPOMR_O2)

      !increments
      !alive
      d_Phy = GrowthPhy-MortPhy-ExcrPhy-GrazPhy
      d_Het = self%Uz*Grazing-MortHet-RespHet
      !d_Baae = -GrazBaae
      !d_Baan = -GrazBaan
      !d_Bhae = -GrazBhae
      !d_Bhan = -GrazBhan
      !nutrients
      d_NH4 = N_fixation&
              +(-dphy_in_m*quota(LimNH4,LimN)&
              +dzoo_resp_in_m&
              +ddoml_o2_in_m&
              +dpoml_o2_in_m)/self%c_to_n
      d_NO2 = -dphy_in_m/self%c_to_n*quota(LimNO3,LimN)*quota(NO2,NO2+NO3)
      d_NO3 = -dphy_in_m/self%c_to_n*quota(LimNO3,LimN)*quota(NO3,NO2+NO3)
      d_PO4 =(-dphy_in_m&
              +dzoo_resp_in_m&
              +ddoml_o2_in_m&
              +dpoml_o2_in_m)/self%c_to_p
      d_Si = (-dphy_in_m&
              +dzoo_resp_in_m&
              +ddoml_o2_in_m&
              +dpoml_o2_in_m)/self%c_to_si
      !d_Sipart = (MortPhy+GrazPhy)*self%r_si_n
      !OM
      d_POML = -Autolysis_L-DcPOML_O2+MortPhy+MortHet&
               +Grazing*(1._rk-self%Uz)*(1._rk-self%Hz)&
               -GrazPOP
      d_DOML = Autolysis_L-DcDOML_O2+ExcrPhy&
               +Grazing*(1._rk-self%Uz)*self%Hz
      d_POMR = DcPOML_O2-DcPOMR_O2-Autolysis_R
      d_DOMR = DcDOML_O2-DcDOMR_O2+Autolysis_R
      !
      d_DIC = -dphy_in_m+dzoo_resp_in_m&
              +ddomr_o2_in_m+dpomr_o2_in_m

      d_O2  =  dphy_in_m-dzoo_resp_in_m&
              -ddomr_o2_in_m-dpomr_o2_in_m
      !Components of temporal derivarives calculated in this module:
      d_alk = -d_PO4-d_NO3-d_NO2+d_NH4
      ! -1 mole per 1 mole of NO3- or NO2- or PO4-
      ! +1 mole per 1 mole of NH4+ (Wollf-Gladrow, Zeebe,.. 2007)

      _SET_ODE_(self%id_POML,d_POML)
      _SET_ODE_(self%id_DOML,d_DOML)
      _SET_ODE_(self%id_POMR,d_POMR)
      _SET_ODE_(self%id_DOMR,d_DOMR)

      _SET_ODE_(self%id_NH4,d_NH4)
      _SET_ODE_(self%id_NO2,d_NO2)
      _SET_ODE_(self%id_NO3,d_NO3)
      _SET_ODE_(self%id_PO4,d_PO4)
      _SET_ODE_(self%id_Si,d_Si)
      !_SET_ODE_(self%id_Sipart,d_Sipart)

      _SET_ODE_(self%id_DIC,d_DIC)
      _SET_ODE_(self%id_O2,d_O2)

      _SET_ODE_(self%id_Phy,d_Phy)
      _SET_ODE_(self%id_Het,d_Het)
      !_SET_ODE_(self%id_Baae,d_Baae)
      !_SET_ODE_(self%id_Baan,d_Baan)
      !_SET_ODE_(self%id_Bhae,d_Bhae)
      !_SET_ODE_(self%id_Bhan,d_Bhan)

      _SET_ODE_(self%id_Alk,d_alk)

      O2_sat = oxygen_saturation_concentration(temp, salt)

      !POMTot=POML+POMR
      !DOMTot=DOML+DOMR

      _SET_DIAGNOSTIC_(self%id_LimNH4,LimNH4)
      _SET_DIAGNOSTIC_(self%id_LimN,LimN)
      _SET_DIAGNOSTIC_(self%id_LimP,LimP)
      _SET_DIAGNOSTIC_(self%id_LimNO3,LimNO3)
      _SET_DIAGNOSTIC_(self%id_LimSi,LimSi)
      _SET_DIAGNOSTIC_(self%id_ChlCratio,ChlC)
      _SET_DIAGNOSTIC_(self%id_Biorate,biorate)
      _SET_DIAGNOSTIC_(self%id_growthrate,growthrate)
      _SET_DIAGNOSTIC_(self%id_GrowthPhy,GrowthPhy)
      _SET_DIAGNOSTIC_(self%id_MortPhy,MortPhy)
      _SET_DIAGNOSTIC_(self%id_ExcrPhy,ExcrPhy)

      _SET_DIAGNOSTIC_(self%id_GrazPhy,GrazPhy)
      !_SET_DIAGNOSTIC_(self%id_GrazBhae,GrazBhae)
      !_SET_DIAGNOSTIC_(self%id_GrazBhan,GrazBhan)
      !_SET_DIAGNOSTIC_(self%id_GrazBaae,GrazBaae)
      !_SET_DIAGNOSTIC_(self%id_GrazBaan,GrazBaan)
      _SET_DIAGNOSTIC_(self%id_GrazPOP,GrazPOP)
      _SET_DIAGNOSTIC_(self%id_Grazing,Grazing)
      !_SET_DIAGNOSTIC_(self%id_GrazBact,GrazBact)
      _SET_DIAGNOSTIC_(self%id_MortHet,MortHet)
      _SET_DIAGNOSTIC_(self%id_RespHet,RespHet)

      _SET_DIAGNOSTIC_(self%id_O2_sat,O2_sat)
      _SET_DIAGNOSTIC_(self%id_O2_rel_sat,max(0.0_rk,100.0_rk*O2/O2_sat))
      _SET_DIAGNOSTIC_(self%id_AOU,(O2_sat-O2))
      _SET_DIAGNOSTIC_(self%id_N_fixation,N_fixation)

      !_SET_DIAGNOSTIC_(self%id_DcTOM_O2,DcTOM_O2)
      !_SET_DIAGNOSTIC_(self%id_DOMTot,DOMTot)
      !_SET_DIAGNOSTIC_(self%id_POMTot,POMTot)
      _SET_DIAGNOSTIC_(self%id_DcPOML_O2,DcPOML_O2)
      _SET_DIAGNOSTIC_(self%id_DcPOMR_O2,DcPOMR_O2)
      _SET_DIAGNOSTIC_(self%id_DcDOMR_O2,DcDOMR_O2)
      _SET_DIAGNOSTIC_(self%id_DcDOML_O2,DcDOML_O2)
    _LOOP_END_
  end subroutine do
  !
  !O2 saturation
  !
  subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
    class (type_niva_brom_bio),intent(in) :: self

    _DECLARE_ARGUMENTS_DO_SURFACE_
    real(rk)                   :: O2, temp, salt, windspeed
    real(rk)                   :: Ox, Oa, TempT, Obe, Q_O2

    _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_O2,O2)
      _GET_(self%id_temp,temp) ! temperature
      _GET_(self%id_salt,salt) ! salinity
      _GET_HORIZONTAL_(self%id_windspeed,windspeed)

      Ox = 1953.4_rk-128._rk*temp+3.9918_rk*temp*temp-&
           0.050091_rk*temp*temp*temp !(Wanninkoff, 1992)
      if (Ox>0._rk) then
        Oa = 0.028_rk*(windspeed**3._rk)*sqrt(400._rk/Ox)
      else
        Oa = 0._rk
      end if

      !Calculation of O2 saturation Obe according to UNESCO, 1986
      TempT = (temp+273.15_rk)/100._rk
      Obe = exp(-173.4292_rk+249.6339_rk/TempT+143.3483_rk*&
            log(TempT)-21.8492_rk*TempT+salt*(-0.033096_rk+&
            0.014259_rk*TempT-0.0017_rk*TempT*TempT)) !O2_sat
      Obe = Obe*1000._rk/22.4_rk !convert from ml/l into uM
      Q_O2 = windspeed*(Obe-O2) !After (Burchard et al., 2005)

      _SET_SURFACE_EXCHANGE_(self%id_O2,Q_O2)
    _HORIZONTAL_LOOP_END_
  end subroutine do_surface
  !
  !temperature limiter (q10 type)
  !
  pure real(rk) function f_t(temperature,q10,treference)
    real(rk),intent(in):: temperature
    real(rk),intent(in):: q10
    real(rk),intent(in):: treference

    f_t = exp((temperature - treference)/10 * log(q10))
  end function f_t
  !
  !adapted from ersem
  !
  pure function oxygen_saturation_concentration(ETW,X1X) result(O2_sat)
    real(rk),                      intent(in) :: ETW,X1X
    real(rk)                                  :: O2_sat

    real(rk),parameter :: A1 = -173.4292_rk
    real(rk),parameter :: A2 = 249.6339_rk
    real(rk),parameter :: A3 = 143.3483_rk
    real(rk),parameter :: A4 = -21.8492_rk
    real(rk),parameter :: B1 = -0.033096_rk
    real(rk),parameter :: B2 = 0.014259_rk
    real(rk),parameter :: B3 = -0.0017_rk
    real(rk),parameter :: R = 8.3145_rk
    real(rk),parameter :: P = 101325_rk
    real(rk),parameter :: T = 273.15_rk

    ! volume of an ideal gas at standard temp (25C) and pressure (1 atm)
    real(rk),parameter :: VIDEAL = (R * 298.15_rk / P) *1000._rk

    real(rk)           :: ABT

    ! calc absolute temperature
    ABT = ETW + T

    ! calc theoretical oxygen saturation for temp + salinity
    ! From WEISS 1970 DEEP SEA RES 17, 721-735.
    ! units of ln(ml(STP)/l)
    O2_sat = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
           + A4 * (ABT/100._rk) &
           + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

    ! convert units to ml(STP)/l then to mMol/m3
    O2_sat = exp( O2_sat )
    O2_sat = O2_sat * 1000._rk / VIDEAL
  end function oxygen_saturation_concentration
  !
  !Chl:C relationship, Cloern et al., 1995
  !
  pure real(rk) function chl_c_ratio(temperature, irradiance, nutrient_limiter)
    real(rk),intent(in):: temperature
    real(rk),intent(in):: irradiance
    real(rk),intent(in):: nutrient_limiter

    real(rk):: a0, a, b, c

    a0 = 0.003 ! minimal Chl:C ratio
    a = 0.0154; b = 0.050; c = 0.059 ! achieved by experiment

    chl_c_ratio = a0+a*exp(b*temperature)&
                 *exp(-1._rk*c*irradiance)*nutrient_limiter
  end function chl_c_ratio
  !
  !Calculates amount of hours in a day
  !
  pure real(rk) function photoperiod(latitude, day)
    real(rk),intent(in):: latitude
    real(rk),intent(in):: day

    real(rk):: pi
    real(rk):: a0, a1, a2, a3, b1, b2, b3
    real(rk):: th0, th02, th03, delta, wh

    pi = 3.141592_rk
    a0 = 0.006918_rk
    a1 =-0.399912_rk
    a2 =-0.006758_rk
    a3 =-0.002697_rk
    b1 = 0.070257_rk
    b2 = 0.000907_rk
    b3 = 0.001480_rk

    th0 = pi*day/182.5_rk
    th02 = 2._rk*th0
    th03 = 3._rk*th0

    delta =(a0&
           + a1*cos(th0)+b1*sin(th0)&
           + a2*cos(th02)+b2*sin(th02)&
           + a3*cos(th03)+b3*sin(th03))

    wh = (2._rk*pi)/24._rk
    photoperiod = (2/wh)*acos(-tan((pi/180._rk)*latitude)*tan(delta))
  end function photoperiod
  !
  !D is photoperiod, hours
  !pbm is the maximum hourly rate of photosynthesis, [mg C (mg Chl a h)-1]
  !alpha is photosynthetic efficiency at low irradiance
  !I is instanteneous irradance, PAR [microM quanta m-2 s-1]!
  !
  pure real(rk) function photosynthetic_rate(D, pbm, alpha, I)
    real(rk),intent(in):: D, pbm, alpha, I

    photosynthetic_rate &
        = D*pbm*(1._rk-exp(-1._rk*I*alpha/pbm))
  end function photosynthetic_rate
  !
  !Coefficiens inside evaluate respiration;
  !biorate is the daily rate of photosynthesis, [mg C (mg Chl a d)-1]
  !ChlCratio, joint nutrients limiter
  !
  pure real(rk) function daily_growth(biorate, ChlCratio)
    real(rk),intent(in):: biorate, ChlCratio

    daily_growth &
        = 0.85_rk*biorate*ChlCratio
  end function daily_growth
  !
  !
  !
  pure real(rk) function mortality(subject, rate)
    real(rk),intent(in):: subject, rate

    mortality = 1._rk-exp(-rate*subject)
  end function mortality
  !
  !
  !
  !real(rk) function graz(self, var, Het)
  !  class (type_niva_brom_bio),intent(in):: self
  !  real(rk),intent(in):: var, Het

  !  graz = self%K_het_bac_gro*Het&
  !       * monod_squared(self%limGrazBac, quota(var, Het))
  !end function graz
end module fabm_niva_brom_bio
