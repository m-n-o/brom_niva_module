#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
   module fabm_niva_brom_hg
!
! !DESCRIPTION:
!
! !USES:

   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !REVISION HISTORY:!
!  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev
!

! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_niva_brom_hg
!     Variable identifiers
      type (type_state_variable_id)        :: id_H2S, id_O2, id_Mn4, id_Fe2,id_Fe3,id_Baae,id_Bhae,id_Baan,id_Bhan
      type (type_state_variable_id)        :: id_Phy, id_Het, id_NH4, id_PON, id_DON
      type (type_dependency_id)            :: id_temp,id_par,id_depth
      type (type_state_variable_id)        :: id_Hg0,id_Hg2,id_MeHg,id_HgS
      type (type_state_variable_id)        :: id_MeHg_biota, id_MeHg_POM, id_MeHg_DOM
      type (type_diagnostic_variable_id)   :: id_MeHg_Fe3, id_MeHg_Mn4, id_MeHg_free, id_MeHg_tot, id_MeHg_tot_diss
      type (type_state_variable_id)        :: id_Hg2_biota, id_Hg2_POM, id_Hg2_DOM
      type (type_diagnostic_variable_id)   :: id_Hg2_Fe3, id_Hg2_Mn4, id_Hg2_free, id_Hg2_tot, id_Hg2_tot_diss
      type (type_diagnostic_variable_id)   :: id_Hg_tot_diss, id_Hg_tot, id_dMeHg
      !---- Hg---------!
      real(rk) :: K_hg2_mehg , K_mehg_hg2, K_hg2_hg0, K_hg0_hg2
      real(rk) :: K_mehg_irr_degr ! photo-degradation
      real(rk) :: K_hg0_irr_ox ! photo-oxidation of hg0    
      real(rk) :: K_hg2_irr_red ! photo-reduction of Hg2      
      real(rk) :: K_HgS, K_hgs_form, K_hgs_ox, K_hgs_diss, O2s_nf      
      real(rk) ::  r_fe3_mehg,  r_mn4_mehg, r_fe3_hg2,  r_mn4_hg2
!     Model parameters
      real(rk) :: Wsed, Wphy, Wm            

   contains
      procedure :: initialize
      procedure :: do
      procedure :: do_surface
   end type
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
   class (type_niva_brom_hg), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
! !REVISION HISTORY:
!  Original author(s): 
!
!EOP
!-----------------------------------------------------------------------
!BOC
!-----------------------------------------------------------------------   
! Parameters, i.e. rate constants  
   !---- Hg---------!  
   call self%get_parameter(self%K_HgS,          'K_HgS',          ' - ', 'HgS equilibrium constant',       default=0.0_rk)     
   call self%get_parameter(self%K_hg2_mehg,     'K_hg2_mehg',     ' - ', 'Coef. of mercury methilation',   default=0.0_rk)  
   call self%get_parameter(self%K_mehg_hg2,     'K_mehg_hg2',     ' - ', 'Coef. of demethylation of MeHg', default=0.0_rk) 
   call self%get_parameter(self%K_hg2_hg0,      'K_hg2_hg0',      ' - ', 'Coef. biotic reduction of Hg2',  default=0.0_rk)    
   call self%get_parameter(self%K_hg0_hg2,      'K_hg0_hg2',      ' - ', 'Coef. of dark oxidation of Hg0', default=0.0_rk)    
   call self%get_parameter(self%K_mehg_irr_degr,'K_mehg_irr_degr',' - ', 'photo-degradation of MeHg',      default=0.0_rk) 
   call self%get_parameter(self%K_hg0_irr_ox,   'K_hg0_irr_ox',   ' - ', 'photo-oxidation of hg0',         default=0.0_rk)    
   call self%get_parameter(self%K_hgs_form,     'K_hgs_form',     ' - ', 'Formation of HgS',               default=0.0_rk) 
   call self%get_parameter(self%K_hgs_ox,       'K_hgs_ox',       ' - ', 'Oxidation of HgS',               default=0.0_rk) 
   call self%get_parameter(self%K_hgs_diss,     'K_hgs_diss',     ' - ', 'Dissolution of HgS',             default=0.0_rk)  
   call self%get_parameter(self%K_hg2_irr_red,  'K_hg2_irr_red',' - ', 'photo-reduction of Hg2 ',        default=0.0_rk)    
   call self%get_parameter(self%r_fe3_mehg, 'r_fe3_mehg', '[-]',' Fe3[uM]/MeHg[uM] partitioning coeff. for Fe3', default=833.0_rk)
   call self%get_parameter(self%r_mn4_mehg,  'r_mn4_mehg',  '[-]',' MnO2[uM]/MeHg[uM] partitioning coeff. for MnO2', default=833.0_rk)
   call self%get_parameter(self%r_fe3_mehg, 'r_fe3_hg2', '[-]',' Fe3[uM]/Hg2[uM] partitioning coeff. for Fe3', default=833.0_rk)
   call self%get_parameter(self%r_mn4_mehg,  'r_mn4_hg2',  '[-]',' MnO2[uM]/Hg2[uM] partitioning coeff. for MnO2', default=833.0_rk)

   call self%get_parameter(self%Wsed, 'Wsed', '[m/day]',  'Rate of sinking of detritus (POP, PON)',       default=5.00_rk)     
   call self%get_parameter(self%Wphy, 'Wphy', '[m/day]',  'Rate of sinking of Phy',                       default=0.10_rk)
   call self%get_parameter(self%Wm,   'Wm','   [m/day]',  'Rate of accelerated sinking of metals',        default=7.0_rk)
   
   call self%get_parameter(self%O2s_nf, 'O2s_nf', '[uM O]','half saturation for nitrification',default=4.488_rk)
      !---- Hg---------!      
   call self%register_state_variable(self%id_Hg0, 'Hg0',  'mmol/m**3','Hg(0)',       minimum=0.0_rk)
   call self%register_state_variable(self%id_Hg2, 'Hg2',  'mmol/m**3','Hg(II)',       minimum=0.0_rk)
   call self%register_state_variable(self%id_MeHg,'MeHg', 'mmol/m**3','MeHg',  minimum=0.0_rk)
   call self%register_state_variable(self%id_HgS, 'HgS',  'mmol/m**3','HgS', minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)
   call self%register_state_variable(self%id_MeHg_biota,'MeHg_biota','mmol/m**3','MeHg_biota', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_MeHg_POM,'MeHg_POM','mmol/m**3','MeHg_POM', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_MeHg_DOM,'MeHg_DOM','mmol/m**3','MeHg_DOM', minimum=0.0_rk)
   call self%register_state_variable(self%id_Hg2_biota,'Hg2_biota','mmol/m**3','Hg2_biota', minimum=0.0_rk,vertical_movement=-self%Wphy/86400._rk)
   call self%register_state_variable(self%id_Hg2_POM,'Hg2_POM','mmol/m**3','Hg2_POM', minimum=0.0_rk,vertical_movement=-self%Wsed/86400._rk)
   call self%register_state_variable(self%id_Hg2_DOM,'Hg2_DOM','mmol/m**3','Hg2_DOM', minimum=0.0_rk)

   call self%register_state_dependency(self%id_O2,  'O2',    'mmol/m**3', 'O2')
   call self%register_state_dependency(self%id_Mn4, 'Mn4',   'mmol/m**3', 'Mn4')
   call self%register_state_dependency(self%id_Fe2, 'Fe2',   'mmol/m**3', 'Fe2')
   call self%register_state_dependency(self%id_Fe3, 'Fe3',   'mmol/m**3', 'Fe3')
   call self%register_state_dependency(self%id_H2S, 'H2S',   'mmol/m**3',  'H2S')
   call self%register_state_dependency(self%id_Phy, 'Phy', 'mmol/m**3','Phy')
   call self%register_state_dependency(self%id_Het, 'Het', 'mmol/m**3','Het')
   call self%register_state_dependency(self%id_Baae, 'Baae', 'mmol/m**3','Aerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhae, 'Bhae', 'mmol/m**3','Aerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_Baan, 'Baan', 'mmol/m**3','Anaerobic Autotrophic Bacteria')
   call self%register_state_dependency(self%id_Bhan, 'Bhan', 'mmol/m**3','Anaerobic Heterotrophic Bacteria')
   call self%register_state_dependency(self%id_PON,'PON','mmol/m**3','particulate organic nitrogen')
   call self%register_state_dependency(self%id_DON,'DON','mmol/m**3','dissolved organic nitrogen')
   
    call self%register_diagnostic_variable(&
         self%id_Hg2_Mn4,'Hg2_Mn4','mmol/m**3',&
         'Hg2 adsorped on Mn4',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Hg2_Fe3,'Hg2_Fe3','mmol/m**3',&
         'Hg2 adsorped on Fe3',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_Hg2_free,'Hg2_free','mmol/m**3',&
         'Hg2+ free + part on biota',&
         output=output_time_step_integrated) ! dissolved inorganic
    call self%register_diagnostic_variable(&
         self%id_Hg2_tot_diss,'Hg2_tot_diss','mmol/m**3',&
         'Hg2+ free dissolved',&
         output=output_time_step_integrated) ! dissolved organic and inorganic
    call self%register_diagnostic_variable(&
         self%id_Hg2_tot,'Hg2_total','mmol/m**3',&
         'Hg2 total',&
         output=output_time_step_integrated) ! dissolved and particulate
    
   call self%register_diagnostic_variable(&
         self%id_MeHg_Mn4,'MeHg_Mn4','mmol/m**3',&
         'MeHg adsorped on Mn4',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MeHg_Fe3,'MeHg_Fe3','mmol/m**3',&
         'MeHg adsorped on FeS',&
         output=output_time_step_integrated)
    call self%register_diagnostic_variable(&
         self%id_MeHg_free,'MeHg_free','mmol/m**3',&
         'MeHg+ free + part on biota',&
         output=output_time_step_integrated) ! dissolved inorganic
    call self%register_diagnostic_variable(&
         self%id_MeHg_tot_diss,'MeHg_tot_diss','mmol/m**3',&
         'MeHg+ free dissolved',&
         output=output_time_step_integrated) ! dissolved organic and inorganic
    call self%register_diagnostic_variable(&
         self%id_MeHg_tot,'MeHg_total','mmol/m**3',&
         'MeHg total',&
         output=output_time_step_integrated) ! dissolved and particulate 
    call self%register_diagnostic_variable(&
         self%id_Hg_tot_diss,'Hg_tot_diss','mmol/m**3',&
         'Hg total dissolved',&
         output=output_time_step_integrated) ! dissolved  
    call self%register_diagnostic_variable(&
         self%id_Hg_tot,'Hg_total','mmol/m**3',&
         'Hg total',&
         output=output_time_step_integrated) ! dissolved and particulate 
    call self%register_diagnostic_variable(&
         self%id_dMeHg,'dMeHg','%',&
         'dMeHg',&
         output=output_time_step_integrated) ! dissolved and particulate 

   call self%register_dependency(self%id_temp, standard_variables%temperature) 
   call self%register_dependency(self%id_depth,standard_variables%pressure)  
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
!   call self%register_dependency(self%id_Izt,'Izt','W/m2','downwelling_photosynthetic_radiative_flux')
! Specify that are rates computed in this module are per day (default: per second)
   self%dt = 86400.

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
! !DESCRIPTION: This module descibes biogeochemical transformation of mercury (Hg) in the seawater
! 
!
! !INPUT PARAMETERS:
   class (type_niva_brom_hg),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !REVISION HISTORY:
!  Original author(s): Jorn Bruggeman, Svetlana Pakhomova, Evgeniy Yakushev
!
! !LOCAL VARIABLES:
   real(rk) ::  temp, O2, Mn4, Fe2, Fe3, depth
   real(rk) ::  Hg0, Hg2, MeHg, HgS, Iz, H2S, Om_HgS, Hg2_flux 
   real(rk) ::  Phy, Het, Bhae, Baae, Bhan, Baan, DON, PON
   real(rk) ::  Hg2_biota,Hg2_POM,Hg2_DOM
   real(rk) ::  Hg2_Mn4, Hg2_Fe3, Hg2_free, Hg2_tot, Hg2_tot_diss
   real(rk) ::  MeHg_biota,MeHg_POM,MeHg_DOM
   real(rk) ::  MeHg_Mn4, MeHg_Fe3, MeHg_free, MeHg_tot, MeHg_tot_diss
   real(rk) ::  Hg_tot, Hg_tot_diss, dMeHg
   real(rk) ::  hgs_form, hgs_diss, hgs_ox, hg2_hg0, hg0_hg2
   real(rk) ::  hg0_irr_ox, hg2_irr_red, mehg_irr_degr, hg2_mehg, mehg_hg2
   real(rk) ::  dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM
!EOP 
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Environment
   _GET_(self%id_temp,temp)              ! temperature
   _GET_(self%id_depth,depth)            ! local photosynthetically active radiation

   _GET_(self%id_par,Iz)              ! local photosynthetically active radiation

    ! Our own state variables
    _GET_(self%id_Hg0,Hg0)
    _GET_(self%id_Hg2,Hg2)
    _GET_(self%id_Hg2_biota,Hg2_biota) 
    _GET_(self%id_Hg2_POM,Hg2_POM) 
    _GET_(self%id_Hg2_DOM,Hg2_DOM)
    _GET_(self%id_MeHg,MeHg) 
    _GET_(self%id_MeHg_biota,MeHg_biota) 
    _GET_(self%id_MeHg_POM,MeHg_POM) 
    _GET_(self%id_MeHg_DOM,MeHg_DOM)
    _GET_(self%id_HgS,HgS)
    _GET_(self%id_Phy,Phy) 
    _GET_(self%id_Het,Het)
    _GET_(self%id_PON,PON)
    _GET_(self%id_DON,DON)
    _GET_(self%id_Baae,Baae)   
    _GET_(self%id_Bhae,Bhae)
    _GET_(self%id_Baan,Baan)   
    _GET_(self%id_Bhan,Bhan)
    _GET_(self%id_H2S,H2S)
    _GET_(self%id_Mn4,Mn4)    
    _GET_(self%id_Fe2,Fe2)
    _GET_(self%id_Fe3,Fe3)
    _GET_(self%id_O2,O2)   
    
    ! Hg species (Knigthes 2008)
!% Hg0 bioreduction  Hg0 -> Hg2+  ()
    hg0_hg2=self%K_hg0_hg2*Hg0         
!% Hg2 biooxydation  Hg2+ + 0.5O2 + 2H+-> Hg0 + H2O   ()
    hg2_hg0=self%K_hg2_hg0*Hg2*0.5*(1.+tanh(o2-self%O2s_nf))  
!% Hg2 methylation Hg2+  -> MeHg   ()
    hg2_mehg=self%K_hg2_mehg*Hg2*Bhan
!% MeHg demethylation MeHg  -> Hg2+   ()
    mehg_hg2=self%K_mehg_hg2*MeHg
     if (H2S>400) mehg_hg2=self%K_mehg_hg2*MeHg*10.
    Om_HgS=H2S*Hg2/(self%K_HgS) 
!% HgS formation Hg2+ + H2S -> HgS + 2H+ ()
    hgs_form=max(0._rk,self%K_hgs_form*max(0._rk,(Om_HgS-1._rk)))
!    if (Hg2<0.000001.or.H2S<0.00001) hgs_form=0._rk
!% HgS dissolution  HgS + 2H+ -> Hg2+ + H2S   ()
    hgs_diss=self%K_hgs_diss*HgS*max(0._rk,(1._rk-Om_HgS))
!    if (HgS<0.000001) hgs_diss=0._rk 
!% HgS oxydation  HgS + 2O2 -> Hg2+ + SO4   ()
    hgs_ox=self%K_hgs_ox*HgS*O2     
!% Hg2 biooxydation  Hg2+ + 0.5O2 + 2H+-> Hg0 + H2O   ()
!    hg2_hg0=self%K_hg2_hg0*Hg2*O2     
!% Hg2 photo reduction  Hg2+ -> Hg0   ()
    hg2_irr_red=self%K_hg2_irr_red*Hg2*Iz/80.
!% Hg0 photo oxydation  Hg0 -> Hg2+   ()
    hg0_irr_ox=self%K_hg0_irr_ox*Hg0*Iz/80.
!% Hg2 methylation Hg2+  -> MeHg   ()
!    hg2_mehg=self%K_hg2_mehg*Hg2*Bhan
!% MeHg demethylation MeHg  -> Hg2+   ()
!    mehg_hg2=self%K_mehg_hg2*MeHg
!% MeHg photo degradation MeHg  -> Hg2+   
    mehg_irr_degr=self%K_mehg_irr_degr*MeHg

    call partit (Hg2, Hg2_biota, Hg2_POM, Hg2_DOM, &
                 Phy, Het, Baae, Bhae, Baan, Bhan, PON, DON, &
                 dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM)
    

    !!  Total Hg2 free concentration from the previous time step
    Hg2_free = Hg2-Hg2_biota-Hg2_POM-Hg2_DOM
!
!! calculate adsorption  (impossible for small Hg2 content), 
!! and Hg2_free remained after adsorption, that will be 
!! available for mineral formation and partitioning with biota.    
    Hg2_Fe3 = min(Fe3/self%r_fe3_hg2, Hg2_free)
    Hg2_free = max (0.0, Hg2_free-Hg2_Fe3)
    Hg2_Mn4 = min(Mn4/self%r_mn4_hg2, Hg2_free)
    Hg2_free = max (0.0, Hg2_free-Hg2_Mn4)



    Hg2_tot=Hg2_free+Hg2_Mn4+Hg2_Fe3+Hg2_biota+Hg2_POM+Hg2_DOM
    Hg2_tot_diss=Hg2_free+Hg2_DOM
    
   _SET_ODE_(self%id_Hg2_biota,dSubst_biota)
   _SET_ODE_(self%id_Hg2_POM,dSubst_POM)
   _SET_ODE_(self%id_Hg2_DOM,dSubst_DOM)
   _SET_ODE_(self%id_Hg2, hg0_hg2-hg2_hg0-hg2_mehg+mehg_hg2-hgs_form+hgs_diss+mehg_irr_degr+hg0_irr_ox-hg2_irr_red+dSubst_dis) !-dSubst_POM-dSubst_DOM-dSubst_biota)

   
   !!  Total MeHg
   
    call partit (MeHg, MeHg_biota, MeHg_POM, MeHg_DOM, &
                 Phy, Het, Baae, Bhae, Baan, Bhan, PON, DON, &
                 dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM)
!
!!  Total MeHg free concentration from the previous time step
    MeHg_free = MeHg-MeHg_biota-MeHg_POM-MeHg_DOM
!
!! calculate adsorption  (impossible for small MeHg content), 
!! and MeHg_free remained after adsorption, that will be 
!! available for mineral formation and partitioning with biota.    
      MeHg_Mn4 = min(Mn4/self%r_mn4_mehg, MeHg_free)
    MeHg_free = max (0.0, MeHg_free-MeHg_Mn4)
      MeHg_Fe3 = min(Fe3/self%r_fe3_mehg, MeHg_free)
    MeHg_free = max (0.0, MeHg_free-MeHg_Fe3)


    MeHg_tot=MeHg_free+MeHg_Mn4+MeHg_Fe3+MeHg_biota+MeHg_POM+MeHg_DOM
    MeHg_tot_diss=MeHg_free+MeHg_DOM
    
    
    Hg_tot_diss=MeHg_tot_diss+Hg2_tot_diss+Hg0
    Hg_tot=MeHg_tot+Hg2_tot+Hg0
    dMeHg=MeHg_tot_diss/Hg_tot_diss


   _SET_ODE_(self%id_MeHg,hg2_mehg-mehg_hg2-mehg_irr_degr+dSubst_dis) !-dSubst_POM-dSubst_DOM-dSubst_biota)
   _SET_ODE_(self%id_MeHg_biota,dSubst_biota)
   _SET_ODE_(self%id_MeHg_POM,dSubst_POM)
   _SET_ODE_(self%id_MeHg_DOM,dSubst_DOM)
   
   _SET_ODE_(self%id_Hg0, -hg0_hg2+hg2_hg0+hgs_ox-hg0_irr_ox+hg2_irr_red)
   _SET_ODE_(self%id_HgS,hgs_form-hgs_diss-hgs_ox)
   _SET_ODE_(self%id_Baan,0.0_rk)   
   _SET_ODE_(self%id_Bhan,0.0_rk)
   _SET_ODE_(self%id_Fe3,0.0_rk)
   _SET_ODE_(self%id_Mn4,0.0_rk)
   _SET_ODE_(self%id_H2S,-hgs_form+hgs_diss+hgs_ox)
   _SET_ODE_(self%id_O2,-hg2_hg0-hgs_ox)
   
   
      _SET_DIAGNOSTIC_(self%id_Hg2_Mn4,Hg2_Mn4)
      _SET_DIAGNOSTIC_(self%id_Hg2_Fe3,Hg2_Fe3)
      _SET_DIAGNOSTIC_(self%id_Hg2_free,Hg2_free)
      _SET_DIAGNOSTIC_(self%id_Hg2_tot,Hg2_tot)
      _SET_DIAGNOSTIC_(self%id_Hg2_tot_diss,Hg2_tot_diss)
      _SET_DIAGNOSTIC_(self%id_MeHg_Mn4,MeHg_Mn4)
      _SET_DIAGNOSTIC_(self%id_MeHg_Fe3,MeHg_Fe3)
      _SET_DIAGNOSTIC_(self%id_MeHg_free,MeHg_free)
      _SET_DIAGNOSTIC_(self%id_MeHg_tot,MeHg_tot)
      _SET_DIAGNOSTIC_(self%id_MeHg_tot_diss,MeHg_tot_diss)
      _SET_DIAGNOSTIC_(self%id_Hg_tot_diss,Hg_tot_diss)
      _SET_DIAGNOSTIC_(self%id_Hg_tot,Hg_tot)
      _SET_DIAGNOSTIC_(self%id_dMeHg,dMeHg)


! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC
!-----------------------------------------------------------------------
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!
! !DESCRIPTION:
! Sea water Hg(0) exchange.   
   
! !INPUT PARAMETERS:
   class (type_niva_brom_hg),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!
! !LOCAL VARIABLES:
   real(rk) :: pCO2w, xk, Ox, Q_Hg0, Hg0
   real(rk) :: temp, Kc0, salt
   real(rk) :: Sc, TK, fwind !PML
   real(rk) :: Hg0_air
   real(rk) :: windspeed

   _HORIZONTAL_LOOP_BEGIN_
      _GET_(self%id_temp,temp)              ! temperature
      _GET_(self%id_Hg0,Hg0)              ! temperature
      
      TK=(temp + 273.15)
      windspeed=5.
      Hg0_air=5.0e-8 ! 0.01/200./1000. !convert from ng/l into mmol/m3

! PML
! calculate the Scmidt number and unit conversions
      Sc = 2073.1_rk-125.62_rk*temp+3.6276_rk*temp**2._rk-0.043219_rk*&
           temp**3.0_rk
      fwind = (0.222_rk*windspeed**2_rk+0.333_rk*windspeed)*&
              (Sc/660._rk)**(-0.5_rk)
      fwind=fwind*24._rk/100._rk !convert to m/day
! flux depends on the difference in partial pressures, wind and henry
! here it is rescaled to mmol/m2/d
!          flux = fwind * HENRY * ( PCO2A - PCO2W ) * dcf      

      Q_Hg0= fwind * (Hg0_air- max(0e0,Hg0))

      _SET_SURFACE_EXCHANGE_(self%id_Hg0,Q_Hg0)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
   
   
   

!-----------------------------------------------------------------------
   subroutine partit (Subst_dis,Subst_biota, Subst_POM, Subst_DOM, &
                      Phy, Het, Baae, Bhae, Baan, Bhan, PON, DON, &
                      dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM)
   ! !LOCAL VARIABLES:
   real(rk) :: Subst_dis, Subst_biota, Subst_POM, Subst_DOM,  Subst_tot 
   real(rk) :: Phy, Het, Baae, Bhae, Baan, Bhan, PON, DON
   real(rk) :: dSubst_dis, dSubst_biota, dSubst_POM, dSubst_DOM, dSubst_tot
   real(rk) :: dSubst_tot_diss, dSubst_tot_part
   real(rk) :: pol_bio  ! pollutant in BIOTA,"ng?"/l
   real(rk) :: pol_dom  ! pollutant in DOM, "ng?"/l
   real(rk) :: pol_pom  ! pollutant in POM, "ng?"/l
   
   real(rk) :: Subst_total   ! total pollutant, "ng?"/l
   real(rk) :: Subst_free     ! pollutant in dissolved INORGANIC,"ng?"/l, i.e. not partitioned  
   
   real(rk) :: sha_bio ! % share of polutant in living organisms
   real(rk) :: sha_pom ! % share of polutant in POM
   real(rk) :: sha_dom ! % share of polutant in DOM
   real(rk) :: sha_free ! % share of 'free' polutant 
   real(rk) :: uMn2lip=0.0084 !coeff.to transfer POM (umol/l N)->(g DryWeight/l) 
   !real(rk) :: rho_FeS= 5.90E7 !    # Density of FeS [mmolFe/m3] (default = 5.90E7 mmolFe/m3)
   !real(rk) :: rho_FeS2=4.17E7 !    # Density of FeS2 [mmolFe/m3] (default = 4.17E7 mmolFe/m3)
   real(rk) :: rho_Mn4= 5.78E7 !    # Density of Mn4 [mmolMn/m3] (default = 5.78E7 mmolMn/m3)   
   real(rk) :: rho_Fe3= 5.90E7 !    # Density of Fe3 [mmolFe/m3] (default = 5.90E7 mmolFe/m3)

!====================================================

!  Ni
  real(rk) :: Kow_bio = 39811.  ! part.coeff. BIO/water (Allisson, 2005, Table 1, in L/kg)
  real(rk) :: Kow_pom = 39811.  ! part.coeff. POM/water (Allisson, 2005, Table 1, in L/kg)
  real(rk) :: Kow_dom = 125893. ! part.coeff. DOM/water (Allisson, 2005, Table 1, in L/kg)
! /Ni   
!====================================================
!EOP 
!---
   
! Let's first assume that all the polutant is dissolved INORGANIC...
        Subst_total = Subst_dis+ Subst_biota+ Subst_POM+ Subst_DOM ! total amount of pollutant 

! We assume that density of organic matter is the same as that of 
!  water, i.e. 1 g=1 ml and operate with weight units to caluclate 
!  the shares of pollutant partitioning medias:
       if((Phy+Het+Baae+Bhae+Baan+Bhan)<=0.) then 
        sha_bio=0. 
       else
        sha_bio=uMn2lip/1000.*(Phy+Het+Baae+Bhae+Baan+Bhan)  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif 
       
       if(PON<=0.) then 
        sha_pom=0. 
       else
        sha_pom=uMn2lip/1000.*PON  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif     
       
       if(DON<=0.) then 
        sha_dom=0. 
       else
        sha_dom=uMn2lip/1000.*DON  ! Volume(weight in kg, g->kg=/1000) of BIO
       endif    

      sha_free = 1.-sha_bio-sha_pom-sha_dom ! i.e Volume(weight in [kg]) of 1l of water minus volumes of org. and part. forms 
!! The free subst.conc. left as dissolved free
      Subst_free = Subst_total* sha_free /(sha_free +Kow_bio*sha_bio              &
     &            +Kow_pom*sha_pom +Kow_dom*sha_dom) ! Kow_water=1. needed for correcn units is excluded
!! subst.conc. partitioned to biota
     pol_bio=max(0.,Kow_bio*Subst_free*sha_bio/sha_free)
!! subst.conc. partitioning to POM
     pol_pom=max(0.,Kow_pom*Subst_free*sha_pom/sha_free)
!! subst.conc. partitioning to DOM
     pol_dom=max(0.,Kow_dom*Subst_free*sha_dom/sha_free)

! difference betweeen new and old state variable, needed for FABM
    dSubst_dis   = -Subst_dis   +Subst_free
    dSubst_biota = -Subst_biota +pol_bio
    dSubst_POM   = -Subst_POM   +pol_pom
    dSubst_DOM   = -Subst_DOM   +pol_dom

   
   end subroutine partit
end module