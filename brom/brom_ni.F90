#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: 
!
! !INTERFACE:
module fabm_niva_brom_ni
    ! 
    ! This module describes the transformation of Ni connceted with 
    ! formation/dissolution of minerals and adsorption on solids.
    ! !USES:
    use fabm_types
    implicit none
    !  default: all is private.
    private
    !
    ! !REVISION HISTORY: version 25.03.2017
    !  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev
    ! !PUBLIC DERIVED TYPES:
    type,extends(type_base_model),public :: type_niva_brom_ni
        type (type_state_variable_id)         :: id_H2S, id_Mn4, id_FeS, id_FeS2
        type (type_state_variable_id)         :: id_Ni, id_NiS, id_Ni_biota, id_Ni_POM, id_Ni_DOM
        type (type_diagnostic_variable_id)    :: id_Ni_FeS, id_Ni_FeS2, id_Ni_Mn4, id_Ni_free, id_Ni_tot
        type (type_diagnostic_variable_id)    :: id_NiS_diss, id_NiS_form, id_Ni_tot_diss

        real(rk) ::  K_NiS, K_NiS_form, K_NiS_diss
        real(rk) ::  r_fes_ni, r_fes2_ni, r_mn4_ni
        real(rk) :: Wsed, Wphy, Wm                  
    contains
        procedure :: initialize
        procedure :: do
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
        ! !DESCRIPTION:
        ! This module described transformation of Ni connceted with 
        ! formation/dissolution of minerals and adsorption on solids.
        ! !INPUT PARAMETERS:
        class (type_niva_brom_ni), intent(inout), target :: self
        integer,                     intent(in)            :: configunit
        !
        ! !REVISION HISTORY:
        !  Original author(s): Svetlana Pakhomova, Evgeniy Yakushev 
        !
        !EOP
        !-----------------------------------------------------------------------
        !BOC
        !-----------------------------------------------------------------------   
        ! Parameters, i.e. rate constants  
        call self%get_parameter(self%Wsed, 'Wsed', '[m/day]', &
             'Rate of sinking of detritus (POP, POM)', &
              default=5.00_rk)     
        call self%get_parameter(self%Wphy, 'Wphy', '[m/day]',&
             'Rate of sinking of Phy',& 
              default=0.10_rk)
        call self%get_parameter(self%Wm,'Wm','[m/day]', &
            'Rate of accelerated sinking of metals', & 
            default=7.0_rk)
        call self%get_parameter(self%K_NiS, 'K_NiS','[uM]',&
             'K_NiS equilibrium constant (Solubility Product Constant)', &
              default=2510.0_rk)
        call self%get_parameter(self%K_NiS_form, 'K_NiS_form', '[1/day]', &
             'Specific rate of precipitation of NiS from Ni with H2S', &
              default=5.e-5_rk)
        call self%get_parameter(self%K_NiS_diss, 'K_NiS_diss', '[1/day]', &
             'Specific rate of dissollution of NiS to Ni and H2S', &
              default=1.e-6_rk)
        call self%get_parameter(self%r_fes_ni,  'r_fes_ni',  '[-]', &
              'FeS[uM]/Ni[uM] partitioning coeff. for FeS', & 
              default=833.0_rk)
        call self%get_parameter(self%r_fes2_ni, 'r_fes2_ni', '[-]', &
            'FeS2[uM]/Ni[uM] partitioning coeff. for FeS2', &
            default=833.0_rk)
        call self%get_parameter(self%r_mn4_ni,  'r_mn4_ni', &
            '[-]',' MnO2[uM]/Ni[uM] partitioning coeff. for MnO2', &
            default=833.0_rk)
        !---- Ni---------!      
        call self%register_state_variable(self%id_Ni, &    
                    'Ni','mmol/m**3','Ni', minimum=0.0_rk)
        call self%register_state_variable(self%id_Ni_biota, &
                    'Ni_biota','mmol/m**3','Ni_biota', minimum=0.0_rk, &
                    vertical_movement=-self%Wphy/86400._rk)
        call self%register_state_variable(self%id_Ni_POM, &
                    'Ni_POM','mmol/m**3','Ni_POM',minimum=0.0_rk, &
                    vertical_movement=-self%Wsed/86400._rk)
        call self%register_state_variable(self%id_Ni_DOM, &
                    'Ni_DOM','mmol/m**3','Ni_DOM', minimum=0.0_rk)
        call self%register_state_variable(self%id_NiS, &
                    'NiS',  'mmol/m**3','NiS',minimum=0.0_rk, &
                    vertical_movement=-self%Wm/86400._rk)
        
        !   call self%register_state_variable(self%id_NiS2,  
        !'NiS2',  'mmol/m**3','NiS2',       
        ! minimum=0.0_rk,vertical_movement=-self%Wm/86400._rk)

        call self%register_state_dependency(self%id_Mn4, 'Mn4',   'mmol/m**3', 'Mn4')
        call self%register_state_dependency(self%id_FeS, 'FeS',   'mmol/m**3', 'FeS')
        call self%register_state_dependency(self%id_FeS2,'FeS2',  'mmol/m**3', 'FeS2')
        call self%register_state_dependency(self%id_H2S, 'H2S',   'mmol/m**3', 'H2S')

        call self%register_diagnostic_variable(&
                self%id_Ni_Mn4,'Ni_Mn4','mmol/m**3',&
                'Ni adsorped on Mn4',&
                output=output_time_step_integrated)
        call self%register_diagnostic_variable(&
                self%id_Ni_FeS,'Ni_FeS','mmol/m**3',&
                'Ni adsorped on FeS',&
                output=output_time_step_integrated)
        call self%register_diagnostic_variable(&
                self%id_Ni_FeS2,'Ni_FeS2','mmol/m**3',&
                'Ni adsorped on FeS2',&
                output=output_time_step_integrated)
        call self%register_diagnostic_variable(&
                self%id_Ni_free,'Ni_free','mmol/m**3',&
                'Ni2+ free dissolved',&
                output=output_time_step_integrated) ! dissolved inorganic
        call self%register_diagnostic_variable(&
                self%id_Ni_tot_diss,'Ni_tot_diss','mmol/m**3',&
                'Ni2+ free dissolved',&
                output=output_time_step_integrated) ! dissolved organic and inorganic
        call self%register_diagnostic_variable(&
                self%id_Ni_tot,'Ni_total','mmol/m**3',&
                'Ni total',&
                output=output_time_step_integrated) ! dissolved and particulate 
        call self%register_diagnostic_variable(&
                self%id_NiS_form,'NiS_form','mmol/m**3',&
                'NiS formation rate',&
                output=output_time_step_integrated)
        call self%register_diagnostic_variable(&
                self%id_NiS_diss,'NiS_diss','mmol/m**3/d',&
                'NiS dissolution rate',&
                output=output_time_step_integrated)

        self%dt = 86400.
   
    end subroutine initialize
    !EOC

    !-----------------------------------------------------------------------
    !BOP
    !IROUTINE: 
    !INTERFACE:
    subroutine do(self,_ARGUMENTS_DO_)
        ! !INPUT PARAMETERS:
        class (type_niva_brom_ni),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
        ! REVISION HISTORY: 
        ! Original author(s): Svetlana Pakhomova, Evgeniy Yakushev 
        !
        ! LOCAL VARIABLES:
        real(rk) ::  Ni, NiS, Ni_biota, Ni_POM, Ni_DOM, d_Ni
        real(rk) ::  Ni_Mn4, Ni_FeS, Ni_FeS2, Ni_free, Ni_tot, Ni_tot_diss
        real(rk) ::  temp, O2, Mn4, FeS, FeS2, depth, Hplus
        real(rk) ::   H2S 
        real(rk) ::  Om_NiS, NiS_form, NiS_diss
        !EOP 
        !-----------------------------------------------------------------------
        !BOC
        ! Enter spatial loops (if any)
        _LOOP_BEGIN_
        _GET_(self%id_Ni_biota,Ni_biota)
        _GET_(self%id_Ni_POM,Ni_POM)
        _GET_(self%id_Ni_DOM,Ni_DOM)    
        _GET_(self%id_Ni,Ni) 
        _GET_(self%id_NiS,NiS) 

        _GET_(self%id_H2S,H2S)
        _GET_(self%id_FeS,FeS)
        _GET_(self%id_FeS2,FeS2)
        _GET_(self%id_Mn4,Mn4)    

        !Total Ni free concentration from the previous time step
        Ni_free = Ni+Ni_biota+Ni_POM+Ni_DOM

        ! calculate adsorption  (impossible for small Ni content), 
        ! and Ni_free remained after adsorption, that will be 
        ! available for mineral formation and partitioning with biota.    
        Ni_Mn4 = min(Mn4/self%r_mn4_ni, Ni_free)
        Ni_free = max(0.0, Ni_free-Ni_Mn4)
        Ni_FeS = min(FeS/self%r_fes_ni, Ni_free)
        Ni_free = max(0.0, Ni_free-Ni_FeS)
        Ni_FeS2= min(FeS2/self%r_fes2_ni, Ni_free)
        Ni_free = max(0.0, Ni_free-Ni_FeS2)

        Ni_tot=Ni_free+Ni_Mn4+Ni_FeS+Ni_FeS2 !+Ni_biota+Ni_POM+Ni_DOM

        !NiS  formation/dissollution (REF1) is calculated for Ni_free excuding aborbed Ni 
        Om_NiS=H2S*Ni_free/(self%K_NiS)
        !% FeS formation Fe2+ + HS- -> FeS + H+ (Bektursunova, 11)
        NiS_form=self%K_NiS_form*max(0._rk,(Om_NiS-1._rk)) !*1.e-20
        !% FeS dissollution FeS + H+ -> Fe2+ + HS (Bektursunova, 11)
        NiS_diss=self%K_NiS_diss*NiS*max(0._rk,(1._rk-Om_NiS)) !*1.e-20
    
        _SET_ODE_(self%id_Ni,NiS_diss-NiS_form)  ! "Ni" includes Ni free and Ni adsorped
        _SET_ODE_(self%id_Ni_biota,0.0_rk)
        _SET_ODE_(self%id_Ni_DOM,0.0_rk)
        !_SET_ODE_(self%id_Ni_POM,0.0_rk)
        _SET_ODE_(self%id_NiS,-NiS_diss+NiS_form)
        _SET_ODE_(self%id_H2S,NiS_diss-NiS_form)

        _SET_DIAGNOSTIC_(self%id_Ni_Mn4,Ni_Mn4)
        _SET_DIAGNOSTIC_(self%id_Ni_FeS,Ni_FeS)
        _SET_DIAGNOSTIC_(self%id_Ni_FeS2,Ni_FeS2)
        _SET_DIAGNOSTIC_(self%id_Ni_free,Ni_free)
        _SET_DIAGNOSTIC_(self%id_Ni_tot,Ni_tot)
        _SET_DIAGNOSTIC_(self%id_Ni_tot_diss,Ni+Ni_DOM)
        _SET_DIAGNOSTIC_(self%id_NiS_form,NiS_form)
        _SET_DIAGNOSTIC_(self%id_NiS_diss,NiS_diss)

        _LOOP_END_

    end subroutine do
end module