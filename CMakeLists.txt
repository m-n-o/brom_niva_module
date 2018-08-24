add_library(fabm_models_niva OBJECT
            niva_model_library.F90
            brom/brom_bio.F90
            brom/brom_carbon.F90
            brom/brom_eq_constants.F90
            brom/brom_main_nutrients.F90
            brom/brom_nitrogen.F90
            brom/brom_pH.F90
            brom/brom_functions.F90
           )

add_dependencies(fabm_models_niva fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
