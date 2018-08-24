module niva_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: niva_model_factory

contains

   subroutine create(self,name,model)

      use fabm_niva_brom_bio
      use fabm_niva_brom_carbon
      use fabm_niva_brom_eq_constants
      use fabm_niva_brom_main_nutrients
      use fabm_niva_brom_nitrogen
      use fabm_niva_brom_pH
      ! Add new NIVA models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('brom_bio');           allocate(type_niva_brom_bio::model)
         case ('brom_carbon');        allocate(type_niva_brom_carbon::model)
         case ('brom_eq_constants');  allocate(type_niva_brom_eq_constants::model)
         case ('brom_main_nutrients');allocate(type_niva_brom_main_nutrients::model)
         case ('brom_nitrogen');      allocate(type_niva_brom_nitrogen::model)
         case ('brom_pH');            allocate(type_niva_brom_pH::model)
         ! Add new NIVA models here
      end select

   end subroutine

end module
