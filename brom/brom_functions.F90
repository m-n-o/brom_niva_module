!-----------------------------------------------------------------------
! brom_functions is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.
!-----------------------------------------------------------------------

module brom_functions
  use fabm_types
contains
  !
  !
  !
  pure real(rk) function n_zero(var)
    real(rk),intent(in):: var

    n_zero = max(var, 1.e-10_rk)
  end function n_zero
  !
  !
  !
  pure real(rk) function quota(var, var2)
    real(rk),intent(in):: var, var2

    quota = var/n_zero(var2)
  end function quota
  !
  ! It is limiter type of function
  ! if concentration equals threshold_value hyper_limiter = 0.5
  ! and growths coef 1 gives hyperbolic function
  ! larger coef gives larger gradients
  !
  pure real(rk) function hyper_limiter(threshold_value, r, coef)
    real(rk), intent(in) :: threshold_value, r, coef

    hyper_limiter = 0.5_rk+0.5_rk*tanh((r-threshold_value)*coef)
  end function hyper_limiter
  !
  ! It is inhibitor type of function
  !
  pure real(rk) function hyper_inhibitor(threshold_value, r, coef)
    ! Threshold value for the reaction
    real(rk), intent(in) :: threshold_value, r, coef

    hyper_inhibitor = 0.5_rk-0.5_rk*tanh((r-threshold_value)*coef)
  end function hyper_inhibitor
  !
  !Monod limiter
  !
  pure real(rk) function monod(ks, r)
      real(rk),intent(in):: ks
      real(rk),intent(in):: r

      monod = r/(r+ks)
  end function monod
  !
  !Monod inhibitor
  !
  pure real(rk) function inhib(ks, r)
    real(rk),intent(in):: ks
    real(rk),intent(in):: r

    inhib = ks/(r+ks)
  end function inhib
  !
  !This is a squared Michaelis-Menten type of limiter
  !
  pure real(rk) function monod_squared(ks,r)
    real(rk),intent(in):: ks,r

    monod_squared = r**2._rk/(ks**2._rk+r**2._rk)
  end function monod_squared
  !
  !This is a squared Michaelis-Menten type of inhibitor
  !
  pure real(rk) function monod_squared_inhibitor(ks,r)
    real(rk),intent(in):: ks,r

    monod_squared_inhibitor = ks**2._rk/(r**2._rk+ks**2._rk)
  end function monod_squared_inhibitor

  pure real(rk) function carbon_g_to_mole(carbon)
    real(rk), intent(in):: carbon

    carbon_g_to_mole = carbon/12.011_rk
  end function carbon_g_to_mole

end module brom_functions
