//! Derivatives of relations as a function of Mach number
#[doc(no_inline)]
use num::Float;

/// Total temperature ratio for given mach number and specific heat ratio
///
pub fn der_mach_to_t0_t<F: Float>(mach: F, gamma: F) -> F {
    (gamma - F::one()) * mach
}

pub fn der_mach_to_p0_p<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    gamma
        * mach
        * (F::one() + half * (gamma - F::one()) * mach.powi(2)).powf(F::one() / (gamma - F::one()))
}

pub fn der_mach_to_rho0_rho<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    mach * (F::one() + half * (gamma - F::one()) * mach.powi(2)).powf(F::one() / (F::one() - gamma))
}

pub fn der_mach_to_a_ac<F: Float>(mach: F, gamma: F) -> F {
    todo!()
}
