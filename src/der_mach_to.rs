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
    if mach.is_zero() {
        return F::neg_infinity();
    }
    let half = F::from(0.5).unwrap();
    let t0_t = F::one() + half * (gamma - F::one()) * mach.powi(2);
    (F::from(2.0).unwrap() / (gamma + F::one()) * t0_t)
        .powf(half * (gamma + F::one()) / (gamma - F::one()))
        * (-mach.powi(-2) + half * (gamma + F::one()) * t0_t.powi(-1))
}

pub fn der_mach_to_v_cpt0<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    (gamma - F::one()).sqrt()
        * (F::one() + half * (gamma - F::one()) * mach.powi(2))
            .sqrt()
            .powi(-3)
}

pub fn der_mach_to_mcpt0_ap0<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    gamma / (gamma - F::one()).sqrt()
        * (F::one()
            - (half * (gamma + F::one()) * mach.powi(2))
                / (F::one() + half * (gamma - F::one()) * mach.powi(2)))
        * (F::one() + half * (gamma - F::one()) * mach.powi(2))
            .powf(-half * (gamma + F::one()) / (gamma - F::one()))
}

pub fn der_mach_to_mcpt0_ap<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    gamma / (gamma - F::one()).sqrt() * (F::one() + (gamma - F::one()) * mach.powi(2))
        / (F::one() + half * (gamma - F::one()) * mach.powi(2)).sqrt()
}

pub fn der_mach_to_f_mcpt0<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    (gamma - F::one()).sqrt() / gamma * (F::one() + gamma * mach.powi(2)) / mach.powi(2)
        * (F::one() + half * (gamma - F::one()) * mach.powi(2))
            .sqrt()
            .powi(-3)
}
