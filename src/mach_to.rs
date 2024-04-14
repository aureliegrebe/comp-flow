//! Collection of functions for isentropic compressible flow.

use num::Float;

/// Prandtl-Meyer angle in radians for a given mach number and specific heat ratio.
///
/// # Examples
///
/// ```
/// use comp_flow::mach_to_pm_angle;
///
/// assert_eq!(mach_to_pm_angle(2.0_f32, 1.4_f32), 0.4604136818474);
/// assert_eq!(mach_to_pm_angle(1.0_f64, 1.4_f64), 0.0);
/// ```
pub fn mach_to_pm_angle<F: Float>(mach: F, gamma: F) -> F {
    ((gamma + F::one()) / (gamma - F::one())).sqrt()
        * ((gamma - F::one()) / (gamma + F::one()) * (mach.powi(2) - F::one()))
            .sqrt()
            .atan()
        - (mach.powi(2) - F::one()).sqrt().atan()
}

/// Mach angle in radians for a given mach number.
///
/// # Examples
///
/// ```
/// use comp_flow::mach_to_mach_angle;
///
/// assert_eq!(mach_to_mach_angle(2.0_f32), 0.5235988);
/// assert_eq!(mach_to_mach_angle(1.0_f64), 1.5707963267948966);
/// ```
pub fn mach_to_mach_angle<F: Float>(mach: F) -> F {
    (F::one() / mach).asin()
}

/// Total temperature ratio for given mach number and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_to_t_t0;
///
/// assert_eq!(mach_to_t_t0(0.0, 1.4), 1.0);
/// assert_eq!(mach_to_t_t0(1.0, 1.4), 0.8333333333333334);
/// assert_eq!(mach_to_t_t0(2.0_f32, 1.4), 0.55555556);
/// ```
pub fn mach_to_t_t0<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    (F::one() + half * (gamma - F::one()) * mach.powi(2)).powi(-1)
}

/// Total pressure ratio for given mach number and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_to_p_p0;
///
/// assert_eq!(mach_to_p_p0(0.0, 1.4), 1.0);
/// assert_eq!(mach_to_p_p0(1.0, 1.4), 0.5282817877171742);
/// assert_eq!(mach_to_p_p0(2.0, 1.4), 0.12780452546295096);
/// ```
pub fn mach_to_p_p0<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    (F::one() + half * (gamma - F::one()) * mach.powi(2)).powf((gamma) / (F::one() - gamma))
}

/// Stagnation density ratio for given mach number and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_to_rho_rho0;
///
/// assert_eq!(mach_to_rho_rho0(0.0, 1.4), 1.0);
/// assert_eq!(mach_to_rho_rho0(1.0, 1.4), 0.633938145260609);
/// assert_eq!(mach_to_rho_rho0(2.0, 1.4), 0.2300481458333117);
/// ```
pub fn mach_to_rho_rho0<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    (F::one() + half * (gamma - F::one()) * mach.powi(2)).powf(F::one() / (F::one() - gamma))
}

/// Critical area ratio for given mach number and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_to_a_ac;
///
/// assert_eq!(mach_to_a_ac(0.1, 1.4), 5.821828750000001);
/// assert_eq!(mach_to_a_ac(1.0, 1.4), 1.0);
/// assert_eq!(mach_to_a_ac(2.0, 1.4), 1.6875000000000002);
/// ```
pub fn mach_to_a_ac<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    F::one() / mach
        * ((F::one() + half * (gamma - F::one()) * mach.powi(2)) / (half * (gamma + F::one())))
            .powf(half * (gamma + F::one()) / (gamma - F::one()))
}
