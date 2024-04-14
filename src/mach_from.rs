//! Collection of fuctions for isentropic compressible flow.

use crate::{mach_to_a_ac, mach_to_pm_angle};
use eqsolver::single_variable::FDNewton;
use num::Float;

/// Prandtl-Meyer angle in radians for a given mach number and specific heat ratio.
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_pm_angle;
///
/// assert_eq!(mach_from_pm_angle(0.4604136818474_f32, 1.4_f32), 2.0);
/// assert_eq!(mach_from_pm_angle(0.0_f64, 1.4_f64),  1.00000022981460310);
/// ```
pub fn mach_from_pm_angle<F: Float>(pm_angle: F, gamma: F) -> F {
    let f = |m| mach_to_pm_angle(m, gamma) - pm_angle;
    let x0 = F::from(2.).unwrap();
    FDNewton::new(f).solve(x0).unwrap()
}

/// Mach number for a given mach angle in radians.
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_mach_angle;
///
/// assert_eq!(mach_from_mach_angle(0.5235988_f32), 2.0);
/// assert_eq!(mach_from_mach_angle(1.5707963267948966_f64), 1.0);
/// ```
pub fn mach_from_mach_angle<F: Float>(mach_angle: F) -> F {
    // TODO check for invalid input i.e. mach_angle > 90 deg
    (F::one()) / mach_angle.sin()
}

/// Total temperature ratio for given mach number and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_t_t0;
///
/// assert_eq!(mach_from_t_t0(1.0, 1.4), 0.0);
/// assert_eq!(mach_from_t_t0(0.8333333333333334, 1.4), 1.0);
/// assert_eq!(mach_from_t_t0(0.55555556_f32, 1.4), 2.0);
/// ```
pub fn mach_from_t_t0<F: Float>(t_t0: F, gamma: F) -> F {
    let two = F::from(2.0).unwrap();
    (two / (gamma - F::one()) * (F::one() / t_t0 - F::one())).sqrt()
}

/// Mach number for given total pressure ratio and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_p_p0;
///
/// assert_eq!(mach_from_p_p0(1.0, 1.4), 0.0);
/// assert_eq!(mach_from_p_p0(0.5282817877171742, 1.4), 1.0);
/// assert_eq!(mach_from_p_p0(0.1278045254629509, 1.4), 2.0);
/// ```
pub fn mach_from_p_p0<F: Float>(p_p0: F, gamma: F) -> F {
    let two = F::from(2.0).unwrap();
    (two / (gamma - F::one()) * (p_p0.powf((F::one() - gamma) / gamma) - F::one())).sqrt()
}

/// Mach number for given stagnation density ratio and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_rho_rho0;
///
/// assert_eq!(mach_from_rho_rho0(1.0, 1.4), 0.0);
/// assert_eq!(mach_from_rho_rho0(0.633938145260609, 1.4), 1.0);
/// assert_eq!(mach_from_rho_rho0(0.2300481458333117, 1.4), 2.0);
/// ```
pub fn mach_from_rho_rho0<F: Float>(rho_rho0: F, gamma: F) -> F {
    let two = F::from(2.0).unwrap();
    (two / (gamma - F::one()) * (rho_rho0.powf(F::one() - gamma) - F::one())).sqrt()
}

/// Critical area ratio for given mach number and specific heat ratio
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_a_ac;
///
/// assert_eq!(mach_from_a_ac(5.821828750000001, 1.4, false), 0.1);
/// assert_eq!(mach_from_a_ac(1.0, 1.4, false), 1.0);
/// assert_eq!(mach_from_a_ac(1.0, 1.4, true), 1.0);
/// assert_eq!(mach_from_a_ac(1.6875000000000002, 1.4, true), 2.0);
/// ```
pub fn mach_from_a_ac<F: Float>(a_ac: F, gamma: F, supersonic: bool) -> F {
    if a_ac.is_one() {
        return F::one();
    }
    let f = |m| mach_to_a_ac(m, gamma) - a_ac;
    let x0: F;
    if supersonic {
        x0 = F::from(1.01).unwrap();
    } else {
        x0 = F::from(0.99).unwrap();
    }
    FDNewton::new(f).solve(x0).unwrap()
}
