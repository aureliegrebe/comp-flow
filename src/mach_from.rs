//! Collection of functions for isentropic compressible flow.

use crate::{
    der_mach_to_f_mcpt0, der_mach_to_mcpt0_ap0, mach_to_a_ac, mach_to_f_mcpt, mach_to_mcpt0_ap0,
    mach_to_pm_angle,
};
use eqsolver::single_variable::{FDNewton, Newton};
use num::Float;

/// Mach number for a given Prandtl-Meyer angle in radians.
///
/// <div class="warning">
///
/// This function uses Newton's method to solve for the Mach number. If this
/// function must be called many times, it may be preferable to make a look up
/// table with `mach_to_pm_angle` and interpolate those values.
///
/// </div>
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

/// Mach number for a given total temperature ratio.
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

/// Mach number for a given total temperature ratio.
///
pub fn mach_from_t0_t<F: Float>(t0_t: F, gamma: F) -> F {
    let two = F::from(2.0).unwrap();
    (two / (gamma - F::one()) * (t0_t - F::one())).sqrt()
}

/// Mach number for a given total pressure ratio.
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

pub fn mach_from_p0_p<F: Float>(p0_p: F, gamma: F) -> F {
    let two = F::from(2.0).unwrap();
    (two / (gamma - F::one()) * (p0_p.powf((gamma - F::one()) / gamma) - F::one())).sqrt()
}

/// Mach number for a given stagnation density ratio.
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

pub fn mach_from_rho0_rho<F: Float>(rho0_rho: F, gamma: F) -> F {
    let two = F::from(2.0).unwrap();
    (two / (gamma - F::one()) * (rho0_rho.powf(gamma - F::one()) - F::one())).sqrt()
}

pub fn mach_from_v_cpt0<F: Float>(v_cpt0: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    (F::one() / (gamma - F::one()) * v_cpt0.powi(2) / (F::one() - half * v_cpt0.powi(2))).sqrt()
}

pub fn mach_from_f_mcpt0<F: Float>(f_mcpt0: F, gamma: F, supersonic: bool) -> F {
    let f = |m| mach_to_f_mcpt(m, gamma);
    let df = |m| der_mach_to_f_mcpt0(m, gamma);
    let x0 = match supersonic {
        true => F::from(1.5).unwrap(),
        false => F::from(0.5).unwrap(),
    };
    Newton::new(f, df).with_itermax(100).solve(x0).unwrap()
}

pub fn mach_from_mcpt0_ap0<F: Float>(f_mcpt0: F, gamma: F, supersonic: bool) -> F {
    let f = |m| mach_to_mcpt0_ap0(m, gamma);
    let df = |m| der_mach_to_mcpt0_ap0(m, gamma);
    let x0 = match supersonic {
        true => F::from(1.5).unwrap(),
        false => F::from(0.5).unwrap(),
    };
    Newton::new(f, df).with_itermax(100).solve(x0).unwrap()
}

pub fn mach_from_mcpt0_ap<F: Float>(f_mcpt0: F, gamma: F) -> F {
    let f = |m| mach_to_mcpt0_ap(m, gamma);
    let df = |m| der_mach_to_mcpt0_ap(m, gamma);
    let x0 = F::from(0.5).unwrap();
    Newton::new(f, df).with_itermax(100).solve(x0).unwrap()
}

/// Mach number for a given critical area ratio.
///
/// <div class="warning">
///
/// This function uses a bisection algorythm to solve for the Mach number. If this
/// function must be called many times, it may be preferable to make a look up
/// table with `mach_to_a_ac` and interpolate those values.
///
/// </div>
///
/// # Examples
///
/// ```
/// use comp_flow::mach_from_a_ac;
///
/// assert_eq!(mach_from_a_ac(1.6875000000000002_f64, 1.4, true), 2.0);
/// assert_eq!(mach_from_a_ac(5.821828750000001_f64, 1.4, false), 0.1);
/// assert_eq!(mach_from_a_ac(6.25, 1.4, false), 0.09307469911759117);
/// assert_eq!(mach_from_a_ac(1.0, 1.4, false), 1.0);
/// assert_eq!(mach_from_a_ac(1.0, 1.4, true), 1.0);
/// ```
pub fn mach_from_a_ac<F: Float>(a_ac: F, gamma: F, supersonic: bool) -> F {
    if a_ac.is_one() {
        return F::one();
    }
    let mut m_max: F;
    let mut m_min: F;
    if supersonic {
        m_max = F::max_value();
        m_min = F::one();
    } else {
        m_max = F::one();
        m_min = F::zero();
    }
    let mut m_try = (m_max + m_min) / F::from(2.).unwrap();
    let mut val: F;
    let mut i: usize = 0;
    let max_iter: usize = 10000;
    while i < max_iter {
        val = mach_to_a_ac(m_try, gamma);
        if val == a_ac {
            return m_try;
        } else if val < a_ac {
            if supersonic {
                m_min = m_try;
            } else {
                m_max = m_try;
            }
        } else if val > a_ac {
            if supersonic {
                m_max = m_try;
            } else {
                m_min = m_try;
            }
        }
        m_try = (m_max + m_min) / F::from(2.).unwrap();
        i += 1;
    }
    return m_try;
}
