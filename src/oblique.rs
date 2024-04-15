//! Weak oblique shock functions

use eqsolver::single_variable::FDNewton;
use num::Float;

/// Wave angle for weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_beta;
///
/// assert_eq!(oblique_beta(2.0_f32, 1.4_f32, 0.1745329_f32), 0.6861576);
/// assert_eq!(oblique_beta(5.0_f64, 1.4_f64, 0.3490659_f64), 0.5201241529003784);
///
///
/// ```
pub fn oblique_beta<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let beta_max: F = oblique_beta_max(mach, gamma);
    let mut x0 = beta_max;
    let two = F::from(2.0).unwrap();

    let f = |x: F| {
        theta.tan()
            - two / x.tan() * (mach.powi(2) * x.sin().powi(2) - F::one())
                / (mach.powi(2) * (gamma + (two * x).cos()) + two)
    };
    let beta_result = FDNewton::new(f).solve(x0);
    let mut beta = match beta_result {
        Ok(x) => x,
        Err(_) => return F::nan(),
    };

    while (beta > beta_max) || (beta < F::zero()) {
        x0 = x0 - F::from(0.1).unwrap();
        let beta_result = FDNewton::new(f).solve(x0);
        beta = match beta_result {
            Ok(x) => x,
            Err(_) => return F::nan(),
        };
    }

    beta
}

/// Maximum oblique shock angle
pub fn oblique_beta_max<F: Float>(mach: F, gamma: F) -> F {
    ((F::one() / (gamma * mach.powi(2))
        * (((gamma + F::one()) / F::from(4.0).unwrap() * mach.powi(2)) - F::one()
            + ((gamma + F::one()) * (gamma + F::one()) / F::from(16.0).unwrap() * mach.powi(4)
                + (gamma - F::one()) / F::from(2.0).unwrap() * mach.powi(2)
                + F::one())
            .sqrt()))
    .sqrt())
    .asin()
}

/// Mach number after weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_mach2;
///
/// assert_eq!(oblique_mach2(2.0_f32, 1.4_f32, 0.1745329_f32), 1.6405221);
/// assert_eq!(oblique_mach2(5.0_f64, 1.4_f64, 0.3490659_f64), 3.022151379742916);
///
/// ```
pub fn oblique_mach2<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let beta = oblique_beta(mach, gamma, theta);
    let two = F::from(2.).unwrap();

    ((F::one() + (gamma - F::one()) / two * mach.powi(2))
        / (gamma * mach.powi(2) * beta.sin().powi(2) - (gamma - F::one()) / two)
        + (mach.powi(2) * beta.cos().powi(2))
            / (F::one() + (gamma - F::one()) / two * mach.powi(2) * beta.sin().powi(2)))
    .sqrt()
}

/// Stagnation pressure ratio across weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_p02_p01;
///
/// assert_eq!(oblique_p02_p01(2.0_f32, 1.4_f32, 0.1745329_f32), 0.98464406);
/// assert_eq!(oblique_p02_p01(5.0_f64, 1.4_f64, 0.3490659_f64), 0.5050701357775494);
///
/// ```
pub fn oblique_p02_p01<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let beta = oblique_beta(mach, gamma, theta);
    let mach1n = mach * beta.sin();
    let two = F::from(2.).unwrap();

    F::one()
        / ((two * gamma / (gamma + F::one()) * mach1n.powi(2)
            - (gamma - F::one()) / (gamma + F::one()))
        .powf(F::one() / (gamma - F::one()))
            * (two / (gamma + F::one()) / mach1n.powi(2) + (gamma - F::one()) / (gamma + F::one()))
                .powf(gamma / (gamma - F::one())))
}

/// Static pressure ratio across weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_p2_p1;
///
/// assert_eq!(oblique_p2_p1(2.0_f32, 1.4_f32, 0.1745329_f32), 1.7065787);
/// assert_eq!(oblique_p2_p1(5.0_f64, 1.4_f64, 0.3490659_f64), 7.037411017501249);
///
/// ```
pub fn oblique_p2_p1<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let beta = oblique_beta(mach, gamma, theta);
    let mach1n = mach * beta.sin();
    F::from(2.).unwrap() * gamma / (gamma + F::one()) * (mach1n.powi(2) - F::one()) + F::one()
}

/// Static density ratio across weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_rho2_rho1;
///
/// assert_eq!(oblique_rho2_rho1(2.0_f32, 1.4_f32, 0.1745329_f32), 1.4584259);
/// assert_eq!(oblique_rho2_rho1(5.0_f64, 1.4_f64, 0.3490659_f64), 3.3154179190165545);
///
/// ```
pub fn oblique_rho2_rho1<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let beta = oblique_beta(mach, gamma, theta);
    let mach1n = mach * beta.sin();
    (gamma + F::one()) * mach1n.powi(2)
        / ((gamma - F::one()) * mach1n.powi(2) + F::from(2.).unwrap())
}

/// Static temperature ratio across weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_t2_t1;
///
/// assert_eq!(oblique_t2_t1(2.0_f32, 1.4_f32, 0.1745329_f32), 1.17015128);
/// assert_eq!(oblique_t2_t1(5.0_f64, 1.4_f64, 0.3490659_f64), 2.122631652901466);
///
/// ```
pub fn oblique_t2_t1<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let two = F::from(2.).unwrap();
    let beta = oblique_beta(mach, gamma, theta);
    let mach1n = mach * beta.sin();
    (two + (gamma - F::one()) * mach1n.powi(2))
        * (two * gamma * mach1n.powi(2) - (gamma - F::one()))
        / ((gamma + F::one()).powi(2) * mach1n.powi(2))
}

/// Speed of sound ratio across weak oblique shock
///
/// # Examples
///
/// ```
/// use comp_flow::oblique_a2_a1;
///
/// assert_eq!(oblique_a2_a1(2.0_f32, 1.4_f32, 0.1745329_f32), 1.08173530);
/// assert_eq!(oblique_a2_a1(5.0_f64, 1.4_f64, 0.3490659_f64), 1.4569254108915342);
///
/// ```
pub fn oblique_a2_a1<F: Float>(mach: F, gamma: F, theta: F) -> F {
    let two = F::from(2.).unwrap();
    let beta = oblique_beta(mach, gamma, theta);
    let mach1n = mach * beta.sin();
    ((two + (gamma - F::one()) * mach1n.powi(2))
        * (two * gamma * mach1n.powi(2) - (gamma - F::one()))
        / ((gamma + F::one()).powi(2) * mach1n.powi(2)))
    .sqrt()
}
