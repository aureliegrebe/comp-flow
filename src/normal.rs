//! Normal Shock relations
use num::Float;

/// Mach number after normal shock
///
/// # Examples
///
/// ```
/// use comp_flow::normal_mach2;
///
/// assert_eq!(normal_mach2(2.0_f32, 1.4_f32), 0.57735026);
/// assert_eq!(normal_mach2(5.0_f64, 1.4_f64), 0.41522739926869984);
///
/// ```
pub fn normal_mach2<F: Float>(mach: F, gamma: F) -> F {
    let two = F::from(2.).unwrap();
    ((F::one() + (gamma - F::one()) / two * mach.powi(2))
        / (gamma * mach.powi(2) - (gamma - F::one()) / two))
        .sqrt()
}

/// Total pressure ratio across normal shock
///
/// # Examples
///
/// ```
/// use comp_flow::normal_p02_p01;
///
/// assert_eq!(normal_p02_p01(2.0_f32, 1.4_f32), 0.7208737);
/// assert_eq!(normal_p02_p01(5.0_f64, 1.4_f64), 0.061716319748617694);
///
/// ```
pub fn normal_p02_p01<F: Float>(mach: F, gamma: F) -> F {
    let two = F::from(2.).unwrap();
    F::one()
        / ((two * gamma / (gamma + F::one()) * mach.powi(2)
            - (gamma - F::one()) / (gamma + F::one()))
        .powf(F::one() / (gamma - F::one()))
            * (two / (gamma + F::one()) / mach.powi(2) + (gamma - F::one()) / (gamma + F::one()))
                .powf(gamma / (gamma - F::one())))
}

/// Static pressure ratio across normal shock
///
/// # Examples
///
/// ```
/// use comp_flow::normal_p2_p1;
///
/// assert_eq!(normal_p2_p1(2.0_f32, 1.4_f32), 4.5);
/// assert_eq!(normal_p2_p1(5.0_f64, 1.4_f64), 29.0);
///
/// ```
pub fn normal_p2_p1<F: Float>(mach: F, gamma: F) -> F {
    F::from(2.).unwrap() * gamma / (gamma + F::one()) * (mach.powi(2) - F::one()) + F::one()
}

/// Static density ratio across normal shock
///
/// # Examples
///
/// ```
/// use comp_flow::normal_rho2_rho1;
///
/// assert_eq!(normal_rho2_rho1(2.0_f32, 1.4_f32), 2.66666666);
/// assert_eq!(normal_rho2_rho1(5.0_f64, 1.4_f64), 5.000000000000001);
///
/// ```
pub fn normal_rho2_rho1<F: Float>(mach: F, gamma: F) -> F {
    (gamma + F::one()) * mach.powi(2) / ((gamma - F::one()) * mach.powi(2) + F::from(2.).unwrap())
}

/// Static temperature ratio across normal shock
///
/// # Examples
///
/// ```
/// use comp_flow::normal_t2_t1;
///
/// assert_eq!(normal_t2_t1(2.0_f32, 1.4_f32), 1.6875);
/// assert_eq!(normal_t2_t1(5.0_f64, 1.4_f64), 5.799999999999999);
///
/// ```
pub fn normal_t2_t1<F: Float>(mach: F, gamma: F) -> F {
    let two = F::from(2.).unwrap();
    (two + (gamma - F::one()) * mach.powi(2)) * (two * gamma * mach.powi(2) - (gamma - F::one()))
        / ((gamma + F::one()).powi(2) * mach.powi(2))
}

/// Speed of sound ratio across normal shock
///
/// # Examples
///
/// ```
/// use comp_flow::normal_a2_a1;
///
/// assert_eq!(normal_a2_a1(2.0_f32, 1.4_f32), 1.29903810);
/// assert_eq!(normal_a2_a1(5.0_f64, 1.4_f64), 2.408318915758459);
///
/// ```
pub fn normal_a2_a1<F: Float>(mach: F, gamma: F) -> F {
    let two = F::from(2.).unwrap();
    ((two + (gamma - F::one()) * mach.powi(2)) * (two * gamma * mach.powi(2) - (gamma - F::one()))
        / ((gamma + F::one()).powi(2) * mach.powi(2)))
    .sqrt()
}
