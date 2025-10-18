use num::Float;

pub fn der_normal_mach2<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    let t0_t = F::one() + half * (gamma - F::one()) * mach.powi(2);
    let a = (gamma + F::one()).powi(2) * mach * half.sqrt();
    let c = gamma * (F::from(2.0).unwrap() * mach.powi(2) - F::one()) + F::one();
    -a * t0_t.sqrt().powi(-1) * c.sqrt().powi(-3)
}

pub fn der_normal_p02_p01<F: Float>(mach: F, gamma: F) -> F {
    let half = F::from(0.5).unwrap();
    let t0_t = F::one() + half * (gamma - F::one()) * mach.powi(2);
    let a = gamma * mach * (mach.powi(2) - F::one()).powi(2) / t0_t.powi(2);
    let b = (gamma + F::one()) * mach.powi(2) / t0_t * half;
    let c = F::from(2.0).unwrap() * gamma / (gamma + F::one()) * mach.powi(2)
        - (gamma - F::one()) / (gamma + F::one());
    -a * b.powf(F::one() / (gamma - F::one())) * c.powf(-gamma / (gamma - F::one()))
}
