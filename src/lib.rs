//! # Compressible Flow
//!
//! `comp-flow` is a collection of functions for basic compressible flow relations.
//!
//! <div class="warning">
//!
//! The included functions have no input checking or error handling whatsoever.
//! Invalid (non-physical) inputs such as mach < 1 for a shock relation or gamma < 1
//! may produce non-sensical outputs.
//!
//! </div>
//!
#![warn(missing_docs)]

pub mod der_mach_to;
pub mod der_normal;
pub mod mach_from;
pub mod mach_to;
pub mod normal;
pub mod oblique;

#[doc(inline)]
pub use der_mach_to::*;
#[doc(inline)]
pub use der_normal::*;
#[doc(inline)]
pub use mach_from::*;
#[doc(inline)]
pub use mach_to::*;
#[doc(inline)]
pub use normal::*;
#[doc(inline)]
pub use oblique::*;
