//#![cfg_attr(feature = "simd", feature(portable_simd))]
#![doc = include_str!("../README.md")]
// Conventions for generics used in this crate:
//    T: general type, eg f32, Vector3d, Vector3df32, Quaternion etc
//    R: real number type ie f32 or f64
//    F: filter type, eg Pt1Filter, BiquadFilter etc
//#![doc(html_math_jax_enabled)]
#![no_std]
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
#![deny(clippy::panic)]
#![deny(unused_must_use)]
#![warn(unused_results)]
#![warn(clippy::pedantic)]
#![allow(clippy::inline_always)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::must_use_candidate)]

mod biquad_filter;
mod circular_buffer;
mod filters;
mod median_filter;
mod median_filter5;
mod moving_average_filter;
mod pt_filters;
mod rolling_buffer;
mod slew_filter;

pub use biquad_filter::{
    BiquadFilter, BiquadFilterVector2df32, BiquadFilterVector2df64, BiquadFilterVector3df32, BiquadFilterVector3df64,
    BiquadFilterf32, BiquadFilterf64,
};

pub use circular_buffer::CircularBuffer;

pub use filters::{SignalFilter, UpdateFilter};

pub use median_filter::{MedianFilter3, MedianFilter3f32, MedianFilter3f64};
pub use median_filter5::{MedianFilter5, MedianFilter5f32, MedianFilter5f64};

pub use moving_average_filter::{MovingAverageFilter, MovingAverageFilterf32, MovingAverageFilterf64};
pub use moving_average_filter::{
    MovingAverageFilterVector2df32, MovingAverageFilterVector2df64, MovingAverageFilterVector3df32,
    MovingAverageFilterVector3df64, MovingAverageFilterVector4df32, MovingAverageFilterVector4df64,
};

pub use pt_filters::{Pt1Filter, Pt1Filterf32, Pt1Filterf64};
pub use pt_filters::{Pt1FilterVector2df32, Pt1FilterVector2df64, Pt1FilterVector3df32, Pt1FilterVector3df64};

pub use pt_filters::{Pt2Filter, Pt2Filterf32, Pt2Filterf64};
pub use pt_filters::{Pt2FilterVector2df32, Pt2FilterVector2df64, Pt2FilterVector3df32, Pt2FilterVector3df64};

pub use pt_filters::{Pt3Filter, Pt3Filterf32, Pt3Filterf64};
pub use pt_filters::{Pt3FilterVector2df32, Pt3FilterVector2df64, Pt3FilterVector3df32, Pt3FilterVector3df64};

pub use rolling_buffer::RollingBuffer;

pub use slew_filter::{LimitSlew, SlewRateLimiter, SlewRateLimiterf32, SlewRateLimiterf64};
