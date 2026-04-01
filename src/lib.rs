//#![feature(portable_simd)]
//use core::simd::f32x4; // A vector of four 32-bit floats#![doc = include_str!("../README.md")]
//#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))]

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

mod biquad_filter;
mod circular_buffer;
mod filters;
mod median_filter;
mod moving_average_filter;
mod pt_filters;
mod pt_filters_simd;
mod rolling_buffer;
mod slew_filter;

pub use biquad_filter::{BiquadFilter, BiquadFilterf32, BiquadFilterf64, BiquadFilterVector2df32, BiquadFilterVector2df64, BiquadFilterVector3df32, BiquadFilterVector3df64};

pub use circular_buffer::CircularBuffer;

pub use filters::{SignalFilter, UpdateFilter};

pub use median_filter::{Median3Filter, Median3Filterf32, MedianFilter, MedianFilterf32};
pub use median_filter::{Median3FilterVector3df32, MedianFilterVector3df32};

pub use moving_average_filter::MovingAverageFilterVector3df32;
pub use moving_average_filter::{MovingAverageFilter, MovingAverageFilterf32};

pub use pt_filters::{
    Pt1Filter, Pt1FilterVector2df32, Pt1FilterVector2df64, Pt1FilterVector3df32, Pt1FilterVector3df64, Pt1Filterf32,
    Pt1Filterf64,
};
pub use pt_filters::{
    Pt2Filter, Pt2FilterVector2df32, Pt2FilterVector2df64, Pt2FilterVector3df32, Pt2FilterVector3df64, Pt2Filterf32,
    Pt2Filterf64,
};
pub use pt_filters::{
    Pt3Filter, Pt3FilterVector2df32, Pt3FilterVector2df64, Pt3FilterVector3df32, Pt3FilterVector3df64, Pt3Filterf32,
    Pt3Filterf64,
};

pub use rolling_buffer::RollingBuffer;

pub use slew_filter::{LimitSlew, SlewRateLimiter, SlewRateLimiterf32, SlewRateLimiterf64};
