// Conventions for generics used in this crate:
//    T: general type, eg f32, Vector3d, Vector3df32, Quaternion etc
//    R: real number type ie f32 or f64
//    F: filter type, eg Pt1Filter, BiquadFilter etc

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
mod rolling_buffer;
mod slew_filter;

pub use biquad_filter::{BiquadFilter, BiquadFilterState, BiquadFilterf32, BiquadFilterf64};
pub use circular_buffer::CircularBuffer;
pub use filters::{ApplyFilter, FilterSignal};
pub use median_filter::{
    Median3Filter, Median3FilterVector3df32, Median3Filterf32, MedianFilter, MedianFilterVector3df32, MedianFilterf32,
};
pub use moving_average_filter::{MovingAverageFilter, MovingAverageFilterVector3df32, MovingAverageFilterf32};
pub use pt_filters::{
    Pt1Filter, Pt1Filterf32, Pt1Filterf64, Pt2Filter, Pt2Filterf32, Pt2Filterf64, Pt3Filter, Pt3Filterf32, Pt3Filterf64,
};
pub use rolling_buffer::RollingBuffer;

pub use slew_filter::{
    LimitSlew, LimitSlewSymmetric, SlewRateLimiter, SlewRateLimiterf32, SlewRateLimiterf64,
    SymmetricSlewRateLimiterf32, SymmetricSlewRateLimiterf64,
};
