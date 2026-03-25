#![no_std]
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
#![deny(clippy::panic)]
#![deny(unused_must_use)]

mod circular_buffer;
mod filters;
mod rolling_buffer;

pub use circular_buffer::CircularBuffer;
pub use filters::{
    BiquadFilter, BiquadFilterState, BiquadFilterf32, BiquadFilterf64, FilterSignal, MovingAverageFilter, Pt1Filter,
    Pt1Filterf32, Pt1Filterf64, Pt2Filter, Pt2Filterf32, Pt2Filterf64, Pt3Filter, Pt3Filterf32, Pt3Filterf64,
};
pub use rolling_buffer::RollingBuffer;
