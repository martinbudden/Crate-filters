# filters Rust Crate ![license](https://img.shields.io/badge/license-MIT-green) [![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) ![open source](https://badgen.net/badge/open/source/blue?icon=github)

This crate contains a collection of filters and utilities for digital signal processing.

This crate is `no_std`, that it does not link to the standard library and so does not depend on an operating system
and uses no allocation. This means it is suitable for embedded system.

## Overview

This crate contains the following filters:

`Pt1Filter`, `Pt2Filter`, and `Pt3Filter` : basic first, second, and third order low-pass filters.

`BiquadFilter`: second order filter that can be used as a low-pass, high-pass, or notch filter.

`MedianFilter3` and `MedianFilter5`: median filters for spike rejection.

`SlewRateLimiter`: filter that limits the rate of change of a signal.

`MovingAverageFilter` : simple moving average filter.

The filters have aliases for `f32` and `f64` versions:

| `f32`                   | `f64`                   |
| ----------------------- | ------------------------|
| `Pt1Filterf32`          | `Pt1Filterf64`          |
| `Pt2Filterf32`          | `Pt2Filterf64`          |
| `BiquadFilterf32`       | `BiquadFilterf64`       |
| `MedianFilter3f32`      | `MedianFilter3f64`      |
| `MedianFilter5f32`      | `MedianFilter5f64`      |
| `SlewRateLimiterf32`    | `SlewRateLimiterf64`    |
| `MovingAverageFilterf32`| `MovingAverageFilterf64`|

Additionally the all the filters except the Median Filters and Slew Rate Limiters have aliases for their vectorized forms for both `f32` and `f64` vectors.

So for `Pt1Filter` we have:

| vector `f32`            | vector `f64`          |
| ----------------------- | ----------------------|
| `Pt1FilterVector2df32`  |`Pt1FilterVector2df64` |
| `Pt1FilterVector3df32`  |`Pt1FilterVector3df64` |
| `Pt1FilterVector4df32`  |`Pt1FilterVector4df64` |

and similarly for the other filters.

## Rolling and Circular buffers

This crate also has basic circular and rolling buffers:
`CircularBuffer<T, const N: usize>` and `RollingBuffer<T, const N: usize>`.

## Motivation

They have been developed for use in stabilized vehicles (self balancing robots and aircraft)and have been used to:

1. Filter gyro and accelerometer output for use in the Attitude and Heading Reference System (AHRS).
2. Filter motor encoder values for use in the motor controller.
3. Filter derivative terms in a PID controller.
4. Filter motor power input values to smooth the motor speed.

## Original implementation

I originally implemented this crate as a C++ library:
[Library-Filter](https://github.com/martinbudden/Library-Filter).

## License

Licensed under either of:

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.
