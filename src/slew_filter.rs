//#![allow(unused)]
//use core::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign};
//use core::default::Default;
use core::ops::{Add, Mul, Neg, Sub};
use num_traits::Zero;

pub type SlewRateLimiterf32 = SlewRateLimiter<f32>;
pub type SlewRateLimiterf64 = SlewRateLimiter<f64>;

pub type SymmetricSlewRateLimiterf32 = SymmetricSlewRateLimiter<f32>;
pub type SymmetricSlewRateLimiterf64 = SymmetricSlewRateLimiter<f64>;

/// Slew rate limiter.
///
/// Useful for motor control: we often want a slow ramp-up to protect the gears and reduce current in-rush,
/// and a fast ramp-down (or "emergency" deceleration) for safety and responsiveness.
///
/// In  control loop, a Slew-Rate Limiter ensures that a motor doesn't try to change speed faster than it physically can.
/// A low-pass filter (eg a Pt1Filter) "rounds off" sharp edges, a slew-rate limiter "slopes" them at a constant, linear rate.
//
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlewRateLimiter<T> {
    last_output: T,
    rise_rate_per_second: T,
    fall_rate_per_second: T,
    dt: T,
}

impl<T> Default for SlewRateLimiter<T>
where
    T: Copy + Zero,
{
    fn default() -> Self {
        Self::new(T::zero(), T::zero(), T::zero())
    }
}

impl<T> SlewRateLimiter<T>
where
    T: Copy + Zero,
{
    pub fn new(rise_rate_per_second: T, fall_rate_per_second: T, dt: T) -> Self {
        Self { last_output: T::zero(), rise_rate_per_second, fall_rate_per_second, dt }
    }
    pub fn reset(&mut self) {
        self.last_output = T::zero();
    }
}

impl<T> SlewRateLimiter<T>
where
    T: Copy + PartialOrd + Zero + Neg<Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>,
{
    pub fn update(&mut self, input: T) -> T {
        let diff = input - self.last_output;

        // Choose the rate based on whether the signal is rising or falling
        let max_change =
            if diff > T::zero() { self.rise_rate_per_second * self.dt } else { self.fall_rate_per_second * self.dt };

        // Clamp the change (fall_rate_per_second is used as a magnitude, so we clamp -max to +max)
        let actual_change: T = if diff < -max_change {
            -max_change
        } else if diff > max_change {
            max_change
        } else {
            diff
        };

        self.last_output = self.last_output + actual_change;
        self.last_output
    }
}

/// Extension trait for slew rate limiter.
pub trait LimitSlew<T> {
    fn limit_slew_using(&mut self, limiter: &mut SlewRateLimiter<T>) -> &mut Self;
}

impl<T> LimitSlew<T> for T
where
    T: Copy + PartialOrd + Zero + Neg<Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>,
{
    fn limit_slew_using(&mut self, limiter: &mut SlewRateLimiter<T>) -> &mut Self {
        *self = limiter.update(*self);
        self
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymmetricSlewRateLimiter<T> {
    last_output: T,
    max_rate_per_second: T, // Units per second (e.g., max change in duty cycle/sec)
    dt: T,
}

impl<T> Default for SymmetricSlewRateLimiter<T>
where
    T: Copy + Zero,
{
    fn default() -> Self {
        Self::new(T::zero(), T::zero())
    }
}

impl<T> SymmetricSlewRateLimiter<T>
where
    T: Copy + Zero,
{
    pub fn new(max_rate_per_second: T, dt: T) -> Self {
        Self { last_output: T::zero(), max_rate_per_second, dt }
    }
}

impl<T> SymmetricSlewRateLimiter<T>
where
    T: Copy + Zero,
{
    pub fn reset(&mut self) {
        self.last_output = T::zero();
    }
}

impl<T> SymmetricSlewRateLimiter<T>
where
    T: Copy + PartialOrd + Neg<Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>,
{
    pub fn update(&mut self, input: T) -> T {
        // Calculate the maximum change allowed for this specific time step
        let max_change = self.max_rate_per_second * self.dt;
        let diff = input - self.last_output;

        // Clamp the difference to the max allowed change
        //let actual_change = diff.clamp(-max_change, max_change);
        let actual_change: T = if diff < -max_change {
            -max_change
        } else if diff > max_change {
            max_change
        } else {
            diff
        };

        self.last_output = self.last_output + actual_change;
        self.last_output
    }
}

// Extension trait for symmetric slew limiter.
pub trait LimitSlewSymmetric<T> {
    fn limit_symmetric_slew_using(&mut self, limiter: &mut SymmetricSlewRateLimiter<T>) -> &mut Self;
}

impl<T> LimitSlewSymmetric<T> for T
where
    T: Copy + PartialOrd + Neg<Output = T> + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>,
{
    fn limit_symmetric_slew_using(&mut self, limiter: &mut SymmetricSlewRateLimiter<T>) -> &mut Self {
        *self = limiter.update(*self);
        self
    }
}

#[cfg(any(debug_assertions, test))]
mod tests {
    #![allow(unused)]
    use crate::ApplyFilter;

    use super::*;

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<SymmetricSlewRateLimiter<f32>>();
        is_full::<SymmetricSlewRateLimiter<f64>>();
        is_full::<SlewRateLimiter<f32>>();
        is_full::<SlewRateLimiter<f64>>();
    }
    #[test]
    fn biquad_filter_f32() {}
    #[test]
    fn slew_rate_asymmetric() {
        // Rise: 10.0 units/sec, Fall: 100.0 units/sec (Fast stop)
        let dt = 0.1;
        let rise_rate_per_second = 10.0;
        let fall_rate_per_second = 100.0;
        let mut limiter = SlewRateLimiterf32::new(rise_rate_per_second, fall_rate_per_second, dt);

        // 1. Test Ramping Up
        // Max rise per step = 10.0 * 0.1 = 1.0 unit
        // Target 10.0, last_output 0.0 -> result should be 1.0
        assert_eq!(1.0, limiter.update(10.0));
        // Next step: 1.0 + 1.0 = 2.0
        assert_eq!(2.0, limiter.update(10.0));

        // 2. Test Ramping Down (Fast)
        // Max fall per step = 100.0 * 0.1 = 10.0 units
        // Current state is 2.0. Target 0.0.
        // Difference is -2.0, which is less than the max fall (10.0).
        // It should reach 0.0 in a single step.
        assert_eq!(0.0, limiter.update(0.0));
    }
    #[test]
    fn slew_reset() {
        let dt = 0.1;
        let rise_rate_per_second = 4.0;
        let fall_rate_per_second = 50.0;
        let mut limiter = SlewRateLimiterf32::new(rise_rate_per_second, fall_rate_per_second, dt);

        let input = 100.0;
        let output = limiter.update(input);
        assert_eq!(0.4, output);

        let input = 100.0;
        let output = limiter.update(input);
        assert_eq!(0.8, output);

        // After reset, last_output is 0.0.
        // With dt=0.1 and rise=100.0, result should be 0.4.
        limiter.reset();
        let input = 100.0;
        let output = limiter.update(input);
        assert_eq!(0.4, output);
    }
    #[test]
    fn slew_clamping() {
        let dt = 0.1;
        let rise_rate_per_second = 4.0;
        let fall_rate_per_second = 50.0;
        let mut limiter = SlewRateLimiterf32::new(rise_rate_per_second, fall_rate_per_second, dt);

        let input = 100.0;
        let output = limiter.update(input);
        assert_eq!(0.4, output);

        let input = 100.0;
        let output = limiter.update(input);
        assert_eq!(0.8, output);

        let input = 1.1;
        let output = limiter.update(input);
        assert_eq!(1.1, output);
    }

    #[test]
    fn extended_trait() {
        let dt = 0.1;
        let rise_rate_per_second = 4.0;
        let fall_rate_per_second = 50.0;
        let mut limiter = SlewRateLimiterf32::new(rise_rate_per_second, fall_rate_per_second, dt);

        let mut value = 100.0;
        value.limit_slew_using(&mut limiter);
        assert_eq!(0.4, value);
        value = 100.0;
        value.limit_slew_using(&mut limiter);
        assert_eq!(0.8, value);

        value = 1.1;
        value.limit_slew_using(&mut limiter);
        assert_eq!(1.1, value);
    }

    #[test]
    fn output_array() {
        /*type Outputs = [f32; 4];
        type OutputSlewRateLimiter = SlewRateLimiter<Outputs, f32>;

        let dt = 0.1;
        let rise_rate_per_second: Outputs = [1.0, 2.0, 3.0, 4.0];
        let fall_rate_per_second: Outputs = [10.0, 20.0, 30.0, 40.0];
        let mut limiter: OutputSlewRateLimiter =
            OutputSlewRateLimiter::new(rise_rate_per_second, fall_rate_per_second, dt);

        let mut outputs: Outputs = [0.0, 0.0, 0.0, 0.0];*/
        //outputs.limit_slew_using(&mut limiter);
        //assert_eq!(0.4, value);
        /*value = 100.0;
        value.limit_slew_using(&mut limiter);
        assert_eq!(0.8, value);

        value = 1.1;
        value.limit_slew_using(&mut limiter);
        assert_eq!(1.1, value);*/
    }
}
