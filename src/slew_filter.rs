use core::ops::{AddAssign, Mul, Neg, Sub};
use num_traits::{ConstZero, Zero};

/// `SlewRateLimiter` for `f32`<br>
pub type SlewRateLimiterf32 = SlewRateLimiter<f32>;
/// `SlewRateLimiter` for `f64`<br><br>
pub type SlewRateLimiterf64 = SlewRateLimiter<f64>;

#[allow(clippy::doc_paragraphs_missing_punctuation)]
/// Limits the maximum rate of change (dV/dt) of a signal.<br>
/// It allows for different rates depending on whether the signal is
/// increasing (rising) or decreasing (falling).
///
/// The algorithm calculates the difference `delta = x{n} - y{n-1}` and
/// clamps it based on the elapsed time `delta_t`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlewRateLimiter<T> {
    last_output: T,
    rise_step: T, // rise_rate_per_second * dt
    fall_step: T, // fall_rate_per_second * dt
}

impl<T> Default for SlewRateLimiter<T>
where
    T: Copy + ConstZero + Mul<T, Output = T>,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> SlewRateLimiter<T>
where
    T: Copy + ConstZero + Mul<T, Output = T>,
{
    pub fn with_rates(rise_rate_per_second: T, fall_rate_per_second: T, dt: T) -> Self {
        Self { last_output: T::ZERO, rise_step: rise_rate_per_second * dt, fall_step: fall_rate_per_second * dt }
    }

    pub const fn new() -> Self {
        Self { last_output: T::ZERO, rise_step: T::ZERO, fall_step: T::ZERO }
    }
}

impl<T> SlewRateLimiter<T>
where
    T: Copy + Zero + Mul<T, Output = T>,
{
    pub fn reset(&mut self) {
        self.last_output = T::zero();
    }
}

impl<T> SlewRateLimiter<T>
where
    T: Copy + Zero + PartialOrd + Neg<Output = T> + Sub<T, Output = T> + AddAssign,
{
    pub fn update(&mut self, input: T) -> T {
        let diff = input - self.last_output;

        // Select pre-calculated limit
        let limit = if diff > T::zero() { self.rise_step } else { self.fall_step };

        // Clamp using min/max logic (often compiles to branchless CMOV or MIN/MAX instructions)
        let actual_change = if diff > limit {
            limit
        } else if diff < -limit {
            -limit
        } else {
            diff
        };

        self.last_output += actual_change;
        self.last_output
    }
}

/// Adds `value.limit_slew_using(&mut slew_rate_limiter)` method call to `SlewRateLimiter`.
pub trait LimitSlew<T> {
    fn limit_slew_using(&mut self, limiter: &mut SlewRateLimiter<T>);
}

impl<T> LimitSlew<T> for T
where
    T: Copy + Zero + PartialOrd + Neg<Output = T> + Sub<T, Output = T> + AddAssign,
{
    fn limit_slew_using(&mut self, limiter: &mut SlewRateLimiter<T>) {
        *self = limiter.update(*self);
    }
}

#[cfg(any(debug_assertions, test))]
mod tests {
    #![allow(clippy::float_cmp)]
    //use core::default;

    //use crate::UpdateFilter;

    #[allow(unused)]
    use super::*;

    #[allow(unused)]
    fn is_normal<T: Sized + Send + Sync + Unpin>() {}
    #[allow(unused)]
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + PartialEq>() {}

    #[test]
    fn normal_types() {
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
        let mut limiter = SlewRateLimiterf32::with_rates(rise_rate_per_second, fall_rate_per_second, dt);

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
        let mut limiter = SlewRateLimiterf32::with_rates(rise_rate_per_second, fall_rate_per_second, dt);

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
        let mut limiter = SlewRateLimiterf32::with_rates(rise_rate_per_second, fall_rate_per_second, dt);

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
        let mut limiter = SlewRateLimiterf32::with_rates(rise_rate_per_second, fall_rate_per_second, dt);

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
        use crate::SlewRateLimiterf32;
        const MAX_MOTOR_COUNT: usize = 8;
        type MotorOutputs = [f32; MAX_MOTOR_COUNT];
        const QUAD_MOTOR_COUNT: usize = 4;
        type QuadOutputs = [f32; QUAD_MOTOR_COUNT];

        let mut motor_outputs = MotorOutputs::default();

        let quad_outputs = QuadOutputs::default();

        let mut output_filters = <[SlewRateLimiterf32; QUAD_MOTOR_COUNT]>::default();

        for ii in 0..QUAD_MOTOR_COUNT {
            // 1. Take raw value from quad_outputs
            // 2. Update it using the corresponding filter
            // 3. Store result directly in the motor_outputs array
            motor_outputs[ii] = output_filters[ii].update(quad_outputs[ii]);
        }
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
