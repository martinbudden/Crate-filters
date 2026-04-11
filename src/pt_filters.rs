use core::ops::{Add, AddAssign, Div, Mul, Sub};
use num_traits::{MulAdd, One, Zero};
use vector_quaternion_matrix::{MathConstants, Vector2d, Vector3d};

use crate::SignalFilter;

pub type Pt1Filterf32 = Pt1Filter<f32, f32>;
pub type Pt1FilterVector2df32 = Pt1Filter<Vector2d<f32>, f32>;
pub type Pt1FilterVector2df64 = Pt1Filter<Vector2d<f64>, f64>;

pub type Pt1Filterf64 = Pt1Filter<f64, f64>;
pub type Pt1FilterVector3df32 = Pt1Filter<Vector3d<f32>, f32>;
pub type Pt1FilterVector3df64 = Pt1Filter<Vector3d<f64>, f64>;

pub type Pt2Filterf32 = Pt2Filter<f32, f32>;
pub type Pt2FilterVector2df32 = Pt2Filter<Vector2d<f32>, f32>;
pub type Pt2FilterVector2df64 = Pt2Filter<Vector2d<f64>, f64>;

pub type Pt2Filterf64 = Pt2Filter<f64, f64>;
pub type Pt2FilterVector3df32 = Pt2Filter<Vector3d<f32>, f32>;
pub type Pt2FilterVector3df64 = Pt2Filter<Vector3d<f64>, f64>;

pub type Pt3Filterf32 = Pt3Filter<f32, f32>;
pub type Pt3FilterVector2df32 = Pt3Filter<Vector2d<f32>, f32>;
pub type Pt3FilterVector2df64 = Pt3Filter<Vector2d<f64>, f64>;

pub type Pt3Filterf64 = Pt3Filter<f64, f64>;
pub type Pt3FilterVector3df32 = Pt3Filter<Vector3d<f32>, f32>;
pub type Pt3FilterVector3df64 = Pt3Filter<Vector3d<f64>, f64>;

/// Discrete-time, first-order low-pass filter (Proportional Time element).<br>
/// It is implemented as a stateful struct that allows for efficient, in-place smoothing of sensor data or motor setpoints."
///
/// The discrete-time transfer function is:
///
/// $$y_{n} = y_{n-1} + k \cdot (x_{n} - y_{n-1})$$
///
/// where $k$ is calculated from the time constant $T$ and sample time $dt$:
/// $k = \frac{dt}{T + dt}$
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt1Filter<T, R> {
    state: T,
    k: R,
}

/// Default is k = 1.0, which is passthrough
impl<T, R> Default for Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T> + AddAssign,
    R: One,
{
    fn default() -> Self {
        Self::new(R::one())
    }
}

impl<T, R> Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T> + AddAssign,
    R: One,
{
    pub fn new(k: R) -> Self {
        Self { state: T::zero(), k }
    }
}
/*
use num_traits::MulAdd;

impl<T, R> SignalFilter<T, R> for Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy,
{
    fn update(&mut self, input: T) -> T {
        // Equation: state = (input - state) * k + state
        // Using mul_add(multiplier, addend)
        self.state = (input - self.state).mul_add(self.k, self.state);
        self.state
    }
}
 */

impl<T, R> SignalFilter<T, R> for Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T> + AddAssign + MulAdd<R, T, Output = T>,
    R: Copy,
{
    fn reset(&mut self) {
        self.state = T::zero();
    }

    fn update(&mut self, input: T) -> T {
        // Equation: state = (input - state) * k + state
        // Using mul_add(multiplier, addend)
        //self.state = (input - self.state).mul_add(self.k, self.state);
        self.state += (input - self.state) * self.k;
        self.state
    }
}

impl<T, R> Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T> + AddAssign + MulAdd<R, T, Output = T>,
    R: Copy + Zero + One,
{
    pub fn set_to_passthrough(&mut self) {
        self.k = R::one();
        self.reset();
    }

    pub fn set_k(&mut self, k: R) {
        self.k = k;
        self.reset();
    }

    // for testing
    #[allow(dead_code)]
    fn k(self) -> R {
        self.k
    }
    #[allow(dead_code)]
    fn state(self) -> T {
        self.state
    }
}

impl<T, R> Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T> + AddAssign + MulAdd<R, T, Output = T>,
    R: Copy + Zero + One + MathConstants + PartialOrd + Div<R, Output = R>,
{
    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: R, delta_t: R) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: R, delta_t: R) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    // Calculates filter gain based on delay (time constant of filter) - time it takes for filter response to reach 63.2% of a step input.
    pub fn gain_from_delay(delay: R, delta_t: R) -> R {
        if delay <= R::zero() {
            // gain of 1.0 means no filtering
            return R::one();
        }
        let omega = delta_t / delay;
        omega / (omega + R::one())
    }

    pub fn gain_from_frequency(cutoff_frequency_hz: R, delta_t: R) -> R {
        let omega = (R::one() + R::one()) * R::PI * cutoff_frequency_hz * delta_t;
        omega / (omega + R::one())
    }
    pub fn gain_from_frequency2(cutoff_frequency_hz: R, delta_t: R) -> R {
        let omega = (R::one() + R::one()) * R::one() * cutoff_frequency_hz * delta_t;
        omega / (omega + R::one())
    }
}

/// Discrete-time, second-order low-pass filter (Proportional Time element).<br>
/// This is equivalent to two cascaded PT1 filters with the same time constant.
///
/// The discrete-time difference equations are:
///
/// $$w_{n} = w_{n-1} + k \cdot (x_{n} - w_{n-1})$$
/// $$y_{n} = y_{n-1} + k \cdot (w_{n} - y_{n-1})$$
///
/// where:
/// - $x_{n}$ is the raw input
/// - $w_{n}$ is the internal state (output of the first stage)
/// - $y_{n}$ is the final filtered output
/// - $k$ is the filter gain
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt2Filter<T, R> {
    state: [T; 2],
    k: R,
}

/// Default is k = 1.0, which is passthrough
impl<T, R> Default for Pt2Filter<T, R>
where
    T: Zero + AddAssign,
    R: One,
{
    fn default() -> Self {
        Self::new(R::one())
    }
}

impl<T, R> Pt2Filter<T, R>
where
    T: Zero + AddAssign,
    R: One,
{
    pub fn new(k: R) -> Self {
        Self { state: [T::zero(), T::zero()], k }
    }
}

impl<T, R> SignalFilter<T, R> for Pt2Filter<T, R>
where
    T: Copy + Zero + AddAssign + Sub<Output = T> + Mul<R, Output = T> + MulAdd<R, T, Output = T>,
    R: Copy,
{
    fn reset(&mut self) {
        self.state = [T::zero(), T::zero()];
    }

    fn update(&mut self, input: T) -> T {
        self.state[1] += (input - self.state[1]) * self.k;
        self.state[0] += (self.state[1] - self.state[0]) * self.k;
        self.state[0]
    }
}

impl<T, R> Pt2Filter<T, R>
where
    T: Copy + Zero + AddAssign + Sub<Output = T> + Mul<R, Output = T> + MulAdd<R, T, Output = T>,
    R: Copy + Zero + One,
{
    pub fn set_to_passthrough(&mut self) {
        self.k = R::one();
        self.reset();
    }

    pub fn set_k(&mut self, k: R) {
        self.k = k;
        self.reset();
    }

    // for testing
    #[allow(dead_code)]
    fn k(self) -> R {
        self.k
    }
    #[allow(dead_code)]
    fn state(self) -> [T; 2] {
        self.state
    }
}

impl<T, R> Pt2Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T> + AddAssign + MulAdd<R, T, Output = T>,
    R: Copy + Zero + One + MathConstants + PartialOrd + Div<R, Output = R>,
{
    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: R, delta_t: R) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: R, delta_t: R) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    pub fn gain_from_delay(delay: R, delta_t: R) -> R {
        Pt1Filter::<T, R>::gain_from_delay(delay * R::FILTER_PT2_CUTOFF_CORRECTION, delta_t)
    }
    pub fn gain_from_frequency(cutoff_frequency_hz: R, delta_t: R) -> R {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        Pt1Filter::<T, R>::gain_from_frequency(cutoff_frequency_hz * R::FILTER_PT2_CUTOFF_CORRECTION, delta_t)
    }
}

/// Discrete-time, third-order low-pass filter (Proportional Time element).<br>
/// This is equivalent to three cascaded PT1 filters. It provides a very steep
/// 60dB/decade roll-off.
///
/// The discrete-time difference equations are:
///
/// $$u_{n} = u_{n-1} + k \cdot (x_{n} - u_{n-1})$$
/// $$v_{n} = v_{n-1} + k \cdot (u_{n} - v_{n-1})$$
/// $$y_{n} = y_{n-1} + k \cdot (v_{n} - y_{n-1})$$
///
/// where $u_{n}$ and $v_{n}$ are internal intermediate states, and $y_{n}$ is the final output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt3Filter<T, R> {
    state: [T; 3],
    k: R,
}

/// Default is k = 1.0, which is passthrough
impl<T, R> Default for Pt3Filter<T, R>
where
    T: Zero,
    R: One,
{
    fn default() -> Self {
        Self::new(R::one())
    }
}

impl<T, R> Pt3Filter<T, R>
where
    T: Zero,
    R: One,
{
    pub fn new(k: R) -> Self {
        Self { state: [T::zero(), T::zero(), T::zero()], k }
    }
}

impl<T, R> SignalFilter<T, R> for Pt3Filter<T, R>
where
    T: Copy + Zero + AddAssign + Sub<Output = T> + Mul<R, Output = T> + MulAdd<R, T, Output = T>,
    R: Copy + Zero,
{
    fn reset(&mut self) {
        self.state = [T::zero(), T::zero(), T::zero()];
    }

    fn update(&mut self, input: T) -> T {
        self.state[2] += (input - self.state[2]) * self.k;
        self.state[1] += (self.state[2] - self.state[1]) * self.k;
        self.state[0] += (self.state[1] - self.state[0]) * self.k;
        self.state[0]
    }
}

impl<T, R> Pt3Filter<T, R>
where
    T: Copy + Zero + AddAssign + Sub<Output = T> + Mul<R, Output = T> + MulAdd<R, T, Output = T>,
    R: Copy + Zero + One,
{
    pub fn set_to_passthrough(&mut self) {
        self.k = R::one();
        self.reset();
    }

    pub fn set_k(&mut self, k: R) {
        self.k = k;
        self.reset();
    }

    // for testing
    #[allow(dead_code)]
    fn k(self) -> R {
        self.k
    }
    #[allow(dead_code)]
    fn state(self) -> [T; 3] {
        self.state
    }
}

impl<T, R> Pt3Filter<T, R>
where
    T: Copy + Zero + AddAssign + Sub<Output = T> + Mul<R, Output = T> + AddAssign + MulAdd<R, T, Output = T>,
    R: Copy + Zero + One + MathConstants + PartialOrd + Div<R, Output = R>,
{
    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: R, delta_t: R) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: R, delta_t: R) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    pub fn gain_from_delay(delay: R, delta_t: R) -> R {
        Pt1Filter::<T, R>::gain_from_delay(delay * R::FILTER_PT3_CUTOFF_CORRECTION, delta_t)
    }

    pub fn gain_from_frequency(cutoff_frequency_hz: R, delta_t: R) -> R {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        Pt1Filter::<T, R>::gain_from_frequency(cutoff_frequency_hz * R::FILTER_PT3_CUTOFF_CORRECTION, delta_t)
    }
}

#[cfg(any(debug_assertions, test))]
mod tests {
    #![allow(unused)]
    #![allow(clippy::float_cmp)]
    #![allow(unused_results)]
    use crate::UpdateFilter;

    use super::*;
    use vector_quaternion_matrix::Vector3df32;
    use vector_quaternion_matrix::Vector3di16;

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<Pt1Filter<f32, f32>>();
        is_full::<Pt1Filterf32>();
        is_full::<Pt2Filter<f32, f32>>();
        is_full::<Pt2Filterf32>();
        is_full::<Pt3Filter<f32, f32>>();
        is_full::<Pt3Filterf32>();
    }
    #[test]
    fn pt1_filter_f32() {
        let mut filter = Pt1Filterf32::new(1.0);

        let mut reading: f32 = 2.7;
        reading = filter.update(reading);
        assert_eq!(2.7, reading);

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(-1.0, filter.update(-1.0));

        filter.reset();
        assert_eq!(0.0, filter.state());
        assert_eq!(4.0, filter.update(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.385_869_56, filter.k());
        assert_eq!(0.385_869_56, filter.update(1.0));
        assert_eq!(1.008_713_4, filter.update(2.0));

        filter.set_k(1.0);
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));

        filter.reset();
        filter.set_k(0.5);
        assert_eq!(0.5, filter.update(1.0));
        assert_eq!(1.25, filter.update(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.385_869_56, filter.update(1.0));
        assert_eq!(1.008_713_4, filter.update(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));
    }
    #[test]
    fn pt1_filter_f32_method_call() {
        let mut filter = Pt1Filterf32::new(0.2);
        assert_eq!(0.2, filter.update(1.0));
        assert_eq!(0.2, filter.update(0.2));

        filter.reset();
        let mut value: f32 = 1.0;
        value.update_using(&mut filter);
        assert_eq!(0.2, value);
        value.update_using(&mut filter);
        assert_eq!(0.2, value);
    }
    #[test]
    fn pt1_filter_vector3df32_method_call() {
        let mut filter = Pt1Filterf32::new(0.25);
        assert_eq!(0.05, filter.update(0.2));
        filter.reset();
        assert_eq!(0.125, filter.update(0.5));
        filter.reset();
        assert_eq!(0.375, filter.update(1.5));

        let mut filter = Pt1FilterVector3df32::new(0.25);
        let value = Vector3df32 { x: 0.2, y: 0.5, z: 1.5 };
        let output = filter.update(value);
        assert_eq!(Vector3df32 { x: 0.05, y: 0.125, z: 0.375 }, output);

        filter.reset();
        let mut value = Vector3df32 { x: 0.2, y: 0.5, z: 1.5 };
        value.update_using(&mut filter);
        assert_eq!(Vector3df32 { x: 0.05, y: 0.125, z: 0.375 }, value);
    }
    #[test]
    fn pt2_filter_f32() {
        let mut filter = Pt2Filterf32::new(1.0);

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(-1.0, filter.update(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.update(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.244_031_07, filter.update(1.0));
        assert_eq!(0.735_024_03, filter.update(2.0));

        filter.set_k(1.0);
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.244_031_07, filter.update(1.0));
        assert_eq!(0.735_024_03, filter.update(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));
    }
    #[test]
    fn pt2_filter_f32_method_call() {
        let mut filter = Pt2Filterf32::new(0.2);
        assert_eq!(0.040_000_003, filter.update(1.0));
        assert_eq!(0.0656, filter.update(0.040_000_003));

        filter.reset();
        let mut value: f32 = 1.0;
        value.update_using(&mut filter);
        assert_eq!(0.040_000_003, value);
        value.update_using(&mut filter);
        assert_eq!(0.0656, value);
    }
    #[test]
    fn pt3_filter_f32() {
        let mut filter = Pt3Filterf32::new(1.0);

        let mut state = filter.state();
        assert_eq!([0.0, 0.0, 0.0], state);

        // test that filter with default settings performs no filtering
        let mut output = filter.update(1.0);
        assert_eq!(1.0, output);
        state = filter.state();
        assert_eq!([1.0, 1.0, 1.0], state);

        output = filter.update(1.0);
        state = filter.state();
        assert_eq!([1.0, 1.0, 1.0], state);
        assert_eq!(1.0, output);

        assert_eq!(-1.0, filter.update(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.update(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.168_247_66, filter.update(1.0));
        assert_eq!(0.562_591_97, filter.update(2.0));

        filter.set_k(1.0);
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.168_247_66, filter.update(1.0));
        assert_eq!(0.562_591_97, filter.update(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));
    }
    #[test]
    fn pt1_filter_vector3df32() {
        let mut filter = Pt1Filter::<Vector3df32, f32>::new(1.0);
        let mut output: Vector3df32;

        // test that filter with default settings performs no filtering
        output = filter.update(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 });
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, output);
        let state = filter.state();
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, state);

        filter.reset();
        let state = filter.state();
        assert_eq!(Vector3df32 { x: 0.0, y: 0.0, z: 0.0 }, state);

        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.385_869_56, filter.update(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(1.008_713_4, filter.update(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);

        filter.set_k(1.0);
        assert_eq!(1.0, filter.update(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.update(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.385_869_56, filter.update(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(1.008_713_4, filter.update(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.update(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.update(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);
    }
    /*#[test]
    fn pt1_filter_vector3df32_i16() {
        let mut filter = Pt1Filter::<Vector3di16, f32>::new(1.0);
        let mut output: Vector3di16;
        let mut state: Vector3di16;

        // test that filter with default settings performs no filtering
        output = filter.update(Vector3di16 { x: 2, y: 3, z: 5 });
        assert_eq!(Vector3di16 { x: 2, y: 3, z: 5 }, output);
        state = filter.state();
        assert_eq!(Vector3di16 { x: 2, y: 3, z: 5 }, state);
    }
    #[test]
    fn pt1_filter_vector3df32_i32() {
        let mut filter = Pt1Filter::<Vector3di32, f32>::new(1.0);
        let mut output: Vector3di32;
        let mut state: Vector3di32;

        // test that filter with default settings performs no filtering
        output = filter.update(Vector3di32 { x: 2, y: 3, z: 5 });
        assert_eq!(2, output.x);
        assert_eq!(3, output.y);
        assert_eq!(5, output.z);
        state = filter.state();
        assert_eq!(2, state.x);
        assert_eq!(3, state.y);
        assert_eq!(5, state.z);
    }*/
}
