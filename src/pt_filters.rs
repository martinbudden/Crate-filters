use core::ops::{Add, Div, Mul, Sub};
use num_traits::{One, Zero};
use vector_quaternion_matrix::MathConstants;

use crate::FilterSignal;

pub type Pt1Filterf32<T> = Pt1Filter<T, f32>;
pub type Pt1Filterf64<T> = Pt1Filter<T, f64>;

pub type Pt2Filterf32<T> = Pt2Filter<T, f32>;
pub type Pt2Filterf64<T> = Pt2Filter<T, f64>;

pub type Pt3Filterf32<T> = Pt3Filter<T, f32>;
pub type Pt3Filterf64<T> = Pt3Filter<T, f64>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt1Filter<T, R> {
    state: T,
    k: R,
}

/// Default is k = 1.0, which is passthrough
impl<T, R> Default for Pt1Filter<T, R>
where
    T: Zero,
    R: One,
{
    fn default() -> Self {
        Self::new(R::one())
    }
}

impl<T, R> Pt1Filter<T, R>
where
    T: Zero,
    R: One,
{
    pub fn new(k: R) -> Self {
        Self { state: T::zero(), k }
    }
}

impl<T, R> FilterSignal<T, R> for Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy,
{
    fn apply(&mut self, input: T) -> T {
        self.state = self.state + (input - self.state) * self.k; // equivalent to self.state = self.k*input + (1.0 - self.k)*self.state;
        self.state
    }

    fn reset(&mut self) {
        self.state = T::zero();
    }
}

impl<T, R> Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
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
    fn state(self) -> T {
        self.state
    }
}

impl<T, R> Pt1Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
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

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt2Filter<T, R> {
    state: [T; 2],
    k: R,
}

/// Default is k = 1.0, which is passthrough
impl<T, R> Default for Pt2Filter<T, R>
where
    T: Zero,
    R: One,
{
    fn default() -> Self {
        Self::new(R::one())
    }
}

impl<T, R> Pt2Filter<T, R>
where
    T: Zero,
    R: One,
{
    pub fn new(k: R) -> Self {
        Self { state: [T::zero(), T::zero()], k }
    }
}

impl<T, R> FilterSignal<T, R> for Pt2Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy,
{
    fn apply(&mut self, input: T) -> T {
        self.state[1] = self.state[1] + (input - self.state[1]) * self.k;
        self.state[0] = self.state[0] + (self.state[1] - self.state[0]) * self.k;
        self.state[0]
    }

    fn reset(&mut self) {
        self.state = [T::zero(), T::zero()];
    }
}

impl<T, R> Pt2Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
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
    fn state(self) -> [T; 2] {
        self.state
    }
}

impl<T, R> Pt2Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
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

impl<T, R> FilterSignal<T, R> for Pt3Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero,
{
    fn apply(&mut self, input: T) -> T {
        self.state[2] = self.state[2] + (input - self.state[2]) * self.k;
        self.state[1] = self.state[1] + (self.state[2] - self.state[1]) * self.k;
        self.state[0] = self.state[0] + (self.state[1] - self.state[0]) * self.k;
        self.state[0]
    }

    fn reset(&mut self) {
        self.state = [T::zero(), T::zero(), T::zero()];
    }
}

impl<T, R> Pt3Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
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
    fn state(self) -> [T; 3] {
        self.state
    }
}

impl<T, R> Pt3Filter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
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
    use super::*;
    use vector_quaternion_matrix::Vector3df32;
    use vector_quaternion_matrix::Vector3di16;
    use vector_quaternion_matrix::Vector3di32;

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<Pt1Filter<f32, f32>>();
        is_full::<Pt1Filterf32<f32>>();
        is_full::<Pt2Filter<f32, f32>>();
        is_full::<Pt2Filterf32<f32>>();
        is_full::<Pt3Filter<f32, f32>>();
        is_full::<Pt3Filterf32<f32>>();
    }
    #[test]
    fn pt1_filter_f32() {
        let mut filter = Pt1Filterf32::<f32>::new(1.0);

        let mut reading: f32 = 2.7;
        reading = filter.apply(reading);
        assert_eq!(2.7, reading);

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(-1.0, filter.apply(-1.0));

        filter.reset();
        assert_eq!(0.0, filter.state());
        assert_eq!(4.0, filter.apply(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.38586956, filter.apply(1.0));
        assert_eq!(1.0087134, filter.apply(2.0));

        filter.set_k(1.0);
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.38586956, filter.apply(1.0));
        assert_eq!(1.0087134, filter.apply(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));
    }
    #[test]
    fn pt2_filter_f32() {
        let mut filter = Pt2Filterf32::<f32>::new(1.0);

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(-1.0, filter.apply(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.apply(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.24403107, filter.apply(1.0));
        assert_eq!(0.73502403, filter.apply(2.0));

        filter.set_k(1.0);
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.24403107, filter.apply(1.0));
        assert_eq!(0.73502403, filter.apply(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));
    }
    #[test]
    fn pt3_filter_f32() {
        let mut filter = Pt3Filterf32::<f32>::new(1.0);

        let mut state = filter.state();
        assert_eq!([0.0, 0.0, 0.0], state);

        // test that filter with default settings performs no filtering
        let mut output = filter.apply(1.0);
        assert_eq!(1.0, output);
        state = filter.state();
        assert_eq!([1.0, 1.0, 1.0], state);

        output = filter.apply(1.0);
        state = filter.state();
        assert_eq!([1.0, 1.0, 1.0], state);
        assert_eq!(1.0, output);

        assert_eq!(-1.0, filter.apply(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.apply(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.16824766, filter.apply(1.0));
        assert_eq!(0.56259197, filter.apply(2.0));

        filter.set_k(1.0);
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.16824766, filter.apply(1.0));
        assert_eq!(0.56259197, filter.apply(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));
    }
    #[test]
    fn pt1_filter_vector3df32() {
        let mut filter = Pt1Filter::<Vector3df32, f32>::new(1.0);
        let mut output: Vector3df32;

        // test that filter with default settings performs no filtering
        output = filter.apply(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 });
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, output);
        let state = filter.state();
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, state);

        filter.reset();
        let state = filter.state();
        assert_eq!(Vector3df32 { x: 0.0, y: 0.0, z: 0.0 }, state);

        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.38586956, filter.apply(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(1.0087134, filter.apply(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);

        filter.set_k(1.0);
        assert_eq!(1.0, filter.apply(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.apply(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.38586956, filter.apply(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(1.0087134, filter.apply(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.apply(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.apply(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);
    }
    #[test]
    fn pt1_filter_vector3df32_i16() {
        let mut filter = Pt1Filter::<Vector3di16, f32>::new(1.0);
        let mut output: Vector3di16;
        let mut state: Vector3di16;

        // test that filter with default settings performs no filtering
        output = filter.apply(Vector3di16 { x: 2, y: 3, z: 5 });
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
        output = filter.apply(Vector3di32 { x: 2, y: 3, z: 5 });
        assert_eq!(2, output.x);
        assert_eq!(3, output.y);
        assert_eq!(5, output.z);
        state = filter.state();
        assert_eq!(2, state.x);
        assert_eq!(3, state.y);
        assert_eq!(5, state.z);
    }
}
