use core::ops::{Add, Div, Mul, Sub};
use num_traits::{One, Zero};
use vector_quaternion_matrix::{MathConstants, MathMethods};

use crate::FilterSignal;

pub type BiquadFilterf32<T> = BiquadFilter<T, f32>;
pub type BiquadFilterf64<T> = BiquadFilter<T, f64>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BiquadFilterState<T> {
    x1: T,
    x2: T,
    y1: T,
    y2: T,
}

impl<T> Default for BiquadFilterState<T>
where
    T: Zero,
{
    fn default() -> Self {
        Self { x1: T::zero(), x2: T::zero(), y1: T::zero(), y2: T::zero() }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BiquadFilter<T, R> {
    state: BiquadFilterState<T>,
    weight: R,
    a1: R,
    a2: R,
    b0: R,
    b1: R,
    b2: R,
    loop_time_seconds: R,
    two_pi_loop_time_seconds: R, // cached value of 2.0 * PI * loop_time_seconds
    q: R,
    one_over_2q: R,
}

impl<T, R> Default for BiquadFilter<T, R>
where
    T: Zero,
    R: Zero + One + Div<R, Output = R>,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Zero,
    R: Zero + One + Div<R, Output = R>,
{
    fn new() -> Self {
        Self {
            state: BiquadFilterState { x1: T::zero(), x2: T::zero(), y1: T::zero(), y2: T::zero() },
            weight: R::one(),
            a1: R::zero(),
            a2: R::zero(),
            b0: R::one(),
            b1: R::zero(),
            b2: R::zero(),
            loop_time_seconds: R::zero(),
            two_pi_loop_time_seconds: R::zero(),
            q: R::one(),
            one_over_2q: R::one() / (R::one() + R::one()),
        }
    }
}

impl<T, R> FilterSignal<T, R> for BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One + Div<R, Output = R>,
{
    fn apply(&mut self, input: T) -> T {
        let output = input * self.b0 + self.state.x1 * self.b1 + self.state.x2 * self.b2
            - self.state.y1 * self.a1
            - self.state.y2 * self.a2;

        self.state.x2 = self.state.x1;
        self.state.x1 = input;
        self.state.y2 = self.state.y1;
        self.state.y1 = output;
        output
    }

    fn reset(&mut self) {
        self.state.x1 = T::zero();
        self.state.x2 = T::zero();
        self.state.y1 = T::zero();
        self.state.y2 = T::zero();
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One + Div<R, Output = R>,
{
    pub fn set_q(&mut self, q: R) {
        self.q = q;
        self.one_over_2q = R::one() / ((R::one() + R::one()) * q); // cache value for faster setting of frequencies
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One,
{
    pub fn set_weight(&mut self, weight: R) {
        self.weight = weight;
    }

    pub fn weight(&self) -> R {
        self.weight
    }

    pub fn set_parameters_and_weight(&mut self, a1: R, a2: R, b0: R, b1: R, b2: R, weight: R) {
        self.weight = weight;
        self.a1 = a1;
        self.a2 = a2;
        self.b0 = b0;
        self.b1 = b1;
        self.b2 = b2;
    }

    pub fn set_parameters(&mut self, a1: R, a2: R, b0: R, b1: R, b2: R) {
        self.set_parameters_and_weight(a1, a2, b0, b1, b2, R::one());
    }

    /// Copy parameters from another Biquad filter
    pub fn set_parameters_from(&mut self, other: &BiquadFilter<T, R>) {
        self.weight = other.weight;
        self.a1 = other.a1;
        self.a2 = other.a2;
        self.b0 = other.b0;
        self.b1 = other.b1;
        self.b2 = other.b2;
    }

    pub fn calculate_omega(&self, frequency: R) -> R {
        frequency * self.two_pi_loop_time_seconds
    }

    pub fn q(&self) -> R {
        self.q
    }

    pub fn loop_time_seconds(&self) -> R {
        self.loop_time_seconds
    }

    // for testing
    #[allow(dead_code)]
    fn state(self) -> BiquadFilterState<T> {
        self.state
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One + MathConstants + MathMethods + Div<R, Output = R> + Sub<R, Output = R>,
{
    pub fn set_to_passthrough(&mut self) {
        self.b0 = R::one();
        self.b1 = R::zero();
        self.b2 = R::zero();
        self.a1 = R::zero();
        self.a2 = R::zero();
        self.weight = R::one();
        self.reset();
    }

    pub fn apply_weighted(&mut self, input: T) -> T {
        let output = self.apply(input);
        // weight of 1.0 gives just output, weight of 0.0 gives just input
        (output - input) * self.weight + input
    }
    pub fn init_low_pass(&mut self, frequency_hz: R, loop_time_seconds: R, q: R) {
        //assert(Q != 0.0 && "Q cannot be zero");
        self.set_loop_time(loop_time_seconds);
        self.set_q(q);
        self.set_low_pass_frequency_assuming_q(frequency_hz);
        self.reset();
    }

    pub fn init_notch(&mut self, frequency_hz: R, loop_time_seconds: R, q: R) {
        //assert(Q != 0.0 && "Q cannot be zero");
        self.set_loop_time(loop_time_seconds);
        self.set_q(q);
        self.set_notch_frequency_assuming_q(frequency_hz);
        self.reset();
    }
    //Note: weight must be in range [0, 1].
    pub fn set_low_pass_frequency_weighted_assuming_q(&mut self, frequency_hz: R, weight: R) {
        self.weight = weight;

        let omega = frequency_hz * self.two_pi_loop_time_seconds;
        let (sin_omega, cos_omega) = omega.sin_cos();
        let alpha = sin_omega * self.one_over_2q;
        let a0_reciprocal = R::one() / (R::one() + alpha);

        self.b1 = (R::one() - cos_omega) * a0_reciprocal;
        self.b0 = self.b1 * (R::one() / (R::one() + R::one()));
        self.b2 = self.b0;
        self.a1 = R::zero() - (R::one() + R::one()) * cos_omega * a0_reciprocal;
        self.a2 = (R::one() - alpha) * a0_reciprocal;
    }

    pub fn set_low_pass_frequency_assuming_q(&mut self, frequency_hz: R) {
        self.set_low_pass_frequency_weighted_assuming_q(frequency_hz, R::one());
    }

    pub fn set_notch_frequency_weighted_from_sin_cos_assuming_q(&mut self, sin_omega: R, cos_omega: R, weight: R) {
        self.weight = weight;

        let alpha = sin_omega * self.one_over_2q;
        let a0reciprocal = R::one() / (R::one() + alpha);

        self.b0 = a0reciprocal;
        self.b2 = a0reciprocal;
        self.b1 = R::zero() - (R::one() + R::one()) * cos_omega * a0reciprocal;
        self.a1 = self.b1;
        self.a2 = (R::one() - alpha) * a0reciprocal;
    }

    pub fn set_notch_frequency_weighted_assuming_q(&mut self, frequency_hz: R, weight: R) {
        let omega = frequency_hz * self.two_pi_loop_time_seconds;
        let (sin_omega, cos_omega) = omega.sin_cos();
        self.set_notch_frequency_weighted_from_sin_cos_assuming_q(sin_omega, cos_omega, weight);
    }

    pub fn set_notch_frequency_assuming_q(&mut self, frequency_hz: R) {
        // assumes Q already set
        self.set_notch_frequency_weighted_assuming_q(frequency_hz, R::one());
    }

    pub fn set_notch_frequency(&mut self, center_frequency_hz: R, lower_cutoff_frequency_hz: R) {
        self.set_q(Self::calculate_q(center_frequency_hz, lower_cutoff_frequency_hz));
        self.set_notch_frequency_assuming_q(center_frequency_hz);
    }

    pub fn calculate_q(center_frequency_hz: R, lower_cutoff_frequency_hz: R) -> R {
        center_frequency_hz * lower_cutoff_frequency_hz
            / (center_frequency_hz * center_frequency_hz - lower_cutoff_frequency_hz * lower_cutoff_frequency_hz)
    }

    pub fn set_q_from_frequencies(&mut self, center_frequency_hz: R, lower_cutoff_frequency_hz: R) {
        self.set_q(Self::calculate_q(center_frequency_hz, lower_cutoff_frequency_hz));
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One + MathConstants,
{
    pub fn set_loop_time(&mut self, loop_time_seconds: R) {
        self.loop_time_seconds = loop_time_seconds;
        self.two_pi_loop_time_seconds = (R::one() + R::one()) * R::PI * loop_time_seconds; // cache value for faster setting of frequencies
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
        is_full::<BiquadFilter<f32, f32>>();
        is_full::<BiquadFilterf32<f32>>();
        is_full::<BiquadFilterState<f32>>();
    }
    #[test]
    fn biquad_filter_f32() {
        let mut filter = BiquadFilterf32::<f32>::default();

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(-1.0, filter.apply(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.apply(4.0));

        filter.set_parameters_and_weight(2.0, 3.0, 5.0, 7.0, 11.0, 13.0);
        filter.set_to_passthrough();
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(2.0, filter.apply(2.0));
        assert_eq!(1.0, filter.apply_weighted(1.0));
        assert_eq!(2.0, filter.apply_weighted(2.0));
    }
    #[test]
    fn biquad_filter_vector3df32() {
        let mut filter = BiquadFilterf32::<Vector3df32>::default();
        let mut output: Vector3df32;
        let mut state: BiquadFilterState<Vector3df32>;

        // test that filter with default settings performs no filtering
        output = filter.apply(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 });
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, output);
        state = filter.state();
        assert_eq!(2.0, state.x1.x);
        assert_eq!(0.0, state.x2.x);
        assert_eq!(2.0, state.y1.x);
        assert_eq!(0.0, state.y2.x);

        filter.reset();
        state = filter.state();
        assert_eq!(0.0, state.x1.x);
        assert_eq!(0.0, state.x2.x);
        assert_eq!(0.0, state.y1.x);
        assert_eq!(0.0, state.y2.x);
        assert_eq!(4.0, filter.apply(Vector3df32 { x: 4.0, y: 0.0, z: 0.0 }).x);

        filter.set_parameters_and_weight(2.0, 3.0, 5.0, 7.0, 11.0, 13.0);
        filter.set_to_passthrough();
        assert_eq!(1.0, filter.apply(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.apply(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(1.0, filter.apply_weighted(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.apply_weighted(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);
    }
}
