use core::ops::{Add, Div, Mul, Sub};
use num_traits::{One, Zero};
use vector_quaternion_matrix::{MathConstants, MathMethods};

pub type Pt1Filterf32<T> = Pt1Filter<T, f32>;
pub type Pt1Filterf64<T> = Pt1Filter<T, f64>;
pub type Pt2Filterf32<T> = Pt2Filter<T, f32>;
pub type Pt2Filterf64<T> = Pt2Filter<T, f64>;
pub type Pt3Filterf32<T> = Pt3Filter<T, f32>;
pub type Pt3Filterf64<T> = Pt3Filter<T, f64>;
pub type BiquadFilterf32<T> = BiquadFilter<T, f32>;
pub type BiquadFilterf64<T> = BiquadFilter<T, f64>;

pub trait FilterSignal<T, F> {
    fn apply(&mut self, input: T) -> T;
    fn reset(&mut self);
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt1Filter<T, F> {
    state: T,
    k: F,
}

/// Default is k = 1.0, which is passthrough
impl<T, F> Default for Pt1Filter<T, F>
where
    T: Zero,
    F: One,
{
    fn default() -> Self {
        Self::new(F::one())
    }
}

impl<T, F> Pt1Filter<T, F>
where
    T: Zero,
    F: One,
{
    pub fn new(k: F) -> Self {
        Self { state: T::zero(), k }
    }
}

impl<T, F> FilterSignal<T, F> for Pt1Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy,
{
    fn apply(&mut self, input: T) -> T {
        self.state = self.state + (input - self.state) * self.k; // equivalent to self.state = self.k*input + (1.0 - self.k)*self.state;
        self.state
    }

    fn reset(&mut self) {
        self.state = T::zero();
    }
}

impl<T, F> Pt1Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One,
{
    pub fn set_to_passthrough(&mut self) {
        self.k = F::one();
        self.reset();
    }

    pub fn set_k(&mut self, k: F) {
        self.k = k;
        self.reset();
    }

    // for testing
    #[allow(dead_code)]
    fn state(self) -> T {
        self.state
    }
}

impl<T, F> Pt1Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + MathConstants + PartialOrd + Div<F, Output = F>,
{
    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: F, delta_t: F) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: F, delta_t: F) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    // Calculates filter gain based on delay (time constant of filter) - time it takes for filter response to reach 63.2% of a step input.
    pub fn gain_from_delay(delay: F, delta_t: F) -> F {
        if delay <= F::zero() {
            // gain of 1.0 means no filtering
            return F::one();
        }
        let omega = delta_t / delay;
        omega / (omega + F::one())
    }

    pub fn gain_from_frequency(cutoff_frequency_hz: F, delta_t: F) -> F {
        let omega = (F::one() + F::one()) * F::PI * cutoff_frequency_hz * delta_t;
        omega / (omega + F::one())
    }
    pub fn gain_from_frequency2(cutoff_frequency_hz: F, delta_t: F) -> F {
        let omega = (F::one() + F::one()) * F::one() * cutoff_frequency_hz * delta_t;
        omega / (omega + F::one())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt2Filter<T, F> {
    state: [T; 2],
    k: F,
}

/// Default is k = 1.0, which is passthrough
impl<T, F> Default for Pt2Filter<T, F>
where
    T: Zero,
    F: One,
{
    fn default() -> Self {
        Self::new(F::one())
    }
}

impl<T, F> Pt2Filter<T, F>
where
    T: Zero,
    F: One,
{
    pub fn new(k: F) -> Self {
        Self { state: [T::zero(), T::zero()], k }
    }
}

impl<T, F> FilterSignal<T, F> for Pt2Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy,
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

impl<T, F> Pt2Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One,
{
    pub fn set_to_passthrough(&mut self) {
        self.k = F::one();
        self.reset();
    }

    pub fn set_k(&mut self, k: F) {
        self.k = k;
        self.reset();
    }

    // for testing
    #[allow(dead_code)]
    fn state(self) -> [T; 2] {
        self.state
    }
}

impl<T, F> Pt2Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + MathConstants + PartialOrd + Div<F, Output = F>,
{
    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: F, delta_t: F) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: F, delta_t: F) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    pub fn gain_from_delay(delay: F, delta_t: F) -> F {
        Pt1Filter::<T, F>::gain_from_delay(delay * F::FILTER_PT2_CUTOFF_CORRECTION, delta_t)
    }
    pub fn gain_from_frequency(cutoff_frequency_hz: F, delta_t: F) -> F {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        Pt1Filter::<T, F>::gain_from_frequency(cutoff_frequency_hz * F::FILTER_PT2_CUTOFF_CORRECTION, delta_t)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pt3Filter<T, F> {
    state: [T; 3],
    k: F,
}

/// Default is k = 1.0, which is passthrough
impl<T, F> Default for Pt3Filter<T, F>
where
    T: Zero,
    F: One,
{
    fn default() -> Self {
        Self::new(F::one())
    }
}

impl<T, F> Pt3Filter<T, F>
where
    T: Zero,
    F: One,
{
    pub fn new(k: F) -> Self {
        Self { state: [T::zero(), T::zero(), T::zero()], k }
    }
}

impl<T, F> FilterSignal<T, F> for Pt3Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero,
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

impl<T, F> Pt3Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One,
{
    pub fn set_to_passthrough(&mut self) {
        self.k = F::one();
        self.reset();
    }

    pub fn set_k(&mut self, k: F) {
        self.k = k;
        self.reset();
    }

    // for testing
    #[allow(dead_code)]
    fn state(self) -> [T; 3] {
        self.state
    }
}

impl<T, F> Pt3Filter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + MathConstants + PartialOrd + Div<F, Output = F>,
{
    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: F, delta_t: F) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: F, delta_t: F) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    pub fn gain_from_delay(delay: F, delta_t: F) -> F {
        Pt1Filter::<T, F>::gain_from_delay(delay * F::FILTER_PT3_CUTOFF_CORRECTION, delta_t)
    }

    pub fn gain_from_frequency(cutoff_frequency_hz: F, delta_t: F) -> F {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        Pt1Filter::<T, F>::gain_from_frequency(cutoff_frequency_hz * F::FILTER_PT3_CUTOFF_CORRECTION, delta_t)
    }
}

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
pub struct BiquadFilter<T, F> {
    state: BiquadFilterState<T>,
    weight: F,
    a1: F,
    a2: F,
    b0: F,
    b1: F,
    b2: F,
    loop_time_seconds: F,
    two_pi_loop_time_seconds: F, // cached value of 2.0 * PI * loop_time_seconds
    q: F,
    one_over_2q: F,
}

impl<T, F> Default for BiquadFilter<T, F>
where
    T: Zero,
    F: Zero + One + Div<F, Output = F>,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, F> BiquadFilter<T, F>
where
    T: Zero,
    F: Zero + One + Div<F, Output = F>,
{
    fn new() -> Self {
        Self {
            state: BiquadFilterState { x1: T::zero(), x2: T::zero(), y1: T::zero(), y2: T::zero() },
            weight: F::one(),
            a1: F::zero(),
            a2: F::zero(),
            b0: F::one(),
            b1: F::zero(),
            b2: F::zero(),
            loop_time_seconds: F::zero(),
            two_pi_loop_time_seconds: F::zero(),
            q: F::one(),
            one_over_2q: F::one() / (F::one() + F::one()),
        }
    }
}

impl<T, F> FilterSignal<T, F> for BiquadFilter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + Div<F, Output = F>,
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

impl<T, F> BiquadFilter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + Div<F, Output = F>,
{
    pub fn set_q(&mut self, q: F) {
        self.q = q;
        self.one_over_2q = F::one() / ((F::one() + F::one()) * q); // cache value for faster setting of frequencies
    }
}

impl<T, F> BiquadFilter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One,
{
    pub fn set_weight(&mut self, weight: F) {
        self.weight = weight;
    }

    pub fn weight(&self) -> F {
        self.weight
    }

    pub fn set_parameters_and_weight(&mut self, a1: F, a2: F, b0: F, b1: F, b2: F, weight: F) {
        self.weight = weight;
        self.a1 = a1;
        self.a2 = a2;
        self.b0 = b0;
        self.b1 = b1;
        self.b2 = b2;
    }

    pub fn set_parameters(&mut self, a1: F, a2: F, b0: F, b1: F, b2: F) {
        self.set_parameters_and_weight(a1, a2, b0, b1, b2, F::one());
    }

    /// Copy parameters from another Biquad filter
    pub fn set_parameters_from(&mut self, other: &BiquadFilter<T, F>) {
        self.weight = other.weight;
        self.a1 = other.a1;
        self.a2 = other.a2;
        self.b0 = other.b0;
        self.b1 = other.b1;
        self.b2 = other.b2;
    }

    pub fn calculate_omega(&self, frequency: F) -> F {
        frequency * self.two_pi_loop_time_seconds
    }

    pub fn q(&self) -> F {
        self.q
    }

    pub fn loop_time_seconds(&self) -> F {
        self.loop_time_seconds
    }

    // for testing
    #[allow(dead_code)]
    fn state(self) -> BiquadFilterState<T> {
        self.state
    }
}

impl<T, F> BiquadFilter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + MathConstants + MathMethods + Div<F, Output = F> + Sub<F, Output = F>,
{
    pub fn set_to_passthrough(&mut self) {
        self.b0 = F::one();
        self.b1 = F::zero();
        self.b2 = F::zero();
        self.a1 = F::zero();
        self.a2 = F::zero();
        self.weight = F::one();
        self.reset();
    }

    pub fn apply_weighted(&mut self, input: T) -> T {
        let output = self.apply(input);
        // weight of 1.0 gives just output, weight of 0.0 gives just input
        (output - input) * self.weight + input
    }
    pub fn init_low_pass(&mut self, frequency_hz: F, loop_time_seconds: F, q: F) {
        //assert(Q != 0.0 && "Q cannot be zero");
        self.set_loop_time(loop_time_seconds);
        self.set_q(q);
        self.set_low_pass_frequency_assuming_q(frequency_hz);
        self.reset();
    }

    pub fn init_notch(&mut self, frequency_hz: F, loop_time_seconds: F, q: F) {
        //assert(Q != 0.0 && "Q cannot be zero");
        self.set_loop_time(loop_time_seconds);
        self.set_q(q);
        self.set_notch_frequency_assuming_q(frequency_hz);
        self.reset();
    }
    //Note: weight must be in range [0, 1].
    pub fn set_low_pass_frequency_weighted_assuming_q(&mut self, frequency_hz: F, weight: F) {
        self.weight = weight;

        let omega = frequency_hz * self.two_pi_loop_time_seconds;
        let (sin_omega, cos_omega) = omega.sin_cos();
        let alpha = sin_omega * self.one_over_2q;
        let a0_reciprocal = F::one() / (F::one() + alpha);

        self.b1 = (F::one() - cos_omega) * a0_reciprocal;
        self.b0 = self.b1 * (F::one() / (F::one() + F::one()));
        self.b2 = self.b0;
        self.a1 = F::zero() - (F::one() + F::one()) * cos_omega * a0_reciprocal;
        self.a2 = (F::one() - alpha) * a0_reciprocal;
    }

    pub fn set_low_pass_frequency_assuming_q(&mut self, frequency_hz: F) {
        self.set_low_pass_frequency_weighted_assuming_q(frequency_hz, F::one());
    }

    pub fn set_notch_frequency_weighted_from_sin_cos_assuming_q(&mut self, sin_omega: F, cos_omega: F, weight: F) {
        self.weight = weight;

        let alpha = sin_omega * self.one_over_2q;
        let a0reciprocal = F::one() / (F::one() + alpha);

        self.b0 = a0reciprocal;
        self.b2 = a0reciprocal;
        self.b1 = F::zero() - (F::one() + F::one()) * cos_omega * a0reciprocal;
        self.a1 = self.b1;
        self.a2 = (F::one() - alpha) * a0reciprocal;
    }

    pub fn set_notch_frequency_weighted_assuming_q(&mut self, frequency_hz: F, weight: F) {
        let omega = frequency_hz * self.two_pi_loop_time_seconds;
        let (sin_omega, cos_omega) = omega.sin_cos();
        self.set_notch_frequency_weighted_from_sin_cos_assuming_q(sin_omega, cos_omega, weight);
    }

    pub fn set_notch_frequency_assuming_q(&mut self, frequency_hz: F) {
        // assumes Q already set
        self.set_notch_frequency_weighted_assuming_q(frequency_hz, F::one());
    }

    pub fn set_notch_frequency(&mut self, center_frequency_hz: F, lower_cutoff_frequency_hz: F) {
        self.set_q(Self::calculate_q(center_frequency_hz, lower_cutoff_frequency_hz));
        self.set_notch_frequency_assuming_q(center_frequency_hz);
    }

    pub fn calculate_q(center_frequency_hz: F, lower_cutoff_frequency_hz: F) -> F {
        center_frequency_hz * lower_cutoff_frequency_hz
            / (center_frequency_hz * center_frequency_hz - lower_cutoff_frequency_hz * lower_cutoff_frequency_hz)
    }

    pub fn set_q_from_frequencies(&mut self, center_frequency_hz: F, lower_cutoff_frequency_hz: F) {
        self.set_q(Self::calculate_q(center_frequency_hz, lower_cutoff_frequency_hz));
    }
}

impl<T, F> BiquadFilter<T, F>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<F, Output = T>,
    F: Copy + Zero + One + MathConstants,
{
    pub fn set_loop_time(&mut self, loop_time_seconds: F) {
        self.loop_time_seconds = loop_time_seconds;
        self.two_pi_loop_time_seconds = (F::one() + F::one()) * F::PI * loop_time_seconds; // cache value for faster setting of frequencies
    }
}

/// Simple moving average filter.
/// See [Moving Average Filter - Theory and Software Implementation - Phil's Lab #21](https://www.youtube.com/watch?v=rttn46_Y3c8).

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingAverageFilter<T, const N: usize> {
    count: usize,
    index: usize,
    sum: T,
    samples: [T; N],
}

impl<T, const N: usize> Default for MovingAverageFilter<T, N>
where
    T: Copy + Zero,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, const N: usize> MovingAverageFilter<T, N>
where
    T: Copy + Zero,
{
    pub fn new() -> Self {
        Self { count: 0, index: 0, sum: T::zero(), samples: [T::zero(); N] }
    }
}

impl<T, const N: usize> MovingAverageFilter<T, N>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<f32, Output = T>,
{
    pub fn apply(&mut self, input: T) -> T {
        self.sum = self.sum + input;
        if self.count < N {
            self.samples[self.index] = input;
            self.index += 1;
            self.count += 1;
            return self.sum * (1.0 / self.count as f32);
        }
        if self.index == N {
            self.index = 0;
        }
        self.sum = self.sum - self.samples[self.index];
        self.samples[self.index] = input;
        self.index += 1;

        self.sum * (1.0 / N as f32)
    }

    pub fn reset(&mut self) {
        self.sum = T::zero();
        self.count = 0;
        self.index = 0;
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
        is_full::<BiquadFilter<f32, f32>>();
        is_full::<BiquadFilterf32<f32>>();
        is_full::<BiquadFilterState<f32>>();
        is_full::<MovingAverageFilter<f32, 2>>();
    }
    #[test]
    fn pt1_filter_f32() {
        let mut filter = Pt1Filterf32::<f32>::new(1.0);

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
    fn moving_average_filter_f32() {
        let mut filter = MovingAverageFilter::<f32, 3>::new();
        assert_eq!(1.0, filter.apply(1.0));
        assert_eq!(1.5, filter.apply(2.0));
        assert_eq!(2.0, filter.apply(3.0));
        assert_eq!(3.0, filter.apply(4.0));
        assert_eq!(4.0, filter.apply(5.0));
        assert_eq!(5.0, filter.apply(6.0));
        assert_eq!(7.0, filter.apply(10.0));

        filter.reset();
        assert_eq!(4.0, filter.apply(4.0));
        assert_eq!(12.0, filter.apply(20.0));
        assert_eq!(5.0, filter.apply(-9.0));
    }
    #[test]
    fn pt1_filter_vector3df32() {
        let mut filter = Pt1Filter::<Vector3df32, f32>::new(1.0);
        let mut output: Vector3df32;
        let mut state: Vector3df32;

        // test that filter with default settings performs no filtering
        output = filter.apply(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 });
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, output);
        state = filter.state();
        assert_eq!(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 }, state);

        filter.reset();
        state = filter.state();
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
    #[test]
    fn moving_average_filter_vector3df32() {
        let mut filter = MovingAverageFilter::<Vector3df32, 4>::new();
        let mut m = filter.apply(Vector3df32 { x: 1.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 1.0, y: 0.0, z: -3.0 }, m);

        m = filter.apply(Vector3df32 { x: 2.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 1.5, y: 0.0, z: -3.0 }, m);

        m = filter.apply(Vector3df32 { x: 3.0, y: 3.0, z: 0.0 });
        assert_eq!(Vector3df32 { x: 2.0, y: 1.0, z: -2.0 }, m);

        m = filter.apply(Vector3df32 { x: 4.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 2.5, y: 1.25, z: -2.25 }, m);

        m = filter.apply(Vector3df32 { x: 5.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 3.5, y: 1.75, z: -2.25 }, m);

        m = filter.apply(Vector3df32 { x: 6.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 4.5, y: 2.25, z: -2.25 }, m);

        m = filter.apply(Vector3df32 { x: 10.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 6.25, y: 2.0, z: -3.0 }, m);

        filter.reset();
        m = filter.apply(Vector3df32 { x: 4.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 4.0, y: 2.0, z: -3.0 }, m);

        m = filter.apply(Vector3df32 { x: 20.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 12.0, y: 1.0, z: -3.0 }, m);

        m = filter.apply(Vector3df32 { x: -9.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 5.0, y: 2.0 / 3.0, z: -3.0 }, m);
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
    fn moving_average_filter_vector3df32_i16() {
        let mut filter = MovingAverageFilter::<Vector3di16, 4>::new();
        let mut m = filter.apply(Vector3di16 { x: 4, y: 0, z: -12 });
        assert_eq!(Vector3di16 { x: 4, y: 0, z: -12 }, m);

        m = filter.apply(Vector3di16 { x: 8, y: 0, z: -12 });
        assert_eq!(Vector3di16 { x: 6, y: 0, z: -12 }, m);

        m = filter.apply(Vector3di16 { x: 12, y: 12, z: 0 });
        assert_eq!(Vector3di16 { x: 8, y: 4, z: -8 }, m);

        m = filter.apply(Vector3di16 { x: 16, y: 8, z: -12 });
        assert_eq!(Vector3di16 { x: 10, y: 5, z: -9 }, m);
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
    #[test]
    fn moving_average_filter_vector3df32_i32() {
        let mut filter = MovingAverageFilter::<Vector3di32, 4>::new();
        let mut m = filter.apply(Vector3di32 { x: 4, y: 0, z: -12 });
        assert_eq!(Vector3di32 { x: 4, y: 0, z: -12 }, m);

        m = filter.apply(Vector3di32 { x: 8, y: 0, z: -12 });
        assert_eq!(Vector3di32 { x: 6, y: 0, z: -12 }, m);

        m = filter.apply(Vector3di32 { x: 12, y: 12, z: 0 });
        assert_eq!(Vector3di32 { x: 8, y: 4, z: -8 }, m);

        m = filter.apply(Vector3di32 { x: 16, y: 8, z: -12 });
        assert_eq!(Vector3di32 { x: 10, y: 5, z: -9 }, m);
    }
}
