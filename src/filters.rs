#![allow(unused)]

use core::f32::consts;
use core::ops::{Add, Mul, Sub};
use num_traits::Zero;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilterPT1<T> {
    state: T,
    k: f32,
}

/// Default is k = 1.0, which is passthrough
impl<T> Default for FilterPT1<T>
where
    T: Zero,
{
    fn default() -> Self {
        Self {
            state: T::zero(),
            k: 1.0,
        }
    }
}

impl<T> FilterPT1<T>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<f32, Output = T>,
{
    pub fn new(k: f32) -> Self {
        Self {
            state: T::zero(),
            k,
        }
    }

    pub fn init(&mut self, k: f32) {
        self.k = k;
        self.reset();
    }

    pub fn reset(&mut self) {
        self.state = T::zero();
    }

    pub fn filter(&mut self, input: T) -> T {
        self.state = self.state + (input - self.state) * self.k; // equivalent to self.state = self.k*input + (1.0 - self.k)*self.state;
        self.state
    }

    pub fn set_to_passthrough(&mut self) {
        self.k = 1.0;
        self.reset();
    }

    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: f32, delta_t: f32) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: f32, delta_t: f32) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    // Calculates filter gain based on delay (time constant of filter) - time it takes for filter response to reach 63.2% of a step input.
    pub fn gain_from_delay(delay: f32, delta_t: f32) -> f32 {
        if delay <= 0.0 {
            // gain of 1.0 means no filtering
            return 1.0;
        }
        let omega = delta_t / delay;
        omega / (omega + 1.0)
    }

    pub fn gain_from_frequency(cutoff_frequency_hz: f32, delta_t: f32) -> f32 {
        let omega = 2.0 * consts::PI * cutoff_frequency_hz * delta_t;
        omega / (omega + 1.0)
    }

    // for testing
    fn state(self) -> T {
        self.state
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilterPT2<T> {
    state: [T; 2],
    k: f32,
}

/// Default is k = 1.0, which is passthrough
impl<T> Default for FilterPT2<T>
where
    T: Zero,
{
    fn default() -> Self {
        Self {
            state: [T::zero(), T::zero()],
            k: 1.0,
        }
    }
}

impl<T> FilterPT2<T>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<f32, Output = T>,
{
    // PowerTransferFilter<n> cutoff correction = 1/sqrt(2^(1/n) - 1)
    #[allow(clippy::excessive_precision)]
    const CUTOFF_CORRECTION: f32 = 1.553773974;
    pub fn new(k: f32) -> Self {
        Self {
            state: [T::zero(), T::zero()],
            k,
        }
    }

    pub fn init(&mut self, k: f32) {
        self.k = k;
        self.reset();
    }

    pub fn reset(&mut self) {
        self.state = [T::zero(), T::zero()];
    }

    pub fn set_to_passthrough(&mut self) {
        self.k = 1.0;
        self.reset();
    }

    pub fn filter(&mut self, input: T) -> T {
        self.state[1] = self.state[1] + (input - self.state[1]) * self.k;
        self.state[0] = self.state[0] + (self.state[1] - self.state[0]) * self.k;
        self.state[0]
    }

    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: f32, delta_t: f32) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: f32, delta_t: f32) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    pub fn gain_from_delay(delay: f32, delta_t: f32) -> f32 {
        FilterPT1::<T>::gain_from_delay(delay * Self::CUTOFF_CORRECTION, delta_t)
    }
    pub fn gain_from_frequency(cutoff_frequency_hz: f32, delta_t: f32) -> f32 {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        FilterPT1::<T>::gain_from_frequency(cutoff_frequency_hz * Self::CUTOFF_CORRECTION, delta_t)
    }
    // for testing
    fn state(self) -> [T; 2] {
        self.state
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilterPT3<T> {
    state: [T; 3],
    k: f32,
}

/// Default is k = 1.0, which is passthrough
impl<T> Default for FilterPT3<T>
where
    T: Zero,
{
    fn default() -> Self {
        Self {
            state: [T::zero(), T::zero(), T::zero()],
            k: 1.0,
        }
    }
}

impl<T> FilterPT3<T>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<f32, Output = T>,
{
    // PowerTransferFilter<n> cutoff correction = 1/sqrt(2^(1/n) - 1)
    #[allow(clippy::excessive_precision)]
    const CUTOFF_CORRECTION: f32 = 1.961459177;
    pub fn new(k: f32) -> Self {
        Self {
            state: [T::zero(), T::zero(), T::zero()],
            k,
        }
    }

    pub fn init(&mut self, k: f32) {
        self.k = k;
        self.reset();
    }

    pub fn reset(&mut self) {
        self.state = [T::zero(), T::zero(), T::zero()];
    }

    pub fn set_to_passthrough(&mut self) {
        self.k = 1.0;
        self.reset();
    }

    pub fn filter(&mut self, input: T) -> T {
        self.state[2] = self.state[2] + (input - self.state[2]) * self.k;
        self.state[1] = self.state[1] + (self.state[2] - self.state[1]) * self.k;
        self.state[0] = self.state[0] + (self.state[1] - self.state[0]) * self.k;
        self.state[0]
    }

    pub fn set_cutoff_frequency(&mut self, cutoff_frequency_hz: f32, delta_t: f32) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
    }

    pub fn set_cutoff_frequency_and_reset(&mut self, cutoff_frequency_hz: f32, delta_t: f32) {
        self.k = Self::gain_from_frequency(cutoff_frequency_hz, delta_t);
        self.reset();
    }

    pub fn gain_from_delay(delay: f32, delta_t: f32) -> f32 {
        FilterPT1::<T>::gain_from_delay(delay * Self::CUTOFF_CORRECTION, delta_t)
    }

    pub fn gain_from_frequency(cutoff_frequency_hz: f32, delta_t: f32) -> f32 {
        // shift cutoffFrequency to satisfy -3dB cutoff condition
        FilterPT1::<T>::gain_from_frequency(cutoff_frequency_hz * Self::CUTOFF_CORRECTION, delta_t)
    }

    // for testing
    fn state(self) -> [T; 3] {
        self.state
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
        Self {
            x1: T::zero(),
            x2: T::zero(),
            y1: T::zero(),
            y2: T::zero(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BiquadFilter<T> {
    state: BiquadFilterState<T>,
    weight: f32,
    a1: f32,
    a2: f32,
    b0: f32,
    b1: f32,
    b2: f32,
    _2_pi_loop_time_seconds: f32,
    _2_q_reciprocal: f32,
}

impl<T> Default for BiquadFilter<T>
where
    T: Zero,
{
    fn default() -> Self {
        Self {
            state: BiquadFilterState {
                x1: T::zero(),
                x2: T::zero(),
                y1: T::zero(),
                y2: T::zero(),
            },
            weight: 1.0,
            a1: 0.0,
            a2: 0.0,
            b0: 1.0,
            b1: 0.0,
            b2: 0.0,
            _2_pi_loop_time_seconds: 0.0,
            _2_q_reciprocal: 0.0,
        }
    }
}

impl<T> BiquadFilter<T>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<f32, Output = T>,
{
    pub fn set_weight(&mut self, weight: f32) {
        self.weight = weight;
    }
    pub fn weight(&self) -> f32 {
        self.weight
    }
    pub fn set_parameters_and_weight(
        &mut self,
        a1: f32,
        a2: f32,
        b0: f32,
        b1: f32,
        b2: f32,
        weight: f32,
    ) {
        self.weight = weight;
        self.a1 = a1;
        self.a2 = a2;
        self.b0 = b0;
        self.b1 = b1;
        self.b2 = b2;
    }
    pub fn set_parameters(&mut self, a1: f32, a2: f32, b0: f32, b1: f32, b2: f32) {
        self.set_parameters_and_weight(a1, a2, b0, b1, b2, 1.0);
    }
    /// Copy parameters from another Biquad filter
    pub fn set_parameters_from(&mut self, other: &BiquadFilter<T>) {
        self.weight = other.weight;
        self.a1 = other.a1;
        self.a2 = other.a2;
        self.b0 = other.b0;
        self.b1 = other.b1;
        self.b2 = other.b2;
    }

    pub fn reset(&mut self) {
        self.state.x1 = T::zero();
        self.state.x2 = T::zero();
        self.state.y1 = T::zero();
        self.state.y2 = T::zero();
    }
    pub fn set_to_passthrough(&mut self) {
        self.b0 = 1.0;
        self.b1 = 0.0;
        self.b2 = 0.0;
        self.a1 = 0.0;
        self.a2 = 0.0;
        self.weight = 1.0;
        self.reset();
    }

    pub fn filter(&mut self, input: T) -> T {
        let output = input * self.b0 + self.state.x1 * self.b1 + self.state.x2 * self.b2
            - self.state.y1 * self.a1
            - self.state.y2 * self.a2;

        self.state.x2 = self.state.x1;
        self.state.x1 = input;
        self.state.y2 = self.state.y1;
        self.state.y1 = output;
        output
    }

    pub fn filter_weighted(&mut self, input: T) -> T {
        let output = self.filter(input);
        // weight of 1.0 gives just output, weight of 0.0 gives just input
        (output - input) * self.weight + input
    }

    pub fn init_low_pass(&mut self, frequency_hz: f32, loop_time_seconds: f32, q: f32) {
        //assert(Q != 0.0 && "Q cannot be zero");
        self.set_loop_time(loop_time_seconds);
        self.set_q(q);
        self.set_low_pass_frequency(frequency_hz);
        self.reset();
    }

    pub fn init_notch(&mut self, frequency_hz: f32, loop_time_seconds: f32, q: f32) {
        //assert(Q != 0.0 && "Q cannot be zero");
        self.set_loop_time(loop_time_seconds);
        self.set_q(q);
        //setNotchFrequency(frequency_hz);
        self.reset();
    }

    pub fn calculate_omega(&self, frequency: f32) -> f32 {
        frequency * self._2_pi_loop_time_seconds
    }

    //Note: weight must be in range [0, 1].
    pub fn set_low_pass_frequency_weighted(&mut self, frequency_hz: f32, weight: f32) {
        self.weight = weight;

        let omega = frequency_hz * self._2_pi_loop_time_seconds;
        let (sin_omega, cos_omega) = omega.sin_cos();
        let alpha = sin_omega * self._2_q_reciprocal;
        let a0_reciprocal = 1.0 / (1.0 + alpha);

        self.b1 = (1.0 - cos_omega) * a0_reciprocal;
        self.b0 = self.b1 * 0.5;
        self.b2 = self.b0;
        self.a1 = -2.0 * cos_omega * a0_reciprocal;
        self.a2 = (1.0 - alpha) * a0_reciprocal;
    }

    pub fn set_low_pass_frequency(&mut self, frequency_hz: f32) {
        self.set_low_pass_frequency_weighted(frequency_hz, 1.0);
    }

    pub fn set_notch_frequency_weighted_from_sin_cos_assuming_q(
        &mut self,
        sin_omega: f32,
        cos_omega: f32,
        weight: f32,
    ) {
        self.weight = weight;

        let alpha = sin_omega * self._2_q_reciprocal;
        let a0reciprocal = 1.0 / (1.0 + alpha);

        self.b0 = a0reciprocal;
        self.b2 = a0reciprocal;
        self.b1 = -2.0 * cos_omega * a0reciprocal;
        self.a1 = self.b1;
        self.a2 = (1.0 - alpha) * a0reciprocal;
    }

    pub fn set_notch_frequency_weighted_assuming_q(&mut self, frequency_hz: f32, weight: f32) {
        let omega = frequency_hz * self._2_pi_loop_time_seconds;
        let (sin_omega, cos_omega) = omega.sin_cos();
        self.set_notch_frequency_weighted_from_sin_cos_assuming_q(sin_omega, cos_omega, weight);
    }

    pub fn set_notch_frequency_assuming_q(&mut self, frequency_hz: f32) {
        // assumes Q already set
        self.set_notch_frequency_weighted_assuming_q(frequency_hz, 1.0);
    }

    pub fn set_notch_frequency(
        &mut self,
        center_frequency_hz: f32,
        lower_cutoff_frequency_hz: f32,
    ) {
        self.set_q(Self::calculate_q(
            center_frequency_hz,
            lower_cutoff_frequency_hz,
        ));
        self.set_notch_frequency_assuming_q(center_frequency_hz);
    }

    pub fn calculate_q(center_frequency_hz: f32, lower_cutoff_frequency_hz: f32) -> f32 {
        center_frequency_hz * lower_cutoff_frequency_hz
            / (center_frequency_hz * center_frequency_hz
                - lower_cutoff_frequency_hz * lower_cutoff_frequency_hz)
    }

    pub fn set_q_from_frequencies(
        &mut self,
        center_frequency_hz: f32,
        lower_cutoff_frequency_hz: f32,
    ) {
        self._2_q_reciprocal =
            1.0 / (2.0 * Self::calculate_q(center_frequency_hz, lower_cutoff_frequency_hz));
    }

    pub fn set_q(&mut self, q: f32) {
        self._2_q_reciprocal = 1.0 / (2.0 * q);
    }

    pub fn q_calculated(&self) -> f32 {
        (1.0 / self._2_q_reciprocal) / 2.0
    }

    pub fn set_loop_time(&mut self, loop_time_seconds: f32) {
        self._2_pi_loop_time_seconds = 2.0 * consts::PI * loop_time_seconds;
    }

    // for testing
    fn state(self) -> BiquadFilterState<T> {
        self.state
    }
}

/// Simple moving average filter.
/// See [Moving Average Filter - Theory and Software Implementation - Phil's Lab #21](https://www.youtube.com/watch?v=rttn46_Y3c8).

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FilterMovingAverage<T, const N: usize> {
    count: usize,
    index: usize,
    sum: T,
    samples: [T; N],
}

impl<T, const N: usize> FilterMovingAverage<T, N>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<f32, Output = T>,
{
    fn new() -> Self {
        Self {
            count: 0,
            index: 0,
            sum: T::zero(),
            samples: [T::zero(); N],
        }
    }
    pub fn reset(&mut self) {
        self.sum = T::zero();
        self.count = 0;
        self.index = 0;
    }

    pub fn filter(&mut self, input: T) -> T {
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
}

#[cfg(any(debug_assertions, test))]
mod tests {
    use super::*;
    use vector_quaternion_matrix::Vector3d;
    use vector_quaternion_matrix::Vector3dI16;
    use vector_quaternion_matrix::Vector3dI32;

    #[test]
    fn filter_pt1_f32() {
        let mut filter = FilterPT1::<f32>::new(1.0);

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(-1.0, filter.filter(-1.0));

        filter.reset();
        assert_eq!(0.0, filter.state());
        assert_eq!(4.0, filter.filter(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.38586956, filter.filter(1.0));
        assert_eq!(1.0087134, filter.filter(2.0));

        filter.init(1.0);
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.38586956, filter.filter(1.0));
        assert_eq!(1.0087134, filter.filter(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));
    }
    #[test]
    fn filter_pt2_f32() {
        let mut filter = FilterPT2::<f32>::new(1.0);

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(-1.0, filter.filter(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.filter(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.24403107, filter.filter(1.0));
        assert_eq!(0.73502403, filter.filter(2.0));

        filter.init(1.0);
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.24403107, filter.filter(1.0));
        assert_eq!(0.73502403, filter.filter(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));
    }
    #[test]
    fn filter_pt3_f32() {
        let mut filter = FilterPT3::<f32>::new(1.0);

        let mut state = filter.state();
        assert_eq!([0.0, 0.0, 0.0], state);

        // test that filter with default settings performs no filtering
        let mut output = filter.filter(1.0);
        assert_eq!(1.0, output);
        state = filter.state();
        assert_eq!([1.0, 1.0, 1.0], state);

        output = filter.filter(1.0);
        state = filter.state();
        assert_eq!([1.0, 1.0, 1.0], state);
        assert_eq!(1.0, output);

        assert_eq!(-1.0, filter.filter(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.filter(4.0));

        filter.reset();
        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(0.16824766, filter.filter(1.0));
        assert_eq!(0.56259197, filter.filter(2.0));

        filter.init(1.0);
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(0.16824766, filter.filter(1.0));
        assert_eq!(0.56259197, filter.filter(2.0));

        filter.set_to_passthrough();
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));
    }
    #[test]
    fn biquad_filter_f32() {
        let mut filter = BiquadFilter::<f32>::default();

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(-1.0, filter.filter(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.filter(4.0));

        filter.set_parameters_and_weight(2.0, 3.0, 5.0, 7.0, 11.0, 13.0);
        filter.set_to_passthrough();
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(2.0, filter.filter(2.0));
        assert_eq!(1.0, filter.filter_weighted(1.0));
        assert_eq!(2.0, filter.filter_weighted(2.0));
    }
    #[test]
    fn moving_average_filter_f32() {
        let mut filter = FilterMovingAverage::<f32, 3>::new();
        assert_eq!(1.0, filter.filter(1.0));
        assert_eq!(1.5, filter.filter(2.0));
        assert_eq!(2.0, filter.filter(3.0));
        assert_eq!(3.0, filter.filter(4.0));
        assert_eq!(4.0, filter.filter(5.0));
        assert_eq!(5.0, filter.filter(6.0));
        assert_eq!(7.0, filter.filter(10.0));

        filter.reset();
        assert_eq!(4.0, filter.filter(4.0));
        assert_eq!(12.0, filter.filter(20.0));
        assert_eq!(5.0, filter.filter(-9.0));
    }
    #[test]
    fn filter_pt1_vector3d() {
        let mut filter = FilterPT1::<Vector3d>::new(1.0);
        let mut output: Vector3d;
        let mut state: Vector3d;

        // test that filter with default settings performs no filtering
        output = filter.filter(Vector3d {
            x: 2.0,
            y: 3.0,
            z: 5.0,
        });
        assert_eq!(
            Vector3d {
                x: 2.0,
                y: 3.0,
                z: 5.0
            },
            output
        );
        state = filter.state();
        assert_eq!(
            Vector3d {
                x: 2.0,
                y: 3.0,
                z: 5.0
            },
            state
        );

        filter.reset();
        state = filter.state();
        assert_eq!(
            Vector3d {
                x: 0.0,
                y: 0.0,
                z: 0.0
            },
            state
        );

        filter.set_cutoff_frequency(100.0, 0.001);
        assert_eq!(
            0.38586956,
            filter
                .filter(Vector3d {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            1.0087134,
            filter
                .filter(Vector3d {
                    x: 2.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );

        filter.init(1.0);
        assert_eq!(
            1.0,
            filter
                .filter(Vector3d {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            2.0,
            filter
                .filter(Vector3d {
                    x: 2.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );

        filter.set_cutoff_frequency_and_reset(100.0, 0.001);
        assert_eq!(
            0.38586956,
            filter
                .filter(Vector3d {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            1.0087134,
            filter
                .filter(Vector3d {
                    x: 2.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );

        filter.set_to_passthrough();
        assert_eq!(
            1.0,
            filter
                .filter(Vector3d {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            2.0,
            filter
                .filter(Vector3d {
                    x: 2.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
    }
    #[test]
    fn biquad_filter_vector3d() {
        let mut filter = BiquadFilter::<Vector3d>::default();
        let mut output: Vector3d;
        let mut state: BiquadFilterState<Vector3d>;

        // test that filter with default settings performs no filtering
        output = filter.filter(Vector3d {
            x: 2.0,
            y: 3.0,
            z: 5.0,
        });
        assert_eq!(
            Vector3d {
                x: 2.0,
                y: 3.0,
                z: 5.0
            },
            output
        );
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
        assert_eq!(
            4.0,
            filter
                .filter(Vector3d {
                    x: 4.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );

        filter.set_parameters_and_weight(2.0, 3.0, 5.0, 7.0, 11.0, 13.0);
        filter.set_to_passthrough();
        assert_eq!(
            1.0,
            filter
                .filter(Vector3d {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            2.0,
            filter
                .filter(Vector3d {
                    x: 2.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            1.0,
            filter
                .filter_weighted(Vector3d {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
        assert_eq!(
            2.0,
            filter
                .filter_weighted(Vector3d {
                    x: 2.0,
                    y: 0.0,
                    z: 0.0
                })
                .x
        );
    }
    #[test]
    fn moving_average_filter_vector3d() {
        let mut filter = FilterMovingAverage::<Vector3d, 4>::new();
        let mut m = filter.filter(Vector3d {
            x: 1.0,
            y: 0.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 1.0,
                y: 0.0,
                z: -3.0
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 2.0,
            y: 0.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 1.5,
                y: 0.0,
                z: -3.0
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 3.0,
            y: 3.0,
            z: 0.0,
        });
        assert_eq!(
            Vector3d {
                x: 2.0,
                y: 1.0,
                z: -2.0
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 4.0,
            y: 2.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 2.5,
                y: 1.25,
                z: -2.25
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 5.0,
            y: 2.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 3.5,
                y: 1.75,
                z: -2.25
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 6.0,
            y: 2.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 4.5,
                y: 2.25,
                z: -2.25
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 10.0,
            y: 2.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 6.25,
                y: 2.0,
                z: -3.0
            },
            m
        );

        filter.reset();
        m = filter.filter(Vector3d {
            x: 4.0,
            y: 2.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 4.0,
                y: 2.0,
                z: -3.0
            },
            m
        );

        m = filter.filter(Vector3d {
            x: 20.0,
            y: 0.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 12.0,
                y: 1.0,
                z: -3.0
            },
            m
        );

        m = filter.filter(Vector3d {
            x: -9.0,
            y: 0.0,
            z: -3.0,
        });
        assert_eq!(
            Vector3d {
                x: 5.0,
                y: 2.0 / 3.0,
                z: -3.0
            },
            m
        );
    }
    #[test]
    fn filter_pt1_vector3d_i16() {
        let mut filter = FilterPT1::<Vector3dI16>::new(1.0);
        let mut output: Vector3dI16;
        let mut state: Vector3dI16;

        // test that filter with default settings performs no filtering
        output = filter.filter(Vector3dI16 { x: 2, y: 3, z: 5 });
        assert_eq!(Vector3dI16 { x: 2, y: 3, z: 5 }, output);
        state = filter.state();
        assert_eq!(Vector3dI16 { x: 2, y: 3, z: 5 }, state);
    }
    #[test]
    fn moving_average_filter_vector3d_i16() {
        let mut filter = FilterMovingAverage::<Vector3dI16, 4>::new();
        let mut m = filter.filter(Vector3dI16 { x: 4, y: 0, z: -12 });
        assert_eq!(Vector3dI16 { x: 4, y: 0, z: -12 }, m);

        m = filter.filter(Vector3dI16 { x: 8, y: 0, z: -12 });
        assert_eq!(Vector3dI16 { x: 6, y: 0, z: -12 }, m);

        m = filter.filter(Vector3dI16 { x: 12, y: 12, z: 0 });
        assert_eq!(Vector3dI16 { x: 8, y: 4, z: -8 }, m);

        m = filter.filter(Vector3dI16 {
            x: 16,
            y: 8,
            z: -12,
        });
        assert_eq!(Vector3dI16 { x: 10, y: 5, z: -9 }, m);
    }
    #[test]
    fn filter_pt1_vector3d_i32() {
        let mut filter = FilterPT1::<Vector3dI32>::new(1.0);
        let mut output: Vector3dI32;
        let mut state: Vector3dI32;

        // test that filter with default settings performs no filtering
        output = filter.filter(Vector3dI32 { x: 2, y: 3, z: 5 });
        assert_eq!(2, output.x);
        assert_eq!(3, output.y);
        assert_eq!(5, output.z);
        state = filter.state();
        assert_eq!(2, state.x);
        assert_eq!(3, state.y);
        assert_eq!(5, state.z);
    }
    #[test]
    fn moving_average_filter_vector3d_i32() {
        let mut filter = FilterMovingAverage::<Vector3dI32, 4>::new();
        let mut m = filter.filter(Vector3dI32 { x: 4, y: 0, z: -12 });
        assert_eq!(Vector3dI32 { x: 4, y: 0, z: -12 }, m);

        m = filter.filter(Vector3dI32 { x: 8, y: 0, z: -12 });
        assert_eq!(Vector3dI32 { x: 6, y: 0, z: -12 }, m);

        m = filter.filter(Vector3dI32 { x: 12, y: 12, z: 0 });
        assert_eq!(Vector3dI32 { x: 8, y: 4, z: -8 }, m);

        m = filter.filter(Vector3dI32 {
            x: 16,
            y: 8,
            z: -12,
        });
        assert_eq!(Vector3dI32 { x: 10, y: 5, z: -9 }, m);
    }
}
