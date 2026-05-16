use core::ops::{Add, Div, Mul, Neg, Sub};
use num_traits::{ConstOne, ConstZero, One, Zero};
use vqm::{MathConstants, TrigonometricMethods, Vector2d, Vector3d, Vector4d};

use crate::SignalFilter;

/// `BiquadFilter` for `f32`<br>
pub type BiquadFilterf32 = BiquadFilter<f32, f32>;
/// `BiquadFilter` for `Vector2df32`<br>
pub type BiquadFilterVector2df32 = BiquadFilter<Vector2d<f32>, f32>;
/// `BiquadFilter` for `Vector3df32`<br>
pub type BiquadFilterVector3df32 = BiquadFilter<Vector3d<f32>, f32>;
/// `BiquadFilter` for `Vector4df32`<br>
pub type BiquadFilterVector4df32 = BiquadFilter<Vector4d<f32>, f32>;

/// `BiquadFilter` for `f64`<br><br>
pub type BiquadFilterf64 = BiquadFilter<f64, f64>;
/// `BiquadFilter` for `Vector2df64`<br><br>
pub type BiquadFilterVector2df64 = BiquadFilter<Vector2d<f64>, f64>;
/// `BiquadFilter` for `Vector3df64`<br><br>
pub type BiquadFilterVector3df64 = BiquadFilter<Vector3d<f64>, f64>;
/// `BiquadFilter` for `Vector4df64`<br><br>
pub type BiquadFilterVector4df64 = BiquadFilter<Vector4d<f64>, f64>;

pub trait ConstHalf {
    const HALF: Self;
}

impl ConstHalf for f32 {
    const HALF: f32 = 0.5;
}

impl ConstHalf for f64 {
    const HALF: f64 = 0.5;
}

#[allow(clippy::doc_paragraphs_missing_punctuation)]
/// Second-order biquad IIR filter.<br><br>
///
/// This implementation uses the Direct Form I structure.
///
/// The difference equation is:
///
/// ```math
/// {n} = b{0} * x{n} + b{1} * x{n-1} + b{2} * x{n-2} - a{1} * y{n-1} - a{2} * y{n-2}
/// ```
///
/// where `x` represents the input signal and `y` represents the filtered output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BiquadFilter<T, R> {
    state: BiquadFilterState<T>,
    coeffs: BiquadFilterCoefficients<R>,
    weight: R,
    loop_time_seconds: R,
    two_pi_loop_time_seconds: R, // cached value of 2.0 * PI * loop_time_seconds
    q: R,
    one_over_2q: R, // cached value of 1.0 / (2.0 * q)
}

impl<T, R> Default for BiquadFilter<T, R>
where
    T: Copy + ConstZero,
    R: Copy + ConstZero + ConstOne + ConstHalf,
{
    fn default() -> Self {
        Self::new()
    }
}
impl<T, R> BiquadFilter<T, R>
where
    T: Copy + ConstZero,
    R: Copy + ConstZero + ConstOne + ConstHalf,
{
    pub const fn new() -> Self {
        Self {
            state: BiquadFilterState::new(),
            coeffs: BiquadFilterCoefficients::new(),
            weight: R::ONE,
            loop_time_seconds: R::ZERO,
            two_pi_loop_time_seconds: R::ZERO,
            q: R::ONE,
            one_over_2q: R::HALF,
        }
    }
    pub fn with_coefficients(coeffs: BiquadFilterCoefficients<R>) -> Self {
        Self {
            state: BiquadFilterState::new(),
            coeffs,
            weight: R::ONE,
            loop_time_seconds: R::ZERO,
            two_pi_loop_time_seconds: R::ZERO,
            q: R::ONE,
            one_over_2q: R::HALF,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BiquadFilterCoefficients<R> {
    pub a1: R,
    pub a2: R,
    pub b0: R,
    pub b1: R,
    pub b2: R,
}

impl<R> BiquadFilterCoefficients<R>
where
    R: Copy + ConstZero + ConstOne,
{
    pub const fn new() -> Self {
        Self { a1: R::ZERO, a2: R::ZERO, b0: R::ONE, b1: R::ZERO, b2: R::ZERO }
    }
}

impl<T> Default for BiquadFilterCoefficients<T>
where
    T: Copy + ConstZero + ConstOne,
{
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct BiquadFilterState<T> {
    // Input history
    x1: T,
    x2: T,
    // Output
    y1: T,
    y2: T,
}

impl<T> BiquadFilterState<T>
where
    T: Copy + ConstZero,
{
    pub const fn new() -> Self {
        Self { x1: T::ZERO, x2: T::ZERO, y1: T::ZERO, y2: T::ZERO }
    }
}

impl<T> Default for BiquadFilterState<T>
where
    T: Copy + ConstZero,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, R> SignalFilter<T, R> for BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy,
{
    fn reset(&mut self) {
        self.state.x1 = T::zero();
        self.state.x2 = T::zero();
        self.state.y1 = T::zero();
        self.state.y2 = T::zero();
    }

    fn update(&mut self, input: T) -> T {
        // 9 operations: 5 multiplications, 4 additions
        let output = input * self.coeffs.b0 + self.state.x1 * self.coeffs.b1 + self.state.x2 * self.coeffs.b2
            - self.state.y1 * self.coeffs.a1
            - self.state.y2 * self.coeffs.a2;

        self.state.x2 = self.state.x1;
        self.state.x1 = input;
        self.state.y2 = self.state.y1;
        self.state.y1 = output;
        output
    }
}

/// Biquad update notch.
///
/// For a notch filter a1 == b1 and b0 == b2, so this optimized form can be used.
impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Zero + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy,
{
    pub fn update_notch(&mut self, input: T) -> T {
        // 8 operations: 3 multiplications, 5 additions
        let output = (input + self.state.x2) * self.coeffs.b0
            + (self.state.x1 + self.state.x1 - self.state.y1) * self.coeffs.a1
            - self.state.y2 * self.coeffs.a2;

        self.state.x2 = self.state.x1;
        self.state.x1 = input;
        self.state.y2 = self.state.y1;
        self.state.y1 = output;
        output
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One + Div<R, Output = R>,
{
    pub fn set_q(&mut self, q: R) {
        self.q = q;
        self.one_over_2q = R::one() / ((R::one() + R::one()) * q); // cache value for faster setting of frequencies
    }
}

impl<T, R> BiquadFilter<T, R>
where
    T: Copy + Add<Output = T> + Sub<Output = T> + Mul<R, Output = T>,
    R: Copy + Zero + One,
{
    pub fn set_weight(&mut self, weight: R) {
        self.weight = weight;
    }

    pub fn weight(&self) -> R {
        self.weight
    }

    pub fn set_coefficients(&mut self, coeffs: BiquadFilterCoefficients<R>) {
        self.coeffs = coeffs;
    }

    pub fn set_coefficients_and_weight(&mut self, coeffs: BiquadFilterCoefficients<R>, weight: R) {
        self.weight = weight;
        self.coeffs = coeffs;
    }

    /// Copy parameters from another Biquad filter.
    pub fn set_parameters_from(&mut self, other: &BiquadFilter<T, R>) {
        self.weight = other.weight;
        self.coeffs = other.coeffs;
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
    R: Copy
        + Zero
        + One
        + Neg<Output = R>
        + MathConstants
        + TrigonometricMethods
        + Div<R, Output = R>
        + Sub<R, Output = R>,
{
    pub fn set_to_passthrough(&mut self) {
        self.coeffs.b0 = R::one();
        self.coeffs.b1 = R::zero();
        self.coeffs.b2 = R::zero();
        self.coeffs.a1 = R::zero();
        self.coeffs.a2 = R::zero();
        self.weight = R::one();
        self.reset();
    }

    pub fn update_weighted(&mut self, input: T) -> T {
        let output = self.update(input);
        // weight of 1.0 gives just output, weight of 0.0 gives just input
        (output - input) * self.weight + input
    }

    pub fn update_notch_weighted(&mut self, input: T) -> T {
        let output = self.update_notch(input);
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

        self.coeffs.b1 = (R::one() - cos_omega) * a0_reciprocal;
        self.coeffs.b0 = self.coeffs.b1 * (R::one() / (R::one() + R::one()));
        self.coeffs.b2 = self.coeffs.b0;
        self.coeffs.a1 = -(R::one() + R::one()) * cos_omega * a0_reciprocal;
        self.coeffs.a2 = (R::one() - alpha) * a0_reciprocal;
    }

    pub fn set_low_pass_frequency_assuming_q(&mut self, frequency_hz: R) {
        self.set_low_pass_frequency_weighted_assuming_q(frequency_hz, R::one());
    }

    pub fn set_notch_frequency_weighted_from_sin_cos_assuming_q(&mut self, sin_omega: R, cos_omega: R, weight: R) {
        self.weight = weight;

        let alpha = sin_omega * self.one_over_2q;
        let a0reciprocal = R::one() / (R::one() + alpha);
        // NOTE: b0 == b2 and a1 == b1 for notch filter
        self.coeffs.b0 = a0reciprocal;
        self.coeffs.b2 = a0reciprocal;
        self.coeffs.b1 = R::zero() - (R::one() + R::one()) * cos_omega * a0reciprocal;
        self.coeffs.a1 = self.coeffs.b1;
        self.coeffs.a2 = (R::one() - alpha) * a0reciprocal;
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
    #![allow(clippy::float_cmp)]
    #[allow(unused)]
    use super::*;

    #[allow(unused)]
    fn is_normal<T: Sized + Send + Sync + Unpin>() {}
    #[allow(unused)]
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<BiquadFilter<f32, f32>>();
        is_full::<BiquadFilterf32>();
        is_full::<BiquadFilterState<f32>>();
    }
    #[test]
    fn biquad_filter_f32() {
        let mut filter = BiquadFilterf32::default();

        // test that filter with default settings performs no filtering
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(-1.0, filter.update(-1.0));

        filter.reset();
        assert_eq!(4.0, filter.update(4.0));

        filter.set_coefficients_and_weight(
            BiquadFilterCoefficients { a1: 2.0, a2: 3.0, b0: 5.0, b1: 7.0, b2: 11.0 },
            13.0,
        );
        filter.set_to_passthrough();
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(2.0, filter.update(2.0));
        assert_eq!(1.0, filter.update_weighted(1.0));
        assert_eq!(2.0, filter.update_weighted(2.0));
    }
    #[test]
    fn biquad_filter_vector3df32() {
        use vqm::Vector3df32;
        let mut filter = BiquadFilterVector3df32::default();
        let mut state: BiquadFilterState<Vector3df32>;

        // test that filter with default settings performs no filtering
        let output = filter.update(Vector3df32 { x: 2.0, y: 3.0, z: 5.0 });
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
        assert_eq!(4.0, filter.update(Vector3df32 { x: 4.0, y: 0.0, z: 0.0 }).x);

        filter.set_coefficients_and_weight(
            BiquadFilterCoefficients { a1: 2.0, a2: 3.0, b0: 5.0, b1: 7.0, b2: 11.0 },
            13.0,
        );
        filter.set_to_passthrough();
        assert_eq!(1.0, filter.update(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.update(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(1.0, filter.update_weighted(Vector3df32 { x: 1.0, y: 0.0, z: 0.0 }).x);
        assert_eq!(2.0, filter.update_weighted(Vector3df32 { x: 2.0, y: 0.0, z: 0.0 }).x);
    }
}
