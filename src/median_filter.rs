#![allow(unused)]
use num_traits::Zero;
use vector_quaternion_matrix::Vector2df32;

use crate::SignalFilter;

pub type Median3Filterf32 = Median3Filter<f32>;
pub type Median3FilterVector3df32<const N: usize> = Median3Filter<Vector2df32>;

pub type MedianFilterf32<const N: usize> = MedianFilter<f32, N>;
pub type MedianFilterVector3df32<const N: usize> = MedianFilter<Vector2df32, N>;

/// Non-linear median-of-3 filter for spike rejection.<br>
/// Maintains a window of the last three samples and returns the median value.
///
/// It is effective at removing single-sample outliers without "smearing"
/// the error into subsequent samples like a linear low-pass filter would.
///
/// The output $y_{n}$ is defined as:
///
/// $$y_{n} = \text{median}(x_{n}, x_{n-1}, x_{n-2})$$
///
/// **Note:** This filter introduces a fixed lag of 1 sample. During the
/// first two samples after a reset, the filter returns the raw input.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Median3Filter<T> {
    buffer: [T; 3],
    index: usize,
    count: usize,
}

impl<T> Default for Median3Filter<T>
where
    T: Copy + Zero,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Median3Filter<T>
where
    T: Copy + Zero,
{
    pub fn new() -> Self {
        Self { buffer: [T::zero(); 3], index: 0, count: 0 }
    }
}

pub struct Median3FilterTF<T, F> {
    buffer: [T; 3],
    index: usize,
    count: usize,
    dummy: F,
}

impl<T, F> Default for Median3FilterTF<T, F>
where
    T: Copy + Zero,
    F: Zero,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, F> Median3FilterTF<T, F>
where
    T: Copy + Zero,
    F: Zero,
{
    pub fn new() -> Self {
        Self { buffer: [T::zero(); 3], index: 0, count: 0, dummy: F::zero() }
    }
}

impl<T, F> SignalFilter<T, F> for Median3FilterTF<T, F>
where
    T: Copy + Zero + PartialOrd,
{
    fn reset(&mut self) {
        self.buffer = [T::zero(); 3];
        self.index = 0;
        self.count = 0;
    }

    fn update(&mut self, input: T) -> T {
        // Store new sample in the ring buffer
        self.buffer[self.index] = input;
        self.index = (self.index + 1) % 3;
        if self.count < 3 {
            self.count += 1;
        }

        // If buffer isn't full, just return the input
        if self.count < 3 {
            return input;
        }

        // Fast sorting network for 3 values (no loops)
        let mut a = self.buffer[0];
        let mut b = self.buffer[1];
        let mut c = self.buffer[2];

        if a > b {
            core::mem::swap(&mut a, &mut b);
        }
        if b > c {
            core::mem::swap(&mut b, &mut c);
        }
        if a > b {
            core::mem::swap(&mut a, &mut b);
        }

        // b is now the median
        b
    }
}

/// General purpose Moving Median (Window of N)
/// Maintains a window of the last `N` samples and returns the median value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MedianFilter<T, const N: usize> {
    buffer: [T; N],
    index: usize,
}

impl<T, const N: usize> Default for MedianFilter<T, N>
where
    T: Copy + Zero,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, const N: usize> MedianFilter<T, N>
where
    T: Copy + Zero,
{
    pub fn new() -> Self {
        Self { buffer: [T::zero(); N], index: 0 }
    }
}

impl<T, F, const N: usize> SignalFilter<T, F> for MedianFilter<T, N>
where
    T: Copy + Zero + PartialOrd,
    F: Copy,
{
    fn reset(&mut self) {
        self.buffer.fill(T::zero());
        self.index = 0;
    }

    fn update(&mut self, input: T) -> T {
        self.buffer[self.index] = input;
        self.index = (self.index + 1) % N;

        // Copy and sort to find the median without disturbing the ring buffer order
        let mut sorted = self.buffer;
        // no_std stable sort (insertion sort for small N)
        // TODO: sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));

        // Return the middle element
        sorted[N / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //use filters::

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<Median3Filter<f32>>();
        is_full::<MedianFilter<f32, 5>>();
    }
    #[test]
    fn test_median3_spike_rejection() {
        //let mut filter = MovingAverageFilter::<f32, 3>::new();
        //let mut filter: Median3Filterf32; // as SignalFilter<f32, f32>>;
        let mut filter = <Median3FilterTF<f32, f32>>::new();
        //let mut filter = Median3Filter<f32>::new();

        // 1. Initial values (filling the buffer)
        // Values: [10.0, 0.0, 0.0] -> Count 1 -> Returns 10.0
        let input = 10.0f32;
        //let output = <crate::median_filter::Median3Filter<f32> as crate::median_filter::SignalFilter<f32, f32>>::update(&mut filter, input);
        //let output = <crate::Median3Filter<f32> as crate::SignalFilter<f32, f32>>::update(&mut filter, input);
        //let output = <Median3Filter<f32> as SignalFilter<f32, f32>>::update(&mut filter, input);

        let output = filter.update(input);
        assert_eq!(input, output);
        // Values: [10.0, 20.0, 0.0] -> Count 2 -> Returns 20.0
        //assert_eq!(filter.update(20.0), 20.0);
        /*
        // 2. The Buffer is now full: [10.0, 20.0, 30.0]
        // Median of {10, 20, 30} is 20.0
        assert_eq!(filter.update(30.0), 20.0);

        // 3. Test a massive outlier "spike"
        // Buffer: [400.0, 20.0, 30.0] (400 replaces 10)
        // Median of {400, 20, 30} is 30.0 (The spike is ignored!)
        assert_eq!(filter.update(400.0), 30.0);

        // 4. Return to normal
        // Buffer: [400.0, 25.0, 30.0] (25 replaces 20)
        // Median of {400, 25, 30} is 30.0
        assert_eq!(filter.update(25.0), 30.0);

        // Buffer: [400.0, 25.0, 22.0] (22 replaces 30)
        // Median of {400, 25, 22} is 25.0
        assert_eq!(filter.update(22.0), 25.0);*/
    }

    /*#[test]
    fn test_median3_identical_values() {
        let mut filter = Median3Filter::new();
        filter.update(5.0);
        filter.update(5.0);
        assert_eq!(filter.update(5.0), 5.0);
        assert_eq!(filter.update(100.0), 5.0); // Spike rejected
    }

    #[test]
    fn test_median3_reset() {
        let mut filter = Median3Filter::new();
        filter.update(100.0);
        filter.update(100.0);
        filter.update(100.0);

        filter.reset();

        // After reset, the first update should return the input directly (count < 3)
        assert_eq!(filter.update(5.0), 5.0);
    }*/
}
