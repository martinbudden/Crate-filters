pub type MedianFilter3f32 = MedianFilter3<f32>;
pub type MedianFilter3f64 = MedianFilter3<f64>;

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
pub struct MedianFilter3<T> {
    buffer: [T; 3],
    index: usize,
    count: usize,
}

impl<T> Default for MedianFilter3<T>
where
    T: Copy + Default,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> MedianFilter3<T>
where
    T: Copy + Default,
{
    pub fn new() -> Self {
        const COUNT:usize = 3;
        Self { buffer: [T::default(); COUNT], index: 0, count: 0 }
    }
}

impl<T> MedianFilter3<T>
where
    T: Copy + Default + PartialOrd,
{
    pub fn reset(&mut self) {
        const COUNT:usize = 3;
        self.buffer = [T::default(); COUNT];
        self.index = 0;
        self.count = 0;
    }

    pub fn update(&mut self, input: T) -> T {
        const COUNT:usize = 3;
        // Store new sample in the ring buffer
        self.buffer[self.index] = input;
        self.index = (self.index + 1) % COUNT;
        if self.count < COUNT {
            self.count += 1;
        }

        // If buffer isn't full, just return the input
        if self.count < COUNT {
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

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]
    #![allow(unused_results)]
    use super::*;

    #[allow(unused)]
    fn is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<MedianFilter3<f32>>();
    }
    #[test]
    fn median3_spike_rejection() {
        //let mut filter = MovingAverageFilter::<f32, 3>::new();
        //let mut filter: MedianFilter3f32; // as SignalFilter<f32, f32>>;
        let mut filter = MedianFilter3f32::new();

        // 1. Initial values (filling the buffer)
        // Values: [10.0, 0.0, 0.0], count = 1: returns 10.0
        let output = filter.update(10.0);
        assert_eq!(10.0, output);
        // Values: [10.0, 20.0, 0.0], count = 2: returns 20.0
        let output = filter.update(20.0);
        assert_eq!(20.0, output);

        // 2. The Buffer is now full: [10.0, 20.0, 30.0]
        // Median is 20.0
        let output = filter.update(30.0);
        assert_eq!(20.0, output);

        // 3. Test a massive outlier "spike"
        // Buffer: [400.0, 20.0, 30.0] (400 replaces 10)
        // Median is 30.0 (The spike is ignored)
        let output = filter.update(400.0);
        assert_eq!(30.0, output);

        // 4. Return to normal
        // Buffer: [400.0, 25.0, 30.0] (25 replaces 20)
        // Median of {400, 25, 30} is 30.0
        let output = filter.update(25.0);
        assert_eq!(30.0, output);

        // Buffer: [400.0, 25.0, 22.0] (22 replaces 30)
        // Median of {400, 25, 22} is 25.0
        let output = filter.update(22.0);
        assert_eq!(25.0, output);
    }

    #[test]
    fn median3_identical_values() {
        let mut filter = MedianFilter3f32::new();
        filter.update(5.0);
        filter.update(5.0);
        let output = filter.update(5.0);
        assert_eq!(5.0, output);
        let output = filter.update(100.0);
        assert_eq!(5.0, output);
    }

    #[test]
    fn median3_reset() {
        let mut filter = MedianFilter3f32::new();
        filter.update(100.0);
        filter.update(100.0);
        filter.update(100.0);

        filter.reset();

        // After reset, the first update should return the input directly.
        let output = filter.update(5.0);
        assert_eq!(5.0, output);
    }
}
