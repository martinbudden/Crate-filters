pub type MedianFilter5f32 = MedianFilter5<f32>;
pub type MedianFilter5f64 = MedianFilter5<f64>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MedianFilter5<T> {
    buffer: [T; 5],
    sorted_buffer: [T; 5],
    index: usize,
}

impl<T> Default for MedianFilter5<T>
where
    T: Copy + Default,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> MedianFilter5<T>
where
    T: Copy + Default,
{
    pub fn new() -> Self {
        const COUNT:usize = 5;
        Self { buffer: [T::default(); COUNT], sorted_buffer: [T::default(); COUNT], index: 0 }
    }
}

impl<T> MedianFilter5<T>
where
    T: Copy + Default + PartialOrd,
{
    pub fn reset(&mut self) {
        const COUNT:usize = 5;
        self.buffer = [T::default(); COUNT];
        self.sorted_buffer = [T::default(); COUNT];
        self.index = 0;
    }

    pub fn update(&mut self, input: T) -> T {
        use core::cmp::Ordering;
        const COUNT:usize = 5;
        const MID:usize = 2;

        let oldest = self.buffer[self.index];
        self.buffer[self.index] = input;
        self.index = (self.index + 1) % COUNT;

        // 1. Find indices without manual loops or range bounds
        // .position() returns the first index where the predicate is true
        let old_pos = self.sorted_buffer.iter().position(|&val| val == oldest).unwrap_or(0); // Should always find a match

        let new_pos = self.sorted_buffer.iter().position(|&val| input < val).unwrap_or(COUNT);

        // 2. Perform the update

        match old_pos.cmp(&new_pos) {
            Ordering::Less => {
                // Shift left
                self.sorted_buffer.copy_within(old_pos + 1..new_pos, old_pos);
                self.sorted_buffer[new_pos - 1] = input;
            }
            Ordering::Greater => {
                // Shift right
                self.sorted_buffer.copy_within(new_pos..old_pos, new_pos + 1);
                self.sorted_buffer[new_pos] = input;
            }
            Ordering::Equal => {
                self.sorted_buffer[old_pos] = input;
            }
        }
        self.sorted_buffer[MID]
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]
    #![allow(unused_results)]
    use super::*;

    #[test]
    fn median5() {
        // Initialize with 0.0 (assuming T = f32)
        let mut filter = MedianFilter5::<f32>::new();

        // Fill the buffer: [0, 0, 0, 0, 10] -> Sorted: [0, 0, 0, 0, 10]
        // Median is index 2 (the middle 0)
        let output = filter.update(10.0);
        assert_eq!(0.0, output);

        // [0, 0, 0, 10, 20] -> Sorted: [0, 0, 0, 10, 20]
        let output = filter.update(20.0);
        assert_eq!(0.0, output);

        // [0, 0, 10, 20, 30] -> Sorted: [0, 0, 10, 20, 30]
        // Now the middle element is 10
        let output = filter.update(30.0);
        assert_eq!(10.0, output);

        // [0, 10, 20, 30, 5] -> Sorted: [0, 5, 10, 20, 30]
        let output = filter.update(5.0);
        assert_eq!(10.0, output);

        // [10, 20, 30, 5, 15] -> Sorted: [5, 10, 15, 20, 30]
        // Median is now 15
        let output = filter.update(15.0);
        assert_eq!(15.0, output);
    }
    #[test]
    fn median5_reset() {
        let mut filter = MedianFilter5::<f32>::new();
        filter.update(100.0);
        filter.update(200.0);

        filter.reset();
        // assert_eq!(filter.update(5.0), 0.0);
    }
}
