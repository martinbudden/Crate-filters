#![allow(clippy::cast_precision_loss)]

use core::ops::{Add, Mul, Sub};
use num_traits::Zero;
use vector_quaternion_matrix::Vector2df32;

pub type MovingAverageFilterf32<const N: usize> = MovingAverageFilter<f32, N>;
pub type MovingAverageFilterVector3df32<const N: usize> = MovingAverageFilter<Vector2df32, N>;

/// Simple moving average filter.<br>
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
    pub fn reset(&mut self) {
        self.sum = T::zero();
        self.count = 0;
        self.index = 0;
    }

    pub fn update(&mut self, input: T) -> T {
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
    #![allow(unused)]
    #![allow(clippy::float_cmp)]
    use super::*;
    use vector_quaternion_matrix::Vector3df32;
    use vector_quaternion_matrix::Vector3di16;

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {
        is_full::<MovingAverageFilter<f32, 2>>();
    }
    #[test]
    fn moving_average_filter_f32() {
        let mut filter = MovingAverageFilter::<f32, 3>::new();
        assert_eq!(1.0, filter.update(1.0));
        assert_eq!(1.5, filter.update(2.0));
        assert_eq!(2.0, filter.update(3.0));
        assert_eq!(3.0, filter.update(4.0));
        assert_eq!(4.0, filter.update(5.0));
        assert_eq!(5.0, filter.update(6.0));
        assert_eq!(7.0, filter.update(10.0));

        filter.reset();
        assert_eq!(4.0, filter.update(4.0));
        assert_eq!(12.0, filter.update(20.0));
        assert_eq!(5.0, filter.update(-9.0));
    }
    #[test]
    fn moving_average_filter_vector3df32() {
        let mut filter = MovingAverageFilter::<Vector3df32, 4>::new();
        let mut m = filter.update(Vector3df32 { x: 1.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 1.0, y: 0.0, z: -3.0 }, m);

        m = filter.update(Vector3df32 { x: 2.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 1.5, y: 0.0, z: -3.0 }, m);

        m = filter.update(Vector3df32 { x: 3.0, y: 3.0, z: 0.0 });
        assert_eq!(Vector3df32 { x: 2.0, y: 1.0, z: -2.0 }, m);

        m = filter.update(Vector3df32 { x: 4.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 2.5, y: 1.25, z: -2.25 }, m);

        m = filter.update(Vector3df32 { x: 5.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 3.5, y: 1.75, z: -2.25 }, m);

        m = filter.update(Vector3df32 { x: 6.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 4.5, y: 2.25, z: -2.25 }, m);

        m = filter.update(Vector3df32 { x: 10.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 6.25, y: 2.0, z: -3.0 }, m);

        filter.reset();
        m = filter.update(Vector3df32 { x: 4.0, y: 2.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 4.0, y: 2.0, z: -3.0 }, m);

        m = filter.update(Vector3df32 { x: 20.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 12.0, y: 1.0, z: -3.0 }, m);

        m = filter.update(Vector3df32 { x: -9.0, y: 0.0, z: -3.0 });
        assert_eq!(Vector3df32 { x: 5.0, y: 2.0 / 3.0, z: -3.0 }, m);
    }
}
