use vector_quaternion_matrix::Vector3df32;

/// The Filter Definition trait
///
/// Note the "filter" function is called "apply"
/// This is because filter.apply() reads better than filter.filter()
/// and it also avoids confusion with the filter function in the Iterator trait.
/// ```
/// use filters::{Pt1Filterf32,FilterSignal};
///
/// let mut filter = Pt1Filterf32::<f32>::new(1.0);
///
/// let mut reading:f32 = 2.7;
/// reading = filter.apply(reading);
///
/// assert_eq!(2.7, reading);
/// ```
pub trait FilterSignal<T, F> {
    fn apply(&mut self, input: T) -> T;
    fn reset(&mut self);
}

/// The Filter Extension Trait (for the data)
///
/// T is the type being filtered, so it might be an f32 or a Vector3df32
/// R is the type of the filter's internal constant. It is either f32 or f64
/// and should correspond to T.
///
/// So:
///     Ff T is f32, R is f32.
///     If T is f64, R is f64.
///     If T is Vector3df32, R is f32
///     If T is Vector3df64, R is f64
pub trait ApplyFilter<T, R> {
    fn apply_using<F: FilterSignal<T, R>>(&mut self, filter: &mut F) -> &mut Self;
}

impl ApplyFilter<f32, f32> for f32 {
    fn apply_using<F: FilterSignal<f32, f32>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is f32, filter.apply takes and returns f32
        *self = filter.apply(*self);
        self
    }
}

impl ApplyFilter<f64, f64> for f64 {
    fn apply_using<F: FilterSignal<f64, f64>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is f64, filter.apply takes and returns f64
        *self = filter.apply(*self);
        self
    }
}

impl ApplyFilter<Vector3df32, f32> for Vector3df32 {
    fn apply_using<F: FilterSignal<Vector3df32, f32>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is Vector3df32, filter.apply handles the whole vector at once
        *self = filter.apply(*self);
        self
    }
}

#[cfg(any(debug_assertions, test))]
mod tests {
    #[allow(unused)]
    use super::*;

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn _is_full<T: Sized + Send + Sync + Unpin + Copy + Clone + Default + PartialEq>() {}

    #[test]
    fn normal_types() {}
    #[test]
    fn default() {}
}
