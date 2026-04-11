/// Note the "filter" function is called "update" rather than "apply" or "filter".
/// This is because update implies the filter has internal state whereas
/// apply sometimes implies a "pure" mathematical function with no memory.
///
/// `filter.update()` reads better than `filter.filter()`
/// and it also avoids confusion with the filter function in the Iterator trait.
use vector_quaternion_matrix::{Vector2df32, Vector2df64, Vector3df32, Vector3df64};

/// Filter Definition trait.
/// `filter.reset()`
/// `filter.update(value)`
/// ```
///
/// use filters::{Pt2Filterf32,SignalFilter};
///
/// let mut filter = Pt2Filterf32::new(0.25);
/// let mut value:f32 = 1.0;
///
/// value = filter.update(value);
///
/// assert_eq!(0.0625, value);
/// ```
pub trait SignalFilter<T, F> {
    fn reset(&mut self);
    fn update(&mut self, input: T) -> T;
}

// `T` is the type being filtered, so it might be an `f32` or a `Vector3df32`
// `R` is the type of the filter's internal constant. It is either `f32` or `f64`
// and should correspond to `T`.
//
// So:
// * If `T` is `f32`, `R` is `f32`.
//
//  * If `T` is `f64`, `R` is `f64`.
// * If `T` is `Vector3df32`, `R` is `f32`
// * If `T` is `Vector3df64`, `R` is `f64`

/// Adds `value.update_using(&mut filter)` method call syntax to `SignalFilter`
/// ```
/// use filters::{Pt2Filterf32,UpdateFilter};
/// let mut filter = Pt2Filterf32::new(0.25);
/// let mut value: f32 = 1.0;
///
/// value.update_using(&mut filter);
///
/// assert_eq!(0.0625, value);
/// ```
pub trait UpdateFilter<T, R> {
    fn update_using<F: SignalFilter<T, R>>(&mut self, filter: &mut F) -> &mut Self;
}

impl UpdateFilter<f32, f32> for f32 {
    fn update_using<F: SignalFilter<f32, f32>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is f32, filter.update takes and returns f32
        *self = filter.update(*self);
        self
    }
}

impl UpdateFilter<f64, f64> for f64 {
    fn update_using<F: SignalFilter<f64, f64>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is f64, filter.update takes and returns f64
        *self = filter.update(*self);
        self
    }
}

impl UpdateFilter<Vector2df32, f32> for Vector2df32 {
    fn update_using<F: SignalFilter<Vector2df32, f32>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is Vector2df32, filter.update handles the whole vector at once
        *self = filter.update(*self);
        self
    }
}

impl UpdateFilter<Vector2df64, f64> for Vector2df64 {
    fn update_using<F: SignalFilter<Vector2df64, f64>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is Vector2df64, filter.update handles the whole vector at once
        *self = filter.update(*self);
        self
    }
}
impl UpdateFilter<Vector3df32, f32> for Vector3df32 {
    fn update_using<F: SignalFilter<Vector3df32, f32>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is Vector3df32, filter.update handles the whole vector at once
        *self = filter.update(*self);
        self
    }
}

impl UpdateFilter<Vector3df64, f64> for Vector3df64 {
    fn update_using<F: SignalFilter<Vector3df64, f64>>(&mut self, filter: &mut F) -> &mut Self {
        // *self is Vector3df64, filter.update handles the whole vector at once
        *self = filter.update(*self);
        self
    }
}
#[cfg(any(debug_assertions, test))]
mod tests {
    #[allow(unused)]
    use super::*;

    fn _is_normal<T: Sized + Send + Sync + Unpin>() {}
    fn _is_full<T: Sized + Send + Sync + Unpin + Clone + Copy + Default + PartialEq>() {}

    #[test]
    fn normal_types() {}
    #[test]
    fn default() {}
}
