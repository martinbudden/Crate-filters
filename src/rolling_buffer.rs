#![allow(unused)]

use num_traits::Zero;

#[derive(Debug)]
pub struct RollingBuffer<T, const N: usize> {
    /// The virtual beginning of the rolling buffer.
    begin: usize,
    /// The virtual end of the rolling buffer (one behind the last element).
    end: usize,
    /// The number of items in the rolling buffer.
    size: usize,
    // need one spare empty cell so we can avoid end == begin when full
    //buffer: [core::mem::MaybeUninit<T>; N], // Use MaybeUninit for safe uninitialized data
    buffer: [T; N],
}

impl<T, const N: usize> Default for RollingBuffer<T, N>
where
    T: Copy + Zero,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T, const N: usize> RollingBuffer<T, N>
where
    T: Copy + Zero,
{
    pub fn new() -> Self {
        // SAFETY: Creating an uninitialized array is safe since we use MaybeUninit
        // let buffer = unsafe { core::mem::MaybeUninit::uninit().assume_init() };
        Self {
            begin: 0,
            end: 0,
            size: 0,
            //buffer: unsafe { core::mem::MaybeUninit::uninit().assume_init() }
            buffer: [T::zero();N],
        }
    }
    // need one spare empty cell so we can avoid end == begin when full
    pub fn capacity(&self) -> usize {
        N - 1
    }
    pub fn size(&self) -> usize {
        self.size
    }
    pub fn is_empty(&self) -> bool {
        self.size == 0
    }
    pub fn is_full(&self) -> bool {
        self.size == self.capacity()
    }
    pub fn front(&self) -> Option<T> {
        if self.size == 0 {
            return None;
        }
        Some(self.buffer[self.begin])
    }
    pub fn back(&self) -> Option<T> {
        if self.size == 0 {
            return None;
        }
        if self.end > 0 {
            Some(self.buffer[self.end - 1])
        } else {
            Some(self.buffer[self.capacity()])
        }
    }
    pub fn begin(&self) -> usize {
        self.begin
    }
    pub fn end(&self) -> usize {
        self.end
    }

    pub fn at(&self, index: usize) -> Option<T> {
        if index >= self.size {
            return None;
        }
        let mut pos = self.begin + index;
        if pos > self.capacity() {
            pos -= self.capacity() + 1;
        }
        Some(self.buffer[pos])
    }

    pub fn push_back(&mut self, value: T) {
        self.buffer[self.end] = value; // buffer size is capacity() + 1, so always OK to store value at self.end
        self.end += 1;

        if self.size >= self.capacity() {
            // buffer is full, so don't increment size, instead drop items off front by incrementing self.begin
            self.begin += 1;
            // wrap self.begin if required
            if self.begin > self.capacity() {
                self.begin = 0;
            }
            // wrap self.end if required
            if self.end > self.capacity() {
                self.end = 0;
            }
        } else {
            self.size += 1;
        }
    }
}

#[cfg(any(debug_assertions, test))]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let rb = RollingBuffer::<f32, 3>::new();
        assert_eq!(true, rb.is_empty());
        assert_eq!(false, rb.is_full());
        assert_eq!(2, rb.capacity());
        assert_eq!(0, rb.size());
        assert_eq!(0, rb.begin());
        assert_eq!(0, rb.end());
        assert_eq!(None, rb.at(0));
    }

    #[test]
    fn rolling_buffer_size() {
        let mut rb = RollingBuffer::<i32, 5>::new();
        assert_eq!(true, rb.is_empty());
        assert_eq!(4, rb.capacity());

        assert_eq!(0, rb.size());
        rb.push_back(10);
        assert_eq!(1, rb.size());

        rb.push_back(11);
        assert_eq!(2, rb.size());

        rb.push_back(12);
        assert_eq!(3, rb.size());

        rb.push_back(13);
        assert_eq!(4, rb.size());

        // the buffer is full, so size will no longer increase
        rb.push_back(14);
        assert_eq!(4, rb.size());

        rb.push_back(15);
        assert_eq!(4, rb.size());
        assert_eq!(4, rb.capacity());
    }

    #[test]
    fn rolling_buffer_front_back() {
        let mut rb = RollingBuffer::<i32, 5>::new();

        rb.push_back(10);
        assert_eq!(Some(10), rb.front());
        assert_eq!(Some(10), rb.back());
        assert_eq!(Some(10), rb.at(0));
        assert_eq!(None, rb.at(1));
        assert_eq!(None, rb.at(2));
        assert_eq!(None, rb.at(3));
        assert_eq!(None, rb.at(4));
        assert_eq!(None, rb.at(5));
        assert_eq!(None, rb.at(6));

        rb.push_back(11);
        assert_eq!(Some(10), rb.front());
        assert_eq!(Some(11), rb.back());
        assert_eq!(Some(10), rb.at(0));
        assert_eq!(Some(11), rb.at(1));
        assert_eq!(None, rb.at(2));
        assert_eq!(None, rb.at(3));
        assert_eq!(None, rb.at(4));
        assert_eq!(None, rb.at(5));
        assert_eq!(None, rb.at(6));

        rb.push_back(12);
        assert_eq!(Some(10), rb.front());
        assert_eq!(Some(12), rb.back());
        assert_eq!(Some(10), rb.at(0));
        assert_eq!(Some(11), rb.at(1));
        assert_eq!(Some(12), rb.at(2));
        assert_eq!(None, rb.at(3));
        assert_eq!(None, rb.at(4));
        assert_eq!(None, rb.at(5));
        assert_eq!(None, rb.at(6));

        rb.push_back(13);
        assert_eq!(Some(10), rb.front());
        assert_eq!(Some(13), rb.back());
        assert_eq!(Some(10), rb.at(0));
        assert_eq!(Some(11), rb.at(1));
        assert_eq!(Some(12), rb.at(2));
        assert_eq!(Some(13), rb.at(3));
        assert_eq!(None, rb.at(4));
        assert_eq!(None, rb.at(5));
        assert_eq!(None, rb.at(6));

        // now items start dropping off the front
        rb.push_back(14);
        assert_eq!(Some(11), rb.front());
        assert_eq!(Some(14), rb.back());
        assert_eq!(Some(11), rb.at(0));
        assert_eq!(Some(12), rb.at(1));
        assert_eq!(Some(13), rb.at(2));
        assert_eq!(Some(14), rb.at(3));
        assert_eq!(None, rb.at(4));
        assert_eq!(None, rb.at(5));
        assert_eq!(None, rb.at(6));

        rb.push_back(15);
        assert_eq!(Some(12), rb.front());
        assert_eq!(Some(15), rb.back());
        assert_eq!(Some(12), rb.at(0));
        assert_eq!(Some(13), rb.at(1));
        assert_eq!(Some(14), rb.at(2));
        assert_eq!(Some(15), rb.at(3));
        assert_eq!(None, rb.at(4));
        assert_eq!(None, rb.at(5));
        assert_eq!(None, rb.at(6));
    }
}
