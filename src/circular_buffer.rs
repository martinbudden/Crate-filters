#![allow(unused)]

pub struct CircularBuffer<T, const N: usize> {
    /// The virtual beginning of the circular buffer.
    begin: usize,
    /// The virtual end of the circular buffer (one behind the last element).
    end: usize,
    /// The number of items in the circular buffer.
    size: usize,
    // need one spare empty cell so we can avoid end == begin when full
    buffer: [T; N],
}

impl<T, const N: usize> CircularBuffer<T, N>
where
    T: Copy,
{
    pub fn new() -> Self {
        // SAFETY: Creating an uninitialized array is safe since we use MaybeUninit
        let buffer = unsafe { core::mem::MaybeUninit::uninit().assume_init() };
        Self {
            begin: 0,
            end: 0,
            size: 0,
            buffer,
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

    pub fn push_back(&mut self, value: T) -> bool {
        if self.is_full() {
            return false;
        }
        self.size += 1;
        self.buffer[self.end] = value; // sizeof(_buffer) = CAPACITY + 1, so always OK to store value at _end
        self.end += 1;
        // wrap _end if required
        if self.end > self.capacity() {
            self.end = 0;
        }
        true
    }
    pub fn pop_front(&mut self) -> Option<T> {
        if self.is_empty() {
            return None;
        }
        self.size -= 1;
        let value = self.buffer[self.begin];
        self.begin += 1;
        // wrap _begin if required
        if self.begin > self.capacity() {
            self.begin = 0;
        }
        Some(value)
    }
}

#[cfg(any(debug_assertions, test))]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let cb = CircularBuffer::<f32, 3>::new();
        assert_eq!(true, cb.is_empty());
        assert_eq!(false, cb.is_full());
        assert_eq!(2, cb.capacity());
        assert_eq!(0, cb.size());
        assert_eq!(0, cb.begin());
        assert_eq!(0, cb.end());
        assert_eq!(None, cb.at(0));
    }

    #[test]
    fn circular_buffer_size() {
        let mut cb = CircularBuffer::<i32, 5>::new();
        assert_eq!(true, cb.is_empty());
        assert_eq!(4, cb.capacity());
        assert_eq!(None, cb.pop_front());

        assert_eq!(0, cb.size());
        cb.push_back(10);
        assert_eq!(1, cb.size());

        cb.push_back(11);
        assert_eq!(2, cb.size());

        cb.push_back(12);
        assert_eq!(3, cb.size());

        cb.push_back(13);
        assert_eq!(4, cb.size());

        // the buffer is full, so size will no longer increase
        cb.push_back(14);
        assert_eq!(4, cb.size());

        cb.push_back(15);
        assert_eq!(4, cb.size());
        assert_eq!(4, cb.capacity());
    }

    #[test]
    fn circular_buffer_front_back() {
        let mut cb = CircularBuffer::<i32, 5>::new();
        let mut success: bool;

        success = cb.push_back(10);
        assert_eq!(true, success);
        assert_eq!(Some(10), cb.front());
        assert_eq!(Some(10), cb.back());
        assert_eq!(Some(10), cb.at(0));
        assert_eq!(None, cb.at(1));
        assert_eq!(None, cb.at(2));
        assert_eq!(None, cb.at(3));
        assert_eq!(None, cb.at(4));
        assert_eq!(None, cb.at(5));
        assert_eq!(None, cb.at(6));

        success = cb.push_back(11);
        assert_eq!(true, success);
        assert_eq!(Some(10), cb.front());
        assert_eq!(Some(11), cb.back());
        assert_eq!(Some(10), cb.at(0));
        assert_eq!(Some(11), cb.at(1));
        assert_eq!(None, cb.at(2));
        assert_eq!(None, cb.at(3));
        assert_eq!(None, cb.at(4));
        assert_eq!(None, cb.at(5));
        assert_eq!(None, cb.at(6));

        success = cb.push_back(12);
        assert_eq!(true, success);
        assert_eq!(Some(10), cb.front());
        assert_eq!(Some(12), cb.back());
        assert_eq!(Some(10), cb.at(0));
        assert_eq!(Some(11), cb.at(1));
        assert_eq!(Some(12), cb.at(2));
        assert_eq!(None, cb.at(3));
        assert_eq!(None, cb.at(4));
        assert_eq!(None, cb.at(5));
        assert_eq!(None, cb.at(6));

        success = cb.push_back(13);
        assert_eq!(true, success);
        assert_eq!(Some(10), cb.front());
        assert_eq!(Some(13), cb.back());
        assert_eq!(Some(10), cb.at(0));
        assert_eq!(Some(11), cb.at(1));
        assert_eq!(Some(12), cb.at(2));
        assert_eq!(Some(13), cb.at(3));
        assert_eq!(None, cb.at(4));
        assert_eq!(None, cb.at(5));
        assert_eq!(None, cb.at(6));

        // now buffer is full, so pushing items will fail
        success = cb.push_back(14);
        assert_eq!(false, success);
        assert_eq!(Some(10), cb.front());
        assert_eq!(Some(13), cb.back());
        assert_eq!(Some(10), cb.at(0));
        assert_eq!(Some(11), cb.at(1));
        assert_eq!(Some(12), cb.at(2));
        assert_eq!(Some(13), cb.at(3));
        assert_eq!(None, cb.at(4));
        assert_eq!(None, cb.at(5));
        assert_eq!(None, cb.at(6));

        assert_eq!(Some(10), cb.pop_front());
        cb.push_back(15);
        assert_eq!(Some(11), cb.front());
        assert_eq!(Some(15), cb.back());
        assert_eq!(Some(11), cb.at(0));
        assert_eq!(Some(12), cb.at(1));
        assert_eq!(Some(13), cb.at(2));
        assert_eq!(Some(15), cb.at(3));
        assert_eq!(None, cb.at(4));
        assert_eq!(None, cb.at(5));
        assert_eq!(None, cb.at(6));
    }
}
