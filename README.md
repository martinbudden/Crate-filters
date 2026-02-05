# filters Rust Crate ![license](https://img.shields.io/badge/license-MIT-green) [![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) ![open source](https://badgen.net/badge/open/source/blue?icon=github)


This crate contains a collection of filter classes. They have been developed for use in stabilized vehicles (self balancing robots and aircraft)and have been used to:

1. Filter gyro and accelerometer output for use in the Attitude and Heading Reference System (AHRS).
2. Filter motor encoder values for use in the motor controller.
3. Filter derivative terms in a PID controller.
4. Filter motor power input values to smooth the motor speed.

## Original implementation

This crate was originally implemented as a c++ library. The 
[original implementation can be found here](https://github.com/martinbudden/Library-Filter).

## License

Licensed under either of

* Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
