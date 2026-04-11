#![warn(clippy::pedantic)]
#![warn(unused_results)]

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use rand::{RngExt, rng};
use std::hint::black_box;

use filters::{SignalFilter, Pt1Filterf32,Pt2Filterf32,Pt3Filterf32,Pt1FilterVector3df32};
use vector_quaternion_matrix::Vector3df32;

// see target/criterion/Matrix%20Math/report/index.html for results

// # Replace 'v3d_bench' with the name defined in your Cargo.toml [[bench]] section
// RUSTFLAGS="-C target-cpu=native" cargo asm --bench vq_bench "mul_add"

fn bench_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("filter");

    let mut pt1_filter = Pt1Filterf32::new(1.0);
    let mut pt2_filter = Pt2Filterf32::new(1.0);
    let mut pt3_filter = Pt3Filterf32::new(1.0);

    let mut pt1_v3_filter = Pt1FilterVector3df32::new(1.0);

    _ = group.throughput(Throughput::Elements(1));

    _ = group.bench_function("pt1", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| pt1_filter.update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("pt2", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| pt2_filter.update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("pt3", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| pt3_filter.update(black_box(v)),
            BatchSize::SmallInput,
        );

    });

    _ = group.bench_function("pt1_v", |b| {
        b.iter_batched(
            || {
                // Setup: Generate two random vectors
                let a:[f32;3]  = rng().random();
                Vector3df32::from(a)
            },
            |v| pt1_v3_filter.update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, bench_filter);
criterion_main!(benches);
