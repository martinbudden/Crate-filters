#![warn(clippy::pedantic)]
#![warn(unused_results)]

use criterion::{BatchSize, Criterion, Throughput, criterion_group, criterion_main};
use rand::{RngExt, rng};
use std::hint::black_box;

use signal_filters::{
    BiquadFilterf32, MedianFilter3f32, MedianFilter5f32, MovingAverageFilterf32, Pt1FilterVector3df32, Pt1Filterf32,
    Pt2Filterf32, Pt3Filterf32, SignalFilter, SlewRateLimiterf32,
};
use vqm::Vector3df32;

#[allow(clippy::too_many_lines)]
fn bench_filter(c: &mut Criterion) {
    type MovingAverageFilter4f32 = MovingAverageFilterf32<4>;

    let mut group = c.benchmark_group("filter");

    let mut pt1_filter = Pt1Filterf32::new();
    let mut pt2_filter = Pt2Filterf32::new();
    let mut pt3_filter = Pt3Filterf32::new();
    let mut biquad_filter = BiquadFilterf32::new();
    let mut median_filter3 = MedianFilter3f32::new();
    let mut median_filter5 = MedianFilter5f32::new();
    let mut ma_filter4 = MovingAverageFilter4f32::new();
    let mut skew_limiter = SlewRateLimiterf32::with_rates(10.0, 100.0, 0.1);

    let mut pt1_v3_filter = Pt1FilterVector3df32::new();

    _ = group.throughput(Throughput::Elements(1));

    _ = group.bench_function("pt1", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut pt1_filter).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("pt2", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut pt2_filter).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("pt3", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut pt3_filter).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("biquad", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut biquad_filter).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("median3", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut median_filter3).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("median3_c", |b| {
        b.iter_batched(|| 1.0, |v| black_box(&mut median_filter3).update(black_box(v)), BatchSize::SmallInput);
    });

    _ = group.bench_function("median5", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut median_filter5).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("ma4", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut ma_filter4).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("skew_limiter", |b| {
        b.iter_batched(
            || {
                let v: f32 = rng().random();
                v
            },
            |v| black_box(&mut skew_limiter).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    _ = group.bench_function("pt1_v3", |b| {
        b.iter_batched(
            || {
                // Setup: Generate two random vectors
                let a: [f32; 3] = rng().random();
                Vector3df32::from(a)
            },
            |v| black_box(&mut pt1_v3_filter).update(black_box(v)),
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, bench_filter);
criterion_main!(benches);
