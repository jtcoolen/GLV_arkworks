use ark_bls12_381::Fr;
use ark_ec::bls12::Bls12Config;
use ark_ec::scalar_mul::glv::GLVConfig;
use ark_ec::short_weierstrass::Projective;
use ark_ff::{MontFp, PrimeField};
use ark_std::UniformRand;
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use glv_rs::{mul, random_point};
use num_bigint::BigUint;
use rand::rngs::OsRng;

pub const BETA: <ark_bls12_381::Config as Bls12Config>::Fp = MontFp!("4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436");
/* Generated with PARI/GP: BETA is an element of order 3 in Fq
? q
4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
? isprime(q)
1
? 1/znprimroot(q)^((q-1)/3)
Mod(4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436, 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787) */

/* LAMBDA is a primitive 3rd root of unity mod r
? r
52435875175126190479447740508185965837690552500527637822603658699938581184513
? znprimroot(r)^((r-1)/3)
Mod(228988810152649578064853576960394133503, 52435875175126190479447740508185965837690552500527637822603658699938581184513)*/
pub const LAMBDA: ark_bls12_381::Fr = MontFp!("228988810152649578064853576960394133503");

pub fn bench_glv(c: &mut Criterion) {
    c.bench_function("BLS12-381 G1 mul (glv-rs crate)", |b| {
        b.iter(|| {
            mul(
                black_box(random_point::<ark_bls12_381::g1::Config>()),
                black_box(&BigUint::from(Fr::rand(&mut OsRng).into_bigint())),
                black_box(&BETA),
                black_box(&LAMBDA),
            )
        })
    });
}

pub fn bench_glv_projective_arkworks(c: &mut Criterion) {
    c.bench_function("BLS12-381 G1 projective mul (arkworks crate)", |b| {
        b.iter(|| {
            GLVConfig::glv_mul_projective(
                black_box(Projective::from(random_point::<ark_bls12_381::g1::Config>())),
                black_box(Fr::rand(&mut OsRng)),
            )
        })
    });
}

pub fn bench_glv_affine_arkworks(c: &mut Criterion) {
    c.bench_function("BLS12-381 G1 affine mul (arkworks crate)", |b| {
        b.iter(|| {
            GLVConfig::glv_mul_affine(
                black_box(random_point::<ark_bls12_381::g1::Config>()),
                black_box(Fr::rand(&mut OsRng)),
            )
        })
    });
}

criterion_group!(
    benches,
    bench_glv,
    bench_glv_affine_arkworks,
    bench_glv_projective_arkworks
);
criterion_main!(benches);
