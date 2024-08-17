use std::time::Instant;
use ark_bls12_381::Fr;
use ark_ec::bls12::Bls12Config;
use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{AdditiveGroup, MontFp, PrimeField};
use ark_std::UniformRand;
use num_bigint::BigUint;
use num_traits::Zero;
use rand::rngs::OsRng;
use glv_rs::{mul, random_point};

// TODO: compute BETA for any elliptic curve defined over a prime field
// here the method doesn't work with root 1/BETA
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
pub const LAMBDA: Fr = MontFp!("228988810152649578064853576960394133503");


pub fn double_and_add_affine<P: SWCurveConfig>(
    base: &Affine<P>,
    scalar: impl AsRef<[u64]>,
) -> Projective<P> {
    let mut res = Projective::<P>::zero();
    let mut n_adds = 0;
    let mut n_doubles = 0;
    for b in ark_ff::BitIteratorBE::without_leading_zeros(scalar) {
        res.double_in_place();
        n_doubles += 1;
        if b {
            res += base;
            n_adds += 1;
        }
    }
    println!("n_doubles={}, n_adds={}", n_doubles, n_adds);

    res
}

fn main() {
    let k = Fr::rand(&mut OsRng);
    let k_uint = BigUint::from(k.into_bigint());
    println!("bitlen of scalar multiplier {}", k_uint.bits());

    let p = random_point::<ark_bls12_381::g1::Config>();


    let now = Instant::now();
    let _ = mul(p, &k_uint, &BETA, &LAMBDA); // the GLV method computes [k]p
    println!("elapsed GLV = {:?}", now.elapsed());

    let now = Instant::now();
    let _ = double_and_add_affine(&p, k_uint.to_u64_digits().as_slice());
    println!("elapsed double and add = {:?}", now.elapsed());
}