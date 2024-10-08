use ark_ec::short_weierstrass::{Affine, Projective, SWCurveConfig};
use ark_ff::{AdditiveGroup, PrimeField, Zero};
use ark_std::UniformRand;
use num_bigint::Sign::Minus;
use num_bigint::{BigInt, BigUint, Sign, ToBigInt};
use num_integer::Integer;
use rand::rngs::OsRng;
use std::ops::Mul;

fn phi_inplace<C: SWCurveConfig>(projective: &mut Projective<C>, beta: &C::BaseField) {
    projective.x *= beta;
}

// Computes the remainders and Bézout coefficients for b of steps m,m+1,m+2 for the greatest r_m>sqrt(a)
pub fn truncated_extended_gcd(a: &BigInt, b: &BigInt) -> ([BigInt; 3], [BigInt; 3]) {
    let a_sqrt = a.sqrt();
    let mut remainders: [BigInt; 3] = [a.clone(), b.clone(), BigInt::zero()]; // intermediate remainders
    let mut vs: [BigInt; 3] = [BigInt::zero(), BigInt::from(1u8), BigInt::zero()]; // intermediate Bézout coefficients for b
    let mut q: BigInt; // intermediate quotients

    loop {
        (q, remainders[0]) = remainders[0].div_rem(&remainders[1]);
        q *= &vs[1];
        vs[0] -= &q;
        remainders.swap(0, 1);
        vs.swap(0, 1);
        if remainders[0] <= a_sqrt {
            break;
        }
    }
    (q, remainders[2]) = remainders[0].div_rem(&remainders[1]);
    q *= &vs[1];
    vs[2] = &vs[0] - &q;
    (remainders, vs)
}

fn norm_squared(a: &BigInt, b: &BigInt) -> BigInt {
    a.mul(a) + b.mul(b)
}

// Returns a short vector in the integer lattice spanned by vectors vec1 and vec2
// TODO example 1 of https://arxiv.org/pdf/1310.5250 (but adapt to BLS12-381)
// b1 = (1/2 t_pi-1, c) and (c, 1-1/2 * t_pi)
// where c^2  = q-(t_pi/2)^2
pub fn short_vector(k: &BigUint, r: &[BigInt; 3], v: &[BigInt; 3]) -> (BigInt, BigInt) {
    let k_bigint = k.to_bigint().unwrap();

    let (r0, v0) = (r[0].clone(), -&v[0]);
    let (r1, v1) = (r[1].clone(), -&v[1]);
    let (r2, v2) = (r[2].clone(), -&v[2]);

    let vec1 = (r1.clone(), v1.clone());
    let vec2 = if norm_squared(&r0, &v0) < norm_squared(&r2, &v2) {
        (r0, v0)
    } else {
        (r2, v2)
    };

    let mut det = &vec1.0 * &vec2.1 - &vec1.1 * &vec2.0;
    assert!(!det.is_zero());

    let two_k = &k_bigint << 1;
    let mut b1: BigInt = &two_k * &vec2.1 + &det;
    let mut b2: BigInt = &det - &two_k * &vec1.1;

    det <<= 1;
    b1 = b1.div_floor(&det);
    b2 = b2.div_floor(&det);

    let ux = k_bigint - (&b1 * &vec1.0 + &b2 * &vec2.0);
    let uy = -(&b1 * &vec1.1 + &b2 * &vec2.1);

    (ux, uy)
}

// Precomputations for Shamir's trick
// TODO use WNAF decomposition
pub fn simultaneous_multiple_scalar_multiplication_create_precomputations<C: SWCurveConfig>(
    beta: &C::BaseField,
    window_width: u64,
    p: &Projective<C>,
    u: Sign,
    v: Sign,
) -> Vec<Projective<C>> {
    assert_eq!(64 % window_width, 0);

    let table_size = 1 << (2 * window_width);
    let mut table = vec![Projective::zero(); table_size];

    let mut ps = Vec::with_capacity(1 << window_width);
    let mut qs = Vec::with_capacity(1 << window_width);

    let mut p_step = Projective::zero();
    for _ in 0..(1 << window_width) {
        ps.push(p_step.clone());
        p_step += p;
    }

    let mut transformed_ps = ps.clone();
    for e in &mut transformed_ps {
        phi_inplace(e, beta);
    }
    qs.extend(transformed_ps);

    if u == Minus {
        ps.iter_mut().for_each(|p| {
            let _ = p.neg_in_place();
        });
    }

    if v == Minus {
        qs.iter_mut().for_each(|q| {
            let _ = q.neg_in_place();
        });
    }

    for (i, p) in ps.iter().enumerate() {
        let offset = i * (1 << window_width);
        for (j, q) in qs.iter().enumerate() {
            table[offset + j] += p + q;
        }
    }

    table
}

// get i-th digit of u in base 2^window_width
// assumes window_width divides 64
#[inline(always)]
fn get_digit(u: &[u64], i: u64, window_width: u64) -> usize {
    let idx = (window_width * i) as usize;
    let e = u[idx / 64];
    let shift = idx % 64;
    ((e >> shift) & ((1 << window_width) - 1)) as usize
}

// Shamir's trick (exponentiation using vector-addition chains)
// windowed method
pub fn simultaneous_multiple_scalar_multiplication<C: SWCurveConfig>(
    window_width: u64,
    u: &BigInt,
    v: &BigInt,
    precomputations: Vec<Projective<C>>,
) -> Projective<C> {
    assert_eq!(64 % window_width, 0);
    let mut r = Projective::zero();
    let t: u64 = u.bits().max(v.bits());
    let d = t.div_ceil(window_width);
    let mut u = u.to_u64_digits().1;
    let mut v = v.to_u64_digits().1;
    let m = u.len().max(v.len());
    u.resize(m, 0_u64);
    v.resize(m, 0_u64);
    let mut ui;
    let mut vi;

    for i in (0..d).rev() {
        for _ in 0..window_width {
            r.double_in_place();
        }

        ui = get_digit(&u, i, window_width);
        vi = get_digit(&v, i, window_width);
        if ui != 0 || vi != 0 {
            r += &precomputations[ui * (1 << window_width) + vi];
        }
    }

    r
}

// returns a point in the prime order subgroup G1
// if G1=<g> then random_point=[n]g for random integer n
pub fn random_point<C: SWCurveConfig>() -> Affine<C> {
    Projective::rand(&mut OsRng).into()
}

pub fn mul<C: SWCurveConfig>(
    point: Affine<C>,
    scalar: &BigUint,
    beta: &C::BaseField,
    gcd: &([BigInt; 3], [BigInt; 3]),
) -> Projective<C> {
    let a: BigUint = C::ScalarField::MODULUS.into();

    let scalar = scalar.mod_floor(&a);

    let (u, v) = short_vector(&scalar, &gcd.0, &gcd.1);

    let window_width: u64 = 2;
    let precomputations = simultaneous_multiple_scalar_multiplication_create_precomputations(
        beta,
        window_width,
        &point.into(),
        u.sign(),
        v.sign(),
    );
    simultaneous_multiple_scalar_multiplication(window_width, &u, &v, precomputations)
}

#[cfg(test)]
pub mod tests {
    use std::ops::{Mul, Rem};

    use ark_bls12_381::Fr;
    use ark_ec::bls12::Bls12Config;
    use ark_ec::CurveGroup;
    use ark_ff::{MontConfig, MontFp};
    use num_bigint::BigUint;
    use num_traits::Signed;
    use proptest::proptest;

    use super::*;

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

    fn phi_inplace<C: SWCurveConfig>(affine: &mut Affine<C>, beta: &C::BaseField) {
        affine.x *= beta;
    }

    #[test]
    fn test_phi() {
        let mut p = random_point::<ark_bls12_381::g1::Config>();
        assert!(p.is_on_curve());
        assert!(p.is_in_correct_subgroup_assuming_on_curve());
        // at this point order(p)=|Fr|=r
        let mul = Affine::mul(p.clone(), LAMBDA);
        phi_inplace(&mut p, &BETA);
        assert_eq!(p, mul.into_affine()); // [lambda]p=phi(p)
    }

    use proptest::prelude::*;

    // Custom strategy to generate random `Fr` elements
    fn arb_fr() -> impl Strategy<Value = Fr> {
        any::<[u8; 32]>().prop_map(|bytes| Fr::from_le_bytes_mod_order(&bytes))
    }

    fn arb_bigint() -> impl Strategy<Value = BigInt> {
        any::<[u8; 32]>().prop_map(|bytes| {
            BigInt::from_bytes_le(Sign::Plus, bytes.as_slice())
                .rem(BigInt::from(ark_bls12_381::FrConfig::MODULUS))
        })
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(10_000))]
        #[test]
        fn test_simultaneous_multiple_scalar_multiplication(k in arb_fr()) {
            let n: BigUint = ark_bls12_381::FrConfig::MODULUS.into();
            let lambda: BigUint = LAMBDA.into();
            let k_uint = BigUint::from(k.into_bigint());

            let p = random_point::<ark_bls12_381::g1::Config>();
            let n_bigint = n.to_bigint().unwrap();
            let lambda_bigint = lambda.to_bigint().unwrap();
            let gcd = truncated_extended_gcd(&n_bigint, &lambda_bigint);

            assert_eq!(mul(p, &k_uint, &BETA, &gcd), p.mul(k)); // the GLV method computes [k]p
        }

        #[test]
        fn test_short_vector(k in arb_bigint()) {
            let a: BigUint = ark_bls12_381::FrConfig::MODULUS.into();
            let b: BigUint = LAMBDA.into();

            let n_bigint = a.to_bigint().unwrap();
            let lambda_bigint = b.to_bigint().unwrap();
            let gcd = truncated_extended_gcd(&n_bigint, &lambda_bigint);

            let (k1, k2) = short_vector(&k.clone().to_biguint().unwrap(), &gcd.0, &gcd.1);
            assert!(k1.abs() < BigInt::from(a.clone()));
            assert!(k2.abs() < BigInt::from(a.clone()));
            let k_ = (&k1 + Into::<BigInt>::into(b) * &k2).mod_floor(&a.to_bigint().unwrap());

            assert_eq!(k_.to_biguint().unwrap(), k.to_biguint().unwrap()); // k = k1 + lambda * k2 mod a
            assert!(BigInt::max(k1, k2) < a.sqrt().into()); // (k1, k2) is a short vector
        }

        #[test]
        fn test_truncated_extended_gcd(a in arb_bigint(), b in arb_bigint()) {
            let (r, v) = truncated_extended_gcd(&a, &b);
            let mut res: BigInt = r[1].clone();
            res *= v[2].abs();
            let mut res2 = r[2].clone();
            res2 *= v[1].abs();
            res += res2;

            assert_eq!(res, a);

            let mut res: BigInt = r[0].clone();
            res *= v[1].abs();
            let mut res2 = r[1].clone();
            res2 *= v[0].abs();
            res += res2;

            assert_eq!(res, a);
        }
    }
}
