use ark_bls12_381::{g1, Fq, Fr};
use ark_ff::{field_new, Zero};
use ark_std::UniformRand;
use num_bigint::{BigInt, BigUint, Sign, ToBigInt};
use num_integer::Integer;
use num_traits::Signed;

// TODO: compute BETA for any elliptic curve defined over a prime field
// here the method doesn't work with root 1/BETA
const BETA: Fq = field_new!(Fq, "4002409555221667392624310435006688643935503118305586438271171395842971157480381377015405980053539358417135540939436");
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
const LAMBDA: Fr = field_new!(Fr, "228988810152649578064853576960394133503");

fn phi_inplace(affine: &mut g1::G1Affine) {
    affine.x *= &BETA;
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
    a.pow(2) + b.pow(2)
}

// Returns a short vector in the integer lattice spanned by vectors vec1 and vec2
// TODO: change inputs to scalar field elements
pub fn short_vector(n: &BigUint, lambda: &BigUint, k: &BigUint) -> (BigInt, BigInt) {
    let (r, v) = truncated_extended_gcd(&n.to_bigint().unwrap(), &lambda.to_bigint().unwrap());
    let (r0, v0) = (r[0].clone(), -&v[0]);
    let (r1, v1) = (r[1].clone(), -&v[1]);
    let (r2, v2) = (r[2].clone(), -&v[2]);
    let vec1 = (r1, v1);
    let vec2: (BigInt, BigInt) = if norm_squared(&r0, &v0) < norm_squared(&r2, &v2) {
        (r0, v0)
    } else {
        (r2, v2)
    };
    // TODO: check that vec1 and vec2 are linearly independent
    let mut det = &vec1.0 * &vec2.1 - &vec1.1 * &vec2.0;
    let mut b1: BigInt = 2 * k.to_bigint().unwrap() * &vec2.1 + &det;
    let mut b2: BigInt = -2 * k.to_bigint().unwrap() * &vec1.1 + &det;
    det *= 2;
    b1 = b1.div_floor(&det);
    b2 = b2.div_floor(&det);
    let ux = k.to_bigint().unwrap() - (&b1 * &vec1.0 + &b2 * &vec2.0);
    let uy = -(&b1 * &vec1.1 + &b2 * &vec2.1);
    (ux, uy)
}

// Precomputations for Shamir's trick
pub fn simultaneous_multiple_scalar_multiplication_create_precomputations(
    window_width: usize,
    p: &g1::G1Projective,
    q: &g1::G1Projective,
    u: Sign,
    v: Sign,
) -> Vec<Vec<g1::G1Projective>> {
    assert!(64 % window_width == 0);
    // TODO: use flat vector
    let mut table = vec![vec![g1::G1Projective::zero(); 1 << window_width]; 1 << window_width];
    let mut pp = g1::G1Projective::zero();
    let mut qq = g1::G1Projective::zero();
    let (p, q) = match (u, v) {
        (Sign::Minus, Sign::Minus) => (-*p, -*q),
        (Sign::Plus, Sign::Minus) => (*p, -*q),
        (Sign::Minus, Sign::Plus) => (-*p, *q),
        _ => (*p, *q),
    };
    for i in 0..=((1 << window_width) - 1) {
        qq.set_zero();
        for j in 0..=((1 << window_width) - 1) {
            table[i][j] = pp + qq;
            qq += q;
        }
        pp += p;
    }
    table
}

// get i-th digit of u in base 2^window_width
// assumes window_width divides 64
#[inline(always)]
fn get_digit(u: &[u64], i: usize, window_width: usize) -> usize {
    let idx = window_width * i;
    let e = u[idx / 64] as usize;
    let idx_e = idx % 64;
    let mask: usize = (1 << window_width) - 1;
    ((mask << idx_e) & e) >> idx_e
}

fn pad_vec(v: &mut Vec<u64>, len: usize) {
    assert!(v.len() <= len);
    while v.len() < len {
        v.push(0_u64);
    }
}

// Shamir's trick (exponentiation using vector-addition chains)
// windowed method
pub fn simultaneous_multiple_scalar_multiplication(
    window_width: usize,
    u: &BigInt,
    v: &BigInt,
    precomputations: Vec<Vec<g1::G1Projective>>,
) -> g1::G1Projective {
    assert!(64 % window_width == 0);
    let mut r = g1::G1Projective::zero();
    // TODO: handle unwrap
    let t: usize = std::cmp::max(u.abs().bits(), v.abs().bits())
        .try_into()
        .unwrap();
    let d = num_integer::Integer::div_ceil(&t, &window_width);
    let mut u = u.abs().to_u64_digits().1;
    let mut v = v.abs().to_u64_digits().1;
    let dd = std::cmp::max(u.len(), v.len());
    pad_vec(&mut u, dd as usize);
    pad_vec(&mut v, dd as usize);
    let mut ui;
    let mut vi;
    for i in (0..=(d - 1)).rev() {
        r *= (1 << window_width).into();
        ui = get_digit(&u, i, window_width);
        vi = get_digit(&v, i, window_width);
        if !(ui == 0 && vi == 0) {
            r += precomputations[ui][vi];
        }
    }
    r
}

// returns a point in the prime order subgroup G1
// if G1=<g> then random_point=[n]g for random integer n
pub fn random_point() -> g1::G1Affine {
    let mut rng = ark_std::test_rng();
    g1::G1Projective::rand(&mut rng).into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::AffineCurve;
    use ark_ff::{BigInteger256, FpParameters};
    use num_bigint::BigUint;
    use num_traits::Signed;
    use std::ops::Rem;

    #[test]
    fn test_truncated_extended_gcd() {
        // TODO: make test with randomized values
        let a: BigUint = ark_bls12_381::FrParameters::MODULUS.into();
        let a: BigInt = a.into();
        let b: BigUint = LAMBDA.into();
        let b: BigInt = b.into();

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

    #[test]
    fn test_short_vector() {
        let a: BigUint = ark_bls12_381::FrParameters::MODULUS.into();
        let b: BigUint = LAMBDA.into();

        let k: BigInteger256 =
            field_new!(Fr, "654589647885213749678705807402976612273608534928").into();

        let (k1, k2) = short_vector(&a, &b, &k.into());
        let kk: BigInt = BigInt::rem(&k1 + Into::<BigInt>::into(b) * &k2, &a.to_bigint().unwrap());
        assert_eq!(kk.to_biguint().unwrap(), k.into()); // k = k1 + lambda * k2 mod a
        assert!(BigInt::max(k1, k2) < a.sqrt().into()); // (k1, k2) is a short vector
    }

    #[test]
    fn test_phi() {
        let mut p = random_point();
        assert!(p.is_on_curve());
        assert!(p.is_in_correct_subgroup_assuming_on_curve());
        // at this point order(p)=|Fr|=r
        let mul = ark_ec::AffineCurve::mul(&p, LAMBDA);
        phi_inplace(&mut p);
        assert_eq!(p.into_projective(), mul); // [lambda]p=phi(p)
    }

    #[test]
    fn test_simultaneous_multiple_scalar_multiplication() {
        let a: BigUint = ark_bls12_381::FrParameters::MODULUS.into();
        let b: BigUint = LAMBDA.into();

        let k: BigInteger256 =
            field_new!(Fr, "654589647885213749678705807976612273608534928").into();

        let p = random_point();
        let mut phi_p = p.clone();
        phi_inplace(&mut phi_p);
        let (u, v) = short_vector(&a, &b, &k.into());

        let window_width: usize = 2;
        let precomputations = simultaneous_multiple_scalar_multiplication_create_precomputations(
            window_width,
            &p.into(),
            &phi_p.into(),
            u.sign(),
            v.sign(),
        );

        assert_eq!(
            simultaneous_multiple_scalar_multiplication(window_width, &u, &v, precomputations),
            p.mul(k)
        ); // the GLV method computes [k]p
    }
}
