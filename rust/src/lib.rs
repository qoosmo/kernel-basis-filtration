//! Core algorithms for the Boolean kernel basis (companion to "The Boolean
//! Kernel Basis and Its Low-Degree Filtration over Arbitrary Fields").
//!
//! Provides naive O(N^2) and fast O(N log N) implementations of the
//! Boolean zeta transform and its Moebius inversion (Section 4.4 of the
//! paper), the exact kernel coefficient formula (Proposition 4.3), and a
//! direct check of the low-degree filtration character condition
//! (Theorem 5.2), so the two can be benchmarked against each other and
//! cross-checked for correctness.
//!
//! Indexing convention: a Boolean vector y in {0,1}^m is represented as an
//! integer in 0..N (N = 2^m) via bit i <-> coordinate y_i, exactly
//! matching |a|_2 = sum_i a_i 2^i from the paper. Complement is the
//! bitwise NOT restricted to m bits, `(N - 1) ^ a`.
//!
//! This crate is a research reference implementation, not a constant-time
//! or audited finite-field library.

#![forbid(unsafe_code)]

use std::ops::{Add, Mul, Neg, Sub};

/// The Goldilocks prime, 2^64 - 2^32 + 1.
pub const GOLDILOCKS_P: u64 = 0xFFFF_FFFF_0000_0001;

/// An element of the prime field GF(p). Multiplication uses a u128
/// intermediate so this is correct (not necessarily fastest) for any
/// prime that fits in a u64.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct F {
    v: u64,
    p: u64,
}

impl F {
    /// Construct a canonical residue modulo `p`.
    ///
    /// The caller is responsible for choosing a prime modulus when field
    /// semantics are required.
    pub fn new(v: i128, p: u64) -> Self {
        assert!(p >= 2, "modulus must be at least 2");
        let modulus = p as i128;
        let r = ((v % modulus) + modulus) % modulus;
        Self { v: r as u64, p }
    }

    pub fn zero(p: u64) -> Self {
        assert!(p >= 2, "modulus must be at least 2");
        Self { v: 0, p }
    }

    pub fn one(p: u64) -> Self {
        assert!(p >= 2, "modulus must be at least 2");
        Self { v: 1, p }
    }

    #[inline]
    pub fn value(self) -> u64 {
        self.v
    }

    #[inline]
    pub fn modulus(self) -> u64 {
        self.p
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.v == 0
    }

    #[inline]
    fn assert_same_modulus(self, other: Self) {
        assert_eq!(
            self.p, other.p,
            "field modulus mismatch: {} != {}",
            self.p, other.p
        );
    }
}

impl Add for F {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        let sum = self.v as u128 + rhs.v as u128;
        Self {
            v: (sum % self.p as u128) as u64,
            p: self.p,
        }
    }
}

impl Sub for F {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        let v = if self.v >= rhs.v {
            self.v - rhs.v
        } else {
            self.p - (rhs.v - self.v)
        };
        Self { v, p: self.p }
    }
}

impl Mul for F {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self.assert_same_modulus(rhs);
        let product = self.v as u128 * rhs.v as u128;
        Self {
            v: (product % self.p as u128) as u64,
            p: self.p,
        }
    }
}

impl Neg for F {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.v == 0 {
            self
        } else {
            Self {
                v: self.p - self.v,
                p: self.p,
            }
        }
    }
}

#[inline]
fn domain_size(m: u32) -> usize {
    assert!(
        m < usize::BITS,
        "m must be smaller than usize::BITS ({})",
        usize::BITS
    );
    1usize << m
}

fn validate_table(values: &[F], m: u32, p: u64, name: &str) {
    let expected = domain_size(m);
    assert_eq!(
        values.len(),
        expected,
        "{name} must have length 2^m = {expected}"
    );
    assert!(
        values.iter().all(|value| value.modulus() == p),
        "{name} contains an element with a different modulus"
    );
}

/// Hamming weight of an m-bit integer (only the low `m` bits are meant to
/// be set, but this just counts whatever bits are present).
#[inline]
pub fn wt(y: usize) -> u32 {
    y.count_ones()
}

/// Boolean complement of `a` within `m` bits.
#[inline]
pub fn comp(a: usize, m: u32) -> usize {
    let n = domain_size(m);
    assert!(a < n, "Boolean index a={a} is outside 0..{n}");
    (n - 1) ^ a
}

/// Coefficient vector (length N = 2^m) of the kernel polynomial
/// K_y(X) = prod_i (X^{2^i} + y_i), built via the recursive even/odd
/// construction of Section 4.3: O(N log N) time.
pub fn kernel_poly_coeffs(y: usize, m: u32, p: u64) -> Vec<F> {
    let n = domain_size(m);
    assert!(y < n, "Boolean index y={y} is outside 0..{n}");
    let mut c = vec![F::zero(p); n];
    c[0] = F::one(p);
    let mut deg = 1usize;
    for i in 0..m {
        let yi = (y >> i) & 1;
        let mut new_c = vec![F::zero(p); n];
        for j in 0..deg {
            if yi == 1 {
                new_c[j] = new_c[j] + c[j];
            }
            new_c[j + (1 << i)] = new_c[j + (1 << i)] + c[j];
        }
        c = new_c;
        deg <<= 1;
    }
    c
}

/// Exact coefficient formula, Proposition 4.3/2.2 (complement-corrected):
/// [K_y]_a = 1 iff y >= comp(a), else 0.
pub fn kernel_coeff_formula(y: usize, a: usize, m: u32, p: u64) -> F {
    let n = domain_size(m);
    assert!(y < n, "Boolean index y={y} is outside 0..{n}");
    assert!(a < n, "Boolean index a={a} is outside 0..{n}");
    let ab = comp(a, m);
    if (y & ab) == ab {
        F::one(p)
    } else {
        F::zero(p)
    }
}

/// Naive O(N^2) zeta transform composed with the complement permutation:
/// u[a] = sum_{y >= comp(a)} lambda[y]. Reference implementation, and the
/// slow half of the benchmark.
pub fn zeta_naive(lambda: &[F], m: u32, p: u64) -> Vec<F> {
    validate_table(lambda, m, p, "lambda");
    let n = domain_size(m);
    let mut u = vec![F::zero(p); n];
    for (a, ua) in u.iter_mut().enumerate() {
        let ab = comp(a, m);
        let mut s = F::zero(p);
        for (y, &ly) in lambda.iter().enumerate() {
            if (y & ab) == ab {
                s = s + ly;
            }
        }
        *ua = s;
    }
    u
}

/// Fast O(N log N) zeta transform composed with the complement
/// permutation, via the standard "sum over supersets" butterfly schedule.
pub fn zeta_fast(lambda: &[F], m: u32, p: u64) -> Vec<F> {
    validate_table(lambda, m, p, "lambda");
    let n = domain_size(m);
    let mut f = lambda.to_vec();
    for i in 0..m {
        let bit = 1usize << i;
        for mask in 0..n {
            if mask & bit == 0 {
                let hi = f[mask | bit];
                f[mask] = f[mask] + hi;
            }
        }
    }
    let mut u = vec![F::zero(p); n];
    for a in 0..n {
        u[a] = f[comp(a, m)];
    }
    u
}

/// Naive O(N^2) Moebius inversion:
/// lambda[y] = sum_{b>=y} (-1)^{wt(b)-wt(y)} u[comp(b)].
pub fn mobius_naive(u: &[F], m: u32, p: u64) -> Vec<F> {
    validate_table(u, m, p, "u");
    let n = domain_size(m);
    let mut lambda = vec![F::zero(p); n];
    for (y, lambda_y) in lambda.iter_mut().enumerate() {
        let mut s = F::zero(p);
        for b in 0..n {
            if (b & y) == y {
                let term = u[comp(b, m)];
                s = if (wt(b) - wt(y)).is_multiple_of(2) {
                    s + term
                } else {
                    s - term
                };
            }
        }
        *lambda_y = s;
    }
    lambda
}

/// Fast O(N log N) Moebius inversion via the inverse butterfly schedule.
pub fn mobius_fast(u: &[F], m: u32, p: u64) -> Vec<F> {
    validate_table(u, m, p, "u");
    let n = domain_size(m);
    let mut h = vec![F::zero(p); n];
    for (b, hb) in h.iter_mut().enumerate() {
        *hb = u[comp(b, m)];
    }
    for i in 0..m {
        let bit = 1usize << i;
        for mask in 0..n {
            if mask & bit == 0 {
                let hi = h[mask | bit];
                h[mask] = h[mask] - hi;
            }
        }
    }
    h
}

/// Degree of a monomial-coefficient vector: highest index with a nonzero
/// coefficient, or `None` for the zero polynomial.
pub fn degree(u: &[F]) -> Option<usize> {
    u.iter().rposition(|c| !c.is_zero())
}

/// Theorem 5.2 check: does `lambda` (full kernel-coefficient table, N =
/// 2^m entries, low `k` bits = x, high `m-k` bits = h) satisfy the
/// character condition for this k?
pub fn satisfies_character_condition(lambda: &[F], m: u32, k: u32, p: u64) -> bool {
    assert!(k <= m, "k must satisfy k <= m");
    validate_table(lambda, m, p, "lambda");
    let q = m - k;
    let low_n = 1usize << k;
    let high_n = 1usize << q;
    for x in 0..low_n {
        let mu_x = lambda[x | ((high_n - 1) << k)]; // h = all-ones (1_q)
        for h in 0..high_n {
            let y = x | (h << k);
            let expected = if (q - wt(h)) % 2 == 1 { -mu_x } else { mu_x };
            if lambda[y] != expected {
                return false;
            }
        }
    }
    true
}

/// Build the kernel-coefficient table satisfying the character condition
/// for the given `mu : {0,1}^k -> F` (i.e. the forward direction of
/// Theorem 5.2: lambda[(x,h)] = (-1)^{q-wt(h)} mu[x]).
pub fn character_table_from_mu(mu: &[F], m: u32, k: u32, p: u64) -> Vec<F> {
    assert!(k <= m, "k must satisfy k <= m");
    let q = m - k;
    let n = domain_size(m);
    let low_n = domain_size(k);
    let high_n = domain_size(q);
    assert_eq!(mu.len(), low_n, "mu must have length 2^k = {low_n}");
    assert!(
        mu.iter().all(|value| value.modulus() == p),
        "mu contains an element with a different modulus"
    );
    let mut lambda = vec![F::zero(p); n];
    for (x, &mux) in mu.iter().enumerate() {
        for h in 0..high_n {
            let y = x | (h << k);
            lambda[y] = if (q - wt(h)) % 2 == 1 { -mux } else { mux };
        }
    }
    lambda
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = GOLDILOCKS_P;

    fn rand_field_vec(n: usize, seed: u64) -> Vec<F> {
        // Small deterministic LCG, good enough for test data.
        let mut state = seed ^ 0x9E3779B97F4A7C15;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                F::new((state >> 1) as i128, P)
            })
            .collect()
    }

    #[test]
    fn kernel_poly_matches_coeff_formula_m4() {
        let m = 4u32;
        let n = 1usize << m;
        for y in 0..n {
            let coeffs = kernel_poly_coeffs(y, m, P);
            for (a, &coeff) in coeffs.iter().enumerate() {
                assert_eq!(
                    coeff,
                    kernel_coeff_formula(y, a, m, P),
                    "mismatch at y={y}, a={a}"
                );
            }
        }
    }

    #[test]
    fn zeta_fast_matches_naive_m_up_to_10() {
        for m in 1..=10u32 {
            let n = 1usize << m;
            let lambda = rand_field_vec(n, 12345 + m as u64);
            assert_eq!(zeta_fast(&lambda, m, P), zeta_naive(&lambda, m, P), "m={m}");
        }
    }

    #[test]
    fn mobius_fast_matches_naive_m_up_to_10() {
        for m in 1..=10u32 {
            let n = 1usize << m;
            let u = rand_field_vec(n, 54321 + m as u64);
            assert_eq!(mobius_fast(&u, m, P), mobius_naive(&u, m, P), "m={m}");
        }
    }

    #[test]
    fn mobius_inverts_zeta_round_trip() {
        for m in 1..=12u32 {
            let n = 1usize << m;
            let lambda = rand_field_vec(n, 999 + m as u64);
            let u = zeta_fast(&lambda, m, P);
            let back = mobius_fast(&u, m, P);
            assert_eq!(back, lambda, "round trip failed at m={m}");
        }
    }

    #[test]
    fn theorem_5_2_forward_direction_all_k() {
        // Build lambda from a character condition for each k, verify the
        // resulting polynomial's monomial degree is < 2^k.
        for m in 1..=8u32 {
            for k in 0..=m {
                let low_n = 1usize << k;
                let mu = rand_field_vec(low_n, 7 + m as u64 * 100 + k as u64);
                let lambda = character_table_from_mu(&mu, m, k, P);
                let u = zeta_fast(&lambda, m, P);
                let deg = degree(&u);
                let bound = 1usize << k;
                match deg {
                    None => {} // zero polynomial trivially has degree < bound
                    Some(d) => assert!(d < bound, "m={m} k={k}: degree {d} not < {bound}"),
                }
            }
        }
    }

    #[test]
    fn theorem_5_2_converse_direction_all_k() {
        // Build a random low-degree U directly in monomial coefficients,
        // Mobius-invert to kernel coordinates, check the character
        // condition holds exactly.
        for m in 1..=8u32 {
            for k in 0..=m {
                let n = 1usize << m;
                let bound = 1usize << k;
                let mut u = vec![F::zero(P); n];
                let rvals = rand_field_vec(bound, 55 + m as u64 * 100 + k as u64);
                u[..bound].copy_from_slice(&rvals);
                let lambda = mobius_fast(&u, m, P);
                assert!(
                    satisfies_character_condition(&lambda, m, k, P),
                    "m={m} k={k}: character condition failed"
                );
            }
        }
    }

    // Deterministic seed sweeps provide additional reproducible coverage
    // without adding external test-framework dependencies.
    #[test]
    fn zeta_mobius_round_trip_many_seeds() {
        let m = 6u32;
        let n = 1usize << m;
        for seed in 0u64..200 {
            let lambda = rand_field_vec(n, seed);
            let u = zeta_fast(&lambda, m, P);
            let back = mobius_fast(&u, m, P);
            assert_eq!(back, lambda, "seed={seed}");
        }
    }

    #[test]
    fn character_condition_forward_many_seeds() {
        let m = 6u32;
        for seed in 0u64..50 {
            for k in 0..=m {
                let low_n = 1usize << k;
                let mu = rand_field_vec(low_n, seed * 1000 + k as u64);
                let lambda = character_table_from_mu(&mu, m, k, P);
                let u = zeta_fast(&lambda, m, P);
                let bound = 1usize << k;
                if let Some(d) = degree(&u) {
                    assert!(d < bound, "seed={seed} k={k}: degree {d} not < {bound}");
                }
            }
        }
    }

    #[test]
    fn zero_dimension_transforms_are_identity() {
        let lambda = vec![F::new(7, P)];
        assert_eq!(zeta_fast(&lambda, 0, P), lambda);
        assert_eq!(mobius_fast(&lambda, 0, P), lambda);
        assert_eq!(kernel_poly_coeffs(0, 0, P), vec![F::one(P)]);
    }

    #[test]
    fn characteristic_two_condition_is_constant_on_high_fibers() {
        let p = 2;
        let m = 3;
        let k = 1;
        let mu = vec![F::zero(p), F::one(p)];
        let lambda = character_table_from_mu(&mu, m, k, p);

        for x in 0..(1usize << k) {
            for h in 0..(1usize << (m - k)) {
                assert_eq!(lambda[x | (h << k)], mu[x]);
            }
        }

        assert!(satisfies_character_condition(&lambda, m, k, p));
        let u = zeta_fast(&lambda, m, p);
        if let Some(d) = degree(&u) {
            assert!(d < (1usize << k));
        }
    }

    #[test]
    fn k_equals_m_accepts_an_arbitrary_kernel_table() {
        let m = 5;
        let lambda = rand_field_vec(1usize << m, 2026);
        assert!(satisfies_character_condition(&lambda, m, m, P));
    }

    #[test]
    #[should_panic(expected = "k must satisfy k <= m")]
    fn rejects_k_larger_than_m() {
        let lambda = vec![F::zero(P); 4];
        let _ = satisfies_character_condition(&lambda, 2, 3, P);
    }

    #[test]
    #[should_panic(expected = "lambda must have length")]
    fn rejects_wrong_transform_table_length() {
        let lambda = vec![F::zero(P); 3];
        let _ = zeta_fast(&lambda, 2, P);
    }

    #[test]
    #[should_panic(expected = "different modulus")]
    fn rejects_mixed_modulus_table() {
        let lambda = vec![F::zero(P), F::zero(17)];
        let _ = zeta_fast(&lambda, 1, P);
    }

    #[test]
    #[should_panic(expected = "field modulus mismatch")]
    fn arithmetic_rejects_cross_modulus_operands() {
        let _ = F::one(17) + F::one(19);
    }
}
