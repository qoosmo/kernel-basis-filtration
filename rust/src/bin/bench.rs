use kernel_basis_bench::*;
use std::time::Instant;

fn rand_field_vec(n: usize, seed: u64, p: u64) -> Vec<F> {
    let mut state = seed ^ 0x9E3779B97F4A7C15;
    (0..n)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            F::new((state >> 1) as i128, p)
        })
        .collect()
}

fn time_it<T>(f: impl FnOnce() -> T) -> (T, f64) {
    let start = Instant::now();
    let out = f();
    (out, start.elapsed().as_secs_f64())
}

fn main() {
    let p = GOLDILOCKS_P;

    println!("Zeta transform: naive O(N^2) vs fast O(N log N), field = Goldilocks (2^64 - 2^32 + 1)");
    println!("{:>3} {:>10} {:>14} {:>14} {:>10}", "m", "N", "naive (s)", "fast (s)", "speedup");
    println!("{}", "-".repeat(56));

    // Naive is O(N^2); keep it to sizes that finish in reasonable time.
    let naive_max_m = 14u32;
    let fast_max_m = 22u32;

    for m in 1..=fast_max_m {
        let n = 1usize << m;
        let lambda = rand_field_vec(n, 42, p);

        let (u_fast, t_fast) = time_it(|| zeta_fast(&lambda, m, p));

        if m <= naive_max_m {
            let (u_naive, t_naive) = time_it(|| zeta_naive(&lambda, m, p));
            assert_eq!(u_naive, u_fast, "naive/fast disagreement at m={m}");
            let speedup = if t_fast > 0.0 { t_naive / t_fast } else { f64::INFINITY };
            println!(
                "{:>3} {:>10} {:>14.6} {:>14.6} {:>9.1}x",
                m, n, t_naive, t_fast, speedup
            );
        } else {
            println!("{:>3} {:>10} {:>14} {:>14.6} {:>10}", m, n, "(skipped)", t_fast, "-");
        }
    }

    println!();
    println!("Mobius inversion: naive O(N^2) vs fast O(N log N)");
    println!("{:>3} {:>10} {:>14} {:>14} {:>10}", "m", "N", "naive (s)", "fast (s)", "speedup");
    println!("{}", "-".repeat(56));

    for m in 1..=fast_max_m {
        let n = 1usize << m;
        let u = rand_field_vec(n, 7, p);

        let (lam_fast, t_fast) = time_it(|| mobius_fast(&u, m, p));

        if m <= naive_max_m {
            let (lam_naive, t_naive) = time_it(|| mobius_naive(&u, m, p));
            assert_eq!(lam_naive, lam_fast, "naive/fast disagreement at m={m}");
            let speedup = if t_fast > 0.0 { t_naive / t_fast } else { f64::INFINITY };
            println!(
                "{:>3} {:>10} {:>14.6} {:>14.6} {:>9.1}x",
                m, n, t_naive, t_fast, speedup
            );
        } else {
            println!("{:>3} {:>10} {:>14} {:>14.6} {:>10}", m, n, "(skipped)", t_fast, "-");
        }
    }

    println!();
    println!("(naive skipped above m={naive_max_m} -- O(N^2) makes it impractical; fast continues to m={fast_max_m}, N={}.)", 1usize << fast_max_m);
}
