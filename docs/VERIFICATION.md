# Verification boundary

This note records exactly what is mathematically proved in the manuscript, what is computationally cross-checked in Rust, and what remains unfinished in Lean.

## Layer 1 — mathematical manuscript

The manuscript is the source of the mathematical proofs. Its core claims include:

1. every kernel polynomial is monic of degree `2^m - 1`;
2. the Boolean-indexed kernel polynomials form a basis of `F[X]_{<2^m}`;
3. kernel-to-monomial coordinates are a Boolean zeta transform composed with complement;
4. Boolean Möbius inversion recovers kernel coordinates;
5. `deg U < 2^k` is equivalent to the high-coordinate sign-character condition on kernel coefficients;
6. in characteristic two, that sign condition reduces to high-fiber constancy.

The manuscript explicitly leaves FRI folding/proximity results for future work.

## Layer 2 — Rust computational checks

The Rust implementation is an executable cross-check over explicit modular arithmetic.

Committed tests cover:

- coefficient construction vs the closed coefficient formula for `m = 4`;
- fast zeta vs naive zeta for `m <= 10`;
- fast Möbius vs naive Möbius for `m <= 10`;
- zeta/Möbius round trips for `m <= 12`;
- the forward filtration direction for every `k <= m`, `m <= 8`;
- the converse filtration direction for every `k <= m`, `m <= 8`;
- deterministic seed sweeps for transform and filtration behavior;
- zero-dimensional transform behavior;
- the `k = m` boundary;
- characteristic-two high-fiber constancy;
- invalid dimension/table/modulus guards.

These tests are reproducible finite computations. They are not a proof for arbitrary `m` or arbitrary fields.

## Layer 3 — Lean 4 formalization

`KernelBasisFiltration.lean` is intentionally marked as **in progress**.

At the current repository stage, four central theorem bodies are admitted with `sorry`:

- `kernelPoly_natDegree`;
- `kernelPoly_monic`;
- `kernelPoly_coeff`;
- `kernelPoly_linearIndependent`.

The basis packaging and full low-degree filtration theorem are not yet formalized.

Therefore:

- `lake build` checks that the current definitions/statements elaborate against the pinned Lean/Mathlib environment;
- it does **not** certify the manuscript;
- the repository must not be described as a completed machine-checked proof while these admissions remain.

## CI interpretation

A green Rust job means formatting, compilation, tests, strict Clippy, and release benchmark compilation succeeded.

A green Lean job means the current Lean development builds in its pinned environment. It does not imply `sorry`-free formal verification.
