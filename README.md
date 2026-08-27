# Boolean Kernel Basis Filtration

[![CI](https://github.com/qoosmo/kernel-basis-filtration/actions/workflows/ci.yml/badge.svg)](https://github.com/qoosmo/kernel-basis-filtration/actions/workflows/ci.yml)
[![Lean 4](https://img.shields.io/badge/formalization-Lean%204-blue.svg)](KernelBasisFiltration.lean)
[![Rust](https://img.shields.io/badge/computational%20checks-Rust-orange.svg)](rust/)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/software%20license-MIT%20OR%20Apache--2.0-blue.svg)](#license)

**Research artifact for _“The Boolean Kernel Basis and Its Low-Degree Filtration over Arbitrary Fields”_ by Ali Mkhida (Algorizk Labs).**

This repository studies a Boolean-indexed polynomial basis in which every basis element has maximal ordinary degree, then characterizes exactly how low-degree polynomials appear in those coordinates.

- **Paper:** [`docs/paper.pdf`](docs/paper.pdf)
- **LaTeX source:** [`docs/paper.tex`](docs/paper.tex)
- **Rust computational checks:** [`rust/`](rust/)
- **Lean 4 formalization:** [`KernelBasisFiltration.lean`](KernelBasisFiltration.lean)
- **Verification notes:** [`docs/VERIFICATION.md`](docs/VERIFICATION.md)

## Main result

Let \(N = 2^m\), let \(F\) be an arbitrary field, and index Boolean vectors by
\(y=(y_0,\ldots,y_{m-1}) \in \{0,1\}^m\). Define

\[
K_y(X)=\prod_{i=0}^{m-1}\left(X^{2^i}+y_i\right).
\]

The paper proves three structural facts.

### 1. Boolean kernel basis

The family

\[
\{K_y : y\in\{0,1\}^m\}
\]

is a basis of \(F[X]_{<N}\), even though every \(K_y\) has degree exactly \(N-1\).

This is deliberately unlike a degree-graded basis: low degree cannot be read off by simply asking which kernel coefficients vanish.

### 2. Exact change of basis

For

\[
U(X)=\sum_y \lambda_y K_y(X)=\sum_a u_a X^{|a|_2},
\]

the monomial coefficients satisfy

\[
u_a=\sum_{y\ge \bar a}\lambda_y.
\]

Thus kernel coordinates and monomial coordinates are related by a Boolean zeta transform composed with the complement permutation. The inverse is Boolean Möbius inversion. Both admit \(O(N\log N)\) butterfly algorithms.

### 3. Low-degree filtration in kernel coordinates

Split a kernel index as \(y=(x,h)\), where \(x\in\{0,1\}^k\) and
\(h\in\{0,1\}^{q}\) with \(q=m-k\). Then

\[
\deg U < 2^k
\quad\Longleftrightarrow\quad
\lambda_{(x,h)}
=
(-1)^{q-\operatorname{wt}(h)}\mu_x
\]

for a unique function \(\mu:\{0,1\}^k\to F\).

In characteristic two, \(-1=1\), so this becomes constancy across each high-coordinate fiber:

\[
\lambda_{(x,h)}=\mu_x.
\]

## Why this artifact exists

The research question is representation-theoretic: **what does an ordinary degree bound look like after changing from monomial coordinates to a basis whose individual elements are all high-degree?**

The result is relevant to later work on FRI-style folding because folding is fundamentally tied to degree reduction. This repository does **not** claim an FRI folding theorem, a proximity theorem, or a polynomial commitment construction. Those are separate research questions.

## Repository map

| Path | Role |
| --- | --- |
| `docs/paper.pdf` | Human-readable manuscript |
| `docs/paper.tex` | Reproducible paper source |
| `docs/references.bib` | Bibliography |
| `rust/src/lib.rs` | Zeta/Möbius transforms, coefficient formula, filtration checks |
| `rust/src/bin/bench.rs` | Naive-vs-fast timing harness |
| `KernelBasisFiltration.lean` | In-progress Lean/Mathlib formalization |
| `docs/VERIFICATION.md` | Exact computational/formal verification boundary |
| `.github/workflows/ci.yml` | Rust and Lean build gates |

## Verification status

### Rust

The Rust artifact contains deterministic executable checks for:

- the exact coefficient formula;
- naive \(O(N^2)\) vs fast \(O(N\log N)\) zeta transforms;
- naive vs fast Möbius inversion;
- zeta/Möbius round trips;
- both directions of the low-degree filtration theorem over the Goldilocks prime field;
- boundary cases, including \(m=0\), \(k=m\), invalid dimensions, and characteristic two.

The committed tests are cross-checks, not a formal proof of the theorem.

### Lean 4

The Lean file currently formalizes the definitions and states central results for:

- kernel polynomial degree;
- monicity;
- the coefficient formula;
- linear independence.

Those central proofs still contain `sorry`, and the full filtration theorem has not yet been formalized. A successful Lean CI build therefore means that the current declarations elaborate and type-check **with admitted proof obligations**; it is not evidence of a completed machine-checked proof.

## Quick start

### Rust checks

```sh
cargo fmt --manifest-path rust/Cargo.toml --all -- --check
cargo test --manifest-path rust/Cargo.toml --all-features
cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings
```

Run the timing harness with:

```sh
cargo run --release --manifest-path rust/Cargo.toml --bin kernel-basis-bench
```

The harness prints measurements for the machine on which it is run. The repository does not treat one machine's timings as universal performance claims.

### Lean 4

The Lean toolchain is pinned in [`lean-toolchain`](lean-toolchain), and the Mathlib revision is pinned in [`lakefile.toml`](lakefile.toml).

```sh
lake update
lake exe cache get
lake build
```

### Paper

With a LaTeX installation containing `elsarticle`, TikZ, and BibTeX:

```sh
cd docs
pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
```

## Rust implementation scope

The Rust code intentionally uses a small modular-arithmetic type so that the transforms remain easy to inspect.

It is **research code**, not a production finite-field library:

- the modulus is assumed to be prime when field semantics are required;
- arithmetic is not constant-time;
- there is no side-channel hardening;
- there has been no security audit;
- the implementation is not a SNARK/STARK or polynomial commitment system.

For production proving-system work, use an audited field library and treat this crate as a reference artifact.

## Reproducibility

The repository pins the Lean toolchain and Mathlib revision used by this artifact. Rust has no third-party dependencies and commits `Cargo.lock`.

CI independently checks:

1. Rust formatting;
2. Rust compilation;
3. Rust tests;
4. warning-free Clippy;
5. release compilation of the benchmark harness;
6. Lean/Mathlib build of the current formalization.

See [`docs/VERIFICATION.md`](docs/VERIFICATION.md) for the exact boundary between mathematical proof, computational testing, and unfinished formalization.

## Research scope

This repository establishes and checks the kernel-basis algebra described above. It does not claim:

- a Fiat–Shamir transcript;
- a polynomial commitment scheme;
- a FRI soundness or folding theorem;
- a SNARK/STARK implementation;
- zero-knowledge masking;
- production cryptographic security;
- a completed Lean proof.

## Author

**Ali Mkhida** — Algorizk Labs
ORCID: [`0009-0009-2101-9070`](https://orcid.org/0009-0009-2101-9070)

## License

The **software and formalization source** in this repository is available under either the MIT License or Apache License 2.0, at your option; see [`LICENSE-MIT`](LICENSE-MIT) and [`LICENSE-APACHE`](LICENSE-APACHE).

The manuscript files under `docs/` are scholarly works and are **not granted under those software licenses**. Copyright and publication rights for the manuscript remain with the author unless a separate license is stated.
