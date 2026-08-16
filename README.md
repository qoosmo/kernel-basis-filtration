# kernel-basis-filtration

Public research artifact for **\"The Boolean Kernel Basis and Its Low-Degree Filtration over Arbitrary Fields\"** by Ali Mkhida (Algorizk Labs).

The paper studies a Boolean-indexed polynomial basis over arbitrary fields and gives an exact characterization of the ordinary degree filtration in that basis.

- **Paper:** [`docs/paper.pdf`](docs/paper.pdf)
- **Source:** [`docs/paper.tex`](docs/paper.tex)
- **Rust implementation and tests:** [`rust/`](rust/)
- **Lean 4 formalization (in progress):** [`KernelBasisFiltration.lean`](KernelBasisFiltration.lean)

## Verification status

- **Rust:** committed executable cross-checks for the coefficient formula, Boolean zeta/Möbius transforms, and both directions of the filtration theorem over the Goldilocks prime field. The repository also contains a benchmark harness.
- **Lean 4:** formalization is in progress; proofs are not yet complete and are not presented as a machine-checked proof of the paper.
- **Paper:** the mathematical arguments are self-contained and do not depend on either computational artifact.

## Building

### Rust

```sh
cd rust
cargo test --release
cargo run --release --bin bench
```

### Lean 4

```sh
lake update
lake build
```

### Paper

```sh
cd docs
pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
```
