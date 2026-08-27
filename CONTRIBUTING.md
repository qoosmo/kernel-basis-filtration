# Contributing

This repository is a research artifact. Changes should preserve the distinction between the mathematical manuscript, executable cross-checks, and incomplete formal verification.

## Before opening a pull request

Run the Rust quality gates:

```sh
cargo fmt --manifest-path rust/Cargo.toml --all -- --check
cargo check --manifest-path rust/Cargo.toml --all-targets --all-features
cargo test --manifest-path rust/Cargo.toml --all-features
cargo clippy --manifest-path rust/Cargo.toml --all-targets --all-features -- -D warnings
cargo build --release --manifest-path rust/Cargo.toml --bin kernel-basis-bench
```

Run the Lean build:

```sh
lake update
lake exe cache get
lake build
```

## Research claims

When changing a theorem statement, transform convention, index ordering, or benchmark claim:

1. update the manuscript/source of truth if appropriate;
2. update the Rust checks that exercise the same statement;
3. update the Lean statement when it has a corresponding declaration;
4. state clearly when a result is experimentally checked rather than proved;
5. do not describe admitted Lean declarations as machine-checked proofs.

## Style

- Rust must pass `rustfmt` and Clippy with `-D warnings`.
- Prefer deterministic tests.
- Keep asymptotic claims separate from machine-specific timings.
- Avoid introducing proving-system claims outside the scope of this artifact.
