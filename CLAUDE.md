# Engine Bench Engineering Standards

Engine Bench, built in Rust: physical simulation and animated 2D visualization of 38 engines, from the aeolipile to rotating detonation and ion thrusters, in real-time and offline modes. Architecture and roadmap: `docs/plan/0001-architecture-and-roadmap.md`. Keep this file in sync as the project grows.

## Project Structure

Cargo workspace with 3 members:

```
engine_bench/                          (workspace root)
├── Cargo.toml                         [workspace] manifest — all dependency versions here
├── docs/plan/                         Approved plans
│
└── crates/
    ├── engine-bench-physics/          (library — units, thermo, chemistry, numerics, solvers)
    │   ├── build.rs                   Parses data/thermo.inp → const species tables (OUT_DIR)
    │   ├── data/thermo.inp            Vendored NASA Glenn 9-coefficient records (from NASA CEA)
    │   ├── benches/kernels.rs         criterion benchmarks of hot kernels (baseline, not a gate)
    │   ├── proptest-regressions/      Saved proptest failure seeds (committed, re-run first)
    │   └── src/
    │       ├── lib.rs                 Module declarations
    │       ├── error.rs               PhysicsError (thiserror), carried inside anyhow::Error
    │       ├── flow/                  Perfect-gas compressible flow (γ frozen upstream)
    │       │   ├── isentropic.rs      Ratios vs Mach, A/A*, M(A/A*), Prandtl–Meyer
    │       │   ├── shock.rs           Normal-shock relations
    │       │   ├── orifice.rs         Reservoir, Orifice: bidirectional, choked, C¹ near Δp = 0
    │       │   └── nozzle/            Nozzle operating point: every quasi-1D regime, thrust, Isp
    │       ├── units.rs               uom SI quantities for public APIs
    │       ├── numerics/
    │       │   ├── roots.rs           Brent, bracketed Newton
    │       │   └── ode/               OdeSystem, rk4_step, DormandPrince45 (dense output, events)
    │       └── thermo/
    │           ├── gases.rs           Generated species constants (N2, O2, …, ALL)
    │           ├── nasa9.rs           TemperaturePowers + 9-coefficient polynomial kernels
    │           ├── species.rs         Species, Element; constant-cp extrapolation below range
    │           ├── species_set.rs     SpeciesSet<N> (compile-time), Air3, Major7
    │           └── mixture/           Mixture<S, N>: composition, properties, T from u / h
    │
    ├── engine-bench-core/             (library — engine models, Engine trait, Scene, recording)
    │   └── src/lib.rs                 Module declarations + re-exports
    │
    └── engine-bench/                  (binary)
        └── src/
            ├── main.rs                Entry point
            └── logging.rs             tracing setup (console, RUST_LOG override)
```

**Dependency boundary:** physics ← core ← engine-bench. Never reverse. Physics knows nothing about specific engines. Core has no binary-only dependencies (logging setup, CLI parsing, UI): engines describe what to draw as a renderer-agnostic `Scene`, and the binary renders it.

**Units:** `uom` quantities (re-exported from `engine_bench_physics::units`) in all public APIs. Hot kernels unwrap to raw `f64` struct-of-arrays and re-wrap on output.

## Planning Process

**Planning is interactive, not autonomous.** When asked to plan a step or feature:

1. Research the codebase and external dependencies. Present findings.
2. **Ask questions** before finalizing the plan — surface ambiguities, trade-offs, API changes, and design decisions that need the user's input.
3. Only finalize and save the plan after the user has reviewed and approved the approach.
4. Save the final plan to `docs/plan/` before implementation begins.

Do not silently make architectural decisions. If the implementation plan document conflicts with the current codebase state (e.g., outdated dependency versions, changed APIs), flag the discrepancy and ask how to proceed.

## Quick Reference

```bash
cargo nextest run --workspace                 # All workspace tests (use nextest, not cargo test)
cargo nextest run -p engine-bench-physics     # Physics-only tests
cargo nextest run -p engine-bench-core        # Core-only tests
cargo clippy --workspace --all-targets        # Zero warnings required (includes tests)
cargo check --workspace --all-targets         # Type-check everything including tests
cargo fmt --all -- --check                    # Format check
cargo build -p engine-bench                   # Build the binary
cargo bench -p engine-bench-physics --profile bench-unwind   # Kernel benchmarks (plain `cargo bench` fails: build-std + immediate-abort)
```

**Always use `cargo nextest run` instead of `cargo test`.** Nextest is the project's test runner.

**Toolchain:** nightly (`rust-toolchain.toml`). The release profile uses `panic = "immediate-abort"` and `.cargo/config.toml` rebuilds std (`build-std`) with `target-cpu=native`.

**`--all-targets` is mandatory** for `clippy` and `check`. Lints must pass in tests and examples — not just lib/bin targets. Test code follows the same quality standards as production code.

## Test-Driven Development

**Strict red-green-refactor. No exceptions.**

1. **Red** — Write a failing test first. Run it. Confirm it fails for the right reason.
2. **Green** — Write the minimum code to make the test pass. Nothing more.
3. **Refactor** — Clean up while all tests stay green.

### Test Requirements

| Context | Requirement |
|---------|------------|
| Every public function | At least one unit test |
| Data processing / numeric logic | Property-based tests (`proptest`, 1000+ cases) |
| Async flows | Integration test with `#[tokio::test]` + `tokio::time::timeout` |
| API client changes | `wiremock`-based integration tests |
| `thiserror` enum variants | Test verifying each variant's Display output |

### Test Conventions

- **Naming:** `test_<unit>_<scenario>` (e.g., `test_parser_rejects_empty_input`)
- **Location:** `#[cfg(test)] mod tests` inline in the source file. Integration tests in `crates/<crate>/tests/`.
- **Assertions:** `assert_eq!`, `assert!`, `prop_assert!`, `assert_relative_eq!` (approx)
- **Errors:** Tests return `anyhow::Result<()>` with `.context()` for diagnostics.
- **Quality:** Test code follows the same lint rules as production code. No `.unwrap()`, `.expect()`, or `panic!()` in tests — use `?` with `anyhow::Result` or `anyhow::bail!`.
- **Proptest config:** `#![proptest_config(ProptestConfig::with_cases(1000))]`
- **Time mocking:** Abstract time behind a `Clock` trait with a mock (`set_time()` / `advance()`) instead of real time.

## Performance Rules

**Runtime performance beats compile time in every trade-off.** Slower builds, heavier generics, more monomorphisation, `lto`/`codegen-units = 1`, and compile-time-heavy crates are all acceptable when they make the program faster or safer at no runtime cost. Never trade runtime speed for faster compiles.

### Mandatory

- **Pre-allocate when size is known.** `String::with_capacity()`, `Vec::with_capacity()`.
- **`rust_decimal` or integer types for financial/percentage math where precision matters.** Use `f64` only where approximate values are acceptable (statistics, display values).
- **No allocation in tight loops.** Data processing and analytics paths should minimize allocations.

### Avoid

- `clone()` when a borrow suffices
- `Box<dyn Trait>` when monomorphic generics work
- `String` for fixed-vocabulary identifiers — use enums or newtypes
- `Arc<Mutex<T>>` when channels or `RwLock` suffice

## Error Handling

- **Public APIs:** Return `anyhow::Result<T>`. Add `.context("msg")` or `.with_context(|| format!(...))` (lazy) on every `?`.
- **Domain errors:** `thiserror` enums for recoverable, matchable errors (e.g. `AppError`).
- **Retryable errors:** Expose `is_retryable()` on error enums that cover transient failures.
- **Forbidden everywhere (including tests):** `.unwrap()`, `.expect()`, `panic!()`, `todo!()` — enforced by lints with `--all-targets`.
- **Fallible constructors:** Return `Result<Self>` not `Self`. Validate inputs at construction.
- **Hot-path kernels are infallible.** Property evaluation, fluxes and ODE derivatives run millions of times per second and take already-validated inputs; they return plain values, not `Result`. Validation happens at the edges: constructors and solver steps (`step()`, `solve()`) return `Result` and check for out-of-range or non-finite state.

## Coding Standards

### Lints (Cargo.toml, applied to all targets)

```
unsafe_code       = "forbid"
unwrap_used       = "deny"
expect_used       = "deny"
panic             = "deny"
todo              = "deny"
print_stdout      = "warn"
print_stderr      = "warn"
clippy::pedantic  = "warn" (priority -1)
```

These lints apply to **all code** — lib, bin, tests. Run `cargo clippy --all-targets` to verify.

### Formatting (`rustfmt.toml`)

- `max_width = 100`
- `imports_granularity = "Crate"`
- `group_imports = "StdExternalCrate"`
- `wrap_comments = true`
- `format_code_in_doc_comments = true`
- Edition 2024

### Import Order

```rust
use std::...;

use anyhow::...;       // external crates
use serde::...;

use crate::...;        // local
use super::...;
```

One blank line between groups.

### Visibility

Minimum necessary. `pub(super)` or `pub(crate)` for internal types. Private fields by default.

### Module Organization

- One responsibility per file.
- Split at ~300 lines. Convert to directory module (`mod.rs` + submodules) at ~500 lines.
- Directory modules: `mod.rs` is the coordinator (struct, trait impl, re-exports). Submodules own specific concerns.
- Binary-specific modules live in the binary crate, shared logic in `engine-bench-core`. No feature gates needed.

### Comments

- `///` doc comments on all public items.
- `// SAFETY:` for correctness reasoning (cancellation safety, invariant guarantees, why a conversion is sound).
- No trivial comments restating what code does. Comments explain *why*.

### Type Safety

- Domain types over raw primitives in public APIs where a domain type exists.
- Trait abstractions for testability (e.g. `Clock`). Always accept trait bounds, not concrete types, in functions that need time or other side effects.

### Dependencies

**Rule 1: Always Disable Default Features.** Always explicitly set `default-features = false` when declaring dependencies. Only enable the specific features required.

**Rule 2: Minimum Features Enabled.** Strictly limit feature opt-ins to the absolute bare minimum required for the code to compile and run. Never use blanket features like `features = ["full"]`. This keeps compile times fast, binary sizes small, and the attack surface minimal.

**Rule 3: Deps live in the narrowest crate.** Binary-only dependencies (logging setup, CLI, UI) belong in the binary crate, never in `engine-bench-core` or `engine-bench-physics`. All versions are defined in the workspace root `[workspace.dependencies]`.

**General:**
- New dependencies require justification: what problem, why this crate, what alternatives were considered.
- Prefer zero-cost abstractions over convenience crates.
- No exact (`=`) version pins: every dependency goes through the upgrade workflow below.


### Upgrading Dependencies

On a feature branch (`chore/upgrade-deps`), never directly on `dev`:

```bash
cargo upgrade --dry-run -vv                   # Every dependency: current, compatible and latest version
cargo upgrade                                 # Compatible upgrades (Cargo.toml and Cargo.lock)
cargo upgrade --incompatible                  # Also major / 0.x-minor upgrades: read their changelogs first
```

Then run the full check (fmt, clippy, check, nextest, as in the Pre-Merge Checklist).

Crates whose output feeds results (numeric, benchmarking or measurement crates) can change results silently, even in patch releases: pin their outputs with golden-value tests, and when one fails, compare old and new versions on the same input before updating the golden values. Put the comparison in the PR.

## Tracing & Observability

- **`#[instrument]`** on all public `async fn`. Use `skip(self)` or `skip_all` with explicit `fields(...)`.
- **Field syntax:** `%value` for Display, `?value` for Debug, bare primitives for Copy types.
- **Levels:** `error!` = failures requiring attention. `warn!` = degraded states. `info!` = lifecycle events. `debug!` = protocol detail.
- **Default filter:** `engine_bench_core=DEBUG,engine_bench=DEBUG` (debug), `INFO` (release).
- **Override:** `RUST_LOG` env var via `EnvFilter::try_from_default_env()`.

## Git Workflow

No CI, no hooks. All verification is the developer's responsibility before merging.

### Branch Structure

```
main          Stable release branch. Always clean: compiles, all tests pass, zero warnings.
  └── dev     Integration branch. Features merge here first via PR. Periodically merged to main.
       ├── feature/<name>   New functionality
       ├── fix/<name>       Bug fixes
       ├── refactor/<name>  Code restructuring
       ├── test/<name>      Test additions/improvements
       └── docs/<name>      Documentation changes
```

### Branch Lifecycle

1. **Create** a feature branch from `dev`:
   ```bash
   git checkout dev
   git pull origin dev
   git checkout -b feature/<name>
   ```

2. **Work** on the branch. Make atomic commits (each compiles, passes tests, zero warnings).

3. **Push** the branch and open a PR targeting `dev`:
   ```bash
   git push -u origin feature/<name>
   ```

4. **Verify** locally before merging (no CI — you must run these):
   ```bash
   cargo fmt --all -- --check
   cargo clippy --workspace --all-targets
   cargo nextest run --workspace
   ```

5. **Merge** the PR into `dev` via merge commit (no squash — preserve commit history).

6. **Delete** the feature branch after merge.

7. **Promote** `dev` to `main` when a stable milestone is reached:
   ```bash
   git checkout main
   git merge dev
   git push origin main
   ```

### Branch Rules

- **Never commit directly to `main` or `dev`.** All changes go through feature branches + PRs.
- **Never force-push** to `main` or `dev`.
- **No WIP commits** on `dev` or `main`. Feature branches may have WIP commits but clean them up before PR.
- **New commits over amends** — preserve history.
- **Delete merged branches** to keep the branch list clean.

### Commit Messages

```
<type>: <imperative summary, max 72 chars>

<optional body: explains WHY, not what. Motivation, trade-offs, context.>

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
```

Types: `feat`, `fix`, `refactor`, `test`, `docs`, `chore`

### Pull Requests

All PRs target `dev` unless it is a `dev` → `main` promotion.

#### PR Format

```
Title: <type>: <short description> (under 72 chars)

## Summary
- Bullet points of what changed and why

## Test Plan
- [ ] Tests written first (red-green-refactor)
- [ ] cargo fmt --all -- --check passes
- [ ] cargo nextest run passes
- [ ] cargo clippy --all-targets clean
- [ ] cargo check --all-targets clean
```

#### Pre-Merge Checklist

Run locally before every merge — there is no CI to catch failures:

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --workspace --all-targets` — zero warnings
- [ ] `cargo nextest run --workspace` — all tests pass
- [ ] `cargo check --workspace --all-targets` — all targets type-check
- [ ] Tests written first (TDD evidence in commit history)
- [ ] No `.unwrap()`, `.expect()`, `panic!()`, `todo!()`
- [ ] Error paths have `.context()`
- [ ] New public items have doc comments
- [ ] Minimum visibility applied
