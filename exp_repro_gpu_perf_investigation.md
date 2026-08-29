# `exp_repro` performance investigation — full results

**Branch:** `port-set_viscous_ML-tile` (all experiment branches below fork from
`exp-repro/baseline`, which is this branch plus the cherry-picked `exp_repro` /
`set_viscous_ML` integration, committed at `a18abae3e`)
**Repo:** `src/MOM6` submodule of `MOM6-examples` (all branches live here, not in the
outer repo)
**Test case:** `ocean_only/benchmark`, 360×180×22, 96 steps (`DAYMAX=1.0`), single GPU
**Hardware/toolchain:** NVIDIA RTX 4060 Laptop GPU (cc8.9, 8GB), nvhpc 26.5, CUDA 13.2
**Correctness gate:** `mom6test --base ocean.stats.exprepro` — sha256 byte-exact match
required on every row unless noted otherwise
**Timing:** median of 3 `mom6test` runs unless noted; read from
`benchmark/logfile.000000.out`, clocks `(Ocean set_viscous_ML)` (routine-grain,
isolates this one subroutine) and `(Ocean vertical viscosity)` (module-grain, includes
`set_viscous_ML` plus `vertvisc_coef`/`vertvisc`/`vertvisc_remnant`/`thickness_to_dz`)
**Static fingerprint:** `cuobjdump -elf` (relocation count for a name → device-call
count) and `cuobjdump -res-usage` (REG/STACK per kernel), on
`ocean_only/build/MOM6`

---

## 1. Starting point

The cherry-picked range `45f6b30..b63251d` added `exp_repro()`, a bitwise-reproducible
`exp()` implemented as a Fortran **submodule** procedure (`MOM_exp.F90`,
`submodule (MOM_intrinsic_functions) MOM_exp`). The (then-uncommitted, now baseline)
change swapped the two `exp()` calls in `set_viscous_ML` for `exp_repro()`
(`MOM_set_viscosity.F90:2430`, `:2783`) and added `!$omp declare target` to `exp_repro`
and its two lookup tables. It reproduces the CPU binary bit-for-bit, but is slow:

| Clock | Baseline (`exp_repro`) |
|---|---|
| `(Ocean set_viscous_ML)` | 0.690 s |
| `(Ocean vertical viscosity)` | 2.175 s |
| `(Ocean set BBL viscosity)` | 0.103 s |
| Total runtime | 32.85 s |

### Static diagnosis of the baseline binary

Before any experiment, `cuobjdump` on the as-cherry-picked binary showed `exp_repro` is
a genuine cross-module RDC device call, not inlined:

- PTX contains `.extern .func (.param .b64 func_retval0) mom_intrinsic_functions_exp_repro_`
  — one relocated `CALL.ABS.NOINC` per kernel (`R_CUDA_ABS47_34`).
- The argument is passed **by reference**; the callee reads the double back **byte by
  byte** (`ld.u8 [%rd1+7]`, `ld.u8 [%rd1+6]`, …) because `transfer(x, int_mold)` on a
  by-reference dummy lowers to individual byte loads.
- The two `exp_repro`-calling kernels are the *only* `set_viscous_ML` kernels with a
  stack frame:

  | Kernel | REG | STACK |
  |---|---|---|
  | u-half, contains `exp_repro` (`F1L2395_12_`) | 96 | 8 |
  | v-half, contains `exp_repro` (`F1L2747_16_`) | 106 | 8 |
  | u-half, no `exp_repro` (`F1L2347_10_`) | 60 | 0 |
  | v-half, no `exp_repro` (`F1L2697_14_`) | 64 | 0 |
  | `exp_repro` itself (standalone `.func`) | 38 | — |

Two candidate causes were ruled out before experimenting: `ENABLE_FAST_RINT` is already
on (`fast_rint`, not `ieee_rint`, is in use), and RDC was already the default before
this change (nothing was newly forced into relocatable-device-code linking).

This diagnosis — an out-of-line, by-reference cross-file call — turned out to explain
only a small fraction of the actual slowdown (see §7). It was a real, measurable effect,
just not the dominant one.

---

## 2. Experiment index

| # | Branch | What changed | `exp_repro` itself modified? |
|---|---|---|---|
| B1 | `exp-repro/baseline` | (nothing — this *is* the starting point) | — |
| B0 | `exp-repro/b0-intrinsic-exp` | both call sites reverted to intrinsic `exp(` | n/a |
| E2 | `exp-repro/e2-value-attr` | `VALUE` attribute on `exp_repro`'s `x` dummy | yes (1 line) |
| E3 | `exp-repro/e3-de-submodule` | move `exp_repro` out of the submodule into an ordinary module procedure | yes (restructure) |
| E4 | `exp-repro/e4-mextract-minline` | hoist `exp_remez_c` to module scope + `declare target` (prerequisite fixes); `-Mextract`/`-Minline` itself tried un-committed | yes (2 touch-ups) |
| E5a | `exp-repro/e5a-manual-inline-fn` | private same-file `exp_repro_inline()`, manually flattened | no (new function, `exp_repro` untouched) |
| E6 | `exp-repro/e6-inline-plus-threadlimit` | E5a + E7 combined (cherry-picked) | no |
| **E7** | `exp-repro/e7-thread-limit-only` | `teams thread_limit(128)` on the two target regions | **no** |
| **E8** | `exp-repro/e8-bind-parallel-teams` | `target teams` + `bind(parallel,teams)` on the loop clauses, no `thread_limit` | **no** |
| E9 | `exp-repro/e9-inline-plus-bind` | E5a + E8 combined (cherry-picked) | no |
| E10 | `exp-repro/e10-value-plus-bind` | E2 + E8 combined (cherry-picked) | yes (E2's 1 line) |

All branches fork directly from `exp-repro/baseline`; E6, E9, and E10 cherry-pick two
prerequisite branches' commits on top of a fresh checkout of baseline.

---

## 3. Results table

| # | ocean.stats | `exp_repro` relocs | REG u/v | STACK | `set_viscous_ML` (s) | `vertical viscosity` (s) |
|---|---|---|---|---|---|---|
| B0 (intrinsic `exp`) | ≈`.base`, roundoff-close but not byte-identical (different nvfortran version than `.base` was built with) — not independently fingerprinted | n/a | — | — | **0.271** | 1.759 |
| B1 (baseline) | =`.exprepro` | 2 | 96/106 | 8 | 0.690 | 2.175 |
| E2 (`VALUE`) | =`.exprepro` | 2 (unchanged) | 90/104 | **0** | 0.687 | 2.173 |
| E3 (de-submodule) | =`.exprepro` | 2 (unchanged) | 96/106 | 8 | 0.689 | 2.177 |
| E4 (`-Mextract`/`-Minline`) | — | — | — | — | **blocked** — see §5 | |
| E5a (manual inline) | =`.exprepro` | **0** | 90/80 | **0** | 0.636 | 2.123 |
| E6 (E5a+E7) | =`.exprepro` | 0 | — | — | 0.271 | 1.757 |
| **E7 (`thread_limit(128)`)** | =`.exprepro` | 2 (unchanged, still out-of-line) | 114/130 | 8 | **0.275** | **1.766** |
| **E8 (`bind(parallel,teams)`)** | =`.exprepro`, re-confirmed 3/3 runs | 2 (expected) | 90/96 | 8, SHARED:0 | **0.269** | **1.755** |
| E9 (E5a+E8) | **FAIL — non-deterministic, 3/3 runs disagree with each other and with `.exprepro`** | 0 | 90/80 | 0, but **SHARED:20** | 0.862 (3× slower) | — |
| E10 (E2+E8) | =`.exprepro`, 3/3 runs | 2 (expected) | 88/90 | 0, SHARED:0 | 0.270 | 1.758 |

All correctness checks: **PASS**, 3/3 runs, except E6 (1 confirmation run after the
proper cherry-pick, matching 3 earlier runs of the identical content: 0.271, 1.763; a
later 3/3-run re-check for §8 also passed) and **E9, which fails outright — see §8**.

---

## 4. Per-experiment detail

### B0 — intrinsic `exp()` (speed floor)
Revert both call sites to `exp(`, leave the `exp_repro` import in place (minimal diff).
Purpose: establish what `set_viscous_ML` costs with *no* reproducible-exp overhead at
all, as the target to match. 3 runs: 0.272274, 0.269621, 0.270528 (median 0.271);
vertical viscosity 1.762819, 1.755322, 1.758848 (median 1.759); total runtime
32.677–32.552 s. Its `ocean.stats` is roundoff-close to but not byte-identical with the
repo's `ocean.stats.base` (differs from step 24 onward around the 11th significant
digit) — expected, since `.base` predates this nvhpc version and native `exp()` isn't
bit-reproducible across compiler versions (the entire reason `exp_repro` exists).

### E2 — `VALUE` attribute
`real, intent(in) :: x` → `real, intent(in), value :: x` on `exp_repro`'s interface —
the *only* line changed (the submodule body inherits dummy attributes from the
interface). Removes the by-reference marshalling and byte-wise reload identified in the
static diagnosis. Result: STACK frame disappears (8→0), REG drops slightly (96→90,
106→104), but wall time is essentially unchanged (0.690→0.687s, −0.4%). **Conclusion:
marshalling cost was real but negligible; the call itself, not its argument passing, is
where the time was hypothesized to go.**

### E3 — de-submodule
Move `exp_repro`, `exp_remez_expm1_estrin_4`, and `fast_rint` out of the submodule
`MOM_exp.F90` into `MOM_intrinsic_functions.F90` as ordinary module procedures; delete
the submodule file. Required renaming `signbit`/`expbit` (already used by
`cuberoot`'s `rescale_cbrt`/`descale`) to `signbit_exp`/`expbit_exp` to avoid a
duplicate-declaration collision now that both live in one module scope. Purpose: test
whether `-Mextract`/`-Minline` needs this restructuring (it does not — see Stage 0 in
§5). On its own, performance-neutral as predicted: 0.690→0.689s, REG/STACK unchanged.

### E4 — `-Mextract`/`-Minline` cross-file — **blocked**
See §5 for the full investigation. Two source touch-ups were committed (hoisting
`exp_remez_expm1_estrin_4`'s coefficient array `c(0:4)` to module scope as
`exp_remez_c`, and adding `!$omp declare target` to that function) because Stage 0
proved they are necessary regardless of submodule structure. But applying the actual
`-Mextract`/`-Minline` flags to the real `MOM_set_viscosity.F90` compile — not
committed, since these are build flags, not source — fails with a reproducible
nvfortran-internal compiler error (`NVFORTRAN-S-1066`) specific to the v-half kernel's
internal sub-compile. This did not reproduce in an isolated microbenchmark with the
same flags, same extracted library, and equivalent kernel structure minus the real
file's nested-loop complexity.

### E5a — manual same-file inline
A private `elemental function exp_repro_inline(x)` added to `MOM_set_viscosity.F90`'s
`contains` section: the full range-reduction + Remez-polynomial arithmetic flattened
into one function (no nested calls), reusing `MOM_exp_data_n128`'s lookup tables
directly (only ~30 lines of arithmetic duplicated, no data duplication). Both call
sites switched to `exp_repro_inline(...)`. This is the mechanism expected to be the
**guaranteed upper bound** for a call-overhead fix — and it worked exactly as
advertised at the static level: `exp_repro` relocation count → 0, STACK → 0, REG down
to 90/80. But wall time only dropped 0.690→0.636s (**13% of the 0.690→0.271s gap to
B0**). **This was the finding that triggered re-profiling instead of continuing down
the call-mechanism path** (E1 `-gpu=lto` was never run, given what followed).

### E6 — E5a + E7 combined
Cherry-picked E5a's commit then E7's commit onto a fresh branch off baseline. No
further improvement over E7 alone (0.271s vs E7's 0.275s, both within run-to-run
noise) — direct confirmation that once the launch configuration is fixed, further
optimizing the call itself buys nothing measurable.

### E7 — `teams thread_limit(128)` — **the recommended fix**
Found via `ncu` profiling of E5a (see §6), not hypothesized. 4-line diff: the two
`!$omp target` lines that enclose `exp_repro`'s call sites become
`!$omp target teams thread_limit(128)`, and their matching `!$omp end target` become
`!$omp end target teams`. `exp_repro` itself is completely untouched — still the
original submodule, still called out-of-line (relocation count still 2, REG actually
*higher* than baseline at 114/130, because there is now real per-warp parallelism to
allocate registers across instead of none). Result: 0.690→0.275s median, matching B0
(0.271s) to within noise. This is the pattern `set_viscous_BBL` already uses correctly
elsewhere in the same file (`!$omp target teams loop collapse(2) thread_limit(128) ...`).

### E8 — `bind(parallel,teams)` (alternative to E7)
Same root fix, different mechanism: keep `target` → `target teams` (required —
`bind` is illegal without an enclosing `teams` region, confirmed by
`NVFORTRAN-S-1195` when tried on a bare `target`), but instead of capping team size
with `thread_limit(128)`, add `bind(parallel,teams)` to all four `!$omp loop` clauses
in the two `exp_repro`-calling regions. Result: matches E7's fix (0.690→0.269s median,
slightly *faster* than E7's 0.275s) with measurably lower register pressure (90/96 vs
E7's 114/130) — giving the compiler real information about the loop's binding lets it
allocate more precisely than an artificial thread cap does. See §7 for the tradeoff
that makes E7 the recommended choice anyway.

---

## 5. The E4 dead end, in detail

**Question tested**: does `-Mextract`/`-Minline` cross-file inlining require `exp_repro`
to first be moved out of its submodule (E3)? **Answer: no**, confirmed empirically, not
assumed — but the actual blocker turned out to be unrelated to submodules entirely.

Stage 0 (scratchpad, before touching any branch): a standalone nvfortran program
reproduced the production symptom exactly (4 relocations to `exp_repro`, `REG:38
STACK:0` on the standalone device function, matching the real binary).

1. `-Mextract=name:exp_repro,lib:exp.il -c MOM_exp.F90` (the **unmodified submodule
   file**) succeeded and produced a real extract entry.
2. `-Minline=lib:exp.il,name:exp_repro` on the driver reported `exp_repro
   cross-inlined, size=32` — **a submodule procedure inlines cross-file with zero
   restructuring.**
3. The real blocker: `exp_repro` calls a private helper, `exp_remez_expm1_estrin_4`,
   which itself failed to inline, producing first `NVFORTRAN-S-1073 ... must have
   'omp declare target' information`, then (once that was fixed) `subprogram not
   inlined -- static variable during crossing files: exp_remez_expm1_estrin_4`.
4. An isolated **plain-module** test (no submodule at all — `mymod2.F90`/`helper2`)
   reproduced the identical `static variable during crossing files` error for a
   function with a local `real, parameter :: c(...)` array, and reproduced *success*
   once that array was hoisted to module scope (`mymod3.F90`/`helper3`). **This is a
   general nvfortran `-Minline` restriction on functions with local static/parameter
   data, unrelated to submodules** — `fast_rint`, `exp_repro`'s *other* helper, has no
   such local array and inlined cross-file without issue, submodule and all.
5. With the array hoisted (`exp_remez_c`, module-scope) and `!$omp declare target`
   added directly to `exp_remez_expm1_estrin_4` (previously implicit only because its
   caller was in the same file — that implicit grant doesn't survive `exp_repro` being
   physically relocated into a foreign TU), the full scratchpad chain compiled,
   linked, and reproduced the exact same result as the unmodified baseline driver.
   `exp_repro`'s own relocation count dropped to 0; one smaller residual relocated
   call to `exp_remez_expm1_estrin_4` remained (it did not itself cross-inline even
   after the fix — a partial-inlining outcome, not full elimination).
6. Applying the same two fixes to the **real** `MOM_intrinsic_functions.F90` /
   `MOM_exp.F90` (committed on the E4 branch) and then running the identical
   `-Mextract`/`-Minline` two-pass against the real `MOM_set_viscosity.o`: compile
   fails with `NVFORTRAN-S-1066-The -cuda flag should be used with CUDA DEVICE
   variable - ..inline`, scoped specifically to the internal sub-compile of the
   v-half kernel (`__nv_mom_set_visc_set_viscous_ml__F1L2747_15`, "2 severes"). Bisected
   down to: `-Minline=lib:exp.il` alone (no `name:` filter, no `-stdpar=gpu`, minimal
   other flags) still triggers it on the real file, while the exact same flag against a
   trivial single-loop kernel (the Stage 0 driver) does not. The u-half kernel's
   sub-compile never even got reached (compilation aborted after the v-half error).

**Conclusion**: `-Mextract`/`-Minline` is a real, working mechanism for simple cases,
but hits a genuine nvfortran-internal limitation on `set_viscous_ML`'s specific kernel
structure (a sequential `do k` recurrence containing *two* separate `!$omp loop`
regions). Not something fixable from the MOM6 side within reasonable effort.

---

## 6. The pivotal finding: `ncu` Launch Statistics

E5a (manual inline, the mechanism expected to be the guaranteed upper bound for a
call-overhead fix) only recovered 13% of the B1→B0 gap
(`(0.690-0.636)/(0.690-0.271) = 12.9%`), despite achieving everything a call-overhead
fix should: zero relocations, zero stack frame, fewer registers than baseline. This
meant the out-of-line call was never the dominant cost.

`ncu --set basic -k "regex:set_viscous_ml__F1L2857" --launch-skip 5 --launch-count 2`
on E5a's binary (`Section: Launch Statistics`):

```
Block Size                                                     1
Grid Size                                                  65160
Registers Per Thread             register/thread              80
Waves Per SM                                              113.12
...
OPT   Est. Speedup: 96.88%
      Threads are executed in groups of 32 threads called warps. This kernel launch
      is configured to execute 1 threads per block. Consequently, some threads in a
      warp are masked off and those hardware resources are unused. ...
```

**Block Size = 1 thread per block.** Every single thread runs alone; 31/32 of every
warp is masked off. A separate `ncu` run on the *non*-`exp_repro` sibling region in the
same file (`F1L2457`, structurally identical `!$omp target` + `!$omp loop collapse(2)`,
just without the `exp_repro` call) showed a completely healthy launch:

```
Block Size                                                   128
Grid Size                                                    508
Duration                         us        46.18
```

46 µs vs. 2.95 ms for the `exp_repro` region — a 64× difference, from launch
configuration alone. Both regions are bare `!$omp target` with no `teams`/
`thread_limit` clause; nvfortran's own block-size heuristic apparently trips into a
degenerate block-size-1 fallback for the higher-register-pressure (`exp_repro`-calling)
region while choosing a sane 128 for the lower-pressure sibling. This — not `exp_repro`'s
call mechanism, marshalling, submodule housing, or even its arithmetic complexity — was
the actual root cause of the disproportionate slowdown.

---

## 7. E7 vs E8: `thread_limit(128)` vs `bind(parallel,teams)`

Both confirmed via `ncu` to produce the identical fix (Block Size 128, Grid ~510), by
different mechanisms:

- **`thread_limit(128)`** is a blunt cap on the `target teams` construct: it hardcodes
  a team size without telling the compiler anything about how the `loop` construct's
  iterations should map onto it.
- **`bind(parallel,teams)`** is a clause on the `loop` construct itself: it tells the
  compiler the loop binds across *both* the teams and parallel levels simultaneously,
  and the compiler derives its own block size (128 — the same value it already picks
  correctly for the non-`exp_repro` sibling region) from that information, rather than
  being told a fixed number.
- **`bind` requires `teams`, full stop**: tried on a bare `!$omp target` (no `teams`),
  nvfortran refuses to compile: `NVFORTRAN-S-1195 ... 'bind(teams)' can be used only
  when 'loop' region is strictly nested inside a 'teams' region.` So E8's diff still
  needed `target` → `target teams` — only the `thread_limit(128)` cap was replaced.
- **Registers, measurably lower** with `bind`: 90/96 (u/v) vs E7's 114/130 — real
  information about the loop's structure lets the compiler allocate more precisely
  than an artificial thread cap does.
- **Timing, marginally faster** with `bind`: 0.269s vs 0.275s median — plausibly
  downstream of the lower register pressure (more room for concurrent warps per SM).
- **Portability — the real tradeoff**: `bind(parallel,teams)` is **not standard
  OpenMP**. The OpenMP 5.x `bind` clause takes exactly one value (`teams`, `parallel`,
  or `thread`); nvfortran silently accepts the two-value form with **zero diagnostic**
  that it is nonstandard — no warning, no `-Minfo` note, nothing. Confirmed working on
  nvhpc 26.5; no guarantee it survives a future release or is portable to any other
  OpenMP-target compiler. `thread_limit(128)` is a standard clause already used
  correctly elsewhere in this exact file (`set_viscous_BBL`), with no such risk.

**Recommendation: ship E7.** The gap between E7 and E8 (0.275s vs 0.269s) is smaller
than the noise between repeated runs of either, and both already match the
intrinsic-`exp()` reference (0.271s) — nowhere near large enough to justify the
portability risk of an undocumented compiler-specific extension. E8 is kept on record
as a working, faster, measured alternative in case a future workload's register
pressure makes the difference matter more. §8 below found a second, sharper reason to
prefer E7.

---

## 8. Follow-up: does combining inlining with `bind(parallel,teams)` help further?

Asked directly: since E8 (`bind(parallel,teams)`, `exp_repro` still out-of-line) and E6
(full manual inline + `thread_limit(128)`) both independently match B0, does *stacking*
manual inlining on top of `bind(parallel,teams)` buy anything more? Two more branches
were built to answer this, both cherry-picked onto a fresh checkout of
`exp-repro/baseline`, and each re-run 3× to check not just a single `ocean.stats` match
but run-to-run determinism of the same binary:

- **E9** (`exp-repro/e9-inline-plus-bind`): E5a's `exp_repro_inline` cherry-picked, then
  E8's `bind(parallel,teams)` commit on top.
- **E10** (`exp-repro/e10-value-plus-bind`): E2's `VALUE` attribute (still an
  out-of-line call, just by-value) cherry-picked, then the same E8 commit on top.

**E10 passes cleanly**: `ocean.stats` matches `.exprepro` on all 3/3 runs, `exp_repro`
relocation count unchanged at 2 (still out-of-line), REG drops slightly further than E8
alone (88/90 vs E8's 90/96 — consistent with E2's by-value marshalling savings stacking
on top of `bind`'s savings), SHARED:0, timing indistinguishable from E8 alone
(0.270s/1.758s vs E8's 0.269s/1.755s). **Conclusion: partial call-mechanism changes
(VALUE attribute) combine safely with `bind(parallel,teams)`, with no measurable
additional benefit beyond `bind` alone.**

**E9 fails outright — not a rounding difference, a real, reproducible bug.** The very
first `mom6test` run already disagreed with `.exprepro` starting at step 24 (En differs
at the 4th significant digit — `0.546574` vs `0.546815` — and the mass/salt/heat
conservation errors `Me`/`Se`/`Te` are six orders of magnitude worse than reference,
`~1e-10` vs `~1e-16`, not the ~11th-digit roundoff seen in the benign B0/`.base`
mismatch in §4). Suspecting a fluke, the identical binary was re-run three more times
with no rebuild: **three different `sha256`s, three different `ocean.stats`, all
different from each other and from `.exprepro`** —

```
run1: cb4b4c63d22b4c7e532fb752ad6bda8a387248055e4dc74e5acdab986debe5d3  En(step24)=0.54657433...
run2: 9124763d47a122c5e728c96635d25d4e3f6f419193063a46cff4649d3c7f0292  En(step24)=0.54657397...
run3: bcd43d9c2d8cc77e6cb564ffdaaed38b0b553f695cceeb33c4d762e6f85a4379  En(step24)=0.54656911...
```

This is a genuine, **runtime** data race, not a compile-time miscompilation — a
miscompilation would produce the same wrong answer every run. `cuobjdump -res-usage`
on the E9 binary pinpointed the mechanism: **`SHARED:20` bytes on exactly the two
`exp_repro`-bearing kernels** (`F1L2505_12_` REG:90 SHARED:20, `F1L2857_16_` REG:80
SHARED:20) — every other kernel in the file, in every other experiment measured in this
whole investigation (B1, E2, E3, E5a, E6, E7, E8, E10), reports `SHARED:0`. Something in
the fully-inlined function body is being placed in team-shared memory once
`bind(parallel,teams)` is also present, and multiple threads within a team appear to
race on it — exactly reproducing the observed non-determinism. E9 is also **~3.2×
slower** than E6/E7/E8/E10 (0.862s vs ~0.27s), consistent with the extra shared-memory
traffic/serialization: this combination is a double loss, not a trade-off.

Two isolation checks pin the trigger down precisely:

- **`bind(parallel,teams)` alone (E8, E10) is solid**: out-of-line calls, whether
  by-reference (E8) or by-value (E10), never show `SHARED>0` and never show
  non-determinism (3/3 identical runs each).
- **Full inlining alone (E5a) or full inlining + `thread_limit` (E6) is solid**: E6 was
  re-run 3× for this follow-up specifically (not just the original single confirmation
  run) — `ocean.stats` matches `.exprepro` bit-for-bit on all three, `SHARED:0` on both
  exp-bearing kernels (REG:96/96, higher than E9's 90/80 — no shared-memory shortcut was
  taken).

So the bug requires **both** ingredients: `exp_repro`'s body must be physically inlined
into the loop (E5a/E9), *and* the loop must carry `bind(parallel,teams)` (E8/E9) rather
than `thread_limit` (E7/E6). Neither ingredient alone is sufficient.

**Answer to "does inlining improve on `bind(parallel,teams)`?": no — it breaks it.**
E8/E10 (bind, call left out-of-line) is the performance ceiling reachable through this
mechanism; pushing further via manual inlining does not yield a faster *correct*
configuration, it yields a slower, non-deterministic one. This is also a second,
independent argument (beyond §7's portability point) for preferring **E7
(`thread_limit(128)`) over E8** as the shipped fix: E7's safety is verified to survive
being combined with full inlining (E6), while E8's does not — a future change that
happens to get `exp_repro` inlined (a smarter compiler, a future manual-inlining PR)
would silently reintroduce this exact race if `bind(parallel,teams)` were the shipped
mechanism instead of `thread_limit`.

---

## 9. What each diff actually touches (verified, not assumed)

Confirmed by direct `git diff exp-repro/baseline exp-repro/e7-thread-limit-only --
.` and the equivalent for E8, across the **entire repo**, not just the one file:

- **E7**: `src/parameterizations/vertical/MOM_set_viscosity.F90`, 4 lines changed, 1
  file touched. `exp_repro`'s own implementation (`MOM_intrinsic_functions.F90`,
  `MOM_exp.F90`, `MOM_exp_data_n128.F90`) is byte-for-byte the pristine cherry-picked
  baseline. Both call sites (`MOM_set_viscosity.F90:2430`, `:2783`) are still literally
  `exp_repro(...)`.
- **E8**: same file, 8 lines changed, 1 file touched. Same guarantee — `exp_repro`
  completely untouched.

Neither fix required any of the E2/E3/E5 call-mechanism work — that entire line of
investigation, while thoroughly measured and worth keeping for the record, was solving
a problem that was not the actual bottleneck.

---

## 10. Branch/commit reference

All in `src/MOM6`, forked from `port-set_viscous_ML-tile` @ `6e63e6778`:

| Branch | Commit | Note |
|---|---|---|
| `exp-repro/baseline` | `a18abae3e` | shared parent for everything below |
| `exp-repro/b0-intrinsic-exp` | `83aa7ec59` | |
| `exp-repro/e2-value-attr` | `b662bfeab` | |
| `exp-repro/e3-de-submodule` | `d4d9c9fd2` | |
| `exp-repro/e4-mextract-minline` | `b089a6133` | source touch-ups only; flags not committed (see §5) |
| `exp-repro/e5a-manual-inline-fn` | `53fda9360` | |
| `exp-repro/e6-inline-plus-threadlimit` | `4b1b04863` (on top of `493bb28ae`) | cherry-picked E5a + E7 |
| `exp-repro/e7-thread-limit-only` | `6244a6edf` | **recommended fix** |
| `exp-repro/e8-bind-parallel-teams` | `3cd1f94c3` | measured alternative |
| `exp-repro/e9-inline-plus-bind` | `9790dc530` (on top of `53fda9360`) | cherry-picked E5a + E8 — **fails: non-deterministic, see §8** |
| `exp-repro/e10-value-plus-bind` | `6a774e7fc` (on top of `b662bfeab`) | cherry-picked E2 + E8 — passes, no gain over E8 alone |

`src/MOM6` is currently checked out to `exp-repro/e7-thread-limit-only`, built, and
verified. A durable note on the underlying `!$omp target` block-size-1 launch bug, and
on the inline+`bind(parallel,teams)` race found in §8, is saved to project memory
(`mom6-omp-target-block-size-1-launch-bug.md`) for future GPU-porting sessions on this
codebase.
