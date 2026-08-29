# Handoff: running the `exp_repro` GPU investigation on an A100

**Purpose**: this investigation (`exp_repro_gpu_perf_investigation.md`, same directory)
was done entirely on a laptop RTX 4060 (cc8.9, no clock-lock, single GPU). The
conclusions — the block-size-1 launch bug, the ~2.5x recovery from `thread_limit(128)`
(E7), the `bind(...)` argument-order correctness bug (§8.2), and the occupancy-vs-noise
finding on inlining (§8.3) — are all real on that hardware, but "how much of this
matters on a real datacenter GPU" is an open question this laptop can't answer. These
are step-by-step instructions to reproduce the core comparisons on an A100.

**You do not need to re-derive anything.** The environment, branches, and tooling
already exist and are proven — this is a transplant, not a fresh investigation.

---

## 0. What's already done for you

- **Container image**: `~/netcdf-nvhpc26.5-el8.sif` (1.4 GB, built 2026-08-09) — a
  Ubuntu 22.04-based Apptainer/Singularity image with nvhpc 26.5 (CUDA 13.2), OpenMPI
  4.1.6 built against it with CUDA + Slurm (PMIx 4.1.2) support, and netcdf-fortran
  4.6.2. Recipe: `netcdf-nvhpc26.5.def` in this directory (`src/MOM6/`), in case it
  ever needs rebuilding (`sudo singularity build <name>.sif netcdf-nvhpc26.5.def` — this
  needs root/fakeroot, which most cluster login nodes don't grant, so **prefer copying
  the already-built `.sif` over rebuilding it on the cluster**).
- **All `exp-repro/*` branches are pushed** to `origin` on the `src/MOM6` submodule,
  which is your own fork (`git@github.com:edoyango/MOM6.git`) — not the shared
  `mom6-gpu` upstream. Nothing further needs pushing; `git clone`/`git fetch` on the
  A100 side will see all of them.
- **`mom6test`/`mom6run`** (`~/.local/bin/`): the wrapper scripts used for every
  build+run+compare cycle in this investigation. Copy both over — they're plain bash,
  no laptop-specific paths baked in beyond defaults you pass with `--dir`/`--container`.
  `mom6run` already has a `--container IMAGE` flag built for exactly this situation
  (wraps `make`/`mpiexec` in `singularity exec --nv --bind "$HOME" IMAGE`).

## 1. Sanity-check the target machine before anything else

- **Confirm you have full, exclusive access to an A100** — not a MIG slice. `nvidia-smi
  -L` should list one whole A100; if MIG is enabled you'll see instance UUIDs instead,
  and timing comparisons against this investigation's numbers (measured on one whole,
  if small, GPU) won't be apples-to-apples.
- **Confirm Apptainer/Singularity is available** on the node that will actually run the
  job (`apptainer --version` or `singularity --version`), and that it's launched with
  `--nv` so the container sees the driver.
- **Build on a node where the GPU is actually visible**, not a CPU-only login node.
  nvfortran auto-detects the target compute capability from the visible driver at
  compile time (that's why `config.mk.sample`'s `-gpu=cc90` override is commented out
  by default) — building without a GPU present, or on a different node than you'll run
  on, risks compiling for the wrong architecture silently. If your cluster's login and
  compute nodes have different GPUs (or none on login), either build inside a GPU-having
  interactive/batch allocation, or force it explicitly with `-gpu=cc80` (A100's compute
  capability) in `config.mk` — see step 4.
- If this cluster uses Slurm (the container was built against a specific host PMIx/Slurm
  version — see `netcdf-nvhpc26.5.def`'s `%post` comments), direct-launch (`srun`) should
  work; if PMIx versions don't line up, fall back to launching `mpiexec` from inside an
  interactive allocation instead of via `srun` directly.

## 2. Transfer what git won't carry

Everything in `src/MOM6` travels via git (branches are already pushed — see §0). Three
things are gitignored and need an explicit copy (`scp`/`rsync`) from this machine:

```bash
# the container image
~/netcdf-nvhpc26.5-el8.sif

# the correctness references this whole investigation gates on
~/MOM6-examples/ocean_only/benchmark/ocean.stats.exprepro
~/MOM6-examples/ocean_only/benchmark/ocean.stats.base

# the wrapper scripts
~/.local/bin/mom6test
~/.local/bin/mom6run
```

`ocean.stats.exprepro` is the one that matters most: every experiment's correctness gate
is `sha256sum ocean.stats` == `sha256sum ocean.stats.exprepro`. **Do not regenerate it on
the A100** — copy the exact file from this machine. (Whether a *fresh* A100-generated
`exp_repro` reference would also match this file bit-for-bit is actually an interesting
question in its own right — `exp_repro` is built entirely from portable IEEE-754
arithmetic and a lookup table, no hardware transcendental unit, specifically so it
*should* be architecture-independent — but that's a side experiment, not something to
rely on for the correctness gate here. Keep the laptop's reference as ground truth.)

## 3. Get the repo onto the A100 machine

```bash
git clone --recursive git@github.com:edoyango/MOM6-examples.git
cd MOM6-examples/src/MOM6
git fetch origin 'refs/heads/exp-repro/*:refs/heads/exp-repro/*'
git branch --list 'exp-repro/*'   # should list all 12: baseline, b0, e2-e11
```

(If `MOM6-examples` is already present via a different remote/method, just make sure
`src/MOM6`'s `origin` points at `edoyango/MOM6.git` and fetch from there — that's where
all the `exp-repro/*` branches live, not the `mom6-gpu` upstream.)

## 4. Configure and build, inside the container

From `ocean_only/`, create `config.mk` (gitignored on this machine too — it was never
actually used here; the equivalent flags were hand-edited into the generated
`build/Makefile` instead, but `config.mk` is the cleaner, documented path for a fresh
setup and matches those flags exactly):

```bash
cd MOM6-examples/ocean_only
cp config.mk.sample config.mk
```

`config.mk.sample` already carries this investigation's exact flags
(`-O4 -mp=gpu -stdpar=gpu -gpu=mem:separate -Mnovect -Mnofma -Minfo=mp,accel,inline`,
plus the `-Minline=name:flux_elem,...` MOM6-side requirement). Leave the
`#FCFLAGS += -gpu=cc90` line commented out to let nvfortran auto-detect the A100 (cc80)
from the driver, unless step 1's node-visibility check says you need to force it —
in that case set it to `-gpu=cc80`, not `-cc90` (that override predates this handoff and
targets a different GPU).

Build FMS and MOM6 together, wrapped in the container, targeting a fresh build
directory name so it's obviously distinct from anything already on this cluster:

```bash
apptainer exec --nv --bind "$HOME" ~/netcdf-nvhpc26.5-el8.sif \
  make BUILD=build-a100 FMS_BUILD=../shared/fms/build-a100 -j6
```

This is the same two-variable pattern (`BUILD=`/`FMS_BUILD=`) already used for the
`build-container` directory on the laptop — bootstraps FMS via autoconf, then MOM6,
both against the container's own spack-built netcdf-fortran/OpenMPI.

## 5. What to actually run

Don't rerun the full 13-branch matrix blind — most of it (E2–E6) was diagnostic
call-mechanism work that turned out not to be the bottleneck (§5/§6 of the main
report). The comparisons that actually answer "how much does this matter on a real
GPU" are:

| Priority | Branch | Why |
|---|---|---|
| 1 | `exp-repro/b0-intrinsic-exp` | speed floor: no `exp_repro` overhead at all |
| 1 | `exp-repro/baseline` | the original bug: block-size-1 launch, `exp_repro` as cherry-picked |
| 1 | `exp-repro/e7-thread-limit-only` | **the recommended, shipped fix** |
| 2 | `exp-repro/e8-bind-parallel-teams` | alternative fix (now `bind(teams,parallel)`, §8.2) |
| 2 | `exp-repro/e9-inline-plus-bind` | fully inlined + safe `bind` — §8.3 asks whether inlining's occupancy gain (real on the 4060, but swamped by wall-clock noise there) shows up more clearly on an A100 with far more SMs and a much steadier clock |

For each branch: `git checkout exp-repro/<name>` in `src/MOM6`, rebuild (§4's `make`
line — only touched files recompile, not a full rebuild each time), then:

```bash
cd ocean_only/benchmark
MOM6TEST_BENCH_DIR=$PWD apptainer exec --nv --bind "$HOME" \
  ~/netcdf-nvhpc26.5-el8.sif \
  mom6run --dir ../build-a100 --base ocean.stats.exprepro
```

(or adapt `mom6test`'s own env-loading preamble — on the laptop it sources spack +
`load_netcdf`; inside this container none of that is needed since the container's
`%environment` already puts nvhpc/spack's netcdf-fortran on `PATH`/`LD_LIBRARY_PATH`,
so `mom6run` alone, run inside `apptainer exec --nv`, is suffient — no spack sourcing
step to port over.)

Take the median of at least 3 runs per branch as before. Read
`benchmark/logfile.000000.out` for `(Ocean set_viscous_ML)` and
`(Ocean vertical viscosity)` — the same two clocks §3 of the main report tracks.
Confirm every run prints `PASS: ocean.stats matches ocean.stats.exprepro` before
trusting its timing.

## 6. Optional: kernel-level profiling (`ncu`)

Only worth doing if you want the §8.3-style occupancy comparison, not just wall-clock
numbers. Two cluster-specific gotchas to check for *before* assuming `ncu` just works:

- **Profiling counter permissions**: most clusters lock these down by default
  (`ERR_NVGPUCTRPERM` / "permission issue"). This needs either an admin-enabled node
  (`NVreg_RestrictProfilingToAdminUsers=0`) or a profiling-specific queue/reservation —
  check your cluster's docs rather than assuming `ncu` "just works" the way it does on
  a personal laptop where you own the driver.
- **MIG, again**: `ncu`'s launch-statistics/occupancy metrics are meaningless on a MIG
  slice for the same reason as §1 — confirm a whole-GPU allocation first.

If both check out, the invocation pattern is identical to what this investigation
already used (see §6/§8.3 of the main report):
`ncu --set full -k "regex:<kernel-name-fragment>" --launch-skip N --launch-count M -o
report mpiexec -n 1 build-a100/MOM6`, then `ncu --import report.ncu-rep --page details
--print-summary per-kernel`. Kernel names embed source line numbers and shift between
branches (inlining adds lines) — re-derive them per branch with `cuobjdump -res-usage`
first (`cuobjdump` ships inside the nvhpc install at
`compilers/../cuda/13.2/bin/cuobjdump` inside the container), don't assume the laptop's
`F1L23xx`/`F1L28xx` line numbers carry over.

## 7. What to bring back

For a clean comparison against this investigation's own tables (§3 of
`exp_repro_gpu_perf_investigation.md`), report at minimum:

- Median `(Ocean set_viscous_ML)` and `(Ocean vertical viscosity)` for B0, baseline,
  and E7 — this is the core "how much did the fix recover, and how big was the bug in
  the first place" comparison, on real hardware.
- Whether E8/E9 still land within noise of E7 and of each other (as they do on the
  4060), or whether the A100's much larger SM count and steadier clock actually
  separates them — particularly whether E9's real occupancy gain (§8.3: +2.5–3.7pp
  compute throughput, ~5% shorter kernel duration on the 4060) becomes visible at the
  routine-grain wall clock here where it didn't on the laptop.
- Correctness: confirm every configuration still gates `PASS` against the copied
  `ocean.stats.exprepro` — a real cross-architecture mismatch here would itself be a
  significant, separate finding about `exp_repro`'s portability claim (see the note in
  §2 above).
