# Progress on `convert_to_submodules.py`

This file tracks the fixes made while iterating on the MOM6 module-to-submodule conversion experiment.

## Overall Goal

`convert_to_submodules.py` converts Fortran module procedure bodies into companion `_s.F90` submodules, leaving the parent module with separate module procedure interfaces. The work so far has focused on making the converter robust enough for MOM6 source patterns and Intel Fortran compilation.

## Major Improvements Made

- Added a `Conversion` record so each successful conversion carries the parent path, submodule path, module name, and submodule name.
- Made the existing `--suffix` option functional instead of hardcoding `_s`.
- Added `--deps-file`, defaulting to `submodule_deps.mk`.
- Added generation of make dependency rules so submodule objects depend on the parent object and the parent `.mod`/`.smod` sidecars. This addresses parallel build races and stale sidecar issues.
- Added `SKIP_FILES` and excluded `MOM_error_handler.F90` because converting it causes a compiler ICE.

## Parser Fixes

- Fixed module-level `contains` detection so type-bound `contains` blocks and interface blocks are not mistaken for the module-level `contains`.
- Added stale submodule cleanup when a source is skipped or has no module-level procedures, avoiding orphan `_s.F90` files.
- Added support for `impure` procedure prefixes, which fixed missed procedures such as `check_capacity_by_id`.
- Improved procedure declaration scanning around preprocessor conditionals:
  - Declaration-only `#if/#ifdef/#ifndef` blocks can remain in the interface.
  - Executable preprocessor blocks stay in the submodule body, preventing split `#ifdef/#else/#endif` errors.
- Fixed `find_procedure_decl_end` to consume continuation lines belonging to a declaration: previously the scanner checked `is_declaration()` on each line independently. A continuation line like `gridLatT => NULL(), & ! comment` (following `real, pointer, dimension(:) :: &`) has no `::` or type keyword, so `is_declaration` returned False and the loop broke — leaving `last_decl` pointing at the dangling `& ` line and the continuation lines mis-classified as body. The fix advances `i` and updates `last_decl` for every continuation line that follows a recognised declaration.
- Fixed function return variable type not visible in `module procedure` body: ifx does not propagate array return-variable dimensions from the parent module interface into a separate module subprogram body (scalar types propagate fine). Fixed in `build_submodule_procedure` by scanning `iface_decls` for the result variable's declaration and prepending it to `local_decls` only when the return type is an array (`dimension` or `allocatable` in the type spec, or array bounds on the variable name).
- Fixed `submodule_preprocessor_lines` to preserve `#ifdef/#else/#endif` structure around `#define` blocks: previously it extracted only `#define` lines, emitting both branches of a conditional unconditionally (causing macro-redefinition warnings). Now it collects entire conditional blocks and emits the whole block when it contains `#define` lines.
- Fixed `ends_with_continuation` to recognise `& ! comment` patterns: previously `ends_with_continuation` only matched lines ending literally with `&` (after `rstrip()`), so lines like `gridLatT => NULL(), & ! label comment` were not recognised as continuations. This caused multi-line pointer declarations to be broken into individual single-line groups, leaving the first line (`real, pointer, dimension(:) :: &`) as an incomplete fragment in the interface block. The fix uses `re.search(r'&\s*(?:!.*)?$', s)` to match `&` optionally followed by whitespace and a comment.
- Fixed mixed dummy-arg/local declarations on the same line (e.g. `integer :: func_name, local_var` where `func_name` is the implicit function return variable). Previously the whole line went to the interface, leaving the local undeclared in the submodule body. Now such lines are split into two synthetic declaration lines — one in the interface (dummy names only) and one in the submodule body (local names only).
- Fixed top-level comma splitting in `split_decls`: replaced `after.split(',')` with `split_top_level(after)` so that commas inside dimension expressions (e.g. `weights(widths(layer+1),widths(layer))`) are not mistaken for variable-name separators. Previously, `widths` from the dimension spec was extracted as a "declared name", matched a dummy arg, and the whole local declaration was misclassified as an interface declaration.
- Grouped continued declaration lines before classifying them as interface or local declarations. This fixed split declarations like:
  ```fortran
  integer(kind=int64), &
           dimension(:,:), intent(inout) :: field
  ```
- Improved handling of old-style declarations without `::`, such as `integer rc`, so locals stay in the submodule body unless they declare dummy arguments.
- Duplicated procedure-local `use` statements into both the parent interface and submodule body when needed. This fixed cases where local declarations used imported kind parameters such as `int64`.

## Preprocessor and Include Handling

- Treated `#include "version_variable.h"` as a local implementation include so it stays with the procedure body that uses `version`.
- Added submodule propagation for selected module-level macro includes, currently `MOM_memory.h`, so submodule bodies can use macros such as `ALLOC_` and `DEALLOC_`.
- Added propagation of module-level `#define` lines into submodules, fixing macros such as `DIAG_ALLOC_CHUNK_SIZE` used by moved procedure bodies.

## Current Generated-Source Repairs Applied During Debugging

These were direct repairs to already-generated files after fixing the converter logic:

- `src/framework/MOM_intrinsic_functions.F90` / `_s.F90`: kept executable `#ifdef __INTEL_COMPILER` logic in the submodule body.
- `src/equation_of_state/MOM_EOS_base_type.F90` / `_s.F90`: reconverted after fixing module-level `contains` detection.
- `src/framework/MOM_file_parser.F90` / `_s.F90`: moved `version_variable.h` into the submodule body using `version`.
- `src/core/MOM_grid_s.F90`: added `#include <MOM_memory.h>` for allocation macros.
- `src/framework/MOM_io.F90` / `_s.F90`: restored `integer :: rc` as a submodule-local declaration.
- `config_src/infra/FMS2/MOM_coms_infra.F90` / `_s.F90`: repaired a split continued declaration.
- `src/framework/MOM_random_s.F90`: added `use iso_fortran_env, only : int64` to the submodule procedure using `int64`.
- `src/framework/MOM_diag_buffers.F90` / `_s.F90`: added the missing `impure module function check_capacity_by_id` interface and body.
- `src/ice_shelf/MOM_ice_shelf_diag_mediator_s.F90`: added `#define DIAG_ALLOC_CHUNK_SIZE 15` for submodule macro use.

## Verification Done

- Repeatedly ran:
  ```bash
  python3 -m py_compile convert_to_submodules.py
  ```
  after script edits.
- Did not run the full MOM6 build or test suite. Compile issues were addressed incrementally from the user-provided compiler diagnostics.

## makedep Bug Fix

- Fixed a bug in `ac/makedep` where the post-scan loop (lines 197-207) processed both `.smod`
  and `.mod` entries from `parent2subobjs`. The loop uses `smod[:-5] + '.mod'` to find the
  parent `.mod` file, assuming a 5-character `.smod` suffix. But `submodule_parent_files`
  also puts `.mod` entries (4-char suffix) into `parent2subobjs`, causing `smod[:-5]` to
  strip one character from the module name (e.g., `mom_dynamics_split_rk2b.mod` → stripped to
  `mom_dynamics_split_rk2` → parent looked up as `mom_dynamics_split_rk2.mod`). This
  overwrote the correct `mod2o` entry, making the submodule depend on the WRONG parent object
  and excluding the correct parent from the link target. Fix: added `if not smod.endswith('.smod'): continue`.

## Remaining Cautions

- The converter is still heuristic, not a full Fortran parser.
- It has been hardened against MOM6 patterns seen so far, but more compile errors may expose additional source patterns.
- After conversion, stale `.o`, `.mod`, and `.smod` files should be cleaned before rebuilding, especially after parent interface changes.
- `submodule_deps.mk` must be included or equivalent dependencies must be represented in the build system for parallel builds.
