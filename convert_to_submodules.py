#!/usr/bin/env python3
"""Convert Fortran modules to module+submodule pairs for compile-speed testing.

For each .F90 file containing a module with a contains section:
  - Rewrites the module file with only the declaration section + module procedure interfaces
  - Writes a companion _s.F90 submodule file with all procedure bodies

Usage:
  python3 convert_to_submodules.py [--dry-run] [--verbose] dir1 [dir2 ...]

The script is conservative: if it cannot confidently parse a file it skips it and
prints a warning rather than producing broken output.
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Optional, NamedTuple


def split_top_level(s, sep=','):
    """Split s by sep, ignoring separators inside parentheses."""
    parts = []
    depth = 0
    current = []
    for ch in s:
        if ch == '(':
            depth += 1
            current.append(ch)
        elif ch == ')':
            depth -= 1
            current.append(ch)
        elif ch == sep and depth == 0:
            parts.append(''.join(current))
            current = []
        else:
            current.append(ch)
    if current:
        parts.append(''.join(current))
    return parts


class Conversion(NamedTuple):
    parent: Path
    submodule: Path
    module_name: str
    submodule_name: str

# ---------------------------------------------------------------------------
# Regex helpers (applied to .strip().lower() of each raw line)
# ---------------------------------------------------------------------------

def _re(pat):
    return re.compile(pat, re.IGNORECASE)

RE_MODULE      = _re(r'^module\s+(\w+)\s*$')
RE_SUBMODULE   = _re(r'^submodule\s*\(')
RE_END_MODULE  = _re(r'^end\s*module\b')
RE_CONTAINS    = _re(r'^contains\s*$')
# type definition: "type [,attrs] :: name" or old-style "type name" (not "type(...")
RE_TYPE_START  = _re(r'^type\b(?!\s*\().*::|^type\s+\w+\s*$')
RE_TYPE_END    = _re(r'^end\s*type\b')
RE_IFACE_START = _re(r'^interface\b')
RE_IFACE_END   = _re(r'^end\s*interface\b')
RE_SUB_START   = _re(r'^((?:(?:module|pure|impure|elemental|recursive)\s+)*)subroutine\s+(\w+)')
# function prefixes: procedure attributes + optional return type
_FPFX = (r'(?:module|pure|impure|elemental|recursive'
         r'|real|integer|logical|complex|character'
         r'|double\s+precision'
         r'|type\s*\([^)]*\)|class\s*\([^)]*\)'
         r'|character\s*\([^)]*\)'
         r')\s+')
RE_FUNC_START  = _re(rf'^(\s*(?:{_FPFX})*)function\s+(\w+)')
RE_SUB_END     = _re(r'^end\s*subroutine\b')
RE_FUNC_END    = _re(r'^end\s*function\b')
RE_INCLUDE     = _re(r'^#\s*include\s+["<]([^">]+)[">]')

LOCAL_INCLUDE_FILES = {
    'version_variable.h',
}

SUBMODULE_INCLUDE_FILES = {
    'MOM_memory.h',
}

SKIP_FILES = {
    'MOM_error_handler.F90',
}

# Declaration-section starters (checked against stripped+lowercased line)
_DECL_STARTS = (
    'use ', 'use,',
    'implicit ',
    'integer', 'real', 'double precision', 'complex',
    'logical', 'character',
    'type(', 'type ::', 'class(', 'class ::',
    'procedure(', 'procedure,', 'procedure ::',
    'external ', 'intrinsic ',
    'public', 'private', 'protected',
    'parameter', 'dimension',
    'include ',
    'namelist',
    'enum,', 'enumerator',
)

# ---------------------------------------------------------------------------
# Low-level line helpers
# ---------------------------------------------------------------------------

def stripped(line: str) -> str:
    return line.strip()

def slow(line: str) -> str:
    return line.strip().lower()

def is_blank_or_comment(line: str) -> bool:
    s = stripped(line)
    return not s or s.startswith('!')

def is_preprocessor(line: str) -> bool:
    return stripped(line).startswith('#')

def preprocessor_kind(line: str) -> Optional[str]:
    """Return the preprocessor directive keyword, or None."""
    m = re.match(r'^\s*#\s*(\w+)', line)
    return m.group(1).lower() if m else None

def find_matching_endif(lines: list[str], start: int) -> Optional[int]:
    """Return the matching #endif for a #if/#ifdef/#ifndef line."""
    depth = 0
    for i in range(start, len(lines)):
        kind = preprocessor_kind(lines[i])
        if kind in ('if', 'ifdef', 'ifndef'):
            depth += 1
        elif kind == 'endif':
            depth -= 1
            if depth == 0:
                return i
    return None

def preprocessor_block_is_declaration(lines: list[str], start: int,
                                      end_limit: int) -> bool:
    """Return true if a conditional preprocessing block only has declarations."""
    end = find_matching_endif(lines, start)
    if end is None or end >= end_limit:
        return False

    for line in lines[start + 1:end]:
        kind = preprocessor_kind(line)
        if kind in ('if', 'ifdef', 'ifndef', 'elif', 'else', 'endif'):
            continue
        if not is_declaration(line):
            return False
    return True

def is_declaration(line: str) -> bool:
    """Heuristic: is this line part of the specification (declaration) section?"""
    s = stripped(line)
    lo = s.lower()
    if not s or s.startswith('!') or s.startswith('#'):
        return True
    # OpenMP/OpenACC that appear in declarations (e.g. !$OMP THREADPRIVATE)
    if s.startswith('!$') and 'threadprivate' in lo:
        return True
    if s.startswith('!$'):
        return False
    for starter in _DECL_STARTS:
        if lo.startswith(starter):
            return True
    # Modern style :: (but not do concurrent which has ::)
    if '::' in s and not lo.startswith('do ') and not lo.startswith('do\t'):
        return True
    # end-of-block markers that appear in declaration sections
    if re.match(r'end\s+(type|interface|enum)\b', lo):
        return True
    return False

def ends_with_continuation(line: str) -> bool:
    s = line.rstrip()
    # Match & optionally followed by whitespace and a comment (e.g. "val, & ! comment")
    return bool(re.search(r'&\s*(?:!.*)?$', s))

# ---------------------------------------------------------------------------
# Module-level parsing
# ---------------------------------------------------------------------------

def find_module_name(lines: list[str]) -> Optional[str]:
    """Return the module name (original case) if this file defines a simple module, else None."""
    for line in lines:
        lo = slow(line)
        if RE_SUBMODULE.match(lo):
            return None   # already a submodule
        m = RE_MODULE.match(lo)
        if m:
            if m.group(1).lower() == 'procedure':
                return None  # 'module procedure' line, not a module definition
            # Recover original-case name from the raw line
            m2 = re.match(r'^\s*module\s+(\w+)\s*$', line, re.IGNORECASE)
            return m2.group(1) if m2 else m.group(1)
    return None

def find_module_contains(lines: list[str]) -> Optional[int]:
    """Return line index of the module-level 'contains', or None."""
    block_stack = []

    for i, line in enumerate(lines):
        lo = slow(line)
        if not lo or lo.startswith('!') or lo.startswith('#'):
            continue

        if RE_TYPE_START.match(lo):
            block_stack.append('type')
            continue
        if RE_IFACE_START.match(lo):
            block_stack.append('interface')
            continue
        if RE_TYPE_END.match(lo):
            if block_stack and block_stack[-1] == 'type':
                block_stack.pop()
            continue
        if RE_IFACE_END.match(lo):
            if block_stack and block_stack[-1] == 'interface':
                block_stack.pop()
            continue
        if not block_stack and RE_CONTAINS.match(lo):
            return i

    return None

def find_end_module(lines: list[str], start: int) -> Optional[int]:
    """Return index of 'end module' starting from 'start'."""
    for i in range(start, len(lines)):
        if RE_END_MODULE.match(slow(lines[i])):
            return i
    return None

def submodule_preprocessor_lines(module_decl_lines: list[str]) -> list[str]:
    """Return module-level preprocessing lines also needed by submodules."""
    result = []
    seen_includes = set()
    n = len(module_decl_lines)
    i = 0
    while i < n:
        line = module_decl_lines[i]
        sl = line.strip()

        m = RE_INCLUDE.match(sl)
        if m:
            include_name = Path(m.group(1)).name
            key = ('include', include_name)
            if include_name in SUBMODULE_INCLUDE_FILES and key not in seen_includes:
                result.append(line)
                seen_includes.add(key)
            i += 1
            continue

        # Preprocessor conditional block: collect the whole block and emit it
        # only if it contains #define lines, preserving the surrounding structure.
        if re.match(r'^#\s*(?:if|ifdef|ifndef)\b', sl):
            block = []
            depth = 0
            j = i
            while j < n:
                bl = module_decl_lines[j].strip()
                block.append(module_decl_lines[j])
                if re.match(r'^#\s*(?:if|ifdef|ifndef)\b', bl):
                    depth += 1
                elif re.match(r'^#\s*endif\b', bl):
                    depth -= 1
                    if depth == 0:
                        j += 1
                        break
                j += 1
            if any(re.match(r'^#\s*define\b', bl.strip()) for bl in block):
                result.extend(block)
            i = j
            continue

        # Standalone #define outside any conditional
        if re.match(r'^#\s*define\b', sl):
            result.append(line)

        i += 1
    return result

# ---------------------------------------------------------------------------
# Procedure-block extraction
# ---------------------------------------------------------------------------

class Procedure:
    def __init__(self, is_function: bool, name: str, prefixes: str,
                 all_lines: list[str], header_end: int, decl_end: int):
        self.is_function = is_function
        self.name = name
        self.prefixes = prefixes.strip()   # pure, elemental, recursive (NOT module)
        self.all_lines = all_lines         # complete procedure, header through end xxx
        self.header_end = header_end       # index into all_lines: last line of the header
        self.decl_end = decl_end           # index into all_lines: last declaration line

    @property
    def header_lines(self):
        return self.all_lines[:self.header_end + 1]

    @property
    def decl_lines(self):
        """All lines from after the header through end of declarations."""
        return self.all_lines[self.header_end + 1: self.decl_end + 1]

    def split_decls(self, dummy_args: list[str]):
        """
        Partition decl_lines into (iface_decls, local_decls).

        iface_decls  – lines that belong in the interface block: USE statements,
                       IMPLICIT, preprocessor directives, blank/comment lines,
                       and declarations of dummy arguments (recognised by the
                       presence of 'intent(' or a match against dummy_args).
        local_decls  – everything else: local variable declarations that must
                       stay in the module-procedure body.
        """
        dummy_set = {a.lower() for a in dummy_args}
        iface, local = [], []
        inside_iface = 0

        groups = []
        i = 0
        while i < len(self.decl_lines):
            group = [self.decl_lines[i]]
            while ends_with_continuation(group[-1]) and i + 1 < len(self.decl_lines):
                i += 1
                group.append(self.decl_lines[i])
            groups.append(group)
            i += 1

        for group in groups:
            line = ''.join(part.rstrip().rstrip('&') + ' ' for part in group)
            lo = slow(line)

            # Blank and comments preserve spacing in the interface.
            if not lo or lo.startswith('!'):
                iface.extend(group)
                continue

            # Some includes declare local data used by the implementation.
            if lo.startswith('#'):
                m = RE_INCLUDE.match(lo)
                if m and m.group(1) in LOCAL_INCLUDE_FILES:
                    local.extend(group)
                else:
                    iface.extend(group)
                continue

            # USE statements may support both dummy declarations in the
            # interface and local declarations in the submodule body.
            if lo.startswith('use ') or lo.startswith('use,'):
                iface.extend(group)
                local.extend(group)
                continue

            if lo.startswith('implicit '):
                iface.extend(group)
                continue

            # Internal interface blocks -> interface
            if RE_IFACE_START.match(lo):
                inside_iface += 1
            if RE_IFACE_END.match(lo) and inside_iface:
                inside_iface -= 1
            if inside_iface:
                iface.extend(group)
                continue

            if '::' in line:
                type_spec_raw, after_raw = line.split('::', 1)
                after = after_raw
                # Strip inline comment first
                if '!' in after:
                    after = after[:after.index('!')]
                # Parse (name, original_part) pairs, splitting only at top-level commas
                parts_and_names = []
                for part in split_top_level(after):
                    name = part.strip().split('(')[0].split('=')[0].strip().lower()
                    if name:
                        parts_and_names.append((name, part))
                names = {n for n, _ in parts_and_names}

                # Dummy-argument declaration: has 'intent(' OR names overlap dummy list
                if 'intent(' in lo or not names:
                    iface.extend(group)
                else:
                    dummy_parts = [p for n, p in parts_and_names if n in dummy_set]
                    local_parts = [p for n, p in parts_and_names if n not in dummy_set]
                    if not local_parts:
                        # All dummy args
                        iface.extend(group)
                    elif not dummy_parts:
                        # All local
                        local.extend(group)
                    else:
                        # Mixed (e.g. "integer :: func_name, local_var"): split into two lines
                        indent = group[0][:len(group[0]) - len(group[0].lstrip())]
                        type_spec = type_spec_raw.strip() + ' ::'
                        iface.append(indent + type_spec + ' ' + ', '.join(p.strip() for p in dummy_parts) + '\n')
                        local.append(indent + type_spec + ' ' + ', '.join(p.strip() for p in local_parts) + '\n')
            else:
                words = re.split(r'[\s,()]+', lo)
                decl_name = words[-1] if words else ''
                decl_names = {decl_name} if decl_name else set()

                # Old-style dummy declarations are rare but must remain in the interface.
                if 'intent(' in lo or (decl_names and decl_names & dummy_set):
                    iface.extend(group)
                elif is_declaration(line):
                    local.extend(group)
                else:
                    iface.extend(group)

        return iface, local

    @property
    def body_lines(self):
        """Lines from after declarations through (but not including) the end statement."""
        return self.all_lines[self.decl_end + 1: -1]

    @property
    def end_line(self):
        return self.all_lines[-1]

    def kind(self):
        return 'function' if self.is_function else 'subroutine'


def extract_dummy_args(header_lines: list[str]) -> list[str]:
    """Extract dummy argument names + return variable from the procedure header."""
    joined = ''
    for line in header_lines:
        s = line.rstrip()
        if s.endswith('&'):
            joined += s[:-1] + ' '
        else:
            joined += s
    # Dummy arguments from the argument list
    m = re.search(r'\(([^)]*)\)', joined)
    if not m:
        return []
    args_str = m.group(1)
    args = [a.strip().lower() for a in args_str.split(',') if a.strip()]
    # Explicit result variable via result(var)
    rm = re.search(r'\bresult\s*\(\s*(\w+)\s*\)', joined, re.IGNORECASE)
    if rm:
        args.append(rm.group(1).lower())
    else:
        # Implicit result variable: for a function with no result() clause the
        # return variable has the same name as the function itself.
        fm = re.search(r'\bfunction\s+(\w+)', joined, re.IGNORECASE)
        if fm:
            args.append(fm.group(1).lower())
    return args


def find_procedure_decl_end(all_lines: list[str], header_end: int,
                             dummy_args: list[str]) -> int:
    """
    Return the index of the last declaration line in the procedure.
    Scans from after the header; stops at the first executable line or at
    an internal 'contains'.
    """
    n = len(all_lines)
    last_decl = header_end
    inside_iface = 0

    i = header_end + 1
    while i < n - 1:  # -1 to exclude end subroutine/function
        lo = slow(all_lines[i])
        if not lo or lo.startswith('!'):
            last_decl = i
            i += 1
            continue
        if lo.startswith('#'):
            kind = preprocessor_kind(all_lines[i])
            if kind in ('if', 'ifdef', 'ifndef'):
                end = find_matching_endif(all_lines, i)
                if preprocessor_block_is_declaration(all_lines, i, n - 1):
                    last_decl = end
                    i = end + 1
                    continue
                break
            last_decl = i
            i += 1
            continue
        # Internal interface block is declarations
        if RE_IFACE_START.match(lo):
            inside_iface += 1
            last_decl = i
            i += 1
            continue
        if RE_IFACE_END.match(lo) and inside_iface:
            inside_iface -= 1
            last_decl = i
            i += 1
            continue
        if inside_iface:
            last_decl = i
            i += 1
            continue
        # Module-level contains inside this procedure = end of spec section
        if RE_CONTAINS.match(lo):
            break
        if is_declaration(all_lines[i]):
            last_decl = i
            # Consume continuation lines — they belong to the same declaration
            # even if they don't look like declarations individually.
            while ends_with_continuation(all_lines[i]) and i + 1 < n - 1:
                i += 1
                last_decl = i
        else:
            break
        i += 1
    return last_decl


def extract_procedures(lines: list[str]) -> list[Procedure]:
    """
    Extract top-level procedures from the lines that follow a module-level 'contains'.
    Handles nested subroutines/functions by tracking depth.
    """
    procedures = []
    n = len(lines)
    i = 0

    while i < n:
        lo = slow(lines[i])
        ms = RE_SUB_START.match(lo)
        mf = RE_FUNC_START.match(lo) if not ms else None

        if ms or mf:
            is_func = bool(mf)
            match = mf if is_func else ms
            prefixes_raw = match.group(1)
            name = match.group(2)          # lowercase (from lo)
            # Recover original-case name
            m_orig = (RE_FUNC_START if is_func else RE_SUB_START).match(lines[i].strip())
            if m_orig:
                name = m_orig.group(2)
            # Strip 'module' from prefixes – we'll add it ourselves
            prefixes = re.sub(r'\bmodule\b', '', prefixes_raw, flags=re.IGNORECASE).strip()

            # Find the end of the procedure header (handle continuation lines)
            header_end = i
            while ends_with_continuation(lines[header_end]) and header_end + 1 < n:
                header_end += 1

            # Collect all lines of the procedure (tracking depth)
            depth = 1
            j = header_end + 1
            proc_lines = list(lines[i:header_end + 1])
            while j < n and depth > 0:
                lj = slow(lines[j])
                if RE_SUB_START.match(lj) or (RE_FUNC_START.match(lj) and not RE_SUB_END.match(lj)):
                    depth += 1
                elif RE_SUB_END.match(lj) or RE_FUNC_END.match(lj):
                    depth -= 1
                proc_lines.append(lines[j])
                j += 1

            if depth != 0:
                # Unmatched: skip
                i = j
                continue

            # Find declarations end
            inner_header_end = header_end - i  # index within proc_lines
            dummy_args = extract_dummy_args(proc_lines[:inner_header_end + 1])
            decl_end = find_procedure_decl_end(proc_lines, inner_header_end, dummy_args)

            procedures.append(Procedure(
                is_function=is_func,
                name=name,
                prefixes=prefixes,
                all_lines=proc_lines,
                header_end=inner_header_end,
                decl_end=decl_end,
            ))
            i = j
        else:
            i += 1

    return procedures


# ---------------------------------------------------------------------------
# Code generation
# ---------------------------------------------------------------------------

def rewrite_procedure_header(header_lines: list[str], is_function: bool,
                              prefixes: str) -> list[str]:
    """
    Rewrite the first line of the procedure header to add 'module' prefix.
    Handles subroutine/function declarations.
    """
    result = list(header_lines)
    first = result[0]
    # Remove existing 'module' prefix if any
    cleaned = re.sub(r'\bmodule\b\s*', '', first, flags=re.IGNORECASE)
    kw = 'function' if is_function else 'subroutine'
    # Find where the keyword is
    m = re.search(kw, cleaned, re.IGNORECASE)
    if not m:
        return result
    before = cleaned[:m.start()]
    after = cleaned[m.start():]
    # Rebuild: preserve leading whitespace + non-module prefixes + 'module ' + rest
    indent = len(first) - len(first.lstrip())
    ws = first[:indent]
    prefix_str = (prefixes + ' ').lstrip() if prefixes else ''
    result[0] = ws + prefix_str + 'module ' + after.lstrip()
    return result


def build_interface_entry(proc: Procedure, dummy_args: list[str]) -> list[str]:
    """Build the module subroutine/function interface block entry."""
    lines = []
    lines.extend(rewrite_procedure_header(proc.header_lines, proc.is_function, proc.prefixes))
    iface_decls, _ = proc.split_decls(dummy_args)
    lines.extend(iface_decls)
    kw = 'function' if proc.is_function else 'subroutine'
    indent = len(proc.all_lines[0]) - len(proc.all_lines[0].lstrip())
    lines.append(' ' * indent + f'end {kw} {proc.name}\n')
    return lines


def build_submodule_procedure(proc: Procedure, dummy_args: list[str]) -> list[str]:
    """Build the 'module procedure' block for the submodule."""
    lines = []
    indent = len(proc.all_lines[0]) - len(proc.all_lines[0].lstrip())
    ws = ' ' * indent
    lines.append(ws + f'module procedure {proc.name}\n')
    iface_decls, local_decls = proc.split_decls(dummy_args)

    # ifx does not propagate array return-variable dimensions from the parent
    # module interface into a module procedure body (scalar types are fine).
    # Re-emit the result declaration locally for array returns only.
    if proc.is_function:
        hdr = ''.join(l.rstrip().rstrip('&') + ' ' for l in proc.header_lines)
        rm = re.search(r'\bresult\s*\(\s*(\w+)\s*\)', hdr, re.IGNORECASE)
        result_name_lo = rm.group(1).lower() if rm else proc.name.lower()
        result_decls = []
        for decl_line in iface_decls:
            lo = decl_line.strip().lower()
            if '::' not in lo:
                continue
            type_part, after = lo.split('::', 1)
            if '!' in after:
                after = after[:after.index('!')]
            names_here = {p.strip().split('(')[0].split('=')[0].strip()
                          for p in split_top_level(after) if p.strip()}
            if result_name_lo not in names_here:
                continue
            # Only re-emit when dimension expressions reference variables:
            # ifx loses dimension info for runtime-sized arrays (e.g.
            # dimension(CS%nk+1)) but handles constant-size and deferred-shape
            # arrays correctly.  Check for alphabetic chars in the dimension
            # content as a proxy for "references a variable".
            dm = re.search(r'\bdimension\s*\(([^)]*)\)', type_part)
            if dm:
                is_array = bool(re.search(r'[a-zA-Z_]', dm.group(1)))
            else:
                # Old-style bounds on variable name: real :: f(n)
                is_array = any(
                    re.search(r'[a-zA-Z_]', re.search(r'\(([^)]*)\)', p).group(1))
                    for p in split_top_level(after)
                    if p.strip().split('(')[0].split('=')[0].strip() == result_name_lo
                    and '(' in p and re.search(r'\(([^)]*)\)', p)
                )
            if is_array:
                result_decls.append(decl_line)
        local_decls = result_decls + local_decls

    lines.extend(local_decls)
    lines.extend(proc.body_lines)
    lines.append(ws + f'end procedure {proc.name}\n')
    return lines


# ---------------------------------------------------------------------------
# File-level transformation
# ---------------------------------------------------------------------------

def remove_stale_submodule(path: Path, suffix: str, dry_run: bool,
                           verbose: bool) -> None:
    """Remove a generated submodule when its parent is not converted."""
    smod_path = path.with_name(path.stem + suffix + path.suffix)
    if not smod_path.exists():
        return

    try:
        first_line = smod_path.read_text(encoding='utf-8', errors='replace').splitlines()[0]
    except IndexError:
        first_line = ''
    except OSError as e:
        print(f'WARNING: cannot read {smod_path}: {e}', file=sys.stderr)
        return

    if not RE_SUBMODULE.match(first_line.strip().lower()):
        return

    if dry_run:
        print(f'  DRY-RUN: would remove stale {smod_path}')
        return

    smod_path.unlink()
    if verbose:
        print(f'  REMOVED stale {smod_path.name}')


def convert_file(path: Path, suffix: str, dry_run: bool,
                 verbose: bool) -> Optional[Conversion]:
    """
    Convert one file. Returns conversion details if a conversion was written
    (or would be in dry_run).
    """
    try:
        text = path.read_text(encoding='utf-8', errors='replace')
    except OSError as e:
        print(f'WARNING: cannot read {path}: {e}', file=sys.stderr)
        return None

    lines = text.splitlines(keepends=True)

    module_name = find_module_name(lines)
    if not module_name:
        return None

    contains_idx = find_module_contains(lines)
    if contains_idx is None:
        if verbose:
            print(f'  SKIP {path.name}: no module-level contains')
        remove_stale_submodule(path, suffix, dry_run, verbose)
        return None

    end_module_idx = find_end_module(lines, contains_idx)
    if end_module_idx is None:
        print(f'WARNING: {path.name}: cannot find end module', file=sys.stderr)
        return None

    # Lines of the module declaration section (before contains)
    module_decl_lines = lines[:contains_idx]
    # Lines between contains and end module (exclusive)
    proc_lines = lines[contains_idx + 1: end_module_idx]
    end_line = lines[end_module_idx]

    # Parse procedures
    procedures = extract_procedures(proc_lines)
    if not procedures:
        if verbose:
            print(f'  SKIP {path.name}: no procedures found after contains')
        remove_stale_submodule(path, suffix, dry_run, verbose)
        return None

    # -----------------------------------------------------------------------
    # Build new parent module
    # -----------------------------------------------------------------------
    new_module_lines = list(module_decl_lines)
    # Collect any non-procedure lines in the proc_lines (comments, blank lines
    # between procedures) that appear before the first procedure – keep them.
    # Actually just add the interface block.
    new_module_lines.append('\n')
    # Use two-space indent for interface block (MOM6 style)
    new_module_lines.append('  interface\n')
    for proc in procedures:
        dummy_args = extract_dummy_args(proc.header_lines)
        new_module_lines.extend(build_interface_entry(proc, dummy_args))
    new_module_lines.append('  end interface\n')
    new_module_lines.append('\n')
    new_module_lines.append(end_line)

    # -----------------------------------------------------------------------
    # Build submodule
    # -----------------------------------------------------------------------
    smod_name = module_name + suffix
    submod_lines = []
    submod_lines.append(f'submodule ({module_name}) {smod_name}\n')
    submod_lines.extend(submodule_preprocessor_lines(module_decl_lines))
    submod_lines.append('  implicit none\n')
    submod_lines.append('contains\n')
    for proc in procedures:
        dummy_args = extract_dummy_args(proc.header_lines)
        submod_lines.extend(build_submodule_procedure(proc, dummy_args))
    submod_lines.append(f'end submodule {smod_name}\n')

    # -----------------------------------------------------------------------
    # Write
    # -----------------------------------------------------------------------
    smod_path = path.with_name(path.stem + suffix + path.suffix)

    if dry_run:
        print(f'  DRY-RUN: would rewrite {path}')
        print(f'  DRY-RUN: would write   {smod_path}')
        return Conversion(path, smod_path, module_name, smod_name)

    path.write_text(''.join(new_module_lines), encoding='utf-8')
    smod_path.write_text(''.join(submod_lines), encoding='utf-8')

    if verbose:
        print(f'  CONVERTED {path.name} -> {path.name} + {smod_path.name}'
              f'  ({len(procedures)} procedures)')
    return Conversion(path, smod_path, module_name, smod_name)


def write_dependency_file(path: Path, conversions: list[Conversion],
                          dry_run: bool) -> None:
    """Write make rules that force submodules to compile after parents."""
    if not conversions:
        return

    lines = [
        '# Generated by convert_to_submodules.py.\n',
        '# Include this from make-based builds after converting sources.\n',
        '# Fortran submodules need fresh ancestor module sidecars from the parent object.\n',
        '\n',
    ]
    for conv in conversions:
        parent_obj = conv.parent.with_suffix('.o').name
        submodule_obj = conv.submodule.with_suffix('.o').name
        ancestor = conv.module_name.lower()
        parent_mod = ancestor + '.mod'
        parent_smod = ancestor + '.smod'
        lines.append(f'# {conv.submodule}: submodule of {conv.parent}\n')
        lines.append(f'{parent_mod} {parent_smod}: {parent_obj}\n')
        lines.append(f'{submodule_obj}: {parent_obj} {parent_mod} {parent_smod}\n')
        lines.append('\n')

    if dry_run:
        print(f'  DRY-RUN: would write   {path}')
        return

    path.write_text(''.join(lines), encoding='utf-8')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('dirs', nargs='+', help='Source directories to process')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print what would be done without writing files')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print one line per converted file')
    parser.add_argument('--suffix', default='_s',
                        help='Suffix for submodule files (default: _s)')
    parser.add_argument('--deps-file', default='submodule_deps.mk',
                        help='Make dependency fragment to write (default: submodule_deps.mk)')
    args = parser.parse_args()

    converted = 0
    skipped = 0
    conversions = []

    for d in args.dirs:
        p = Path(d)
        if not p.exists():
            print(f'ERROR: {d} does not exist', file=sys.stderr)
            continue
        files = sorted(p.rglob('*.F90'))
        print(f'Processing {len(files)} files in {d} ...')
        for f in files:
            if f.name in SKIP_FILES:
                if args.verbose:
                    print(f'  SKIP {f.name}: excluded from conversion')
                skipped += 1
                continue
            # Skip files that are already submodule files (if re-running)
            if f.stem.endswith(args.suffix):
                continue
            conversion = convert_file(f, args.suffix, args.dry_run, args.verbose)
            if conversion:
                conversions.append(conversion)
                converted += 1
            else:
                skipped += 1

    write_dependency_file(Path(args.deps_file), conversions, args.dry_run)

    print(f'\nDone: {converted} converted, {skipped} skipped.')
    if conversions:
        print(f'Wrote submodule compile dependencies to {args.deps_file}.')


if __name__ == '__main__':
    main()
