#!/usr/bin/env python3
"""Remove .false. gate from APOLLO 11-12um BTD thin cirrus test in all Fortran files.

Strategy:
1. Remove .false. .and. from the APOLLO entry condition
2. Remove .false. .and. from the secant calculation (force schi=99.0 → static threshold)
3. Remove .false. .and. from tri-spectral bit-setting (ocean files)
4. Remove .false. .and. from polar night snow guard

This effectively restores the 11-12um thin cirrus test using only the static threshold,
without relying on the MODIS-specific APOLLO lookup table.
"""

import re
from pathlib import Path

CLOUDMASK_DIR = Path('/home/liusy2020/yuq/cloudmask/fy3_cloudmask/src/fortran/cloudmask')


def fix_file(filepath: Path) -> bool:
    """Fix one Fortran file. Returns True if changed."""
    with open(filepath, 'r') as f:
        content = f.read()
    original = content

    # Pattern 1: APOLLO entry condition
    # ".false. .and. nint(masir11)..."  →  "nint(masir11)..."
    # Variants include optional spaces around .and., different line continuations
    content = re.sub(
        r'\.false\.\s*\.and\.\s*('
        r'nint\(masir11\)|'
        r'\.\s*not\.\s*\(antarctic|'
        r'masdf2\.lt\.tri_thres'
        r')',
        r'\1',
        content
    )

    # Pattern 2: Secant calculation guard
    # ".false. .and. abs(cosvza).gt.Rel_equality_EPS"  →  "abs(cosvza).gt.Rel_equality_EPS"
    # But we want to FORCE static threshold, so we comment out the secant branch
    content = re.sub(
        r'\.false\.\s*\.and\.\s*abs\(cosvza\)\.gt\.Rel_equality_EPS',
        r'.false.',
        content
    )

    if content != original:
        with open(filepath, 'w') as f:
            f.write(content)
        return True
    return False


def main():
    f90_files = sorted(CLOUDMASK_DIR.glob('*.f90'))
    changed = []
    for fp in f90_files:
        if fix_file(fp):
            changed.append(fp.name)
            print(f"  Fixed: {fp.name}")

    print(f"\nChanged {len(changed)} files:")
    for name in changed:
        print(f"  - {name}")


if __name__ == '__main__':
    main()
