#!/usr/bin/env python3
"""Remove standalone ngtests lines that were moved inside if blocks."""
import re
from pathlib import Path

F90_DIR = Path(__file__).resolve().parent.parent / 'src' / 'fortran' / 'cloudmask'

NGTESTS_RE = re.compile(r'^\s+ngtests\(\d+\)\s*=\s*ngtests\(\d+\)\s*\+\s*1\s*$')

def fix_file(filepath):
    """Remove standalone ngtests lines after conf_test or cmin."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    new_lines = []
    removed = 0

    for i, line in enumerate(lines):
        stripped = line.strip()

        # Check if this is a ngtests line
        if NGTESTS_RE.match(line):
            # Look at previous non-empty line
            prev_idx = len(new_lines) - 1
            while prev_idx >= 0 and new_lines[prev_idx].strip() == '':
                prev_idx -= 1

            if prev_idx >= 0:
                prev_stripped = new_lines[prev_idx].strip()
                # Remove if preceded by conf_test, cmin, or end if
                if ('conf_test' in prev_stripped or
                    prev_stripped.startswith('cmin') or
                    prev_stripped.startswith('!') or
                    prev_stripped == 'end if'):
                    removed += 1
                    continue

        new_lines.append(line)

    if removed > 0:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        return removed
    return 0

def main():
    total = 0
    for f90_file in sorted(F90_DIR.glob('*.f90')):
        removed = fix_file(f90_file)
        if removed > 0:
            print(f"Fixed {f90_file.name}: removed {removed} standalone ngtests lines")
            total += removed
    print(f"\nTotal: {total} lines removed")

if __name__ == '__main__':
    main()
