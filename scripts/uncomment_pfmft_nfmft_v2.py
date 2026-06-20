#!/usr/bin/env python3
"""Carefully uncomment PFMFT/NFMFT blocks, adding missing endif statements.

This script:
1. Finds blocks between "=== PFMFT test disabled" and "=== PFMFT test disabled end ==="
2. Uncomments each line (removes leading '!' added for the disable)
3. Adds 'endif' before the end marker to close the outer if
4. Removes start/end markers
"""

import re
from pathlib import Path

CLOUDMASK_DIR = Path('/home/liusy2020/yuq/cloudmask/fy3_cloudmask/src/fortran/cloudmask')


def uncomment_block(lines, start_idx, end_idx, test_name):
    """Uncomment lines between start_idx (inclusive) and end_idx (exclusive).

    Adds endif before end marker. Returns list of replacement lines.
    """
    result = []
    block_lines = lines[start_idx:end_idx]

    for line in block_lines:
        stripped = line.lstrip()
        if stripped.startswith('!!'):
            # Pre-existing comment (had ! before disable) → keep as single comment
            idx = line.index('!!')
            result.append(line[:idx] + '!' + line[idx+2:])
        elif stripped.startswith('!'):
            # Disable comment → remove the leading '!'
            # Find the position of '!' (may have whitespace before it)
            idx = line.index('!')
            # Remove just the '!' character, preserving rest of indentation
            result.append(line[:idx] + ' ' + line[idx+1:])
        else:
            result.append(line)

    # Add endif before what would have been the end marker
    result.append('      endif\n')

    return result


def fix_file(filepath):
    with open(filepath) as f:
        lines = f.readlines()

    original = ''.join(lines)
    new_lines = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # Check for PFMFT block start
        if '=== PFMFT test disabled (btclr requires NWP RTM) ===' in line:
            # Find end marker
            end_idx = None
            for j in range(i+1, len(lines)):
                if '=== PFMFT test disabled end ===' in lines[j]:
                    end_idx = j
                    break

            if end_idx is None:
                new_lines.append(line)
                i += 1
                continue

            # Replace with uncommented block
            new_lines.append(line.replace(
                '=== PFMFT test disabled (btclr requires NWP RTM) ===',
                '=== PFMFT test (btclr from NWP sfctmp) ==='))
            uncommented = uncomment_block(lines, i+1, end_idx, 'PFMFT')
            new_lines.extend(uncommented)
            i = end_idx + 1  # Skip past end marker
            continue

        # Check for NFMFT block start
        if '=== NFMFT test disabled (btclr requires NWP RTM) ===' in line:
            end_idx = None
            for j in range(i+1, len(lines)):
                if '=== NFMFT test disabled end ===' in lines[j]:
                    end_idx = j
                    break

            if end_idx is None:
                new_lines.append(line)
                i += 1
                continue

            new_lines.append(line.replace(
                '=== NFMFT test disabled (btclr requires NWP RTM) ===',
                '=== NFMFT test (btclr from NWP sfctmp) ==='))
            uncommented = uncomment_block(lines, i+1, end_idx, 'NFMFT')
            new_lines.extend(uncommented)
            i = end_idx + 1
            continue

        new_lines.append(line)
        i += 1

    new_content = ''.join(new_lines)
    if new_content != original:
        with open(filepath, 'w') as f:
            f.write(new_content)
        return True
    return False


def main():
    f90_files = sorted(CLOUDMASK_DIR.glob('*.f90'))
    changed = []
    for fp in f90_files:
        if fix_file(fp):
            changed.append(fp.name)
            print(f"  Uncommented: {fp.name}")

    print(f"\nChanged {len(changed)} files")


if __name__ == '__main__':
    main()
