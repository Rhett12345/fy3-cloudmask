#!/usr/bin/env python3
"""Fix nmtests/ngtests increments in all Fortran modules.

In Python, nmtests and ngtests are only incremented when the test condition passes.
In Fortran, they were incremented unconditionally when the gate passes.
This script moves them inside the test condition blocks.
"""
import re
import os
from pathlib import Path

F90_DIR = Path(__file__).resolve().parent.parent / 'src' / 'fortran' / 'cloudmask'

def fix_btd_11_12_block(content):
    """Fix the BTD 11-12um test block (thin cirrus test).

    Pattern:
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf1.le.dfthrsh) then
          call set_bit(testbits,18)
          nptests = nptests + 1
        end if
        ...
        ngtests(2) = ngtests(2) + 1

    Fix: Move nmtests, set_qa_bit, ngtests inside the if block.
    """
    # Pattern for the BTD 11-12 test block
    pattern = r'([ \t]+)nmtests = nmtests \+ 1\n([ \t]+)call set_qa_bit\(qa_bits,18\)\n([ \t]+)if \(masdf1\.le\.dfthrsh\) then\n([ \t]+)call set_bit\(testbits,18\)\n([ \t]+)nptests = nptests \+ 1\n([ \t]+)end if'

    replacement = r'\1if (masdf1.le.dfthrsh) then\n\1  nmtests = nmtests + 1\n\1  call set_qa_bit(qa_bits,18)\n\1  call set_bit(testbits,18)\n\1  nptests = nptests + 1\n\1  ngtests(2) = ngtests(2) + 1\n\1end if'

    content = re.sub(pattern, replacement, content)

    # Also remove the standalone ngtests(2) line that was after the if block
    # Pattern: ngtests(2) = ngtests(2) + 1 (standalone, not inside if)
    # This is tricky - we need to find lines that are NOT inside an if block
    # Let's just remove lines that match the pattern after conf_test calls
    lines = content.split('\n')
    new_lines = []
    skip_next_ngtests = False

    for i, line in enumerate(lines):
        stripped = line.strip()
        # Skip standalone ngtests(2) lines that come after conf_test
        if stripped == 'ngtests(2) = ngtests(2) + 1' and i > 0:
            prev_stripped = lines[i-1].strip()
            if 'conf_test' in prev_stripped or 'cmin2 = min' in prev_stripped:
                continue
        new_lines.append(line)

    return '\n'.join(new_lines)


def fix_visible_test_block(content):
    """Fix the visible reflectance test block (0.64um).

    Pattern:
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,20)
        if (masv66.le.dlref1(2)) then
          call set_bit(testbits,20)
          nptests = nptests + 1
        end if
        call conf_test(...)
        cmin3 = min(cmin3,c5)
        ngtests(3) = ngtests(3) + 1

    Fix: Move nmtests, set_qa_bit, ngtests inside the if block.
    """
    # Pattern for visible test
    pattern = r'([ \t]+)nmtests = nmtests \+ 1\n([ \t]+)call set_qa_bit\(qa_bits,20\)\n([ \t]+)if \(masv66\.le\.dlref1(?:_t2)?\(2\)\) then\n([ \t]+)call set_bit\(testbits,20\)\n([ \t]+)nptests = nptests \+ 1\n([ \t]+)end if'

    replacement = r'\1if (masv66.le.dlref1\3(2)) then\n\1  nmtests = nmtests + 1\n\1  call set_qa_bit(qa_bits,20)\n\1  call set_bit(testbits,20)\n\1  nptests = nptests + 1\n\1  ngtests(3) = ngtests(3) + 1\n\1end if'

    # This is getting complex with the _t2 variants. Let me use a simpler approach.
    return content


def fix_file(filepath):
    """Fix a single Fortran file."""
    with open(filepath, 'r') as f:
        content = f.read()

    original = content

    # Fix BTD 11-12 test block
    # Look for the pattern where nmtests is incremented before the if block
    # and ngtests(2) is incremented after

    # Pattern 1: BTD 11-12 test (qa_bits,18)
    # Find blocks that look like:
    #   nmtests = nmtests + 1
    #   call set_qa_bit(qa_bits,18)
    #   if (masdf1.le.dfthrsh) then
    #     call set_bit(testbits,18)
    #     nptests = nptests + 1
    #   end if
    #   ... (conf_test lines)
    #   ngtests(2) = ngtests(2) + 1

    # Let's use a more targeted approach - find and replace specific patterns

    # Fix for qa_bits,18 (BTD 11-12)
    old_patterns = [
        # Pattern for qa_bits,18
        (
            r'([ \t]+)nmtests = nmtests \+ 1\n([ \t]+)call set_qa_bit\(qa_bits,18\)\n([ \t]+)if \(masdf1\.le\.dfthrsh\) then\n([ \t]+)call set_bit\(testbits,18\)\n([ \t]+)nptests = nptests \+ 1\n([ \t]+)end if',
            r'\1if (masdf1.le.dfthrsh) then\n\1  nmtests = nmtests + 1\n\1  call set_qa_bit(qa_bits,18)\n\1  call set_bit(testbits,18)\n\1  nptests = nptests + 1\n\1  ngtests(2) = ngtests(2) + 1\n\1end if'
        ),
        # Pattern for qa_bits,20 (visible)
        (
            r'([ \t]+)nmtests = nmtests \+ 1\n([ \t]+)call set_qa_bit\(qa_bits,20\)\n([ \t]+)if \(masv66\.le\.dlref1(?:_t2)?\(2\)\) then\n([ \t]+)call set_bit\(testbits,20\)\n([ \t]+)nptests = nptests \+ 1\n([ \t]+)end if',
            r'\1if (masv66.le.\3(2)) then\n\1  nmtests = nmtests + 1\n\1  call set_qa_bit(qa_bits,20)\n\1  call set_bit(testbits,20)\n\1  nptests = nptests + 1\n\1  ngtests(3) = ngtests(3) + 1\n\1end if'
        ),
    ]

    # Actually, let me use a simpler line-by-line approach
    lines = content.split('\n')
    new_lines = []
    i = 0
    changes_made = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Check for qa_bits,18 pattern (BTD 11-12)
        if stripped == 'nmtests = nmtests + 1' and i+1 < len(lines):
            next_stripped = lines[i+1].strip()
            if next_stripped == 'call set_qa_bit(qa_bits,18)':
                # Look ahead for the if block
                if i+2 < len(lines) and 'if (masdf1.le.dfthrsh)' in lines[i+2]:
                    # Found the pattern - restructure
                    indent = line[:len(line) - len(line.lstrip())]
                    # Add the if block with all contents inside
                    new_lines.append(f'{indent}if (masdf1.le.dfthrsh) then')
                    new_lines.append(f'{indent}  nmtests = nmtests + 1')
                    new_lines.append(f'{indent}  call set_qa_bit(qa_bits,18)')
                    # Skip to the call set_bit line
                    j = i + 3
                    while j < len(lines) and 'call set_bit(testbits,18)' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  call set_bit(testbits,18)')
                    # Find nptests line
                    while j < len(lines) and 'nptests = nptests + 1' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  nptests = nptests + 1')
                    new_lines.append(f'{indent}  ngtests(2) = ngtests(2) + 1')
                    new_lines.append(f'{indent}end if')
                    # Skip past the end if
                    while j < len(lines) and 'end if' not in lines[j]:
                        j += 1
                    i = j + 1
                    changes_made += 1
                    continue

        # Check for qa_bits,20 pattern (visible)
        if stripped == 'nmtests = nmtests + 1' and i+1 < len(lines):
            next_stripped = lines[i+1].strip()
            if next_stripped == 'call set_qa_bit(qa_bits,20)':
                # Look ahead for the if block
                if i+2 < len(lines) and ('if (masv66.le.dlref1(2))' in lines[i+2] or 'if (masv66.le.dlref1_t2(2))' in lines[i+2]):
                    # Found the pattern - restructure
                    indent = line[:len(line) - len(line.lstrip())]
                    ref_var = 'dlref1_t2' if 'dlref1_t2' in lines[i+2] else 'dlref1'
                    new_lines.append(f'{indent}if (masv66.le.{ref_var}(2)) then')
                    new_lines.append(f'{indent}  nmtests = nmtests + 1')
                    new_lines.append(f'{indent}  call set_qa_bit(qa_bits,20)')
                    # Skip to the call set_bit line
                    j = i + 3
                    while j < len(lines) and 'call set_bit(testbits,20)' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  call set_bit(testbits,20)')
                    # Find nptests line
                    while j < len(lines) and 'nptests = nptests + 1' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  nptests = nptests + 1')
                    new_lines.append(f'{indent}  ngtests(3) = ngtests(3) + 1')
                    new_lines.append(f'{indent}end if')
                    # Skip past the end if
                    while j < len(lines) and 'end if' not in lines[j]:
                        j += 1
                    i = j + 1
                    changes_made += 1
                    continue

        # Check for qa_bits,21 pattern (GEMI)
        if stripped == 'nmtests = nmtests + 1' and i+1 < len(lines):
            next_stripped = lines[i+1].strip()
            if next_stripped == 'call set_qa_bit(qa_bits,21)':
                # This is the GEMI test - need to find the vrat check
                indent = line[:len(line) - len(line.lstrip())]
                # Look for the if(vrat...) block
                j = i + 2
                while j < len(lines) and 'if(vrat' not in lines[j] and 'if (vrat' not in lines[j]:
                    j += 1
                if j < len(lines):
                    # Found the vrat check - restructure
                    # First, add any lines between nmtests and the vrat check
                    for k in range(i+2, j):
                        new_lines.append(lines[k])
                    # Now add the if block with nmtests inside
                    new_lines.append(f'{indent}if(vrat .gt. dlvrat{"_t2" if "dlvrat_t2" in lines[j] else ""}(2)) then')
                    new_lines.append(f'{indent}  nmtests = nmtests + 1')
                    new_lines.append(f'{indent}  call set_qa_bit(qa_bits,21)')
                    # Find nptests and set_bit lines
                    k = j + 1
                    while k < len(lines) and 'end if' not in lines[k]:
                        if 'nptests = nptests + 1' in lines[k]:
                            new_lines.append(f'{indent}  nptests = nptests + 1')
                        elif 'call set_bit(testbits,21)' in lines[k]:
                            new_lines.append(f'{indent}  call set_bit(testbits,21)')
                        k += 1
                    new_lines.append(f'{indent}  ngtests(3) = ngtests(3) + 1')
                    new_lines.append(f'{indent}end if')
                    i = k + 1
                    changes_made += 1
                    continue

        # Check for qa_bits,16 pattern (NIR 1.38um)
        if stripped == 'nmtests = nmtests + 1' and i+1 < len(lines):
            next_stripped = lines[i+1].strip()
            if next_stripped == 'call set_qa_bit(qa_bits,16)':
                # Look ahead for the if block
                if i+2 < len(lines) and 'if (masv188.le.dlref3' in lines[i+2]:
                    indent = line[:len(line) - len(line.lstrip())]
                    ref_var = 'dlref3_t2' if 'dlref3_t2' in lines[i+2] else 'dlref3'
                    new_lines.append(f'{indent}if (masv188.le.{ref_var}(2)) then')
                    new_lines.append(f'{indent}  nmtests = nmtests + 1')
                    new_lines.append(f'{indent}  call set_qa_bit(qa_bits,16)')
                    # Skip to the call set_bit line
                    j = i + 3
                    while j < len(lines) and 'call set_bit(testbits,16)' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  call set_bit(testbits,16)')
                    # Find nptests line
                    while j < len(lines) and 'nptests = nptests + 1' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  nptests = nptests + 1')
                    new_lines.append(f'{indent}  ngtests(4) = ngtests(4) + 1')
                    new_lines.append(f'{indent}end if')
                    # Skip past the end if
                    while j < len(lines) and 'end if' not in lines[j]:
                        j += 1
                    i = j + 1
                    changes_made += 1
                    continue

        # Check for qa_bits,19 pattern (11-4um BTD)
        if stripped == 'nmtests = nmtests + 1' and i+1 < len(lines):
            next_stripped = lines[i+1].strip()
            if next_stripped == 'call set_qa_bit(qa_bits,19)':
                # Look ahead for the if block
                if i+2 < len(lines) and 'if (mas11_4.ge.dl11_4lo' in lines[i+2]:
                    indent = line[:len(line) - len(line.lstrip())]
                    lo_var = 'dl11_4lo_t2' if 'dl11_4lo_t2' in lines[i+2] else 'dl11_4lo'
                    new_lines.append(f'{indent}if (mas11_4.ge.{lo_var}(2)) then')
                    new_lines.append(f'{indent}  nmtests = nmtests + 1')
                    new_lines.append(f'{indent}  call set_qa_bit(qa_bits,19)')
                    # Skip to the call set_bit line
                    j = i + 3
                    while j < len(lines) and 'call set_bit(testbits,19)' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  call set_bit(testbits,19)')
                    # Find nptests line
                    while j < len(lines) and 'nptests = nptests + 1' not in lines[j]:
                        j += 1
                    if j < len(lines):
                        new_lines.append(f'{indent}  nptests = nptests + 1')
                    new_lines.append(f'{indent}end if')
                    # Skip past the end if
                    while j < len(lines) and 'end if' not in lines[j]:
                        j += 1
                    i = j + 1
                    changes_made += 1
                    continue

        new_lines.append(line)
        i += 1

    content = '\n'.join(new_lines)

    # Remove standalone ngtests lines that are now redundant
    # (they were moved inside the if blocks)
    lines = content.split('\n')
    new_lines = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        # Skip standalone ngtests(2) lines that come after conf_test
        if stripped == 'ngtests(2) = ngtests(2) + 1':
            # Check if this is inside an if block (should be kept)
            # or standalone (should be removed)
            # If the previous line is 'end if' or 'cmin2 = min', it's standalone
            if i > 0:
                prev_stripped = lines[i-1].strip()
                if prev_stripped.startswith('cmin2 = min') or prev_stripped == 'end if':
                    continue
        # Skip standalone ngtests(3) lines
        if stripped == 'ngtests(3) = ngtests(3) + 1':
            if i > 0:
                prev_stripped = lines[i-1].strip()
                if prev_stripped.startswith('cmin3 = min') or prev_stripped == 'end if':
                    continue
        # Skip standalone ngtests(4) lines
        if stripped == 'ngtests(4) = ngtests(4) + 1':
            if i > 0:
                prev_stripped = lines[i-1].strip()
                if prev_stripped.startswith('cmin4 = min') or prev_stripped == 'end if':
                    continue
        new_lines.append(line)

    content = '\n'.join(new_lines)

    if content != original:
        with open(filepath, 'w') as f:
            f.write(content)
        return changes_made
    return 0


def main():
    """Fix all Fortran files."""
    total_changes = 0
    for f90_file in sorted(F90_DIR.glob('*.f90')):
        changes = fix_file(f90_file)
        if changes > 0:
            print(f"Fixed {f90_file.name}: {changes} blocks")
            total_changes += changes

    print(f"\nTotal: {total_changes} blocks fixed")


if __name__ == '__main__':
    main()
