#!/usr/bin/env python3

import os
import re
import sys

# Usage: python find_missing.py /path/to/folder

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <folder>")
    sys.exit(1)

folder = sys.argv[1]

pattern = re.compile(r"^tnpNtuple_Reference_(\d+)\.root$")

found = set()

for filename in os.listdir(folder):
    m = pattern.match(filename)
    if m:
        found.add(int(m.group(1)))

expected = set(range(1002))  # 0 through 1001 inclusive
missing = sorted(expected - found)
extra = sorted(found - expected)

if missing:
    print("Missing files:")
    for x in missing:
        print(f"  tnpNtuple_Reference_{x}.root")
else:
    print("No missing files.")

if extra:
    print("\nUnexpected file numbers:")
    for x in extra:
        print(f"  {x}")
