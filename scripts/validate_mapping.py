#!/usr/bin/env python3
import sys, json

if len(sys.argv) < 2:
    print("Usage: validate_mapping.py <mapping_file>", file=sys.stderr)
    sys.exit(1)

p = sys.argv[1]
count = 0
try:
    with open(p, 'r') as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                json.loads(line)
                count += 1
            except Exception as e:
                print(f"Malformed JSON at line {i}: {e}", file=sys.stderr)
                sys.exit(2)
except FileNotFoundError:
    print(f"File not found: {p}", file=sys.stderr)
    sys.exit(3)

print(f"Valid JSON lines: {count}")
