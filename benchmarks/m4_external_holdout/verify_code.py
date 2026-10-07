#!/usr/bin/env python3
"""Reject edits to the delivered validation code before using cached stages."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def verify():
    freeze = ROOT / 'CODE_FREEZE.json'
    if not freeze.is_file():
        raise SystemExit('Missing CODE_FREEZE.json; use the complete delivered package.')
    for name, expected in json.loads(freeze.read_text())['files'].items():
        path = ROOT / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit(f'Frozen validation file changed: {name}. Keep this experiment unchanged; use a separate package and new output for revised rules.')
    print('Validation code and protocol hashes verified.')


if __name__ == '__main__':
    verify()
