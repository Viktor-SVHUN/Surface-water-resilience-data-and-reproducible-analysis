"""Verify all distributed files against MANIFEST.sha256 using the standard library."""
from pathlib import Path
import hashlib
import sys


def verify(root):
    root = Path(root).resolve()
    failures = []
    count = 0
    for line in (root / 'MANIFEST.sha256').read_text(encoding='utf-8').splitlines():
        expected, relative = line.split('  ', 1)
        path = (root / relative).resolve()
        if root not in path.parents or not path.is_file():
            failures.append(f'Missing or invalid path: {relative}')
            continue
        with path.open('rb') as stream:
            actual = hashlib.file_digest(stream, 'sha256').hexdigest()
        if actual != expected:
            failures.append(f'Checksum mismatch: {relative}')
        count += 1
    return count, failures


if __name__ == '__main__':
    count, failures = verify(Path(__file__).resolve().parents[1])
    if failures:
        print('\n'.join(failures), file=sys.stderr)
        raise SystemExit(1)
    print(f'Integrity verified: {count} distributed files.')
