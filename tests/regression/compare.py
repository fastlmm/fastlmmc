"""Compare a regression run's outputs with the expected outputs, allowing floating-point noise.

Usage: python3 compare.py <expected_dir> <actual_dir> [--rel 1e-4] [--abs 1e-8]

Text files are compared token by token: non-numeric tokens must match exactly; numbers must
satisfy |a - b| <= max(abs, rel * max(|a|, |b|)). When rows differ only in order (near-tied
p-values sorted differently), rows are re-aligned by their non-numeric tokens.

*.bin files are float64 arrays. U.bin (eigenvectors, one column per entry of the sibling
S.bin) is compared column by column up to sign, skipping columns whose eigenvalue is repeated,
since any orthonormal basis of a repeated eigenvalue's eigenspace is equally correct.
"""
import argparse
import fnmatch
import math
import os
import struct
import sys


def is_close(x, y, rel, abs_tol):
    return abs(x - y) <= max(abs_tol, rel * max(abs(x), abs(y)))


def to_float(token):
    try:
        return float(token)
    except ValueError:
        return None


def compare_tokens(row_a, row_b, rel, abs_tol):
    """Return None if the rows match, else a short description of the first difference."""
    if len(row_a) != len(row_b):
        return f"token count {len(row_a)} vs {len(row_b)}"
    for token_a, token_b in zip(row_a, row_b):
        if token_a == token_b:
            continue
        x, y = to_float(token_a), to_float(token_b)
        if x is None or y is None:
            return f"{token_a!r} vs {token_b!r}"
        if math.isnan(x) or math.isnan(y) or not is_close(x, y, rel, abs_tol):
            return f"{token_a} vs {token_b}"
    return None


def row_key(row):
    return tuple(token for token in row if to_float(token) is None)


def compare_text(path_a, path_b, rel, abs_tol):
    with open(path_a, encoding="latin-1") as file_a, open(path_b, encoding="latin-1") as file_b:
        rows_a = [line.split() for line in file_a]
        rows_b = [line.split() for line in file_b]
    if len(rows_a) != len(rows_b):
        return [f"line count {len(rows_a)} vs {len(rows_b)}"]
    problems = []
    for line_number, (row_a, row_b) in enumerate(zip(rows_a, rows_b), 1):
        problem = compare_tokens(row_a, row_b, rel, abs_tol)
        if problem:
            problems.append(f"line {line_number}: {problem}")
    if not problems:
        return []
    # Retry with rows aligned by their non-numeric tokens (e.g. SNP ids), if those are unique.
    keyed_b = {row_key(row): row for row in rows_b}
    if len(keyed_b) != len(rows_b) or len({row_key(row) for row in rows_a}) != len(rows_a):
        return problems
    realigned = []
    for row_a in rows_a:
        row_b = keyed_b.get(row_key(row_a))
        problem = "missing row" if row_b is None else compare_tokens(row_a, row_b, rel, abs_tol)
        if problem:
            realigned.append(f"row {' '.join(row_key(row_a))}: {problem}")
    return realigned


def read_doubles(path):
    with open(path, "rb") as file:
        data = file.read()
    if len(data) % 8:
        raise ValueError(f"{path}: size {len(data)} is not a multiple of 8")
    return struct.unpack(f"<{len(data) // 8}d", data)


def compare_bin(path_a, path_b, rel, abs_tol):
    values_a, values_b = read_doubles(path_a), read_doubles(path_b)
    if len(values_a) != len(values_b):
        return [f"length {len(values_a)} vs {len(values_b)}"]
    eigenvalues_path = os.path.join(os.path.dirname(path_a), "S.bin")
    if os.path.basename(path_a) == "U.bin" and os.path.exists(eigenvalues_path):
        return compare_eigenvectors(values_a, values_b, read_doubles(eigenvalues_path), rel, abs_tol)
    return [f"index {index}: {x} vs {y}" for index, (x, y) in enumerate(zip(values_a, values_b))
            if not is_close(x, y, rel, abs_tol)][:5]


def compare_eigenvectors(values_a, values_b, eigenvalues, rel, abs_tol):
    column_count = len(eigenvalues)
    assert len(values_a) % column_count == 0, "U.bin size is not a multiple of len(S.bin)"
    row_count = len(values_a) // column_count
    scale = max(abs(value) for value in eigenvalues)
    problems = []
    for column in range(column_count):
        repeated = any(other != column and abs(eigenvalues[other] - eigenvalues[column]) <= 1e-8 * scale
                       for other in range(column_count))
        if repeated:
            continue
        start = column * row_count  # column-major
        column_a = values_a[start:start + row_count]
        column_b = values_b[start:start + row_count]
        same = all(is_close(x, y, rel, abs_tol) for x, y in zip(column_a, column_b))
        flipped = all(is_close(x, -y, rel, abs_tol) for x, y in zip(column_a, column_b))
        if not (same or flipped):
            problems.append(f"eigenvector column {column} (eigenvalue {eigenvalues[column]:.6g}) differs")
    return problems[:5]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("expected_dir")
    parser.add_argument("actual_dir")
    parser.add_argument("--rel", type=float, default=1e-4)
    parser.add_argument("--abs", type=float, default=1e-8)
    parser.add_argument("--ignore", action="append", default=[], metavar="GLOB",
                        help="skip outputs whose relative path matches GLOB (repeatable)")
    args = parser.parse_args()

    expected_files, actual_files = set(), set()
    for root_dir, found in ((args.expected_dir, expected_files), (args.actual_dir, actual_files)):
        for directory, _, names in os.walk(root_dir):
            for name in names:
                found.add(os.path.relpath(os.path.join(directory, name), root_dir))

    failures = 0
    for relative in sorted(expected_files | actual_files):
        if any(fnmatch.fnmatch(relative, pattern) for pattern in args.ignore):
            print(f"skip {relative}")
            continue
        if relative not in actual_files:
            problems = ["missing from this run"]
        elif relative not in expected_files:
            problems = ["not in expected outputs (run with --update if it is new)"]
        else:
            path_a = os.path.join(args.expected_dir, relative)
            path_b = os.path.join(args.actual_dir, relative)
            compare = compare_bin if relative.endswith(".bin") else compare_text
            problems = compare(path_a, path_b, args.rel, args.abs)
        if problems:
            failures += 1
            print(f"FAIL {relative}")
            for problem in problems[:5]:
                print(f"       {problem}")
        else:
            print(f"ok   {relative}")
    skipped = sum(any(fnmatch.fnmatch(relative, pattern) for pattern in args.ignore)
                  for relative in expected_files | actual_files)
    print(f"\n{len(expected_files | actual_files) - failures - skipped} ok, {failures} failed, {skipped} skipped")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
