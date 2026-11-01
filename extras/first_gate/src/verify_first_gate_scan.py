"""Exhaustive First-Gate verification for every even N with 122 <= N < 10007^2
(Computation comp:finite): all 1,225 blocks q^2 < N < r^2 with 11 <= q <= 9973.

For N = 2C in the block of q, put W = q^2 - C.  An offset d survives if neither C - d
nor C + d has a prime divisor below q, i.e. gcd(C^2 - d^2, P_q) = 1 with P_q the product
of the primes below q.  For 0 <= d < W both endpoints lie in [2, q^2), where having no
prime divisor below q is equivalent to being a prime >= q.  Hence N satisfies First-Gate
(Proposition prop:full-FG-window) if and only if some d < W survives.  The scan records
the least surviving offset lambda(N) from a prime table; a failure is a target with no
survivor below W.  For two targets per block, lambda(N) is recomputed from the gcd
definition with exact integers, independently of the prime table.

Run from this directory:  python verify_first_gate_scan.py
"""
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import flint
import numpy as np
from flint import fmpz

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANUSCRIPT = ROOT / "first_gate_prime_square_blocks.tex"
CERTIFICATE = ROOT / "certificates" / "first_gate_scan_q9973.json"
LOG = ROOT / "logs" / "first_gate_scan_q9973.log"

Q_MIN, Q_MAX, R_LAST = 11, 9973, 10007
LIMIT = 10 ** 8
PI_CHECKPOINTS = {10: 4, 10 ** 2: 25, 10 ** 3: 168, 10 ** 4: 1229, 10 ** 5: 9592,
                  10 ** 6: 78498, 10 ** 7: 664579, 10 ** 8: 5761455}
WINDOW = 200000


class Log:
    def __init__(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = path.open("w", encoding="utf-8", newline="\n")

    def __call__(self, msg=""):
        print(msg, flush=True)
        self.fh.write(msg + "\n")
        self.fh.flush()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prime_table(n):
    isp = np.ones(n + 1, dtype=bool)
    isp[:2] = False
    isp[4::2] = False
    for p in range(3, math.isqrt(n) + 1, 2):
        if isp[p]:
            isp[p * p::2 * p] = False
    return isp


def validate_table(isp, log):
    pi = {x: int(np.count_nonzero(isp[:x + 1])) for x in PI_CHECKPOINTS}
    pi_ok = pi == PI_CHECKPOINTS
    log(f"pi(x) checkpoints {pi}: {'OK' if pi_ok else 'FAIL'}")
    windows = [(0, WINDOW), (LIMIT - WINDOW, LIMIT)]
    mismatches = 0
    for lo, hi in windows:
        mismatches += sum(bool(isp[n]) != fmpz(n).is_prime() for n in range(lo, hi + 1))
    log(f"prime table vs FLINT is_prime on {windows}: {mismatches} mismatches")
    return {"pi_checkpoints": {str(k): v for k, v in pi.items()}, "pi_checkpoints_pass": pi_ok,
            "flint_windows": [list(w) for w in windows], "flint_mismatches": mismatches,
            "pass": pi_ok and mismatches == 0}


def scan_block(isp, q, r):
    X = q * q
    N = np.arange(X + 1, r * r, 2, dtype=np.int64)
    C = N // 2
    W = X - C
    lam = np.full(N.size, -1, dtype=np.int64)
    d = (C + 1) & 1
    active = np.arange(N.size)
    failures = []
    while active.size:
        da, Ca = d[active], C[active]
        lo = Ca - da
        # In the window d < W both endpoints lie in [2, q^2): q-rough means prime >= q.
        hit = (da < W[active]) & isp[lo] & isp[Ca + da] & (lo >= q)
        lam[active[hit]] = da[hit]
        active = active[~hit]
        d[active] += 2
        out = d[active] >= W[active]
        if out.any():
            failures.extend(N[active[out]].tolist())
            active = active[~out]
    return N, C, W, lam, failures


def first_survivor_gcd(N, P, bound):
    C = N // 2
    for d in range(bound + 1):
        if math.gcd(C * C - d * d, P) == 1:
            return d
    return None


def main():
    t0 = time.time()
    log = Log(LOG)
    log("Exhaustive First-Gate scan, Computation comp:finite")
    log(f"started {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    env = {"python": sys.version.split()[0], "python-flint": flint.__version__,
           "FLINT": flint.__FLINT_VERSION__, "numpy": np.__version__, "platform": platform.platform()}
    sources = {Path(__file__).name: sha256(Path(__file__).resolve())}
    log(f"environment {env}")
    log(f"manuscript sha256 {sha256(MANUSCRIPT)}")
    log(f"source {Path(__file__).name} sha256 {sources[Path(__file__).name]}")

    isp = prime_table(LIMIT)
    log(f"prime table to {LIMIT} built ({time.time() - t0:.0f}s)")
    table = validate_table(isp, log)

    primes = np.nonzero(isp[:R_LAST + 1])[0].tolist()
    blocks = [(q, r) for q, r in zip(primes, primes[1:]) if Q_MIN <= q <= Q_MAX]
    geometry = {
        "blocks": len(blocks),
        "first_q": blocks[0][0], "last_q": blocks[-1][0], "last_r": blocks[-1][1],
        "consecutive": all(b[1] == c[0] for b, c in zip(blocks, blocks[1:])),
        "covers_122_to_10007^2": blocks[0][0] ** 2 + 1 == 122 and blocks[-1][1] == R_LAST,
        "largest_square_below_limit": Q_MAX ** 2 < LIMIT,
    }
    geometry["pass"] = (geometry["blocks"] == 1225 and geometry["first_q"] == Q_MIN
                        and geometry["last_q"] == Q_MAX and geometry["last_r"] == R_LAST
                        and geometry["consecutive"] and geometry["covers_122_to_10007^2"]
                        and geometry["largest_square_below_limit"])
    log(f"block geometry {geometry}")

    digest = hashlib.sha256()
    targets = 0
    failures = []
    next_N = 122
    coverage_ok = True
    worst = (-1, 0, 0)
    samples = []
    for i, (q, r) in enumerate(blocks, 1):
        N, C, W, lam, fails = scan_block(isp, q, r)
        coverage_ok &= bool(N[0] == next_N and N[-1] == r * r - 1 and N.size == (r * r - q * q) // 2
                            and np.all(N % 2 == 0) and np.all(W > 0))
        next_N = r * r + 1
        targets += N.size
        failures += fails
        digest.update(lam.astype("<i8").tobytes())
        j = int(lam.argmax())
        if lam[j] > worst[0]:
            worst = (int(lam[j]), int(N[j]), q)
        k = 1 + int(lam[1:].argmax())
        samples.append((q, int(N[0]), int(lam[0])))
        samples.append((q, int(N[k]), int(lam[k])))
        if i % 100 == 0 or i == len(blocks):
            log(f"  block {i:4d}/{len(blocks)} q={q:5d} r={r:5d} targets so far {targets:9d} "
                f"failures {len(failures)} max lambda {worst[0]} ({time.time() - t0:.0f}s)")
    expected = (R_LAST ** 2 - 121) // 2
    coverage_ok &= next_N == R_LAST ** 2 + 1 and targets == expected
    log(f"even targets scanned {targets} (expected {expected}); coverage {'OK' if coverage_ok else 'FAIL'}")
    log(f"First-Gate failures: {len(failures)}")
    log(f"largest first survivor lambda = {worst[0]} at N = {worst[1]} (q = {worst[2]})")
    log(f"sha256 of first-survivor array: {digest.hexdigest()}")

    mismatches = []
    P = 1
    p_index = 0
    small = [p for p in primes if p < Q_MAX]
    for q, N, lam in samples:
        while p_index < len(small) and small[p_index] < q:
            P *= small[p_index]
            p_index += 1
        C = N // 2
        d = first_survivor_gcd(N, P, lam)
        if d != lam or not (fmpz(C - d).is_prime() and fmpz(C + d).is_prime()
                            and q <= C - d and C + d < q * q):
            mismatches.append({"q": q, "N": N, "lambda_table": lam, "lambda_gcd": d})
    log(f"gcd(C^2-d^2, P_q) cross-checks: {len(samples)} targets "
        f"(first target and largest-lambda target of each block), mismatches {len(mismatches)}")

    status = "PASS" if (table["pass"] and geometry["pass"] and coverage_ok and not failures
                        and not mismatches and len(samples) == 2450) else "FAIL"
    cert = {
        "certificate": "first_gate_scan_q9973",
        "manuscript_claims": ["Computation comp:finite", "Proposition prop:full-FG-window"],
        "manuscript": {"file": MANUSCRIPT.name, "sha256": sha256(MANUSCRIPT)},
        "sources": sources,
        "environment": env,
        "range": {"q": [Q_MIN, Q_MAX], "r_last": R_LAST, "N": [122, R_LAST ** 2 - 1], "even_only": True},
        "prime_table": dict(limit=LIMIT, **table),
        "block_geometry": geometry,
        "coverage": {"even_targets": targets, "expected": expected, "pass": coverage_ok},
        "survivor_predicate": "d in [0, W) with C-d, C+d both prime >= q; equivalent to "
                              "gcd(C^2-d^2, prod_{p<q} p) = 1 for d < W",
        "first_gate_failures": failures,
        "max_first_survivor": {"lambda": worst[0], "N": worst[1], "q": worst[2]},
        "first_survivor_sha256": digest.hexdigest(),
        "first_survivor_encoding": "int64 little-endian lambda(N) for all targets in increasing N",
        "gcd_cross_check": {"targets": len(samples),
                            "selection": "per block: N = q^2+1 and the largest-lambda target among the rest",
                            "mismatches": mismatches},
        "status": status,
    }
    CERTIFICATE.parent.mkdir(parents=True, exist_ok=True)
    CERTIFICATE.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8", newline="\n")
    log(f"certificate written: {CERTIFICATE.relative_to(ROOT).as_posix()}")
    log(f"STATUS {status}")
    log(f"elapsed {time.time() - t0:.0f}s")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
