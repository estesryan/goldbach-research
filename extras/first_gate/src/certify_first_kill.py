"""Certify the first-kill inequality (sp:eq:certified-margin) behind Theorem thm:two-sided-0548.

    f(L/a) > E(a, b, eps),   E = 2 a^2 / omega(1/a) * int_a^b omega((1-al)/al) F(L/al) al^-3 dal,

with a = 2/25, b = 137/500, eps = 1/10^4, L = 1/2 - eps.  With beta = 1/al,

    int_a^b omega((1-al)/al) F(L/al) al^-3 dal = int_{1/b}^{1/a} omega(beta - 1) F(L beta) beta dbeta,

and on each segment between consecutive integer crossings of beta - 1 and L beta the
integrand is a product of two fixed Taylor models and beta, i.e. a polynomial plus a
bounded remainder, which is integrated exactly in ball arithmetic.

Run from this directory:  python certify_first_kill.py
"""
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import flint
import numpy as np
from flint import acb, arb, arb_poly, ctx, fmpq, fmpz

import sieve_functions as sfm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MANUSCRIPT = ROOT / "first_gate_prime_square_blocks.tex"
CERTIFICATE = ROOT / "certificates" / "first_kill_0548.json"
LOG = ROOT / "logs" / "first_kill_0548.log"

A = fmpq(2, 25)
B = fmpq(137, 500)
EPS = fmpq(1, 10 ** 4)
L = fmpq(1, 2) - EPS

REQUIRED = [
    ("f(L/a)", fmpq(99994969, 10 ** 8), fmpq(99994970, 10 ** 8)),
    ("E(a,b,eps)", fmpq(99952912, 10 ** 8), fmpq(99952914, 10 ** 8)),
    ("f(L/a)-E(a,b,eps)", fmpq(42055, 10 ** 8), None),
]


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


def exact(x):
    m, e = x.man_exp()
    return fmpq(m * fmpz(2) ** e) if e >= 0 else fmpq(m, fmpz(2) ** (-e))


def decimal(q, digits, up):
    scaled = q * fmpz(10) ** digits
    n = scaled.ceil() if up else scaled.floor()
    sign = "-" if n < 0 else ""
    s = str(abs(n)).rjust(digits + 1, "0")
    return f"{sign}{s[:-digits]}.{s[-digits:]}"


def enclosure(x):
    lo, hi = exact(x.lower()), exact(x.upper())
    m_lo, e_lo = x.lower().man_exp()
    m_hi, e_hi = x.upper().man_exp()
    return {
        "lower": f"{m_lo}*2^{e_lo}",
        "upper": f"{m_hi}*2^{e_hi}",
        "lower_decimal_outward": decimal(lo, 30, up=False),
        "upper_decimal_outward": decimal(hi, 30, up=True),
        "radius": float(x.rad()),
    }


def closed_form_tests(sf):
    eg2 = 2 * arb.const_euler().exp()

    def J(v):
        """int_2^v log(t-1)/t dt."""
        return acb.integral(lambda t, an: (t - 1).log(analytic=an) / t, 2, arb(v)).real

    def omega_23(u):
        return (1 + arb(u - 1).log()) / arb(u)

    def omega_34(u):
        return (1 + arb(u - 1).log() + J(u - 1)) / arb(u)

    def f_24(s):
        return eg2 * arb(s - 1).log() / arb(s)

    def F_35(s):
        return eg2 * (1 + J(s - 1)) / arb(s)

    h = fmpq(1, 2)
    cases = [
        ("omega", 2, fmpq(5, 2), omega_23), ("omega", 2, 3, omega_23), ("omega", 3, 3, omega_23),
        ("omega", 3, 3 + h, omega_34), ("omega", 3, 4, omega_34), ("omega", 4, 4, omega_34),
        ("f", 2, fmpq(5, 2), f_24), ("f", 2, 3, f_24), ("f", 3, 3, f_24),
        ("f", 3, 3 + h, f_24), ("f", 3, 4, f_24), ("f", 4, 4, f_24),
        ("F", 2, 3, lambda s: eg2 / arb(s)), ("F", 3, 3, lambda s: eg2 / arb(s)),
        ("F", 3, 3 + h, F_35), ("F", 3, 4, F_35), ("F", 4, 4, F_35),
        ("F", 4, 4 + h, F_35), ("F", 4, 5, F_35), ("F", 5, 5, F_35),
    ]
    out = []
    for name, k, x, ref in cases:
        val = getattr(sf, name)(x, k=k)
        r = ref(fmpq(x))
        ok = bool(val.overlaps(r) and val.rad() < arb("1e-30") and r.rad() < arb("1e-30"))
        out.append({"function": name, "piece": k, "argument": str(fmpq(x)), "pass": ok,
                    "model": val.str(25), "closed_form": r.str(25)})
    return out


def segments(lo, hi):
    cuts = set()
    for m in range(int(lo.floor()) + 1, int(hi.ceil())):
        cuts.add(fmpq(m))
    for m in range(int((L * lo).floor()) + 1, int((L * hi).ceil())):
        cuts.add(fmpq(m) / L)
    nodes = [lo] + sorted(c for c in cuts if lo < c < hi) + [hi]
    return list(zip(nodes[:-1], nodes[1:]))


def first_kill_integral(sf, log):
    """Enclosure of int_{1/b}^{1/a} omega(beta - 1) F(L beta) beta dbeta."""
    total = arb(0)
    segs = segments(1 / B, 1 / A)
    rows = []
    half = fmpq(1, 2)
    for b0, b1 in segs:
        mid = (b0 + b1) / 2
        om = sf.omega_pieces[int((mid - 1).floor())]
        Fp = sf.F_pieces[int((L * mid).floor())]
        shift = 1 + om.center
        offsets = [b0 - shift, b1 - shift, L * b0 - Fp.center, L * b1 - Fp.center]
        if not all(-half <= t <= half for t in offsets):
            raise RuntimeError("segment leaves the validated interval of its pieces")
        d = L * shift - Fp.center
        poly = om.poly * Fp.poly(arb_poly([arb(d), arb(L)])) * arb_poly([arb(shift), 1])
        P = poly.integral()
        main = P(arb(b1 - shift)) - P(arb(b0 - shift))
        err = (arb(om.tail) * Fp.sup_poly() + om.sup_poly() * arb(Fp.tail)
               + arb(om.tail) * arb(Fp.tail)) * arb((b1 * b1 - b0 * b0) / 2)
        piece = main + arb(0, err.upper())
        total += piece
        rows.append({"beta": [str(b0), str(b1)], "omega_piece": om.k, "F_piece": Fp.k,
                     "remainder_bound": float(err.upper()), "value": piece.str(20)})
        log(f"  beta in [{float(b0):.9f}, {float(b1):.9f}]  omega piece {om.k:2d}  F piece {Fp.k}  "
            f"value {piece.str(20)}  remainder <= {float(err.upper()):.2e}")
    return total, rows


def float_cross_check(H=20000):
    """Independent double-precision solution of the delay equations (non-rigorous)."""
    eg2 = 2 * np.exp(np.euler_gamma)

    def advance(values, grid, k, same, driver):
        i0 = (k - 1) * H
        prev = driver[i0 - H:i0 + 1]
        cum = np.concatenate(([0.0], np.cumsum((prev[1:] + prev[:-1]) / 2) / H))
        values[i0:i0 + H + 1] = (k * same[i0] + cum) / grid[i0:i0 + H + 1]

    u = 1 + np.arange(12 * H + 1) / H
    om = np.zeros_like(u)
    om[:H + 1] = 1 / u[:H + 1]
    for k in range(2, 13):
        advance(om, u, k, om, om)
    s = 1 + np.arange(6 * H + 1) / H
    F = np.where(s <= 3, eg2 / s, 0.0)
    f = np.zeros_like(s)
    advance(f, s, 2, f, F)
    for k in range(3, 7):
        advance(F, s, k, F, f)
        advance(f, s, k, f, F)
    total = 0.0
    for b0, b1 in segments(1 / B, 1 / A):
        beta = np.linspace(float(b0), float(b1), 200001)
        g = np.interp(beta - 1, u, om) * np.interp(float(L) * beta, s, F) * beta
        total += float(np.sum((g[1:] + g[:-1]) / 2 * np.diff(beta)))
    fL = float(np.interp(float(L / A), s, f))
    om_a = float(om[int(11.5 * H)])
    E = 2 * float(A) ** 2 / om_a * total
    return {"f(L/a)": fL, "omega(1/a)": om_a, "E(a,b,eps)": E, "margin": fL - E}


def main():
    t0 = time.time()
    log = Log(LOG)
    sfm.set_precision()
    log("First-kill certificate for Theorem thm:two-sided-0548 (eq. sp:eq:certified-margin)")
    log(f"started {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    env = {"python": sys.version.split()[0], "python-flint": flint.__version__,
           "FLINT": flint.__FLINT_VERSION__, "numpy": np.__version__, "platform": platform.platform()}
    sources = {p.name: sha256(p) for p in (HERE / "sieve_functions.py", Path(__file__).resolve())}
    log(f"environment {env}")
    log(f"manuscript sha256 {sha256(MANUSCRIPT)}")
    for k, v in sources.items():
        log(f"source {k} sha256 {v}")
    log(f"Arb precision {ctx.prec} bits, Taylor degree {sfm.DEG}, Cauchy radius {sfm.R}, "
        f"evaluation radius {sfm.RHO}")
    log(f"parameters a={A} b={B} eps={EPS} L={L} L/a={L / A} 2b={2 * B}")

    hyp = {
        "0<eps<1/6": bool(0 < EPS < fmpq(1, 6)),
        "0<a<b<1/3": bool(0 < A < B < fmpq(1, 3)),
        "a<=(1/2-eps)/2": bool(A <= L / 2),
        "a<1/4": bool(A < fmpq(1, 4)),
        "X^b=q^0.548": bool(2 * B == fmpq(548, 1000)),
        "L/a=6.24875": bool(L / A == fmpq(624875, 100000)),
    }
    log(f"hypotheses of Theorem sp:thm:criterion: {hyp}")

    sf = sfm.SieveFunctions(omega_pieces=12, linear_pieces=6)
    log(f"Taylor models: omega pieces 1..12, F pieces 1..6, f pieces 2..6 ({time.time() - t0:.1f}s)")
    for name, pieces in (("omega", sf.omega_pieces), ("F", sf.F_pieces), ("f", sf.f_pieces)):
        log(f"  {name} max truncation tail {max(float(p.tail) for p in pieces.values()):.2e}")

    tests = closed_form_tests(sf)
    for t in tests:
        log(f"  closed form {t['function']}({t['argument']}) piece {t['piece']}: "
            f"{t['model']}  {'OK' if t['pass'] else 'FAIL'}")
    tests_ok = all(t["pass"] for t in tests)

    fL = sf.f(L / A)
    om_a = sf.omega(1 / A)
    log(f"f(L/a)      = {fL.str(30)}")
    log(f"omega(1/a)  = {om_a.str(30)}")
    log("integral over beta = 1/alpha:")
    I, rows = first_kill_integral(sf, log)
    E = 2 * arb(A) ** 2 / om_a * I
    margin = fL - E
    log(f"segments    = {len(rows)}")
    log(f"integral    = {I.str(30)}")
    log(f"E(a,b,eps)  = {E.str(30)}")
    log(f"margin      = {margin.str(30)}")

    values = {"f(L/a)": fL, "E(a,b,eps)": E, "f(L/a)-E(a,b,eps)": margin}
    checks = []
    for name, lo, hi in REQUIRED:
        v = values[name]
        ok = bool(v > arb(lo)) if hi is None else bool(arb(lo) < v < arb(hi))
        stmt = f"{name} > {decimal(lo, 8, False)}" if hi is None else \
            f"{decimal(lo, 8, False)} < {name} < {decimal(hi, 8, True)}"
        checks.append({"statement": stmt, "pass": ok})
        log(f"required {stmt}: {'PASS' if ok else 'FAIL'}")

    fc = float_cross_check()
    rigorous = {"f(L/a)": fL, "omega(1/a)": om_a, "E(a,b,eps)": E, "margin": margin}
    deviation = {k: abs(fc[k] - float(rigorous[k].mid())) for k in fc}
    fc_ok = max(deviation.values()) < 1e-7
    log(f"double-precision cross-check (independent solver, non-rigorous): {fc}")
    log(f"  max deviation from enclosure midpoints {max(deviation.values()):.2e}: "
        f"{'OK' if fc_ok else 'FAIL'}")

    status = "PASS" if (all(hyp.values()) and tests_ok and fc_ok and len(rows) == 15
                        and all(c["pass"] for c in checks)) else "FAIL"
    cert = {
        "certificate": "first_kill_0548",
        "manuscript_claims": ["Theorem thm:two-sided-0548", "Theorem sp:thm:criterion",
                              "eq. sp:eq:certified-margin (Appendix A, Section sp:sec:numerics)"],
        "manuscript": {"file": MANUSCRIPT.name, "sha256": sha256(MANUSCRIPT)},
        "sources": sources,
        "environment": env,
        "parameters": {"a": str(A), "b": str(B), "eps": str(EPS), "L": str(L), "L/a": str(L / A)},
        "hypotheses": hyp,
        "method": {
            "arithmetic": "Arb ball arithmetic (python-flint)",
            "precision_bits": ctx.prec,
            "taylor_degree": sfm.DEG,
            "taylor_centres": "k+1/2 on each [k,k+1]",
            "cauchy_radius": str(sfm.R),
            "evaluation_radius": str(sfm.RHO),
            "integration": "exact polynomial integration in beta=1/alpha with explicit remainder bounds",
            "segments": len(rows),
        },
        "segments": rows,
        "enclosures": {"f(L/a)": enclosure(fL), "omega(1/a)": enclosure(om_a),
                       "integral": enclosure(I), "E(a,b,eps)": enclosure(E),
                       "f(L/a)-E(a,b,eps)": enclosure(margin)},
        "required_bounds": checks,
        "validation": {"closed_form_tests": tests,
                       "double_precision_cross_check": {"values": fc, "abs_deviation": deviation,
                                                        "pass": fc_ok}},
        "scope": "Numerical inequality only; the analytic reduction in Appendix A is not certified here.",
        "status": status,
    }
    CERTIFICATE.parent.mkdir(parents=True, exist_ok=True)
    CERTIFICATE.write_text(json.dumps(cert, indent=2) + "\n", encoding="utf-8", newline="\n")
    log(f"certificate written: {CERTIFICATE.relative_to(ROOT).as_posix()}")
    log(f"STATUS {status}")
    log(f"elapsed {time.time() - t0:.1f}s")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
