# First-Gate computational supplement

Computational evidence for *Goldbach Representations in Consecutive Prime-Square Blocks*
(`first_gate_prime_square_blocks.tex`, SHA-256
`8ca6fde20856f0b86e74303951eeb5e4d3d0a90a6b1bedaa8ffc65b95b4ab9b8`).

| Manuscript claim | Script | Certificate | Log |
|---|---|---|---|
| Theorem `thm:two-sided-0548`, Appendix A, inequality `sp:eq:certified-margin` | `src/certify_first_kill.py` (uses `src/sieve_functions.py`) | `certificates/first_kill_0548.json` | `logs/first_kill_0548.log` |
| Computation `comp:finite` | `src/verify_first_gate_scan.py` | `certificates/first_gate_scan_q9973.json` | `logs/first_gate_scan_q9973.log` |

## Requirements

Python 3.12 and the pinned packages in `requirements.txt` (python-flint 0.9.0 bundles FLINT 3.6.0).
The recorded runs used Python 3.12.10 on Windows 11 (x86-64).

```sh
python -m pip install -r requirements.txt
```

## Reproduction

```sh
cd src
python certify_first_kill.py
python verify_first_gate_scan.py
cd ..
sha256sum -c SHA256SUMS.txt
```

Each script regenerates its certificate and log and exits with status 0 only if every check passes.
The certificates contain no timestamps, so a rerun in the same environment reproduces them byte for byte.

## First-kill inequality

With exact rationals a = 2/25, b = 137/500, ε = 1/10⁴ and L = 1/2 − ε, the script encloses
f(L/a), ω(1/a) and E(a, b, ε) = 2a²/ω(1/a) ∫ₐᵇ ω((1−α)/α) F(L/α) α⁻³ dα, and checks the three stated bounds.

- On each [k, k+1], ω, F and f are degree-110 Taylor models about k + 1/2, computed in 192-bit Arb
  arithmetic from the delay equations. Their coefficients enclose the exact Taylor coefficients.
  Every piece is holomorphic on |x| < 3/2, and truncation is bounded by Cauchy's estimate on |x| ≤ 5/4.
- With β = 1/α the integral equals ∫ ω(β−1) F(Lβ) β dβ over [1/b, 1/a]. At the 14 integer crossings of
  β − 1 and Lβ it splits into 15 segments. On each segment the integrand is a polynomial, which is
  integrated exactly, plus an explicitly bounded remainder.
- The models are validated against closed forms on [2, 5]. An independent double-precision solver
  agrees to 1.1·10⁻¹⁰; this check is not rigorous.

Certified: f(L/a) − E(a, b, ε) ∈ [0.000420563830407454, 0.000420563830407455].

## Finite First-Gate scan

All 1,225 blocks q² < N < r² with 11 ≤ q ≤ 9973 (r = 10007 for the last block), covering every even N with
122 ≤ N < 10007². That is 50,069,964 targets.

- For N = 2C and W = q² − C, an offset d survives when gcd(C² − d², ∏_{p<q} p) = 1. When d < W,
  both C ± d lie in [2, q²), so this is equivalent to C ± d being primes ≥ q. Therefore First-Gate holds at
  N (Proposition `prop:full-FG-window`) if and only if some d < W survives.
- The least surviving offset is found from a sieve of Eratosthenes to 10⁸. The sieve is validated
  against π(10ᵏ) for k = 1…8 and against FLINT primality tests on [0, 2·10⁵] and [10⁸ − 2·10⁵, 10⁸].
- In every block, the least survivor of the first target and of the hardest target is recomputed from
  the gcd definition with exact integers: 2,450 cross-checks, with no prime table.

Result: no First-Gate failures and no cross-check mismatches. The largest least survivor is 3429, at N = 92,661,736.

## Scope

The certificates establish numerical statements only. Theorem `thm:two-sided-0548` also depends on the
analytic argument of Appendix A: Jurkat–Richert, Bombieri–Vinogradov and Siegel–Walfisz,
the bilinear Lemma `sp:lem:BBV`, and the Buchstab asymptotics. None of that is verified by computation here.

## Citation

This paper is an unrefereed research manuscript. Suggested citation:

> Ryan S. Estes, *Goldbach Representations in Consecutive Prime-Square Blocks: A Multidimensional
> Factorization and Sieve Analysis*, unrefereed research manuscript,
> https://github.com/estesryan/goldbach-visualizations.

The manuscript is [`first_gate_prime_square_blocks.tex`](first_gate_prime_square_blocks.tex)
([PDF](first_gate_prime_square_blocks.pdf)). Citation metadata is in [`CITATION.cff`](CITATION.cff).
