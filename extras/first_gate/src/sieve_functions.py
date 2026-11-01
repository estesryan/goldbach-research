"""Validated enclosures of Buchstab's function omega and the linear-sieve functions F, f.

    omega(u) = 1/u            (1 <= u <= 2),   (u omega(u))' = omega(u - 1)   (u > 2),
    F(s)     = 2 e^gamma / s  (0 < s <= 3),    (s F(s))'     = f(s - 1)       (s > 3),
    f(s)     = 0              (0 < s <= 2),    (s f(s))'     = F(s - 1)       (s > 2).

On each unit interval [k, k + 1] a function is represented by its Taylor expansion
about c = k + 1/2 in the offset x = u - c.  The coefficients of degree <= DEG are Arb
enclosures of the exact Taylor coefficients; the truncation tail is bounded by
Cauchy's estimate on the disc |x| <= R.
"""
from flint import arb, arb_poly, ctx, fmpq

PREC = 192
DEG = 110
R = fmpq(5, 4)
RHO = fmpq(1, 2)


def set_precision():
    ctx.prec = PREC


class Piece:
    """Taylor model of a function on [k, k + 1] about k + 1/2, valid for |x| <= RHO."""

    def __init__(self, k, poly, bound):
        self.k = k
        self.center = fmpq(2 * k + 1, 2)
        self.poly = poly
        self.bound = bound
        ratio = arb(RHO / R)
        self.tail = (bound * ratio ** (DEG + 1) / (1 - ratio)).upper()

    def __call__(self, x):
        x = arb(x)
        if not abs(x).upper() <= arb(RHO):
            raise ValueError("offset outside the validated interval")
        return self.poly(x) + arb(0, self.tail)

    def sup_poly(self):
        """Upper bound for |poly(x)| on |x| <= RHO."""
        s = arb(0)
        for n, c in enumerate(self.poly.coeffs()):
            s += abs(c).upper() * arb(RHO) ** n
        return s.upper()


def _reciprocal(c):
    return arb_poly([arb(-1) ** n / arb(c) ** (n + 1) for n in range(DEG + 1)])


def _truncate(p):
    return arb_poly(p.coeffs()[: DEG + 1])


def _reciprocal_piece(k, scale):
    c = fmpq(2 * k + 1, 2)
    return Piece(k, _reciprocal(c) * scale, (scale / arb(c - R)).upper())


def _log_piece(scale):
    """scale * log(s - 1) / s on [2, 3]."""
    c = fmpq(5, 2)
    a = arb(c - 1)
    log_series = arb_poly([a.log()] + [arb(-1) ** (n + 1) / (n * a ** n) for n in range(1, DEG + 1)])
    # On |x| <= 5/4: 1/4 <= |3/2 + x| <= 11/4 and Re(3/2 + x) > 0, so |log(3/2 + x)| <= log 4 + pi/2.
    bound = scale * (arb(4).log() + arb.pi() / 2) / arb(c - R)
    return Piece(2, _truncate(log_series * _reciprocal(c)) * scale, bound.upper())


def _delay_step(k, same, driver):
    """Piece k of h with (u h(u))' = d(u - 1), from piece k - 1 of h and of d.

    G(x) = (c + x) h(c + x) = k h(k) + int_{-1/2}^{x} d_{k-1}(t) dt, where the
    driver d_{k-1} uses the same offset about c - 1.
    """
    c = fmpq(2 * k + 1, 2)
    left = k * same(RHO)
    antider = driver.poly.integral()
    g0 = left - antider(-arb(RHO)) + arb(0, (arb(RHO) * driver.tail).upper())
    g = arb_poly([g0] + antider.coeffs()[1: DEG + 1])
    bound = (abs(left).upper() + arb(R + RHO) * driver.bound) / arb(c - R)
    return Piece(k, _truncate(g * _reciprocal(c)), bound.upper())


class SieveFunctions:
    def __init__(self, omega_pieces, linear_pieces):
        two_eg = 2 * arb.const_euler().exp()
        self.omega_pieces = {1: _reciprocal_piece(1, arb(1))}
        for k in range(2, omega_pieces + 1):
            self.omega_pieces[k] = _delay_step(k, self.omega_pieces[k - 1], self.omega_pieces[k - 1])
        self.F_pieces = {1: _reciprocal_piece(1, two_eg), 2: _reciprocal_piece(2, two_eg)}
        self.f_pieces = {2: _log_piece(two_eg)}
        for k in range(3, linear_pieces + 1):
            self.F_pieces[k] = _delay_step(k, self.F_pieces[k - 1], self.f_pieces[k - 1])
            self.f_pieces[k] = _delay_step(k, self.f_pieces[k - 1], self.F_pieces[k - 1])

    @staticmethod
    def _eval(pieces, u, k):
        u = fmpq(u)
        if k is None:
            k = int(u.floor())
        return pieces[k](arb(u - pieces[k].center))

    def omega(self, u, k=None):
        return self._eval(self.omega_pieces, u, k)

    def F(self, s, k=None):
        return self._eval(self.F_pieces, s, k)

    def f(self, s, k=None):
        return self._eval(self.f_pieces, s, k)
