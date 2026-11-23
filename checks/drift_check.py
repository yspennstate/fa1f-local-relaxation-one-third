"""Exact check of the drift of W = D + u b for the marked infection (Lemma 3.2 and Remark 6.1).

The mark sits at y, with y^+ the neighbour towards x and y^- the other one. The rules of Definition 3.1 are applied
transition by transition to every local configuration of y^{--}, y^-, y, y^+ (each infected, healthy, or a frozen
infected exterior site where the geometry allows it), with D >= 2. The script prints, for each q = j/300 with
100 < j < 300, the maximum drift over configurations and checks that it equals -gamma = -q(3q-1)/(1+q) for
u = p/(1+q). For the remark it checks, on a grid of u in [0, 1], that the maximum drift is >= 0 for every u when
q <= 1/3, and that the two drifts there (the Case 1 and Case 2 values) are equal at u = p/(1+q).

Run: python drift_check.py   (exact rational arithmetic, about six minutes)
"""
from fractions import Fraction as Fr
from itertools import product

INF, HEA, FRZ = 'inf', 'hea', 'frz'   # infected in the box, healthy, frozen infected exterior site


def infected(s):
    return s in (INF, FRZ)


def b_of(mark_frozen, rear):
    """b = 1 if the mark is in the box and its rear neighbour is infected (frozen or not)."""
    return 0 if mark_frozen else int(infected(rear))


def configs():
    """Yield (mark_frozen, s_mm, s_m, s_p): states of y^{--}, y^-, y^+ (y itself is infected)."""
    # mark frozen: y is an exterior site, nothing behind it, y^+ is in the box
    for s_p in (INF, HEA):
        yield True, None, None, s_p
    # rear neighbour frozen: nothing behind it
    for s_p in (INF, HEA):
        yield False, None, FRZ, s_p
    # rear neighbour in the box: y^{--} may be in the box or frozen
    for s_mm, s_m, s_p in product((INF, HEA, FRZ), (INF, HEA), (INF, HEA)):
        yield False, s_mm, s_m, s_p


def drift(q, u, cfg):
    p = 1 - q
    mark_frozen, s_mm, s_m, s_p = cfg
    W0 = u * b_of(mark_frozen, s_m)          # D is measured relative to the current D
    total = Fr(0)
    # rule (1): clock of y^+ rings with spin 0 (rate q); mark moves to y^+, its rear neighbour is y (infected)
    total += q * ((-1 + u * 1) - W0)
    if mark_frozen:
        return total
    # clock of y^- (rate 1): spin 0 infects it, spin 1 heals it (allowed: y is infected); frozen never changes
    if s_m != FRZ:
        total += q * (u * b_of(False, INF) - W0)
        total += p * (u * b_of(False, HEA) - W0)
    # clock of y with spin 1 (rate p): y heals if a neighbour is infected
    if infected(s_m) or infected(s_p):
        if infected(s_p):
            new = -1 + u * 0                  # mark to y^+, rear is y, now healthy
        elif s_m == FRZ:
            new = +1 + 0                      # mark to the frozen site, b = 0
        else:
            new = +1 + u * b_of(False, s_mm)  # mark to y^-, its rear is y^{--}
        total += p * (new - W0)
    # clock of y^+ with spin 1 heals y^+, clock of y with spin 0 does nothing: no change of D or b
    return total


def max_drift(q, u):
    return max(drift(q, u, c) for c in configs())


def main():
    bad = 0
    for j in range(101, 300):
        q = Fr(j, 300)
        p = 1 - q
        u = p / (1 + q)
        gamma = q * (3 * q - 1) / (1 + q)
        m = max_drift(q, u)
        if m != -gamma:
            bad += 1
            print('MISMATCH q =', q, 'max drift', m, '-gamma', -gamma)
    print('Lemma 3.2: q = j/300, 100 < j < 300:', 'all equal -gamma' if bad == 0 else '%d mismatches' % bad)
    grid = [Fr(i, 600) for i in range(601)]
    bad2 = 0
    for j in range(1, 101):
        q = Fr(j, 300)
        if min(max_drift(q, u) for u in grid) < 0:
            bad2 += 1
            print('negative max drift at q =', q)
    print('Remark 6.1: q = j/300, j <= 100, u in {i/600}:', 'max drift >= 0 for every u' if bad2 == 0 else
          '%d failures' % bad2)
    bad3 = 0
    for j in range(101, 300):
        q = Fr(j, 300)
        p = 1 - q
        u = p / (1 + q)
        d1 = q * (-1 + 2 * u)
        d2 = -q - p * u + p
        if d1 != d2:
            bad3 += 1
    print('Remark 6.1: Case 1 and Case 2 drifts equal at u = p/(1+q):', 'yes' if bad3 == 0 else '%d failures' % bad3)
    return bad + bad2 + bad3


if __name__ == '__main__':
    raise SystemExit(1 if main() else 0)
