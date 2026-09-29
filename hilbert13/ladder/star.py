"""Isotypic components of the divisor behind condition (*):  Phi_g - 3 Phi_c  on C."""
from cusp import *
tau = a
mu = cyc((0, 2), (1, 3)); tmu = mul(tau, mu)
goods = [cyc((0, 1), (i, j)) for (i, j) in [(4, 5), (4, 6), (5, 6)]]
goods += [mul(tau, g) for g in goods]
Phi_c = fix_divisor(mu) + fix_divisor(tmu)
Phi_g = sum(fix_divisor(g) for g in goods)
print("deg Phi_c", Phi_c.sum(), " deg Phi_g", Phi_g.sum())
v = Phi_g - 3 * Phi_c
print("Phi_g - 3Phi_c :", {k: round(x, 4) for k, x in isotypic_norms(v).items()})
# bad-involution fixed points t_b and the free orbit t_0, for later
