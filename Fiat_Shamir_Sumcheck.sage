"""
Fiat–Shamir variant of sumcheck over GF(101).
- Polynomial: g(x1, x2) = x1*x2 + x1 + x2 (multilinear)
- Challenges r_i are derived deterministically from the transcript via SHA-256,
  so anyone can reproduce the verifier's randomness.
"""

import hashlib
from itertools import product

# Field and ring
F = GF(101)
R.<x1, x2> = PolynomialRing(F)

# Example polynomial (already multilinear)
g = x1 * x2 + x1 + x2


def hypercube_sum(poly, vars_tuple):
    """Sum poly over all points in {0,1}^n for the given vars."""
    total = F(0)
    for bits in product([0, 1], repeat=len(vars_tuple)):
        subs = {v: F(b) for v, b in zip(vars_tuple, bits)}
        total += poly.subs(subs)
    return total


def hash_to_field(data_bytes):
    """Hash bytes to a field element in F using SHA-256, then reduce mod |F|."""
    h = hashlib.sha256(data_bytes).digest()
    return F(Integer(int.from_bytes(h, "big") % F.order()))


class FSProver:
    def __init__(self, poly, vars_tuple):
        self.f = poly
        self.vars = vars_tuple

    def total_sum(self):
        return hypercube_sum(self.f, self.vars)

    def send_univariate(self, i, prefix_vals):
        """
        Return S_i(X_i) = sum_{suffix in {0,1}^{ell-i}} f(prefix, X_i, suffix).
        i is 1-based; prefix_vals = [r1,...,r_{i-1}].
        """
        Xi = self.vars[i - 1]
        suffix_vars = self.vars[i:]
        defaults = {v: val for v, val in zip(self.vars[: i - 1], prefix_vals)}
        # interpolate degree-1 polynomial from Xi=0,1
        s0 = hypercube_sum(self.f.subs(defaults).subs({Xi: F(0)}), suffix_vars)
        s1 = hypercube_sum(self.f.subs(defaults).subs({Xi: F(1)}), suffix_vars)
        return s0 + (s1 - s0) * Xi


class FSVerifier:
    def __init__(self, poly, vars_tuple):
        self.f = poly
        self.vars = vars_tuple
        self.r_values = []
        self.claimed = None

    def derive_challenge(self, S_i, round_idx):
        """Derive r_i from transcript (round, poly repr, prior r's, claimed)."""
        transcript = f"{round_idx}|{S_i}|{self.r_values}|{self.claimed}".encode()
        return hash_to_field(transcript)

    def run(self, prover, claimed_sum=None):
        ell = len(self.vars)
        self.claimed = prover.total_sum() if claimed_sum is None else claimed_sum
        print(f"\nFiat–Shamir sumcheck claim: C = {self.claimed}")
        for i in range(1, ell + 1):
            S_i = prover.send_univariate(i, self.r_values)
            Xi = self.vars[i - 1]
            val0 = S_i.subs({Xi: F(0)})
            val1 = S_i.subs({Xi: F(1)})
            if val0 + val1 != self.claimed:
                print(f"Round {i}: reject (S_i(0)+S_i(1) != claimed)")
                return False
            r_i = self.derive_challenge(S_i, i)
            self.r_values.append(r_i)
            self.claimed = S_i.subs({Xi: r_i})
            print(f"Round {i}: r_{i}={r_i}, S_{i}(r_{i})={self.claimed}")
        final_val = self.f(*self.r_values)
        if final_val == self.claimed:
            print(f"Final check ok: g(r1,r2)={final_val}")
            return True
        print(f"Final check failed: g(...)={final_val} vs {self.claimed}")
        return False


if __name__ == "__main__":
    prover = FSProver(g, (x1, x2))
    verifier = FSVerifier(g, (x1, x2))
    verifier.run(prover)
