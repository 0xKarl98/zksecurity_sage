"""
Basic setup for sumcheck examples.
Field: GF(101)
Ring: polynomial ring in three variables x1, x2, x3 over GF(101)
"""

F = GF(101)
R.<x1, x2, x3> = PolynomialRing(F)

# Example exercise: define a 2-variable polynomial (ignoring x3), e.g. f = x1*x2 + x1 + 1
f = x1 * x2 + x1 + 1
print("Field F =", F)
print("Ring R =", R)
print("Sample polynomial f =", f)

# Step 2: sum over the Boolean hypercube {0,1}^n
from itertools import product

def hypercube_eval(poly, vars_tuple, mult=False, defaults=None):
    """
    Evaluate 'poly' on all points in {0,1}^n for the given vars,
    and either sum (default) or multiply (if mult=True) the results.
    'vars_tuple' selects which variables are iterated; others can be fixed via defaults.
    """
    base_subs = {} if defaults is None else defaults
    total = F(1) if mult else F(0)
    for bits in product([0, 1], repeat=len(vars_tuple)):
        subs = dict(base_subs)
        subs.update({v: F(b) for v, b in zip(vars_tuple, bits)})
        val = poly.subs(subs)
        total = total * val if mult else total + val
    return total

# demo: sum f(x1,x2) over {0,1}^2 (ignoring x3)
print("Sum_{(x1,x2) in {0,1}^2} f(x1,x2) =", hypercube_eval(f, (x1, x2), defaults={x3: F(0)}))
# demo: multiply f(x1,x2) over {0,1}^2 (ignoring x3)
print("Prod_{(x1,x2) in {0,1}^2} f(x1,x2) =", hypercube_eval(f, (x1, x2), mult=True, defaults={x3: F(0)}))

# Step 3: minimal prover/verifier for sumcheck (assuming multilinear f)
class Prover:
    def __init__(self, poly, vars_tuple):
        self.f = poly
        self.vars = vars_tuple

    def total_sum(self):
        return hypercube_eval(self.f, self.vars)

    def send_univariate(self, i, prefix_vals):
        """
        i is 1-based index of current variable.
        prefix_vals is a list of r_1,...,r_{i-1}.
        Returns univariate S_i(X_i) as a polynomial in R.
        """
        assert i >= 1 and i <= len(self.vars)
        assert len(prefix_vals) == i - 1
        Xi = self.vars[i - 1]
        suffix_vars = self.vars[i:]  # vars after Xi
        defaults = {v: val for v, val in zip(self.vars[: i - 1], prefix_vals)}
        # degree-1, so interpolate from Xi=0 and Xi=1
        s0 = hypercube_eval(self.f, suffix_vars, defaults={**defaults, Xi: F(0)})
        s1 = hypercube_eval(self.f, suffix_vars, defaults={**defaults, Xi: F(1)})
        return s0 + (s1 - s0) * Xi


class Verifier:
    def __init__(self, poly, vars_tuple):
        self.f = poly
        self.vars = vars_tuple
        self.r_values = []

    def run(self, prover):
        ell = len(self.vars)
        claim = prover.total_sum()
        print(f"\nProver claims total sum C = {claim}")
        prev = claim
        for i in range(1, ell + 1):
            S_i = prover.send_univariate(i, self.r_values)
            Xi = self.vars[i - 1]
            val0 = S_i.subs({Xi: F(0)})
            val1 = S_i.subs({Xi: F(1)})
            if val0 + val1 != prev:
                print(f"Round {i}: reject (consistency check failed)")
                return False
            r_i = F.random_element()
            self.r_values.append(r_i)
            prev = S_i.subs({Xi: r_i})
            print(f"Round {i}: received S_{i}(X), sampled r_{i}={r_i}, S_{i}(r_{i})={prev}")
        final_val = self.f(*self.r_values)
        if final_val == prev:
            print(f"Final check: g(r_1,...,r_ell)={final_val} matches S_ell(r_ell); accept")
            return True
        print(f"Final check failed: g(...)={final_val} vs {prev}; reject")
        return False


class SumcheckVerifier:
    """
    Modular verifier à la the solution snippet:
    - check_round_poly stores the poly and checks the 0/1 sum
    - update_claimed_sum applies challenge r_i and stores it
    - final_check compares against f at all sampled r's
    """
    def __init__(self, claimed_sum, poly, vars_tuple):
        self.claimed_sum = claimed_sum
        self.f = poly
        self.vars = vars_tuple
        self.r_values = []
        self.poly = None

    def check_round_poly(self, poly):
        Xi = self.vars[len(self.r_values)]
        self.poly = poly
        val0 = poly.subs({Xi: F(0)})
        val1 = poly.subs({Xi: F(1)})
        return val0 + val1 == self.claimed_sum

    def update_claimed_sum(self, r):
        Xi = self.vars[len(self.r_values)]
        self.claimed_sum = self.poly.subs({Xi: r})
        self.r_values.append(r)

    def final_check(self):
        return self.claimed_sum == self.f(*self.r_values)


def run_sumcheck(prover, claimed_sum=None):
    claim = prover.total_sum() if claimed_sum is None else claimed_sum
    verifier = SumcheckVerifier(claim, prover.f, prover.vars)
    print(f"\nProver claims total sum C = {claim}")
    for i in range(len(prover.vars)):
        S_i = prover.send_univariate(i + 1, verifier.r_values)
        if not verifier.check_round_poly(S_i):
            print(f"Round {i+1}: reject (consistency check failed)")
            return False
        r_i = F.random_element()
        verifier.update_claimed_sum(r_i)
        print(f"Round {i+1}: received S_{i+1}(X), sampled r_{i+1}={r_i}, S_{i+1}(r_{i+1})={verifier.claimed_sum}")
    if verifier.final_check():
        print(f"Final check: g(r_1,...,r_ell)={verifier.claimed_sum}; accept")
        return True
    print("Final check failed; reject")
    return False


if __name__ == "__main__":
    prover = Prover(f, (x1, x2, x3))
    # run with modular verifier flavor (correct claim)
    run_sumcheck(prover)
    # run again with an intentionally wrong claim to see rejection
    run_sumcheck(prover, claimed_sum=prover.total_sum() + 1)
