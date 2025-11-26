"""
MLE on 3 bits for g:{0,1}^3 -> F:
    g(0,0,0) = 0, and g(b) = 1 otherwise.
We build its multilinear extension  1 - (1 - x1)(1 - x2)(1 - x3)
and verify it matches g on the Boolean cube.
"""

F = GF(101)  # any field containing {0,1} works
R.<x1, x2, x3> = PolynomialRing(F)

# multilinear extension from the closed form in the note
g_tilde = 1 - (1 - x1) * (1 - x2) * (1 - x3)
print("g_tilde =", g_tilde)

# original Boolean function
def g(b1, b2, b3):
    return F(0) if (b1, b2, b3) == (0, 0, 0) else F(1)

# check the MLE on all Boolean inputs
print("\nCheck on {0,1}^3:")
for b1 in (0, 1):
    for b2 in (0, 1):
        for b3 in (0, 1):
            lhs = g_tilde(F(b1), F(b2), F(b3))
            rhs = g(b1, b2, b3)
            print(f"{(b1, b2, b3)} -> g_tilde={lhs}, g={rhs}")
