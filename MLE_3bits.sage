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

# evaluate g_tilde at (0.2,0.8,0.5) && (1,0,0)
print("\nEvaluate g_tilde at two points:")
xs1, xs2, xs3 = SR.var("x1 x2 x3")
g_tilde_sym = 1 - (1 - xs1) * (1 - xs2) * (1 - xs3) 
for pt in [(0.2, 0.8, 0.5), (1, 0, 0)]:
    val = g_tilde_sym(x1=pt[0], x2=pt[1], x3=pt[2])
    print(f"g_tilde{pt} = {val}")

# Schwartz-Zippel demo: take h as the difference between two degree-2 polynomials.
# Think of g_tilde as the "true" MLE and g_tilde_wrong as a (wrong) alternative
# that still has total degree <= 2. h = g_tilde - g_tilde_wrong is the
# polynomial Schwartz–Zippel reasons about (the difference between two MLEs).
print("\nSchwartz–Zippel: h = g_tilde - g_tilde_wrong")
g_tilde_wrong = g_tilde + x1 * x2  # inject a wrong term so they disagree
h = g_tilde - g_tilde_wrong
d = h.total_degree()
bound = d / F.order()
print(f"g_tilde_wrong = {g_tilde_wrong}")
print(f"h = g_tilde - g_tilde_wrong = {h}")
print(f"total degree d = {d}, bound Pr[h(r)=0] <= d/|F| = {bound}")

# empirical probability that h(r)=0 for random r in F^3 (use many samples)
trials = 20000
hits = 0
for _ in range(trials):
    r = (F.random_element(), F.random_element(), F.random_element())
    if h(*r) == 0:
        hits += 1
print(f"Empirical Pr[h(r)=0] ~ {hits/trials} over {trials} samples")
