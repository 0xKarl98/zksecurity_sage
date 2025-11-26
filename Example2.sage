# Define a variable in Sage
a = var('a')  # Creates a symbolic variable and assigns it to 'a'

# Define a polynomial ring
R.<x> = PolynomialRing(QQ)  # Creates a polynomial ring R in variable x over the rational numbers QQ

# Define a polynomial
f = x^2 + 3*x + 2  # Creates a polynomial f in the ring R

# Evaluate the polynomial at a specific value
result = f(1)  # Evaluates f at x=1
print(result)  # Output: 6