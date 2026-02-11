// This file certifies that the sextic K3 surface given by f2, f3 below has geometric picard rank 2 over F_47, and also checks that the chosen lift contains no lines over Qbar.

p := 47;

P4<x0, x1, x2, x3, x4> := ProjectiveSpace(FiniteField(p), 4);

f2_string := "x0^2 + 44*x0*x1 + 3*x1^2 + 5*x0*x2 + 4*x1*x2 + 5*x2^2 + 46*x0*x3 + 45*x1*x3
+ 44*x2*x3 + 42*x0*x4 + 5*x1*x4";
f3_string := "2*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3 + 46*x0*x1*x2 + 44*x1^2*x2 +
4*x0*x2^2 + 43*x1*x2^2 + 5*x2^3 + 4*x0^2*x3 + x0*x1*x3 + 5*x1^2*x3 + 4*x0*x2*x3
+ 4*x1*x2*x3 + 44*x2^2*x3 + 4*x1*x3^2 + 46*x2*x3^2 + 5*x0^2*x4 + 43*x0*x1*x4 +
2*x1^2*x4 + x0*x2*x4 + 4*x1*x2*x4 + 45*x2^2*x4 + 4*x0*x3*x4 + 44*x2*x3*x4 +
46*x0*x4^2 + 46*x1*x4^2 + 5*x2*x4^2";

f2 := eval f2_string;
f3 := eval f3_string;

Xp := Scheme(P4, [f2, f3]);

printf "X is smooth over F_%o: %o\n", p, IsNonsingular(Xp);

// COMPUTE f6

M := JacobianMatrix([f2, f3]);
MSub := Submatrix(M, 1, 4, 2, 2); // Submatrix of determinant computing derivatives with respect to x3, x4 only.

q := Determinant(Submatrix(M, 1, 4, 2, 2)); // ramification locus of rational map P4 - - - > P2
Z := Scheme(P4, [f2, f3, q]); // Ramification locus of restriction X --> P2

P2<y0, y1, y2> := ProjectiveSpace(FiniteField(p), 2);
f := map<P4 -> P2 | [x0, x1, x2] >; // Projection map P4 -> P2

f6 := Generators(Ideal(f(Z)))[1]; // Generator of branch locus

printf "f6 defining branch locus (a smooth sextic curve): %o\n", f6;

printf "sextic curve defined by f6 is smooth (over F_%o): %o\n", p, IsNonsingular(f(Z));

// COMPUTE WEIL POLYNOMIAL WITH RESPECT TO f6

printf "computing Weil polynomial using degree 2 model...\n";

wp1, wp2 := WeilPolynomialOfDegree2K3Surface(f6);

// COMPUTE WEIL POLYNOMIAL OF TWISTED OPTION

printf "computing Weil polynomial using degree 2 model with a quadratic twist...\n";

wpt1, wpt2 := WeilPolynomialOfDegree2K3Surface(-1*f6);

R<t> := PolynomialRing(Integers());

wp1 := Evaluate(wp1, t);
wp2 := Evaluate(wp2, t);
wpt1 := Evaluate(wpt1, t);
wpt2 := Evaluate(wpt2, t);

printf "Weil polynomial without twist: %o\n", Factorization(wp1 * wp2);
printf "Weil polynomial with twist: %o\n", Factorization(wpt1 * wpt2);

// LIFT Xp TO X OVER Q

P4Q<x0, x1, x2, x3, x4> := ProjectiveSpace(Rationals(), 4);

f2_Q := eval f2_string;
f3_Q := eval f3_string;

f2_L := f2_Q + p*x3^2 + p * x4^2;
f3_L := f3_Q;

printf "f2_L: %o\n", f2_L;
printf "f3_L: %o\n", f3_L;


X := Scheme(P4Q, [f2_L, f3_L]);

printf "X is smooth over Q: %o\n", IsNonsingular(X);

C := CoordinateRing(P4Q);

R_<y1, y2, y3, y4, y5, y6> := PolynomialRing(Rationals(), 6);
S<u, v> := PolynomialRing(R_, 2);

function getLineScheme(v1, v2)
    hom_cs := hom<C -> S | [u * v1[i] + v * v2[i]: i in [1..5]]>;
    coeffs := Coefficients(hom_cs(f2_L)) cat Coefficients(hom_cs(f3_L));
    return IsProper(ideal<R_ | coeffs>);
end function;

print "Check properness of ideal for each schubert cell (should be false in each case)";

getLineScheme(
    [1, 0, y1, y2, y3],
    [0, 1, y4, y5, y6]);
getLineScheme(
    [1, y1, 0, y2, y3],
    [0, y4, 1, y5, y6]);
getLineScheme(
    [1, y1, y2, 0, y3],
    [0, y4, y5, 1, y6]);
getLineScheme(
    [1, y1, y2, y3, 0],
    [0, y4, y5, y6, 1]);
getLineScheme(
    [0, 1, 0, y1, y2],
    [0, 0, 1, y4, y5]);
getLineScheme(
    [0, 1, y1, 0, y2],
    [0, 0, y4, 1, y5]);
getLineScheme(
    [0, 1, y1, y2, 0],
    [0, 0, y4, y5, 1]);
getLineScheme(
    [0, 0, 1, 0, y1],
    [0, 0, 0, 1, y4]);
getLineScheme(
    [0, 0, 1, y1, 0],
    [0, 0, 0, y4, 1]);
getLineScheme(
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1]);
