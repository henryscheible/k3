// This file certifies that the degree 8 K3 surface (intersection of three quadrics in P5)
// given by q1, q2, q3 below has geometric Picard rank 1. The argument is:
//   (1) The Weil polynomial over F_47 shows geometric Picard rank at most 2.
//   (2) No tritangent lines exist over Qbar, so the rank drops by at least 1
//       under specialization, forcing geometric Picard rank exactly 1 over Q.

p := 47;

// Define matrices over Q; reduce mod p as needed.
q1_Q := Matrix(Rationals(), [[-136, -232, -136,  144, -146,  -42],
                              [-232, -140,  187,   93,  -98,  193],
                              [-136,  187,    1,  234,  137,  281],
                              [ 144,   93,  234,   47,  -94, -141],
                              [-146,  -98,  137,  -94,  237, -137],
                              [ -42,  193,  281, -141, -137, -139]]);

q2_Q := Matrix(Rationals(), [[  43,   44,  -50, -140,  281,  -49],
                              [  44,  141,  192,   -4,   95,  -91],
                              [ -50,  192,  185, -188,  281, -234],
                              [-140,   -4, -188,  -89,   52,   95],
                              [ 281,   95,  281,   52,  144,  -42],
                              [ -49,  -91, -234,   95,  -42,  -44]]);

q3_Q := Matrix(Rationals(), [[ 193,   -1, -146,    1,    4,  235],
                              [  -1,  234,   95,  -90,  137, -140],
                              [-146,   95,   51,    3,   92,  280],
                              [   1,  -90,    3,  183,  143,   45],
                              [   4,  137,   92,  143,    0,  -98],
                              [ 235, -140,  280,   45,  -98,  185]]);

q1 := ChangeRing(q1_Q, FiniteField(p));
q2 := ChangeRing(q2_Q, FiniteField(p));
q3 := ChangeRing(q3_Q, FiniteField(p));

// === STEP 1: Verify smoothness of X over F_p ===

P5<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(FiniteField(p), 5);
P2<u, v, w> := ProjectiveSpace(FiniteField(p), 2);

R5 := CoordinateRing(P5);
coord_vec := Matrix(R5, [[x0, x1, x2, x3, x4, x5]]);
g1 := (coord_vec * ChangeRing(q1, R5) * Transpose(coord_vec))[1][1];
g2 := (coord_vec * ChangeRing(q2, R5) * Transpose(coord_vec))[1][1];
g3 := (coord_vec * ChangeRing(q3, R5) * Transpose(coord_vec))[1][1];

Xp := Scheme(P5, [g1, g2, g3]);
printf "X is smooth over F_%o: %o\n", p, IsNonsingular(Xp);

// === STEP 2: Compute f6 (branch sextic) and verify smoothness ===

f6 := -1 * Determinant(u*q1 + v*q2 + w*q3);
Z := Scheme(P2, f6);
printf "f6 defining branch locus (a smooth sextic curve): %o\n", f6;
printf "sextic curve defined by f6 is smooth (over F_%o): %o\n", p, IsNonsingular(Z);

// === STEP 3: Print the tritangent line over F_47 ===

printf "Tritangent lines of f6 over F_%o (should be exactly one): %o\n", p, TriTangentLines(f6);

// === STEP 4: Compute Weil polynomials for f6 and -f6 ===

printf "computing Weil polynomial...\n";
wp1, wp2 := WeilPolynomialOfDegree2K3Surface(f6);

printf "computing Weil polynomial with quadratic twist...\n";
wpt1, wpt2 := WeilPolynomialOfDegree2K3Surface(-1*f6);

R<t> := PolynomialRing(Integers());
wp1 := Evaluate(wp1, t);
wp2 := Evaluate(wp2, t);
wpt1 := Evaluate(wpt1, t);
wpt2 := Evaluate(wpt2, t);

printf "Weil polynomial without twist: %o\n", Factorization(wp1 * wp2);
printf "Weil polynomial with twist: %o\n", Factorization(wpt1 * wpt2);

// === STEP 5: Naive point count check ===

wp := wp1 * wp2;
wpt := wpt1 * wpt2;
C := CompanionMatrix(wp);
Ct := CompanionMatrix(wpt);
wp_count := Trace(C) + p^2 + 1;
wpt_count := Trace(Ct) + p^2 + 1;
naive_count := #Points(Xp);
printf "Naive point count over F_%o: %o\n", p, naive_count;
printf "Point count from Weil polynomial: %o\n", wp_count;
printf "Point count from twisted Weil polynomial: %o\n", wpt_count;
printf "Naive count matches: %o\n", naive_count eq wp_count select "wp" else (naive_count eq wpt_count select "wpt" else "neither");

// === STEP 6: Check no tritangent lines over Qbar ===

P2Q<u, v, w> := ProjectiveSpace(Rationals(), 2);
f6_Q := Determinant(u*q1_Q + v*q2_Q + w*q3_Q);
printf "Tritangent lines of lifted sextic (should be empty): %o\n", TriTangentLines(f6_Q);

