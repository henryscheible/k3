logs := Read("k3_wp_logs.csv");

R<x> := PolynomialRing(Rationals());

data := Split(Split(logs, "\n")[25], ",");


// NTS not cyclotomic

p := data[1];
f2 := data[3];
f3 := data[4];
q := data[5];
f6 := data[6];
wp1 := data[7];
wp2 := data[8];

wp1 := eval wp1;
Factorization(wp1);
weilpoly := Evaluate(wp1, x/47);
polys := [Evaluate(CyclotomicPolynomial(i), x) : i in [1..20]];
[weilpoly in ideal<R | p> : p in polys];

P4<x0, x1, x2, x3, x4> := ProjectiveSpace(Rationals(), 4);
f2 := eval f2;
f3 := eval f3;
p := eval p;
X := Scheme(P4, [f2, f3]);
X;
IsNonsingular(X);
IsNonsingular(ChangeRing(X, FiniteField(47)));
// Dimension(X);
// Degree(X);

f2_L := f2 + p*x3^2 + p * x4^2;
f3_L := f3;

Y := Scheme(P4, [f2_L, f3_L]);

C := CoordinateRing(P4);

R<y1, y2, y3, y4, y5, y6> := PolynomialRing(Rationals(), 6);
S<u, v> := PolynomialRing(R, 2);

function getLineScheme(v1, v2)
    hom_cs := hom<C -> S | [u * v1[i] + v * v2[i]: i in [1..5]]>;
    hom_cs(f2_L);
    Coefficients(hom_cs(f2_L));
    coeffs := Coefficients(hom_cs(f2_L)) cat Coefficients(hom_cs(f3_L));
    return IsProper(ideal<R | coeffs>);
end function;

print "Check Ideals";

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

f2;
f3;
