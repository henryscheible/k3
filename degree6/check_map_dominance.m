logs := Read("k3_wp_logs.csv");

R<x0,x1,x2,x3,x4> := PolynomialRing(Rationals(),5);
A<x0A,x1A,x2A> := PolynomialRing(Rationals(),3);
S<x3S,x4S> := PolynomialRing(A,2);
phi := hom< R -> S | x0A, x1A, x2A, x3S, x4S >;

i := 25;
f2 := Split(Split(logs, "\n")[i], ",")[3];
f3 := Split(Split(logs, "\n")[i], ",")[4];
f2 := eval f2;
f3 := eval f3;

b := MonomialCoefficient(phi(f2), x3S);
c := MonomialCoefficient(phi(f2), x4S);
e := MonomialCoefficient(phi(f3), x3S^2);
f := MonomialCoefficient(phi(f3), x4S^2);
g := MonomialCoefficient(phi(f3), x3S*x4S);
e*c^2 + f*b^2 - g*c*b;
