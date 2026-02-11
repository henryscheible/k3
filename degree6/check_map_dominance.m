logs := Read("k3_wp_logs.csv");

R<x0,x1,x2,x3,x4> := PolynomialRing(Rationals(),5);
A<x0A,x1A,x2A> := PolynomialRing(Rationals(),3);
S<y0,y1> := PolynomialRing(A,2);
phi := hom< R -> S | x0A, x1A, x2A, y0, y1 >;

i := 25;
f2 := Split(Split(logs, "\n")[i], ",")[3];
f3 := Split(Split(logs, "\n")[i], ",")[4];

printf "f2: %o\n\n", f2;
printf "f3: %o\n\n", f3;
f2 := eval f2;
f3 := eval f3;

l0 := MonomialCoefficient(phi(f2), y0);
l1 := MonomialCoefficient(phi(f2), y1);
l00 := MonomialCoefficient(phi(f3), y0^2);
l01 := MonomialCoefficient(phi(f3), y0*y1);
l11 := MonomialCoefficient(phi(f3), y1^2);
g := MonomialCoefficient(phi(f3), x3S*x4S);
e*c^2 + f*b^2 - g*c*b;
