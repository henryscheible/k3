p := 47;

R<x0,x1,x2,x3,x4> := PolynomialRing(GF(p),5);
A<x0A,x1A,x2A> := PolynomialRing(GF(p),3);
S<y0,y1> := PolynomialRing(A,2);
phi := hom< R -> S | x0A, x1A, x2A, y0, y1 >;

f2_string := "x0^2 + 44*x0*x1 + 3*x1^2 + 5*x0*x2 + 4*x1*x2 + 5*x2^2 + 46*x0*x3 + 45*x1*x3
+ 44*x2*x3 + 42*x0*x4 + 5*x1*x4";
f3_string := "2*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + x1^3 + 46*x0*x1*x2 + 44*x1^2*x2 +
4*x0*x2^2 + 43*x1*x2^2 + 5*x2^3 + 4*x0^2*x3 + x0*x1*x3 + 5*x1^2*x3 + 4*x0*x2*x3
+ 4*x1*x2*x3 + 44*x2^2*x3 + 4*x1*x3^2 + 46*x2*x3^2 + 5*x0^2*x4 + 43*x0*x1*x4 +
2*x1^2*x4 + x0*x2*x4 + 4*x1*x2*x4 + 45*x2^2*x4 + 4*x0*x3*x4 + 44*x2*x3*x4 +
46*x0*x4^2 + 46*x1*x4^2 + 5*x2*x4^2";

f2 := eval f2_string;
f3 := eval f3_string;

l0 := MonomialCoefficient(phi(f2), y0);
l1 := MonomialCoefficient(phi(f2), y1);
q := MonomialCoefficient(phi(f2), [0, 0]);

l00 := MonomialCoefficient(phi(f3), y0^2);
l01 := MonomialCoefficient(phi(f3), y0*y1);
l11 := MonomialCoefficient(phi(f3), y1^2);
q0 := MonomialCoefficient(phi(f3), y0);
q1 := MonomialCoefficient(phi(f3), y1);
c := MonomialCoefficient(phi(f3), [0, 0]);

// We need to show that the x0, x1, x2 for which the fiber y0, y1 has no solutions is a Zariski closed (or possibly empty) set in P2
// After considering f2 and f3 in S and homogenizing, we have a line and a quadric curve in P2, which meet in exactly two points over Qbar.
// To show that there are solutions away from the line at infinity (values of y0, y1) satisfying the equations, we need to set the third homogenizing coordinate to 0 to check that there are no solutions on the line at infinity.

// Thus, we get equations l0 y0 + l1 y1 = 0 and l00 y0^2 + l01 y0 y1 + l11 y1^2 = 0. These are two zero dimensional subvarieties of P1
// The only point satisfying the first equation is [-l1 : l0]. Plugging this into the second equation, we get
// l00*l1^2 -l01*l0*l1 + l11*l0 = 0.
// This is a cubic equation in x0, x1, x2, and its solutions are points where the fiber
// over [x0 : x1 : x2] could potentially have no solutions. Thus, as long as it is not the zero polynomial,
// the map X -> P2 is dominant
l00*l1^2 -l01*l0*l1 + l11*l0;
