F := Rationals();
C<l0, l1, l00, l11, l01, q0, q1, q, c> := PolynomialRing(F, 9);
R<y0, y1> := PolynomialRing(C, 2);

f2 := l0*y0 + l1*y1 + q;
f3 := l00*y0^2 + l11*y1^2 + l01*y0*y1 + q0*y0 + q1*y1 + c;

J := JacobianMatrix([f2, f3]);
D := Determinant(J);
