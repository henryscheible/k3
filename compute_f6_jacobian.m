load "sextic_k3_util.m";


p := 3;

for i := 1 to 5 do

    T := Time();

    P4<x0, x1, x2, x3, x4> := ProjectiveSpace(FiniteField(p), 4);

    I := ideal<CoordinateRing(P4) | x0, x1, x2>;

    X, f2, f3 := random_sextic_k3_containing_ideal(P4, I, p);

    M := JacobianMatrix([f2, f3]);
    MSub := Submatrix(M, 1, 4, 2, 2);

    q := Determinant(Submatrix(M, 1, 4, 2, 2));
    Z := Scheme(P4, [f2, f3, q]);

    P2<y0, y1, y2> := ProjectiveSpace(FiniteField(p), 2);
    R2 := CoordinateRing(P2);

    f := map<P4 -> P2 | [x0, x1, x2] >;

    f6 := Generators(Ideal(f(Z)))[1];
    SetVerbose("Degree2K3", 2);

    wp1, wp2 := WeilPolynomialOfDegree2K3Surface(f6);
    R<x> := PolynomialRing(Rationals());

    // Check that scheme defined by f6 is smooth

    fprintf "k3_wp_logs.csv", "%o, %o, %o, %o, %o, %o, %o, %o\n", p, Time(T), f2, f3, q, f6, Evaluate(wp1, x), Evaluate(wp2, x);

end for;
// GensOverP2 := [Evaluate(p, [y0, y1, y2, 0, 0]): p in Generators(Ideal(Z))];
// Y := Scheme(P2, GensOverP2);

// Dimension(Y);
// Degree(Y);


// F<a, b, c, d, e, f> := RationalFunctionField(FiniteField(p), 6);

// P4<x0, x1, x2, x3, x4> := ProjectiveSpace(F, 4);

// h1 := a * x0 + b * x1 + c * x2;
// h2 := d * x0 + e * x1 + f * x2;

// I := ideal<CoordinateRing(P4) | x0, x1, x2>;

// X, f2, f3 := reduce_to_indeterminate_ring(X, f2, f3, F);
// P4 = ChangeRing(P4, FiniteField(p));

// R4 := CoordinateRing(P4);

// P2<y0, y1, y2> := ProjectiveSpace(FiniteField(p), 2);
// R2 := CoordinateRing(P2);

// IntIdeal := ideal<R4 | h1, h2, f2, f3>;
// // if Evaluate(h1, [1, 0, 0, 0, 0]) = 0 then 


// S := EliminationIdeal(IntIdeal, 2);
// Gens := Generators(S);

// GensOverP2 := [Evaluate(p, [0, 0, y0, y1, y2]): p in Gens];
// IdealP2 := Ideal(GensOverP2);
// IdealP2;

// residual_portion := ColonIdeal(IdealP2, ideal<R2 | y0>);
// residual_portion_eliminated := EliminationIdeal(residual_portion, 1);
// P1<z0, z1> := ProjectiveSpace(FiniteField(p), 1);
// residual_portion;
// residual_portion_eliminated;
// p1_eq := Evaluate(Generators(residual_portion_eliminated)[1], [0, z0, z1]);
// p1_eq;
// IsNonsingular(Scheme(P1, [p1_eq]));
