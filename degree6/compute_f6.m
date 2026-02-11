load "sextic_k3_util.m";

p := 5;

h1 := Random(1, p-1) * x0 + Random(1, p-1) * x1 + Random(1, p-1) * x2;
h2 := Random(1, p-1) * x0 + Random(1, p-1) * x1 + Random(1, p-1) * x2;

// h1 := Random(1, 10) * x1 + Random(1, 10) * x2;
// h2 := Random(1, 10) * x1 + Random(1, 10) * x2;

I := LineIdeal(P4);

X, f2, f3 := random_sextic_k3_containing_ideal(P4, I, p);
X, f2, f3 := reduce_mod_p(X, f2, f3, p);
P4 := ChangeRing(P4, FiniteField(p));

R4 := CoordinateRing(P4);

P2<y0, y1, y2> := ProjectiveSpace(FiniteField(p), 2);
R2 := CoordinateRing(P2);

IntIdeal := ideal<R4 | h1, h2, f2, f3>;
// if Evaluate(h1, [1, 0, 0, 0, 0]) = 0 then 


S := EliminationIdeal(IntIdeal, 2);
Gens := Generators(S);

GensOverP2 := [Evaluate(p, [0, 0, y0, y1, y2]): p in Gens];
IdealP2 := Ideal(GensOverP2);
IdealP2;

residual_portion := ColonIdeal(IdealP2, ideal<R2 | y0>);
residual_portion_eliminated := EliminationIdeal(residual_portion, 1);
P1<z0, z1> := ProjectiveSpace(FiniteField(p), 1);
residual_portion;
residual_portion_eliminated;
p1_eq := Evaluate(Generators(residual_portion_eliminated)[1], [0, z0, z1]);
p1_eq;
IsNonsingular(Scheme(P1, [p1_eq]));
