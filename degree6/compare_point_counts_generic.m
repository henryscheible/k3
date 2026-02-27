// Setup rings over Q
P4<x0, x1, x2, x3, x4> := ProjectiveSpace(Rationals(), 4);
R3<b0, b1, b2> := PolynomialRing(Rationals(), 3);
Rfib<y0, y1> := PolynomialRing(R3, 2);
S := CoordinateRing(P4);
I := ideal<S | x0, x1, x2>;
Runi<t> := PolynomialRing(Rationals());

// Random homogeneous polynomial helpers
function RandomHomog(R, k, b)
    return &+[Random(-b,b)*m : m in MonomialsOfDegree(R,k)];
end function;

function RandomHomogInIdeal(R, II, k, b)
    return &+[Random(-b,b)*m : m in [m : m in MonomialsOfDegree(R,k) | m in II]];
end function;

// Helper: extract coefficient of a y-monomial from a polynomial in Rfib
function FiberCoeff(f, mon)
    coeffs, mons := CoefficientsAndMonomials(f);
    for i in [1..#mons] do
        if mons[i] eq mon then
            return coeffs[i];
        end if;
    end for;
    return R3!0;
end function;

// Random prime between 40 and 70
function RandomPrime()
    primes := [41, 43, 47, 53, 59, 61, 67];
    return primes[Random(1, #primes)];
end function;

// Find a non-square in F_p
function NonSquare(p)
    F := FiniteField(p);
    for i in [2..p-1] do
        if not IsSquare(F!i) then
            return i;
        end if;
    end for;
end function;

phi := hom<S -> Rfib | b0, b1, b2, y0, y1>;

iteration := 0;
while true do
    iteration +:= 1;
    T := Time();

    p := RandomPrime();
    printf "=== Iteration %o (p = %o) ===\n", iteration, p;

    // Generate random smooth K3 surface with f2, f3 in (x0, x1, x2)
    is_viable := false;
    while not is_viable do
        f2 := RandomHomogInIdeal(S, I, 2, 5);
        f3 := RandomHomogInIdeal(S, I, 3, 5);
        X := Scheme(P4, [f2, f3]);
        is_viable := IsNonsingular(X) and IsNonsingular(ChangeRing(X, FiniteField(p)));
    end while;

    printf "f2 = %o\n", f2;
    printf "f3 = %o\n", f3;

    // Decompose into fiber form
    F2 := phi(f2);
    F3 := phi(f3);

    l0 := FiberCoeff(F2, y0);
    l1 := FiberCoeff(F2, y1);
    q := FiberCoeff(F2, Rfib!1);

    l00 := FiberCoeff(F3, y0^2);
    l01 := FiberCoeff(F3, y0*y1);
    l11 := FiberCoeff(F3, y1^2);
    q0 := FiberCoeff(F3, y0);
    q1 := FiberCoeff(F3, y1);
    c := FiberCoeff(F3, Rfib!1);

    g6 := -1 * l0^2 * q1^2 - l1^2*q0^2 - l01^2*q^2 + 4*l00*l11*q^2
        + 2*l0*l1*q0*q1 + 2*l0*l01*q*q1 + 2*l1*l01*q*q0 - 4*l1*l00*q*q1-4*l0*q*l11*q0
        + 4*l1^2*l00*c + 4*l0^2*l11*c - 4*l0*l1*l01*c;

    // Reduce g6 to F_p and compute Weil polynomial
    P2<z0, z1, z2> := ProjectiveSpace(FiniteField(p), 2);
    psi := hom<R3 -> CoordinateRing(P2) | z0, z1, z2>;
    g6_p := psi(g6);

    ns := NonSquare(p);
    wp1, wp2 := WeilPolynomialOfDegree2K3Surface(g6_p);
    wpt1, wpt2 := WeilPolynomialOfDegree2K3Surface(ns*g6_p);

    wp := wp1 * wp2;
    wpt := wpt1 * wpt2;

    elapsed := Time(T);
    printf "Done in %o seconds\n", elapsed;
    printf "wp = %o\n", wp;
    printf "wpt = %o\n\n", wpt;

    // Append to CSV
    fprintf "compare_point_counts_generic.csv", "%o, %o, %o, %o, %o, %o\n",
        p, f2, f3, Evaluate(wp, t), Evaluate(wpt, t), elapsed;

end while;
