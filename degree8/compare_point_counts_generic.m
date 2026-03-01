// Setup rings over Q
P5<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(Rationals(), 5);
R5 := CoordinateRing(P5);
Runi<t> := PolynomialRing(Rationals());

// Random symmetric 6x6 matrix with entries in [-b, b]
function RandomSymmetricMatrix(b)
    M := ZeroMatrix(Rationals(), 6, 6);
    for i in [1..6] do
        for j in [i..6] do
            a := Random(-b, b);
            M[i][j] := a;
            M[j][i] := a;
        end for;
    end for;
    return M;
end function;

// Build quadratic form x^T * Q * x in CoordinateRing(P5)
function QuadraticForm(Q)
    v := Matrix(R5, 1, 6, [x0, x1, x2, x3, x4, x5]);
    return (v * ChangeRing(Q, R5) * Transpose(v))[1][1];
end function;

// Random prime between 41 and 67
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

// Compact upper-triangle string for a symmetric matrix (single line, CSV-safe)
function MatrixStr(M)
    n := NumberOfRows(M);
    s := "[";
    sep := "";
    for i in [1..n] do
        for j in [i..n] do
            s := s * sep * Sprint(M[i][j]);
            sep := " ";
        end for;
    end for;
    return s * "]";
end function;

iteration := 0;
while true do
    iteration +:= 1;
    T := Time();

    p := RandomPrime();
    printf "=== Iteration %o (p = %o) ===\n", iteration, p;

    // Generate random smooth degree-8 K3 (complete intersection of three quadrics in P5)
    is_viable := false;
    while not is_viable do
        Q1 := RandomSymmetricMatrix(5);
        Q2 := RandomSymmetricMatrix(5);
        Q3 := RandomSymmetricMatrix(5);
        f1 := QuadraticForm(Q1);
        f2 := QuadraticForm(Q2);
        f3 := QuadraticForm(Q3);
        X := Scheme(P5, [f1, f2, f3]);
        is_viable := IsNonsingular(X) and IsNonsingular(ChangeRing(X, FiniteField(p)));
    end while;

    printf "Q1 = %o\n", Q1;
    printf "Q2 = %o\n", Q2;
    printf "Q3 = %o\n", Q3;

    // Reduce matrices mod p and compute Weil polynomial via sextic discriminant
    p_field := FiniteField(p);
    P2<u, v, w> := ProjectiveSpace(p_field, 2);
    Q1_p := ChangeRing(Q1, p_field);
    Q2_p := ChangeRing(Q2, p_field);
    Q3_p := ChangeRing(Q3, p_field);

    f6_p := Determinant(u*Q1_p + v*Q2_p + w*Q3_p);

    ns := NonSquare(p);
    wp1, wp2 := WeilPolynomialOfDegree2K3Surface(f6_p);
    wpt1, wpt2 := WeilPolynomialOfDegree2K3Surface(ns * f6_p);

    wp := wp1 * wp2;
    wpt := wpt1 * wpt2;

    elapsed := Time(T);
    printf "Done in %o seconds\n", elapsed;
    printf "wp = %o\n", wp;
    printf "wpt = %o\n\n", wpt;

    // Append to CSV: p, Q1 (upper tri), Q2 (upper tri), Q3 (upper tri), wp(t), wpt(t), elapsed
    fprintf "compare_point_counts_generic.csv", "%o, %o, %o, %o, %o, %o, %o\n",
        p, MatrixStr(Q1), MatrixStr(Q2), MatrixStr(Q3), Evaluate(wp, t), Evaluate(wpt, t), elapsed;

end while;
