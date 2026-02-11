P5<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(GF(47), 5);
RandomSymmetricMatrix := function(n)
    Q := GF(47);
    M := ZeroMatrix(Q, n, n);
    for i in [1..n] do
        for j in [i..n] do
            a := Q!Random([-5..5]);
            M[i][j] := a;
            M[j][i] := a;
        end for;
    end for;
    return M;
end function;

// Projective plane for the linear combination variables
P2<u, v, w> := ProjectiveSpace(GF(47), 2);

// Optional: set maxIter to a positive integer to limit attempts. 0 means unlimited.
maxIter := 0;

i := 0;
while true do
    i := i + 1;
    q1 := RandomSymmetricMatrix(6);
    q2 := RandomSymmetricMatrix(6);
    q3 := RandomSymmetricMatrix(6);
    f6 := Determinant(u*q1 + v*q2 + w*q3);
    lines := TriTangentLines(f6);
    n := #lines;
    printf "Iteration %o: %o tri-tangent lines\n", i, n;
    if n eq 1 then
        printf "Found example on iteration %o\n", i;
        print lines;
        printf "q1 := %o;\n\n", q1;
        printf "q2 := %o;\n\n", q2;
        printf "q3 := %o;\n\n", q3;
        wp1, wp2 := WeilPolynomialOfDegree2K3Surface(f6);
        printf "Factored Weil Polynomial: %o", Factorization(wp1 * wp2);
        R<x> := PolynomialRing(Rationals());
        fprintf "k3_wp_logs.csv", "%o, %o, %o, %o, %o, %o\n", q1, q2, q3, f6, Evaluate(wp1, x), Evaluate(wp2, x);
    end if;
    if maxIter ne 0 and i ge maxIter then
        printf "Stopped after reaching maxIter = %o without finding deg 2 K3 with only one tritangent line.\n", maxIter;
        break;
    end if;
end while;