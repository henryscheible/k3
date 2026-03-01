// Read results from compare_point_counts_generic.csv and verify point counts

function PointCountsFromWP(wp_poly, p, n)
    C := CompanionMatrix(wp_poly);
    return [ Trace(C^k) + p^(2*k) + 1 : k in [1..n] ];
end function;

// Parse upper-triangle string like "[-5 -5 -2 1 5 ...]" into a 6x6 symmetric matrix over F
function ParseMatrixStr(s, F)
    s := Trim(s);
    s := s[2..#s-1];  // strip [ and ]
    parts := Split(s, " ");
    entries := [F!StringToInteger(Trim(e)) : e in parts | Trim(e) ne ""];
    M := ZeroMatrix(F, 6, 6);
    k := 1;
    for i in [1..6] do
        for j in [i..6] do
            M[i][j] := entries[k];
            M[j][i] := entries[k];
            k +:= 1;
        end for;
    end for;
    return M;
end function;

// Build quadratic form x^T * Q * x in CoordinateRing(P5)
function QuadraticForm(Q, R5, vars)
    v := Matrix(R5, 1, 6, vars);
    return (v * ChangeRing(Q, R5) * Transpose(v))[1][1];
end function;

F_file := Open("compare_point_counts_generic.csv", "r");
row := 0;

while true do
    line := Gets(F_file);
    if IsEof(line) then
        break;
    end if;
    row +:= 1;

    parts := Split(line, ",");
    // parts: p, Q1_str, Q2_str, Q3_str, wp, wpt, elapsed
    p := StringToInteger(Trim(parts[1]));
    Fp := FiniteField(p);

    // Reconstruct symmetric matrices over F_p
    Q1 := ParseMatrixStr(parts[2], Fp);
    Q2 := ParseMatrixStr(parts[3], Fp);
    Q3 := ParseMatrixStr(parts[4], Fp);

    // Build K3 in P5 over F_p and count points naively
    P5<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(Fp, 5);
    R5 := CoordinateRing(P5);
    vars := [x0, x1, x2, x3, x4, x5];
    f1 := QuadraticForm(Q1, R5, vars);
    f2 := QuadraticForm(Q2, R5, vars);
    f3 := QuadraticForm(Q3, R5, vars);
    Xp := Scheme(P5, [f1, f2, f3]);
    naive_1 := #Points(Xp);

    // Extract Weil polynomials and compute predicted counts for r=1..22
    R<t> := PolynomialRing(Integers());
    wp := eval Trim(parts[5]);
    wpt := eval Trim(parts[6]);

    wp_counts := PointCountsFromWP(wp, p, 22);
    wpt_counts := PointCountsFromWP(wpt, p, 22);

    printf "=== Row %o (p=%o) ===\n", row, p;
    printf "r=1: naive=%o, wp=%o, wpt=%o\n", naive_1, wp_counts[1], wpt_counts[1];
    printf "wp  counts (r=1..22): %o\n", wp_counts;
    printf "wpt counts (r=1..22): %o\n\n", wpt_counts;
end while;

delete F_file;
