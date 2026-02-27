// Read results from compare_point_counts_generic.csv and verify point counts

function PointCountsFromWP(wp_poly, p, n)
    C := CompanionMatrix(wp_poly);
    return [ Trace(C^k) + p^(2*k) + 1 : k in [1..n] ];
end function;

F := Open("compare_point_counts_generic.csv", "r");
row := 0;

while true do
    line := Gets(F);
    if IsEof(line) then
        break;
    end if;
    row +:= 1;

    parts := Split(line, ",");
    // parts: p, f2, f3, wp, wpt, elapsed
    p := StringToInteger(Trim(parts[1]));

    // Setup P4 over F_p and eval f2, f3
    P4<x0, x1, x2, x3, x4> := ProjectiveSpace(FiniteField(p), 4);
    S := CoordinateRing(P4);
    f2 := eval Trim(parts[2]);
    f3 := eval Trim(parts[3]);

    // Naive point counts over F_p and F_p^3
    Xp := Scheme(P4, [f2, f3]);
    naive_1 := #Points(Xp);

    // Extract point counts from Weil polynomials for r=1..22
    R<t> := PolynomialRing(Integers());
    wp := eval Trim(parts[4]);
    wpt := eval Trim(parts[5]);

    wp_counts := PointCountsFromWP(wp, p, 22);
    wpt_counts := PointCountsFromWP(wpt, p, 22);

    printf "=== Row %o (p=%o) ===\n", row, p;
    printf "r=1: naive=%o, wp=%o, wpt=%o\n", naive_1, wp_counts[1], wpt_counts[1];
    printf "wp  counts (r=1..22): %o\n", wp_counts;
    printf "wpt counts (r=1..22): %o\n\n", wpt_counts;
end while;

delete F;
