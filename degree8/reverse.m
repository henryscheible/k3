Q := GF(5);

// Number of unknown symmetric entries per 6x6 matrix: 6*7/2 = 21
numSym := 21;

// Total variables: u,v,w plus 3*21 matrix variables
totalVars := 6+1*numSym;
S := PolynomialRing(Q, totalVars);
SS<u, v, w> := PolynomialRing(S, 3);
gens := [S.i : i in [1..totalVars]];

// Helper: allocate symmetric matrix M over S using successive generators
makeSymmetricMatrixFromGens := function(startIndex)
    M := ZeroMatrix(S, 6, 6);
    pos := startIndex;
    for i in [1..6] do
        for j in [i..6] do
            M[i][j] := gens[pos];
            M[j][i] := gens[pos];
            pos +:= 1;
        end for;
    end for;
    return M, pos;
end function;

// Build symbolic Q1,Q2,Q3 using the polynomial ring generators
pos := 7; // S.1..S.3 are u,v,w
Q1 := IdentityMatrix(GF(5), 6);
Q2 := DiagonalMatrix([S.i: i in [1..6]]);
Q3, pos := makeSymmetricMatrixFromGens(pos);

// User input / example: define a target sextic f6 in u,v,w.
// Replace the random construction below with your concrete f6 if available.
RandomSymmetricMatrixOverS := function()
    M := ZeroMatrix(S,6,6);
    for i in [1..6] do
        for j in [i..6] do
            // small integer constants embedded in S
            a := S!(Random([-3..3]));
            M[i][j] := a;
            M[j][i] := a;
        end for;
    end for;
    return M;
end function;

// By default construct a random instance and compute the corresponding f6.
// The user can override f6 by editing the line below.
Q1_true := IdentityMatrix(GF(5), 6);
Q2_true := DiagonalMatrix([S!Random([-3..3]): i in [1..6]]);
Q3_true := RandomSymmetricMatrixOverS();
f6 := Determinant(u*Q1_true + v*Q2_true + w*Q3_true);

print "Built test f6 from random symmetric matrices (edit script to supply f6 directly): ";
f6;

Q1; Q2; Q3;

// Build symbolic determinant
M := u*Q1 + v*Q2 + w*Q3;
detM := Determinant(M);

// Collect degree-6 monomials in u,v,w (there are 28 of them)
mons := [];
for a in [0..6] do
    for b in [0..6-a] do
        c := 6 - a - b;
        Append(~mons, u^a * v^b * w^c);
    end for;
end for;

// Form equations by equating coefficients of corresponding monomials
eqs := [];
// Helper: safe coefficient extractor for multivariate polynomials
getCoeff := function(f, mon)
    monsF := Monomials(f);
    coeffsF := Coefficients(f);
    for i in [1..#monsF] do
        if monsF[i] eq mon then
            return coeffsF[i];
        end if;
    end for;
    return SS!0;
end function;
for m in mons do
    c_det := getCoeff(detM, m);
    c_target := getCoeff(f6, m);
    eqn := c_det - c_target;
    if not IsZero(eqn) then
        Append(~eqs, eqn);
    end if;
end for;

print "Number of unknown matrix variables:", totalVars;
print "Number of nonzero coefficient equations:", #eqs;

// Display a few example equations
// print "Example equations (first 8):";
// for i in [1..Min(#eqs,8)] do
//     print i, ":", eqs[i];
// end for;

// Package into an ideal (over S)
I := ideal< S | eqs >;

// Optionally compute a Groebner basis. This is expensive for realistic
// instances. Set to true only if you know what you're doing and have
// plenty of time / memory.
COMPUTE_GB := true;
if COMPUTE_GB then
    print "Computing Groebner basis...";
    G := GroebnerBasis(I);
    print "Groebner basis length:", #G;
    for i in [1..#G] do
        print "G[", i, "] = ", G[i];
    end for;
end if;

// Save equations to a file for external processing if desired
out := Open("reverse_system.txt", "w");
for eqn in eqs do
    Put(out, Sprint(eqn) * "\n");
end for;
Flush(out);
// Close(out); // Close is not available/declared in this Magma environment; leaving file to be closed on program exit
print "Wrote polynomial system to reverse_system.txt";