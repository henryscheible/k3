Qbar := Rationals(); // AlgebraicClosure(Rationals());

P5<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(Qbar, 5);
P2<u, v, w> := ProjectiveSpace(Qbar, 2);

R := CoordinateRing(P5);

RandomSymmetricMatrix := function(n)
    Qb := Qbar;
    M := ZeroMatrix(Qb, n, n);
    for i in [1..n] do
        for j in [i..n] do
            a := Qb!Random([-5..5]);
            M[i][j] := a;
            M[j][i] := a;
        end for;
    end for;
    return M;
end function;

q1 := Matrix(Qbar, [[ 5,  3,  5,  3, 42,  5],
[ 3,  1, 46, 46, 43,  5],
[ 5, 46,  1, 46, 43, 46],
[ 3, 46, 46,  0,  0,  0],
[42, 43, 43,  0,  2,  4],
[ 5,  5, 46,  0,  4,  2]]);
q1 := q1 + 47 * RandomSymmetricMatrix(6);

q2 := Matrix(Qbar, [[43, 44, 44,  1, 46, 45],
[44,  0 , 4, 43,  1,  3],
[44,  4, 44,  0, 46,  1],
[ 1, 43,  0,  5,  5,  1],
[46,  1, 46,  5,  3,  5],
[45 , 3 , 1 , 1 , 5 , 3]]);
q2 := q2 + 47 * RandomSymmetricMatrix(6);

q3 := Matrix(Qbar, [[ 5, 46, 42,  1,  4,  0],
[46, 46,  1,  4, 43,  1],
[42,  1,  4,  3, 45, 45],
[ 1,  4,  3, 42,  2, 45],
[ 4, 43, 45,  2,  0, 43],
[ 0,  1, 45, 45, 43, 44]]);
q3 := q3 + 47 * RandomSymmetricMatrix(6);

f6 := Determinant(u*q1 + v*q2 + w*q3);
lines := TriTangentLines(f6);
print lines;

coord_vec := Matrix(CoordinateRing(P5), [[x0, x1, x2, x3, x4, x5]]);

coord_vec * ChangeRing(q1, R) * Transpose(coord_vec);
coord_vec * ChangeRing(q2, R) * Transpose(coord_vec);
coord_vec * ChangeRing(q3, R) * Transpose(coord_vec);

