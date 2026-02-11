Qbar := GF(47); // AlgebraicClosure(Rationals());

P5<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(Qbar, 5);
P2<u, v, w> := ProjectiveSpace(Qbar, 2);

R := CoordinateRing(P5);

q1 := Matrix(R, [[ 5,  3,  5,  3, 42,  5],
[ 3,  1, 46, 46, 43,  5],
[ 5, 46,  1, 46, 43, 46],
[ 3, 46, 46,  0,  0,  0],
[42, 43, 43,  0,  2,  4],
[ 5,  5, 46,  0,  4,  2]]);

q2 := Matrix(R, [[43, 44, 44,  1, 46, 45],
[44,  0 , 4, 43,  1,  3],
[44,  4, 44,  0, 46,  1],
[ 1, 43,  0,  5,  5,  1],
[46,  1, 46,  5,  3,  5],
[45 , 3 , 1 , 1 , 5 , 3]]);

q3 := Matrix(R, [[ 5, 46, 42,  1,  4,  0],
[46, 46,  1,  4, 43,  1],
[42,  1,  4,  3, 45, 45],
[ 1,  4,  3, 42,  2, 45],
[ 4, 43, 45,  2,  0, 43],
[ 0,  1, 45, 45, 43, 44]]);

coord_vec := Matrix(CoordinateRing(P5), [[x0, x1, x2, x3, x4, x5]]);

coord_vec * q1 * Transpose(coord_vec);
coord_vec * q2 * Transpose(coord_vec);
coord_vec * q3 * Transpose(coord_vec);
