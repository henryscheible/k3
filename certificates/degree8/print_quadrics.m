p := 47;

q1_Q := Matrix(Rationals(), [[-136, -232, -136,  144, -146,  -42],
                              [-232, -140,  187,   93,  -98,  193],
                              [-136,  187,    1,  234,  137,  281],
                              [ 144,   93,  234,   47,  -94, -141],
                              [-146,  -98,  137,  -94,  237, -137],
                              [ -42,  193,  281, -141, -137, -139]]);

q2_Q := Matrix(Rationals(), [[  43,   44,  -50, -140,  281,  -49],
                              [  44,  141,  192,   -4,   95,  -91],
                              [ -50,  192,  185, -188,  281, -234],
                              [-140,   -4, -188,  -89,   52,   95],
                              [ 281,   95,  281,   52,  144,  -42],
                              [ -49,  -91, -234,   95,  -42,  -44]]);

q3_Q := Matrix(Rationals(), [[ 193,   -1, -146,    1,    4,  235],
                              [  -1,  234,   95,  -90,  137, -140],
                              [-146,   95,   51,    3,   92,  280],
                              [   1,  -90,    3,  183,  143,   45],
                              [   4,  137,   92,  143,    0,  -98],
                              [ 235, -140,  280,   45,  -98,  185]]);

// Print as polynomials over Q
P5Q<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(Rationals(), 5);
R5Q := CoordinateRing(P5Q);
cv_Q := Matrix(R5Q, [[x0, x1, x2, x3, x4, x5]]);
g1_Q := (cv_Q * ChangeRing(q1_Q, R5Q) * Transpose(cv_Q))[1][1];
g2_Q := (cv_Q * ChangeRing(q2_Q, R5Q) * Transpose(cv_Q))[1][1];
g3_Q := (cv_Q * ChangeRing(q3_Q, R5Q) * Transpose(cv_Q))[1][1];

printf "Over Q:\n";
printf "g1 = %o\n", g1_Q;
printf "g2 = %o\n", g2_Q;
printf "g3 = %o\n\n", g3_Q;

// Print as polynomials over F_p
P5p<x0, x1, x2, x3, x4, x5> := ProjectiveSpace(FiniteField(p), 5);
R5p := CoordinateRing(P5p);
cv_p := Matrix(R5p, [[x0, x1, x2, x3, x4, x5]]);
g1_p := (cv_p * ChangeRing(ChangeRing(q1_Q, FiniteField(p)), R5p) * Transpose(cv_p))[1][1];
g2_p := (cv_p * ChangeRing(ChangeRing(q2_Q, FiniteField(p)), R5p) * Transpose(cv_p))[1][1];
g3_p := (cv_p * ChangeRing(ChangeRing(q3_Q, FiniteField(p)), R5p) * Transpose(cv_p))[1][1];

printf "Over F_%o:\n", p;
printf "g1 = %o\n", g1_p;
printf "g2 = %o\n", g2_p;
printf "g3 = %o\n", g3_p;
