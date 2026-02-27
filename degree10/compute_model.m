// compute_model.m
// Degree 10 K3 surface as a linear section of the Grassmannian G(2,5).
//
// Construction: S = G(2,5) ∩ H_1 ∩ H_2 ∩ H_3 ∩ Q  in P^9
//   where H_1, H_2, H_3 are hyperplanes cutting out a P^6,
//   and Q is a quadric hypersurface.
//
// G(2,5) has dimension 6 and degree 5 in P^9.
// Intersecting with 3 hyperplanes gives a 3-fold of degree 5.
// Intersecting with Q gives a surface of degree 10.

load "grassmannian.m";

function RandomLinear(P, b)
    R := CoordinateRing(P);
    return &+[Random(-b, b) * m : m in MonomialsOfDegree(R, 1)];
end function;

function RandomQuadric(P, b)
    R := CoordinateRing(P);
    return &+[Random(-b, b) * m : m in MonomialsOfDegree(R, 2)];
end function;

function random_degree10_k3(P9)
    G := Grassmannian25(P9);
    is_viable := false;
    while not is_viable do
        L1 := RandomLinear(P9, 5);
        L2 := RandomLinear(P9, 5);
        L3 := RandomLinear(P9, 5);
        Q  := RandomQuadric(P9, 5);
        S  := G meet Scheme(P9, [L1, L2, L3, Q]);
        is_viable := Dimension(S) eq 2 and IsNonsingular(S);
    end while;
    return S;
end function;

K := Rationals();
P9<[x]> := ProjectiveSpace(K, 9);

S := random_degree10_k3(P9);
print "Dimension:", Dimension(S);
print "Degree:", Degree(S);
