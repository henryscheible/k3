function RandomHomog(R, k, b)
    return &+[Random(-b,b)*m : m in MonomialsOfDegree(R,k)];
end function;

function RandomHomogInIdeal(R, I, k, b)
    return &+[Random(-b,b)*m : m in [m : m in MonomialsOfDegree(R,k) | m in I]];
end function;

function P4();
    return ProjectiveSpace(Rationals(), 4);
end function;

function LineIdeal(P4);
    return ideal<CoordinateRing(P4)|P4.1, P4.2, P4.3>;
end function;

function random_sextic_k3(P4)
    is_viable := false;
    while not is_viable do 
        R := CoordinateRing(P4);
        f2 := RandomHomog(R, 2, 5);
        f3 := RandomHomog(R, 3, 5);
        X := Scheme(P4, [f2, f3]);
        is_viable := IsNonsingular(X);
    end while;
    return X, f2, f3;
end function;

function random_sextic_k3_containing_ideal(P4, I, p)
    is_viable := false;
    while not is_viable do 
        R := CoordinateRing(P4);
        f2 := RandomHomogInIdeal(R, I, 2, 5);
        f3 := RandomHomogInIdeal(R, I, 3, 5);
        X := Scheme(P4, [f2, f3]);
        is_viable := IsNonsingular(X) and IsNonSingular(ChangeRing(X, FiniteField(p)));
    end while;
    return X, f2, f3;
end function;

function reduce_mod_p(X, f2, f3, p)
    return ChangeRing(X, FiniteField(p)), ChangeRing(f2, FiniteField(p)), ChangeRing(f3, FiniteField(p));
end function;

function reduce_to_indeterminate_ring(X, f2, f3, F)
    return ChangeRing(X, F), ChangeRing(f2, F), ChangeRing(f3, F);
end function;
