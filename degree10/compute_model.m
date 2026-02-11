Q := Rationals();
P9<[x]> := ProjectiveSpace(Q, 9);

RandomLinearEquation := function()
    coeffs := [ Random([-10..10]) : i in [1..10] ];
    while &and[ c eq 0 : c in coeffs ] do
        coeffs := [ Random([-10..10]) : i in [1..10] ]; // Ensure not all zero
    end while;
    return &+[ coeffs[i]*x[i] : i in [1..10] ];
end function;


RandomP6 := function()
    return Scheme(P9, [RandomLinearEquation(), RandomLinearEquation(), RandomLinearEquation()]);
end function;

RP6 := RandomP6();
G := Grassmanian(2, 5);
