logs := Read("k3_wp_logs.csv");

R<x> := PolynomialRing(Rationals());

// for i := 1 to #Split(logs, "\n") do
//     wp1 := Split(Split(logs, "\n")[i], ",")[7];
//     poly := eval wp1;
//     print(i);
//     Factorization(poly);
//     print "==============================";
// end for;

    // NTS that irreducible part doesn't have a factor which is cyclotomic, cyclotomic must have constant. Also check that deg 6 k3 is smooth mod 47.

i := 25;
wp1 := Split(Split(logs, "\n")[i], ",")[7];
poly := eval wp1;
print(i);
Factorization(poly);
