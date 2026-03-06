# k3

## Overview

This repository contains **Magma** scripts for experimental work on K3 surfaces over finite fields.  The code implements algorithms for generating random K3 surfaces of degree 6, 8 and 10, computing their Weil polynomials, counting points, and certifying Picard‑rank statements.

## Directory structure

- `degree6/` – Scripts dealing with sextic K3 surfaces (double covers of \(\mathbb{P}^2\) branched along a sextic).  Includes generators, point‑count comparison, and certificates.
- `degree8/` – Scripts for degree‑8 K3 surfaces (complete intersections of three quadrics in \(\mathbb{P}^5\)).  Contains utilities for tritangent‑line searches and verification of point counts.
- `degree10/` – Scripts for degree‑10 K3 surfaces (linear sections of the Grassmannian \(G(2,5)\)).
- `certificates/` – Pre‑computed certificates proving Picard‑rank statements for specific examples.
- `general/` – Placeholder for shared utilities (currently minimal).

## Prerequisites

- **Magma Computational Algebra System** (https://magma.maths.usyd.edu.au/magma/).  The scripts use native Magma syntax (`PolynomialRing`, `Factorization`, etc.) and rely on a recent version (≥ 2.27‑...).
- (Optional) **Python 3** – some auxiliary scripts such as `degree8/test.py` are written in Python.

## Quick start

```bash
# Launch Magma
magma

# Load a script, e.g. to compute a random degree‑6 K3 surface and its Weil polynomial:
Attach("degree6/compute_f6.m");
```

Typical workflow:
1. Generate a random surface using one of the `compute_f*.m` scripts.
2. Compute its Weil polynomial with `WeilPolynomialOfDegree2K3Surface`.
3. Verify point counts with `compare_point_counts_generic.m` (or the “small” variant).
4. Use the certificate scripts in `certificates/` to certify Picard rank.

## Example (degree 6)

```magma
R<x,y,z> := PolynomialRing(Rationals());
// Random sextic defining a K3 surface as a double cover of P^2
f2 := Random(R, 2);
f3 := Random(R, 3);
f6 := f2^3 + f3^2; // sextic polynomial
wp1, wp2 := WeilPolynomialOfDegree2K3Surface(f6);
print wp1, wp2;
```


