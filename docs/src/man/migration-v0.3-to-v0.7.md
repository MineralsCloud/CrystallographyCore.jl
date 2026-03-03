# Migration Guide: `v0.3.x` to `v0.7.x`

```@contents
Pages = ["migration-v0.3-to-v0.7.md"]
Depth = 2:3
```

## Summary

- `Lattice` and `ReciprocalLattice` are now immutable, static, and no longer `AbstractMatrix` subtypes.
- `eachbasisvector` was removed; use `basisvectors`.
- Matrix conversion and matrix-action behavior changed (`convert`, `3x3` multiplication rules).
- Transformation support expanded (`inv`, callable lattices, change-of-basis types).
- Unit handling is now extension-based (`UnitfulExt`, `UnitfulLinearAlgebraExt`).
- A few behaviors are now explicitly rejected (`transpose(lattice)`, `Number / lattice`).

## Immutable lattice updates

`Lattice` is immutable in `v0.7+`, so index mutation with `setindex!` is not supported.

Use a functional update pattern. If you want index-like updates, prefer `@set` from `Accessors.jl`:

```julia
using Accessors

lat2 = @set lat.data[1, 1] = 2.0
```

## `eachbasisvector` removed

`eachbasisvector(::Lattice)` and `eachbasisvector(::ReciprocalLattice)` were removed.

Use:

```julia
a, b, c = basisvectors(lat)
```

## `AbstractLattice` no longer subtypes `AbstractMatrix`

Dispatch like `f(A::AbstractMatrix)` no longer accepts lattice objects directly.

Convert explicitly:

```julia
f(convert(Matrix, lat))
# or
f(parent(lat))
```

## `Matrix(lattice)` replaced by `convert`

Lattice-to-matrix conversion now follows:

```julia
convert(::Type{T}, lattice) where {T<:AbstractMatrix}
```

Example:

```julia
A = convert(Matrix, lat)
S = convert(SMatrix{3,3,Float64}, lat)
```

## Six-constant lattice constructor added

`v0.7+` includes:

```julia
Lattice(a, b, c, α, β, γ; axis=:a)
```

`axis=:c` is also supported.

## Matrix-left action on `Lattice`

`R * lattice` is supported and treated as a left (active) action.

Only `3x3` matrices are accepted.

## Matrix-right action on `Lattice`

`lattice * P` is supported and treated as a right (basis-change) action.

Only `3x3` matrices are accepted.

## Non-`3x3` matrix multiplication rejected

Matrix-lattice multiplication now throws `DimensionMismatch` unless matrix shape is exactly `3x3`.

## `transpose(lattice)` rejected

`transpose(lattice)` now throws a `DomainError`.

Use:

```julia
permutedims(lattice)
```

## `Number / Lattice` rejected

`x / lattice` is explicitly rejected (including broadcasted forms like `x ./ lattice`).

Use `lattice / x` when appropriate.

## `Number / ReciprocalLattice` rejected

`x / reciprocal_lattice` is explicitly rejected (including `x ./ reciprocal_lattice`).

Use `reciprocal_lattice / x` when appropriate.

## Callable coordinate transforms

Lattices are callable for reduced-to-Cartesian transforms, and `inv(lattice)` is callable for Cartesian-to-reduced transforms:

```julia
cart = lattice(reduced)
reduced = inv(lattice)(cart)
```

## Inverse-lattice multiplication added

`inv(lat1) * lat2` and `lat2 * inv(lat1)` are supported for relative transform extraction.

## Wrapper-type check for inverse multiplication

Inverse multiplication requires matching wrappers (`Lattice` with `Lattice`, `ReciprocalLattice` with `ReciprocalLattice`), otherwise `ArgumentError` is thrown.

## Primitive/Standardized change-of-basis types added

The following types are available in core:

- `PrimitiveFromStandardized`
- `StandardizedFromPrimitive`
- `PrimitiveToStandardized` (alias)
- `StandardizedToPrimitive` (alias)

## Unitful conversion in extension

`UnitfulExt` provides unit-aware conversions and helpers for lattice types (`convert`, `uconvert`, `ustrip`).

Ensure `Unitful` is present in your environment to load this extension.

## Unitful linear solve in extension

`UnitfulLinearAlgebraExt` provides unit-aware inverse solve behavior for `inv(lattice)(vector_with_units)`.

Ensure `UnitfulLinearAlgebra` is present in your environment when you need this behavior.

## Cross `isapprox` with matrices added

`isapprox` now supports mixed comparisons:

```julia
isapprox(lattice, A)
isapprox(A, lattice)
```
