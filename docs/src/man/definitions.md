# Definitions and conventions

```@contents
Pages = ["definitions.md"]
Depth = 2
```

## Basis vectors

In this package, basis vectors are represented by three-column vectors:

```math
\mathbf{a} = \begin{bmatrix}
    a_x \\
    a_y \\
    a_z
\end{bmatrix},
\quad
\mathbf{b} = \begin{bmatrix}
    b_x \\
    b_y \\
    b_z
\end{bmatrix},
\quad
\mathbf{c} = \begin{bmatrix}
    c_x \\
    c_y \\
    c_z
\end{bmatrix},
```

in Cartesian coordinates. Depending on the situation,
``\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}``
is used instead of
``\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix}``.

Therefore, a lattice is represented as

```math
\mathrm{A} =
\begin{bmatrix} \mathbf{a} & \mathbf{b} & \mathbf{c} \end{bmatrix} =
\begin{bmatrix}
    a_x & b_x & c_x \\
    a_y & b_y & c_y \\
    a_z & b_z & c_z
\end{bmatrix}.
```

## [Left and right matrix actions on a lattice](@id lattice_matrix_actions)

With the convention above, the two operator overloads
`Base.:*(::AbstractMatrix, ::Lattice)` and
`Base.:*(::Lattice, ::AbstractMatrix)` represent different transformations:

- Left action (`R * lattice`): ``\bar{\mathrm{A}} = \mathrm{R}\mathrm{A}``.
  This is an active (reverse/alibi) transformation, e.g. a rigid Cartesian
  rotation that rotates the crystal basis vectors themselves.
- Right action (`lattice * P`): ``\mathrm{A}' = \mathrm{A}\mathrm{P}``.
  This is a passive (forward/alias) transformation, i.e. a change of basis.

In both cases, only ``3\times 3`` matrices are valid, because
`Lattice` is explicitly a 3D object.
Supplying any other matrix shape throws `DimensionMismatch`.

These actions can be composed:

```math
\bar{\mathrm{A}}' = \mathrm{R}\,\mathrm{A}\,\mathrm{P}.
```

In general, left and right actions have different meanings and do not commute.

## [Transformation to the primitive cell](@id primitive)

In the standardized unit cells, there are five different centring
types available, base centrings of A and C, rhombohedral (R), body-centred (I),
and face-centred (F). The transformation is applied to the
standardized unit cell by

```math
\begin{bmatrix} \mathbf{a}_p & \mathbf{b}_p & \mathbf{c}_p \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_s & \mathbf{b}_s & \mathbf{c}_s \end{bmatrix}
\mathrm{P}
```

where ``\mathbf{a}_p``, ``\mathbf{b}_p``, and ``\mathbf{c}_p``
are the basis vectors of the primitive cell and ``\mathrm{P}`` is the
transformation matrix from the standardized unit cell to the primitive
cell. Matrices ``\mathrm{P}`` for different centring types are given as follows:

```math
\mathrm{P}_\text{A} = \begin{bmatrix}
    1 & 0 & 0 \\
    0 & \dfrac{1}{2} & \dfrac{-1}{2} \\
    0 & \dfrac{1}{2} & \dfrac{1}{2}
\end{bmatrix},
\quad
\mathrm{P}_\text{C} = \begin{bmatrix}
    \dfrac{1}{2} & \dfrac{1}{2} & 0 \\
    \dfrac{-1}{2} & \dfrac{1}{2} & 0 \\
    0 & 0 & 1
\end{bmatrix},
\quad
\mathrm{P}_\text{R} = \begin{bmatrix}
    \dfrac{2}{3} & \dfrac{-1}{3} & \dfrac{-1}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{\bar{2}}{3} \\
    \dfrac{1}{3} & \dfrac{1}{3} & \dfrac{1}{3}
\end{bmatrix},
\quad
\mathrm{P}_\text{I} = \begin{bmatrix}
    \dfrac{-1}{2} & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{-1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & \dfrac{-1}{2}
\end{bmatrix},
\quad
\mathrm{P}_\text{F} = \begin{bmatrix}
    0 & \dfrac{1}{2} & \dfrac{1}{2} \\
    \dfrac{1}{2} & 0 & \dfrac{1}{2} \\
    \dfrac{1}{2} & \dfrac{1}{2} & 0
\end{bmatrix}.
```

The choice of transformation matrix depends on the purpose.

For rhombohedral lattice systems with the H setting (hexagonal lattice),
``\mathrm{P}_\text{R}`` is applied to obtain
primitive basis vectors. However, with the R setting (rhombohedral lattice),
no transformation matrix is used because it is already a primitive cell.

## Supercell generation

The basis vectors of unit cell
``\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}``
can be transformed to basis vectors of supercell
``\begin{bmatrix} \mathbf{a}'_1 & \mathbf{a}'_2 & \mathbf{a}'_3 \end{bmatrix}``
by linear transformation

```math
\begin{bmatrix} \mathbf{a}'_1 & \mathbf{a}'_2 & \mathbf{a}'_3 \end{bmatrix} =
\begin{bmatrix} \mathbf{a}_1 & \mathbf{a}_2 & \mathbf{a}_3 \end{bmatrix}
\mathrm{P},
```

just as in [primitive cell transformations](@ref primitive).

Usually, all elements ``P_{ij}`` should be integers satisfying ``i \ne 0``
so that ``\det(\mathrm{P}) \ge 1``.
When ``\det(\mathrm{P}) = 1``, the transformation preserves volume.
However, sometimes we could have ``\det(\mathrm{P}) < 0`` if the transformation
changes the handedness of the lattice.
See [this post](https://gitlab.com/ase/ase/-/issues/938) for more information.

As stated above, a new (super)cell is produced by replication of the initial one over
cell vectors ``\mathbf{a}``, ``\mathbf{b}``, and ``\mathbf{c}``.
The new cell will have a size of
``l \mathbf{a} \times m \mathbf{b} \times n \mathbf{c}``.
Each crystallographic site in the initial cell with Cartesian
coordinates ``\mathbf{r}`` will have a total of ``l m n``
images in the supercell with Cartesian coordinates

```math
\mathbf{r}^{i, j, k} = \mathbf{r} + i \mathbf{a} + j \mathbf{b} + k \mathbf{c},
```

where ``i = 0, 1, \ldots, l - 1``, ``m = 0, 1, \ldots, m - 1``,
and ``k = 0, 1, \ldots, n - 1``
(See [this paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0129-3)).

It is important to keep in mind that this supercell expansion approach
is a special case: the simplest one.
It does not allow, for example, transformations of a primitive cell into a
conventional (super)cell or the opposite. A more general approach exists, which
creates supercell vectors on the basis of linear combinations of the initial cell vectors,
a desirable improvement that will be considered for a future version of the code.
