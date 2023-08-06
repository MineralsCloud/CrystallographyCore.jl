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
