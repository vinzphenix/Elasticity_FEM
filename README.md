# Linear elasticity simulations - 2D FEM

This work is a Finite Element sovler for 2D isotropic linear elasticty:

$$\begin{aligned}
0 &= \nabla \cdot \boldsymbol\sigma(\mathbf{u}) + \mathbf{f}
&& \text{in } \Omega\\
\boldsymbol\sigma &= 2\mu \boldsymbol\epsilon + \lambda \mathrm{tr}(\boldsymbol\epsilon) \boldsymbol\delta\\
\boldsymbol \epsilon &= 
(\nabla \mathbf{u}) + (\nabla \mathbf{u})^* \\
\mathbf{\mathbf{u}} &= \mathbf{u}_D && \text{on } \partial\Omega_D\\
\boldsymbol{\sigma}\cdot\mathbf{n} &= \mathbf{t}_N && \text{on } \partial\Omega_N\\
\boldsymbol{\sigma}\cdot\mathbf{n} + k \mathbf{u} &= \mathbf{t}_R && \text{on } \partial\Omega_R
\end{aligned}$$

The deformation problem reduces to a linear system $K u = b$ when the above is discretized with finite elements.

The elasto-dynamic problem is given below for the continuum and the discrete problem

$$\begin{aligned}
\rho \partial_t^2 \mathbf{u} &= \nabla \cdot \boldsymbol\sigma(\mathbf{u}) + \mathbf{f}\\
M \ddot{u} + Ku &= b
\end{aligned}$$

Using a fourier expansion $U=\Phi \exp(\imath \omega t)$, we find a generalized eigenvalue problem whose solutions are the natural frequencies and shapes of the structure:
$$K \Phi = \omega^2 M \Phi.$$

This solver is capable of solving both deformation and modal analysis problems.

## Launch a simulation

To compile the files, just type ```make```. 

To execute it:
```sh
./deformation <model_name> <nb> <mesh_size_factor>
```
where 
- `<model_name>` must be one of the models in `models/`, 
- `<nb>` is either 
  - $0$ for the static case or
  - $k>0$ to compute the first $k$ eigenmodes, and
- `<mesh_size_factor>` indicates the mesh size target for gmsh

## Structure

    ├── include           <- Header files .h
    ├── models            <- Models with geometry + b.c.
    ├── src               <- Source files .c
    └── tests             <- Test files + benchmark files

## Solver capabilities
- [x] P1/Q1 finite elements
- [x] FEM matrices assembled in packed band storage
- [x] 2D plane stress / 2D plane strain / 3D axisymmetric
- [x] Dirichlet - Neumann - Robin boundary conditions
- [x] $x$ / $y$ or normal / tangent boundary conditions
- [x] RCMK renumbering of nodes
- [x] Band $LDL^*$ linear solver
- [x] Iterative sparse PCG (Jacobi, SSOR, ILU0, ILU1)
- [ ] Assemble in COO storage
- [ ] Impose b.c. in local matrices (better for COO)
- [ ] Write `.msh` files and parse them to avoid `gmsh` memory leaks

The solver was validated against analytic solutions
- in 2D with a clamped rectangular beam
- in 2D and axisymmetric with a hollow cylinder

## Post-processing
The post-processing is done with `Gmsh`. Multiple views are available to 
visualize the deformation:
- Displacement $u$
- Boundary forces
- Stress (Von Mises scalar)
- Stress (Min/Max eigenvalue)
- Stress components (cartesian / polar)

## Dependencies
- `gmsh` to generate the mesh and to realize the post-processing
- `openblas` for linear algebra
- `lapack` for testing

## Write a new model

Write a file `models/<name>.c` with at least 5 functions

1. Geometry instructions for `gmsh`

    ```mesh_<name>(double lc, int e_type)```

    - `lc` is a mesh size factor
    - `e_type` is the element type (triangle / quad)

2. Bulk source term
  
    `set_bk_source_<name>(double rho, double xy[2], double f[2])`

    - `rho` is a scalar for density
    - `xy` is the position
    - `f` contains the bulk force

3. Boundary force term (Neumann/Robin)
    
    ```set_bd_force_<name>(int e, char d, double xy[2], double g[2])```

    - `e` is the entity tag of the boundary curve
    - `d` is the direction (`x`, `y`, `n` or `t`)
    - `xy` is the position
    - `g[0]` contains the stiffness coefficient $k$ of a Robin b.c.
    - `g[1]` contains the force term $\mathbf{t}_R\cdot \mathbf{d}$

4. Imposed boundary displacement (Dirichlet)

    ```set_bd_disp_<name>(int e, char d, double xy[2], double u[1])```
    - `e` is the entity tag of the boundary curve
    - `d` is the direction (`x`, `y`, `n` or `t`)
    - `xy` is the position
    - `u[0]` contains the displacement $\mathbf{u}\cdot \mathbf{d}$

5. Set physical constants

    ```set_physics_<name>(double params[4], Model2D *type)```
    - `params[0]` contains the Young's modulus $[\textrm{Pa}]$
    - `params[1]` contains the Poisson coefficient $[-]$
    - `params[2]` contains the density $[\textrm{kg} \: \textrm{m}^{-3}]$
    - `params[3]` contains a reference length $[\textrm{m}]$
    - `type[0]` indicates plane stress / plane strain / axisymmetric
   
Add your functions in `src/model.c` so that they can be selected at run time.

