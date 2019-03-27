# High-Order Curved Mesh Generation - On the use of Radial Basis Functions (RBF) for mesh movement.

## Objective: 
The code has the objective of:
1. Read a .cgns file as a mesh input (TODO: specify part of the code where it happens)
2. Read a .iges file as a geometric input (TODO: same as above)
3. Combine both into high-order (curved) mesh. 
4. Write down the mesh in a recognized pattern (currently GMSH)
5. Use Radial Basis Functions for completing the task of moving the curved mesh in an efficient way. 
- It makes use of the LAPACK library for optimal speed.
