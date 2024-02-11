# FEA_1D_BAR_ELEMENT
Meshing - To solve the problem numerically by using the finite element method, we first need to MESH the 1D structure into several nodes and elements. I have meshed it such that it has 4 elements and 5 nodes with both end nodes fixed(using the Dirichlet boundary condition) and load applied at the 2nd element 3rd node(Using Neumann boundary condition).

Code Structure:

1. PreProcessing: Read input files
2. Loop over elements to construct the global stiffness
3. Loop over nodes to construct the global force vector
4. Apply boundary condition 
5. Solve system for nodal displacement
6. PostProcessing: plot stress distribution, strain field, etc
