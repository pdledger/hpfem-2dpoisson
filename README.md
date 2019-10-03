# hpfem-2dpoisson

This is a simple 2D finite element program that employs a hp-finite element discretisation for the solution of the
Poisson equation. Discretisations consisting of triangular or quadrilateral meshes are possible with elements of 
arbitrary order.

The program has been developed as a MATLAB teaching tool to allow students to play with different types of discretisations
and to experiment with h and p refinements for problems with smooth solutions and for problems with singularities associated
with edges and corners. A number of benchmark problems are included for which an analytical solution is known and the program
output the error measured in the L2, H1 and energy norms.

For triangular meshes the code calls the Mesh2d matlab mesh generation tool that has been developed by D. Engwirda

D. Engwirda, Locally-optimal Delaunay-refinement and optimisation-based mesh generation, Ph.D. Thesis, School of Mathematics
and Statistics, The University of Sydney, http://hdl.handle.net/2123/13148, 2014.

D. Engwirda, Unstructured mesh methods for the Navier-Stokes equations, Honours Thesis, School of Aerospace, Mechanical and
Mechatronic Engineering, The University of Sydney, 2005.

https://www.mathworks.com/matlabcentral/fileexchange/25555-mesh2d-delaunay-based-unstructured-mesh-generation

The hierarchic H1 conforming finite element basis functions are based on those proposed in

B. Szabo, I. Babuska Finite Element Analysis, Wiley 1991

S. Zaglmayr High Order Finite Elements for Electromagnetic Field Computation, PhD Thesis, Johannes Kepler University Linz,
Austria, 2006 https://www.numa.uni-linz.ac.at/Teaching/PhD/Finished/zaglmayr
