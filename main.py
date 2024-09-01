import numpy as np
from loadfile import*
from matrixfunc import*

from plot import*
# Example
    # Truss_1, Truss_2 are folders that can be used for input
nodes, elements, bc_dr, bc_dv, bc_f=load_all("Truss_1")
noe = 2 # number of node at each element
nldof = 3 # number of degree of freedom (DOF) at each node
ngdof = nldof*len(nodes) # Total number of DOF

er=element_end_rections(nodes, elements,bc_dr,bc_dv, bc_f, noe, nldof, ngdof, 10)
# ouput is in ouput folder
plot_truss(nodes,elements,bc_dr,bc_f)
