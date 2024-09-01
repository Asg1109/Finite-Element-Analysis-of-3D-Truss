import numpy as np
from tabulate import tabulate
from loadfile import *
import copy

def printM(Matrix):
    M=Matrix.tolist()
    print(tabulate(M,tablefmt="plain"))
    print("\n")
def printL(List2D):
    print(tabulate(List2D,tablefmt="plain"))
    print("\n")

def Km_local(A,E,L):
    size = 2 # size of stiffness matrix
    K = np.matrix(np.zeros((size, size), dtype=float))
    K[0,0] = 1
    K[0,1] = -1
    K[1,0] = -1
    K[1,1] = 1
    K=(A*E/L)*K
    # printM(K)
    return K

def T_matrix(x1,y1,z1, x2,y2,z2,L):
    # nx1 = unit vector at node1 in nodal x direction
    # ny1 = unit vector at node1 in nodal y direction
    # nx2 = unit vector at node2 in nodal x direction
    # ny2 = unit vector at node2 in nodal y direction
    lambda_x = (x2-x1)/L
    lambda_y = (y2-y1)/L
    lambda_z = (z2-z1)/L

    T = np.matrix(np.zeros((2,6), dtype=float))
    T[0,0] = lambda_x
    T[0,1] = lambda_y
    T[0,2] = lambda_z

    T[1,3] = lambda_x
    T[1,4] = lambda_y
    T[1,5] = lambda_z
    return T

def Km_global(A,E, x1,y1,z1, x2,y2,z2):
    L = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    # KL local stiffness matrix
    Kle = Km_local(A,E,L)
    T  = T_matrix( x1,y1,z1, x2,y2,z2, L)
    # Kg global stiffness matrix
    Kge = T.T * Kle * T
    return Kge

def Kms_matrix(nodes, elements,noe, nldof, ngdof):
    Ks=(np.matrix(np.zeros((ngdof,ngdof), dtype=float)))
    for element in elements:
        node1 = (element[0]-1)
        node2 = (element[1]-1)
        x1 = nodes[node1][0]
        y1 = nodes[node1][1]
        z1 = nodes[node1][2]

        x2 = nodes[node2][0]
        y2 = nodes[node2][1]
        z2 = nodes[node2][2]
        
        # Ae Ee Ie : area elasticity and inertia of element
        Ae, Ee= element[2],element[3]
        # Ke element stiffness matrix global
        Kge = Km_global(Ae,Ee, x1,y1,z1, x2,y2,z2)
        # n = 2, number of node in each element 
        # ngdof = number of total global dof in structure 
        Ge = gather_operator(element,noe,nldof,ngdof)
        Ks=Ks+Ge.T*Kge*Ge
    

    export("Output/Kms.csv",range(1,len(Ks)+1), Ks.tolist())
    return Ks

def gather_operator(element,noe,nldof, ngdof):
    # n number of nodes in an element
    # nldof number of local dof at any node
    # ngdof number of global dof
    # nodei = element[i]
    
    G = np.matrix(np.zeros((noe*nldof,ngdof), dtype=int))
    for k in range(noe):
        node_k = element[k]
        for ith_ldof in range(nldof): # range = (0,1,2,3,4,5) 
            jth_gdof=int(nldof*(node_k-1)+ith_ldof)
            G[k*nldof+ith_ldof, jth_gdof] = 1

    return G
        
def force_vector(bc_f,nldof, ngdof):
    # f= force vector
    f = np.matrix(np.zeros(ngdof)).T
    if not bc_f[0]:
        return f
    else:
        for force in bc_f:
            node = force[0]
            for ith_dof in range(nldof):
                if force[ith_dof+1]!=0:
                    j = nldof*(node-1) + ith_dof 
                    f[j] = force[ith_dof+1]
    return f

def disp_vector(bc_dv,nldof, ngdof):
    # u  displacement vector
    u = np.matrix(np.zeros(ngdof)).T

    if not bc_dv[0]:
        return u
    else:
        for displacment in bc_dv:
            node = displacment[0]
            for ith_dof in range(nldof):
                j = nldof*(node-1) + ith_dof 
                u[j] = displacment[ith_dof+1]
    return u

def deletionapproach(bc_dr, nldof, ngdof ):
    # K Structures stiffess matrix
    # bc_dr Boundary condition displacement restraint (Support)
    # bc_dv Boundary condition displacement values (Settlement)
    # bc_f  Boundary condition force (force-vector)
    I=np.matrix(np.eye(ngdof))
    D1=np.matrix(np.eye(ngdof))
    D2=np.matrix(np.eye(ngdof))
    for support in bc_dr:
        node = support[0]
        for ith_dof, restraint in zip(range(3) ,support[1:]):
            if restraint == 1:
                j = nldof*(node-1) + ith_dof 
                D1[j,j] = 0
    D2=I-D1
    
    return D1, D2

def element_end_rections(nodes, elements,bc_dr,bc_dv, bc_f, noe, nldof, ngdof, iter):
    # element end reactions in local coordinate system
    # R Reaction
    # 1 and 2 implies local node
    # element-er = [[R1xt, R1yt, R1zr, R2xt, R2yt, R2zr]]
    element_er=[]
    f=force_vector(bc_f,nldof,ngdof)
    # u = displacement vector
    u=disp_vector(bc_dv,nldof,ngdof)
    ug,rg = u_and_r_all(nodes, elements,bc_dr,u,f, noe, nldof, ngdof) 
    # ug,rg = iteration(nodes, elements,bc_dr,bc_dv, bc_f, noe, nldof, ngdof,iter)
    for element in elements:
        node1 = (element[0]-1)
        node2 = (element[1]-1)
        x1 = nodes[node1][0]
        y1 = nodes[node1][1]
        z1 = nodes[node1][2]

        x2 = nodes[node2][0]
        y2 = nodes[node2][1]
        z2 = nodes[node2][2]

        # length of element
        Le = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        Ae, Ee = element[2],element[3]
        Kle = Km_local(Ae,Ee,Le)
        # Ge gather operator
        Ge = gather_operator(element,noe,nldof,ngdof)
        T = T_matrix(x1,y1,z1, x2,y2,z2, Le)
        # ue nodal displacement at end of elements in global coordinate system
        uge = Ge*ug

        # q memeber end reactions in local coordinate system
        q = Kle*T*uge 
        element_er.append((q.T).tolist()[0])
        # printL(element_er)
    export('Output\Element end Reaction.csv',["F1x","F2x"],element_er)
    return element_er

def u_and_r_all(nodes, elements,bc_dr, u, f, noe, nldof, ngdof):
    Ks = Kms_matrix(nodes, elements, noe, nldof, ngdof)

    D1, D2= deletionapproach(bc_dr, nldof,ngdof)
    # printM(D1)

    K11=D1*Ks*D1+D2
    K12=Ks*D2-D2*Ks*D2 + D2
    # f = force vector
    f1 = D1*f
    u2 = D2*u
    # u3 new displacement vector
    u3=np.linalg.solve(K11,(f1-K12*u2))
    r =Ks*u3-D2*f
    
    export(r"Output\U.csv",["Displacement"],u3.tolist())
    export(r"Output\F.csv",["Force"],r.tolist())
    return u3, r
