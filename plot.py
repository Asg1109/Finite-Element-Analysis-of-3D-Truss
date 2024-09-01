import numpy as np
import matplotlib.pyplot as plt
# works only for 2D frame with 3 nldof
def plot_truss(nodes,elements, bc_dr, bc_f):

    # Creating the figure and the 3D axis
    # fig = plt.figure()
    fig = plt.figure() 
    ax = fig.add_subplot(111, projection='3d')

    # finding smallest length of element
    le = 10e9
    for element in elements:
        node1 = (element[0]-1)
        node2 = (element[1]-1)
        x1 = nodes[node1][0]
        y1 = nodes[node1][1]
        z1 = nodes[node1][2]

        x2 = nodes[node2][0]
        y2 = nodes[node2][1]
        z2 = nodes[node2][2]
          
        Le = np.sqrt((x2-x1)**2+(y2-y1)**2 + (z2-z1)**2)
        if Le < le:
            le=Le

    # plotting original nodes
    nodes = np.array(nodes)
    x,y,z = nodes[:,0], nodes[:,1], nodes[:,2]
    ax.scatter(x,y,z, c='r', zorder=5)

    # plotting displacement boundary conditions  (restraint)
    
    if bc_dr[0]:
        for support in bc_dr:
            node = support[0]-1
            p=0.1
            if support[1]==1: 
                ax.plot([nodes[node][0],nodes[node][0]+ p*le], [nodes[node][1],nodes[node][1]],[nodes[node][2],nodes[node][2]], c='black',linewidth='4', zorder=5)
            if support[2]==1:
                ax.plot([nodes[node][0],nodes[node][0]], [nodes[node][1],nodes[node][1]+ p*le],[nodes[node][2],nodes[node][2]], c='black',linewidth='4', zorder=5)
            if support[3]==1:
                ax.plot([nodes[node][0],nodes[node][0]], [nodes[node][1],nodes[node][1]],[nodes[node][2],nodes[node][2]+ p*le], c='black',linewidth='4', zorder=5)
    

   # Plotting Elements
    for element in elements:
        node1 = (element[0]-1)
        node2 = (element[1]-1)
        x1 = nodes[node1][0]
        y1 = nodes[node1][1]
        z1 = nodes[node1][2]

        x2 = nodes[node2][0]
        y2 = nodes[node2][1]
        z2 = nodes[node2][2]

        ax.plot( [x1,x2], [y1,y2], [z1,z2], c='skyblue', zorder=1)
            
    plt.tight_layout()
    plt.axis('equal') 
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()


