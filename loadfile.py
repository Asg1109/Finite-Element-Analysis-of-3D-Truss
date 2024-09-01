import numpy as np
import csv

# writing to csv file
def format_number(value, band):
    if value < band[0]:
        return '{:.2e}'.format(value)  # Format in scientific notation with 2 decimal places
    if value > band[1]:
        return '{:.2e}'.format(value)  # Format in scientific notation with 2 decimal places

    else:
        return str(value)

def export(filename, fields, data):
    with open(filename, 'w',newline='') as csvfile:
        # creating a csv writer object
        csvwriter = csv.writer(csvfile)
    
        # writing the fields
        csvwriter.writerow(fields)
    
        # writing the data rows
        # csvwriter.writerows(rows)
        for row in data:
            formatted_row = [format_number(value, [-1e-3, 1e+4]) for value in row]
            csvwriter.writerow(formatted_row)

def load_nodes( path, delimiter=',', skiprows= 1  ):
    data = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows)
    data=data.tolist()
    if ((np.matrix(data)).shape)[0] == 1:
        return [data]
    return data


def load_elements( path, delimiter=',', skiprows= 1  ):
    data=[]
    dtype = [('column1', int), ('column2', int), ('column3', float),('column4', float)]
    data = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows, dtype=dtype)
    data=data.tolist()
    if ((np.matrix(data)).shape)[0] == 1:
        return [data]
    return data

def load_bc_disp_restrain( path, delimiter=',', skiprows= 1  ):
    # functon to loadd boundary condition for displacement
    # t translation
    # r rotation
    # xt translation in x direction
    # yt translation in y direction 
    # zt translation in z direction 
    data = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows, dtype=int)
    data=data.tolist()
    if ((np.matrix(data)).shape)[0] == 1:
        return [data]
    return data

def load_bc_disp_value( path, delimiter=',', skiprows= 1  ):
    dtype = [('column1', int), ('column2', float), ('column3', float),('column4', float)]
    data = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows, dtype=dtype)

    data=data.tolist()
    if ((np.matrix(data)).shape)[0] == 1:
        return [data]
    return data

def load_bc_force(path, delimiter=',', skiprows= 1  ):
    dtype = [('column1', int), ('column2', float), ('column3', float),('column4', float)]
    data = np.loadtxt(path, delimiter=delimiter, skiprows=skiprows, dtype=dtype)
    data=data.tolist()
    if ((np.matrix(data)).shape)[0] == 1:
        return [data]
    return data

def load_all(truss_name):

    nodes = load_nodes(f"{truss_name}\{'nodes.csv'}")
    elements = load_elements(f"{truss_name}\{'elements.csv'}")
    bc_dr = load_bc_disp_restrain(f"{truss_name}\{'bc_disp_restraint.csv'}") # displacement Boundary Condition for restraint
    bc_dv = load_bc_disp_value(f"{truss_name}\{'bc_disp_value.csv'}")
    bc_f = load_bc_force(f"{truss_name}\{'bc_force_nodal.csv'}")
    return nodes, elements, bc_dr, bc_dv, bc_f