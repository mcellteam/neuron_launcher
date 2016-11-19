import os

from math import *

# Import an MDL file geometry as a list
def import_mdl(fname, obj_name):

    vert_list = []
    el_list = []
    reg_dict = {}

    # Flag whether we are reading data for the correct object in the file
    READ_OBJ = False
    # Flag to read: (vertices, element connections, surface regions)
    READ_FLAGS = (False, False, False)
    # Current region name
    CURR_REG_NAME = ""

    f = open(fname, 'r')
    for line in f:
        line_s = line.split()
        if len(line_s) > 0:
    
            # Check that it's the correct object
            if len(line_s) == 2 and line_s[1] == 'POLYGON_LIST':
                if line_s[0] == obj_name:
                    READ_OBJ = True
                    continue
                else:
                    READ_OBJ = False
                    continue
            elif len(line_s) == 1 and line_s[0] == 'VERTEX_LIST':
                READ_FLAGS = (True, False, False)
                continue
            elif len(line_s) == 1 and line_s[0] == 'ELEMENT_CONNECTIONS':
                READ_FLAGS = (False, True, False)
                continue
            elif len(line_s) == 1 and line_s[0] == 'DEFINE_SURFACE_REGIONS':
                READ_FLAGS = (False, False, True)
                continue

            # Do the reading
            
            # Verts
            if READ_FLAGS[0] == True:
                if len(line_s) == 5:
                    vert_list.append((float(line_s[1][:-1]),float(line_s[2][:-1]),float(line_s[3])))
                continue
            
            # Element connections
            if READ_FLAGS[1] == True:
                if len(line_s) == 5:
                    el_list.append((int(line_s[1][:-1]),int(line_s[2][:-1]),int(line_s[3])))
                continue

            # Region defs
            if READ_FLAGS[2] == True:
                if len(line_s) == 1 and line_s[0] != '{' and line_s[0] != '}':
                    CURR_REG_NAME = line_s[0]
                elif len(line_s) > 0 and line_s[0] == 'ELEMENT_LIST':
                    vals = []
                    vals.append(int(line_s[2][1:-1]))
                    vals += [int(line_s[i][:-1]) for i in range(3,len(line_s))]
                    reg_dict[CURR_REG_NAME] = vals
                    
    f.close()
                
    return vert_list, el_list, reg_dict

