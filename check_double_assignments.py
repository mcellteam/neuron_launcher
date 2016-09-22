import bpy, bmesh
from mathutils import Vector
import math

# Main

def f_check_double_assignments(context):
    
    print("> Running: f_check_double_assignments")

    # Flag to fix doubles
    FIX_DOUBLES = False
    
    # Check list
    f_check = []
    double_ctr = 0

    # Get the active object
    ob_list = context.selected_objects
    
    if len(ob_list) == 1:
        ob = ob_list[0]
    
        # Get the region list
        reg_list = ob.mcell.regions.region_list

        # Make sure everything is de-selected before we start
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        # Go through the regions
        for reg in reg_list:
        
            print("Checking region: " + str(reg.name))

            # Get the faces
            f_list = reg.get_region_faces(ob.data)

            # Duplicate list if planning to fix things
            if FIX_DOUBLES:
                f_list_new = list(f_list)

            # Go through all faces
            for f in f_list:
                if not f in f_check:
                    f_check.append(f)
                else:
                    print("Face: " + str(f) + " is double - one of the regions is: " + str(reg.name))
                    print("(the other I did not search out)")
                    ob.data.polygons[f].select = True
                    double_ctr += 1
    
                    # Fix?
                    if FIX_DOUBLES:
                        
                        # De-assign
                        del f_list_new[f_list_new.index(f)]
                        reg.set_region_faces(ob.data,f_list_new)
                        
                        print("Fixed double for face: " + str(f) + " is now de-assigned from: " + str(reg.name))

        # Return in edit mode to show any doubles
        bpy.ops.object.mode_set(mode='EDIT')

        print("Finished - number of doubles: " + str(double_ctr))

    print("> Finished: f_check_double_assignments")



