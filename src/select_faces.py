import bpy, bmesh

# Main

if __name__ == "__main__":

    print("> Running")
    
    fs = [1956-1,2110-1]

    # Get the active object
    ob_list = bpy.context.selected_objects
    
    if len(ob_list) == 1:
        ob = ob_list[0]

        # Deelect all edges
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.context.tool_settings.mesh_select_mode = (False, False, True)
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        # Select
        ob.data.polygons[fs[0]].select = True
        ob.data.polygons[fs[1]].select = True

        # Return in edit mode
        bpy.ops.object.mode_set(mode='EDIT')