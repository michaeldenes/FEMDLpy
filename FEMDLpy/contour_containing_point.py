import numpy as np

def contour_containing_point(contours, max_point=None):
    """
    Compute the closed contour containing a point. Useful in cases where there
    are two peaks to a function, causing disconnected contour sets.
    """
    if max_point is None:
        contours_positions = []
        for i in range(len(contours.collections)):
            current_contours_positions = []
            for j in range(len(contours.collections[i].get_paths())):
                current_contours_positions.append(contours.collections[i].get_paths()[j].vertices)
            contours_positions.append(current_contours_positions)

    else:
        contours_positions = []
        for i in range(len(contours.collections)):
            for j in range(len(contours.collections[i].get_paths())):
                path = mpltPath.Path(contours.collections[i].get_paths()[j].vertices)
                if max_point is None:
                    if path.contains_point([i_particles[0,max_particle_id], i_particles[1,max_particle_id]]):
                        contour = contours.collections[i].get_paths()[j].vertices
                        contours_positions.append(contour)
                        break # Breaks the j loop
                elif path.contains_point(max_point):
                    contour = contours.collections[i].get_paths()[j].vertices
                    contours_positions.append(contour)
                    break # Breaks the j loop
                elif j == len(contours.collections[i].get_paths())-1:
                    contours_positions.append(np.array([])) # Pass back an empty array if the max particle doesn't exist in a given contour.


    return contours_positions
