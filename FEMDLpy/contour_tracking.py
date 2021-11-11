import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath

def contour_tracking(initial = None, final = None, contours_to_track = None, connected=False, max_point_initial=None, max_point_final=None, return_particles=False):
    """
    Given an initial dataset and a final dataset, find the contour value in the
    final dataset that best matches the given contours in contours_to_track.
    The best matching contour value is the line integral along the given contour
    divided by it's length.
    """

    # Initial data
    i_normed_eigenfunction = initial['eigenfunction']
    i_particles = initial['particles']
    i_domain = initial['domain']
    if initial['triangulation'] is None:
        i_tri = Triangulation(i_particles[0,i_domain], i_particles[1,i_domain])
    else:
        i_tri = initial['triangulation']

    # Final data
    f_normed_eigenfunction = final['eigenfunction']
    f_particles = final['particles']
    f_domain = final['domain']

    if final['triangulation'] is None:
        f_tri = Triangulation(f_particles[0,f_domain], f_particles[1,f_domain])
    else:
        f_tri = final['triangulation']

    # Refresh figures
    plt.close('all')

    # Find maximum particle
    #max_particle_id = i_domain[np.argmax(i_normed_eigenfunction)]

    # Create contours
    contours = plt.tricontour(i_tri, i_normed_eigenfunction, levels=contours_to_track)
    if max_point_initial is not None:
        contours_positions = contour_containing_point(contours, max_point_initial)#[i_particles[0,max_particle_id], i_particles[1,max_particle_id]])
        single_contour=True
    else:
        contours_positions = contour_containing_point(contours) # No max point given
        single_contour=False

    # Pull back the final eigenvector to the initial particle positions
    tri_pullback = Triangulation(i_particles[0,f_domain], i_particles[1,f_domain])
    function_evaluator = LinearTriInterpolator(tri_pullback, f_normed_eigenfunction)

    # Interpolate the pulled back eigenvector onto to the contour
    # Also determine particle ids inside contour (not simple value > contour due to possible multiple contours)
    tracked_contours = []
    for i in range(len(contours_to_track)):
        #Compute distance between contour points and new contour value
        if single_contour:
            if len(contours_positions[i]==0):
                new_contour_value = 1 # i.e. there is no curve that satisfies.
            else:
                interpolated_contour_points = function_evaluator(contours_positions[i][:,0], contours_positions[i][:,1])
                distances = distance(contours_positions[i][:-1,0],contours_positions[i][:-1,1],contours_positions[i][1:,0],contours_positions[i][1:,1])
                line_integral = np.sum(moving_average(interpolated_contour_points)*distances)
                new_contour_value = line_integral/np.sum(distances)
        else:
            distances_total = 0
            line_integral = 0
            for j in range(len(contours_positions[i])):
                interpolated_contour_points = function_evaluator(contours_positions[i][j][:,0], contours_positions[i][j][:,1])
                distances = distance(contours_positions[i][j][:-1,0],contours_positions[i][j][:-1,1],contours_positions[i][j][1:,0],contours_positions[i][j][1:,1])
                line_integral += np.sum(moving_average(interpolated_contour_points)*distances)
                distances_total += np.sum(distances)
            new_contour_value = line_integral/distances_total

        # Save new contour value
        tracked_contours.append(new_contour_value)

    if return_particles:
        # Compute the new set of touched particles in the final contours
        touched_particles = []
        contours = plt.tricontour(f_tri, f_normed_eigenfunction, levels=tracked_contours)
        if max_point_final is not None:
            contours_positions = contour_containing_point(contours, max_point_final)#[i_particles[0,max_particle_id], i_particles[1,max_particle_id]])
            single_contour=True
        else:
            contours_positions = []
            for i in range(len(contours.collections)):
                for j in range(len(contours.collections[i].get_paths())):
                    contour = contours.collections[i].get_paths()[j].vertices
                    contours_positions.append(contour)

            single_contour=False

        for i in range(len(tracked_contours)):
            # Compute particles touched
            if single_contour:
                domain_particles = f_domain[f_normed_eigenfunction>=tracked_contours[i]]
                contour_path = mpltPath.Path(contours_positions[i])
                f_touched_particles = domain_particles[contour_path.contains_points(np.array([f_particles[0,domain_particles], f_particles[1,domain_particles]]).T)]
                touched_particles.append(f_touched_particles)
            else:
                domain_particles = f_domain[f_normed_eigenfunction>=tracked_contours[i]]
                f_touched_particles = np.array([])
                for j in range(len(contours_positions[i])):
                    if len(contours_positions[i][j]) < 3:
                        continue
                    contour_path = mpltPath.Path(contours_positions[i][j])
                    current_touched_particles = domain_particles[contour_path.contains_points(np.array([f_particles[0,domain_particles], f_particles[1,domain_particles]]).T)]
                    f_touched_particles = np.hstack(f_touched_particles, current_touched_particles)
                touched_particles.append(f_touched_particles)
        return tracked_contours, touched_particles
    else:
        return tracked_contours
