import numpy as np

#From src/gradbasis.py
def spherical_gradbasis(particles_spherical,triangulation_spherical):

    # Create a congruent triangular element
    v = np.empty((3,2,triangulation_spherical.shape[0]))

    ab = particles_spherical[:,triangulation_spherical[:,1]] - particles_spherical[:,triangulation_spherical[:,0]]
    ac = particles_spherical[:,triangulation_spherical[:,2]] - particles_spherical[:,triangulation_spherical[:,0]]
    ab_norm = np.linalg.norm(ab, axis=0)
    ac_norm = np.linalg.norm(ac, axis=0)

    cosA = np.abs(np.sum(ab * ac, axis=0) / (ab_norm*ac_norm))
    sinA = np.abs(np.linalg.norm(np.cross(ab, ac, axis=0), axis=0) / (ab_norm*ac_norm))
    l1 = ab_norm
    l2 = ac_norm*cosA
    l3 = ac_norm*sinA

    ## Original triangle
    # A - 0th vertex
    # B - 1st vertex
    # C - 2nd vertex

    ## Congruent triangles
    # A (0,0)
    # B (l1,0) -> l1 = |AB|
    # C (l2 cos(A),l2 sin(A)) -> l2 = |AC|, A = angle(BAC)

    v[0,:,:] = [l2-l1,l3]
    v[1,:,:] = [-l2, -l3]
    v[2,:,:] = [l1, np.zeros(triangulation_spherical.shape[0])]

    #v[0,:,:] = [l1-l2,-l3]
    #v[1,:,:] = [-l1, np.zeros(triangulation_spherical.shape[0])]
    #v[2,:,:] = [l2, l3]

    # Compute area and gradients in \mathbb{R}^{2} of the congruent triangular element
    area = 0.5*(-v[2,0,:]*v[1,1,:] + v[2,1,:]*v[1,0,:])

    dphi = np.empty((2,3,triangulation_spherical.shape[0]))
    dphi[0] = -v[:,1,:]/(2*area)
    dphi[1] = v[:,0,:]/(2*area)
    area = abs(area)

    return dphi, area
