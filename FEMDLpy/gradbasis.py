import numpy as np

#From src/gradbasis.py
def gradbasis(particles,triangulation):

    v = np.empty((3,2,triangulation.shape[0]))
    v[0] = particles[:,triangulation[:,2]] - particles[:,triangulation[:,1]]
    v[1] = particles[:,triangulation[:,0]] - particles[:,triangulation[:,2]]
    v[2] = particles[:,triangulation[:,1]] - particles[:,triangulation[:,0]]

    area = 0.5*(-v[2,0,:]*v[1,1,:] + v[2,1,:]*v[1,0,:])

    dphi = np.empty((2,3,triangulation.shape[0]))
    dphi[0] = -v[:,1,:]/(2*area)
    dphi[1] = v[:,0,:]/(2*area)
    area = abs(area)

    return dphi, area
