import numpy as np
from scipy.sparse import lil_matrix, coo_matrix


#From src/assemble.py
def assemble(p,t,pb=None,G=None, spherical=False):

    nparticles = p.shape[1] ## Number of particles
    nt = t.shape[0] ## Number of elements

    if spherical == True:
        dphi, area = spherical_gradbasis(p,t)
    else:
        dphi, area = gradbasis(p,t)

    PG = np.transpose(G[:,:,:], (2,1,0))

    D = lil_matrix((nparticles, nparticles))
    M = lil_matrix((nparticles, nparticles))


    for i in range(3):
        for j in range(i,3):
            Dij = -area*(dphi[0,i,:] * PG[:,0,0] * dphi[0,j,:] +
                         dphi[0,i,:] * PG[:,0,1] * dphi[1,j,:] +
                         dphi[1,i,:] * PG[:,1,0] * dphi[0,j,:] +
                         dphi[1,i,:] * PG[:,1,1] * dphi[1,j,:])
            Mij = area/12.

            I = t[:,i]
            J = t[:,j]

            if (j == i):
                D = D + coo_matrix((Dij,(I,J)), shape=(nparticles,nparticles))
                M = M + coo_matrix((Mij + area/12.,(I,J)), shape=(nparticles,nparticles))
            else:
                D = D + coo_matrix((Dij,(I,J)), shape=(nparticles,nparticles))
                D = D + coo_matrix((Dij,(J,I)), shape=(nparticles,nparticles))

                M = M + coo_matrix((Mij,(I,J)), shape=(nparticles,nparticles))
                M = M + coo_matrix((Mij,(J,I)), shape=(nparticles,nparticles))

    return [D.tocsr(),M.tocsr()]
