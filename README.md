# FEMDLpy
 A python package for a finite element based discretization of dynamic Laplacians

## Usage
You will require the following packages:
```Python
 import numpy as np
 import matplotlib.pyplot as plt
 import scipy as sp
 from scipy.spatial import Delaunay
 from scipy.sparse import csc_matrix
 from scipy.sparse.linalg import eigsh
 from matplotlib.tri import Triangulation
```

 Suppose you have trajectory data $p$ stored as an numpy array in the format $p[i,j,k] =$ value (in degrees), where $i$ is the index for the timestep (usually in days), $j$ is the index for longitude or latitude (0 for longitude, 1 for latitude), and $k$ is the index for the particle number.



 Then you are ready to construct the stiffness and mass matrix:
 ```Python
 n = p.shape[2]
 nt = p.shape[0]

 #Delaunay options
 options = 'Qt Qbb Qc'

 t = [] #stores our Delaunay simplices

 D = sp.sparse.csc_matrix((n, n))
 M = sp.sparse.csc_matrix((n, n))

 for k in range(0,nt):
     t.insert(k, Delaunay(np.transpose(p[k]), qhull_options=options).vertices)

     # Want to remove simplices with zero volume
     keep = np.ones(len(t[k]), dtype = bool)
     for i, z in enumerate(t[k]):
         if abs(np.linalg.det(np.hstack((np.transpose(p[k])[z], np.ones([1,3]).T)))) < 1E-15:
             keep[i] = False # Point is coplanar, we don't want to keep it
     t[k] = t[k][keep]

     A = np.kron([1, 0, 1], np.ones((len(t[k]),1)))

     [Dt, Mt] = assemble(p[k], t[k], pb, A)

     D = D + Dt
     M = M + Mt
 ```
 Next solve the corresponding eigenproblem
 ```Python
 L,V = eigsh(D, 20, M, sigma=0, which='LM')

 pos = (-L).argsort()
 lam = L[pos]
 ord = V[:,pos]
 ```

 The eigenvectors in $ord$ contain geometric information regarding the finite-time coherent sets.

 ## Reference
 G. Froyland and O. Junge. *Robust FEM-based extraction of finite-time
 coherent sets using scattered, sparse, and incomplete trajectories*.
 arXiv, 2017.
