import numpy as np
import scipy as sp

def SEBA(V, Rinit=None, sort=True, maxiter=5000, accur=1e-14):
    """
    Compute the Sparse Eigenbasis approximation of a set of linearly
    independent vectors.
    """
    V,_ = sp.linalg.qr(V, mode='economic') #maybe transpose?
    p,r = V.shape
    mu = 0.99/np.sqrt(p)

    S = np.zeros((V.shape))

    #Perturb near-constant vectors
    for j in range(r):
        if np.max(V[:,j]) - np.min(V[:,j]) < 1e-14:
            V[:,j] = V[:,j] + (np.random.rand(p,1)-0.5)*1e-12


    #Initialise rotation
    if Rinit is None:
        Rnew = np.eye(r)
    else:
        P,_,Q = sp.linalg.svd(Rinit)
        #Rinit = np.matmul(P,Q)
        Rinit = P@Q
        Rnew = Rinit
    R = 0

    iteration=0

    while np.linalg.norm(Rnew - R,2) > accur and iteration < maxiter:
        iteration += 1
        R = Rnew
        Z = V@R.T

        #Threshold to solve sparse approximation problem
        for i in range(r):
            X = softthresh(Z[:,i], mu)
            S[:,i] = X/np.linalg.norm(X,2)

        #Polar decomposition to solve Procrustes problem
        #P,_,Q = sp.linalg.svd(S.T@V)
        #Rnew = P@Q
        Rnew = sp.linalg.polar(S.T@V)[0]


    #Choose correct parity of vectors and scale so largest value is 1
    for i in range(r):
        S[:,i] = S[:,i]*np.sign(np.sum(S[:,i]))
        S[:,i] = S[:,i]/np.max(S[:,i])

    #Sort eigenvectors
    if sort == True:
        rel = np.amin(S, 0);
        relord = np.argsort(rel)[::-1]
        #relord = np.argsort(rel, kind='stable')[::-1] # requires numpy 1.15+
        S = S[:,relord]

    Su = np.maximum(S,0)

    return S, Su, Rnew
