import numpy as np
import ModCluster as mc
import ModClusterExplicit as mce

def func():
    D = np.genfromtxt("test.txt", delimiter=',')
    A = np.dot(D.T, D)
    modMat = mce.makeInitModMap(A)
    degrees = mce.degSum(A)
 
    clust1 = mce.splitCluster(modMat)
    clust2 = mc.splitCluster(D)
    
    groupMat = mce.makeGroupMat(modMat, clust1[0])
    vals, vects = np.linalg.eig(groupMat) #get eigenvalues/vectors of the modularity matrix
    eig_index = np.argmax(vals) #find the index of the largest eigenvalue
    eigs = vects[:, eig_index] #get principle eigenvector

    eigs2 = mc.modEig(D, clust1[0])
    
    print "Vals\n"
    print vals
    print eig_index
   
    print "REAL ", eigs
    print "Implicit", eigs2



if __name__ == "__main__":
    func()
