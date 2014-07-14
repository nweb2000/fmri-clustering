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
    print "vals", np.linalg.norm(groupMat, 1)
    vals, vects = np.linalg.eig(groupMat) #get eigenvalues/vectors of the modularity matrix
    eig_index = np.argmax(vals) #find the index of the largest eigenvalue
    eigs = vects[:, eig_index] #get principle eigenvector

    eigs2 = mc.modEig(D, clust1[0])
    print vals[eig_index]

   # print vals[eig_index]
   # print vals[eig_index] / eigs2[1]

    

    
    print "REAL ", eigs
    print "Implicit", eigs2
    print (eigs / abs(eigs)) / (eigs2[0] / abs(eigs2[0]))

    print "TEST\n\n"
    print np.dot(groupMat, eigs2[0]) / (eigs2[0] * eigs2[1])
    print eigs2[0] * eigs2[1]

    



if __name__ == "__main__":
    func()
