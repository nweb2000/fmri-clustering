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

    clust3 = mce.splitCluster(modMat, clust1[0])
    clust4 = mce.splitCluster(modMat, clust1[1])
    clust5 = mc.splitCluster(D, clust2[0])
    clust6 = mc.splitCluster(D, clust2[1])
    
    clust7 = mce.splitCluster(modMat, clust3[0])
    clust8 = mce.splitCluster(modMat, clust3[1])
    clust9 = mc.splitCluster(D, clust5[0])
    clust10 = mc.splitCluster(D, clust5[1])

    
  
if __name__ == "__main__":
    func()
