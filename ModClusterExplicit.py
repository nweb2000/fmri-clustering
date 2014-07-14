import numpy as np

def degSum(A):
    """
    Find the degree sums of adjacecy matrix A (this is equal to the row sums)
    """
    degrees = np.zeros(A.shape[0]) #initialize a array to hold the degree of each vertex
    for i in np.arange(A.shape[0]):
        degrees[i] = np.sum(A[i, :]) #sum the ith row to get degree of ith vertex

    return degrees

def modularity(modMat, clust):
    """
    Calculate the modularity value obtained when the clustering described in clust is applied
    on the graph whose modularity is modMat
    clust should be a vector with values of either 1 or -1, where the sign indicates which
    cluster a element is in
    """
    #calculate clust * modMat * clust_T
    return np.dot(np.dot(clust, modMat), clust.T)
    
def makeInitModMap(A, degrees=None):
    """
    Create a modularity matrix using the values from Adjacency Matrix
        A  (A should be represented by a 2D numpy array)
        @INPUTS 
        (2D numpy array) A - the adjacency matrix of the graph
        (numpy array) degrees - a vector whose ith element holds the sum of the degrees of the 
        ith row of the adjacency matrix, optional
        @OUTPUTS
            A modularity matrix with the same dimensions as A     
    """
    
    if degrees == None:
        degrees = degSum(A)
            
    modMat = np.zeros(A.shape)  #initialize the mod mat as the same size as D,    
    m = np.sum(degrees) #m is equal to twice the number of edges
    
    for i in np.arange(A.shape[0]):
        for j in np.arange(A.shape[1]):
            modMat[i, j] = A[i, j] - ((degrees[i] * degrees[j]) / m)

    return modMat

def makeGroupMat(modMat, clustList):
    """
    Create a modularity matrix for a group of graph nodes specified in clustList
    @INPUTS
    (2D numpy array) modMat - the modularity matrix for the entire graph
    (numpy array) clustList - a list of indices for the group of graph nodes
    """
    groupMat = np.zeros((clustList.size, clustList.size))

    for i in np.arange(clustList.size):
        for j in np.arange(clustList.size):
            if i == j:
                b = 0
                for k in np.arange(clustList.size):
                   b = b + modMat[clustList[i], clustList[k]] #see Newman's modularity paper for details
                   
                groupMat[i, j] = modMat[clustList[i], clustList[j]] - b
            else:
                groupMat[i, j] = modMat[clustList[i], clustList[j]]

    return groupMat


def splitCluster(modMat=None, clustList=None, A=None):
    """
    Split the group of graph nodes (whose indices in the adjacency matrix are listed in clustList), into
    two distinct groups
    NOTE: Passing in only modMat implies a first time split
    NOTE: If only the adjacency matrix is passed in, it will generate a modularity matrix
    @INPUT
    (numpy array) clustList - A list of indices for the group of graph nodes to split
    (2D numpy array) modMat - the modularity matrix for the graph
    (2D numpy array) A - the adjacency matrix
    @OUTPUT
        Two lists of graph node indices which show how clustList was split.
        These are stored in a python tuple.
        TODO: If no optimal split can be found, returns None.

    TODO: Need to check clustList to make sure that it is not empty, and other general error checking
    """
    if modMat==None and A != None:
        modMat = makeInitModMat(A)
    
    if clustList == None:  #if this is a first time split
        clustList = np.arange(modMat.shape[0])
        B = modMat
    else:
        B = makeGroupMat(modMat, clustList)
       
    vals, vects = np.linalg.eig(B) #get eigenvalues/vectors of the modularity matrix
    eig_index = np.argmax(vals) #find the index of the largest eigenvalue
    p_eig = vects[:, eig_index] #get principle eigenvector
   
    #print modularity(B, p_eig / np.abs(p_eig))

    clust_1 = np.array([clustList[x] for x in np.arange(p_eig.size) if p_eig[x] > 0])
    clust_2 = np.array([clustList[x] for x in np.arange(p_eig.size) if p_eig[x] < 0])
    clusters = (clust_1, clust_2)

    
    return clusters
    
if __name__ == "__main__":
    A = np.genfromtxt("test.txt", delimiter=',')
    modMat = makeInitModMap(A)
    degrees = degSum(A)

    clust1 = splitCluster(modMat)

    clust3 = splitCluster(modMat, clust1[0])
    clust4 = splitCluster(modMat, clust1[1])

    print clust1

    print "clust3:\n", clust3
    print "clust4:\n", clust4
    
    
