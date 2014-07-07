import numpy as np

def makeModMat(D):
    """Create a modularity matrix using the values from data matrix,
       where the columns of D represent different variables, and each row
       corresponds to a certain piece of data for a specific variable
       (D should be represented by a 2D numpy array)"""
    modMat = np.zeros(D.shape)  #initialize the mod mat as the same size as D,
    #In actual code we will need to calculate the values of D on the fly
    degrees = np.zeros(D.shape[0]) #initialize a array to hold the degree of each vertex
    for i in np.arange(D.shape[0]):
        degrees[i] = np.sum(D[i, :]) #sum the ith row to get degree of ith vertex
        
    m = np.sum(degrees) #m is equal to twice the number of edges
    for i in np.arange(D.shape[0]):
        for j in np.arange(D.shape[1]):
            modMat[i, j] = D[i, j] - ((degrees[i] * degrees[j]) / m)

    return modMat


def split_cluster(D):
    """Using modularity measures, split the data matrix D into two different
       clusters and return a vector whose ith element gives the cluster#(either 1 or -1)
       that the corresponding vertex on the graph belongs to"""

    modMat = makeModMat(D) #create the modularity matrix
    vals, vects = np.linalg.eig(modMat) #get eigenvalues/vectors of the modularity matrix
    eig_index = np.argmax(vals) #find the index of the largest eigenvalue
    cluster = vects[:, eig_index] #get principle eigenvector
    #cluster = cluster / np.absolute(cluster)
    return cluster
    
if __name__ == "__main__":
    D = np.genfromtxt("test.txt", delimiter=',')
    B = np.dot(D.T, D)
    cluster = split_cluster(B)
    print cluster
        
