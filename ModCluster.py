import numpy as np

def modMat(D, i, j, ki, kj, m):
    """THIS FUNCTION IS FOR CONVENIENCE, AND IS NOT ACTUALL NEEDED
       Get the (i, j) index of the modularity matrix created by data matrix D.
       @INPUT: D (a data matrix),
               i, j (Index of modularity matrix to access)
               ki, kj (degrees of vertex i, and vertex j)
               m (the sum of degrees for every vertex, which equals 2 * number of graph edges)"""
    
    return np.dot(D[:, i], D[:, j]) - ((ki * kj) / m)

def implicitDegSum(D):
    """Finds the sums of the degrees of each node in the graph of the implicitly represented adjacency matrix of data matrix D"""
    degrees = np.zeros(D.shape[1])
    rowsums = np.zeros(D.shape[0]) 
    for i in np.arange(D.shape[0]): #find the sum of each row and save it 
        rowsums[i] = np.sum(D[i, :])
        
    for i in np.arange(D.shape[1]):
        degrees[i] = np.dot(D[:, i], rowsums)

    return degrees

def modEig(D, maxIter=50, start=None):
    """Find the principle eigenvector for the implicitly represented modularity matrix of
       data matrix D using the power method
       @INPUT: D (data matrix where the columns of D represent different nodes and the 
                  rows represent a piece of time series data for the corresponding node)
               maxIter (how much iterations of power method to go through, this is currently just for testing
                        we will use some kind of tolerence level later on...)
               start (a vector to start the power method with, optional) """
       
    #set x to the first estimate of the principle eigen vector
    if start == None:
        x = np.zeros(D.shape[1])
        x[0] = 1 #initialize x_0 as a unit vector, with all elements = to zero except for the first
    else:
        x = start 
   
    degrees = implicitDegSum(D)    
    m = np.sum(degrees) #m is 2 * the number of edges (or the sum of all edge weights)        
    iterations = 0
    
    while iterations < maxIter:  #do the power method for maxIter iterations
        dx = np.dot(D, x) #get D . X
        eigEst = np.zeros(D.shape[1])
        
        for i in np.arange(D.shape[1]):  #calculate eigEst
            alpha = (degrees[i] / m) * degrees #get alpha vector
            eigEst[i] = np.dot(D[:, i], dx) - np.dot(alpha, x)  #calc ith element of eigEst
        #all of eigEst has been calculated, end outer for loop    

        valEst = np.amax(x) #get infinity norm of x (just the max element in this case)
        if valEst != 0:
            x = eigEst / valEst  #set the next value of x and scale it
        else:
            x = eigEst
        iterations = iterations + 1

    return (x, valEst)

       
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
    eig = modEig(D, 30)
    print eig[0]
    print "Value = "
    print eig[1]

     
        
