import numpy as np

def modMat(D, i, j, ki, kj, m):
    """Get the (i, j) index of the modularity matrix created by data matrix D.
       @INPUT: D (a data matrix),
    i, j (Index of modularity matrix to access)
                ki, kj (degrees of vertex i, and vertex j)
                m (the sum of degrees for every vertex, which equals 2 * number of graph edges)
       NOTE: The modularity matrix for D will not actually be created, this function will
       be used to represent the modularity matrix implicitly
       (D should be represented by a 2D numpy array)"""
    
    return np.dot(D[:, i], D[:, j]) - ((ki * kj) / m)
                                       
def modEig(D, maxIter, start=None):
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
   else
       x = start 

   degrees = np.zeros(D.shape[1])     
   for i in np.arange(D.shape[1]): #find the degree sums of each node 
       deg_sum = 0
       for j in np.arange(D.shape[1]):
           deg_sum = deg_sum + npdot(D[:, i], D[:, j]) #get the (i, j) index of (D * D-Transpose) and add it to sum

       degrees[i] = deg_sum #save the degree sum for the ith element



   m = np.sum(degrees) #m is 2 * the number of edges (or the sum of all edge weights)        
   iterations = 0
   while iterations < maxIter
       dx = np.dot(D, x) #get D . X
       eigEst = np.zeros(D.shape[1]) 
       for i in np.arange(D.shape[1]):  #calculate eigEst

           alpha = np.zeros(D.shape[1])

           for j in np.arange(D.shape[1]):   #find alpha vector
               alpha[j] = (degrees[i] * degrees[j]) / m

           eigEst[i] = np.dot(D(:, i), dx) - dot(alpha, x)  #ith element of eig index
        #all of eigEst has been calculated, end outer for loop    

        valEst = np.amax(x) #get infinity norm of x (just the max element in this case)
        x = eigEst / valEst 
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
    cluster = split_cluster(D)
    print cluster
        
