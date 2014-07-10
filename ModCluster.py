import numpy as np

def modMat(D, i, j, ki, kj, m):
    """THIS FUNCTION IS JUST TO SHOW HOW MODULARITY INDEX IS CALCULATED, AND IS NOT ACTUALL NEEDED
       Get the (i, j) index of the modularity matrix created by data matrix D.
       @INPUT: D (data matrix where the columns of D represent different nodes and the 
    rows represent a piece of time series data for the corresponding node)D (a data matrix),
               i, j (Index of modularity matrix to access)
               ki, kj (degrees of vertex i, and vertex j)
               m (the sum of degrees for every vertex, which equals 2 * number of graph edges)"""
    
    return np.dot(D[:, i], D[:, j]) - ((ki * kj) / m)

def implicitDegSum(D):
    """Finds the sums of the degrees of each node in the graph of the implicitly represented adjacency matrix of data matrix D"""
    degrees = np.zeros(D.shape[1])
    rowsums = np.zeros(D.shape[0]) 
    for i in np.arange(D.shape[0]): #find the sum of each row of D and save it 
        rowsums[i] = np.sum(D[i, :])
        
    for i in np.arange(D.shape[1]):   #find the sum of each row in the implicit adjacency matrix
        degrees[i] = np.dot(D[:, i], rowsums)

    return degrees

def modEig(D, clustList=None, degrees=None, maxIter=30, start=None):
    """
    Find the principle eigenvector for the implicitly represented modularity matrix of
    data matrix D using the power method
    NOTE: For now pass in the entire data matrix to D, this function will take care of pulling out
    the correct columns, do the same for degrees (just pass the entire thing in
    NOTE: If clustList is not passed in, then it is implied that the entire graph is being split 
    @INPUT:
    (2D numpy array) D - data matrix where the columns of D represent different nodes and the 
    rows represent a piece of time series data for the corresponding node
    (numpy array) clustList - A list of indices for the data matrix columns which specify what graph
    nodes (columns in the data matrix D) are in the cluster that is to be split
    (numpy array) degrees - a array of degree sums for each graph node(row sums of the implicit
     adjaceny matrix
    (int) maxIter - how much iterations of power method to go through, this is currently just for testing
     we will use some kind of tolerence level later on...
    (numpy array) start - a vector to start the power method with, optional
    
    
    TODO: Replace the iterations parameter with some sort of normalized tolerence level 
    """
    if clustList == None:
        clustList = np.arange(D.shape[1])  #the cluster is just the entire graph
       
    #set x to the first estimate of the principle eigen vector
    if start == None:
        x = np.zeros(clustList.size)
        x[0] = 1 #initialize x_0 as a unit vector, with all elements = to zero except for the first
    else:
        x = start 

    if degrees == None:    
        degrees = implicitDegSum(D)
        
    m = np.sum(degrees) #m is 2 * the number of edges (or the sum of all edge weights)        
    iterations = 0
    
    if clustList.size != D.shape[1]: 
        Dg = np.array([D[:, i] for i in clustList]) #pull out the columns of the group
        Dg = Dg.T
        Kg = np.array([degrees[i] for i in clustList]) #get the degree sum vector for the group
        
    else:
        Dg = D
        Kg = degrees

    Dg_sum = np.sum(Dg, 1)
    Kg_sum = np.sum(Kg)
                      
    while iterations < maxIter:  #do the power method for maxIter iterations
        Dgx = np.dot(Dg, x) #get Dg . X
        eigEst = np.zeros(Dg.shape[1])
        
        for i in np.arange(Dg.shape[1]):  #calculate eigEst
            term1 = np.dot(Dg[:, i], Dgx - (Dg_sum * x[i]))
            term2 = (Kg[i] / m) * ((Kg_sum * x[i]) - np.dot(Kg, x))
            eigEst[i] = term1 + term2  #calc ith element of eigEst
        #all of eigEst has been calculated, end outer for loop    

        valEst = np.amax(x) #get infinity norm of x (just the max element in this case)
        if valEst != 0:
            x = eigEst / valEst  #set the next value of x and scale it
        else:
            x = eigEst
        iterations = iterations + 1

    return (x, valEst)

       
def splitCluster(D, clustList=None, degrees=None):
   """
   Using modularity measures split the cluster listed in clustList into two different clusters
   @INPUTS:
      (2D numpy array) D - data matrix whose columns represent different nodes of the graph while the 
      rows represent a piece of time series data for the corresponding node
      (numpy array) clustList - A list of indices for the data matrix columns which specify what graph
      nodes are in the cluster that is to be split
      (numpy array) degrees - a array of degree sums for each graph node(row sums of the implicit
      adjaceny matrix
   @OUPUT:
      Two lists of graph node indeces which show how clustList was split,
      these are stored in a python tuple
   """
   if clustList==None:
       clustList = np.arange(D.shape[1])
       
   p_eig = modEig(D, clustList, degrees)[0]
   clust_1 = np.array([clustList[x] for x in np.arange(p_eig.size) if p_eig[x] > 0])
   clust_2 = np.array([clustList[x] for x in np.arange(p_eig.size) if p_eig[x] < 0])
   clusters = (clust_1, clust_2)
   return clusters
    
if __name__ == "__main__":
    D = np.genfromtxt("test.txt", delimiter=',')
    
    clusters = splitCluster(D)

  
    c1 = splitCluster(D, clusters[0])
    c2 = splitCluster(D, clusters[1])

    b1 = splitCluster(D, c2[0])
    b2 = splitCluster(D, c2[1])
    
    print clusters[0], " " , clusters[1]
    print "\nc1\n"
    print c1[0], " " ,c1[1]
    print "\nc2\n"
    print c2[0], " " ,c2[1]
    

    print "\nb1\n"
    print b1[0], " " ,b1[1]
    print "\nc2\n"
    print b2[0], " " ,b2[1]
     
