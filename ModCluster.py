import numpy as np
np.set_printoptions(threshold=np.nan)

def modularity(D, split, clustList=None, degrees=None):
    """
    Calculate the modularity increase value obtained when the clustering described in split is applied
    on the graph of data matrix/group matrix Dg
    @INPUTS:
    (2d numpy array) D - data matrix of the entire graph
    (numpy array) split - A array with values of either 1 or -1, where the sign indicates which group the ith element in Dg belongs to
    (numpy array) clustList - A list of indices for the data matrix columns which specify what graph
    nodes (columns in the data matrix D) are in the cluster we want to find the modularity of
    (numpy array)  - an array whose ith entry is the degree of the graph node indicated by the ith column in Dg
    @OUTPUTS:
    The unnormalized modularity value
    """
    if degrees == None: #if degrees not supplied calculate them    
        degrees = implicitDegSum(D)

    if clustList != None and clustList.size != D.shape[1]: #if clustList exists and is a subset of the entire graph
        sg = getSubGraph(D, clustList, degrees) #get back a tuple containing group subgraph and group degrees
        Dg = sg[0] 
        Kg = sg[1] 
    else:
        Dg = D
        Kg = degrees

    m = np.sum(degrees)
    
    #Find the product of the modularity matrix of Dg, (denoted as Mg) and split, then find the dot product of the result and split
     
    modifier = np.dot(Dg.T, np.sum(Dg, 1)) - ((Kg * np.sum(Kg)) / m)  #calculate the modifier values to apply to diagnols of Mg
    Dgs = np.dot(Dg, split) 
    Kgs = np.dot(Kg, split) 
    Mgs = np.dot(Dg.T, Dgs) - ((Kg * Kgs) / m) - (modifier * split) #calculate product of Dg and split
    return np.dot(split, Mgs)  #calculate and return the modularity, divide this value by 2*m to normalize

def implicitDegSum(D):
    """Finds the sums of the degrees of each node in the graph of the implicitly represented adjacency matrix of data matrix D"""
    degrees = np.zeros(D.shape[1])
    rowsums = np.zeros(D.shape[0]) 
    for i in np.arange(D.shape[0]): #find the sum of each row of D and save it 
        rowsums[i] = np.sum(D[i, :])
        
    for i in np.arange(D.shape[1]):   #find the sum of each row in the implicit adjacency matrix
        degrees[i] = np.dot(D[:, i], rowsums)

    return degrees

def getSubGraph(D, clustList, degrees=None):
    """
    Get the subgraph of data matrix D containing only the nodes indicated in clustList
    @INPUTS
    D - The entire data matrix
    clustList - a list of graph node ids whose columns should be included in Dg
    degrees - a list of the degrees of the entire graph, if degrees is not passed in only the subgraph is returned,
    if it is passed in, a tuple of both the group subgraph and its respective degrees are returned
    @OUTPUTS
    A subgraph of D which contains the columns of D indicated in clustList, or a tuple containing the subgraph of D
    and and array of the respective degrees of that subgraph
    """
    if clustList == None:
        Dg = D
    else:
        Dg = np.zeros((D.shape[0], clustList.size))  #create a new data matrix made up of the columns in the group
        for i in np.arange(clustList.size):
            Dg[:, i] = D[:, clustList[i]]

    if degrees != None:
        if clustList != None:
            Kg = np.array([degrees[i] for i in clustList])
        else:
            Kg = degrees
            
        out = (Dg, Kg)
    else:
        out = Dg
        
    return out


def powerMethod(Dg, Kg, start, m, modifier, precision, shift=0):
    """
    Does the powerMethod on Dg to find the dominent eigen vector
    @INPUTS:
    Dg - The data matrix, or a subgraph of the data matrix
    Kg - A array whose ith entry contains the degree of the graph node represented by the ith
    column in Dg
    start - the initial vector to start the power method with
    m - a constant value equal to 2 * (the sum of the degrees of the nodes in the entire graph)
    precision - how good of an estimate is needed, this value should be between (0, 1)
    shift - the shift value to use, shifting the matrix helps find the most positive eigen vector
    @OUTPUTS:
    A python tuple, the first element contains the dominent eigenvector estimate, while the second element
    contains the estimated dominent eigenvalue
    """
    converged = False #whether or not we have converged to desired precision
    x = start
    while not converged:  #do the power method until we converge
        Dgx = np.dot(Dg, x) #get Dg . X
        Kgx = np.dot(Kg, x) #get Kg . X
        eigEst = np.dot(Dg.T, Dgx) - ((Kg * Kgx) / m) - (modifier * x) + (shift * x) #calculate eigenvector estimate 
        
        valEst = np.dot(x, eigEst) / np.dot(x, x) #estimate eigen value using rayleigh quotient

        if np.linalg.norm(eigEst - (valEst * x)) <= precision: #if our desired precision has been reached
            converged = True #we hast converged unto the eigen vector, hazaaaah!

        n = np.linalg.norm(eigEst)
        if n != 0:
            x = eigEst / np.linalg.norm(eigEst)   #set the next value of x and scale it
        else:
            x = eigEst
       
    return (eigEst, valEst)

    
def modEig(D, clustList=None, degrees=None, precision=0.95, start=None):
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
    (numpy array) degrees - a array of degree sums for the ENTIRE graph (row sums of the implicit
     adjaceny matrix)
    (float) precision - how much precision is needed with the power method, should be a value between (0, 1)
    (numpy array) start - a vector to start the power method with, optional
    @OUTPUT:
    A python tuple, the first element contains the leadin eigenvector estimate, while the second element
    contains the estimated leading eigenvalue
    """
    #Set initial (the first guess for the eigenvector in the power method
    if start == None:
        if clustList != None:
            initial = np.random.random_integers(1, 100, clustList.size)#initial = np.zeros(clustList.size)
            initial = initial / np.linalg.norm(initial)
        else:
            initial = np.random.random_integers(1, 100, D.shape[1])#np.zeros(D.shape[1])
            initial = initial / np.linalg.norm(initial)
    else:
        initial = start

    if degrees == None: #if degrees not supplied calculate them    
        degrees = implicitDegSum(D)
        
    m = np.sum(degrees) #m is 2 * the number of edges (or the sum of all edge weights)        
    iterations = 0 #this is just for testing
    
    if clustList != None and clustList.size != D.shape[1]: #if clustList exists and is a subset of the entire graph
        sg = getSubGraph(D, clustList, degrees) #get back a tuple containing group subgraph and group degrees
        Dg = sg[0] 
        Kg = sg[1] #make a new array containing the degrees sums of the nodes in cluster
        modifier = np.dot(Dg.T, np.sum(Dg, 1)) - ((Kg * np.sum(Kg)) / m)  #calculate the modifier values to apply to diagnols of group mod mat
    else: #we are splitting the entire graph
        Dg = D
        Kg = degrees
        modifier = 0

    vect, val = powerMethod(Dg, Kg, initial, m, modifier, precision) #use power method to get dominent eigen vector
    if val < 0: #if the dominint eigen vector is negative, do the power method again, shifting the matrix by val/2
        vect, val = powerMethod(Dg, Kg, initial, m, modifier, precision, np.abs(val)/2) #this will get us the positive leading eigenvector
            
    return (vect, val)
       
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
    Two lists of graph node indices which show how clustList was split,
    these are stored in a python tuple
    """
    if clustList != None and clustList.size == 0:
        print "Cannot split an empty partition!"
    
    if degrees == None:
        degrees = implicitDegSum(D)

    if clustList == None:
        clustList = np.arange(D.shape[1])

    p_eig = modEig(D, clustList, degrees)[0]
      # split = p_eig / np.abs(p_eig)
    
    clust_1 = np.array([clustList[x] for x in np.arange(p_eig.size) if p_eig[x] > 0])
    clust_2 = np.array([clustList[x] for x in np.arange(p_eig.size) if p_eig[x] < 0])
    clusters = (clust_1, clust_2)

    return clusters
    
if __name__ == "__main__":
    D = np.genfromtxt("test.txt", delimiter=',')
    c = splitCluster(D)
    clust1 = c[0]
    clust2 = c[1]

    cl = splitCluster(D, clust1)
    cl2= splitCluster(D, clust2)
    a = cl[0]
    b = cl[1]
    c = cl2[0]
    d = cl2[1]
 
    e = splitCluster(D, a)
    f = splitCluster(D, b)
    g = splitCluster(D, c)
    h = splitCluster(D, d)

    print "CLUST A\n\n", e
    print "\n\nCLUST B\n\n", f
    print "\n\nCLUST C\n\n", g
    print "\n\nCLUST D\n\n", h

