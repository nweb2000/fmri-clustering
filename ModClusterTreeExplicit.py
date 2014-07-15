import numpy as np
from heapq import *
from ModCluster import *
from ModClusterExplicit import *

class ModClusterNode:
	"""
	This represents a Cluster Node that is split based on modularity. 
	
	Data Members:
	cluster - An array of nodes in the cluster. Nodes are represented by their node ids.
	modClusterId - The cluster id for the current cluster.
	size - The number of nodes in the cluster.
	left/ right - The left and right children of the current cluster. The cluster is usually split
				  into two based on the modularity.
	parent - Pointer to the parent cluster.
	
	Member Functions:
	constructor - Creates an instance of a modClusterNode.
	< operator overloaded - This is used in the heapq data structure to compare two instances of 
							modClusterNode. Here a cluster with higher size is given higher priority.
	string operator - Returns a string that represents a modClusterNode.
	"""
	def __init__(self, clus, clId = None, lChild = None, rChild = None, par = None):
		"""
		Summary:
		Constructor - requires clus, an array of node ids in the cluster.
		
		Input:
		clId - A designated cluster Id, optional
		lChild/ rChild/ par - A designated left/ right child and parent, optional
		"""
		self.cluster = clus
		self.modClusterId = clId
		self.size = clus.shape[0]
		self.left = lChild
		self.right = rChild
		self.parent = par
		
	def __lt__(self, other):
		return self.size >= other.size
		
	def __str__(self):
		return "ClusterId " + str(self.modClusterId) + \
				 "\nCluster Size " + str(self.size) + \
				 "\nCluster Nodes " + str(self.cluster)

class ModClusterTree:
	"""
	This class represents a cluster tree. Each cluster is split based on the modularity into 
	two new clusters which become ites children.
	
	Data Members:
	root - A modClusterNode that points to the root cluster in the tree.
	nodeIdLabeler - This is used to label the nodes as they are split. This is a unique id 
					to each cluster.
	clusterMap - A dictionary that maps cluster ids to actual modClusterNodes. This is used in 
				extracting all clusters from a certain depth.
				
	Functions:
	Constructor - Creates an instance od modClusterTree.
	buildModClusterTree - builds the entire tree from a data matrix D.
	extractClusters - Returns all clusters as a list of modClusterNodes that are greater than a 
					  given depth. 
	disp - Displays the entire tree. This is a wrapper function that calls displayTree.
	displayTree - Recursively prints the entire contents of the tree.
	"""
	def __init__(self):
		self.root = None
		self.nodeIdLabeler = 1
		self.clusterMap = {}
		
	def buildModClusterTree (self, D):
		"""
		Summary:
		This function builds the entire cluster tree defined by modularity. 
		
		Input:
		D - A data matrix in which each voxel's data is represented along the columns.
		"""
		(lenData, numNodes) = D.shape					# Get the length of each data and the total number of nodes.
		self.root = ModClusterNode(np.arange(numNodes))	# Create a root node that represents the entire graph.
		
		clusterQueue = []								# A priority queue, used to select the highest priority cluster 
														# to split.
		heappush(clusterQueue, self.root)
		
		A = D.transpose().dot(D)						# Get the corresponding adjacency matrix.
		modMat = makeInitModMap(A)						# Get the modularity matrix.
		degrees = degSum(A)								# Get the degree vector.
		
		while (clusterQueue):							# Do this as long as there is something in the queue.
			clusterToSplit = heappop(clusterQueue)		# Remove the highest priority cluster.
			clusterToSplit.modClusterId = self.nodeIdLabeler		# Note that the clusterId is assigned during the split.
			self.clusterMap[self.nodeIdLabeler] = clusterToSplit	# Update the clusterMap.
			self.nodeIdLabeler += 1
			
			clusterList = clusterToSplit.cluster		# Get the actual array of clusters
			(clusterAlist, clusterBlist) = splitCluster(modMat, clusterList)	# Split the cluster into two based on modularity.
			
			# Insert valid clusters into the queue and let them be children of ths split cluster.
			if not np.array_equal(clusterAlist, clusterList) and len(clusterAlist) != 0:
				clusterA = ModClusterNode(clusterAlist, par = clusterToSplit)
				clusterToSplit.left = clusterA
				heappush(clusterQueue, clusterA)
			if not np.array_equal(clusterBlist, clusterList) and len(clusterBlist) != 0:
				clusterB = ModClusterNode(clusterBlist, par = clusterToSplit)
				clusterToSplit.right = clusterB
				heappush(clusterQueue, clusterB)
				
	def extractClusters(self, depth):
		"""
		Summary:
		Extracts all clusters that are greater than a given depth.
		
		Input:
		depth - The depth (or resolution) of the clusters need to be returned. The way the clusters are split,
				depth also corresponds to modClusterId.
				
		Output:
		resultClusters - A list of modClusterNodes that are at depths bigger than depth. Note that only the largest 
						cluster greater than depth is returned. Smaller clusters that are part of a larger cluster
						is represented within the larger cluster.
		"""
		resultClusterList = []		# A list of modClusterIds.
		resultClusters = []			# A lsit of modClusterNodes.
		
		for i in range(depth, self.nodeIdLabeler):		# Iterate from depth till the size od clusters.
			cluster = self.clusterMap[i]				# Get the cluster (modClusterNode) that represents a certain depth.
			
			if np.any(np.intersect1d(np.array(resultClusterList), cluster.cluster)):	# If it is already included, ignore.
				continue
			resultClusterList += cluster.cluster.tolist()	# First occurrence of a cluster, include into the lists.
			resultClusters.append(cluster)
			
		return resultClusters
		
	def displayTree(self, node, char):
		if node == None:
			return
		print char,
		print(node)
		self.displayTree(node.left, char + ' ')
		self.displayTree(node.right, char + ' ')
	
	def disp(self):
		self.displayTree(self.root, ' ')

if __name__ == "__main__":
	modTree = ModClusterTree()
	D = np.genfromtxt("test.txt", delimiter=',')
	modTree.buildModClusterTree(D)
	modTree.disp()
	
	listClusters = modTree.extractClusters(4)
	print ("\nCluster List")	
	for clus in listClusters:
		print(clus)
