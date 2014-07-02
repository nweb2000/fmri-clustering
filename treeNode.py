class treeNode:
    """
    A base node for the cluster tree which holds the data as well as pointers to the
    left and right children of the node.
    """
    def __init__(self, data, left=None, right=None):
        """
        Init the node's data and it's left and right children
        """
        self.left = left
        self.right = right
        self.data = data


    
