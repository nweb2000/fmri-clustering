import sys

class FmriNode:
	"""
	This class represents an induvidual voxel.
	
	Data Members:
	i - Voxel's x dimension.
	j - Voxel's y dimension.
	k - Voxel's z dimension.
	timeSeries - The time series associated with the voxel.
	"""
	def __init__(self, i, j, k, timeSeries):
		self.i = i
		self.j = j
		self.k = k
		self.timeSeries = timeSeries
		
class FmriArray:
	"""
	This class represents a list of FmriNodes.

	Data Members:
	header - The complete header of the nifti subject image .
	affine - The affine matrix associated with the nifti image.
	extra - Additional extra information associated with the nifti image.
	file_map - Additional file map associated with the nifti image.
	fmriNodeList - A list of FmriNodes. This contains all the intersting voxels. A voxel in a nifti image is 
				  interesting if it contains grey matter and a time series. 
	imgGrey - The entire grey matter file is also stored.
	shape - The dimensions of the nifti image. Note that the dimensions of the nifty subject image should match the 
			dimensions of the grey matter image.
	
	Methods:
	Constructor - Creates an instance of FmriArray from nifti and grey images.
	addFmriNode - Creates and adds a FmriNode into a fmriNodeList.
	sizeFmriArray - Returns the number of voxels (FmriNode) stored.
	"""

	def __init__(self, imgNifti, imgGrey):
		"""
		Summary:
		Constructor - Creates an instance of FmriArray from nifti and grey images.
		
		Input:
		imgNifti - The nifti image of the subject.
		imgGrey - The nifti image that contains probabilities of grey matter for each voxel.
		"""
		self.header = imgNifti.get_header()
		self.affine = imgNifti.get_affine()
		self.extra =  imgNifti.extra
		self.file_map = imgNifti.file_map
		self.fmriNodeList = []
		self.imgGrey = imgGrey
		self.shape = self.header.get_data_shape()
		
	def addFmriNode(self, i, j, k, timeSeries):
		"""
		Summary:
		addFmriNode - Creates and adds a FmriNode into a fmriNodeList.
		
		Input:
		i - Voxel's x dimension.
		j - Voxel's y dimension.
		k - Voxel's z dimension.
		timeSeries - The time series associated with the voxel.
		"""
		self.fmriNodeList.append(FmriNode(i, j, k, timeSeries))
		
	def sizeFmriArray(self):
		"""
		Summary:
		sizeFmriArray - Returns the number of voxels (FmriNode) stored.
		"""
		return len(self.fmriNodeList)
		
