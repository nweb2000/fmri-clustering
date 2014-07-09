import nibabel as nib
from FmriArray import *

def filterGreyMatter (niftifile, greyfile, prob = 0):
	"""
	Summary:
	This function takes in 2 nifti files, 1 that represents subject data and another that represents grey matter
	and extracts all the interesting voxels and returns a data structure that contains an array of voxels.
	 
	Input:
	niftifile - The name of a nifti file that represents subject data in either .nii or .nii.gz formats
	greyfile - The name of a nifti file that represents grey matter in either .nii or .nii.gz formats. This file contains 
	a value for each voxel that gives the probability that the voxel contains grey matter.
	prob - The probability threshold that represents interesting voxels. All voxels are interesting with associate grey 
	matter probability > prob.
	
	Output:
	fmriArr - A FmriArray data structure. This is a list of interesting voxels.
	 
	The x-y-z dimensions of the nifti image must match the grey image.
	"""
	
	# Load the Nifti files.
	imgNifti = nib.nifti1.load(niftifile)
	imgGrey = nib.nifti1.load(greyfile)
	
	# Check if the dimensions match.
	imgNiftiShape = imgNifti.get_header().get_data_shape()
	imgGreyShape = imgGrey.get_header().get_data_shape()
	if imgNiftiShape[0] != imgGreyShape[0] or \
	   imgNiftiShape[1] != imgGreyShape[1] or \
	   imgNiftiShape[2] != imgGreyShape[2]:
		print 'The x-y-z dimensions of the nifti and grey images must be the same'
		print 'The dimensions of the nifti image', self.shape
		print 'The dimensions of the grey image', imgGreyShape
		sys.exit()
	
	# Create a FmriArray object.
	fmriArr = FmriArray(imgNifti, imgGrey)
	
	# Get the nifti image data.
	dataNifti = imgNifti.get_data()
	dataGrey = imgGrey.get_data()

	# Extract the interesting voxels.
	for i in range(imgNiftiShape[0]):
		for j in range(imgNiftiShape[1]):
			for k in range(imgNiftiShape[2]):
				if dataGrey[i][j][k] > prob and not (dataNifti[i][j][k].min() == 0 and dataNifti[i][j][k].max() == 0):
					fmriArr.addFmriNode(i, j, k, dataNifti[i][j][k])

	return fmriArr;


