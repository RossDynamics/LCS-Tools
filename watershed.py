import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import h5py as hp
from scipy import ndimage as ndi
from skimage.morphology import watershed, disk
from skimage.feature import peak_local_max
from skimage.filters import rank
f = hp.File('attFTLEOutput.mat','r')
attdata = np.transpose(f['F'][0,:,:])
f.close()
print np.max(rank.gradient(attdata,disk(1)))
print np.min(rank.gradient(attdata,disk(1)))
markers = rank.gradient(attdata,disk(1)) < 1
markers = ndi.label(markers)[0]
gradient = rank.gradient(attdata, disk(1))
labels = watershed(attdata, markers)

# display results
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8), sharex=True, sharey=True, subplot_kw={'adjustable':'box-forced'})
ax = axes.ravel()

ax[0].imshow(attdata, cmap=plt.cm.gray, interpolation='nearest')
ax[0].set_title("Original")

ax[1].imshow(gradient, cmap=plt.cm.spectral, interpolation='nearest')
ax[1].set_title("Local Gradient")

ax[2].imshow(markers, cmap=plt.cm.spectral, interpolation='nearest')
ax[2].set_title("Markers")

ax[3].imshow(labels, cmap=plt.cm.spectral, interpolation='nearest', alpha=.7)
ax[3].set_title("Segmented")

plt.savefig('watershed.tif',tranparent=False,bbox_inches='tight')

