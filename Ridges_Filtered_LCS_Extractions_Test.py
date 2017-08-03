#Domain [-300 300]^2 km
#Origin (41.3209371228N, 289.46309961W)
#Projection Lambert
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import h5py as hp
import scipy.ndimage.filters as filter

plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
#origin = 36.7964N,-120.822E
origin = [36.7964, -120.822]

#f=hp.File('attFTLEOutput.mat','r')
f=hp.File('repFTLEOutput.mat','r')
attdata = f['F'][:,:,:]
f.close()
dim = attdata.shape 


dx, dy = np.gradient(attdata[0,:,:],0.3,0.3)
dxdx, dydx = np.gradient(dx,0.3,0.3)
dxdy, dydy = np.gradient(dy,0.3,0.3) 

dirdiv = np.empty([dim[1],dim[2]])
mineig = np.empty([dim[1],dim[2]])

for i in range(dim[1]):
    for j in range(dim[2]):
        eig = np.linalg.eig([[dxdx[i,j],dxdy[i,j]],[dydx[i,j],dydy[i,j]]])
        eigmin =  np.argmin(eig[0])
        dirdiv[i,j] = np.dot(eig[1][:,eigmin],[dx[i,j],dy[i,j]])
        mineig[i,j] = eig[0][eigmin]

#attquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0,:,:])),shading='gouraud')#,cmap='CyanOrange')
#cbar = plt.colorbar()

sd = [0,1,2,3]
for s in sd:
    filtereddirdiv = filter.gaussian_filter(dirdiv,s,s)
    filteredmineig = filter.gaussian_filter(mineig,s,s)
    #attdata = filter.gaussian_filter(attdata,sd,sd)
    tol =[0,-1e-4,-1e-3,-1e-2,-1e-1,-0.25,-0.5,-1]
    index = 0
    for t in tol:
        m = Basemap(width=2000000,height=2000000,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='c',area_thresh=10000.,projection='lcc',\
            lat_1=35.,lat_0=origin[0],lon_0=origin[1])
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        x = np.linspace(0, m.urcrnrx, dim[1])
        y = np.linspace(0, m.urcrnry, dim[2])
        xx, yy = np.meshgrid(x, y)
        attquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0,:,:])),shading='gouraud')
        potridge = np.ma.masked_where(filteredmineig>=t,filtereddirdiv)
        ridge = m.contour(xx, yy, np.transpose(potridge),levels =[0],colors='white')
        ttl = plt.title("FTLE, 7-09-16 0400 UTC, Guass Filtered Ridges, sd = {0:d}, Threshold = {1:e}".format(s, t))

        plt.savefig('owens_GuassFilterRidges_sd{0}_tol_{1}.tif'.format(s,index), transparent=False, bbox_inches='tight')
        #plt.savefig('owens_GuassFilterRidges_{0:d}.tif'.format(int(index)), transparent=False, bbox_inches='tight')
        plt.clf()
        index += 1
