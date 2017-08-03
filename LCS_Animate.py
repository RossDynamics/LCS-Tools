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

f=hp.File('attFTLEOutput.mat','r')
attdata = f['F'][:,:,:]
f.close()
dim = attdata.shape 
f=hp.File('repFTLEOutput.mat','r')
repdata = f['F'][:,:,:]
f.close()


atol =-0.25
rtol =-0.175
#rtol =-0.15
for t in range(dim[0]):
    print t
    adx, ady = np.gradient(attdata[t,:,:],0.3,0.3)
    adxdx, adydx = np.gradient(adx,0.3,0.3)
    adxdy, adydy = np.gradient(ady,0.3,0.3) 
    
    rdx, rdy = np.gradient(repdata[t,:,:],0.3,0.3)
    rdxdx, rdydx = np.gradient(adx,0.3,0.3)
    rdxdy, rdydy = np.gradient(ady,0.3,0.3) 

    adirdiv = np.empty([dim[1],dim[2]])
    amineig = np.empty([dim[1],dim[2]])

    rdirdiv = np.empty([dim[1],dim[2]])
    rmineig = np.empty([dim[1],dim[2]])

    for i in range(dim[1]):
        for j in range(dim[2]):
            aeig = np.linalg.eig([[adxdx[i,j],adxdy[i,j]],[adydx[i,j],adydy[i,j]]])
            aeigmin =  np.argmin(aeig[0])
            adirdiv[i,j] = np.dot(aeig[1][:,aeigmin],[adx[i,j],ady[i,j]])
            amineig[i,j] = aeig[0][aeigmin]
            
            reig = np.linalg.eig([[rdxdx[i,j],rdxdy[i,j]],[rdydx[i,j],rdydy[i,j]]])
            reigmin =  np.argmin(reig[0])
            rdirdiv[i,j] = np.dot(reig[1][:,reigmin],[rdx[i,j],rdy[i,j]])
            rmineig[i,j] = reig[0][reigmin]
    m = Basemap(width=2000000,height=2000000,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='h',area_thresh=100.,projection='lcc',\
        lat_1=35.,lat_0=origin[0],lon_0=origin[1])
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(np.arange(20,50,2.5),labels=[True,False,False,False])
    m.drawmeridians(np.arange(-135,-105,2.5),labels=[False,False,False,True])
    m.drawstates()
    x = np.linspace(0, m.urcrnrx, dim[1])
    y = np.linspace(0, m.urcrnry, dim[2])
    xx, yy = np.meshgrid(x, y)
    #plt.hold(True)
    apotridge = np.ma.masked_where(amineig>=atol,adirdiv)
    aridge = m.contour(xx, yy, np.transpose(apotridge),levels =[0],colors='blue')
    rpotridge = np.ma.masked_where(rmineig>=rtol,rdirdiv)
    rridge = m.contour(xx, yy, np.transpose(rpotridge),levels =[0],colors='red')
    minute = 15*t
    h, minute = divmod(minute,60)
    x, y = m(-119.3,36.2167)
    print x,y
    m.scatter(x,y,marker='*',color='g',s=20*16)
    plt.annotate('Visalia',xy=(x-0.05*x,y+0.03*y),size=15)
    plt.title("FTLE, 6-{0}-2016 {1:02d}{2:02d} UTC".format(9+(4+h)//24, (4+h)%24, minute))
    plt.savefig('owens_lcs_{0:04d}.tif'.format(t), transparent=False, bbox_inches='tight')
    plt.clf()
