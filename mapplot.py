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

cdict = {'red':  [(0.0, 0.0000, 0.0000),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 1.0000, 1.0000)],
        'green': [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.5450, 0.5450)],
        'blue':  [(0.0, 0.5450, 0.5450),
                  (0.5, 1.0000, 1.0000),
                  (1.0, 0.0000, 0.0000)]}
plt.register_cmap(name='CyanOrange', data=cdict)
plt.figure(num=None, figsize=(16, 16), dpi=80, facecolor='w', edgecolor='k')
#origin = 36.7964N,-120.822E
origin = [36.7964, -120.822]
# setup lambert conformal basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel (defaults to lat_1).
# lon_0,lat_0 is central point.
# rsphere=(6378137.00,6356752.3142) specifies WGS84 ellipsoid
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(width=2000000,height=2000000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='c',area_thresh=1000.,projection='lcc',\
            lat_1=35.,lat_0=origin[0],lon_0=origin[1])
m.drawcoastlines()
m.drawparallels(np.arange(20,50,2.5),labels=[True,False,False,False])
m.drawmeridians(np.arange(-135,-105,2.5),labels=[False,False,False,True])
m.drawcountries()
m.drawstates()
f=hp.File('attFTLEOutput.mat','r')
attdata = f['F'][:,:,:]
f.close()
dim = attdata.shape 
'''
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

print dx.shape
#print dim
attmean = np.mean(attdata)
#print attmean
"""
f=hp.File('RepFTLEOutput.mat','r')
repdata = f['F'][:,:,:]
f.close
#print repdata.shape
repmean = np.mean(repdata)
#print repmean
"""
attdata = np.ma.masked_less(attdata,3*attmean)
#repdata = -1*np.ma.masked_less(repdata,3*repmean)
attmax = np.max(attdata)
#repmin = np.min(repdata)

#for f in range(3):
dirdiv = filter.gaussian_filter(dirdiv,3,3)
mineig = filter.gaussian_filter(mineig,3,3)
'''

x = np.linspace(0, m.urcrnrx, dim[1])
y = np.linspace(0, m.urcrnry, dim[2])

xx, yy = np.meshgrid(x, y)
#tol = -1e-4
attquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0,:,:])),shading='gouraud')
#attquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(attdata[0,:,:])),shading='gouraud',vmin=0,vmax=attmax)#,cmap='CyanOrange')
#potridge = np.ma.masked_where(mineig>=tol,dirdiv)
#attquad = m.pcolormesh(xx, yy, np.transpose(potridge),shading='gouraud')#,cmap='CyanOrange')
#potridge = np.ma.masked_where(mineig>=tol,dirdiv)
#attquad = m.contour(xx, yy, np.transpose(potridge),levels =[0],colors='blue')#,cmap='CyanOrange')
#v = attquad.collections[0].get_paths()[0].vertices

#repquad = m.pcolormesh(xx, yy, np.transpose(np.squeeze(repdata[0,:,:])),shading='gouraud',vmin=repmin,vmax=attmax,cmap='CyanOrange')
#cbar = plt.colorbar()
#cbar.set_label('Attractiveness, s^-1')
#ttl = plt.title("FTLE, 7-21-17 1000 UTC/6am EDT")
#plt.savefig('owens.tif', transparent=False, bbox_inches='tight')

#def init():
#attquad.set_array([])
#repquad.set_array([])
#ttl.set_text("")
#return attquad, ttl
#def animate(t, *args)
for t in range(dim[0]):
    attquad.set_array(np.ravel(np.transpose(np.squeeze(attdata[t,:,:]))))
    #repquad.set_array(np.ravel(np.transpose(np.squeeze(args[3][t,:,:]))))
    ampm = ['am', 'pm']
    #h = int(t)
    #minute = int(t%1*15)
    m = 15*t
    h, minute = divmod(m,60)
    plt.title("FTLE, 6-{0}-2016 {1:02d}{2:02d} UTC".format(9+(4+h)//24, (4+h)%24, minute))
    plt.savefig('owens-frame-{:04d}.tif'.format(int(t)), transparent=False, bbox_inches='tight')
    #return attquad, ttl
    

#myargs = (xx,yy,attdata,repdata)
#myargs = (xx,yy,attdata)
#anim = animation.FuncAnimation(plt.gcf(),animate,fargs=myargs,frames=dim[0],repeat=False)
#plt.savefig('day3-frame-{:04d}.tif'.format(int(t*4+1)), transparent=False, bbox_inches='tight')
#plt.close()
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=4)
#anim.save('FTLE.mp4')#, writer=writer)
#plt.show()

