# python 3
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from pytictoc import TicToc

t = TicToc()

# Set meta info

# Generate month grouped list of annual URLs
urlroot='https://rsg.pml.ac.uk/thredds/dodsC/cci/v4.0-release/geographic/monthly/all_products/'
yrange=list(range(1998,2019))
pref='ESACCI-OC-L3S-OC_PRODUCTS-MERGED-1M_MONTHLY_4km_GEO_PML_OCx_QAA-'
mrange=list(range(1,13))
urlmd={}
for yy in yrange:
 for mm in mrange:
  mt='%2.2d'%(mm)
  mty='%2.2d-y'%(mm)
  url='%s%4.4d/%s%4.4d%2.2d-fv4.0.nc'%(urlroot,yy,pref,yy,mm)
  if mt in urlmd:
   urlmd[mt].append(url)
   urlmd[mty].append(yy)
  else:
   urlmd[mt]=[url]
   urlmd[mty]=[yy]
# print(urlmd)

# Get common info from first URL (assume same for all)
dsname=urlmd['01'][0]
dataset = netCDF4.Dataset(dsname)
# List of all variables
vlist=dataset.variables.keys()
print(vlist)
# Check there are some things there we are expecting
if 'water_class1' not in vlist:
 exit()
if 'lat' not in vlist:
 exit()
if 'lon' not in vlist:
 exit()
if 'chlor_a' not in vlist:
 exit()
# Get latitude grid
lats=dataset.variables['lat'][:].data[:]
# Get longitude grid
lons=dataset.variables['lon'][:].data[:]
# Get chunking (three numbers 0 - time, 1 - x, 2- y)
# CCI limits download to 50MB per request, so choose whole chunks that are less than 50MB total.
chk=dataset.variables['chlor_a']._ChunkSizes
chki=chk[1]*16
chkj=chk[2]*8
# Get size (three numbers 0 - time, 1 - x, 2- y)
siz=dataset.variables['chlor_a'].shape
ni=siz[1]
nj=siz[2]
dataset.close()

# Read some data
urlmdcur={}
# urlmdcur['01']=urlmd['01']
# urlmdcur['01-y']=urlmd['01-y']
# urlmdcur['02']=urlmd['02']
urlmdcur=urlmd
mlo=0
mhi=21
vl=['chlor_a']

print(urlmdcur)
# exit()

# Create arrays to average into
favgs={}
for m in urlmdcur:
 favgs[m]={}
 for vv in vl:
  if vv not in vlist:
   exit() 
  favgs[m][vv]={}
  favgs[m][vv]['sum']=np.zeros( (ni,nj) )
  favgs[m][vv]['npt']=np.zeros( (ni,nj) )

# Get data and accumulate averages
mlist=['01','02','03','04','05','06','07','08','09','10','11','12']
for m in mlist:
 for nu,u in enumerate(urlmdcur[m][mlo:mhi]):
  print(u)
  ky='%s-y'%(m)
  yn=urlmdcur[ky][nu]
  ncds=netCDF4.Dataset(u)
  for vv in vl:
   if vv not in vlist:
    exit()
   print("Getting ",vv,u)
   # Temporary array to load into
   fval=np.zeros( (ni,nj) )
   # Download data in using data set chunk size
   nbj=(int)(nj/chkj)
   nbi=(int)(ni/chki)
   nbs=nbi*nbj
   nbcur=0
   for j in range(nbj):
    for i in range(nbi):
     nbcur=nbcur+1
     ilo=i*chki; ihi=ilo+chki-1;
     jlo=j*chkj; jhi=jlo+chkj-1;
     print("Getting chunk",ilo,ihi,jlo,jhi,nbcur/nbs*100.)
     t.tic()
     fval[ilo:ihi,jlo:jhi]=ncds.variables[vv][0,ilo:ihi,jlo:jhi]
     mval=fval[0,0]
     phi=np.where(fval!=mval,fval,0)
     phiC=np.where(fval!=0,1,0)
     favgs[m][vv]['sum']=favgs[m][vv]['sum']+phi
     favgs[m][vv]['npt']=favgs[m][vv]['npt']+phiC
     t.toc()
  ncds.close()
  # Save intermediate per year field to local netcdf
  fname='%s-%s-%4.4d.nc'%(vv,m,yn)
  print(fname)
  ncdsw=netCDF4.Dataset(fname,"w",format="NETCDF4")
  nclat=ncdsw.createDimension("lat",len(lats))
  latitudes=ncdsw.createVariable( "lat","f4",("lat",) )
  latitudes[:]=lats
  nclon=ncdsw.createDimension("lon",len(lons))
  longitudes=ncdsw.createVariable( "lon","f4",("lon",) )
  longitudes[:]=lons
  vals=ncdsw.createVariable( vv,"f4",("lat","lon") )
  vals[:,:]=fval
  ncdsw.close()

print('DONE')

# Calculate averages
for m in favgs:
 for v in favgs[m]:
  phiSum=favgs[m][v]['sum']
  phiNpt=favgs[m][v]['npt']
  phiNpt=np.where(phiNpt!=0,phiNpt,1)
  phiAvg=phiSum/phiNpt
  favgs[m][v]['avg']=phiAvg
  # Save averaged values to local netcdf
  fname='%s-%s-avg.nc'%(v,m)
  print(fname)
  ncdsw=netCDF4.Dataset(fname,"w",format="NETCDF4")
  nclat=ncdsw.createDimension("lat",len(lats))
  latitudes=ncdsw.createVariable( "lat","f4",("lat",) )
  latitudes[:]=lats
  nclon=ncdsw.createDimension("lon",len(lons))
  longitudes=ncdsw.createVariable( "lon","f4",("lon",) )
  longitudes[:]=lons
  vals=ncdsw.createVariable( v,"f4",("lat","lon") )
  vals[:,:]=favgs[m][v]['avg']
  ncdsw.close()

# Make a plot
# mval=fval[0,0]
# phi=np.where(fval!=mval,fval,0)
phi=favgs['01']['chlor_a']['avg']
# phi=np.flip(phi,0)
fig, ax1 = plt.subplots(1, 1, figsize=(48, 18))
vals=np.arange(0,phi.mean()+phi.std()*2,phi.mean()/50.)
plt.imshow(phi[:,:],vmin=0,vmax=phi.mean()+phi.std()*2)
plt.savefig('foo.png')
