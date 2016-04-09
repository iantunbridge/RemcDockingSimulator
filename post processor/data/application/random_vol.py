#!/usr/bin/env python

import math,random,sys

def crd(line):
	x = float(line[30:38])
	y = float(line[38:46])
	z = float(line[46:54])
	return x,y,z

def minmax(crds):
	print crds[0]
	xmin,ymin,zmin = crds[0]
	xmax,ymax,zmax = crds[0]
	for r in crds:
		x,y,z=r
		if x<xmin:
			xmin = x
		if x>xmax:
			xmax = x
		if y<ymin:
			ymin = y
		if y>ymax:
			ymax = y
		if z<zmin:
			zmin = z
		if z>zmax:
			zmax = z
	return xmin,xmax,ymin,ymax,zmin,zmax



for i in range(10):
	print random.random()

pdblines = filter(lambda x: x[0:6]=="ATOM  ", open(sys.argv[1]).readlines())

rcut = 6.0
xyz = map(crd,pdblines)
natom = len(xyz)
print natom
xmin,xmax,ymin,ymax,zmin,zmax = minmax(xyz)
xmin -= rcut
ymin -= rcut
zmin -= rcut
xmax += rcut
ymax += rcut
zmax += rcut
xbox = xmax-xmin
ybox = ymax-ymin
zbox = zmax-zmin

V = xbox*ybox*zbox


nmc =  10000
nblock = 10

sum = 0.0
sum2 = 0.0
for b in range(nblock):
	goodpoints = 0
	for i in range(nmc):
		#print random.random()
		
		x=xmin+xbox*random.random()
		y=ymin+ybox*random.random()
		z=zmin+zbox*random.random()
		good = 0
		for t in range(natom):
			xt,yt,zt = xyz[t]
			dx = xt-x
			dy = yt-y
			dz = zt-z
			dr = ( dx*dx + dy*dy + dz*dz ) ** 0.5
	#		print dr
			if dr < rcut:
				good = 1
				break
		goodpoints += good

#	print goodpoints,nmc,V
	Vprot = float(goodpoints)/float(nmc)*V
	sum += Vprot
	sum2 += Vprot**2

meanV = sum / nblock
sdevV = (sum2/nblock-meanV**2)**0.5
sdevMeanV = sdevV/(float(nblock))**0.5

print meanV, sdevMeanV
#print Vprot
#print xbox*ybox*zbox
