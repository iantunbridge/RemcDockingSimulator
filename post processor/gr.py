#!/usr/bin/env python

import sys,os,math

aalist = { 'ALA': 'A', 'GLY': 'G', 'THR': 'T', 'TYR': 'Y',
	   'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'TRP': 'W',
	   'GLU': 'E', 'ASP': 'D', 'SER': 'S', 'ASN': 'N',
	   'GLN': 'Q', 'PRO': 'P', 'PHE': 'F', 'ARG': 'R',
	   'CYS': 'C', 'HIS': 'H', 'LYS': 'K', 'MET': 'M',
	   'CGU': 'E', 'DBY': 'Y', 'GLZ': 'G', 'GLQ': 'E',
	   'HSD': 'H', 'HEM': 'X', 'ABA': 'B', 'CSO': 'C',
	   'ASPP': 'D'}

def setup(pdb):
	blob = {}
	atoms = filter(lambda x: x[0:6]=="ATOM  ", open(pdb).readlines())
	for a in atoms:
		restype = a[17:20]
		oneletter = aalist[restype]
		chain = a[21]
		if chain not in blob.keys():
			blob[chain] = { "seq": "", "xyz": [], "idx": [] }
		blob[chain]["seq"] += oneletter
		X = float(a[29:37])
		Y = float(a[37:45])
		Z = float(a[45:53])
		resnum = int(a[22:26])
		blob[chain]["xyz"].append((X,Y,Z))
		blob[chain]["idx"].append(resnum)
	chains = blob.keys()
	chains.sort()
	types = {}
	com = {}
	for c in chains:
		found = 0
		for k in types.keys():
			if blob[c]["seq"] == blob[k]["seq"]:
				types[k].append(c)
				found =1 
				break
		if not found:
			types[c] = [c]
			#
			c_atoms = blob[c]["xyz"]
			resnums = blob[c]["idx"]
			nc = len(c_atoms)
			Xav,Yav,Zav=0.,0.,0.
			# here we assume all atoms in one molecule are in 
			# the same periodic image
			for i in range(nc):
				Xav+=c_atoms[i][0]
				Yav+=c_atoms[i][1]
				Zav+=c_atoms[i][2]
			Xav/=float(nc)
			Yav/=float(nc)
			Zav/=float(nc)
			mind = 1000.
			for i in range(nc):
				Xdiff=c_atoms[i][0]-Xav
				Ydiff=c_atoms[i][1]-Yav
				Zdiff=c_atoms[i][2]-Zav
				dist = (Xdiff**2+Ydiff**2+Zdiff**2)**0.5
				if dist<mind:
					mind = dist
					idx = resnums[i]
			com[c] = idx
	return types,com

def parse_com(f,types,com):
	atoms = filter(lambda x: x[0:6]=="ATOM  ", open(f).readlines())
	pchain = "@"
	coms = {}
	for a in atoms:
		#restype = a[17:20]
		#oneletter = aalist[restype]
		chain = a[21]
		if chain != pchain:
			for k in types.keys():
				if chain in types[k]:
					comres = com[k]
			pchain = chain
		resnum = int(a[22:26])
		if resnum  == comres:
			X = float(a[29:37])
			Y = float(a[37:45])
			Z = float(a[45:53])
			coms[chain] = (X,Y,Z)
	return coms


	

boxlen = float(sys.argv[1])
volume = boxlen**3

pdb_files = sys.argv[2:]

types,com = setup(pdb_files[0])

lo,nbin = 0.,75
hi = boxlen/2.0 - 1.0
dr = (hi-lo)/float(nbin)

typekeys = types.keys()
ntype = len(typekeys)

gr = {}
for i in range(ntype):
	gr[typekeys[i]] = {}
	for j in range(i,ntype):
		gr[typekeys[i]][typekeys[j]] = []
		for b in range(nbin):
			gr[typekeys[i]][typekeys[j]].append(0.0)

nframe = 0 
for f in pdb_files:
	nframe += 1
	com_crd = parse_com(f,types,com)
	#print com_crd
	for i in range(ntype):
		typeI = typekeys[i]
		for j in range(i,ntype):
			typeJ = typekeys[j]
			for Ichain in types[typeI]:
				icom = com_crd[Ichain]
				for Jchain in types[typeJ]:
					if Jchain == Ichain:
						continue
					jcom = com_crd[Jchain]
					dX = icom[0]-jcom[0]
					dY = icom[1]-jcom[1]
					dZ = icom[2]-jcom[2]
					dX -= boxlen*round(dX/boxlen)
					dY -= boxlen*round(dY/boxlen)
					dZ -= boxlen*round(dZ/boxlen)
					dR = (dX*dX+dY*dY+dZ*dZ)**0.5
					bin = int(math.floor((dR-lo)/dr))
					if bin<0 or bin>=nbin:
						continue
					gr[typeI][typeJ][bin] += 1.0


for i in range(ntype):
	typeI = typekeys[i]
	nI = len(types[typeI])
	for j in range(i,ntype):
		typeJ = typekeys[j]
		nJ = len(types[typeJ])
		outp = open("gr_%s_%s.dat"%(typeI,typeJ),"w")
		for b in range(nbin):
			r = lo+(0.5+float(b))*dr
			normfac = float(nframe)*4.*math.pi*r*r*dr*nI*nJ/volume
			gr_norm = gr[typeI][typeJ][b]/normfac
			outp.write("%8.3f %8.3f\n"%(r,gr_norm))
		outp.close()




