import itertools as it
import numpy as np
import sys
import os.path

dim = 3
L = 4
liz= 3 #must be odd, only restriction

### tuples come with (z,y,x) or (y,x) etc.
def dFlatten(dcoord):
	dindex = 0
	count =  0
	for i in dcoord:
		dindex =  dindex + i * L**(count)
		count = count + 1
	return dindex

latticeShape = tuple([L for i in range(dim)])
latticeValues = range(L)
translationRange = np.array(range(liz)) - (liz-1)/2

#build translations of liz center
lizShape = tuple([liz for i in range(dim)])
tdimlist = [translationRange for i in range(dim)]
translations = it.product(*tdimlist)

#Make the lattice
latticeCoords = np.ndindex(*[ L for i in range(dim)])
fLatticeCoords = np.zeros(shape=latticeShape)

for index in latticeCoords:
	fLatticeCoords[index] = int(dFlatten(index))

with open("./lizList",'w') as fp:
	fp.write("dim,%d,L,%d,liz,%d\n" % (dim,L,liz)  )
	latticeCoords = np.ndindex(*[ L for i in range(dim)])
	for center in latticeCoords:
		centerCoord = fLatticeCoords[center]
		fp.write("%d  ,   %s  ,  center \n" % (centerCoord + 1,center))
		fp.write("%d   ,  %s\n" % (centerCoord + 1, center) )
		translations = it.product(*tdimlist)
		for i in translations:
			ind = tuple((np.array(center) + np.array(i) + L) % L)
			fLatticeValue = fLatticeCoords[ind]
			if centerCoord == fLatticeValue :
				continue
			else:
				fp.write("%d   ,  %s\n" % (fLatticeValue+1, ind) )


