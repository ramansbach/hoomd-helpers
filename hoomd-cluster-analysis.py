# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 08:43:11 2017

@author: rachael
Perform a set of analyses on a GSD file, eventually with some inputs
"""
import gsd.hoomd
import numpy as np
import imp

try: 
    mk = imp.load_source('/home/mansbac2/coarsegraining/code/markov/markov.py')
except IOError:
    mk = imp.load_source('/home/rachael/home/mansbac2/coarsegraining/code/markov/markov.py')
    
#sizeprof requires a numpy 1D array of positions and a knowledge of the number of atoms in a molecule
#return peps from a function that assumes a gsd trajectory already exists
def getPos(t,traj):
    snap = traj[t]
    pos = snap.particles.position
    sz = np.shape(pos)
    peps = np.reshape(pos,[1,sz[0]*sz[1]])[0]
    return peps

#Returns a numpy array with the cluster size that each peptide is part of (eg. [2, 2, 3, 3, 1, 3, ...] means peptides 0 and 1 are part of dimers, peptides 2, 3, and 5 are part of trimers, peptide 4 is a monomer, etc...). t is the time frame of trajectory xtc to look at from run tpr, outgro is the temporary gro file to write to when getting atom positions, cutoff is the minimum distance between two atoms in a peptide to define the peptides as being part of the same cluster, and ats is the number of atoms in a peptide.
def sizeprof(t, traj, cutoff, ats):
	peps = getPos(t,traj)
	pots = range(len(peps)/3/ats)
	sizes = np.zeros(len(peps)/3/ats)
	while len(pots) > 0:
		init = pots[0]
		pots.remove(init)
		clust = mk.getClust(init, cutoff, peps, pots, ats, False) + [init]
		for at in clust:
			sizes[at] = len(clust)
			
	return sizes
#Returns a list of clusters sizes at each timestep in tlist. Each clustersize is a list of sizes 
def getsizes(tlist, traj, cutoff, ats):
	return [sizeprof(tlist[i], traj, cutoff, ats) for i in range(len(tlist))]
 
#Saves a list of cluster sizes for all time steps in tlist from run trajectory xtc.
def sizerun(tlist, traj, cutoff, ats, fnm='out.txt'):

    mk.savedclusts([getsizes(tlist, traj, cutoff, ats)], fnm=fnm)

 
if __name__ == "__main__":
    #do a little sanity checking
    #basically savesizes.py
    print "under construction"
    tlist = range(100)
    ats = 17
    cutoff = 0.5
    outputFilename = "outtest.dat"
    inputfname = "rigid5.gsd"
    traj = gsd.hoomd.open(inputfname)
    sizerun(tlist,inputfname,cutoff,ats,fnm=outputFilename)
    