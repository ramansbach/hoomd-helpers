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
    mk = imp.load_source('mk','/home/mansbac2/coarsegraining/code/markov/markov.py')
except IOError:
    mk = imp.load_source('mk','/home/rachael/cluster/home/mansbac2/coarsegraining/code/markov/markov.py')
    
#sizeprof requires a numpy 1D array of positions and a knowledge of the number of atoms in a molecule
#return peps from a function that assumes a gsd trajectory already exists
def getPos(t,traj,subset='all'):

    snap = traj[t]

    pos = snap.particles.position
    #rearrange so they're all segregated according to the body they're in
    binds = np.argsort(snap.particles.body)
    pos = pos[binds]
    #get a subset of the positions of a particular type and return those
    if subset != 'all':
        typeid = snap.particles.types.index(subset)
        types = snap.particles.typeid[binds]
        pos = pos[np.where(types==typeid)[0]]
    
    sz = np.shape(pos)
    peps = np.reshape(pos,[1,sz[0]*sz[1]])[0]
    return peps

#Returns a numpy array with the cluster size that each peptide is part of (eg. [2, 2, 3, 3, 1, 3, ...] means peptides 0 and 1 are part of dimers, peptides 2, 3, and 5 are part of trimers, peptide 4 is a monomer, etc...). t is the time frame of trajectory xtc to look at from run tpr, outgro is the temporary gro file to write to when getting atom positions, cutoff is the minimum distance between two atoms in a peptide to define the peptides as being part of the same cluster, and ats is the number of atoms in a peptide.
def sizeprof(t, traj, cutoff, ats, metric,types='all'):
    if metric == 'contact':
        peps = getPos(t,traj)
        pots = range(len(peps)/3/ats)
        sizes = np.zeros(len(peps)/3/ats)
        while len(pots) > 0:
            init = pots[0]
            pots.remove(init)
            clust = mk.getClust(init, cutoff, peps, pots, ats, False, metric) + [init]
            for at in clust:
                sizes[at] = len(clust)
    elif metric == 'optical':
        peps = getPos(t,traj,types)
        pots = range(len(peps)/3/ats)
        sizes = np.zeros(len(peps)/3/ats)
        while len(pots) > 0:
            init = pots[0]
            pots.remove(init)
            clust = mk.getClust(init, cutoff, peps, pots, ats, False, metric) + [init]
            for at in clust:
                sizes[at] = len(clust)  
    return sizes
#Returns a list of clusters sizes at each timestep in tlist. Each clustersize is a list of sizes 
def getsizes(tlist, traj, cutoff, ats, metric,types='all'):
	return [sizeprof(tlist[i], traj, cutoff, ats, metric,types) for i in range(len(tlist))]
 
#Saves a list of cluster sizes for all time steps in tlist from run trajectory .
def sizerun(tlist, traj, cutoff, ats, metric, fnm='out.txt',types='all'):
    if metric == 'contact':
        mk.savedclusts([getsizes(tlist, traj, cutoff, ats, metric)], fnm=fnm)
    elif metric == 'optical':
        mk.savedclusts([getsizes(tlist, traj, cutoff, ats, metric,types)], fnm=fnm)
    else:
        print "This is not a recognized metric."
    
 
def clustMorph(t,traj,cutoff,ats):
    #return the index of which clusters each peptide is in
    #also the eigenvectors and eigenvalues of the gyration tensor of each cluster
    #also the hydrodynamic radius
    #start = time.clock();
    snap = traj[t]
    box_length = snap.configuration.box[0:3]
    peps = getPos(t,traj)
    #print box_length

    pots = range(len(peps)/3/ats)
    inds = np.zeros(len(peps)/3/ats)
    ind = 1
    eigvals = np.array([])
    eigvecs = np.array([])
    ms = np.array([])
    Rhs = np.array([])
    Rgs = np.array([])
    while len(pots) > 0:
        init = pots[0]
        pots.remove(init)
        clusts = mk.getClust(init,cutoff,peps,pots,ats,False) + [init]
        #clusts is a list of peptides that are found in the cluster
        #each index in clusts corresponds to the indices index*ats*3:(index+1)*ats*3 in peps
        pepList = np.zeros(len(clusts)*ats*3)
        curr = 0
        #mass = len(clusts);
        for clust in clusts:
            inds[clust] = ind
            pepList[curr*ats*3:(curr+1)*ats*3]=peps[clust*ats*3:(clust+1)*ats*3]
            curr+=1
        pepList = mk.fixPBC(pepList,box_length,ats,cutoff)
        #writeGro('mol_'+str(ind)+'.gro',pepList,box_length)
        gyrationTensor = mk.gyrTens(pepList,box_length)
        eigstuff = np.linalg.eig(gyrationTensor)
        ind+=1
        eigMorph = np.zeros(5)
        eigval = np.sort(eigstuff[0])
        eigMorph[0] = eigval[0]
        eigMorph[1] = eigval[1]
        eigMorph[2] = eigval[2]
        eigMorph[3] = eigMorph[0]+eigMorph[1]+eigMorph[2] #Rg^2
        eigMorph[4] = 1.5*eigMorph[2]-0.5*eigMorph[3]
        eigvals = np.append(eigvals,eigMorph)
        eigOrder = np.argsort(eigstuff[0])
        eigvec = eigstuff[1]
        eigvect = np.append(eigvec[:,eigOrder[2]],eigvec[:,eigOrder[1]])
        eigvect = np.append(eigvect,eigvec[:,eigOrder[0]])
        eigvecs = np.append(eigvecs,eigvect)
        #eigvecs = np.append(eigvecs,eigstuff[1])
        Rh = mk.hydroRad(pepList)
        Rhs = np.append(Rhs,Rh)
        #mstuff = np.zeros(3)
        mass = float(len(clusts))
        #print Rh
        #print mass
        Rgs = np.append(Rgs,eigMorph[3])
        ms = np.append(ms,mass)
        #ms = np.append(ms,np.array([mass,eigMorph[3],Rh]))
    #end = time.clock()
    #t = end - start
    #print "time: ", t
    return (inds,eigvals,eigvecs,Rhs,ms,Rgs)

 
if __name__ == "__main__":
    #do a little sanity checking
    #basically savesizes.py
    
    tlist = range(100)
    ats = 17
    cutoff = 0.5
    outputFilename = "contactsizes.dat"
    inputfname = "/home/rachael/coarsegraining/hoomd/patchytest/analysis_testing/rigid5.gsd"
    traj = gsd.hoomd.open(inputfname)
    #contact clusters
    #sizerun(tlist,traj,cutoff,ats,'contact',fnm=outputFilename)
    #optical clusters
    optats = 12
    outputFname = "opticalsizes.dat"
    aromType = u'LS'
    sizerun(tlist,traj,cutoff,optats,'optical',outputFname,aromType)
    #testClusts.py
    indsfname = "clustInds.dat"
    scalarfname = "clustScalars.dat"
    vecfname = "clustVecs.dat"
    hydrfname = "clustRh.dat"
    mRfname = "mass_v_Rg.dat"
    indsf = open(indsfname,'w')
    scalarf = open(scalarfname,'w')
    vecf = open(vecfname,'w')
    hydrf = open(hydrfname,'w')
    mRf = open(mRfname,'w')
    for t in tlist:
        (inds,eigvals,eigvecs,Rhs,ms,Rgs) = clustMorph(t,traj,cutoff,ats)
        #print inds
        #print eigvals
        #print eigvecs
        #print Rhs
        linei = ''
        lines = ''
        linev = ''
        liner = ''
        linem = ''
        for i in range(len(inds)):
            linei+=str(t)+' '
            for j in range(ats):
                #linds[j+i*ats] = inds[i]
                linei += str(inds[i])+' '
            linei+='\n'
        for i in range(0,len(eigvals),5):
            #for j in range(5):
            lines += str(t)+' '+str(eigvals[i])+' '+ str(eigvals[i+1])+' ' + str(eigvals[i+2])+' '+str(eigvals[i+3])+' '+str(eigvals[i+4])+'\n'
        for i in range(len(ms)):
            linem += str(ms[i])+' '+str(Rhs[i])+' '+str(Rgs[i])+'\n'
            #linem += '\n'		
        for i in range(len(Rhs)):
            liner += str(t)+' '+str(Rhs[i])+'\n'
        for i in range(0,len(eigvecs),9):
            linev += str(t) + ' ' + str(eigvecs[i]) + ' ' + str(eigvecs[i+1]) + ' ' + str(eigvecs[i+2]) + ' ' + str(eigvecs[i+3]) + ' ' + str(eigvecs[i+4]) + ' ' + str(eigvecs[i+5]) + ' ' + str(eigvecs[i+6]) + ' ' + str(eigvecs[i+7]) + ' ' + str(eigvecs[8]) + '\n'
        indsf.write(linei)
        scalarf.write(lines)
        vecf.write(linev)
        hydrf.write(liner)
        mRf.write(linem)	
    hydrf.close()
    indsf.close()
    scalarf.close()
    vecf.close()