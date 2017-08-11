# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 11:31:54 2017

@author: rachael
Helper functions for HOOMD visualization
"""
import hoomd
import gsd.hoomd
import numpy as np
#given a matrix of orientation quaternions and an axis, return a set of vectors in 3D space
#that come from rotating the axis by the quaternion
#r' = rvec + 2 q * (qvec x rvec) + 2 qvec x (qvec x rvec)
def qrotmat(qmat,ax):
    rmat = np.zeros([np.shape(qmat)[0],3])
    for qi in range(np.shape(qmat)[0]):
        rmat[qi,:] = ax + 2 * qmat[qi,0] * np.cross(qmat[qi,1:4],ax) + 2.0 * np.cross(qmat[qi,1:4],np.cross(qmat[qi,1:4],ax))
    return rmat    
    
#given a trajectory name, produce and write a new trajectory with unwrapped coordinates for visualization purposes
def unwrap_write(trajname,outname):
    print "under debugging"
    traj = gsd.hoomd.open(trajname)
    trajunwrap = gsd.hoomd.open(outname,mode='wb')
    #trajunwrap.extend([unwrap(traj[i],i) for i in range(len(traj))])
    for i in range(len(traj)):
        newsnap = unwrap(traj[i],i)
        trajunwrap.append(newsnap)

#takes a snapshot and returns it with positions unwrapped
def unwrap(snap,i):
    print "under debugging"
    lx = snap.configuration.box[0]
    ly = snap.configuration.box[1]
    lz = snap.configuration.box[2]
    ls = np.array([[lx],[ly],[lz]])
    pos = snap.particles.position
    img = snap.particles.image
    if np.sum(np.sum(img)):
        print "outside box images"
    posun = pos.copy()
    for k in range(3):
        posun[:,k] = pos[:,k] + img[:,k] * ls[k]
    posun[:,1] = posun[:,1] + pos[:,1] * pos[:,2] * img[:,2] * lz
    posun[:,0] = posun[:,0] + pos[:,0] * pos[:,1] * img[:,1] * ly + pos[:,0] * pos[:,2] * img[:,2] * lz
    newsnap = gsd.hoomd.Snapshot()
    newsnap.configuration.step = snap.configuration.step
    newsnap.particles.N = snap.particles.N
    newsnap.particles.position = posun
    newsnap.configuration.box = snap.configuration.box
    newsnap.configuration.dimensions = snap.configuration.dimensions
    newsnap.particles.diameter = snap.particles.diameter
    newsnap.particles.image = np.zeros(np.shape(snap.particles.image))
    newsnap.particles.orientation = snap.particles.orientation
    newsnap.particles.types = snap.particles.types
    newsnap.particles.typeid = snap.particles.typeid
    return newsnap

#takes a snapshot and generates a new snapshot with twice as many particles
#also has positions, set the positions to the original ones and also to slightly along the
#patch axis which can be extracted from the orientation for each one
#set diameters to 2/3 of large particle size
#set location of new particles to a distance Delta along the patch axis (a1) according to the law of cosines
#phi' = phi R/r & Delta^2 = R^2 + r^2 - 2 R r cos(phi'/2 - phi/2)
#keep also orientations (the same) and diameters
#phi is subtended patch angle
def add_patch(snap,phi,patch_ax):
    newtypes = np.concatenate((snap.particles.types,[u'P']))
    Ds = snap.particles.diameter
    ds = (5./6.) * Ds
    Rs = Ds/2.0
    rs = ds/2.0
    Deltas = np.sqrt(Rs * Rs + rs * Rs - 2.0 * Rs * rs * np.cos(phi*(Rs/rs)/2.0 - phi/2.0))
    
    newpos = snap.particles.position + np.reshape(Deltas,[len(Deltas),1]) * qrotmat(snap.particles.orientation,patch_ax)
    
    newsnap = gsd.hoomd.Snapshot()
    newsnap.configuration.step = snap.configuration.step
    newsnap.particles.N = 2*snap.particles.N
    newsnap.configuration.box = snap.configuration.box
    newsnap.configuration.dimensions = snap.configuration.dimensions
    newsnap.particles.diameter = np.concatenate((Ds,ds))
    newsnap.particles.orientation = np.concatenate((snap.particles.orientation,snap.particles.orientation))
    newsnap.particles.image = np.concatenate((snap.particles.image,snap.particles.image))
    newsnap.particles.types = newtypes
    newsnap.particles.typeid = np.concatenate((snap.particles.typeid,snap.particles.typeid+1))
    newsnap.particles.position = np.concatenate((snap.particles.position,newpos))
    return newsnap
#function that takes in a particular trajectory and writes out a new trajectory
#with dummy particles added at each timestep to be read into Ovito for visualization
def ovito_write(trajname,outname,phi,patch_ax):
    print "under construction"
    traj = gsd.hoomd.open(trajname)
    trajovito = gsd.hoomd.open(outname,mode='wb')
    for i in range(len(traj)):
        newsnap = add_patch(traj[i],phi,patch_ax)
        trajovito.append(newsnap)


    
if __name__ == "__main__":
    #do a simple test
    trajinput = "/home/rachael/coarsegraining/hoomd/patchytest/analysis_testing/rigid5.gsd"
    trajoutput = "/home/rachael/coarsegraining/hoomd/patchytest/analysis_testing/rigid5_unwrap.gsd"
    #unwrap_write(trajinput,trajoutput)
    traj = gsd.hoomd.open(trajinput)
    s1 = traj[0]
    ax = np.array([0.,0.,1.])
    os = s1.particles.orientation
    print qrotmat(os,ax)