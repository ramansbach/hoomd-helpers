# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 14:16:11 2017

@author: rachael
HOOMD run/debug helper functions
"""
import hoomd
from hoomd import *
from hoomd import md
import numpy as np
#function to determine the moment of inertia of a sphere at a distance of r from the origin, with a mass of m and a radius R
def isphere(m,R,r):
    I = (2./5.)*m*R*R+m*r*r
    return I
    
#function to determine the moment of inertia of a line of spheres spaced symmetrically about the origin, with masses m and radii R, and a spacing r
#n is the number of spheres to each side of the central sphere, so the total number is 2n+1
def isphereline(m,R,r,n):
    I = isphere(m,R,0)
    for i in range(n):
        I += 2*isphere(m,R,(i+1)*r)
    return I    
#make a quaternion of a particular angle, mainly a reminder function
#quatnernion is a rotation of theta degrees about the omega axis
def make_quat(theta,omega):
    q = np.zeros(4)
    q[0] = np.cos(theta/2.)
    q[1:4] = np.sin(theta/2.)*omega
    return q
beadMass = 1
beadR = 0.5

#run but write out certain data at every run step--very very slow and only good for debuggin
def runwrite(steps,fid,system,step=1):
    hoomd.option.set_notice_level(0)
    for i in range(int(steps)):
        run(step)
        for pind in range(len(system.particles)):
            p1 = system.particles[pind]
    
            u1 = p1.net_energy
            #u2 = p2.net_energy
            if u1!=0:
                #f.write("Displacement: {0}\n".format(np.linalg.norm(np.array(p2.position)-np.array(p1.position))))
                fid.write("{0},{1}: pos1: {2}\n orient1: {3}\nU1: {4}\nF1: {5}\nTau1: {6}\n\n".format(i,pind,p1.position,p1.orientation,p1.net_energy,p1.net_force,p1.net_torque))
    hoomd.option.set_notice_level(3)
    

#initialize five-bead rigid body
def rigid5init(system,beadMass,beadR,lbeadR,sticky_theta):
    for p in system.particles:
        #set moment of inertia based on overlapping spheres and radius = 0.5 with rigid beads laid out along x axis
        p.moment_inertia = (5*isphere(beadMass,beadR,0),isphereline(beadMass,beadR,beadR,2),isphereline(beadMass,beadR,beadR,2))
    #add secondary particle types to simulation
    system.particles.types.add('E2')
    system.particles.types.add('LB')
    system.particles.types.add('LS')
     #LB is large LJ sticky particle, LS is small LJ sticky particle
    #create rigid particles
    rigid = md.constrain.rigid()
    rigid.set_param('E',positions=[(-beadR,0,0),(beadR,0,0),(-2*beadR,0,0),(2*beadR,0,0),(0,beadR*np.cos(sticky_theta),beadR*np.sin(sticky_theta)),(0,-beadR*np.cos(sticky_theta),beadR*np.sin(sticky_theta)),(beadR,beadR*np.cos(sticky_theta),beadR*np.sin(sticky_theta)),(beadR,-beadR*np.cos(sticky_theta),beadR*np.sin(sticky_theta)),(-beadR,beadR*np.cos(sticky_theta),beadR*np.sin(sticky_theta)),(-beadR,-beadR*np.cos(sticky_theta),beadR*np.sin(sticky_theta)),(0,beadR*np.cos(sticky_theta),-beadR*np.sin(sticky_theta)),(0,-beadR*np.cos(sticky_theta),-beadR*np.sin(sticky_theta)),(beadR,beadR*np.cos(sticky_theta),-beadR*np.sin(sticky_theta)),(beadR,-beadR*np.cos(sticky_theta),-beadR*np.sin(sticky_theta)),(-beadR,beadR*np.cos(sticky_theta),-beadR*np.sin(sticky_theta)),(-beadR,-beadR*np.cos(sticky_theta),-beadR*np.sin(sticky_theta))],types=['E2','E2','LB','LB','LS','LS','LS','LS','LS','LS','LS','LS','LS','LS','LS','LS']) #type of particle has to be different from central particle because otherwise they are defined recursively
    rigid.create_bodies()
    
    #create different groups
    groupR = group.rigid_center()
    groupLB = group.type('LB')
    groupLS = group.type('LS')
    
    #set masses and diameters of beads
    for p in groupR: 
        p.mass = beadMass
        p.diameter = 2*beadR
        
    for p in groupLB:
        p.mass = beadMass
        p.diameter = 2*beadR
    
    for p in groupLS:
        p.mass = 0.0
        p.diameter = 2*lbeadR
    return (rigid,groupR,groupLB,groupLS)
    
#takes a list of types and a dictionary of keys (typeA,typeB) for epsilon values and for sigma values
#sets up the corresponding slj style pair potential
def setupSLJ(types,epsDict,sigDict,r_cut,nl,dmax):
    slj = md.pair.slj(r_cut = r_cut,nlist = nl, d_max = dmax)
    setpairs(slj,types,epsDict,sigDict)
    slj.set_params(mode='shift')
    return slj

#as previous, but with LJ instead
def setupLJ(types,epsDict,sigDict,r_cut,nl):
    lj = md.pair.lj(r_cut=r_cut,nlist=nl)
    setpairs(lj,types,epsDict,sigDict)
    lj.set_params(mode='shift')
    return lj
    
#set pair coefficients for given potential from given dictionaries
def setpairs(pot,types,epsDict,sigDict):
    for i in range(len(types)):
        for j in range(i+1):
            typeA = types[i]
            typeB = types[j]
            if (typeA,typeB) in epsDict.keys():
                epsVal = epsDict[(typeA,typeB)]
            elif (typeB,typeA) in epsDict.keys():
                epsVal = epsDict[(typeB,typeA)]
            else:
                epsVal = 0.0
            if (typeA,typeB) in sigDict.keys():
                sigVal = sigDict[(typeA,typeB)]
            elif (typeB,typeA) in sigDict.keys():
                sigVal = sigDict[(typeB,typeA)]
            else:
                sigVal = 1.0
            pot.pair_coeff.set(typeA,typeB,epsilon=epsVal,sigma=sigVal)
    