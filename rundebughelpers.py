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
def runwrite(steps,fid,system):
    hoomd.option.set_notice_level(0)
    for i in range(int(steps)):
        run(1)
        for pind in range(len(system.particles)):
            p1 = system.particles[pind]
    
            u1 = p1.net_energy
            #u2 = p2.net_energy
            if u1!=0:
                #f.write("Displacement: {0}\n".format(np.linalg.norm(np.array(p2.position)-np.array(p1.position))))
                f.write("{0},{1}: pos1: {2}\n orient1: {3}\nU1: {4}\nF1: {5}\nTau1: {6}\n\n".format(i,pind,p1.position,p1.orientation,p1.net_energy,p1.net_force,p1.net_torque))
    hoomd.option.set_notice_level(2)