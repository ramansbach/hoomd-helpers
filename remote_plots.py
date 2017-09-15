# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 07:58:01 2017

@author: rachael
Functions using matplotlib to open a set of data files with the same base filename 
and opening them to plot the average with error bars of column versus column
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
plt.ioff()
#function that opens a set of files with the same base name and extracts data for plotting column X vs the avg of cols Y with errorbars in Y
#rows is an array [start end] or an array of arrays [[start1 end1], [start1 end1]] for subsampling of the data or
#starting at different times in different files
#if colX is None just use 1:length
def dataInFiles(bname,ext,nofiles,colX,colY,rows):
    #find the length of the data
    if type(rows[0]) == int:
        dlen = rows[1] - rows[0] 
    else:
        dlen = rows[0][1] - rows[0][0] 
    #initialize arrays
    if colX:
        xs = np.zeros(dlen)
    else:
        xs = np.array(range(dlen))
    
    Ys = np.zeros([nofiles,dlen])
    #iterate through each file and go through all lines, pulling the desired
    #data and putting it into xs and Ys
    for fi in range(nofiles):
        if type(rows[0]) == int:
            start = rows[0]
            end = rows[1]
        else:
            start = rows[fi][0]
            end = rows[fi][1]
        
        f = open(bname+str(fi+1)+ext)
        li = 0
        for line in f:
            if li >= start and li <= end:
                spline = line.split()
                if colX:
                    xs[li-start] = float(spline[colX])
                Ys[fi][li-start] = float(spline[colY])
            li+=1
        f.close()
    return (xs,Ys)
    
def plotMeanErrs(xs,Ys,labelX,labelY,figname):
    ys = np.mean(Ys,axis=0)
    yerrs = np.std(Ys,axis=0)
    plt.figure()
    plt.errorbar(xs,ys,yerr=yerrs)
    plt.xlabel(labelX)
    plt.ylabel(labelY)
    plt.savefig(figname)
    plt.close()
    
#function that plots column X vs the cols Y
def plotMultiData(xs,Ys,labelX,labelY,figname):
    nofiles = np.shape(Ys)[0]

    plt.figure()
    for fi in range(nofiles):
        plt.plot(xs,Ys[fi,:])
    plt.xlabel(labelX)
    plt.ylabel(labelY)
    plt.savefig(figname)
    plt.close()
    
