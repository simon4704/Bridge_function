# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 16:25:15 2022

@author: simon
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt

os.chdir('C:\\Users\\simon\\OneDrive - Danmarks Tekniske Universitet\\Min mappe - oneD\\Fritid\\Programmer\\BridgeFunction')


#%% Point generator

def GeneratePoints(N):
    points = np.zeros([N,2])
    for i in range(0, N):
        points[i, 0] = np.random.random()
        points[i, 1] = np.random.random()
    
    return points



#%% Plotter

def PlotPoints(points):
    plt.figure()
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid()
    plt.plot(points[:,0], points[:,1], color='black', marker='o', linestyle='')
    plt.show()



#%% Bridge basic

# Basic definition calculating all distances

def BridgeBasic(points, r):
    rsq = r**2
    dim = np.size(points,0)
    distanceMatrix = np.zeros([dim, dim])
    overlapMatrix = np.zeros([dim, dim]).astype(int)

    for i in range(0, dim):
        xi = points[i,0]
        yi = points[i,1]
        for j in range(0, dim):
            xj = points[j,0]
            yj = points[j,1]
            
            thisDistsq = (xj - xi)**2 + (yj - yi)**2
            distanceMatrix[i,j] = thisDistsq
            
            if thisDistsq <= rsq:
                overlapMatrix[i,j] = 1
    
    return distanceMatrix, overlapMatrix
            


#%% Bridge2

# Uses that the distance matrix is symmetric

def Bridge2(points, r):
    rsq = r**2
    dim = np.size(points,0)
    distanceMatrix = np.zeros([dim, dim])
    overlapMatrix = np.zeros([dim, dim]).astype(int)

    jrange = 1
    for i in range(0, dim):
        xi = points[i,0]
        yi = points[i,1]
        for j in range(0, jrange):
            xj = points[j,0]
            yj = points[j,1]
            
            thisDistsq = (xj - xi)**2 + (yj - yi)**2
            distanceMatrix[i,j] = thisDistsq
            
            if thisDistsq <= rsq:
                overlapMatrix[i,j] = 1
        
        jrange += 1
    
    distanceMatrix = distanceMatrix + np.transpose(distanceMatrix)
    overlapMatrix = overlapMatrix + np.transpose(overlapMatrix)
    return distanceMatrix, overlapMatrix



#%% Point struct

class PointStruct:
    def __init__(self):
        self.points = np.zeros([1,2]) # Point coordinates
        self.indexes = [0] # Index of point from original list
    
    def AddPoint(self, point, index):
        self.points = np.vstack((self.points, point))
        self.indexes.append(index)



#%%

# Uses symmetry and splits area into cells
# Requires points to be organized into new arrays, this is done with
    # "OrganizePointsInCells()"
# ncells is the number of cells per row and col, so if ncells = 3 it will
# create a 3x3 grid of cells, resulting in 9 cells!

def OrganizePointsInCells(inputPoints, r, ncells):
    b = 2*r
    a = -(b*ncells - b - 1)/ncells
    
    print("a=",a," b=",b)
    
    if a < 0:
        print("r or ncells is too large")
        return
    
    pointStructs = []
    for i in range(0, ncells**2 + 1):       
        pointStructs.append(PointStruct())
    
    
    for i, inputPoint in enumerate(inputPoints):
        pointx = inputPoint[0]
        pointy = inputPoint[1]
        
        cellIndexx = ncells
        cellIndexy = ncells
        
        for j in range(0, ncells):
            lowerBound = (a+b)*j
            upperBound = a + (a+b)*j
            
            
            if lowerBound < pointx and pointx < upperBound:
                cellIndexx = j
            
            if lowerBound < pointy and pointy < upperBound:
                cellIndexy = j
        
        if cellIndexx == ncells or cellIndexy == ncells:
            cellIndexTotal = ncells**2
        else:
            cellIndexTotal = ncells * cellIndexy + cellIndexx
        
        pointStructs[cellIndexTotal].AddPoint(inputPoint, i)
            
    return pointStructs




def Bridge3(inputPoints, r, ncells):
    rsq = r**2
    dim = np.size(inputPoints,0)
    distanceMatrix = np.zeros([dim, dim])
    overlapMatrix = np.zeros([dim, dim]).astype(int)
    
    pointStructs = OrganizePointsInCells(inputPoints, r, ncells)
    
    for pointStruct in pointStructs[0:ncells**2]:
        jrange = 2
        for i in range(1, len(pointStruct.indexes)):
            xi = pointStruct.points[i,0]
            yi = pointStruct.points[i,1]
            for j in range(1, jrange):
                xj = pointStruct.points[j,0]
                yj = pointStruct.points[j,1]
                
                thisDistsq = (xj - xi)**2 + (yj - yi)**2
                
                thisIndexi = pointStruct.indexes[i]
                thisIndexj = pointStruct.indexes[j]
                
                distanceMatrix[thisIndexi,thisIndexj] = thisDistsq
                
                if thisDistsq <= rsq:
                    overlapMatrix[thisIndexi, thisIndexj] = 1
            
            jrange += 1
    
    
    lastPointStruct = pointStructs[-1]
    for i in range(1, len(lastPointStruct.indexes)):
        xi = lastPointStruct.points[i,0]
        yi = lastPointStruct.points[i,1]
        for pointStruct in pointStructs[0:ncells*2]:
            for j in range(1, len(pointStruct.indexes)):
                xj = pointStruct.points[j,0]
                yj = pointStruct.points[j,1]
                
                
                thisDistsq = (xj - xi)**2 + (yj - yi)**2
                
                thisIndexi = lastPointStruct.indexes[i]
                thisIndexj = pointStruct.indexes[j]
                
                distanceMatrix[thisIndexi,thisIndexj] = thisDistsq
                
                if thisDistsq <= rsq:
                    overlapMatrix[thisIndexi, thisIndexj] = 1
    
    
    jrange = 2
    for i in range(1, len(lastPointStruct.indexes)):
        xi = lastPointStruct.points[i,0]
        yi = lastPointStruct.points[i,1]
        
        for j in range(1, jrange):
            xj = lastPointStruct.points[j,0]
            yj = lastPointStruct.points[j,1]
            
            thisDistsq = (xj - xi)**2 + (yj - yi)**2
            
            thisIndexi = lastPointStruct.indexes[i]
            thisIndexj = lastPointStruct.indexes[j]
            
            distanceMatrix[thisIndexi,thisIndexj] = thisDistsq
            
            if thisDistsq <= rsq:
                overlapMatrix[thisIndexi, thisIndexj] = 1
        
        jrange += 1
        
    
    distanceMatrix = distanceMatrix + np.transpose(distanceMatrix)
    overlapMatrix = overlapMatrix + np.transpose(overlapMatrix)
    return distanceMatrix, overlapMatrix



#%% Parameter definitions

N = 4000

r = 0.02

points = GeneratePoints(N)



#%% Test of "BridgeBasic"

t0 = time.time()
M, OM = BridgeBasic(points, r)
t1 = time.time()

dt = t1 - t0
print(dt)



#%% Test of "Bridge2"

t0 = time.time()
M, OM = Bridge2(points, r)
t1 = time.time()

dt = t1 - t0
print(dt)


#%% Test of "Bridge3"

t0 = time.time()
M, OM = Bridge3(points, r, 4)
t1 = time.time()

dt = t1 - t0
print(dt)


