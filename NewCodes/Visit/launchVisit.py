# -*- coding: utf-8 -*-
from numpy import *
from scipy import *
import scipy.io as sio
import time
import sys
sys.path.append("/home/aaron/Software/visit2_7_1.linux-x86_64/2.7.1/linux-x86_64/lib/site-packages")
import vtk
#reader = vtk.vtkDataSetReader()
reader = vtk.vtkUnstructuredGridReader()
#import visit
#visit.Launch()
from visit import *

""" --------------------------------------------------
INPUT SETTINGS
-------------------------------------------------- """
"""
#Looking into injector mouth, 129530
#"""
#vScale = 3e-6 # Vector scaling
#nVectors = 2000 # number of vectors to plot
#pMin = -5e4 # Pseudocolor velocity min
#pMax = 5e4 # Pseudocolor velocity max
#
#plotNim = 1
#nimPath = "fin_beta_ramp/dump_*_b.vtk"
#nimTimes = (20100, 21100)
#
#plotTet = 0
#tetPath = "aaronData_140106/out_*.xmf" # First Dynamic PSI-TET run
#tetTimes = (63, 17)
#
#plotIds1 = 1
#idsShot1 = "12979310" # string
#idsTP1 = "p" # "t" for toroidal fiber, "p" for poloidal fiber
#
#plotIds2 = 1
#idsShot2 = "12979310" # string
#idsTP2 = "t" # "t" for toroidal fiber, "p" for poloidal fiber
#
#idsTimes = (199, 203) # shift down by 1 from vtk file number
#idsvectorColor = (255, 0, 0, 0) # "firebrick" red
#idslineWidth = 0
#
## Window Information
#x0 = 2600
#y0 = 50
#width = 2400
#height = 1000
#
## SliceAtttributes
#torPlane = 0 # true for midplane slice
#polPlane = 1 # true for poloidal plane slice
#theta = 135.
#xAxistitle = "R [m]"
#yAxistitle = "Z [m]"
#boxAtt = (-0.6, 0.6, -0.6, 0.6, -0.45, 0.45)
#viewportCoords = (0.25, 0.95, 0.15, 0.85)
#windowCoords = (-0.6, 0.6, -0.45, 0.45)
#
## Line Out Settings
#makeLineOut = 0
#nSamples = 1000
#lineOutV = 1 # make Velocity LineOut
#lineOutn = 1 # make density LineOut
#lineOutT = 1 # make temperature LineOut
## Geometry settings, related to 'makeLineOutGeom' Matlab code
#tor = 1
#config = 3
#
#databaseInfoFlag = 1
#timeInfoFlag = 1

"""
Toroidal midplane, 129499
"""
vScale = 3e-6 # Vector scaling
nVectors = 2000 # number of vectors to plot
pMin = -5e4 # Pseudocolor velocity min
pMax = 5e4 # Pseudocolor velocity max

plotNim = 0
nimPath = "fin_beta_ramp"
nimTimes = range(110, 294)
#nimTimes = range(0, 294)

plotTet = 1
tetPath = "aaronData_141022" # First Dynamic PSI-TET run
#tetTimes = tuple(map(tuple, range(1, 75))) # all times from this 
#tetTimes = range(193, 209)
tetTimes = range(0, 209)

plotIds1 = 0
idsShot1 = "12949910" # string
idsTP1 = "t" # "t" for toroidal fiber, "p" for poloidal fiber

plotIds2 = 0
idsShot2 = "12949910" # string
idsTP2 = "p" # "t" for toroidal fiber, "p" for poloidal fiber

idsTimes = (167, 170) # shift down by 1 from vtk file number
idsvectorColor = (255, 0, 0, 0) # "firebrick" red
idslineWidth = 0

# Window Information
x0 = 2600
y0 = 50
width = 2400
height = 1000

# SliceAtttributes
torPlane = 1 # true for midplane slice
polPlane = 0 # true for poloidal plane slice
theta = 135
boxAtt = (-0.6, 0.6, -0.6, 0.6, -0.6, 0.6)
viewportCoords = (0.25, 0.95, 0.15, 0.85)
windowCoords = (-0.6, 0.7, -0.6, 0.7)

# Line Out Settings
makeLineOut = 1
#nSamples = 1000
lineOutV = 1 # make Velocity LineOut
lineOutn = 1 # make density LineOut
lineOutT = 1 # make temperature LineOut
# Geometry settings, related to 'makeLineOutGeom' Matlab code
tor = 1
config = 2

databaseInfoFlag = 0
timeInfoFlag = 1

""" --------------------------------------------------
Begin actual Scrip
-------------------------------------------------- """

def idsVectors(idsShot, idsTP):
    OpenDatabase("~/IDS/Visit/IDSvtk/ids_" + idsShot + idsTP + "_*.vtk database")
    idsTS = GetActiveTimeSlider()
    SetTimeSliderState(idsTimes[0])
    
    # Get number of points
    reader.SetFileName("/home/aaron/IDS/Visit/IDSvtk/ids_" + idsShot + idsTP + "_1.vtk")
    reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    reader.Update()
    data = reader.GetOutput()
    
    AddPlot("Vector", "V")
    v = VectorAttributes()
    v.nVectors = data.GetNumberOfPoints() / 2
    v.useLegend = 0
    v.autoScale = 0
    v.scale = vScale
    v.colorByMag = 0
    v.glyphLocation = 0
    v.vectorColor = idsvectorColor
    v.lineWidth = idslineWidth
    SetPlotOptions(v)
    SetDefaultPlotOptions(v) # ensures future instances use these settings
    return idsTS


Launch()

SetWindowArea(x0, y0, width, height)
a = GetAnnotationAttributes()
a.userInfoFlag = 0 # disable user info
a.databaseInfoFlag = databaseInfoFlag
a.timeInfoFlag = timeInfoFlag
if polPlane:
    a.axes2D.xAxis.title.userTitle = 1
    a.axes2D.xAxis.title.title = xAxistitle
    a.axes2D.yAxis.title.userTitle = 1
    a.axes2D.yAxis.title.title = yAxistitle
SetAnnotationAttributes(a)

if plotNim:
    OpenDatabase("~/IDS/NimrodData/" + nimPath + "/DUMP/dump_*_b.vtk database")
    nimTS = GetActiveTimeSlider()
    SetTimeSliderState(nimTimes[0])
    simTimes = nimTimes # used for saving LineOut files
    simCode = 'nim'
    simFile = '/home/aaron/IDS/NimrodData/' + nimPath + '/LO/'
    DefineScalarExpression("Vx", "ve[0]")
    DefineScalarExpression("Vy", "ve[1]")
    DefineScalarExpression("Vz", "ve[2]")
    DefineVectorExpression("V", "ve")
    if makeLineOut and lineOutn:
        DefineScalarExpression("n", "nd")
    if makeLineOut and lineOutT:
        DefineScalarExpression("T", "tion")
    
if plotTet:
    OpenDatabase("~/IDS/PsitetData/" + tetPath + "/DUMP/out_*.xmf database")
    tetTS = GetActiveTimeSlider()
    SetTimeSliderState(tetTimes[0])
    simTimes = tetTimes # used for saving LineOut files
    simCode = 'tet'
    simFile = '/home/aaron/IDS/PsitetData/' + tetPath + '/LO/'
    DefineScalarExpression("Vx", "V[0]")
    DefineScalarExpression("Vy", "V[1]")
    DefineScalarExpression("Vz", "V[2]")
    if makeLineOut and lineOutn:
        DefineScalarExpression("n", "N")
    # 'T' already exists as a variable in PSI-TET
    
if plotNim or plotTet:
    # Calculate component of V normal to slice
    if torPlane:
        DefineScalarExpression("Vn", "Vz")
    if polPlane:
         thetaRad = deg2rad(theta)
         xHat = -sin(thetaRad)
         yHat = cos(thetaRad)
         DefineScalarExpression("Vn", "V[0] * " + str(xHat) + " + V[1] * " + str(yHat))
        
    AddPlot("Pseudocolor", "Vn")
    p = PseudocolorAttributes()
    p.minFlag = 1
    p.min = pMin
    p.maxFlag = 1
    p.max = pMax
    SetPlotOptions(p)
    SetDefaultPlotOptions(p) # ensures future instances use these settings
    
    plotName = GetPlotList().GetPlots(0).plotName
    l = GetAnnotationObject(plotName)
    l.xScale = 0.4
    l.yScale = 0.7
    l.fontHeight = 0.05
    l.numberFormat = "%+5.f"
    l.drawTitle = 0
    l.drawMinMax = 0
    
    AddPlot("Vector", "V")
    v = VectorAttributes()
    v.nVectors = nVectors
    v.useLegend = 0
    v.autoScale = 0
    v.scale = vScale
    v.colorByMag = 0
    v.glyphLocation = 1
    SetPlotOptions(v)

    if makeLineOut:
        # load start and end points for LineOut
        contents = sio.loadmat('/home/aaron/IDS/Geometry/coords' + str(tor) + str(config))
        originArr = contents['origin']
        origin = tuple(map(tuple, originArr)) # must convert numpy array to tuple
        origin = origin[0] # this is so stupid
        originArr = originArr[0] # this is so stupid
        ptsArr = contents['pts']
        pts = tuple(map(tuple, ptsArr))
        chan_range = contents['chan_range'] # can be left as an array 
        
        # Calculate length of each chord
        dl = zeros((ptsArr.shape[0], 1)) # preallocate
        nSamples = zeros((1, ptsArr.shape[0])) # preallocate
        for (m, ptsx) in enumerate(ptsArr[:, 0]):
            dl[m, 0] = 1e3 * math.sqrt((ptsArr[m, 0] - originArr[0]) ** 2 +
                (ptsArr[m, 1] - originArr[1]) ** 2 +
                (ptsArr[m, 2] - originArr[2]) ** 2)
            nSamples[0, m] = dl[m, 0]
        nSamples = nSamples.astype(int64)
       
        DrawPlots()
        if lineOutV and not (lineOutn and lineOutT):
            lineMode = "V"
            nVars = 3
            Lineout(origin, pts[0], ("Vx", "Vy", "Vz"), nSamples[0,0])
        if lineOutV and lineOutn and not lineOutT:
            lineMode = "Vn"
            nVars = 4
            Lineout(origin, pts[0], ("Vx", "Vy", "Vz", "n"), nSamples[0,0])
        if lineOutV and lineOutT and not lineOutn:
            lineMode = "VT"
            nVars = 4
            Lineout(origin, pts[0], ("Vx", "Vy", "Vz", "T"), nSamples[0,0])
        if lineOutV and lineOutn and lineOutT:
            lineMode = "VnT"
            nVars = 5
            Lineout(origin, pts[0], ("Vx", "Vy", "Vz", "n", "T"), nSamples[0,0])
        

if plotIds1:
    idsTS1 = idsVectors(idsShot1, idsTP1)
    
if plotIds2:
    idsTS2 = idsVectors(idsShot2, idsTP2)
        
if not makeLineOut: # refine the plots for display
        
    AddOperator("Slice", 1) # "1" applies to all plots
    s = SliceAttributes()
    s.project2d = 1
    #s.originPoint = (0, 0, 1e-4)
    s.originIntercept = 1e-6
    if torPlane:
        s.normal = (0, 0, 1)
    if polPlane:
        s.theta = theta
        s.axisType = 4
    SetOperatorOptions(s, 0, 1)
    SetDefaultOperatorOptions(s) # ensures future instances use these settings
    
    va = View2DAttributes()
    va.viewportCoords = viewportCoords
    va.windowCoords = windowCoords
    SetView2D(va)
    
    AddOperator("Box", 1)
    b = BoxAttributes()
    b.amount = 1
    b.minx = boxAtt[0]
    b.maxx = boxAtt[1]
    b.miny = boxAtt[2]
    b.maxy = boxAtt[3]
    b.minz = boxAtt[4]
    b.maxz = boxAtt[5]
    SetOperatorOptions(b, 0, 1)

for m in range(0, len(simTimes)): # Time Advance
    
    if plotNim:
        SetActiveTimeSlider(nimTS)
        SetTimeSliderState(nimTimes[m])
        
    if plotTet:
        SetActiveTimeSlider(tetTS)
        SetTimeSliderState(tetTimes[m])
        
    if plotIds1:
        SetActiveTimeSlider(idsTS1)
        SetTimeSliderState(idsTimes[m])
        
    if plotIds2:
        SetActiveTimeSlider(idsTS2)
        SetTimeSliderState(idsTimes[m])
        
    if makeLineOut:
        SetActiveWindow(2)
        lineVals = []
        for (i, pt) in enumerate(pts):
            for k in range(nVars):
                SetActivePlots(k)
                t = GetOperatorOptions(0)
                t.point2 = pt
                t.numberOfSamplePoints = nSamples[0, i]
                SetOperatorOptions(t)
 
            SetActivePlots(0)
            vx = GetPlotInformation()["Curve"]
            SetActivePlots(1)
            vy = GetPlotInformation()["Curve"]
            SetActivePlots(2)
            vz = GetPlotInformation()["Curve"]
            if lineMode is "V":
                lineVals.append({'vx': vx, 'vy': vy, 'vz': vz})
            if lineMode is "Vn":
                SetActivePlots(3)
                n = GetPlotInformation()["Curve"]
                lineVals.append({'vx': vx, 'vy': vy, 'vz': vz, 'n': n})
            if lineMode is "VT":
                SetActivePlots(3)
                T = GetPlotInformation()["Curve"]
                lineVals.append({'vx': vx, 'vy': vy, 'vz': vz, 'T': T})
            if lineMode is "VnT":
                SetActivePlots(3)
                n = GetPlotInformation()["Curve"]
                SetActivePlots(4)
                T = GetPlotInformation()["Curve"]
                lineVals.append({'vx': vx, 'vy': vy, 'vz': vz, 'n': n, 'T': T})
                
        sio.savemat(simFile + 'LO' + simCode + str(tor) + '_' + str(config) 
                + '_' + str(simTimes[m]).zfill(4) + '.mat', 
                {'lineVals': lineVals, 'chan_range': chan_range})
        
    DrawPlots()
    time.sleep(5)
#




