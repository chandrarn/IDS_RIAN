# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 10:33:29 2014

@author: aaron
"""

class idsVectors:
    
    
    
    def __init__(self, idsShot, idsTP):
        from visit import *
        OpenDatabase("~/IDS/Visit/IDSvtk/ids_" + idsShot + idsTP + "_*.vtk database")
        idsTS = GetActiveTimeSlider()
        SetTimeSliderState(idsTimes[0])
        
    #    # Get number of points
    #    reader.SetFileName("~/IDS/Visit/IDSvtk/ids_" + idsShot1 + idsTP1 + "_1.vtk")
    #    reader.ReadAllScalarsOn()  # Activate the reading of all scalars
    #    reader.Update()
    #    data = reader.GetOutput()
        
        AddPlot("Vector", "V")
        v = VectorAttributes()
    #    v.nVectors = data.GetNumberOfPoints()
        v.useLegend = 0
        v.autoScale = 0
        v.scale = vScale
        v.colorByMag = 0
        v.glyphLocation = 0
        SetPlotOptions(v)
        SetDefaultPlotOptions(v) # ensures future instances use these settings
        return idsTS
