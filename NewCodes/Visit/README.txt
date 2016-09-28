Basic description of IDS - Simulation - VisIt scripts

--- Process_IDS_2.m ---
Pulls in processed IDS data (dat<shotnum>.mat) and converts to a .vtk database for display in VisIt.

    --- findIDSforVisit.m ---
    subfunction ALSO CALLED BY 'makeLineOutGeom_1.m' which does most of the math for calculating chords. 

    --- pts2grid.m ---
    takes in 'pts' and 'vec' for making vtk files and makes a false grid by offsetting them to either side of their original plane.
    this is necessary for making a vtk grid so that IDS vectors can be plotted with 2D slicing in VisIt. (2D slices "miss" point data in VisIt).

--- Lineout_2_Movie.m ---
This pulls in line-out files (ie: LOnim1_2_100.mat) and turns them into movie files ("shot<fake number>.mat" and "t<fake number>.mat")

--- launchVisit.py ---
This does lots of things.
1) Visualization: read in NIMROD vtk database OR Psi-TET HDF5 database and up to 2 IDS .vtk databases for simultaneous display.
2) Line Out: load in simulation database and pull lineout data, save to "LOnim1_2_100.mat" type file.
