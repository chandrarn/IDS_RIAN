# Psi-TET:
# time = 27 # static max
# time = 38 # static min
# time = 63 # dynamic max
# time = 17 # dynamic min 
# NIMROD:
# time = 33450 # static max
# time = 34050 # static min
# time = 20100 # dynamic max
time = 21100 # dynamic min
n_samples = 100
chord_start = (-0.54666, 0.20548, 0.)
# chord_ends = ((0.5055, -0.1715, 0.),(0.5211, -0.1161, 0.),(0.5305, -0.0593, 0.),(0.5339, -0.0018, 0.), # all
# (0.5309, 0.0558, 0.),(0.5219, 0.1128, 0.),(0.5065, 0.1685, 0.),(0.4853, 0.2224, 0.),(0.4583, 0.2738, 0.),
# (0.4257, 0.3221, 0.),(0.3880, 0.3668, 0.),(0.3453, 0.4072, 0.),(0.2982, 0.4428, 0.),(0.2471, 0.4732, 0.),
# (0.1347, 0.5166, 0.),(0.0744, 0.5286, 0.))
# chord_ends = ((0.5055, -0.1715, 0.),(0.5211, -0.1161, 0.),(0.5305, -0.0593, 0.),(0.5339, -0.0018, 0.), # nimrod 1
# (0.5309, 0.0558, 0.),(0.5219, 0.1128, 0.),(0.5065, 0.1685, 0.),(0.4853, 0.2224, 0.))
chord_ends = ((0.4583, 0.2738, 0.),(0.4257, 0.3221, 0.),(0.3880, 0.3668, 0.),(0.3453, 0.4072, 0.),# nimrod 2
(0.2982, 0.4428, 0.),(0.2471, 0.4732, 0.),(0.1347, 0.5166, 0.),(0.0744, 0.5286, 0.))
#
# chord_ends = ((0.3880, 0.3668, 0.),(0.3453, 0.4072, 0.),(0.2982, 0.4428, 0.)) # For testing
#
# Psi-TET:
# OpenDatabase("S:\MC_IDS\Psi TET Comparisons\\aaronData_131029\out_00" + str(time) + ".xmf") # static
# OpenDatabase("S:\MC_IDS\Psi TET Comparisons\\aaronData_140106\PSITET_14kHz_xMHD\out_00" + str(time) + ".xmf") # dynamic
# NIMROD:
# OpenDatabase("S:\MC_IDS\Matlab Code\NIMROD\\vtk_nimrod\\0beta\dump_" + str(time) + "_b.vtk") # static
OpenDatabase("S:\MC_IDS\Matlab Code\NIMROD\\vtk_nimrod\\fin_beta\dump_" + str(time) + "_b.vtk") # dynamic
#
# Psi-TET:
# DefineScalarExpression("Vx", "V[0]")
# DefineScalarExpression("Vy", "V[1]")
# DefineScalarExpression("Vz", "V[2]")
# NIMROD:
DefineScalarExpression("Vx", "ve[0]")
DefineScalarExpression("Vy", "ve[1]")
DefineScalarExpression("Vz", "ve[2]")
#
# SetTimeSliderState(time)
AddPlot("Pseudocolor", "Vx")
DrawPlots()
Lineout(chord_start, chord_ends[0], n_samples)
SetActiveWindow(2)
#
# Psi-TET:
# f = open('S:\Line_tet_static_max.txt', 'w') # static max
# f = open('S:\Line_tet_static_min.txt', 'w') # static min
# f = open('S:\Line_tet_dynamic_max.txt', 'w') # dynamic max
# f = open('S:\Line_tet_dynamic_min.txt', 'w') # dynamic min
# NIMROD:
# f = open('S:\Line_nim_static_max1.txt', 'w') # static max 1
# f = open('S:\Line_nim_static_max2.txt', 'w') # static max 2
# f = open('S:\Line_nim_static_min1.txt', 'w') # static min 1
# f = open('S:\Line_nim_static_min2.txt', 'w') # static min 2
# f = open('S:\Line_nim_dynamic_max1.txt', 'w') # dynamic max 1
# f = open('S:\Line_nim_dynamic_max2.txt', 'w') # dynamic max 2
# f = open('S:\Line_nim_dynamic_min1.txt', 'w') # dynamic min 1
f = open('S:\Line_nim_dynamic_min2.txt', 'w') # dynamic min 2
print f
for (i,chord_end) in enumerate(chord_ends):
#
	t = GetOperatorOptions(0)
	t.point2 = chord_end
	SetOperatorOptions(t)
	#
	ChangeActivePlotsVar("Vx")
	Query("Average Value")
	vx = GetQueryOutputValue()
	ChangeActivePlotsVar("Vy")
	Query("Average Value")
	vy = GetQueryOutputValue()
	ChangeActivePlotsVar("Vz")
	Query("Average Value")
	vz = GetQueryOutputValue()
	print i+1, vx, vy, vz
	s = str(vx) + '\t' + str(vy) + '\t' + str(vz) + '\n'
	f.write(s)
	#
f.close()