import os
import numpy as np
import gdspy
import jhComponents as comps
import copy

#print 'Using gdspy module version ' + str(gdspy.__version__)

## Create a cell structure called 'CELL0' to hold some objects.
chip_name = 'Grating2_V1'
myCell = gdspy.Cell('MAIN')
### Ring Params ###pi
radius = 80
wgw = 0.8 #um waveguide width

###Grating Params###
#NORMAL grating parameters
pitch = 0.78 #um
width_grating = 35#um
nteeth = int(np.ceil(width_grating/pitch))
normal_GC = {'period': pitch,
             'number_of_teeth':nteeth,
             'fill_frac':0.37,
             'width':width_grating,
             'position': (0,0), #dummy for now,
             'direction': (0,0), #dummy
             'lda': 1.6,
             'sin_theta': np.sin(np.pi*-12/180),
             'focus_distance': 25,
             'focus_width': wgw,
             'layer': 0, #dummy
             'datatype':0} #dummy
#SPD grating parameters
pitch = 0.805 #um
width_grating = 14#um
nteeth = int(np.ceil(width_grating/pitch))
spd_GC = {'period': pitch,
             'number_of_teeth':nteeth,
             'fill_frac':0.3,
             'width':width_grating,
             'position': (0,0), #dummy for now,
             'direction': (0,0), #dummy
             'lda': 1.55,
             'sin_theta': np.sin(np.pi*-8/180),
             'focus_distance': 12,
             'focus_width': wgw,
             'layer': 0, #dummy
             'datatype':0} #dummy

#MIT grating parameters
pitch = 0.805 #um
width_grating = 12#um
nteeth = int(np.ceil(width_grating/pitch))
MIT_GC = {'period': pitch,
             'number_of_teeth':nteeth,
             'fill_frac':0.3,
             'width':width_grating,
             'position': (0,0), #dummy for now,
             'direction': (0,0), #dummy
             'lda': 1.55,
             'sin_theta': np.sin(np.pi*-8/180),
             'focus_distance': 10,
             'focus_width': wgw,
             'layer': 0, #dummy
             'datatype':0} #dummy


###Useful cells
cnormal250 = gdspy.Cell('NORMALGC250')
comps.add_grating_io(cnormal250,[0, 0],50,250,normal_GC, label='N250', layer=0, grating_layer=2)
cnormal127 = gdspy.Cell('NORMALGC127')
comps.add_grating_io(cnormal127,[0, 0], 50, 127,normal_GC, label='N127', layer=0, grating_layer=2)
cSPD = gdspy.Cell('SPDGC')
comps.add_grating_io(cSPD,[0, 0], 50, 127,spd_GC, label='S127', layer=0, grating_layer=2)
cMIT = gdspy.Cell('MITGC')
comps.add_grating_io(cMIT,[0, 0], 50, 127,MIT_GC, label='M127', layer=0, grating_layer=2)



period_x = 150
period_y = 300
start_x = 0
start_y = 0

def place_normal_tester_250(x,y):
      myCell.add(gdspy.CellReference(cnormal250,[x,y]))

def place_SPD_tester_127(x,y):
	myCell.add(gdspy.CellReference(cSPD,[x,y]))

def place_normal_tester_127(x,y):
	myCell.add(gdspy.CellReference(cnormal127,[x,y]))

def place_MIT_tester_127(x,y):
	myCell.add(gdspy.CellReference(cMIT,[x,y]))


place_normal_tester_250(start_x,start_y)
place_SPD_tester_127(start_x+period_x,start_y)
place_normal_tester_127(start_x+2*period_x,start_y)
place_MIT_tester_127(start_x+3*period_x,start_y)

place_normal_tester_250(start_x+3*period_x,start_y + period_y)
place_SPD_tester_127(start_x+period_x,start_y + period_y)
place_normal_tester_127(start_x+2*period_x,start_y + period_y)
place_MIT_tester_127(start_x,start_y + period_y)

place_normal_tester_250(start_x,start_y + 2*period_y)
place_SPD_tester_127(start_x+period_x,start_y + 2*period_y)
place_normal_tester_127(start_x+2*period_x,start_y + 2*period_y)
place_MIT_tester_127(start_x+3*period_x,start_y + 2*period_y)




##Draw device

## Create a file to save the GDSII stream.
fname = os.path.abspath(os.path.dirname(os.sys.argv[0])) + os.sep + chip_name + '.gds'
out = open(fname, 'wb')
## Write the GDSII stream in the file with all created cells (by default).
## The units we used are set to micrometers and the precision to nanometers.
gdspy.write_gds(out, unit=1.0e-6, precision=1.0e-9, cells=(myCell,cnormal250,cnormal127, cSPD, cMIT))

## Close the file.
out.close()

#Cell area

textarea = 'Area = %d um^2' % (myCell.area())
print('Sample gds file saved: ' + fname)
print('Area: ' + textarea)


## Plot all created cells (by default) using the module matplotlib, if present in your system.

gdspy.LayoutViewer()

print('DONE!')
