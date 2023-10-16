python
import sys
from glob import glob
fname = sys.argv[1]+'-vib-mode*.pdb'
lst = glob(fname)
lst.sort()
for i in lst: cmd.load(i)
# http://www.zonums.com/online/color_ramp/
my_colors_rgb = [
    [0,0,255], [10,0,244], [ 21,0,233], [31,0,223], [42,0,212],
    [53,0,201], [63,0,191], [74,0,180], [85,0,170], [95,0,159],
    [106,0,148], [116,0,138], [ 127,0,127], [138,0,116], [148,0,106],
    [ 159,0,95], [170,0,85], [180,0,74], [191,0,63], [201,0,53],
    [ 212,0,42], [223,0,31], [233,0,21], [244,0,10], [255,0,0]
]
my_colors_name = []
for i in range( len(my_colors_rgb)):
    my_colors_name.append( "cc" + str(i))
    cmd.set_color( my_colors_name[i], my_colors_rgb[i])
color_ramp = []
for i in my_colors_name:
    color_ramp.append("white_" + i)
for i in cmd.get_object_list():
    cmd.cartoon('putty', i)
for i, c in zip( cmd.get_object_list(), color_ramp):
    cmd.spectrum('b', c, i + ' and name CA or (organic and not hydrogens) or metals')
python end
