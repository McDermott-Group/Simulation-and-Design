from gdsCAD import *
from sys import *
from numpy import *
import argparse
import datetime

now = datetime.datetime.now()
parser = argparse.ArgumentParser(add_help=False,description='Generates the gds for the SQUID and Zcontrol.\nRefer to the pdf file for all the parametrizations')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Shows the default value of each parameters')
parser.add_argument('-m', action="store",  dest="m", default=45, type = int, help = "45")
parser.add_argument('-w', action="store", dest="w", default=6,  type=int, help = "6")
parser.add_argument('-lx', action="store", dest="lx", default=4, type=int, help = "4")
parser.add_argument('-g', action="store",  dest="g", default=4, type = int, help = "4")
parser.add_argument('-g2', action="store", dest="g2", default=20,  type=int, help = "20")
parser.add_argument('-x', action="store", dest="x", default=50, type=int, help = "50")
parser.add_argument('-y', action="store",  dest="y", default=20, type = int, help = "20")
parser.add_argument('-z', action="store", dest="z", default=150,  type=int, help = "150")
parser.add_argument('-s', action="store", dest="s", default=8, type=int, help = "8")
parser.add_argument('-w2', action="store",  dest="w2", default=4, type = int, help = "6")
parser.add_argument('-l', action="store",  dest="l", default=18, type = int, help = "18")
parser.add_argument('-main', dest="main", action="store_true", default=False)
parser.add_argument('-x_ref', action="store",  dest="x_ref", default=0, type = int, help = "4")
parser.add_argument('-y_ref', action="store",  dest="y_ref", default=0, type = int, help = "4")
args = parser.parse_args()

cell=core.Cell('TOP4')
main = args.main
s = args.s
m = args.m
w = args.w
g = args.g
lx = args.lx
x_ref = args.x_ref
y_ref = args.y_ref - w/2 - g -lx

points1=[(x_ref-m/2,y_ref-w/2), (x_ref-m/2,y_ref + w/2), (x_ref + m/2,y_ref+w/2), (x_ref +m/2,y_ref-w/2)]
points2=[(x_ref-m/2-g-lx, y_ref-w/2), (x_ref-m/2-g, y_ref-w/2), (x_ref-m/2-g, y_ref + w/2 + g), (x_ref +m/2 + g, y_ref + w/2+g), (x_ref +m/2 + g, y_ref-w/2), (x_ref +m/2+g+lx, y_ref-w/2), (x_ref +m/2 + g +lx, y_ref + w/2 + g + lx), (x_ref-m/2 - g - lx, y_ref + w/2 + g + lx), (x_ref-m/2-g-lx, y_ref-w/2)]    
points3 = [(x_ref + s/2,y_ref+ w/2), (x_ref -s/2, y_ref+w/2), (x_ref-s/2,y_ref+ w/2 + g), (x_ref+s/2,y_ref+ w/2 + g)]
bdy = core.Boundary(points1)
bdy2 = core.Boundary(points2)
bdy3 = core.Boundary(points3)
cell.add(bdy)
cell.add(bdy2)
cell.add(bdy3)
    #_________________________________________________________________________________________________________________________________________
x_ref = x_ref - (m/2+3*g/2+lx)
y_ref -= (w/2+g)
l = args.l
g = args.g
g2 = args.g2
w2 = args.w2
x = args.x
y = args.y
z = args.z
g3 = (g2-g)/2
points1=[(x_ref + g/2, y_ref), (x_ref + g/2 + l, y_ref), (x_ref + g/2 +l,y_ref -w2), (x_ref + g/2 + w2, y_ref - w2), (x_ref + g/2+w2,y_ref -x), (x_ref +  g2/2 + w2, y_ref -x-y), (x_ref + g2/2 + w2,y_ref -x-y-z), (x_ref + g2/2,y_ref -x-y-z), (x_ref + g2/2,y_ref -x-y), (x_ref + g/2,y_ref -x), (x_ref + g/2, y_ref)]
points2 = [(x_ref - g/2, y_ref), (x_ref - g/2 - w2, y_ref), (x_ref - g/2 - w2, y_ref -x), (x_ref - g2/2 - w2, y_ref - x - y), (x_ref - g2/2 - w2, y_ref - x - y - z), (x_ref - g2/2, y_ref -x - y - z), (x_ref - g2/2, y_ref - x - y), (x_ref - g/2, y_ref -x), (x_ref - g/2, y_ref)]
bdy = core.Boundary(points1)
bdy2 = core.Boundary(points2)
cell.add(bdy)
cell.add(bdy2)

if main: 
    name = 'XMON_'+now.strftime('%H_%M_%S')+'.gds'
    layout = core.GdsImport('temp.gds')
    cell.add(layout['TOP3'], origin=(0,0))
    print('File created : ' + name)
else:
    name = 'SQUID_'+now.strftime('%H_%M_%S')+'.gds'
    layout = core.Layout('LIBRARY')
    
layout.add(cell)
layout.save(name)

print('\n\nFile created : ' + name)
