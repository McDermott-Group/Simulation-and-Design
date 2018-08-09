from gdsCAD import *
from sys import *
from numpy import *
import argparse
import datetime

now = datetime.datetime.now()
parser = argparse.ArgumentParser(add_help=False,description='Generates the gdb for the XYcontrol.\nRefer to the pdf file for all the parametrizations')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Shows the default value of each parameters')
parser.add_argument('-c', action="store",  dest="c", default=18, type = int, help = "18")
parser.add_argument('-w1', action="store", dest="w1", default=3,  type=int, help = "3")
parser.add_argument('-w2', action="store", dest="w2", default=6, type=int, help = "6")
parser.add_argument('-g', action="store",  dest="g", default=4, type = int, help = "4")
parser.add_argument('-g2', action="store", dest="g2", default=20,  type=int, help = "20")
parser.add_argument('-x', action="store", dest="x", default=60, type=int, help = "60")
parser.add_argument('-y', action="store",  dest="y", default=30, type = int, help = "30")
parser.add_argument('-z', action="store", dest="z", default=150,  type=int, help = "150")
parser.add_argument('-s', action="store", dest="s", default=8, type=int, help = "8")
parser.add_argument('-w', action="store",  dest="w", default=4, type = int, help = "4")
args = parser.parse_args()
cell=core.Cell('TOP')
c = args.c
w1 = args.w1
w2 = args.w2
g = args.g
g2 = args.g2
x = args.x
y = args.y
z = args.z
s = args.s
w = args.w
y_ref = 0
x_ref = -c


points1=[(x_ref , y_ref + g/2 + w1),(x_ref - x , y_ref + g/2 + w1), (x_ref - x - y, y_ref + g2/2 + w2), (x_ref - x -y-z,  y_ref + g2/2 + w2), (x_ref -x -y -z, y_ref + g2/2), (x_ref - x - 0.8*y, y_ref + g2/2), (x_ref - x, y_ref +g/2), (x_ref - w1, y_ref + g/2), (x_ref - w1, y_ref - g/2), (x_ref - x, y_ref - g/2), (x_ref - x - 0.8 * y, y_ref - g2/2), (x_ref - x - y -z, y_ref - g2/2), (x_ref - x - y -z, y_ref - g2/2 - w2), (x_ref-x-y, y_ref - g2/2 - w2), (x_ref - x, y_ref - g/2 - w1), (x_ref, y_ref-g/2 - w1), (x_ref, y_ref+g/2 + w1)]
bdy = core.Boundary(points1)
cell.add(bdy)
layout = core.Layout('LIBRARY')
layout.add(cell)
layout.save('XYcontrol' + now.strftime('%H_%M_%S')+'.gds')
