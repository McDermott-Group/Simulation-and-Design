from gdsCAD import *
from sys import *
from argparse import *
import argparse
import datetime

now = datetime.datetime.now()
parser = argparse.ArgumentParser(add_help=False,description='Generates the gds for the capacitor.\nRefer to the pdf file for all the parametrizations')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='shows the default value of each parameters')
parser.add_argument('-s', action="store",  dest="s", default=8, type = int, help="8")
parser.add_argument('-w', action="store", dest="w", default=4,  type=int, help="4")
parser.add_argument('-w2', action="store", dest="w2", default=3, type=int, help="3")
parser.add_argument('-g', action="store",  dest="g", default=1, type = int, help="1")
parser.add_argument('-g2', action="store", dest="g2", default=5,  type=int, help="5")
parser.add_argument('-w1', action="store", dest="w1", default=1, type=int, help="1")
parser.add_argument('-x', action="store", dest="x", default=50, type=int, help = "50")
parser.add_argument('-y', action="store",  dest="y", default=20, type = int, help = "20")
parser.add_argument('-z', action="store", dest="z", default=150,  type=int, help = "150")
args = parser.parse_args()
cell=core.Cell('top')
y_ref = 0
x_ref = 0
s = args.s
w = args.w
w2 = args.w2
g = args.g
g2 = args.g2
w1 = args.w1
x = args.x
y = args.y
z = args.z
g3 = (g2-g)/2
f= g*2
points = [(x_ref - s/2, y_ref), (x_ref + s/2, y_ref), (x_ref + s/2, y_ref - w1), (x_ref + g + w1, y_ref - w1), (x_ref + g + w1, y_ref - w1 -x), (x_ref + g2 +w2,  y_ref -w1 - x - y), (x_ref + g2+w2, y_ref-w1-x-y-z), (x_ref+g2, y_ref-w1-x-y-z), (x_ref+g2, y_ref-w1-x-0.95*y), (x_ref+g, y_ref -w1-x), (x_ref+g, y_ref-w1), (x_ref - s/2, y_ref - w1)]
points1 = [(x_ref - g- 2*w1, y_ref - w1-f), (x_ref-g-2*w1, y_ref - 2*w1-f), (x_ref -g -w1, y_ref - 2*w1-f), (x_ref -g-w1, y_ref - w1 -x), (x_ref - g2-w2, y_ref-w1-x-y), (x_ref - g2 - w2, y_ref - w1-x-y-z), (x_ref - g2, y_ref-w1-x-y-z), (x_ref-g2, y_ref-w1-x-0.95*y), (x_ref - g, y_ref - w1-x), (x_ref -g, y_ref-w1-f)]
bdy = core.Boundary(points)
bdy2 = core.Boundary(points1)
cell.add(bdy)
cell.add(bdy2)
layout = core.Layout('LIBRARY')
layout.add(cell)
layout.save('Kelly_SQUID_'+now.strftime('%H_%M_%S')+'.gds')
