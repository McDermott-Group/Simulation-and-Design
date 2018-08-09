from gdsCAD import *
from sys import *
from numpy import *
import argparse
import datetime

now = datetime.datetime.now()
parser = argparse.ArgumentParser(add_help=False,description='Generates the gds for the quantum bus.\nRefer to the pdf file for all the parametrizations')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Shows the default value of each parameters')
parser.add_argument('-c', action="store",  dest="c", default=18, type = int, help = "18")
parser.add_argument('-d', action="store", dest="d", default=70,  type=int, help = "70")
parser.add_argument('-w1', action="store", dest="w1", default=4, type=int, help = "4")
parser.add_argument('-g', action="store",  dest="g", default=4, type = int, help = "4")
parser.add_argument('-g2', action="store", dest="g2", default=20,  type=int, help = "20")
parser.add_argument('-g3', action="store", dest="g3", default=16, type=int, help = "16")
parser.add_argument('-g4', action="store",  dest="g4", default=14, type = int, help = "14")
parser.add_argument('-l1', action="store", dest="l1", default=350,  type=int, help = "350")
parser.add_argument('-l2', action="store", dest="l2", default=50,  type=int, help = "50")
parser.add_argument('-x', action="store", dest="x", default=60,  type=int, help = "60")
parser.add_argument('-s', action="store", dest="s", default=8, type=int, help = "8")
parser.add_argument('-w', action="store",  dest="w", default=4, type = int, help = "4")
args = parser.parse_args()
cell=core.Cell('TOP')
layout = core.Layout('LIBRARY')

w = args.w
s = args.s
y_ref = 0
x_ref = 0
d = args.d
c = args.c
w1 = args.w1
g = args.g
g2 = args.g2
g3 = args.g3
g4 = args.g4
l1 = args.l1
l2 = args.l2
x = args.x
t = y_ref+s/2+w+g
points = [(x_ref + g + 2*w1 + g3 + l2, y_ref - g4/2 - w1), (x_ref + g + 2*w1 + g3, y_ref - g4/2 - w1), (x_ref + g + 2*w1 + g3, y_ref - s/2 - w - g -2*w1-g2), (x_ref - x, y_ref - s/2 - w - g -2*w1-g2), (x_ref - x, y_ref - s/2 - w - g), (x_ref + g, y_ref - s/2 - w - g), (x_ref + g, y_ref + s/2 + w +g), (x_ref - x, y_ref + s/2 + w +g), (x_ref-x, y_ref + s/2 + w + g + 2*w1 + g2), (x_ref + g + 2*w1 + g3, y_ref + s/2 + w + g + 2*w1 + g2), (x_ref + g + 2*w1 + g3, y_ref + g4/2 + w1), (x_ref + g + 2*w1 + g3 + l2, y_ref+ g4/2 + w1)] 
points1 = [(x_ref + g + 2*w1 + g3 + l2, y_ref - g4/2 ), (x_ref + g + w1 + g3, y_ref - g4/2 ), (x_ref + g + w1 + g3, y_ref - s/2 - w - g -w1-g2), (x_ref - x+w1, y_ref - s/2 - w - g -w1-g2), (x_ref - x+w1, y_ref - s/2 - w - g-w1), (x_ref + g+w1, y_ref - s/2 - w - g-w1), (x_ref + g+w1, y_ref + s/2 + w +g+w1), (x_ref - x+w1, y_ref + s/2 + w +g+w1), (x_ref-x+w1, y_ref + s/2 + w + g + w1 + g2), (x_ref + g + w1 + g3, y_ref + s/2 + w + g + w1 + g2), (x_ref + g + w1 + g3, y_ref + g4/2 ), (x_ref + g + 2*w1 + g3 + l2, y_ref+ g4/2 )]    
yy = y_ref - d/2 + g4/2 + w1/2    
ref1 = x_ref + g + 2*w1 + g3 + l2 + d/2 +w1/2
points2=[(ref1, yy), (ref1, yy-l1), (ref1-w1, yy-l1), (ref1-w1, yy)]    
ref2 = ref1 - w1 - g4
points3=[(ref2, yy), (ref2, yy-l1), (ref2-w1, yy-l1), (ref2-w1, yy)]    
xx = x_ref + g + 2*w1 + g3 + l2    
circ_segment = shapes.Circle((xx,yy), d/2, w1, initial_angle=0, final_angle=90, layer=1)
circ2 = shapes.Circle((xx,yy), (d/2 - w1 - g4), w1, initial_angle = 0, final_angle=90,layer=1)    
bdy = core.Boundary(points)
bdy2= core.Boundary(points1)
bdy3= core.Boundary(points2)
bdy4= core.Boundary(points3)
cell.add(circ2)
cell.add(bdy3)
cell.add(bdy4)
cell.add(bdy2)
cell.add(bdy)
cell.add(circ_segment)
layout.add(cell)
layout.save('Qbus'+now.strftime('%H_%M_%S')+'.gds')
