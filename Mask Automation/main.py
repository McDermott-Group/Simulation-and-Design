from gdsCAD import *
from sys import *
import argparse
import sys

s = int(sys.argv[1])
w = int(sys.argv[2])
l = int(sys.argv[3])

c = 18
layout = core.Layout('LIBRARY')
cell=core.Cell('Capacitor')
def capacitor(s, w, l):
    s1 = [(l+s/2, s/2), (s/2,s/2), (s/2, l+s/2), (-s/2, l+s/2), (-s/2, s/2), (-(l+s/2), s/2), (-(l+s/2), -s/2), (-s/2, -s/2), (-s/2, -(l+s/2)), (s/2, -(l+s/2)), (s/2, -s/2), (l+s/2, -s/2), (l+s/2, s/2)]
    s2 = [(l+(s/2 + w), (s/2 + w)), ((s/2 + w),(s/2 + w)), ((s/2 + w), l+(s/2 + w)), (-(s/2 + w), l+(s/2 + w)), (-(s/2 + w), (s/2 + w)), (-(l+(s/2 + w)), (s/2 + w)), (-(l+(s/2 + w)), -(s/2 + w)), (-(s/2 + w), -(s/2 + w)), (-(s/2 + w), -(l+(s/2 + w))), ((s/2 + w), -(l+(s/2 + w))), ((s/2 + w), -(s/2 + w)), (l+(s/2 + w), -(s/2 + w)), (l+(s/2 + w), (s/2 + w))] 
    x1 = core.Path(s1, 0.01)
    x2 = core.Path(s2, 0.01)
    cell.add(x1)
    cell.add(x2)

def XYcontrol(x_ref, y_ref):
    w1 = 3
    w2 = 6
    g = 4
    g2 = 20
    x = 60
    y = 30
    z = 200
    points1=[(x_ref , y_ref + g/2 + w1),(x_ref - x , y_ref + g/2 + w1), (x_ref - x - y, y_ref + g2/2 + w2), (x_ref - x -y-z,  y_ref + g2/2 + w2), (x_ref -x -y -z, y_ref + g2/2), (x_ref - x - 0.8*y, y_ref + g2/2), (x_ref - x, y_ref +g/2), (x_ref - w1, y_ref + g/2), (x_ref - w1, y_ref - g/2), (x_ref - x, y_ref - g/2), (x_ref - x - 0.8 * y, y_ref - g2/2), (x_ref - x - y -z, y_ref - g2/2), (x_ref - x - y -z, y_ref - g2/2 - w2), (x_ref-x-y, y_ref - g2/2 - w2), (x_ref - x, y_ref - g/2 - w1), (x_ref, y_ref-g/2 - w1), (x_ref, y_ref+g/2 + w1)]
    bdy = core.Boundary(points1)
    cell.add(bdy)

def Qbus(x_ref, y_ref):
    d = 70
    c = 18
    w1 = 4
    g = 4
    g2 = 20
    g3 = 16
    g4 = 14
    l1 = 350
    l3 = 50
    x = 60
    t = y_ref+s/2+w+g
    points = [(x_ref + g + 2*w1 + g3 + l3, y_ref - g4/2 - w1), (x_ref + g + 2*w1 + g3, y_ref - g4/2 - w1), (x_ref + g + 2*w1 + g3, y_ref - s/2 - w - g -2*w1-g2), (x_ref - x, y_ref - s/2 - w - g -2*w1-g2), (x_ref - x, y_ref - s/2 - w - g), (x_ref + g, y_ref - s/2 - w - g), (x_ref + g, y_ref + s/2 + w +g), (x_ref - x, y_ref + s/2 + w +g), (x_ref-x, y_ref + s/2 + w + g + 2*w1 + g2), (x_ref + g + 2*w1 + g3, y_ref + s/2 + w + g + 2*w1 + g2), (x_ref + g + 2*w1 + g3, y_ref + g4/2 + w1), (x_ref + g + 2*w1 + g3 + l3, y_ref+ g4/2 + w1)] 
    points1 = [(x_ref + g + 2*w1 + g3 + l3, y_ref - g4/2 ), (x_ref + g + w1 + g3, y_ref - g4/2 ), (x_ref + g + w1 + g3, y_ref - s/2 - w - g -w1-g2), (x_ref - x+w1, y_ref - s/2 - w - g -w1-g2), (x_ref - x+w1, y_ref - s/2 - w - g-w1), (x_ref + g+w1, y_ref - s/2 - w - g-w1), (x_ref + g+w1, y_ref + s/2 + w +g+w1), (x_ref - x+w1, y_ref + s/2 + w +g+w1), (x_ref-x+w1, y_ref + s/2 + w + g + w1 + g2), (x_ref + g + w1 + g3, y_ref + s/2 + w + g + w1 + g2), (x_ref + g + w1 + g3, y_ref + g4/2 ), (x_ref + g + 2*w1 + g3 + l3, y_ref+ g4/2 )]
    
    yy = y_ref - d/2 + g4/2 + w1/2
    
    ref1 = x_ref + g + 2*w1 + g3 + l3 + d/2 +w1/2
    points2=[(ref1, yy), (ref1, yy-l1), (ref1-w1, yy-l1), (ref1-w1, yy)]
    
    ref2 = ref1 - w1 - g4
    points3=[(ref2, yy), (ref2, yy-l1), (ref2-w1, yy-l1), (ref2-w1, yy)]
    
    xx = x_ref + g + 2*w1 + g3 + l3
    
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

def SQUID(x_ref, y_ref):
    w2 = 3
    g = 1
    g2 = 5
    w1 = 1
    x = 50
    y = 20
    z = 150
    g3 = (g2-g)/2
    f= g*2
    points = [(x_ref - s/2, y_ref), (x_ref + s/2, y_ref), (x_ref + s/2, y_ref - w1), (x_ref + g + w1, y_ref - w1), (x_ref + g + w1, y_ref - w1 -x), (x_ref + g2 +w2,  y_ref -w1 - x - y), (x_ref + g2+w2, y_ref-w1-x-y-z), (x_ref+g2, y_ref-w1-x-y-z), (x_ref+g2, y_ref-w1-x-0.95*y), (x_ref+g, y_ref -w1-x), (x_ref+g, y_ref-w1), (x_ref - s/2, y_ref - w1)]
    points1 = [(x_ref - g- 2*w1, y_ref - w1-f), (x_ref-g-2*w1, y_ref - 2*w1-f), (x_ref -g -w1, y_ref - 2*w1-f), (x_ref -g-w1, y_ref - w1 -x), (x_ref - g2-w2, y_ref-w1-x-y), (x_ref - g2 - w2, y_ref - w1-x-y-z), (x_ref - g2, y_ref-w1-x-y-z), (x_ref-g2, y_ref-w1-x-0.95*y), (x_ref - g, y_ref - w1-x), (x_ref -g, y_ref-w1-f)]
    bdy = core.Boundary(points)
    bdy2 = core.Boundary(points1)
    cell.add(bdy)
    cell.add(bdy2)

def resonator(x_ref, y_ref): 
    g4 = 4
    g2 = 20
    g3 = 16
    g = g4 + s + w
    x = 60
    points = [(x_ref-g-2*w-g2, y_ref - x), (x_ref-g-2*w-g2, y_ref+g4+2*w+g3), (x_ref+g+2*w+g2, y_ref+g4+2*w+g3), (x_ref+g+2*w+g2, y_ref -x), (x_ref+g, y_ref -x), (x_ref+g, y_ref +g4), (x_ref-g, y_ref+g4), (x_ref-g, y_ref -x), (x_ref-g-2*w-g2, y_ref - x)]
    pth = core.Path(points, 4)
    cell.add(pth)

points = [(500,500), (-500,500), (-500,-500), (500,-500), (500,500)]
boundary = core.Path(points, 0.01)
cell.add(boundary)
layout.add(cell)

capacitor(s,w,l)
XYcontrol( -l-w-c, 0)
Qbus(l+w+4, 0)
SQUID(0, -l-w-6)
resonator(0, l+w+6)
layout.save('xmon.gds')
