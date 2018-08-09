from gdsCAD import *
from sys import *
from argparse import *
import argparse
import datetime

now = datetime.datetime.now()
parser = argparse.ArgumentParser(description='generates the gdb for the capacitor')

parser.add_argument('-s', action="store",  dest="s", default=8, type = int)
parser.add_argument('-w', action="store", dest="w", default=4,  type=int)
parser.add_argument('-l', action="store", dest="l", default=135, type=int)
args = parser.parse_args()
l = args.l
s = args.s
w = args.w
cell=core.Cell('TOP')
s1 = [(l+s/2, s/2), (s/2,s/2), (s/2, l+s/2), (-s/2, l+s/2), (-s/2, s/2), (-(l+s/2), s/2), (-(l+s/2), -s/2), (-s/2, -s/2), (-s/2, -(l+s/2)), (s/2, -(l+s/2)), (s/2, -s/2), (l+s/2, -s/2), (l+s/2, s/2)]
s2 = [(l+(s/2 + w), (s/2 + w)), ((s/2 + w),(s/2 + w)), ((s/2 + w), l+(s/2 + w)), (-(s/2 + w), l+(s/2 + w)), (-(s/2 + w), (s/2 + w)), (-(l+(s/2 + w)), (s/2 + w)), (-(l+(s/2 + w)), -(s/2 + w)), (-(s/2 + w), -(s/2 + w)), (-(s/2 + w), -(l+(s/2 + w))), ((s/2 + w), -(l+(s/2 + w))), ((s/2 + w), -(s/2 + w)), (l+(s/2 + w), -(s/2 + w)), (l+(s/2 + w), (s/2 + w))] 
x1 = core.Path(s1, 0.01)
x2 = core.Path(s2, 0.01)
cell.add(x1)
cell.add(x2)


layout = core.Layout('LIBRARY')
layout.add(cell)


layout.save('capacitor'+now.strftime('%H_%M_%S')+'.gds')
~                                                                                                                                                     
~                                                                                                                                                     
~                                                                                                                                                     
~                                                   
