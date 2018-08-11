from gdsCAD import *
from sys import *
from argparse import *
import argparse
import datetime

now = datetime.datetime.now()
parser = argparse.ArgumentParser(add_help=False,description='Generates the gds for the Capacitor.\nRefer to the pdf file for all the parametrizations')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Shows the default value of each parameters')
parser.add_argument('-s', action="store",  dest="s", default=8, type = int, help="8")
parser.add_argument('-w', action="store", dest="w", default=4,  type=int, help="4")
parser.add_argument('-l', action="store", dest="l", default=135, type=int, help="135")
parser.add_argument('-main', action="store_true", default=False)
args = parser.parse_args()
l = args.l
s = args.s
w = args.w
main = args.main
cell=core.Cell('TOP')
s1 = [(l, s/2), (s/2,s/2), (s/2, l), (-s/2, l), (-s/2, s/2), (-(l), s/2), (-(l), -s/2), (-s/2, -s/2), (-s/2, -(l)), (s/2, -(l)), (s/2, -s/2), (l, -s/2), (l, s/2)]
s2 = [(l+(w), (s/2 + w)), ((s/2 + w),(s/2 + w)), ((s/2 + w), l+( w)), (-(s/2 + w), l+(w)), (-(s/2 + w), (s/2 + w)), (-(l+(w)), (s/2 + w)), (-(l+(w)), -(s/2 + w)), (-(s/2+w), -(s/2+w)), (-(s/2+w), -(l+(w))), ((s/2 + w), -(l+(w))), ((s/2 + w), -(s/2+w)), (l+(w), -(s/2+w)), (l+( w), (s/2 + w))]
s = s1+s2
x = core.Boundary(s)
cell.add(x)

if main:
    name = 'temp.gds'
else:
    name = 'capacitor'+now.strftime('%H_%M_%S')+'.gds'
    
    
layout = core.Layout('LIBRARY')
layout.add(cell)
layout.save(name)
