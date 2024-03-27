import sys
import os
import shlex
import argparse

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__

cap = input("Enter parameters for the capacitor : (Enter blank for default values)\n")
XY = input("Enter parameters for the XY Control : (Enter blank for default values)\n")
typ = input("Enter type of SQUID : (kelly/default)\n")
Z = input("Enter parameters for the SQUID and Z control : (Enter blank for default values)\n")
Qbus = input("Enter parameters for the Quantum Bus : (Enter blank for default values)\n")
resonator = input("Enter parameters for the Resonator : (Enter blank for default values)\n")

parser = argparse.ArgumentParser(add_help=False,description='Generates the gds for the Capacitor.\nRefer to the pdf file for all the parametrizations')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='Shows the default value of each parameters')
parser.add_argument('-s', action="store",  dest="s", default=8, type = int, help="8")
parser.add_argument('-w', action="store", dest="w", default=4,  type=int, help="4")
parser.add_argument('-l', action="store", dest="l", default=135, type=int, help="135")
args = parser.parse_args(shlex.split(cap))

s = args.s
w = args.w
l = args.l

blockPrint()
os.system(('python capacitor.py -main '+cap))
os.system(('python XYcontrol.py -main -s '+str(s)+' -w '+str(w)+' -x_ref '+str((-l-w))+' '+XY))
os.system(('python Qbus.py -main -s '+str(s)+' -w '+str(w)+' -x_ref '+str((l+w))+' '+Qbus))
os.system(('python resonator.py -main -s '+str(s)+' -w '+str(w)+' -y_ref '+str((+l+w))+' '+resonator))
enablePrint()
if(typ == 'kelly'):
    os.system(('python kelly.py -main -s '+str(s)+' -w '+str(w)+' -y_ref '+str((-l-w))+' '+Z))
else:
    os.system(('python SQUID.py -main -s '+str(s)+' -w '+str(w)+' -y_ref '+str((-l))+' '+Z))

os.remove('temp.gds')
