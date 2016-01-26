import ImportGDSII as gds

filebase = "C:\\Users\\5107-1\\Desktop\\Inductor Test\\"

g = gds.GDSII()

f1 = "rfSquidGradiometer.gds"
g.GDS2FH(filebase+f1, 0.1, 'Big gradiometer and flux bias coil for rf Squid readout')
