import ImportGDSII as gds

g = gds.GDSII()

filebase = r"C:\FolderOfStuff"

fname = os.path.join(filebase, 'rfSQUID.gds')
thickness = 0.1
g.GDS2FH(filebase+f1, thickness, 'Comments')
