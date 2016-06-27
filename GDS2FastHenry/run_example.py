import ImportGDSII as gds

filebase = "C:\\GDSii-File-Folder"

g = gds.GDSII()

f1 = "Intersting_GDS_File.gds"
thickness = 0.1
g.GDS2FH(filebase+f1, thickness, 'Comments')
