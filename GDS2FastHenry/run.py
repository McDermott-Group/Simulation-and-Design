import ImportGDSII as gds

filebase = "Z:\\mcdermott-group\\MaskDesign\\50 Ohm JPM\\InductorSimulations\\"

g = gds.GDSII()

f1 = "Inductor_300pH_v2.gds"
g.GDS2FH(filebase+f1, 0.1, '300 pH matching network inductor')
