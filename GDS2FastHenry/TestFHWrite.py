import FastHenryFile as fh

fhfile = fh.FastHenryFile('testfile.imp')

fhfile.WriteDefaults()

P1 = fh.point(1,2,3)
P2 = fh.point(1,2,3)
P3 = fh.point(1,2,3)

GP = fh.GroundPoly(P1,P2,P3)

fhfile.WriteGroundPlanes([GP])

N1 = fh.point(1, 3)
N2 = fh.point(4, 5)
N3 = fh.point(5, 6)

fhfile.WriteNodeList([N1, N2, N3], name='test 1')

S1 = fh.segment(1,2)
S2 = fh.segment(2,3)

fhfile.WriteSegmentList([S1,S2], name='test segments')

fhfile.WriteExternals([(1,3)])

fhfile.WriteTerminator()
