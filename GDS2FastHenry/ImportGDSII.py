#!/usr/bin/python

##ImportGDSII.py
##Guilhem Ribeill -- 2/14/14

##Imports a GDSII file and converts it to the internal representation appropriate to write to a FastHenry file.

import FastHenryFile as fhf
from gdsii.library import Library
from gdsii.elements import *
import warnings
import os

SCALE = 1000.

class GDSII:
    
    '''
    GDSII class -- contains all functions needed to convert a gdsII file to an internal representation suited 
    to writing to a Fast Henry .imp file.
    '''
    
    def Load(self,filename,thick=0.1):
        '''
        Imports a gdsii file.
        Inputs:
            filename: GDSII file to be imported. Must be able to open in binary format.
            thick: Thickness of all elements in GDSII file. Defaults to 100nm.
        '''
        try:
            print(filename)
            stream = open(filename,'rb')
            self.gds_lib = Library.load(stream)
        except:
            raise Exception('Could not import GDSII file!!')
        
        self.thick=thick
        
######################################################
        
    def ConvertToInternal(self, skipLayers=[]):
        '''
        Converts GDSII file to internal representation.
        Here's how this works:
            1. All rectangles on GDS layer 0 are assumed to be ground planes
            2. All wires on separate layers are assumed to be connected segments
            3. Ports will be placed at the start and end of all wires.
        The GDSII cell should only have one top level cell (i.e. one cell that is not referenced by any other cell),
        as this is the one that will be convereted to the internal representation. Layers in skipLayers will be ignored.
        
        Inputs:
            skipLayers - List of layer numbers to skip over. Defaults to empty.
        Outpus:
            None.
        
        '''
        
        try:
            self.gds_lib
        except:
            raise Exception('GDSII file not imported!')
        
        if len(self.gds_lib) > 1:
            raise Exception('More than one cell in GDSII library! I don''t know how to deal with this!')
        if len(self.gds_lib) == 0:
            raise Exception('No cells in GDSII library!')
        
        cell = self.gds_lib[0]
        
        self.GroundPlaneList = []
        self.NodeListList = []
        self.SegmentListList = []
        self.ExternalList = []
        
        #counters
        nodeCounter = 1
        startNode = 0
        
        #loop through the elements in the cell
        
        for elem in cell:
            
            #print elem.xy
            
            if type(elem) is Boundary:
                
                #the element is a ground plane!
                
                if len(elem.xy) != 5:
                    warnings.warn('Found a polygon with more than 4 vertices! Skipping...')
                    continue
                
                if (elem.layer not in skipLayers) and (elem.layer == 0):
                    
                    P1 = fhf.point(elem.xy[0][0]/SCALE,elem.xy[0][1]/SCALE)
                    P2 = fhf.point(elem.xy[1][0]/SCALE,elem.xy[1][1]/SCALE)
                    P3 = fhf.point(elem.xy[2][0]/SCALE,elem.xy[2][1]/SCALE)
                    
                    self.GroundPlaneList.append(fhf.GroundPoly(P1,P2,P3, self.thick))
                
            elif type(elem) is Path:
                
                #the element is a wire
                
                if elem.layer not in skipLayers:
                    
                    startNode = nodeCounter
                    nodeList = []
                    segList = []
                    
                    for vertex in elem.xy:
                        
                        nodeList.append(fhf.point(vertex[0]/SCALE, vertex[1]/SCALE, elem.layer*self.thick))
                        nodeCounter = nodeCounter+1
                        
                    for idx in range(len(elem.xy)-1):
                        
                        segList.append(fhf.segment(startNode+idx,startNode+idx+1,elem.width/SCALE,self.thick))
                    
                    
                    self.NodeListList.append(nodeList)
                    self.SegmentListList.append(segList)
                    self.ExternalList.append((startNode,startNode+len(elem.xy)-1))        
                                    
            else:
                #don't know what to do with this, print warning
                warnStr = 'Found unknown GDSII element: '+str(type(elem))+' in cell. Skipping...'
                warnings.warn(warnStr)
                                
######################################################

    def GDS2FH(self,filename,thick=0.1,comments='',skipLayers=[]):
        '''
        Import a GDSII file and export a FastHenry file.
        
        Inputs:
            filename - string path/file to GDSII file to be imported
            comments - string of comments to include in the FastHenry file
            skipLayers - any layers to be skipped in the GDSII file
            thick - thickness of the GDSII file
        Outputs:
            None.
        '''
        
        path = os.path.splitext(filename)
        
        if path[1] != '.gds':
            raise Exception('File '+filename+' is not a GDSII file!!')
        
        print('Importing GDSII file {0}.'.format(filename))
        
        self.Load(filename,thick)
        
        print('Parsing GDSII File...')
        
        self.ConvertToInternal(skipLayers)
        
                
        fh_fn = path[0]+'.inp'
        
        fastHenry = fhf.FastHenryFile(fh_fn, comments)
        
        print('Exporting FastHenry file...')
        
        fastHenry.WriteDefaults()
        fastHenry.WriteGroundPlanes(self.GroundPlaneList)
        
        for nodes in self.NodeListList:
            fastHenry.WriteNodeList(nodes)
            
        for segs in self.SegmentListList:
            fastHenry.WriteSegmentList(segs)
        
        fastHenry.WriteExternals(self.ExternalList)
        
        fastHenry.WriteTerminator()
        
        print('Wrote file: '+fh_fn)
        print('All done!')
                
                
        
        
        
       
          
    
