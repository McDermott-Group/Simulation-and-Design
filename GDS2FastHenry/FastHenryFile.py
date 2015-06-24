#!/usr/bin/python

##FastHenryFile.py
##Guilhem Ribeill -- 2/14/13

##Contains class to write data to a fast henry file, as well as classes that represent internal groundplane, point and 
##segment structures

import time

#################### 

class point:
    '''
    Point object that stores x,y,z values.
    '''
    def __init__(self,x=0.,y=0.,z=0.):
        self.x=x
        self.y=y
        self.z=z

####################        

class segment:
    '''
    Segment for FastHenry. Contains a starting node number, an end node number, a width, and a height.
    '''
    def __init__(self,start,end,w=3,h=0.1):
        self.start = start
        self.end = end
        self.w=w
        self.h=h
        
#################### 

class GroundPoly:
    '''
    Object that stores the attributes of a ground polygon object. Variables are
    three points P1,P2,P3 that define the corners of the ground plane, thickness that
    defines the thickness of the ground plane, and xseg and yseg that define the segmentation
    of the ground plane
    '''
    def __init__(self,P1,P2,P3,thickness=0.1,xseg=51,yseg=51):
        self.P1 = P1
        self.P2 = P2
        self.P3 = P3
        self.thick = thickness
        self.xseg=xseg
        self.yseg=yseg
    
    def WriteToFile(self,fileObj,num):
        '''
        Writes this ground plane object to the file given by fileObj.
        
        Inputs:
            fileObj -- python file object (file pointer should be ready to accept ground plane
            definitions
            num -- number to identify the ground plane
        Outputs:
            None.
        '''
        try:
            
            fileObj.write('* GROUND PLANE %d \n'%num)
            fileObj.write('g%d x1 = %f\ty1 = %f\tz1 = %f \n'%(num, self.P1.x, self.P1.y, self.P1.z))
            fileObj.write('+  x2 = %f\ty2 = %f\tz2 = %f \n'%(self.P2.x, self.P2.y, self.P2.z))
            fileObj.write('+  x3 = %f\ty3 = %f\tz3 = %f \n'%(self.P3.x, self.P3.y, self.P3.z))
            fileObj.write('*    thickness:\n')
            fileObj.write('+ thick = %f \n'%self.thick)
            fileObj.write('*    discretization:\n')
            fileObj.write('+ seg1 = %d seg2 = %d \n'%(self.xseg, self.yseg))
            
        except:
            raise Exception('Could not write Ground Plane #%d!!!'%num)

#################### 

class FastHenryFile:
    
    def __init__(self,filename,comments=''):
        '''
        Creates a text file with the given filename in the current path and stores the 
        file object. Also writes a simple header with any comments.
        
        Inputs: filename -- string name of current file
                commetns -- string with comments. should not contain any newlines
        Outputs: none
        '''
        
        try:
            
            self.file = open(filename, 'w')
            
        except:
            
            raise Exception('Could not create .inp file in current folder!')
        
        self.file.write('* '+filename+'\n')
        self.file.write('* Created by GDS2FastHenry '+time.strftime('%x %X')+'\n')
        self.file.write('* '+comments+'\n\n')
        
        #number of ground planes, nodes, and segments we've written to the file
        self.nGP = 1
        self.nN  = 1
        self.nS  = 1
        self.nNL = 1
        self.nSL = 1
        
        
    def WriteDefaults(self,z=0.,pdepth=0.1,nwinc=30,nhinc=10,units='um'):
        '''
        Write the default unit settings to the FastHenry file. 
        
        Inputs:
            z -- default z height of all elements in the file
            pdepth -- penetration depth of superconductors
            nwinc -- number of filaments in x-y plane
            nhinc -- number of filaments in z direction
            units -- default units of the file
        Outputs: 
            None.
            '''
        self.file.write('.Units '+units+'\n')
        self.file.write('.default nwinc=%d nhinc=%d'%(nwinc, nhinc)+'\n')
        self.file.write('.default lambda=%f'%pdepth+'\n')
        self.file.write('.default z=%f'%z+'\n')
        self.file.write('\n')
        
    def WriteTerminator(self):
        '''
        Writes the stuff you need at the end of the FastHenry file, and closes the file.
        
        Inputs:
            None.
        Outputs:
            None.
        '''
        
        self.file.write('.freq fmin=0.159154943 fmax=0.159154944 ndec=1 \n')
        self.file.write('.end\n')
        
        self.file.close()
        
        
    def WriteGroundPlanes(self, GroundPolyList):
        '''
        Write the ground plane definitions to the FastHenry file. 
        
        Inputs:
            GroundPolyList -- list of GroundPoly class objects that defines the ground planes.
        Outputs:
            none
        '''
        
        for gp in GroundPolyList:
            gp.WriteToFile(self.file, self.nGP)
            self.nGP = self.nGP + 1
        
    def WriteNodeList(self, NodeList, name=''):
        '''
        Writes a list of nodes to the file. Each node is a point object -- the z value will be ignored.
        
        Inputs:
            NodeList: list of point objects that constitute a set of nodes
            name: optional node name that will be written as a comment in the FastHenry file
        Outputs:
            None.
        '''
        
        try:
            
            self.file.write('* NODE LIST %d '%self.nNL+name+'\n')
            self.nNL = self.nNL + 1
            
            for node in NodeList:
                
                self.file.write('N%02d x=%f y=%f z=%f\n'%(self.nN, node.x, node.y, node.z))
                self.nN = self.nN + 1
                
        except:
            raise Exception('Could not write node list '+name+'!')
                
    def WriteSegmentList(self, SegmentList, name=''):
        '''
        Write a list of connecting segments stored in a list of segment objects to the file.
        
        Inputs:
            SegmentList: list of segment objects that connect nodes
            name: optional name that will be written as a comment
        Outputs:
            None.
        '''
        try:
            
            self.file.write('* SEGMENT LIST %d '%self.nSL+name+'\n')
            self.nSL = self.nSL + 1
            
            for seg in SegmentList:
                
                self.file.write('E%02d N%02d N%02d w=%f h=%f \n'%(self.nS, seg.start, seg.end, seg.w, seg.h))
                self.nS = self.nS + 1
                
        except:
            raise Exception('Could not write node list '+name+'!')
        
    def WriteExternals(self, ExtList):
        '''
        Writes externals (ports across which to compute inductance) to the file.
        
        Inputs:
            ExtList: list of tuples that define each external port
        Outputs:
            None.
        '''
        
        try:
            self.file.write('* External Ports\n')
        except:
            raise Exception('Could not write external port list!')
        
        for ext in ExtList:
            try:
                self.file.write('.external N%02d N%02d \n'%ext)
            except:
                raise Exception('Could not write external port list!')
        
            
        
        
    
            
            
        
        
        