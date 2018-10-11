#########################################################################
#       (C) 2017 Department of Petroleum Engineering,                   # 
#       Univeristy of Louisiana at Lafayette, Lafayette, US.            #
#                                                                       #
# This code is released under the terms of the BSD license, and thus    #
# free for commercial and research use. Feel free to use the code into  #
# your own project with a PROPER REFERENCE.                             #
#                                                                       #
# PyGRDECL Code                                                         #
# Author: Bin Wang                                                      # 
# Email: binwang.0213@gmail.com                                         # 
#########################################################################

import operator
import collections
import numpy as np

#Find intersection of lines
#https://github.com/ideasman42/isect_segments-bentley_ottmann
#from poly_point_isect import *


import matplotlib.pyplot as plt
#Line plot
from matplotlib.path import Path
import matplotlib.patches as Patches
#Polygon
import matplotlib
from matplotlib.collections import PatchCollection

#Shapely Geometry Computation Engine
try:
    from shapely.geometry import LineString,Point,MultiLineString,MultiPoint,Polygon
    from shapely.ops import split,polygonize,nearest_points
except ImportError:
    warnings.warn("No shapely module loaded.")


class FaultProcess:
    def __init__(self,GRDECL=[]):
        """Fault Process module (2D fault lines detect and process) 
        (assuming fault will penetrate all layers in Z direction)
        
        All fault lines are reference coordinates in terms of grid, (0,0) - (GRDECL_Data.NX,GRDECL_Data.NY)
        Shapely library will be used for fast and robust 2D geometry computation
        https://shapely.readthedocs.io/en/stable/manual.html

        Arguments
        ---------
        GRDECL_Data     -- Petrel geoglogy model class
        BoundaryLines   -- Polylines along boundary (may cut by fault)
        FaultLines      -- Internal Fault Lines (may cut by each other)
        IntersectPts    -- Intersection Points for all lines (including Boundary and Fault)

        [Optional]
        BoundaryLines_Split   -- Polylines along boundary, cutted by extended fault line
        FaultLines_Split      -- Internal Fault Lines, extended until hit the boundary
        IntersectPts_Split    -- Intersection Points, more intersection point added
        SplitPolygons         -- Splitted polygon divided by extended Faultlines
                                 When extend fault lines to the boundary:
                                 1. More intersection points added
                                 2. boundary line split by these additional intersection point
                                 3. More fault lines added if extended fault line intersectes

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        self.NumFaultLines=0
        self.GRDECL_Data=GRDECL
        self.BoundaryLines=[]
        self.FaultLines=[]
        self.IntersectPts=[]

        self.SplitPolygons=[]

    def findFaultLines(self):
        #1. Find internal fault lines
        NX,NY=self.GRDECL_Data.NX,self.GRDECL_Data.NY
        
        RawFaultVerts=[]
        for j in range(NY):
            for i in range(NX):
                CellFault=self.GRDECL_Data.findCellFault([i,j,0])
                BdMarker,BDFaces=self.GRDECL_Data.isBoundaryCell([i,j,1])
                FaultIndex=sum(CellFault)-BdMarker
                if(FaultIndex>0): #We find a internal fault
                    #print('(%d,%d)'%(i,j),CellFault,FaultIndex,'BoundaryCell',BdMarker,BDFaces)
                    vert=deriveFaultLoc(self.GRDECL_Data,i,j,CellFault,BdMarker,BDFaces)
                    for vi in vert:
                        RawFaultVerts.append(vi)
        
        self.FaultLines=computeInternalFaultLine(self.GRDECL_Data,RawFaultVerts)
        self.NumFaultLines=len(self.FaultLines)
        #2. Intersection point are the end points of fault lines
        self.IntersectPts=[]
        for line in self.FaultLines:
            self.IntersectPts.append(line[0])
            self.IntersectPts.append(line[-1])
        print('[FaultProcess] Found %d faults.'%(self.NumFaultLines))
        
        #debug
        #self.plotLines(self.FaultLines)

    def findBoundaryLines(self):
        #Find the bounary line which will cut by faults lines
        #ACTNUM is not considered here
        NX,NY=self.GRDECL_Data.NX,self.GRDECL_Data.NY

        if(self.NumFaultLines==0):
            print("Please find the Fault lines first! Boundary Lines cutted by it")
            return

        #Raw boundary lines (currently only support straight boundary)
        Edge1=((0, 0), (NX, 0))
        Edge2=((NX, 0), (NX, NY))
        Edge3=((NX, NY), (0, NY))
        Edge4=((0, NY), (0, 0))
        RawBoundaryLines=MultiLineString([Edge1,Edge2,Edge3,Edge4]) #Shapely object

        #Split boundary lines
        unique_intersectPts=list(set(self.IntersectPts))
        #print(unique_intersectPts)
        pts=MultiPoint(unique_intersectPts)#Shapely object
        result=split(RawBoundaryLines, pts)#Shapely function

        #Convert shapely object to list of tuple
        self.BoundaryLines=Shapely2List_MultiLineString(result)
        #print(self.BoundaryLines)
        #self.plotLines(self.FaultLines,self.BoundaryLines,self.IntersectPts)

    
    def SplitDomainByFault(self):
        # Extend fault line to remove the hanging point
        # hanging point is not intersect with other fault and boundary
        #https://gis.stackexchange.com/questions/232771/splitting-polygon-by-linestring-in-geodjango
        #https://gis.stackexchange.com/questions/283352/most-efficient-way-to-split-a-polygon-with-lines-c-api

        #Extend Fault lines
        BoundaryLine_Splitted,FaultLine_Extend, NewIntersectPts=self.extendFaultLines()
        
        #Find the 2D sub domain
        results = polygonize(MultiLineString(BoundaryLine_Splitted+FaultLine_Extend))
        self.SplitPolygons=Shapely2List_MultiPolygon(results)

        for i,p in enumerate(self.SplitPolygons):#Remove unnecessary points in polygon
            self.SplitPolygons[i]=simplify_Polygon(p)
        
        for i,poly in enumerate(self.SplitPolygons):#Convert float to int
            poly=list(reversed(poly)) #Reverse to anti-clock wise order
            for j,node in enumerate(poly):
                self.SplitPolygons[i][j]=(int(node[0]),int(node[1]))

        print('[FaultProcess] Domain is splitted as %d polygons.'%(len(self.SplitPolygons)))

        #debug
        #self.plotLines(BoundaryLine_Splitted,FaultLine_Extend,NewIntersectPts)
        
        self.plotSplittedDomain()

    def extendFaultLines(self):
        """Extend Fault lines
            When extend fault lines to the boundary:
                1. More intersection points added
                2. boundary line split by these additional intersection point
                3. More fault lines added if extended fault line intersectes
        Arguments
        ---------
        FaultLines   -- [dict] Unique fault line data [Verts][LocID] 
        IntersectPts -- intersection points (end points) for fault lines

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        debug=0
        OldLines=MultiLineString(self.BoundaryLines+self.FaultLines)
        FaultLine_Extend=self.FaultLines[:]
        BoundaryLine_Splitted=self.BoundaryLines[:]

        #Step 1. Extend Faults Lines
        ExtendLineIDs=[]
        ExtendLines=[]
        NewIntersectPts=[]
        for i,Line in enumerate(self.FaultLines):
            #if(i>25):
            #    continue
            flag=0
            StartPoint,EndPoint=Line[0],Line[-1]
            countSP=self.IntersectPts.count(StartPoint)
            countEP=self.IntersectPts.count(EndPoint)
            
            NewLine=Line[:]
            NewEndPoint=[]
            if(debug): print('Before',NewLine,countSP,countEP)
            if(countSP==1 and isBoundaryVert(self.GRDECL_Data,StartPoint)==False):
                #if(debug):print('SV ',StartPoint,'is a hanging vert')
                NewEndPoint=extend_FaultLines(self.GRDECL_Data,Line, FaultLine_Extend,'StartPoint')
                NewLine=NewEndPoint + NewLine #+NewLine[1:] 
                NewIntersectPts.append(NewEndPoint[0])
                flag=1
            if(countEP==1 and isBoundaryVert(self.GRDECL_Data,EndPoint)==False):
                #if(debug): print('EV ',EndPoint,'is a hanging vert')
                NewEndPoint=extend_FaultLines(self.GRDECL_Data,Line, FaultLine_Extend,'EndPoint')
                NewLine=NewLine + NewEndPoint#NewLine[:-1]+NewEndPoint#
                NewIntersectPts.append(NewEndPoint[0])
                flag=1
            if(flag==1):
                if(debug): print('After',NewLine)
                ExtendLines.append(NewLine)
                ExtendLineIDs.append(i)
                FaultLine_Extend[i]=NewLine

        
        if(debug): 
            print('Added EndPoint',sorted(NewIntersectPts))
            print('Extended Lines',ExtendLineIDs)

        if(len(ExtendLines)>0): #We have extenable lines    
            #Step 2. Find the intersection points between newly extended lines
            NewLines=MultiLineString(ExtendLines)
            PossibileIntersectPts=[]
            for i,line_i in enumerate(NewLines):
                for j,line_j in enumerate(NewLines):
                    if(j>i): 
                        result=line_i.intersection(line_j)
                        if(result.geom_type in ['LineString','Point']):
                            #print('--------',result.geom_type)
                            result=list(result.coords)
                        else:
                            if(len(result)>0): 
                                #print('--------',result.geom_type)
                                if(result.geom_type=='MultiPoint'):
                                    result=Shapely2List_MultiPoint(result)
                                else:
                                    print("!!!!!!!!!!!!!!May have problem...Check extendFaultLines!!")
                        if(len(result)>0): 
                            #print(result)
                            #print(i,j,line_i,line_j)
                            PossibileIntersectPts+=result
            print('Added %d new intersection pts'%(len(PossibileIntersectPts)))
            NewIntersectPts+=PossibileIntersectPts
            
            #Step 3. Split the old line in terms of new intersection point
            if(len(NewIntersectPts)>0):
                result=split(MultiLineString(FaultLine_Extend), MultiPoint(NewIntersectPts))
                FaultLine_Extend=Shapely2List_MultiLineString(result)
                result=split(MultiLineString(self.BoundaryLines), MultiPoint(NewIntersectPts))
                BoundaryLine_Splitted=Shapely2List_MultiLineString(result)

            #debug
            #self.plotLines(BoundaryLine_Splitted,FaultLine_Extend)
        

        return BoundaryLine_Splitted,FaultLine_Extend,self.IntersectPts+NewIntersectPts



    def plotLines(self,bdlines=[],faultlines=[],endpoints=[]):
        
        #Plot the fault line map
        if(len(bdlines)+len(faultlines)==0):
            BoundaryLabels=['Edge'+str(i) for i in range(len(self.BoundaryLines))]
            FaultLabels=['Fault'+str(i) for i in range(len(self.FaultLines))]
            Lines=self.BoundaryLines+self.FaultLines
        else:
            BoundaryLabels=['Edge'+str(i) for i in range(len(bdlines))]
            FaultLabels=['Fault'+str(i) for i in range(len(faultlines))]
            Lines=bdlines+faultlines

        Labels=BoundaryLabels+FaultLabels
        DrawPath(Lines,Labels,endpoints)
        #print(Lines,Labels)
    
    def plotSplittedDomain(self):
        DrawPolygons(self.SplitPolygons)


def computeInternalFaultLine(GRDECL_Data,RawFaultVerts):
    """Connect fault vertex based on its frequence 

    FaultVerts
    Fault_coords    LocID (Counter)
    (3,5)             4
    (2,1)             3
    (0,1)             2
    (6,4)             1

    Arguments
    ---------
    FaultVerts         -- [dict] Unique fault line data [Verts][LocID] 
    LocID              -- Special ID shows the location of the vertex and 
                          the number faults shared by this vertex (only when LocalID>2)
                            Counter>2   start/end vertex on fault line, shared by [Counter] fault lines
                            Counter=2   normal vertex on fault line, shared by [1] fault line
                            Counter=1   start/end vertex on fault line, shared by only [1] fault line
    StartEndVerts      -- [list] a array of start/end vertices
    SearchMarker       -- [list] Search marker used to mark the if this vertices has been 
                                 picked as a fault line vertex, 

    Author:Bin Wang(binwang.0213@gmail.com)
    Date: Sep. 2018
    """
    debug=0

    #Collect the unqiue fault verts and save its number of shared faults
    FaultVerts=collections.Counter(RawFaultVerts)
    FaultVerts={vert:int(FaultVerts[vert]/2) for vert in FaultVerts}    
    FaultVerts=collections.OrderedDict(sorted(FaultVerts.items()))
    Verts=list(FaultVerts.keys())
    LocID=list(FaultVerts.values())

    #print(Verts)


    NumVerts=len(Verts)
    #Collect Start and End Points and corrsponding markers
    StartEndVerts=[]
    SearchMarker=np.array([1 for i in range(NumVerts)])
    for i in range(NumVerts):
        if(LocID[i]==1 or LocID[i]>2):#This is Start vertx or end vertx of a fault line
            StartEndVerts.append(Verts[i])
            SearchMarker[i]=LocID[i]

    if(debug): print('EndVerts',StartEndVerts)
    #if(debug): print("All verts",Verts)

    def calcRelativeDist(vert1,vert2):
        #Calc releative distance between two vert
        #return abs(vert1[0]-vert2[0])+abs(vert1[1]-vert2[1])
        return abs(vert1[0]-vert2[0])+abs(vert1[1]-vert2[1])
    
    def calcMinOffset(vert1,vert2):
        #Calc the min offset in x and y direction
        return min(abs(vert1[0]-vert2[0]),abs(vert1[1]-vert2[1]))

    def countNumBoundaryVert(StartIDs):
        #Prefer to start at the boundary to avoid bug
        count=0
        for id in StartIDs:
            if(isBoundaryVert(GRDECL_Data,Verts[id])==True):
                count+=1
        return count



    FaultLines=[]
    for line_i in range(100):
        if(len(np.nonzero(SearchMarker)[0])>0):#The first nonzero element index
            StartIDs=np.nonzero(SearchMarker)[0]
            #StartID=np.nonzero(SearchMarker)[0][0]
        else:
            if(debug): print("Searching Complete!")
            break
        
        #for StartID in StartIDs:#Search the next start Vert
        NumBoundaryVert=countNumBoundaryVert(StartIDs)
        for StartID in StartIDs:
            if(Verts[StartID] in StartEndVerts):
                if(isBoundaryVert(GRDECL_Data,Verts[StartID])==False and NumBoundaryVert>0): #Prefer to start on the boundary
                    continue
                verts=[Verts[StartID]]
                SearchMarker[StartID]-=1 #NumOfLine-=1
                StartSearchID=StartID
                if(debug): print("Searching Line",line_i,'Start@',Verts[StartID],SearchMarker[StartID])
                break
            if(StartID==StartIDs[-1]):
                print('\n[Error] Can not find start/End Point')
                print(StartID,Verts[StartID])
                print(SearchMarker)
                startIDs=np.nonzero(SearchMarker)[0]
                print([Verts[i] for i in startIDs])
                break
        EndID=-1
        loopcount=0
        while EndID==-1:#Keep searching until found all connected lines
            if(loopcount>0):#Start Shrink the search range after the first loop
                StartSearchID=np.nonzero(SearchMarker)[0][0]
                #EndSearchID=SearchIDRange[1]
            for i in range(StartSearchID,NumVerts):
                if(SearchMarker[i]>0):
                    dist=calcRelativeDist(Verts[i],verts[-1])
                    #print('StartSearchID',StartSearchID)
                    #print("Checking ",i,Verts[i],SearchMarker[i],'Last Vert',verts[-1],dist,'LoopCount',loopcount)
                    if(dist==1 and Verts[i] not in verts):#One Line can not pass start node two times
                        #print(i,dist,Verts[i],verts[-1])
                        if(isFaultEdge(GRDECL_Data,(Verts[i],verts[-1]))==False):
                            if(debug): print("!!!This is not a fault edge!!!",(Verts[i],verts[-1]))
                            continue
                        SearchMarker[i]-=1
                        verts.append(Verts[i])
                        if(verts[-1] in StartEndVerts):
                            EndID=i
                            break
                    if(calcMinOffset(Verts[i],verts[-1])>max(GRDECL_Data.NX,GRDECL_Data.NY)):#We are far away
                        if(debug): print("NewSearchLoop!",Verts[i],verts[-1])
                        loopcount+=1
                        break
            loopcount+=1
            #Can not find in the last loop,reset start search id at first available vert
            if(loopcount>(GRDECL_Data.NX+GRDECL_Data.NY)):#Something bad happend
                #print('Bad Line!!')
                #print('LineID',line_i,'StartPoint',StartID,Verts[StartID])
                break 

        #print(SearchMarker)
        if(EndID!=-1):
            FaultLines.append(verts)
        if(debug):
            print("Line%d Start@"%(line_i),StartID,Verts[StartID],' - End@',EndID,Verts[EndID])
            print(verts)

    #Simplify the Faultlines, striaght line will be simplifed as two vert
    NumFaultLines=len(FaultLines) 
    for i in range(NumFaultLines):
        NumVerts=len(FaultLines[i])
        Length=int(calcDist(FaultLines[i][0],FaultLines[i][-1]))
        if(NumVerts==Length+1):
            FaultLines[i]=[FaultLines[i][0],FaultLines[i][-1]]
        if(debug): print("Line%d NumVerts=%d Length=%d"%(i,NumVerts,Length))
    if(debug): 
        print(FaultLines)

    return FaultLines
    

def deriveFaultLoc(GRDECL_Data,i,j,CellFault,BdMarker,BDFaces):
    #Derive the Fault Coords based on input Fault condition
    debug=0
    #if(i==1 and j==10): debug=1
    
    #vert=set()
    vert=[]
    #Find the real fault face
    if(BdMarker):#This is a boundary cell
        if('X-' in BDFaces): CellFault[0]=False
        if('X+' in BDFaces): CellFault[1]=False
        if('Y-' in BDFaces): CellFault[2]=False
        if('Y+' in BDFaces): CellFault[3]=False
    
    if(CellFault[0]==True): 
        #vert.add((i-1+1,j))
        #vert.add((i-1+1,j+1))
        vert.append((i-1+1,j))
        vert.append((i-1+1,j+1))
    if(CellFault[1]==True): 
        #vert.add((i+1,j))
        #vert.add((i+1,j+1))
        vert.append((i+1,j))
        vert.append((i+1,j+1))
    if(CellFault[2]==True): 
        #vert.add((i,j-1+1))
        #vert.add((i+1,j-1+1))
        vert.append((i,j-1+1))
        vert.append((i+1,j-1+1))
    if(CellFault[3]==True):
        #vert.add((i,j+1))
        #vert.add((i+1,j+1))
        vert.append((i,j+1))
        vert.append((i+1,j+1))
    
    #if((10,3)in vert):
    #   debug=1

    if(debug):
        print("ij(%d,%d)"%(i,j))
        print('Modified Fault Face Marker',CellFault)
        print('Fault Coord=',vert)

    return list(vert)

        
#--------------------Auxilary Func-------------------
def isFaultEdge(GRDECL_Data,edge):
    #test if a cell edge (shared by two cell) is a fault edge
    p1,p2=edge[0],edge[1]

    if(abs(p1[0]-p2[0])<1e-10):#Vertical Line
        #print('Vertical Line')
        Cell_left=(p1[0]-1,min(p1[1],p2[1]),0)
        return GRDECL_Data.findCellFault(Cell_left)[1] #[X-,X+,Y-,Y+]
    else:
        #print('Horizontal Line')
        Cell_down=(min(p1[0],p2[0]),p1[1]-1,0)
        return GRDECL_Data.findCellFault(Cell_down)[3] #[X-,X+,Y-,Y+]
    
def isBoundaryVert(GRDECL_Data,vert):
    #test if a vert is a boundary vert
    #print(0,GRDECL_Data.NX,0,GRDECL_Data.NY)
    if(vert[0]>0 and vert[0]<GRDECL_Data.NX and
       vert[1]>0 and vert[1]<GRDECL_Data.NY):
       return False
    return True

def extend_FaultLines(GRDECL_Data,line,OldFaults,startfrom='StartPoint or EndPoint'):
    #extend a line from its start point or end point
    #the end point must on the boundary
    debug=0
    
    if(startfrom=='StartPoint'):
        p1,p2=line[0],line[1]
    if(startfrom=='EndPoint'):
        p1,p2=line[-1],line[-2]

    if(abs(p1[0]-p2[0])<1e-10):
        if(debug): print("Line along Y direction")
        if(p1[1]-p2[1]<0):#Y- direction
            NewEndPoint=(p1[0],0)
            NextPoint=(p1[0],p1[1]-0.00001)
        if(p1[1]-p2[1]>0):#Y- direction
            NewEndPoint=(p1[0],GRDECL_Data.NY)
            NextPoint=(p1[0],p1[1]+0.00001)
    if(abs(p1[1]-p2[1])<1e-10):
        if(debug): print("Line along X direction")
        if(p1[0]-p2[0]<0):#X- direction
            NewEndPoint=(0,p1[1])
            NextPoint=(p1[0]-0.00001,p1[1])
        if(p1[0]-p2[0]>0):#X+ direction
            NewEndPoint=(GRDECL_Data.NX,p1[1])
            NextPoint=(p1[0]+0.00001,p1[1])
    
    #Check is hit old fault line and upadte the new endpoint
    if(debug): print('P2P1',(p2,p1),'ExtendSeg',NextPoint,NewEndPoint)
    ExtendedSegment=LineString([NextPoint,NewEndPoint])
    OldFaults=MultiLineString(OldFaults)
    objects=ExtendedSegment.intersection(OldFaults)
    
    if(objects.is_empty==False):#We have hit point
        #print('HitGeometry',objects,objects.geom_type)
        if(objects.geom_type in ['LineString','Point']):
            pts=nearest_points(Point(p1),objects)
            pts=Shapely2List_MultiPoint(pts)[1]
            #print('NearestPoint',pts)
            #pts=sorted(list(objects.coords))[0]
        elif(objects.geom_type in ['MultiLineString','MultiPoint','GeometryCollection']):
            pts=nearest_points(Point(p1),objects)
            pts=Shapely2List_MultiPoint(pts)[1]
            #print('NearestPts',pts)
            #pts=Shapely2List_MultiLineString(objects)
            #pts=sorted([j for i in pts for j in i])[0]
        else:
            print('Unkonwn shapely type',objects.geom_type,objects)
        pts=(int(pts[0]),int(pts[1]))
        if(debug): print('HitPoints',pts)
        NewEndPoint=pts

    #NewLine=sorted(line+[NewEndPoint])    
    return [NewEndPoint]

def simplify_Polygon(polygon):
    # Simplify polygon by merging stright line points into two points

    if(polygon[0]==polygon[-1]):#Remove the round end connect polygon
        test_polygon=polygon[:-1]
    else:
        test_polygon=polygon[:]
    #print(test_polygon)
    area=Polygon(test_polygon).area #Shapely polygon
    NumNodes=len(test_polygon)
    
    RemoveablePts=[]
    for i in range(NumNodes):
        temp=Polygon(test_polygon[:i] + test_polygon[(i + 1):])
        #print(temp,area,temp.area)
        if(abs(temp.area-area)<1e-10):#This node will no change the area
            #print(polygon[i],'is removeable',area)
            RemoveablePts.append(polygon[i])

    #Remove these points
    polygon=list(polygon)#convert tuple to list if applicable
    for pts in RemoveablePts:
        polygon.remove(pts)

    return polygon

def deriveFaultCellSide(edge,poly):
    '''Derive the cell location and side of a edge

    A edge is always shared by two cell, 
    e.g Edge15 shared by Y+ or Y- cells of (1,1,0) and (1,0,0) respectivly
    --------------
    | Y+ (1,1,0) |
    1---.---.----5
    | Y- (1,0,0) |
    --------------

    Testing of cell center is within the subdomain poly

    Author:Bin Wang(binwang.0213@gmail.com)
    Date: Sep. 2018
    '''
    debug=0
    p1,p2=edge[0],edge[1]

    if(abs(p1[0]-p2[0])<1e-10):#Vertical Line
        if(debug): print("Vertical Line")
        maxY=max(p1[1],p2[1])
        CellLeft=(p1[0]-0.5,maxY-0.5)
        CellRight=(p1[0]+0.5,maxY-0.5)
        if(debug): print(CellLeft,CellRight)
        if(point_in_polygon(CellLeft,poly)):
            return 'X-'
        else:
            return 'X+'
    else:#Horizontal Line
        if(debug): print('Horizontal Line')
        maxX=max(p1[0],p2[0])
        CellUp=(maxX-0.5,p1[1]+0.5)
        CellDown=(maxX-0.5,p2[1]-0.5)
        if(debug): print(CellUp,CellDown)
        if(point_in_polygon(CellUp,poly)):
            return 'Y+'
        else:
            return 'Y-'

def deriveFaultCells(CellSide,Edge,k):
    '''Derive all cell location along this edge

    k is the designed value of k

    Author:Bin Wang(binwang.0213@gmail.com)
    Date: Sep. 2018
    '''
    CellLocs=[]
    StartPos,EndPos=Edge[0],Edge[1]
    step=1
    offset=0

    if(CellSide=='X-' or CellSide=='X+'):#Horizontal Line
        NumEdgeNodes=int(StartPos[1]-EndPos[1])
        if(NumEdgeNodes>0):# j change from large val to low val
            step=-1
            NumEdgeNodes=abs(NumEdgeNodes)
            offset=1
        for j in range(StartPos[1],EndPos[1],step):
            if(CellSide=='X-'): CellLocs.append((StartPos[0]-1,j-offset,k))
            if(CellSide=='X+'): CellLocs.append((StartPos[0],j-offset,k))
    else:#Vertical Line
        NumEdgeNodes=int(StartPos[0]-EndPos[0])

        if(NumEdgeNodes>0):# i change from large val to low val
            step=-1
            NumEdgeNodes=abs(NumEdgeNodes)
            offset=1
        for i in range(StartPos[0],EndPos[0],step):
            if(CellSide=='Y-'): CellLocs.append((i-offset,StartPos[1]-1,k))
            if(CellSide=='Y+'): CellLocs.append((i-offset,StartPos[1],k))

    return CellLocs

def isFaultOnBoundaryEdge(GRDECL_Data,fault):
    '''Determine if a fault is a boundary edge

    Fault edge (0,5)->(0,15), the constant axis is 0 which is on the boundary

    Author:Bin Wang(binwang.0213@gmail.com)
    Date: Sep. 2018
    '''
    BoundaryEdge='InternalFault'
    NX,NY=GRDECL_Data.NX,GRDECL_Data.NY

    if(fault[0][0]==fault[1][0]): #Constant X coord
        X_const=fault[0][0]
        if(X_const==0):  BoundaryEdge='X-'
        if(X_const==NX): BoundaryEdge='X+'
    if(fault[0][1]==fault[1][1]): #Constant X coord
        y_const=fault[0][1]
        if(y_const==0):  BoundaryEdge='Y-'
        if(y_const==NY): BoundaryEdge='Y+'
    
    return BoundaryEdge


##------------------Geometry--------------------
def Shapely2List_MultiLineString(lines):
    LinesList=[]
    for line in lines:
        temp=list(line.coords)
        LinesList.append(tuple(temp))
    return LinesList

def Shapely2List_MultiPoint(points):
    PointList=[]
    for p in points:
        if(p.geom_type=='Point'):
            temp=list(p.coords)[0]
            PointList.append(tuple([int(temp[0]),int(temp[1])]))
    return PointList

def Shapely2List_MultiPolygon(polygons):
    PolygonList=[]
    for p in polygons:
        temp=list(p.exterior.coords)
        PolygonList.append(tuple(temp))
    return PolygonList

def calcDist(Pts0=(0,0),Pts1=(1,1)):
    '''Calculating distance of two points
    '''
    return np.sqrt((Pts1[0]-Pts0[0])**2+(Pts1[1]-Pts0[1])**2)

def point_in_polygon(pts,polygon):
    #Test of point is in a 2D polygon
    #https://stackoverflow.com/questions/21339448/how-to-get-list-of-points-inside-a-polygon-in-python
    #pts has to be numpy array with size of Nx2
    #polygon is list of tuple, eg. [(0,0), (0, 1), (1, 1), (1, 0),(0,0)]

    #Create polygon 
    p = Path(polygon)
    return p.contains_point(pts)

def points_in_polygon(pts,polygon,flag=1):
    #Test of point is in a 2D polygon
    #https://stackoverflow.com/questions/21339448/how-to-get-list-of-points-inside-a-polygon-in-python
    #pts has to be numpy array with size of Nx2
    #polygon is list of tuple, eg. [(0,0), (0, 1), (1, 1), (1, 0),(0,0)]

    #Create polygon 
    p = Path(polygon)
    pts=np.array(pts) # convert list of tuple into Nx2 array
    return p.contains_points(pts)*flag


def point_in_line(pts,A,B):
    #Test of point(pts) lies on line segment (AB)
    #https://stackoverflow.com/questions/328107/how-can-you-determine-a-point-is-between-two-other-points-on-a-line-segment?noredirect=1&lq=1
    
    epsilon=0.0000000001
    squaredlengthba=(A[0] - B[0])**2 + (A[1] - B[1])**2
    
    crossproduct = (pts[1] - A[1]) * (B[0] - A[0]) - (pts[0] - A[0]) * (B[1] - A[1])
    if abs(crossproduct) > epsilon : return False   # (or != 0 if using integers)

    dotproduct = (pts[0] - A[0]) * (B[0] - A[0]) + (pts[1] - A[1])*(B[1] - A[1])
    if dotproduct < 0 : return False

    squaredlengthba = (B[0] - A[0])*(B[0] - A[0]) + (B[1] - A[1])*(B[1] - A[1])
    if dotproduct > squaredlengthba: return False

    return True

def DrawPolygons(polygons):
    #https://stackoverflow.com/questions/32141476/how-to-fill-polygons-with-colors-based-on-a-variable-in-matplotlib
    
    font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 14,
                }

    fig,ax = plt.subplots(figsize=(6, 6), dpi=80, facecolor='w', edgecolor='k')

    patches=[]
    for p in polygons:
        patches.append(Patches.Polygon(np.array(p),True))

    p = PatchCollection(patches,cmap=matplotlib.cm.rainbow, alpha=0.8)
    p.set_edgecolor('k')
    colors = 10*np.random.random(len(patches))
    p.set_array(np.array(colors))
    ax.add_collection(p)
    fig.colorbar(p, ax=ax)

    plt.axis('equal')
    plt.gca().invert_xaxis()
    plt.title('Domain Decomposition Map (%d domains)' %(len(patches)), fontdict=font)
    plt.xlabel('X', fontdict=font)
    plt.ylabel('Y', fontdict=font)
    plt.show()


def DrawPath(lines,labels=[],endpoints=[]):
    #https://matplotlib.org/users/path_tutorial.html

    plt.figure(num=None, figsize=(10, 10), dpi=80, facecolor='w', edgecolor='k')
    font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 16,
                }

    for i,verts in enumerate(lines):
        #print(i,verts)
        plt.plot(*zip(*verts),label=labels[i])
        #plt.scatter(*zip(*verts),alpha=0.7,s=20)
    
    if(len(endpoints)>0):
        plt.scatter(*zip(*endpoints),facecolors='none', edgecolors='r',alpha=0.7,s=20)

    #ax.set_xlim(-2,2)
    #ax.set_ylim(-2,2)
    #plt.legend(bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.,fontsize=10)
    plt.axis('equal')
    plt.gca().invert_xaxis()
    plt.title('X-Y Plane Fault Map', fontdict=font)
    plt.xlabel('X', fontdict=font)
    plt.ylabel('Y', fontdict=font)
    plt.grid()
    plt.show()