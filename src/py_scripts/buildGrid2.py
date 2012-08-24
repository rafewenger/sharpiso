#!/usr/bin/python
import sys
import struct

# the version 2 of the buildgrid
# additional information
# take as input the largest dim of the model and
# the scale


# example run python buildGrid2.py temp/test-cube-1-dc-1-0.5.dcf longest_dim scale oct_tree_depth



positions=[[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]
edges= [[0,4],[1,5],[2,6],[3,7],[0,2],[1,3],[4,6],[5,7],[0,1],[2,3],[4,5],[6,7]]

locations=[]

counter=10
level = 0
size_of_grid = 0
scale = 0
width_of_grid = 0
longest_dim_of_model = 0
loc = [0,0,0]
fi=open ('vert.off','width_of_grid')
fi2=open('scalar.txt', 'w')

#MAIN
def main():
    l=[]
    with open(sys.argv[1], "rb") as f:
        byte = f.read(1)
        while byte != "":
            l.append(byte)
            byte = f.read(1)
    global counter
    print counter
    X = convert("i", l,   4)
    Y = convert("i", l,   4)
    Z = convert("i", l,   4)
    print X, Y,Z, counter
    #node
    node_type = convert ("i", l, 4)
    print 'node_type',node_type

    global longest_dim_of_model
    global scale
    global width_of_grid
    longest_dim_of_model = float(sys.argv[2])
    scale = float(sys.argv[3])
    octtree_depth = sys.argv[4]
    width_of_grid =  longest_dim_of_model / scale
    
    CreateEmptyGrid(width_of_grid, octtree_depth )
    readNode(l)




def readNode(l):
    global counter
    global level
    global loc
    global curr_width
    # initialize both the location to loc
    parent_location = loc
    level = level + 1
    global width_of_grid
    curr_width = width_of_grid*1.0/2**(level)
    print 'width_of_grid in readNode', width_of_grid
    print 'curr width',curr_width, 'Level[',level,']'

    #node
    for n in range (8):
        #print '----- N ',n,'key positions ',positions[n],
        child_location=[0,0,0]
        #print 'Parent location',parent_location,
        computeOriginOfTheChildLocation(curr_width, n, parent_location, child_location)
        #print' Child location ', child_location
        node_type=convert("i", l,  4)
        #print 'node_type',node_type

        if node_type==1:
            isempty=convert("h", l,  2)
            print 'isempty',isempty,

        if node_type==2:
            sign_corner=[]
            for i in range(8):
                sign_corner.append(convert("h", l,  2))
            #Store all the intersection points and then compute the centroid
            set_of_intersection_points=[]
            
            
            for edge_number in range(12):
                ne=convert("i",l,4)
                #print 'ne[',ne,']'
                intersection_data=[]
                offset=0
                # If there is an intersection, figure out the edge points and the intersection point on the edge
                if ne>0:
            #Compute the end points of the edge which is intersected
            #DEBUG also compute the actual point of intersection
                    endPts = ComputeEdgePoints(edge_number, child_location, curr_width)
                    print >>fi2, endPts[0] , ' ', 0
                    print >>fi2, endPts[1], ' ', curr_width
                    
                    #print 'edge Number', edge_number, 'endpoints', endPts
                    #display the endpoints
                    for i in range(2):
                        for j in range(3):
                            print >>fi,endPts[i][j],
                        print >>fi, ''
                    vecDir = [(endPts[1][x]-endPts[0][x]) for x in range(3)]
                    vecDirNorm=[x/curr_width for x in vecDir]
                    
                    
                    offset=convert("f", l, 4)
                    n0=convert("f", l, 4)
                    n1=convert("f", l, 4)
                    n2=convert("f", l, 4)
                    intersection_data.append(offset)
                    intersection_data.append(n0)
                    intersection_data.append(n1)
                    intersection_data.append(n2)
                    intersection_point = [endPts[0][x]+vecDirNorm[x]*curr_width*offset for x in range(3)]
                    
                    set_of_intersection_points.append(intersection_point)
                if ne >=2:
                    print 'More than 2 intersections'
                    sys.exit()
                    
                                        
            print 'intersection points ', set_of_intersection_points
            #CentroidSum=[0,0,0]
            #for k in range(len(set_of_intersection_points)):
            #    CentroidSum[0] = CentroidSum[0] + set_of_intersection_points[k][0]
            #    CentroidSum[1] = CentroidSum[1] + set_of_intersection_points[k][1]
            #    CentroidSum[2] = CentroidSum[2] + set_of_intersection_points[k][2]
            #print 'centroid sum', CentroidSum
            #print 'centroid ' ,[CentroidSum[x]/len(set_of_intersection_points) for x in range(3)] 
        if node_type==0:
            print 'calling readNode again'
            loc = temploc
            readNode(l)
            level = level-1
            loc = origloc



# find the end points of the edges
def ComputeEdgePoints(edge_number, child_location, curr_width):
    endPts=[]
    for j in range(2):
        local_coord = [curr_width*x for x in positions[edges[edge_number][j]] ]
        # DEBUG
        #print 'localcoord',local_coord
        TempEndPt=[0,0,0]
        for m in range (3):
            TempEndPt[m]=child_location[m]+local_coord[m]
        endPts.append(TempEndPt)
    return endPts



def convert(typ,l, num):
    global counter
    s=counter
    counter = counter + num
    return   struct.unpack("<"+typ, ''.join(l[s:counter]))[0]


#
def computeOriginOfTheChildLocation(curr_width, n, parent_location, child_location):
    pos = positions[n] # is a 3 tuple
    child_location[0]=parent_location[0]+pos[0]*curr_width
    child_location[1]=parent_location[1]+pos[1]*curr_width
    child_location[2]=parent_location[2]+pos[2]*curr_width

# create the empty scalar grid
def CreateEmptyGrid(width_of_grid, octtree_depth ):
    for i in range (2**octtree_depth+1):
    for j in range (2**octtree_depth+1):
        for k in range (2**octtree_depth+1):    
            tempLoc =[ i * width_of_grid/(2**octtree_depth), j * width_of_grid/(2**octtree_depth),
                       k * width_of_grid/(2**octtree_depth)]
            locations.append (tempLoc)

# call the main function
if __name__ == "__main__":
    main()
