'''
Convert a line file into a vect file for geomview 
Arindam bhattacharya
'''
import sys
def main ():
    if (len(sys.argv) != 2):
        print 'enter only one argument'
    else:
        f=open (sys.argv[1],'r')
        points=[]
        edges=[]
        num_points=0
        num_edges=0
        for index,line in enumerate(f):
            if index==0:
                if line.split()[0] != 'LINEC':
                    print 'not a line file '
            elif index==2:
                vals=line.split()
                num_points=int(vals[0])
                num_edges=int(vals[1])
            else:
                if  index>=3 and index < 3+int(num_points):
                    temp = line.split()
                    points=points+temp[:]
                if index > 3+int(num_points):
                    temp = line.split()
                    edges=edges+temp[:]
        setup_vect(points,edges, num_points, num_edges)
            
def setup_vect(points, edges, num_points, num_edges):
    
    outname=sys.argv[1].split('.')[0]+'.vect'
    print 'printing the',str(outname),'file'
    f=open(outname,'w')
    print >>f,'VECT'
    print >>f,num_edges, 2*num_edges, 1
    li=[str(2)]*num_edges
    print >>f,(' ').join(li)
    li=[str(0)]*(num_edges-1)
    print >>f,'1',(' ').join(li),'\n'
    for i in range(num_edges):
        pt1=int(edges[2*i])
        pt2=int(edges[2*i+1]) 
        print >>f, points[3*pt1],points[3*pt1+1],points[3*pt1+2],' ',
        print >>f, points[3*pt2],points[3*pt2+1],points[3*pt2+2]
    print >>f,'\n1 0 0 1'
    print 'DONE'
    
if __name__ == "__main__":
    main()
    
