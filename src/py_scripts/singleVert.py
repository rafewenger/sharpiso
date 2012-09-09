'''
This will check if there are single vertex per cube
'''
import sys
import numpy


depth = 1
cubePtsList =[]
cube=numpy.zeros((2**depth,2**depth,2**depth))
unitLength =1

class  CubePtsType:
	numOfPoints = 0
	pointsLocation=[]
	
def main ():
	print 'this is main '
	readOff()
	#test()


def test():
	print 'in test'
	global depth
	global cubePtsList
	for i in range (2**depth):
		for j in range (2**depth):
			for k in range (2**depth):
				temp = CubePtsType()
				temp.numOfPoints = 2;
				a= [i,j,k]
				temp.pointsLocation = a
				cubePtsList.append(temp)
	
	for i in range (len (cubePtsList)):
		temp = CubePtsType()
		temp = cubePtsList[i]
		print 'n = ', temp.numOfPoints,'numv', i, ' loc ', temp.pointsLocation
		
	for i in range (2**depth):
		for j in range (2**depth):
			for k in range (2**depth):
				print  i,' ',j,' ',k,' --  ', (2**depth)*(2**depth)*i + 2**depth*j + k
	



def setList():
	global depth
	global cube
	global unitLength
	global cubePtsList
	
	for i in range (2**depth):
		for j in range (2**depth ):
			for k in range (2**depth):
				temp=CubePtsType()
				temp.numOfPoints=0
				temp.pointsLocation =[]
				cubePtsList.append(temp)

'''
Read from a off file 
'''
def readOff():
	global depth
	global cube
	global unitLength
	global cubePtsList
	global scale 
	fileName = sys.argv[1]
	depth = int(sys.argv[2])
	scale = float(sys.argv[3])
	unitLength =(1.0/scale)/2**depth
	
	print 'File name is ', fileName
	f = open (fileName, 'r')
	data = f.readlines()
	metadata=data[1].split(" ")
	numv= int(metadata[0])
	print 'numv',numv
	
	# set up the List 
	setList()
	
	for i in range (numv):
		lineNo = 2+i
		vertLoc=data[lineNo].split(" ")	
		# for each Pt 
		for m in range(2**depth):
			for n in range (2**depth):
				for p in range (2**depth):
					xMin = m*unitLength
					xMax = m*unitLength + unitLength
					
					yMin = n*unitLength
					yMax = n*unitLength + unitLength
					
					zMin = p*unitLength
					zMax = p*unitLength + unitLength
					if float(vertLoc[0]) > xMin and float(vertLoc[0]) < xMax:
						if float(vertLoc[1]) >yMin and float(vertLoc[1])<yMax:
							if float(vertLoc[2]) >zMin and float(vertLoc[2])<zMax:
								ind = (2**depth)*(2**depth)*m + 2**depth*n + p
								#increase the number of points 
								cubePtsList[ind].numOfPoints = cubePtsList[ind].numOfPoints+1
								#add the corresponding locations
								cubePtsList[ind].pointsLocation.append([float(vertLoc[0]),float (vertLoc[1]), float (vertLoc[2])])
	print 'length of List    ', len(cubePtsList)
	for i in range (len(cubePtsList)):
		if cubePtsList[i].numOfPoints  > 1 :
			print ' num points ',' -- ',cubePtsList[i].numOfPoints ,
			for j in range (len((cubePtsList[i].pointsLocation))):
				print cubePtsList[i].pointsLocation[j],
				if j%3 ==0 :
					print ''
			
				

		
			
	
	
	
	
	
# call the main function
if __name__ == "__main__":
    main()