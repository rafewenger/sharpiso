#!/usr/bin/python
'''
Merge Test 
'''
import subprocess as proc
import sys as sys
import os

#configurations

#setup isovalues
isoval_base = 3.0
isoval_offset = [0.0, 0.23,0.5,0.8,0.9]

#set up location of file locations
fread = open ('./file-names.txt','r') # contains the names of the files on which the test is run
fcountdegree= open ('./edge-count.txt','w') # stores the edge count values
ftestdetails = open ('./test-details.txt','w') # stores details of the test runs



############
## mergeSharp  runs 
## set up the test you want to run.
## write the test as 'mergesharp' and then append it to the mergesharpTests, as shown below
############
mergesharpTests=[]

######################################################################
#test1
mergesharp = ['./mergesharp', '-trimesh']
mergesharpTests.append(mergesharp)
'''
#test2
mergesharp = ['./mergesharp', '-no_merge_sharp','-position','gradEC','-trimesh' ,'-dist2centroid','-allow_conflict',\
            '-multi_isov','-max_dist','0','-clamp_far','-sharp_edgeI','-lindstrom_fast']
mergesharpTests.append(mergesharp)

#test3 
mergesharp = ['./mergesharp', '-no_merge_sharp','-position','gradEC','-trimesh' ,'-dist2centroid','-clamp_conflict',\
            '-multi_isov','-sharp_edgeI','-lindstrom_fast']
mergesharpTests.append(mergesharp)


#test4
mergesharp = ['./isodual3D', '-merge_sharp','-position','gradEC','-trimesh' ,'-dist2centroid','-allow_conflict',\
            '-multi_isov','-sharp_edgeI','-lindstrom_fast']
mergesharpTests.append(mergesharp)
'''
#######################################################################
                        
findsharp = ['findsharp','140']
countdegree =['countdegree', '-fshort']



'''
for each file run isodual3D
'''
def test(n,iso):    
    for f in fread:
        fname_with_nrrd = f.split("/")[len(f.split("/"))-1]
        fname = fname_with_nrrd.split(".")[0]
        fgradnoisename = fname+".grad.nrrd"
        for j in isoval_offset:
                #iso_temp = isodual3D[:]
                i = j+isoval_base;
                iso_temp = iso[:]        
                test_details = str(n)+","+fname+","+str(i)
                n = n+1
                #fname_line_name = loc+"output_offs/"+fname+"-iso-"+str(i)+".line"
                iso_temp.append('-o')
                iso_temp.append('out.off')
                iso_temp.append('-s')
                #iso_temp.append('-gradient')
                #iso_temp.append(fgradnoisename)
                iso_temp.append(str(i))
                iso_temp.append(f.strip())
                print >>ftestdetails,test_details
		
                procced = proc.check_call(iso_temp)
                
                if procced==0:
                    if procced == 0:
                        findsharpTemp = findsharp[:]
                        findsharpTemp.append("out.off")
                        procced = proc.check_call(findsharpTemp)
                        if procced ==0:
                            countdegreeTemp = countdegree[:]
                            countdegreeTemp.append("out.line")
                            
                            procced=proc.check_call(countdegreeTemp,stdout=fcountdegree)
                        else:
                            print 'error in degree count'
                            sys.exit()
                    else:
                        print 'error in meshconvert'
                        sys.exit()
                                
                else:
                    print "error"
                    sys.exit();
		
                
                
def main():
    print "status: start test"
    n=1
    for iso in mergesharpTests:
       # print iso
       print "status: in test\n ", iso
       print >>ftestdetails, iso  
       test (n,iso)
if __name__ == "__main__":
    main()            
            
            
            
    
