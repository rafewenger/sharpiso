'''
HERE WE DEFINE THE DEFAULT PARAMETERS OF THE SCRIPT 
'''
fi=open ('testData/data_sheet.txt','r')
flist=[]
for files in fi:
  a=files.split()
  flist.append(a[0])

POS_DEF = ['gradEC','gradCD','gradNS','gradCDdup','gradCS','gradN',
              'gradIE','gradIES','gradC','gradNIE','gradNIES',
              'gradES','cube_center','centroid']              

ISO_DEF = ['4.1']
data_loc='testData/'
temp_loc='tempFiles/' # all the temporary files[*.off and *.line are shifted to this folder]
DATA_DEF =['annulus1.nrrd', 'flange2.nrrd' ]
DATA_DEF =flist[:]
 

POS_DEF_LONG = ['gradEC','gradNS']
ISO_DEF_LONG = ['4.1', '4.2']
data_loc='testData/'
temp_loc='tempFiles/' # all the temporary files[*.off and *.line are shifted to this folder]
DATA_DEF_LONG =a[:] 


