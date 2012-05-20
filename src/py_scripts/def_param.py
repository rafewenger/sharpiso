'''
HERE WE DEFINE THE DEFAULT PARAMETERS OF THE SCRIPT 
'''
fi=open ('testData/data_sheet.txt','r')
flist=[]
for files in fi:
  a=files.split()
  flist.append(a[0])

# IF YOU CHANGE THE  ORDER, THEN CHANGE THE NAMES ACCORDINGLY
'''
POS_DEF = ['gradEC','gradCD','gradNS','gradCDdup','gradCS','gradN',
              'gradIE','gradIES','gradC','gradNIE','gradNIES',
              'gradES','cube_center','centroid'] 
'''
POS_DEF = ['gradEC','gradCD','gradNS','gradCDdup','gradCS','gradN',
              'gradIE','gradIES','gradC','gradNIE','gradNIES',
              'gradES','cube_center','centroid']
POS_DEF_NAMES = ['EC','CD','NS','CDdup','CS','N',
              'IE','IES','C','NIE','NIES',
              'ES','cc','cntrd']                              

ISO_DEF = ['10.1','10.5','10.6']

data_loc='testData/'
temp_loc='tempFiles/' # all the temporary files[*.off and *.line are shifted to this folder]
DATA_DEF =flist[1:30]
tests=['default','-reposition','lindstrom','-allow_conflict','centroid_conflict'\
,'clamp_far','reselctg','reselectg and allow_conflict']
