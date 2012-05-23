'''
HERE WE DEFINE THE DEFAULT PARAMETERS OF THE SCRIPT 
'''


# IF YOU CHANGE THE  ORDER, THEN CHANGE THE NAMES ACCORDINGLY

POS_DEF = ['gradEC','gradCD','gradNS','gradCDdup','gradCS','gradN',
              'gradIE','gradIEDir','gradIES','gradC','gradNIE','gradNIES',
              'gradES','cube_center','centroid',]
POS_DEF_NAMES = ['EC','CD','NS','CDdup','CS','N',
              'IE','IES','C','NIE','NIES',
              'ES','cc','cntrd']                              

ISO_DEF = ['15.5','15.8']

data_loc='testData/'
temp_loc='tempFiles/' # all the temporary files[*.off and *.line are shifted to this folder]
DATA_DEF =[]


TEST_HEADER_STRING1=\
'fname    dir    iso    clmp   repo  lin   allw    cntrd   clmp  rslctg  alwcnf clmpcnf'
TEST_HEADER_STRING2=\
'                       cnflct             cnflct  cnflct   far          repos   repos  '
