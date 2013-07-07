#!/usr/bin/python
'''
Code to generate data
OUTPUT
* data-set-info.txt
* file-names.txt
'''
import subprocess as proc

#standard arguments to ijkgenscalar2
standard_args = ['ijkgenscalar2', '-grad', '-dim', '3', '-asize', '100']
cgradient = ['cgradient']


center_opts = ['50.0 50.0 50.0','50.6 53.1 49.7']

dir_opts = ['1 0 0','1 1 0','1 1 1', '2 2 5 ','2 5 7','1 1 4','1 1 7' ]
dir_names = ['100','110','111', '225','257','114','117']

side_dir_opts = ['1 0 0', '1 1 0', '1 1 1']
side_dir_names = ['100','110','111']


datainfo = open ('./data-set-info.txt','w')
fname = open ('./file-names.txt','w')



type = 'cube'
n=0
for cen in center_opts:
    temp_args = standard_args[:]
    temp_args.append('-field')
    temp_args.append(type)
    temp_args.append('-n')
    temp_args.append('3')
    temp_args.append('-stack')
    temp_args.append('-center')
    temp_args.append(cen)
    dirNum=0
    for dir in dir_opts:
        temp_args_center_fixed = temp_args[:]
        temp_args_center_fixed.append('-dir')
        temp_args_center_fixed.append(dir)
        temp_args_center_fixed.append('-translate')
	temp_args_center_fixed.append('7 8 9')  
        side_dir=0
        for sid in side_dir_opts:
            if dir != sid :
                n=n+1
                temp_args_dir_fixed = temp_args_center_fixed[:]
                temp_args_dir_fixed.append("-side_dir")
                temp_args_dir_fixed.append(sid)
                fname_nrrd_name = type+ '-dir-' +dir_names[dirNum]+"-sd-"+side_dir_names[side_dir]+"-"+str(n)+".nrrd"
                temp_args_dir_fixed.append(fname_nrrd_name)
                print >> datainfo, fname_nrrd_name, " center ",cen, ' dir ', dir
                print >>fname, fname_nrrd_name
		print temp_args_dir_fixed
                proc.call(temp_args_dir_fixed)
                
                side_dir = side_dir + 1
                print " create data number  " , n , "done"
        #print 'return code ', proc1
        dirNum=dirNum+1 
       
print "ALL TESTS RAN "
