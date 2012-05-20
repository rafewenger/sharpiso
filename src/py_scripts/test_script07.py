#run the test_script06 over a bunch of test data sets
import subprocess as sp
args=['python','test_script06.py']
fi=open ('compile_annulus.txt','w')
print 'default parameters run...'
p1 = sp.call (args)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Lindstrom----'
print 'lindstrom parameters run...'
temp.append ('-lindstrom')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Reposition----'
print '----Reposition----'
temp.append ('-reposition')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----allow_conflict----'
print '----allow_conflict----'
temp.append ('-allow_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----clamp_conflict----'
print '----clamp_conflict----'
temp.append ('-clamp_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----reselectg----'
print '----reselectg----'
temp.append ('-reselectg')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----clamp_far----'
print '----clamp_far----'
temp.append ('-clamp_far')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,	

#--------------with central gradients-----------------
print >>fi, '----with cgradients----'
print '----with cgradients----'
args=['python','test_script06.py','-cgradient']
p1 = sp.call (args)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Lindstrom with cgradient----'
print '----Lindstrom with cgradient----'
temp.append ('-lindstrom')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Reposition with cgradient----'
print '----Reposition with cgradient----'
temp.append ('-reposition')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----allow_conflict with cgradient----'
print '----allow_conflict with cgradient----'
temp.append ('-allow_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----clamp_conflict with cgradient----'
print '----clamp_conflict with cgradient----'
temp.append ('-clamp_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----reselectg----'
print '----reselectg----'
temp.append ('-reselectg')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
	
temp=args[:];
print >>fi,'----clamp_far with cgradient----'
temp.append ('-clamp_far')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line, 

	
