#run the test_script06 over a bunch of test data sets
import subprocess as sp
args=['python','test_script06.py']
fi=open ('compile.txt','w')
p1 = sp.call (args)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Lindstrom----'
temp.append ('-lindstrom')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Reposition----'
temp.append ('-reposition')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----allow_conflict----'
temp.append ('-allow_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----clamp_conflict----'
temp.append ('-clamp_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----clamp_far----'
temp.append ('-clamp_far')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,	

#--------------with central gradients-----------------

args=['python','test_script06.py','-cgradient']
p1 = sp.call (args)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Lindstrom with cgradient----'
temp.append ('-lindstrom')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,

temp=args[:];
print >>fi,'----Reposition with cgradient----'
temp.append ('-reposition')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----allow_conflict with cgradient----'
temp.append ('-allow_conflict')
p1 = sp.call (temp)
#read the output.txt
f1 = open ('output.txt','r')
for line in f1:
	print >>fi,line,
	
temp=args[:];
print >>fi,'----clamp_conflict with cgradient----'
temp.append ('-clamp_conflict')
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

	
