#!/usr/bin/python
#test 2 versions of religrad for differences
import subprocess
import sys
from subprocess import Popen, PIPE



testFiles=['data/cube.A.nrrd','data/cube.B.nrrd','data/annulus.dir111.A.nrrd','data/flange.dir321.A.nrrd','data/twocubes.dir111.A.nrrd']

fnull = open('/dev/null','w') # null file
def main():
  #print str(sys.argv)
  if (len(sys.argv) > 2):
    religradNew=sys.argv[1]
    religradOld=sys.argv[2]
    
    for fname in testFiles:
      print fname
      p=subprocess.check_call([religradNew, "-angle_based","-reliable_scalar_pred_dist","2", fname, "outNew.nrrd"],stderr=fnull,stdout=fnull)
      p=subprocess.check_call([religradOld, "-angle_based", "-reliable_scalar_pred_dist","2",fname, "outOld.nrrd"],stderr=fnull,stdout=fnull)
      p = Popen(["diff","outNew.nrrd","outOld.nrrd"], stdout=PIPE)
      output = p.communicate()[0]
      if (p.returncode):
        print "differnces detected"
  else:
    print "not enough arguments"
    

if __name__ == "__main__":
    main()  
