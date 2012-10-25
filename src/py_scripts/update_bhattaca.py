#!/usr/bin/python 

import subprocess as sp 
import sys as s

sp.call(['sudo' ,'apt-get', 'update'])
sp.call(['sudo', 'apt-get', 'dist-upgrade'])
sp.call(['sudo', 'reboot'])
