#! /usr/bin/env python
# encoding: utf-8

import os, shutil, re

print "Install Integration - Ensta ParisTech"
print "Clean old version"
do=0
if os.path.isfile('src/integrate/ibex_integrate.h'):
    shutil.rmtree('src/integrate')
    do=1

src = 'tmp_install/integrate'
dst = 'src/integrate'
shutil.copytree(src, dst, symlinks=False)
    
shutil.copy('tmp_install/example_integrate.cpp','examples/example_integrate.cpp')
    
print 'You can now re-configure and re-compile Ibex'
  
  
  
  
  