from distutils.core import setup, Extension
from Cython.Build import cythonize

import fnmatch
import os
import sys
sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), "SuperGeneTrees")))


libraries = []
for root, dirnames, filenames in os.walk('SuperGeneTrees'):
    for filename in fnmatch.filter(filenames, '*.cpp'):
    	if 'main.' not in filename:
        	libraries.append(os.path.join(root, filename))

setup(ext_modules=cythonize(Extension(
	"minSGT", 
	sources=["minSGT.pyx"]+libraries,
	language="c++",
	extra_compile_args=["-std=c++0x"],
    extra_link_args=["-std=c++0x"]
	)))