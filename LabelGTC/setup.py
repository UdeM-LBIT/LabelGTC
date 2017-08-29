from distutils.core import setup, Extension
from Cython.Build import cythonize

import fnmatch
import os
import sys
sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), "src/SuperGeneTrees")))


libraries = []
for root, dirnames, filenames in os.walk('src/SuperGeneTrees'):
    for filename in fnmatch.filter(filenames, '*.cpp'):
        if 'main.' not in filename:
            libraries.append(os.path.join(root, filename))

setup(
    name = "labelgtc",
    author = "Kylian Berrebi",
    version = "1.0.0rc",
    classifiers= [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: Education',
        ],
    packages=['lib', 'lib.TreeLib', 'lib.PolyRes', 'lib.SGT', 'lib.LabelGTC'],
    py_modules=[],
    #  scripts=['bin/labelgtc'], # a labelgtc script should be defined
    install_requires=['ete3', 'numpy', 'cython'],
    ext_modules=cythonize(Extension(
    "lib.SGT.minSGT", 
    sources=["src/minSGT.pyx"]+libraries,
    language="c++",
    extra_compile_args=["-std=c++0x"],
    extra_link_args=["-std=c++0x"]
    )))