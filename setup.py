import os
from sys import platform as _platform

if _platform == "darwin":
    os.environ["CC"] = "gcc-5"
    os.environ["CXX"] = "g++-5"

try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize
import numpy

sourcefiles = ['src/gropt.pyx']

if _platform == "darwin":
    extensions = [Extension("gropt",
                    sourcefiles,
                    language="c++",
                    include_dirs=[".",  "./src", "/usr/local/include/", numpy.get_include()],
                    library_dirs=[".", "./src", "/usr/local/lib/"],
                    extra_compile_args=['-std=c++11'],
                   )]
elif _platform == "win32":
    extensions = [Extension("gropt",
                    sourcefiles,
                    language="c++",
                    include_dirs=[".", "./src", numpy.get_include()],
                    library_dirs=[".", "./src"]
                   )]

setup(
    ext_modules = cythonize(extensions)
)