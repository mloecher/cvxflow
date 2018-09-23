# CVXFlow

This repository contains sample code to use for generating convex optimized waveforms for 4D-Flow MRI

## Installation

All C++ code is in the src/ directory.  The optimization source file can be compiled directly, or a module for use with python can be compiled with the included setup.py.

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) needs to be in your path somewhere.

gropt_proc.cpp and gropt_proc_spoiler.cpp can be compiled to standalone programs that were used to generate waveforms on the scanner.

## Python

Compile the module with 'python setup.py build_ext --inplace'

An example python notebook is included.
