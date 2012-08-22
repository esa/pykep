"""
This module contains several examples on the use of PyKEP, try one of the
run_example() methods!!!
"""

from PyKEP import __extensions__
from _ex2 import run_example2

if __extensions__['pygmo']:
	from _ex1 import run_example1
	from _ex3 import run_example3
	from _ex4 import run_example4
	from _ex5 import run_example5
