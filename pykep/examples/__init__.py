"""
This module contains several examples on the use of pykep, try one of the
run_example() methods!!!
"""

from pykep import __extensions__
from ._ex2 import run_example2

if __extensions__['pygmo']:
    #from ._ex1 import run_example1
    from ._ex3 import run_example3
    from ._ex4 import run_example4
    from ._ex5 import run_example5
    from ._ex6 import run_example6
    from ._ex7 import run_example7
    from ._ex8 import run_example8
    from ._ex9 import run_example9
    from ._ex10 import run_example10
    from ._ex_utilities import add_gradient, algo_factory
