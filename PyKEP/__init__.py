import sims_flanagan
import examples

__doc__ = 'PyKEP is the answer'
__all__ = ['core', 'sims_flanagan', 'examples']

# For convenience, bring all core classes into the root namespace when importing *.
from core import *
__all__ += filter(lambda name: not name.startswith('_'),dir(core))
